/* メインアプリ (classic script)
 * WebGL描画 + ノードUI + 入出力 + レシピ + 実在モード + 書き出し
 */
(function () {
  'use strict';
  var LOOP = window.CAUSTICS_SHADER.LOOP_SECONDS; // 8
  var $ = function (id) { return document.getElementById(id); };

  var state = window.CausticsRecipe.defaultState();
  var params = null; // resolveEnvironment の結果キャッシュ
  // FIX-2: 水の状態(慣性つきスペクトル)。refresh() が目標値を作り、
  // render() が毎フレーム追従させる。プリセット切替・レシピ読込では snap。
  var waterState = window.CausticsEnv.createWaterState(null);
  var waterTarget = { specL: 1.0, specAmp: 0.5, dirX: 0, dirY: -1, flowX: 0, flowY: 0, fine: 0.5 };
  var lastT = -1;
  function waterTargetFromParams(p) {
    return {
      specL: p.specL, specAmp: p.specAmp,
      dirX: p.waveDir[0], dirY: p.waveDir[1],
      flowX: p.flowPerSec[0], flowY: p.flowPerSec[1],
      fine: p.fine
    };
  }

  // ---------- WebGL ----------
  var canvas = $('gl');
  var gl = canvas.getContext('webgl', { preserveDrawingBuffer: true, antialias: false });
  var glOK = !!gl;
  var prog = null, uni = {}, tex = null;
  var quadBuf = null, quadLoc = -1;
  // FIX-7: 前方蓄積パス(スプラット描画+コースティクスバッファ)。
  // FIX-13: コースティクスバッファは描画解像度と同等(上限 1024)。型は半浮動小数点→
  // 浮動小数点→8bit加算近似の順に試し、すべて不可なら後方ヤコビアンに
  // フォールバックする(uUseForward=0)。水の状態の生成は現行のまま。
  var splatProg = null, splatUni = {}, splatLoc = -1, gridBuf = null, idxBuf = null, idxCount = 0, gridN = 256;
  var caus = { fbo: null, tex: null, w: 0, h: 0, mode: 'none', filter: 'nearest' }; // mode: half|float|byte|none, filter: linear|nearest(FIX-12)
  // FIX-15: 性能パス(適応品質ラダー)。変えるのは描画量(コースティクスバッファの
  // 解像度・スプラットグリッドの細かさ)のみ。数式・較正値・定数・プリセットには触れない。
  // ラダー: high(バッファ倍率1.0・グリッド256)/mid(0.75・192)/low(0.5・160)。
  // バッファ倍率は initCausticsTarget 呼び出し寸法への乗算(上限1024は維持)。
  var FIX15_QUALITY_ORDER = ['high', 'mid', 'low'];
  var FIX15_QUALITY_LADDER = {
    high: { buffer: 1.0, grid: 256 },
    mid: { buffer: 0.75, grid: 192 },
    low: { buffer: 0.5, grid: 160 }
  };
  var quality = 'high'; // 初期段は high
  var qualityMode = 'auto'; // 'auto' | 'high' | 'mid' | 'low'(手動選択値)
  var autoQuality = true; // 既定で自動制御ON。手動選択で止まり、「自動」に戻すと再開。
  var frameTimes16 = []; // FIX-16: 直近2秒のフレーム描画時刻(秒・昇順)。fpsは窓内数/実効秒。
  var windowFps16 = 0; // FIX-16: 窓ベースFPSの内部状態(0.5秒ごとに再計算)
  var lastFpsUpdateSec = -1e9; // FIX-16: fps値の最終更新時刻(秒)
  var lastOverlaySec16 = -1e9; // FIX-16: オーバーレイ表示の最終更新時刻(秒)
  var lastQualityEvalSec = -1e9; // 自動制御の最終評価時刻(秒)
  var lastQualitySwitchSec = -1e9; // 品質段の最終切替時刻(秒)
  var upStreak = 0; // FPS>57が続いた評価回数
  // FIX-15: 自動制御の純粋判定(cur: 現在段、fps: 平均FPS、streak: 上振れ連続回数、
  // sinceSwitchSec: 前回切替からの経過秒)。verify-node.js が取出し評価できるよう
  // 副作用なしにする。閾値: 1秒毎評価・45未満で1段下げ・57超が3回続けば1段上げ・
  // 切替は最短2秒間隔・下限low/上限high。
  function fix15Decide(cur, fps, streak, sinceSwitchSec) {
    var idx = FIX15_QUALITY_ORDER.indexOf(cur);
    if (idx < 0) idx = 0;
    if (fps < 45) {
      var down = idx < FIX15_QUALITY_ORDER.length - 1 ? FIX15_QUALITY_ORDER[idx + 1] : cur;
      return { next: down, streakNext: 0, doSwitch: down !== cur && sinceSwitchSec >= 2 };
    }
    if (fps > 57) {
      var s = streak + 1;
      if (s >= 3 && idx > 0 && sinceSwitchSec >= 2) {
        return { next: FIX15_QUALITY_ORDER[idx - 1], streakNext: 0, doSwitch: true };
      }
      return { next: cur, streakNext: s, doSwitch: false };
    }
    return { next: cur, streakNext: 0, doSwitch: false };
  }
  // FIX-16: 窓ベースFPS算出の純粋関数(frameTimes: 昇順の描画時刻配列(秒)、
  // nowSec: 現在時刻(秒))。verify-node.js が取出し評価できるよう副作用なし。
  // 直近2秒より古い時刻は窓外として数えない。fps=窓内フレーム数/窓の実効秒と
  // し、実効秒=min(2, now-最古時刻)とする(起動直後の過小評価を避けるため)。
  // 標本が2未満では0を返す(起動直後の不定値を出さない)。
  function fix16WindowFps(frameTimes, nowSec) {
    var n = 0, oldest = 0;
    for (var i = 0; i < frameTimes.length; i++) {
      if (frameTimes[i] > nowSec - 2) {
        if (n === 0) oldest = frameTimes[i];
        n++;
      }
    }
    if (n < 2) return 0;
    var span = nowSec - oldest;
    if (!(span > 1e-6)) return 0;
    if (span > 2) span = 2;
    return n / span;
  }
  // FIX-15: 現在の品質段に応じたコースティクスバッファ寸法(上限1024は維持)。
  function qualityBufferDims() {
    var q = FIX15_QUALITY_LADDER[quality] || FIX15_QUALITY_LADDER.high;
    return {
      w: Math.min(1024, Math.max(2, Math.round(canvas.width * q.buffer))),
      h: Math.min(1024, Math.max(2, Math.round(canvas.height * q.buffer)))
    };
  }
  // FIX-15: 品質段の切替。バッファ再生成とグリッド再構築を既存経路で行う。
  // canvas自体の寸法・プリセット・レシピ・書き出し(PNG/WebM)の結果寸法は変えない。
  function applyQuality(next, reason) {
    if (!FIX15_QUALITY_LADDER[next] || next === quality) return false;
    quality = next;
    var d = qualityBufferDims();
    initCausticsTarget(d.w, d.h);
    if (splatProg) buildSplatGrid(FIX15_QUALITY_LADDER[next].grid);
    if (window.__caustics) {
      window.__caustics.quality = quality;
      window.__caustics.caustics = { mode: caus.mode, w: caus.w, h: caus.h, grid: gridN, filter: caus.filter };
    }
    void reason;
    return true;
  }
  var mediaAspect = 16 / 9;
  var mediaKind = 'none'; // 'image' | 'video' | 'sample'
  var imgEl = null;
  var vidEl = $('vid');
  var sampleCanvas = null;
  var frameCount = 0;
  var shaderLog = '';

  function compileShader(type, src) {
    var s = gl.createShader(type);
    gl.shaderSource(s, src);
    gl.compileShader(s);
    var ok = gl.getShaderParameter(s, gl.COMPILE_STATUS);
    var log = gl.getShaderInfoLog(s) || '';
    shaderLog += (type === gl.VERTEX_SHADER ? '[vertex]' : '[fragment]') + (ok ? 'OK ' : 'FAIL ') + log + '\n';
    if (!ok) { gl.deleteShader(s); return null; }
    return s;
  }

  function initGL() {
    if (!gl) {
      setStatus('このブラウザではWebGLが使えません。');
      window.__caustics = { shaderOK: false, log: 'no webgl', frameCount: 0, fps: 0, quality: 'high' };
      return false;
    }
    try { gl.getExtension('OES_standard_derivatives'); } catch (e) { /* fwidth fadeはshader側でフォールバック */ }
    var vs = compileShader(gl.VERTEX_SHADER, window.CAUSTICS_SHADER.VERT_SRC);
    var fs = compileShader(gl.FRAGMENT_SHADER, window.CAUSTICS_SHADER.FRAG_SRC);
    if (!vs || !fs) {
      setStatus('シェーダーのコンパイルに失敗しました。');
      window.__caustics = { shaderOK: false, log: shaderLog, frameCount: 0, fps: 0, quality: 'high' };
      return false;
    }
    prog = gl.createProgram();
    gl.attachShader(prog, vs);
    gl.attachShader(prog, fs);
    gl.linkProgram(prog);
    if (!gl.getProgramParameter(prog, gl.LINK_STATUS)) {
      shaderLog += '[link]FAIL ' + (gl.getProgramInfoLog(prog) || '') + '\n';
      window.__caustics = { shaderOK: false, log: shaderLog, frameCount: 0, fps: 0, quality: 'high' };
      return false;
    }
    gl.useProgram(prog);
    quadBuf = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, quadBuf);
    gl.bufferData(gl.ARRAY_BUFFER, new Float32Array([-1, -1, 3, -1, -1, 3]), gl.STATIC_DRAW);
    quadLoc = gl.getAttribLocation(prog, 'aPos');
    gl.enableVertexAttribArray(quadLoc);
    gl.vertexAttribPointer(quadLoc, 2, gl.FLOAT, false, 0, 0);
    var names = ['uTex', 'uMediaAspect', 'uRes', 'uTime', 'uSunDir', 'uSunColor',
      'uSunI', 'uDiffuse', 'uSpecL', 'uSpecAmp', 'uWindDir', 'uFine', 'uFlowPerSec',
      'uDepthM', 'uRefract', 'uChroma', 'uCausticsI', 'uGlitterI', 'uAbsorb', 'uScatter',
      'uDepthShade', 'uTurbid', 'uExposure', 'uKc', 'uFocusM', 'uDgCap',
      'uCaustics', 'uCausticsAvg', 'uUseForward', 'uCausticsTexel'];
    names.forEach(function (n) { uni[n] = gl.getUniformLocation(prog, n); });
    tex = gl.createTexture();
    gl.bindTexture(gl.TEXTURE_2D, tex);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
    // FIX-7: 前方パスのスプラットプログラムをコンパイル。失敗時は後方フォールバックに倒す。
    try {
      var svs = compileShader(gl.VERTEX_SHADER, window.CAUSTICS_SHADER.SPLAT_VERT_SRC);
      var sfs = compileShader(gl.FRAGMENT_SHADER, window.CAUSTICS_SHADER.SPLAT_FRAG_SRC);
      if (svs && sfs) {
        splatProg = gl.createProgram();
        gl.attachShader(splatProg, svs);
        gl.attachShader(splatProg, sfs);
        gl.linkProgram(splatProg);
        if (!gl.getProgramParameter(splatProg, gl.LINK_STATUS)) {
          shaderLog += '[splat-link]FAIL ' + (gl.getProgramInfoLog(splatProg) || '') + '\n';
          splatProg = null;
        }
      }
    } catch (eSplat) { shaderLog += '[splat]FAIL ' + (eSplat && eSplat.message ? eSplat.message : eSplat) + '\n'; splatProg = null; }
    if (splatProg) {
      var snames = ['uRes', 'uTime', 'uSunDir', 'uSpecL', 'uSpecAmp', 'uWindDir',
        'uFine', 'uFlowPerSec', 'uDepthM', 'uRefract', 'uDgCap',
        'uAspect', 'uMeshN', 'uDispK', 'uDelta2Fwd',
        'uDispGL', 'uDispGM', 'uDispGH'];
      snames.forEach(function (n) { splatUni[n] = gl.getUniformLocation(splatProg, n); });
      splatLoc = gl.getAttribLocation(splatProg, 'aGrid');
      if (splatLoc >= 0) gl.enableVertexAttribArray(splatLoc);
      try {
        gridN = window.CausticsEnv.FIX9_MESH_GRID ? window.CausticsEnv.FIX9_MESH_GRID : 256;
      } catch (eG) { gridN = 256; }
      buildSplatGrid(gridN);
    }
    setStatus('シェーダー準備OK。画像・動画を開くか「サンプルで試す」を押してください。');
    window.__caustics = { shaderOK: true, log: shaderLog, frameCount: 0, w: 0, h: 0, fps: 0, quality: quality };
    return true;
  }

  function uploadSource() {
    if (!gl || !prog) return;
    gl.bindTexture(gl.TEXTURE_2D, tex);
    try {
      if (mediaKind === 'video' && vidEl.videoWidth) {
        gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, gl.RGBA, gl.UNSIGNED_BYTE, vidEl);
      } else if (mediaKind === 'image' && imgEl && imgEl.naturalWidth) {
        gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, gl.RGBA, gl.UNSIGNED_BYTE, imgEl);
      } else if (mediaKind === 'sample' && sampleCanvas) {
        gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, gl.RGBA, gl.UNSIGNED_BYTE, sampleCanvas);
      } else {
        return;
      }
    } catch (e) { setStatus('テクスチャの転送に失敗: ' + e.message); }
  }

  function sizeCanvasToMedia() {
    var w = 640, h = 400;
    if (mediaKind === 'video' && vidEl.videoWidth) { w = vidEl.videoWidth; h = vidEl.videoHeight; }
    else if (mediaKind === 'image' && imgEl && imgEl.naturalWidth) { w = imgEl.naturalWidth; h = imgEl.naturalHeight; }
    else if (mediaKind === 'sample' && sampleCanvas) { w = sampleCanvas.width; h = sampleCanvas.height; }
    var scale = Math.min(1, 1280 / w);
    canvas.width = Math.max(2, Math.round(w * scale));
    canvas.height = Math.max(2, Math.round(h * scale));
    mediaAspect = w / h;
    if (window.__caustics) { window.__caustics.w = canvas.width; window.__caustics.h = canvas.height; }
    // FIX-7: コースティクスバッファは描画解像度と同等(上限 1024×1024 程度)
    // FIX-8: 上限 512→1024(描画解像度と同等まで可)。フットプリント半径の下限式は不変。
    // FIX-13: バッファ全解像度化。>>1 をやめ canvas.width/height をそのまま使う(上限1024は維持)。
    // FIX-15: 品質段のバッファ倍率を乗算(上限1024は維持)。グリッドは現在の段に合わせる。
    var qd15 = qualityBufferDims();
    initCausticsTarget(qd15.w, qd15.h);
    var qg15 = (FIX15_QUALITY_LADDER[quality] || FIX15_QUALITY_LADDER.high).grid;
    if (splatProg && gridN !== qg15) buildSplatGrid(qg15);
    if (window.__caustics) {
      window.__caustics.caustics = { mode: caus.mode, w: caus.w, h: caus.h, grid: gridN, filter: caus.filter };
    }
  }

  // FIX-9: メッシュ用の水面グリッド頂点配列(aGrid: 0..1)と
  // インデックスバッファ(セル2三角形)を生成。256×256頂点=65536で
  // 最大index 65535のため Uint16 で足りる。
  function buildSplatGrid(n) {
    gridN = n;
    var arr = new Float32Array(n * n * 2);
    var k = 0;
    for (var iy = 0; iy < n; iy++) {
      for (var ix = 0; ix < n; ix++) {
        arr[k++] = ix / Math.max(n - 1, 1);
        arr[k++] = iy / Math.max(n - 1, 1);
      }
    }
    if (!gridBuf) gridBuf = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, gridBuf);
    gl.bufferData(gl.ARRAY_BUFFER, arr, gl.STATIC_DRAW);
    var quads = (n - 1) * (n - 1);
    var idx = new Uint16Array(quads * 6);
    var q = 0;
    for (var jy = 0; jy < n - 1; jy++) {
      for (var jx = 0; jx < n - 1; jx++) {
        var a = jy * n + jx, b = a + 1, c = a + n + 1, d = a + n;
        idx[q++] = a; idx[q++] = b; idx[q++] = c;
        idx[q++] = a; idx[q++] = c; idx[q++] = d;
      }
    }
    if (!idxBuf) idxBuf = gl.createBuffer();
    gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, idxBuf);
    gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, idx, gl.STATIC_DRAW);
    idxCount = quads * 6;
    // 頂点バッファを切り替えたので合成パスのポインタを戻す
    if (quadBuf && quadLoc >= 0 && prog) {
      gl.useProgram(prog);
      gl.bindBuffer(gl.ARRAY_BUFFER, quadBuf);
      gl.vertexAttribPointer(quadLoc, 2, gl.FLOAT, false, 0, 0);
    }
  }

  // FIX-7: コースティクスバッファ(浮動小数/半浮動小数、不可なら8bit加算にフォールバック)を確保。
  // いずれも不可なら caus.mode='none' とし、合成パスは後方ヤコビアンにフォールバックする。
  function initCausticsTarget(bw, bh) {
    if (!gl || !prog) return;
    caus.w = bw; caus.h = bh; caus.mode = 'none'; caus.filter = 'nearest';
    if (caus.fbo) { gl.deleteFramebuffer(caus.fbo); caus.fbo = null; }
    if (caus.tex) { gl.deleteTexture(caus.tex); caus.tex = null; }
    var halfExt = null, floatExt = null;
    try { halfExt = gl.getExtension('OES_texture_half_float'); } catch (eH) { halfExt = null; }
    try { floatExt = gl.getExtension('OES_texture_float'); } catch (eF) { floatExt = null; }
    // FIX-12: タイプ別フィルタ選択。half/float の LINEAR サンプリングには
    // 各別のリニア拡張が必要(OES_texture_half_float_linear /
    // OES_texture_float_linear)。無い環境で LINEAR を指定するとテクスチャが
    // 不完全扱いになり texture2D が (0,0,0,1) を返すため NEAREST に落とす。
    // byte(UNSIGNED_BYTE)は常に LINEAR。候補順位(half→float→byte)と
    // FBO 完成チェックは現行どおり維持。
    var halfLinearExt = null, floatLinearExt = null;
    try { halfLinearExt = gl.getExtension('OES_texture_half_float_linear'); } catch (eHL) { halfLinearExt = null; }
    try { floatLinearExt = gl.getExtension('OES_texture_float_linear'); } catch (eFL) { floatLinearExt = null; }
    var cands = [];
    if (halfExt) cands.push({ mode: 'half', type: halfExt.HALF_FLOAT_OES, filter: halfLinearExt ? gl.LINEAR : gl.NEAREST, filterName: halfLinearExt ? 'linear' : 'nearest' });
    if (floatExt) cands.push({ mode: 'float', type: gl.FLOAT, filter: floatLinearExt ? gl.LINEAR : gl.NEAREST, filterName: floatLinearExt ? 'linear' : 'nearest' });
    cands.push({ mode: 'byte', type: gl.UNSIGNED_BYTE, filter: gl.LINEAR, filterName: 'linear' });
    gl.activeTexture(gl.TEXTURE0);
    for (var i = 0; i < cands.length; i++) {
      var t = gl.createTexture();
      gl.bindTexture(gl.TEXTURE_2D, t);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, cands[i].filter);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, cands[i].filter);
      var allocOK = true;
      try {
        gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, bw, bh, 0, gl.RGBA, cands[i].type, null);
      } catch (eA) { allocOK = false; }
      if (!allocOK) { gl.deleteTexture(t); continue; }
      var fb = gl.createFramebuffer();
      gl.bindFramebuffer(gl.FRAMEBUFFER, fb);
      gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, t, 0);
      var ok = false;
      try { ok = (gl.checkFramebufferStatus(gl.FRAMEBUFFER) === gl.FRAMEBUFFER_COMPLETE); } catch (eC) { ok = false; }
      gl.bindFramebuffer(gl.FRAMEBUFFER, null);
      if (ok) { caus.fbo = fb; caus.tex = t; caus.mode = cands[i].mode; caus.filter = cands[i].filterName; break; }
      gl.deleteFramebuffer(fb);
      gl.deleteTexture(t);
    }
    // 入力テクスチャの束縛を戻す
    gl.activeTexture(gl.TEXTURE0);
    if (tex) gl.bindTexture(gl.TEXTURE_2D, tex);
  }

  function render(tSec) {
    if (!gl || !prog || mediaKind === 'none') return;
    if (mediaKind === 'video' && vidEl.readyState >= 2 && !vidEl.paused && !vidEl.ended) uploadSource();
    // FIX-2: 波の見えは水状態(慣性追従後)からのみ決まる
    var now = performance.now() / 1000;
    var dt = lastT < 0 ? 0 : Math.min(now - lastT, 0.25);
    lastT = now;
    // FIX-16: FPS計測(直近2秒窓の実フレーム数ベース。EMAは廃止)。
    // 毎フレーム描画時刻を記録し、2秒より古い履歴は押し出す(自然押し出し)。
    // fps値の再計算は0.5秒ごと。自動制御・露出・表示はこの値を使う。
    frameTimes16.push(now);
    while (frameTimes16.length && frameTimes16[0] <= now - 2) frameTimes16.shift();
    if (now - lastFpsUpdateSec >= 0.5) {
      windowFps16 = fix16WindowFps(frameTimes16, now);
      lastFpsUpdateSec = now;
    }
    // FIX-15: 自動制御(既定ON)。1秒毎に評価し、45未満で1段下げ、
    // 57超が3回続けば1段上げ。切替は最短2秒間隔(判定は fix15Decide)。
    // FIX-16: 評価に使うfpsは窓ベースの windowFps16。閾値・間隔・回数は不変。
    if (autoQuality && qualityMode === 'auto' && now - lastQualityEvalSec >= 1) {
      var dec15 = fix15Decide(quality, windowFps16, upStreak, now - lastQualitySwitchSec);
      upStreak = dec15.streakNext;
      lastQualityEvalSec = now;
      if (dec15.doSwitch && applyQuality(dec15.next, 'auto')) {
        lastQualitySwitchSec = now;
        upStreak = 0;
      }
    }
    var ws = waterState.update(waterTarget, dt);
    var st = state.render.strength == null ? 1 : state.render.strength;
    // FIX-7: 前方パスが使えるか(バッファ確保済み+スプラットプログラム正常)
    var useForward = (caus.mode !== 'none' && caus.fbo && splatProg) ? 1 : 0;
    var causAvg = 1.0;
    if (useForward === 1) {
      // --- FIX-9 前方パス: 水面グリッドの三角形メッシュを屈折着地点へ ---
      // 面積比つきで加算描画(ブレンド ONE, ONE)。原子正規化はエネルギー保存
      // から解析的に 1.0(三角形は水面積を過不足なく敷き詰める)。
      // 変位ゲインは uDispK(FIX11_DISP_K、単一箇所・FIX-11再較正値)、ピーク正則化は
      // uDelta2Fwd(=1/VCAP^2、FIX9_VCAP_FWD由来・フラグメント側で使用)。帯域別変位ゲインは
      // uDispGL/GM/GH(FIX10_DISP_GL/GM/GH、変位のみ・法線等は全帯域のまま)。
      var dispK = 2.0, delta2Fwd = 0.04, dispGL = 1.0, dispGM = 0.55, dispGH = 0.25;
      try {
        if (window.CausticsEnv.FIX11_DISP_K) dispK = window.CausticsEnv.FIX11_DISP_K;
        else if (window.CausticsEnv.FIX10_DISP_K != null) dispK = window.CausticsEnv.FIX10_DISP_K;
        else if (window.CausticsEnv.FIX9_DISP_K) dispK = window.CausticsEnv.FIX9_DISP_K;
        if (window.CausticsEnv.FIX9_VCAP_FWD) delta2Fwd = 1 / (window.CausticsEnv.FIX9_VCAP_FWD * window.CausticsEnv.FIX9_VCAP_FWD);
        if (window.CausticsEnv.FIX10_DISP_GL != null) dispGL = window.CausticsEnv.FIX10_DISP_GL;
        if (window.CausticsEnv.FIX10_DISP_GM != null) dispGM = window.CausticsEnv.FIX10_DISP_GM;
        if (window.CausticsEnv.FIX10_DISP_GH != null) dispGH = window.CausticsEnv.FIX10_DISP_GH;
      } catch (eK) { /* 既定のまま */ }
      causAvg = 1.0;
      gl.bindFramebuffer(gl.FRAMEBUFFER, caus.fbo);
      gl.viewport(0, 0, caus.w, caus.h);
      gl.useProgram(splatProg);
      if (splatLoc >= 0) {
        gl.bindBuffer(gl.ARRAY_BUFFER, gridBuf);
        gl.vertexAttribPointer(splatLoc, 2, gl.FLOAT, false, 0, 0);
      }
      gl.uniform2f(splatUni.uRes, caus.w, caus.h);
      gl.uniform1f(splatUni.uTime, tSec % LOOP);
      gl.uniform3f(splatUni.uSunDir, params.sunDir[0], params.sunDir[1], params.sunDir[2]);
      gl.uniform1f(splatUni.uSpecL, ws.specL);
      gl.uniform1f(splatUni.uSpecAmp, ws.specAmp);
      gl.uniform2f(splatUni.uWindDir, ws.dirX, ws.dirY);
      gl.uniform1f(splatUni.uFine, ws.fine);
      gl.uniform2f(splatUni.uFlowPerSec, ws.flowX, ws.flowY);
      gl.uniform1f(splatUni.uDepthM, params.depthM);
      gl.uniform1f(splatUni.uRefract, params.refract * st);
      if (splatUni.uDgCap) gl.uniform1f(splatUni.uDgCap, params.dgCap == null ? 4.2 : params.dgCap);
      gl.uniform1f(splatUni.uAspect, canvas.width / Math.max(canvas.height, 1));
      gl.uniform1f(splatUni.uMeshN, gridN);
      if (splatUni.uDispK) gl.uniform1f(splatUni.uDispK, dispK);
      if (splatUni.uDelta2Fwd) gl.uniform1f(splatUni.uDelta2Fwd, delta2Fwd);
      if (splatUni.uDispGL) gl.uniform1f(splatUni.uDispGL, dispGL);
      if (splatUni.uDispGM) gl.uniform1f(splatUni.uDispGM, dispGM);
      if (splatUni.uDispGH) gl.uniform1f(splatUni.uDispGH, dispGH);
      gl.clearColor(0, 0, 0, 0);
      gl.clear(gl.COLOR_BUFFER_BIT);
      gl.enable(gl.BLEND);
      gl.blendFunc(gl.ONE, gl.ONE);
      gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, idxBuf);
      gl.drawElements(gl.TRIANGLES, idxCount, gl.UNSIGNED_SHORT, 0);
      gl.disable(gl.BLEND);
      gl.bindFramebuffer(gl.FRAMEBUFFER, null);
    }
    // --- 合成パス(現行HDR乗算。光量は前方バッファ由来、不可時は後方フォールバック) ---
    gl.viewport(0, 0, canvas.width, canvas.height);
    gl.useProgram(prog);
    if (quadLoc >= 0) {
      gl.bindBuffer(gl.ARRAY_BUFFER, quadBuf);
      gl.vertexAttribPointer(quadLoc, 2, gl.FLOAT, false, 0, 0);
    }
    gl.activeTexture(gl.TEXTURE0);
    gl.bindTexture(gl.TEXTURE_2D, tex);
    gl.uniform1i(uni.uTex, 0);
    gl.uniform1f(uni.uMediaAspect, mediaAspect);
    gl.uniform2f(uni.uRes, canvas.width, canvas.height);
    gl.uniform1f(uni.uTime, tSec % LOOP);
    gl.uniform3f(uni.uSunDir, params.sunDir[0], params.sunDir[1], params.sunDir[2]);
    gl.uniform3f(uni.uSunColor, params.sunColor[0], params.sunColor[1], params.sunColor[2]);
    gl.uniform1f(uni.uSunI, params.sunI);
    gl.uniform1f(uni.uDiffuse, params.diffuse);
    // FIX-2: 波の見えは水状態(慣性追従後)からのみ決まる(本関数の先頭で追従済みの ws を使う)
    gl.uniform1f(uni.uSpecL, ws.specL);
    gl.uniform1f(uni.uSpecAmp, ws.specAmp);
    gl.uniform2f(uni.uWindDir, ws.dirX, ws.dirY);
    gl.uniform1f(uni.uFine, ws.fine);
    gl.uniform2f(uni.uFlowPerSec, ws.flowX, ws.flowY);
    gl.uniform1f(uni.uDepthM, params.depthM);
    gl.uniform1f(uni.uRefract, params.refract * st);
    gl.uniform1f(uni.uChroma, params.chroma);
    gl.uniform1f(uni.uCausticsI, params.causticsI * st);
    gl.uniform1f(uni.uGlitterI, params.glitterI * st);
    gl.uniform3f(uni.uAbsorb, params.absorb[0], params.absorb[1], params.absorb[2]);
    gl.uniform3f(uni.uScatter, params.scatter[0], params.scatter[1], params.scatter[2]);
    gl.uniform1f(uni.uDepthShade, params.depthShade);
    gl.uniform1f(uni.uTurbid, params.turbid);
    gl.uniform1f(uni.uExposure, state.render.exposure == null ? 1 : state.render.exposure);
    // FIX-5: HDR振幅 kc と焦点深度。uniformが無い旧シェーダでも落ちないようガード。
    if (uni.uKc) gl.uniform1f(uni.uKc, params.causticsKc == null ? 0.8 : params.causticsKc);
    if (uni.uFocusM) gl.uniform1f(uni.uFocusM, params.focusM == null ? 0.8 : params.focusM);
    // FIX-6: dgL圧縮上限。uniformが無い旧シェーダでも落ちないようガード。
    if (uni.uDgCap) gl.uniform1f(uni.uDgCap, params.dgCap == null ? 4.2 : params.dgCap);
    // FIX-7: 前方蓄積バッファの束縛(テクスチャユニット1)と正規化・分岐uniform。
    if (uni.uUseForward) gl.uniform1f(uni.uUseForward, useForward);
    if (uni.uCausticsAvg) gl.uniform1f(uni.uCausticsAvg, causAvg);
    if (uni.uCausticsTexel) gl.uniform2f(uni.uCausticsTexel,
      1 / Math.max(caus.w, 1), 1 / Math.max(caus.h, 1));
    if (uni.uCaustics) {
      gl.activeTexture(gl.TEXTURE1);
      if (useForward === 1 && caus.tex) gl.bindTexture(gl.TEXTURE_2D, caus.tex);
      gl.uniform1i(uni.uCaustics, 1);
      gl.activeTexture(gl.TEXTURE0);
      gl.bindTexture(gl.TEXTURE_2D, tex);
    }
    gl.drawArrays(gl.TRIANGLES, 0, 3);
    frameCount++;
    if (window.__caustics) {
      window.__caustics.frameCount = frameCount;
      window.__caustics.forward = useForward;
      window.__caustics.caustics = { mode: caus.mode, w: caus.w, h: caus.h, grid: gridN, filter: caus.filter };
      // FIX-15: FPS(小数1桁)と品質段を露出する。
      // FIX-16: 露出するfpsは窓ベースの windowFps16(0.5秒ごとに更新)。
      window.__caustics.fps = Math.round(windowFps16 * 10) / 10;
      window.__caustics.quality = quality;
    }
    // FIX-16: キャンバス右下に小さく薄い文字で「FPS・品質段」を重ねて表示する。
    // 更新は0.5秒ごと(FIX-15の15フレーム間引きは廃止。0.5秒間隔で十分軽い)。
    if (now - lastOverlaySec16 >= 0.5 || frameCount < 2) {
      var ov15 = document.getElementById('fps-overlay');
      if (ov15) ov15.textContent = (Math.round(windowFps16 * 10) / 10).toFixed(1) + ' fps ' + quality;
      lastOverlaySec16 = now;
    }
    // headless検証用フック: DOMに状態を露出(dump-domで回収できる)
    var vo = document.getElementById('verify-out');
    if (vo) {
      vo.textContent = JSON.stringify({
        shaderOK: !!(window.__caustics && window.__caustics.shaderOK),
        frameCount: frameCount, w: canvas.width, h: canvas.height,
        forward: useForward, causticsMode: caus.mode,
        fps: Math.round(windowFps16 * 10) / 10, quality: quality
      });
    }
  }

  function loop() {
    var t = (performance.now() / 1000) % LOOP;
    render(t);
    requestAnimationFrame(loop);
  }

  // ---------- 状態 ⇄ UI ----------
  function readUI() {
    state.sun.date = $('in-date').value || state.sun.date;
    state.sun.time = $('in-time').value || state.sun.time;
    state.sun.latitudeDeg = parseFloat($('in-lat').value);
    state.sun.longitudeDeg = parseFloat($('in-lon').value);
    state.water.depthM = parseFloat($('in-depth').value);
    state.water.transparency = parseFloat($('in-transp').value);
    state.water.turbidity = parseFloat($('in-turbid').value);
    state.water.pigment = parseFloat($('in-pigment').value);
    state.windWave.windSpeedMS = parseFloat($('in-wind').value);
    state.windWave.windDirDeg = parseFloat($('in-windir').value);
    state.windWave.fineness = parseFloat($('in-fine').value);
    state.windWave.flowSpeed = parseFloat($('in-flow').value);
    state.weather.cloudPct = parseFloat($('in-cloud').value);
    state.weather.clarity = parseFloat($('in-clarity').value);
    state.render.strength = parseFloat($('in-strength').value);
    state.render.exposure = parseFloat($('in-exposure').value);
    state.name = $('recipe-name').value || '無題の環境';
  }

  function applyUI() {
    $('in-date').value = state.sun.date;
    $('in-time').value = state.sun.time;
    $('in-lat').value = state.sun.latitudeDeg;
    $('in-lon').value = state.sun.longitudeDeg;
    $('in-depth').value = state.water.depthM;
    $('in-transp').value = state.water.transparency;
    $('in-turbid').value = state.water.turbidity;
    $('in-pigment').value = state.water.pigment;
    $('in-wind').value = state.windWave.windSpeedMS;
    $('in-windir').value = state.windWave.windDirDeg;
    $('in-fine').value = state.windWave.fineness;
    $('in-flow').value = state.windWave.flowSpeed;
    $('in-cloud').value = state.weather.cloudPct;
    $('in-clarity').value = state.weather.clarity;
    $('in-strength').value = state.render.strength;
    $('in-exposure').value = state.render.exposure;
    $('recipe-name').value = state.name || '無題の環境';
    buildCouplingEditor();
    refresh();
    // 離散的な状態変化(プリセット・レシピ読込・実在モード反映)は水状態に即時反映し、
    // 検証・書き出しの決定性を保つ。スライダー操作中は render 側の追従で慣性が出る。
    waterState.snap(waterTarget);
  }

  function fmt(x, d) { return (Math.round(x * Math.pow(10, d)) / Math.pow(10, d)).toFixed(d); }

  function refresh() {
    // スライダー表示
    $('v-strength').textContent = fmt(state.render.strength, 2);
    $('v-exposure').textContent = fmt(state.render.exposure, 2);
    $('v-transp').textContent = fmt(state.water.transparency, 2);
    $('v-turbid').textContent = fmt(state.water.turbidity, 2);
    $('v-pigment').textContent = fmt(state.water.pigment, 2);
    $('v-fine').textContent = fmt(state.windWave.fineness, 2);
    $('v-flow').textContent = fmt(state.windWave.flowSpeed, 2);
    $('v-clarity').textContent = fmt(state.weather.clarity, 2);
    // 物理解決
    params = window.CausticsEnv.resolveEnvironment(state);
    waterTarget = waterTargetFromParams(params);
    $('sun-alt').textContent = fmt(params.sun.altitudeDeg, 1);
    $('sun-az').textContent = fmt(params.sun.azimuthDeg, 1);
    $('wave-h').textContent = fmt(params.waveHeightM, 2);
    $('wx-sun').textContent = fmt(params.sunI, 2);
    $('wx-dif').textContent = fmt(params.diffuse, 2);
    $('d-cloud').textContent = state.weather.cloudPct;
    $('d-wind').textContent = state.windWave.windSpeedMS;
    $('d-wave').textContent = fmt(params.waveHeightM, 2);
    $('d-alt').textContent = fmt(params.sun.altitudeDeg, 1);
    $('d-caus').textContent = fmt(params.causticsI, 2);
    $('d-depth').textContent = state.water.depthM;
    $('d-turb').textContent = fmt(state.water.turbidity, 2);
    $('env-readout').textContent =
      '太陽高度 ' + fmt(params.sun.altitudeDeg, 1) + '° / 方位 ' + fmt(params.sun.azimuthDeg, 1) +
      '° / 直射 ' + fmt(params.sunI, 2) + ' / 波高目安 ' + fmt(params.waveHeightM, 2) +
      'm / 集光 ' + fmt(params.causticsI, 2) + ' / グリッター ' + fmt(params.glitterI, 2);
  }

  function setStatus(msg) { $('gl-status').textContent = msg; }

  // ---------- 結合エディタ ----------
  function buildCouplingEditor() {
    var ul = $('coupling-list');
    ul.innerHTML = '';
    state.couplings.forEach(function (c, idx) {
      var li = document.createElement('li');
      var lab = document.createElement('label');
      var cb = document.createElement('input');
      cb.type = 'checkbox'; cb.checked = !!c.enabled;
      cb.setAttribute('aria-label', c.label + 'のON/OFF');
      cb.addEventListener('change', function () { state.couplings[idx].enabled = cb.checked; refresh(); });
      lab.appendChild(cb);
      lab.appendChild(document.createTextNode(c.label + ' (' + c.from + '→' + c.to + ')'));
      var rg = document.createElement('input');
      rg.type = 'range'; rg.min = '0'; rg.max = '2'; rg.step = '0.05'; rg.value = c.weight;
      rg.setAttribute('aria-label', c.label + 'の重み');
      var sp = document.createElement('span');
      sp.textContent = fmt(c.weight, 2);
      rg.addEventListener('input', function () {
        state.couplings[idx].weight = parseFloat(rg.value);
        sp.textContent = fmt(state.couplings[idx].weight, 2);
        refresh();
      });
      li.appendChild(lab); li.appendChild(rg); li.appendChild(sp);
      ul.appendChild(li);
    });
  }

  // ---------- メディア ----------
  function makeSample() {
    var c = document.createElement('canvas');
    c.width = 640; c.height = 400;
    // ?seed=N で決定論的サンプル(検証用)。通常はMath.random。
    var seed = NaN;
    try { seed = parseInt(new URLSearchParams(location.search).get('seed'), 10); } catch (e) { seed = NaN; }
    var st = isNaN(seed) ? 0 : (seed >>> 0) || 1;
    var R = function () {
      if (!st) return Math.random();
      st ^= st << 13; st >>>= 0; st ^= st >> 17; st ^= st << 5; st >>>= 0;
      return (st % 100000) / 100000;
    };
    var g = c.getContext('2d');
    var sky = g.createLinearGradient(0, 0, 0, 400);
    sky.addColorStop(0, '#7fb6c9'); sky.addColorStop(0.55, '#c8d8b8'); sky.addColorStop(1, '#3f6b4f');
    g.fillStyle = sky; g.fillRect(0, 0, 640, 400);
    g.fillStyle = '#fff6d8';
    g.beginPath(); g.arc(470, 90, 34, 0, Math.PI * 2); g.fill();
    // 川底の石っぽい斑点
    var i, x, y, r;
    for (i = 0; i < 260; i++) {
      x = R() * 640; y = 200 + R() * 200; r = 3 + R() * 14;
      g.fillStyle = 'rgba(' + (40 + R() * 60 | 0) + ',' + (90 + R() * 60 | 0) + ',' + (70 + R() * 40 | 0) + ',0.55)';
      g.beginPath(); g.arc(x, y, r, 0, Math.PI * 2); g.fill();
    }
    // 花(型ガラス越しの花束のような色点)
    var cols = ['#e05575', '#f2c14e', '#e8e8e8', '#c65b8f'];
    for (i = 0; i < 24; i++) {
      x = 60 + R() * 520; y = 60 + R() * 180; r = 8 + R() * 16;
      g.fillStyle = cols[i % cols.length];
      g.beginPath(); g.arc(x, y, r, 0, Math.PI * 2); g.fill();
    }
    sampleCanvas = c;
    mediaKind = 'sample';
    sizeCanvasToMedia(); uploadSource();
    setStatus('サンプル画像でプレビュー中。お手持ちの画像・動画も開けます。');
  }

  function openMediaFile(f) {
    if (!f) return;
    if (f.type.indexOf('video') === 0) {
      if (imgEl) imgEl = null;
      vidEl.src = URL.createObjectURL(f);
      vidEl.hidden = false;
      vidEl.play().catch(function () { /* 自動再生不可時は手動 */ });
      var onMeta = function () {
        mediaKind = 'video';
        sizeCanvasToMedia(); uploadSource();
        setStatus('動画を読み込みました(再生しながら適用します)。音声はミュートです。');
      };
      vidEl.onloadedmetadata = onMeta;
      if (vidEl.readyState >= 1) onMeta();
    } else {
      vidEl.pause(); vidEl.removeAttribute('src'); vidEl.hidden = true;
      var url = URL.createObjectURL(f);
      var im = new Image();
      im.onload = function () {
        imgEl = im; mediaKind = 'image';
        sizeCanvasToMedia(); uploadSource();
        setStatus('画像を読み込みました。');
        URL.revokeObjectURL(url);
      };
      im.onerror = function () { setStatus('画像の読み込みに失敗しました。'); };
      im.src = url;
    }
  }

  // ---------- 書き出し ----------
  function savePNG() {
    if (mediaKind === 'none') { setStatus('先に画像・動画を開くか「サンプルで試す」を押してください。'); return; }
    render((performance.now() / 1000) % LOOP);
    canvas.toBlob(function (b) {
      if (!b) { setStatus('PNGの保存に失敗しました。'); return; }
      var a = document.createElement('a');
      a.href = URL.createObjectURL(b);
      a.download = (state.name || 'caustics') + '.png';
      a.click();
      setTimeout(function () { URL.revokeObjectURL(a.href); }, 5000);
      setStatus('静止画PNGを保存しました。');
    }, 'image/png');
  }

  function pickMime() {
    var cands = ['video/webm;codecs=vp9', 'video/webm;codecs=vp8', 'video/webm'];
    for (var i = 0; i < cands.length; i++) {
      try { if (window.MediaRecorder && MediaRecorder.isTypeSupported(cands[i])) return cands[i]; } catch (e) { /* noop */ }
    }
    return '';
  }

  function saveWebM() {
    if (mediaKind === 'none') { setStatus('先に画像・動画を開くか「サンプルで試す」を押してください。'); return; }
    if (!window.MediaRecorder || !canvas.captureStream) {
      setStatus('このブラウザでは動画書き出し(MediaRecorder/captureStream)に未対応です。');
      return;
    }
    var mime = pickMime();
    var stream = canvas.captureStream(30);
    var rec;
    try {
      rec = mime ? new MediaRecorder(stream, { mimeType: mime, videoBitsPerSecond: 8000000 }) : new MediaRecorder(stream);
    } catch (e) { setStatus('MediaRecorderの初期化に失敗: ' + e.message); return; }
    var chunks = [];
    rec.ondataavailable = function (e) { if (e.data && e.data.size) chunks.push(e.data); };
    rec.onstop = function () {
      var b = new Blob(chunks, { type: 'video/webm' });
      var a = document.createElement('a');
      a.href = URL.createObjectURL(b);
      a.download = (state.name || 'caustics') + '-loop8s.webm';
      a.click();
      setTimeout(function () { URL.revokeObjectURL(a.href); }, 5000);
      $('webm-status').textContent = 'ループ動画WebM(8秒)を保存しました。';
      $('btn-webm').disabled = false;
    };
    $('btn-webm').disabled = true;
    $('webm-status').textContent = '録画中…8秒お待ちください(ちょうど1ループ分です)。';
    rec.start(500);
    setTimeout(function () { try { rec.stop(); } catch (e) { $('btn-webm').disabled = false; } }, 8000);
  }

  // ---------- レシピUI ----------
  function refreshLsList() {
    var sel = $('ls-list');
    sel.innerHTML = '';
    var names = window.CausticsRecipe.lsList();
    if (!names.length) {
      var o = document.createElement('option');
      o.value = ''; o.textContent = '(保存なし)';
      sel.appendChild(o);
      return;
    }
    names.forEach(function (n) {
      var o = document.createElement('option');
      o.value = n; o.textContent = n;
      sel.appendChild(o);
    });
  }

  function buildPresets() {
    var row = $('preset-row');
    row.innerHTML = '';
    window.CausticsRecipe.presets().forEach(function (p) {
      var b = document.createElement('button');
      b.type = 'button';
      b.textContent = p.name;
      b.addEventListener('click', function () {
        state = JSON.parse(JSON.stringify(p));
        applyUI();
        setStatus('プリセット「' + p.name + '」を適用しました。');
      });
      row.appendChild(b);
    });
  }

  function downloadFile(filename, text, mime) {
    var b = new Blob([text], { type: mime || 'application/json' });
    var a = document.createElement('a');
    a.href = URL.createObjectURL(b);
    a.download = filename;
    a.click();
    setTimeout(function () { URL.revokeObjectURL(a.href); }, 5000);
  }

  // ---------- 実在モード ----------
  function runRealMode() {
    var place = $('in-place').value.trim();
    var d = $('in-rdate').value, t = $('in-rtime').value;
    var st = $('real-status');
    if (!place) { st.textContent = '地名を入力してください。'; return; }
    st.textContent = '取得中…(ジオコーディング→気象データ)';
    window.CausticsEnv.fetchRealEnvironment(place, d, t).then(function (res) {
      state.sun.date = d;
      state.sun.time = t;
      state.sun.latitudeDeg = Math.round(res.latitude * 1000) / 1000;
      state.sun.longitudeDeg = Math.round(res.longitude * 1000) / 1000;
      state.place = { label: res.label, latitude: res.latitude, longitude: res.longitude, datetime: d + ' ' + t };
      if (res.weather) {
        state.weather.cloudPct = Math.round(res.weather.cloudPct);
        state.windWave.windSpeedMS = Math.round(res.weather.windSpeedMS * 10) / 10;
        state.windWave.windDirDeg = Math.round(res.weather.windDirDeg);
        st.textContent = '取り込み完了: ' + res.label + ' / 雲' + state.weather.cloudPct +
          '% / 風' + state.windWave.windSpeedMS + 'm/s(' + state.windWave.windDirDeg + '°)。手で微調整できます。';
      } else {
        st.textContent = '緯度経度のみ取り込みました(' + res.label + ')。' +
          '天気の取得に失敗したため太陽計算だけ成立・天気は手動値のままです' +
          (res.weatherError ? '(' + res.weatherError + ')' : '') + '。';
      }
      applyUI();
    }).catch(function (e) {
      st.textContent = '取り込み失敗: ' + (e && e.message ? e.message : e) + '。入力値を確認してください。';
    });
  }

  // ---------- 配線 ----------
  function bind() {
    var ids = ['in-date', 'in-time', 'in-lat', 'in-lon', 'in-depth', 'in-transp',
      'in-turbid', 'in-pigment', 'in-wind', 'in-windir', 'in-fine', 'in-flow',
      'in-cloud', 'in-clarity', 'in-strength', 'in-exposure'];
    ids.forEach(function (id) {
      $(id).addEventListener('input', function () { readUI(); refresh(); });
      $(id).addEventListener('change', function () { readUI(); refresh(); });
    });
    var modes = document.getElementsByName('editmode');
    for (var i = 0; i < modes.length; i++) {
      modes[i].addEventListener('change', function () {
        var adv = document.querySelector('input[name="editmode"]:checked').value === 'advanced';
        $('coupling-editor').hidden = !adv;
      });
    }
    $('btn-open').addEventListener('click', function () { $('file-media').click(); });
    $('file-media').addEventListener('change', function (e) { openMediaFile(e.target.files[0]); e.target.value = ''; });
    $('btn-sample').addEventListener('click', makeSample);
    $('btn-png').addEventListener('click', savePNG);
    $('btn-webm').addEventListener('click', saveWebM);
    $('btn-real').addEventListener('click', runRealMode);
    $('btn-save-ls').addEventListener('click', function () {
      readUI(); refresh();
      try {
        window.CausticsRecipe.lsSave(state.name, state);
        refreshLsList();
        setStatus('ブラウザに保存しました: ' + state.name);
      } catch (e) { setStatus('保存に失敗: ' + e.message); }
    });
    $('btn-load-ls').addEventListener('click', function () {
      var n = $('ls-list').value;
      if (!n) return;
      try { state = window.CausticsRecipe.lsLoad(n); applyUI(); setStatus('読み込みました: ' + n); }
      catch (e) { setStatus('読み込み失敗: ' + e.message); }
    });
    $('btn-del-ls').addEventListener('click', function () {
      var n = $('ls-list').value;
      if (!n) return;
      window.CausticsRecipe.lsDelete(n);
      refreshLsList();
      setStatus('削除しました: ' + n);
    });
    $('btn-download').addEventListener('click', function () {
      readUI(); refresh();
      downloadFile((state.name || 'caustics') + '.json', window.CausticsRecipe.serialize(state));
      setStatus('レシピJSONをファイルに保存しました。');
    });
    $('btn-upload').addEventListener('click', function () { $('file-recipe').click(); });
    $('file-recipe').addEventListener('change', function (e) {
      var f = e.target.files[0];
      if (!f) return;
      var r = new FileReader();
      r.onload = function () {
        try { state = window.CausticsRecipe.parse(r.result); applyUI(); setStatus('レシピを読み込みました: ' + state.name); }
        catch (err) { setStatus('レシピの読み込み失敗: ' + err.message); }
      };
      r.readAsText(f);
      e.target.value = '';
    });
    // FIX-15: 描画品質セレクト(自動/高/中/低)。手動を選んだら自動制御は止まり、
    // 「自動」に戻すと再開する。
    var qsel = $('in-quality');
    if (qsel) {
      qsel.value = qualityMode;
      qsel.addEventListener('change', function () {
        var v = qsel.value;
        if (v === 'auto') {
          qualityMode = 'auto'; autoQuality = true;
          lastQualityEvalSec = -1e9; upStreak = 0;
        } else if (FIX15_QUALITY_LADDER[v]) {
          qualityMode = v; autoQuality = false;
          if (applyQuality(v, 'manual')) lastQualitySwitchSec = performance.now() / 1000;
        }
      });
    }
    // 動画テクスチャ更新用
    vidEl.addEventListener('play', function () { sizeCanvasToMedia(); });
  }

  // ---------- 起動 ----------
  function boot() {
    var ok = initGL();
    bind();
    buildPresets();
    refreshLsList();
    // 深リンク: ?preset=N でN番目のプリセットを初期適用(検証・共有用)
    try {
      var q = new URLSearchParams(location.search).get('preset');
      if (q !== null && q !== '') {
        var idx = parseInt(q, 10);
        var ps = window.CausticsRecipe.presets();
        if (!isNaN(idx) && ps[idx]) state = JSON.parse(JSON.stringify(ps[idx]));
      }
    } catch (e) { /* 無視して既定状態 */ }
    applyUI();
    if (ok) {
      // 自動でサンプルを出して初回描画(受け入れ確認しやすくする)
      makeSample();
      requestAnimationFrame(loop);
    }
  }
  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', boot);
  } else { boot(); }
})();
