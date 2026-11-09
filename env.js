/* 環境ノードの物理解決 + 実在モード(Open-Meteo) (classic script)
 * window.CausticsEnv = { COUPLINGS, solarPosition, sunColorFromAltitude,
 *   resolveEnvironment, fetchRealEnvironment }
 */
(function () {
  'use strict';

  // 結合の一覧(上級モードで ON/OFF・重み 0-2 を編集)
  var COUPLINGS = [
    { id: 'cloud-sun',   from: '天気',   to: '太陽',       label: '雲→太陽減衰',       enabled: true, weight: 1.0 },
    { id: 'wind-wave',   from: '風と波', to: '風と波',      label: '風→波高',           enabled: true, weight: 1.0 },
    { id: 'sun-caustics', from: '太陽',  to: 'コースティクス', label: '太陽→コースティクス強度', enabled: true, weight: 1.0 },
    { id: 'turbid-caustics', from: '水', to: 'コースティクス', label: '濁り→コースティクス減衰', enabled: true, weight: 1.0 },
    { id: 'depth-tint',  from: '水',     to: '水の色',       label: '水深→吸収ティント', enabled: true, weight: 1.0 },
    { id: 'cloud-glitter', from: '天気', to: 'グリッター',   label: '雲→グリッター拡散', enabled: true, weight: 1.0 }
  ];

  function couplingById(list, id) {
    for (var i = 0; i < list.length; i++) {
      if (list[i].id === id) return list[i];
    }
    return { enabled: false, weight: 0 };
  }
  function cw(list, id) { // 結合の効き具合 (OFFなら0、ONなら重み)
    var c = couplingById(list, id);
    return c.enabled ? Math.max(0, Math.min(2, c.weight)) : 0;
  }

  function toRad(d) { return d * Math.PI / 180; }
  function toDeg(r) { return r * 180 / Math.PI; }
  function clamp(x, a, b) { return Math.max(a, Math.min(b, x)); }

  // 太陽位置の天文計算(近似・NOAA系)。
  // dateStr 'YYYY-MM-DD', timeStr 'HH:MM' は「その場所の地方時」とみなし、
  // UTCオフセットは経度からの推定(15度=1hに丸め)を使う。作品用途の近似。
  // 戻り値 { altitudeDeg, azimuthDeg } (方位は北=0・東回り)
  function solarPosition(dateStr, timeStr, latDeg, lonDeg) {
    var dp = String(dateStr || '').split('-');
    var tp = String(timeStr || '').split(':');
    var Y = parseInt(dp[0], 10) || 2000;
    var M = parseInt(dp[1], 10) || 1;
    var D = parseInt(dp[2], 10) || 1;
    var hh = parseInt(tp[0], 10) || 0;
    var mm = parseInt(tp[1], 10) || 0;
    var zoneH = Math.round((lonDeg || 0) / 15); // 推定標準時
    // UTCに直してユリウス日
    var utcH = hh - zoneH + mm / 60;
    var a = Math.floor((14 - M) / 12);
    var y = Y + 4800 - a;
    var m = M + 12 * a - 3;
    var jdn = D + Math.floor((153 * m + 2) / 5) + 365 * y +
      Math.floor(y / 4) - Math.floor(y / 100) + Math.floor(y / 400) - 32045;
    var jd = jdn + (utcH - 12) / 24;
    var n = jd - 2451545.0;
    var L = (280.460 + 0.9856474 * n) % 360;
    var g = toRad((357.528 + 0.9856003 * n) % 360);
    var lambda = toRad(L + 1.915 * Math.sin(g) + 0.020 * Math.sin(2 * g));
    var eps = toRad(23.439 - 0.0000004 * n);
    var alpha = Math.atan2(Math.cos(eps) * Math.sin(lambda), Math.cos(lambda));
    var delta = Math.asin(Math.sin(eps) * Math.sin(lambda));
    // グリニッジ恒星時→時角
    var gmst = (18.697374558 + 24.06570982441908 * n) % 24;
    if (gmst < 0) gmst += 24;
    var lst = (gmst + (lonDeg || 0) / 15) % 24;
    var H = toRad(lst * 15) - alpha;
    var lat = toRad(latDeg || 0);
    var sinAlt = Math.sin(lat) * Math.sin(delta) + Math.cos(lat) * Math.cos(delta) * Math.cos(H);
    var alt = Math.asin(clamp(sinAlt, -1, 1));
    var az = Math.atan2(Math.sin(H), Math.cos(H) * Math.sin(lat) - Math.tan(delta) * Math.cos(lat));
    var azDeg = (toDeg(az) + 180 + 360) % 360; // 北=0
    return { altitudeDeg: toDeg(alt), azimuthDeg: azDeg };
  }

  // 太陽高度から光色。低い=橙。地平線下は薄暗い青。
  function sunColorFromAltitude(altDeg) {
    var dayAmt = clamp((altDeg + 2) / 37, 0, 1); // -2°→0, 35°→1
    dayAmt = dayAmt * dayAmt * (3 - 2 * dayAmt); // smoothstep
    var low = [1.0, 0.42, 0.13];
    var high = [1.0, 0.97, 0.92];
    var c = [
      low[0] + (high[0] - low[0]) * dayAmt,
      low[1] + (high[1] - low[1]) * dayAmt,
      low[2] + (high[2] - low[2]) * dayAmt
    ];
    if (altDeg < -2) { // 夜: ほのかな青だけ
      var k = clamp((altDeg + 12) / 10, 0, 1);
      c = [0.10 * k + 0.02, 0.14 * k + 0.02, 0.22 * k + 0.03];
    }
    return c;
  }

  // state(レシピ) → シェーダuniform群 + 表示用派生値
  function resolveEnvironment(state) {
    var couplings = state.couplings || COUPLINGS;
    var cloud = clamp((state.weather.cloudPct || 0) / 100, 0, 1);
    var clarity = clamp(state.weather.clarity == null ? 0.8 : state.weather.clarity, 0, 1);
    var wind = Math.max(0, state.windWave.windSpeedMS || 0);
    var windDir = toRad(state.windWave.windDirDeg || 0);
    var fineness = clamp(state.windWave.fineness == null ? 0.5 : state.windWave.fineness, 0, 1);
    var flowSpd = clamp(state.windWave.flowSpeed == null ? 0.5 : state.windWave.flowSpeed, 0, 1);
    var depth = Math.max(0.1, state.water.depthM || 1);
    var transp = clamp(state.water.transparency == null ? 0.7 : state.water.transparency, 0, 1);
    var turbid = clamp(state.water.turbidity == null ? 0.3 : state.water.turbidity, 0, 1);
    var pigment = clamp(state.water.pigment == null ? 0.3 : state.water.pigment, 0, 1);

    var sun = solarPosition(state.sun.date, state.sun.time, state.sun.latitudeDeg, state.sun.longitudeDeg);
    var sunColor = sunColorFromAltitude(sun.altitudeDeg);
    var up = clamp(Math.sin(toRad(Math.max(sun.altitudeDeg, 0))), 0, 1); // 高度による減衰
    var nightK = sun.altitudeDeg < -2 ? clamp((sun.altitudeDeg + 12) / 10, 0.05, 1) : 1;

    // 雲→太陽減衰 (+大気の澄み)
    var cloudW = cw(couplings, 'cloud-sun');
    var sunI = (1 - cloud * 0.88 * Math.min(cloudW, 1.5)) * (0.35 + 0.65 * clarity);
    sunI *= (0.15 + 0.85 * up);
    sunI *= nightK;
    sunI = clamp(sunI, 0, 1.6);

    // 風→波高 (m/s → 振幅)。無風でも微小なうねりを残す
    var windW = cw(couplings, 'wind-wave');
    var waveAmp = (0.12 + Math.pow(Math.min(wind, 20) / 20, 0.8) * 1.5) * (0.25 + 0.75 * Math.min(windW, 1.5));
    if (windW === 0) waveAmp = 0.35; // 結合OFF時は固定の中波

    // FIX-2: 風スペクトルの駆動パラメータ。Phillipsピーク尺度 L = V^2 / g (m)
    // を場面単位に写像し、シェーダ側で Phillips型 A(k)∝exp(-1/(kL)^2)/k^4 の
    // 重み・分散関係 omega=sqrt(g·k)・量子化位相を組み立てる。
    // 風が強いほど大きな波が立ち、弱いほど細かいさざめきだけが残る。
    var Lphys = (wind * wind) / 9.8;
    var specL = clamp(0.30 + Lphys * 0.30, 0.30, 4.0);
    var specAmp = (0.10 + Math.pow(Math.min(wind, 20) / 20, 0.8) * 1.30) * (0.25 + 0.75 * Math.min(windW, 1.5));
    if (windW === 0) specAmp = 0.30;

    // 太陽→コースティクス強度
    var scW = cw(couplings, 'sun-caustics');
    var causticsI = sunI * (0.35 + 0.65 * up) * 0.9 * Math.min(Math.max(scW, 0), 2);
    // 濁り→コースティクス減衰
    var tbW = cw(couplings, 'turbid-caustics');
    causticsI *= (1 - turbid * 0.75 * Math.min(tbW, 1.5));
    causticsI = clamp(causticsI, 0, 2.5);

    // 雲→グリッター拡散
    var glW = cw(couplings, 'cloud-glitter');
    var glitterI = (0.4 + sunI * 1.6) * nightK;
    var diffuse = clamp(cloud * Math.min(Math.max(glW, 0), 1.5), 0, 1);
    if (glW === 0) diffuse = cloud * 0.5;

    // 水深→吸収ティント。水深が深いほど効く。透明度が高いと弱まる。
    var dtW = cw(couplings, 'depth-tint');
    var depthShade = clamp(Math.log10(1 + depth) / Math.log10(31), 0, 1); // 0.1m→~0, 30m→1
    depthShade *= (0.25 + 0.75 * (1 - transp * 0.6)) * Math.min(Math.max(dtW, 0), 2);
    if (dtW === 0) depthShade = 0.25;
    depthShade = clamp(depthShade, 0, 1);
    // FIX-5: 吸収係数は純水開始値 β=(R 0.34, G 0.057, B 0.0092)/m とし、
    // 濁り・色素・透明度ノードでスケールする現行構造を守る。
    // シェーダ側で路程 (0.25+水深m)×(1+濁り×1.5)×0.22 を掛けて
    // trans=exp(-absorb×路程) とする。水の青緑 + 色素(藻)の緑を混ぜる。
    var absorbScale = (0.55 + 0.45 * (1 - transp)) * (1.0 + turbid * 1.2);
    var absorb = [
      0.34 * absorbScale + pigment * 0.10,
      0.057 * absorbScale + pigment * 0.12,
      0.0092 * absorbScale + pigment * 0.03 + (1 - transp) * 0.02
    ];
    var scatter = [0.03 + pigment * 0.10, 0.06 + pigment * 0.12, 0.05 + pigment * 0.03];
    // FIX-6: HDR合成の振幅 kc(既定0.8、プリセット側で0.5〜1.1の範囲で調整可)と
    // 焦点深度(既定0.8m、レシピの水深uDepthMに連動してシェーダ側で係数化)。
    // プリセット0(真夏の沖縄・正午)のみ面のコントラストで強さを出すため kc=1.0。
    // 他プリセットは据え置き0.8。
    var isPreset0 = state.name === '真夏の沖縄・正午';
    var causticsKc = isPreset0 ? 1.0 : 0.8;
    var focusM = 0.8;
    // FIX-6: 伝播ゲイン dgL の圧縮上限。既定4.20、プリセット0のみ3.40(粒の密度では
    // なく面のコントラストで強さを出す。kc=1.0との組合せ)。シェーダ側は uDgCap。
    var dgCap = isPreset0 ? 3.40 : 4.20;

    // 屈折量: IOR 1.333相当を基準に波高で微変調
    var refract = 0.028 * (0.5 + waveAmp * 0.7);

    // 流れ: 風向に流す。*8秒が整数格子になるよう量子化(継ぎ目なし用)
    var flowSpeed = 0.05 + flowSpd * 0.5 + wind * 0.008;
    var rawFx = Math.sin(windDir) * flowSpeed * 0.35;
    var rawFy = -Math.cos(windDir) * flowSpeed * 0.35;
    var flowPerSec = [Math.round(rawFx * 8) / 8, Math.round(rawFy * 8) / 8];

    // 太陽方向ベクトル(画面=天頂見下ろし近似: zが高度)
    var altR = toRad(sun.altitudeDeg);
    var azR = toRad(sun.azimuthDeg);
    var sunDir = [Math.sin(azR) * Math.cos(altR), Math.cos(azR) * Math.cos(altR), Math.sin(altR)];
    if (sunDir[2] < 0.02) sunDir[2] = 0.02;
    var sl = Math.sqrt(sunDir[0] * sunDir[0] + sunDir[1] * sunDir[1] + sunDir[2] * sunDir[2]);
    sunDir = [sunDir[0] / sl, sunDir[1] / sl, sunDir[2] / sl];

    return {
      sun: sun,
      sunDir: sunDir,
      sunColor: sunColor,
      sunI: sunI,
      diffuse: diffuse,
      waveAmp: waveAmp,
      waveFreq: 5.0 + fineness * 14.0,
      waveDir: [Math.sin(windDir), -Math.cos(windDir)],
      // FIX-2: シミュレーションの駆動パラメータ(スペクトルの L・振幅・水深・細かさ)
      specL: specL,
      specAmp: specAmp,
      depthM: depth,
      fine: fineness,
      flowPerSec: flowPerSec,
      refract: refract,
      chroma: 0.10 + (1 - clarity) * 0.10,
      causticsI: causticsI,
      causticsKc: causticsKc,
      dgCap: dgCap,
      focusM: focusM,
      glitterI: clamp(glitterI, 0, 2.5),
      absorb: absorb,
      scatter: scatter,
      depthShade: depthShade,
      turbid: turbid,
      waveHeightM: waveAmp * 0.12
    };
  }

  // FIX-2: 水の状態(慣性つき)。風の入力が変わっても波のスペクトルは
  // 即座に切り替わらず、時定数 TAU_W(約2.5秒)で追従する。シェーダが毎
  // フレーム評価する「状態」がこれであり、単なる時刻関数ではないことの核。
  // target = { specL, specAmp, dirX, dirY, flowX, flowY, fine }。
  var TAU_W = 2.5;
  function createWaterState(initial) {
    var cur = {
      specL: 1.0, specAmp: 0.5, dirX: 0, dirY: -1, flowX: 0, flowY: 0, fine: 0.5
    };
    if (initial) snap(initial);
    function snap(t) {
      if (t.specL != null) cur.specL = t.specL;
      if (t.specAmp != null) cur.specAmp = t.specAmp;
      if (t.dirX != null && t.dirY != null) {
        var l = Math.sqrt(t.dirX * t.dirX + t.dirY * t.dirY) || 1;
        cur.dirX = t.dirX / l; cur.dirY = t.dirY / l;
      }
      if (t.flowX != null) cur.flowX = t.flowX;
      if (t.flowY != null) cur.flowY = t.flowY;
      if (t.fine != null) cur.fine = t.fine;
    }
    function update(t, dt) {
      var h = Math.max(0, Math.min(dt || 0, 1));
      var k = 1 - Math.exp(-h / TAU_W); // 0に近いほど旧状態を保つ(慣性)
      if (t.specL != null) cur.specL += (t.specL - cur.specL) * k;
      if (t.specAmp != null) cur.specAmp += (t.specAmp - cur.specAmp) * k;
      if (t.dirX != null && t.dirY != null) {
        var l = Math.sqrt(t.dirX * t.dirX + t.dirY * t.dirY) || 1;
        cur.dirX += (t.dirX / l - cur.dirX) * k;
        cur.dirY += (t.dirY / l - cur.dirY) * k;
        var cl = Math.sqrt(cur.dirX * cur.dirX + cur.dirY * cur.dirY) || 1;
        cur.dirX /= cl; cur.dirY /= cl;
      }
      if (t.flowX != null) cur.flowX += (t.flowX - cur.flowX) * k;
      if (t.flowY != null) cur.flowY += (t.flowY - cur.flowY) * k;
      if (t.fine != null) cur.fine += (t.fine - cur.fine) * k;
      return { specL: cur.specL, specAmp: cur.specAmp,
        dirX: cur.dirX, dirY: cur.dirY,
        flowX: cur.flowX, flowY: cur.flowY, fine: cur.fine };
    }
    return { update: update, snap: snap, TAU_W: TAU_W };
  }

  // 実在モード: 地名→緯度経度→太陽計算+気象取得
  function fetchRealEnvironment(placeName, dateStr, timeStr) {
    var geoUrl = 'https://geocoding-api.open-meteo.com/v1/search?name=' +
      encodeURIComponent(placeName) + '&count=1&language=ja&format=json';
    return fetch(geoUrl).then(function (r) {
      if (!r.ok) throw new Error('geocode HTTP ' + r.status);
      return r.json();
    }).then(function (gj) {
      if (!gj.results || !gj.results.length) throw new Error('地名が見つかりません');
      var g = gj.results[0];
      var lat = g.latitude, lon = g.longitude;
      var label = g.name + (g.country ? ' (' + g.country + ')' : '');
      var wUrl = 'https://api.open-meteo.com/v1/forecast?latitude=' + lat +
        '&longitude=' + lon +
        '&hourly=temperature_2m,cloud_cover,wind_speed_10m,wind_direction_10m' +
        '&wind_speed_unit=ms&timezone=auto&forecast_days=16';
      // 日時指定が過去/未来でも拾えるよう一旦予報APIを叩く(失敗しても太陽計算は返す)
      return fetch(wUrl).then(function (r) {
        if (!r.ok) throw new Error('weather HTTP ' + r.status);
        return r.json();
      }).then(function (wj) {
        var out = { latitude: lat, longitude: lon, label: label, weather: null };
        try {
          var target = new Date(dateStr + 'T' + (timeStr || '12:00') + ':00');
          var times = (wj.hourly && wj.hourly.time) || [];
          var best = -1, bestDiff = Infinity;
          for (var i = 0; i < times.length; i++) {
            var d = Math.abs(new Date(times[i]).getTime() - target.getTime());
            if (d < bestDiff) { bestDiff = d; best = i; }
          }
          if (best >= 0 && bestDiff < 1000 * 3600 * 24 * 20) {
            out.weather = {
              cloudPct: wj.hourly.cloud_cover[best],
              windSpeedMS: wj.hourly.wind_speed_10m[best],
              windDirDeg: wj.hourly.wind_direction_10m[best],
              tempC: wj.hourly.temperature_2m[best]
            };
          }
        } catch (e) { out.weather = null; }
        return out;
      }).catch(function (werr) {
        // BRIEF指定: API失敗時は太陽計算だけでも成立させ、天気は手動値のまま
        return { latitude: lat, longitude: lon, label: label, weather: null, weatherError: String(werr && werr.message || werr) };
      });
    });
  }

  window.CausticsEnv = {
    COUPLINGS: COUPLINGS,
    FIX5_BANDS: { L: 12, M: 12, H: 12, micro: 8, total: 44 },
    FIX5_DELTA: 0.0833,
    FIX5_KC: 0.8,
    FIX5_FOCUS_M: 0.8,
    FIX6_NYQUIST: 0.65,
    FIX6_FINE_EPS_RATIO: 0.8,
    FIX6_HBAND_GAIN: 0.5,
    FIX6_LAYER_BROAD: 0.55,
    FIX6_LAYER_FINE: 0.45,
    FIX6_BROAD_POW: 1.25,
    FIX6_KC_PRESET0: 1.0,
    FIX6_DGCAP: 4.20,
    FIX6_DGCAP_PRESET0: 3.40,
    // FIX-7: 前方スプラッティングの解像度・正規化定数(水の状態の生成は不変。
    // 光輸送の評価だけを後方ヤコビアンから前方蓄積に替える)。
    // グリッド N×N 頂点(既定384: FIX-8で192→384に細分化)、
    // バッファは描画の1/2で上限1024(目安384は据え置き)。
    // FOOTPRINT_K はガウスフットプリント w=exp(-2.2r^2) の実効面積係数
    // (連続積分 πS^2/8.8≈0.357を起点に、離散化+端部損失ぶんをverify-nodeで較正して0.31)。
    // SUN_DISC は FIX-6 の太陽ぼかし半径式(0.0093×水深)の流用。
    FIX7_SPLAT_GRID: 384,
    FIX7_CAUSTICS_SIZE: 384,
    FIX7_CAUSTICS_MAX: 1024,
    FIX7_FOOTPRINT_K: 0.31,
    FIX7_SUN_DISC: 0.0093,
    FIX7_POINT_MAX: 64,
    FIX7_FORWARD: true,
    // FIX-9: 前方パスの点描画(GL_POINTS)を三角形メッシュに置換。
    // (a) MESH_GRID: 水面グリッド解像度。既定256(厳しければ192可)。
    // (b) DISP_K: 着地点変位スケールの物理較正(単一箇所・全プリセット共通、
    // プリセット別調整・露出調整なし)。波座標→着地点写像の換算が過小だった
    // ぶんを補う大域係数。後方フォールバック分岐・像の揺らぎは従来式のまま。
    // (c) VCAP_FWD: 前方メッシュ用のピーク正則化上限(単一箇所)。
    // メッシュ離散化+折り畳み重なりのエネルギー較正で verify-node により決定
    // (FIX-7 の FOOTPRINT_K 較正と同等の位置づけ)。後方の δ(0.0833)は不変。
    FIX9_MESH_GRID: 256,
    FIX9_DISP_K: 5.0,
    FIX9_VCAP_FWD: 5.0,
    // FIX-10: 帯域別変位ゲイン(単一箇所・全プリセット共通、プリセット別/露出調整なし)。
    // 着地点変位に使う勾配のみ L/M/H 帯域別に重みづけし、変位は主に長波長で作る。
    // L=1.0/M=0.55/H=0.25 が初期値(微細波はH扱いでGHを掛ける)。法線・反射など
    // 変位以外の勾配用途は現行どおり全帯域。DISP_K は単一較正を維持し、
    // 1.5〜3.5 を優先掃引したが全プリセットの変位std(0.02〜0.15)を満たさず、
    // 上限を広げて 5.0 に再較正した(理由は RESULT.md FIX-10節)。
    FIX10_DISP_GL: 1.0,
    FIX10_DISP_GM: 0.55,
    FIX10_DISP_GH: 0.25,
    FIX10_DISP_K: 5.0,
    // FIX-11: DISP_K の単一再較正。輝度が頂点解析密度の補間からフラグメント微分
    // による幾何面積比へ変わったため K=1.5〜5.0 を新評価で掃引し、直接指標セット
    // (平均1±0.15・最小≤0.5・max6〜12・相関≥0.5かつFIX-5水準以上・CV≥0.25)を
    // 最もよく満たす K=2.0 を採用(理由は RESULT.md FIX-11節)。新方式では折り畳み
    // 稜線が頂点位置に依らず線として明るくなるため、旧方式ほどの幾何学的折り畳み
    // (大K)は不要になる。変位 std は今回から診断値でありハードゲートではない。
    // FIX9/FIX10 の値は履歴として残す。δ^2(=uDelta2Fwd=1/VCAP^2)は現行定数と整合。
    FIX11_DISP_K: 2.0,
    solarPosition: solarPosition,
    sunColorFromAltitude: sunColorFromAltitude,
    resolveEnvironment: resolveEnvironment,
    createWaterState: createWaterState,
    TAU_W: TAU_W,
    fetchRealEnvironment: fetchRealEnvironment
  };
})();
