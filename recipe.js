/* レシピ(JSON schema v1)+プリセット+保存/読込 (classic script)
 * window.CausticsRecipe = { SCHEMA, defaultState, PRESETS, serialize, parse, ... }
 */
(function () {
  'use strict';
  var SCHEMA = 'caustics-recipe/v1';
  var LS_PREFIX = 'caustics-recipe:';

  function clone(o) { return JSON.parse(JSON.stringify(o)); }

  function defaultCouplings() {
    return clone(window.CausticsEnv.COUPLINGS);
  }

  function defaultState() {
    return {
      schema: SCHEMA,
      name: '無題の環境',
      place: null, // { label, latitude, longitude, datetime }
      sun: { date: '2026-08-15', time: '12:00', latitudeDeg: 26.3, longitudeDeg: 127.8 },
      water: { depthM: 2.0, transparency: 0.7, turbidity: 0.3, pigment: 0.3 },
      windWave: { windSpeedMS: 3.0, windDirDeg: 135, fineness: 0.5, flowSpeed: 0.5 },
      weather: { cloudPct: 20, clarity: 0.8 },
      couplings: defaultCouplings(),
      render: { strength: 1.0, exposure: 1.0 }
    };
  }

  // refsの質感に寄せた5プリセット
  function presets() {
    var base = defaultState;
    var p1 = base();
    p1.name = '真夏の沖縄・正午';
    p1.sun = { date: '2026-08-15', time: '12:00', latitudeDeg: 26.3, longitudeDeg: 127.8 };
    p1.water = { depthM: 8.0, transparency: 0.9, turbidity: 0.08, pigment: 0.12 };
    p1.windWave = { windSpeedMS: 4.0, windDirDeg: 135, fineness: 0.55, flowSpeed: 0.55 };
    p1.weather = { cloudPct: 5, clarity: 0.95 };
    p1.render = { strength: 1.0, exposure: 0.9 }; // FIX-3(欠陥3): 最強プリセットの白飛び抑制

    var p2 = base();
    p2.name = '緑の川の浅瀬・午前';
    p2.sun = { date: '2026-06-10', time: '09:30', latitudeDeg: 35.7, longitudeDeg: 139.7 };
    p2.water = { depthM: 0.6, transparency: 0.55, turbidity: 0.35, pigment: 0.8 };
    p2.windWave = { windSpeedMS: 1.5, windDirDeg: 90, fineness: 0.7, flowSpeed: 0.6 };
    p2.weather = { cloudPct: 20, clarity: 0.75 };
    p2.render = { strength: 1.1, exposure: 1.0 };

    var p3 = base();
    p3.name = '夕方の湖';
    p3.sun = { date: '2026-10-05', time: '17:20', latitudeDeg: 35.7, longitudeDeg: 139.7 };
    p3.water = { depthM: 3.0, transparency: 0.6, turbidity: 0.25, pigment: 0.35 };
    p3.windWave = { windSpeedMS: 1.0, windDirDeg: 270, fineness: 0.4, flowSpeed: 0.3 };
    p3.weather = { cloudPct: 15, clarity: 0.7 };
    p3.render = { strength: 1.0, exposure: 1.0 };

    var p4 = base();
    p4.name = '曇りの日';
    p4.sun = { date: '2026-09-20', time: '11:00', latitudeDeg: 35.7, longitudeDeg: 139.7 };
    p4.water = { depthM: 2.0, transparency: 0.4, turbidity: 0.55, pigment: 0.4 };
    p4.windWave = { windSpeedMS: 5.0, windDirDeg: 0, fineness: 0.5, flowSpeed: 0.5 };
    p4.weather = { cloudPct: 92, clarity: 0.35 };
    p4.render = { strength: 0.9, exposure: 0.95 };

    var p5 = base();
    p5.name = '深い海';
    p5.sun = { date: '2026-07-20', time: '10:00', latitudeDeg: 26.3, longitudeDeg: 127.8 };
    p5.water = { depthM: 20.0, transparency: 0.9, turbidity: 0.1, pigment: 0.15 };
    p5.windWave = { windSpeedMS: 3.5, windDirDeg: 180, fineness: 0.45, flowSpeed: 0.45 };
    p5.weather = { cloudPct: 10, clarity: 0.9 };
    p5.render = { strength: 1.0, exposure: 1.0 };

    return [p1, p2, p3, p4, p5];
  }

  function serialize(state) {
    var s = clone(state);
    s.schema = SCHEMA;
    return JSON.stringify(s, null, 2);
  }

  function num(x, dflt, lo, hi) {
    var v = Number(x);
    if (!isFinite(v)) v = dflt;
    return Math.max(lo, Math.min(hi, v));
  }

  function parse(jsonStr) {
    var o = JSON.parse(jsonStr);
    if (!o || o.schema !== SCHEMA) throw new Error('schemaが ' + SCHEMA + ' ではありません');
    var s = defaultState();
    s.name = String(o.name || '無題の環境').slice(0, 80);
    if (o.place) s.place = o.place;
    if (o.sun) {
      s.sun.date = String(o.sun.date || s.sun.date);
      s.sun.time = String(o.sun.time || s.sun.time);
      s.sun.latitudeDeg = num(o.sun.latitudeDeg, 35.7, -90, 90);
      s.sun.longitudeDeg = num(o.sun.longitudeDeg, 139.7, -180, 180);
    }
    if (o.water) {
      s.water.depthM = num(o.water.depthM, 2, 0.1, 60);
      s.water.transparency = num(o.water.transparency, 0.7, 0, 1);
      s.water.turbidity = num(o.water.turbidity, 0.3, 0, 1);
      s.water.pigment = num(o.water.pigment, 0.3, 0, 1);
    }
    if (o.windWave) {
      s.windWave.windSpeedMS = num(o.windWave.windSpeedMS, 3, 0, 40);
      s.windWave.windDirDeg = num(o.windWave.windDirDeg, 135, 0, 360);
      s.windWave.fineness = num(o.windWave.fineness, 0.5, 0, 1);
      s.windWave.flowSpeed = num(o.windWave.flowSpeed, 0.5, 0, 1);
    }
    if (o.weather) {
      s.weather.cloudPct = num(o.weather.cloudPct, 20, 0, 100);
      s.weather.clarity = num(o.weather.clarity, 0.8, 0, 1);
    }
    if (o.couplings && Array.isArray(o.couplings)) {
      var defs = defaultCouplings();
      s.couplings = defs.map(function (d) {
        for (var i = 0; i < o.couplings.length; i++) {
          if (o.couplings[i].id === d.id) {
            return {
              id: d.id, from: d.from, to: d.to, label: d.label,
              enabled: !!o.couplings[i].enabled,
              weight: num(o.couplings[i].weight, 1, 0, 2)
            };
          }
        }
        return d;
      });
    }
    if (o.render) {
      s.render.strength = num(o.render.strength, 1, 0, 2);
      s.render.exposure = num(o.render.exposure, 1, 0.2, 2.5);
    }
    return s;
  }

  // localStorage 一覧保存
  function lsList() {
    var out = [];
    try {
      for (var i = 0; i < localStorage.length; i++) {
        var k = localStorage.key(i);
        if (k && k.indexOf(LS_PREFIX) === 0) out.push(k.slice(LS_PREFIX.length));
      }
    } catch (e) { /* プライベートモード等では使えない */ }
    return out.sort();
  }
  function lsSave(name, state) {
    var s = clone(state); s.name = name;
    localStorage.setItem(LS_PREFIX + name, serialize(s));
  }
  function lsLoad(name) { return parse(localStorage.getItem(LS_PREFIX + name)); }
  function lsDelete(name) { localStorage.removeItem(LS_PREFIX + name); }

  window.CausticsRecipe = {
    SCHEMA: SCHEMA,
    defaultState: defaultState,
    presets: presets,
    serialize: serialize,
    parse: parse
,
    lsList: lsList,
    lsSave: lsSave,
    lsLoad: lsLoad,
    lsDelete: lsDelete
  };
})();
