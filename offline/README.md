# 高忠実度オフラインパイプライン v1

WebGLの枠を外した物理経路の最終レンダリング。`caustics-recipe/v1` JSONを入力に、
静止画1枚＋16フレームのテストクリップ＋数値レポートを出す。

## 構成

- `env.py` — `web/env.js` のポート（太陽天文計算・高度別光色・結合解決）。
  オフラインとWeb定義の整合用。
- `state.py` — 水の状態。Phillipsスペクトルを2048×2048グリッドでFFT生成
  （パッチ S=4·specL、kmin〜Nyquistで帯域3桁、方向分布cos型、微小波減衰長
  l=0.005m）。Choppy水平変位 λ=0.85。角周波数 ω=√(gk) を 2π/8s の整数倍に
  量子化して8秒厳密周期化。位相はレシピ由来シードで決定的に固定。
  高度・勾配・変位は周期双三次（Catmull-Rom）補間で任意点参照。
- `transport.py` — 光の輸送。毎フレーム400万フォトンをパッチ上へ一様
  ジッターサンプリング。入射方向＝太陽方向＋太陽円盤ジッター。
  水面でSnell屈折（IOR帯域別 R 1.33132 / G 1.33434 / B 1.33925、
  フォトンはRGB 3サブレイ）。底面交点を双線形スプラットで 2048×2048 HDR
  バッファへ蓄積（周期ラップでエネルギー保存）、σ0.75pxガウシアン平滑、
  物理正規化（平均＝解析的期待値）。
- `image.py` — 像の形成。同じ状態で写真の揺らぎサンプル＋Beer-Lambert減衰
  （路程 (0.25+水深)(1+濁り1.5)×0.22、帯域別）＋水中散乱アンビエント＋
  Fresnel空反射（F0=0.02）。ACES→sRGB。
- `run.py` — 編成。レシピJSONパスとプリセット番号を引数に取る。

## 実行

```bat
cd caustics-filter\offline
python run.py --recipe ../blender/recipe-1.json --preset 1
```

出力（`out/`）：

- `preset1-still.png`（735×539、t=0フレーム）
- `preset1-test16.mp4`（736×540・20fps・16フレーム、t=0〜0.75sを0.05s刻み）
- `preset1-numbers.json`（状態・輸送・場統計・ref1統計・ループ検証・実行時間）

依存: Python 3.10 + NumPy + Pillow + opencv-python（mp4書出）。
スモークテスト: `python run.py --recipe ../blender/recipe-1.json --preset 1
--grid 256 --photons 20000 --frames 2 --out out-smoke`

## パイプライン v2（主経路転換・OUTPUT-SPEC/PIPELINE-V2準拠、プリセット1/t=0）

- `det.py` — 決定論的forward輸送（主経路）。1頂点1レイ×Q方向、
  Snell屈折（代表帯域G）＋双線形スプラット＋σ0.75px FFTガウシアン
  （フォトン経路と同一ピクセルフィルタ）。出力は無次元倍率
  C = E_direct / E_base,dir（解析的除数のみ、無クランプ・無正規化）。
  `flow_backward_master` で逆変位（G帯、マスター格子px、2ch）を返す。
- `image.py render_frame_v2` — §3分離式（環境光はCで変調しない）。
- `run_v2.py` — 編成＋ゲート試験。`diag_v2.py` — ゲート分析専用。

```bat
cd caustics-filter\offline
python run_v2.py --recipe ../blender/recipe-1.json --preset 1
```

出力（`out/v2/`）：`C-t0.npy`・`flow-t0.npy`・`metadata.json`・
`numbers-v2.json`（＋`numbers-v2-diag.json`）・`preset1-still-v2.png`。
スモーク: `--grid 256 --photons 20000 --qlevels 5,7 --out out/v2-smoke
--no-baseline`。SciPy必須（v2のMS-SSIM・連結成分・面積平均用）。

## パイプライン v3（合成調整ラウンド・PIPELINE-V3準拠、プリセット1/t=0）

- 状態生成と輸送は凍結（`out/v2/C-t0.npy` を再利用し det/transport は
  呼ばない）。調整は `image.py render_frame_v3` の合成段のみ。
- `run_v3.py` — 編成・6バリアント計測。A=×1.0/×1.7/×2.6（E_amb と
  カメラ側加算項のみ、直接光項不変）× B=full/low（状態σ3.0px平滑）。
  線形平均・潰れ率・飛び率・破砕メトリクス（Laplacian比）を記録。

```bat
cd caustics-filter\offline
python run_v3.py --recipe ../blender/recipe-1.json --preset 1
```

出力（`out/v3/`）：`still-A{A}-B{B}.png` ×6・`numbers-v3.json`。
t=0静止画のみ（全フレーム本番は範囲外）。

## パイプライン v5 チャンク1（水と底の変数・PIPELINE-V5準拠、プリセット1/t=0）

追加ノードは report.md §1・§2・§5 の式に基づく導出量のみ（気分スライダーは無し）。
`state.py` / `transport.py` / `det.py` / `image.py` / `run_v2,v3,prod` は無変更。

- `v5_physics.py` — 物理ノード。Pope & Fry 吸収表（≈650/550/470nm 補間）、
  b_w(550)=0.0018・λ^-4.32、濁度 Δa+Δb、Beer–Lambert（下り s↓=z/cosθ_w、
  戻り s↑=z/cosθ_view）、Smith G₁、Cox–Munk σ²/異方性、基質アルベド、散乱ブラー。
- `v5_state.py` — `V5WaterState`（`state.WaterState` のサブクラス、H0 を再スケール
  するだけ）。H_s=4√m0 の読み出し／高さ駆動モード（Stokes H/λ≤1/7 キャップ付き）、
  実現 σ² の読み出しと高波数側の Cox–Munk 較正 g(k)=1/(1+(k/kc)^n)。
- `v5_composite.py` — §3 分離式を維持した合成。直接光項に Smith G₁ と基質ゲイン、
  物理 Beer–Lambert（下り・戻り）、C への濁度散乱ブラー。A/B は v3 最終の定数、
  C 層露出 E は表示レベル正規化（下記）。
- `recipe_v5.py` — `caustics-recipe/v2` の追加フィールド解析。v1 は新フィールド無しで
  従来挙動（基質ゲイン恒等・風駆動）。
- `run_v5.py` — 編成＋数値ゲート＋検証スチル4組。
- `diag_v5_wind.py` — 診断（検証スチル外）。風速を上げるとコースティクスが戻ることの確認。

```bat
cd caustics-filter\offline
python run_v5.py --recipe ../blender/recipe-1.json --preset 1
```

出力（`out/v5/`）：`numbers-v5.json`、`still-*.png` 8枚、`pair1..pair4-*.png`、
`C-baseline-v5.npy` / `C-depth-*.npy` / `C-wave-fixed.npy`、`recipes/*-v2.json`。
診断: `python diag_v5_wind.py --winds 5,8` → `diag-wind*.png` / `diag-wind.json`。

### v5 の近似・較正（申告）

- **Smith G₁ は統計近似**（Beckmann 斜面分布の Walter et al. 2007 有理近似）。
  レイ厳密版（フットプリント再交差テスト）は実装していない。両者は併用しない
  （report §1.3 / §7）。本プリセットでは α=0.1085, θ=31.0° で a=15.3≥1.6 → G₁=1
  （影は low-sun・steep-wave の効果であり、この条件では出ない）。
- **グリッター**: 異方性 (σ_a², σ_c²) のサンプラは実装・記録するが、オフライン v3
  合成にはグリッター項が無い（Fresnel 空反射のみ）ため検証スチルには適用しない。
  等方版は σ_a=σ_c の特殊ケースとして記録。
- **導出でない較正定数**（numbers-v5.json の `calibrations_not_derived` に列挙）:
  濁度→Δa+Δb のスケール 0.30 /m（散乱分率0.7）、散乱ブラー比例係数 40 px@2048²、
  b_w(550)=0.0018 /m のアンカー、表示レベル正規化 E_v5=E_v3·(L_v3/L_v5)。
- **表示レベル正規化は物理ではない**。物理 Beer–Lambert は v3 の ad-hoc 水深減衰より
  暗いため、平坦水面の直接光輝度が v3 最終と一致するよう C 層露出 E を一度だけ
  解き直している（E=2.2→2.708）。E=2.2 固定の診断値も numbers に併記。
- **既定プリセットの帰結**: U=1.5m/s で Cox–Munk に合わせるとスペクトル σ² は
  0.136→0.0107 rad²（実現値 0.0118、rms 勾配 20.3°→6.2°）、m0 の65%が除かれ
  H_s 0.148→0.084m（実現 0.089m）。v3 の状態が 1.5m/s 風としては過度に粗かった
  ことの修正であり、既定プリセット（z=0.6m）のコースティクスは大幅に穏やかに
  なる（C p99 4.28→1.25）。ただし同じ水面でも z=3m ではレンズ幾何で C p99 6.03
  まで戻る。風を上げても戻る（U=8 で C p99 2.59, max 16.4）。強いコースティクスが
  要る場合は風ノードか水深を上げる（どちらも既存のシーン入力）。

## パイプライン v6 チャンク2（雨・無風の凪・PIPELINE-V6準拠、プリセット1/t=0）

追加ノードは report.md §3（雨）・§5 の式に基づく導出量のみ。新しいシーン入力は
`rain.rate_mm_h` ただ1つ（既定 0）。`state.py` / `transport.py` / `det.py` /
`image.py` / `v5_*` / `recipe_v5.py` / `run_v2,v3,v5` は無変更。

- `v6_physics.py` — 雨の統計鎖。Marshall–Palmer Λ(R)=4.1R^-0.21、Atlas–Ulbrich
  v(D)=9.65−10.3exp(−0.6D)、衝突頻度 F(R)=∫Nv dD（3変種: 物理的クランプ／生の
  適合／厳密0.5–5mm域）、流量加重平均衝突エネルギー、D₀=3.67/Λ、輪波数
  k_ring(R)（D₀→単調減少、観測バンド1–2cmにアンカー）、重力–毛管分散
  ω²=gk+(σ/ρ)k³、輪バンド分散 σ_ring²=η·F·Ē·τ/(ρg+σk²)、第二効果（低k減衰）。
  較正定数は **η（RING_ETA=0.10）ただ1つ**。
- `v6_state.py` — `V6WaterState`（`V5WaterState` のサブクラス）。凪の
  Cox–Munk 下限 σ²≈0.003 の上に、等方な輪バンドを分散ターゲット一致で付加。
  ω は雨ありのとき重力–毛管分散（2π/8s 倍数へ量子化、8秒ループ厳密閉止）。
  R=0 は早期リターンで v5 とビット一致。
- `v6_composite` は新設せず `v5_composite.render_frame_v5` をそのまま使う
  （このチャンクは合成段を変えない）。
- `recipe_v6.py` — `rain.rate_mm_h` の解析（v5 解析の上に追加、v1 は従来どおり）。
- `run_v6.py` — 編成＋数値ゲート＋検証スチル（R=0/2/10/50）＋並置パネル＋
  拡大クロップ。風 0 の凪シーンを明示的に構成する。

```bat
cd caustics-filter\offline
python run_v6.py --recipe ../blender/recipe-1.json --preset 1
```

出力（`out/v6/`）：`numbers-v6.json`、`still-R00/R02/R10/R50.png`、
`panel-R00-R02-R10-R50.png`（2×2）、`crop-zoom-R00-R02-R10-R50.png`
（2行×4列：上段=納品スチルの3倍拡大、下段=同じ物理窓の C マスター（輸送出力）
を1:1・固定階調で。輪の粒状感がサンプリング由来でないことの証拠）、
`still-v5-calm.png`（R=0 の参照）、`C-rain-R*.npy`、
`C-baseline-v5-calm.npy`、`recipes/rain-R*-v2.json`。
`--crop-only` で輸送を再計算せずクロップ/パネルのみ再生成できる。
スモーク: `--grid 512 --q 9 --out out/v6-smoke`（約30秒）。

### v6 の近似・較正（申告）

- **輪の置き方 k_ring(R)**: report §3.4 は観測バンド（λ=1–2cm）と「大きい滴ほど
  長い輪」のみを与え、周波数の閉形式は無い。よって
  λ_ring = 10mm + 10mm·(D₀−D₀(R=1))/(D₀(R=50)−D₀(R=1)) という**単調配置則**
  （両端は観測バンドの限界、補間変数は質量中央径 D₀）を使う。R=10 で
  k=422 rad/m（設計の ≈420 と 0.5% 一致）、R∈[1,50] で 310–630 内、単調減少。
  これは較正定数ではなく配置則（numbers の
  `approximations_and_design_choices` に明記）。
- **輪バンド分散**は深水線形のエネルギー密度から σ²=E/(ρg+σk²)。バンド中心
  k_ring で評価（帯は狭いので中心近似）。**τ=150ms** は report §3.4 の実測
  プロセス時間で、粘性減衰はここに畳み込む（β は独立定数にしない）。
- **第二効果（低k減衰）**: 設計は「機構だけ・係数は1に近い近似」を要求。形は
  gain(k)=1/(1+γ·τ·w(k)), γ=F·Ē/(ρg·m0_lowk)、w=1/(1+(k/310)⁴) で**調整
  定数ゼロ**。凪シーンでは掛け先が Cox–Munk 下限しかないため実効 gain は
  R=2/10/50 で 0.998/0.987/0.926（＝ほぼ no-op）。8秒ループ分の値
  exp(−γ·8) は記録のみで適用しない。
- **バンドは k のガウス（σ_k=(630−310)/6）でハードクリップしない**。R=1/50 では
  中心がバンド端に載るため裾が 310/630 を越える（バンド内分散比は R=2 で 0.89、
  R=10 で 0.99、R=50 で 0.59）。ガウスの対称性から、バンドの**平均**波数は
  k_ring(R) に厳密一致するので、中心の単調性ゲートは厳密。
- **凪の背景**は既存パイプラインの cos² ローブ形状をそのまま残す（R=0 の
  ビット一致に必要）。等方性を持つのは雨の輪バンド側。
- **パッチ**: 風0では specL が 0.30 にクランプされパッチは 1.20m（設計の
  1.48m ではなく、輪波長1.5cm は ≈26px とより細かく解像）。
- **DSD の体積整合は ±20%**: 同じ Marshall–Palmer 分布で体積フラックス
  （(π/6)∫N D³v dD）を積分して雨量へ戻すと R=1/10/25/50 に対し
  1.18/11.6/28.2/54.2 mm/h（+18/+16/+13/+8%）。設計の「Marshall–Palmer は
  平均であって法則ではない（±factor 2）」の範囲内であり、**N₀ を再規格化
  して R を合わせる操作はしていない**（隠れた正規化を入れない）。numbers の
  `flux.rate_back_check_mm_h` に各値あり。
- **継承した v5 の較正**（濁度スケール・散乱ブラー比例・b_w(550)・表示レベル
  正規化）は不変更。

## 近似点（PIPELINE.mdに対する申告）

- 入射時の水面反射損失は無視し、全エネルギーを水中に保存する
  （平均保存の検証で確認）。
- 写真は底面アルベドとして扱う（撮影時の照明を含む近似）。
- フォトン側の法線参照も像形成側と同じ双三次補間・同じ状態。
- ωは深水近似 √(gk)（PIPELINE指定どおり）。
