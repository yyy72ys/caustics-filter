# BRIEF-PJ4-TRIAL-ASSETS — PJ4 試し描き用アセット3点の切り出し

Workdir: D:\yohey2026\dsh\muse\caustics-filter\offline
Model: opencode-go/deepseek-v4.1-flash (this repo's standing rule)

## 目的

PJ4（物理絵画）の試し描き用に、PJ1 の本番レシピ（v19 系）で作った静止画1枚と、
それに画素対応する C（コースティクス照度マップ）と底アルベドを切り出す。
既存ファイルは一切変更しない。新規スクリプト1本と新規出力ディレクトリだけを作る。

## 手順

1. 新規スクリプト `run_pj4_trial.py` を offline/ に作る。
   - `run_v19_prod.py` を import し、`build_scene()` と `load_photo()` をそのまま使う
     （シーン条件・シード・光学は v19 本番と同一になる）。
   - フレーム i=96（t=4.0 s）を対象にする。`render_one` は dbg を捨てるので使わず、
     その中身（fields_at → det_frame → V16.render_frame_v16(..., return_debug=True)）を
     このスクリプト内に複製して1回だけ実行し、C・合成画像・dbg["alb"] を全部保持する。
   - 合成条件は render_one と完全同一（A_FROZEN / E_FROZEN / B_SIGMA_FROZEN /
     substrate_gain (1,1,1) / g1 は render_one と同じ求め方）。
2. 出力先 `out/pj4-trial/` に以下を保存する。
   - `still.png` — 合成済み静止画 (736x540, uint8)
   - `C-master.npy` — det_frame の返す C を無加工で float32 保存 (2048x2048)
   - `C-frame.npy` — C を出力画素グリッドへサンプリングした float32 (540,736)。
     サンプリングは render_frame_v2 の master ratio サンプルと同式
     （u=(ix+0.5)/736, v=(iy+0.5)/540, texel = u*2048-0.5 / v*2048-0.5, 周期 wrap の
     バイリニア）ただし屈折シフト（shx/shy）は入れない。この定義はレポートに明記する。
   - `albedo.png` — dbg["alb"]（線形 RGB float）を linear→sRGB 変換し clip して uint8 化。
     image.py の linear_to_srgb を使う。
   - `numbers-pj4-trial.json` — 各ファイルの shape/dtype、C-master と C-frame の
     mean/min/max/p50/p99、still の sha256、still と既存
     `out/v19/frames/frame-0096.png` の画素比較（max abs diff、一致画素の割合）、
     albedo の band 別 mean、所要時間。
   - `REPORT-PJ4-TRIAL.md` — 上記の要約と、各ファイルの定義（特に C-frame の
     サンプリング定義とシフトなしの断り）を書く。
3. 検証（スクリプト内で assert する）:
   - still が frame-0096.png と完全一致すること。一致しなければ numbers に
     差分統計を記録し、レポートに大きく明記する（黙って成功にしない）。
   - C-master の mean が v19 numbers（numbers-v19-prod.json のフレーム96の C_mean）と
     1e-6 以内で一致すること。
   - C-frame の shape が still の (540,736) と一致すること。
4. 冪等性: 出力一式が揃っていれば再計算せず skip とログする。

## 禁止事項

- 既存ファイルの変更・削除（run_v19_prod.py 含む）。
- 数値の創作。測れないものは測れないと書く。
- out/v19/ 配下への書き込み。
