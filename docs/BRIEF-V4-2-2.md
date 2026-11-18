# BRIEF-V4-2-2 — 床の作り直し＋監査の較正＋輸送の再確認

Workdir: `D:\yohey2026\dsh\muse\caustics-filter\offline`
Design: `DESIGN-V4-2-2.md`（本 brief と一緒に読むこと）
前作: `REPORT-V4-2-1.md`（out/v4-2-1/）の §1・§6 が切り分けの出発点。

## STAGE 0 — 絶対パス事前確認（1つでも欠けたら停止して報告、他ツリーを探さない）

- `D:\yohey2026\dsh\muse\caustics-filter\offline\bed_gen.py`
- `D:\yohey2026\dsh\muse\caustics-filter\offline\run_v4_2_1.py`
- `D:\yohey2026\dsh\muse\caustics-filter\offline\stream_surface2.py`
- `D:\yohey2026\dsh\muse\caustics-filter\offline\out\v4-2-1\REPORT-V4-2-1.md`
- 既存 .py 全本の sha256 を stage-0 記録として保存（G-LEGACY の基準）

## 実装

新規ファイルのみ: `bed_gen2.py` と `run_v4_2_2.py`。既存ファイルは一切変更しない。
場面定数（床パッチ 0.6 m、D50=0.020、σg=2.1、装甲 2.0、岩 0.15/0.22 m、
水面 H=0.25/U=0.6/S0=0.003/D84=0.042、プリセット5太陽、clear、η=1、E=2.2、
シード 3255045276382562994）は v4-2.1 から verbatim に取る。

生成の変更点（DESIGN-V4-2-2 の生成側 1–4）:
- 突出量 p の非対称分布（大半は埋没、一部が最大 ~0.5·D まで突出）
- σz(matrix) を実測反復（最大4回）で 0.40·D50_surface を中心に閉じる。全試行を記録
- 空隙率 0.30–0.45、クラスター被覆率 0.10 中心になるよう構成パラメータを置く

## ゲート

### G-BED2（監査は DESIGN の監査側 1–3 の定義で測る）
| 項目 | 帯 |
| --- | --- |
| 歪度 (matrix) | 0.3–0.9 |
| σz / D50_surface（正確な監査） | 0.32–0.46 |
| 構造関数 H（折れ点より大スケール側） | 0.7–0.85 |
| D84/D50（正確な監査） | 1.9–2.3 |
| 装甲比（正確な監査） | 1.7–2.3 |
| 傾斜角 中央値（上流向き） | 25–35° |
| クラスター被覆率 | 0.05–0.15 |
| 空隙率 | 0.30–0.45 |
併記（合否なし）: 復元演算子の生値と参照床での偏り、折れ点の位置（×D50）、
小スケール側の H、σz/入力D50、ファセット傾斜の分布（>60° / >80° の割合）。

### G-TRANS（新しい床で v4-2.1 の輸送を再実行）
- (a) エネルギー保存: 平坦 vs 高さ場の直接フラックス差 ≤ 5%
- (b) 岩の影: コア/照射 < 0.15

### G-LEGACY / G-OUT
- 既存 py の sha256 が stage-0 記録と一致
- 動画 48 frames / 2.000 s / 736×540 / 24 fps（ffprobe）

## 成果物（全て out/v4-2-2/）
- REPORT-V4-2-2.md（σz 反復の全試行と、v4-2.1 との項目別対比表を含む）
- numbers-v4-2-2.json、bed-stats.json
- **heightfield.npy**（統計グリッドの高さ場、float32、単位 m）
- **transects.png**（流れ方向・横断方向の断面図。mm 単位、D50 スケールバー付き）
- depth-map.png、frame-0024.png、piece-2s-v4-2-2.mp4
- checkpoint.json、stage0-py-hashes.json

## 実行規律
- 全ステージ foreground。jobkit checkpoint を `out/v4-2-2/checkpoint.json` に
  stage ごとに原子的にコミット
- 判定は数値統計のみ（目視なし）。ゲートの帯を後から動かすのは禁止 —
  外れたら外れたまま、理由の切り分けを付けて報告する
