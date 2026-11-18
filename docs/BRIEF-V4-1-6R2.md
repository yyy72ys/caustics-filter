# BRIEF-V4-1-6R2 — 検証メインの深さ 0.6 m で出す（分離形）

Workdir: `D:\yohey2026\dsh\muse\caustics-filter\offline`

## 背景

ユーザー決定（2026-11-13 05:01）: **検証のメインの深さは 0.6 m**
（v2 の揺らぎが出た深さ）。4 m（プリセット5）は時々使う。
v4-1.6R（深さ 0.35 m）はこの決定で差し替え。本作がその実施。

## STAGE 0 — 絶対パス事前確認（欠けたら停止して報告）

- `run_v4_1_2.py`、`stream_surface3.py`、`stream_surface4.py`、
  `run_v4_1_5.py`、`out\v4-1-5\REPORT-V4-1-5.md`
- 既存 .py 全本の sha256 を stage-0 記録として保存

## 作業

新規ファイルのみ: `run_v4_1_6r2.py`。既存ファイルは一切変更しない。

1. **波の場の生成**（川の水、深さ 0.6 m）: stream_surface3 の会計で、
   環境は自然シーン（run_v4_1_2.py のシーンB）を土台に **H=0.60 m**
   とする。他は verbatim: U=0.80 m/s, U0=0.92, S0=0.002,
   D84=0.09 m, seed=20261112, 媒体 algae_pond。
   生成パッチは L_gen = 10 × λ0（λ0 はこの場面の定常根から取る）。
   底応答の底は元写真を L_gen に写像。
   輸送段は v4-1.5 の折り返し規則（|ds/di| のラダーから帯域制限幅を
   その場面で導く）。
2. **ベース画像**: 元写真を L_img = **3.68 m**、736×540、
   屈折の水深 z = **0.60 m**。合成は凍結形 A=2.6 / Blow σ=3 /
   E=2.2 / Q=25。
3. **引き伸ばし合成**: 生成した場を L_gen → L_img へ写す。
   倍率と描かれる波長（λ0 × 倍率）はレポート冒頭と numbers に明記。

## ゲート

- G-LINE: frame 24 の軸平行周期（5–15 px 帯）が元写真以下
  （v4-1.5 と同じ測定法・同じ表形式）
- G-SPEC: 生成場の実測スペクトルピーク（annulus）が解析 k0 の
  0.75–1.33 倍以内
- G-MOVE: ユニークフレーム 48/48、フレーム間平均絶対差（緑） > 0.5 諧調
- G-WARP: frame 24 の歪み（緑・平均）と C の標準偏差を記録（合否なし）
- G-LEGACY: 既存 py の sha256 が stage-0 記録と一致
  （並行レーンの外部差分は帰属を分けて報告）
- G-OUT: ffprobe で 48 frames / 2.000 s / 736×540 / 24 fps

## 成果物（全て out/v4-1-6r2/）
REPORT-V4-1-6R2.md、numbers-v4-1-6r2.json、piece-2s-v4-1-6r2.mp4、
frame-0024.png、checkpoint.json、stage0-py-hashes.json

## 実行規律
- 全ステージ foreground。jobkit checkpoint を stage ごとに原子的にコミット
- 判定は数値統計のみ（目視なし）。後付け調整は禁止
