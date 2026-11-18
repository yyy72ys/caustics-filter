# BRIEF-V4-1-4 — 予算会計の修正と両シーン再レンダ

Workdir: `D:\yohey2026\dsh\muse\caustics-filter\offline`
Design: `DESIGN-V4-1-4.md`（本 brief と一緒に読むこと）

## STAGE 0 — 絶対パス事前確認（1つでも欠けたら停止して報告、他ツリーを探さない）

- `D:\yohey2026\dsh\muse\caustics-filter\offline\stream_surface2.py`
- `D:\yohey2026\dsh\muse\caustics-filter\offline\run_v4_1_2.py`
- `D:\yohey2026\dsh\muse\caustics-filter\offline\run_v4_1_3.py`
- `D:\yohey2026\dsh\muse\caustics-filter\offline\test-inputs` は無いので入力写真は
  `D:\mutaro-work\moth-on-bark.jpg`
- 既存 .py 全本の sha256 を stage-0 記録として保存（G4 の基準）

## 実装

新規ファイルのみ: `stream_surface3.py` と `run_v4_1_4.py`。既存ファイルは一切変更しない。

`stream_surface3.py`: `stream_surface2.py` を import して再利用し、予算配分だけを
DESIGN-V4-1-4 の会計に差し替える:

- `var(ii) = φ · σ_reg²`、`var(iii) = (1−φ) · σ_reg²`
- (i) 底応答は従来どおり計算し、合計分散は `σ_reg² + var_i`（二乗和）
- 単体確認: var_i を強制的に 0 にしたとき、(ii)(iii) の配分と生成面の σ が
  stream_surface2 と一致すること（数値で報告）

## シーン（定数は既存 run ファイルから verbatim に取ること。新しい値を発明しない）

- **水路検証**: v4-1.2 と同じ13条件の σ 比（numbers のみ）
- **シーンB（自然シーン）**: `run_v4_1_2.py` の場面定数をそのまま使い 2 秒
  （48 frames, 736×540）。ファイル名 `piece-2s-v4-1-4-natscene.mp4`、静止画 `frame-natscene-0024.png`
- **シーンC（実寸シーン）**: `run_v4_1_3.py` の場面定数をそのまま使い 2 秒
  （48 frames, 736×540）。ファイル名 `piece-2s-v4-1-4-realscale.mp4`、静止画 `frame-realscale-0024.png`

## ゲート

- G1: 13条件の σ 比が全て 0.5–2.0（v4-1.2 と同じ帯。in-sample である旨を開示）
- G2: シーンC のユニークフレーム = 48/48、かつフレーム間平均絶対差（緑）が
  0.5 諧調を超えること。リングの実測スペクトルピーク（rad/m と mm）を報告
- G3: シーンB の歪み（緑・平均 px）を v4-1.2 の 1.747 px と並べて記録（合否なし）
- G4: 既存 py の sha256 が stage-0 記録と一致
- G5: 両動画とも 48 frames / 2.000 s / 736×540 / 24 fps（ffprobe）

## 実行規律

- 全ステージ foreground 実行。jobkit checkpoint を `out/v4-1-4/checkpoint.json` に
  stage ごとに原子的にコミット
- 判定は数値統計のみ（目視なし）。後付け調整は禁止 — 帯を外れたら外れたまま報告
- 出力は全て `out/v4-1-4/`： REPORT-V4-1-4.md、numbers-v4-1-4.json、上記 mp4/png、
  checkpoint.json、stage0-py-hashes.json
- REPORT には新会計での σ 内訳（σ_reg / var_i / σ_total）を両シーン分、表で出すこと
