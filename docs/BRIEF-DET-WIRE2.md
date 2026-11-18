# BRIEF-DET-WIRE2 — det_fast2（radix+fused）の本番配線

出典: `out/det-binsort/REPORT-DET-BINSORT.md`（cf-0062、判定 PASS・ビット一致＋高速）。配線の型は cf-0056（`out/det-wire/REPORT-DET-WIRE.md`、caustics-det-wire）の踏襲。
Workdir: `D:\yohey2026\dsh\muse\caustics-filter\offline`
成果物は `out/det-wire2/` のみ新規。変更してよい既存ファイルは **`run_prod_v4.py` 1 本だけ**。

## 目的

本番入口 `run_prod_v4.py` の det 既定実装を、現行の `det_fast.det_frame_fast` から `det_fast2.det_frame_fast2`（`sort_mode="radix"`, `deposit="fused"`、row_block=256 固定）へ切り替える。出力は現行とビット一致しなければならない（cf-0062 で frames 0・24 の C・PNG 一致、frames 0–5 のバイト一致を実証済み）。

## 配線の形（cf-0056 と同じ）

- det 呼び出しは既存の `build_det()` 経由のまま。既定の実装を fast2 に変える。
- 戻し口を 2 段用意する: 環境変数 `DET_IMPL` ／ CLI `--det-impl` で `fast2`（新既定）／ `fast`（現行 det_fast）／ `base`（元の det）を選べること。cf-0056 で入れた `base` の戻し口は壊さない。
- `det_fast2.py` 自体は変更しない（cf-0062 の納品物のまま。SHA256 `4fd299001b4b713e14457ef4c0555c77cd7331f4a416440879667cdb034c1ec0`）。

## 必須ゲート（1 つでも落ちたら配線を自動で元に戻し、FAIL として報告）

1. **本番入口ゲート**: `run_prod_v4.py` 経由で frame 0 と 24 を生成し、C・flow・PNG の **6/6 が `out/prod/` の既存成果物と SHA256 一致**。
2. **戻し口ゲート**: `DET_IMPL=fast` と `DET_IMPL=base` の両方で frame 0・24 が同じく 6/6 一致（各切替で 6/6）。
3. **48 フレーム全数**: 本番入口の既定（fast2）で 48 フレームを回し、**144/144 が `out/prod/` と一致**。壁時計も測る（cf-0056 の実測 594.85 秒＝12.39 秒/フレームが比較対象。負荷条件は正直に記録）。
4. **変更範囲**: 変更は `run_prod_v4.py` 1 本のみ。実行前後に offline 配下の全 `*.py` の SHA256 を取り、他ファイルの変更が 0 件であることを示す（cf-0062 の hash_py の型）。`out/prod/` への書き込み禁止（読み取りのみ）。

## 注意（既知の仕様）

- ランナーの seam 検査は `--out` 先の `C-0000.npy` を読む既存仕様で、フレーム分割実行時に落ちることがある（cf-0056 §7 記録済み）。これは配線の欠陥ではない。48 フレームは分割せず 1 回で回すこと。
- JIT ウォームアップ（約 5 秒/プロセス）は測定から除外するか、その旨を明記。
- 実行は `python -B`、Numba cache=False（既存の作法どおり）。

## 成果物（`out/det-wire2/`）

- `REPORT-DET-WIRE2.md`（配線内容の diff 要約、ゲート結果、48 フレーム壁時計、変更範囲の検証）
- `numbers-det-wire2.json`（ゲートの一致数、ハッシュ、時間）
- `checkpoint.json`（jobkit 形式、工程ごとに更新）

## やらないこと

- `det_fast2.py`・`det_fast.py`・`det.py` の変更。
- ゲート失敗時の原因の勝手な修正（戻して FAIL 報告までが仕事）。
- 192 フレーム本走（48 フレームまで）。
