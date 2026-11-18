# BRIEF-DET-WIRE — 高速版 det_fast を本番経路へ配線する（ゲート付き）

Workdir: `D:\yohey2026\dsh\muse\caustics-filter\offline`
背景: `det_fast.det_frame_fast` は等価性ゲート合格・ビット一致（REPORT-DET-SPEEDUP §3：frame 0・24 の C／PNG が本番と完全一致、ファイル SHA256 も一致）。本ジョブはそれを本番ランナーから使えるようにし、本番入口経由でもビット一致を再確認する。

## 手順

1. `run_prod_v4.py`（と、もし `det.det_frame` を呼んでいれば `make_video_v4.py`）の det 呼び出し箇所を特定する。
2. 変更は上記のランナーファイルだけ。det の呼び出しを `det_fast.det_frame_fast`（`row_block=256` 固定）に差し替え、環境変数または引数で従来版へ戻せるスイッチを残す（既定は高速版）。`det.py` と `det_fast.py` 自体は無変更。
3. **ゲート（必須）**：本番入口から frame 0 と 24 をスクラッチ出力先（`out/det-wire/` 配下、`out/prod/` へは書かない）で生成し、C／flow／PNG の SHA256 を `out/prod/` の本番成果物と比較する。**全部一致しなければ配線を元に戻し（revert）、失敗として報告する。**
4. ゲート合格後、48 フレームを 3 ワーカー×4 スレッドで本番経路のまま生成し、壁時計を実測する。`out/prod/` に同フレームの成果物があるものは全て SHA256 を比較する（フレーム 0〜47）。マシンは共用なので、壁時計は「共用負荷下の実測」と明記する。

## 制約・成果物

- `out/prod/` は読み取りのみ。他プロセスは止めない。
- 変更したファイルは手順2のランナーだけ（実行前後の SHA256 で、変更がその範囲に収まることを確認）。
- `out/det-wire/REPORT-DET-WIRE.md`：配線の差分要約、ゲート結果（フレームごとのハッシュ一致表）、48 フレーム実測（壁時計・スループット）、失敗時は revert 済みであることの明記。
- `out/det-wire/numbers-det-wire.json`：生の数値。
