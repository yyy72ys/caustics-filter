# BRIEF-V4-1-12 — G-LINE-stat 診断の実装と初回測定（v4-1.12）

## 目的
採用済み設計 `DESIGN-GLINE-STAT.md`（改訂版・地位＝「分布診断」）を実装し、既存の stream シーンで初回診断を行う。**本ジョブは測定と診断のみ。合否判定は出さない**（閾値 1.05/1.25 は暫定配置値であり、本運用は新規独立 seed での検証後。レポートでは参考線としてのみ描く）。

## STAGE 0 — 前提確認（絶対パス。無ければ STOP して報告。別ツリーを探さない）
作業ツリー: `D:\yohey2026\dsh\muse\caustics-filter\offline`
- `run_v4_1_10.py`（warp-only 再生経路の正本）
- `out\v4-1-10\numbers-v4-1-10.json`（§5 の 5 realization の seed と frame 24 の窓別 depth が記録されている）
- `stream_surface2.py` / `v5_composite.py` / `run_v16_water.py`（1.10 と同じ warp 経路）
- 写真側の数値: `out\v4-1-10\numbers-v4-1-10.json` 内の photo 4 窓 depth（63.7/61.6/47.4/39.3‰）と配置感度 sd

## STAGE 1 — 行動再現ゲート（通らなければ STOP）
1.10 と同一経路で production seed・frame 24 の warp-only 4 窓 depth を再計算し、`numbers-v4-1-10.json` の値と**完全一致**すること。あわせて G-WARP（mean |d|）= 10.067452 px を再現すること。一致しなければ理由を報告して止まる（自分の都合で経路を変えない）。

## STAGE 2 — 実装（新規ファイル `run_v4_1_12.py` のみ。既存ファイルは変更禁止）
測定対象: 1.10 §5 と同じ 5 realization（production + 代替 4 件、seed は numbers から読む）× フレーム {0, 16, 32, 24} の warp-only 画像。
1. **窓別 depth**: G-LINE と同一の定義（5–15 px 帯、4 窓 sx/sy/bx/by、写真 63.7/61.6/47.4/39.3‰ との比）。
2. **集計は realization 単位**: 各 realization で先にフレーム平均を取り、分布は R=5 の realization 平均の列から作る。報告量 = 窓ごとの平均比・sd・SE（realization 平均から）・クラスタ・ブートストラップ（realization を再抽出）90% CI。**n=15 という表記は使わない**（実効 n = R = 5 と明記）。
3. **フレーム相関 ρ**: realization 内で窓 depth のフレーム間相関（0–16 / 16–32 / 0–32）を窓別に測り、プールした値も出す。相関を仮定せず実測値として報告。
4. **prom2d との相関（G-LINE2 独立性チェック）**: DESIGN-GLINE2 と同一の定義（2-D FFT リング 5–15 px の max/median、256² 窓、Hann、HP σ=3.0）で frame 24 の realization 別 prom2d を測り、窓 depth（frame 24）との realization 間相関を窓別＋全体で出す。n=5 の相関であることを明記し、独立の主張はしない（相関の実測値を出すまでが仕事）。
5. **写真側**: 1.10 の配置感度 sd を「配置感度の約 ○ 倍」という言い方でのみ併記する（σ による有意差の言い方は禁止）。
6. 振幅正規化（Sn）版の窓 depth も 1.10 と同様に併記する（実現パターンの散らばりか振幅かの切り分け用）。

## STAGE 3 — レポート
`out\v4-1-12\REPORT-V4-1-12.md` と `out\v4-1-12\numbers-v4-1-12.json` を書く。
- 窓ごとの分布表（平均比・sd・SE・90% CI・配置感度比）
- フレーム相関 ρ の実測値と実効 n の扱い
- prom2d 相関の実測値（独立性の判定はしない、数だけ出す）
- bark-y の系統的超過（1.10 で +12%）が分布で見ても残るか
- **合否は書かない**。閾値 1.05/1.25 は「暫定配置値（未検証）」として参考表示のみ
- 限界の明記: R=5 は薄い、写真は 1 枚、seed の選び方（+1/+2/+3/+7）は無作為抽出ではない

## 実行規律
- jobkit checkpoint を stage ごとに書く（`out\v4-1-12\checkpoint.json`）
- 前景で実行し、成果物が書けるまでセッションを終えない
- 既存ファイルの変更はゼロ（開始時と終了時に主要ファイルの sha256 を記録し、レポートに差分なしを明記）
- 画像の目視は使えない環境なので、判断はすべて数値で行う
