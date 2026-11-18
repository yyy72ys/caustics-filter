# BRIEF-ABERLE-H — Aberle DEM 13面の H_loc 測定（cf-0081）

## 目的

受領・保管済みの Aberle DEM 13面について、凍結の H 物差しで H_loc を測り、cf-0059（底ラインの次設計）の判断材料にする。Aberle DEM が再建の参照床として使えるかどうかは、ここでは決めない（判断は cf-0059 で本人がする）。測るのは数値の並びだけである。

## 背景

- キラーテスト（cf-0058、REPORT-V4-2-10.md）は Hodge DEM（D50 25–34 mm）で現行の粒生成アーキテクチャを閉じた。実測床の H_loc ≈ 0.3 は測定の産物ではなく床の性質で、帯 4–36 mm で自己アフィンでない。mean 埋めレーンは検出力なし、no-fill レーンに検出力あり。
- Aberle DEM（cf-0014 で受領）は 13面、d50 9.09–25.3 mm と Hodge より細かい。格子 DEM（2383×180、dx 約1 mm・dy 4 mm）がある。
- 問いは「Aberle の床も Hodge と同じ H_loc ≈ 0.3 の性質を持つか、それとも違うか」である。

## 対象と条件

| 項目 | 値 |
| --- | --- |
| データ | PC `D:\yohey2026\dsh\muse\caustics-filter\data\dem-alt\aberle\` の .mat 13面（I_Qc_120／I_Qc_180／I_Qc_220／I_Qc250a／I_Qc250b／II_Qc120／II_Qc180／II_Qc220／III_Qc120／III_Qc180／III_Qc240／IV_Qc120／IV_Qc180） |
| 各面 | 変数 `xyzdaten`（N×3 float64 の点群）。格子 2383×180、dx 約1 mm、dy 4 mm、x 0–2385 mm、y 0–716 mm |
| 粒径 | `General_Data.xls` の d50（面別 9.09–25.3 mm） |
| 参照値 | Hodge DEM 参照の H_loc 0.289（キラーテストの DEM 参照値、REPORT-V4-2-10.md §3） |

## 方法

1. **物差しは凍結のものを使う。** `run_v4_2_5.py` の H 関数を import して verbatim 使用する（式は写し直さない。キラーテストが `run_v4_2_10.py` で同じ方法を取った）。H_loc(ℓ) は同じラダー上の 0.5·dlnS/dlnd（np.gradient）。
2. **レーンは no-fill のみ。** 両端とも有効なペアだけを使う。mean 埋めは検出力がないことがキラーテストで確定しているので使わない。
3. **方向は分ける。** 格子が異方性（dx 約1 mm・dy 4 mm）なので、x 方向と y 方向のラダーを別々に取る。y 方向のラダーは 4 mm ラグから始まることを記録する。
4. **帯はキラーテストと同じ。** 4–36 mm 帯の平均 H_loc を面ごとに報告する。帯ごとの H_loc(ℓ) の表も出す。
5. **同一性ゲート（測る前）。** 先に Hodge DEM 参照の H_loc 0.289 を同じスクリプトで再現する。Hodge DEM は PC `data/dem-alt/` にある（cf-0053 の実データ）。再現値が 0.289 と一致しなければ（許容はキラーテストの記録と同じ桁、差 0.01 以内）測定に進まない。
6. **d50 との対応。** 各面の H_loc（4–36 mm 帯平均）と General_Data.xls の d50 を並べた表を出す。d50 と H_loc の関係が読めるかどうかも一言記録する。

## 判定の見方（事前固定）

このジョブは測定であり、合否ゲートではない。報告する数値は次の3点である。

- 13面それぞれの H_loc（4–36 mm 帯平均、x／y 方向別）
- Hodge 参照 0.289 との差（面ごと）
- d50 との対応表

「Aberle DEM を再建の参照床として使えるか」の判断はしない。判断は cf-0059 で本人がする。測定の結論は数値の並びそのものとする。

## していないこと（禁止事項）

- 生成器・粒生成コードを動かさない。新しい床を生成しない。既存 `*.py` の変更は 0 件（新しいスクリプトだけ `out/aberle-h/` に置く）。
- キラーテストの再発注ではない。Hodge の再測定は同一性ゲートの範囲でのみ行う。
- 本番配線・設計変更・ビルド・動画化なし。

## 成果物

`out/aberle-h/` に置く。

- REPORT-ABERLE-H.md（本レポート）
- numbers-aberle-h.json（全数値：面×方向の H_loc ラダー、帯平均、d50 対応、同一性ゲート）
- checkpoint.json
- run_aberle_h.py（診断スクリプト）

## 再現手順

```
python out/aberle-h/run_aberle_h.py --stage0
python out/aberle-h/run_aberle_h.py --stage1
python out/aberle-h/run_aberle_h.py --stage2
python out/aberle-h/run_aberle_h.py --stage3
```
