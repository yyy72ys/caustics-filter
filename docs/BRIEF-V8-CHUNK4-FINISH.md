# BRIEF — v3 チャンク4（3.4）完了専用セッション

前任セッションは transport エラーで落ちた。到達点は実ファイルで確認済み：
- 作成済み： offline/v8_physics.py、v8_prereq.py、v8_state.py
- out/v8/ には C-baseline-v33.npy と still-preset1-v33.png のみ（前提タスクの途中）
- 未了： run 系の編成、ゲート一式、numbers-v8.json、REPORT-v8-gates.md

やること： 正本ブリーフ D:\mutaro-work\rally-in\BRIEF-V8-CHUNK4.md と設計 PIPELINE-V8.md を読み直し、既存の v8 コードを点検してから、残りだけを完了させる（前提タスクの判定・Q収束・MS-SSIM 分解、チャンク4 のゲート1–6、成果物一式を out/v8/ へ）。既に書かれたコードは活かし、書き直しは最小限に。

厳守： この PC では別の作業（proto5）が並行して走っている。`taskkill /F /IM python.exe` のような一括終了は絶対にしないこと。自分のレンダーが固まったら、その PID だけを終了する。
