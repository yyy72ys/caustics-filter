# BRIEF — v3 チャンク4（3.4）最終完了セッション（残ステージのみ）

状況（実測済み）：
- フル解像度ランは [M] フォトン検算の段で死亡した（python.exe プロセス消滅、out/v8 は 23:21 以降無変化）。
- out/v8 に在る物： C-baseline-v33.npy、C-v8-t0.npy（フル解像度）、still-preset1-v33.png、still-preset1-v8.png。
- out/v8-smoke2 にスモーク規模の全ゲート結果（numbers-v8.json / numbers-v8-prereq.json / panel / ablate・clamp・shallow の C 群）が揃っている。
- ログ out/v8-run.log の到達点： [E] Smith 監査完了、v8 状態（sigma2_spec=0.020954, realised=0.022135, Hs=0.1140)、[Q] 収束完了（Q=9/25/49 で p99 1.5097/1.5098/1.5095)。

やること（これだけ）：
1. 既存の C-v8-t0.npy を再利用し、不足している段だけを**このセッションのフォアグラウンドで**実行する（バックグラウンド化・セッション外プロセス化は禁止 — セッションが終わる前に全ファイルが書かれていること）：
   - [M] フォトン検算＋MS-SSIM 分解（フォトン数は v2 検算と同条件。メモリ不足で死んだ疑いがあるので、ピークメモリを抑える実装にし、ピーク値を記録する）
   - フル解像度でまだ無いゲート項目（ablation 群・浅水検証・クランプ検証 — スモーク版が在るものはフル版を作成）
2. out/v8/numbers-v8.json と out/v8/REPORT-v8-gates.md を書く（前提タスク判定 [E] の結論と根拠、[Q]、[M]、ゲート1–6 の合否と数値、死亡したランの扱いと再実行範囲の明記を含む）。
3. 各段の完了ごとに D:\mutaro-work\jobkit 形式の checkpoint（進捗ファイル）を out/v8/checkpoint.json に追記する（途中で死んでも到達点が残るように）。

厳守： `taskkill /F /IM python.exe` 禁止。PC では別作業（proto5）が並行中。モデルは deepseek-v4.1-flash 固定（呼出側で指定済み）。
