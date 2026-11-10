# RESULT-TWO — Wan 2.2 蛾2匹 1秒テスト（2026-11-10）

結論: 3テイク生成完了・全て検証済み。振り付けは take1・take2 が成立、take3 は不成立（右も飛んだ）。

## 環境・手順
- Python: butterfly-wan\venv（torch 2.14.1+cu130）、ComfyUI 0.44.0（既に起動していたものを再利用、ポート8188）
- 生成スクリプト: scripts\gen_wan2_two.py（前任セッション作成）。組立 ffmpeg = C:\Program Files\SteelSeries\GG\apps\moments\ffmpeg.exe（h264_nvenc）
- 入力: moth-two-832x480.jpg（832x480、141,071 B、MD5 b6db9d8f67f04e6405cf66fc0742c774）。2匹合成の開始フレーム（左=反転・0.55倍・+5°、右=0.60倍・-3°、背景は拡散インペイント+樹皮パッチキルティングで再構成）。マスター PNG は VM 側 delivery/butterfly-wan-two/moth-two-832x480.png
- 条件: 832x480、25フレーム、steps 30、モデルは前回と同じ fp8 セット
- プロンプト: "Two green moths resting side by side on lichen-covered tree bark. The moth on the right flutters its wings gently but stays perched on the bark the whole time. The moth on the left opens its wings, flutters them rapidly, takes off and flies up and out of the frame. Locked-off static camera, photorealistic, natural insect motion."
- ネガティブ: "morphing, extra wings, extra moths, distorted bodies, moving background, watermark"

## テイク結果（gen-wan-two.log の実測）
| take | seed | 生成時間 | MP4 | サイズ | QC |
|---|---|---|---|---|---|
| two1 | 11 | 233s（05:49:19–05:53:13） | out\two1\two1.mp4 | 1,689,365 B | 初回 probe が fps=null で ok=False 表示。再 probe で fps 24.0 を確認し、主担当の ffprobe でも 25f/832x480/24fps/1.041s を検証済み（probe のパース問題であり動画は正常） |
| two2 | 12 | 130s（05:54:08–05:56:19） | out\two2\two2.mp4 | 1,521,216 B | ok=True |
| two3 | 13 | 分離不能（下記） | out\two3\two3.mp4 | 1,676,138 B | ok=True（組立後 probe） |

take3 の経緯: 05:56:20 に生成開始、ドライバ（エージェントセッション）がセッション再起動で死亡。ComfyUI 側で生成は完了し25フレームが output に残ったため、scripts\deliver_two3.py が 06:00:36 に組立・スチル抽出・QC・checkpoint 記録を実施。純粋な生成時間はログから分離できない（開始から組立まで ≤256s）。

## GPU（logs\gpu-monitor3.csv 196サンプル 05:48:54–06:05:36 + gpu-monitor4.csv）
- 最大温度 **68℃**、最大電力 **195.4 W**、最大 VRAM **7,860 MiB** / 8,192 MiB
- **OOM なし**（3テイクとも生成完了）

## 振り付け判定（主担当が f00/f12/f24 を目視）
- **two1: 成立**。f12 で左はモーションブラーで離陸中・右は鮮明に静止。f24 で左は画面外（左上に抜ける途中の残像）・右は静止のまま
- **two2: 成立（離陸は遅め）**。f12 では両方静止、f24 で左のみ画面外・右は静止
- **two3: 不成立**。f12 で右が手前に飛び出す大きなブラー＋左も上方へ退出中。f24 では両方とも不在（樹皮のみ）

## 備考
- エージェントセッションが transport 切断で2回終了した。checkpoint-two.json の completed_ids に seed13_done の重複が残っている（Checkpoint の競合書き込み）が、成果物は全て atomic commit 済みで個別検証済みのため影響なし。
- 合成背景は、蛾がいなくなったフレーム（two3 f24）では中央がやや滑らかな暗い樹皮として見える。蛾が覆っている間は目立たない。
