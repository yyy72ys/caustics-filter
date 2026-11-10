# @LexnLin（Leon Lin） — ESCAPEMENT：全コード駆動の Blender 短編映画

対象投稿： https://x.com/LexnLin/status/2119558994352304202 （2026-11-08、6投稿スレッド）
一次情報： 公開リポジトリ https://github.com/Leonxlnx/escapement （README 読了、2026-11-10）
調査日： 2026-11-10。

## 工程（リポジトリ README より）
- 36秒の短編映画を Claude Code（Opus 5.5）が Blender 5.2 headless を Python で駆動して制作。統合GPUのノートPC（Radeon 780M）で実行。
- **全オブジェクト・マテリアル・ライト・カメラが Python 生成**。ダウンロードした3Dモデル／画像テクスチャはゼロ。
- **BRIEF.md（production bible）＋ shotlist.json**（全ショット： id・スクリプト・カット尺・前後ハンドル）を先に置く。
- **幕ごとに builder エージェントを並行稼働 → director エージェントが批評 → 2巡目で修正**。反復ログは notes/act*.md に残す。
- レンダーは tools/rb.py 経由に一本化し、**GPU ロック**（gpu_guard.py）で同時実行を禁止（iGPU のブルースクリーン 0x119 を踏んだための恒常対策）。
- **"Measured, not eyeballed"**： exposure_check.py が各スチルを「白飛ばし／黒浮き／色の dull さ」で採点。コンポジタはハイライトを白に届く前にロールオフ。
- 編集・グレーディング・音響ミックスもコード（numpy＋ffmpeg、音楽の120BPMグリッドに同期）。
- 正直な但し書きも README に明記： 雲の作り直しはGPUを凍らせて不採用、カットしたショット、速い動きは83%解像度・16サンプルで描いてアップスケール。

## 本人の訂正（スレッド6件目）
- 「every frame is Python」は正確には「Python が全シーンを構築し、描画は Blender の EEVEE」。事後の訂正を自分で出す姿勢も含めて記録。

## 周囲
- 車シーンの参照： @Ror_Fly の工程（Three.js でブロッキング → Midjourney でシネマトグラフィ参照 → Seedance で本描画、GPT-6 Astra の多角ブロッキング）。「安い媒体で構図を固めてから高い描画に行く」分業。

## 蛾への処方として効く点
- スチルQCの自動採点（露出系）→ 蛾では「シルエット重なり率・翅色の写真実測との差・影の差分量」を採点項目に置き換える（SYNTHESIS (A) 参照）。
- BRIEF＋shotlist の先行固定、director 批評パスの分離、手順と反復ログのファイル化 — proto3 の進め方にそのまま移植可能。
