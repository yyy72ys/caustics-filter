# ROADMAP-V4 — コースティクス研究と3段階の拡張（ユーザー指示 2026-11-10）

v4 の主題は「コースティクスそのものの研究」。光の出方を決める要因を洗い出し、物理として正しく取り込む。

## ユーザーが挙げた要因（2026-11-10）
- 光の入射角でコースティクスの出方は変わるはず
- 水面と水中の水の流れ
- 水面のしぶき
- 揺れによる水中の空気（気泡）で反射が変わる
- ボトムに光が当たる／ボトムの凹凸があれば影になる複雑さ
- 影も大事な要素
- **視点（カメラ）の位置・距離・角度**（2026-11-10 19:43 追加）： 見る側の幾何学も変数に入れる。物理的な効き方： 視線角で水面の Fresnel 反射/透過の比率が変わる（斜めから見るほど水面が鏡寄りになり、水中が見えにくくなる）、屈折の変位量と実効光路長が変わる、2nd の立体ボトムでは視差と遮蔽が生まれる。v3 の網膜（観察者の生理、合成後の後段処理）とは別物 — こちらは光輸送の幾何そのもの。なお v3 不足ラリーで GPT が挙げた「カメラ撮像特性（露出・分光感度・MTF）」は露出・記録側の話で、こことは区別して扱う。
- **視点は「目」と「カメラ」の2モード切替**（2026-11-10 19:51 追加、ユーザー）: 観察者を人間の目にするかカメラにするかを切り替えられるようにする。目モードでは v3.3 の網膜研究が視点側に接続される — 視点の位置・方向・角度が「目に入る光量」を決め、それが瞳孔の調整・順応・網膜と脳への刺激に連動する（v3.3 では表示輝度を前提に置いた順応輝度 L_a が、ここでは視点を通して実際のシーンから決まる）。カメラモードでは位置・距離・角度に撮像特性（露出・分光感度など）が付く。同じ放射輝度場を、どちらの観察者で受けるかの切替。
- **知覚・サンプリングのノード群**（2026-11-11 01:15 追加、ユーザーの気づき）: 解像度は世界の性質ではなく観察者のサンプリングの性質であり、脳は入力をわざと間引いてメモリを節約している — その仕組み自体をノード群にする。目モード側： 中心窩の高解像と周辺への解像度フォールオフ、注視点／注視の経路（サッカード）、周辺の脱飽和と時間解像度の低下、変わらない場所は再サンプリングせず予測で埋める（補完）、そして総サンプル予算（メモリ上限）を注意で配分するパラメータ群。カメラモード側の対応物： センサー解像度・画素ピッチ・MTF・シャッターのサンプリング。位置づけ： v3.3 網膜（瞳孔・順応・ボケ・CSF＝目の光学と感度）とは別の層で、こちらは「どこをどれだけ細かく取るか」の配分の層。同じ放射輝度場を、予算つきの観察者がどう間引いて受け取るか — 現行パイプラインの一様 2048² は「予算無限のカメラ」の仮定であり、このノード群がその仮定を外す。

## 3段階の流れ（ユーザー確定）
- **1st（現在）**： 2D メディア＋水面シミュレーション（v3.x 系。写真・動画に水面とコースティクス）
- **2nd**： 水面の前にボトム（立体）を作り、そこに物体を配置する有機的なシミュレーション（ボトムの凹凸・物体の影・立体へのコースティクス投影）
- **3rd**： 立体水面（水そのものを立体として扱う段階。Blender の「場で保たれる水オブジェクト」の parked idea と合流する）

## 進め方（v3 と同じ型）
調査（チャンク0相当）→ 統合（機構への接続）→ 段階ごとに設計→実装→数値ゲートと目視で検証。

## 状態
- 2026-11-10： v4 調査を deep research で開始（入射角／流れ／しぶき／気泡／ボトム凹凸と影／影の一級要素化、の6問）。成果は research-v4-caustics/ に統合する。
- v3 の残り（チャンク3 網膜）は未完のまま保持。扱いの順序は調査結果と合わせて決める。

## Viewpoint — user reference photos (2026-11-11)
User brought two reference photos and asked whether viewpoint changes how caustics appear even in shallow water:
- `refs/ref-shallow-caustics-web.jpg` — looking THROUGH the surface (near-normal): crisp caustic web projected on a shallow pebbly bottom.
- `refs/ref-grazing-green-surface.jpg` — looking AT the surface at a grazing angle: reflection dominates, bottom gone, chromatic (dispersion) fringes on the ripple highlights.
Recorded physics: viewing angle sets the Fresnel reflect/transmit split (R ~2% near-normal → 21% at 75° → mirror at grazing), so the same water reads as "bottom caustics" or as "surface mirror" purely by viewpoint; sun glitter is defined around the specular direction (viewpoint-bound); the bottom caustic WEB is strong even in shallows because focusing is set by short-scale surface curvature, while depth mainly scales the lateral warp of the submerged image (consistent with the user's z=4 m favourite). This is direct evidence for the v4 viewpoint branch (EYE/CAMERA + angle as first-class variables).

## v4 最初の設計対象（ユーザー決定 2026-11-12）

**流れのある水面**——流れの強い川の上流の水面から始める。経緯： v20 の対照実験で「v2 の迫力の主因は古い強すぎる風波だった」と確定したあと、ユーザーが「逆に v3 で 23px 作る環境はあるか」と問い、深い澄んだ湖（深さで稼ぐ線）と並んで「流れの強い川の上流」を挙げ、**v4 の最初はそれで**と決定した。

物理の見立て（着手前の整理）：
- 上流の水面の荒れは風ではなく流れが作る：水が底の岩を越えるときの盛り上がりと落ち込み（定在波・ボイル）。傾きは風波より急になりうる。歪みは「深さ×傾き」なので、浅くても傾きで稼げる可能性が本命の仮説。
- 現行エンジン（v3.x）の水面生成は風スペクトルのみで、この水面は作れない。流れで水面がどう変形するかの計算が新しい核になる。
- 進め方は v3 と同じ型：調査（チャンク0相当）→ 機構への接続 → 設計 → 実装 → 数値ゲートと目視。
- なお v3 は 2026-11-12 に ALPHA として push 済み（本決定と同日）。雨（3.2）の再導出も v4 時代の課題として残る。
