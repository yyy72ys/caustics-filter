# FIX-12: コースティクスバッファのフィルタリング・フォールバック（検証環境で実証された実バグの修正）

## バグの内容（メイン側の切り分け実験で確定）
- `initCausticsTarget` は half-float / float のレンダーターゲットにも無条件で LINEAR フィルタを設定している。
- WebGL1 では half-float テクスチャの LINEAR サンプリングには別拡張 `OES_texture_half_float_linear`、float には `OES_texture_float_linear` が必要。無い環境ではテクスチャが不完全扱いになり、`texture2D` は (0,0,0,1) を返す → **コースティクスが無言で消える**。
- 実証（SwiftShader 環境）: スプラットの出力を定数 0.7 に強制しても合成側の Cn は 0 のまま。ターゲットのフィルタを NEAREST に変えただけで場（fold 線・セル・細部）が出現した。
- 実GPU環境では概ね拡張があるため動くが、無いクライアントでは全滅する実バグ。加えて、メイン側の検証経路（SwiftShader）では FIX-7 以降ずっと前方バッファが見えていなかったことになる。

## 修正要件（app.js の initCausticsTarget のみ）
1. 候補タイプごとにフィルタを決める:
   - half: `OES_texture_half_float_linear` があれば LINEAR、無ければ NEAREST
   - float: `OES_texture_float_linear` があれば LINEAR、無ければ NEAREST
   - byte (UNSIGNED_BYTE): 常に LINEAR
2. 候補の優先順位（half → float → byte）と FBO 完成チェックは現行どおり維持。
3. 検証用に、選択されたフィルタを `window.__caustics.caustics.filter`（'linear' / 'nearest'）として露出する（既存の mode/w/h/grid と並べて）。

## 凍結事項（重要）
- スペクトル、DISP_K（=2.0）、帯域ゲイン、合成式、グリッター、その他すべての数式・定数は変更しない。今回の変更はターゲット生成のフィルタ選択と検証フックのみ。
- 理由: バッファが見えない状態で下した視覚判断は無効であり、まず正しく見える状態を作ってから再評価する。

## 完了条件
- node --check 全 JS OK。verify-node.js にフィルタ選択ロジックの文字列検査を追加し全項目 PASS。
- RESULT.md に FIX-12 節（変更点のみ簡潔に）、STATUS.md に完了記録。
- ブラウザ実測は不要（メイン側が SwiftShader / 実GPU の両経路で実施する）。
- 作業範囲は caustics-filter/ 内のみ。未達があれば STATUS.md に未完を明記し完了扱いにしない。
