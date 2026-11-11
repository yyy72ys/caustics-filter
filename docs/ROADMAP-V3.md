# ROADMAP v3 — 変数の拡張と磨き（ユーザー指示 2026-11-10）

進め方（ユーザー決定）：全部一気にはやらない。**リサーチ → バージョン分け → チャンク単位で更新**。
各変数は「数学的に正しいもの」に磨く（恣意的なスライダーではなく、物理の式と既存変数への連動で定義する）。

## ユーザーが挙げた不足変数

### A. 水面まわり
- 波の高さ（波高）
- 穏やかさ（波高・風との連動量として定義する想定）
- 波の影（水面自身の陰影／底に落ちる波の陰）

### B. 底まわり
- 水面とボトムの光の反射具合（表面反射率／底の反射率＝アルベド）
- ボトムのレンジ（底までの距離＝水深レンジと、光の届き方の変化）

### C. 人間の網膜（新ノード群）
- 網膜に見えている波・コースティクスとの関係。
- 網膜側の変数もノードで継ぎ接ぎできるようにする（他のノード群と同じ流儀で）。
- 候補：瞳孔径（明るさで変わる）、順応（明暗への慣れ）、網膜像のぼけ（ピント・被写界深度）、コントラストの見え方。

### D. 変数そのものの磨き
- 既存変数も含め、数学的に正しい定義へ寄せる（例：波高は風速から式で決まる、穏やかさは派生量として定義、など）。

### E. 新シーン
- **無風に雨**：風のない穏やかな水面に雨が降るシーン（雨滴が作る輪波の世界）も作ってみたい。

### F. 環境の連動そのもの（ユーザー 2026-11-10 03:56 の着想）
「結局、環境のバタフライエフェクト的に波が揺らいでる」――波の場は環境変数すべての連動の結果だという描像。リサーチ対象に追加：
- 重力：水面波は重力波。分散関係 ω² = gk·tanh(kh)（＋表面張力項）に g が入る。重力を変数にできるか、毛管長（重力⇔表面張力の境目）も含めて。
- 磁場：通常の水の波に効くのか、正直に評価する（反磁性。効かないなら「効かない」と明記）。
- 吹走距離（フェッチ）と水深：小さな池は風でもさざなみまでしか育たない（fetch/depth-limited growth）。「水の量・体積に規定される揺らぎ」の物理量への翻訳（フェッチ・水深・水域サイズ・副振動セイシュ）。
- 波が波を作る：非線形の四波相互作用（スペクトル内エネルギー移動）。
- 川の流れ：流速のドップラー効果と移流。流れは常にあり速度もある、という川の素性。
- 地震・プレート（ユーザー 2026-11-10 03:58 の着想）：海底の変位→津波。池や湖の水は地震動で副振動（セイシュ）として揺れる（固有周期は水域サイズと水深で決まる）。遠地地震でもセイシュが起きる実例があり、「地球の裏側の地震で池が揺れる」＝環境連動の最たる例。フィルターでは連続変数ではなくイベント（インパルス）型の入力として扱うのが自然。

## チャンク計画（案）

- **チャンク0（リサーチ）**：A〜E 全変数の数学的土台を調べる。各変数について「定義・式・現実の値域・既存パイプライン変数との連動・実装メモ」を1本のレポートにまとめる。正直に式にできないものは、その旨を明記する。
- **チャンク1（v3 第一弾）**：水と底の変数（波高／穏やかさ／反射具合／ボトムレンジ／波の影）を物理連動で実装し、数値検証して納品。
- **チャンク2（雨）**：無風＋雨モード（雨滴の輪波）を実装。
- **チャンク3（網膜）**：網膜ノード群を実装（環境ノード群と継ぎ接ぎ可能に）。

順番はユーザー確定（2026-11-10）：**チャンク1 水と底 → チャンク2 雨 → チャンク3 網膜**。

## V6（少し先・ユーザー指定 2026-11-10）：メタ変数マップ

- 「地球の水は何の影響を受けているのか」を全体として調べる。それがそのまま変数の全体像になる。
- メタ的視点が必須：水の近くの変数だけ拾うのではバタフライエフェクト（遠くの・間接の連動）を拾えない。影響のグラフ全体を先に描いてから、どれをノードにするかを選ぶ。
- 候補の射程（例）：太陽と月の引力・潮汐、風と気圧、気温と蒸発、地球の自転（コリオリ）、雨と河川の流入、地形と水域の形、地震とプレート、生物、人間の活動……。
- これは v3 のチャンクとは別枠。**今はやらない**（ユーザーが「少し先」と明言）。着手は V6 のタイミングで。

## 再訪課題：雨の計算（ユーザー指示 2026-11-11）
- チャンク2（v3.2）として実装・ゲート通過済みだが、ユーザーは「雨の計算はまだ怪しい気がする」と判定。**今の作業が片付いたら、v4 の前に雨をもう一度やる**（再導出・再検証）。
- それまで雨は本編・実験の表舞台では使わない（2026-11-11 の強風池スティルでも雨は入れない決定）。

## 既知の不整合：勾配アンカー vs SMB波高（2026-11-11 強風池スティルで発見）
- 風6m/s・フェッチ100m・z=1.0 のスティルで、実測 σ²=0.0556（アンカー目標通り）に対し実測 H_s=0.332 m。SMBフェッチ則の H_s は 0.031 m で**10倍の不一致**。現行スペクトル形状では「正しい傾き」と「正しい波高」が両立しない（実在の池は毛細波の裾が傾き分散を担い波高をほぼ増やさない）。スペクトル形状の分離（重力波ピーク＋毛細裾）が今後の課題。雨の再訪課題と同族（スペクトル決定 P0-C の残件）。
- 副次発見：パッチ実寸が風で広がる（U=6 で 5.61 m、プリセット1は 1.48 m）。写真への割当てが風と連動して「カメラが引く」ため、メートル単位のワープが昨日並み（平均変位 48.5 mm）でも画素単位では穏やかに出る。パッチ↔写真の対応は場面の実寸で固定すべきか、設計論点として残す。
- 進め方（ユーザー事前承認 2026-11-11）：v3.7 は診断→リサーチ→ラリー→設計（PIPELINE-V9）→**承認待ちなしで実装まで進めてよい**（事前承認済み）。結果で報告する。
- 診断確定（2026-11-11, delivery/spectrum-diag/）：不一致の原因は2つ独立。(1) 実測 H_s = 4×waveHeightM（風だけの振幅カーブ、フェッチ・水深なし）— スペクトル形状は波高に入っていなかった。(2) 毛細裾の欠落＋テーパーアンカーは「k_c より上を切って目標に合わせる」仕掛けで、ピークも解像されていなかった（λ_p = パッチ寸法）。U=10 ではアンカーが無言で無効化。
- v3.7 設計確定（PIPELINE-V9、ラリー3者統合）：加算型2レジーム（JONSWAP を分散のヤコビアンで厳密変換＋ECKV B_h を加算、接続なし・再正規化なし）／σ² 目標は廃止し出力として3参照と並記／受け入れは G0–G10 の電池式ゲート／深さ補正はレジームゲート式／カスケードは単一スペクトルの分割＋共通量子化クロック／低風速 α_m フロアは明示＋記録。実装は同日 07:26 に deepseek-v4.1-flash で発注済み（台帳 caustics-v37-chunk7）。

## v3.8（チャンク8）：JONSWAP裾の打ち切り — ユーザー決定 2026-11-11 14:27
- 経緯：v3.7 の未達2件（1cm未満帯域の割合 0.056–0.137 vs 文献目安 0.25–0.40／G9 p99.9 バンド 4.72・2.73 vs 6–12）は同根 — JONSWAP 平衡裾が 50–200 rad/m の傾き分散の87%を占め毛細帯域を二重数えしている（v3.7 レポート §1・§7-1 実測）。
- ユーザー判断の軸：「本当の自然をつくるなら改修」— 見た目を盛る改修ではなく「数え方を正す」改修として実施する。目安数字は海の実測由来で池スケールの実測検証は存在しない点は本人に明示済み（改修で言えるのは「最も信頼できるモデルの作法への一致」まで）。
- 進め方（v3.7 と同じ流れ）：打ち切りルール候補のリサーチ → 三者ラリー（既存 v3-gaps スレッド継続）→ 設計（PIPELINE-V10）→ 実装（PC・deepseek-v4.1-flash）。自由パラメータのあるブレンド/接続は失格（ラリー D2 の確定事項）。
- リサーチ発注済み（2026-11-11 14:28、research_notes/gravity-tail-termination-spectrum-20261111-0528、agent 835025c0）：候補 = ECKV B_l 本体（若い海 Ω_c≈6.3 での妥当性を含む）／Kudryavtsev モデル／Hwang 2005／Donelan–Banner–Plant・Apel／k_t=√(ρg/γ) アンカー／Phillips 飽和限界。目安 0.25–0.40 自体の堅さも検証対象。
- 受け入れの芯：G2 帯域配分と G9 バンドの再測＋G8 回帰電池（凍結コンポジットのビット安定ほか）が壊れないこと。G10 風池の再描画で配分の変化を確認する。
- リサーチ完了（2026-11-11、research_notes/gravity-tail-termination-spectrum-20261111-0528/report.md）：診断が組み替わった。欠陥は「裾の位置」ではなく「50–200 rad/m での水準」— 変換後の JONSWAP 裾は曲率スペクトルで B_g≈α_J/2 の plateau（0.006–0.023）に張り付き、実測飽和水準 B̂≈0.007–0.008 以上。分散変換自体が 370 超では裾を殺している（係数 0.25@370、0.02@1000）ので、k_t=370 での打ち切りでは 87% は直らない。
- 目安の訂正：k>370 で 0.25–0.40 という v3.7 の検証エンベロープは文献が支持しない（Cox–Munk の slick 分割は k≈21、Laxague の純毛細 k>1173 は ≈5%）。正しい目標は「同じ入力での ECKV 自身の >370 割合」（推定 0.10–0.25）。
- 第1候補（Rank 1）：S_g を Elfouhaily B_l に丸ごと置換（S=k^-3(B_l+B_h)、継ぎ目なし・C∞）。k_p/Ω_c は当方のフェッチ則から。H_s は入力から予測に変わる（セルごとに報告、再スケール禁止＝自由パラメータ化を避ける）。Ω_c=6.3 は ECKV 公表域（≤5）外の外挿である点が露出箇所。
- Kudryavtsev は第2候補として保留（α(k/k_g) 等の式が公開源で未確認、総量が数十%動く恐れ）。三者ラリー（D1–D5）を既存 v3-gaps スレッドへ発注済み（2026-11-11）。

## v3.8 implementation outcome (2026-11-11, REPORT-V11-GATES.md)
MERGE WITHHELD — default stays `jonswap`. A1/A2 PASS (diagnosis algebra verified), K1 no flags (worst dHs 24.7%), K4/K5/K6 PASS. Blockers: (1) C1 fails in the 4 low-wind cells — inherited v3.7 low-wind alpha_m floor zeroes B_h at U=1.5 (band 100% gravity by construction), NOT a B_l defect; (2) C2 fails all cells because its premise was false — v3.7's "totals within 1-2% of full-ECKV" was circular (v10 compared against an integral of its own spectrum); against true matched-input ECKV, v3.7 totals are 0.81-3.32x the reference, and the new spectrum matches ECKV to 0.12% and is necessarily 40-70% lower. C1 and C2 cannot both hold. K3 windy: warp 42.0->23.7 mm, fine texture -27.6% (absolute >370 variance -54% despite share rising), caustic peaks stronger (p99.9 2.73->3.84, max 4.05->6.57); legacy still byte-identical to v3.7's. Next lever if the merge is wanted: the low-wind floor itself (a v3.7 inheritance), plus the known parasitic-capillary gap.

## v3.9 (chunk 9) — low-wind floor, rally-corrected (2026-11-11)
Research reframed the C1 failure: (1) the alpha_m hard zero is an implementation invention (ECKV's log fit crosses zero at u*=c_m/e, U≈2.72; ECKV prescribes nothing below); (2) gate C1 is unsatisfiable at U=1.5 under matched definitions. The research's proposed fix (floor 0.003 + re-scope C1) was REJECTED by all three rally AIs as goalpost-moving; R3 weakened 3/3 (the 0.003 is an ocean residual, not a measured pond capillary level). Design PIPELINE-V11.md implements the rally consensus instead: hard zero removed; below-domain level = LABELLED provisional prior (sigma2_h target 0.003, provenance printed); three-state gates (low-wind C1 = UNRESOLVED, disclosed); Grok's minimum-capillary test C1b for all cells; sensitivity sweep of the prior mandatory; C2 retired with reason. Implementation dispatched 2026-11-11 (run_v13*, out/v13). Merge re-evaluation per the design's rule. Full record: RALLY-V3-GAPS.md v3.9 round; report research_notes/low-wind-capillary-floor-20261111/report.md.

## v3.6 (chunk 6) — IMPLEMENTED + MERGED 2026-11-11
G1–G4 all PASS (delivery/offline-v12/REPORT-V12-GATES.md). The a/b split + HG single-scattering path is now the default water model (WATER_MODEL_DEFAULT=ab); legacy recipes keep the pre-3.6 lane byte-identical. Measured headline: the derived scattering halo is ~9.3× wider than the old empirical blur yet carries only ~8% of the light — the old model was wrong in both directions at once. Same-c falsification pair delivered (absorption split keeps 1.26× the contrast of the scattering split; green level 1.47× higher under scattering).

## v3.9 (chunk 9) — IMPLEMENTED 2026-11-11, MERGE WITHHELD (U=3 boundary)
Delivery: delivery/offline-v13/REPORT-V13-GATES.md. The labelled prior works where active (U=1.5: gravity share 1.000→0.840, sigma2_h delivered = 0.0030, cells UNRESOLVED as designed). But the merge is withheld again, on a NEW and precise blocker: U=3 sits just INSIDE the ECKV fit domain (u*/u*_fit = 1.113), so the prior does not act there, raw alpha_m = 0.001067 leaves the 50–200 band 85–89% gravity → C1 FAIL (both U3 cells) and C1b FAIL (U3/F100: share 0.106, sigma2_h 0.001305). This is rally claim R4 made concrete: at U=3 the suspect is ECKV's level, and no design on the table may invent a floor above the crossing. A1/A2, K1, K4–K6 PASS; sensitivity sweep delivered. Next lever, if pursued: the domain-boundary treatment itself (e.g. Gemini's viscous-transition form) — a new design + rally round, not a patch.

## v3.9 redo (chunk 10) — MEASURED 2026-11-11, MERGE WITHHELD (B-envelope clause)
Delivery: delivery/offline-v14/REPORT-V14-GATES.md. Three arms (RAW = v13 as shipped; B = prior extended to every cell whose raw sigma2_h < 0.003; D = RAW numbers judged at low wind by Gemini's absolute criteria) over the same 8-cell matrix. Construction gates ALL PASS in all arms; regression GREEN (v13 lanes byte-identical, v12 default `ab`, frozen modules unchanged). Low-wind cells: UNRESOLVED in every arm, no FAIL. At U=3, Arm B delivers sigma2_h = 0.003 (capillary share 0.214 at F=100 / 0.259 at F=1000; B max 0.0118 / 0.0045, inside the absolute envelope); Arm D shows RAW already satisfies the absolute criteria there (B <= 0.012, B_l band max <= 0.008 in 8/8 cells for every arm). Optical consequence of B vs RAW at U=3/F=1000 is negligible (MS-SSIM 0.99999, max pixel diff 6/255, own-shift +1%). Decision under PIPELINE-V12's pre-registered rule (applied by Muse as decision-maker): the arms tie on every independent constraint; the winner reading is Arm D (fewest invented inputs among the tied arms — B adds an invented input at U=3 that buys no constraint or optical benefit; D's stated win condition, B-envelope passing at U=3, is met). Going forward inside the elfouhaily-form family: keep v13's numbers, judge low wind by the absolute criteria, low-wind cells stay UNRESOLVED and labelled. MERGE WITHHELD anyway, for every arm: 3 high-wind cells FAIL the B-envelope clause imported from the Gemini battery (max B(k) <= 0.012 over [10,1500]): U=6/F=100 B=0.0183 at k=14.9 (gravity peak flank), U=10/F=100 B=0.0250 at k=10.5 (gravity peak), U=10/F=1000 B=0.01249 at k=368.2 (capillary peak). Whether 0.012 is a legitimate universal cap — including over a young sea's gravity peak, where Elfouhaily's own B_l peak runs higher — is itself unverified; that clause is now the blocker. Next lever: verify the B-envelope clause against literature (Elfouhaily B_l peak vs inverse wave age; measured saturation ranges), NOT a patch and NOT a quiet gate edit. Default spectrum remains `jonswap` (v3.7, defects disclosed).

## v3.10 (chunk 11, run v15) — MERGED + VERIFIED 2026-11-11/12
Delivery: delivery/offline-v15/REPORT-V15-MERGE.md (PC run r4, exit=0, report landed 2026-11-11 20:47). The B-envelope clause was retired after the literature check + three-way rally (0.012 is Phillips' equilibrium-range mean β, never a total-spectrum cap; the three "failing" values equal Elfouhaily's own formula outputs). Default spectrum is now `elfouhaily-form` under the Arm-D rule set; `jonswap` retained as rollback (regression item 5: byte-identical lane). Region-wise two-column battery replaces the clause: conformance (a)/(c) PASS in all 8 cells; observation (b) reported per cell (PASS 3, FAIL_OUT_OF_BAND 2 at low wind — reported, not gated there — UNRESOLVED band-contaminated 3 at high wind, accepted per the brief's rule). Verdicts: PASS 4 / FAIL 0 / UNRESOLVED 4 (all low-wind, by design). Transport through the merged stack: all three gates PASS; labels travel (T2 carries UNRESOLVED, T3 carries warn_extrapolated_wave_age=true). Regression all pass; no v11–v14 file modified (hash-verified). DISCLOSED JUDGEMENT (§8.2 of the report): U=10/F=100 fitted β = 0.00584 sits 2.6% below the 0.006 floor AND is band-contaminated; the brief's literal "contaminated = UNRESOLVED, never silent PASS/FAIL" rule was applied, so the cell passes the high-wind rule — under the stricter reading it would FAIL and the merge would not stand. Closeout verified by the owning chat 2026-11-12: pulled deliverables' sha256 match the report's §7 hashes exactly.

## v3.6 correction (chunk 12, run v16) — COMPLETED + CLOSEOUT VERIFIED 2026-11-12
Delivery: delivery/offline-v16/REPORT-V16-SCATTER.md (job caustics-v16-scatter-r2, exit=0). The validation research (research_notes/scattering-validation-20261112) REFUTED the v3.6 phase function: HG g=0.85 scores 0/5 against published measurements (Petzold average-particle ⟨cos θ⟩=0.924, B=0.0183; HG puts 2× too much light backward and 4.6× too much at 90°). Scorecard (pre-registered, β̃(90°) reference declared = Petzold flat-back 0.00291 sr⁻¹): FF-Petzold (n=1.10, μ=3.5835) 5/5 (β90 ratio ×1.442), HG-0.924 3/5 (shape clauses fail — the defect is the FORM, not just g). Halo invariants: H bit-identical to v12 in all channels and both lanes (phase acts on the kernel only); σ scales exactly as √(⟨θ²⟩ ratio) (algae ×0.694); halo share unchanged. Kernel shape: no Gaussian is within 0.05 encircled energy of the FF kernel at r = 8.5/32/79 px → the FF angular CDF is used directly as the kernel. Regression all green (legacy HG+Gaussian lane bit-identical incl. G1 sha cfb086b22a0e; v15 spectrum cell identical; 21 frozen hashes; no v11–v15 file modified). **MERGE RULE FIRED: PHASE_DEFAULT = "ff", KERNEL_DEFAULT = "ff-cdf" (run_v16_water.py).** Sediment preset stays a labelled ASSUMPTION (nominal FF 1.15/4.0, envelope printed, phase_assumption=true). Disclosed caveats: H held at the v12 budget by the brief's pre-registered reading (the alternative ρ-detour reading would move H by +0.54%, report §6); the Jonasz turbid FF fit was not available in the checkout — the clear fit is printed for the record, not gated. Closeout verified by the owning chat 2026-11-12: all 5 hash-listed deliverables pulled and sha256-matched against the report.

## v3 closing plan — "v3 ALPHA" (user decision 2026-11-12)
Sequence fixed by the user: (1) close the two remaining known limits — v3.6 scattering validated against published measured data (Petzold VSF / phase-function literature), and the low-wind cells upgraded from UNRESOLVED to ADOPTED_WITH_BASIS (fresh three-way rally in the existing v3-gaps threads; the basis wording must also live in the GitHub repo docs); (2) produce the 8-second piece with the v3.10 standard stack (base v2 frames, η=1.0, frozen composite A2.6+Blow+E2.2) — that closes v3 as **v3 ALPHA**; (3) rain (v3.2) stays ON HOLD past the close, carried into the v4 era as a known deferred item, not a blocker; (4) then v4 starts. Delivery gate: the finished 8 s video goes to the user; his next action is watching it and approving the GitHub push. No push before that approval.

Outcome note (2026-11-12): (1) Scattering — CLOSED, see the chunk-12 section above. (2) Low wind — rally v42 ran TWO-WAY (Grok skipped by the user's one-time exception, recorded in RALLY-V3-GAPS.md): GPT REJECTED the ADOPTED_WITH_BASIS claim (the cited basis constrains total MSS, not the capillary-band component; the interval is a sensitivity interval, not an evidence bracket); Gemini ADOPT WITH CONDITIONS (rename the identifier, UNRESOLVED-first headline, provenance test). Both agree: no number, gate, or printed status changes, and the prior stays in operational use. The final status identifier is the user's choice (asked 2026-11-12); once settled, the wording + the `low_wind_prior` provenance test go into code and the GitHub repo docs, and the 8 s piece proceeds. (Settled 2026-11-12: the user chose GPT's wording — UNRESOLVED — OPERATIONAL PRIOR RETAINED. Implementation dispatched as chunk 13 / BRIEF-V17-LOWWIND-STATUS.md; repo doc written to the freeze staging as docs/LOW-WIND-STATUS.md.)
