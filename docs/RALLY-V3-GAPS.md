# RALLY — v3 知見の不足チェック（GPT / Grok / Gemini(AI Studio)）

目的： 「やりたいこと（地球環境をコードで作り視覚化する）」に対して、v3 の知見・変数セットに足らない部分が無いかを、外部AI 3者に攻撃的に検証させる。開始： 2026-11-10 18:55（ユーザー指示）。

スレッド（トピックごとに1本、ここへ記録して継続使用する）：
- GPT: https://chatgpt.com/c/6af2ea4a-c8d0-83e8-8cc1-5b232f4f0ead （初回応答 2026-11-10、要旨は rally-v3-gaps-gpt.md）
- Grok: https://grok.com/c/c03959aa-48d5-46e8-8cdb-eb289de9a297 （初回応答 2026-11-10 18:57、全文は rally-v3-gaps-grok.md）
- Gemini (AI Studio): 第1ラウンドは Temporary chat で未保存。第2ラウンドで保存済みスレッドを取得： https://aistudio.google.com/prompts/1qoJkIHSoA67BrMjDktzEsFUE59oKUjmy （要旨は rally-v3-gaps-gemini.md / rally-v3-gaps-gemini-r2.md）。今後の継続はこのスレッドへ。
※既存の caustics-output スレッド（research-gpt.md / research-grok.md）とは別トピックなので混ぜない。

## 初回プロンプト（3者共通、設計済み）

あなたは物理ベースの水面コースティクス・レンダラーの設計を攻撃するレビュアーです。擁護や同意ではなく、不足と誤りを探すのが仕事です。

「下記の v3 変数セットで、『地球の環境そのものをノード条件として、水面の揺らぎとコースティクスを物理から生成する』という作品目的には十分である」


- 目的： 写真/動画に、水面屈折による揺らぎ＋コースティクスの反射・煌めきを、環境変数から物理的に生成する。水の状態（スペクトル波形）→光輸送（全 preimage のヤコビアンによる集光）→像、という構造。手続き模様の重ね描きは禁止。
- v3 で実装済みの変数： 太陽（日時・緯度経度→天文学計算、実在モードは Open-Meteo の実天気）、水深、濁度、風と波（波高 H_s の読み出しと固定指定、Cox–Munk 較正の穏やかさ σ²、Smith の波の影）、底の基質アルベドと深度のスペクトル減衰（Beer–Lambert)、雨（雨量→輪波スペクトル、Marshall–Palmer)、人間の網膜（瞳孔・順応・焦点ぼけ・CSF、表示輝度を前提）。
- 既に否定済み： 磁場の影響（約18桁不足）、水圧の光学効果（Δn≈8×10⁻⁶、導出量）、気圧の平均値。既に別バージョン行き： 水温・表面膜・CDOM は v5、立体ボトムと物体は v4 の 2nd、宇宙レベルは v8 で無影響判定。
- 検証の物差しは実在感・奥行き・複雑さ。対象シーンは池・川の浅瀬（深さ 0.3–5m、数 m スケール）が主。


1. 不足している変数・機構を1つずつ： [名前／なぜ物理的に効くか／大きさの見積もり（桁で）／v3 に入れるべきか後のバージョンか／入れないと何がどう嘘になるか]
2. 変数でなく「結合（カップリング）」の欠落があれば別枠で。
3. 実装済みの中に、物理として間違っている・二重計上・前提が崩れているものがあれば指摘。
4. 最後に「v3 が完成と言える条件」を3行で。

遠慮なく反論し、失敗条件を明示してください。一般論でなく、この変数セットとシーン規模に対して答えてください。

## 統合（Muse 評価、2026-11-10）— 3者の攻撃を実装実態と突き合わせた判定

### 3者一致・本物の不足
1. **天空光（直達/拡散の分離）** — 全員一致で最上位。現実装は直達＋一律の環境光項で、天気→直達率→コースティクスのコントラストが結合していない（v4調査の「フラッシュには直達≳30%」と同根）。
2. **吸収と散乱の分離＋前方散乱** — 全員一致。現実装は Beer–Lambert 減衰＋濁度ブラー σ∝√(b·z) の近似止まり。a(λ)/b(λ) 分離と位相関数は無い。
3. **有限水深の分散 tanh(kh) とスペクトル決定の完成** — H_s は SMB フェッチ制限済みだが分散は深水式のまま。ピークも Phillips の風速依存止まり。
4. **流れ（流速ベクトル）** — GPT/Grok。川シーン用。現実装なし。

### ブリーフ未記載ゆえに攻撃されたが実装済み（不足ではない）
- 太陽の有限視直径＋求積 Q=25（Gemini の筆頭攻撃）— オフライン実装に入っている。GPT は前回合意を引いて「維持せよ」と扱った。
- FFT 実現場（Grok「統計だけでは」）— 実現場を持つ。統計は較正アンカーのみ。
- 雨による長波減衰の近似（Gemini の結合指摘）— 3.2 に実装済み。

### 部分的に当たり・要確認
- Cox–Munk: 二重計上ではない（σ² は較正目標のみ）が、池スケールで海洋式を目標に置く妥当性は残る。フェッチ制限した傾斜統計への置換が筋。
- Smith 影： 直達項のみの統計的影として設計したが、実現場とのスケール分離がコード上で明示的か要確認（GPT/Gemini 両指摘）。
- 砕波限界 H≤0.78h: 未実装（現行は Stokes H/λ≤1/7)。浅いプリセットで効きうる安い追加。
- カメラ撮像モデル（GPT): 網膜は後段に正しく配置済み。ただし写真合成側の露出・分光感度は無く、入力が写真の作品としては本物の論点。

## 第2ラウンド（2026-11-10 21:52 発注、ユーザー決定）
反証（実装の実態5点： 有限太陽円盤 Q=25 求積／FFT 実現場＋Cox–Munk は較正アンカーのみ／雨減衰の近似あり／網膜は後段／濁度ブラー反映済み）と、統合評価 A–F を3者に提示し、「崩れる指摘は撤回、残る指摘は再攻撃」形式で継続中。GPT・Grok は同一スレッドへ継続投稿。Gemini は前回が未保存のため、保存される形で開き直し＋前回指摘の要約を前提として同梱。


## 第2ラウンド統合 — 確定版（Muse 評価、2026-11-10）

撤回で確定した「不足ではない」もの： 有限太陽円盤（Q=25 求積）、FFT 実現場、Cox–Munk の二重加算なし、網膜の後段配置、雨の長波減衰。3者とも反証を受けて撤回した。

### 確定した不足（3者の合意が厚い）
- **P0-A 天空光の分離**： 天気→直達率 η=DNI/(DNI+DHI)→コースティクス・コントラストの結合が無い。3者一致で最大の欠陥。
- **P0-B 吸収・散乱の分離**： a(λ)/b(λ)＋位相関数（HG g≈0.85）。ブラー近似は代替にならない。3者一致。
- **P0-C スペクトル決定の完成**： 有限水深分散 tanh(kh)＋表面張力項＋方向分布＋**Cox–Munk アンカーの短フェッチ置換**（Gemini の新攻撃： F̃≈39 の池で海洋 σ² に合わせると高周波傾斜を 1/2–1/3 に潰す）。毛管波がカスプの鋭さを作る点も Gemini/GPT/Grok で一致。
- **検証 P0**： Q=25 求積の収束試験（9/25/49 点、ピークで見る）、パッチ 1.48m＝最大波長の明示と低波数切断の扱い（GPT 新規）。

### 意見が割れた／コードで決着させるもの
- **E Smith 影**： GPT「コード確認を P0 に」、Gemini「決定論追跡なら全廃せよ」。我々の輸送は前方ヤコビアン写像で遮蔽テストを内蔵しないため Gemini の前提は強すぎる。実装確認（直達項のみ・統計的影として入っている）をもって判断する。
- **F 砕波限界**： GPT「基準条件ではほぼ発動しない（優先度下げ可）」vs Gemini「浅端 h=0.1m では H>h の数学的破綻を防ぐ必須クランプ」。両立する → 安いクランプ H_max=min(0.142λ, 0.78h) を P1 で追加するのが落とし所。

### 優先度低で確定
- **D 流れ**： 池主体では二次的（3者一致）。川シーン拡張時／v4 の流れと一緒に扱う。
- 雨滴の詳細結合・クレーター波、カメラ撮像モデル（v4 のカメラモードへ接続済み）、底面アルベド 2D 分布・BRDF は P1–P2。

---

## v3.7 round (2026-11-11) — D1–D5 design under test, same threads

New evidence: the spectrum diagnosis (delivery/spectrum-diag/) —
TWO independent defects (H_s = 4×waveHeightM wind-only curve;
capillary tail absent + taper anchor cuts to target + lambda_p =
patch artefact + U=10 silent no-op).
Verdicts: GPT — D1-replace SURVIVES, D2 WEAKENED, D3 taper-retire
SURVIVES / 0.034-swap WEAKENED, D4 SURVIVES/WEAKENED, D5 dispersion
SURVIVES / quantisation WEAKENED. Grok — D1 WEAKENED, D2 DEAD
(splice continuity unspecified), D3 WEAKENED, D4/D5 SURVIVES.
Gemini — D1 WEAKENED, D2 DEAD, D3 DEAD (ECKV's own B_h integrates
above 0.034; taper was the only closed-loop check), D4/D5 SURVIVES,
plus an eta=1 Hessian/quadrature breakage claim.
Synthesis adopted in PIPELINE-V9.md: additive two-regime spectrum
(JONSWAP via exact Jacobian conversion + ECKV B_h ADDED — no
splice, no renormalisation); NO MSS target (MSS is an output,
reported vs CM / retired target / full-ECKV); an acceptance
battery G0–G10 replaces the retired anchor; depth correction
regime-gated; cascades partition ONE spectrum on a shared
quantised clock; low-wind alpha_m floor explicit + logged.
Digests: rally-v37-{gpt,grok,gemini}.md.

## v3.8 round — gravity-tail termination (2026-11-11)
Claims D1–D5 from research_notes/gravity-tail-termination-spectrum-20261111-0528. Digests: rally-v38-gpt.md / rally-v38-grok.md / rally-v38-gemini.md.
- D1 (level, not cutoff): SURVIVES (Gemini, Grok) / WEAKENED (GPT). Finite-depth attack dead at all three (kh≥15 in band). GPT's condition adopted as design gate A1: reproduce B_g from the implemented transform under the code's convention (deep + finite depth, ratios reported).
- D2 (envelope reset): WEAKENED at all three — circularity. Resolution: ECKV fraction is a conformance reference only; independent anchors (Cox–Munk, Hwang, Laxague) carry validation. Gemini adds: ECKV lacks the young-sea parasitic-capillary source, so even conformance may under-represent a real Ω_c≈6.3 pond — recorded as a known missing source.
- D3 (B_l swap): SURVIVES (Gemini) / WEAKENED (GPT) / DEAD-as-pure-Elfouhaily, WEAKENED-as-hybrid (Grok). Consensus path: proceed as a challenger, honestly labelled "Elfouhaily-form spectrum with fetch-derived sea-state parameters", behind a compatibility gate (unchanged renderer, old vs new: H_s, height quantiles, slope/curvature, caustic stats). ΔH_s demoted from rejection to report+flag (Gemini's scaling-exponent argument: a 15% gate rejects valid cells systematically).
- D4 (defer Kudryavtsev): SURVIVES at all three.
- D5 (battery): rewritten into three classes — conformance / validation / compatibility; bit-identity scoped to identical buffers (Gemini's contrary reading was of an imprecise prompt sentence, clarified in the design).

## v3.9 round — low-wind capillary floor (2026-11-11)
Claims R1–R6 + decision (a)+(d), from research_notes/low-wind-capillary-floor-20261111. Digests: rally-v39-{gpt,grok,gemini}.md.
- R1 (alpha_m zero is a fit artefact): SURVIVES at GPT (scope-limited), Grok, Gemini.
- R2 (C1 unsatisfiable at U=1.5): WEAKENED at GPT (matched-definition care), SURVIVES at Grok, Gemini.
- R3 (floor at sigma2_h=0.003): WEAKENED at all three — not a measured pond capillary level; labelled provisional prior only.
- R4 (suspect flips at U=3): WEAKENED at GPT, SURVIVES at Grok, Gemini.
- R5 (Kudryavtsev no rescue): WEAKENED at GPT (implementation vs mechanism), SURVIVES at Grok, Gemini.
- R6 (parasitic increment unquantified): WEAKENED at GPT, Grok; SURVIVES at Gemini (direction only above the crest-steepness threshold).
- Decision (a)+(d): REJECTED at all three (goalpost-moving). Converged requirements: three-state gates with low-wind C1 = UNRESOLVED (GPT); a wind-independent minimum-capillary test for all cells (Grok); absolute-B envelope reported + viscous-transition physics noted (Gemini); sensitivity sweep of the prior (GPT).
- All three name the same decisive measurement: a lab/field slope spectrum at U=1.5–3 resolving k>200 rad/m. Not available to us — design carries UNRESOLVED openly.
Design: PIPELINE-V11.md (chunk 9 = v3.9). Implementation brief: BRIEF-V13-FLOORFIX.md.

## Audit round 2026-11-11 (3.1–3.9 full audit + 3.9 boundary) — digests rally-v40-{gpt,grok,gemini}.md
Consensus: 3.1 STANDS; 3.2 WEAKENED→flagged REJECTED_INVALID_DERIVATION (Gemini label); 3.3 STANDS; 3.4 STANDS as machinery; 3.5 WEAKENED (η=1 = scope decision, not a repair); 3.6 WEAKENED (merged as best-available, not validated — battery compares the model to itself); 3.7 WRONG as a validated default (3/3) → operational fallback with defects disclosed; 3.8 STANDS as challenger; 3.9 WEAKENED. Boundary ranking: (b) prior extension rank 1 at GPT+Grok; (d) absolute-saturation gate redesign rank 1 at Gemini, on a new physical argument (viscous dissipation ∝4νk² makes capillary dominance at U≤3 possibly unattainable; the 50% share gate is a high-wind intuition). GPT/Grok admit (d) only with independent justification — Gemini's B-envelope test is that candidate justification. Settlement (PIPELINE-V12): three-arm matched measurement (RAW / B / D), winner by independent constraints only, low-wind cells UNRESOLVED at best (3/3), no universal-default merge on a prior. Gemini's flat-surface runtime fallback NOT adopted (project scenes are low-wind ponds; labelled uncertainty, not a deleted surface).

## B-envelope adjudication round 2026-11-11 — digests rally-v41-{gpt,grok,gemini}.md; research: research_notes/benvelope-cap-literature-20261111-1014
Question: the chunk-10 merge blocker — the B-envelope clause (max B(k) <= 0.012 over [10,1500], from Gemini's audit battery) failed 3 of 4 high-wind cells identically in all arms. Literature verdict: the clause is a CATEGORY ERROR — 0.012 is Phillips' equilibrium-range coefficient beta, a MEAN (Juszko 1995: mean 0.012, range 0.006-0.024), valid only in 1.3-3 f_p; no source caps the total spectrum with it. Elfouhaily's own B_l(k_p) = 0.1433*alpha_p*gamma equals 0.0120 @ Omega 4.9, 0.01765 @ 7.3, 0.02411 @ 10.3 (flank x1.04 -> 0.0183 / 0.0250 = our failing cells exactly); B_h(k_m) = (1/2)*alpha_m = 0.01249 at U=10 (implementation alpha_m = 0.024970, Eq. 44 high branch — the research's alpha_m flag resolved: my mission text's 0.016 was my own error, the implementation was already Eq. 44).
Rally: C1 SURVIVES 3/3. C2 SURVIVES (GPT/Gemini conditional: beyond Omega<=5 it is model-conformance, not empirical validation). C3 SURVIVES at Grok/Gemini, WEAKENED at GPT (JONSWAP sub-km values extrapolate the alpha fetch law). C4 SURVIVES (GPT conditional: consistency check, not independent physics). Retirement = principled correction, NOT goalpost-moving (3/3). Replacement: region-wise gate — (a) gravity peak/flank vs Elfouhaily's own value x1.5 tolerance + Omega>5 FLAG; (b) equilibrium fitted beta in [0.006, 0.024]; (c) capillary vs (1/2)*alpha_m. GPT condition: report conformance vs observation as TWO columns; for Omega>5 cells also print the fetch-anchored (JONSWAP) alternative peak. Grok conditions: warning travels with every output; regression stays green; H_s as prediction with a FORWARD 15% tripwire (adopted forward; on this matrix |dHs| is -21%/-22% at the two F=100 cells — annotated, not retroactively re-gated; the pre-registered 35% K1 stands). Gemini condition: WARN_EXTRAPOLATED_WAVE_AGE tag; UNRESOLVED flags persist.
MERGE APPROVED 3/3 with those conditions -> chunk 11 (BRIEF-V15-MERGE.md, run_v15): elfouhaily-form (Arm-D rules) becomes the default spectrum (v3.10); jonswap stays as rollback. GPT restriction recorded: the new battery is a conformance + region-wise observation battery, never "physical truth".

## v42 round — low-wind adoption (2026-11-12)
Question: adopt the σ²_h = 0.003 low-wind prior as ADOPTED_WITH_BASIS (claims L1–L5; research: research_notes/low-wind-adoption-20261112; prompt: rally-v42-lowwind-prompt.md). Digests: rally-v42-gemini.md / rally-v42-gpt.md.
- **Grok: SKIPPED by the user's one-time exception** (2026-11-12 01:51: 「grokで止まるのもったいないからパスでいいよ」「今回だけ特例」). The prompt was posted in the thread but rate-limited; the 04:00 collection cron was cancelled the same night (automation registry updated). This round's record is a two-way synthesis; the exception is for this round only.
- Gemini (AI Studio, 3.8 Flash): **ADOPT WITH CONDITIONS**. L1 SURVIVES, L2 WEAKENED, L3 SURVIVES, L4 WEAKENED, L5 SURVIVES. Conditions: (1) rename the identifier to PROVISIONAL_UNRESOLVED_CONVENTION with the opening line stating Gate C1 NOT PASSED (counterargument #9 "valid and fatal to the draft's headline"); (2) explicit pond-slicks disclaimer (true capillary MSS may be orders below the bracket floor on sheltered/film ponds; convention = intentional upper-bound clean-water); (3) automated test that `low_wind_prior` provenance is present in every low-wind output. #1 answered: not v3.9 goalpost-moving (gates untouched). #7: rain fails criterion 1, the 3.6 halo fails criterion 2 — locked out.
- GPT: **REJECT** the adoption claim (retaining the 0.003 operational floor is NOT rejected). L1 WEAKENED, L2 WEAKENED, L3 DEAD, L4 DEAD, L5 WEAKENED. Core objection: the same-observable failure — Cox–Munk/Hu constrain TOTAL MSS, not the capillary-band component assigned 0.003, so criterion 1 ("zero positively excluded") and the "bracket" do not hold for the component; [0.0015, 0.006] is a declared sensitivity interval, not an evidence bracket; a Kudryavtsev-2013 implementation is a model, not a measurement (supersession = model-comparison trigger only). Tightened 5-clause admission rule recorded in rally-v42-gpt.md (low-wind fails clause 2; rain fails clause 1; halo not automatically eligible). Required replacement status: **UNRESOLVED — OPERATIONAL PRIOR RETAINED (not validated; not a gate pass)**.
- Convergence (2/2): no number, gate outcome, or printed status changes; the sensitivity sweep ships; provenance travels; the headline must lead with UNRESOLVED and must not read "ADOPTED". Difference: Gemini permits a provisional-convention adoption framing; GPT rejects the adoption claim itself.
- Synthesis (Muse): GPT's wording satisfies Gemini's headline condition and the user's closure intent — the implementation decision is closed, the physical question stays open, and v4 may depend on the prior. The final identifier is the user's choice (asked 2026-11-12, choice-format); implementation of the settled wording + the provenance test follows in code and repo docs.
- **RESOLVED 2026-11-12 02:24: the user chose GPT's wording** — UNRESOLVED — OPERATIONAL PRIOR RETAINED. Implementation dispatched the same night as chunk 13 (BRIEF-V17-LOWWIND-STATUS.md: new `lowwind_status.py` + `test_lowwind_provenance.py`, no existing file modified); the repo-facing document was written to the v3 ALPHA freeze staging as `docs/LOW-WIND-STATUS.md`.
