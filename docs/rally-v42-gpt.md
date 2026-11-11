# Rally v42 — GPT digest (low-wind adoption)

- Thread: https://chatgpt.com/c/6af2ea4a-c8d0-83e8-8cc1-5b232f4f0ead (existing v3-gaps thread, continued)
- Prompt: `rally-v42-lowwind-prompt.md` posted in full 2026-11-12 01:39 JST; answer completed in the same session.
- Model caveat: the answering model name was not shown on screen (display read only "ChatGPT"); the task did not change the model. Recorded as-is.
- **Final ruling: REJECT** the proposed ADOPTED_WITH_BASIS status. (Explicitly NOT a rejection of retaining the 0.003 operational floor.)

## Per-claim verdicts

- **L1 — WEAKENED.** Candidate A introduces fewer unconstrained choices than B–E, but "least invented" is a comparison among candidates, not evidence that 0.003 is a valid capillary-band anchor. Calling it the centre of an interval narrows the commitment only if the interval is genuinely supported; the supplied MSS constraints do not establish that.
- **L2 — WEAKENED.** Not changing any number/gate/status is a real distinction from the v3.9 manoeuvre. But a governance relabel still moves the epistemic boundary — it changes what the project claims grounds to adopt — and criteria 1–4 do not all hold, so the proposed protection does not apply.
- **L3 — DEAD.** The criteria do not establish a capillary-band floor: nonzero wave inception and nonzero total MSS do not prove a sustained capillary-band variance of 0.0015–0.006. Nor do the criteria as written reliably exclude the unvalidated halo; they need a same-observable, same-regime evidence requirement.
- **L4 — DEAD.** The caveat is present but the status line is not honest enough for the stated reading environment: ADOPTED_WITH_BASIS in a "physically correct water" repo will reasonably be read as scientific endorsement, and "not validated; not a gate pass" does not undo that first impression — especially as the claimed basis itself overstates what the measurements constrain.
- **L5 — WEAKENED.** A named supersession target and automatic expiry are good governance controls but do not repair an inadequate present basis. A faithful Kudryavtsev implementation is a model implementation, not a measurement; it may supersede the parameterisation without independently validating the low-wind capillary level.

## Counterarguments answered

- **#1 (goalpost-moving):** Not exactly the v3.9 manoeuvre — C1 is untouched, nothing turns greener. But the scientific defect from v3.9 remains (magnitude and mechanism unmeasured), and the proposal adds a stronger assertion — that measurements positively exclude the zero alternative — which is NOT established for the sustained capillary-band component. Not yet a legitimate adoption with basis.
- **#2 (name one greener thing):** None — that is a strength of the proposal, not a scientific justification for it. Administrative closure is legitimate only if it is explicitly closure of an operational convention, not of the physical question.
- **#7 (criteria excluding rain/halo):** The four criteria are insufficient as written. Rain is excluded by its positively contradicted derivation. The halo is NOT inherently excluded by the four criteria. GPT's tightened admission rule (all five required):
  1. No positively contradicted derivation.
  2. **The basis must constrain the same physical quantity** (same observable/component, relevant regime; a total that mixes components is not enough).
  3. Evidence and modelling assumptions explicitly separated (borrowed scalar, assumed shape, extrapolation each identified).
  4. No gate laundering (no gate result, number, or printed status improves).
  5. A concrete expiry condition whose evidence can actually resolve the specific uncertainty (not merely swap one model for another).
  Under these: rain fails 1; the halo is not automatically eligible (must show its basis independently constrains the claimed optical effect); **the low-wind prior currently fails 2**.
- **#9 (README reader):** Yes, some will infer endorsement. The strongest phrase (ADOPTED_WITH_BASIS) is prominent while the key qualification is buried. The status itself must lead with the unresolved physical state and describe the adoption as operational retention.

## Core objection (the critical error, in GPT's words paraphrased)

The inference from nonzero inception + total-MSS envelopes to a nonzero sustained capillary-band floor is not established. Cox–Munk and Hu ranges show the floored TOTAL MSS is not obviously outside broad low-wind totals; they do not bracket the capillary COMPONENT assigned 0.003. Also: do not call [0.0015, 0.006] an evidence-bracketed interval — it is a declared sensitivity interval around an operational centre.

## Required replacement wording (GPT's draft, key lines)

- Status line: **"UNRESOLVED — OPERATIONAL PRIOR RETAINED (not validated; not a gate pass)."**
- The floor is described as "a modelling prior, not a measured low-wind capillary spectrum"; affected results must report the sensitivity interval σ²_h ∈ [0.0015, 0.006].
- Zero-alternative paragraph weakened: inception/nonzero-roughness measurements argue against a perfectly smooth surface; they do not establish a nonzero sustained capillary-band variance of 0.003, nor rule out a very small capillary contribution under sheltered/film-covered/below-threshold conditions.
- Supersession: reconsider when a published slope spectrum at U = 1–3 m/s resolves k > 200 rad/m AND distinguishes the capillary contribution from total MSS; a Kudryavtsev-2013 implementation triggers a model-comparison review only, not observational validation.
- Governance point preserved for v4: the renderer may use the prior, and v4 may depend on it, without the project claiming the low-wind capillary spectrum has been adopted on a sufficient empirical basis.

## Relation to the Gemini answer

Gemini: ADOPT WITH CONDITIONS, headline condition = rename the identifier to PROVISIONAL_UNRESOLVED_CONVENTION (lead with non-validation). GPT goes one step further: the adoption claim itself is rejected, but the operational retention + UNRESOLVED-first wording converges with Gemini's naming demand in substance. Both reviewers in so far agree the headline must NOT read "ADOPTED". Grok outstanding (collection cron 2026-11-12 04:00 JST).
