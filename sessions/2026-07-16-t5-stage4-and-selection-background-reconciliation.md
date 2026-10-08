# 2026-07-16 — T5 stage-4 completed, and the effect re-homed to the selection background

*Private session log (sessions/, gitignored). Saved at J.B.'s request. One session
spanning execution (finishing the T5 whitening on data already on disk), a conceptual
refinement from J.B. about where the vacuum actually acts, and the manuscript edit that
refinement forced. Formal artifacts are in the repo; this is the narrative record.*

## 1. Execution: T5 stage 4 completed on the v2 data (commit df3a9f9)

The 8-band re-stream on disk had quietly grown past the 48.3 d documented in
PILOT_RESULTS.md to **85.9 d / 39,653 chunks**, but was never analyzed and sat
uncommitted (would have vanished on a clean checkout). Completed stage 4 on it with **no
new streaming**:

- New script `whiten_fold.py` combines the two half-steps that existed separately —
  `whiten_test.py`'s whitening (regress quantum band b7 on auxiliary bands b1–b6) and
  `fold_diurnal.py`'s sidereal + solar/lunar-tidal fit — so whitening is measured at the
  level of the **κξ bound**, not just N_eff.
- Whitening removes **66%** of b7's classical variance; τ **3.5 h → 2.0 h**; N_eff ×1.7.
- **The new result:** whitening tightens the vector bound **×2.1 (κξ 34 → 16 of natural)**
  on matched 85.9-d data — *more* than the ×1.3 the √N_eff error-shrink predicts, because
  it also **halves the sidereal amplitude** (1.67 → 0.86 × 10⁻²): classical power was
  sitting *in* the sidereal band, and b1–b6 regression removes it while preserving any
  genuine quantum modulation (which lives in b7/b8, correlated +0.83).
- **Honest caveat:** at 85.9 d the design condition is **653** (solar/sidereal/P1 still
  collinear), so κξ < 16 does **not** supersede the definitive 256.8-d v1 bound
  (condition 5.4, **κξ < 5.3**). This is a *method* result; projecting the ×2 whitening
  gain onto the full baseline points to κξ ~ 2–3 — better, but still above natural (κξ=1).

## 2. Conceptual refinement (J.B.): the vacuum acts at the boundary, not in the dissipation

J.B. drew a distinction he felt my T5 summary had flattened: the quantum vacuum is nearly
inert in **free flight**; its significant contribution to the wavefunction's evolution
happens **at the boundary** between near-vacuum and condensed matter (his example: a photon
at a prism surface, or a photon meeting an electron bound to an atom). That boundary effect
is **non-dissipative** — weak, and primarily **adds randomness**. *Measurement* is the
**further, dissipative** evolution as the system locks to the neighborhood of atoms;
dissipation arises from the **mass asymmetry** (condensed matter ≫ electron).

This is the **fluctuation/dissipation split**, and it is already the paper's own (§3.1),
which I had blurred: §3.1 states three things are distinct — *selection* of an outcome
("the role of the background configuration, **not of the dissipation itself**"),
*decoherence*, and *erasure* — and makes **Stage-3 registration equivariant**,
"contributes no statistics of its own." So:

- **Vacuum = fluctuation.** The zero-point background, sampled at the boundary where the
  wave's amplitude is deformed (§2.5), supplies the stochastic, outcome-selecting kick —
  weak, zero-mean, non-dissipative. This is **Stage 2 (selection)**.
- **Bulk = dissipation.** The mass-asymmetric coupling to the heavy bath is the
  irreversible closure — **Stage 3 (registration)** — and is frame-independent.
- **The prism nails it:** refraction is a *coherent* boundary interaction (a phase shift,
  no record) that does **not** measure — proof that boundary interaction ≠ dissipative
  measurement. §2.5 calls this a coherent (vs recording) boundary.

## 3. The reconciliation this forced (commit c8424ca)

If Stage-3 registration "contributes no statistics of its own," a **sidereal statistic
cannot originate there** — so T5 was mis-filed. Edited `LIGO_SIDEREAL_TEST_T5.md`
throughout to re-home the frame-dependent effect from the **Stage-3 dissipative
registration** channel to the **Stage-2 selection background**:

- **T2 and T5 both read the Stage-2 selection background**, at different moments: T2 its
  **mean** (a coherent phase, GHZ-amplified); T5 its **variance** (the anisotropic
  quantum-noise floor). Stage-3 dissipation is frame-independent (mass-set) and hosts
  neither.
- Observable (sidereal PSD modulation), reach table, pilot results, and vector/tensor
  verdicts are **unchanged** — only the theoretical seat moved. Added a dated reconciliation
  note; retitled "Quantum-Noise Decoherence" → "Quantum-Noise **Anisotropy**."
- **Bonus for the hardest linchpin:** "dormant in free flight" becomes cleaner — the vacuum
  fluctuation is always present but **nets to zero except where a boundary breaks the
  symmetry**, rather than a coupling switching on/off.
- **Caveat kept honest:** at pilot sensitivity the data cannot distinguish anisotropic
  *added noise* from anisotropic *loss* (both modulate the floor); the re-identification is
  theoretical — decisive for the §7 derivations, not adjudicated by the current bound.

## 4. Open items

- §7 derivations recast onto the selection background: boundary-only selectivity (now the
  zero-mean-except-at-boundary argument), vector-vs-tensor symmetry of the selection
  coupling, and sidereal coherence of the selection background across the band.
- The full 256.8-d 8-band re-stream (to get a whitened *full-baseline* bound) remains
  blocked on unmetered bandwidth, or needs the collaboration's internal aux channels.
- `t5_pilot/smoke.csv` left untracked (a reducer smoke-test artifact, not stage-4 data).
- The manuscript companion `SIDEREAL_DECOHERENCE_PAPER.md` likely still carries the old
  Stage-3 framing — not yet reconciled; flagged for a later pass.
- A tracked `discussions/` entry could distill §2–3 for the public archive if wanted; this
  private log is the raw record.
