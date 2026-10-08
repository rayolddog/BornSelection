# Superposition_Equations.md
## Working notes — synchronization interpretation of spin superposition (Case 1: electron spin ↑/↓)

### 0. Scope and relation to prior work

This document sets up the **first concrete case** of the synchronization-based
account of superposition proposed alongside the Born-selection / measurement
microscopic simulation work: the ordinary spin-up/down superposition of a
free electron, modeled via the positive- and negative-energy components of
the Dirac spinor rather than as an abstract two-level system.

This is **Interpretation (B)** from the prior discussion: synchronization
between an observer's internal clock and the system's internal clock(s) is
a physical coupling that does measurement-like work, not merely a change of
reference frame.

**Reconciliation note:** the coupling here is reconstructed from the
free-particle Dirac Hamiltonian (α, β structure) rather than copied from
your existing K = m identification in PAPER_UNIFIED.md / EQUATIONS.md,
since I don't have those exact expressions in front of me. Section 3
in particular should be checked against your established notation before
this is folded into the main framework — the physics should agree, but
notation and normalization conventions may not.

---

### 1. Rest-frame Dirac spinor basis and the spin-flip subtlety

Dirac representation, γ⁰ = diag(1, 1, −1, −1), free Hamiltonian:

$$ H = c\,\boldsymbol\alpha\cdot\mathbf p + \beta m c^2 $$

Plane-wave solutions labeled by momentum **p** and spin index s = 1, 2:

$$ \psi^{(+)}_{p,s}(x) = u_s(p)\, e^{-i(E_p t - \mathbf p\cdot\mathbf x)/\hbar}, \qquad E_p = +\sqrt{p^2c^2+m^2c^4} $$
$$ \psi^{(-)}_{p,s}(x) = v_s(p)\, e^{+i(E_p t - \mathbf p\cdot\mathbf x)/\hbar} $$

At rest (**p** = 0), the four independent solutions are:

$$ u_1(0)=N\begin{pmatrix}1\\0\\0\\0\end{pmatrix},\quad
   u_2(0)=N\begin{pmatrix}0\\1\\0\\0\end{pmatrix},\quad
   v_1(0)=N\begin{pmatrix}0\\0\\1\\0\end{pmatrix},\quad
   v_2(0)=N\begin{pmatrix}0\\0\\0\\1\end{pmatrix} $$

with u₁, u₂ at E = +mc² and v₁, v₂ at the E = −mc² eigenvalue.

**Important subtlety to carry forward:** under hole-theory / charge-conjugation,
the *physical* spin of the positron state built from vₛ is opposite to the
label s. If channel "↑" is meant to denote a single physically consistent
spin axis, it must be built from u₁ paired with the v-spinor whose physical
spin is also ↑ — not v₁ by label. This needs to be nailed down explicitly
before Section 2's pairing is trusted; I've assumed the "same physical spin"
pairing below and flagged it inline.

---

### 2. The internal beat: one "clock" per spin channel

Define the spin-↑ internal state as a superposition of the positive- and
negative-energy solutions sharing the same *physical* spin axis:

$$ \Phi_\uparrow(t) = a\, u_\uparrow\, e^{-imc^2t/\hbar} + b\, v_\uparrow\, e^{+imc^2t/\hbar}, \qquad |a|^2+|b|^2=1 $$

and analogously for spin-↓ with amplitudes c, d:

$$ \Phi_\downarrow(t) = c\, u_\downarrow\, e^{-imc^2t/\hbar} + d\, v_\downarrow\, e^{+imc^2t/\hbar}, \qquad |c|^2+|d|^2=1 $$

The relative phase between the positive- and negative-energy amplitudes
within a channel evolves at the **Compton (zitterbewegung) frequency**:

$$ \omega_C = \frac{2mc^2}{\hbar} $$

This relative phase — not the overall phase — is what I'm calling the
**internal clock** of a channel:

$$ \varphi_\uparrow(t) = \omega_C t, \qquad \varphi_\downarrow(t) = \omega_C t $$

At rest and in zero external field, the two clocks tick at the *same*
frequency and are dynamically independent — nothing in the free Hamiltonian
mixes the ↑ and ↓ sectors. Something external has to break that symmetry
before "synchronization" can mean anything selective.

**A concrete, computable handle:** turn on a uniform field **B** along z.
The non-relativistic (Pauli/Foldy–Wouthuysen) reduction of the Dirac
equation gives a Zeeman splitting of the two channels' energies by
∓μ_B B, which detunes the two clocks:

$$ \omega_\uparrow \approx \omega_C + \frac{2\mu_B B}{\hbar}, \qquad
   \omega_\downarrow \approx \omega_C - \frac{2\mu_B B}{\hbar} $$

*(This should be re-derived carefully from your FW reduction rather than
taken on faith — the sign convention and the factor of 2 depend on exactly
how "clock frequency" is defined relative to the beat term.)*

This gives the model something real to lock to: two nearby but distinct
frequencies, rather than a purely formal distinction between channels.

---

### 3. Full spin superposition as a two-clock superposition

The ordinary spin state:

$$ |\chi\rangle = \alpha|\uparrow\rangle + \beta|\downarrow\rangle, \qquad |\alpha|^2+|\beta|^2=1 $$

is represented, in this picture, as:

$$ \Psi(t) = \alpha\, \Phi_\uparrow(t) + \beta\, \Phi_\downarrow(t) $$

i.e. a superposition not just of two spin labels but of two independently
beating pos/neg-energy clock pairs, each running at its own (possibly
Zeeman-shifted) frequency.

---

### 4. Observer clock and synchronization-selection (Conjecture B)

Introduce the observer's own internal phase θ_O(t), running at natural
frequency ω_O, coupled to *both* channel clocks. Reusing Kuramoto-style
phase coupling — not for the system's own free evolution (which you've
moved away from) but specifically at the **observer–system interface**:

$$ \dot\theta_O = \omega_O + \kappa_\uparrow \sin(\varphi_\uparrow - \theta_O) + \kappa_\downarrow \sin(\varphi_\downarrow - \theta_O) $$

This is the Adler-equation structure from your earlier work, redeployed
as a measurement-coupling rather than an internal-dynamics assumption.

**Working claim:** if the coupling is such that θ_O locks into phase with
φ_↑ (or φ_↓), the observer registers a definite spin value. If ω_O sits
outside the locking range of both — or the coupling strengths κ_↑, κ_↓ are
too weak or too symmetric to produce a lock — the observer's clock remains
unsynchronized to either channel, and (per the original proposal) that
observer perceives superposition rather than a collapsed outcome.

---

### 5. The actual target: Born weights from locking-basin size

This is the open problem the rest of this file exists to set up cleanly,
not something derived yet.

**Claim to test:** the fraction of (κ_↑, κ_↓, ω_O, initial-θ_O) parameter
space for which θ_O locks to φ_↑ should equal |α|², and likewise the
↓-locking basin should equal |β|², for the amplitudes appearing in Ψ(t)
in Section 3.

This is structurally the same kind of claim as the basin-of-attraction
arguments in your earlier Kuramoto work — so if it holds here, it's a real
continuity, not just borrowed language.

**The hard, unresolved part:** what physical process fixes κ_↑ and κ_↓ in
the first place. Without an independent (non-circular) answer to that, the
model can be tuned to reproduce |α|², |β|² after the fact, which is curve
fitting, not a derivation. This is the crux question for the next stage of
work, not a detail to defer indefinitely.

---

### 6. Next steps for a verification script

1. Re-derive ω_↑, ω_↓ (Section 2) properly from the Foldy–Wouthuysen or
   Pauli reduction in a field, rather than assuming the leading-order form.
2. Numerically integrate the coupled system (θ_O, φ_↑, φ_↓) over a grid of
   (κ_↑, κ_↓, ω_O, initial θ_O) and classify each trajectory as locked-↑,
   locked-↓, or unlocked.
3. Compute locking-basin fractions and compare against |α|², |β|² for
   several test values of (α, β) — start with the symmetric case
   α = β = 1/√2 as the sharpest test, since basin symmetry should be
   easiest to check or falsify there.
4. Only after (1)–(3) hold up: attempt an independent physical argument for
   what sets κ_↑, κ_↓ (candidates: coupling strength to the field that
   produced the Zeeman splitting; thermal/environmental coupling magnitude;
   something derived from the measuring apparatus's own internal clock
   structure). This is where the derivation either becomes real or the
   model needs to be revised.

---

### 7. Proposed experimental realization (working choice)

**Microscopic source: Cooper-pair splitting.** A superconductor's Cooper
pairs are spin-singlet entangled electron pairs by construction. Splitting
them via Andreev reflection at a superconductor–normal-conductor interface
(e.g. into two quantum dots on either side, as in carbon-nanotube Cooper-pair
beamsplitters) gives a physically concrete singlet source: two electrons
with a definite, shared phase relationship prior to separation. This is the
natural real-world instantiation of Ψ(t) in Section 3 — the singlet is the
α = 1/√2, β = −1/√2 (or equivalent sign convention) case, making it the
sharpest test case flagged in Section 6, step 3.

**Spin readout: energy-selective spin-to-charge conversion (E-RO), read out
via a QPC electrometer.** After splitting, each electron sits in its own
quantum dot. A magnetic field splits the local spin-↑/spin-↓ levels by the
Zeeman energy, and the dot's chemical potential is tuned so that only the
higher-energy spin state can tunnel out to the reservoir; whether that
tunneling event occurs is registered by a quantum point contact placed next
to the dot, whose conductance is highly sensitive to the dot's charge state.
This choice is what ties the detector directly into the existing framework:
**the E-RO threshold is set relative to the same Zeeman-split frequencies
ω↑, ω↓ from Section 2** — so the detector's own tuning literally is a choice
of which system frequency it is built to resolve. This is the most direct
candidate for "the observer clock" of Section 4 to be identified with a real
piece of apparatus, rather than a free parameter.

**Pair-level entanglement witness: Pauli spin blockade**, as a cross-check
independent of the single-dot E-RO readout. Whether the two electrons can
jointly tunnel between the two dots depends on whether their combined state
is singlet or triplet, so blockade behavior gives a direct entanglement
signature without needing individual spin resolution first — useful as an
independent verification that the source is producing genuine singlets
before trusting per-electron E-RO statistics against the |α|², |β|²
prediction of Section 5.

**Open item carried forward:** E-RO is a destructive, single-shot
measurement (the electron leaves the dot), which is the correct kind of
measurement to test against for this proposal, but it forecloses repeated
weak measurement on the same pair. If a repeated-measurement or weak-locking
version of this experiment is wanted later, non-destructive Kerr-rotation
readout is the natural substitute to revisit — noted here rather than
adopted, since E-RO is the more direct match to the "definite outcome vs.
superposition" framing of Section 4.

---

### 8. Numerical verification: results so far

Scripts: `basin_locking_sim.py` (deterministic), `basin_locking_noisy_sim.py`
(observer clock with added noise). Both implement the model of Section 4
with the ansatz of Section 6 step 4 under test:

$$ \kappa_\uparrow = \kappa_0|\alpha|, \qquad \kappa_\downarrow = \kappa_0|\beta| $$

**8.1 Deterministic result — FAILS.** Sweeping the observer's initial phase
θ(0) does not produce basin fractions matching |α|², |β|². Instead the
system behaves as a classic two-tone injection-locking oscillator: whichever
channel has even slightly larger coupling captures essentially *every*
initial phase (a sharp threshold near α² = β² = 0.5, jumping directly from
0% to 100% locked, with no graded region in between). Randomizing θ(0) is
not the right source of randomness for this dynamical system — there is no
multistability of the kind the original Kuramoto basin-of-attraction
argument relied on. Verified: pure states (α=1,β=0 and α=0,β=1) lock
correctly to their respective channel in both scripts, so the coupling and
classification machinery itself is working as intended.

**8.2 Adding observer-clock noise — partial, and reveals a second problem.**
A Langevin term (Euler–Maruyama, noise strength σ) does turn the sharp
threshold into a smoothed transition rather than a hard step. But at
σ = 0.20 (κ₀ = 0.20, Zeeman split = 0.15, same units), the resulting curve
still does not match |α|²: pure-state endpoints (α²=0, α²=1) are recovered
correctly, but the transition is steep and off-center, with a large
"unlocked" fraction dominating near the symmetric point rather than a clean
graded split (e.g. at α²=0.6: P(up)≈0.07, P(down)≈0, unlocked≈0.93 — nowhere
near the target P(up)=0.6).

**8.3 Deeper issue found: "unlocked" is not a stable outcome, it's a
window-length artifact.** Extending the total integration time T_MAX at
fixed α²=0.6 does not shrink the unlocked fraction toward the target split —
it *grows* (0.97 → 0.99 → 1.00 as T_MAX increased from 2500 to 12000, same
units as before). Diagnosis: noise causes intermittent phase slips — a
trajectory can reside near ω↑ for a while, get kicked loose, drift, reside
near ω↓, and so on. The classification method used (net average frequency
over the whole run) averages over these slips and, given long enough
windows, converges toward neither target frequency. This means "unlocked,"
as currently defined, measures how long you watched rather than a genuine
third physical outcome.

**8.4 Proposed fix, not yet implemented.** Redefine "locked" as the
*fraction of time* spent near each channel's frequency within a finite
window, and set that window's length to a physically motivated timescale
rather than an arbitrary simulation duration — specifically, the detector's
own response/integration time from the E-RO readout scheme in Section 7.
This turns a numerical nuisance parameter into a physical one: the
classification window becomes "how long the actual apparatus integrates
before registering a tunneling event," not an arbitrary choice.

**Status:** the κ ∝ |amplitude| ansatz has not yet produced a match to Born
statistics under either the deterministic or the naively-noisy version of
the model. Next step is the windowed-classification rebuild described in
8.4, not yet built.

---

### 9. Windowed classification: results

Script: `basin_locking_windowed_sim.py`. Implements the fix from 8.4:
classification now uses the average frequency over a finite window T_det
following a settling period, rather than a net average over the whole run.

**9.1 Window-length trade-off, confirmed empirically.** Scanning T_det at
the symmetric point (α²=β²=0.5) shows a genuine sweet spot near
**T_det ≈ 1 × the natural resolving time** (2π/Zeeman split). Shorter
windows can't yet distinguish ω↑ from ω↓ at all (ambiguous fraction stays
high, ~0.65–0.80). Longer windows give noise time to produce a phase slip
mid-measurement, so ambiguous grows again past that point too. Both
constraints from Section 8.4 turned out to matter — the "too long" direction
in particular contradicted the original guess (8 × resolving time), which
was corrected after running the scan.

**9.2 At the optimal window, the symmetric case gives a genuine three-way
split**, not a 50/50 up/down split: P(up) ≈ 0.33, P(down) ≈ 0.33,
ambiguous ≈ 0.34. Conditioning on a definite outcome (excluding ambiguous
trials) recovers P(up|resolved) ≈ 0.5 — a clean match to Born statistics at
this one point.

**9.3 Away from the symmetric point, the conditioned curve is a steep
sigmoid, not the required straight line.** Full sweep at T_det = 1×resolving
time, κ ansatz unchanged (κ↑=κ₀|α|, κ↓=κ₀|β|):

| α² | predicted | P(up) | P(down) | ambiguous | P(up\|resolved) |
|---|---|---|---|---|---|
| 0.00 | 0.00 | 0.000 | 0.993 | 0.007 | 0.00 |
| 0.25 | 0.25 | 0.028 | 0.860 | 0.113 | 0.03 |
| 0.40 | 0.40 | 0.170 | 0.537 | 0.292 | 0.24 |
| 0.50 | 0.50 | 0.333 | 0.333 | 0.335 | 0.50 |
| 0.60 | 0.60 | 0.527 | 0.117 | 0.355 | 0.82 |
| 0.75 | 0.75 | 0.880 | 0.022 | 0.098 | 0.98 |
| 1.00 | 1.00 | 0.988 | 0.000 | 0.013 | 1.00 |

Matches only at the two endpoints and the midpoint; overshoots badly
everywhere else (α²=0.6 predicts 0.60, resolves to 0.82). This is the
signature of too much positive feedback in the locking competition: small
coupling asymmetries are amplified into large outcome-probability swings.
The windowing fix worked methodologically (turned a meaningless, growing
"unlocked" fraction into a stable, interpretable three-way outcome); it also
cleanly exposed that the κ(α,β) ansatz has the wrong functional shape.

---

### 10. Isolating the sigmoid: intrinsic response curve

**10.1 Ruling out the varying-total-coupling explanation.** The Section 6
ansatz lets total coupling κ↑+κ↓ vary with (α,β) — largest near the
symmetric point (κ₀√2), smallest at the pure states (κ₀) — a second
uncontrolled variable tangled up with the sigmoid. To isolate this, the
total coupling budget was fixed at κ_total = κ₀ and only its *split* was
varied: κ↑ = κ_total·x, κ↓ = κ_total·(1−x), x ∈ [0,1]. Result: **the sigmoid
persists, and is if anything sharper** (P(up|resolved) goes from 0.10 at
x=0.4 to 0.98 at x=0.7). This rules out the varying-total-coupling ansatz as
the cause — the steep sigmoidal response is intrinsic to the competitive
two-tone locking mechanism itself, not an artifact of that particular
choice.

**10.2 Noise-strength sweep: an approximately linear regime exists, at a
cost.** Repeating the fixed-budget response curve across observer-noise
strengths σ = 0.2, 0.5, 1.0, 2.0 (κ_total = 0.2 throughout):

- σ = 0.2 (Section 9's value): sharp sigmoid, as above.
- **σ = 0.5: approximately linear**, but with slope compressed to roughly
  0.68 rather than 1 (P(up|resolved) ≈ 0.15, 0.24, 0.35, 0.52, 0.74, 0.83 at
  x = 0, 0.2, 0.4, 0.6, 0.8, 1.0 — no longer sigmoidal in shape, but shifted
  toward 0.5 at both ends rather than reaching 0 and 1). Ambiguous fraction
  is large here, roughly 0.39–0.48 across the sweep.
- σ = 1.0: response is washed out almost flat (~0.50 ± 0.02 regardless of
  x) — noise dominates the coupling asymmetry entirely, no discriminating
  information left.
- σ = 2.0: fully washed out, and ambiguous fraction is very large (~0.86).

**Interpretation:** there is a real transition from sigmoid (low σ) through
an approximately linear regime (σ ≈ 0.5) to a fully noise-dominated flat
response (σ ≳ 1). The σ ≈ 0.5 regime is the first result in this
investigation that is linear in shape rather than merely tuned at isolated
points — but its slope is below 1, and its ambiguous fraction is large
(~40–50%).

**Suggestive, unconfirmed note:** that ambiguous fraction sits in roughly
the same range as the ~65% single-shot readout visibility reported for real
quantum-dot E-RO spin measurements (i.e. ~35% infidelity) — recorded here as
a coincidence worth keeping in view, not as evidence of anything yet.

**Open items for the next stage:**
1. Fine-tune (σ, κ_total) jointly to find where the response slope is
   exactly 1, not just approximately linear — the σ=0.5 point above was a
   coarse grid sample, not a fit.
2. Once a slope-1 combination is found, check whether its ambiguous
   fraction is a fixed constant (independent of α,β) or itself varies with
   α,β — Born's rule requires |α|²+|β|²=1 to exhaustively cover the
   probability, so a α,β-dependent ambiguous fraction would need to be
   folded back into a renormalization, not treated as a separate leak.
3. Revisit whether the ambiguous-fraction/readout-infidelity echo in 10.2
   is worth pursuing as a real correspondence, or is coincidental at this
   sample size.

---

### 11. Searching for the exact slope-1 point: a genuine impasse

**11.1 A slope-1 fit was found, but it is a least-squares statement, not a
pointwise one.** Fine-scanning σ at fixed κ_total=0.2, T_det=1×resolving
time locates **σ ≈ 0.40** as the point where a linear regression of
P(up|resolved) against the coupling fraction x gives slope = 0.998,
intercept = −0.004 (smallest deviation from the identity line in the scan).

**11.2 That fit does not survive contact with the actual physical ansatz.**
Plugging in the calibration this implies — κ↑ = κ_total·α², κ↓ =
κ_total·β² (coupling proportional to *probability*, not amplitude; see 11.5
for why this itself is a separate concern) — the middle of the curve tracks
α² reasonably (α²=0.6 → P(up|resolved)≈0.66, α²=0.75 → ≈0.84), but the pure
states no longer resolve with near-certainty: at α²=1.0 (β=0, zero coupling
to the down channel at all), P(up) is only ≈0.67–0.69, not ≈0.98 as at lower
noise. A definite state must give a definite, certain outcome — this is
non-negotiable, more fundamental than a straight middle section, and it is
failing at the noise level the linear fit required.

**11.3 Scaling κ_total and σ together (fixed ratio) makes it worse, not
better.** Tested at κ_total/σ = (0.2, 0.4), (0.4, 0.8), (0.8, 1.6), (1.6,
3.2), holding T_det fixed at 1× the Zeeman resolving time throughout:
ambiguous fraction climbs monotonically with scale (≈30% → 52% → 80% → 91%)
instead of improving. Diagnosis: the locking timescale set by κ_total
shrinks as κ_total grows, so a T_det fixed at the Zeeman-based resolving
time increasingly spans many locking cycles at larger κ — more time within
one measurement window for noise to produce a phase slip. T_det cannot be
held fixed while κ_total, σ scale; it has to be retuned jointly with them.

**11.4 Retuning T_det alone at fixed (κ_total=0.2, σ=0.40) does not recover
it either.** Scanned T_det from 0.1× to 1.0× the resolving time: shorter
windows fail via the original resolving-time constraint (ambiguous climbs
back to 76–83% at 0.1×), and T_det = 1× resolving time remains the best
point in this set, at which pure-state fidelity is still only ≈0.68–0.93,
not the required ≈1.

**11.5 Conclusion — this is a genuine impasse under the current coupling
form, not an unfinished search.** All three available knobs were checked —
σ alone, σ and κ_total together, T_det alone — and each one traded one
failure mode for the other: whatever softens the middle of the curve into a
line also erodes certainty at the pure-state edges, and whatever restores
edge certainty steepens the middle back toward the original sigmoid. No
single (κ_total, σ, T_det) triple found satisfies both requirements
simultaneously under the Adler-type sin(φ−θ) coupling used throughout
Sections 4–10. This points at the coupling *functional form itself* — not
further parameter tuning — as the next thing to revisit: sin-coupling is
what produces both the winner-take-all competition (needed for endpoint
certainty) and the noise-softened sigmoid (needed for a graded middle), and
it appears unable to deliver both at once.

**11.6 Separate, still-open concern: the κ∝α² ansatz risks circularity.**
Even where the fit above partially worked, it required setting the coupling
strength proportional to the Born probability itself (α², not |α|), which
is what the model is supposed to derive, not assume. This is not
automatically illegitimate — physical couplings arising from an interaction
Hamiltonian's expectation value are generically quadratic in state
amplitudes, so a ⟨Ψ|H_int|Ψ⟩-type origin for this scaling is plausible in
principle — but as used here it was fit to match the target, not derived
from an explicit interaction term. This needs to be addressed independently
of 11.5 even if a coupling form is found that resolves the impasse.

**Next step, not yet started:** reconsider the coupling functional form
itself (Section 4) — e.g. a coupling that saturates rather than following a
pure sine, or two structurally different mechanisms for "capture" (drives
certainty at the pure-state limit) versus "competition" (drives the graded
middle) — rather than continuing to search the (κ, σ, T_det) space of the
current form.

---

### 12. Replacing the coupling form: linear response + separate race process

Script: `linear_response_race_sim.py`. Motivated directly by 11.5/11.6, and
by the parallel to abandoning the analogous Kuramoto sine-coupling in the
Born-selection/measurement work in favor of the linear Dirac spinor
(mass-term) coupling.

**12.1 Structural diagnosis.** sin(φ−θ) is nonlinear in the phase
difference. When two such couplings act on one observer simultaneously, the
nonlinear cross-term makes them *compete* — this is the direct cause of the
winner-take-all thresholding found throughout Sections 8–11. A *linear*
coupling has no such cross-term: two simultaneous linear drives at
different frequencies contribute additively to a linear system's
steady-state response, since their cross-term time-averages to zero for
different frequencies. This is the same structural mechanism behind Fermi's
Golden Rule and other amplitude-squared results in physics — it does not
require assuming Born's rule to get amplitude-squared scaling.

**12.2 New model.** Replace the Adler equation entirely with a damped,
linearly driven oscillator for the observer:

$$ \dot z = -i\omega_O z - \gamma z + \kappa_\uparrow e^{-i\omega_\uparrow t} + \kappa_\downarrow e^{-i\omega_\downarrow t} $$

Closed-form steady-state power at each drive frequency:

$$ P_\uparrow = \frac{\kappa_\uparrow^2}{(\omega_O-\omega_\uparrow)^2+\gamma^2}, \qquad P_\downarrow = \frac{\kappa_\downarrow^2}{(\omega_O-\omega_\downarrow)^2+\gamma^2} $$

With ω_O set to the midpoint (ω↑+ω↓)/2, the two denominators are equal, so
**P↑ : P↓ = κ↑² : κ↓² exactly**. Using κ↑ = κ₀|α|, κ↓ = κ₀|β| (amplitude,
not probability) is now non-circular: the square comes from power in a
linear response, a well-established mechanism, not from the ansatz.

**12.3 Outcome generation is now a separate, explicit layer.** Rather than
folding "definite outcome" and "graded competition" into one dynamical
equation (which Section 11 showed cannot do both), outcome generation is a
race between two independent Poisson click processes with rates
proportional to P↑, P↓. Whichever fires first within a detector window
T_window is the outcome; if neither fires, the trial is genuinely
inconclusive — a separate, independently tunable parameter (the overall
rate scale) rather than an artifact entangled with the coupling dynamics.

**12.4 Numerical verification.** A real discretization bug was caught and
fixed here: the first integration pass used forward-Euler stepping, which
showed a systematic ~5–8% amplitude bias against the closed-form formula.
Confirmed as a first-order discretization artifact (bias shrank
proportionally with step size: 7.8% at dt=0.01, 2.25% at dt=0.003, 0.74% at
dt=0.001), then eliminated properly with an exact exponential-integrator
step (closed-form solution of the linear ODE per step, not an approximation
regardless of step size). After the fix, numeric and analytic power agree
to 4+ decimal places (e.g. 0.64003 vs 0.64000 at α²=0.25).

**12.5 Result: matches Born statistics across the full range, cleanly.**

| α² | predicted | P(up) | P(down) | inconclusive |
|---|---|---|---|---|
| 0.10 | 0.10 | 0.097 | 0.903 | 0.000 |
| 0.25 | 0.25 | 0.236 | 0.763 | 0.000 |
| 0.40 | 0.40 | 0.387 | 0.613 | 0.000 |
| 0.50 | 0.50 | 0.492 | 0.508 | 0.000 |
| 0.60 | 0.60 | 0.605 | 0.395 | 0.000 |
| 0.75 | 0.75 | 0.754 | 0.246 | 0.000 |
| 0.90 | 0.90 | 0.902 | 0.098 | 0.000 |

Deviations are consistent with Monte Carlo sampling noise at 2000 trials per
point (expected std ≈ 0.01 near the middle), not a systematic shape error —
this is the first result in the whole investigation with the *right shape*
across the entire range, not a locally-tuned point.

**12.6 Open items, not yet resolved:**
1. **Why ω_O = midpoint?** Tuning the observer's natural frequency exactly
   to (ω↑+ω↓)/2 is what makes the two resonance denominators equal. This
   was chosen by design here, not derived — needs a physical reason the
   observer/detector would sit at that specific frequency rather than
   anywhere else.
2. **κ ∝ |amplitude| still needs a derivation, not just non-circularity.**
   12.1–12.2 show the ansatz is no longer circular in the way the Section
   11.6 concern flagged, but "non-circular" is not the same as "derived." A
   concrete next step: obtain κ↑, κ↓ from an explicit interaction
   Hamiltonian term (something like ⟨Ψ|H_int|Ψ⟩) rather than positing the
   proportionality directly.
3. **γ (observer damping) and the rate-scale/T_window pair are free
   parameters with no physical origin yet assigned.** Section 7's proposed
   apparatus (Cooper-pair source + E-RO/QPC readout) should be revisited to
   see whether any of these three map onto real, measurable properties of
   that setup (γ ↔ some decoherence rate; rate-scale/T_window ↔ the E-RO
   tunneling rate and detector integration time) rather than remaining
   free-floating simulation knobs.
4. Demonstrate the inconclusive fraction is a genuinely independent knob by
   rerunning at a smaller rate-scale·T_window product and checking that the
   resolved-outcome ratio still tracks α² even as the inconclusive fraction
   grows — not yet done, but the model's structure (12.3) predicts it
   should hold.

---

### 13. Parallel branch: bilocal Dirac current instead of an observer clock

**This is a separate line of investigation running alongside Sections
4–12, not a replacement for them.** Sections 4–12 keep an explicit observer
clock racing against two system clocks in some global time. This branch
asks what happens if "observer" is replaced entirely by the two interacting
particles themselves, with no third clock and no global time-ordering
assumed at all.

**13.0 Motivating critique.** The Section 12 linear-response fix corrected
the *coupling* (linear, so no competition between channels) but the
*readout* still used |amplitude|² — a modulus, which is orientation-blind
and discards exactly the relative-phase information a real projection
retains. The Dirac equation already has the right kind of real, projected
object built in: the conserved probability current ψ†γ^μψ. This branch
replaces both the observer-clock coupling and the modulus-based readout
with constructions built from that current instead.

**13.1 Single-point transition current (starting point).** For two Dirac
spinor solutions ψ₁, ψ₂ of the same mass, evaluated at the same spacetime
point, define

$$ J^\mu_{\text{sync}} \equiv \mathrm{Re}[\bar\psi_1\gamma^\mu\psi_2] = \tfrac12(\bar\psi_1\gamma^\mu\psi_2 + \bar\psi_2\gamma^\mu\psi_1) $$

Real by construction (via $(\bar\psi_1\gamma^\mu\psi_2)^*=\bar\psi_2\gamma^\mu\psi_1$),
and conserved by the same cancellation argument as the ordinary Dirac
current: using $i\gamma^\mu\partial_\mu\psi_2=(mc/\hbar)\psi_2$ and the
adjoint relation for $\bar\psi_1$, the mass terms cancel exactly:

$$ \partial_\mu J^\mu_{\text{sync}} = 0 $$

**13.2 Extension to genuinely separated particles.** The single-point
construction breaks down once the two electrons separate after Cooper-pair
splitting (Section 7) — $\psi_1$ and $\psi_2$ live at different spacetime
points. Resolved via standard two-particle QFT machinery rather than
forcing the single-point object to apply: let $\Psi(x_1,x_2)$ be the joint
bispinor obeying the free two-particle Dirac equation (tensor product of
each particle's own Dirac operator, acting on its own coordinate), and
define the bilocal current

$$ J^{\mu\nu}_{\text{sync}}(x_1,x_2) \equiv \mathrm{Re}\big[\bar\Psi(x_1,x_2)(\gamma^\mu\otimes\gamma^\nu)\Psi(x_1,x_2)\big] $$

well-defined at $x_1\ne x_2$ since it is a correlation function, not a
local density. It satisfies two separate conservation laws, one per
particle's coordinate:

$$ \partial_{1\mu}J^{\mu\nu}_{\text{sync}}=0, \qquad \partial_{2\nu}J^{\mu\nu}_{\text{sync}}=0 $$

**13.3 This is the precise version of "no disagreement between the
particles."** By Gauss's law on the first index alone, $Q^\nu\equiv\int
d^3x_1\,J^{0\nu}_{\text{sync}}(x_1,x_2)$ is the same on any spacelike slice
through particle 1's worldline, regardless of which frame drew that slice —
and likewise for particle 2 via the second index. Different outside frames
can disagree about which detection event happened first; $Q^\nu$ does not
depend on that choice, because it is built from a 4-current conserved
separately in each argument. This is standard QFT machinery applied to make
the synchronization claim rigorous, not a new conjecture.

**13.4 Detector dispersal, in the same language (qualitative so far).** As
particle 1 enters the detector and interacts with detector constituents
$a_1, a_2, \dots$, each interaction is governed by the same pairwise rule,
generating new currents $J^{\mu\nu}_{1,a_1}$, then correlations extending to
$a_2$, and so on. The original 1–2 correlation is not destroyed (nothing
here is non-conserving) but becomes encoded across a rapidly growing joint
object involving many detector particles. Any current an experimenter can
access locally only involves a handful of degrees of freedom — recovering
the original clean phase relationship would require access to the full
many-body object, which is not available. This is decoherence, reconstructed
in the current formalism instead of the density-matrix formalism — a
translation, not new physics.

**13.5 What this branch does and does not achieve, stated plainly.** It
gives a rigorous, frame-independent account of *why* the two particles'
outcomes are correlated with no disagreement about ordering (13.3), and a
believable qualitative mechanism for why detector outcomes look definite
and classical (13.4). It does **not** derive the Born-rule weights |α|²,
|β|² — 13.3 explains why the correlation is frame-independent, not what its
statistics are; 13.4 explains apparent definiteness, not why one channel is
favored in the right proportion. This is the same gap every serious
decoherence program runs into (it is why Zurek's program needed envariance
as a separate argument on top of ordinary decoherence) — relocated to a
more rigorous setting here, not yet closed.

**Next step in this branch, not yet started:** build a concrete toy model
of the dispersal step in 13.4 — e.g. a short chain of sequential pairwise
interactions (particle 1 with detector particle $a_1$, then $a_1$ with
$a_2$, etc.), tracking numerically how quickly the original two-particle
current becomes unrecoverable from any local subset as the number of
participating detector particles grows.

---

### 13.6 Toy dispersal-chain model

Script: `dispersal_chain_toy.py`. Makes 13.4's qualitative dispersal claim
quantitative and checkable, within an explicitly simplified toy setting.

**Simplification, stated plainly.** This is not a bispinor/QFT simulation
— not tractable to run directly for a many-body detector chain. Qubits
stand in for the two-level (up/down) degree of freedom, and the pairwise
interaction is a hopping-type ("XY") gate, chosen for one specific
structural reason: it conserves total excitation number exactly, the
discrete analog of the conserved charge $Q^\nu$ from 13.3. That is the one
feature of the real physics this toy is built to preserve faithfully. The
spinor structure and the actual Dirac current construction are set aside
here — this checks the *logic* of 13.4 ("dispersed, not destroyed"), not
the specific construction of 13.1–13.3, and the connection between them is
an analogy, not a derivation.

**Model.** Particles 1, 2 (qubits 0, 1) start in a singlet. Particle 1
interacts sequentially, one at a time, with detector qubits $a_1,\dots,
a_k$ (qubits 2, 3, ...), via the same XY gate each time (coupling angle θ),
modeling it entering a detector made of many constituents.

**Result 1 — nothing is destroyed, to machine precision.** Total
excitation (summed over *all* qubits, the full system) stays at exactly
1.000000000000 through all 8 sequential interactions, for every θ tested —
confirming the toy model's one required structural feature actually holds
in the simulation, not just by construction on paper.

**Result 2 — local recoverability decays, and by an exact closed-form
law.** The reduced two-particle correlation, tracked two ways
(coherence $|\rho_{12}[01,10]|$ and concurrence), decays geometrically with
the number of detector interactions $k$. Concurrence in particular follows

$$ C(k) = \cos^k(\theta) $$

confirmed exactly against simulation at $k=8$ across four coupling
strengths (θ = 30°, 45°, 60°, 75° → predicted 0.316406, 0.062500, 0.003906,
0.000020 — matching simulated values to displayed precision in every case).
Any nonzero per-step coupling (θ < 90°) drives the original 1–2 correlation
to zero as more detector constituents participate — this is the concrete,
quantitative version of "encoded across a growing joint object... no longer
recoverable from any local subset."

**Caveat on purity.** Purity of the reduced (1,2) state does *not*
monotonically approach 0 — it plateaus near 0.5. This is not evidence
against dispersal; purity alone isn't the right diagnostic here, since a
mixed two-qubit state can retain substantial purity while holding zero
entanglement or coherence with its original partner. Coherence and
concurrence are the correct diagnostics for this claim, and both cleanly
vanish.

**What this does and does not establish.** It gives a concrete, checkable
demonstration that a conserved global quantity can coexist with exact
exponential loss of local recoverability — the specific logical structure
13.4 needed. It does **not** bridge to the actual Dirac bilocal current of
13.1–13.3 (that connection remains an analogy), and — consistent with
13.5 — it says nothing about Born-rule statistics; this branch still
doesn't touch that question.

---

### 13.7 Testing the toy model against a genuine two-particle Dirac equation

Script: `dirac_dispersal_chain.py`. Directly checks whether the qubit toy's
cos^k(θ) law survives when the interaction is built from actual Dirac
bispinor structure rather than an abstract two-level analog.

**13.7.0 A modeling choice made explicit.** 13.1–13.3 establish a conserved
*current*, not a specific interaction *Hamiltonian* for coupling two
separate particles. Two candidates were built to bracket the question:

- **(A) Spin-exchange**, built from the Dirac spin operator Σ = diag(σ,σ),
  block-diagonal in the positive/negative-energy (u/v) index — never mixes
  u and v.
- **(B) γ-matrix coupling**, built from the spatial matrices γ¹, γ² —
  off-diagonal in the standard Dirac representation, genuinely mixes u and v.

**13.7.1 (A) reduces exactly to the qubit toy, as required.** Starting
purely in the u-sector, Σ cannot populate the v-sector at all, so the
dynamics stay confined to an effective two-level subspace per particle —
concurrence matches cos^k(θ) to 6 decimal places at every k tested. This
confirms the toy model is not an arbitrary abstraction; it is the exact
restriction of a genuine (if spin-only) Dirac coupling to the positive-energy
sector.

**13.7.2 (B) does NOT reproduce cos^k(θ) — genuinely new physics appears.**
At θ=60° (matching 13.6's example): concurrence(k=1) = 0.800, not the
predicted 0.500, and 37.5% of the population has already leaked into the
negative-energy sector after a *single* interaction — leakage that has no
analog in the qubit toy at all. Leakage then saturates toward 50% as k
grows, and concurrence within the surviving u-sector population does
eventually settle into roughly the same asymptotic per-step decay rate as
(A) — but only after this large initial transient the toy model completely
missed.

**13.7.3 A θ→90°−θ symmetry, confirmed exactly, and explained.** θ=30°
and θ=60° produce *identical* results at every k (not approximately —
matching to 6 decimals throughout); confirmed across four further angle
pairs (10°/80°, 20°/70°, 35°/55°, 40°/50°), always exact. Diagnosed via the
interaction Hamiltonian's spectrum: eigenvalues are exactly {−2, 0, +2}
(multiplicities 4, 8, 4) — a clean discrete spectrum. Dynamics therefore
depend on θ only through phases e^{±2iθ}, so observable probabilities
depend on trigonometric functions of 2θ, which are exactly invariant under
θ→90°−θ. This has no analog in the spin-exchange/qubit picture, whose decay
depends on cos(θ) directly, with no such symmetry.

**13.7.4 A qualitatively new phenomenon: instantaneous complete dispersal
at θ=45°.** At exactly θ=45° (2θ=90°, a Rabi π/2-pulse point), a *single*
interaction drives concurrence to exactly zero and leakage to exactly 50%
— and **every subsequent interaction does nothing at all**: the state has
hit an exact fixed point of the single-particle-plus-fresh-environment map.
This is fundamentally different from the toy model's gradual geometric
decay: here, dispersal can complete in one step rather than accumulating
over many, depending sensitively on coupling strength.

**13.7.5 The free-evolution check found no effect — but this is a
diagnostic limitation, not a physics finding.** Inserting free evolution
(the Compton-frequency phase rotation) between interactions changed nothing
in the concurrence/leakage numbers. Reason: the rest-frame free Hamiltonian
mc²γ⁰ is spin-degenerate (u↑, u↓ share the same eigenvalue, as do v↑, v↓),
so it acts as a pure global phase on any state confined to a single u/v
sector — invisible to a concurrence diagnostic computed within that sector.
Free evolution *could* still matter for phase coherence between the u and v
sectors once leakage has occurred, but the diagnostic used here (projection
onto the u-sector block) cannot see that. **Open item:** redesign the
diagnostic to track u–v cross-coherences directly before concluding timing
is irrelevant.

**13.7.6 Bottom line.** The qubit toy's cos^k(θ) law is the correct
description only for the spin-exchange-restricted, positive-energy-confined
special case (13.7.1) — not for a genuinely faithful bispinor coupling.
The real Dirac structure introduces immediate (not gradual) negative-energy
leakage, a θ→90°−θ symmetry with a clean spectral explanation, and a
special angle where dispersal completes in a single interaction. This is a
real answer, not a null result: the toy model in Section 13.6 captured the
right *logical structure* (conservation + local irrecoverability) but not
the right *functional form* once the actual bispinor degrees of freedom are
included.

---

### 13.8 A u-v coherence diagnostic, and the mass-drag hypothesis

Script: `uv_coherence_diagnostic.py`. Built to fix the limitation flagged
in 13.7.5 (the concurrence diagnostic is blind to u-v coherence), and to
test a specific physical claim from the Born-selection/MathX work: that the
Dirac mass term "drags" the v-component along with u dynamically, in a way
a massless particle (photon) has no counterpart for.

**Diagnostic design, deliberately a real projection, not a modulus.** For a
particle's reduced density matrix, define $c_{uv} \equiv \langle u\!\uparrow|\rho|v\!\uparrow\rangle + \langle u\!\downarrow|\rho|v\!\downarrow\rangle$
and track $\mathrm{Re}(c_{uv})$ and $\mathrm{Im}(c_{uv})$ *separately* —
the point being that $|c_{uv}|$ can stay perfectly constant while the phase
carries all the physically relevant information, exactly the failure mode
flagged earlier when the coupling model used a modulus-based readout.

**13.8.1 First attempt revealed a genuine, informative false start.**
Testing $c_{uv}$ on the leaked state from Section 13.7 (generated via the
entangling γ-interaction with a detector particle) gave $c_{uv}=0$
*exactly* — not a bug: the reduced single-particle density matrix came out
exactly diagonal, all coherence gone, not just between u and v but
entirely. Diagnosis: the γ-interaction is an *entangling* interaction, and
it pushed all of particle 0's local coherence into its correlation with the
detector particle rather than destroying it — the same "dispersed, not
destroyed, but locally unrecoverable" logic as Sections 13.4/13.6, now
shown to apply to the u-v degree of freedom as well. This meant the leaked
state was the wrong test case for isolating the mass-drag question
specifically, since entanglement alone already erases all local coherence
regardless of mass.

**13.8.2 Corrected test: an unentangled local superposition.** A single
particle prepared directly in $(|u\!\uparrow\rangle+|v\!\uparrow\rangle)/\sqrt2$,
with no other particles involved, isolates the question cleanly.

**13.8.3 Result: exact confirmation of the mass-drag hypothesis, closed
form.** For mc²=1 (electron-like): $c_{uv}(t)$ oscillates as detailed
below, while $|c_{uv}(t)|=0.5$ stays exactly constant throughout — the
modulus is *blind* to the entire effect, confirmed numerically, not just
argued. For mc²=0 (photon-proxy): $c_{uv}(t)=0.5$ for all t, frozen solid,
no oscillation at all.

The oscillation matches an exact closed form:

$$ c_{uv}(t) = \tfrac12 e^{-2imc^2t} $$

confirmed to 6 decimal places against the simulation at every t tested
(e.g. t=0.5: predicted Re=0.270151, Im=−0.420735; simulated: identical).
This is exactly the Compton frequency $\omega_C=2mc^2/\hbar$ established in
Section 2 — the mass term genuinely drives u and v apart in phase at
precisely that rate, with no analog when mc²=0.

**13.8.4 Interpretation, connecting back to the Born-selection/MathX
separation.** This gives a clean, quantitative version of "the second
spinor is dragged along in the electron, but not necessarily in the
photon": for a massive particle, any local u-v superposition is
*dynamically forced* to rotate at the Compton frequency purely by free
evolution — it cannot sit still. A massless particle has no such forcing
term at all. This supports treating the coupling/synchronization as
acting on the positive-energy component specifically (as in the MathX
work), with the negative-energy component's behavior for a massive
particle following automatically from the mass term itself, rather than
needing to be independently specified — while for a massless case, no such
automatic entrainment exists, and the two-component relationship (if any)
would have to be established some other way.

**Open items:** (1) this uses mc²=0 as a rough proxy for "photon-like" —
the photon is spin-1, not a massless Dirac fermion, so the analogy is
suggestive, not literal; a genuine treatment would need the actual photon
polarization/chiral-clock structure noted elsewhere in this research
program, not this Dirac-equation massless limit. (2) 13.8.1's finding — that
entanglement erases local u-v coherence just as it erases local spin
coherence — has not yet been connected back to the 13.7.5 open question of
whether free evolution *timing* affects the outcome of the dispersal chain;
that connection (does inserting free evolution between chain interactions,
now measured via $c_{uv}$ rather than concurrence, actually reveal a timing
dependence in the multi-particle chain) is still untested.

---

### 13.9 Settling 13.7.5: timing, entanglement structure, and genuine thermal noise

Script: `uv_coherence_diagnostic.py` (extended). This section fully closes
13.7.5, but the honest answer is more structural than originally expected,
and it required distinguishing two very different kinds of "noise."

**13.9.1 Why concurrence/leakage are provably timing-blind — not just
empirically null.** u↑ and u↓ share the identical γ⁰ eigenvalue (degenerate),
so free evolution restricted to either sector is *always* a pure global
phase. Any function of a density matrix (concurrence, populations) is
invariant under a global phase by construction. This is a proof, not an
observation: no choice of free-evolution time can ever change these
quantities, for this Hamiltonian, full stop.

**13.9.2 A second, independent reason, found by direct inspection.**
Printing the actual nonzero structure of the (particle0,particle1) reduced
density matrix after one γ-interaction showed the leaked branch
(particle0=v_down) has *exactly zero* coherence with either original
singlet branch. Mechanism: γ⁰ anticommutes with γ¹, γ² individually, but
γ⁰⊗γ⁰ *commutes* with γ¹⊗γ¹+γ²⊗γ² (two sign flips cancel) — making
γ⁰(particle0)⊗γ⁰(current detector) an exactly conserved quantity through
each interaction. Consequently, particle0 flipping u→v necessarily forces
the *current detector particle* to flip in lockstep — leaking and
entangling-with-the-detector are the same event, not two separable things.
Any coherence between leaked and unleaked branches is therefore carried
entirely by the (traced-out) detector correlation, invisible to *any*
diagnostic built from the (particle0,particle1) reduced state, regardless
of what phase is inserted between interactions. Two attempted fixes (a
phase kick on the fresh detector's pure state; a phase kick via inserted
free evolution on particle0) were both tried and both null, for exactly
this reason — confirmed, not merely observed, since 13.9.1 and this
mechanism both independently forbid it.

**13.9.3 Genuine thermal noise — modeled correctly this time — has a real,
substantial effect.** The two null results above are specific to *phase*
perturbations. Thermal disorder in a detector atom's collision strength
(from velocity/position spread at the moment of interaction) is physically
a different kind of randomness: it changes the actual coupling angle θ_k
per collision, not a phase on an already-fixed transition. Modeled as
θ_k ~ Normal(θ_mean, θ_std) per interaction, independently drawn per
realization, with the experimenter's state being the *ensemble average*
over unknown θ_k (200 realizations per point):

| θ_std | k=1 | k=2 | k=3 | k=4 |
|---|---|---|---|---|
| 0.00 | 0.800 | 0.471 | 0.246 | 0.125 |
| 0.15 | 0.736 | 0.416 | 0.215 | 0.110 |
| 0.30 | 0.595 | 0.309 | 0.169 | 0.095 |
| 0.50 | 0.403 | 0.168 | 0.136 | 0.095 |

Larger θ_std produces faster loss of recoverable concurrence at *every*
step — at k=1, θ_std=0.5 roughly halves concurrence relative to the
noiseless case. This is not subject to the phase-invariance arguments of
13.9.1–13.9.2: randomizing the actual transition strength changes real
population/coherence structure, not an unobservable global phase, so
averaging over unknown θ_k genuinely produces a more mixed state than any
single deterministic trajectory would. (Minor non-monotonicity at large k,
small θ_std — e.g. θ_std=0.5's k=5→k=6 uptick — is Monte Carlo noise at
200 realizations on an already-small concurrence value, not a real effect.)

**13.9.4 Resolution, stated plainly.** Timing between interactions and
phase-type noise are *provably* irrelevant to dispersal in this model —
two independent structural reasons, not an unlucky parameter choice.
Genuine thermal noise in the coupling strength itself is a categorically
different kind of perturbation, and it matters a great deal, in exactly the
direction physical intuition suggests: more thermal disorder means faster,
not slower, loss of locally recoverable correlation.

**13.9.5 A concrete quantum-computer experiment.** This maps directly onto
a buildable near-term circuit: a fixed two-qubit entangling gate
(approximating the γ-matrix interaction) applied repeatedly between one
"signal" qubit and a fresh "detector" qubit each round, with the gate angle
θ_k drawn independently at random each shot (classically randomized
compiling — directly supported on current hardware) rather than held fixed.
Running many shots at several noise levels and reconstructing the reduced
two-qubit state of the original signal pair via state tomography after
each round would directly measure the concurrence-decay curves above on
real hardware — and would additionally pick up the device's own intrinsic
thermal/decoherence noise as an uncontrolled but physically relevant bonus
contribution, which could itself be compared against the deliberately
injected θ-randomization to see whether they act similarly.

---

### 14. Detector geometry as the axis-selecting factor (new direction, placeholder pending Born-selection paper content)

**14.0 A clarification that motivated this section: which current is
complex.** The *ordinary* single-particle Dirac current $J^\mu=\bar\psi\gamma^\mu\psi$
is manifestly real (γ⁰γ^μ is Hermitian, so ψ†γ⁰γ^μψ is real for any ψ) — it
literally is the electron's own probability/charge flow, matching the
classical "flow of electrons" intuition directly, no complex quantity
involved. What *is* generically complex is the **transition current
between two different spinors**, $\bar\psi_1\gamma^\mu\psi_2$ — since
ψ₁≠ψ₂, the Hermiticity argument forcing reality no longer applies. This is
exactly why Section 13.1 took its real part explicitly:
$J^\mu_{\text{sync}}=\mathrm{Re}[\bar\psi_1\gamma^\mu\psi_2]$ is a real,
conserved quantity extracted from something that is not automatically real.
It is already real by construction — projecting it onto any axis below
stays real, with no further step needed.

**14.1 The proposal.** $\vec J_{\text{sync}}$ (the spatial part of
$J^\mu_{\text{sync}}$) is a real 3-vector, not a scalar — so it has a
*direction*, which can be projected onto whichever axis n̂ the detector
geometrically defines:

$$ J_{\text{sync}}(\hat n) \equiv \vec J_{\text{sync}} \cdot \hat n $$

exactly as ordinary spin measurement always means measuring $\vec\sigma\cdot\hat n$
for an experimenter-chosen axis, with different choices of n̂ giving
different up/down statistics for the *same* state. The proposal: n̂ is set
by the detector's geometry, and "up" vs "down" is defined relative to that
choice, not as an absolute property of the electron.

**14.2 This is not a new, free ingredient — it reuses an axis already in
the model.** The external B field from Section 2 already defines the
Zeeman-splitting axis (which energies count as ω↑ vs ω↓), and reappears as
the detector's field orientation in the Cooper-pair/E-RO proposal (Section
7). The proposal here is that this *same* geometric axis is what n̂ should
be — unifying "which axis splits the energies" and "which axis the
detector projects onto" into one geometric fact, rather than two
independent assumptions.

**14.3 What this does and does not yet solve.** Projecting onto n̂ answers
"up or down *relative to what*" — a genuinely missing piece before now. It
does **not** by itself answer how a continuous projected value becomes one
discrete click; that is the same gap Sections 5 and 12 were built to
address (via the race process), and it remains open here as a separate
question, not automatically solved by adding geometric structure.

**14.4 Proposed resolution: discretization via orbital resonance
(PLACEHOLDER — pending the actual Born-selection paper content).** The
working hypothesis is that the measurement process detailed in the latest
Born-selection paper (MathX) resolves 14.3: orbital resonances discretize
what would otherwise be a continuous distribution of $J_{\text{sync}}(\hat n)$
values, in a manner structurally analogous to how resonance/boundary
conditions produce discrete stable configurations elsewhere in physics
(e.g. Bohr-Sommerfeld quantization, normal-mode selection). **I do not have
the actual mechanism from that paper in front of me** — only the general
project summary (Many Clocks Interpretation, coupling moved to the Dirac
spinor form). This subsection is a placeholder naming the hypothesis, not
a derivation. To proceed rigorously rather than invent a plausible-sounding
substitute, the next step needs the actual resonance condition from the
paper — paste or upload the relevant section and this can be aligned
properly, including checking whether it operates directly on
$J_{\text{sync}}(\hat n)$ as proposed here or on a different quantity
entirely.

**Open items:** (1) obtain and align the actual orbital-resonance
discretization mechanism (14.4). (2) Once aligned, check whether it
naturally produces a binary (two-outcome) selection from a continuous
$J_{\text{sync}}(\hat n)$, or whether it needs supplementing with something
race-process-like (Section 12) — these are logically separate questions
even if the same underlying mechanism turns out to answer both. (3) This
branch still does not touch Born-rule statistics (consistent with 13.5,
13.9) — resolving 14.3/14.4 would need to additionally reproduce |α|²,
|β|² weighting, not just discreteness, to fully close the gap.

---

### 14.5 Correction, against the actual Born-selection paper (PAPER_DRAFT.md v0.5)

14.1–14.4 above are kept as written, not deleted, following the paper's
own practice of keeping misses in the record with a dated note rather than
silently rewriting them.

**There is no orbital-resonance mechanism in the paper.** 14.4's placeholder
guessed wrong. The actual discretization mechanism (the paper's clauses 2
and 5) is a **race among exponential, memoryless clocks**, one per
candidate site/mode, with rates proportional to ordinary Born intensity
(|ψ|², clause 1), and exclusivity via "first fire ends every other clock,
a null renormalizes the survivors" (the normalized conditional state,
clause 5). This is explicitly named as the quantum-jump unravelling of
photodetection (Dalibard–Castin–Mølmer / Carmichael / Wiseman–Milburn),
identical to Gillespie's exact stochastic simulation algorithm. The paper
is explicit throughout that this is *located*, not derived: "It does not
derive the square: clause 1 is the rate law. It does not derive
exclusivity: clause 5 is the conditional-state update."

**Validating convergence with Section 12.** This is structurally identical
to the race process built in Section 12 (independent Poisson-type
processes, rates ∝ power, first to fire wins, null renormalizes) — arrived
at there from a linear-response/Fermi's-Golden-Rule argument, entirely
independently of this paper. Two different routes landing on the same
mechanism is a real point in its favor.

**14.4's specific proposal has a demonstrated failure mode, not just an
unproven one.** The paper's §3 terminology section is explicit: *intensity*
means modulus-squared, never the modulus itself, "since the two readings
that differ by that power give Born and its square root." Its
`vacuum_race/` run confirms this concretely: when selection runs on a
quantity linear in amplitude (structurally analogous to a real projection
like $J_{\text{sync}}\cdot\hat n$) rather than on modulus-squared, the
model reproduces $\sqrt r/(\sqrt r+1)$ — the square root of Born — not
Born itself. Section 14's proposal to drive outcome selection from the
projected real current risks exactly this error. This is a demonstrated
failure mode from a real run, not a hypothetical concern.

**A second, independent problem with the premise underlying Sections
1–12.** The paper's §9 (`spad_surface_contact/`) solves the Dirac equation
for a bound electron directly and finds the lower (second) spinor is
*slaved* to the upper one — fixed at 10⁻⁸ degrees, identical across
identical atoms, negligible free content. Conclusion: "the model's phase
race has nothing to select... the phase mechanism leaves no trace in a
passive detector," and clause 7's runs explicitly recommend **dropping**
the second-spinor phase-mismatch idea as a selection mechanism. This bears
on the premise behind Sections 1–12 of this file (spin superposition
modeled via a phase relationship between spinor components as the
selection variable) — the sponsor's own more careful bound-state
calculation already tested a structurally similar idea and moved past it.

**Where the detector's geometric axis actually belongs.** The paper's §2
Stern-Gerlach account resolves this cleanly, and matches the original
intuition ("detector geometry is a factor in up vs down") without the
sqrt(Born) failure mode: the magnet's gradient (the geometric axis) sets up
*which two ordinary complex modes exist* — it spatially separates the
state into two packets, each an ordinary Hilbert-space projection (a
complex amplitude, whose modulus-squared is the ordinary intensity). "The
magnet chooses what the alternatives are... it does not choose which one
happens. The screen does that." The race (clauses 1/2/5) then selects
between *those* ports using ordinary Born intensities. Geometry sets the
*basis*, not an additional real-projected rate term.

**Revised status of Section 14's proposal:** the "project $J_{\text{sync}}$
onto detector axis n̂ and use that for selection" idea, as stated in
14.1–14.4, is withdrawn as a rate-selection mechanism, for the reason
given above. What survives: $J_{\text{sync}}$ and its frame-independence
properties (Sections 13.1–13.3) remain a potentially useful account of
*why correlations don't disagree across frames*, entirely separate from,
and not a substitute for, the ordinary modulus-squared/race mechanism that
actually does the selecting.
