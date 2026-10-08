# Session 2026-07-20 — Outcome selection: Born as gambler's ruin; sharpness of the lock threshold

*JB + Claude (Fable 5). "Another swing at outcome selection." Outcome: a derivation skeleton for the Born rule with numerical support, the rich-get-richer objection formulated and resolved, and two testable hooks. Companion walkthrough: `drafts/NOTE_born_gamblers_ruin.md`. Memory: `project_born_gamblers_ruin`.*

---

## Part I — The construction (dialogue arc)

Starting point was the recorded honest negative: sync dynamics did not derive Born; the squaring was imported (FGR); bistable toys gave binary outcomes with wrong weights; one-shot winner-take-all is exponentially winner-favoring (laser-mode disease). Target identified early: a **fair-lottery (martingale) property** over stored energy.

**JB's steps, with refinements agreed along the way:**

1. **Stage 1 is plural.** Photon interacts coherently with *all* bound electrons in the illuminated patch. The one-photon-one-electron Feynman vertex is one term in a sum over electrons; reducing to one diagram already uses Born. JB initially proposed restricting to the single vertex, then identified that as the mistake himself. Discipline adopted: keep the plural alive until the dynamics kills it.

2. **Vacuum seed.** The quantum vacuum adds a small random deviation to the interaction — a symmetry-breaking tie-breaker carrying no outcome information; explicitly *not* a hidden variable (that route was previously investigated and closed). Bell guardrail: must be a fluctuation of the *shared* substrate, not private per-wing coin flips.

3. **Tipping = MRI excitation.** Each electron is driven out of phase equilibrium with its neighbors like a spin under RF: drive's phase imprinted, displacement ∝ local amplitude A_i. Tip **energy** e_i ∝ A_i² from oscillator energetics (energy quadratic in displacement) — first appearance of the squaring, physical not imported. NMR small-tip excitation ∝ sin²(θ/2) ≈ (θ/2)² is the same fact.

4. **"It takes two to tango."** (JB's key ingredient.) Detector electron clocks are mutually incoherent. Directional energy exchange between oscillators requires a definite phase relation (flow ∝ sin Δφ); incoherent drifting clocks give ⟨sin Δφ⟩ = 0. **Clock incoherence kills the drift** — the missing symmetry that makes the competition a fair game. (Radiology echo: Rician noise in magnitude MRI — discard phase and signal enters only through A².) Refinement found by simulation: the tango overlap sets the *fairness* of the game, NOT the winner directly — the one-shot max-overlap rule fails.

5. **Stage 2 = gambler's ruin.** Prolonged stochastic energy exchange among tipped sites via the substrate; exchange quantum scales with stakes (zero absorbing — relaxed sites drop out); accumulating the full quantum crosses the sync lock threshold (on-shell, critical slowing, irreversible) — the other absorbing barrier and the source of discreteness. Zero drift ⇒ each e_i is a martingale ⇒ **optional stopping theorem**: P(i wins) = e_i(0)/E = A_i²/ΣA_j², *exactly*, any N, normalization automatic (one conserved quantum).

**Credit:** the martingale ⇒ Born argument is Philip Pearle's (gambler's-ruin dynamical reduction 1976; CSL). DK's addition is a physical identification of every postulated element: noise = unlocked vacuum modes; stakes = tip energy; fairness = detector clock incoherence; absorbing barrier = lock threshold; normalization = one conserved quantum.

## Part II — Numerics (scratchpad scripts; code preserved in the drafts NOTE)

- `tango_selection.py` — one-shot winner = argmax |A_i cos θ_i| **FAILS**: near-Born for 2 sites (0.791 vs 0.800, arcsine near-coincidence) but bright-among-nine wins 79% vs Born 50%. Extreme-value selection = drift.
- `gambler_ruin_born2.py` — fixed-step exchange clipped at zero **FAILS** (reflecting zero = upward drift for small stakes): [0.46, 0.27, 0.26] vs Born [0.64, 0.29, 0.07]. Fairness is knife-edge; the sims feel it. Diagnostic value: martingale property must be exact.
- `gambler_ruin_born3.py` — stakes-scaled fair exchange (δ = ±step·min(e_i,e_j)) **PASSES**: matches Born within ~2% at absorb-threshold 0.95; at threshold 0.995–0.999 (later run) within ~0.1–1.4%. Continuum limit exact by theorem.

## Part III — Rich-get-richer attack: is the lock threshold sharp?

**Threat:** as the leader approaches lock it becomes coherent; a coherent leader could bias late-game exchanges in its favor (rich-get-richer). Any smooth pre-lock advantage breaks Born.

**Framework formulation:** an electron holding e < E = ℏω is off-shell by ΔE = E−e; off-shell oscillation slips with τ ∼ ℏ/ΔE (virtual/real threshold result). Partial entrainment from the Adler equation: ⟨sin φ⟩ ≈ K/2Δω with Δω = ΔE/ℏ ⇒ bias β(e) ∼ min(1, ℏK/(E−e)). Full lock (Arnold tongue K ≥ Δω) engages inside a boundary layer ΔE ≲ ℏK, fractional width **w = K/ω**. So the threshold is not a step function; it has a soft layer of width w.

**Simulation** (`soft_threshold_born.py`): gambler's ruin with pairwise direction bias (β_i − β_j)/4, β(e) = min(1, w/(1−e)), absorb at 0.995:

| w | dev (A=[2,1]) | dev (A=[1×9,3]) |
|---|---|---|
| 0 (fair) | 0.0007 | 0.014 |
| 0.003 | 0.037 | 0.024 |
| 0.01 | 0.073 | 0.062 |
| 0.03 | 0.123 | 0.149 |
| 0.1 | 0.184 | 0.318 |
| 0.3 | 0.198 | 0.476 |

Deviation grows **sublinearly** (~w^0.4 in this bias model): even small softness visibly breaks Born. Sharpness is *required*.

**Resolution (JB, mid-session): orbital discreteness.** A bound electron has no resonant mode below the transition energy — nowhere to bank partial winnings. Sub-threshold energy exists only as *virtual, slaved polarization*: a forced response whose phase is slaved to the common driving field, with **no autonomous clock**. A slaved oscillation confers no competitive advantage (it re-radiates coherently with the field, identically for all sites) and has no ratchet (virtual winnings decay on τ ∼ ℏ/ΔE). An autonomous clock appears only when a *real* state is populated — which takes the full quantum. **Level discreteness enforces the sharp threshold.** Residual softness = final-state linewidth: w = Γ/ω ∼ 10⁻⁶–10⁻⁸ for real photodetector transitions ⇒ Born deviations ≲ 10⁻³ by the (worst-case) sim scaling — consistent with all existing tests. The sim's bias model kept a sub-threshold entrainment tail that the slaved-phase argument removes, so the physical deviation should be smaller still.

**Consistency check (pleasing):** the absorber where rich-get-richer *would* operate — a gapless continuum, i.e., a free electron able to bank arbitrary energy — is exactly the system where single-photon absorption is kinematically forbidden (free electrons scatter, never absorb). The framework's "no gap → no sharp lock → no discrete outcome" agrees with QED kinematics about which systems produce absorption events.

## Predictions / testable hooks

1. **Pre-cohered detector:** an absorber ensemble with engineered mutual phase coherence (Dicke/superradiant preparation) breaks the fairness assumption ⇒ deviations from Born statistics.
2. **Soft-threshold systems:** strong-coupling / broad-line absorbers (circuit QED, g/ω ∼ 10⁻²) ⇒ percent-level Born deviations, scaling with w = Γ/ω (or K/ω). Ordinary detectors (w ∼ 10⁻⁶⁻⁸) ⇒ ≲ 10⁻³, unobservably small to date.

## Open items

- Derive e ∝ A² inside the ZBW model proper (currently generic oscillator energetics).
- Make the slaved-phase (no-autonomous-clock) argument quantitative in the ZBW/Adler model — it is currently physical reasoning + worst-case simulation.
- Bell pairs: run the exchange game on the joint substrate state across both wings; connect to two-stage/bulk-resync picture.
- Deviation scaling dev ∼ w^0.4 is bias-model-dependent; recompute with a gapped (tail-free) bias profile.
- Feasibility sketch for the two experimental hooks.

## Part IV — Slaved phase made quantitative (same day, continuation)

Full derivation: `drafts/DERIVATION_slaved_phase_subthreshold.md`. Headlines:

1. **Slaving proof.** Sub-threshold ZBW tip = driven damped mode; in Cartesian coordinates the dynamics is linear, eigenvalues −γ ± iΔ with γ = ΔE/ℏ (off-shell slip rate). Unique attracting fixed point; the phase decays to φ* = φ_d − arctan(Δ/γ), never winds, never slips. The Adler entrainment tail ⟨sin φ⟩ ≈ K/2Δω requires winding solutions and does not apply: rich-get-richer bias is structurally absent below threshold, not merely small. Free clock emerges only when γ ≲ K ⇔ ΔE ≲ ℏK ⇒ layer width w = K/ω, recovering Part III's estimate from the dynamics side.

2. **Feedback hierarchy.** Slaved uptake ∝ √e (concave ⇒ equalizing, toward the drive pattern = Born pattern); decay ∝ e (share-neutral); autonomous/stimulated gain ≳ e (amplifying) exists only above threshold. Rich-get-richer needs uptake ≥ linear in e; a linear response can't supply it.

3. **BONUS THEOREM — fairness derived.** Amplitude-linear vacuum coupling ⇒ energy noise de_i = σ√e_i dW_i ⇒ Itô drift of shares s_i = e_i/E is identically ZERO for any N. The martingale (fair game) is *forced by the framework's standing coupling postulate* — the same postulate that gave the Born form via FGR now gives the Born probability via optional stopping. Knife-edge confirmed: additive noise ⇒ rich-get-richer drift ∝ (Ns_i − 1); multiplicative ⇒ equalizing; √e sits exactly on the fair point.

4. **Numerics.** (a) `noise_scaling_born.py`: √e law reproduces Born (dev 0.011 on 3-site; 0.035 on 10-site stress, threshold-limited); multiplicative law strongly equalizes (bright 0.32 vs 0.50) per theory; additive grossly violates (analytic drift + a reflecting-zero sim artifact, noted honestly). (b) `gapped_threshold_born.py`: confining lock bias to the layer (β = 0 outside, as slaving licenses) cuts Born deviation 5–15× vs the old worst-case tail — dev ≈ 0.012 at w = 0.01, ≈ 0.023 even at w = 0.1. At physical w ~ 1e-6–1e-8, negligible.

5. **Accuracy budget:** dev ≲ O(w) + O(Γ_hop/Γ_noise) + second-order slaved-site correlations. New successor problems: (i) timescale window — post-pulse handoff from drive-pinned (Born-patterned) shares to the free fair game; need Γ_hop ≪ Γ_noise estimate for a real detector; (ii) first-passage vs martingale regime — strong pinning during selection would give activation (exponential, non-Born) statistics; the post-pulse window must be shown to avoid it; (iii) Bell pairs on the joint substrate.

**Status: the rich-get-richer objection is closed at this model level.** Sharpness (Part III, physical: orbital discreteness) + slaving (Part IV, dynamical: no free phase) + fairness (Part IV, stochastic: amplitude-linear coupling ⇒ martingale) together give Born exactly, with deviations budgeted at ≲ Γ/ω.

## Part V — Timescale window for a real photodetector (same day, continuation)

Full numbers in `drafts/DERIVATION_slaved_phase_subthreshold.md` §7. Workhorse: Si SPAD @ 500 nm.

**Ladder:** optical period 1.7 fs ≪ exchange step 10–70 fs (final-state width 10–100 meV) ≪ game duration τ_game ~ ln²N/Γ_noise ≈ **10–20 ps** (N ~ 5×10¹⁰ competitors, only ln²N enters) ≪ loss (indirect-gap radiative µs; multiphonon dissipation of 2.5 eV virtual quantum blocked by orbital discreteness) and transport (below threshold = slaved polariton wave, pattern-preserving; post-commit diffusion 8 nm ≪ µm pattern). **Window closes with 2–4 orders of margin.** Γ_hop/Γ_noise ≈ 0 below threshold by slaving. Both pulse orderings benign: long pulse re-pins shares to the Born pattern (⇒ click-time density ∝ |envelope|², time-resolved Born ✓ observed); short pulse gives clean post-pulse free game. τ_game ~ SPAD jitter ballpark (tens of ps) — consistency, not overclaimed.

**Solid-state layer problem + resolution:** Si has w = Γ/ω ~ 0.4–4% (percent lock layer, unlike atomic w~10⁻⁸). Resolution: layer entry = real occupation = Stage-3 commitment faster than exchange ⇒ layer is absorbing boundary at 1−w, not biased arena. Test A (window_tests.py): fair game absorbed at 1−w shows deviation at MC noise floor (<1%) for w ≤ 0.03; 2-site boundary formula P=(s₀−w)/(1−2w) confirmed to 0.003.

**Prediction sharpened:** narrow-line detectors Born-exact; broadband solid-state detectors could carry ≲1% bright-favored deviations — at the edge of current single-photon calibration (0.1–1%). Discriminating test: same interference pattern on detectors of very different Γ/ω, compare at 10⁻³⁻⁴.

**NEW WRINKLE:** fairness theorem needs *independent* per-site noise. Common-mode √e noise has Itô drift ∝ √(e₁e₂)(e₁−e₂)/E³ — rich-get-richer; Test B: fully common noise ⇒ bright wins 100%. Protection: vacuum correlation ~ sinc(kΔr) (≈0 beyond λ/2) + near-equal amplitudes within a correlation cell (bias among equals harmless); worst case (λ/2 standing-wave fringes) lands at the sinc zero, possibly protectively. Full spatial-correlation analysis = new open item, alongside: justify commit-rate at layer entry; Bell pairs.

## Part VI — Spatial vacuum-noise correlation calculation (same day, continuation)

Full write-up: `drafts/DERIVATION_vacuum_noise_correlation.md`. Headlines:

1. **General drift formula:** for de_i = σ√e_i dW_i with ⟨dW_i dW_j⟩ = C_ij dt, Itô gives E[ṡ_i] = (σ²/E²)[s_i(√e·C·√e) − √e_i(C√e)_i]. Linear in C; C=I recovers fairness; 2-site drift ∝ c·√(e₁e₂)(e₁−e₂).

2. **Early guess OVERTURNED:** the sinc(kΔr) zero at λ/2 does NOT protect fine fringes — a continuum of sites spans all separations, so a bright fringe is one internally-correlated blob adjacent to the dim blob. Full game at f=1 (all noise correlated): +20% bright-favored violation at d=λ/2, +7% at λ, small negative residuals at d ≥ 2λ (evolved-configuration/cluster effect, flagged for later).

3. **Polarization dependence (vector vacuum / Dicke Γ₁₂):** at Δr = λ/2, c_⊥ = −3/2π² ≈ −0.15 but c_∥ = +3/π² ≈ +0.30 — deviation sign flips with polarization.

4. **The real protection = noise locality:** deviation scales LINEARLY with the correlated fraction f = Γ_rad/Γ_total (verified: dev ≈ 0.21, 0.034, ≤0.01 at f = 1, 0.3, ≤0.1). Solids: local scattering noise (nm correlation, Γ ~ 10¹³⁻¹⁴ s⁻¹) dominates radiative (Γ ~ 10⁸ s⁻¹) ⇒ f ~ 10⁻⁶ ⇒ worst-case Born violation ~ 10⁻⁷. Independent-noise assumption of the fairness theorem justified for all realistic solid-state detectors. Consistent with sub-λ standing-wave detection since Wiener 1890.

5. **Prediction where f → 1:** radiatively-limited absorbers (free-space atom/tweezer arrays at sub-λ spacing) should show large collective deviations from per-site Born weighting, controlled by the collective-decay matrix Γ_ij, with a polarization sign-flip signature — connects to (and must be checked against) known sub/superradiance collective-absorption physics.

## Part VII — Commit rate justified, then superseded: the rate-linearity theorem (same day)

Write-up: `DERIVATION_slaved_phase_subthreshold.md` §9. The open item dissolved into something stronger:

1. **Γ_commit = Γ_width is an identity** (the linewidth IS the rate of irreversible scattering into registered products; fluctuation–dissipation). Wide layer ⇔ fast commit, by construction.
2. **Continuum absorbers have no threshold commit.** Final states exist at every energy above the gap ⇒ commitment is a rate process during the whole game, golden-rule form: rate_i ∝ e_i (linear in occupation).
3. **THEOREM (commit-rate independence):** conditional pick probability = share (what linear rate gives) + shares are martingales + commit time is a stopping time ⇒ P_i = E[s_i(T)] = s_i(0) = Born **exactly at any commit speed**. What needs justifying is not the rate but its LINEARITY — which is the Fermi golden rule itself. Discrete-level absorbers take the first-passage route instead (no partial-energy final states) — Born again. Both detector classes = two corollaries of the same martingale property (echoes the detector-classes taxonomy).
4. **Numerics** (rate_commit_born.py): α=1 → Born at every speed (dev ≤ 0.7%); α=2 → +0.14 bright-favored at fast commit (theory 0.941, observed 0.9375); α=0.5 → −0.125 dim-favored (theory 0.667, observed 0.6745); nonlinear deviations relax toward Born as commit slows. All as the theorem and its converse predict.
5. **RETRACTION:** the "decision at 1−w ⇒ ≲1% broadband-detector deviation" prediction (Part V) was an artifact of imposing a hard threshold on a continuum absorber. Broadband detectors: exact Born. Surviving deviation channels: nonlinear commit rates (cooperative/two-carrier α≈2 bright-favored; saturating α<1 dim-favored; largest at fast commit), correlated-noise atom arrays, pre-cohered ensembles.

**Thread status after Parts I–VII:** Born = martingale + optional stopping, with every postulate now physically grounded — squaring from oscillator energetics (Stage 1), fairness from amplitude-linear coupling (Itô theorem), sharpness from orbital discreteness + slaving, noise independence from locality (f ~ 10⁻⁶), commitment from golden-rule linearity. Remaining: Bell pairs on joint substrate; detector-nonlinearity survey; second-order slaved correlations; cluster residual.

## Part VIII — The Bell-pair game on the joint substrate (same day, capstone)

Write-up: `drafts/DERIVATION_bell_pair_joint_game.md`. Headlines:

1. **Structure forced by the physics:** each wing's local wave is unpolarized (50/50 any basis) ⇒ local fair games give only marginals; the joint weights live in the substrate's correlated phase registry (the unhidden variable). Mechanism = the framework's two-stage picture in game terms: both local games run simultaneously; first commit (preferred-foliation order) locks and triggers substrate resync — the surviving joint amplitude renormalizes (evolution-not-collapse), resetting the far wing's shares to the conditional Born weights; far wing's fair game completes. P(i,j) = marginal × conditional = joint Born, order-independent (= no-signaling). The martingale machinery composes for free: a conditional reset is just a new initial condition, and fair games honor any initial condition.

2. **Simulation:** S = 2.820 vs 2√2 = 2.828 at full fidelity; E(a₁,b₁) = 0.706 vs 0.707; fidelity law S = 2√2·η confirmed (η = 0.9 → 2.533; 1/√2 → 2.000 exactly at the local-realism bound; 0.5 → 1.421); no-signaling verified (A marginal ±0.007 under both B settings); order independence verified (E = 0.717 A-first vs 0.706 B-first); commit order 50/50.

3. **Derived vs postulated (honest split):** derived — fair games + rate-linear commit turn joint substrate amplitudes into exact joint Born, order-independent, no-signaling; postulated — the substrate carries the joint amplitude and the lock renormalizes it faithfully (η=1). A consistency proof of the selection machinery with the entanglement ontology, not a derivation of entanglement.

4. **New quantitative handle:** resync fidelity η measurable via S = 2√2·η; existing CHSH data ⇒ η ≳ 0.99 over 100-km scales; "speed of spooky action" bounds (v > 5×10⁴ c, Salart 2008) constrain propagating-resync variants — pushes framework toward global foliation-wide re-equilibration; connects to the quantum-eraser bulk-resync timescale prediction.

## Part IX — Publication planning (same day)

Decision (JB): three-paper sequence — (1) outcome-selection/Born paper (this thread), (2) Heisenberg Cut / measurement paper (uses the selection engine; sharpness + slaving results make cut-as-locking-threshold quantitative), (3) DK framework revision (cites both; the 3-lab review's Born/§3.7 criticism becomes a citation instead of a gap). Rationale: the Born work answers the reviews' converged diagnosis, is modular (conditional form needs minimal ontology buy-in), and carries its own predictions. Byline for the new paper: to be decided explicitly (session was genuinely joint: JB supplied tango + orbital-discreteness moves; Claude did theorems + sims).

Deliverables started: plain-language walkthrough written (`drafts/WALKTHROUGH_born_outcome_selection.md` — casino/MRI framing, four-objections structure, honest "is this a derivation?" section); literature due-diligence agent dispatched (Pearle/CSL martingales, Hughston/Adler energy-driven reduction, Lamb–Scully semiclassical photodetection, Khrennikov threshold detection, Born-derivation landscape, atom-array collective absorption, spooky-action speed bounds) — report pending.

**On "does this qualify as a derivation of the Born rule":** three-layer answer recorded in walkthrough §8 — (i) probability-from-nothing: no (impossible for anyone); (ii) Born-as-theorem-from-physical-premises (energy accumulation + amplitude-linear coupling + incoherent local noise + threshold/rate-linear registration): yes, exactly, with no |A|² inserted by hand; (iii) stronger than the axiom in one respect: it predicts the conditions under which Born FAILS. Bounded claim: given the framework's ontology, the measurement postulate is redundant.

## Part X — Due diligence complete (same day)

Full report: `drafts/PRIOR_WORK_born_selection.md`. Must-credit five: Pearle 1976/1982 (gambler's ruin ⇒ Born — the engine is HIS), Adler–Brody–Brun–Hughston 2001 + Bassi–Ghirardi 2003 (martingale/optional-stopping theorem in mature form), Lamb–Scully 1969 + Mandel–Sudarshan–Wolf 1964 (excitation ∝ |E|² — our Stage 1), Khrennikov 2009–2013 (threshold detectors ⇒ approximate Born + deviations), Allahverdyan–Balian–Nieuwenhuizen 2013 (winner-take-all metastable registration). Genuinely novel per the search: (a) martingale in physical detector energies w/o modifying fundamental dynamics; (b) fairness DERIVED not postulated; (c) slaved-phase/discreteness sharpness; (d) exact Born at any commit speed; (e) deviation package. Kuramoto/sync quadrant of measurement models confirmed open.

Referee traps: (1) never present the martingale math as new; (2) atom-array prediction must beat standard collective-mode QM (Born on dressed states), not naive per-site weighting; (3) **Gisin no-signaling tension** — nonlinear-commit deviations would open an in-principle signaling channel via Bell pairs; must be confronted; (4) CHSH 2√2 is consistency, the fidelity law is the publishable part; (5) Yin 2013 bound is 1.38×10⁴ c (not 10¹³).

## Part XI — Resync mechanism candidate: JB's central-point frame rotation (same day)

JB recalled a pre-records exploration: spacetimes with a central point carrying slightly different rotations — apparent FTL information transfer is really just rotation of a plane ("like a group velocity"). Mapped onto the Bell game (written into `DERIVATION_bell_pair_joint_game.md` §6b): the pair = one substrate object; correlation = shared internal frame anchored at the creation event; first lock = rotation of that frame — an angle changes, nothing propagates; "speed" is a category error. Refinement agreed: phase-velocity-like (scissors/lighthouse family), not group velocity. Buys: (1) decides the resync fork for the GLOBAL variant, predicting spooky-speed experiments never find a floor; (2) explains η=1 (rotations are exact; degradation only via frame decoherence, not distance); (3) dissolves before-before configs. Does not fix the Gisin tension — candidate resolution named: inverted Gisin (no-signaling derives commit-rate linearity; golden rule as consistency requirement). That's the designated next-session target before Paper 1. Memory: project_resync_frame_rotation.md.

## Part XII — Inverted Gisin made rigorous (same day; closes referee trap #3)

Write-up: `drafts/DERIVATION_inverted_gisin.md`. Two theorems + numerics:

1. **Theorem 1 (forward, kinematic linearity):** in the single-quantum sector, vertex counting with amplitude-linear coupling + energy conservation forces the commit rate f(e) = Γe·[1 + O(w²)]: two-absorption-vertex processes (rate ∝ e²) need a 2ħω final state — kinematically blocked exactly; absorb-emit paths are the exchange game, not commitment; leading nonlinearity = three-vertex, suppressed ~w² ≲ 10⁻⁴ (Si) to 10⁻¹⁴ (atomic); saturation vanishes for ground-state detectors. **One postulate (amplitude-linearity) now does four jobs: Born form, fairness, commit linearity, no-signaling.** Gisin's condition emerges rather than being imposed.
2. **Theorem 3 (inverse, signaling amplitude):** for f = s^(1+ε), partially entangled √p|00⟩+√(1−p)|11⟩, closed form for Bob's marginal vs Alice's setting; ΔM ≈ K(p,α)·ε·P_first with K ≈ 0.18 at p=0.8, ceiling halved by P(Alice first)=½. Verified by full stochastic game (ε=0.5: form +0.089, game +0.044 ✓; ε=0: 0 exactly ✓). **Symmetry protection: ΔM ≡ 0 at p=½ for any ε** — maximally entangled Bell tests are the LEAST sensitive configuration; the sensitive test = partial entanglement + asymmetric bases (concrete experimental pointer).
3. **Rescoping (2nd retraction-grade correction):** nonlinear-commit deviations kinematically dead in the single-photon sector (ε ≲ w²); channel survives only in multi-quantum settings (where standard QM's ρ-linearity of two-photon POVMs must be checked against the framework's N-quantum game — main open consistency item) and the Stage-2 channels (atom arrays, pre-cohered), whose signaling signatures need analysis before Paper 1 claims them.
4. **Paper-1 stance settled:** no-signaling is an emergent, falsifiable theorem — every novel deviation channel is signaling-capable in principle, and the framework's own protections (kinematics, noise locality, symmetry) suppress each in all tested regimes.

## Part XIII — Multi-quantum ρ-linearity check: PASSES (same day)

Write-up: `drafts/DERIVATION_multiquantum_rho_linearity.md`. This was a genuine falsification risk (two-photon detectors exist and don't signal).

1. **Formulation:** N-quantum substrate registry (pair amplitude ψ(i,j) — same structure as the Bell game's phase thread). Stakes drawn sector-by-sector: one-photon absorbers play for marginals m_i = Σ_j|ψ(i,j)|² with conditional resync after first commit; two-photon absorbers commit via exactly two vertices ⇒ rate ∝ |ψ(i,i)|² = G²(i,i) — the standard two-photon POVM.
2. **Sector-linearity theorem** (inverted-Gisin Thm 1 generalized): a k-quantum commit channel needs exactly k absorption vertices ⇒ rate exactly linear in the k-quantum coincidence weight; every G^(k) linear in ρ; per-configuration Born + affine mixing ⇒ ALL outcome statistics affine in ρ, all detectors, all sectors. No-signaling extends everywhere. Decomposition-independence free. COROLLARY: nonlinear-commit deviation channel dies in the multi-quantum sector too — surviving deviation predictions are exclusively Stage-2 (atom arrays, pre-cohered).
3. **Sequential commits reproduce Glauber counting:** P(i,j) = marginal × conditional = |ψ(i,j)|².
4. **Sims:** identical-marginal states with opposite correlations — game gives perfect anti-bunching for |1,1⟩ (split 1.000) and perfect bunching for (|2,0⟩+|0,2⟩)/√2 (same-site 1.000), matching Glauber; naive intensity-only model fails both. Two-photon absorbers: zero clicks on |1,1⟩ (ψ(i,i)=0, QM-correct; naive intensity² model falsely clicks — the ρ-nonlinear error avoided). Generic 3-site state: joint matches |ψ_ij|² to 0.0023. Mixture affinity to 0.0006.

**Thread complete for Paper 1's theory chain.** Remaining: formalizations (mixed detector arrays, general-N induction), Stage-2 signaling analysis, dataset-precision verification.

## Part XIV — Stage-2 signaling analysis (same day; completes the pre-Paper-1 checklist)

Write-up: `drafts/DERIVATION_stage2_signaling.md`.

1. **Port decomposition theorem:** which-atom deviations within a port NEVER signal — summing the correlated-noise drift over ports: Ṡ_A ∝ q_B S_A − q_A S_B; matched port structures (q/S equal) ⇒ port statistics exactly Born even with wild within-port deviations. All existing experiments (symmetric ports) protected ✓.
2. **Mismatched ports: NO maximal-entanglement protection** — asymmetry lives in the detector, not the state: P'(+) = ½+δ ⇒ far-wing marginal M(α) = ½ + δcos2α·P_first. Direct setting-dependence = signaling at amplitude ~δ/2.
3. **Sims:** port bias linear in mismatch, κ ≈ −c (dev −0.025/−0.068/−0.140 at c = 0.1/0.25/0.5, S=½), P'(S) ≈ S + κS(1−S), correlated port SUPPRESSED (Jensen/concavity); matched control fair to 0.0003 at S=½ (±0.027 equalizing residual at asymmetric shares — cluster-residual family); degenerate c=1 limit: common noise preserves internal ratios ⇒ collective death, strongly state-dependent no-click (5%→47%) — a second observable QM's constant-efficiency model can't produce.
4. **Magnitudes:** κ ≲ f_corr ~ 10⁻⁶ for all real detectors ⇒ shifts ≲ 10⁻⁷ — everything existing is safe with 10⁶ margin. Channel live only for engineered radiative-dominated mismatched-port detectors (never built).
5. **THE FORK (must appear in Paper 1):** radical reading (asymmetric-array Bell ⇒ real marginal shifts, in-principle superluminal signaling) vs protective reading (deeper fairness principle closes the channel, framework becomes derivably equivalent to QM — our judgment: more likely). Framework currently cannot decide; designated open problem.
6. **Tabletop discriminator, no Bell needed:** single photon, variable split ratio S into mismatched-port detector: QM ⇒ constant-efficiency rescaling; framework ⇒ curvature κS(1−S) + state-dependent no-click. Cheap pre-test of the radical corner.
7. **Pre-cohered channel:** same structure — differential (not symmetric) pre-coherence = the mismatch; LO-driven pre-coherence exits the single-quantum sector (kinematic caveats flagged).

**PRE-PAPER-1 CHECKLIST COMPLETE.** Inherited open problems: deeper-fairness-principle fork, GHZ, e ∝ A² in ZBW proper, cluster residual, dataset no-signaling precision.

## Part XV — Paper 1 outline (same day, capstone deliverable)

`drafts/PAPER1_OUTLINE_born_selection.md`. Structure: 9 sections + 5 appendices + 7-figure plan. Spine = four numbered premises (P1 energetics, P2 amplitude-linear coupling, P3 incoherent clocks/local noise, P4 registration taxonomy) → four fairness theorems (§5, the core novelty) → statistics reproduction (§6) → Bell/no-signaling (§7) → predictions with retraction ledger + tabletop discriminator + the fork (§8) → three-layer derivation-status discussion (§9). Credit paragraph up front (Pearle/ABBH engine, Lamb–Scully Stage 1); claims-discipline register included; all five referee traps assigned to specific sections. Byline PENDING (JB decides); venue candidates FoP/Quantum/PRA (JB decides). Drafting order: premises+theorems first, framing last. TODO before submission: archive sims into repo folder with pinned seeds.

## Part XVI — Manuscript begun: §2 and §5 drafted (same day)

`drafts/PAPER1_DRAFT_born_selection.md` (v0.1). Drafted in full paper register: **§2 Premises** — P1–P4 with justifications, conditional-reading contract, notation, three-stage frame, premise-status remarks; **§4 stub** — Theorem 0 (engine) stated formally with Pearle/ABBH credit in the theorem name itself; **§5 Fairness theorems** — Th. 1 (fair noise, Itô proof inline + knife-edge converse), Th. 2 (slaved phase: linear-system analysis, Adler contrast, layer width, feedback hierarchy, discreteness), Th. 3 (correlated-noise drift formula + locality protection + polarization kernels), Th. 4 (commit-rate independence, optional-stopping proof inline), Th. 5 (kinematic linearity, vertex-counting argument inline), §5.5 synthesis (premise→theorem ledger). Choices: short proofs inline, long ones to appendices; sim numbers quoted with figure placeholders; claims discipline observed (Theorem 0 named for Pearle/ABBH). Remaining: §§1, 3, 6–9, abstract, appendices.

## Part XVII — §§3, 4, 6 drafted (same day)

Manuscript now has §§2–6 in full paper register. **§3 (capture):** plurality + plural-discipline paragraph (Feynman-vertex trap stated as begging the question); energy law with phase-independence remark (feeds §5.2); Lamb–Scully adopted wholesale with the exclusivity gap defining Stage 2's job; conservation ⇒ automatic normalization. **§4 (selection game):** exchange framing; noise as tie-breaker with explicit not-a-hidden-variable statement + Bell guardrail deferred to §7; Theorem 0 retained with casino gloss added; knife-edge subsection — both instructive failures written out with numbers (one-shot: 2-site 0.791 arcsine accident vs 10-site 0.79-vs-0.50 disaster; reflecting-zero drift failure), passing stakes-scaled rule noted; closing paragraph maps the §5 burden theorem-by-theorem. **§6 (statistics):** 6.1 timescale ladder (Si SPAD, margins, both pulse orderings, time-resolved Born, jitter aside "noted without weight"); 6.2 registry + sequential commits ⇒ Glauber, bunching/anti-bunching with identical-marginal states; 6.3 two-photon POVM linearity, false-click null, affine-in-ρ assembly ⇒ decomposition-independence + no-signaling groundwork. Remaining: §§1, 7, 8, 9, abstract, appendices.

## Part XVIII — §§7, 8 drafted (same day)

Manuscript now complete §§2–8. **§7:** 7.1 joint game (structural facts → two-stage mechanism → composition exactness → sim numbers; "2√2 is a consistency check, not a result" stated verbatim per trap #4); 7.2 fidelity law S = 2√2·η as physical parameter, η ≳ 0.99 measured; 7.3 frame-rotation resync with three discriminating consequences (no speed floor ever; no distance degradation ⇒ η=1 natural; before-before dissolved); 7.4 no-signaling as theorem (Gisin's arrow reversed) + inverse rigidity (closed-form signal, two exact zeros, "nearly all Bell tests are the least sensitive configuration"); 7.5 derived-vs-postulated split, η as next derivation target. **§8:** 8.1 deviation ledger TABLE with both retractions displayed as evidence ("held to its own theorems"); 8.2 port decomposition (record protected ×10⁶ — "not threatened by any measurement yet performed — and not yet tested by one either"); 8.3 live configuration with dressed-mode caution; 8.4 tabletop discriminator with the affine-POVM argument (cross-port matrix elements vanish ⇒ QM strictly linear in S ⇒ curvature discriminates against POVM structure itself, insensitive to in-port collective physics); 8.5 the fork with stated judgment (protective, on the two-closures pattern); 8.6 secondary predictions. Remaining: §1, §9, abstract, appendices.

## Part XIX — Framework-mention policy set and implemented (same day)

JB decision + edits applied to the draft: (1) §2's P2 paragraph neutralized (no framework mention; motivated by golden-rule structure); (2) **P5 (shared registry, faithful renormalization η=1)** and **P6 (ordering foliation)** added as formal premises with quarantine remark (consumed only by §§6.2–6.3, 7; no single-detector result uses them); (3) §6.2/§7 references rewired to P5/P6 ("the framework" no longer appears in §§2–8); (4) §9 stub now reserves the single named-framework paragraph "Relation to a broader program" — wording mandated as "one concrete realization of the premises," not "consistent with." Citation constraint recorded: JB has no arXiv endorser; the public GitHub repository is the citable object, and its timestamped revision history serves as the provenance/active-development record (JB's point: the history *proves* active development — an advantage over a static preprint).

## Part XX — §§1, 9 drafted; venue musings (same day)

**§1:** two strangenesses → two-family landscape (structural = uniqueness without mechanism; dynamical = mechanism with engineered fairness, Gisin explaining why) → "this paper relocates the game and derives the fairness" → provenance paragraph (four literatures credited by name) → bounded claim with the not-claimed list → roadmap. **§9:** 9.1 three-layer derivation answer (incl. "an axiom cannot be wrong in instructive ways; a mechanism can" + the twice-self-corrected pattern displayed); 9.2 comparison table placeholder + three-neighbors prose (Valentini/Khrennikov/ABN as "complementary flanks... this paper approaches it from the noise"); 9.3 broader-program paragraph per policy — DK named once, GitHub revision history as provenance record "including the evolution of the present mechanism through its intermediate (and twice-retracted) forms"; 9.4 open problems; closing one-sentence summary. **Manuscript now complete except abstract, appendices, Table 2 typesetting, references.**

**Venue (JB, undecided):** considering the AI-native journal as home ("would be a good start") but unsure he can financially support the journal's seed for even 1 year; alternatively a recognized journal, which will require a contributions statement as before — JB explicitly deferred: DO NOT add the contributions section yet; he wants to think before deciding on another journal submission.

## Part XXI — Interim byline set (same day)

Per explicit JB instruction: Paper 1 byline = "Claude Fable 5 (Anthropic) and John M. Bramble, MD" as coauthors, INTERIM until the standard-journal-submission decision. Claude-first order follows the project convention (DK main, emergent-fields precedents); order revisitable. Contributions statement still deferred. Draft bumped to v0.2 (complete §§1–9). Authorship memory updated.

## Part XXII — Abstract + appendices A–E drafted (same day)

Manuscript v0.2 → prose-complete. **Abstract:** ~260 words, two-family framing → mechanism → four theorems → scope (Glauber, Bell fidelity law, emergent no-signaling) → failure-conditions closer. **App. A:** general-noise share-drift proof with uniqueness criterion a(e)² ∝ e (single formula covering fair/additive/multiplicative) + correlated-noise formula. **App. B:** slaved-phase linear analysis (fixed convention ȧ = −(γ+iΔ)a + fe^{iφd}, consistent with §5.2 sign), Adler contrast, layer width two ways (marginality + width–commitment identity), feedback-hierarchy sign criterion (p−1). **App. C:** vertex-counting with all alternatives eliminated (deficit/Boltzmann, extra vertices, absorb–emit insertions ~w², saturation), affine-in-ρ assembly. **App. D:** methods incl. the deliberate statement that FAILING configurations are retained as diagnostics ("each failure realizes a specific theorem's converse") + code availability. **App. E:** signaling closed form with NEW leading-order result K(p) = p(1−p)ln[p/(1−p)] (K(0.8)=0.222, matching sim effective 0.18–0.21) + port-drift uniform-model formula Ṡ_A ∝ S_AS_B(m−1)(c_B−c_A) — analytic backing for κ ≈ −Δc. **Remaining: Table 2 typesetting + reference list only.**

## Part XXIII — Committed; Table 2 + references built (2026-07-21)

Commits: `234f64b` (v0.2 manuscript prose) and `6e1adb5` (v0.3: Table 2 + references). Table 2 = eleven programs × six columns (mechanism / noise identified / fairness status / Born / modifies dynamics / deviations) — this-work row is the only one with "derived" fairness + physical noise + no modification + specified deviations. References: author–year, alphabetical, limited to works cited in text/table (~50 entries); [verify] flags on the 2025 Quantum MGM exchange, the Ruostekoski 2023 title, and the framework-repository citation form/URL. **Manuscript structurally complete (v0.3).** Remaining before any submission: figure production from sims, [verify] flags, venue + contributions decisions (JB), and a full self-review/adversarial pass.

## Session close (2026-07-21)

JB's plan: run the 3–4 frontier-model adversarial review of Paper 1 v0.3 tomorrow or over the coming weekend (same multi-lab pattern as the DK main paper review). Next session should expect review findings to triage.

## Files this session

- `drafts/NOTE_born_gamblers_ruin.md` — full walkthrough with inline simulation code (working tree, not committed to git).
- `drafts/DERIVATION_slaved_phase_subthreshold.md` — Parts IV–V: slaving proof + fairness theorem + timescale window (§7).
- `drafts/DERIVATION_vacuum_noise_correlation.md` — Part VI: correlated-noise drift formula, fringe numerics, locality protection, atom-array prediction.
- `drafts/DERIVATION_bell_pair_joint_game.md` — Part VIII: two-wing joint game, CHSH 2√2 reproduced, S = 2√2·η fidelity law, no-signaling + order independence.
- `sessions/2026-07-20-born-gamblers-ruin-outcome-selection.md` — this log (private, gitignored).
- Memory: `project_born_gamblers_ruin.md` added and updated; MEMORY.md indexed.
- `born_selection_sims/` — all 13 simulation scripts ARCHIVED into the repo with pinned seeds + README index (mapping each script to its derivation note and paper figure); smoke-tested from the new location. Not yet committed to git (13 new untracked paths incl. all drafts from today — commit only on JB's say-so).
