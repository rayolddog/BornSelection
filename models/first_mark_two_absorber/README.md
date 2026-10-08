# first_mark_two_absorber — the 576-state two-absorber model (no-bath slice)

This directory holds a small, self-contained quantum simulation. It is isolated
from every manuscript simulation in this repository and nothing outside this
directory reads from it or writes to it.

## What this is, in plain English

Imagine one photon arriving at a beamsplitter and travelling down two paths at
once. At the end of each path sits an absorber that can catch it. Each absorber
also has two small "bath" companions: one that keeps a record of the local
background phase, and one that would carry away the energy if the absorber ever
made a permanent mark.

Write down every combination of those things and you get 576 possibilities. This
code builds all 576 of them explicitly, writes down a single energy operator (a
Hamiltonian) that says how amplitude flows between them, and evolves the state
forward in time by the ordinary Schrödinger equation.

That is the whole model. There is **no collapse, no random outcome, no jump, and
no measurement** anywhere in it. It is the plain unitary baseline that any later
claim about how one outcome gets selected would have to be compared against.

This first slice runs the model in its **no-bath setting**: the phase-history and
mark/record couplings are set to zero, so the only thing happening is a photon
being coherently caught by the absorbers and coherently released again. That is
the calibration step, and it is the step the numbers below actually verify.

## The 576 states

The subsystem order is fixed once in `model.py` and is used everywhere:

```
field, absorber_A, absorber_B, phase_A, record_A, phase_B, record_B
  4         3           3         2         2         2         2     = 576
```

| Subsystem | Dimension | Basis order |
| --- | --- | --- |
| Field | 4 | vacuum, photon in arm A, photon in arm B, collected output |
| Absorber A | 3 | ground, reversible excitation, candidate mark |
| Absorber B | 3 | ground, reversible excitation, candidate mark |
| Phase A | 2 | phase-tag 0, phase-tag 1 |
| Record A | 2 | blank record, excited record |
| Phase B | 2 | phase-tag 0, phase-tag 1 |
| Record B | 2 | blank record, excited record |

Two energy relations are enforced by construction rather than fitted:

- the photon energy equals the reversible excitation energy, `E_e = ħω`;
- forming a mark pays the record qubit exactly the difference,
  `E_e = E_m + ħΩ_r`.

## The Hamiltonian

`H(t)` is the sum of seven named pieces, each with one job:

| Term | Job |
| --- | --- |
| `H_field` | The photon carries the same energy in arm A, arm B, or the output mode |
| `H_absorber_i` | Absorber energies, plus a controllable detuning `Δ_i(t)` |
| `H_bath` | Free energies of the phase and record qubits |
| `H_capture` | Photon in arm `i` ↔ reversible excitation of absorber `i`, strength `K_i(t)` |
| `H_phase` | The phase qubit shifts the phase of the excitation without moving population, strength `χ_i·u_i(t)` |
| `H_mark` | Reversible excitation + blank record ↔ candidate mark + excited record, strength `G_i` |
| `H_route` | A controllable recombiner sending arm amplitude into the shared output mode, strength `J_i(t)`, phase `θ_i` |

The **no-bath setting is `χ_A = χ_B = 0` and `G_A = G_B = 0`**, which switches off
the last two rows' bath involvement and leaves capture, detuning and routing.

`χ` and `G` live in `Params` (physical constants), while `K`, `Δ`, `u`, `J` and
`θ` live in `Controls` (things an experimenter sets). That split is deliberate:
the plan states that operational reversal is *not* permitted to flip `G_i`, so
`G_i` is not exposed as a control knob at all. Only the exact global inverse
reverses it.

That split has a consequence for the return helper. Reversing the sign of `K`
inverts the capture term and nothing else, so a detuning or a bath coupling left
running through the pulse is never undone and the photon does not come back.
`capture_and_return_schedule` is therefore a **resonant, no-bath helper, and it
enforces that** — asking it for a detuned or bath-coupled return raises
`ResonantReturnError` rather than handing back a schedule that would silently
fail to return the photon.

Enforcing it requires knowing the model, which is why the helper's **first
argument is a required `Params` or `Model`**:

```python
capture_and_return_schedule(params_or_model, K, t_dark=...)
```

`Δ` and `u` live in `Controls`, so a caller reveals them just by calling. `G` does
not — it lives in `Params`, where nothing about the call exposes it. While that
context was optional, a caller could build the schedule without it, evolve under
a non-zero `G`, and get no error at all even though the return fails badly
(at `K = 0.85`, `t_dark = 1.9`: `G = 0.1` returns 0.830 of the photon, `G = 0.3`
returns 0.087, `G = 0.5` returns 0.254). Making the context mandatory closes that
path: there is no public call that constructs the advertised no-bath return
schedule without supplying what is needed to check it.

A non-zero `χ` is *not* rejected, and that is deliberate rather than an oversight.
The phase term is `ħ·(χ·u/2)·P_e·Z`, so forcing `u = 0` makes it vanish
identically whatever `χ` is; the suite verifies that a `χ ≠ 0`, `G = 0` model
still returns the photon exactly (2.6e-16). `capture_schedule` keeps its detuning,
because it makes no return claim.

## Time evolution and the exact global inverse

Controls are piecewise constant. Inside one segment `H` does not change, so the
exact solution is a single matrix exponential, computed from a Hermitian
eigendecomposition. There is no integrator step size to tune and no accumulated
truncation error.

The exact global inverse is exactly what the plan specifies:

```
H_global_reverse(s) = -H_forward(τ - s)
```

which for piecewise-constant controls means *same segments, reverse order, whole
Hamiltonian negated*. The minus sign is applied to the assembled matrix, so every
term reverses together — capture, detuning, phase history, mark coupling,
record-qubit energy and routing. This is the mathematical Loschmidt reversal. It
is not claimed to be experimentally available, and the *operational* reversal
generator (which reverses only permitted controls) is **not** implemented in this
slice.

The reversed propagator is built by independently diagonalising the negated
Hamiltonian, not by reusing or inverting the forward propagator. The recovery
result below is therefore a real numerical check rather than an identity that
holds by construction.

## Files

| File | Contents |
| --- | --- |
| `model.py` | Subsystem order, all 576×576 operators, `H(t)`, initial state, A↔B swap |
| `controls.py` | `Controls`, `Segment`, schedule builders, the exact global inverse |
| `simulate.py` | Unitary propagation, propagator construction, phase-aligned state comparison |
| `observables.py` | Populations, record excitation, phase-qubit `⟨Z⟩`, energy ledger, visibility |
| `verify.py` | The command-line verification script |
| `README.md` | This file |

Only NumPy and SciPy are used, both already present in the repository
environment (verified against NumPy 2.4.3 / SciPy 1.16.3 on Python 3.13).

## How to run it

From the repository root:

```bash
python3 -m first_mark_two_absorber.verify              # the checks
python3 -m first_mark_two_absorber.verify --verbose    # plus a detail line each
python3 -m first_mark_two_absorber.verify --demo       # plus a capture/return table
```

The script exits `0` only when every check passes. Each check prints the actual
numerical residual it measured and the tolerance it was held to, so a "PASS" can
be audited instead of taken on trust. The whole run takes a few seconds.

## What the checks confirm

All 18 pass. Measured residuals from a clean run:

| Check | Residual | Tolerance |
| --- | --- | --- |
| Operator shapes and tensor order | 0 | 1e-12 |
| Resonance relations | 0 | 1e-14 |
| Hermiticity of every generator term | 0 | 1e-13 |
| Propagator unitarity and `expm` cross-check | 6.7e-14 | 1e-11 |
| Norm conservation | 9.9e-14 | 1e-12 |
| Excitation number conserved, and `N=2` states exist | 1.4e-13 | 1e-11 |
| Fixed-Hamiltonian energy conservation | 2.1e-12 | 1e-11 |
| No-bath resonant capture vs analytic Rabi | 4.4e-14 | 1e-11 |
| No-bath detuned capture vs generalized Rabi | 1.0e-13 | 1e-11 |
| Single-arm recombiner routes into the output mode | 1.1e-13 | 1e-11 |
| Two-arm recombiner: bright/dark and `√2 J` rate | 6.1e-14 | 1e-11 |
| Resonant return helper refuses non-resonant use | 2.6e-16 | 1e-10 |
| No-bath return contract cannot be bypassed by omission | 2.6e-16 | 1e-10 |
| Mark/record exchange frequency | 5.0e-14 | 1e-11 |
| No-bath capture and coherent return | 6.7e-16 | 1e-10 |
| A/B symmetry | 1.9e-13 | 1e-12 |
| Exact global inverse, no bath | 1.1e-14 | 1e-11 |
| Exact global inverse, all couplings on | 3.4e-14 | 1e-11 |

Five of them are worth a sentence each.

**The analytic capture check.** On resonance, "photon in arm `i`" and "absorber
`i` excited" form an exactly degenerate two-level pair coupled with strength
`ħK`, and the two arms never talk to each other during capture. So the excited
population of each absorber must be `½ sin²(Kt)` — the one-half being the input
split. The simulation matches that to 4e-14 over two full cycles. A detuned
version matches the generalized Rabi formula to 1e-13, which tests the detuning
term as well.

**Capture and return.** Running the capture for a π/2 pulse area puts *all* the
amplitude into the absorbers; reversing the sign of `K` for the same pulse area
brings it back. The photon returns with population 1.000000000000, the prepared
A/B relative phase of 0.7 comes back as 0.70000000000000, and the fringe
visibility returns to 1. Between those, while the photon is inside the absorbers,
the arm-field visibility is reported as `nan` rather than 0 — it is undefined
there, and printing zero would read as lost coherence when nothing has been lost.

**The two-arm recombiner.** Both arms drive the *same* output mode, so routing is
interferometric and the recombination phases have to be tested as phases. With
`J_A = J_B = J`, the bright combination couples to the output with strength
`√2 ħJ` — enhanced over the single-arm rate — while the orthogonal dark
combination is annihilated by `H_route` and never leaves the arms. The check
scans the input phase across 13 values and matches
`P_out(t) = cos²(D/2)·sin²(√2 J t)` with `D = θ_in + θ_B − θ_A` to 6.1e-14, then
confirms the specific points: constructive input reaches output population
1.000000000000 at `t = π/(2√2 J)`, destructive input reaches 7.2e-28, and a pulse
timed at the *single-arm* rate `J` undershoots to 0.633 (it would be 1.0 if the
bright rate were `J`). The declared phases are deliberately distinct and
non-zero, so the test discriminates: mutating the code to conjugate the
recombination phases, to ignore them, to flip the sign of one arm, or to swap
`θ_A` and `θ_B` all make this check fail, while the older single-arm check passes
every one of those mutations.

**The resonant return helper refuses non-resonant use.** Flipping the sign of `K`
inverts the capture term and nothing else. A detuning, a phase-history coupling
or a mark coupling keeps running forward through capture, hold and return, and is
simply not undone. So `capture_and_return_schedule` raises `ResonantReturnError`
rather than building a schedule that its own name misdescribes. The check
confirms the guard fires, and then hand-builds the refused schedule to show the
guard is protecting against something real: at `K = 0.85`, `t_dark = 1.9`,
`Δ = 0.5` the photon comes back only to population 0.8226530516, with a
phase-aligned state error of 0.4312719942. Detuning belongs on
`capture_schedule`, which makes no return claim.

**...and the contract cannot be dodged by saying nothing.** A guard that only
fires when the caller volunteers the relevant information is not enforcement. The
second check covers the omission path: it confirms that calling the helper with
no model context is a `TypeError` rather than an unguarded schedule, that
non-zero `G` raises whether the context is a `Params`, a `Model`, or asymmetric
(`G_A = 0`, `G_B = 0.2`), and that a wrong first argument — a float, `None`, a
dict — fails loudly instead of accidentally satisfying the requirement. It then
sweeps `G` to quantify what is being prevented: 0.1 → 0.8298890013 returned,
0.3 → 0.0870779969, 0.5 → 0.2540097674. Eight contract probes, all holding.

**A/B symmetry, and one thing that is not an asymmetry.** With equal input and
identical parameters, the complete 576-component state is invariant under
exchanging the A and B halves at every sampled time (1.9e-13). With an unequal
input *phase*, per-arm populations still match exactly (3.0e-14) — provided the
recombiner is off. With the recombiner on they differ by about 2.3e-2, because
`J_A` and `J_B` feed the *same* output mode, so the two arm amplitudes interfere
there and per-arm populations depend on the relative phase by construction. That
is the interference `H_route` exists to produce, not a broken symmetry. The check
reports that number explicitly rather than hiding it.

## What this does NOT show

This is the part that matters most, and it is stated plainly.

- **Nothing here is a Born-rule result.** No probability is selected, no outcome
  is sampled, no branch is weighted. The model is unconditioned unitary
  evolution from start to finish.
- **Nothing here says one world is real.** It would be wrong to say the state
  space forbids a second excitation — it does not. The field and absorber
  factors are independent, so `|vac, e_A, e_B⟩` and `|photon-A, e_A, g_B⟩` are
  ordinary vectors in the 576-dimensional space, both with excitation number
  `N = 2`, and the verification suite builds them to show it. What keeps the
  calibration runs at `N = 1` is *dynamical*: the prepared state starts in the
  `N = 1` sector and the Hamiltonian conserves `N` (`[H, N] = 0`, verified), so
  the trajectory cannot leave it. That is a conservation law about excitation
  counting — not evidence that only one alternative is realised. The four-state
  field truncation is a separate matter: it limits which field modes exist, not
  how many excitations the tensor space can hold.
- **Nothing here demonstrates irreversibility or a real mark.** The mark/record
  check confirms only that the `H_mark` operator does what its definition says;
  a single record qubit oscillates back, and the exact global inverse recovers
  everything. Durable commitment is not in this model.
- **Nothing here is a phase-history or decoherence result.** The headline runs
  are the no-bath setting, `χ = G = 0`. The runs with those couplings switched on
  exist solely to exercise operator construction and the global inverse against a
  fully populated generator, and no physical conclusion is drawn from them.
- **No echo, random-reference, recorded-twin, or operational-reversal protocol is
  implemented in this slice.**

## Known limitations of this slice

- Only the exact global inverse is implemented. The operational reversal
  generator — the one that reverses only permitted controls and deliberately does
  *not* reverse `G_i` or the record qubit — is future work, and it is the branch
  where a meaningful recovery failure could appear. Nothing here tests it.
- Parameter values are calibration choices, not fitted physics. `Ω_r = 0.35ω` was
  picked only so the record gap is a clearly distinct scale; `K`, `J` and `G` are
  chosen per check so the relevant oscillation completes within the run. No
  parameter was tuned to make any effect appear.
- The unequal-path input (`weight_A`) is implemented in `initial_state` but is
  not exercised by the checks, since the plan puts the asymmetric case after the
  symmetric protocol passes.
- Finite-bath recurrence is not measured here, because the bath is off.

No manuscript claim should cite this model.
