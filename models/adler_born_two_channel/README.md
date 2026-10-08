# adler_born_two_channel — deterministic calibration and the analytic control

This package is the foundation of a test of one sharply bounded question:

> Does direct noisy Adler phase dynamics, applied to two competing populations of
> absorber clocks while a photon-pulse envelope opens and closes their Arnold
> tongues, produce first-commitment frequencies proportional to the **square** of
> their amplitude-linear peak couplings?

**None of that question is answered here.** What is in this package is the part
that has to be right *before* the question can be asked honestly. Five things:
the deterministic clock dynamics, calibrated against closed-form Adler theory;
kept structurally separate from them, the arithmetic of the prediction they will
later be tested against; the separation that will keep the future event
generator from ever reaching that prediction, described in
[The raw-process boundary](#the-raw-process-boundary); the supply of random
numbers that generator runs on, in
[The keyed noise source](#the-keyed-noise-source-and-the-brownian-tree); and —
new, and the subject of
[One clock, moving](#one-clock-moving-evolution-handoffs-and-the-dwell-clock) —
**one clock** evolving under that noise, with the exact tongue handoffs, the
timestamp dwell clock, and the two mechanism-negative controls.

Because of the third part, importing the package no longer imports the analytic
prediction. Nothing else about how you use the package has changed.

**A first clock is a model stopping event, not an outcome.** There is a
population now: `N` clocks race on one immutable grid and the channel stops at
the earliest dwell completion. That stop is an experimental endpoint defined by
this one-channel first-passage study. It is not a detector click, an absorption,
a measurement outcome, a unique actuality, a microscopic bath, one side of a
two-channel competition, or a Born-rule result, and no physical exclusivity is
claimed or implied by it. There is still no survival estimator, no commitment
law and no fitted exponent, so nothing here produces a quantity that could agree
or disagree with Born's rule.

---

## Read this first: what these controls do not show

**These controls do not demonstrate Born probabilities, stochastic selection,
commitment, exclusivity, or energy routing.**

There is now a phase that moves, a lock band, a dwell clock and a thing this
package calls a *commitment*. None of those is an outcome.

**A commitment is one clock staying inside a declared band around its moving
target for a declared time. It is not a detector click, an absorption, a
measurement outcome, or a Born event.** The word is the name of a dwell rule and
carries no physical interpretation beyond one; whether anything like it has to
do with quantum measurement is the question the programme is trying to ask, and
naming cannot answer it.

**Endpoint accounting can overestimate uninterrupted dwell: a path can leave the
band and return between two samples. The stationary killed-diffusion oracle
measures that error and meets its reference tolerances; the diagnostic
moving-band audit measures what it costs a dwell and currently returns a
`numerical_no_result` at its reference sample. Neither makes a pulsed result
here interpretable as physics: the population and the ledger now exist, and
the production numerical budget that ticket 07 froze for them is **not met** —
its disposition is `numerical_no_result` and no matrix was feasible among the
options it searched.**

**The killed-diffusion solver is the continuum authority and the moving-band
Brownian bridge is a frozen-drift diagnostic beside it, never a correction to the
primary ledger and never a way to choose a timestep or an exponent.**

**Failure to meet a frozen numerical budget is a numerical no-result. It blocks
interpretation; it falsifies nothing physical, and it is never a reason to widen
the budget.**

**There is a survival curve, a hazard and a fitted coupling exponent now, and
every scaling verdict this package can reach is `machinery_only`.** The
estimator recovers declared exponents from synthetic counts, and that is
evidence about the estimator. It is not evidence about the physics: every run
written here records `numerical_gate = "diagnostic_only"`, and
`compare.scaling_verdict` returns `machinery_only` for such a run **however good
the fit is**. There is no argument to any function in the analysis layer that
can bypass that gate.

**A statistically adequate exponent is not a Born-rule result.** It is the
coupling dependence of a stopping rule this study imposed on a declared model
bath, over one frozen population and one frozen pulse. An exponent near two is
not by itself support for the width-times-rate mechanism either: that
attribution requires the central-clock and width-only controls to *fail* to
reproduce it, and the strongest verdict `compare.causal_decision` can return is
`not_blocked` — there is no verdict spelled `supported`.

**Events divided by risk is a discrete conditional failure probability, not a
rate.** Only the width-divided hazard has units of inverse time. Every binned
row carries its risk count, its event count, its censoring count and its width,
so the two cannot be confused by omission.

**The bare relaxation-rate sum on the exact production clocks is the
finite-population comparator; the spacing-weighted continuum flux is a continuum
control and never a finite-population prediction.** They are different
quantities with different units, and a curve that claimed otherwise could not be
constructed.

A run writes what happened and stops. The analysis layer reads it afterwards and
writes nothing, and nothing in this package yet turns those records into an
outcome *frequency*, so nothing in it can agree or disagree with Born's rule.

**A run written by this package is diagnostic and is not production-cleared.**
Every manifest says so in its own fields. The stationary killed-diffusion gate
passed at its reference tolerances; the moving-band audit returned a frozen
`numerical_no_result`, pooled and in fourteen of fifteen per-regime cells; and
the production numerical budget ticket 07 froze is **not met** — its disposition
is a numerical no-result and no matrix was feasible among the options searched. The quadratic that appears in
`analytic.py` is a statement about an integral of the Adler relaxation rate over
a flat detuning grid — it is arithmetic about a synchronization model, and by
itself it is not evidence about quantum measurement.

The size of that error is worth seeing. **The stationary killed-diffusion
oracle that bounds it has passed its reference budget; the moving-band
diagnostic beside it has not, and no production budget exists for either.** For
a Brownian path with an upper boundary
at `a = 0.1`, endpoints both at `0.09`, `D = 0.01` and `dt = 0.01`, the
conditional probability that the path crossed the boundary in between is
`exp(-(a-x)(a-y)/(D·dt)) = exp(-1) = 0.368`. Endpoint accounting calls every
such segment continuously inside.

Even when the later stages are built and even if they reproduce the Born curve
exactly, the strongest available conclusion would be that *a two-channel noisy
synchronization race can reproduce the Born frequency curve under its stated
clock, bath, and stopping assumptions*. **Exclusivity, energy routing, and
Born's rule remain unproved.** The model contains no photon field, no
energy-bearing absorber levels, no amplification, no detector bias supply, and no
norm-preserving route for losing amplitude; its eventual "first winner stops the
trial" rule is bookkeeping imposed by hand and must never be described as derived
physical exclusivity.

---

## The physics, in plain English

A **clock** is an absorber whose phase we track relative to a reference. Left
alone it runs at its own natural frequency, a little fast or a little slow; that
offset is its **detuning**, `Δ`. Coupling it to the reference with strength `K`
gives its phase error the Adler equation

```
dφ/dt = Δ − K·sin(φ)
```

Everything else follows from that one line.

**The Arnold tongue.** If `K > |Δ|` the sine term can cancel the detuning and the
phase error stops moving: the clock has a stable phase to sit at, `φ* =
arcsin(Δ/K)`, and it is **eligible** to lock. If `K < |Δ|` no such point exists;
the phase never settles and the clock **slips**, running through a full turn
every `2π/√(Δ² − K²)`. The set of `(K, Δ)` where locking is possible is
wedge-shaped — the tongue — and it gets wider in proportion to `K`.

A point worth stating because it is easy to get backwards: **a clock does not
leave a stationary tongue just because its phase rotates.** The tongue boundary
lives in coupling-versus-detuning space, not in phase-angle space. Leaving it
requires `K` or `Δ` to actually change. Here that change comes from the photon
pulse.

**How fast is the lock?** Near `φ*` the error decays exponentially at the local
relaxation rate

```
λ(K, Δ) = √(K² − Δ²)
```

This is fastest at the centre of the tongue and falls to zero at its edge —
**critical slowing**. So the tongue is not uniform: a clock near the boundary is
technically eligible but relaxes so slowly it may never arrive anywhere.

**Why anyone thinks this might be a Born candidate.** Sum `λ` over a flat
spectrum of detunings and the two effects — a tongue that is *wider* in
proportion to `K`, and an interior that is *faster* in proportion to `K` —
multiply:

```
∫ over |Δ| < K of √(K² − Δ²) dΔ  =  (π/2)·K²
```

A universal constant times the **square** of an amplitude-linear coupling. That
is the entire reason the rate-weighted tongue is interesting. Note that tongue
width *alone* would give an amplitude-**linear** law; the second factor comes
from the interior speed, and the later stochastic stages exist to find out
whether a real race actually harvests it.

**The photon pulse.** A passing wave packet makes `K` rise from zero, peak, and
fall back to zero — a raised cosine, `K(t) = K_peak·½(1 + cos(2π(t−t₀)/T))`. Each
clock's tongue therefore *opens and then closes*:

```
before pulse        rising edge          near peak          falling edge       after pulse
no eligible clocks  tongue admits clocks widest tongue      clocks drop out    no eligible clocks
```

Now finite time bites. A clock near the boundary is admitted late, given a stable
phase to chase, and relaxes towards it at a crawl — and then the falling edge
takes the fixed point away before it arrives. The useful bookkeeping quantity is
the number of e-folds the pulse **offers** it:

```
A = ∫ over the eligibility window of √(K(t)² − Δ²) dt
```

`A` of order one or less means most release phases are too far away to be
reached in the time on offer, so the *probability* of arriving is small — not
that arrival is impossible. A clock released close to where the stable phase
appears can already be inside the lock band at the moment its tongue opens.
Finite exposure shifts the odds; it does not forbid the outcome for every
initial phase.

That distinction is not pedantry, and the suite measures it rather than
asserting it. The three outcomes are counted separately — inside the lock band
at the exact tongue-entry time, entered the band later, never entered it — and
the *causal* claim is tested by reversibility: give the identical release phases
a longer pulse and the same near-edge clocks come back (see the results below).
This is the mechanism most likely to *break* the quadratic when the later stages
run.

---

## Vocabulary the later stages will need

These categories are defined here and are not simulated yet. They are written
down now so nobody redefines them later:

| Term | Meaning |
| --- | --- |
| **eligible** | A stable Adler fixed point exists right now: `K(t) > \|Δ\|`. |
| **lock-failed** | Eligible for a while, and **never** reached the lock band. |
| **dwell-failed** | Reached the lock band, but lost it — or lost eligibility on the falling edge — before the fixed commitment dwell elapsed. |
| **unresolved** | The pulse ended with no commitment, for any of the above reasons. |

`lock-failed` is reserved for clocks that *never* reached the band. A clock that
did reach it may go on to commit or to be dwell-failed, and which of those
happens is a question for the dwell state machine — which does not exist yet.
This stage must not answer it on that machine's behalf, so reaching the band is
reported as reaching the band and nothing more.

The distinction matters because conditioning only on successful commitments
would manufacture a clean-looking channel ratio while hiding a strong pulse
duration dependence in the overall detection efficiency. Unresolved trials are
to be reported, never discarded.

---

## What is in each file

| File | Responsibility |
| --- | --- |
| `__init__.py` | Nothing. Deliberately: it contains no import statement at all. See [The raw-process boundary](#the-raw-process-boundary). |
| `model.py` | Phase wrapping and circular distance; eligibility; stable, unstable and boundary fixed points; the raised-cosine envelope; deterministic Adler stepping (RK4 and Euler). |
| `simulate.py` | Trajectories — stationary and pulsed — propagation to an exact time, and the measurements taken *from* them: fitted relaxation rate, eligibility entry/exit, phase-slip frequency. |
| `validation.py` | Argument checks, and no physics. Shared by every other module so that neither side of the isolation boundary has to import the other. |
| `analytic.py` | The prediction, and nothing else: local relaxation rate, fixed flat detuning grid, finite-grid sum, unconstrained exponent fit, refinement study. |
| `raw_config.py` | What a raw event run may be told: the clock population, the pulse, the noise strength, the lock band, the dwell, the step, the trial count, the stream namespace, the pilot/production label. Numbers and labels only. |
| `raw_ledger.py` | The close marker, the frozen three-field manifest, and the digest a raw ledger must carry before any comparison process may read it. The protocol, not a ledger. |
| `stochastic.py` | The keyed physical white-noise source: counter-keyed Gaussian phase kicks streamed in bounded blocks, the nested Brownian tree that splits one of them at an exact off-grid crossing, and the irregular elementary mesh both refinement levels read. Numbers, not dynamics. |
| `dynamics.py` | The stepper: Euler–Maruyama over the elementary mesh, the pulse train and its exact tongue schedule, the three eligibility states, the lifted-target coordinate handoffs, and the fixed-contraction width-only control. |
| `commitment.py` | The lock predicate — proximity **and** strict local contraction — and the timestamp dwell clock, with the four ledger categories. No phase, no coupling, no noise. |
| `raw_experiments.py` | Experiments S1 (deterministic limit), S2 (synthetic dwell sequences) and S4 (one stochastic clock), and the in-memory record each produces. One clock at a time. |
| `raw_runner.py` | The raw entry point: the input boundary, the configuration-to-noise-stream, criterion and path factories, and the honest list of what is not implemented. It does not run a race. |
| `killed_diffusion.py` | **Outside the raw graph.** The independent stationary killed-diffusion oracle — survival, exit-time quantiles and edge-resolved exit counts inside a fixed band — plus the closed-form limiting case it is calibrated against, and the frozen hash-checked budget machinery every numerical validation reports through. It imports `validation.py` and nothing else from this package. |
| `moving_band_audit.py` | **Outside the raw graph.** The diagnostic piecewise-linear Brownian-bridge replay for a *moving* band: admissible geometry, circular local representatives, the capped one-sided crossing probabilities, separately keyed auxiliary uniforms, paired coarse-level aggregation, and the reset-only replay whose commitments are a pathwise subset of the endpoint ones. |
| `observables.py` | **Outside the raw graph.** Risk sets, survival, discrete conditional failure, the width-divided time-rate hazard, cumulative hazard, the unresolved fraction with its three causes, and the winner mismatch with every co-completer retained. Three frozen conventions — event, censoring, tie — and a display record that always carries risk, events, censoring and width. |
| `analysis.py` | **Outside the raw graph.** The binomial complementary-log-log exponent fitted directly to committed/unresolved counts in every frozen coupling cell, its curvature and lack-of-fit tests and its curvature-identifiability gate, the closed vocabulary of `no valid exponent` reasons, provenance-bound master-trial identity, and the paired resampling that is the authoritative uncertainty — beside the deliberately wrong independent-cell control. |
| `compare.py` | **Outside the raw graph.** The comparison process: it opens closed hash-matching runs read-only, carries their fingerprints, and assembles them into a `RunSweep` whose frozen fields are checked as a set. It holds the four separately typed comparators, the frozen rising/falling and winner-mismatch holdout computed inside the decision path, the paired full-minus-control exponent contrast, the role-bound causal defeat rule, and the manifest-bound numerical gate that keeps every scaling verdict `machinery_only`. It is the only module allowed to hold a closed ledger and `analytic.py` at the same time. |
| `experiments.py` | **Outside the raw graph.** The feasibility, power and freeze layer: measured benchmarks of the real streamed race, the design-based power calculation for the binomial cloglog exponent and the paired full-minus-control contrast, the production numerical budget tied to planned statistical uncertainty, the finite primary matrix proposal and its cost, the pilot firewall with its count-only selection rule, and the dry-run production manifest proposal — which is structurally unable to clear or sign itself. |
| `verify.py` | The registered checks, with residuals and tolerances, and the CLI. |

### The isolation rule

`model.py` and `simulate.py` **never import, name, or call anything in
`analytic.py`.** The dynamics must not be able to see the answer they are meant
to be tested against, or a later "confirmation" would be circular.

This is enforced structurally, not by convention: `verify.py` parses both files
with `ast` and fails if they import anything outside an allowlist or mention any
analytic name anywhere in their syntax trees. `analytic.py` is checked in the
other direction too — it imports no dynamics module, which is what makes it an
independent prediction rather than a restatement of the simulation.

One consequence: `analytic.py` re-derives the eligibility predicate `K > |Δ|`
locally instead of importing it. That duplication is deliberate, and it is made
safe by a check that compares the two implementations over 200,000 random
`(K, Δ)` pairs plus the exact boundary cases.

The two sides share exactly one module, `validation.py`, which holds the
argument checks. That exemption is *earned* rather than granted: `verify.py`
parses `validation.py` and requires it to be physics-free — no identifier
anywhere in its syntax tree may name a coupling, a detuning, a phase, a tongue,
a rate or a pulse. A shared module that cannot express a physical quantity
cannot carry one between the prediction and the simulation.

The dynamics also receive a **coupling**, never an amplitude squared, an
intensity, or a precomputed hazard. `verify.py` checks the source for those too.

### The raw-process boundary

The rule above names two files. That was enough while the package contained two
dynamics modules. It is not enough for what comes next: the stochastic stage
adds a noise generator, a dwell state machine, a ledger writer, a
killed-diffusion oracle, a bridge audit, an analysis layer and a comparison
layer, and only the last two may see the prediction. A check that names a fixed
pair of files cannot notice a third file importing a fourth.

So the rule is now stated as a graph, and one thing had to change first.

**What changed.** `__init__.py` used to say `from . import analytic, model,
simulate, validation`. That meant importing *any* submodule — including a raw
event module whose entire purpose is not to know the answer — initialised the
package with the prediction already loaded. The structural check parsed the
dynamics and found nothing wrong, because nothing was wrong with the dynamics.
The leak was in the package root.

`__init__.py` now contains no import statement at all.

**What still works.** All three explicit forms, unchanged:

```python
from adler_born_two_channel import analytic
from adler_born_two_channel.analytic import tongue_rate_sum
import adler_born_two_channel.analytic
```

Each of those is a statement that *names* `analytic`, so each is visible to a
reader, to a reviewer, and to the syntax-tree scan that forbids the raw modules
from writing one.

**What no longer works**, and this is the whole point:

```python
import adler_born_two_channel
adler_born_two_channel.analytic          # AttributeError
```

Attribute access after a bare package import used to succeed *because* of the
eager import, and it is the one route to the prediction that leaves no trace in
the file that takes it. Nothing outside the package used it.

**The graph.** `verify.py` declares three roots — `__init__.py`,
`raw_runner.py`, `raw_ledger.py` — and walks every permitted relative import
until it runs out of modules. Everything reached is scanned; anything named
`raw_*.py` that the walk cannot reach is a failure in itself, so a new raw
module cannot escape by not being imported yet.

A raw module may import `__future__`, `dataclasses`, `typing`, `numpy`, `math`
`hashlib` and `json`, and may import `validation`, `model` and its sibling raw modules.
It may not name the prediction, the killed-diffusion oracle, the moving-band
audit, the analysis, comparison or reporting layers, an intensity, an amplitude
squared or a prescribed hazard. `importlib`, `import_module`, `__import__`,
`eval` and `exec` are on the ban list too, because the scan is over import
*syntax* and those have none.

**Every binding position, not just uses.** The scan checks names where they are
*bound* as well as where they are used: definitions, arguments, keyword
arguments, `global`/`nonlocal` declarations, `except … as` targets, the object
names and aliases of every import, and the binding fields of structural
patterns — `case … as hazard`, `case [1, *hazard]`, `case {…, **hazard}` and
`case Point(intensity=v)`. An independent review defeated the earlier use-only
version in one line — define `def hazard(...)` in an *allowed* module and write
`from .raw_config import hazard as event_rule` — and the walk reported nothing;
a later review did the same through a mapping pattern's `**rest`. Match patterns
bind through plain string fields on the pattern node, so nothing there is an
`ast.Name` and a use-site scan sees a binding only if the name is later read.

Allowed **module** traversal and forbidden **object** checking are two separate
questions: `from . import validation` names a module and is fine,
`from .raw_config import hazard` names an object and is not, whatever module it
came from and whatever it is renamed to on the way in.

**The square rule, and its limits.** A raw module may not write a square:

- no exponentiation, as an expression or an augmented assignment;
- no value multiplied by itself — `coupling * coupling`, `coupling *= coupling`;
- none of the named power helpers `pow`, `square`, `power`, `float_power`; and
- none of the named multiplications `multiply`, `matmul`, `mul`, `dot`, `inner`,
  `outer` applied to two structurally identical operands, which is what
  `np.multiply(coupling, coupling)` is. These are not banned outright, because
  multiplying two *different* things is ordinary arithmetic.

The two multiplicands are identified before they are compared, so the option
arguments a call may also carry — `out`, `where`, `dtype`, `casting` — are
ignored, and the operands may be passed by name. `np.multiply(x, x, where=m)`,
`np.multiply(x1=x, x2=x)` and `np.multiply(x, x2=x)` are all caught;
`np.multiply(x, y)`, `np.multiply(x1=x, x2=y, where=m)` and
`np.multiply(x, x2=y)` are all fine. Where the operands *cannot* be picked out
unambiguously — a starred argument, a `**` mapping, or two spellings of the same
operand — nothing is claimed at all. Four of the mutation battery's six control
packages are valid unequal multiplications, so a rule widened until it flagged
ordinary arithmetic would fail the suite rather than pass it.

The same rule applies to `model.py` and `simulate.py`.

The first version banned only the `**` operator and claimed this meant a raw
module could not form a squared coupling. That was false, and two review rounds
demonstrated it — first with `coupling * coupling`, then with
`np.multiply(coupling, coupling)`. The rule is wider now and the claim is
narrower: it is a **syntactic** guard against *writing* a square, not a proof
that no squared value can arise. `a = coupling; b = coupling; a * b` still
passes, as would a product accumulated in a loop, an operand reached through an
alias, or a helper called with its arguments spread from a list. **No data-flow
or semantic analysis is performed or implied.** The two-factor rate-weighted
tongue is the hypothesis under test, so a raw module forming a square deserves
to be noticed — nothing here can promise it always will be.

**The runtime half.** A syntax-tree scan cannot see an import that happens at
run time, and cannot see the package root executing on the way to a submodule.
So `verify.py` also starts a fresh interpreter and asks it what got loaded:
`import adler_born_two_channel` must leave exactly one package module in
`sys.modules` and bind no public attribute, and `import
adler_born_two_channel.raw_runner` must load only modules on the same
allowlist. Then it imports `analytic` explicitly, and requires that to *appear*
— which is what shows the two lists before it were measuring something.

**The caller's half.** Structural isolation is worth nothing if the caller can
carry the prediction in by hand. A configuration object with arbitrary
attributes lets a well-meaning user hand a raw run the very curve it is meant to
produce independently, after which the "measured" law is a restatement of its
own input. So `raw_config.py` accepts numbers, labels and its own clock grid,
and refuses everything else: callables by kind, arrays by kind (an intensity
weight, an amplitude squared, a per-clock rate and a prescribed hazard all
arrive as an array of numbers), and 42 specific names — `hazard`, `predictor`,
`intensity`, `tongue_rate_sum`, `oracle`, `expected_exponent` and the rest — by
name, with an error that says why rather than "unexpected keyword argument".

**Exactly the declared type, not merely an instance of it.** The entry point
`validate_raw_config` requires the type to be *exactly* `RawEventConfig`, and
rebuilds one by naming its eleven fields individually. Both halves matter, and
both exist because the review got a predictor through the boundary without
either.

It used to check `isinstance` and rebuild with `dataclasses.replace`, which
preserves the dynamic class. A frozen subclass declaring `predictor: object`
therefore passed, was rebuilt as the same subclass, and arrived at the raw
process still holding an arbitrary opaque object — the one input refused by
name, by kind *and* as an unknown key, reaching the code that exists
specifically to not have one. Subclassing was the way round all of it. The same
exact-type rule applies to the nested `RawClockGrid`, which a subclass could
otherwise override to return any detunings at all.

The rebuild is also what stops a *stale* configuration: these dataclasses are
frozen against ordinary assignment and not against `object.__setattr__`, so the
entry point re-runs every field check and every cross-field rule at the moment
the run begins. A NaN peak coupling edited in after construction would otherwise
make `K > |Δ|` false everywhere and produce a clean, physical-looking, entirely
fictitious "no clock was ever admitted" result.

**The clock grid is duplicated on purpose.** `raw_config.RawClockGrid`
re-derives the symmetric midpoint grid rather than importing
`analytic.FlatGrid`, for the same reason `analytic.py` re-derives `K > |Δ|`.
The duplication is made safe the same way: `verify.py` compares the two, at a
tolerance of **zero**, over even and odd counts from 2 to 4000 clocks across six
decades of support. One divergence is deliberate and is asserted rather than
hidden — `FlatGrid` refuses a single point because one point is not a
quadrature, while `RawClockGrid` allows one because a single clock at zero
detuning is the mandatory central-clock negative control.

**Time, as well as space.** A comparison process that opens a ledger while the
raw run is still writing it becomes part of the experiment, silently: unwritten
trials read as missing rather than unresolved, an interrupted sweep reads as a
completed one with fewer commitments, and a rerun that appended mixes two random
streams into one apparent result. So `raw_ledger.py` defines a close marker —
schema version, a per-table row count and content digest, and the manifest
digest — and a read gate:

```python
require_closed_ledger(marker, manifest, ledger, clocks, shadow) -> ClosedLedger
```

Five inputs, because a closed run is five files: the marker, the manifest bytes,
and the three tables. `ledger` is one row per trial, `clocks` the authoritative
trial-by-clock grid the trial rows are *derived* from, and `shadow` the
predeclared shadow trials. What comes back is a `ClosedLedger` — the rebuilt
marker, the parsed manifest and the three parsed tables — so a caller never
holds a reference to something the gate did not validate, and never re-parses
bytes the gate has not vouched for.

It checks **eight** things, in this order, and each is here because skipping it
lets a specific wrong number look right. This is a schema-v3 contract: the
package's earliest revision took a single payload, checked four things and
returned a `CloseMarker`, and that shape no longer exists.

1. **The marker is exactly a `CloseMarker`**, not an instance of one. A subclass
   could override anything the seven checks below read.
2. **The marker is rebuilt here**, which re-runs every field check. A marker is
   frozen against ordinary assignment and not against `object.__setattr__`, and
   one edited afterwards used to sail straight through — the gate would accept
   and return a marker simultaneously asserting `schema_version=999`,
   `row_count=0` and `manifest_digest="not-a-digest"`.
3. **The manifest bytes** hash to `manifest_digest`, and are verified *before*
   they are parsed, so everything read out of them is read out of bytes this
   marker has already vouched for. Until this existed nothing bound the two
   together, and a ledger could be compared against a prediction computed from a
   *different* run's parameters.
4. **Every table's bytes** hash to that table's own declared digest. Truncation,
   extension, emptiness, a hand edit and one flipped byte each raise.
5. **Every table is parsed**, header and all. A renamed or reordered column, a
   row of the wrong width, a number a spreadsheet rewrote, or a category outside
   the declared four fails here.
6. **The counts are reconciled three ways**: what the manifest froze, what the
   marker declares, and how many rows are actually present. This is what
   replaced the ticket-01 manifest-level-only check, which could establish only
   that a marker agreed with a manifest about a number neither had counted.
7. **The row identities are exactly the declared ones.** Trials `0 .. trials-1`,
   each once; the clock table exactly the full trial-by-clock grid, each pair
   once and in order; shadow rows exactly the predeclared shadow grid. A ledger
   holding one trial twice and another not at all has the right row count and is
   not the run's ledger.
8. **Every row is internally possible, and the tables agree with each other.** A
   cell can be the right kind and still be a lie: one word changed from
   `dwell_failed` to `committed` passes the digests — a forger recomputes them —
   the schema and every count, and turns an unresolved trial into a commitment.
   The trial rows are restated from the authoritative clock rows here, so a
   reader establishes the writer's own invariants over the bytes in front of it.

None of the eight is an argument. A previous revision took the expected count
from the caller, which closed nothing — a caller holding a marker that claimed
999 trials passed 999, and the two agreed with each other about a number neither
had any right to. There is now no value a caller can supply that changes the
outcome.

The manifest is therefore a type rather than an opaque blob. `RawManifest` takes
one mapping and requires **exactly** the schema's 45 keys — no more and no
fewer — of exactly the declared kinds, and stores them as an ordered tuple of
normalized pairs. Among them:

| Field | Meaning |
| --- | --- |
| `schema_version`, `ledger_schema_version` | Must equal `MANIFEST_SCHEMA_VERSION` and `LEDGER_SCHEMA_VERSION`. Versioned separately, since the manifest and the table layout will not grow at the same rate. |
| `run_label` | `"pilot"` or `"production"`, imported from `raw_config` rather than restated, so the firewall cannot come to mean different things on the two sides of it. |
| `trials`, `shadow_trials`, `clocks` | The authoritative counts. The gate reads them from here and reconciles them against the marker and the rows. |
| `ledger_columns`, `clock_columns`, `shadow_columns` | The three table layouts, pinned to exactly what this code writes and reads. |
| `model`, `drive`, `peak_coupling`, `detunings`, `fixed_contraction_rate`, … | The physical configuration a comparison would compute a prediction from. |
| `numerical_gate`, `validation_note`, … | The status fields that decide whether a run may be read as physics at all. |

An unknown key and a missing key are both refused, which is what forces a schema
change to be a deliberate version bump rather than a reader quietly ignoring
fields — reading the keys it recognises and skipping the rest is how a run gets
compared against parameters that are not its own.

`canonical_manifest_bytes` gives the writer one spelling (sorted keys, no
insignificant whitespace, ASCII). `parse_manifest` is strict about the schema
and relaxed about the formatting, because the digest is over the bytes as they
were actually written and reformatting them would break it.

It also refuses a manifest that names a member **twice**, and it has to do that
*during* parsing: `dict(pairs)` silently keeps the last value for a repeated
key, so once the parser has returned, `{"trials": 2, "trials": 999}` and
`{"trials": 999}` are the same object. A manifest declaring `trials` twice would
otherwise bind its digest to bytes saying one thing and its parsed count to
another, and which count the gate enforced would depend on the order the two
appeared in — the reader and anyone looking at the file could disagree about the
same run, silently. Whether the repeated values happen to *agree* is not
consulted: a file that names a member twice is malformed however it is read.

There is a second door, for a snapshot that is already open.
`require_unmutated_ledger(ledger) -> ClosedLedger` takes a parsed `ClosedLedger`,
re-emits canonical manifest and table bytes from it, puts *those* through
`require_closed_ledger` above, and returns what that gate parsed. A
`ClosedLedger` is frozen only in the reference: its rows are ordinary
dictionaries, and an earlier guard that checked a snapshot against a
caller-supplied marker proved only self-consistency — four rehashed forgeries
passed it. Running the full gate is what closed that, and the analysis layer
binds the object the gate *returns* rather than the one it was handed.

None of this is authentication. The digests establish that the bytes a reader
holds are the bytes a writer finished with; anyone able to rewrite the ledger can
recompute them. It is protection against truncation, interleaving, mismatched
pairs and accidental edit, which is what actually goes wrong. **It is not a
signature.**

**And it is not a sandbox.** All of this is protection against *accidental
manufacture* — against the prediction becoming reachable from an event generator
without anybody deciding it should. Nothing here treats the caller as an
adversary, and a determined one has `importlib`.

### The keyed noise source and the Brownian tree

The stochastic stage needs random numbers, and *which* random numbers is a
design decision rather than a detail. A sequential generator — draw for each
live clock as the loop reaches it — makes a clock's noise depend on how many
other clocks were still running, on the batch size, on which trials stopped
early, and on the order the loop happened to take. Every one of those is a
property of the *program*, and none of them is physics. Two runs that are
supposed to differ only in coupling would then differ in their noise as well,
and a paired comparison would be measuring both at once.

So a kick belongs to a **key**, not to a run:

```
dk-phase-noise/v3/physical/<dataset>/mesh=<origin>:<step>/trial=<t>/clock=<c>/step=<j>
```

and a split normal — see below — belongs to a longer one:

```
dk-phase-noise/v3/physical/<dataset>/mesh=<origin>:<step>/trial=<t>/clock=<c>
        /step=<j>/node=<L/R ancestry>/crossing=<clock>:<edge>:<pulse>:<time>/split
```

Both are ASCII, and every field in them is exact: the mesh origin and step and
the crossing time are written as **hexadecimal floats**, which round-trip
bit for bit, because decimal formatting would let two distinguishable meshes or
two distinguishable crossings serialize to one key. Every integer in a key —
trial, clock, step — obeys one grammar, `[0, 2⁶⁴ − 1]`, whether it arrives on
its own or inside an array, with no conversion through a narrower type in
either path.

**The address and the seed are two different strings, and v1 conflated them.**
Leaves are generated 1024 at a time, and v1 published the *block seed*
(`.../block=<j // 1024>`) as "the key of one finest-grid kick". All 1024 leaves
of a block shared it, so a manifest storing such a key could not say which
value it meant. In v2 `leaf_key` is the per-leaf address ending in `step=<j>`,
`leaf_from_key` turns that address back into the leaf it names, and the block
seed is private. Blocking is still how a value is *computed* and never what it
depends on, so a caller asking for three values and a caller asking for three
million still get the same numbers.

The version has moved twice, and for the same reason each time: a change to
what a key *means* has to be visible as a change of dataset rather than
silently producing different numbers under the same names. **v2** separated the
address from the block seed; **v3** changed the duration law and the way the
Brownian scale is computed (both below). Every keyed value changes with each
bump. No run has been recorded under any version, so nothing needed migrating;
the statistical residuals in the suite move with the bump and are reported
below.

The public surface, in full: `leaf_key` and `leaf_from_key` address and
reconstruct one kick; `stream_prefix` is the part of an address that names the
dataset and mesh; `standard_normals` and `leaf_kicks` are the unscaled and
scaled root paths; `stream_leaf_blocks` streams the latter in bounded windows;
`step_bounds`, `step_containing`, `time_at` and `identity` are the mesh;
`step_scale` is one step's realized Brownian scale; `elementary_leaves`,
`stream_elementary_leaves` and `coarse_kicks` are the tree. A check requires
every one of those names to appear in this file, so a new public method cannot
arrive undocumented.

A key becomes a number by one fixed route:

```
key -> SHA-256 -> four 64-bit little-endian words -> SeedSequence
    -> Philox counter-based bit generator -> Generator.standard_normal()
```

and the kick is that standard normal scaled by `sqrt(2·D_phase·h)`. **Diffusion
strength and duration scale a keyed normal; they do not select one.** That is
why zero diffusion is exactly zero and consumes no key at all, and why the same
keys at two noise strengths are the same realization at two noise levels.

What is *absent* from the key is the point of it: no coupling, no detuning, no
amplitude, no lock tolerance, no dwell, no batch size, no iteration order, no
clock count, and **no model label**. The full Adler run and the fixed-contraction
width-only control therefore share one physical path by construction rather than
by convention — there is nowhere in the API to tell the noise which model is
asking. The diagnostic moving-band audit, when it is written, gets its own
namespace and cannot reach these values; the `physical` segment is reserved in
the schema now, before the code that needs the separation exists.

**One tree, read two ways.** Coarser refinement levels are *sums* of the same
finest leaves, never separate draws, so a refinement ladder measures timestep
error and not timestep error plus a fresh sample. The summation order is fixed
(strictly left to right in time within each group), so a coarse increment is
reproducible bit for bit, and a partial trailing group is refused rather than
dropped or kept — both choices would make two levels disagree for a reason that
is not physics.

**Where a step ends is where the next one starts.** A step's interval is
`(time_at(i), time_at(i+1))` — never `start + step`, computed separately. The
two are different expressions for one instant and different floats for it: on
`FinestMesh(1.0, 0.1)`, `time_at(1) + step` is `1.2000000000000002` while
`time_at(2)` is `1.2`, so step 1 claimed to end after step 2 began. A crossing
at exactly `1.2` then satisfied *both* steps' containment tests, was split
twice, and produced a spurious child `2.2e-16` wide — seven leaves where six
were due. Over 200 steps of that mesh there were 31 such overlaps and 36 gaps
in the other direction, where a crossing would have belonged to no step at all.
Taking both boundaries from `time_at` makes overlaps and gaps impossible rather
than rare, and `step_containing` is then the single half-open routing rule: a
crossing exactly on a boundary belongs to the step that *starts* there.

**And the realized width is the physical duration.** A step's width can differ
from the declared `step` — on `FinestMesh(1e15, 0.2)` the grid times snap to
multiples of `ulp(1e15) = 0.125`, so the nominal 0.2 realizes as an alternating
**0.25 / 0.125**, neither of them 0.2 and their mean exactly 0.2. The tree used
that realized width for `TreeLeaf.duration`, for the split fraction and for the
coarse aggregation, and the *nominal* step for the kick variance: an interval
recorded as lasting 0.25 carried a kick drawn for 0.2. Averaged over many steps
the two agree, which is precisely why nothing noticed.

`step_scale(index)` now takes the duration from `step_bounds`, so one duration
governs the variance, the geometry and the record. It is a method rather than
the property it used to be, because a mesh does not have *one* leaf scale. The
suite groups kicks by their own step's realized width and holds each group to
`2·D·h` for that width; the nominal law fails the same test by 60–100%.

**The scale is factored, because the product is not survivable.** `2·D·h`
underflows to exactly zero long before its square root does: at `D = 0.2,
h = ulp(0)` the product is `0.0` while the scale `1.4e-162` is an ordinary
double, and at `D = h = 1e-200` the product is `0.0` while the scale is
`1.4e-200`. A zero scale is not a small scale — it is the deterministic limit —
so the old spelling silently replaced a physical bath with no bath, and the
exact-zero tests that protect that limit *agreed with it*. So

```
sqrt(2·D·h)  is computed as  sqrt(2) · sqrt(D) · sqrt(h)
```

which can only overflow when the true scale is genuinely unrepresentable, since
`sqrt(x)` is in range whenever `x` is. The cost is one extra rounding, so
results differ from the product spelling in the last place — part of why v3 is
a new dataset. A nonzero scale that still underflows to zero is **refused**,
not returned as the deterministic limit: those are different physical
statements.

**An empty window is refused before any randomness.** The generator addresses
steps by integer index and never looked at how long they were, so on
`FinestMesh(1e16, 0.1)` — where consecutive grid times are the *same double* —
it returned plausible kicks for a window describing no time at all. Every
public path (`leaf_kicks`, `standard_normals`, `stream_leaf_blocks`,
`elementary_leaves`, `stream_elementary_leaves`, `coarse_kicks`, `step_scale`)
now validates every requested step's width first, in full and in bounded
chunks, so a window whose *later* step collapses fails at the request rather
than after four steps of randomness have been drawn. The two generator-returning
paths raise at the call, not at the first `next()`.

**The elementary mesh is the uniform grid *union* the exact crossings**, and
`coarse_kicks` sums it. Its crossing schedule is a required argument, not an
option — an empty schedule is the uniform case — so there is no public
signature through which a crossing-bearing window can be aggregated from
unsplit roots. `elementary_leaves` and `stream_elementary_leaves` expose the
same sequence directly, and the coarse path holds one group of leaves at a time.

That argument is mandatory because the failure it prevents is invisible: root
sums and leaf sums agree numerically, frequently **bit for bit** and at worst by
about one unit in the last place. No tolerance can separate them, so a
refinement study could sum roots while the fine dynamics consumed crossing
leaves and no numeric check would notice. What separates the two is that
summing the tree *derives the split normals* — the suite counts them, and the
root-summing path, kept private as `_coarse_kicks_uniform`, derives none.

**The tree owns split identity.** There is no public split at all: not
`split_kick`, not `split_node`, not `split_key`, not `split_leaf`. A node path
is half of a split key and the node's geometry is not in that key, so any public
method taking a `node_path` lets one ancestry be paired with two different
intervals — an independent review split step 2 at `node_path='R'` over
`[t₂, t₂+h]` and again over `[t₂, t₂+0.8h]`, got one key, and drew one normal
for two different conditional laws. Validation cannot fix that, because nothing
in the arguments says which node `'R'` is meant to be. Ownership can: the walk
in `elementary_leaves` produces every node path it uses, and it takes the
clock's whole frozen schedule rather than a per-step subset a caller chose.

**The split at an exact crossing.** A tongue entry or exit lands inside an
assigned finest step whose kick has already been drawn. The two sub-intervals
need kicks with the variances their own durations imply, *and* they must still
add up to the parent. Dividing the parent in proportion to the two durations
does the second and fails the first — it gives the left child variance
`α²·2Dh` where the Brownian path requires `α·2Dh`. The construction conditions
on the parent instead:

```
left  = α·parent + sqrt(2·D·h·α·(1−α))·split_normal
right = parent − left
```

so `Var(left) = α·2Dh`, `Var(right) = (1−α)·2Dh`, `Cov(left, right) = 0` — and,
since the pair is jointly Gaussian, zero covariance is unconditional
independence — while `left + right = parent` by construction. The second child
is always the **residual**; sampling it independently gives the right marginals
and destroys the parent. Multiple crossings inside one step recurse
chronologically: the earliest crossing splits the step, the left child is then
finished, and the rest recurse into the right child. The node path records that
`L`/`R` ancestry and is part of the key, so a step split four times uses four
different split normals and never reuses one.

**Geometry is derived, never supplied.** The split works its fraction out from
`crossing.time` and the node's exact endpoints; nothing accepts a fraction
alongside a crossing that says something else. An earlier public primitive did,
and an independent review split the step spanning `[0.04, 0.05]` with a crossing
at `t = 99.0`, and split at three quarters with an identity that said one
quarter — both consuming a split normal keyed to the inconsistent identity.
Every geometry check runs *before* any key is derived, at every level of the
tree including nested node intervals, and the suite asserts that a refused call
consumed nothing — with the unchecked private primitive as the control that
shows the counter was watching.

**An address has one spelling.** `leaf_from_key` reconstructs the canonical key
from whatever it parsed and requires byte-for-byte equality with its input.
`int()` is far more generous than a canonical grammar: `step=007`, `step=+7`,
`trial=01` and `step= 7` were four aliases for one address, and `step=7_0`
parsed as **seventy** and silently returned a different leaf than the caller had
named.

**Conservation, stated honestly.** `left + right = parent` is exact in real
arithmetic and *not* bit-exact in floating point: `fl(left + fl(parent − left))`
can differ from `parent` in the last place or two. So the contract is a
scale-aware **4-unit-in-the-last-place** bound at the largest magnitude
involved, which is a deterministic bound and never a statistical tolerance — a
split that conserves the parent on average is a broken split. The measured
worst case over 40,000 splits at α = ¼ and ½, plus 256 each at seven fractions
from 1e−12 to 1 − 1e−12, is **0.25 of that bound**.

The bound comes from `math.ulp` and is required to be finite before it is
returned. `np.spacing` computes the *forward* gap and overflows at the last
binade, so `conservation_tolerance(DBL_MAX)` used to be `inf` — a bound that
accepts every possible error is not a bound, and the suite keeps `np.spacing`
around as the control that shows the difference.

**The accepted scale domain is
`[2.2250738585072014e-308, 1.3407807929942596e+154]`** — the normal binary64
positives, capped at **`√DBL_MAX`** — plus the exact `0.0` of the deterministic
limit, and one policy enforces it in the constructor, the root window, the tree
node, the split bridge, the stream and the coarse path alike.

The upper end is `√DBL_MAX` and not `DBL_MAX`, and that is the difference
between a finite kick and an infinite one. A kick is `scale·z`, so admitting
scales up to `DBL_MAX` admits `D = h = DBL_MAX/2`, whose scale is `1.27e308`;
any `|z| > 1.41` then overflows to `inf`, silently, and a reproduction returns
a finite number only if that particular draw happened to be small. Capped at
`1.34e154` the product is finite even for a normal as large as the scale
itself, and real standard normals do not exceed about 8 — a margin of 153
orders of magnitude. `DBL_MAX × ¼` and `1e308 × 0.1` still accept; `DBL_MAX ×
DBL_MAX`, `(DBL_MAX/2)²` and `1e200 × 1e200` do not.

The lower end is where the law breaks in the other direction: a **subnormal** has
fewer than 53 significant bits — one, at the very bottom — so multiplying a
standard normal by one does not give a Gaussian, it gives a handful of
quantized magnitudes. At `D = h = ulp(0)` the required scale is `6.987e-324`,
the nearest double is `4.94e-324` — **29% low**, variance halved — the kicks
take exactly **two** magnitudes over 64 draws, and the conditional split's
bridge term `scale·√α·√(1−α)` rounds to **exactly zero** at `α = ¼`, silently
becoming the proportional subdivision the mutation contract exists to reject.
That is not a small numerical error; it is a different law, and it is refused
before any key is derived. The suite keeps the measurement rather than just the
refusal.

A subnormal *input* is fine when the output scale is normal: `D = ulp(0)` with
`h = 1` gives `3.14e-162`, which agrees with a 60-digit `decimal` reference to
the last bit. The policy is about the scale, not the operands.

**The boundary decision follows the mathematics, not the rounding.** The
factored scale is accurate to an ulp or two — plenty for generating, not enough
for *classifying*. At `D = ulp(0)`, `h = nextafter(2⁻⁹⁷¹, 0)` the exact target
lies below the smallest normal double while the computation rounds **up** to
exactly it, so a policy reading its own arithmetic would admit a target the
mathematics rejects. Classification is therefore exact integer arithmetic:
significands multiplied, exponents added, and `sqrt(S) < B` settled with
`math.isqrt` — the unstable product is never formed, and nothing is squared
(the isolation scan forbids a raw module writing a square, and a constant is
not worth an exception).

Admission needs **both** the exact target and the computed scale to be in
range, and the two are held to different intervals. The **exact mathematical
target** must lie in `[MIN_NORMAL, √DBL_MAX)` — lower end closed, **upper end
open**, so a target exactly at `√DBL_MAX` is refused, the single point where
the mathematical product equals `DBL_MAX`. The **computed binary64 scale** must
separately lie in the closed representable interval
`[2.2250738585072014e-308, 1.3407807929942596e+154]`, because that is the
number a run actually multiplies by.

Requiring both is deliberately conservative at the boundaries: a target inside
the domain whose computed scale rounds one ulp outside is refused, and so is a
target outside whose computation rounds inside.

**Every scale a schedule needs is settled before any key is drawn.** The window
check covers the root scales; a split needs the node's own scale and the
conditional bridge, which depend on the sub-interval a crossing carves out. So
the tree's geometry — a function of the crossing times alone — is walked with
no randomness at all before anything is generated. Without it an inadmissible
bridge was found *after* the parent kick had been drawn, and on the streaming
path not until the first `next()`, long after the call that should have refused
it.

The domain is checked by computing the scale, not by asking whether `2·D·h` is
finite. That older question was wrong in both directions: left-associated, the
product overflows at `D = DBL_MAX` before the small duration is ever applied, so
`DBL_MAX` at a quarter step was refused although its scale `9.48e153` is an
ordinary double — as was `1e308` at a tenth — while a product that merely
*underflowed* passed and delivered a scale that no longer carried a Gaussian.
Both directions are now probed: `DBL_MAX × ¼` and `1e308 × 0.1` accept,
`DBL_MAX × DBL_MAX` and the two subnormal-scale cases refuse.

Two cases are bit-exact rather than tolerant, and are asserted as such: a
zero-duration child is exactly `0.0`, and a zero-diffusion tree is exactly zero
throughout — parent, children and every nested descendant — with **zero key
derivations**, counted rather than inferred. An endpoint crossing (`α = 0` or
`α = 1`) likewise draws nothing: it reveals no new information about the path
inside the interval. A nonzero kick where the variance is zero is refused, not
divided, because it means two trees have been mixed.

**White noise here is a minimal effective bath: no microscopic or
quantum-vacuum origin is claimed for it.** It is the smallest assumption that
makes an Adler clock stochastic at all. Because it is additive, the
Itô-versus-Stratonovich distinction does not change this model — a statement
about the equation, not a further physical conclusion. It is also an
infinite-bandwidth idealization whose numerical cutoff is the finest timestep,
which is why timestep convergence is a scientific requirement of the experiment
and not merely a software check.

**What the reproducibility claim covers.** NumPy guarantees its bit generators'
streams and does not guarantee that `standard_normal` will turn a stream into
the same normals in a future major version. So: identical keys give identical
kicks **in a recorded (Python, NumPy) environment**, which is what a run
manifest exists to record — and identical *structure*, meaning independence,
invariance to batching and ordering, and parent conservation, in any environment
at all.

**Memory.** Nothing materializes the trial-by-clock-by-time cube. Streaming 64
clocks × 200,000 steps in windows of 4096 peaks at **6.3 MB** against a 102 MB
cube, and streaming 20,000 steps of the irregular leaf sequence peaks at
**0.04 MB**; both are measured with `tracemalloc` in the suite rather than
argued for. No file is written by any of this, and `verify.py` checks both —
the package directory is unchanged after the run, and the module names none of
`open`, `write`, `savetxt`, `save`, `tofile` anywhere in its syntax tree.

### One clock, moving: evolution, handoffs, and the dwell clock

`dynamics.py`, `commitment.py` and `raw_experiments.py` are where a phase first
moves. Everything before them supplies numbers; these three turn numbers into a
trajectory, a trajectory into a lock history, and a lock history into a record.

#### The step

At every elementary interval, with the drift frozen at the interval's start:

```
phase_next = phase + (Δ − K(t_start)·sin(phase))·h + kick
```

That is Euler–Maruyama, and it is written by calling `model._step_euler` rather
than by re-deriving it, so **at zero diffusion — where every kick is exactly
`0.0` — the walk is bit-for-bit deterministic Euler**, not close to it. The
suite asserts equality rather than a tolerance, and separately measures the
convergence order against a fine RK4 reference: 0.997, 0.999 at timesteps
0.02/0.01/0.005, which is first order as Euler–Maruyama must be.

The intervals are the **elementary mesh** — the finest uniform grid *union* the
exact tongue crossings — not the uniform grid. An entry or exit inside a step
subdivides it, the Brownian tree divides that step's already-drawn kick
conditionally between the halves, and the walk sees two intervals where the grid
has one. That is what "no dwell may begin before exact eligibility" means
operationally: an eligibility change always happens at an interval boundary, so
no interval is ever half-eligible.

An interval's right endpoint is *supplied* to the stepper rather than recomputed
as `start_time + duration`. Those are different floats — a `TreeLeaf` carries a
duration, and `leaf.start_time + leaf.duration` need not be the boundary the
tree was split at — and every eligibility test here is an exact equality against
a crossing time. A recomputed endpoint would put the last step of a window a
rounding away from its exit and send it to `stable_phase` instead of
`boundary_phase`. The two are cross-checked bit for bit, so a caller that
supplies the wrong endpoint is refused rather than integrated.

#### One authoritative phase

Each clock has exactly one state variable: an **unwrapped absolute phase on the
real line**. It is never wrapped and never has its winding discarded.
Periodicity enters in three places and nowhere else — `sin` in the drift, `cos`
in the contraction, and the circular distance used to measure lock error — and
each is a periodic *function of* the authoritative phase, not a substitute for
it.

#### Three eligibility states, and what happens at the boundary

Eligibility is strict, `K > |Δ|`, so a time is in one of three states and the
package names all three:

| State | Where | Target | Dwell |
| --- | --- | --- | --- |
| `outside` | before entry, after exit | none (NaN) | none |
| `boundary` | **exactly** at a crossing, `K = \|Δ\|` | the finite `boundary_phase(Δ) = sign(Δ)·π/2` | **none** |
| `interior` | strictly inside | `stable_phase(K(t), Δ)` | possible |

`stable_phase` returns NaN at equality by design — the attracting and repelling
fixed points have merged there — so it is called *only* in the strict interior,
and the equality state uses the finite limit instead. A mutation that calls
`stable_phase` at equality fails the suite rather than seeding a NaN lift that
would surface many steps later as a clock that had quietly stopped existing.

Equality is a **zero-duration coordinate handoff**: no time passes, no noise is
consumed (the tree emits a leaf of zero duration and exactly zero kick), the
phase is returned bit-identical, and **no dwell is earned**. The last is not a
special case — the endpoint is not strictly eligible, and the merged fixed point
has zero contraction anyway, which the strict guard rejects on its own.

Eligibility is decided by comparing *times* against the exact calibrated
crossings, not by re-testing `K > |Δ|` at a sample. The two spellings disagree
within a rounding of a crossing, and the one that can disagree with the noise
tree — which was split *at* the crossing time — is the one that would put a
step's target on the wrong side of its own subdivision.

#### The lifted target, and why winding is data

While the width-only control is eligible it carries an auxiliary pair: a
continuously **lifted** target and an **unwrapped** error, with

```
phase = target_lift + error
```

exact at every endpoint, because an eligible step *assigns* the phase as that
sum. At entry the finite `boundary_phase(Δ)` is lifted onto the branch nearest
the incoming phase, using the half-open `[-π, π)` convention so that an exact
half-turn tie resolves deterministically to `-π`; the phase itself is left
untouched. While eligible the target follows its own lifted branch and the error
is never wrapped. At exit the pair is discarded and the phase kept. At re-entry
a new lift is built from that unchanged phase.

Without this, a clock sitting just across the `±π` cut gains or loses a whole
turn at a tongue entry or exit **even with the noise switched off** — a run that
teleports by `2π` at an eligibility crossing is not the run anybody meant to
integrate. The suite runs entry/exit/re-entry cycles on a two-pulse train and
requires that no step, handoff included, move the phase by more than 0.5 rad.

The error and the lock error are two different quantities: the error keeps its
winding and may pass a half turn, while the band is measured against its
circular image in `[-π, π)`. A clock outside its tongue slips whole turns and
the record reports them.

#### The lock band: proximity **and** contraction

The lock band is proximity **and** strict local contraction, so the repelling
Adler branch never earns dwell. Both comparisons are strict:

```
inside_now  ⟺  |lock_error| < tolerance  and  contraction < 0
```

Why the second half is not optional: inside the tongue the attracting and
repelling phases approach one another and merge at the edge, so any fixed
tolerance eventually contains both. At `K = 1`, `Δ = 0.999` they are **0.08945
rad apart**, and a 0.1 rad band declares the *repelling* fixed point locked — in
the region where the experiment expects critical slowing to suppress
commitment, which is to say it manufactures events exactly where the hypothesis
predicts none. No noise is needed for the mistake.

Strictness at both boundaries is load-bearing for the same reason: the merged
fixed point at the tongue edge has contraction exactly zero, and a `<=` would
admit it.

The full model supplies `-K·cos(phase)`; the width-only control supplies
`-rate·cos(error)`. Same logical rule, so the control cannot quietly degenerate
into a proximity-only band.

#### The dwell clock is two timestamps

**Dwell is a pair of timestamps, never a count of samples.**

```
first qualifying endpoint  →  inside_since = t          (no retroactive credit)
commit                     →  first later t with t − inside_since ≥ dwell
anything that breaks lock  →  inside_since cleared outright
```

The rejected alternative — add one timestep per inside sample — commits from a
*single* observation when `dwell = dt`, with no elapsed time inside the band at
all, and silently substitutes `ceil(dwell/dt)·dt` for the declared dwell
whenever the ratio is not an integer, so every refinement level would be testing
a different physical rule. The timestamp rule is one subtraction and is
identical for `dwell < dt`, `dwell = dt` and a non-integral ratio; the suite
checks six ratios and predicts each commit index from `ceil(dwell/dt)` rather
than reading it off the implementation.

Three things reset it, all completely rather than partially: the phase leaves
the band, the contraction is lost, or the falling tongue edge removes
eligibility. One exit costs a full fresh dwell, not a step.

Event times are interval-censored to one elementary interval.

#### The two mechanism-negative controls

**Central-clock control.** One clock at exactly zero mismatch, where tongue
widening is absent. It is the full Adler model, not a third set of dynamics —
giving it its own label would be the beginning of it not being a control.

**Fixed-contraction width-only control.** **The width-only control is the full
calibrated Adler model whenever the clock is ineligible; only its eligible
contraction is replaced, by one manifest-frozen rate that never scales with
coupling.**

That first clause is the causal contract, and an earlier draft of the plan broke
it by giving the control free mismatch drift outside the tongue. A marginal
clock on the rising edge feels the subthreshold Adler pull, so removing it
changes the phase *delivered to* tongue entry, and a paired full-minus-control
difference would then confound "we removed interior relaxation" with "we changed
the entry distribution". The two processes here share the initial phase and the
Brownian tree and are the **same run**, bit for bit, until the first eligible
interval.

The rate is frozen by a *rule* rather than chosen: it is the central clock's
Adler relaxation rate at the geometric midpoint of a declared reference
peak-coupling range. For a zero-mismatch clock the stable phase is `arcsin(0) =
0`, so that rate is `K_ref·cos(0) = K_ref = sqrt(lower·upper)` — written out in
`fixed_contraction_rate` in that form rather than collapsed to the square root,
because the collapse is a coincidence of the central clock. Geometric rather
than arithmetic midpoint because the sweep is logarithmic and the estimand is a
log-log exponent.

`REFERENCE_PEAK_COUPLING_RANGE` in `dynamics.py` is **provisional**. The
production range is chosen by a later ticket's labelled range-only pilot under
frozen count rules, and the number written into the production manifest is
`fixed_contraction_rate(*that range)`. What is frozen *now*, before any
stochastic behaviour has been looked at, is the rule: the pilot may choose the
range and may not choose the rate.

#### One Brownian tree, and no coupling anywhere near it

The full model and the control consume the *same* elementary leaves — the same
kicks, durations and node paths — because the split key excludes the model
label. That is what makes the pairing real rather than nominal.

Nothing that could smuggle a coupling dependence into the event generator has
one. `raw_noise_stream` and `raw_lock_criterion` have no coupling, peak or
mismatch argument; the tolerance and the dwell are one pair of numbers for the
whole sweep; the contraction rate is one number from a declared reference range;
and across peak couplings 0.5, 1 and 2 the elementary kicks on every step no
crossing touches are identical. Different couplings *do* move the crossing
times, so the leaves a crossing subdivides differ between runs by construction —
that is the tongue geometry, and the suite compares the uniform steps rather
than pretending otherwise.

#### An example record

`raw_experiments.run_one_clock` returns an `OneClockRecord`. A full-Adler clock
at `Δ = 0.4` under a unit-peak pulse of duration 4, `D = 0.02`, band 0.25 rad,
dwell 0.30, timestep 0.01:

```
label='S4/interior/intermediate/full'  model='full'  drive='pulsed'
detuning=0.4  peak_coupling=1.0  phase_diffusion=0.02  fixed_rate=0.0
initial_phase=0.7  trial=0  clock=0  timestep=0.01
tolerance=0.25  dwell=0.3
category='committed'  committed_at=0.93  band_entries=3  dwell_resets=2
endpoints=294  eligible_endpoints=206  boundary_endpoints=1
final_phase=0.667  winding=-0.00526  lift_residual=0.0
```

That clock entered its band three times and lost the dwell twice before the
third stay lasted 0.30 in a row. The tongue window is
`(−1.128, 1.128)`, so `committed_at = 0.93` is on the *falling* side: noise
broke two earlier attempts and the third survived. Reading anything more into
one path than that would be reading a single sample as a statistic.

`band_entries` and `dwell_resets` are kept apart so that noise-assisted band
entry and noise-broken dwell stay distinguishable; non-monotonic behaviour in
noise strength — weak noise helping a clock in, strong noise breaking the dwell
— is expected and is reported rather than smoothed. `boundary_endpoints` counts
endpoints sitting exactly on a crossing; it is normally small and normally not
zero.

The record schema contains no prediction, weight, probability or hazard, and the
suite fails if it ever does.

Every cross-field relation a record can express is *enforced* rather than merely
observed, and the bar is **reachability**, not internal consistency — a ledger
can be individually valid in every field and still be a history this machine
cannot produce. There is one clock, so there is at most one live dwell, so:

```
band_entries == dwell_resets + int(a dwell is live)      # in DwellState
band_entries − dwell_resets ∈ {0, 1}                     # in the record
ever_inside  == (band_entries > 0)
```

A `committed` record ends on the dwell that completed, so its difference is
exactly one; `never_eligible` and `lock_failed` never reached the band, so both
counts are zero; a state with no last observation is the initial state and
carries no history at all. On top of that: `band_entries ≤ eligible_endpoints`,
the interior and boundary subtotals cannot exceed the endpoint count, the
commitment time agrees with the category, a **stationary** record has exactly
zero boundary endpoints because `OpenSchedule` has no boundary to sit on, and
the full model carries exactly `0.0` while the width-only control carries
exactly `FIXED_CONTRACTION_RATE`.

#### One rate, at three boundaries

The width-only control is defined by having *one* contraction rate. Closing the
raw factory to per-run rates was necessary and not sufficient: `ClockPath` is
exported and `run_one_clock` accepts it, so an independent review ran three
different width-only trajectories under one ordinary
`S4/stationary/interior/zero/width_only` label at rates `0.1`,
`1.0000000000000002` and `7.0`, each producing an ordinary record. The control's
causal claim is that *only* coupling-dependent interior contraction was removed;
a rate that varies per run puts the factor back under the control's own name.

The requirement is enforced at all three public boundaries a bypass would have
to cross — **construction** (`ClockPath`), **execution** (`run_one_clock`, which
re-checks rather than trusting the path, because a frozen dataclass is not
frozen against `object.__setattr__`) and **recording** (`OneClockRecord`).

No sensitivity API replaces it. Ticket 03 does not need one, and an unlabelled
escape hatch would be the same bypass with a longer signature; a real one needs
its own type and a record identity that survives into the ledger, against a
frozen production manifest. The one exported name carrying the word "rate",
`fixed_contraction_rate`, takes a coupling *range* and derives the scalar — it
cannot be handed one.

#### A trusted snapshot at the execution boundary

`ClockPath` derives its tongue schedule once, in the constructor, and caches
it — and that schedule is load-bearing in a way no other field is. It owns
where the Brownian tree splits a step, where the entry and exit handoffs
happen, which endpoints are strictly eligible, which endpoint target is
`boundary_phase` rather than `stable_phase`, and therefore where dwell can start
and where it must reset.

A frozen dataclass is not frozen against `object.__setattr__`, and re-checking
only the contraction rate left the rest of that geometry stale:

```python
path = ClockPath(0.9, train, "full", 0.0)
object.__setattr__(path, "detuning", 0.1)
run_one_clock(path, ...)          # an ordinary committed record
```

The record said `detuning = 0.1`; the run used the `0.9` clock's window
`(1.5903, 2.4097)` instead of `(0.4097, 3.5903)`. Between `0.4097` and `1.5903`
the real clock is strictly eligible and the accepted path treated it as outside,
so the transition branch, the target, the contraction and the dwell geometry
were all the wrong ones.

`run_one_clock` and `synthetic_dwell_history` now take a **trusted snapshot** of
their stateful inputs
before anything else happens — `dynamics.validate_clock_path` and
`commitment.validate_lock_criterion`, the execution-boundary companions of
`raw_runner.validate_raw_config`. The algorithm is: reconstruct by **exact type
and named field, one level at a time**, so a mutation cannot move down a level;
rebuild every `RaisedCosinePulse` from its own `peak`, `duration` and `center`,
then the `PulseTrain` from those pulses, then the `ClockPath` from the detuning,
the rebuilt drive, the model and the rate; **derive a fresh schedule** and
require it to equal the cached one by type and by value; then run the rebuilt
object or reject. Fifteen seams are covered — the top-level fields, the drive as
a whole, a drive swapped for the other kind, each of a nested pulse's own fields,
a malformed nested member, both schedule types planted wholesale, and the
stationary drive's own peak.

The criterion is rebuilt at **both** of the state machine's public execution
boundaries, not one. `synthetic_dwell_history` is the other door, and it used to
return a *committed* history from a criterion edited to `tolerance = pi` — a band
that admits every phase — on a lock error of 3.0 rad. A contract enforced at one
door and not the other is not a contract.

#### Validated wrappers over private trusted kernels

The two helper boundaries were not the only public doors. `LockCriterion.advance`
*is* the state machine and `ClockPath`'s stepping methods *are* the stepper; both
classes are exported, and both accepted objects the helpers refused:

```python
criterion = LockCriterion(0.25, 0.3)
object.__setattr__(criterion, "tolerance", math.pi)
criterion.advance(DwellState(), 0.0, True, 3.0, -1.0)
criterion.advance(...,          0.3, True, 3.0, -1.0)
# -> committed at t = 0.3, from a 3.0 rad lock error
```

`3.0 < pi`, so the edited object admitted an error outside the largest band the
class can express. On the stepper side a width-only path whose detuning had been
edited from `0.9` to `0.1` kept the `0.9` tongue window, and a direct public walk
reported `outside` with a NaN target where the real clock is strictly
`interior` — not a reporting defect, because the stale branch takes the full
Adler transition where the fixed-contraction one was due, so **the phase update
itself differs**.

The repair is the private-kernel convention this package already uses for
`model._wrap_phase`: **public wrappers validate, private kernels do not.**

| Public door | Delegates to |
| --- | --- |
| `LockCriterion.inside_now` | `_inside_now_trusted` |
| `LockCriterion.advance` | `_advance_trusted` (module-level) |
| `ClockPath.schedule` | `_cached_schedule` |
| `ClockPath.crossings` | `_crossings_trusted` |
| `ClockPath.initial_state` | `_initial_state_trusted` |
| `ClockPath.entry_handoff` / `exit_handoff` | `_entry_handoff_trusted` / `_exit_handoff_trusted` |
| `ClockPath.advance` | `_advance_trusted` |
| `ClockPath.observe` | `_observe_trusted` |

Each wrapper rebuilds **every** frozen input — the criterion, the path, and any
supplied `DwellState`, `PathState` or `TreeLeaf` — and then calls the kernel.
Supplied values are rebuilt too, through their own classes' constructors: a
`PathState` given a `target_lift` after construction kept its NaN `error`, and
`observe` read the half-live pair as a live target.

The leaf matters as much as the rest, and less obviously. Four of its six fields
were refused only *incidentally* — a non-finite kick surfaced later as a
non-finite phase, an edited `duration` or `start_time` broke the span and start
cross-checks — which is a failure partway through a transition rather than
before it. The other two were not refused at all: `step_index` and `node_path`
are the leaf's **address**, saying which uniform step this kick belongs to and
which node of the Brownian tree carved it out, and a leaf with a rewritten
address stepped the phase without complaint. `TreeLeaf` is Ticket 02's type, so
the rebuild goes through *its* constructor; this module restates none of those
rules.

Rebuilding a `PathState` through its own constructor is necessary and **not
sufficient**, and that is the subtlest part of the boundary. That constructor
knows one state and no path, so it cannot know whose target the lift is a lift
*of*. On a stationary width-only clock at mismatch `0.35`:

```python
clean  = path.initial_state(0.0, 0.2)
forged = PathState(clean.time, clean.phase,
                   clean.target_lift + 1.0, clean.error - 1.0)
```

`forged` is constructor-valid and its reconstruction still holds to a rounding —
`(lift + error) − phase` is `−5.6e−17` — yet its lift is a whole radian from the
target this path actually has. The compensating error changes `sin(error)`, so
the width-only transition changes: the next authoritative phase came out `−0.708`
where the clean state gives `0.216`.

So supplied states are validated **against their owning path**, before any
kernel:

| Rule | What it refuses |
| --- | --- |
| Reconstruction identity within `conservation_tolerance` | `error += 1`, which carried `lift_residual == 1.0` |
| Pair ownership: only `width_only` has one | a full-model state observed with that pair as its target |
| Lifecycle: required strictly inside, forbidden outside, either at a crossing | an interior control state observed through the full-Adler branch |
| Congruence of the lift to `stable_phase` / `boundary_phase` within `LIFT_RESOLUTION` | the forged pair above |

The **unwrapped error keeps its semantics**: nothing bounds it, and only the
*lift* must be congruent. A five-turn error is legal and its lock error is still
the circular image in `[−π, π)`.

The handoffs are **directed**. `_require_crossing` used to test only that an
instant was *a* boundary, so with a window `(0.41, 3.59)` an entry handoff was
accepted at the exit and an exit handoff at the entry — each performing the
opposite coordinate operation at the wrong edge, each preserving the phase, so
nothing downstream complained. An entry must now land on a window's first
endpoint and an exit on its second, checked in the wrapper so a wrong-edge call
refuses before the kernel.

**One declared limitation.** Rebuilding a `TreeLeaf` proves every field lies in
Ticket 02's constructor domain — that the address is *canonical* — and no more.
A `ClockPath` has no mesh or stream context, so it cannot authenticate a
different yet syntactically valid `step_index` or `node_path` as the one this
kick was actually drawn under. Provenance is the high-level runner's guarantee:
`run_one_clock` accepts leaves only from its own validated `PhaseNoiseStream`,
which is where a leaf's address and its value are bound together. Closing it at
the primitive would mean giving the stepper a stream, which is a different
ownership boundary than this package draws.

The kernels are underscore-prefixed, absent from every `__all__`, and documented
as accepting only rebuilt objects. `run_one_clock` and `synthetic_dwell_history`
validate once and call the kernels in their loops, so nothing rebuilds a pulse
train per leaf. The suite asserts that clean walks through the public doors are
bit-identical to the kernels and to the high-level runs, across pulsed and
stationary × full and width-only.

Every refusal happens **before any key is drawn**, which the suite asserts by
counting calls into `_normals_from_key`: a rejected run cannot have changed the
realization it was rejected from. And the rebuild is exact for an unmutated
path — identical floats give an identical schedule and an identical manifest
digest, therefore identical noise keys — so clean pulsed and stationary,
full and width-only runs remain bit-identical, checked against an independent
walk that does not go through the boundary.

#### Stationary and pulsed

Experiment S4 runs one clock under **stationary and pulsed** coupling, and the
two are different objects rather than a parameter choice. A `StationaryDrive`
has no rising edge, no falling edge, no finite eligibility window and — the part
that matters for the machinery — **no tongue crossings at all**. Its schedule is
an `OpenSchedule`: strictly eligible at every instant, with nothing to enter,
nothing to leave, and no equality instant where `stable_phase` is undefined. A
clock at or above the stationary threshold gets the empty `TongueSchedule`
instead, because eligibility is strict.

Faking a stationary run as a very wide pulse would put two fictitious crossing
times into the Brownian tree's key schedule, split two steps for an eligibility
change that never happens, and describe edges the run does not have. The suite
contains that mutation and requires it to fail.

Every record carries `drive='pulsed'` or `drive='stationary'` as a field rather
than a convention inside its label, so the two halves cannot be aggregated by
accident: a stationary clock's band entries and resets are not comparable with a
pulsed clock's.

What this is **not**: no killed-diffusion oracle, no stationary-hazard
experiment, no burn-in or left-truncation semantics, no exponential-survival
test. Those are experiments S3 and S6 and remain later tickets.

#### Two declared numerical domains

Both are boundaries this package refuses to cross rather than approximations it
makes quietly.

**Angular resolution: `|phase| < 2097152`.** A lift is only a lift if it is
congruent to its principal target modulo a turn, and on binary64 that stops
being true long before the phase stops being finite. At `phase = 1e16` and
target `π/2`, `nearest_lift` returns a value **0.504 rad** from congruent while
`lift + error == phase` still holds exactly — the identity passes and the
auxiliary target has become the wrong target. There is no representation trick
that fixes this inside one float, so the domain is declared: `LIFT_RESOLUTION`
is a nanoradian, and the accepted domain is the largest binade where four ulps
stay inside it. That is about 333,772 turns of accumulated winding; the runs
here wind single digits. The refusal happens at the public boundary —
`nearest_lift`, `continuous_lift`, `PathState`, `ClockPath.initial_state` — and
therefore **before any noise is drawn**, which the suite asserts by counting
keys.

**Dwell equality: a scale-aware few-ulp allowance.** The dwell rule is a
subtraction of two timestamps compared against a duration, and binary64 can
round an exact physical equality the wrong way: with `inside_since = 0.1` and
`moment = 0.3`, `0.3 − 0.1` is `0.19999999999999998` and falls short of a `0.2`
dwell, so a dwell that physically completed at `0.3` was reported at `0.4`. The
comparison now forgives `elapsed_allowance(...)` — four ulps at the timeline's
own scale — because the error is an error *of the timeline*, not of the
difference. It cannot become slack: the allowance is required in code to stay
under one part in a million of the dwell, and a timeline that cannot meet that
(timestamps near `1e12` against a `1e-3` dwell, where one ulp is `1e-4`) raises
rather than committing on the last bit of a subtraction.

`DwellState` also carries the last endpoint it saw, so `LockCriterion.advance`
refuses a duplicate or a backwards observation itself rather than relying on the
sequence helper, and a state that committed before its own `inside_since` — or
at any instant other than its last observation — cannot be constructed.

#### API changes in the review fix-up

| Change | Why |
| --- | --- |
| `raw_one_clock_path(config, detuning, model, **drive**)` — the `fixed_rate` argument is **removed** | The same configuration could be handed `0.1`, `1.0` and `10.0` and returned three different width-only controls, which is the per-run tuning the control exists to exclude. The rate is now derived in one place, `dynamics.FIXED_CONTRACTION_RATE`. No sensitivity factory replaces it: an unlabelled one would be the same bypass under a different name. |
| `OneClockRecord` gains a `drive` field | Stationary and pulsed records must not be aggregated by accident. |
| `DwellState` gains a `last_time` field (appended, so positional construction is unchanged) | Duplicate and backwards endpoints are refused at the public transition. |
| `ClockPath.train` accepts `PulseTrain` **or** `StationaryDrive`, by exact type | The stationary half of S4. |
| `RawEventConfig` rejects `lock_tolerance == pi` | It used to accept a configuration that `raw_lock_criterion` then refused. |
| `nearest_lift`, `continuous_lift`, `PathState`, `ClockPath.initial_state` reject `|phase| >= 2097152` | See the angular-resolution domain above. |
| **`ClockPath(…, "width_only", rate)` requires `rate == FIXED_CONTRACTION_RATE` exactly**; `run_one_clock` and `OneClockRecord` re-check it | The exported constructor was the remaining public route to a per-run rate. Pass `dynamics.FIXED_CONTRACTION_RATE`. |
| `DwellState` and `OneClockRecord` enforce entry/reset conservation and the `ever_inside` latch | Individually valid fields were still assembling histories the state machine cannot reach. |
| `OneClockRecord` requires `boundary_endpoints == 0` when `drive == "stationary"` | `OpenSchedule` has no boundary. |
| **Supplied `PathState` values are validated against their owning path**, and the handoffs are edge-specific | A constructor-valid state can still be wrong for its path; a handoff at the opposite edge performs the opposite operation. Clean callers unaffected. |
| **Every public `ClockPath` and `LockCriterion` transition rebuilds its frozen inputs and delegates to a private trusted kernel** | The exported primitives accepted stale objects the helper boundaries refused. Clean callers are unaffected and bit-identical. |
| **`run_one_clock` rebuilds its `ClockPath` and `LockCriterion`, and `synthetic_dwell_history` rebuilds its `LockCriterion`; both run only the rebuilt copies** (`dynamics.validate_clock_path`, `commitment.validate_lock_criterion`, both newly exported) | A post-construction edit combined new fields with old cached tongue geometry. Clean runs are unaffected and bit-identical. |

One numerical value did move. `fixed_contraction_rate` no longer forms
`lower · upper`, so the provisional rate is now `sqrt(0.5) · sqrt(2.0) =
1.0000000000000002` rather than `sqrt(0.5 · 2.0) = 1.0` — one ulp, traded for a
derivation that is finite and positive across the whole binary64 domain instead
of returning NaN at `1e-300` and raising on an infinity at `1e200`. The rate is
explicitly provisional, no production run has been recorded against it, and
nothing needs migrating.

The v3 noise key schema, every keyed value, and all 56 deterministic and
Brownian-tree residuals are unchanged by all of this, as are the nine ticket-03
check residuals.

### The independent-oracle rule

A test that takes both sides of its comparison from the code under test proves
only that the code agrees with itself. An independent review found two places
where this package did exactly that, and both are now closed by oracles written
out inside `verify.py` and never imported:

- **The pulse.** Pinning endpoints, peak, symmetry, bounds and area does *not*
  pin a raised cosine — `E(u) = (1 − u²)^p` with `p = 2.38175…` satisfies every
  one of those and is a different pulse, differing by 4e-3 at quarter support
  and giving different tongue crossings. The verifier now compares interior
  samples against a transcribed `0.5·(1 + cos(2πu))`, and checks crossings
  against both an independently derived identity,
  `T·acos(√(|Δ|/K_peak))/π`, and a bisection that uses no closed form at all.
  That counterexample envelope is kept in the suite as a live regression: the
  checks assert that it *fails* them.
- **The normalization.** Measuring the flux against `analytic.SEMICIRCLE_CONSTANT`
  cannot detect a common factor error — doubling the sum and changing the
  constant from `π/2` to `π` left the old checks passing with identical
  residuals. The verifier now owns `math.pi / 2`, asserts the exported constant
  equals it exactly, and compares the flux against `0.5·π·K²` computed locally.
  It also demonstrates it can see the error: the same comparison applied to a
  deliberately doubled flux is required to deviate by ~1.0, not by ~tolerance.

### Input validation

Invalid numbers are rejected where they enter, with `TypeError` for the wrong
kind of thing and `ValueError` for an impossible value. This is not decoration:
a NaN peak coupling used to make `K > |Δ|` false everywhere, so a run would
report that no clock was ever admitted to the tongue — a clean, physical-looking,
entirely fictitious result. A negative `dt` used to be silently collapsed to a
single forward step and reported as a positive step that was never taken.

The contract covers the **whole exported surface**, not a sample of it: every
callable and dataclass in the seven modules' `__all__` lists, *and every public
method and property of every exported class* — evaluation arguments
(`envelope(t)`, `coupling(t)`), the eligibility and fixed-point predicates, the
drift, the public steppers including their **callback's return value**, the
analytic helpers, and trajectory queries such as `index_at(t)`.

Booleans and strings are rejected as deliberately as NaN — `bool` is a subclass
of `int` and `float("1")` succeeds, so both would otherwise become silent
physical parameters. The frozen dataclasses write their validated floats back,
so an instance can never retain a value that merely converted once in the
constructor.

Three domain rules are worth stating, because each closes a specific *silent*
failure rather than a loud one:

- **Couplings must be non-negative.** A negative coupling is not a weak one; it
  moves the stable fixed point by π. Unchecked it propagates through `K > |Δ|`
  as a plain `False`, through `√(K² − Δ²)` as NaN, and through anything that
  *squares* it as a positive, physical-looking number —
  `efolds_available([0,1], [-1,-1], 0.5)` used to return a perfectly plausible
  e-fold budget of 0.866.
- **Log-fit inputs are rejected, not filtered.** `fit_log_log_exponent([1,-2,3],
  [1,4,9])` used to drop the middle point and report ≈2.0 for a sweep the caller
  believed had three points. Silent post-selection inside a fitter is worse than
  an error.
- **Eligibility flags must be a genuine boolean array.** `['False','True',
  'False']` coerces to *all true*, and `[0,2,0]` to `[False,True,False]`. Both
  are plausible; neither is what the caller meant. Coercion is refused.

`verify.py` drives that surface from a table of **530 public callables** and
**1380 invalid calls**, then asserts completeness in two senses: over the public
surface *including methods and properties*, and over **every parameter** —
**1274 of them** — of those callables: each one must carry at least one declared
probe or a visible exemption. Adding an unprobed parameter fails the suite by
name. These three numbers, and the check count below, are **pinned**: a check
reads them back out of this file and compares them with what the table
actually holds, because they had gone stale once — the README still said 70,
245 and 141 two tickets after the table had grown past them, and nothing could
tell, since no check read the prose. The raw-boundary modules are in that table on the same terms as the rest,
subclass probes included.

The same check also **reads the prose for schema claims that have stopped being
true**. Three sentences across this file and `killed_diffusion.py` still
described the live validation report as an earlier version several rounds after
it had moved on, and each was found by a human reading rather than by a check.
The first attempt at mechanizing that was a list of forbidden literal strings,
which caught the fully-qualified schema name and missed every abbreviated form
— the three sentences named the version without naming the schema, and all
three passed. It is now a small set of regular expressions matched against
this README, `killed_diffusion.py` and `moving_band_audit.py`, and they match
the *construction* rather than the string: a present-tense assertion that the
live report, schema or tag is an earlier version.

That distinction is the whole point. This package's documentation explains why
each schema version was superseded, so prose about earlier versions is not only
legal but load-bearing; a check that banned mentioning them would delete the
history it should be protecting. Historical wording — a version that *carried*
something, a version that *changed* something, a consumer written against an
older tag that must refuse this shape — is past-tense or third-party and does
not match. No marker comment is required to make a migration note legal, since
a marker is one more thing to remember to write. Ten fixtures pin both
directions: the three original stale sentences and two constructed variants
must be caught, and five pieces of genuine historical prose must pass. Without
those the patterns could rot into matching nothing, and the check would pass
more loudly every round while catching less.

### What that completeness claim does and does not mean

It is coverage of **parameters**, not of **value classes**. Each parameter is
required to have at least one declared probe; nothing derives the set of ways a
parameter could be invalid. A parameter probed only for shape can still be
unprobed for sign.

That is not hypothetical, and the honest version of this section has to say so.
`Trajectory.coupling` was probed for length mismatch and not for sign, so a
same-shape **negative** coupling history was accepted — and the trajectory then
reported `eligible = False` and `stable = NaN` at every sample. That is a clean
"never admitted, never locked" reading, which is precisely the physical result
this package exists to measure, manufactured from input the model does not
define. The probe exists now, but it exists because an independent reviewer
found the gap, not because the mechanism demanded it.

So: read the guarantee as *every parameter has a declared probe*, never as
*every invalid value is rejected*. Name-level completeness is not parameter
completeness, and parameter completeness is not value-class completeness; the
suite enforces the first two and the second was itself only reached after a
table that was exhaustive over `__all__` left all three domain rules above
unenforced.

There is exactly **one** documented exception, and it is narrow on purpose.
`measured_relaxation_rate`'s `deviation` argument may contain `NaN`, because
that is the sentinel `Trajectory.tracking_error()` emits wherever a clock was
ineligible — genuinely missing data that a fit over the eligible band must skip.
It may **not** contain infinities: an infinity is a numerical blow-up, and
quietly filtering it would hide the very failure it signals. Booleans, strings
and complex arrays are refused there too. The carve-out is for NaN-as-missing,
and for nothing else.

Validation costs a scan of its input, which is unaffordable four times per
Runge-Kutta step, so each validated public function is a thin boundary over an
unvalidated private kernel (`_wrap_phase`, `_is_eligible`, `_adler_drift`,
`_raised_cosine`, `_step_rk4`, …). The integrator calls the kernels; everything
reachable from outside the package calls the boundary. The suite runs in about
20 seconds, so the boundary costs nothing that matters.

---

## Running it

```
python3 -m adler_born_two_channel.verify
python3 -m adler_born_two_channel.verify --verbose
python3 adler_born_two_channel/verify.py
```

All three work; the third is the direct-script path, which inserts the repository
root on `sys.path` so it behaves identically. Every check prints a measured
residual and the tolerance it was held to, so a pass is auditable rather than a
bare assertion. `--verbose` adds the full detail line for every check.

The exit code is 0 only if every check passes. To see the failure path end to
end:

```
python3 -m adler_born_two_channel.verify --prove-failure-exit ; echo "exit=$?"
```

This appends one deliberately failing check and exits 1.

---

## Results

**127 checks**, all passing, in a few minutes on one core — itself a
machine-dependent snapshot (Python 3.13.11,
NumPy 2.4.3). Twenty-six of them are the deterministic and analytic calibration
and their residuals are unchanged; six are the raw-process boundary and are
structural, so their residual is a problem count and their tolerance is zero;
eighteen are the keyed noise source and the Brownian tree, also unchanged; nine
are the one-clock stochastic evolution, the exact tongue handoffs, the timestamp
dwell clock, the contraction guard, the width-only control, the
coupling-blindness of every criterion, the S4 records, and a mutation battery
against all of them; six close the public-contract defects an independent
review found — the dwell equality rule, the stationary half of S4, the frozen
contraction rate, the lift-congruence domain, the lock-band domain, and the
ledger invariants; two close the alias bypasses its closure re-review
found — a contract enforced at one public boundary and reachable past another;
one closes the stale-object seam its final re-review found;
and one closes the same seam at the exported *primitive* transitions. The
remaining ten are the continuum audit: the killed-diffusion oracle against its
closed-form limit, the per-position paired endpoint refinement ladder against
the oracle, the frozen contract and sampling design with its numerical
no-result, the bridge expression, the admissible geometry, the reset-only
pathwise-subset replay with its survival and commitment gates, the canonical
disjoint audit address with the frozen reduced matrix reported per regime, a
mutation battery against all of it, the isolation of both new modules from the
raw import graph, and the suite's own bounded resident set. The last nine are
the finite one-channel race and the run on disk: the stop measured against an
independent walk of every clock, the retained same-endpoint co-completers and
the absence of any way to pick one, the four exhaustive categories and their
per-clock reconciliation, the censored exposure and the no-stop shadow sample
that exposes it, the byte-identical files across four live-block sizes, the
measured memory bound, the survival curve reconstructed from raw commit times,
the run-directory gate, a mutation battery against wrong races and forged runs,
and the agreement of the ledger's duplicated vocabularies with the modules they
were copied from. Four more close an independent review's findings: the
authoritative per-`(trial, clock)` table and the aggregate derived from it, the
semantic forgeries that must refuse after every hash is recomputed, the declared
population the race is bound to, and the package-scoped output contract with its
cleanup. A tenth closes the second round: every row's derived observation
interval, the causal ordering of its recorded instants, and the same-key
shadow-prefix relation. An eleventh closes the third: reader and writer I/O
anchored to pinned directory descriptors. The final ten are the separate
read-only analysis layer: the display record that keeps a conditional
probability from being read as a rate, survival and cumulative hazard against
hand-computed censored truth, the three unresolved causes and the retained
co-completers, the complementary-log-log estimator recovering declared exponents
while keeping its zero-event and all-event cells, the refusals — curved,
non-power, saturated, separated and uninformative — each with a named reason,
the paired master-trial bootstrap measured against the wrong independent-cell
one, the four separately typed comparators, the causal defeat rule fired one
defeater at a time, the read-only proof with its before-and-after digests, and
the isolation of all three modules from the raw graph together with the
manifest-bound numerical gate. Two more close an independent review's findings:
master-trial identity binding the stream provenance so two datasets numbered
alike cannot be paired, and control roles and holdout evidence bound to the runs
they claim. A thirteenth closes its re-review: general lack of fit decided by a
nuisance-adjusted cluster-robust score test with a master-trial multiplier
bootstrap, whose size and power are measured over frozen ensembles at four
dependence strengths. A fourteenth closes its own re-review: every
authoritative factory computes from the object it re-read from disk, so a
handle carrying forged rows cannot reach a fit. The last nine are the
feasibility, power and freeze layer: the versioned experiments module and its
closed vocabularies, the measured benchmark of the real streamed race with its
two controlled quadruplings, the power calculation checked against a simulation
of itself, the production numerical budget and the empty trial window it
reports, the finite matrix with its parity, staircase fallback and scaled cost,
the pilot firewall attacked at every forbidden field, the count-only range rule
on five frozen count patterns, the dry-run proposal with every one of its
declared fields omitted in turn, and the checkpoint-1 handoff that derives its
own feasibility verdict. **There is no ticket-08 production gate in this
package.** `production.py` does not exist, no registered check covers a
production gate, the acceptance criteria stop at 121, and ticket 08 has not
begun. Highlights, with the numbers
the suite reports:

**The deterministic limit.** At zero diffusion the Euler–Maruyama walk
reproduces `model.step_euler` **bit for bit** over the whole pulse, at negative,
zero and positive mismatch — an identity, not a tolerance. Against a 1e−5 RK4
reference the measured convergence orders are **0.996 and 1.001** at timesteps
0.02/0.01/0.005. With the peak coupling set to zero, no clock is ever eligible
and none commits at dwells shorter than, equal to, longer than or non-integral
in the timestep.

**The coordinate handoffs.** Thirty entry handoffs across three mismatches and
ten phases — both sides of the ±π cut, seven turns above the principal branch,
five below, and the exact half-turn tie — each return the phase and its sine bit
for bit and put the error in `[-π, π)`; the paired exit handoffs do the same.
The lift identity `phase = target_lift + error` holds to **7.1e−15** against a
4-ulp bound at that scale. Eight entry/exit/re-entry cycles on a two-pulse
train — both mismatch signs, both models, with and without diffusion — move the
phase by at most **0.149 rad** in any one step, so no handoff gained or
lost a turn, and a slipping clock's **3.84 turns** are retained rather than
wrapped away.

**The dwell edges.** Six dwell/timestep ratios — 0.5, 1, 1.5, 2, 2.5 and 3.7 —
all commit at `ceil(dwell/dt)` steps after the first inside sample. One exit
costs a full fresh dwell rather than a step; an exit on the prospective commit
step commits nothing; falling-edge eligibility loss clears the dwell with the
phase still in the band; and two phases 2e−6 rad apart across the ±π cut commit
under the circular distance while a plain subtraction never does.

**The contraction guard.** At `K = 1`, `|Δ| = 0.999` the stable and repelling
phases are **0.08945 rad** apart, inside a 0.1 rad band, and the guard admits
the stable point and the contracting side while refusing the repelling point and
the expanding side, symmetrically in the sign of the mismatch. As `|Δ|/K → 1`
the contraction at the stable phase falls monotonically through 4.36e−1,
1.41e−1, 4.47e−2, 1.41e−2 — critical slowing. At equality itself the two fixed
points have merged, `stable_phase` is NaN, and what refuses the merged point is
that equality is not eligible at all.

**The pairing.** The full model and the width-only control consume the same
**404** elementary leaves from one Brownian tree; their phases agree bit for bit
at all **88** endpoints before the exact tongue entry at `t = −1.1282` and first
differ at `t = −1.1200`; **178** ineligible steps give identical phases for
identical state and noise.

**Free diffusion.** With no coupling and no mismatch the unwrapped phase
variance over 4000 trials is within **1.8%** of the predicted `2·D·t`.

**The dwell equality rule.** `LockCriterion(0.2, 0.2)` with inside samples at
`0.1, 0.2, 0.3, 0.4` commits at **0.3**, not 0.4 — `0.3 − 0.1` is
`0.19999999999999998` in binary64, and the elapsed-time comparison forgives
**2.220e−16** of rounding at that scale. Four other nonzero decimal origins and
a nonuniform mesh agree. Runs reaching 0.9 of a dwell, or falling 1e−9 short of
one, still do not commit.

**The stationary half of S4.** **24** stationary records beside the 24 pulsed
ones, at constant peak coupling 1.0 over 400 steps of 0.01; **19** committed.
Every stationary clock's crossing list is empty and its boundary-endpoint count
is zero.

**The frozen rate.** Across 3 peak couplings × 3 mismatches × 2 drives, with the
tolerance, dwell and diffusion varied underneath, every width-only path carries
the single derived rate **1.0000000000000002** and every full-model path exactly
**0.0**. The
geometric midpoint is finite and positive across seven extremes from
`(5e−324, 5e−324)` to `(DBL_MAX, DBL_MAX)`.

**The lift domain.** `|phase| < 2097152`: four ulps there is **1.863e−9** rad,
outside the 1e−9 budget, while one ulp below it is **9.313e−10** and inside.
Twelve accepted magnitudes — the circular cut, the exact half-turn tie, 37 turns
either way — are congruent to their principal target within **1.850e−10** rad,
and a run at an unresolvable initial phase raises having drawn **0** keys.

**The frozen rate, at every boundary.** The reviewer's three width-only rates
`0.1`, `1.0000000000000002` and `7.0` no longer give three trajectories under
one label: **seven** off-rate values — including both floats adjacent to the
frozen one — are refused by `ClockPath`, a rate smuggled past that constructor
with `object.__setattr__` is refused by `run_one_clock`, and four off-rate rows
plus a full-model row carrying the control's rate are refused by
`OneClockRecord`.

**The stale-path snapshot.** **15** stale-path seams are refused — the
reviewer's detuning edit, the drive as a whole, a drive swapped for the other
kind, each of a nested pulse's own fields, both schedule types planted
wholesale, and the stationary drive's own peak — every one of them with **0**
keys drawn, alongside 4 criteria edited past their own constructor. Four clean
runs reproduce their phase histories exactly against an independent walk.

**Reachability.** All four of the reviewer's impossible `DwellState` objects are
refused against five reachable controls; the conservation law is verified at
every one of **9** endpoints of a scripted history containing two resets, a
falling-edge eligibility loss and a commitment; **7** impossible records are
refused beside 6 reachable ones; and all **48** generated pulsed and stationary
records satisfy the same relations.

**The mutation battery.** Forty-six deliberately wrong implementations — a
wrapped authoritative phase, a dropped target lift, `stable_phase` called at
equality, free mismatch drift outside the tongue, a drift frozen at the wrong
end of the interval, an endpoint recomputed as `start + duration`, a split
stream keyed with the model label, coupling-scaled diffusion, coupling-scaled
contraction rate, coupling-scaled dwell, sample-counted dwell, retroactive
credit, a dwell clock that never resets, a removed contraction guard and
non-strict boundaries, an elapsed-time comparison with no roundoff allowance, a
roundoff allowance wide enough to be a slack budget, a dwell state that forgets
its last observed endpoint, a lift with no angular-resolution domain, a raw
factory that takes the contraction rate from its caller, and a stationary drive
faked as a very wide pulse, the frozen-rate requirement reduced to a sign
convention, the state conservation and latch relation removed, the record
conservation relation removed, the stationary zero-boundary rule removed, and
the execution snapshot reduced to a type check past four different stale seams,
the criterion snapshot reduced to one past a half-turn band and a negative
dwell, each of the seven public primitive wrappers reduced to its old
unchecked body, `advance` without its leaf rebuild, the state rebuild reduced to
its constructor past four path-inconsistent states, and each handoff
precondition reduced to "any boundary" — each caught by the property it breaks, with the unmutated implementation passing
the same property in every case.

**The noise, statistically.** 524,288 kicks — 64 clocks × 8192 finest steps —
have mean **−0.04σ**, variance **−0.64σ** of the declared `2·D·h`, and a
Gaussian-shape χ² of **41.5 on 41 degrees of freedom** (+0.05σ). The worst
autocorrelation over lags 1–8 is **1.71σ**, the worst of 2016 cross-clock pairs
is **3.46σ**, and 3584 pairs straddling a generated block boundary give
**−0.71σ**. Every statistic is a standardized deviate held to a predeclared five
standard errors, and the sample sizes and thresholds were fixed before the
numbers were looked at. Because the streams are keyed rather than seeded, these
are stable numbers rather than a draw: re-running reproduces them exactly.

**The split, deterministically.** Over 40,000 conditional splits at α = ¼ and ½,
the child means, the duration-scaled child variances and the sibling covariance
are all within **1.17σ**, and every single sample conserves its parent to
**0.25 of the 4-ulp bound** — including 256 samples at each of seven fractions
from 1e−12 to 1 − 1e−12. A four-crossing step yields five chronological leaves
with ancestry `L, RL, RRL, RRRL, RRRR`, four *distinct* split normals, and both
kicks and durations summing to the parent **exactly**.

**The elementary mesh.** Six crossings spread across eight finest steps — two
inside one step, two a ten-thousandth of a step apart, one at 0.999 — give 14
chronological leaves. Every coarse increment at strides 1, 2, 4 and 8 equals
that leaf sequence accumulated in time order, **bit for bit**, and the coarse
path derives exactly six split normals where the root-summing control derives
none while agreeing numerically in every group. Windows of 1, 3, 5 and 8 finest
steps give one sequence; 20,000 steps stream at **0.04 MB**.

**Boundary routing.** Over 500 steps of five meshes — the review's
`FinestMesh(1.0, 0.1)`, one rounding the other way, a negative origin, a
large-index mesh, and a dyadic control — there are **no overlaps and no gaps**,
every step start routes to its own step, every step end opens the next, and one
ulp below a start routes to the previous. The review's crossing at `1.2` gives
six leaves owned by one step with no sliver child, and four crossings — three
of them exactly on boundaries — give one leaf sequence under windows of 1, 2, 3
and 5 and under reversed adjacent windows.

**The duration law.** `FinestMesh(1e15, 0.2)` realizes widths **0.125 and
0.25** over 4000 steps — mean exactly the nominal 0.2, and neither of them 0.2.
Grouped by their own width, the kicks match the realized law at **−1.39σ** and
**+0.01σ**; the nominal law, computed beside them as the control, fails the same
test at **11σ**. A midpoint split of a 0.25-wide step has child means,
duration-scaled child variances and sibling covariance within **1.02σ**.

**The stable scale.** Across ten regimes every representable scale is positive,
finite, and within **5.0e−14 relative** of an independently computed
logarithmic reference — including three where the product spelling gives
exactly `0.0`: `duration = ulp(0)` (1.4e−162), both operands `1e−200`
(1.4e−200), and both subnormal (4.9e−324). In the `ulp(0)` regime the
conditional split is still conditional: **+0.45σ** from its law, where
proportional subdivision is **−16σ** away.

**Empty windows.** Fourteen public generation calls across two collapsed meshes
each raise `ValueError` having derived **zero keys**; a window collapsing at
step 100 of 500 fails at the request while its resolvable prefix still works.

**Addresses, again.** Eleven noncanonical spellings — leading zeros on the step
and on the trial, an explicit plus, inner and trailing whitespace, an
underscore separator, a Unicode digit, a trailing separator, a truncated and an
extended path, and a negative step — are each refused, while all 36 keys the
emitter produces round-trip and re-emit themselves exactly.

**Addresses.** All 3082 public leaf addresses over a full block, its boundaries,
three (trial, clock) pairs and a step past 70,000 are distinct — where the
private block seed, kept as the control, still groups 1024 leaves at a time —
and `leaf_from_key` returns the elementary leaf bit for bit at five addresses
including both sides of a block boundary.

**The mutation battery.** Seven wrong implementations — six named in the plan's
mutation contract and one an independent review found in the code — are run
through the same statistics, with the real implementation beside them as the
control: proportional subdivision fails the variance test at **−34σ** (control
+1.03σ); an independently sampled right child breaks conservation by **6.0e15
times** the ulp bound; a model-keyed split separates the paired full and
width-only models by 2.5e−3 where the real one separates them by **0**; one
shared bath for every clock fails cross-clock independence at **45σ** (control
2.68σ); a split that reads a default diffusion instead of the run's leaves a
**2.2e−17** descendant where the deterministic limit requires exactly zero — a
difference no tolerance would have caught; a block addressed by where the
*request* began changes value when the caller splits its loop; and a fraction
supplied independently of its crossing moves the child under one identity,
which the public boundary now cannot be given at all.

**Stationary calibration.** Clocks released from 32 initial phases spread around
the whole circle arrive at `arcsin(Δ/K)` to **4.4e-14 rad**, at eight detunings
from tongue centre to `|Δ|/K = 0.99`; the drift vanishes at both fixed points to
5.6e-16. The *other* fixed point is confirmed to repel, at `+λ`, to 3.1e-6
relative — without which the word "unstable" would be decoration.

**Relaxation rate.** The measured near-lock decay matches `√(K² − Δ²)` to
**4.1e-6 relative, worst case** over `|Δ|/K ≤ 0.9` at two couplings (K = 1.3 and
0.4). Critical slowing is confirmed out to `|Δ|/K = 0.9999`, where the rate has
collapsed by a factor of 30.8 and still matches theory to 1.0e-5.

**Outside the tongue.** No eligibility, no stable phase (NaN, not a number), and
a slip frequency matching `√(Δ² − K²)` to **3.8e-6 relative, worst case** at both
signs of the detuning. A clock just *inside* the tongue is run through the
identical code path as a discriminating control: it does not slip, and settles to
within 1.1e-14 rad of its stable phase.

**The envelope.** Exactly 1 at the peak, exactly 0 at both endpoints and off
support — checked as equalities, not tolerances. Symmetric to 4.4e-16; the area
matches `peak·T/2` exactly. Two channels sharing one envelope keep a
time-independent coupling ratio (to 5.6e-17), so polarization cannot leak into
the pulse *shape*.

**Tongue entry and exit.** Simulated eligibility opens and closes at the
closed-form crossing times of `K(t) = |Δ|`, every residual within **half a step**
— the honest bound, since the crossing is bracketed between the last ineligible
sample and the first eligible one — and falling cleanly as the step is refined
(1.00e-2 → 5.00e-3 → 2.50e-3 → 1.25e-3 across an eightfold refinement, over six
detunings including one at `Δ/K = 0.999`). Eligibility is a single contiguous
interval straddling the peak, not a flicker.

**Finite-time failure.** With `K_peak = 1`, `T = 40` and 16 release phases per
detuning, as `Δ/K` runs 0.2 → 0.995. The "reached the band" total is split into
three exhaustive categories named for exactly what is measured:

| `Δ/K` | e-folds offered | mean lag at peak (rad) | inside at entry | entered later | never entered |
| --- | --- | --- | --- | --- | --- |
| 0.20 | 18.06 | 3.2e-3 | 0.0625 | 0.9375 | 0.0000 |
| 0.60 | 10.35 | 2.3e-2 | 0.0625 | 0.9375 | 0.0000 |
| 0.90 | 2.77 | 2.7e-1 | 0.0000 | 0.6875 | 0.3125 |
| 0.98 | 0.56 | 4.6e-1 | 0.0625 | 0.2500 | 0.6875 |
| 0.995 | 0.14 | 5.0e-1 | 0.0625 | 0.1250 | 0.8125 |

Every one of these clocks *is* eligible at the pulse peak, so this is
admitted-then-lost, not never-admitted. Only the **never entered** column is
`lock-failed`; the other two reached the band, and whether they would hold it
long enough to commit is for the dwell stage to decide.

### Where "inside at entry" is measured, and why it matters

A tongue opens at a time the sample grid does not contain. Asking instead
whether a clock was in band at the first *stored eligible sample* answers a
subtly different question, because a clock can enter the band during the
interval between the true crossing and that sample — and with `store_every` set,
that interval is several steps wide.

So each trajectory is propagated to the **analytic entry time** and compared
there, against the boundary phase `sign(Δ)·π/2`. (Exactly at `K = |Δ|` the
stable and unstable fixed points have merged and `stable_phase` is NaN by
construction, but the limit of `arcsin(Δ/K)` is well defined, and it is the only
sensible target at that instant.)

This is not a cosmetic change. At `Δ/K = 0.90, T = 40` the first-stored-sample
rule reports **1/16** in band; exact-entry classification reports **0/16**. That
clock's true distance at entry is 0.0538 rad — outside the 0.05 band — so it
entered *after* the crossing and the sampled rule had miscounted it.

The classification is step-independent, and by a margin rather than a
coincidence: the entry phase agrees with a 64× finer reference to **1.5e-10
rad**, while no clock lies closer than **5.3e-4** to the lock-band edge. The
boundary is some 4.3 million times further away than the numerical error, so
refining the step cannot move a clock across it.

The category is called *inside the band at entry* and nothing more. It records
where a clock was at one instant; it does not establish what did or did not
happen before that instant, so it is not described as lucky or as a
non-synchronisation.

**The causal control.** A single duration cannot distinguish "ran out of time"
from "near-edge response is simply broken" — both predict the same failures. So
the identical 16 release phases are given longer pulses at fixed peak, with
nothing else changed. The same clocks come back:

| `Δ/K` | `T = 40` | `T = 160` | `T = 640` |
| --- | --- | --- | --- |
| 0.90 | 0.6875 | 1.0000 | 1.0000 |
| 0.98 | 0.3125 | 0.7500 | 1.0000 |
| 0.995 | 0.1875 | 0.4375 | 0.8750 |

(fraction ever reaching the lock band). Recovery is monotone, `Δ/K = 0.98`
returns to complete success by `T = 640`, and the fractions are **identical**
across a fourfold time-step refinement (dt = 0.008 / 0.004 / 0.002) — so this is
exposure time, not sampling resolution or integrator error. That reversibility
is what makes finite synchronization time the *cause* rather than a coincident
observation.

Broken out by category, the recovery is visible in both directions — the
never-entered population drains as the pulse lengthens:

| `Δ/K` | `T = 40` (entry / later / never) | `T = 640` (entry / later / never) |
| --- | --- | --- |
| 0.90 | 0.0000 / 0.6875 / 0.3125 | 0.1875 / 0.8125 / 0.0000 |
| 0.98 | 0.0625 / 0.2500 / 0.6875 | 0.1875 / 0.8125 / 0.0000 |
| 0.995 | 0.0625 / 0.1250 / 0.8125 | 0.2500 / 0.6250 / 0.1250 |

**Pulse duration.** At fixed peak the eligibility window is strictly
*proportional* to duration — `width/T` is constant to exactly zero across seven
durations and six detunings — so a longer pulse can never give a clock less
tongue time. A **fixed-area** sweep is run beside it and
deliberately fails the same monotonicity: at fixed area a longer pulse is a
weaker one, `K_peak = 2A/T` falls, and clocks drop out of the tongue entirely
(windows 1.96 → 3.43 → 5.32 → 4.10 → 0 → 0 → 0). The two sweeps answer different
physical questions and must be stored and labelled separately.

**Adiabatic tie-back.** A long slow pulse reproduces the stationary calibration
over its plateau, with the residual falling like `1/T²` (6.9e-4 → 4.3e-5 →
2.7e-6 for `T` = 100, 400, 1600), which identifies the agreement as adiabatic
convergence rather than three numbers that happen to be small.

**The analytic control.** Over 25 couplings spanning the full decade `K` = 1 → 10
with support `±80`, refining the grid from 2000 to 64000 points drives the worst
relative deviation of `flux/K²` from `π/2` monotonically down:

| grid points | spacing `h` | worst `\|flux/K² ÷ (π/2) − 1\|` | fitted exponent | `\|exponent − 2\|` |
| --- | --- | --- | --- | --- |
| 2000 | 8.0e-2 | 8.45e-3 | 2.00065603 | 6.6e-4 |
| 4000 | 4.0e-2 | 1.11e-3 | 2.00017300 | 1.7e-4 |
| 8000 | 2.0e-2 | 3.77e-4 | 1.99992271 | 7.7e-5 |
| 16000 | 1.0e-2 | 1.26e-4 | 1.99998460 | 1.5e-5 |
| 32000 | 5.0e-3 | 4.30e-5 | 1.99999753 | 2.5e-6 |
| 64000 | 2.5e-3 | **2.10e-5** | **2.00000086** | **8.6e-7** |

Both columns fall monotonically. Enlarging the support eightfold at fixed spacing
changes the flux by 2.8e-14 — the extra detunings are ineligible and contribute
exactly zero — so grid refinement and support size are cleanly separated effects.

The exponent fit is deliberately unconstrained, and the same fitting routine is
pointed at two quantities that are genuinely *linear* — the tongue width (fits
1.00000000) and the eligible-clock count (fits 1.000101) — where it returns 1,
not 2. A fitter that manufactured quadratics would pass the first half of that
check and fail this one.

### The declared safe range

`analytic.py` declares, before any fitting, that the flat-grid flux is
trustworthy only when the tongue spans at least 200 grid cells and stays within
a quarter of the spectral support. At the finest grid that range is
`K ∈ [0.5, 20]`, and the fitted decade `K ∈ [1, 10]` lies strictly inside it, so
a good exponent cannot have been obtained by cropping the sweep after the fact.

The *coarse* rungs of the ladder are deliberately outside that range — at 2000
points the safe range is only `K ∈ [16, 20]` — and they are reported anyway. They
are the evidence of convergence: a grid too coarse to resolve the tongue gives a
visibly worse exponent. Hiding them would remove the only reason to believe the
fine one.

**The race, against an independent walk.** For every trial of a sixteen-clock
population, the stop the race reports and the whole winner list it retains are
compared against a walk of every clock with no stop at all, reduced by an
independently transcribed minimum rule. They agree exactly, trial by trial —
not to a tolerance, since both are the same endpoint of the same mesh. Seven
synthetic completion-time sets, including a three-way tie and a tie that is not
at the minimum, reproduce the declared rule.

**Ties are real and are kept.** In a pinned twenty-four-clock configuration
several of eight trials end with two clocks completing on the same endpoint.
Both survive with their own mismatch and phase. A winner list trimmed to either
end, reordered, or extended by a clock that never completed is refused, because
the list is derived from the per-clock ledger rather than recorded beside it.
`raw_race.py` names none of twelve ways to draw or pick a value — no generator,
no `choice`, no `argmin` — so the module that decides who won has no way to
choose.

**Censoring, measured.** In one 64-trial diagnostic run the shadow sample
records **39** commitments over its eight trials where the race retained
**8**: thirty-one commitments that the channel's own stop made unobservable.
Every censored count is at most its uncensored one and every clock that
committed under the stop committed at the *identical* time without it, so the
censored walk is a prefix of the free one rather than a second realization.

**Byte-identical.** One twelve-trial configuration written at
`(clock_block, step_window)` of (1, 1), (3, 7), (16, 4096) and (5, 13) produces
identical manifests, ledgers, clock tables, shadow samples and close markers —
compared as bytes. The prepared ensemble survives an eight-fold timestep change
and a 400-fold diffusion change unchanged, so a refinement ladder reuses one
ensemble rather than comparing two.

**Memory.** A 64-trial, 16-clock, 80-step run peaks at **1.05 MB** and takes
26 seconds. Quadrupling the step count or the trial count does not move the
peak. The trials-by-clocks-by-steps array this replaces would be 0.66 MB at
that small diagnostic size and grows with all three; the live block is 16
clocks by 64 steps, 8 kB.

**The gate.** A real three-table run is refused after a trial row is deleted,
duplicated, or shuffled; after a clock row is appended or the clock rows are
shuffled; after the shadow sample is truncated; after one manifest field or one
marker digest is changed; when a second run is written over an existing one;
and when the manifest comes from a different run.

**Semantic forgeries, after every hash is recomputed.** An independent review
edited one cell of a closed zero-coupling run, recomputed every table digest and
every marker field, and was handed a `ClosedLedger` three times: a
`never_eligible` trial relabelled `lock_failed` with zero eligible and zero
inside clocks; a manifest detuning of `0.0` bound to clock and shadow detunings
of `0.75` and `-0.75`; and one eligible trial with **999 units of exposure inside
a 0.4-unit window**. All three are well typed, correctly counted and perfectly
hashed, and none describes the run the manifest froze.

All three now refuse, and the reason is that the quantities are *derived* rather
than read: the channel category comes from the authoritative clock rows, every
indexed detuning must equal the manifest's detuning for that clock exactly, and
exposure is recomputed from the manifest window and the channel stop through the
same expression the writer used. Hashes were never authentication; what closes
this is arithmetic.

**Time, which is the quietest thing to move.** A second review round moved
instants instead of counts, and got four more records accepted: a no-stop
exposure changed from `4.0` to `2.0`; a censored loser's first eligibility moved
to the manifest close, long after its trial had stopped; a same-key shadow first
eligibility moved off the value the censored walk had already observed on that
path; and a shadow row inside a band it had not yet entered.

Every row now carries a derived **observation end** — its own commitment, or the
channel stop for a censored row, or the window close for a no-stop row. Exposure
is that end minus the window's opening, exactly. Every recorded instant must lie
inside that row's own interval, not merely inside the run's. Events must be
causally ordered: eligible before inside, inside before a reset, a reset before
the commitment that followed it. And where the censored and no-stop walks of one
`(trial, clock)` cover the same interval they are *the same walk*, so their rows
must be identical field for field; where the censored one is shorter, its first
events may not move and any extra reset or commitment must fall after it stopped
looking.

A third round closed the rest of that family. A commitment already observed
under the stop must recur at the identical instant without it — checked before
the observation ends are compared, since deriving a moved commitment's exposure
also moves its end. Opportunity and endpoint counts are *jointly* monotone
rather than independently so: every counted endpoint closes an interval of
strictly positive duration, so zero eligible endpoints means zero opportunity
and vice versa, a strictly longer walk sees at least one more endpoint, an
unchanged eligible count requires unchanged opportunity, and added eligible
endpoints require added opportunity.

**The population is declared, not inferred.** The same review handed the race the
first path of a declared two-clock grid and got an ordinary outcome back — one
clock, identifiers `[0]`, category `never_eligible`. A trailing omission leaves a
*complete range* of clock identifiers, of a smaller population, so
`range(len(clocks))` could never catch it. The race now takes a
`ClockPopulation` carrying a declared size, clock order and mismatch list taken
from the frozen grid rule, `ChannelOutcome` retains that identity, and right-end
truncation, front omission, reorder, a repeated clock, a same-length detuning
substitution and a different declared drive all refuse before an outcome
exists.

---

## The continuum audit

The primary event process looks at the lock band only at integration endpoints.
A Brownian path can leave the band and come back between two of them, and the
endpoint scheme counts that segment as uninterrupted dwell. This is not a small
effect and it is not a guess: for a bridge between two endpoints at 0.09 with an
edge at 0.10, `D = 0.01` and `dt = 0.01`, the conditional probability that the
path crossed is `exp(−1) = 0.3679`. Adjacent timestep halvings can look like a
plateau while that error is still present, so "the answers stopped moving" is
not evidence of convergence.

Two independent things measure it, and they are deliberately unequal in
authority.

### The oracle: `killed_diffusion.py`

A phase diffusing inside a **fixed** band, killed the instant it touches either
edge. The backward Kolmogorov equation

```
dS/dt = (Δ − K sin φ) dS/dφ + D d²S/dφ²,   S(φ,0) = 1,  S(edge, t) = 0
```

is discretized with second-order central differences and a Crank–Nicolson step.
The initial condition is discontinuous at the corners, so the opening steps are
taken as damped implicit-Euler half-steps (Rannacher startup) before
Crank–Nicolson takes over.

It answers three questions: **survival** through a declared time, **exit-time
quantiles**, and **exit counts resolved by edge**. The third deserves a note,
because the obvious reading of "reset count" has no continuum limit at all: a
Brownian path that touches a boundary touches it infinitely often afterwards, so
the number of band exits an endpoint scheme observes *diverges* as the timestep
shrinks rather than converging. Killing at the first exit is what makes the
question well posed, and once the process is killed the reset count is the
absorption indicator — split by edge, so that it carries information survival
does not.

The oracle is calibrated against a case it does not compute. With `K = 0` the
drift is constant and the problem has a closed-form eigenfunction expansion,
written separately in the same file and sharing no code with the solver. The
suite runs both at three grid sizes: the worst relative departure is
**4.5 × 10⁻⁴**, falling by at least **3.77×** per halving of both grids —
second order, measured rather than asserted — and alive + absorbed-below +
absorbed-above closes to better than 10⁻¹².

**An analytic tail bound is not an error bound**, and the expansion is
certified against both halves. The tail is what truncation discards; the second
half bounds the floating evaluation, on a sum whose exponentially tilted terms
can cancel by a factor of 10¹⁰. An
independent review found the previous version accepting a value whose exact
reflection partner, `S_μ(x,t) = S_₋μ(L−x,t)`, differed by **1.63 × 10⁻⁵** at a
requested tolerance of 10⁻⁶ — while both evaluations reported tails of
5.04 × 10⁻²¹⁰. The analytic statement was true and the number was wrong.

Three things fixed it, in two rounds of review.

**The band coordinate is exact.** `u = (x − lower)/(upper − lower)` is carried
as a double-double — both the edge distance *and* the band width — and `n·u` is
reduced exactly to an integer plus a remainder in `[−½, ½]`. Everything upstream
of `1 − u` is magnified by `1/(1 − u)`, so at a start 1.1 × 10⁻⁶ from an edge a
single rounding of `x − lower` is 2 × 10⁻¹⁰ *relative*, and rounding
`upper − lower` costs the same again. An independent review found an accepted
value wrong by **1.33 × 10⁻¹⁴** — 6.6× its requested tolerance and 7.4× its
certificate — from exactly this. Carrying both exactly takes that same
evaluation's error to **9.5 × 10⁻²⁰**.

**Every derived quantity is scale-safe and checked.** The tilt is
`drift / diffusion / 2` and never `drift / (2 · diffusion)`: the second forms
`2 · diffusion` first, which at `D = 10³⁰⁸` is an infinity, so the tilt of a
heavily drifted process silently became **zero** and the series was evaluated
and certified as if there were no drift at all. Every derived width, split,
mode, squared mode and decay exponent is now validated finite **in the same
grouping the evaluation uses** — a guard that proves `D·h·h·t` safe while the
evaluation computes `D·(h·h)` proves nothing — and an unrepresentable one raises
a named domain error rather than warning, overflowing or returning a NaN.

**The certificate walks the operation graph.** Not a constant: each factor is
charged for what it does — the weight's arithmetic, each exponential for its
own argument magnitude, the amplitude's *absolute* error (which does not shrink
when its two exponentials cancel), the shape's reduction and sine, the
classical summation bound, and a final result floor. `series_error_certificate`
reports `tail`, `roundoff`, `total` and the cancellation factor;
`series_survival` refuses when the *total* exceeds the requested tolerance.

**The certificate is validated, not asserted.** A frozen adversarial sweep — 121
cells across three bands, five edge fractions, four drifts, two times and both
the 120- and 400-term regimes — compares every accepted value against an
independently written 90-digit reference. The worst ratio of true error to
certified error is **5.1 × 10⁻³**: the certificate dominates by at least 200×
everywhere it certifies.

The result floor makes this a **declared conservative refusal domain**: no
tolerance below about 1.8 × 10⁻¹⁵ can be certified for any evaluation, however
accurate. Refusing an accurate value is the safe direction; the alternative
would certify where the bound's own modelling assumptions are the largest term,
and a review has already falsified one such model here.

On the earlier reflection pair the symmetry violation is **8.7 × 10⁻¹¹**. Over a
four-tilt reflection sweep, 190 certified cells agree to 9.3 × 10⁻¹⁰ and 26
ill-conditioned cells are refused outright.

**The oracle imports `validation.py` and nothing else from this package.** Not
the stepper, not the noise source, not the dwell machine. It writes its own
Adler drift, and the verifier compares that expression against `model.py`'s
(they agree exactly over 401 phases). An oracle that imported the thing it is
the reference for would check an implementation against itself.

### What the endpoint scheme actually costs

One keyed Brownian tree of 2048 finest steps, read at three timesteps —
0.0078, 0.0039, 0.0020 — by 6000 paired walkers from three declared starts
inside the band. Every coarse increment is the **bit-for-bit left-to-right sum**
of the same fine leaves the finest level consumes, so the levels are one
realization read three times rather than three realizations.

Survival is **above** the oracle at every level and every start, which is the
only direction an endpoint scheme that misses between-step exits can be wrong
in.

**Every start is its own frozen ladder.** Four observables × three positions is
twelve budgets, twelve convergence requirements, and twelve results that each
have to pass on their own. An earlier revision collapsed the three starts to a
per-timestep worst case, which compares a different path at each level and can
hide a reversal at one start behind a worse value at another; an independent
review called that mutation-fragile and was right.

The measured refinement-error envelopes are reported **per unit, with the
identity of what produced them**:

| Unit | Absolute | Relative | Worst at |
| --- | ---: | ---: | --- |
| probability | 0.0285 | 0.078 | `survival`, start `x = 0.5844` |
| time | 0.1113 | 0.139 | `exit_quantile_p35`, start `x = 0.1844` |

There is deliberately **no** single scalar across units. The first draft had one,
and its headline value was an exit-time quantile error *in time* sitting in the
same maximum as three probabilities. A maximum over a time and a probability is
not a quantity, and a later production budget cannot compare it with anything.

The oracle's own discretization error, measured by refining it, is
**1.2 × 10⁻⁴** — a factor 37 below the smallest endpoint error it is the
reference for, and the check asserts that ratio rather than a bare threshold.

The envelopes are exported, not just printed: `ValidationVerdict.as_dict()`
returns a JSON-serializable report under schema **`dk-numerical-validation/v4`**
carrying the verdict, the digest of the *contract* and the digest of the frozen
*data*, the per-level measured and reference values, each level's unit, its
dependence-aware standard error, its paired standard error and its cluster
count, and the per-unit envelopes with their observable and position.

The live schema is **v4**. It reached that in three steps, each forced by a
shape whose meaning had changed while its tag had not: `v1` carried a single
scalar `envelope` mixing a probability with a time, and no position,
standard-error, cluster or dataset-digest semantics at all; `v2` added those;
`v3` added the frozen contract, budgets and blocking codes; `v4` made the
parser recompute from the frozen dataset. A consumer written against any
earlier tag must refuse this shape rather than read it under the older
meanings.

`parse_validation_report` is the **trust boundary**, and it *recomputes the
whole report from the frozen data* rather than reading it.

That is the v4 change, and it required a v3 step first. Two digests were not enough: a genuine
`numerical_no_result` could be rewritten **coherently** to a pass — verdict,
`passed` and `reasons` all changed together, every level and both digests
byte-identical — and the parser had nothing to check it against, because the
report carried no caps, no floors, no `require_decrease` settings and no
decision evidence. v3 carries the frozen contract, the frozen budgets, a
machine-readable `blocking` list of reason codes, and a per-level `span_error`.
The parser rebuilds the contract and budgets, checks they hash to the declared
budget digest, and re-runs **the same decision function `compare_refinement`
uses** over the levels.

v3 stopped short of the root. It recomputed the verdict from the report's *own*
uncertainty summaries, which were bound to nothing: zeroing every level's
standard, paired and span error, flipping the verdict and emptying the reasons
produced a parser-approved pass with **both digests byte-for-byte unchanged**,
because a self-declared hex string is not a proof. In v4
`parse_validation_report(payload, dataset)` takes the `ValidationDataset` as a
**required argument**, checks its digest is the one the report declares, and
recomputes the entire report — estimates, folded bootstrap, all three
uncertainties, ladders, envelopes, reason codes, reason text, verdict — and
refuses on any difference at all. Recomputing everything rather than checking a
list of invariants is deliberate: it is the only formulation that cannot be
incomplete.

`inspect_validation_report(payload)` remains for shape-and-local-consistency
checking, and is documented as **not** establishing the verdict.

The evidence travels as an object rather than inside the payload because of
size: the stationary S3 dataset holds twelve samples of six thousand clusters
at three levels — 216,000 observations, several megabytes of JSON for a report
whose purpose is to be read. The guarantee is the same, since the digest is
checked against the dataset and everything else is recomputed from it.

It also validates meaning rather than shape. A downstream consumer calls it and then acts on what it
returns, so an independent review was right to treat a lenient parser as the
place where a scientific block gets lost: the previous version accepted a report
carrying `verdict = "numerical_no_result"` beside `passed = true`, which is a
one-line route past a blocked result. It now enforces `passed == (verdict ==
"pass")`, reason/verdict consistency, two distinct well-formed 64-hex digests, a
closed top-level key set, closed level and envelope key sets, declared units,
string positions, finite measurements, non-negative errors, at least two
clusters, `absolute_error == |measured − reference|`, no duplicated
observable/position/timestep, one shared schedule across every ladder, one unit
per ladder, a zero paired error at every ladder's coarsest level, and envelopes
**rebuilt field-for-field from the levels** — exactly one per unit present, each
naming the row that actually attains its maximum.

### The frozen contract, and what freezing can mean

A tolerance supplied after the result it judges is a description of the result.
There is no way inside one process to prove a number was written down first, so
the mechanism aims lower and hits: `FrozenBudgets` hashes its own contents, and
`compare_refinement` refuses to produce a verdict at all unless the caller hands
back the digest it recorded earlier. A widened cap hashes differently and the
comparison raises.

**What is inside the digest is the part that matters, and it grew.** An
independent review pointed out that the caps were hashed while the two other
load-bearing thresholds were not: the required level count was a default
argument a caller could override to two, and the refinement schedule was
nowhere at all, so three levels at 0.040, 0.039 and 0.038 passed as a refinement
ladder. Both now live in a `RefinementContract` that is hashed with the caps:

| Frozen and hashed | Effect |
| --- | --- |
| the exact timestep schedule | a ladder that did not run on it is a no-result |
| the refinement ratio | a non-halving schedule cannot be **constructed**, let alone frozen |
| the level count | three minimum, enforced where no caller can reach it; the override is gone |
| the sampling design | cluster unit, cluster count, auxiliary replications, uncertainty method, resamples, coverage factor |
| each observable's unit and position | a maximum is taken over one quantity, at one start |

The same review showed that a bare `noise_floor` was **permission for systematic
growth**: a ladder of 0.010, 0.020, 0.030 under a floor of 0.011 grew at every
refinement and passed with no reason recorded, because each increment was under
the floor. Three rules replace that:

1. the error at the finest timestep is no larger than at the coarsest, within
   **one** allowance for the whole ladder — not one per step;
2. at most one adjacent pair may go the wrong way; and
3. that one must fit inside the declared resolution **plus** the frozen sampling
   design's own cluster-bootstrap standard errors on the two levels involved.

Then the uncertainty itself moved inside. `compare_refinement` used to take
`RefinementLevel` objects carrying **caller-computed** standard errors, and
under one unchanged contract and one unchanged digest the ladder
`0.010, 0.040, 0.009` returned `numerical_no_result` with those errors set to
0.0 and `pass` with them set to 0.1. Nothing measured had changed. The
signature now takes a `ValidationDataset` of frozen cluster-level observations
— the per-walker or per-master-trial numbers themselves, hashed — and computes
the statistic, the bootstrap and the verdict here. There is no argument through
which an uncertainty can arrive, and the old signature raises.

The bootstrap is **paired**, and the distinction is not pedantic. One resample
of clusters is evaluated at *every* level, so a bootstrap replicate of a
level-to-level change is the change of two numbers that saw the same clusters.
The previous arithmetic took `hypot` of two marginal errors — the
independent-level formula — which inflates the allowance exactly when the levels
are dependent, and they always are, being one realization read three times. On a
fixture with cluster-common scatter whose errors stay one-signed, the two differ
by orders of magnitude. The resampling key is frozen in the sampling design with
everything else.

It resamples **the statistic the gate consumes**: the folded absolute-error
contrast, not the signed estimate. An earlier revision resampled the signed
estimate and argued conservatism from `||a|−|b|| ≤ |a−b|`; an independent review
pointed out that a pointwise inequality between two random variables says
nothing about their standard deviations, and produced a two-cluster sample whose
signed adjacent contrast has standard error **exactly zero** while the folded
contrast the gate actually uses has **0.699**. The suite carries that sample as
a regression.

The frozen data carries **explicit physical identities**: a master-trial
identifier per cluster, a clock-and-replication identifier per member, hashed
into the sample digest. Without them a sample of twenty-four identical rows —
ordinary for a binary observable — is indistinguishable from one trial copied
twenty-four times, and no digest over the observations alone can tell. Duplicate
identities, missing ones, relabelled ones, reordered ones, and two observables
measured over different clusters under one sampling contract are each refused.

The frozen data is an **immutable snapshot**: every array is canonicalized to
little-endian binary64, serialized to a `bytes`, and reconstructed with
`np.frombuffer`, so the whole base chain reads *ndarray → ndarray → bytes* and
there is no writable owner anywhere in it to find. Rebuilding a dataset re-runs
that snapshot.

A read-only NumPy view was not enough. NumPy lets the *owner* of an array
restore its write flag, and an independent review walked one link out — 
`sample.values.base.setflags(write=True)` — wrote zeros through it, turned a
genuine no-result into a pass, and got the forged pass accepted by wrapping the
payload in a mapping subclass whose `items()` performed the write during the
parser's serialization. `bytes` closes the route at the root rather than at each
flag. The parser takes its trusted dataset snapshot **before** any
untrusted payload method runs, then snapshots the payload to an owned plain
dict, captures the digest once, and re-checks it after recomputation.

The manifest encoding is **injective**: every label is length-prefixed under its
field name before hashing. One Unicode policy covers every hashed field —
identities, members, baseline members, observables, positions — refusing the
categories `Cc`, `Cf`, `Cs`, `Co` and `Cn`, so U+0085, U+202E, U+200D, lone
surrogates and private-use code points are all rejected while ordinary
non-ASCII names are not. **No normalization is applied**, so a composed and a
decomposed spelling are distinct identities and hash differently. The
previous encoding joined labels with NUL, and an independent review showed
`("a\0b", "c")` and `("a", "b\0c")` producing identical sample *and* dataset
digests with identical observations — SHA-256 had not collided, two different
physical manifests had been flattened to one byte string before hashing.

And the absolute cap is applied to the **uncertainty-aware** finest-level bound
— the point estimate plus two bootstrap standard errors — not to the point
estimate alone: 0.0112 against a cap of 0.10 reads like a comfortable pass until
you notice its standard error is 0.0829.

Read the digest as **"a moved threshold cannot be silent"**, never as "a
threshold cannot be moved". The suite asserts both halves: the old digest
against a widened budget raises, and the widened budget under *its own* digest
passes. The second assertion is there so nobody reads the mechanism as stronger
than it is.

Failure produces `verdict = "numerical_no_result"` with the reasons attached. It
never produces a widened budget, a dropped level, or anything resembling a
physical finding. Eight deliberate failures are exercised: an overlarge
discrepancy, an error growing under refinement, **repeated sub-floor growth**, a
ladder off the frozen schedule, a two-level ladder, a level from a different
sample, an observable measured in the wrong unit, and a budgeted observable that
was never measured. A waived observable is separately shown to be unable to
rescue a convergent one.

### The audit: `moving_band_audit.py`

During a pulse the band **moves** — the stable target sweeps as the coupling
rises and falls, fastest near tongue entry and exit, which is exactly where the
error is largest. The stationary oracle does not cover that case, so a diagnostic
does.

On every wholly-eligible elementary interval the replay freezes the drift,
represents both admissible edges as straight lines between their endpoint values,
and computes

```
p_upper = exp(−(upper_start − φ_start)(upper_end − φ_end) / (D · dt))
p_lower = exp(−(φ_start − lower_start)(φ_end − lower_end) / (D · dt))
```

combined as `min(p_lower + p_upper, 1)`. Subtracting a straight moving boundary
from a conditioned path leaves a bridge against a *constant* boundary, which is
the entire reason the edges are represented linearly and the reason a curved
"better" representation is a defect. The verifier checks the expression three
ways: against the worked example above, against an independently derived
reflection-principle expression (agreement to 6.7 × 10⁻¹⁶ over five fixed and
moving geometries), and against direct simulation of the bridge, which must
approach the closed form from below and does — 0.1676 at 60 substeps, 0.2126 at
2400, closed form 0.2231.

**Admissible geometry** is the proximity band intersected with the strictly
contracting side of the flow — the same rule `commitment.py` applies pointwise,
restated as a set so the bridge has edges to be linear between. At the pressure
test's own `K = 1`, `Δ = 0.999`, where the stable phase (1.5261) and the
repelling one (1.6155) are 0.0895 apart and a 0.1 tolerance contains both, the
admissible interval contains the first and not the second.

Where the circular representation is ambiguous, where the interval disappears, or
where its topology changes inside the step, the replay **resets conservatively**
rather than choosing the branch that keeps the clock locked. Choosing the
favourable branch is a named mutation and is required to fail.

**Three rules make "diagnostic" mechanical rather than aspirational:**

1. **Reset-only.** The replay may add a dwell reset. It may never remove one,
   move a commitment earlier, or turn an unresolved run into a committed one.
   `AuditedRun` refuses to be constructed if any of those happened.
2. **Pathwise subset.** Because resets only delay, an audited commitment is one
   the endpoint scheme also made, at a time no earlier.
3. **Disjoint randomness.** The auxiliary uniforms live under
   `dk-band-audit/v2/audit`; the physical tree lives under
   `dk-phase-noise/v3/physical`. Neither string is a prefix of the other, and
   256 physical against 512 audit keys under a *deliberately shared dataset
   name* collide nowhere. Changing the auxiliary replicate leaves the primary
   history bit-identical, which is what proves the audit reads no physical key.
4. **Canonical interval identity.** An audit key names one physical interval.
   In v1 it did not: the address was `(finest step, piece ordinal)`, the coarse
   refinement levels are built on the same finest mesh, and an independent
   review found **265 keys shared across strides 1, 2 and 4, every one of them
   describing a different bridge over a different duration**. The key now
   carries the interval's exact endpoints in `float.hex` — the same canonical
   encoding a crossing identity uses — so equal keys imply equal intervals
   across stride, window partition, pulse ordering and refinement. The suite
   runs the review's own reproduction in both forms: the v1 address collides 265
   times, the canonical one zero. Different refinement levels therefore draw
   *independent* auxiliary uniforms, which is the right default — they are
   different bridges over different durations, the physical realization is what
   is paired between them, and the audit's own Monte Carlo variation is resolved
   by replication rather than by coupling.

One accounting subtlety is worth stating, because getting it wrong looks like a
result. The audited walk **outlives** the primary one: its dwell was delayed, so
it keeps going past the endpoint scheme's own stop. Resets it accumulates after
that point are ordinary endpoint resets on exposure the primary never had, not
resets the audit created. The reset-only comparison is therefore frozen at the
shared stopping time, and `AuditedRun.added_resets` reports that difference and
not the difference of two walks of unequal length.

### The reduced matrix, and what converges

Fifteen cells, frozen in source before any audit number existed: central,
interior and near-edge clocks; both mismatch signs off centre; weak, intermediate
and strong noise. Both pulse edges are covered by construction — every cell runs
the whole pulse — and that is *measured*, by requiring every cell to report
eligible intervals on both the rising and the falling side of its window.

Refinement is paired here too: `paired_leaves` reads one keyed tree at three
resolutions by summing its leaves, with the coarse elementary mesh built as the
coarse grid **union the same exact crossings**, so no coarse interval ever
straddles an eligibility change.

Every cell is **run, gated and reported**, not merely reached: twelve master
trials at each of the three paired timesteps, at that cell's own noise strength
and on its own clock identifier, with its own frozen ladders for every
convergent observable — **against the same caps as the pooled ladder**. An
earlier revision ran one initial phase per cell and gated nothing, and an
independent review found a cell at 0.125 against the pooled 0.10 cap passing
with residual zero. The correction is the gate, not a wider cap for the cells
that fail it.

The pooled ladder is **40 master trials × 3 independently keyed clocks × 3
paired timesteps × 4 auxiliary replications**, and it is an *aggregate display*
beside the per-regime gate rather than a substitute for it. The three clocks
shared one clock identifier in an earlier revision, which made their physical
kicks and their audit uniforms identical and the ensemble far more dependent
than its sample size suggested; they are independent now, and the master trial —
the initial phase, common to all three — is the cluster the bootstrap resamples.

Five observables are gated:

| Gated observable | 0.032 | 0.016 | 0.008 |
| --- | ---: | ---: | ---: |
| commitment-probability shift | 0.0063 | 0.0042 | 0.0083 |
| survival shift at 45% of the pulse | 0.0104 | 0.0083 | 0.0000 |
| survival shift at 60% | 0.0104 | 0.0104 | 0.0000 |
| survival shift at 80% | 0.0229 | 0.0125 | 0.0146 |
| commit-time quantile (p20) shift | 0.1344 | 0.1152 | 0.0112 |

Every audited commitment is a delayed copy of an endpoint one at every level and
in every cell, and each upward step is decided against the paired bootstrap's
own standard errors rather than a bare floor.

### What the audit currently concludes: `numerical_no_result`

**The reference audit does not pass, and reporting that is the result.**

The pooled ladder is blocked by its own uncertainty: the p20 commit-time shift
has a point estimate of 0.0112 at the finest timestep and a cluster-bootstrap
standard error of 0.0829, so its two-sigma upper bound is 0.177 against a 0.10
cap. The point estimate says where the estimate landed; it does not say the
quantity is below the cap, and this reference sample cannot tell them apart.

The per-regime ladders are blocked in fourteen of fifteen cells, mostly on the
same uncertainty-aware bound at twelve master trials per cell.

Both are `numerical_no_result`. That is a scientific outcome, not a broken
check: the sample is a *reference* sample sized to exercise the machinery, and
the frozen digest exists precisely so that the alternative — widening a cap
after seeing which cells fail — is a visible act rather than an available one.
The verifier's own check passes, because what it asserts is that the machinery
is sound and the outcome is reported honestly; it separately refuses a
no-result whose reasons are about the machinery rather than the numbers.

The added-reset count does **not** fall — it rises, 0.48 → 0.61 → 0.98 — and its
budget says so explicitly, with `require_decrease = false` inside the digest,
rather than quietly excusing it. That is the same mathematical fact the oracle's
edge-resolved exit counts exist to work around: the endpoint reset count has no
continuum limit, so its error cannot converge. What could converge is the effect
on the dwell — whether a commitment happened, when, and how much survival moved
— and those five are the observables held to a shrinking error. The suite
asserts separately that a waived observable cannot rescue a convergent one.

### Staying inside the probability domain

Two guards and a domain check, because an "oracle" that returns a negative
survival is worse than no oracle. The cell Péclet number `max|f| dx / (2D)` must
be below one, or central differencing oscillates in space; the advective Courant
number `max|f| dt / dx` must be at most one, or the drift is unresolved in time.
Whatever survives both is then **required to be a probability**, and is refused
rather than clipped if it is not — clipping would present a grid that solved a
different problem as one that solved this one badly.

The Courant condition is the one an independent review's counterexample
violates: `D = 0.001`, band `(−1, 1)`, drift `0.05`, eighty cells and five steps
to a horizon of a hundred has a cell Péclet of 0.625 — accepted — and a Courant
number of **40**, and it returned a survival of `−0.0043` beside an edge-exit
probability of `1.0043` while all three fields still summed to one. Closure is
not a domain check: `lower_exit` is defined as the complement, so two fields can
be impossible while the three of them close exactly.

The strict Crank–Nicolson positivity condition `D dt / dx² ≤ 1` would also have
refused that grid — and would have refused every grid this package uses, since
the S3 reference grid sits at 94. It is a monotonicity condition, not an accuracy
one, and enforcing it would have traded a real defect for a false one.

The closed-form series has its own version of the same problem. Its public API
used to accept any term count at any time, and at `band = (0,1)`, `D = 0.1`,
`mu = 2`, `t = 10⁻¹⁰`, four hundred terms return **−2.11** at `x = 10⁻⁴` and
**−15.75** at `x = 10⁻³`, where the true survival tends to one. It now
**refuses** an evaluation whose total certified error — the rigorous analytic
tail plus the operation-graph bound on the floating evaluation — exceeds the
declared budget, with the nested truncation at `terms - 2` reported alongside as
corroborating evidence. The short-time regime needs a term count that grows
without bound, so there is no count the function could silently choose instead. Both absorbing edges are pinned to exact
zero at every time.

### Bounded memory

The first implementation of this ticket peaked at **4.65 GB** resident, measured
by an independent review with `/usr/bin/time -l`, and it was right to call that
a machine requirement rather than a scientific one. Two paths were responsible:
the endpoint ladder materialized a 6000 × 2048 kick array, and the direct bridge
simulation built four 40000 × 2400 arrays at once.

Both are streamed now. The ladder processes 1500 walkers × 1024 finest steps at
a time, with every stride consuming the same block before it is discarded, so the
levels stay paired while nothing is kept but per-walker outcomes. The bridge
simulation runs 1000 paths at a time and keeps an integer hit count. The largest
array anything here allocates is **19.2 MB**; the unstreamed forms would have
been 3.07 GB and 98 MB.

| | Before | After |
| --- | ---: | ---: |
| peak resident set, whole suite | 4.65 GB | **0.63 GB** |
| largest single allocation | ~1.5 GB | **19.2 MB** |
| wall time | 53 s | 65 s |

`check_bounded_memory` asserts both: the declared block sizes against a
per-allocation bound, and the measured peak resident set against a whole-suite
bound. It is registered late in the file on purpose, since `ru_maxrss` can only
report on the checks that have already run — and it is the one row in the suite
whose residual is a measurement of the machine rather than of the code, so it is
the one row that is not byte-identical between two runs.

### The mutations

Each defect the plan or its reviews name, run against the check that must catch
it: a halved series coefficient stops matching the solver; a proximity-only band
admits the repelling phase the real band excludes; a midpoint-curved boundary
departs from the reflection oracle by 0.263; favourable topology selection turns
a conservative reset into a keep and is refused by the nearest-centre rule; the
zero-diffusion-before-clearance ordering reports certain absorption as certain
survival, while the shipped expression does not; audit-created commitments,
earlier commitments and removed primary resets are refused where the record is
created; a key forged into the physical schema does not resolve there; two
intervals of different duration no longer share one canonical audit address or
one uniform; a mismatched budget digest raises rather than producing a verdict; a
caller-supplied level count is a `TypeError`; and a level reporting an undeclared
cluster count is refused rather than compared.

### What this does not say

What passes is the **stationary** gate: the endpoint scheme resolves a fixed
band's exit law to its declared reference tolerances at every declared start,
with the caps applied to uncertainty-aware bounds.

What does **not** pass is the moving-band audit. Its pooled ladder and its
fifteen per-regime ladders both return `numerical_no_result` at this reference
sample. The correct reading is "this sample cannot bound the moving-band
discrepancy below these caps", not "the discrepancy is large" and certainly not
"the discrepancy is small".

Nothing here says anything about any coupling exponent, any population, any
hazard law, any detector outcome, any measurement, or Born's rule.

The budgets exercised are **declared reference values** chosen to test the
freezing mechanism and to record the envelope this resolution actually delivers.
The production budgets come from a later ticket's power calculation, and no
result here may be read as having passed them. The audit's added-reset count is a
reported diagnostic with no continuum limit; it is not evidence of convergence
and cannot substitute for any gate.

---

## The one-channel race and the run on disk

Everything before this section walked **one** clock. This is where `N` of them
race, and where what happened gets written down in a form that cannot be read
until it is finished.

### What a first clock is, and what it is not

**A first clock is a model stopping event.** The trial ends because this
one-channel first-passage study defines it to end there. It is not a detector
click, an absorption, a measurement outcome, a unique actuality, a microscopic
bath, one side of a two-channel competition, or a Born-rule result. Nothing
physical forbids the other clocks from having continued — and the no-stop shadow
sample below walks trials in which they do.

`raw_race.RACE_NON_CLAIMS` carries the full list, `raw_runner.RAW_NON_CLAIMS`
joins it to the one-clock list, and every manifest a run writes contains all of
them, so a ledger read on its own still says what it is not.

### The race: `raw_race.py`

`N` clocks on one immutable exact grid — declared support, midpoint origin,
frozen parity, `n_points` nodes — each with its own counter-keyed noise stream
and its own prepared phase. The channel stops at the earliest endpoint on which
any clock completes its fixed dwell.

**Every clock that completes on that endpoint is retained.** Their identifiers,
mismatches, phases, band-entry times and tongue-entry times all reach the
ledger, and nothing anywhere picks one of them. A same-endpoint tie is a fact
about the timestep, not a coin to be tossed: co-commit multiplicity has to
shrink or reach a stable negligible level under refinement before the
winner-mismatch distribution means anything, and a race that broke its own ties
would have destroyed the evidence for that gate. `raw_race.py` contains no
random source at all — no generator, no `choice`, no `argmin` — and the verifier
scans it to say so.

**The race walks twice, and the second walk is the point.** The stop time is not
known until every clock has been asked, and the per-clock diagnostics have to be
censored *at* the stop. A single pass would have to keep every clock's
time-resolved history in order to truncate it afterwards, which is the
trial-by-clock-by-time cube the plan forbids. So `first_completion` walks each
clock and keeps one number — that clock's own completion time — and
`censored_exposure` walks them again, stopping at the channel stop, and keeps
the diagnostics. Because a leaf belongs to its key and not to the walk that
reached it, the second pass re-derives bit-for-bit what the first one saw. The
cost is named rather than hidden: up to twice the clock-steps of a single-pass
implementation that did not censor honestly.

### Four categories, and nothing conditioned away

| Category | Meaning |
| --- | --- |
| `committed` | at least one clock completed an uninterrupted dwell |
| `never_eligible` | no clock ever entered its tongue |
| `lock_failed` | some clock was eligible; none ever reached the lock band |
| `dwell_failed` | some clock reached the band; none completed the dwell |

Mutually exclusive, exhaustive, and derived from the per-clock ledger rather
than recorded beside it: `ChannelOutcome` refuses a category that contradicts
its own clocks. The unresolved fraction is a **result**. A race that reported
only its commitments would produce a commitment probability of one and a
survival curve of zero length, and both would look like measurements.

### Censoring, exposure, and the no-stop shadow sample

Per-clock statistics in a first-passage race are *informatively* censored. Once
one clock stops the channel, every other clock's later opportunities are
unobserved — so a mismatch that produced no commitments and a mismatch that
never got the chance are the same number in a raw count, and telling them apart
is exactly what a winner-mismatch distribution is for.

Two things are recorded against that. Each clock carries its **exposure**, the
observed time it actually had, and its **eligible time**, the part of that
exposure spent strictly inside its tongue — the opportunity rather than the wall
clock, and the denominator an exposure-normalized statistic needs. A clock whose
tongue never opened has zero of it, which is the honest reading of its zero
commitments.

And a **predeclared no-stop shadow sample**: the first `shadow_trials` trials,
re-walked on *the same physical keys* with the stop switched off. Same
namespace, same trial and clock identifiers, same mesh, therefore the same
prepared phases and the same Brownian leaves — so the censored walk is a
*prefix* of the shadow walk and the comparison is pathwise rather than
statistical. The size is a configuration field, frozen at the input boundary and
written into the manifest before the run starts, because choosing it after
seeing which mismatches looked quiet would make it useless for the one thing it
is for. It is diagnostic: nothing in it is a channel outcome and nothing in it
enters a category count. A commitment recorded there is one that in the actual
race *did not happen*, because the channel had already stopped.

### Reproducibility, and what a block size may not do

A run is executed in **live blocks**: `clock_block` clocks held at once,
`step_window` finest steps of leaves generated at a time, and a losing clock's
walk abandoned as soon as it cannot win. All three are optimizations and all
three are ways to accidentally change a result, so none of them is allowed to.
The same configuration written at `(1, 1)`, at `(N, 4096)` and at anything
between produces **byte-identical** manifests, ledgers, clock tables, shadow
samples and close markers.

That is only possible because the manifest freezes no block size and no
timestamp. Both are deliberately absent: a number that cannot change a result
has no business being frozen into the immutable description of a run as though
it could, and a timestamp would make two runs of one frozen configuration differ
while being the same run. The block sizes are reported by `RawRunReport`
instead, which is not part of the record.

The prepared ensemble is keyed **without the mesh** —
`PhaseNoiseStream.ensemble_prefix`, `initial_phase_key` and `initial_phases`
name a phase by dataset, trial and clock and by nothing else. Halving the
timestep therefore reuses the same initial phases, which is what a paired
refinement ladder needs: otherwise it would be comparing two ensembles as well
as two timesteps.

### Bounded memory

Nothing here allocates an array indexed by trials, and nothing allocates one
indexed by both clocks and steps. Leaves arrive in windows and are dropped as
they are consumed; a trial's result is a row, written and dropped. Peak
allocation follows the live block: quadrupling the step count or the trial count
does not move it, and the verifier measures that rather than asserting it.

Nothing is retained between trials at all: the authoritative per-`(trial, clock)`
rows are written and dropped as they are produced, and the per-clock summary is
derived at *read* time from those rows rather than accumulated during the run.

### The run on disk

```
<directory>/
  manifest.json   the frozen description of the run
  ledger.csv      one row per trial
  clocks.csv      one authoritative row per (trial, clock), censored
  shadow.csv      one row per (shadow trial, clock), uncensored
  CLOSED.json     the close marker, written last
```

Written under a package-scoped results directory that `.gitignore` already
excludes, so generated output never lands beside source.

Written by `write_raw_run(config, run_name)`, which takes a **run name** and not
a directory, and read by `open_raw_run(run_name)`. There is no argument through
which an absolute path or a parent traversal could arrive.

**Scope is a capability, not a check.** Resolving a name gives you another name,
and a review defeated that four ways: replacing the run directory after
validation redirected both boundaries; five filenames that were symbolic links
to an external run let the reader return twelve rows from outside the package;
and five *broken* links — which `Path.exists()` reports as absent — made the
writer create all five files externally.

So the package directory, fixed at import from `__file__`, is the only name
resolved by string. The results root and the run directory are opened beneath it
as directory descriptors with `O_NOFOLLOW`, and every one of the five files is
opened relative to the pinned run descriptor: read-only and required to be a
regular file, or `O_CREAT | O_EXCL` on write. A link standing in for the run
directory or for any single table is refused rather than followed, a broken link
refuses on `O_EXCL`, and nothing that happens after the directory is opened can
redirect anything, because there is no name left to replace. The code fails
closed on a platform without those primitives.

**The marker is written last, and that ordering is the protocol.** A run that
was interrupted, killed, or is still going has no marker at all, so a reader
knows the run is unfinished before it hashes anything. Each table is hashed a
line at a time *while it is being written*, so the digest is over exactly the
bytes that reached the disk and never over a reconstruction of them.

### The authoritative table

`clocks.csv` holds **one row per `(trial, clock)`** — the censored history of
every clock in every trial, not a summary of it:

```
trial,clock,detuning,category,committed_at,first_eligible_at,first_inside_at,
last_reset_at,eligibility_entries,band_entries,dwell_resets,endpoints,
eligible_endpoints,inside_endpoints,exposure_time,eligible_time,final_phase,
won,co_completed
```

`shadow.csv` carries the same fields without `won` and `co_completed`, since a
no-stop walk has no channel stop and therefore nothing to win.

An earlier revision wrote an aggregate instead — one row per clock, totals over
the whole run. Totals are not a record: two runs with identical totals and
different per-trial histories support different conclusions about first-winner
censoring, and a closed run that kept only the totals could never afterwards be
asked. The winner's last-reset time was lost entirely, and the plan's committed
record asks for it by name.

So the summary is **derived, never stored**. `ClosedLedger.clock_totals()`
computes it from rows the gate has already reconciled, which means there is one
place each number lives and it is the place that was checked.

An example ledger row set, from a twelve-clock diagnostic run:

```
trial,category,stop_time,commit_time,censored,co_completion,winners,winner_detunings,...
0,committed,-0.2999999999999998,-0.2999999999999998,true,1,4,-0.875,...
1,committed,-0.8499999999999999,-0.8499999999999999,true,1,9,0.375,...
45,committed,-0.75,-0.75,true,2,5|6,-0.625|-0.375,...
```

```
clock,detuning,trials_eligible,trials_inside,wins,co_completions,band_entries,dwell_resets,exposure_time,eligible_time
0,-1.875,0,0,0,0,0,0,83.39999999999999,0.0
9,0.375,64,18,9,0,21,3,83.39999999999999,50.6...
```

`winners` is a list because ties are real: the pipe-separated cells
(`winners`, `winner_detunings`, `winner_phases`, `winner_band_entry`,
`winner_first_eligible`) are all indexed by the same clock order, so zipping
them zips one clock's fields together. Floats are written through `repr`, the
shortest decimal that reads back as the identical binary64 — a rounded cell
would make the ledger a lossy record of its own run, and the digest would then
certify bytes that no longer say what the run said.

### The gate a consumer must pass

`raw_ledger.require_closed_ledger` is the door, and `raw_runner.open_raw_run`
is the filesystem convenience in front of it. Both return **parsed rows**, never
bytes, so a caller cannot reparse them with a looser parser than the one that
was checked. Seven things are established before a single row is handed over:

1. the marker is exactly a `CloseMarker`, not a subclass;
2. it is rebuilt, so a marker edited after construction fails;
3. the manifest bytes hash to the digest the marker states, and are verified
   *before* they are parsed;
4. every table's bytes hash to its own declared digest;
5. every table parses — header exactly the declared columns in the declared
   order, every row the declared width, every cell the canonical spelling of its
   value, every category one of the declared four;
6. the manifest's count, the marker's count and the rows actually present are
   one number, for all three tables; and
7. the row **identities** are exactly the declared ones — trials `0..T-1`,
   clocks `0..N-1`, shadow rows the predeclared trial-by-clock grid, each once
   and in order.

Step 7 is what catches a table with the right number of rows and the wrong rows
in it, which is what two runs' output interleaved into one file looks like.

Ticket 01's gate could only establish that a marker agreed with a manifest about
a number neither had counted, and its own docstring said the row-level form
should replace it once a parser existed. It has. Both schema versions moved to
2, and a version-1 marker or manifest is now refused rather than partially read.
None of this is authentication: anyone able to rewrite a table can recompute its
digest. It is protection against truncation, interleaving, mismatched pairs and
accidental edit, which is what actually goes wrong.

### What a run records about its own numerical status

Three manifest fields, so that a ledger cannot later be read as though it had
passed a gate that was never applied to it:

| Field | Value |
| --- | --- |
| `endpoint_status` | `passed at the declared reference tolerances` |
| `moving_band_status` | `numerical_no_result pooled and in fourteen of fifteen per-regime cells at the frozen reference sample` |
| `numerical_gate` | `diagnostic_only` |

Ticket 04's implementation is closed and its machinery checks pass. Its
*scientific comparison* is a frozen `numerical_no_result`: the pooled
moving-band ladder and fourteen of its fifteen per-regime cells cannot bound the
discrepancy below their frozen caps at the reference sample. The correct reading
is "this sample cannot bound the moving-band discrepancy" — not that the
discrepancy is small, and not that it is large.

So a run written here exercises the race and the ledger and nothing else. It has
not passed a production numerical budget: ticket 07 froze one and it is **not
met**, with a `numerical_no_result` disposition and no feasible matrix among the
options searched. No survival curve, no commitment law, no
coupling exponent and no scaling claim may be drawn from these tables.

---

## The separate analysis layer

Three modules — `observables.py`, `analysis.py`, `compare.py` — read a closed
run and produce survival, hazard and a coupling exponent. They are a *separate
layer* in the sense the plan means: the raw import graph cannot reach any of
them, the verifier starts a fresh interpreter and confirms that importing
`raw_runner` loads none of the three, and nothing in the layer opens a file for
writing.

### The consumer boundary, and the hash on both sides

`compare.open_for_comparison(run_name)` is the only door. It goes through
`raw_runner.open_raw_run`, which reads the five files once behind pinned
directory descriptors and hands the bytes to
`raw_ledger.require_closed_ledger`; a run with no close marker is *open* and is
refused before anything is hashed, and a run whose manifest, tables, row counts,
row identities, schema, digests or row semantics disagree is refused after
everything is.

What it returns carries a `LedgerFingerprint`: the manifest digest and all three
table digests, as the close marker declares them and as the gate has already
recomputed them from the bytes actually read. `ReadOnlyRun.reverify()` re-opens
the run and requires an identical fingerprint, so a comparison can be bracketed
by a before-and-after hash rather than by a promise. The verifier does exactly
that, on the bytes on disk rather than through the object that was supposed not
to touch them, around a full analysis — counts, risk table, survival, cumulative
hazard, winner mismatch, all four comparators, the rising/falling split, a
discrepancy score, an exponent fit and a scaling verdict. All five digests are
unchanged.

It also scans the three modules' syntax trees for any name that could write —
`open`, `write_text`, `os`, `shutil`, `tempfile`, `subprocess` and the write-mode
flags among them — at every binding and call position. None appears.

### Risk, events, censoring, width — always all four

`events / at_risk` is a **discrete conditional failure probability**. It is
dimensionless and it depends on how wide the bin is.
`events / (at_risk × width)` has units of inverse time and is the only one of
the two that may be compared against a relaxation rate, a flux or a constant
hazard. Both are small positive numbers that rise and fall together, so the
confusion is invisible unless the width is on the row.

So it always is. `observables.DISPLAY_COLUMNS` is the complete list — `left`,
`right`, `width`, `at_risk`, `events`, `censored`, `exits`,
`conditional_failure`, `hazard_rate`, `survival`, `cumulative_hazard` — and
`RiskTable.rows()` emits all of it. The suite measures the difference the
division makes: sixty thousand exponential trials at rate 2 binned ten ways and
five ways recover their own known conditional probability
`1 − exp(−λw)` and time rate `(1 − exp(−λw))/w`; between the two binnings the
**probability changes by 1.82×** while the **rate changes by 0.91×**. Reporting
one as the other would misstate the answer by eighty percent, purely from a
choice of bin.

Three conventions are frozen as module constants rather than chosen per call:
bins are half-open on the left (`left < t ≤ right`); a trial that never commits
is administratively right-censored at the manifest's `window_close` and is never
dropped; and several clocks completing on one endpoint is **one** channel event
with every co-completer retained, so the winner counts sum to co-completers
rather than to committed trials, and both totals are reported.

A risk table must span the run's **whole** observation window. A table that
began late would silently condition on surviving to its first edge and one that
ended early would silently discard the trials still at risk — and either is a
window chosen after seeing the data. Both are refused, and `uniform_edges`
reproduces the window's endpoints exactly so that changing the resolution cannot
become changing the span.

### The exponent, and the refusal

The primary estimand is the plan's, unchanged:

```text
log(-log(1 - q_K)) = amplitude + exponent * log(K)
```

fitted as a binomial complementary-log-log likelihood **directly to the
committed and unresolved counts** of every frozen coupling cell, by Fisher
scoring with step halving. There is no cell filter anywhere in `analysis.py`;
every report states `cells_declared` and `cells_fitted` and they are always
equal.

That matters because the obvious procedure cannot represent the cells that carry
the most information about the slope. A cell with zero events has `q = 0` and no
first logarithm; a cell with every event has `q = 1` and no second one. The
package's own unweighted fitter, `analytic.fit_log_log_exponent`, *refuses* such
a point rather than dropping it, and the suite exercises that refusal so the
naive route cannot be taken by accident. The likelihood keeps both: a
zero-event cell contributes `−n·H` and an all-event cell contributes
`n·log q`, and a steep synthetic sweep with two zero-event and one all-event
cell fits all eight of them. Deleting the boundary cells and refitting gives a
*different* exponent, so the retention is load-bearing rather than bookkeeping.

`no_valid_exponent` is a result, not an exception. It comes back with the counts,
the estimate and a tuple of reasons from a closed vocabulary — `saturated`,
`separation`, `flat`, `ill_conditioned`, `non_convergence`, `lack_of_fit`,
`curvature`, `curvature_unidentified`, `uninformative`, `inadequate_resampling`
— and every threshold behind those reasons is a module constant, so a failing
fit cannot be made to pass by widening a budget at the call site. The suite
fires each one: every cell zero-event, every cell all-event, perfect
separation, a logistic response, a response curved in log-coupling, four trials
per cell, and a coupling grid identical to fifteen digits. Beside them an
ordinary power law of the same shape returns a valid exponent, so the refusals
are discriminating rather than free.

`flat` means a flat *likelihood*, not a flat *response*. A sweep whose
commitment probability genuinely does not vary with coupling, and which is
adequately informed, returns a **valid exponent near zero** — and it should: the
declared model contains `p = 0`, and refusing it would erase a real negative
scaling result rather than report one. An independent review reached the same
disposition.

`curvature_unidentified` closes a defect that review found. The curvature
extension is a *required* diagnostic, so it has to be estimable before the power
model it tests can be called valid. On a six-point grid spanning `1.000` to
`1.010` with two hundred thousand trials per cell and exactly half events, the
power exponent is adequately informed and the curvature information is singular:
the likelihood ratio is zero, its p-value is one, and the earlier implementation
read that as the curvature test *passing*. A test that cannot fail is not a test
that passed. The narrow grid is now refused and a wide grid of the same size
still returns its exponent.

### Common random numbers

A coupling sweep on one keyed noise hierarchy shares master trials: trial 7 at
`K = 0.4` and trial 7 at `K = 0.8` are the same Brownian path at two couplings.
Those cells are correlated, and a bootstrap that resamples inside each cell
independently treats them as though they were not.

The resampling unit is therefore a declared property of the design and is
**checked**: a `master_trial` design must hold identical trial identities in
every cell, an `independent_namespace` design must hold pairwise disjoint ones
under distinct namespaces, and a design that claims one and looks like the other
is refused. The wrong procedure exists under a name that says so —
`cell_independent_control` — because a control that does not exist cannot
discriminate against anything. Over six paired sweeps of 250 replicates each,
built from one shared latent uniform per trial, it returned a **smaller**
exponent spread every time.

**The paired record is the authority, not the Wald interval.** An independent
review measured the three on one 600-trial common-random-number sweep:

```text
expected-information Wald stderr       0.06112
paired master-trial bootstrap stderr   0.07932
wrong independent-cell stderr          0.06219
```

The marginal Wald interval that `fit_exponent` reports tracks the procedure this
package implements *as a control for being wrong*. So it may be displayed and it
may not decide anything. `analysis.calibrate_exponent` returns a
`CalibratedExponent` whose interval, standard error and curvature diagnostic all
come from the paired resampling record, and that is the only object `compare.py`
will accept.

**A replicate counts when its exponent is estimable, and refusals are carried by
reason.** The earlier implementation appended any finite exponent, so a
four-trial-per-cell design whose every replicate was refused reported `50 / 50`
valid with a percentile interval beside it. That reproduction now reports zero.
Adequacy has a frozen floor — `MIN_RESAMPLE_REPLICATES` draws and
`MIN_VALID_RESAMPLE_FRACTION` of them estimable — and an inadequate record is
`no_valid_exponent` under `inadequate_resampling` rather than a narrower
interval computed from the survivors.

*Estimable*, not *valid*, and the distinction is deliberate.
`ESTIMABILITY_REASONS` are the six that mean the estimate does not exist.
`GOODNESS_REASONS` — `lack_of_fit` and `curvature` — are model-adequacy verdicts
from the marginal deviance and likelihood ratio, and those are exactly the two
statistics common random numbers invalidate. A paired bootstrap resamples whole
master trials, so its replicates are clustered and both tests reject far above
their nominal level: at two thousand trials a cell, a perfectly specified power
law has sixty percent of its replicates rejected for lack of fit. Counting those
as failures would make every genuine CRN design "inadequate". So they are
counted and reported but do not disqualify a replicate — and both are replaced
by paired statistics for the design's own verdict.

### Model adequacy, under the pairing

Curvature is decided by the **same one-dimensional cluster score test**, with
the quadratic column W-residualized against `[1, log K]`.

It used to be the percentile interval of bootstrap curvature estimates. A review
held one dataset fixed, changed only the public analysis seed, and got "excludes
zero" for **58 of 100 seeds**; on interior exact nulls that interval excluded
zero 6.7 % of the time at a nominal 5. A percentile interval of point estimates
is not a test, and it was carrying an authority it could not support. The
interval survives on `ResamplingReport` as a **descriptive** figure; the decision
is the test.

General lack of fit is decided by `analysis.cluster_lack_of_fit`: a
**nuisance-adjusted cluster-robust score test** on the `cells − 2` lack-of-fit
contrasts, with a **Rademacher multiplier bootstrap** p-value that re-signs
whole master trials.

Write `s_t` for master trial *t*'s score contribution to the saturated model at
the fitted null, `X = [1, log K]`, and `W = diag(n_i H_i g(H_i))` for the
expected information in `eta`. Let `C` be an orthonormal basis of the Euclidean
complement of `X`. The maximum-likelihood conditions make the **total** score
`S = Σ s_t` exactly orthogonal to `X` — but an individual *cluster's*
contribution is not, and it is the clusters that form the sandwich meat and that
the multiplier re-signs. A review put it exactly:

> The total score is orthogonal to the nuisance design at the restricted MLE,
> but each *cluster contribution* is not. … Re-signing raw `C' s_t` therefore
> does not reproduce the score's distribution after fitting the two nuisance
> parameters.

So the contribution is the efficient score adjusted for the estimated nuisance,

```text
e_t = C' s_t − C' W X (X' W X)⁻¹ X' s_t          (equivalently  e_tᵀ = s_tᵀ B)
B   = C − X (X' W X)⁻¹ X' W C,        with  X' W B = 0  exactly
```

and the statistic is `(Σ e_t)' (Σ e_t e_t')⁻¹ (Σ e_t)`. The **numerator is
unchanged** — `Σ e_t = C'S` because `X'S = 0` at the maximum — so only the meat
and the bootstrap draws move. `verify.py` checks `X'WB = 0` to machine precision
and checks that the unadjusted basis violates it, so the adjustment is
demonstrably doing work. A design with independent namespaces has one
observation per cluster, so the same formula reduces to the ordinary
independent score test — one test covers both units.

**A boundary cell has no calibration.** A cell in which every trial did the same
thing contributes a *constant* cluster score: a systematic offset with no
sampling variation, which the sandwich divides by a variance the bootstrap
invents. The review's counterexample was an exactly specified power law with
three saturated cells, rejected **43 times in 120**. Such a design is now
`adequate=False` with the reason stated and `calibrate_exponent` returns
`uncalibrated_fit` — with **every cell still in the likelihood and the report**.
A calibration refusal is not cell deletion.

Size and power, over ensembles frozen after the method and before the table,
across four dependence structures generated by a Gaussian copula rather than the
monotone shared-uniform special case:

| ensemble (300 nulls each) | lack of fit | curvature | marginal deviance |
| --- | --- | --- | --- |
| null interior, independent cells | 7.0 % | 4.0 % | 6.7 % |
| null interior, copula ρ = 0.0 | 5.7 % | 7.0 % | 6.0 % |
| null interior, copula ρ = 0.4 | 4.0 % | 4.3 % | 1.3 % |
| null interior, copula ρ = 0.8 | 5.0 % | 3.7 % | 0.7 % |
| curved, `bend = 0.5`, ρ = 0.8 | 97 % | 99 % | 97 % |
| curved, `bend = 1.2`, ρ = 0.8 | 100 % | 100 % | 100 % |
| null saturated, ρ = 0.8 | *both decline, 120/120* | | — |

The nominal level is 5 %. The marginal deviance grows monotonically more
conservative as the dependence strengthens — 0.7 % at ρ = 0.8, seven times too
conservative, and a test that under-rejects *misses* non-power responses — while
both adjusted cluster tests hold their level.

The null ensembles are judged by the **exact two-sided binomial acceptance
region** at a guard level of 1e−3, not by a hand-picked percentage band. A
review pointed out that the previous 0.5–13 % guard would have accepted
`19/150 = 12.67 %`, whose one-sided probability under a 5 % null is about
1.9e−4; the suite now asserts that the region for 150 datasets excludes 19 and
admits an ordinary 8.

Nothing about either test is tunable: neither takes a seed or a replicate count,
the count is frozen at `LACK_OF_FIT_REPLICATES = 499`, and the seed is derived
from the design digest, so the same data always give the same p-value. Fewer
than `MIN_LACK_OF_FIT_CLUSTERS` clusters, or any cell without within-cell
variation, is `uncalibrated_fit` / `uncalibrated_curvature` — a *no-result*, and
in particular **no positive curvature claim is made from a design whose
curvature calibration is inadequate**. The marginal deviance remains on the
report as a display quantity.

**A trial number is not an identity.** `analysis.MasterTrial` carries a
namespace as well as a number, with a length-prefixed injective key, and the
namespace of a real design is `compare.crn_identity` — the digest of the
manifest fields that decide which Brownian path trial *i* is. It deliberately
excludes the peak coupling (a sweep varies it) and the model label (the full and
width-only runs share one tree by construction), so a sweep pairs across its
couplings and a full/control pair pairs with each other, while two unrelated
datasets both numbered `0..n-1` are refused. Reordered, duplicated, mixed and
substituted identities are each refused with a valid control beside them.

### Four comparators, four different quantities

| kind | what it is | finite population? |
| --- | --- | --- |
| `eligible_count` | how many of the exact production clocks are inside the tongue | yes |
| `bare_rate_sum` | the sum of their local relaxation rates — the **primary** comparator | yes |
| `continuum_flux` | that sum times the grid spacing: a midpoint quadrature | **no** |
| `constant_hazard` | a declared constant rate, retained as a falsifiable control | control |

Three of the four rise and fall together over a pulse and mean entirely
different things. `ComparatorCurve` derives its unit and its
finite-population flag from its `kind` and refuses to be built with either one
wrong, so handing the continuum quadrature to a caller expecting the
finite-population predictor is a construction error rather than a reading
mistake. `ComparatorSet` requires exactly one curve of each kind on one time
axis, so a set cannot quietly omit the continuum control.

Every quantity is computed on the manifest's **exact production detunings** — not
a refined grid, not a rebuilt one — so the staircase a finite population produces
at weak coupling is reported rather than smoothed away. Predictor curves may be
shifted by the manifest's own dwell duration, and the only alternative is no
shift: **no free time offset is fitted anywhere**, and `comparator_curves` has no
parameter through which one could arrive.

### The causal defeat rule

A quadratic population result supports the width-times-rate mechanism only if
the negative controls do not also look quadratic. `compare.causal_decision`
takes three `ControlEvidence` objects and a computed `HoldoutEvidence`, and
returns `blocked`, with named reasons, whenever:

- the full model has **no valid exponent**;
- **either negative control** has no valid exponent;
- the **central single-clock** control's paired interval still contains two;
- the **width-only** control's paired interval still contains two;
- the paired full-minus-control contrast is **inadequate**, or its lower
  endpoint does not exceed the frozen minimum mechanistic difference; or
- the holdout evidence is **incomplete**, or the control matches or outperforms
  the full model on the frozen out-of-sample rising/falling survival and
  winner-mismatch comparison.

The second and last of those closed defects an independent review found. The
earlier rule read only the controls' intervals, on the premise that an
inadequate control has an unbounded one — which is false for `lack_of_fit` and
`curvature`, both of which keep a finite and often narrow interval. The review
produced `not_blocked` with *both* controls reporting `no_valid_exponent`, one
of them at exponent 2.99 with an interval of (2.92, 3.06). And the holdout used
to be a bare boolean, so "no holdout was computed" and "the control was worse"
were the same argument.

**Roles are facts about sweeps, not labels at a call site.**
`compare.control_evidence(role, sweep, replicates, seed)` re-hashes the sweep,
builds its design from its own ledgers, **computes** the calibration, and
requires the role to be true of it: `central` is the ONE zero-detuning clock
under the full model, `width_only` is the fixed-contraction model, `full` is the
multi-clock population. The suite offers each of the three real sweeps in each
of the three roles and only the diagonal is accepted. A review had passed "an
honest synthetic calibration with counts `[(54,120), …]` as evidence for a
closed sweep holding `[(2,6), …]`"; there is no argument through which a
calibration can arrive any more.

**The holdout is computed inside the decision path.** `compare.holdout_evidence`
takes three opened runs, re-hashes them, requires them to agree on every frozen
field outside `EXPERIMENT_VARYING_FIELDS`, requires the two model roles, and
scores both candidates against the same declared reference on `frozen_edges` — a
binning that puts the pulse centre exactly on a boundary. There is no
caller-writable `complete` bit: a computed record and a `missing_holdout(reason)`
are two different `source` values, and a missing one blocks. A table whose bin
straddles the centre is refused rather than assigned to one side, and two
differently binned splits cannot be differenced at all: an equal number of rows
is not an equal binning.

**Everything is joined before anything is decided.** `causal_decision` requires
the contrast's `full_digest` and `control_digest` to be the full and width-only
evidence's own design digests, **its two ordered fingerprint tuples to be
theirs exactly**, its coupling grid and CRN identity to be theirs, and the
holdout's reference, full and control fingerprints to **belong to** the two
sweeps. Digests alone were not enough: deterministic rows give a gate-rewritten
copy of a sweep the *same* design digest and different fingerprints, and a
review used one to substitute a contrast.

All three roles must also share one `experiment_family`. The central control
legitimately has a different population and therefore a different CRN token, so
equality of *that* token is the wrong check for it; the family digest is the
right one.

It is built by **subtraction**, and that direction is the fix. An earlier
revision derived it from `CRN_IDENTITY_FIELDS` — an allowlist assembled for a
different purpose — and a review found `drive` missing from it: a genuine
one-clock sweep whose manifest said `stationary` joined a `pulsed` causal
experiment with no provenance refusal at all. `EXPERIMENT_FAMILY_FIELDS` is now
the whole manifest minus `EXPERIMENT_FAMILY_EXCLUSIONS` — the four sweep-level
differences (`peak_coupling`, `pulse_digest`, `model`,
`fixed_contraction_rate`) and the seven population fields — so a field nobody
thought about is *included* by default. The suite mutates **every one of the 45
manifest fields** in turn and requires each either to change the family or to be
on that reviewed exclusion list. Each control also carries its sweep's status
blockers, so a diagnostic-only experiment blocks under
`numerical_status_not_cleared`. A review combined public constructors into three
unrelated controls, a caller-authored contrast and a holdout belonging to
neither, and got `not_blocked ()` over calibrations that all said
`no_valid_exponent`; that combination now blocks, and the records it was built
from cannot be constructed at all.

`compare.paired_exponent_contrast(full_sweep, control_sweep, …)` likewise takes
sweeps: it re-hashes both, requires them to be the two models on one
common-random-number hierarchy and one coupling grid, and carries the digests,
grid, identity and fingerprints forward. A replicate counts only when both
exponents are estimable. The out-of-sample score is an unweighted sum of three
non-negative discrepancies and a tie counts as a defeat.

The other verdict is `not_blocked`. It is the absence of one named defeater and
nothing more; there is no verdict spelled `supported`, and the suite asserts
that there is not.

### The gate

`compare.scaling_verdict(sweep)` takes a `RunSweep` — **not a tuple of
strings, not a calibration, and not a seed.** Three independent reviews
attacked this door in turn:

```python
scaling_verdict(good_fit, ("production_cleared",))   # -> permitted, no ledger
forged = CalibratedExponent(..., "valid_exponent", (), 1.0, ...)
scaling_verdict(forged, cleared_sweep)               # -> permitted
ScalingVerdict("physical_verdict_permitted", "closed_runs", (),
               ("production_cleared",), ("not-a-fingerprint",))   # -> accepted
```

All three now refuse. The gate string is a `TypeError`; the calibration is no
longer an argument at all, because `scaling_verdict` **recomputes it** from the
re-opened sweep's own ledgers; and every derived record in this layer is
factory-only.

**No authoritative door takes a seed or a replicate count.** A third review
held one dataset fixed and changed only the public bootstrap seed, and got
"curvature excludes zero" for 58 of 100 seeds. A decision that can be rerolled
is a decision that can be selected, so `calibrate_exponent(design)`,
`scaling_verdict(sweep)`, `control_evidence(role, sweep)` and
`paired_exponent_contrast(full, control)` take neither: the count is frozen at
`RESAMPLE_REPLICATES = 499` and the seed is derived from the design digest.
`resample_exponent` keeps both arguments and is the clearly-exploratory door —
it is what the wrong-procedure control is measured with, and it feeds no
decision.

**A `RunSweep` is one experiment.** It can only be built from `ReadOnlyRun`
objects, which are sealed and can only come from `open_for_comparison`.
Construction requires strictly increasing distinct couplings, distinct
fingerprints, **one model, one fixed contraction rate**, and every frozen
manifest field outside `SWEEP_VARYING_FIELDS` identical across the runs. An
alternating full / width-only mixture is refused rather than fitted.

**A handle is a locator; the disk is the authority.** A review found the gap
between re-hashing and re-*reading*: `reverify()` re-opened the run, compared
fingerprints and returned the fingerprint — and the caller went on computing
from the in-memory ledger it already held. Because the fingerprint derives from
the close *marker*, a `ClosedLedger` carrying the genuine marker beside forged
rows and a forged manifest passed intact, and a fit reported counts
`(6/6, 6/6, 6/6, 6/6)` over a disk holding `(2/6, 4/6, 5/6, 6/6)`.

So `ReadOnlyRun.reopen()` and `RunSweep.reopen()` return the **fresh** object,
parsed from the bytes on disk, and every authoritative factory —
`scaling_verdict`, `control_evidence`, `paired_exponent_contrast`,
`holdout_evidence` — calls one of them immediately before building a design and
computes only from what it returns. The public `design_from_sweep` reopens too.
The supplied handle's name says which run, its fingerprint is compared, and its
data are discarded.

**A parsed row is a mutable dictionary, so a verified snapshot goes through the
gate again.** `ClosedLedger` is a frozen dataclass, but frozen freezes the
*reference*: a review edited one parsed row of a verified ledger — no file
touched, no private name used — and watched the public design factory return a
different digest.

The first answer was a digest check: re-encode the rows, re-emit the manifest,
compare against the marker the object carries. A further review showed that is
not enough. A caller can edit the rows, re-encode all three tables with
`encode_table`, rebuild a `CloseMarker` over the new digests and construct a new
`ClosedLedger` — which then agrees with *its own* marker perfectly. Four
all-public variants were accepted, and the observables answered inconsistently:
a tally reporting five commitments beside a winner table holding four.

So `raw_ledger.require_unmutated_ledger` now reconstructs the canonical manifest
and table bytes and puts them, with the snapshot's marker, through the whole of
`require_closed_ledger` — the same eight checks a run read from disk passes,
including row counts, row identities, ordering and cross-table reconciliation.
**It returns the object that gate parsed**, and every caller binds it: the three
observables, the sweep's own counts, and the design builder. There is no
check-then-use gap, because the thing checked is the thing used.

It needs no disk access, so it protects the *display* routes — which have a
ledger and no run name — as well as the decision ones. The reopen catches a
changed file; this catches a changed object; neither substitutes for the other.

**The design is bound to the observations, not the totals.**
`ScalingDesign.digest` covers the ordered couplings, the canonical master-trial
identities *and every outcome*. Two different outcome matrices with the same
cell totals have different paired bootstraps, so comparing totals was not
enough; the fit, the resampling record and the lack-of-fit test all carry the
same digest and `CalibratedExponent` refuses to hold three records about three
different designs.

**Schema v3 fails closed.** A review pointed out that this suite's own positive
control rewrote `numerical_gate` to `production_cleared` and re-hashed — while
the *same hashed manifest* still said `moving_band_status: numerical_no_result`
and `validation_note: … no scaling claim may be drawn from these tables`. A
free-text field rewritten is not a production budget. So all four status fields
are read jointly, and above them sits `PRODUCTION_STATUS_SCHEMA = 4`, a manifest
schema that **does not exist**: `ScalingVerdict` refuses to hold
`physical_verdict_permitted` under v3 at all. A future permission needs a schema
revision carrying a closed-vocabulary production-budget status, not a rewritten
sentence.

The gate still reads its input, which is what stops it being a constant
function: rewriting one run's `numerical_gate` removes exactly that blocker from
the sweep's status list and leaves the schema, moving-band and validation-note
blockers standing. The suite asserts precisely that set difference.

Exercising the estimator on data that never came from a run is a different
activity with a different **type**: `compare.synthetic_machinery_verdict`
returns a `MachineryReport`, which has no verdict field, answers `False` to
`is_physical`, and carries its own disclaimer in `summary`.

---

## Feasibility, power and the pilot firewall

`experiments.py` answers three questions before any production run is
scheduled: **what does it cost**, **how many trials does the frozen precision
target need**, and **is the numerical error small enough against that precision
to mean anything**. It sits outside the raw import graph in both directions.

It runs no pilot and signs nothing. This section is the closure round of an
independent review; the numbers below replace an earlier set that the review
correctly found under-measured.

### What was measured

The benchmark calls `raw_runner.write_raw_run` at the declared block sizes,
with a discarded warmup and **two** timed repeats each, all recorded. The
`intended` fixture runs at the matrix's **own 64 clocks and its own
`1.953125e-3` timestep**, its pulse shortened so the workload stays bounded;
it commits nothing, so its first pass does not truncate and its rate is the
honest one.

The work unit is a **walked clock-step**: `2 × trials × clocks × steps` for the
two censoring passes, plus `shadow_trials × clocks × steps` for the rewalk.
Quoting per *nominal* step and then multiplying by a shadow factor was double
counting; this unit normalizes once.

| Fixture | Clocks | Steps | Walked | Repeats (s) | Rate | Peak |
| --- | --- | --- | --- | --- | --- | --- |
| intended | 64 | 256 | 49 152 | 11.63, 11.62 | **4 228/s** | 173 kB |
| base | 16 | 80 | 16 640 | 3.51, 3.51 | 4 742/s | 104 kB |
| trials ×4 | 16 | 80 | 62 720 | 12.69, 12.88 | 4 870/s | 114 kB |
| steps ×4 | 16 | 320 | 66 560 | 12.19, 12.17 | 5 460/s | 99 kB |
| shadow ×6 | 16 | 80 | 23 040 | 4.97, 4.98 | 4 624/s | 86 kB |

Costing uses the **slowest**, which is the intended-resolution fixture at
**4 228 walked clock-steps/s**. Across fixtures four times apart in clock count
and twenty-five times apart in timestep the rates span a factor 1.29 — that
agreement is the independent spot check, run inside the suite. The shadow
rewalk's *incremental* cost is measured separately and lands at the same rate,
so charging it at the base rate is a measurement rather than an assumption.

Storage is **91–106 bytes per durable row** with the fixed per-run overhead
(manifest, marker, three headers ≈ 7.1–7.8 kB) separated out.

The no-cube claim rests on flatness: quadrupling trials or steps quadruples the
`trials × clocks × steps` float64 array the plan forbids, while the traced peak
moves by 1.10× and 0.96×. `tracemalloc` tracks NumPy buffers, so a materialized
array would appear in the number.

### How the trial counts were derived

Every arm declares a conservative **probability envelope** before any pilot and
is sized from the least informative point in it, with its own declared
invalid-fit price. Sizing one curve and applying it to all three arms was the
review's finding; a one-clock central control is not as informative as a
sixty-four-clock population.

| Arm | Clocks | Envelope | Weight | Trials/cell |
| --- | --- | --- | --- | --- |
| Full Adler | 64 | 0.10 – 0.95 | 0.0999 | **2 406** |
| Width-only control | 64 | 0.05 – 0.95 | 0.0500 | **2 406** |
| Central-clock control | 1 | 0.02 – 0.95 | 0.0200 | **6 346** |

The two paired arms share a count because their contrast is paired; the central
control does not pair and keeps its own. The analytic independent-cell
requirement alone is **458**; the rest is the declared 2.0 pairing inflation,
the per-arm failure inflation, and a declared 1.25 conservative factor — three
separately reported numbers, never one.

`cloglog_slope_variance` is checked against a simulation of itself: over 150
replicates the closed form gives 0.12525 and the measured independent spread is
0.11894 (ratio 0.95), across 142/150 **valid** fits. Refused fits are excluded
rather than averaged in.

The recovery ensemble draws its shared master-trial latent from the package's
own counter-keyed physical stream and mixes it with a per-cell component at
correlation 0.85, so realizations are **not** monotone across coupling.
Measured pairing inflation **1.29** against the declared 2.0. A control arm is
generated beside the full one at a declared **negative** between-arm
correlation, and the contrast's measured variance inflation is 0.70 against the
declared 1.25 — the earlier claim that pairing "can only reduce" the contrast
variance was false and is gone.

A frozen mechanical **post-range join** re-sizes the design from a count-only
selection and the counts it came from, and nothing else; a selection that did
not come from those counts is refused.

### The numerical budget, and why it blocks

A measured discrepancy, inflated by **the budget's own frozen coverage sigma**,
must sit under one quarter of the planned production 95 % half-width, with one
intended timestep added for a commit-time quantile. Seventeen authoritative
ticket-04 rows are carried — every observable at every position at each
ladder's finest level — and the limiting row is chosen by the full bound inside
its own unit.

| Limiting row | Unit | Bound | Allowance |
| --- | --- | --- | --- |
| moving-band commit-time p20 shift | time | 0.17700 | 0.02195 |
| stationary exit quantile p35 | time | 0.15290 | 0.02195 |
| stationary survival | probability | 0.04125 | 0.00499 |

Disposition: **`numerical_no_result`**, with six blockers. The moving-band
ladder's own verdict is carried through unchanged. Two results are reported
**separately**, because they behave differently:

* the **probability** window caps the trial count at **35 per cell** against the
  2 406 the target needs; and
* both **time** rows fail at *every* trial count, because the time allowance
  does not depend on `n`. Reading 35 as an overall maximum was the misreading
  the review flagged; the overall admissible count is **0**.

### What has actually been ruled out

A finite declared search over **12 options** — refinements to a sixteenth of
the timestep, sample enlargements to 1024×, reduced trial counts, their
combinations, and one option that additionally measures at the intended
configuration — projects each evidence row with **bias and sampling error
scaled by their own levers**. Scaling the whole bound by `√dt` made refinement
look like an escape when it is not.

None of the twelve fits. The closest misses by **1.33×**. The verdict is
therefore `no_feasible_matrix_among_searched`, and it names that set — not
every conceivable design. A claim that no timestep, sample or configuration
could work needs a search or a bound this implementation has not performed.

### The proposed matrix, and what it would cost

Finite, with an **exact execution graph**: the ladder refines named couplings
on a named arm, and each arm carries its own grid.

| | |
| --- | --- |
| Coupling nodes | 6, log-spaced, 0.5 … 2.0 |
| Grids | full and width-only 64 clocks on `[-3, 3]`, spacing 0.09375, even parity; central control **1 clock** |
| Eligible at the weakest coupling | 10, counted on the exact grid, minimum 8 |
| Timestep / ladder | `1.953125e-3`, halved twice |
| Refinement | couplings 0.6598 and 1.1487, `full` arm, 200 trials/cell |
| Shadow | 5 %, rounded up **per run** |

| Line | Runs | Physical trials | Shadow rewalks | Rows | Hours |
| --- | --- | --- | --- | --- | --- |
| primary/full | 6 | 14 436 | 726 | 984 804 | 382 |
| primary/central_control | 6 | 38 076 | 1 908 | 78 060 | 16 |
| primary/width_only_control | 6 | 14 436 | 726 | 984 804 | 382 |
| refinement/level1/full | 2 | 400 | 20 | 27 280 | 21 |
| refinement/level2/full | 2 | 400 | 20 | 27 280 | 42 |
| **total** | **22** | **67 748** | **3 400** | **2 102 228** | **844** |

≈ **794 MB** under the intended-configuration row bound. The range pilot is
priced as a **separate** line and is not inside that total. Physical trial executions and shadow rewalks are counted apart:
a rewalk is validation work over an existing master trial, not a new trial.

**Every wall-time figure on this page is a machine-dependent snapshot**, not a
property of the experiment. Run counts, trial counts, row counts and storage
are reproducible; hours are a measurement of whichever core ran the benchmark,
and an independent review measuring 908 h against an unlabelled 844 h in this
table was right to call that out. The hours above were taken on the environment
recorded with the benchmark digest and proposal design digest in the execution
notes for the round that produced them, and the live numbers are
**re-derived by the verifier on the machine that runs it** — `--verbose` prints the current
hours beside the current digests. Quote the row, storage and run counts; take
the hours from your own run.

### The pilot firewall

Five separations. The two the review broke are closed at the level it asked
for:

* the count door validates the **run name** against the pilot and production
  directory families, so a pilot-labelled, pilot-keyed run in a
  `prod-…` directory is refused; and
* a pilot run is **written directly into a quarantine directory**
  `results/pilot_quarantine/<name>` that `raw_runner.open_raw_run` cannot
  address — its name is two segments, which that door refuses — so the
  unrestricted ledger reader is no longer an alternate path to forbidden pilot
  data. There is no write-then-move: the production results root never holds a
  pilot run, not for an instant. Quarantine costs no auditability: the count
  projection re-reads the five files through the same eight-check closure gate,
  and a mutated quarantined table is refused.

A third separation closes what the first two could not. The writer
authenticated *where* a run was written and not *which pilot was approved*, so
a width-only two-clock `dt=1` two-trial diagnostic could be written to the
right place and then projected into authoritative range evidence. A frozen
**`PilotPlan`** now names the role, the model, the ordered candidate couplings,
the whole 64-clock grid, the intended timestep and the window it implies, the
pulse, noise, band and dwell, 200 trials per cell and its shadow count, the
execution blocking, the reserved stream and run-name family, the count-only
selection rule and its three thresholds, and the source and schema versions it
was frozen against. `write_pilot_run(plan, firewall, coupling, run_name)` has
**no configuration parameter and no model parameter**; every check runs before
any directory is created, so a refusal leaves no residue. The plan is
re-checked at the count door as well: every quarantined run's whole 45-key raw
manifest is rebuilt from the plan and compared field for field, so a directory
placed in the quarantine by any other route is refused there too. The proposal
binds `pilot_plan_digest` inside its design freeze.

`PilotCellCounts` carries exactly coupling, trials, committed and survivors —
four values and no provenance field. `reference_rate`, `select_range`,
`post_range_power` and the manifest factory are **authoritative doors**: each
takes the frozen plan, the firewall it was frozen against and the ordered
quarantined run names, re-opens and re-counts every run through the closure
gate and the plan's manifest rebuild in one snapshot, and recomputes the
selection and the rate from it. There is no selection parameter and no rate
parameter, so a caller-built range cannot authorize its own contraction rate.

### The proposal is not an approval

`ProposedProductionManifest` is **factory-only** and validates all **86**
declared fields by type, enum, range and cross-field relationship. The review's
complete `"wrong-but-json-valid"` record is refused twice over, unsealed and
sealed; every single field is separately refused when given a wrong value with
the seal present; and `PROPOSAL_STATUSES` has exactly one member, so
`checkpoint_one_approved` is not a string this module can spell.

The freeze split in two: **`design_digest`** covers 81 of the 86 declared
fields — everything except the five declared exclusions — and no wall-time
measurement, so re-running the same source reproduces it, while
the whole-proposal digest additionally binds the benchmark evidence and moves
with the machine — which is correct for a record of what a machine did and
wrong for a scientific freeze.

Current design digest
recomputed each round and recorded in the execution notes, status
`checkpoint_blocked`, feasibility `no_feasible_matrix_among_searched`.

**There is no post-range production proposal, and ticket 08 matrix derivation
has not begun.** A pilot's counts would have to be re-derived into a post-range
sampling target, primary matrix and resource estimate before a proposal could
carry them; those records are not built here, and the proposal factory refuses
by name — `post_range_matrix_not_derived` — rather than serializing a
post-range power estimate beside the pre-pilot matrix. The checkpoint states
this as two derived fields, `ticket_08_derivation_status = not_started` and
`ticket_08_authorization_status = not_authorized`, each with a one-member
vocabulary and neither with a factory parameter. They are different facts:
work can be unauthorized and also not begun, and neither implies the other, so
the summary carries both rather than leaving a reader to infer one from a
sentence. Both appear in `Checkpoint1Handoff.summary()` beside
`post_range_production_proposal = None`.

The checkpoint's presented decision text is derived too, from one canonical
constant with no parameter a caller could supply, and it states all three
invariants explicitly: **There is no post-range production proposal. Ticket 08
matrix derivation has not begun. Ticket 08 is not authorized.** The last two
state what exists and what has not been done; neither is a request to authorize
ticket 08.

### One blocked option and one machinery-only diagnostic

Neither is "run the pilot", and they are **not two costed alternatives**:
Option A carries no price at all.

| Option | What it buys | Cost |
| --- | --- | --- |
| **A — intended-configuration validation outline** | a seven-stage outline refining to dt/256 with space refined separately from time, and moving stages separating timestep bias from sampling error | **unpriced** — the two timed kernels omit the moving replay and the stationary endpoint generation, so no range built on them is an estimate or a bound |
| **B — 18-run primary-only diagnostic** | the six-cell three-arm matrix at 35 trials per cell, the largest count the present probability bound admits; **excludes** the four refinement runs, the range pilot and all numerical validation | **12.8 MB** over exactly 18 runs and 29 268 rows; its hours are a machine-dependent snapshot, re-derived per run |

Option B supports a **machinery and diagnostic statement only**: at 35 trials
per cell the exponent interval is far wider than the frozen target, both time
bounds still fail at every trial count, and the moving-band
`numerical_no_result` still blocks. It delimits what the current evidence can
carry; it does not restore a scaling claim.

## Not implemented here

A first clock is a model stopping event, not an outcome. There is a population,
a first-winner stop, a ledger on disk, and now a survival curve, a hazard and a
fitted exponent computed from it — and no physical scaling claim of any kind is
made here, because every run records `numerical_gate = "diagnostic_only"` and
every scaling verdict this package can reach is therefore `machinery_only`.

That is worth stating plainly, because the package now writes a directory of
CSV files with a word like *commitment* in the column headers and can fit an
exponent to them, and a reader could reasonably assume there is a result in it.
There is not. What exists is a record of which clock stopped which trial and
what every other clock was doing at the time, and an estimator that has been
shown to behave correctly on synthetic data with known answers.
`raw_runner.py` carries the list in `UNIMPLEMENTED_RAW_STAGES`, and
it is the same list as this one — with one distinction. The first four rows
below are implemented, and are implemented *somewhere `raw_runner.py` cannot
reach*. From the raw boundary's own point of view they are still absent, which
is the point: an oracle a run could import is not an oracle, and an estimator a
run could import is not an independent estimator.

| Stage | Status |
| --- | --- |
| The stationary killed-diffusion oracle and the endpoint convergence gate | implemented **outside** the raw graph (`killed_diffusion.py`) |
| The diagnostic moving-band Brownian-bridge audit | implemented **outside** the raw graph (`moving_band_audit.py`) |
| Survival, hazard, and the binomial complementary-log-log estimator | implemented **outside** the raw graph (`observables.py`, `analysis.py`) |
| Every isolated comparison curve: eligible count, bare relaxation-rate sum, continuum control, constant-rate control | implemented **outside** the raw graph (`compare.py`) |
| The throughput and memory benchmark, the power calculation, and the frozen production run matrix | implemented **outside** the raw graph (`experiments.py`); the matrix is a **proposal**, not a schedule |
| The pilot firewall, the pilot quarantine and the count-only range rule | implemented **outside** the raw graph (`experiments.py`); **no pilot has been run** |
| The signed production manifest and the production coupling sweep | not implemented — a proposal exists at status `checkpoint_blocked`, and signing needs two named human checkpoints |
| The ticket-08 production gate | **not implemented**: `production.py` does not exist, no registered check covers it, and ticket 08 has not begun; **no production, pilot, campaign or sensitivity run exists** |
| The production numerical budget — nothing written here has passed one | frozen (`experiments.py`), evaluated, and **not met**: the disposition is `numerical_no_result` |
| Validated raw configuration boundary | implemented (`raw_config.py`) |
| Closed-ledger manifest, table schema, marker and digest protocol | implemented (`raw_ledger.py`) |
| Counter-keyed streamed phase kicks; the nested split at exact crossings; the keyed prepared ensemble | implemented (`stochastic.py`) |
| Euler–Maruyama one-clock stepping; exact tongue entry/exit/re-entry handoffs | implemented (`dynamics.py`) |
| The fixed-contraction width-only control, on one shared Brownian tree | implemented (`dynamics.py`) |
| The contraction-guarded lock predicate and the timestamp dwell clock | implemented (`commitment.py`) |
| Experiments S1, S2 and S4 as in-memory one-clock records | implemented (`raw_experiments.py`) |
| The finite-population race, first-completion stopping, retained co-completers, censored exposure, the no-stop shadow sample | implemented (`raw_race.py`) |
| The immutable hashed run on disk and the gate a consumer must pass to read one | implemented (`raw_runner.py`, `raw_ledger.py`) |
| Transitive structural and runtime isolation of the raw import graph | implemented (`verify.py`, checks 21–26) |

Two stages above are **gates**, not merely missing work, and both now exist. The
stationary killed-diffusion oracle measures how much uninterrupted dwell endpoint
accounting overestimates, and the moving-band audit measures what that costs a
dwell during a pulse. Both are held to tolerances frozen and hashed before their
results are opened, and failure of either is a *numerical no-result* — it blocks
interpretation and falsifies nothing. Neither gate has yet been run against a
*production* budget: the budgets exercised here are declared reference values,
and the production numbers come from a later ticket's power calculation.

Also deliberately absent, and belonging to later stages of the plan: tie
handling, two-channel competition, the polarization-angle sweep,
Gaussian/Lorentzian/structured spectral controls, and the fixed-peak versus
fixed-area duration sweeps as *experiments* (the distinction itself is already
enforced above). No manuscript text has been changed.

A note on what the later stages must not do: the commitment dwell is a **fixed
physical time**, shared by every channel and every amplitude. An inverse-coupling
dwell would preserve Adler scale similarity and make a quadratic result much
easier to obtain — and would also make the selection criterion
amplitude-dependent, which is exactly the thing being tested. It may be run later
as a labelled positive control, never as the primary result.
