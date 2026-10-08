"""
adler_born_two_channel — deterministic Adler calibration, the analytic
flat-spectrum control, and the isolated raw-event boundary.

**Importing this package imports nothing else.**  There is not one import
statement below this docstring, and that absence is a contract rather than an
oversight.

Why the package root is empty
-----------------------------
The experiment being built asks whether noisy Adler clocks commit in a pattern
matching an independently computed analytic prediction.  That question is only
meaningful if the code generating the events cannot see the prediction.  An
earlier revision of this file said:

    from . import analytic, model, simulate, validation

which meant that importing *any* submodule — including a raw event module whose
whole purpose is not to know the answer — initialised the package with
``analytic`` already loaded and sitting in ``sys.modules``.  The structural
check parsed the dynamics and found nothing wrong, because nothing was wrong
with the dynamics: the leak was in the package root.

So the root now loads nothing, and ``verify.py`` runs a subprocess that imports
the raw entry point and asserts that ``adler_born_two_channel.analytic`` is
absent from ``sys.modules`` afterwards.

How to reach ``analytic`` on purpose
------------------------------------
Comparison and reporting work still needs it, and all three ordinary explicit
forms still work exactly as before::

    from adler_born_two_channel import analytic
    from adler_born_two_channel.analytic import tongue_rate_sum
    import adler_born_two_channel.analytic

Each of those is a statement that names ``analytic``, so each is visible to a
reader, to a reviewer, and to the syntax-tree scan that forbids the raw modules
from writing one.  The single form that no longer works is the implicit one::

    import adler_born_two_channel
    adler_born_two_channel.analytic          # AttributeError

That is deliberate and is the whole change: attribute access after a bare
package import used to succeed *because* of the eager import above, and it is
the one route to the prediction that leaves no trace in the importing file.

Note that ``from adler_born_two_channel import *`` does load every submodule
named in ``__all__``, ``analytic`` included.  A star import is an explicit act;
the raw modules are forbidden from writing one, since a relative star import is
rejected by the same scan.

What this is and is not
-----------------------
This is protection against *accidental manufacture* — against a prediction
becoming reachable from an event generator without anybody deciding that it
should.  It is not a security sandbox, and nothing here is written as though a
determined caller were an adversary; a determined caller has ``importlib``.

Stages implemented in this package
----------------------------------
1. stationary one-clock Adler calibration (stable phase, near-lock relaxation,
   critical slowing, outside-tongue phase slip);
2. raised-cosine photon-envelope calibration (tongue entry and exit, the moving
   stable phase, finite synchronization time);
3. the analytic finite-grid flat-spectrum relaxation-rate sum and its
   coupling-squared convergence;
4. the raw-event *isolation boundary* — validated configuration types, the
   closed-ledger protocol, and the checks that enforce the separation;
5. the counter-keyed physical white-noise source and its nested Brownian tree;
6. **one-clock** stochastic evolution: Euler-Maruyama over the exact elementary
   mesh, the exact tongue entry/exit/re-entry handoffs, the timestamp dwell
   clock with its contraction guard, and the fixed-contraction width-only
   control that shares the full model's noise; and
7. the **continuum audit** of that endpoint scheme: an independent stationary
   killed-diffusion oracle, a paired endpoint refinement ladder measured
   against it, frozen hash-checked numerical budgets with a machine-readable
   no-result, and a diagnostic moving-band Brownian-bridge reset replay; and
8. the **finite one-channel race**: a declared population of actual competing
   clocks bound to one immutable grid identity, a stop at the earliest dwell
   completion with every same-endpoint co-completer retained, one authoritative
   censored row per (trial, clock) with its exposure and opportunity time, a
   predeclared no-stop shadow sample, and the immutable hashed run on disk
   whose every channel-level number a consumer re-derives from those rows
   before believing any of it; and
9. the **separate read-only analysis layer**: risk sets, survival, the discrete
   conditional failure probability and the width-divided time-rate hazard, the
   unresolved fraction with its three causes and the winner mismatch with every
   co-completer retained; the binomial complementary-log-log exponent fitted
   directly to committed/unresolved counts in every frozen coupling cell, its
   curvature and lack-of-fit tests and its explicit ``no valid exponent``
   verdict; paired master-trial resampling beside the deliberately wrong
   independent-cell control; and the comparison process that opens a closed
   hash-matching run read-only, keeps the four comparators separately named and
   typed, and applies the causal defeat rule.

10. the **feasibility, power and freeze layer**: measured throughput, peak
   memory and storage for the real streamed race with a demonstration that no
   trial-by-clock-by-time noise cube is materialized; the design-based power
   calculation for the binomial complementary-log-log exponent, checked against
   a simulation of itself and priced for common-random-number pairing; the
   production numerical budget tied to that planned uncertainty, which ticket
   04's measured envelopes do **not** meet; one finite primary matrix proposal
   with its parity, staircase fallback, stop and failure rules and its scaled
   cost; the pilot firewall whose count-only door returns four numbers and has
   no field for anything else; and a dry-run production manifest **proposal**
   that is structurally unable to clear or sign itself.

Stages 7, 9 and 10 sit deliberately *outside* the raw event graph.
``killed_diffusion``, ``moving_band_audit``, ``observables``, ``analysis``,
``compare`` and ``experiments`` are absent from ``__all__`` below and are absent
from the transitive imports of every raw module: the oracle a run is measured
against must not be reachable from the code that produces the run, the audit may
only be read after the fact, an estimator a run could import is not an
independent estimator, and a cost model and a power calculation are things an
event generator has no business seeing.  All six are ordinary explicit imports
for anybody who wants them —

::

    from adler_born_two_channel import killed_diffusion, moving_band_audit
    from adler_born_two_channel import observables, analysis, compare
    from adler_born_two_channel import experiments

— and ``verify.py`` walks the raw import graph and asserts that no raw module
writes one, then starts a fresh interpreter and asserts that importing
``raw_runner`` loads none of them.

Stage 8's first clock is a **model stopping event**: the trial ends because a
one-channel first-passage study defines it to end there.  It is not a detector
click, an absorption, a measurement outcome, a unique actuality, a microscopic
bath, one side of a two-channel competition, or a Born-rule result.  Stage 9
now supplies a survival curve, a hazard and a fitted exponent, and each of them
is a statement about the *estimator*: every run written here records itself as
*diagnostic* rather than production-cleared, because the production numerical
budget ticket 07 froze is **not met** -- the disposition is
``numerical_no_result`` -- so
``compare.scaling_verdict`` returns ``machinery_only`` however good a fit is,
and no commitment law or physical scaling claim exists.  ``raw_runner.py``
names what is missing explicitly.  The stage-7 audit measures a
*discretization*: it makes no physical claim, and a failed numerical budget is
a no-result rather than a finding.  Nothing in this package demonstrates Born
probabilities.  See ``README.md``.
"""

__all__ = [
    # Deterministic calibration, and the prediction it will be tested against.
    "analytic",
    "model",
    "simulate",
    "validation",
    # The raw-event boundary.  Start at raw_runner.py.
    "commitment",
    "dynamics",
    "raw_config",
    "raw_experiments",
    "raw_ledger",
    "raw_race",
    "raw_runner",
    "stochastic",
]
__version__ = "0.9.0"
