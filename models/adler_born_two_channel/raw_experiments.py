#!/usr/bin/env python3
"""
raw_experiments.py — the one-clock stochastic experiments S1, S2 and S4, and
the labelled records they produce.

What this module is
-------------------
The place where the three parts of the raw side are put together: the keyed
noise of ``stochastic.py``, the stepper of ``dynamics.py``, and the dwell
machine of ``commitment.py``.  It walks one clock, at one mismatch, through one
pulse train, and returns a record of what happened.

Three experiments, and each answers a different kind of question:

**S1 — the deterministic limit.**  Set the diffusion to zero and require that
the stochastic walk becomes deterministic Euler, term for term.  This is not a
statistical statement: at zero diffusion every leaf kick is exactly ``0.0``, so
the assertion is bit equality against ``model.step_euler`` on the same
sequence of intervals, and the convergence claim is the ordinary first-order
one measured against a fourth-order reference.

**S2 — the synthetic dwell machine.**  Feed declared eligibility, error and
contraction sequences straight into :class:`~.commitment.LockCriterion` and
assert on the whole state history.  Nothing here integrates anything.  The
point is to visit edges that a trajectory would only reach by luck: dwell
shorter than a timestep, dwell exactly equal to one, a non-integral ratio, an
initial inside sample, an exit on the very step that would have committed, the
falling edge, the wrap boundary, the repelling fixed point, the expanding side
of the separatrix, and the stable/unstable coalescence at the tongue edge.

**S4 — one stochastic clock.**  Run a single clock under a pulse and record
what it did.  Central, interior and near-edge mismatches; zero, weak,
intermediate and strong noise; the full Adler model, the central-clock control
and the fixed-contraction width-only control.

What these records are not
--------------------------
They are **not** a population, not a race, not a first-winner ledger, not a
survival curve, not a hazard, not an exponent, and not evidence about Born's
rule.  There is exactly one clock in every record here, so there is nothing for
it to win against.  A commitment is one clock having stayed inside a declared
band around its moving target for a declared time — see
:data:`~.commitment.TRIAL_CATEGORIES` — and nothing in this file licenses
calling that a detector click, an absorption, or an outcome.

Nothing here is written to disk either.  These records live in memory for the
length of a check.  The immutable hashed ledger, the production manifest and
the run directory belong to later tickets, and building half of one now would
produce a file that looked like a result.

Non-monotonic behaviour in noise strength is expected and is *reported*, not
smoothed: weak noise can help a clock into the band while strong noise breaks
the dwell, and the record carries the entry and reset counts that distinguish
those two effects.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional, Tuple

from .commitment import (TRIAL_CATEGORIES, DwellState, _advance_trusted,
                         validate_lock_criterion)
from .dynamics import (DRIVE_LABELS, FIXED_CONTRACTION_RATE, MODEL_LABELS,
                       validate_clock_path)
from .stochastic import FinestMesh, PhaseNoiseStream
from .validation import (require_exact_type, require_finite,
                         require_positive_int)

__all__ = [
    "NON_CLAIMS",
    "OneClockRecord",
    "elementary_segments",
    "run_one_clock",
    "synthetic_dwell_history",
]

#: What a reader of these records is not entitled to conclude.  Pinned here as
#: text, and pinned again in the README, because the disclaimers are the part
#: of a result most likely to be dropped in a summary.
NON_CLAIMS = (
    "one clock is not a population: nothing here is a race, and no clock "
    "wins",
    "a commitment is a dwell in a declared band, not a detector click, an "
    "absorption, a measurement outcome or a Born event",
    "no analytic prediction, amplitude square or prescribed hazard reaches "
    "these dynamics or this state machine",
    "no scaling claim of any kind is made: these records are not fitted, and "
    "the coupling exponent is a later ticket's estimand",
    "the white noise is a minimal effective bath with no microscopic origin",
    "endpoint accounting can overestimate uninterrupted dwell; the continuum "
    "killed-diffusion oracle that measures that error lives outside this "
    "module and cannot be reached from it, and a record here carries no "
    "verdict from it",
)


def elementary_segments(mesh, start_index, count,
                        crossings) -> Tuple[Tuple[float, float], ...]:
    """The exact endpoints of every elementary interval of a step window.

    The elementary physical mesh is the finest uniform grid *union* the exact
    deterministic crossings, and this is that mesh written as ``(start, end)``
    pairs, in the same order ``PhaseNoiseStream.elementary_leaves`` emits its
    leaves.

    It exists because a :class:`~.stochastic.TreeLeaf` carries a ``duration``
    rather than an end time, and ``start_time + duration`` is a *different
    float* from the boundary the tree was split at.  Every eligibility test in
    the stepper is an exact equality against a crossing time, so a recomputed
    endpoint would put the last step of a window a rounding away from the exit
    and send it to ``stable_phase`` instead of ``boundary_phase``.  Deriving
    the endpoints here from the same ``mesh.step_bounds`` and the same crossing
    times the tree used, and cross-checking each against ``leaf.duration`` bit
    for bit in :meth:`~.dynamics.ClockPath.advance`, makes the two meshes one
    mesh or makes the run fail.

    Half-open routing, taken from ``FinestMesh.step_containing`` rather than
    re-derived: a crossing exactly on a step boundary belongs to the step that
    *starts* there, which is why an entry landing on a grid point produces a
    zero-duration handoff interval and not two overlapping ones.
    """
    require_exact_type(mesh, FinestMesh, "mesh")
    start_index = require_positive_int(start_index, "start_index", minimum=0)
    count = require_positive_int(count, "count", minimum=1)

    buckets: dict = {}
    for position, crossing in enumerate(crossings):
        moment = require_finite(getattr(crossing, "time", crossing),
                                f"crossings[{position}].time")
        buckets.setdefault(mesh.step_containing(moment), []).append(moment)

    segments = []
    for index in range(start_index, start_index + count):
        start, end = mesh.step_bounds(index)
        bounds = [start] + sorted(buckets.get(index, ())) + [end]
        for first, second in zip(bounds, bounds[1:]):
            segments.append((first, second))
    return tuple(segments)


@dataclass(frozen=True)
class OneClockRecord:
    """What one clock did, under one model, on one keyed noise realization.

    Everything needed to reproduce the run and everything measured by it, in
    one immutable object.  It is a *record*, not a result: it carries no
    prediction to be matched, no exponent, no probability, and no comparison.

    ``label``
        Free text naming the experiment cell — ``"S4/near-edge/strong"`` and
        the like.  Present so that a table of records reads without a legend.
    ``model``
        ``"full"`` or ``"width_only"``.  The central-clock control appears here
        as ``"full"`` with ``detuning`` exactly zero, because it *is* the full
        model; a separate label would imply separate dynamics.
    ``drive``
        ``"pulsed"`` or ``"stationary"`` — which half of experiment S4 this
        record belongs to.  A separate field rather than a convention inside
        ``label``, so that the two halves cannot be aggregated by accident: a
        stationary clock has no rising or falling edge and no finite
        eligibility window, so its band entries and resets are not comparable
        with a pulsed clock's.
    ``category``
        One of :data:`~.commitment.TRIAL_CATEGORIES`.
    ``committed_at``
        The endpoint at which the dwell completed, or ``None``.
        Interval-censored to one elementary interval by construction.
    ``winding``
        Net whole turns the authoritative unwrapped phase travelled.  Kept
        because it is the quantity a wrapped state variable would destroy, and
        because a control that gained a turn at a coordinate handoff would show
        up here and nowhere else.
    ``lift_residual``
        The largest ``(target_lift + error) - phase`` seen at any endpoint.
        Expected to be exactly zero, and reported anyway: every eligible step
        *assigns* the phase as that sum, so a nonzero value here would mean
        some other path had written a phase and left the pair describing a
        different one.  The one place the identity is not bit-exact is the
        entry handoff, where the phase is deliberately left untouched and the
        pair derived from it; that instant is a coordinate change rather than
        an endpoint and is checked directly rather than through this field.
    ``eligible_endpoints`` / ``boundary_endpoints``
        How many endpoints were strictly inside the tongue and how many were
        exactly on a crossing.  The second is normally small and is *not*
        normally zero: an exact entry or exit landing on a grid point is a
        real, tested case.
    """

    label: str
    model: str
    drive: str
    detuning: float
    peak_coupling: float
    phase_diffusion: float
    fixed_rate: float
    initial_phase: float
    trial: int
    clock: int
    timestep: float
    tolerance: float
    dwell: float
    category: str
    committed_at: Optional[float]
    band_entries: int
    dwell_resets: int
    endpoints: int
    eligible_endpoints: int
    boundary_endpoints: int
    final_phase: float
    winding: float
    lift_residual: float

    def __post_init__(self):
        # Every field is checked, including the ones this package fills in
        # itself.  A record is the thing a later ticket's ledger is built from
        # and the thing a reader quotes, so an unvalidated NaN sitting in one
        # would be a physical-looking number nobody had ever looked at.
        if not isinstance(self.label, str) or not self.label:
            raise TypeError(f"label must be a non-empty string, got "
                            f"{type(self.label).__name__} {self.label!r}")
        if self.model not in MODEL_LABELS:
            raise ValueError(f"model must be one of {list(MODEL_LABELS)}, got "
                             f"{self.model!r}")
        if self.drive not in DRIVE_LABELS:
            raise ValueError(f"drive must be one of {list(DRIVE_LABELS)}, got "
                             f"{self.drive!r}")
        if self.category not in TRIAL_CATEGORIES:
            raise ValueError(f"category must be one of "
                             f"{list(TRIAL_CATEGORIES)}, got "
                             f"{self.category!r}")
        for name in ("detuning", "initial_phase", "final_phase", "winding"):
            object.__setattr__(self, name,
                               require_finite(getattr(self, name), name))
        for name in ("peak_coupling", "phase_diffusion", "fixed_rate",
                     "lift_residual"):
            value = require_finite(getattr(self, name), name)
            if value < 0.0:
                raise ValueError(f"{name} must be non-negative, got {value!r}")
            object.__setattr__(self, name, value)
        for name in ("timestep", "tolerance", "dwell"):
            value = require_finite(getattr(self, name), name)
            if not value > 0.0:
                raise ValueError(f"{name} must be strictly positive, got "
                                 f"{value!r}")
            object.__setattr__(self, name, value)
        for name in ("trial", "clock", "band_entries", "dwell_resets",
                     "endpoints", "eligible_endpoints", "boundary_endpoints"):
            object.__setattr__(
                self, name,
                require_positive_int(getattr(self, name), name, minimum=0))
        if self.committed_at is not None:
            object.__setattr__(self, "committed_at",
                               require_finite(self.committed_at,
                                              "committed_at"))
        if (self.committed_at is None) != (self.category != "committed"):
            raise ValueError(
                f"category {self.category!r} and committed_at "
                f"{self.committed_at!r} disagree about whether this clock "
                f"committed; the category is derived from the dwell state and "
                f"cannot be set independently of the time it names")

        # Cross-field relations.  An independent review pointed out that this
        # constructor said "every field is checked" and checked every field
        # *individually*, so it accepted a record that was internally
        # impossible -- ``category="never_eligible"`` alongside three band
        # entries, nine resets, one endpoint, five eligible endpoints and seven
        # boundary endpoints.  Every generated record obeyed these relations;
        # nothing made them obey.  A record is the handoff to a later ticket's
        # ledger, so a contradiction should fail where it is created rather
        # than become a plausible raw row somebody later reconciles.
        if self.dwell_resets > self.band_entries:
            raise ValueError(
                f"dwell_resets={self.dwell_resets} exceeds band_entries="
                f"{self.band_entries}; a reset ends a dwell some entry had to "
                f"start, and a cleared clock needs a fresh entry before it can "
                f"be reset again")
        if self.band_entries > self.eligible_endpoints:
            raise ValueError(
                f"band_entries={self.band_entries} exceeds eligible_endpoints="
                f"{self.eligible_endpoints}; a dwell can only start at a "
                f"strictly eligible endpoint, so there cannot be more entries "
                f"than there were chances to enter")
        if self.eligible_endpoints + self.boundary_endpoints > self.endpoints:
            raise ValueError(
                f"eligible_endpoints={self.eligible_endpoints} and "
                f"boundary_endpoints={self.boundary_endpoints} sum to more "
                f"than the {self.endpoints} endpoints they are subtotals of; "
                f"interior and boundary are disjoint states of the same "
                f"endpoint list")
        if self.category == "never_eligible" and self.eligible_endpoints:
            raise ValueError(
                f"category 'never_eligible' with "
                f"{self.eligible_endpoints} strictly eligible endpoints; the "
                f"category means the clock never entered its tongue")
        if self.category in ("committed", "dwell_failed") \
                and self.band_entries < 1:
            raise ValueError(
                f"category {self.category!r} needs at least one band entry: "
                f"both categories mean the clock reached the band, and "
                f"'lock_failed' is the one reserved for a clock that never "
                f"did")
        if self.category == "lock_failed" and self.band_entries:
            raise ValueError(
                f"category 'lock_failed' with {self.band_entries} band "
                f"entries; a clock that reached the band and lost it is "
                f"'dwell_failed'")
        if self.category != "never_eligible" and self.eligible_endpoints < 1:
            raise ValueError(
                f"category {self.category!r} means the clock was eligible at "
                f"some endpoint, but eligible_endpoints is "
                f"{self.eligible_endpoints}")
        expected_rate = (FIXED_CONTRACTION_RATE if self.model == "width_only"
                         else 0.0)
        if self.fixed_rate != expected_rate:
            raise ValueError(
                f"a {self.model!r} record must carry a fixed contraction rate "
                f"of exactly {expected_rate!r}, got {self.fixed_rate!r}; an "
                f"ordinary record is a record of the one frozen control, and a "
                f"row carrying some other rate would enter a later ledger "
                f"indistinguishable from one that did not")

        # Conservation of band entries.  At every endpoint the state machine
        # can reach, each past entry is either the one currently live or has
        # already been reset, so ``entries - resets`` is 0 (nothing live) or 1
        # (one live).  Five independent live dwells are not a rare state in a
        # one-clock machine; they are an impossible one, and the constructor
        # used to accept them.
        live = self.band_entries - self.dwell_resets
        if live not in (0, 1):
            raise ValueError(
                f"band_entries={self.band_entries} minus dwell_resets="
                f"{self.dwell_resets} is {live}; one clock has at most one "
                f"live dwell, so every entry is either that one or has been "
                f"reset and the difference can only be 0 or 1")
        if self.category == "committed" and live != 1:
            raise ValueError(
                f"a committed record ends on the live dwell that completed, so "
                f"band_entries must be exactly one more than dwell_resets; got "
                f"{self.band_entries} and {self.dwell_resets}")
        if self.category in ("never_eligible", "lock_failed") \
                and (self.band_entries or self.dwell_resets):
            raise ValueError(
                f"category {self.category!r} means the clock never reached the "
                f"band, so it has no entries to count and none to reset; got "
                f"{self.band_entries} entries and {self.dwell_resets} resets")
        if self.drive == "stationary" and self.boundary_endpoints:
            raise ValueError(
                f"a stationary record reports {self.boundary_endpoints} "
                f"boundary endpoints, but a stationary drive has no tongue "
                f"crossing for an endpoint to sit on: its schedule is an "
                f"OpenSchedule with no boundary at all")

    @property
    def committed(self) -> bool:
        """Did this clock complete its dwell before the run ended?"""
        return self.committed_at is not None


def synthetic_dwell_history(criterion, times, eligible, lock_errors,
                            contractions) -> Tuple[DwellState, ...]:
    """Drive the dwell machine with declared sequences (experiment S2).

    Four equal-length sequences, one entry per endpoint, in time order: the
    endpoint times, whether the clock was strictly eligible there, its signed
    circular lock error, and the local contraction.  Returns the state after
    every endpoint, starting from the initial empty state, so a test can assert
    on the *history* rather than only on the last answer — which is the
    difference between "it committed" and "it committed at the right sample,
    having reset exactly once on the way".

    The lock error and contraction are ignored wherever ``eligible`` is false
    and may be NaN there, which is the same missing-target sentinel the rest of
    the package uses.  The walk stops at the first commitment and the returned
    history is truncated there, because the machine refuses to advance a
    committed state and a run that continued past its own stop would be
    inventing time.

    Nothing in here evolves a phase.  That is the point: every edge of the
    state machine is reachable from a sequence a test chose, including the
    edges a trajectory would visit only by luck.
    """
    # The same trusted snapshot ``run_one_clock`` takes, for the same reason
    # and before the same line.  This is the state machine's *other* public
    # execution boundary, and a contract enforced at one door and not the other
    # is not a contract.  ``LockCriterion`` is frozen against ordinary
    # assignment and not against ``object.__setattr__``, so a criterion edited
    # after construction arrived here carrying a band its own constructor would
    # have refused: at ``tolerance = pi`` a lock error of 3.0 rad — outside any
    # real band — produced a *committed* history.  The rebuilt copy is what
    # drives the sequence; the original is discarded.
    criterion = validate_lock_criterion(criterion)
    moments = [require_finite(t, f"times[{i}]")
               for i, t in enumerate(times)]
    flags = list(eligible)
    errors = list(lock_errors)
    slopes = list(contractions)
    lengths = {len(moments), len(flags), len(errors), len(slopes)}
    if len(lengths) != 1:
        raise ValueError(
            f"times, eligible, lock_errors and contractions must be the same "
            f"length, got {len(moments)}, {len(flags)}, {len(errors)} and "
            f"{len(slopes)}")
    if not moments:
        raise ValueError("a synthetic dwell history needs at least one "
                         "endpoint")
    for first, second in zip(moments, moments[1:]):
        if not second > first:
            raise ValueError(f"endpoint times must strictly increase, got "
                             f"{first!r} then {second!r}")

    # The criterion was rebuilt once, above; the loop uses the trusted kernel
    # so the rebuild does not repeat once per endpoint.  The states are this
    # function's own — the first is constructed here and every later one comes
    # out of the kernel — so they need no revalidation either.
    state = DwellState()
    history = [state]
    for moment, flag, error, slope in zip(moments, flags, errors, slopes):
        state = _advance_trusted(criterion, state, moment, flag, error, slope)
        history.append(state)
        if state.committed:
            break
    return tuple(history)


def run_one_clock(path, criterion, stream, trial, clock, initial_phase,
                  start_index, steps, label) -> Tuple[OneClockRecord,
                                                      Tuple[float, ...]]:
    """Walk one clock over one keyed noise realization (experiments S1 and S4).

    Returns the record and the authoritative unwrapped phase at every endpoint,
    the second so that a check can assert on the trajectory rather than only on
    the summary.

    The ordering is the plan's, once per elementary interval and in this order:

    1. evolve the phase to the interval's exact end;
    2. evaluate the coupling and the eligibility state *there*;
    3. inside the tongue, take the moving target and the local contraction;
    4. measure the circular error to that target; and
    5. update the dwell machine, which stops the run on commitment.

    The intervals are the elementary physical mesh, not the uniform grid: an
    exact tongue entry or exit inside a step subdivides it, the noise tree
    divides that step's already-drawn kick conditionally between the two
    halves, and the walk sees two intervals where the grid has one.  That is
    what "no dwell may begin before exact eligibility" means operationally —
    the eligibility change happens at an interval boundary, so no interval is
    ever half-eligible.

    ``stream`` supplies the leaves and knows nothing about tongues; the
    crossing schedule comes from ``path``.  The two are keyed to the same
    mesh and the same crossing identities, which is what makes the full model
    and the width-only control run on one realization.
    """
    # The trusted snapshot, first and before anything else.  Both objects are
    # rebuilt from their own fields rather than accepted on the strength of
    # their type, and the rebuilt copies are what this run uses; see
    # ``dynamics.validate_clock_path`` for the stale-detuning trajectory that
    # made this necessary.  Doing it here, at the top, is what makes every
    # refusal free: nothing below has run, so no key has been derived and no
    # Brownian leaf has been generated, and a rejected run cannot have changed
    # the realization it was rejected from.
    path = validate_clock_path(path)
    criterion = validate_lock_criterion(criterion)
    require_exact_type(stream, PhaseNoiseStream, "stream")
    trial = require_positive_int(trial, "trial", minimum=0)
    clock = require_positive_int(clock, "clock", minimum=0)
    phase = require_finite(initial_phase, "initial_phase")
    start_index = require_positive_int(start_index, "start_index", minimum=0)
    steps = require_positive_int(steps, "steps", minimum=1)
    if not isinstance(label, str) or not label:
        raise TypeError(f"label must be a non-empty string, got "
                        f"{type(label).__name__} {label!r}")

    # The frozen rate is already enforced, by the rebuild above: reconstructing
    # the path runs ``ClockPath.__post_init__``, which requires exactly 0.0 for
    # the full model and exactly FIXED_CONTRACTION_RATE for the control.  It is
    # asserted again here as a post-condition rather than a check, because the
    # rate is the one field of this record a later ledger cannot re-derive.
    expected_rate = (FIXED_CONTRACTION_RATE if path.model == "width_only"
                     else 0.0)
    if path.fixed_rate != expected_rate:
        raise ValueError(
            f"the rebuilt {path.model!r} path carries a fixed contraction rate "
            f"of {path.fixed_rate!r}, not the required {expected_rate!r}; this "
            f"is unreachable through ClockPath's own constructor and means the "
            f"snapshot above stopped rebuilding the field")

    # The starting state is built *before* any randomness is drawn.  That
    # ordering is the contract, not an accident of layout: an initial phase too
    # coarse to carry a faithful lifted target is refused by
    # ``ClockPath.initial_state``, and a run refused after consuming keys has
    # already changed the realization it was refused from.
    # The path was rebuilt once, above.  Everything below reaches for the
    # trusted kernels and the private cached schedule, so the pulse train and
    # its geometry are not re-derived once per elementary leaf.
    state = path._initial_state_trusted(stream.mesh.time_at(start_index),
                                        phase)

    crossings = path._cached_schedule.crossings(clock, path.train.digest)
    leaves = stream.elementary_leaves(trial, clock, start_index, steps,
                                      crossings)
    segments = elementary_segments(stream.mesh, start_index, steps, crossings)
    if len(leaves) != len(segments):
        raise ValueError(
            f"the noise tree produced {len(leaves)} elementary leaves and the "
            f"schedule {len(segments)} elementary intervals; the phase and the "
            f"noise are not being walked on one mesh")

    dwell = DwellState()
    phases = [state.phase]
    eligible_endpoints = 0
    boundary_endpoints = 0
    worst_residual = 0.0

    for leaf, (first, second) in zip(leaves, segments):
        state = path._advance_trusted(state, leaf, second)
        worst_residual = max(worst_residual, abs(state.lift_residual))
        if second == first:
            continue                      # zero-duration handoff: not an endpoint
        seen = path._observe_trusted(state)
        if seen.state == "interior":
            eligible_endpoints += 1
        elif seen.state == "boundary":
            boundary_endpoints += 1
        phases.append(state.phase)
        dwell = _advance_trusted(criterion, dwell, seen.time, seen.eligible,
                                 seen.lock_error, seen.contraction)
        if dwell.committed:
            break

    record = OneClockRecord(
        label=label,
        model=path.model,
        drive=path.train.drive_label,
        detuning=path.detuning,
        peak_coupling=path.train.peak_coupling,
        phase_diffusion=stream.phase_diffusion,
        fixed_rate=path.fixed_rate,
        initial_phase=phase,
        trial=trial,
        clock=clock,
        timestep=stream.mesh.step,
        tolerance=criterion.tolerance,
        dwell=criterion.dwell,
        category=dwell.category,
        committed_at=dwell.committed_at,
        band_entries=dwell.band_entries,
        dwell_resets=dwell.dwell_resets,
        endpoints=len(phases) - 1,
        eligible_endpoints=eligible_endpoints,
        boundary_endpoints=boundary_endpoints,
        final_phase=state.phase,
        winding=(state.phase - phase) / (2.0 * math.pi),
        lift_residual=worst_residual)
    return record, tuple(phases)
