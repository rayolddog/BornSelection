#!/usr/bin/env python3
"""
commitment.py — the lock predicate and the timestamp dwell clock, and nothing
else.

What this module is
-------------------
A state machine over *timestamps*.  It is handed, once per integration
endpoint, four things — the time, whether the clock is strictly eligible, its
signed circular error from the moving target, and the local contraction of the
deterministic flow at that phase — and it answers one question: has this clock
now stayed continuously locked for the fixed physical dwell?

There is no phase here, no coupling, no mismatch, no noise, no integrator and
no pulse.  This file cannot evolve anything.  That separation is the reason it
can be tested against declared synthetic sequences (experiment S2) rather than
only against trajectories that happen to arise: a state machine driven by
numbers a test chooses is a state machine whose edges can be visited on
purpose.

Why the dwell is a timestamp and not a counter
----------------------------------------------
The original plan added one timestep to a counter at every inside sample and
committed when the counter reached the dwell.  The pressure test rejected that,
and the counterexample is one line long: with ``dwell == dt`` and the very
first evaluated sample inside the band, the counter reaches the dwell on that
same sample and commits — from a single observation, with no elapsed time
inside the band at all.  Worse, ``ceil(dwell / dt) * dt`` silently becomes the
*actual* dwell whenever ``dwell / dt`` is not an integer, so every refinement
level would be testing a slightly different commitment rule and a timestep
convergence study would be comparing two different physical experiments.

So dwell is defined by two timestamps and a subtraction:

* the first sample that qualifies records ``inside_since = t``, and receives no
  credit for the interval before it — the clock may have entered the band at
  any point since the previous sample and this scheme deliberately assumes the
  latest such moment;
* commitment happens at the first later sample with
  ``t - inside_since >= dwell``; and
* anything that breaks the lock clears ``inside_since`` outright, so the next
  entry starts a fresh dwell from zero rather than resuming.

That rule is identical for ``dwell < dt``, ``dwell == dt`` and a non-integral
``dwell / dt``.  The consequence is that with a positive dwell the first inside
sample can never commit, which is what experiment S2 requires, and that the
event time is interval-censored to one timestep.

Why proximity is not enough: the contraction guard
--------------------------------------------------
Inside an Arnold tongue the attracting and repelling fixed phases approach one
another as the tongue edge is approached, and at the edge they merge.  Any
*fixed* angular tolerance therefore eventually contains both of them.  With
``K = 1`` and ``Delta = 0.999`` the stable phase is ``1.5261`` rad and the
repelling one is ``1.6155`` rad — ``0.0894`` rad apart — so a lock band of
``0.1`` rad declares the *repelling* fixed point locked.  No noise is needed
for that mistake, and the slow flow near the saddle-node can let the false lock
survive a short dwell.  This is precisely the region where the experiment
expects critical slowing to *suppress* commitment, so a proximity-only band
would manufacture events exactly where the hypothesis predicts their absence.

The guard is the smallest mathematically consistent repair: a phase earns dwell
only when it is close to the target **and** the local deterministic flow is
strictly contracting there, meaning a small displacement is pulled smaller
rather than pushed larger.  Both halves are strict.  Equality at either
boundary — a phase exactly at the tolerance, or a contraction of exactly zero —
earns nothing, because the merged fixed point at the tongue edge has exactly
zero contraction and calling it locked would be the original error rewritten.

This is not a claim about absorbers.  It is what "locally locking" has to mean
in a phase-only model, and it changes the accepted region near the tongue edge
in a way that any predictor compared against these events must also use.

What a commitment is not
------------------------
It is not a detector click, not an absorption, not a measurement outcome and
not a Born event.  It is the name this experiment gives to one clock having
stayed inside a declared band around its moving target for a declared time.
Whether that has anything to do with quantum measurement is the question the
whole programme is trying to ask, and it cannot be answered by naming.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional

from .validation import require_finite

__all__ = [
    "DWELL_ALLOWANCE_FRACTION",
    "DWELL_EQUALITY_ULPS",
    "TRIAL_CATEGORIES",
    "DwellState",
    "LockCriterion",
    "elapsed_allowance",
    "validate_lock_criterion",
]

#: The four mutually exclusive, exhaustive outcomes of one clock's exposure.
#: Named here rather than at the population layer because the vocabulary has to
#: be settled before anything uses it, and because a per-clock record carries
#: exactly the same four possibilities as the channel-level ledger will.
#:
#: ``lock_failed`` is reserved for a clock that was eligible and *never* reached
#: the band; a clock that reached it and then lost it is ``dwell_failed``.  The
#: distinction is load-bearing: it separates "the tongue never captured this
#: clock" from "noise broke the dwell", which are different physical statements
#: about the same unresolved trial.
TRIAL_CATEGORIES = ("committed", "never_eligible", "lock_failed",
                    "dwell_failed")

#: How many units in the last place of the *timeline* the elapsed-time
#: comparison may forgive.  See :func:`elapsed_allowance`.
DWELL_EQUALITY_ULPS = 4

#: The largest share of the declared dwell that roundoff allowance may occupy
#: before the timeline is declared unable to resolve its own dwell.  One part
#: in a million: far above any rounding at ordinary scales, far below anything
#: that could shorten a dwell in a way an experiment would notice.
DWELL_ALLOWANCE_FRACTION = 1e-6


def elapsed_allowance(moment, inside_since, dwell) -> float:
    """The roundoff the elapsed-time threshold may forgive, at this scale.

    The dwell rule is a subtraction of two timestamps compared against a
    duration, and binary64 subtraction can round an *exact physical equality*
    the wrong way.  The case that found this is ordinary rather than exotic —
    a decimal mesh from a nonzero origin:

    ::

        inside_since = 0.1,  moment = 0.3,  dwell = 0.2
        0.3 - 0.1 == 0.19999999999999998   ->  0.19999999999999998 >= 0.2 is False

    so a dwell that physically completes at ``0.3`` was reported at ``0.4``,
    one whole sample late.  The event time then depends on how the timeline
    happens to be spelled in binary rather than on elapsed time, which is
    exactly what a *fixed physical* dwell rule is supposed to exclude.  (The
    earlier verifier avoided it by choosing a binary-exact ``dt = 0.125``,
    which is why the suite did not see it.)

    The allowance is scale-aware because the error is: ``moment`` and
    ``inside_since`` are each already rounded representations of the physical
    instants they name, so the difference can be off by roughly an ulp *of the
    timeline*, not an ulp of the difference.  The bound is therefore
    :data:`DWELL_EQUALITY_ULPS` ulps at the largest magnitude involved.  At a
    timeline of order one that is an allowance of about ``2e-16``, against a
    dwell of order ``0.1``.

    **It cannot grant a materially early dwell**, and that is enforced rather
    than argued: if the allowance is not smaller than
    :data:`DWELL_ALLOWANCE_FRACTION` of the dwell, the timeline cannot resolve
    its own commitment rule and this raises.  That is the large-origin case —
    timestamps near ``1e12`` with a dwell of ``1e-3`` have an ulp of ``1e-4``,
    which is not a rounding to forgive but a mesh that has lost the ability to
    measure the interval it is being asked about.  A loud failure there is the
    only honest answer; silently committing on a coin-flip of the last bit is
    not.
    """
    later = require_finite(moment, "moment")
    earlier = require_finite(inside_since, "inside_since")
    span = _require_positive(dwell, "dwell")
    scale = max(abs(later), abs(earlier), span)
    allowance = DWELL_EQUALITY_ULPS * math.ulp(scale) if scale else 0.0
    if not allowance < DWELL_ALLOWANCE_FRACTION * span:
        raise ValueError(
            f"at timestamps of magnitude {scale!r} the binary64 rounding of an "
            f"elapsed time is {allowance!r}, which is not small compared with "
            f"the declared dwell {span!r}; this timeline cannot resolve its "
            f"own commitment rule, so the run is refused rather than committed "
            f"on the last bit of a subtraction")
    return allowance


def _require_positive(value, name: str) -> float:
    """Return ``value`` as a strictly positive float, or raise naming it."""
    out = require_finite(value, name)
    if not out > 0.0:
        raise ValueError(f"{name} must be strictly positive, got {out!r}")
    return out


def _require_bool(value, name: str) -> bool:
    """Return ``value`` if it is a genuine ``bool``, else raise naming it.

    Coercion is refused rather than performed, for the reason
    ``validation.require_boolean_array`` gives at the array boundary: ``0`` and
    ``"False"`` both become a plausible flag and neither is one.  Eligibility
    is a predicate and has to arrive as a predicate.
    """
    if not isinstance(value, bool):
        raise TypeError(f"{name} must be a bool, got "
                        f"{type(value).__name__} {value!r}")
    return value


def _require_count(value, name: str) -> int:
    """Return ``value`` as a non-negative plain int, or raise naming it."""
    if isinstance(value, bool):
        raise TypeError(f"{name} must be an integer, got bool {value!r}")
    if not isinstance(value, int):
        raise TypeError(f"{name} must be an integer, got "
                        f"{type(value).__name__} {value!r}")
    if value < 0:
        raise ValueError(f"{name} must be non-negative, got {value!r}")
    return value


def _require_time_or_none(value, name: str) -> Optional[float]:
    """Return ``None`` unchanged, or ``value`` as a finite float."""
    if value is None:
        return None
    return require_finite(value, name)


def _require_measurement(value, name: str, eligible: bool) -> float:
    """Return a lock measurement, allowing NaN only while ineligible.

    NaN is this package's documented sentinel for "no target existed at this
    sample" — ``model.stable_phase`` emits it wherever a clock is outside its
    tongue — so an ineligible endpoint legitimately has no error and no
    contraction to report.  While *eligible* the same NaN is a numerical
    failure wearing the sentinel's clothes, and is refused: a NaN error fails
    ``abs(error) < tolerance`` and so would silently read as "outside the
    band", which is a physical-looking conclusion manufactured from a broken
    number.
    """
    if eligible:
        return require_finite(value, name)
    if isinstance(value, (bool, str, bytes, bytearray, complex)):
        raise TypeError(f"{name} must be a real number, got "
                        f"{type(value).__name__} {value!r}")
    try:
        out = float(value)
    except (TypeError, ValueError):
        raise TypeError(f"{name} must be a real number, got "
                        f"{type(value).__name__} {value!r}") from None
    if math.isinf(out):
        raise ValueError(f"{name} may be NaN while the clock is ineligible — "
                         f"that is the missing-target sentinel — but never "
                         f"infinite, which is a numerical blow-up rather than "
                         f"a missing sample")
    return out


@dataclass(frozen=True)
class DwellState:
    """One clock's lock bookkeeping, as of the endpoint that produced it.

    Immutable, and rebuilt rather than mutated at every endpoint, so a caller
    holding an earlier state still holds what was true then.  That matters for
    the synthetic sequences of experiment S2, which assert on the whole history
    rather than on a final answer.

    ``inside_since``
        The timestamp of the first sample of the current uninterrupted inside
        run, or ``None`` when the clock is not currently inside.  ``None`` and
        not a sentinel number: a dwell that has not started is not a dwell that
        started at zero, and on a mesh whose times straddle zero the two would
        be indistinguishable.
    ``band_entries``
        How many times a fresh ``inside_since`` was started.  Distinguishes
        noise-assisted band entry from noise-broken dwell when read alongside
        ``dwell_resets``.
    ``dwell_resets``
        How many started dwells were cleared before completing.
    ``committed_at``
        The sampled time of commitment, or ``None``.  Interval-censored to one
        timestep by construction: it is an endpoint, and the crossing that
        earned it happened somewhere in the preceding step.
    ``ever_eligible`` / ``ever_inside``
        Latching flags, which are what make the unresolved categories
        distinguishable after the fact.
    ``last_time``
        The timestamp of the most recent endpoint this state has seen, or
        ``None`` before the first.  It exists so that :meth:`advance` can
        refuse a duplicate or a backwards observation.  Without it the public
        state transition accepted both — a dwell started at ``t = 1.0`` would
        take another call at ``1.0``, and a call at ``0.9``, and return a
        perfectly ordinary-looking state each time — while only the
        ``synthetic_dwell_history`` helper checked the ordering.  The helper is
        not the state machine; this is.
    """

    inside_since: Optional[float] = None
    band_entries: int = 0
    dwell_resets: int = 0
    committed_at: Optional[float] = None
    ever_eligible: bool = False
    ever_inside: bool = False
    last_time: Optional[float] = None

    def __post_init__(self):
        inside_since = _require_time_or_none(self.inside_since, "inside_since")
        band_entries = _require_count(self.band_entries, "band_entries")
        dwell_resets = _require_count(self.dwell_resets, "dwell_resets")
        committed_at = _require_time_or_none(self.committed_at, "committed_at")
        ever_eligible = _require_bool(self.ever_eligible, "ever_eligible")
        ever_inside = _require_bool(self.ever_inside, "ever_inside")
        last_time = _require_time_or_none(self.last_time, "last_time")

        # Cross-field rules.  Each rules out a state that is arithmetically
        # representable and physically impossible, and a state machine whose
        # impossible states are unrepresentable is one fewer thing to test for
        # downstream.  The time-ordering rules were added after an independent
        # review built a state that had committed *before* its own dwell began.
        # Conservation of band entries, and the latch that follows from it.
        #
        # At every endpoint this machine can reach, each past entry is either
        # the one currently live or has already been reset — there is one
        # clock, so there is at most one live dwell.  Therefore
        #
        #     band_entries == dwell_resets + int(inside_since is not None)
        #     ever_inside  == (band_entries > 0)
        #
        # exactly, and every proof of that is one line of `advance`: the
        # not-inside branch adds a reset exactly when it clears a live dwell,
        # the fresh-entry branch adds an entry exactly when it starts one, and
        # the continuing branch touches neither.
        #
        # An independent review built four states that satisfied every earlier
        # rule and none of these: a historical band entry with no eligibility
        # (`band_entries=1` alone), a `dwell_failed` category with zero entries
        # (`ever_inside=True` alone), an entry simultaneously live and reset,
        # and two ended entries with no reset.  All four are now refused, which
        # matters because these objects are the handoff to a later ledger and a
        # plausible-looking impossible history is worse than a missing one.
        live = 1 if inside_since is not None else 0
        if band_entries != dwell_resets + live:
            raise ValueError(
                f"band_entries={band_entries} does not equal dwell_resets="
                f"{dwell_resets} plus the {live} currently live dwell; one "
                f"clock has at most one live dwell, so every entry it ever made "
                f"is either that one or has been reset, and no other "
                f"combination is reachable")
        if ever_inside != (band_entries > 0):
            raise ValueError(
                f"ever_inside={ever_inside} disagrees with band_entries="
                f"{band_entries}; the latch is set exactly when the first entry "
                f"is counted and never cleared, so the two cannot differ")
        if committed_at is not None and inside_since is None:
            raise ValueError("committed_at is set with no running dwell; a "
                             "commitment is the end of an uninterrupted dwell, "
                             "not a state reachable without one")
        if ever_inside and not ever_eligible:
            raise ValueError("ever_inside is True but ever_eligible is False; "
                             "an ineligible clock has no target to be near")
        if last_time is None and (band_entries or dwell_resets or ever_eligible
                                  or ever_inside or inside_since is not None
                                  or committed_at is not None):
            raise ValueError(
                f"last_time is unset, so this state has never processed an "
                f"endpoint, yet it carries history: band_entries="
                f"{band_entries}, dwell_resets={dwell_resets}, ever_eligible="
                f"{ever_eligible}, ever_inside={ever_inside}, inside_since="
                f"{inside_since!r}, committed_at={committed_at!r}.  The only "
                f"state with no last observation is the initial one")
        if (inside_since is not None and last_time is not None
                and last_time < inside_since):
            raise ValueError(
                f"last_time={last_time!r} precedes inside_since="
                f"{inside_since!r}; a dwell cannot have started after the most "
                f"recent endpoint that was observed")
        if committed_at is not None:
            if committed_at < inside_since:
                raise ValueError(
                    f"committed_at={committed_at!r} precedes inside_since="
                    f"{inside_since!r}; a commitment is the *end* of an "
                    f"uninterrupted dwell and cannot happen before it began")
            if committed_at != last_time:
                raise ValueError(
                    f"committed_at={committed_at!r} is not the most recent "
                    f"endpoint {last_time!r}; a commitment is declared at the "
                    f"endpoint being processed, so a state that committed "
                    f"somewhere other than its own last observation was never "
                    f"produced by this machine")

        object.__setattr__(self, "inside_since", inside_since)
        object.__setattr__(self, "band_entries", band_entries)
        object.__setattr__(self, "dwell_resets", dwell_resets)
        object.__setattr__(self, "committed_at", committed_at)
        object.__setattr__(self, "ever_eligible", ever_eligible)
        object.__setattr__(self, "ever_inside", ever_inside)
        object.__setattr__(self, "last_time", last_time)

    @property
    def committed(self) -> bool:
        """Has this clock completed its dwell?"""
        return self.committed_at is not None

    @property
    def category(self) -> str:
        """Which of :data:`TRIAL_CATEGORIES` this clock's exposure ended in.

        Read it only at the end of the exposure: mid-run it reports the
        category the clock *would* land in if the run stopped now, which is a
        useful thing and not the same thing.

        The three unresolved categories are mutually exclusive and exhaustive
        by construction here, rather than by a downstream convention: they are
        decided by two latching flags whose only three combinations are the
        three categories.
        """
        if self.committed_at is not None:
            return "committed"
        if not self.ever_eligible:
            return "never_eligible"
        if not self.ever_inside:
            return "lock_failed"
        return "dwell_failed"


@dataclass(frozen=True)
class LockCriterion:
    """The fixed angular band and the fixed physical dwell.

    ``tolerance``
        The half-width of the lock band, in radians, measured as a *circular*
        distance from the moving target.  Strictly less than half a turn: a
        tolerance of ``pi`` or more contains every phase on the circle, so it
        is not a band and the machine would degenerate into "eligible and
        contracting" without saying so.
    ``dwell``
        The uninterrupted physical time inside the band that earns a
        commitment.  Strictly positive: a zero dwell would commit on the first
        qualifying sample, which is the sample-counting rule the pressure test
        rejected, reached by a different route.

    **Neither number may ever depend on coupling.**  They are fixed once, in
    one reference time and angle system, for an entire coupling sweep.  A
    tolerance or a dwell that shrank as the coupling grew would manufacture the
    coupling dependence the experiment exists to measure, which is why the
    verifier's mutation battery contains exactly that mutation and requires it
    to fail.
    """

    tolerance: float
    dwell: float

    def __post_init__(self):
        tolerance = _require_positive(self.tolerance, "lock tolerance")
        dwell = _require_positive(self.dwell, "dwell time")
        if not tolerance < math.pi:
            raise ValueError(
                f"lock tolerance must be less than half a turn, got "
                f"{tolerance!r}; a circular distance never exceeds pi, so a "
                f"tolerance of pi or more admits every phase and is not a band")
        object.__setattr__(self, "tolerance", tolerance)
        object.__setattr__(self, "dwell", dwell)

    def _inside_now_trusted(self, lock_error, contraction) -> bool:
        """Unvalidated kernel behind :meth:`inside_now`.

        **Accepts only a rebuilt criterion.**  ``self`` must be an object that
        :func:`validate_lock_criterion` has just returned, or one this module
        constructed itself; nothing here re-checks the band or the dwell, and a
        criterion edited past its own constructor would be used exactly as
        found.  The scalars are still checked, because they come from a caller
        rather than from a validated object.

        The same private-kernel convention ``model.py`` uses for
        ``_wrap_phase`` and ``_adler_drift``: the public boundary validates,
        the inner loop does not, and the rule is that a kernel may only be
        called with values some boundary has already checked.
        """
        error = require_finite(lock_error, "lock_error")
        slope = require_finite(contraction, "contraction")
        return abs(error) < self.tolerance and slope < 0.0

    def inside_now(self, lock_error, contraction) -> bool:
        """Is this endpoint inside the band *and* locally contracting?

        ``lock_error`` is the signed circular distance from the phase to its
        moving target, in ``[-pi, pi)``; ``contraction`` is the derivative of
        the deterministic drift with respect to phase at that point.  Both
        comparisons are strict:

        * ``abs(lock_error) < tolerance`` — a phase exactly on the boundary is
          outside, so the band is an open interval and no dwell can be earned
          by touching its edge; and
        * ``contraction < 0`` — a flow that neither pulls nor pushes is not
          locking.  At the tongue edge the attracting and repelling phases have
          merged and the contraction is exactly zero, which is the case this
          strictness exists for.

        Nothing here knows what produced either number.  The full Adler model
        supplies ``-K cos(phase)`` and the fixed-contraction width-only control
        supplies ``-rate cos(error)``; the guard is the same logical rule for
        both, which is what stops the control from quietly degenerating into a
        proximity-only band.

        **The criterion is rebuilt before the comparison.**  This is a public
        predicate on an exported frozen dataclass, and a frozen dataclass is
        not frozen against ``object.__setattr__``: a criterion edited to
        ``tolerance = pi`` answered ``True`` for a lock error of ``3.0`` rad,
        which is outside the largest band this class can express.  The two
        helper boundaries refused that object and this door did not, so the
        rebuild happens here too.  Callers in a loop should validate once and
        use the private kernel instead.
        """
        return validate_lock_criterion(self)._inside_now_trusted(lock_error,
                                                                 contraction)

    def advance(self, state, time, eligible, lock_error,
                contraction) -> DwellState:
        """Return the state after one integration endpoint.

        The ordering is the plan's, literally and in this order: the phase has
        already been evolved to ``time`` by the caller, eligibility and the
        target were evaluated *there*, and this is the state update.  Nothing
        in this method looks backwards at the interval just traversed, which is
        exactly why the first inside sample gets no retroactive credit.

        Five branches, and they are the whole machine:

        1. an ineligible endpoint is never inside, whatever its error says —
           there is no target to be near, so the falling tongue edge clears a
           running dwell on the same step it removes eligibility;
        2. an eligible endpoint outside the band, or on the wrong side of the
           contraction boundary, clears a running dwell and counts a reset;
        3. an inside endpoint with no running dwell starts one at ``time`` and
           counts a band entry;
        4. an inside endpoint with a running dwell keeps the original
           ``inside_since`` — the dwell is measured from where it began, not
           from where it was last observed; and
        5. an inside endpoint whose elapsed time has reached the dwell commits.

        Branches 3 and 5 are evaluated in that order on the same endpoint, so a
        freshly started dwell is checked for completion immediately and, with a
        strictly positive dwell, cannot succeed: ``time - time`` is zero.  That
        is the S2 requirement that an initial inside sample starts the clock
        and never commits instantly, and it falls out of the arithmetic rather
        than being special-cased.

        Advancing a state that has already committed raises: the caller's run
        loop stops at commitment, and a machine that quietly kept going would
        turn a missing stop into a plausible later event time.

        **Time must strictly increase.**  A duplicate endpoint and a backwards
        endpoint are both refused here, in the public transition, rather than
        only in the ``synthetic_dwell_history`` helper that happened to check
        it.  Re-observing an instant would let a dwell be measured against a
        timestamp the clock had already passed, and going backwards would let
        it be measured against one it had not reached.

        The elapsed-time comparison forgives :func:`elapsed_allowance` of
        binary64 rounding at the timeline's own scale, so an *exact physical
        equality* on an ordinary decimal mesh commits when it should rather
        than one sample later.  The allowance is bounded to one part in a
        million of the dwell and the bound is enforced, so this cannot become
        a materially early commitment; see that function for the failure it
        repairs and the failure it refuses to hide.

        **Both frozen objects are rebuilt before the transition.**  This is the
        state machine — earlier reviews treated *this* method rather than any
        helper as the thing that has to hold the line — and neither the
        criterion nor the incoming state was being revalidated here.  A
        criterion edited to ``tolerance = pi`` committed from a ``3.0`` rad
        error; a state edited past ``DwellState``'s own reachability rules was
        refused only incidentally, when the *next* state happened to be
        constructed, which is after the transition rather than before it.  Both
        are rebuilt first now, so a refusal costs no state change at all.

        A caller stepping a loop should validate once and use the private
        kernel; :func:`synthetic_dwell_history` and
        ``raw_experiments.run_one_clock`` both do.
        """
        return _advance_trusted(validate_lock_criterion(self),
                                _validated_dwell_state(state), time, eligible,
                                lock_error, contraction)


def _validated_dwell_state(state) -> DwellState:
    """Return a freshly rebuilt, fully revalidated copy of ``state``.

    Named fields, one at a time, through the real constructor, so every
    reachability, conservation, latch and timing rule in ``DwellState`` runs
    again on an object that may have been edited since it last ran.  Private
    because it is the input half of a transition rather than something a caller
    needs; the public door is :meth:`LockCriterion.advance`.
    """
    if type(state) is not DwellState:
        raise TypeError(f"state must be exactly DwellState, got "
                        f"{type(state).__name__}")
    return DwellState(inside_since=state.inside_since,
                      band_entries=state.band_entries,
                      dwell_resets=state.dwell_resets,
                      committed_at=state.committed_at,
                      ever_eligible=state.ever_eligible,
                      ever_inside=state.ever_inside,
                      last_time=state.last_time)


def _advance_trusted(criterion, state, time, eligible, lock_error,
                     contraction) -> DwellState:
    """Unvalidated kernel behind :meth:`LockCriterion.advance`.

    **Accepts only rebuilt objects.**  ``criterion`` must be what
    :func:`validate_lock_criterion` returned and ``state`` what
    :func:`_validated_dwell_state` returned; neither is re-checked here.  The
    scalars are checked, because they arrive from a caller rather than from a
    validated object.

    This is the whole state machine, and the branches it documents are
    described on the public wrapper.
    """
    if state.committed_at is not None:
        raise ValueError(
            f"this clock committed at t={state.committed_at!r}; advancing "
            f"a committed state would continue a run that has already "
            f"ended, so the stop is enforced here rather than assumed")
    moment = require_finite(time, "time")
    if state.last_time is not None and not moment > state.last_time:
        raise ValueError(
            f"this state was last observed at t={state.last_time!r} and is "
            f"being advanced to t={moment!r}; endpoints must strictly "
            f"increase, so a repeated instant and a step backwards are "
            f"both refused rather than folded into the dwell")
    is_eligible = _require_bool(eligible, "eligible")
    error = _require_measurement(lock_error, "lock_error", is_eligible)
    slope = _require_measurement(contraction, "contraction", is_eligible)

    ever_eligible = state.ever_eligible or is_eligible
    inside = is_eligible and criterion._inside_now_trusted(error, slope)

    if not inside:
        return DwellState(
            inside_since=None,
            band_entries=state.band_entries,
            dwell_resets=state.dwell_resets
            + (1 if state.inside_since is not None else 0),
            committed_at=None,
            ever_eligible=ever_eligible,
            ever_inside=state.ever_inside,
            last_time=moment)

    if state.inside_since is None:
        inside_since = moment
        band_entries = state.band_entries + 1
    else:
        inside_since = state.inside_since
        band_entries = state.band_entries

    elapsed = moment - inside_since
    allowance = elapsed_allowance(moment, inside_since, criterion.dwell)
    reached = elapsed >= criterion.dwell or criterion.dwell - elapsed <= allowance
    return DwellState(
        inside_since=inside_since,
        band_entries=band_entries,
        dwell_resets=state.dwell_resets,
        committed_at=moment if reached else None,
        ever_eligible=ever_eligible,
        ever_inside=True,
        last_time=moment)


def validate_lock_criterion(criterion) -> LockCriterion:
    """Return a freshly rebuilt, fully revalidated copy of ``criterion``.

    The state machine's half of the same contract
    ``dynamics.validate_clock_path`` states for the path, applied at the same
    public execution boundary and for the same reason: ``LockCriterion`` is
    frozen against ordinary assignment and not against ``object.__setattr__``,
    so a criterion validated at construction and edited afterwards would arrive
    at a run carrying a band or a dwell its own constructor would have
    refused — a tolerance of ``pi`` or more, which admits every phase and is
    not a band, or a non-positive dwell, which commits on sight.

    The rebuild names the two fields one at a time and re-runs the constructor
    over them, so a criterion that no longer satisfies its own domain fails
    here rather than quietly setting the commitment rule for a whole run.

    **Semantics are unchanged.**  For a criterion that was never edited this
    returns the same two floats and therefore the same rule; the rebuild is
    exact, not a normalization, and no comparison, allowance or boundary moves.
    """
    if type(criterion) is not LockCriterion:
        actual = type(criterion).__name__
        raise TypeError(f"criterion must be exactly LockCriterion, got "
                        f"{actual}")
    return LockCriterion(criterion.tolerance, criterion.dwell)
