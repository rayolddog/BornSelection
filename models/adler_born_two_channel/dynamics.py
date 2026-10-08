#!/usr/bin/env python3
"""
dynamics.py — Euler-Maruyama phase evolution over the exact elementary mesh,
and the fixed-contraction width-only control that shares its noise.

What this module is
-------------------
The stepper.  ``stochastic.py`` supplies numbers and knows nothing about
phases; ``commitment.py`` runs a state machine and knows nothing about
dynamics; this is the file where a phase actually moves.  It takes one clock's
frozen noise leaves and walks them, applying

    dphi = (Delta - K(t) sin(phi)) dt + noise increment

with the drift frozen at the start of each elementary interval, which is
Euler-Maruyama.  At zero diffusion every leaf kick is exactly zero and the walk
reduces, term for term, to ``model._step_euler`` — the deterministic limit is
recovered exactly rather than approximately, and the verifier asserts bit
equality rather than a tolerance.

One authoritative phase
-----------------------
Each clock has exactly one authoritative state variable: an **unwrapped
absolute phase on the real line**.  It is never wrapped, never has its winding
discarded, and never replaced by a principal value.  Periodicity enters in
exactly three places and nowhere else — ``sin`` in the drift, ``cos`` in the
contraction, and the circular distance used to measure lock error — and each of
those is a periodic *function of* the authoritative phase rather than a
substitute for it.

The reason is not tidiness.  The width-only control below writes the phase as
``target_lift + error``, and the target moves; if the phase were wrapped, or if
the target were re-derived from a principal value at each entry, then a clock
sitting just across the ``+-pi`` cut would gain or lose a whole turn at a tongue
entry or exit even with the noise switched off.  A run that silently teleports
by ``2 pi`` at an eligibility crossing is not the run anybody meant to
integrate, and the accumulated winding is data.

The three eligibility states
----------------------------
The calibrated package makes eligibility **strict**: ``K > |Delta|``.  At
``K == |Delta|`` the attracting and repelling fixed points have merged,
``model.stable_phase`` deliberately returns NaN there, and ``boundary_phase``
supplies the finite limit ``sign(Delta) * pi/2``.  So a time is in one of three
states, and this module names all three rather than folding equality into one
of its neighbours:

``outside``
    Before entry or after exit.  No target exists, no dwell can accrue, and the
    full calibrated Adler drift still acts — eligibility controls whether there
    is something to lock onto, not whether the coupling term is switched on.
``boundary``
    Exactly at an entry or an exit.  The target is the finite
    ``boundary_phase(Delta)``; ``stable_phase`` is never called here.  This is a
    zero-duration coordinate handoff: no time passes, no noise is consumed, no
    phase changes, and **no dwell is earned**, because the merged fixed point
    has exactly zero contraction and the guard in ``commitment.py`` is strict.
``interior``
    Strictly inside.  ``stable_phase(K(t), Delta)`` is finite and is the only
    place it is evaluated.

Crossings are taken from the calibrated pulse rather than detected by watching
``K`` cross ``|Delta|`` numerically, and eligibility is decided by comparing
*times* against those exact crossings rather than by re-evaluating the
inequality at a sample.  Those two spellings disagree within a rounding of the
crossing, and the one that can disagree with the noise tree — which was split
at the exact crossing time — is the one that would put a step's target on the
wrong side of its own subdivision.

The width-only control, exactly
-------------------------------
The control exists to answer a specific causal question: if the *only* thing
removed is Adler's coupling-dependent interior contraction, does the coupling
dependence of commitment survive?  For that contrast to isolate anything, every
other thing must be identical — and an earlier draft of the plan failed this,
because it replaced the ineligible drift with free mismatch drift.  A marginal
clock on the rising edge feels the subthreshold Adler pull, so removing it
changes the phase *delivered to* tongue entry, and a paired difference would
then confound "we removed interior relaxation" with "we changed the entry
distribution".

So the contract implemented here is:

* while ineligible, the control is the full calibrated model, step for step,
  bit for bit, on the same noise leaf;
* at exact entry it stores a lifted copy of ``boundary_phase(Delta)`` and the
  unwrapped error from it, changing the authoritative phase by nothing;
* while strictly eligible it evolves that error with ``-rate sin(error) dt``
  plus the same leaf kick, continues the target onto its own lifted branch, and
  reconstructs ``phase = lift + error`` so that the target's motion is retained
  in absolute coordinates;
* at exact exit it keeps the reconstructed phase, discards only the auxiliary
  pair, and resumes the full Adler drift; and
* at re-entry it builds a new lift from the unchanged authoritative phase.

``rate`` is one positive number, common to every clock, coupling, pulse time
and trial.  It is frozen by a *rule*, not chosen: it is the central-clock Adler
relaxation rate at the geometric midpoint of a declared reference peak-coupling
range.  See :func:`fixed_contraction_rate`.  A rate that scaled with coupling
would reintroduce the very factor the control removes, which is why the
verifier's mutation battery contains that mutation and requires it to fail.

What is not here
----------------
No population, no race, no first-winner stop, no ledger written to disk, no
survival curve, no hazard, no fitted exponent, and no analytic prediction of
any kind.  This module cannot import the prediction and the verifier walks its
imports to say so.  It also contains no commitment rule: it reports an
observation at each endpoint and ``commitment.py`` decides what that means.
"""
from __future__ import annotations

import hashlib
import math
from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np

from .model import (RaisedCosinePulse, _circular_distance, _step_euler,
                    _wrap_phase, boundary_phase, stable_phase)
from .stochastic import Crossing, TreeLeaf, conservation_tolerance
from .stochastic import _require_digest as _require_digest_token
from .validation import (require_exact_type, require_finite,
                         require_finite_array, require_positive_int)

__all__ = [
    "DRIVE_LABELS",
    "DRIVE_TYPES",
    "ELIGIBILITY_STATES",
    "FIXED_CONTRACTION_RATE",
    "LIFT_RESOLUTION",
    "LIFT_ULPS",
    "MAX_LIFT_MAGNITUDE",
    "MODEL_LABELS",
    "PULSE_DIGEST_CHARS",
    "REFERENCE_PEAK_COUPLING_RANGE",
    "ClockPath",
    "Observation",
    "OpenSchedule",
    "PathState",
    "PulseTrain",
    "StationaryDrive",
    "TongueSchedule",
    "continuous_lift",
    "fixed_contraction_rate",
    "nearest_lift",
    "validate_clock_path",
]

#: The two data-generating processes.  ``full`` is the calibrated Adler model.
#: ``width_only`` is the mechanism-negative control described above.  The
#: *central-clock* control is not a third entry: it is the full model at zero
#: mismatch, and giving it its own label here would suggest it had its own
#: dynamics, which is exactly the thing that would make it a bad control.
MODEL_LABELS = ("full", "width_only")

#: The three states a time can be in relative to one clock's tongue.  Equality
#: is a state of its own and not a rounding of either neighbour; see the module
#: docstring.
ELIGIBILITY_STATES = ("outside", "boundary", "interior")

#: The two coupling histories experiment S4 runs a clock under.  A stationary
#: drive is not a very long pulse: it has no edges, no finite eligibility
#: window and no tongue crossings, so it is a different object rather than a
#: parameter choice.  Records carry the label so that the two halves of S4 can
#: never be read as one matrix.
DRIVE_LABELS = ("pulsed", "stationary")

#: Characters and length of the pulse-manifest digest carried in a crossing
#: identity.  Thirty-two hexadecimal characters of SHA-256: long enough that
#: two distinguishable pulse manifests will not collide, short enough to read
#: in a key.
PULSE_DIGEST_CHARS = 32

#: The declared reference peak-coupling range whose geometric midpoint fixes
#: the width-only control's contraction rate.
#:
#: **This is not the production range.**  The production range is chosen by the
#: labelled range-only pilot of a later ticket, under its frozen count rules,
#: and the number written into the production manifest is
#: ``fixed_contraction_rate(*that_range)``.  What is frozen *now*, before any
#: stochastic behaviour has been looked at, is the rule: the pilot may choose
#: the range and may not choose the rate.  This provisional range exists so
#: that the one-clock control of experiment S4 has a definite rate to run at,
#: and every record it produces carries that rate explicitly.
REFERENCE_PEAK_COUPLING_RANGE = (0.5, 2.0)

#: The angular resolution the lifted representation guarantees, in radians.
#:
#: A lift is only a lift if it is *congruent to its principal target modulo a
#: turn*, and on binary64 that stops being true long before the phase stops
#: being finite.  At ``phase = 1e16`` and principal target ``pi/2``:
#:
#: ::
#:
#:     nearest_lift(pi/2, 1e16) = 9999999999999998.0
#:     phase - lift             = 2.0
#:     wrap_phase(lift - pi/2)  = -0.5043502211754891
#:
#: — half a radian from the target it claims to lift.  At ``1e17`` the stored
#: error collapses to zero while the congruence error is order one.  The
#: reconstruction ``lift + error == phase`` still holds exactly, so the
#: identity check passes while the auxiliary target has quietly become the
#: wrong target, and the next contraction and lock error are computed against
#: it.
#:
#: There is no representation trick that fixes this inside one float: the phase
#: and its principal residual cannot both be resolved once the phase's ulp is
#: comparable to an angle.  So the domain is *declared and enforced* instead.
#: A nanoradian is far finer than any tolerance in this package — the
#: deterministic suite's tightest angular tolerance is ``1e-12``, and the lock
#: bands are ``1e-1`` — and far coarser than the ``1e-16``-scale roundoff of an
#: ordinary run.
LIFT_RESOLUTION = 1e-9

#: How many ulps of the phase the congruence error can reach.  ``phase -
#: principal``, the modulo, and ``phase - error`` each contribute up to half an
#: ulp of the phase; four is the margin.
LIFT_ULPS = 4

#: The strict bound on ``|phase|`` and ``|target_lift|``.
#:
#: ``ulp(x)`` is ``2**(k-52)`` throughout the binade ``[2**k, 2**(k+1))``, so
#: ``LIFT_ULPS * ulp(x) <= LIFT_RESOLUTION`` requires ``2**k <= 1e-9 * 2**50 =
#: 1125899.9...``.  The largest such binade starts at ``2**20``, so every
#: magnitude strictly below ``2**21`` is resolvable and ``2**21`` itself is
#: not: ``4 * ulp(2**21) = 1.86e-9`` exceeds the budget while ``4 *
#: ulp(2**21 - 1) = 9.31e-10`` does not.  ``verify.py`` re-derives both sides
#: of that boundary rather than trusting this comment.
#:
#: In physical terms this is about 333,772 turns of accumulated winding.  The
#: runs in this package wind single-digit turns, so the domain is a documented
#: guarantee rather than a restriction anything here comes near.
MAX_LIFT_MAGNITUDE = 2097152.0

#: Text prefix of the pulse manifest a digest is taken over.  Versioned for the
#: same reason ``stochastic.KEY_SCHEMA`` is: the digest is part of a crossing
#: identity and therefore part of a noise key, so a change to what is hashed is
#: a change of dataset and has to be visible as one.
_PULSE_SCHEMA = "dk-pulse-train/v1"

#: The same, for a stationary drive.  A distinct prefix, so that a stationary
#: run and a pulsed run can never hash to one manifest digest.
_STATIONARY_SCHEMA = "dk-stationary-drive/v1"


# ---------------------------------------------------------------------------
# 1. Lifts: continuing a principal angle onto the real line
# ---------------------------------------------------------------------------
def _require_resolvable(value, name: str) -> float:
    """Return ``value`` as a finite angle inside the resolvable domain.

    The rejection is deliberately at the *boundary* rather than deep inside
    the evolution: a phase this coarse is refused by the public lift
    functions, by :class:`PathState`, and by
    :meth:`ClockPath.initial_state` — all of which run before any noise is
    drawn — so a run that could not carry a faithful lift never consumes a
    key or takes a step.  See :data:`MAX_LIFT_MAGNITUDE`.
    """
    angle = require_finite(value, name)
    if not abs(angle) < MAX_LIFT_MAGNITUDE:
        raise ValueError(
            f"{name}={angle!r} has magnitude at or beyond "
            f"{MAX_LIFT_MAGNITUDE!r}, where one unit in the last place is "
            f"{math.ulp(abs(angle))!r} radians and a lifted target can no "
            f"longer be held congruent to its principal value to the declared "
            f"{LIFT_RESOLUTION!r} rad resolution; this package refuses such a "
            f"phase rather than returning a lift of the wrong target")
    return angle


def nearest_lift(principal, phase) -> float:
    """Return the lift of ``principal`` nearest ``phase``.

    The returned value is ``principal + 2 pi m`` for some integer ``m``, chosen
    so that ``phase - lift`` lies in the half-open interval ``[-pi, pi)``.

    Half-open, and that is the whole tie rule.  An exact half-turn separation
    is a real case — it happens whenever a clock sits exactly opposite its
    target — and "nearest" does not resolve it, so the convention has to.  It
    is the same ``[-pi, pi)`` convention ``model.wrap_phase`` uses, taken from
    the same function rather than re-improvised here, so the package makes this
    choice exactly once.

    Used at an exact tongue entry, where ``principal`` is the finite
    ``boundary_phase(Delta)`` and ``phase`` is the authoritative unwrapped
    phase.  The point of the operation is what it does *not* do: it does not
    move the phase.  It selects which copy of the target to measure against.

    ``phase`` is required to lie inside :data:`MAX_LIFT_MAGNITUDE`, because
    beyond it the returned value is no longer congruent to ``principal``
    modulo a turn and so is not a lift of anything.
    """
    target = require_finite(principal, "principal")
    absolute = _require_resolvable(phase, "phase")
    return absolute - float(_wrap_phase(absolute - target))


def _lift_and_error(principal: float, phase: float) -> Tuple[float, float]:
    """The lifted target and the unwrapped error, in that order.

    Private, and the *only* way either of them is constructed, because the two
    have to be built in the right order to keep the half-open convention exact.
    Taking the lift first and then subtracting — ``lift = nearest_lift(...)``,
    ``error = phase - lift`` — recovers the error only to within a rounding,
    and a rounding is enough to put it a fraction of an ulp outside
    ``[-pi, pi)``: an entry at a phase exactly a half turn from its target
    produced ``-3.1415926535897936``, which is below ``-pi`` and so not in the
    interval the convention names.

    So the error is taken from ``model._wrap_phase`` directly, where it is in
    ``[-pi, pi)`` by construction, and the lift is defined as ``phase - error``
    — the same expression :func:`nearest_lift` computes, so the two agree bit
    for bit while the error is now exact rather than reconstructed.
    """
    return nearest_lift(principal, phase), float(_wrap_phase(phase - principal))


def continuous_lift(principal, previous_lift) -> float:
    """Return the lift of ``principal`` nearest ``previous_lift``.

    The same construction as :func:`nearest_lift` with a different reference:
    here the branch is chosen to continue an existing lifted target rather than
    to sit near the phase.  That is what makes the target's motion *continuous*
    across a step, including a step during which the principal value crosses
    the ``+-pi`` cut.

    Without it, a target that wanders past the cut jumps by a full turn, the
    reconstructed phase jumps with it, and a clock is reported to have moved a
    whole turn in one timestep with no drift and no noise responsible for it.
    """
    target = require_finite(principal, "principal")
    previous = _require_resolvable(previous_lift, "previous_lift")
    return previous + float(_wrap_phase(target - previous))


def fixed_contraction_rate(lower_peak, upper_peak) -> float:
    """The width-only control's one contraction rate, from the frozen rule.

    The rule, in words: take the geometric midpoint of the declared reference
    peak-coupling range, and use the *central clock's* Adler relaxation rate at
    that coupling.

    In symbols, with ``K_ref = sqrt(lower * upper)``, the local relaxation rate
    of a clock at its stable phase is ``K cos(phi_stable)``, and for a
    zero-mismatch clock ``phi_stable = arcsin(0) = 0``, so the rate is
    ``K_ref``.  It is written out that way below rather than collapsed to the
    square root, because the collapsed form is a coincidence of the central
    clock and reading it as the definition is how a later edit would end up
    with a rate that means something else.

    Three properties matter more than the value:

    * **It is a rule, not a fit.**  The pilot of a later ticket selects the
      coupling *range* under count-only criteria and never inspects stochastic
      behaviour; this function then supplies the rate deterministically.  There
      is no path by which an event result can influence it.
    * **It is one number.**  Common to every clock, every coupling, every pulse
      time and every trial.  The control's whole purpose is to have removed the
      coupling-dependent interior contraction, so a rate that varied with
      coupling would be the full model wearing a different name.
    * **Lower and upper choices are sensitivities.**  They may be run later and
      labelled as such; they may never be substituted for the primary control.

    Geometric rather than arithmetic midpoint because the sweep is logarithmic
    in coupling and the estimand is a log-log exponent: the arithmetic midpoint
    of a decade sits four fifths of the way up it.

    **The midpoint is never formed as ``sqrt(lower * upper)``.**  That product
    is the intermediate this expression does not need and cannot afford:
    ``fixed_contraction_rate(1e-300, 1e-300)`` underflowed the product to
    exactly zero, ``stable_phase(0, 0)`` is NaN by design, and the function
    returned a NaN *rate* — a non-finite number arriving at the one parameter
    the control is supposed to hold fixed.  Large finite pairs overflowed the
    other way and raised on an infinity.  ``sqrt(lower) * sqrt(upper)`` is the
    same real number and cannot do either: each square root is at most about
    ``1.3e154`` and at least about ``2.2e-162``, so the product stays in range
    across the whole binary64 domain.

    The result is then required to be finite and strictly positive before it
    is returned, so no path here can hand a non-finite or non-positive rate to
    a control that is defined by having exactly one positive one.
    """
    lower = require_finite(lower_peak, "lower_peak")
    upper = require_finite(upper_peak, "upper_peak")
    if not lower > 0.0:
        raise ValueError(f"lower_peak must be strictly positive, got {lower!r}")
    if not upper >= lower:
        raise ValueError(f"upper_peak must be at least lower_peak, got "
                         f"lower_peak={lower!r}, upper_peak={upper!r}")
    reference = math.sqrt(lower) * math.sqrt(upper)
    if not (math.isfinite(reference) and reference > 0.0):
        raise ValueError(
            f"the geometric midpoint of [{lower!r}, {upper!r}] is "
            f"{reference!r}, which is not a usable reference coupling")
    central = float(stable_phase(reference, 0.0))
    if not math.isfinite(central):
        raise ValueError(
            f"the central clock has no stable phase at the reference coupling "
            f"{reference!r}, so no relaxation rate can be derived from it")
    rate = reference * math.cos(central)
    if not (math.isfinite(rate) and rate > 0.0):
        raise ValueError(
            f"the derived contraction rate is {rate!r}; the width-only control "
            f"is defined by having exactly one strictly positive rate, so a "
            f"non-finite or non-positive derivation is refused here rather "
            f"than becoming a control that contracts nowhere")
    return rate


#: The rate itself, computed from the rule above at import time so that it is a
#: constant of the module rather than something a experiment passes in.
#: Provisional in the same sense as :data:`REFERENCE_PEAK_COUPLING_RANGE`.
FIXED_CONTRACTION_RATE = fixed_contraction_rate(*REFERENCE_PEAK_COUPLING_RANGE)


# ---------------------------------------------------------------------------
# 2. The coupling source and the exact tongue schedule
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class PulseTrain:
    """One or more non-overlapping raised-cosine pulses, in time order.

    A single-element train is the ordinary case and behaves exactly like the
    calibrated :class:`~.model.RaisedCosinePulse` it contains.  Two or more
    exist for one reason: **re-entry**.  A clock that enters a tongue, exits it,
    and enters again is the case where the width-only control's coordinate
    handoff can silently lose or gain a turn, and a single pulse cannot produce
    one.  The supports are required to be disjoint so that the coupling is
    unambiguous and the eligibility windows of any clock are separated by a
    genuinely ineligible gap.

    Because the supports are disjoint and the raised-cosine envelope is exactly
    zero outside its own support, the train's coupling is the *sum* of its
    pulses' couplings — at most one term is ever nonzero.  Summing rather than
    selecting keeps the expression vectorized and keeps the endpoints exact.
    """

    pulses: Tuple[RaisedCosinePulse, ...]

    def __post_init__(self):
        if isinstance(self.pulses, RaisedCosinePulse):
            raise TypeError("pulses must be a sequence of RaisedCosinePulse "
                            "objects, not a single pulse")
        try:
            ordered = tuple(self.pulses)
        except TypeError:
            raise TypeError(f"pulses must be a sequence of RaisedCosinePulse "
                            f"objects, got "
                            f"{type(self.pulses).__name__}") from None
        if not ordered:
            raise ValueError("a pulse train must contain at least one pulse")
        for position, pulse in enumerate(ordered):
            require_exact_type(pulse, RaisedCosinePulse, f"pulses[{position}]")
        for first, second in zip(ordered, ordered[1:]):
            if not first.support()[1] <= second.support()[0]:
                raise ValueError(
                    f"pulse supports must be disjoint and in time order: "
                    f"{first.support()!r} overlaps or precedes "
                    f"{second.support()!r}; overlapping supports would make "
                    f"the coupling at a shared instant depend on which pulse "
                    f"was asked")
        object.__setattr__(self, "pulses", ordered)

    def _coupling_unchecked(self, t):
        """Unvalidated kernel behind :meth:`coupling`; the stepping hot path.

        Named to match the calibrated coupling sources so that a pulse train
        can be handed to ``simulate.integrate`` and to ``model._step_euler``
        wherever one of those is used.
        """
        total = np.zeros_like(np.asarray(t, dtype=float))
        for pulse in self.pulses:
            total = total + pulse._coupling_unchecked(t)
        return total

    def coupling(self, t):
        """The coupling at time(s) ``t``: zero outside every pulse's support."""
        return self._coupling_unchecked(require_finite_array(t, "t"))[()]

    def support(self) -> Tuple[float, float]:
        """The union interval, from the first pulse's start to the last's end."""
        return (self.pulses[0].support()[0], self.pulses[-1].support()[1])

    @property
    def peak_coupling(self) -> float:
        """The largest coupling this drive ever reaches."""
        return max(pulse.peak for pulse in self.pulses)

    @property
    def drive_label(self) -> str:
        """Which of :data:`DRIVE_LABELS` this is, for a record to carry."""
        return "pulsed"

    @property
    def digest(self) -> str:
        """The pulse manifest's hash, as it appears inside a crossing identity.

        Taken over the schema tag and every pulse's peak, duration and centre
        written as hexadecimal floats, which round-trip exactly.  Decimal
        formatting would let two distinguishable trains hash to one digest, and
        a crossing identity that cannot distinguish two pulse geometries would
        let two different experiments share a split normal.
        """
        parts = [_PULSE_SCHEMA]
        for pulse in self.pulses:
            parts.append(f"{pulse.peak.hex()}:{pulse.duration.hex()}:"
                         f"{pulse.center.hex()}")
        payload = "|".join(parts).encode("ascii")
        return hashlib.sha256(payload).hexdigest()[:PULSE_DIGEST_CHARS]

    def schedule(self, detuning) -> "TongueSchedule":
        """The exact eligibility windows of a clock at this mismatch.

        One window per pulse that admits the clock at all, taken from the
        calibrated ``RaisedCosinePulse.eligibility_window`` — the continuous
        solution of ``K(t) = |Delta|``, not a sampled detection.  A clock the
        train never admits gets an empty schedule, which is a valid schedule
        and evolves under the full Adler drift throughout.
        """
        mismatch = require_finite(detuning, "detuning")
        windows = []
        for pulse in self.pulses:
            window = pulse.eligibility_window(mismatch)
            if window is not None:
                windows.append(window)
        return TongueSchedule(tuple(windows))


@dataclass(frozen=True)
class StationaryDrive:
    """Constant coupling for the whole run: the tongue never opens or closes.

    The other half of experiment S4, which the governing plan states as "run
    one clock under **stationary and pulsed** coupling".  A stationary clock is
    not a very long pulse: it has no rising edge, no falling edge, no finite
    eligibility window and — the part that matters for the machinery here — **no
    tongue crossings at all**.  Its eligibility is decided once, by
    ``peak > |Delta|``, and never changes.

    That is why this is a separate drive rather than a `PulseTrain` with a wide
    support.  Manufacturing two finite crossing times far outside the run would
    put two fictitious entries in the Brownian tree's key schedule, split two
    steps that nothing physically crosses, and describe an eligibility change
    that does not occur.  A stationary clock's schedule is
    :class:`OpenSchedule`, which has no crossings to key and no exit to hand
    off at.

    What is *not* here: any stationary-hazard analysis, burn-in semantics,
    left-truncation, exponential-survival test or killed-diffusion oracle.
    Those are experiments S3 and S6 and belong to later tickets.  This is one
    clock, one constant coupling, one record.
    """

    peak: float

    def __post_init__(self):
        peak = require_finite(self.peak, "peak")
        if peak < 0.0:
            raise ValueError(f"peak coupling must be non-negative, got "
                             f"{peak!r}")
        object.__setattr__(self, "peak", peak)

    def _coupling_unchecked(self, t):
        """Unvalidated kernel behind :meth:`coupling`; the stepping hot path."""
        return np.broadcast_to(np.float64(self.peak),
                               np.shape(np.asarray(t, dtype=float)))

    def coupling(self, t):
        """The coupling at time(s) ``t``; constant, but ``t`` is still checked."""
        return self._coupling_unchecked(require_finite_array(t, "t"))[()]

    @property
    def peak_coupling(self) -> float:
        """The largest coupling this drive ever reaches."""
        return self.peak

    @property
    def drive_label(self) -> str:
        """Which of :data:`DRIVE_LABELS` this is, for a record to carry."""
        return "stationary"

    @property
    def digest(self) -> str:
        """The drive manifest's hash, on the same scheme as a pulse train's.

        A stationary drive contributes no crossing and therefore no key, so
        this never reaches the Brownian tree.  It exists so that a record
        naming a stationary run identifies which stationary run.
        """
        payload = f"{_STATIONARY_SCHEMA}|{self.peak.hex()}".encode("ascii")
        return hashlib.sha256(payload).hexdigest()[:PULSE_DIGEST_CHARS]

    def schedule(self, detuning) -> "TongueSchedule":
        """Eligible for all time, or for none of it.

        ``OpenSchedule()`` when ``peak > |Delta|`` and the empty
        ``TongueSchedule(())`` otherwise.  There is no third possibility and no
        boundary state to reach: the equality ``peak == |Delta|`` is a clock
        that is never strictly eligible, which the empty schedule already says.
        """
        mismatch = require_finite(detuning, "detuning")
        if self.peak > abs(mismatch):
            return OpenSchedule()
        return TongueSchedule(())


@dataclass(frozen=True)
class OpenSchedule:
    """A clock that is strictly eligible at every instant of its run.

    The stationary counterpart of :class:`TongueSchedule`, exposing the same
    three questions.  It answers all of them without a window, because it has
    none: there is nothing to enter, nothing to leave, and no equality instant
    where ``stable_phase`` would be undefined.

    It is a separate type rather than a flag on ``TongueSchedule`` so that
    "eligible everywhere" cannot be spelled as a window with infinite
    endpoints — which would then be routed into the mesh, keyed, and split.
    """

    def state_at(self, time) -> str:
        """Always ``"interior"``: a stationary tongue has no boundary."""
        require_finite(time, "time")
        return "interior"

    def interval_state(self, start, end) -> str:
        """Always ``"interior"``; no interval can straddle a crossing."""
        first = require_finite(start, "start")
        last = require_finite(end, "end")
        if not last > first:
            raise ValueError(f"an interval state needs positive width, got "
                             f"[{first!r}, {last!r}]")
        return "interior"

    def ends_at_exit(self, time) -> bool:
        """Always ``False``: there is no exit for an interval to end at."""
        require_finite(time, "time")
        return False

    def crossings(self, clock, pulse_digest) -> Tuple[Crossing, ...]:
        """No crossings at all — the point of the type.

        The arguments are still validated, so that a caller cannot discover
        later that a malformed clock identifier was silently accepted by the
        stationary path and rejected by the pulsed one.
        """
        require_positive_int(clock, "clock", minimum=0)
        _require_digest_token(pulse_digest, "pulse_digest")
        return ()


@dataclass(frozen=True)
class TongueSchedule:
    """One clock's exact tongue entry and exit times, in order.

    The windows are open intervals: ``entry`` and ``exit`` themselves are the
    equality state, where ``stable_phase`` is NaN and ``boundary_phase`` is the
    finite limit.  Windows must be strictly separated, so that "leaves and
    re-enters" is a statement about two windows with an ineligible gap rather
    than about a single window with a seam in it.
    """

    windows: Tuple[Tuple[float, float], ...]

    def __post_init__(self):
        try:
            ordered = tuple(tuple(window) for window in self.windows)
        except TypeError:
            raise TypeError(f"windows must be a sequence of (entry, exit) "
                            f"pairs, got "
                            f"{type(self.windows).__name__}") from None
        checked = []
        for position, window in enumerate(ordered):
            if len(window) != 2:
                raise ValueError(f"windows[{position}] must be an "
                                 f"(entry, exit) pair, got {len(window)} "
                                 f"values")
            entry = require_finite(window[0], f"windows[{position}] entry")
            exit_time = require_finite(window[1], f"windows[{position}] exit")
            if not exit_time > entry:
                raise ValueError(f"windows[{position}] spans "
                                 f"[{entry!r}, {exit_time!r}], which has no "
                                 f"width; a tongue window a clock is in for no "
                                 f"time is not a window")
            checked.append((entry, exit_time))
        for first, second in zip(checked, checked[1:]):
            if not first[1] < second[0]:
                raise ValueError(
                    f"tongue windows must be separated by an ineligible gap: "
                    f"{first!r} touches or overlaps {second!r}")
        object.__setattr__(self, "windows", tuple(checked))

    def state_at(self, time) -> str:
        """Which of :data:`ELIGIBILITY_STATES` holds at ``time``.

        Equality is tested against the exact crossing times, before the strict
        interior test, so an instant that *is* a crossing is reported as
        ``boundary`` and never as an interior point one rounding inside.
        """
        moment = require_finite(time, "time")
        for entry, exit_time in self.windows:
            if moment == entry or moment == exit_time:
                return "boundary"
        for entry, exit_time in self.windows:
            if entry < moment < exit_time:
                return "interior"
        return "outside"

    def _window_of(self, start: float, end: float) -> Optional[int]:
        """Index of the window containing the open interval, or ``None``.

        Raises if the interval straddles a crossing.  That is the enforcement
        point of "no step spans an eligibility change": the noise tree splits
        every step at its exact crossings, so an interval that still straddles
        one means the schedule handed to the stepper is not the schedule the
        tree was built from, and evolving it would apply one drift law across a
        boundary the noise had already been divided at.
        """
        for index, (entry, exit_time) in enumerate(self.windows):
            if end <= entry or start >= exit_time:
                continue
            if entry <= start and end <= exit_time:
                return index
            raise ValueError(
                f"the interval [{start!r}, {end!r}] straddles the tongue "
                f"window [{entry!r}, {exit_time!r}]; every step must be split "
                f"at its exact eligibility crossing before it is evolved")
        return None

    def interval_state(self, start, end) -> str:
        """``"interior"`` or ``"outside"`` for a whole interval of positive width.

        Never ``"boundary"``: an interval is not an instant.  A crossing at one
        of its endpoints leaves the interval wholly on one side, which is what
        makes an exactly split step unambiguous to evolve.
        """
        first = require_finite(start, "start")
        last = require_finite(end, "end")
        if not last > first:
            raise ValueError(f"an interval state needs positive width, got "
                             f"[{first!r}, {last!r}]; a zero-duration handoff "
                             f"is an instant and uses state_at")
        return "outside" if self._window_of(first, last) is None else "interior"

    def ends_at_exit(self, time) -> bool:
        """Does an eligible interval ending here end at an exact tongue exit?

        The one question the endpoint-target branch asks, and the reason it is
        a method on the schedule rather than an index lookup in the stepper:
        :class:`OpenSchedule` answers it too, with ``False``, and a stationary
        clock therefore needs no special case in
        :meth:`ClockPath._endpoint_target`.
        """
        moment = require_finite(time, "time")
        return any(moment == exit_time for _, exit_time in self.windows)

    def crossings(self, clock, pulse_digest) -> Tuple[Crossing, ...]:
        """The frozen crossing schedule this clock hands to the noise tree.

        One ``rising`` crossing per entry and one ``falling`` crossing per
        exit, carrying the pulse digest so that two different pulse geometries
        cannot key one split normal.  The list is exactly what
        ``PhaseNoiseStream.elementary_leaves`` needs, and passing anything else
        would mean the mesh the phase is evolved on and the mesh the noise was
        divided on are two meshes.
        """
        index = require_positive_int(clock, "clock", minimum=0)
        digest = pulse_digest
        out = []
        for entry, exit_time in self.windows:
            out.append(Crossing(index, "rising", digest, entry))
            out.append(Crossing(index, "falling", digest, exit_time))
        return tuple(out)


#: The coupling histories :class:`ClockPath` accepts, checked by *exact* type.
#: Declared here rather than inline so that adding a third drive is one visible
#: edit and cannot be done by subclassing one of these two.
DRIVE_TYPES = (PulseTrain, StationaryDrive)


# ---------------------------------------------------------------------------
# 3. The evolving state, and what an endpoint looks like
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class PathState:
    """One clock's state at one endpoint.

    ``phase`` is the authoritative unwrapped absolute phase and is always
    present.  ``target_lift`` and ``error`` are the width-only control's
    auxiliary pair and are present only while it is strictly eligible; the pair
    is NaN together, which is the documented "no auxiliary target" sentinel and
    the same convention ``model.stable_phase`` uses for "no fixed point".

    They are validated as a pair rather than individually, because a half-live
    pair is the failure this representation exists to prevent: a lift without
    an error, or an error without a lift, is a state from which the
    reconstruction ``phase = lift + error`` cannot be performed and which would
    otherwise be discovered several steps later as a NaN phase.
    """

    time: float
    phase: float
    target_lift: float = math.nan
    error: float = math.nan

    def __post_init__(self):
        moment = require_finite(self.time, "time")
        absolute = _require_resolvable(self.phase, "phase")
        for name in ("target_lift", "error"):
            raw = getattr(self, name)
            if isinstance(raw, (bool, str, bytes, bytearray, complex)):
                raise TypeError(f"{name} must be a real number, got "
                                f"{type(raw).__name__} {raw!r}")
        try:
            lift = float(self.target_lift)
            error = float(self.error)
        except (TypeError, ValueError):
            raise TypeError("target_lift and error must be real numbers, got "
                            f"{type(self.target_lift).__name__} and "
                            f"{type(self.error).__name__}") from None
        if math.isinf(lift) or math.isinf(error):
            raise ValueError("target_lift and error may be NaN together — the "
                             "no-auxiliary-target sentinel — but never "
                             "infinite")
        if math.isnan(lift) != math.isnan(error):
            raise ValueError(
                f"target_lift and error must be live together or absent "
                f"together, got target_lift={lift!r}, error={error!r}; a "
                f"half-live pair cannot reconstruct a phase")
        if not math.isnan(lift):
            _require_resolvable(lift, "target_lift")
        object.__setattr__(self, "time", moment)
        object.__setattr__(self, "phase", absolute)
        object.__setattr__(self, "target_lift", lift)
        object.__setattr__(self, "error", error)

    @property
    def has_target(self) -> bool:
        """Is the auxiliary lifted target live at this endpoint?"""
        return not math.isnan(self.target_lift)

    @property
    def lift_residual(self) -> float:
        """``(target_lift + error) - phase``, or ``0.0`` with no target.

        The identity ``phase = target_lift + error`` is exact in real
        arithmetic and is *asserted* rather than assumed, so it needs a number.
        It is bit-exact wherever the phase was assigned as that sum — every
        eligible step reconstructs it that way — and off by at most a rounding
        at an entry handoff, where the phase is deliberately left untouched and
        the pair is derived from it instead.  Compare against
        ``stochastic.conservation_tolerance``, which is the same kind of claim
        about the same kind of arithmetic.
        """
        if not self.has_target:
            return 0.0
        return (self.target_lift + self.error) - self.phase


@dataclass(frozen=True)
class Observation:
    """What one endpoint looks like to the dwell state machine.

    Deliberately not a phase.  ``commitment.py`` receives only ``eligible``,
    ``lock_error`` and ``contraction`` from this, which is why it can be driven
    by declared synthetic sequences: everything it needs is here, and nothing
    it needs is a trajectory.

    ``target``, ``lock_error`` and ``contraction`` are NaN wherever no target
    exists.  ``eligible`` is true only in the strict interior — an endpoint
    exactly at a crossing reports a finite boundary target and an ineligible
    flag, so the equality state is visible in the record and still earns no
    dwell.
    """

    time: float
    state: str
    coupling: float
    target: float
    lock_error: float
    contraction: float
    eligible: bool

    def __post_init__(self):
        if self.state not in ELIGIBILITY_STATES:
            raise ValueError(f"state must be one of {list(ELIGIBILITY_STATES)}, "
                             f"got {self.state!r}")
        if not isinstance(self.eligible, bool):
            raise TypeError(f"eligible must be a bool, got "
                            f"{type(self.eligible).__name__}")
        moment = require_finite(self.time, "time")
        strength = require_finite(self.coupling, "coupling")
        if strength < 0.0:
            raise ValueError(f"coupling must be non-negative, got {strength!r}")
        values = {}
        for name in ("target", "lock_error", "contraction"):
            raw = getattr(self, name)
            if isinstance(raw, (bool, str, bytes, bytearray, complex)):
                raise TypeError(f"{name} must be a real number, got "
                                f"{type(raw).__name__}")
            try:
                value = float(raw)
            except (TypeError, ValueError):
                raise TypeError(f"{name} must be a real number, got "
                                f"{type(raw).__name__}") from None
            if math.isinf(value):
                raise ValueError(f"{name} may be NaN where no target exists, "
                                 f"but never infinite")
            values[name] = value
        if self.eligible and any(math.isnan(v) for v in values.values()):
            raise ValueError("an eligible endpoint has a finite target, error "
                             "and contraction; NaN there is a numerical "
                             "failure wearing the missing-target sentinel")
        object.__setattr__(self, "time", moment)
        object.__setattr__(self, "coupling", strength)
        for name, value in values.items():
            object.__setattr__(self, name, value)


# ---------------------------------------------------------------------------
# 4. The stepper
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class ClockPath:
    """One clock's data-generating process: mismatch, drive, model, rate.

    Frozen, and the schedule is derived once in the constructor, because the
    exact crossings are what the noise tree was built from and re-deriving them
    per step is an invitation for the two to drift apart.

    ``train`` is either a :class:`PulseTrain` or a :class:`StationaryDrive` —
    the two halves of experiment S4 — and the type must be *exactly* one of
    them.  Not an ``isinstance`` check: a subclass of either could override
    ``_coupling_unchecked`` or ``schedule`` to return something computed from
    the prediction, which is the same hole a subclassed configuration would be.
    The field keeps its original name so that existing callers and records are
    unchanged.

    ``model`` selects between the calibrated full Adler process and the
    fixed-contraction width-only control.

    **``fixed_rate`` has exactly one legal value per model**: exactly ``0.0``
    for ``"full"`` and exactly :data:`FIXED_CONTRACTION_RATE` for
    ``"width_only"``.  It stays an argument rather than becoming derived so
    that a run configuration says out loud which rate it used, but there is
    nothing to choose.

    That is stricter than it was, and the strictness is the point.  Closing
    ``raw_runner.raw_one_clock_path`` to per-run rates was necessary and not
    sufficient: this class is exported and ``raw_experiments.run_one_clock``
    accepts it, so an independent review ran three different width-only
    trajectories under one ordinary ``S4/stationary/interior/zero/width_only``
    label at rates ``0.1``, ``1.0000000000000002`` and ``7.0``, each producing
    an ordinary record.  The control's whole causal claim is that *only*
    coupling-dependent interior contraction was removed; a rate that varies per
    run reintroduces the factor under the control's own name.

    No sensitivity API replaces it.  Lower and upper fixed-rate choices are a
    later ticket's labelled sensitivity against a frozen production manifest,
    and would need their own type and a record identity that survives into the
    ledger; an unlabelled escape hatch now would be the same bypass with a
    longer signature.
    """

    detuning: float
    train: object
    model: str = "full"
    fixed_rate: float = 0.0

    def __post_init__(self):
        mismatch = require_finite(self.detuning, "detuning")
        if type(self.train) not in DRIVE_TYPES:
            raise TypeError(
                f"train must be exactly one of "
                f"{[t.__name__ for t in DRIVE_TYPES]}, got "
                f"{type(self.train).__name__}")
        if self.model not in MODEL_LABELS:
            raise ValueError(f"model must be one of {list(MODEL_LABELS)}, got "
                             f"{self.model!r}")
        rate = require_finite(self.fixed_rate, "fixed_rate")
        expected = (FIXED_CONTRACTION_RATE if self.model == "width_only"
                    else 0.0)
        if rate != expected:
            raise ValueError(
                f"the {self.model!r} process must be given a fixed contraction "
                f"rate of exactly {expected!r}, got {rate!r}.  "
                f"The width-only control is *defined* by having one rate — "
                f"common to every clock, coupling, mismatch, trial, dwell, "
                f"tolerance and noise strength — so a path carrying any other "
                f"positive value is a different experiment wearing the "
                f"control's label, and the full Adler model has no such rate "
                f"at all.  Closing the raw factory alone was not enough: this "
                f"constructor is exported, so it was the remaining public "
                f"route to a per-run rate.  Pass "
                f"dynamics.FIXED_CONTRACTION_RATE, which is derived once by "
                f"fixed_contraction_rate() from the declared reference range")
        object.__setattr__(self, "detuning", mismatch)
        object.__setattr__(self, "model", self.model)
        object.__setattr__(self, "fixed_rate", rate)
        object.__setattr__(self, "_schedule", self.train.schedule(mismatch))

    @property
    def _cached_schedule(self):
        """The raw cached schedule, unvalidated.

        Private, and the *only* accessor the trusted kernels use.  The public
        :attr:`schedule` rebuilds the whole path before handing one back, so an
        internal caller reaching through it would rebuild the pulse train and
        re-derive the geometry on every step — and
        :func:`validate_clock_path`, which reads this to compare against the
        schedule it derives, would recurse forever.
        """
        try:
            return self.__dict__["_schedule"]
        except KeyError:
            raise ValueError(
                "this path has no cached tongue schedule, so it was never "
                "built by ClockPath's own constructor") from None

    @property
    def schedule(self) -> TongueSchedule:
        """This clock's exact tongue windows, derived once at construction.

        Validated: the path is rebuilt and its cached geometry checked against
        what its own fields derive before the schedule is returned, so a
        caller cannot read the tongue windows of one clock off an object whose
        detuning says another.  See :func:`validate_clock_path`.
        """
        return validate_clock_path(self)._cached_schedule

    def _crossings_trusted(self) -> Tuple[Crossing, ...]:
        """Unvalidated kernel behind :meth:`crossings`."""
        return self._cached_schedule.crossings(0, self.train.digest)

    def crossings(self) -> Tuple[Crossing, ...]:
        """The crossing schedule for clock zero under this train's digest.

        A convenience for the single-clock experiments: a population run keys
        each clock by its own identifier, and passes that identifier to
        ``TongueSchedule.crossings`` directly.

        Validated, because this is the first door of a direct walk: the
        crossings are what the Brownian tree is split at, so a stale one would
        divide the noise at another clock's tongue edges.
        """
        return validate_clock_path(self)._crossings_trusted()

    # -- targets ---------------------------------------------------------
    def _interior_target(self, time: float) -> float:
        """``stable_phase`` at a strictly interior time, or a hard failure.

        The strict-eligibility guard is re-asserted here against the *coupling*
        even though the caller has already established strict interiority from
        the exact crossing times.  The two can disagree within a rounding of a
        crossing, and when they do, ``stable_phase`` returns NaN — which would
        seed a NaN target lift, propagate into the phase, and surface later as
        a clock that had quietly stopped existing.  Failing here names the
        cause.
        """
        strength = float(self.train._coupling_unchecked(time))
        if not strength > abs(self.detuning):
            raise ValueError(
                f"t={time!r} was classified strictly inside the tongue from "
                f"the exact crossings, but the coupling there is "
                f"{strength!r} <= |detuning|={abs(self.detuning)!r}, so "
                f"stable_phase is not defined; the pulse and the crossing "
                f"schedule disagree at this instant")
        return float(stable_phase(strength, self.detuning))

    def _endpoint_target(self, time: float) -> float:
        """The principal target at the end of an eligible interval.

        ``boundary_phase`` if the interval ends exactly at a tongue exit,
        ``stable_phase`` otherwise.  This is the branch the plan requires and
        the reason the equality state is never handed to ``stable_phase``: the
        last eligible step of a window ends *on* the crossing, and asking for a
        stable phase there returns NaN by design.

        The question goes to the schedule rather than to a window index, so a
        stationary clock — whose :class:`OpenSchedule` has no exit — takes the
        interior branch every time without a special case here.
        """
        if self._cached_schedule.ends_at_exit(time):
            return float(boundary_phase(self.detuning))
        return self._interior_target(time)

    # -- state -----------------------------------------------------------
    def initial_state(self, time, phase) -> PathState:
        """The state a run starts from.

        Validated wrapper over :meth:`_initial_state_trusted`: the path is
        rebuilt before its cached geometry decides whether this start is
        inside a window.
        """
        return validate_clock_path(self)._initial_state_trusted(time,
                                                                phase)

    def _initial_state_trusted(self, time, phase) -> PathState:
        """The state a run starts from.

        Validated: the path is rebuilt before its cached geometry decides
        whether this start is inside a window.  See
        :func:`validate_clock_path`.


        The auxiliary pair is established here only when the run begins
        *strictly inside* a window and the model is the width-only control —
        that is a prepared initial condition, not a handoff, and it is the only
        way the pair can appear without an exact entry crossing.  A run that
        starts outside, or exactly at a crossing, starts with no pair and picks
        one up at its first eligible interval.
        """
        moment = require_finite(time, "time")
        absolute = _require_resolvable(phase, "phase")
        if self.model != "width_only":
            return PathState(moment, absolute)
        if self._cached_schedule.state_at(moment) != "interior":
            return PathState(moment, absolute)
        lift, error = _lift_and_error(self._interior_target(moment), absolute)
        return PathState(moment, absolute, lift, error)

    def entry_handoff(self, state) -> PathState:
        """The zero-duration coordinate handoff at an exact tongue entry.

        Validated wrapper over :meth:`_entry_handoff_trusted`.
        """
        path = validate_clock_path(self)
        checked = _validated_path_state(path, state)
        # The edge precondition is asserted here, in the wrapper, so that a
        # wrong-edge call refuses *before* the trusted kernel runs.  The
        # kernel keeps its own copy for the internal handoffs
        # ``_advance_trusted`` performs on states it has just built.
        path._require_crossing(checked.time, "entry")
        return path._entry_handoff_trusted(checked)

    def _entry_handoff_trusted(self, state) -> PathState:
        """The zero-duration coordinate handoff at an exact tongue entry.

        Validated at both ends: the path is rebuilt, and the supplied state
        is rebuilt through ``PathState``'s own constructor, so neither a
        stale tongue geometry nor a half-live auxiliary pair can reach the
        handoff.


        No time passes, no noise is consumed, and **the authoritative phase is
        returned unchanged** — bit for bit, not to a tolerance.  All that
        happens is that a copy of the finite ``boundary_phase(Delta)`` is
        lifted onto the branch nearest the phase, and the unwrapped error from
        it is stored.

        ``boundary_phase`` and never ``stable_phase``: at the entry instant the
        coupling equals the mismatch exactly, the attracting and repelling
        fixed points have merged, and ``stable_phase`` returns NaN there by
        design.  A NaN seeded into the lift here would reappear several steps
        later as a phase that had quietly stopped existing, which is why the
        verifier mutates this call and requires the mutation to fail.

        The reconstruction ``lift + error`` equals the phase to within one
        rounding rather than bit-exactly, because the phase is deliberately
        *not* reassigned here.  Every subsequent eligible step assigns the
        phase as that sum, so the identity is exact from the first step
        onwards; see :attr:`PathState.lift_residual`.
        """
        if self.model != "width_only":
            raise ValueError(f"only the width-only control carries an "
                             f"auxiliary target; this path is {self.model!r}")
        if state.has_target:
            raise ValueError("this state already carries a lifted target; a "
                             "second entry handoff would discard a live error "
                             "and its accumulated winding")
        self._require_crossing(state.time, "entry")
        lift, error = _lift_and_error(float(boundary_phase(self.detuning)),
                                      state.phase)
        return PathState(state.time, state.phase, lift, error)

    def exit_handoff(self, state) -> PathState:
        """The zero-duration coordinate handoff at an exact tongue exit.

        Validated wrapper over :meth:`_exit_handoff_trusted`.
        """
        path = validate_clock_path(self)
        checked = _validated_path_state(path, state)
        path._require_crossing(checked.time, "exit")
        return path._exit_handoff_trusted(checked)

    def _exit_handoff_trusted(self, state) -> PathState:
        """The zero-duration coordinate handoff at an exact tongue exit.

        Validated at both ends, as :meth:`entry_handoff` is.


        The mirror image, and deliberately asymmetric in what it touches: the
        eligible segment that ended here has *already* reconstructed the
        authoritative phase against the lifted boundary target, so there is
        nothing left to convert.  This discards the auxiliary pair and returns
        the same phase, again bit for bit.

        There is no coordinate reset, no wrap, and no call to ``stable_phase``.
        A re-entry later builds a fresh lift from this unchanged phase, which
        is what makes an entry/exit/re-entry cycle preserve the accumulated
        winding instead of quietly returning the clock to its principal value.
        """
        if not state.has_target:
            raise ValueError("this state carries no lifted target, so there is "
                             "nothing to hand back")
        self._require_crossing(state.time, "exit")
        return PathState(state.time, state.phase)

    def advance(self, state, leaf, end_time) -> PathState:
        """Evolve one elementary leaf and return the endpoint state.

        Validated wrapper.  This is the transition itself, and a stale
        path does not merely mis-*report* here — it takes the wrong
        branch.  A width-only path whose detuning had been edited from
        ``0.9`` to ``0.1`` kept the ``0.9`` tongue window, so an interval
        the real clock is strictly eligible on was stepped with the full
        Adler drift instead of the fixed-contraction one, and the phase
        update itself was different.  The path is therefore rebuilt and
        the supplied state revalidated before any of that is decided.

        All three frozen inputs are rebuilt: the path, the supplied state,
        and the leaf itself.  The leaf carries the kick *and its address*
        — which uniform step it belongs to and which node of the Brownian
        tree carved it out — and an edited address used to step the phase
        without complaint; see :func:`_validated_tree_leaf`.

        A caller stepping a whole run should validate once and use the
        private kernel; ``raw_experiments.run_one_clock`` does, so it does
        not rebuild the pulse train, re-derive the schedule, or rebuild a
        leaf its own validated stream produced.
        See :meth:`_advance_trusted` for the four cases.
        """
        path = validate_clock_path(self)
        return path._advance_trusted(_validated_path_state(path, state),
                                     _validated_tree_leaf(leaf), end_time)

    def _advance_trusted(self, state, leaf, end_time) -> PathState:
        """Evolve one elementary leaf and return the endpoint state.

        ``leaf`` is a :class:`~.stochastic.TreeLeaf` — one interval of the
        elementary physical mesh, carrying the kick that belongs to it.
        ``end_time`` is that interval's exact right endpoint, supplied by the
        caller because ``leaf.start_time + leaf.duration`` is a *different
        float* from the boundary the tree was actually built at, and every
        equality test in this method is against an exact crossing.  The two are
        cross-checked: ``end_time - start_time`` must equal ``leaf.duration``
        bit for bit, which is the same subtraction the tree performed, so a
        caller that supplies the wrong endpoint is refused rather than
        integrated.

        Four cases:

        1. **Zero duration.**  The handoff leaf an exact crossing produces when
           it lands on a step boundary.  No time passes, so no phase changes
           and no noise is consumed; its kick is required to be exactly zero.
        2. **Ineligible interval, either model.**  The full calibrated
           Euler-Maruyama step, drift frozen at the interval start.  The
           control takes this branch too, and takes it through the same
           function, which is what makes "identical while ineligible" a fact
           about shared code rather than about two expressions that look alike.
        3. **Eligible interval, full model.**  The same step.  Eligibility
           changes what a *target* is, not what the drift is.
        4. **Eligible interval, width-only control.**  Entry handoff if the
           pair is not yet live, then the fixed-rate restoring step on the
           unwrapped error, the target continued onto its lifted branch, and
           the phase reconstructed as their sum.

        The exit handoff happens at the start of the first ineligible interval
        after a window, which is where the plan puts it: the eligible segment
        that *ends* at the exit has already reconstructed the phase against the
        lifted boundary target, so the handoff discards the pair and changes
        nothing else.
        """
        end = require_finite(end_time, "end_time")
        if leaf.start_time != state.time:
            raise ValueError(
                f"the leaf starts at t={leaf.start_time!r} and the state is at "
                f"t={state.time!r}; a walk that skipped or repeated an "
                f"elementary interval would consume the wrong kick")
        span = end - leaf.start_time
        if span < 0.0:
            raise ValueError(f"end_time {end!r} precedes the leaf start "
                             f"{leaf.start_time!r}")
        if span != leaf.duration:
            raise ValueError(
                f"end_time {end!r} gives a span of {span!r} for a leaf of "
                f"duration {leaf.duration!r}; the endpoint must be the exact "
                f"boundary the Brownian tree was built at, not a recomputed "
                f"one")
        if self.model == "full" and state.has_target:
            raise ValueError("the full Adler model carries no auxiliary "
                             "target; this state belongs to the width-only "
                             "control")

        if span == 0.0:
            if leaf.kick != 0.0:
                raise ValueError(
                    f"a zero-duration handoff leaf carries no variance, so its "
                    f"kick must be exactly zero, got {leaf.kick!r}")
            return state

        eligible = self._cached_schedule.interval_state(leaf.start_time,
                                                       end) == "interior"

        if not eligible or self.model == "full":
            if state.has_target:
                state = self._exit_handoff_trusted(state)
            return PathState(end, self._adler_step(state.phase,
                                                   leaf.start_time, span,
                                                   leaf.kick))

        if not state.has_target:
            state = self._entry_handoff_trusted(state)
        lift, error = state.target_lift, state.error

        lift_next = continuous_lift(self._endpoint_target(end), lift)
        error_next = (error - self.fixed_rate * math.sin(error) * span
                      + leaf.kick)
        return PathState(end, lift_next + error_next, lift_next, error_next)

    def _require_crossing(self, time: float, what: str) -> None:
        """A coordinate handoff may only happen at its **own** exact crossing.

        ``what`` is ``"entry"`` or ``"exit"``, and it is load-bearing rather
        than decorative.  This used to test only that the instant was *a*
        boundary, so with a single window ``(0.41, 3.59)`` an entry handoff was
        accepted at the exit and an exit handoff at the entry — each performing
        the opposite coordinate operation at the wrong edge of the tongue.
        Both preserved the phase, which is why nothing downstream complained,
        and both left the auxiliary pair in the state the other edge should
        have produced.

        A window is ``(entry, exit)``; an entry handoff must land on a first
        endpoint and an exit handoff on a second.  A schedule with no windows
        at all — a stationary :class:`OpenSchedule`, which is eligible
        everywhere and has nothing to cross — accepts neither, which is
        correct: a stationary clock's pair is prepared by
        :meth:`_initial_state_trusted` and never handed off.
        """
        edge = 0 if what == "entry" else 1
        windows = getattr(self._cached_schedule, "windows", ())
        if not any(time == window[edge] for window in windows):
            wanted = [window[edge] for window in windows]
            raise ValueError(
                f"a tongue {what} handoff was reached at t={time!r}, which is "
                f"not one of this clock's exact {what} crossings {wanted}; the "
                f"handoff is directed, so the opposite edge is refused as "
                f"firmly as an ordinary interior instant, and no dwell or "
                f"coordinate change may begin before exact eligibility")

    def _adler_step(self, phase: float, time: float, span: float,
                    kick: float) -> float:
        """One Euler-Maruyama step of the full calibrated Adler drift.

        ``model._step_euler`` supplies the deterministic part, called rather
        than re-derived, so that "zero noise is exactly deterministic Euler" is
        true by construction: with ``kick == 0.0`` this returns the calibrated
        stepper's own float, unmodified.
        """
        deterministic = float(_step_euler(phase, time, span, self.detuning,
                                          self.train._coupling_unchecked))
        return deterministic + kick

    # -- observation -----------------------------------------------------
    def observe(self, state) -> Observation:
        """What the dwell machine is told about this endpoint.

        Validated wrapper: the path is rebuilt and the supplied state is
        rebuilt through ``PathState``'s own constructor before the cached
        tongue geometry is consulted.  A stale path reported ``outside``
        and a NaN target for an endpoint the clock's own detuning puts
        strictly inside its tongue, and a half-live auxiliary pair was
        read as a live target.  See :meth:`_observe_trusted` for what is
        computed and :func:`validate_clock_path` for why.
        """
        path = validate_clock_path(self)
        return path._observe_trusted(_validated_path_state(path, state))

    def _observe_trusted(self, state) -> Observation:
        """What the dwell machine is told about this endpoint.

        Two evaluations, chosen by whether the auxiliary pair is live:

        * **auxiliary pair live** — the width-only control inside a window.
          The lock error is the *circular* image of the unwrapped error, so the
          error itself keeps its winding while the band is measured on the
          circle; the contraction is ``-rate cos(error)``, which is the
          derivative of the control's own restoring flow and is strict in
          exactly the way the full model's guard is.
        * **no pair** — the full model everywhere, and the control before entry
          and after exit, where it *is* the full model.  The target is
          ``stable_phase`` strictly inside, ``boundary_phase`` exactly at a
          crossing, and NaN outside; the contraction is ``-K cos(phase)``.

        ``eligible`` is strict interiority in both branches, so an endpoint
        sitting exactly on a crossing reports its finite boundary target and
        still earns no dwell.
        """
        where = self._cached_schedule.state_at(state.time)
        strength = float(self.train._coupling_unchecked(state.time))
        eligible = where == "interior"

        if state.has_target:
            target = state.target_lift
            lock_error = float(_wrap_phase(state.error))
            contraction = -self.fixed_rate * math.cos(state.error)
            return Observation(state.time, where, strength, target, lock_error,
                               contraction, eligible)

        if where == "interior":
            target = self._interior_target(state.time)
        elif where == "boundary":
            target = float(boundary_phase(self.detuning))
        else:
            return Observation(state.time, where, strength, math.nan, math.nan,
                               math.nan, eligible)
        lock_error = float(_circular_distance(state.phase, target))
        contraction = -strength * math.cos(state.phase)
        return Observation(state.time, where, strength, target, lock_error,
                           contraction, eligible)


# ---------------------------------------------------------------------------
# 5. The execution boundary: a trusted snapshot of a path
# ---------------------------------------------------------------------------
def _validated_tree_leaf(leaf) -> TreeLeaf:
    """Return a freshly rebuilt, fully revalidated copy of ``leaf``.

    A :class:`~.stochastic.TreeLeaf` is a frozen dataclass supplied by a
    caller, so it is under the same ``object.__setattr__`` threat model as the
    path and the state, and leaving it out would make the wrapper contract
    incomplete.  Rebuilt by exact type and every named field through Ticket
    02's own constructor, which is where a leaf's domain is defined and where
    it must be re-checked — this module does not restate those rules.

    Four of the six fields were already refused *incidentally*: a non-finite
    kick surfaced as a non-finite phase in ``PathState``, and an edited
    ``duration`` or ``start_time`` broke the span and start cross-checks in
    :meth:`ClockPath._advance_trusted`.  Incidental is not the same as
    refused-at-the-boundary — each of those failed partway through a
    transition rather than before it.  The other two were not refused at all:
    ``step_index`` and ``node_path`` are the leaf's *identity* — which uniform
    step this kick belongs to, and which node of the Brownian tree carved it
    out — and a leaf edited on either one stepped the phase happily under a
    rewritten address.

    Private: the public door is :meth:`ClockPath.advance`.
    ``raw_experiments.run_one_clock`` stays on the leaves its validated stream
    produced and does not rebuild per leaf.

    **Declared limitation.**  Rebuilding proves that every field lies in Ticket
    02's constructor domain — that the address is *canonical* — and no more.  A
    ``ClockPath`` has no mesh or stream context, so it cannot authenticate a
    different yet syntactically valid ``step_index`` or ``node_path`` as the
    one this kick was actually drawn under.  Provenance is the high-level
    runner's guarantee: ``run_one_clock`` accepts leaves only from its own
    validated ``PhaseNoiseStream``, which is where a leaf's address and its
    value are bound together.  This is recorded rather than repaired here;
    closing it at the primitive would mean giving the stepper a stream, which
    is a different ownership boundary than the one this package draws.
    """
    require_exact_type(leaf, TreeLeaf, "leaf")
    return TreeLeaf(leaf.start_time, leaf.duration, leaf.kick, leaf.node_path,
                    leaf.step_index)


def _validated_path_state(path, state) -> PathState:
    """Return a rebuilt copy of ``state``, validated **against ``path``**.

    Rebuilding through ``PathState``'s own constructor is necessary and was
    not sufficient.  That constructor knows about one state and no path: it
    checks finiteness, the resolvable-magnitude domain, and that the auxiliary
    pair is live-together or absent-together.  It cannot know whose target the
    lift is supposed to be a lift *of*, and an independent review used exactly
    that gap.

    On a stationary width-only path at mismatch ``0.35``::

        clean  = path.initial_state(0.0, 0.2)
        forged = PathState(clean.time, clean.phase,
                           clean.target_lift + 1.0, clean.error - 1.0)

    ``forged`` is constructor-valid and its reconstruction still holds to a
    rounding — ``(lift + error) - phase`` is ``-5.6e-17`` — yet its lift is a
    whole radian from the target this path actually has.  The compensating
    error changes ``sin(error)``, so the width-only transition changes: the
    next authoritative phase came out ``-0.708`` where the clean state gives
    ``0.216``.  A plainer edit, ``error += 1``, was accepted with a
    ``lift_residual`` of exactly ``1.0``.  Neither reached a key, and both
    reached and executed the transition.

    Four path-aware rules close it, checked here so that they run before any
    trusted kernel:

    1. **The reconstruction identity, to a declared bound.**  With a live pair,
       ``(target_lift + error) - phase`` must be within
       ``stochastic.conservation_tolerance`` of the magnitudes involved — the
       same scale-aware few-ulp rule the Brownian tree uses for parent
       conservation, and the honest bound: an eligible step *assigns* the phase
       as that sum and is exact, while an entry handoff deliberately leaves the
       phase untouched and is off by at most one rounding.
    2. **Pair ownership.**  Only the width-only control has an auxiliary pair.
       A full-model state carrying one was observed with that pair as its
       target, which is the control's coordinate system wearing the full
       model's label.
    3. **Lifecycle.**  A width-only state is *required* to carry the pair at a
       strictly interior time and *forbidden* to carry one outside.  At a
       crossing either is legitimate, because entry-before and exit-after are
       pair-less while entry-after and exit-before are live.  A strictly
       interior state without its pair was observed through the full-model
       branch, using ``-K cos(phase)`` where the control's
       ``-rate cos(error)`` was due.
    4. **Congruence.**  A live lift must be congruent, modulo a turn and within
       :data:`LIFT_RESOLUTION`, to the target the schedule says applies at that
       instant: ``stable_phase`` strictly inside, ``boundary_phase`` at a
       crossing.  This is what refuses the forged pair.

    **The unwrapped error keeps its semantics.**  Nothing here bounds it: the
    error may carry any number of whole turns, and only the *lift* is required
    to be congruent.  Wrapping it would destroy the winding these states exist
    to preserve.

    Private: it is the input half of a transition rather than something a
    caller needs.  The public doors are :meth:`ClockPath.advance`,
    :meth:`ClockPath.observe` and the two handoffs, and each of them rebuilds
    the path first and passes it here.
    """
    require_exact_type(state, PathState, "state")
    rebuilt = PathState(state.time, state.phase, state.target_lift,
                        state.error)
    if not rebuilt.has_target:
        if path.model == "width_only" and \
                path._cached_schedule.state_at(rebuilt.time) == "interior":
            raise ValueError(
                f"a width-only state at t={rebuilt.time!r} is strictly inside "
                f"its tongue and carries no lifted target; the control's "
                f"contraction is -rate*cos(error), so a state with no error "
                f"would be stepped and observed through the full Adler branch "
                f"instead")
        return rebuilt

    if path.model != "width_only":
        raise ValueError(
            f"a {path.model!r} state carries the width-only control's "
            f"auxiliary target and error; only the control has that coordinate "
            f"pair, and the full Adler model measures its lock error against "
            f"the moving stable phase directly")

    residual = rebuilt.lift_residual
    bound = conservation_tolerance(rebuilt.phase, rebuilt.target_lift,
                                   rebuilt.error)
    if not abs(residual) <= bound:
        raise ValueError(
            f"this state's target lift and error reconstruct to a phase "
            f"{residual!r} away from the one it carries, against a bound of "
            f"{bound!r}; the authoritative phase and the pair that is supposed "
            f"to describe it have come apart")

    where = path._cached_schedule.state_at(rebuilt.time)
    if where == "outside":
        raise ValueError(
            f"a width-only state at t={rebuilt.time!r} is outside its tongue "
            f"and carries a lifted target; the auxiliary pair is discarded at "
            f"the exact exit and rebuilt at the next exact entry, so it cannot "
            f"exist here")
    principal = (path._interior_target(rebuilt.time) if where == "interior"
                 else float(boundary_phase(path.detuning)))
    congruence = abs(float(_wrap_phase(rebuilt.target_lift - principal)))
    if congruence > LIFT_RESOLUTION:
        raise ValueError(
            f"this state's target lift is {congruence!r} rad from being "
            f"congruent to the {where} target {principal!r} of the clock it "
            f"belongs to, against the declared {LIFT_RESOLUTION!r} rad "
            f"resolution; a lift of the wrong target still reconstructs the "
            f"phase, and still changes sin(error) and therefore the next step")
    return rebuilt


def validate_clock_path(path) -> ClockPath:
    """Return a freshly rebuilt, fully revalidated copy of ``path``.

    The companion of ``raw_runner.validate_raw_config``, at the other end of
    the run, closing the same class of hole for the same reason.

    **Nothing stale.**  ``ClockPath`` derives its tongue schedule once, in the
    constructor, and caches it — and that schedule is load-bearing in a way no
    other field is.  It owns where the Brownian tree splits a step, where the
    entry and exit handoffs happen, which endpoints are strictly eligible,
    which endpoint target is ``boundary_phase`` rather than ``stable_phase``,
    and therefore where dwell can start and where it must reset.

    A frozen dataclass is frozen against ordinary assignment and not against
    ``object.__setattr__``.  An independent review used exactly that to change
    the *detuning* of a validated path after its schedule had been cached::

        path = ClockPath(0.9, train, "full", 0.0)
        object.__setattr__(path, "detuning", 0.1)

    The runner accepted it, recorded ``detuning=0.1``, and ran the whole
    trajectory against the tongue window of the ``0.9`` clock —
    ``(1.5903, 2.4097)`` instead of ``(0.4097, 3.5903)``.  Between ``0.4097``
    and ``1.5903`` the real clock is strictly eligible and the accepted path
    treated it as outside, so the transition branch, the target, the
    contraction and the dwell geometry were all the wrong ones, and the result
    was an ordinary committed record.  An earlier fix-up had adopted this exact
    threat model for ``fixed_rate`` alone; ``detuning`` and the drive are at
    least as semantic, because the schedule is derived from them.

    **Nothing trusted, at any depth.**  The rebuild does not accept a
    previously validated frozen container on the strength of its type: it
    reconstructs by exact type and named field, one level at a time, so that a
    mutation cannot simply move down a level.  Every ``RaisedCosinePulse``
    inside a ``PulseTrain`` is rebuilt from its own ``peak``, ``duration`` and
    ``center``, which re-runs that class's constructor checks; the train is
    rebuilt from those pulses, which re-runs the disjointness and ordering
    rules; and the path is rebuilt from the detuning, the rebuilt drive, the
    model and the rate, which re-runs the model/rate identity and derives a
    *fresh* schedule.

    **Then the two schedules must agree exactly.**  The freshly derived
    schedule is compared with the cached one by type and by value.  That is
    what catches a stale cache directly, including a schedule planted
    wholesale, and it is what makes "the geometry this run uses is the geometry
    its fields describe" a fact rather than an assumption.

    The rebuilt object is what the caller should run; the original is not
    returned.  For an unmutated path the rebuild is exact — identical floats
    give an identical schedule, an identical manifest digest and therefore
    identical noise keys — so a clean run is bit-for-bit what it was.

    ``TypeError`` if this is not exactly a ``ClockPath``, or if its drive or a
    nested member is not exactly one of the declared types.  ``ValueError`` if
    it is one but no longer holds consistent values.
    """
    require_exact_type(path, ClockPath, "path")

    drive = path.train
    if type(drive) is PulseTrain:
        pulses = []
        for position, pulse in enumerate(drive.pulses):
            require_exact_type(pulse, RaisedCosinePulse,
                               f"path.train.pulses[{position}]")
            pulses.append(RaisedCosinePulse(pulse.peak, pulse.duration,
                                            pulse.center))
        rebuilt_drive = PulseTrain(tuple(pulses))
    elif type(drive) is StationaryDrive:
        rebuilt_drive = StationaryDrive(drive.peak)
    else:
        raise TypeError(
            f"path.train must be exactly one of "
            f"{[t.__name__ for t in DRIVE_TYPES]}, got "
            f"{type(drive).__name__}")

    rebuilt = ClockPath(path.detuning, rebuilt_drive, path.model,
                        path.fixed_rate)

    cached = path._cached_schedule
    fresh = rebuilt._cached_schedule
    if type(cached) is not type(fresh) or cached != fresh:
        raise ValueError(
            f"this path's cached tongue schedule is {cached!r}, but its own "
            f"fields derive {fresh!r}; the schedule owns the Brownian crossing "
            f"split, the entry and exit handoffs, the eligibility states and "
            f"the dwell geometry, so a path whose geometry disagrees with its "
            f"detuning and drive would run one clock's trajectory and record "
            f"another's")
    return rebuilt
