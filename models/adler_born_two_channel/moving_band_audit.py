"""
moving_band_audit.py — the diagnostic piecewise-linear Brownian-bridge replay
that measures how much uninterrupted dwell the endpoint scheme overcounts on a
*moving* lock band.

What this is for
----------------
The primary event process evaluates the lock band at integration endpoints.  A
Brownian path can leave the band and come back between two endpoints, and the
endpoint scheme calls that segment continuously inside.  The stationary
killed-diffusion solver in ``killed_diffusion.py`` measures that error exactly,
in continuous time, for a band whose edges do not move.  During a pulse the
edges *do* move — the stable target sweeps as the coupling rises and falls, and
it sweeps fastest near tongue entry and exit, which is exactly where the error
is largest.

This module is the diagnostic for that case.  On every wholly-eligible
elementary interval it freezes the drift, represents both admissible-band edges
as straight lines between their endpoint values, and computes the conditional
probability that the Brownian bridge between the two recorded endpoints touched
either line.  A separately keyed auxiliary uniform decides whether a *hidden
exit* is inserted.

What this is NOT
----------------
It is **not** the continuum authority, and it is not allowed to become one.
Conditioned on its two endpoints, the one-sided bridge calculation is exact for
a frozen-drift Euler-Maruyama interpolation against a straight boundary.  It is
an *approximation* to the continuously varying nonlinear Adler drift and to the
genuinely moving lock interval, and the capped two-boundary combination is
conservative for that discretized bridge and not a rigorous bound on the
original continuous nonlinear process.  ``killed_diffusion.py`` remains the
oracle; this is a diagnostic beside it.

It is also **not** a correction to the primary ledger.  Three rules make that
mechanical rather than aspirational:

1. **Reset-only.**  The replay may add a dwell reset.  It may never remove one,
   never move a commitment earlier, and never turn an unresolved run into a
   committed one.  :class:`AuditedRun` refuses to be constructed if any of
   those happened, so a broken replay fails where it is recorded rather than
   producing a plausible row.
2. **Pathwise subset.**  Because resets only delay, an audited commitment is
   always a commitment the endpoint scheme also made, at a time no earlier.
   The audited committed set is therefore a subset of the endpoint committed
   set on the same paths, not a differently-sampled population.
3. **Disjoint randomness.**  The auxiliary uniforms come from
   :data:`AUDIT_KEY_SCHEMA`, which shares no prefix with
   ``stochastic.KEY_SCHEMA``.  The audit cannot read, consume or perturb a
   single value of the physical Brownian tree, so running it does not change
   the trajectory it is auditing.

The bridge calculation
----------------------
With the package's noise convention — a phase kick over duration ``h`` has
variance ``2 * D_phase * h`` — the conditional probability that a Brownian
bridge from ``phase_start`` to ``phase_end`` touched a straight boundary is

::

    p_upper = exp(-(upper_start - phase_start) * (upper_end - phase_end)
                  / (D_phase * duration))

    p_lower = exp(-(phase_start - lower_start) * (phase_end - lower_end)
                  / (D_phase * duration))

and the two are combined as ``min(p_lower + p_upper, 1)``.  Subtracting a
straight moving boundary from the conditioned path is what lets the fixed-boundary
one-sided result apply unchanged to a linearly moving one; that is the entire
reason the boundaries are represented linearly and the reason a curved
representation is a defect rather than a refinement.

The union cap is deliberately conservative: when both crossings are individually
likely the sum overstates the chance of touching either.  Overstating produces
*more* resets, which is the safe direction for a diagnostic whose only permitted
effect is to remove dwell.

Geometry, and when to give up
-----------------------------
The admissible lock interval is the fixed proximity band intersected with the
strictly contracting side of the flow.  Equality at either the proximity or the
contraction boundary is **outside**, matching ``commitment.LockCriterion``
exactly.  The phase and both edges are given continuous local representatives on
the circle.  Where that representation is ambiguous, where the admissible
interval disappears, or where its topology changes inside the step, the replay
**resets conservatively** instead of choosing whichever branch keeps the clock
locked.  Choosing the favourable branch is a named mutation in ``verify.py`` and
it is required to fail.

Non-claims
----------
See :data:`AUDIT_NON_CLAIMS`.  In particular this module computes no exponent,
selects no timestep, and cannot be used to rescue one.
"""
from __future__ import annotations

import hashlib
import math
from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np

from .commitment import DwellState, validate_lock_criterion, _advance_trusted
from .dynamics import nearest_lift, validate_clock_path
from .model import TWO_PI, boundary_phase, stable_phase
from .raw_experiments import elementary_segments
from .stochastic import FinestMesh, PhaseNoiseStream, TreeLeaf
from .validation import (require_exact_type, require_finite,
                         require_nonnegative, require_positive_int)

__all__ = [
    "AUDIT_KEY_SCHEMA",
    "AUDIT_NON_CLAIMS",
    "AUDIT_OUTCOMES",
    "GEOMETRY_STATES",
    "REDUCED_MATRIX",
    "AuditCell",
    "AuditUniformStream",
    "AuditedRun",
    "BridgeStep",
    "admissible_band",
    "audit_outcome",
    "bridge_crossing_probability",
    "bridge_step",
    "local_representative",
    "paired_leaves",
    "replay_pulse",
]

#: The root of every auxiliary audit key.  **Disjoint from
#: ``stochastic.KEY_SCHEMA`` by construction**, and the disjointness is the
#: point: the physical tree is ``dk-phase-noise/v3/physical`` and this is
#: ``dk-band-audit/v1/audit``, so no audit key can ever be a physical key and
#: no physical key can ever be an audit key.  ``verify.py`` asserts that no
#: string this module produces has the physical prefix, and mutates the schema
#: to the physical one to show the check bites.
#:
#: **v2.**  The key gained the audited interval's exact endpoints.  v1 addressed
#: a uniform by ``(finest step, piece ordinal)``, which is not an address of an
#: interval: coarse refinement levels are built on the same finest mesh, so one
#: v1 key named bridges of three different durations at three different strides.
#: A change to the key grammar is a change of dataset, so the version moves and
#: every auxiliary uniform changes with it.  No audit result has been recorded
#: under v1; nothing needs migrating.
AUDIT_KEY_SCHEMA = "dk-band-audit/v2/audit"

#: What a reader of an audit number is not entitled to conclude.
AUDIT_NON_CLAIMS = (
    "this is a frozen-drift discretized bridge, not the continuous nonlinear "
    "Adler process; the stationary killed-diffusion solver is the continuum "
    "oracle and this is a diagnostic beside it",
    "the capped two-boundary union is conservative for the discretized "
    "bridge only, and is not a rigorous bound on the original process",
    "the audit may only remove dwell; it is not a correction to the primary "
    "ledger and cannot replace, repair or reweight it",
    "an audit difference is a numerical measurement, not evidence about any "
    "mechanism, exponent or hazard law",
    "no timestep and no fit window may be selected by looking at these "
    "numbers; failure means refine or report a numerical no-result",
)

#: What the geometry of one interval turned out to be.  Only ``usable``
#: reaches the bridge expression; the other three reset conservatively.
GEOMETRY_STATES = ("usable", "ambiguous", "disappearing", "topology")

#: What the replay decided about one wholly-eligible elementary interval.
#:
#: ``primary_reset``
#:     An endpoint was on or outside the admissible interval.  The primary
#:     endpoint machine already resets here, so the audit adds nothing.
#: ``conservative_reset``
#:     The geometry was not usable, so dwell is cleared rather than a
#:     favourable representation chosen.
#: ``bridge_reset``
#:     Both endpoints were strictly inside and the auxiliary uniform fell below
#:     the capped crossing probability: a hidden exit is inserted.
#: ``no_reset``
#:     Both endpoints strictly inside and the bridge did not fire.
AUDIT_OUTCOMES = ("primary_reset", "conservative_reset", "bridge_reset",
                  "no_reset")

#: Longest permitted audit dataset namespace, and the characters it may hold.
#: A namespace is serialized into an ASCII key, so it has to round-trip
#: exactly; the same rule ``stochastic.py`` applies to its own namespaces.
_NAMESPACE_CHARS = frozenset(
    "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-")
_MAX_NAMESPACE = 64

#: The strict half-turn beyond which an Adler phase stops contracting.  The
#: full model's local contraction is ``-K cos(phi)``, which is negative exactly
#: on ``(-pi/2, pi/2)`` modulo a turn; the width-only control's is
#: ``-rate cos(error)``, negative on the same open interval in its own
#: coordinate.  One constant, both models, because it is one rule.
CONTRACTION_HALF_TURN = math.pi / 2.0


# ---------------------------------------------------------------------------
# Private validation helpers
# ---------------------------------------------------------------------------
def _require_token(value, name: str) -> str:
    """Return ``value`` as a key-safe namespace token, or raise naming it."""
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, str):
        raise TypeError(f"{name} must be a string, got "
                        f"{type(value).__name__} {value!r}")
    if not value:
        raise ValueError(f"{name} must not be empty")
    if len(value) > _MAX_NAMESPACE:
        raise ValueError(f"{name} must be at most {_MAX_NAMESPACE} "
                         f"characters, got {len(value)}")
    bad = sorted(set(value) - _NAMESPACE_CHARS)
    if bad:
        raise ValueError(f"{name} may only contain letters, digits, "
                         f"underscore and hyphen; {name}={value!r} contains "
                         f"{bad}")
    return value


def _require_label(value, name: str) -> str:
    """Return ``value`` as a non-empty string, or raise naming it."""
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, str):
        raise TypeError(f"{name} must be a string, got "
                        f"{type(value).__name__} {value!r}")
    if not value:
        raise ValueError(f"{name} must not be empty")
    return value


def _require_positive(value, name: str) -> float:
    """Return ``value`` as a strictly positive finite float."""
    number = require_finite(value, name)
    if not number > 0.0:
        raise ValueError(f"{name} must be strictly positive, got {number!r}")
    return number


def _require_time_or_none(value, name: str) -> Optional[float]:
    """Return ``value`` as a finite float, or ``None`` unchanged."""
    return None if value is None else require_finite(value, name)


# ---------------------------------------------------------------------------
# 1. The auxiliary uniform stream
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class AuditUniformStream:
    """Keyed uniforms for the bridge decision, in a namespace of their own.

    ``namespace``
        The audit dataset.  Separate from the physical noise namespace even
        when the two happen to be spelled the same way, because the schema
        prefixes differ; spelling them the same is nonetheless discouraged and
        the verifier reports the two prefixes side by side.
    ``mesh``
        The finest uniform grid the replay walks.  Carried into the key so that
        an audit uniform belongs to one mesh, exactly as a physical kick does.
    ``replicate``
        Which auxiliary replication this is.  The audit adds Monte Carlo
        variation of its own on top of the physical realization, and the plan
        requires that variation to be resolved rather than absorbed into the
        result.  Replicating means re-drawing *only* these uniforms: the phase
        path, the physical tree and the primary ledger are identical across
        replicates by construction.

    The derivation is written here rather than borrowed from
    ``stochastic.py``.  That is deliberate duplication: an audit stream that
    called the physical stream's key builder could be pointed at a physical key
    by a one-character edit, and the whole contract is that it cannot reach one.
    """

    namespace: str
    mesh: FinestMesh
    replicate: int = 0

    def __post_init__(self):
        namespace = _require_token(self.namespace, "namespace")
        mesh = require_exact_type(self.mesh, FinestMesh, "mesh")
        replicate = require_positive_int(self.replicate, "replicate",
                                         minimum=0)
        object.__setattr__(self, "namespace", namespace)
        object.__setattr__(self, "mesh", mesh)
        object.__setattr__(self, "replicate", replicate)

    @property
    def stream_prefix(self) -> str:
        """Everything a key holds before the trial: schema, dataset, mesh."""
        return (f"{AUDIT_KEY_SCHEMA}/{self.namespace}/"
                f"mesh={self.mesh.identity}/replicate={self.replicate}")

    def key(self, trial, clock, index, piece=0, span=None) -> str:
        """The canonical address of one auxiliary uniform.

        ``index`` and ``piece`` are the plan's declared address: the finest
        step, and the ordinal of the elementary interval inside it.  The
        ordinal is a strict refinement of the plan's ``audit / dataset / trial
        / clock / finest-step``, and it is necessary rather than decorative — a
        finest step can contain more than one wholly-eligible elementary
        interval, and a key without the ordinal would hand both the same
        uniform.

        ``span`` is ``(start_time, end_time)``, the exact endpoints of the
        physical interval the uniform decides, and it is what makes the address
        **canonical**.  Equal keys imply equal intervals because the endpoints
        are *in* the key, encoded by ``float.hex`` exactly as
        ``stochastic.Crossing.identity`` encodes a crossing time: two floats
        that differ in the last bit are two different addresses, and no
        rounding of a recomputed endpoint can alias one onto another.

        Why it had to be added: an independent review instrumented the declared
        pulsed path and found 73 keys shared between stride 1 and stride 2, and
        36 between stride 1 and stride 4, **every one of which described a
        different physical interval**.  The ordinal is recomputed from the
        current segmentation, and the coarse levels are built on the same
        finest mesh, so ``step=260/piece=0`` named a 0.008-long bridge at one
        refinement level and a 0.032-long bridge at another.  Reusing a uniform
        across refinement levels could have been a deliberate common-random-number
        coupling, but it was not declared as one and the key did not name what
        it addressed.

        Different refinement levels now draw independent auxiliary uniforms,
        which is the right default: they are different bridges over different
        durations, the physical realization is what is paired between them, and
        the audit's own Monte Carlo variation is resolved by replication rather
        than by coupling.

        ``span`` may be omitted only where there is no interval — the
        distributional probes in ``verify.py`` — and that case is written into
        the key as ``span=none`` rather than left ambiguous.
        """
        trial = require_positive_int(trial, "trial", minimum=0)
        clock = require_positive_int(clock, "clock", minimum=0)
        index = require_positive_int(index, "index", minimum=0)
        piece = require_positive_int(piece, "piece", minimum=0)
        if span is None:
            extent = "none"
        else:
            try:
                start, end = span
            except (TypeError, ValueError):
                raise TypeError(
                    f"span must be a (start_time, end_time) pair or None, got "
                    f"{span!r}") from None
            start = require_finite(start, "span[0]")
            end = require_finite(end, "span[1]")
            if not end > start:
                raise ValueError(
                    f"span must be a strictly positive interval, got "
                    f"({start!r}, {end!r}); a zero-duration handoff carries no "
                    f"bridge and has no auxiliary uniform to draw")
            extent = f"{start.hex()}:{end.hex()}"
        return (f"{self.stream_prefix}/trial={trial}/clock={clock}/"
                f"step={index}/piece={piece}/span={extent}/bridge")

    def uniform(self, trial, clock, index, piece=0, span=None) -> float:
        """The uniform on ``[0, 1)`` that :meth:`key` names.

        SHA-256 of the key, the leading eight bytes read as an unsigned 64-bit
        integer, divided by ``2**64``.  A pure function of the key: no
        generator state, no ordering dependence, and no way for batching, early
        stopping, window partition, pulse ordering or replicate ordering to
        change which value an interval gets.
        """
        digest = hashlib.sha256(
            self.key(trial, clock, index, piece, span).encode("ascii")).digest()
        return int.from_bytes(digest[:8], "big") / 18446744073709551616.0

    def uniforms(self, trial, clock, start_index, count) -> np.ndarray:
        """``count`` consecutive spanless first-piece uniforms, as an array.

        Convenience for the checks that need the stream's distributional
        properties; the replay itself asks for uniforms one interval at a time,
        with that interval's span, because which intervals are eligible is not
        known until the walk reaches them.
        """
        start_index = require_positive_int(start_index, "start_index",
                                           minimum=0)
        count = require_positive_int(count, "count", minimum=1)
        return np.array([self.uniform(trial, clock, start_index + offset)
                         for offset in range(count)], dtype=float)


# ---------------------------------------------------------------------------
# 2. Admissible geometry
# ---------------------------------------------------------------------------
def admissible_band(coupling, detuning, tolerance) -> Tuple[float, float]:
    """The principal admissible lock interval of the full Adler model.

    Returns ``(lower, upper)`` as an **open** interval of principal phase: the
    fixed proximity band around the stable target, intersected with the
    strictly contracting side of the flow.

    ``coupling > |detuning|`` is required — strict interiority, the same
    condition under which ``model.stable_phase`` is finite.  At equality the
    attracting and repelling fixed points have merged, the contraction is
    exactly zero, and there is no admissible interval at all; that is an
    eligibility boundary rather than a band and the replay never bridges
    across one.

    The intersection is a single arc, which is worth stating because on a
    circle it need not have been.  The proximity arc is centred on
    ``phi_star``, which lies strictly inside ``(-pi/2, pi/2)``; the arc the
    proximity band *excludes* is centred half a turn away from ``phi_star``,
    which is never inside ``(-pi/2, pi/2)``.  So the excluded arc cannot fall
    strictly inside the contracting arc and cut it in two, and the intersection
    of the two arcs is the plain interval computed here.  A degenerate or
    inverted result is returned as-is and classified ``disappearing`` by
    :func:`bridge_step` rather than silently repaired.
    """
    coupling = require_nonnegative(coupling, "coupling")
    detuning = require_finite(detuning, "detuning")
    tolerance = _require_positive(tolerance, "tolerance")
    if not tolerance < math.pi:
        raise ValueError(
            f"tolerance must be less than half a turn, got {tolerance!r}; a "
            f"circular distance never exceeds pi, so a wider tolerance admits "
            f"every phase and is not a band.  This is "
            f"commitment.LockCriterion's domain, restated where the geometry "
            f"is built")
    if not coupling > abs(detuning):
        raise ValueError(
            f"an admissible band exists only strictly inside the tongue: got "
            f"coupling={coupling!r} and |detuning|={abs(detuning)!r}.  At "
            f"equality stable_phase is NaN by design, the contraction is "
            f"exactly zero, and boundary_phase supplies the finite limiting "
            f"target for the coordinate handoff -- not a band")
    target = float(stable_phase(coupling, detuning))
    return (max(target - tolerance, -CONTRACTION_HALF_TURN),
            min(target + tolerance, CONTRACTION_HALF_TURN))


def local_representative(phase, lower, upper) -> Tuple[float, float, int, bool]:
    """Lift ``(lower, upper)`` onto the branch nearest ``phase``.

    Returns ``(lifted_lower, lifted_upper, winding, ambiguous)``.

    The band is an interval of principal phase; the phase is authoritative and
    unwrapped and may have accumulated any amount of winding.  The continuous
    local representative is the copy of the band shifted by whole turns so that
    its centre is nearest the phase.

    ``ambiguous`` is true when the phase sits exactly half a turn from the
    band's centre, so two whole-turn shifts are equally near.  There is no
    right answer there, and :func:`bridge_step` resets rather than picking one:
    a tie broken by convention would be a choice of geometry made by roundoff.
    """
    phase = require_finite(phase, "phase")
    lower = require_finite(lower, "lower")
    upper = require_finite(upper, "upper")
    centre = 0.5 * (lower + upper)
    offset = (phase - centre) / TWO_PI
    winding = math.floor(offset + 0.5)
    residual = phase - (centre + winding * TWO_PI)
    ambiguous = abs(abs(residual) - math.pi) == 0.0
    shift = winding * TWO_PI
    return (lower + shift, upper + shift, int(winding), ambiguous)


@dataclass(frozen=True)
class BridgeStep:
    """One wholly-eligible elementary interval, already in one local frame.

    Every field is in the same continuous coordinate: the authoritative phase
    at both endpoints and the two admissible edges at both endpoints, all
    lifted onto the branch the phase is on.  Building it is
    :func:`bridge_step`'s job; this type is what the bridge expression consumes,
    and it exists separately so that the expression can be exercised on declared
    numbers rather than only on a trajectory.

    ``geometry``
        One of :data:`GEOMETRY_STATES`.  Only ``usable`` reaches the bridge
        expression.
    """

    duration: float
    phase_start: float
    phase_end: float
    lower_start: float
    lower_end: float
    upper_start: float
    upper_end: float
    geometry: str = "usable"

    def __post_init__(self):
        duration = _require_positive(self.duration, "duration")
        if self.geometry not in GEOMETRY_STATES:
            raise ValueError(f"geometry must be one of "
                             f"{list(GEOMETRY_STATES)}, got {self.geometry!r}")
        for name in ("phase_start", "phase_end", "lower_start", "lower_end",
                     "upper_start", "upper_end"):
            object.__setattr__(self, name,
                               require_finite(getattr(self, name), name))
        object.__setattr__(self, "duration", duration)

    @property
    def clearances(self) -> Tuple[float, float, float, float]:
        """``(lower_start, lower_end, upper_start, upper_end)`` clearances.

        Each is the signed distance from the phase to that edge at that
        endpoint, oriented so that a *positive* value means strictly inside.
        Zero is outside, matching ``commitment.LockCriterion``'s strict band.
        """
        return (self.phase_start - self.lower_start,
                self.phase_end - self.lower_end,
                self.upper_start - self.phase_start,
                self.upper_end - self.phase_end)

    @property
    def strictly_inside(self) -> bool:
        """Are both endpoints strictly inside the admissible interval?"""
        return all(gap > 0.0 for gap in self.clearances)


def bridge_step(duration, phase_start, phase_end, band_start,
                band_end) -> BridgeStep:
    """Put one interval into a single local frame and classify its geometry.

    ``band_start`` and ``band_end`` are the principal admissible intervals at
    the interval's two endpoints, as :func:`admissible_band` returns them.  The
    edges are then represented by straight lines between their lifted endpoint
    values, which is the piecewise-linear representation the bridge expression
    requires and the only representation it is valid for.

    Three ways the geometry can fail, each producing a conservative reset:

    * ``disappearing`` — the admissible interval is empty or inverted at either
      endpoint.  There is no band to stay inside, so the dwell cannot continue.
    * ``ambiguous`` — a phase exactly half a turn from a band centre, so the
      local representative is a coin flip.
    * ``topology`` — the two endpoints resolve onto *different* whole-turn
      branches.  The clock changed which copy of the band it is near inside one
      step, and a straight line between edges on two different branches is not
      a representation of anything.
    """
    duration = _require_positive(duration, "duration")
    phase_start = require_finite(phase_start, "phase_start")
    phase_end = require_finite(phase_end, "phase_end")
    lower_start, upper_start = (require_finite(band_start[0],
                                              "band_start[0]"),
                                require_finite(band_start[1],
                                               "band_start[1]"))
    lower_end, upper_end = (require_finite(band_end[0], "band_end[0]"),
                            require_finite(band_end[1], "band_end[1]"))

    if not (upper_start > lower_start and upper_end > lower_end):
        return BridgeStep(duration, phase_start, phase_end, lower_start,
                          lower_end, upper_start, upper_end, "disappearing")

    low_s, up_s, wind_s, tie_s = local_representative(phase_start, lower_start,
                                                      upper_start)
    low_e, up_e, wind_e, tie_e = local_representative(phase_end, lower_end,
                                                      upper_end)
    if tie_s or tie_e:
        state = "ambiguous"
    elif wind_s != wind_e:
        state = "topology"
    else:
        state = "usable"
    return BridgeStep(duration, phase_start, phase_end, low_s, low_e, up_s,
                      up_e, state)


# ---------------------------------------------------------------------------
# 3. The bridge expression
# ---------------------------------------------------------------------------
def bridge_crossing_probability(step, diffusion) -> float:
    """The capped two-boundary bridge crossing probability of ``step``.

    Pure arithmetic on a :class:`BridgeStep` and a diffusion strength; it does
    not know what a tongue is and does not consult the geometry classification.
    The replay checks the geometry before calling this.

    Three special cases, in this order:

    1. **A zero or negative endpoint clearance returns exactly one**, before
       anything else is considered.  The phase is on or outside an edge at an
       endpoint, so it left the band with certainty rather than with a
       probability.  The plan calls this an immediate reset instead of an
       evaluation of the expression, and that is what this is.

       This clause is deliberately *first*.  It used to come second, after the
       zero-diffusion branch, and an independent review found the consequence:
       a step whose endpoint sat exactly on an edge with ``diffusion = 0``
       returned a crossing probability of **zero** — certain absorption
       reported as certain survival.  The replay itself never saw it, because
       ``audit_outcome`` checks interiority first, but this is an exported
       arithmetic API and a caller reaching for it directly got the wrong
       answer.  Ordering the two clauses the other way is a named mutation in
       ``verify.py``.
    2. **Zero diffusion then returns exactly zero.**  A path with no noise is
       the straight frozen-drift interpolation between its endpoints, and a
       straight line between two *strictly interior* points does not leave the
       interval.  The expression is not evaluated: ``D * duration`` is zero and
       the quotient does not exist.
    3. Otherwise the two one-sided probabilities are computed and their sum is
       capped at one.

    The cap is a union bound and it overstates the chance of touching either
    edge when both crossings are individually likely.  That is the conservative
    direction for a diagnostic whose only permitted effect is to remove dwell.
    """
    require_exact_type(step, BridgeStep, "step")
    # Rebuilt from its own fields: a step edited past its constructor -- a
    # negative duration, a NaN edge -- would otherwise be evaluated as found.
    step = BridgeStep(step.duration, step.phase_start, step.phase_end,
                      step.lower_start, step.lower_end, step.upper_start,
                      step.upper_end, step.geometry)
    diffusion = require_nonnegative(diffusion, "diffusion")

    # Certain absorption first, whatever the diffusion.  See the docstring for
    # the defect the reverse order produced.
    low_s, low_e, up_s, up_e = step.clearances
    if min(low_s, low_e, up_s, up_e) <= 0.0:
        return 1.0

    if diffusion == 0.0:
        return 0.0

    variance = diffusion * step.duration
    lower = math.exp(-(low_s * low_e) / variance)
    upper = math.exp(-(up_s * up_e) / variance)
    return min(lower + upper, 1.0)


def audit_outcome(step, diffusion, uniform) -> str:
    """What the replay does with one interval; one of :data:`AUDIT_OUTCOMES`.

    The decision order is the plan's, literally:

    1. unusable geometry resets conservatively, *before* any clearance or
       probability is looked at, so a favourable-looking bridge on an
       ambiguous frame is never even computed;
    2. an endpoint on or outside the admissible interval is left to the primary
       endpoint machine, which already resets there — the audit adds nothing
       and must not double-count a reset the ledger already has; and
    3. two strictly interior endpoints are decided by the capped probability
       against the auxiliary uniform, strictly: ``uniform < probability``
       fires, and equality does not.
    """
    require_exact_type(step, BridgeStep, "step")
    diffusion = require_nonnegative(diffusion, "diffusion")
    draw = require_finite(uniform, "uniform")
    if not 0.0 <= draw < 1.0:
        raise ValueError(f"uniform must lie in [0, 1), got {draw!r}")

    if step.geometry != "usable":
        return "conservative_reset"
    if not step.strictly_inside:
        return "primary_reset"
    return ("bridge_reset"
            if draw < bridge_crossing_probability(step, diffusion)
            else "no_reset")


# ---------------------------------------------------------------------------
# 4. The frozen reduced matrix
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class AuditCell:
    """One cell of the frozen reduced audit matrix.

    ``position``
        ``"central"``, ``"interior"`` or ``"near_edge"`` — where the clock sits
        in the tongue at peak coupling.  Named rather than derived, so a table
        of results reads without arithmetic.
    ``noise``
        ``"weak"``, ``"intermediate"`` or ``"strong"``.  Same reason.
    """

    label: str
    position: str
    noise: str
    detuning: float
    phase_diffusion: float

    def __post_init__(self):
        label = _require_label(self.label, "label")
        if self.position not in ("central", "interior", "near_edge"):
            raise ValueError(f"position must be central, interior or "
                             f"near_edge, got {self.position!r}")
        if self.noise not in ("weak", "intermediate", "strong"):
            raise ValueError(f"noise must be weak, intermediate or strong, "
                             f"got {self.noise!r}")
        detuning = require_finite(self.detuning, "detuning")
        diffusion = require_nonnegative(self.phase_diffusion,
                                        "phase_diffusion")
        if (self.position == "central") != (detuning == 0.0):
            raise ValueError(
                f"a central cell is the zero-mismatch clock and nothing else: "
                f"position={self.position!r} with detuning={detuning!r}")
        object.__setattr__(self, "label", label)
        object.__setattr__(self, "detuning", detuning)
        object.__setattr__(self, "phase_diffusion", diffusion)

    @property
    def sign(self) -> str:
        """``"zero"``, ``"positive"`` or ``"negative"`` mismatch sign."""
        if self.detuning == 0.0:
            return "zero"
        return "positive" if self.detuning > 0.0 else "negative"


def _matrix():
    """Build :data:`REDUCED_MATRIX` once, at import."""
    peak = 1.0
    places = (("central", 0.0), ("interior", 0.45 * peak),
              ("interior", -0.45 * peak), ("near_edge", 0.88 * peak),
              ("near_edge", -0.88 * peak))
    noises = (("weak", 0.002), ("intermediate", 0.02), ("strong", 0.08))
    cells = []
    for position, detuning in places:
        for noise, diffusion in noises:
            cells.append(AuditCell(f"S3b/{position}/{noise}/"
                                   f"detuning={detuning:+.3f}",
                                   position, noise, detuning, diffusion))
    return tuple(cells)


#: The frozen reduced matrix the audit runs on.  Fifteen cells: central,
#: interior and near-edge clocks; both mismatch signs at each off-centre
#: position; and weak, intermediate and strong noise at each.  Both pulse edges
#: are covered by construction rather than by a separate axis — every cell runs
#: the whole pulse, so every cell traverses the rising and falling segments of
#: its own eligibility window, and :class:`AuditedRun` counts the eligible
#: intervals on each side so that the coverage is measured rather than assumed.
#:
#: Frozen here, in source, before any audit number exists.  A matrix chosen
#: after seeing which cells behave well would be a selection rather than a
#: sample.
REDUCED_MATRIX = _matrix()


# ---------------------------------------------------------------------------
# 5. The replay
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class AuditedRun:
    """One replicate's audited outcome, recorded beside the primary one.

    **The reset-only contract is enforced in this constructor**, not asserted
    in prose somewhere else:

    * the audited run cannot have fewer resets than the primary;
    * a primary run that never committed cannot produce an audited commitment;
      and
    * an audited commitment cannot be earlier than the primary one.

    A replay that violated any of them would be a bug in the replay, and a
    record is the object a later reader quotes, so the violation fails here
    rather than becoming a plausible row.
    """

    label: str
    model: str
    detuning: float
    peak_coupling: float
    phase_diffusion: float
    timestep: float
    trial: int
    clock: int
    replicate: int
    primary_committed_at: Optional[float]
    audited_committed_at: Optional[float]
    primary_resets: int
    audited_resets: int
    shared_audited_resets: int
    bridge_resets: int
    conservative_resets: int
    shared_hidden_exits: int
    primary_reset_intervals: int
    evaluated_intervals: int
    rising_intervals: int
    falling_intervals: int

    def __post_init__(self):
        object.__setattr__(self, "label", _require_label(self.label, "label"))
        object.__setattr__(self, "model", _require_label(self.model, "model"))
        for name in ("detuning",):
            object.__setattr__(self, name,
                               require_finite(getattr(self, name), name))
        for name in ("peak_coupling", "phase_diffusion"):
            object.__setattr__(self, name,
                               require_nonnegative(getattr(self, name), name))
        object.__setattr__(self, "timestep",
                           _require_positive(self.timestep, "timestep"))
        for name in ("trial", "clock", "replicate", "primary_resets",
                     "audited_resets", "shared_audited_resets",
                     "bridge_resets", "conservative_resets",
                     "shared_hidden_exits", "primary_reset_intervals",
                     "evaluated_intervals", "rising_intervals",
                     "falling_intervals"):
            object.__setattr__(
                self, name,
                require_positive_int(getattr(self, name), name, minimum=0))
        first = _require_time_or_none(self.primary_committed_at,
                                      "primary_committed_at")
        second = _require_time_or_none(self.audited_committed_at,
                                       "audited_committed_at")
        object.__setattr__(self, "primary_committed_at", first)
        object.__setattr__(self, "audited_committed_at", second)

        if self.shared_audited_resets > self.audited_resets:
            raise ValueError(
                f"{self.shared_audited_resets} audited resets as of the "
                f"primary's stopping endpoint exceed the "
                f"{self.audited_resets} the whole audited walk holds; the "
                f"first is a prefix count of the second")
        if self.shared_hidden_exits > self.bridge_resets \
                + self.conservative_resets:
            raise ValueError(
                f"{self.shared_hidden_exits} hidden exits as of the primary's "
                f"stopping endpoint exceed the {self.bridge_resets} bridge "
                f"and {self.conservative_resets} conservative ones the whole "
                f"walk inserted")
        if self.shared_audited_resets < self.primary_resets:
            raise ValueError(
                f"over the exposure the two histories share, the audit "
                f"recorded {self.shared_audited_resets} resets against the "
                f"primary's {self.primary_resets}.  The replay may only *add* "
                f"resets; removing one would make the audited history a "
                f"correction to the ledger rather than a diagnostic beside it")
        if self.shared_audited_resets - self.primary_resets \
                > self.shared_hidden_exits:
            raise ValueError(
                f"over the shared exposure the audited history holds "
                f"{self.shared_audited_resets} resets against the primary's "
                f"{self.primary_resets}, an increase of "
                f"{self.shared_audited_resets - self.primary_resets}, but "
                f"only {self.shared_hidden_exits} hidden exits were inserted "
                f"in that span.  Every added reset must be one of those; the "
                f"reverse inequality is expected and allowed, since a hidden "
                f"exit on an interval with no running dwell interrupts "
                f"nothing")
        if first is None and second is not None:
            raise ValueError(
                f"the primary run never committed and the audited one "
                f"committed at t={second!r}.  Audited commitments are a "
                f"pathwise subset of endpoint commitments; the replay removes "
                f"dwell and can never create it")
        if first is not None and second is not None and second < first:
            raise ValueError(
                f"the audited commitment at t={second!r} precedes the primary "
                f"one at t={first!r}.  Inserting a hidden exit can only delay "
                f"a dwell, so an earlier audited commitment means the reset "
                f"was applied in the wrong direction")
        if self.evaluated_intervals < self.bridge_resets:
            raise ValueError(
                f"{self.bridge_resets} bridge resets from "
                f"{self.evaluated_intervals} evaluated intervals; a bridge "
                f"reset requires an interval whose geometry was usable and "
                f"whose endpoints were strictly inside")

    @property
    def added_resets(self) -> int:
        """How many resets the audit added over the **shared** exposure.

        The comparison is deliberately not ``audited_resets -
        primary_resets``.  The primary walk stops at its own commitment and the
        audited walk, whose dwell was delayed, keeps going; the resets it
        accumulates after that point are ordinary endpoint resets on an
        exposure the primary never had, not resets the audit created.  Counting
        them would make the audit look like it had added dwell interruptions
        that it did not, and would grow with the length of the run rather than
        with the numerical error being measured.
        """
        return self.shared_audited_resets - self.primary_resets

    @property
    def subset_holds(self) -> bool:
        """Is the audited commitment a delayed copy of the primary one?

        True when the audit committed nowhere the primary did not, and no
        earlier.  It is a property rather than an assertion because a reader
        of a table wants to see the column; the constructor has already refused
        the rows where it would be false.
        """
        if self.audited_committed_at is None:
            return True
        return (self.primary_committed_at is not None
                and self.audited_committed_at >= self.primary_committed_at)


def paired_leaves(stream, trial, clock, start_index, steps, stride,
                  crossings) -> Tuple[Tuple, Tuple]:
    """Return ``(leaves, segments)`` for one refinement level of one tree.

    ``stride`` finest steps per coarse step.  ``stride = 1`` is the finest
    level and returns exactly what ``PhaseNoiseStream.elementary_leaves`` and
    ``raw_experiments.elementary_segments`` return; a larger stride returns the
    *same realization* read at a coarser resolution.

    This is what makes the audit's refinement ladder **paired**.  A coarser
    level is not a second run on a coarser mesh — that would key a different
    tree and compare two realizations — it is the identical finest tree with
    its leaves summed into longer intervals, exactly as the plan requires:
    "obtain every coarser increment by exact summation of the corresponding
    fine increments".

    The coarse elementary mesh is the coarse uniform grid **union the same
    exact crossings**, built by the same half-open rule
    ``raw_experiments.elementary_segments`` uses, so no coarse interval ever
    straddles an eligibility change and no interval is ever half-eligible.
    Every fine segment therefore lies wholly inside exactly one coarse segment,
    and a coarse kick is the left-to-right sum of the fine kicks inside it.

    ``duration`` is computed as ``end - start`` from the two boundary floats,
    which is the same subtraction ``ClockPath._advance_trusted`` performs on
    the endpoint it is handed; a duration recomputed any other way is refused
    there rather than integrated.
    """
    require_exact_type(stream, PhaseNoiseStream, "stream")
    trial = require_positive_int(trial, "trial", minimum=0)
    clock = require_positive_int(clock, "clock", minimum=0)
    start_index = require_positive_int(start_index, "start_index", minimum=0)
    steps = require_positive_int(steps, "steps", minimum=1)
    stride = require_positive_int(stride, "stride", minimum=1)
    if steps % stride:
        raise ValueError(
            f"steps must be a whole number of coarse steps: steps={steps} is "
            f"not divisible by stride={stride}.  A partial trailing group "
            f"would be a coarse interval of the wrong duration, and both "
            f"keeping it and dropping it make two refinement levels disagree "
            f"for a reason that is not physics")

    fine_leaves = stream.elementary_leaves(trial, clock, start_index, steps,
                                           crossings)
    fine_segments = elementary_segments(stream.mesh, start_index, steps,
                                        crossings)
    if stride == 1:
        return fine_leaves, fine_segments

    moments = sorted({require_finite(getattr(crossing, "time", crossing),
                                     "crossing time")
                      for crossing in crossings})
    coarse_segments = []
    for group in range(steps // stride):
        first = stream.mesh.time_at(start_index + group * stride)
        second = stream.mesh.time_at(start_index + (group + 1) * stride)
        inner = [moment for moment in moments if first <= moment < second]
        bounds = [first] + inner + [second]
        coarse_segments.extend(zip(bounds, bounds[1:]))

    leaves = []
    position = 0
    for first, second in coarse_segments:
        total = 0.0
        started = False
        span_start = position
        while position < len(fine_segments) and fine_segments[position][1] <= second:
            # Seeded from the first kick rather than from 0.0, so the
            # accumulation is the same sequence of additions
            # ``PhaseNoiseStream.coarse_kicks`` performs and the two agree bit
            # for bit.
            total = (fine_leaves[position].kick if not started
                     else total + fine_leaves[position].kick)
            started = True
            position += 1
        if not started and second != first:
            raise ValueError(
                f"the coarse interval ({first!r}, {second!r}) contains no "
                f"fine elementary interval; the two meshes are not built from "
                f"one crossing schedule")
        leaves.append(TreeLeaf(first, second - first, total if started else 0.0,
                               "root",
                               stream.mesh.step_containing(first)
                               if second != first else
                               fine_leaves[max(span_start - 1, 0)].step_index))
    if position != len(fine_segments):
        raise ValueError(
            f"{len(fine_segments) - position} fine elementary intervals were "
            f"not consumed by any coarse interval; a coarse increment that is "
            f"not the exact sum of its fine leaves is a second realization "
            f"wearing the first one's name")
    return tuple(leaves), tuple(coarse_segments)


def _hidden_exit(state) -> DwellState:
    """Clear a running dwell, counting the reset; leave everything else alone.

    This is the *whole* effect the audit is permitted to have on a history.
    A live dwell is cleared and counted as a reset; a clock that was not inside
    is returned unchanged, since there is nothing to interrupt.  Band entries,
    latching flags, the commitment and the last-seen timestamp are untouched:
    the hidden exit happened *inside* the interval, so the endpoint that
    follows is still a fresh observation and may start a new dwell there.

    Built through ``DwellState``'s own constructor, so the reachability and
    conservation rules that class enforces run on the audited history exactly
    as they run on the primary one.

    A committed state is returned unchanged and is refused rather than cleared:
    the run ended at that commitment, so there is no later dwell for a hidden
    exit to interrupt, and clearing the dwell that produced the commitment
    would leave a ``committed_at`` with nothing behind it.
    """
    if state.committed:
        raise ValueError(
            f"a hidden exit was offered to a history that already committed "
            f"at t={state.committed_at!r}; the run ends at commitment, so "
            f"this is a walk that failed to stop rather than a reset")
    if state.inside_since is None:
        return state
    return DwellState(inside_since=None,
                      band_entries=state.band_entries,
                      dwell_resets=state.dwell_resets + 1,
                      committed_at=state.committed_at,
                      ever_eligible=state.ever_eligible,
                      ever_inside=state.ever_inside,
                      last_time=state.last_time)


def replay_pulse(path, criterion, stream, audit, trial, clock, initial_phase,
                 start_index, steps, label, replicates=1, stride=1):
    """Walk one clock once and audit it under ``replicates`` auxiliary draws.

    Returns a tuple of :class:`AuditedRun`, one per replicate, all describing
    the **same** physical trajectory: the phase path, the physical Brownian
    tree and the primary dwell history are computed once and shared, and only
    the auxiliary uniforms differ between replicates.  That is what makes the
    spread across the returned records a measurement of the audit's own Monte
    Carlo variation rather than of anything physical.

    The walk is the primary walk, interval for interval: the same elementary
    physical mesh, the same trusted kernels, the same end-of-step ordering.  It
    is re-derived here rather than imported from ``raw_experiments.run_one_clock``
    for one reason — that function stops at the primary commitment, and the
    audited history has to keep going past it, since a delayed dwell may commit
    later or not at all.  Everything else is the same call sequence, which is
    what makes "the audit sees the endpoint scheme's own history" a fact rather
    than a claim.

    ``audit`` supplies the auxiliary uniforms and is keyed in a namespace
    disjoint from ``stream``'s; nothing in this function reads one through the
    other.
    """
    path = validate_clock_path(path)
    criterion = validate_lock_criterion(criterion)
    require_exact_type(stream, PhaseNoiseStream, "stream")
    require_exact_type(audit, AuditUniformStream, "audit")
    trial = require_positive_int(trial, "trial", minimum=0)
    clock = require_positive_int(clock, "clock", minimum=0)
    phase = require_finite(initial_phase, "initial_phase")
    start_index = require_positive_int(start_index, "start_index", minimum=0)
    steps = require_positive_int(steps, "steps", minimum=1)
    label = _require_label(label, "label")
    replicates = require_positive_int(replicates, "replicates", minimum=1)
    stride = require_positive_int(stride, "stride", minimum=1)
    if audit.mesh != stream.mesh:
        raise ValueError(
            f"the audit stream is keyed to mesh {audit.mesh.identity!r} and "
            f"the physical stream to {stream.mesh.identity!r}; an audit "
            f"uniform belongs to one interval of one mesh, and two meshes "
            f"would put the decision on a different interval from the bridge")

    streams = tuple(AuditUniformStream(audit.namespace, audit.mesh,
                                       audit.replicate + offset)
                    for offset in range(replicates))

    state = path._initial_state_trusted(stream.mesh.time_at(start_index),
                                        phase)
    crossings = path._cached_schedule.crossings(clock, path.train.digest)
    leaves, segments = paired_leaves(stream, trial, clock, start_index, steps,
                                     stride, crossings)
    if len(leaves) != len(segments):
        raise ValueError(
            f"the noise tree produced {len(leaves)} elementary leaves and the "
            f"schedule {len(segments)} elementary intervals; the phase and "
            f"the noise are not being walked on one mesh")

    primary = DwellState()
    audited = [DwellState() for _ in streams]
    counts = [{"bridge": 0, "conservative": 0} for _ in streams]
    evaluated = 0
    primary_reset_intervals = 0
    rising = 0
    falling = 0
    previous_index = None
    piece = 0
    shared = None

    for leaf, (first, second) in zip(leaves, segments):
        index = stream.mesh.step_containing(first)
        piece = 0 if index != previous_index else piece + 1
        previous_index = index

        start_state = state
        state = path._advance_trusted(state, leaf, second)
        if second == first:
            continue                  # zero-duration handoff: not an endpoint

        eligible_interval = (path._cached_schedule.interval_state(first, second)
                             == "interior")
        if eligible_interval:
            # Which side of the envelope the interval sits on, taken from the
            # coupling itself rather than from a pulse centre, so a train with
            # more than one pulse classifies each of its windows correctly.
            if path.train.coupling(second) >= path.train.coupling(first):
                rising += 1
            else:
                falling += 1
            outcomes = _audit_interval(path, criterion, stream, streams,
                                       start_state, state, first, second,
                                       trial, clock, index, piece)
            if outcomes[0] == "primary_reset":
                primary_reset_intervals += 1
            elif outcomes[0] != "unbridgeable":
                evaluated += 1
            for position, verdict in enumerate(outcomes[1]):
                # A replicate that has already committed has ended its run, so
                # nothing later can interrupt a dwell it no longer has.  The
                # primary walk stops at its own commitment for the same reason;
                # the two stop independently because the audited one can only
                # be later.
                if audited[position].committed:
                    continue
                if verdict == "bridge_reset":
                    audited[position] = _hidden_exit(audited[position])
                    counts[position]["bridge"] += 1
                elif verdict == "conservative_reset":
                    audited[position] = _hidden_exit(audited[position])
                    counts[position]["conservative"] += 1

        seen = path._observe_trusted(state)
        if not primary.committed:
            primary = _advance_trusted(criterion, primary, seen.time,
                                       seen.eligible, seen.lock_error,
                                       seen.contraction)
        for position, current in enumerate(audited):
            if not current.committed:
                audited[position] = _advance_trusted(
                    criterion, current, seen.time, seen.eligible,
                    seen.lock_error, seen.contraction)

        # The shared-exposure snapshot.  The primary run ends here if it just
        # committed; everything the audited walk does after this endpoint
        # happens on exposure the primary never had, so the reset-only
        # comparison is frozen at this instant rather than at the end of the
        # longer walk.  See ``AuditedRun.added_resets``.
        if shared is None and primary.committed:
            shared = [(current.dwell_resets,
                       counts[position]["bridge"]
                       + counts[position]["conservative"])
                      for position, current in enumerate(audited)]

    if shared is None:
        # The primary never committed, so the two walks share the whole window.
        shared = [(current.dwell_resets,
                   counts[position]["bridge"]
                   + counts[position]["conservative"])
                  for position, current in enumerate(audited)]

    return tuple(
        AuditedRun(
            label=label,
            model=path.model,
            detuning=path.detuning,
            peak_coupling=path.train.peak_coupling,
            phase_diffusion=stream.phase_diffusion,
            timestep=stream.mesh.step * stride,
            trial=trial,
            clock=clock,
            replicate=streams[position].replicate,
            primary_committed_at=primary.committed_at,
            audited_committed_at=audited[position].committed_at,
            primary_resets=primary.dwell_resets,
            audited_resets=audited[position].dwell_resets,
            shared_audited_resets=shared[position][0],
            bridge_resets=counts[position]["bridge"],
            conservative_resets=counts[position]["conservative"],
            shared_hidden_exits=shared[position][1],
            primary_reset_intervals=primary_reset_intervals,
            evaluated_intervals=evaluated,
            rising_intervals=rising,
            falling_intervals=falling)
        for position in range(len(streams)))


def _audit_interval(path, criterion, stream, streams, start_state, end_state,
                    first, second, trial, clock, index, piece):
    """Classify one wholly-eligible interval and decide it for every replicate.

    Returns ``(shared_outcome, per_replicate_outcomes)``.  ``shared_outcome`` is
    the part of the decision that does not depend on the auxiliary draw —
    ``"primary_reset"``, ``"unbridgeable"`` for an unusable geometry, or
    ``"bridged"`` — and is what the caller counts intervals by.

    The width-only control is audited in its own coordinate.  Its admissible
    interval is the proximity band intersected with the contracting side of its
    *fixed-rate* restoring flow, which is ``(-b, b)`` about the lifted target
    with ``b = min(tolerance, pi/2)``; the target moves, so the edges move with
    it and the same linear representation applies.  The full model is audited in
    absolute phase against :func:`admissible_band`.
    """
    duration = second - first
    if duration <= 0.0:
        return ("unbridgeable", tuple("conservative_reset"
                                      for _ in streams))

    half = min(criterion.tolerance, CONTRACTION_HALF_TURN)
    if path.model == "width_only":
        # The control is audited in its own lifted-target coordinate: its
        # contraction is ``-rate cos(error)``, negative exactly on a half turn
        # about the target, so the admissible interval is the proximity band
        # intersected with that, which is ``target_lift +- min(tol, pi/2)``.
        # The target moves, so the edges move with it and the same linear
        # representation applies.
        #
        # On the interval that performs the entry handoff the incoming state
        # has no pair yet -- ``_advance_trusted`` establishes it inside the
        # step -- so the starting lift is re-derived here by the same rule the
        # handoff uses: the finite ``boundary_phase`` at eligibility equality,
        # lifted onto the branch nearest the authoritative phase.  Never
        # ``stable_phase``, which is NaN there by design.
        lift_start = (start_state.target_lift if start_state.has_target
                      else nearest_lift(float(boundary_phase(path.detuning)),
                                        start_state.phase))
        if not end_state.has_target:
            return ("unbridgeable",
                    tuple("conservative_reset" for _ in streams))
        band_start = (lift_start - half, lift_start + half)
        band_end = (end_state.target_lift - half,
                    end_state.target_lift + half)
    else:
        try:
            band_start = admissible_band(float(path.train.coupling(first)),
                                         path.detuning, criterion.tolerance)
            band_end = admissible_band(float(path.train.coupling(second)),
                                       path.detuning, criterion.tolerance)
        except ValueError:
            # The interval touches eligibility equality after all: no band
            # exists there, so there is nothing to bridge and the dwell is
            # cleared rather than bridged across a boundary.
            return ("unbridgeable",
                    tuple("conservative_reset" for _ in streams))

    step = bridge_step(duration, start_state.phase, end_state.phase,
                       band_start, band_end)
    if step.geometry != "usable":
        return ("unbridgeable",
                tuple("conservative_reset" for _ in streams))
    if not step.strictly_inside:
        return ("primary_reset", tuple("primary_reset" for _ in streams))

    probability = bridge_crossing_probability(step, stream.phase_diffusion)
    return ("bridged",
            tuple("bridge_reset"
                  if draw.uniform(trial, clock, index, piece,
                                  (first, second)) < probability
                  else "no_reset"
                  for draw in streams))
