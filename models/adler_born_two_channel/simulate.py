#!/usr/bin/env python3
"""
simulate.py — deterministic trajectories and the measurements taken from them.

Two experiments live here, both with the noise switched off:

* **stationary** — one clock, constant coupling.  Used to confirm that the
  attracting fixed point is where ``model.stable_phase`` says it is, that the
  approach to it is exponential, that the approach slows to a crawl at the
  tongue edge, and that a clock outside the tongue keeps slipping instead.

* **pulsed** — one clock, raised-cosine photon envelope.  Used to confirm that
  eligibility opens and closes at the times the instantaneous tongue boundary
  predicts, and to expose the finite-time failure mode: a clock near the
  boundary is admitted, but relaxes so slowly that the falling edge takes the
  fixed point away again before it has arrived.

Everything measured here is measured *from the trajectory* — a fitted slope, a
first-crossing time, a tracking error — never read off an analytic formula.
This module does not import ``analytic.py``, and ``verify.py`` enforces that by
parsing this file.  The comparisons between the two live in ``verify.py``, which
is the only place allowed to see both.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np

from .model import (STEPPERS, TWO_PI, RaisedCosinePulse, StationaryCoupling,
                    _circular_distance, _is_eligible, _stable_phase,
                    _wrap_phase)
from .validation import (require_boolean_array, require_finite,
                         require_finite_array, require_nonnegative,
                         require_nonnegative_array, require_positive_int,
                         require_real_allowing_nan)

__all__ = [
    "Trajectory",
    "integrate",
    "run_stationary",
    "run_pulsed",
    "propagate_to",
    "measured_relaxation_rate",
    "measured_eligibility_window",
    "measured_slip_frequency",
]


@dataclass
class Trajectory:
    """A deterministic run, plus everything derivable from it.

    ``phase`` is the *continuous* integrated phase error, deliberately not
    wrapped: unwrapping after the fact is lossy, and phase slips are counted by
    watching this quantity run away.  Use :attr:`wrapped` for the physical angle
    and :meth:`tracking_error` for distances, both of which go through the
    circular arithmetic in ``model.py``.
    """

    times: np.ndarray            # (T,)
    phase: np.ndarray            # (T, N) continuous
    detunings: np.ndarray        # (N,)
    coupling: np.ndarray         # (T,)
    dt: float
    method: str

    def __post_init__(self):
        self.times = require_finite_array(self.times, "times")
        self.phase = require_finite_array(self.phase, "phase")
        self.detunings = require_finite_array(self.detunings, "detunings")
        # Non-negative, not merely finite and correctly shaped.  A negative
        # coupling *history* of the right length passes every structural test
        # and then produces a completely plausible reading: `eligible` is False
        # at every sample, `stable` and `tracking_error` are NaN throughout —
        # i.e. "no clock was ever admitted to the tongue, and no lock was ever
        # defined".  That is exactly the physical result this package exists to
        # measure, arrived at from input the model does not define.
        self.coupling = require_nonnegative_array(self.coupling, "coupling")
        self.dt = require_finite(self.dt, "dt")
        if not self.dt > 0.0:
            raise ValueError(f"dt must be positive, got {self.dt!r}")
        if self.method not in STEPPERS:
            raise ValueError(f"unknown method {self.method!r}; expected one of "
                             f"{sorted(STEPPERS)}")
        if self.phase.shape != (self.times.size, self.detunings.size):
            raise ValueError(
                f"phase must have shape (n_times, n_clocks) = "
                f"({self.times.size}, {self.detunings.size}), got "
                f"{self.phase.shape}")
        if self.coupling.shape != self.times.shape:
            raise ValueError(f"coupling must match times, got "
                             f"{self.coupling.shape} and {self.times.shape}")

    # These properties use the unvalidated model kernels on purpose.  Their
    # inputs were validated when the trajectory was integrated, and `stable`
    # deliberately *produces* NaN for ineligible clocks — feeding that back
    # through a validating boundary would reject the module's own correct
    # output.  A NaN here means "no fixed point exists", never "bad input".
    @property
    def wrapped(self) -> np.ndarray:
        return _wrap_phase(self.phase)

    @property
    def eligible(self) -> np.ndarray:
        """(T, N) boolean: does a stable fixed point exist at this instant?"""
        return _is_eligible(self.coupling[:, None], self.detunings[None, :])

    @property
    def stable(self) -> np.ndarray:
        """(T, N) moving stable phase, NaN where the clock is ineligible."""
        return np.asarray(_stable_phase(self.coupling[:, None],
                                        self.detunings[None, :]), dtype=float)

    def tracking_error(self) -> np.ndarray:
        """(T, N) |phase - moving stable phase|, NaN where ineligible.

        NaN, not a large number: while a clock is ineligible there is no target
        to be far from, and filling in a distance would invent one.
        """
        target = self.stable
        err = np.abs(_circular_distance(self.phase, np.nan_to_num(target)))
        return np.where(np.isnan(target), np.nan, err)

    def index_at(self, t: float) -> int:
        """Index of the stored sample nearest ``t``.

        ``t`` is validated.  Unchecked, ``index_at(NaN)`` and ``index_at(+inf)``
        both return ``0``, because ``argmin`` over an all-NaN comparison array
        picks the first element — a wrong answer that looks like a legitimate
        "first sample" and would silently misplace whatever was being read off
        the trajectory.
        """
        return int(np.argmin(np.abs(self.times - require_finite(t, "t"))))


def integrate(detunings, source, phi0, t_start, t_end, dt,
              method: str = "rk4", store_every: int = 1) -> Trajectory:
    """Integrate ``dphi/dt = Delta - K(t) sin(phi)`` with the noise off.

    Vectorized over clocks, so the same code path will carry an ensemble when
    the stochastic stage arrives.  ``source`` is any object exposing
    ``coupling(t)`` — :class:`~.model.StationaryCoupling` or
    :class:`~.model.RaisedCosinePulse`.

    Every argument is validated here, at the public boundary.  ``dt`` in
    particular: a negative step used to be accepted silently, because
    ``max(1, round(span/dt))`` collapses it to a single forward step and the
    run then reports a positive effective step it never took.  A step size is
    not the kind of thing that should be quietly reinterpreted.
    """
    if method not in STEPPERS:
        raise ValueError(f"unknown method {method!r}; expected one of "
                         f"{sorted(STEPPERS)}")
    step = STEPPERS[method]

    D = np.atleast_1d(require_finite_array(detunings, "detunings"))
    phi = np.atleast_1d(require_finite_array(phi0, "phi0")) + np.zeros_like(D)
    if phi.shape != D.shape:
        raise ValueError("phi0 must be scalar or match the detuning shape")

    t_start = require_finite(t_start, "t_start")
    t_end = require_finite(t_end, "t_end")
    dt = require_finite(dt, "dt")
    if not dt > 0.0:
        raise ValueError(f"dt must be positive, got {dt!r}")

    # bool is a subclass of int, and `store_every=True` is a mistake, not a 1.
    store_every = require_positive_int(store_every, "store_every")

    if not hasattr(source, "_coupling_unchecked"):
        raise TypeError(f"source must be a coupling source exposing "
                        f"_coupling_unchecked(t), got "
                        f"{type(source).__name__}")
    coupling_of_t = source._coupling_unchecked

    span = t_end - t_start
    if not span > 0.0:
        raise ValueError(f"t_end must exceed t_start, got "
                         f"t_start={t_start!r}, t_end={t_end!r}")
    n_steps = max(1, int(round(span / dt)))
    dt_exact = span / n_steps

    times = [t_start]
    phases = [phi.copy()]
    t = t_start
    for i in range(n_steps):
        phi = step(phi, t, dt_exact, D, coupling_of_t)
        t = t_start + (i + 1) * dt_exact
        if (i + 1) % store_every == 0 or i == n_steps - 1:
            times.append(t)
            phases.append(phi.copy())

    times_arr = np.asarray(times, dtype=float)
    return Trajectory(
        times=times_arr,
        phase=np.asarray(phases, dtype=float),
        detunings=D,
        coupling=np.atleast_1d(np.asarray(coupling_of_t(times_arr),
                                          dtype=float)),
        dt=dt_exact,
        method=method,
    )


def propagate_to(detunings, source, phi0, t_start: float, t_target: float,
                 dt: float, method: str = "rk4") -> np.ndarray:
    """Phase of each clock at *exactly* ``t_target``, by direct integration.

    This is the primitive the exact-entry classification needs.  A trajectory
    stored on a fixed sample grid does not generally contain the moment a clock
    enters its tongue: the continuous entry time falls between samples, and with
    ``store_every`` set the gap can be several steps wide.  Asking "was this
    clock inside the lock band when it entered?" of the first *stored* eligible
    sample therefore answers a slightly different question — a clock can enter
    the band during the interval between the true crossing and that sample.

    Integrating to the target time directly removes the ambiguity.
    :func:`integrate` chooses ``n_steps = round(span/dt)`` and then uses
    ``dt_exact = span/n_steps``, so the final sample lands on ``t_target``
    exactly rather than near it, whatever ``dt`` was asked for.  The remaining
    error is the integrator's own, which ``verify.py`` measures by refinement
    and finds to be around 1e-11 rad — nine orders below the lock band, so the
    classification it feeds cannot be step-sensitive.

    Only the endpoint is retained.  The intermediate samples are what make a
    long integration expensive, and none of them are wanted here.
    """
    t_start = require_finite(t_start, "t_start")
    t_target = require_finite(t_target, "t_target")
    dt = require_finite(dt, "dt")
    if not dt > 0.0:
        raise ValueError(f"dt must be positive, got {dt!r}")
    n_steps = max(1, int(round((t_target - t_start) / dt)))
    traj = integrate(detunings, source, phi0, t_start, t_target, dt,
                     method=method, store_every=n_steps)
    return traj.phase[-1]


def run_stationary(detunings, coupling: float, phi0, t_end: float,
                   dt: float, method: str = "rk4",
                   store_every: int = 1) -> Trajectory:
    """Constant-coupling calibration run starting at ``t = 0``.

    ``coupling`` is handed to :class:`~.model.StationaryCoupling` unconverted,
    so that its own validation sees what the caller actually passed.  Coercing
    it here first would let ``True`` through as ``1.0``, and would strip the
    non-negativity check that class applies.
    """
    return integrate(detunings, StationaryCoupling(coupling), phi0,
                     0.0, require_finite(t_end, "t_end"), dt, method=method,
                     store_every=store_every)


def run_pulsed(detunings, pulse: RaisedCosinePulse, phi0, dt: float,
               margin: float = 0.0, method: str = "rk4",
               store_every: int = 1) -> Trajectory:
    """Photon-pulse run spanning the whole support, plus optional dead margin.

    The margin exists so the run demonstrably begins and ends with the clock
    ineligible, rather than starting at the instant of the crossing.
    """
    if not isinstance(pulse, RaisedCosinePulse):
        raise TypeError(f"pulse must be a RaisedCosinePulse, got "
                        f"{type(pulse).__name__}")
    margin = require_finite(margin, "margin")
    if margin < 0.0:
        raise ValueError(f"margin must be non-negative, got {margin!r}")
    t0, t1 = pulse.support()
    return integrate(detunings, pulse, phi0, t0 - margin, t1 + margin, dt,
                     method=method, store_every=store_every)


# ---------------------------------------------------------------------------
# Measurements taken from a trajectory
# ---------------------------------------------------------------------------
def measured_relaxation_rate(times: np.ndarray, deviation: np.ndarray,
                             lo: float, hi: float) -> Tuple[float, int]:
    """Fit ``deviation ~ exp(-rate * t)`` over the band ``lo <= dev <= hi``.

    The band matters in both directions.  Above ``hi`` the decay is not yet
    linear, so the quadratic term in the Adler expansion biases the slope; below
    ``lo`` the deviation is down among the integrator's own round-off and the
    logarithm turns into noise.  Returns ``(rate, n_points_used)``.
    """
    t = np.atleast_1d(require_finite_array(times, "times"))
    lo = require_finite(lo, "lo")
    hi = require_finite(hi, "hi")
    if not 0.0 < lo < hi:
        raise ValueError(f"require 0 < lo < hi, got lo={lo!r}, hi={hi!r}")
    # `deviation` is the one argument allowed to carry a non-finite entry, and
    # it is allowed exactly one kind: NaN, the documented sentinel for "no
    # target existed at this sample", which `Trajectory.tracking_error()` emits
    # wherever a clock was ineligible.  Infinities are rejected — an infinity
    # is a numerical blow-up, and quietly filtering it would hide the failure it
    # signals — as are booleans, strings and other wrong-kind arrays.
    dev = np.atleast_1d(require_real_allowing_nan(deviation, "deviation"))
    if dev.shape != t.shape:
        raise ValueError(f"times and deviation must have the same shape, got "
                         f"{t.shape} and {dev.shape}")
    band = np.isfinite(dev) & (dev >= lo) & (dev <= hi)
    n = int(np.count_nonzero(band))
    if n < 3:
        return float("nan"), n
    slope, _ = np.polyfit(t[band], np.log(dev[band]), 1)
    return float(-slope), n


def measured_eligibility_window(times: np.ndarray, eligible: np.ndarray
                                ) -> Optional[Tuple[float, float]]:
    """First and last eligible sample times, refined to the sub-step midpoint.

    The true crossing lies between the last ineligible sample and the first
    eligible one, so quoting the midpoint of that bracket halves the sampling
    error and makes the residual scale cleanly with the step size.  Returns
    ``None`` if the clock was never eligible.
    """
    t = np.atleast_1d(require_finite_array(times, "times"))
    # A genuine boolean array, not something coercible to one: an array of the
    # strings 'False'/'True' coerces to all-true, and [0, 2, 0] coerces to
    # [False, True, False].  Both would yield a plausible window from input the
    # caller never meant as a predicate.
    elig = require_boolean_array(eligible, "eligible").ravel()
    if elig.size != t.size:
        raise ValueError(f"times and eligible must have the same length, got "
                         f"{t.size} and {elig.size}")
    idx = np.flatnonzero(elig)
    if idx.size == 0:
        return None
    first, last = int(idx[0]), int(idx[-1])
    t_in = t[first] if first == 0 else 0.5 * (t[first - 1] + t[first])
    t_out = (t[last] if last == t.size - 1
             else 0.5 * (t[last] + t[last + 1]))
    return (float(t_in), float(t_out))


def measured_slip_frequency(times: np.ndarray, phase: np.ndarray) -> float:
    """Mean phase-slip frequency, measured from successive 2*pi crossings.

    For a clock outside the tongue the phase never settles: it runs, slowly near
    ``phi = pi/2`` and quickly elsewhere, and completes a full turn every slip
    period.  That staircase is exactly why a least-squares slope is the wrong
    estimator here — it is biased by whatever fraction of a period the run
    happens to end on.  Timing whole turns removes the bias: the crossing times
    are interpolated linearly, so the estimate is second-order accurate in the
    step and needs only a handful of slips.

    Raises ``ValueError`` if the phase is not monotone (i.e. the clock is not
    slipping) or if fewer than two full turns were completed.
    """
    t = np.atleast_1d(require_finite_array(times, "times"))
    p = np.atleast_1d(require_finite_array(phase, "phase")).ravel()
    if p.size != t.size:
        raise ValueError(f"times and phase must have the same length, got "
                         f"{t.size} and {p.size}")
    total = p[-1] - p[0]
    direction = np.sign(total)
    if direction == 0.0:
        raise ValueError("phase did not advance; this clock is not slipping")

    turns = direction * (p - p[0])                 # monotone rising if slipping
    if np.any(np.diff(turns) < 0.0):
        raise ValueError("phase is not monotone; this clock is not slipping")
    n_turns = int(np.floor(turns[-1] / TWO_PI))
    if n_turns < 2:
        raise ValueError(f"need at least two full slips, saw {n_turns}")

    def crossing_time(level: float) -> float:
        i = int(np.searchsorted(turns, level))
        lo, hi = turns[i - 1], turns[i]
        return float(t[i - 1] + (t[i] - t[i - 1]) * (level - lo) / (hi - lo))

    first = crossing_time(TWO_PI)
    last = crossing_time(n_turns * TWO_PI)
    return float(direction * TWO_PI * (n_turns - 1) / (last - first))
