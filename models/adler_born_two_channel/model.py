#!/usr/bin/env python3
"""
model.py — phases, Arnold tongues, photon-pulse envelopes, and deterministic
Adler dynamics for the two-channel selection test.

Scope of this module
--------------------
This is the *dynamics* layer.  It knows how a clock's phase error moves:

    dphi/dt = Delta - K(t) * sin(phi)

where ``Delta`` is the clock's fixed natural-frequency detuning and ``K(t)`` is
the instantaneous, pulse-scaled coupling to its channel reference.  That is the
Adler equation, written for the phase error rather than the phase itself.

Deliberate boundaries (these are load-bearing, not stylistic)
-------------------------------------------------------------
* This module receives a *coupling*.  It never receives an amplitude squared,
  an intensity, a probability weight, or a precomputed commitment hazard.
* This module never imports, mentions, or calls ``analytic.py``.  The analytic
  relaxation-rate prediction that the whole experiment is trying to test must
  not be able to leak into the thing generating the events.  ``verify.py``
  enforces that structurally, by parsing this file.
* The eligibility predicate ``K > |Delta|`` is implemented here *and*
  independently in ``analytic.py``.  That duplication is on purpose: the two
  copies are cross-checked against each other in ``verify.py``, which is a
  stronger guarantee than an import would give, and it keeps the analytic
  module standalone.

Public API versus private kernels
---------------------------------
Every public callable here validates its arguments and raises on anything that
is not a finite real number — including booleans, strings and NaN.  That is a
correctness requirement rather than defensive habit: a NaN coupling makes
``K > |Delta|`` false everywhere, so an invalid run *looks exactly like* the
physical result "no clock was ever admitted to the tongue".  An invalid input
must never be able to impersonate a physical conclusion.

Validation costs a scan of the input, which is unaffordable four times per
Runge-Kutta step, so each validated public function is a thin boundary over an
unvalidated private kernel (``_wrap_phase``, ``_is_eligible``, ``_adler_drift``,
``_raised_cosine``, ...).  The integrator and the steppers call the kernels;
everything reachable from outside the package calls the boundary.  The rule is
that a kernel may only be called with values some public boundary has already
checked.

What is *not* here yet
----------------------
Noise, the lock-band dwell counter, commitment, channel competition, and
polarization winner frequencies all belong to later stages.  The vocabulary is
kept in the docstrings so that the later stages can slot in without renaming
anything:

* **eligible** — a stable Adler fixed point exists right now (``K > |Delta|``).
  An ineligible clock may still drift; it simply has nothing to lock onto.
* **lock-failed** — eligible for a while, and never reached the lock band.
* **dwell-failed** — reached the lock band, but lost it (or lost eligibility on
  the falling edge) before the fixed commitment dwell had elapsed.
* **unresolved** — the pulse ended with no commitment, for any of those reasons.

Only *eligible* is simulated at this stage.  The other three are documented
here so they are not silently redefined later.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np

from .validation import (require_boolean_array, require_callable,
                         require_finite, require_finite_array,
                         require_nonnegative, require_nonnegative_array,
                         require_positive_int)

TWO_PI = 2.0 * np.pi

# The validators are imported for internal use but deliberately NOT re-exported:
# they are `validation.py`'s public API, and listing them here too would give
# the same contract two owners.
__all__ = [
    "TWO_PI",
    "boundary_phase",
    "wrap_phase",
    "circular_distance",
    "is_eligible",
    "stable_phase",
    "unstable_phase",
    "adler_drift",
    "raised_cosine_envelope",
    "StationaryCoupling",
    "RaisedCosinePulse",
    "step_euler",
    "step_rk4",
    "STEPPERS",
]


# ---------------------------------------------------------------------------
# 0. Input validation
# ---------------------------------------------------------------------------
# The checks themselves live in `validation.py`, which contains no physics, so
# that `analytic.py` can share them without importing anything from the
# dynamics.  See that module for why they are strict about NaN, bool and str.
#
# Every public callable below validates; every private `_`-prefixed kernel does
# not, because the Runge-Kutta inner loop calls the kernels four times a step
# with values a public boundary has already checked.


# ---------------------------------------------------------------------------
# 1. Phase algebra
# ---------------------------------------------------------------------------
def _wrap_phase(phi):
    """Unvalidated kernel behind :func:`wrap_phase`."""
    return np.mod(np.asarray(phi, dtype=float) + np.pi, TWO_PI) - np.pi


def _circular_distance(a, b):
    """Unvalidated kernel behind :func:`circular_distance`."""
    return _wrap_phase(np.asarray(a, dtype=float) - np.asarray(b, dtype=float))


def wrap_phase(phi):
    """Wrap an angle into the half-open interval [-pi, pi).

    The half-open convention matters.  ``+pi`` and ``-pi`` are the same angle,
    so exactly one of them has to be chosen as the representative; here it is
    ``-pi``.  Every comparison in this package goes through :func:`wrap_phase`
    or :func:`circular_distance`, so the choice is made once rather than being
    re-improvised at each call site.
    """
    return _wrap_phase(require_finite_array(phi, "phi"))[()]


def circular_distance(a, b):
    """Signed shortest angular separation ``a - b``, in [-pi, pi).

    This is the only correct way to ask "how far is this clock from its stable
    phase?"  A plain subtraction reports a separation of nearly 2*pi for two
    phases that sit either side of the wrap boundary and are in fact adjacent.
    """
    return _circular_distance(require_finite_array(a, "a"),
                              require_finite_array(b, "b"))[()]


# ---------------------------------------------------------------------------
# 2. Arnold tongue: eligibility and fixed points
# ---------------------------------------------------------------------------
def _is_eligible(coupling, detuning):
    """Unvalidated kernel behind :func:`is_eligible`."""
    K = np.asarray(coupling, dtype=float)
    D = np.asarray(detuning, dtype=float)
    return K > np.abs(D)


def _stable_phase(coupling, detuning):
    """Unvalidated kernel behind :func:`stable_phase`."""
    K = np.asarray(coupling, dtype=float)
    D = np.asarray(detuning, dtype=float)
    elig = _is_eligible(K, D)
    safe_K = np.where(elig, K, 1.0)
    ratio = np.clip(np.where(elig, D / safe_K, 0.0), -1.0, 1.0)
    return np.where(elig, np.arcsin(ratio), np.nan)


def is_eligible(coupling, detuning):
    """True where a stable Adler fixed point exists: ``K > |Delta|``.

    Strict inequality.  At ``K == |Delta|`` the stable and unstable fixed points
    have merged into a single half-stable point at ``phi = +-pi/2``, the local
    relaxation rate is exactly zero, and calling that clock "locked" would be
    an overstatement.  ``K = 0`` therefore admits nothing at all, including a
    clock with ``Delta = 0``.

    The tongue boundary lives in coupling-versus-detuning space, not in
    phase-angle space: a clock does not leave a stationary tongue merely because
    its phase rotates.  Leaving requires ``K`` or ``Delta`` to change.  In this
    model the falling edge of the photon pulse supplies that change.

    Both arguments must be finite, and the coupling must be non-negative.  A
    NaN coupling would otherwise make this predicate false everywhere and so
    imitate a physical no-admission result; a *negative* one would do the same
    thing while looking even more ordinary.  Detunings are signed and may be
    negative — that is the physical convention — but a coupling may not be.
    """
    return _is_eligible(require_nonnegative_array(coupling, "coupling"),
                        require_finite_array(detuning, "detuning"))[()]


def stable_phase(coupling, detuning):
    """The attracting Adler fixed point ``arcsin(Delta / K)``, else NaN.

    Of the two fixed points of ``Delta - K sin(phi)`` the attracting one is the
    branch with ``cos(phi) > 0``, i.e. ``arcsin(Delta/K)`` in [-pi/2, pi/2].
    Where no fixed point exists the answer is NaN rather than a number, so that
    an ineligible clock cannot be accidentally treated as having a target.

    That NaN is an *output* meaning "no fixed point exists here".  A NaN
    *input* is a different thing entirely and is rejected, so the two can never
    be confused — as is a negative coupling, which would silently move the
    attracting fixed point by pi rather than weaken it.
    """
    return _stable_phase(require_nonnegative_array(coupling, "coupling"),
                         require_finite_array(detuning, "detuning"))[()]


def unstable_phase(coupling, detuning):
    """The repelling fixed point, ``pi - arcsin(Delta/K)`` wrapped, else NaN.

    Not used by the dynamics; it exists so ``verify.py`` can confirm that the
    fixed point this module calls "stable" really is the attracting one and the
    other really is not.
    """
    phi_s = np.asarray(stable_phase(coupling, detuning), dtype=float)
    return np.where(np.isnan(phi_s), np.nan, _wrap_phase(np.pi - phi_s))[()]


def boundary_phase(detuning):
    """The limiting fixed phase at the tongue boundary, ``sign(Delta)*pi/2``.

    Exactly at ``K = |Delta|`` the stable and unstable fixed points have merged,
    :func:`stable_phase` is NaN by construction, and yet the *limit* of
    ``arcsin(Delta/K)`` as ``K`` falls to ``|Delta|`` is perfectly well defined.
    That limit is the only sensible target to measure a clock against at the
    instant it enters the tongue, which is what the exact-entry classification
    in ``verify.py`` needs.

    ``Delta = 0`` is the degenerate case: the tongue opens at ``K = 0`` where
    there is no pull at all, and the limit is taken as ``0.0``.
    """
    D = require_finite_array(detuning, "detuning")
    return (np.sign(D) * (0.5 * np.pi))[()]


# ---------------------------------------------------------------------------
# 3. Deterministic Adler drift
# ---------------------------------------------------------------------------
def _adler_drift(phi, coupling, detuning):
    """Unvalidated kernel behind :func:`adler_drift`; the RK4 hot path."""
    return (np.asarray(detuning, dtype=float)
            - np.asarray(coupling, dtype=float) * np.sin(phi))


def adler_drift(phi, coupling, detuning):
    """``dphi/dt = Delta - K sin(phi)``.

    Note what is absent: no ``K**2``, no intensity, no hazard, no noise.  The
    stochastic kick is added by a later stage, outside this function.
    """
    return _adler_drift(require_finite_array(phi, "phi"),
                        require_nonnegative_array(coupling, "coupling"),
                        require_finite_array(detuning, "detuning"))[()]


# ---------------------------------------------------------------------------
# 4. Photon-pulse envelope and coupling sources
# ---------------------------------------------------------------------------
def _raised_cosine(t, center: float, duration: float):
    """Unvalidated kernel behind :func:`raised_cosine_envelope`.

    Called four times per Runge-Kutta step through an already-validated
    :class:`RaisedCosinePulse`, so it must not re-check constants that cannot
    have changed since the constructor checked them.
    """
    s = (np.asarray(t, dtype=float) - center) / duration
    inside = np.abs(s) < 0.5
    return np.where(inside, 0.5 * (1.0 + np.cos(TWO_PI * s)), 0.0)


def raised_cosine_envelope(t, center, duration):
    """Normalized raised-cosine (Hann) envelope on a compact support.

        E(t) = 0.5 * (1 + cos(2*pi*(t - center)/duration))   on the support
        E(t) = 0                                             outside

    It rises from exactly zero, reaches exactly one at ``t = center``, and
    returns to exactly zero at ``center +- duration/2``.  It is symmetric about
    the centre and never negative.

    This is the *shared* envelope: in the two-channel stage both channels are
    multiplied by this same E(t) and differ only in their peak coupling, so the
    polarization angle cannot enter through the pulse shape.

    The exact interior shape is load-bearing, not decorative: several other
    envelopes share these endpoints, this peak, this symmetry and even this
    area while differing through the interior, and they give different tongue
    entry and exit times.  ``verify.py`` therefore compares this function
    against an independently written copy of the formula above rather than
    against itself.
    """
    center = require_finite(center, "envelope centre")
    duration = require_finite(duration, "envelope duration")
    if not duration > 0.0:
        raise ValueError(f"envelope duration must be positive, got {duration!r}")
    return _raised_cosine(require_finite_array(t, "envelope time"),
                          center, duration)[()]


@dataclass(frozen=True)
class StationaryCoupling:
    """Constant coupling ``K(t) = peak`` for all time.

    The Experiment-0 calibration case, and the control that a very long pulse
    must reproduce over its central plateau.
    """

    peak: float

    def __post_init__(self):
        peak = require_finite(self.peak, "peak coupling")
        if peak < 0.0:
            raise ValueError(f"peak coupling must be non-negative, got {peak!r}")
        # Write the normalized float back, so a frozen instance can never hold
        # something that merely *converted* to a float once in the constructor.
        object.__setattr__(self, "peak", peak)

    def _coupling_unchecked(self, t):
        return np.broadcast_to(np.float64(self.peak),
                               np.shape(np.asarray(t, dtype=float)))

    def coupling(self, t):
        """The coupling at time(s) ``t``; constant, but ``t`` is still checked."""
        return self._coupling_unchecked(
            require_finite_array(t, "time"))[()]

    def eligibility_window(self, detuning) -> Optional[Tuple[float, float]]:
        """``(-inf, +inf)`` if the clock is inside the tongue, else ``None``."""
        detuning = require_finite(detuning, "detuning")
        if not bool(_is_eligible(self.peak, detuning)):
            return None
        return (-np.inf, np.inf)

    def support(self) -> Tuple[float, float]:
        return (-np.inf, np.inf)


@dataclass(frozen=True)
class RaisedCosinePulse:
    """A passing photon wave packet: ``K(t) = peak * E(t)``.

    ``peak`` is a *coupling*, amplitude-linear in the field.  ``duration`` is
    the full compact support, so the pulse is non-zero on
    ``(center - duration/2, center + duration/2)``.

    Pulse area is **not** renormalized when the duration changes.  Fixed-peak
    and fixed-area duration sweeps answer different physical questions and must
    be run and labelled separately; silently holding the area fixed would
    conflate them.
    """

    peak: float
    duration: float
    center: float = 0.0

    def __post_init__(self):
        peak = require_finite(self.peak, "peak coupling")
        duration = require_finite(self.duration, "pulse duration")
        center = require_finite(self.center, "pulse centre")
        if peak < 0.0:
            raise ValueError(f"peak coupling must be non-negative, got {peak!r}")
        if not duration > 0.0:
            raise ValueError(f"pulse duration must be positive, got {duration!r}")
        # Store the normalized floats, not whatever converted to them.  Without
        # this a frozen instance built from strings keeps the strings and fails
        # later, somewhere else, with an exception that names the wrong thing.
        object.__setattr__(self, "peak", peak)
        object.__setattr__(self, "duration", duration)
        object.__setattr__(self, "center", center)

    def _envelope_unchecked(self, t):
        return _raised_cosine(t, self.center, self.duration)

    def _coupling_unchecked(self, t):
        return self.peak * _raised_cosine(t, self.center, self.duration)

    def envelope(self, t):
        """The normalized envelope at time(s) ``t``, in [0, 1]."""
        return self._envelope_unchecked(
            require_finite_array(t, "time"))[()]

    def coupling(self, t):
        """The coupling at time(s) ``t``: ``peak`` times the envelope."""
        return self._coupling_unchecked(
            require_finite_array(t, "time"))[()]

    def area(self) -> float:
        """Exact integral of ``K(t) dt`` = ``peak * duration / 2``."""
        return 0.5 * self.peak * self.duration

    def eligibility_window(self, detuning) -> Optional[Tuple[float, float]]:
        """Exact times at which ``K(t)`` crosses ``|Delta|``, or ``None``.

        Solving ``0.5*(1 + cos(2*pi*s)) = |Delta|/peak`` for ``s = (t-center)/T``
        gives ``s = +- arccos(2|Delta|/peak - 1) / (2*pi)``, so the window is

            width = duration * arccos(2|Delta|/peak - 1) / pi

        which is strictly proportional to ``duration`` at fixed peak and fixed
        detuning.  Shortening the pulse can therefore never buy a clock more
        tongue time — the finite-time falsification hinges on that.

        These are *continuous* times, not sampled ones.  The classification of
        a clock at the instant it enters the tongue is done against the entry
        time returned here rather than against the first stored sample after
        it, because those are not the same moment.
        """
        a = abs(require_finite(detuning, "detuning"))
        if not (self.peak > a):
            return None
        s = np.arccos(2.0 * a / self.peak - 1.0) / TWO_PI
        half = float(s) * self.duration
        return (self.center - half, self.center + half)

    def eligibility_duration(self, detuning) -> float:
        """Length of the eligibility window; ``0.0`` if never eligible."""
        window = self.eligibility_window(detuning)
        return 0.0 if window is None else window[1] - window[0]

    def support(self) -> Tuple[float, float]:
        return (self.center - 0.5 * self.duration,
                self.center + 0.5 * self.duration)


# ---------------------------------------------------------------------------
# 5. Deterministic steppers
# ---------------------------------------------------------------------------
# The steppers take an already-validated coupling callable and already-validated
# arrays.  They are the inner loop, so they use the unvalidated kernels.
def _step_euler(phi, t, dt, detuning, coupling_of_t):
    """Unvalidated kernel behind :func:`step_euler`; called once per step."""
    return phi + dt * _adler_drift(phi, coupling_of_t(t), detuning)


def _step_rk4(phi, t, dt, detuning, coupling_of_t):
    """Unvalidated kernel behind :func:`step_rk4`; called once per step."""
    h = 0.5 * dt
    k1 = _adler_drift(phi, coupling_of_t(t), detuning)
    k2 = _adler_drift(phi + h * k1, coupling_of_t(t + h), detuning)
    k3 = _adler_drift(phi + h * k2, coupling_of_t(t + h), detuning)
    k4 = _adler_drift(phi + dt * k3, coupling_of_t(t + dt), detuning)
    return phi + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)


def _validated_step(phi, t, dt, detuning, coupling_of_t):
    """Shared boundary check for the two public steppers.

    The callback's *output* is checked as well as its callability.  A source
    that returns NaN — an unvalidated pulse, a user lambda with a bug — would
    otherwise poison the step and hand back a NaN phase that looks like a
    result.  The wrapper below re-checks at every Runge-Kutta stage time, which
    is affordable precisely because these public steppers are not the inner
    loop; the integrator calls the private kernels instead.
    """
    require_callable(coupling_of_t, "coupling_of_t")

    def checked(t_stage):
        return require_nonnegative_array(coupling_of_t(t_stage),
                                         "coupling_of_t(t)")

    return (require_finite_array(phi, "phi"), require_finite(t, "t"),
            require_finite(dt, "dt"),
            require_finite_array(detuning, "detuning"), checked)


def step_euler(phi, t, dt, detuning, coupling_of_t):
    """Explicit Euler.

    Kept because the later stochastic stage integrates with Euler-Maruyama, and
    it is useful to be able to run the deterministic limit of the *same* scheme.
    It is not the default here: for a calibration whose job is to measure a
    relaxation rate to five digits, first-order accuracy is the wrong tool.
    """
    phi, t, dt, detuning, checked = _validated_step(
        phi, t, dt, detuning, coupling_of_t)
    return _step_euler(phi, t, dt, detuning, checked)[()]


def step_rk4(phi, t, dt, detuning, coupling_of_t):
    """Classical fourth-order Runge-Kutta.

    The coupling is re-evaluated at the interior stage times, so a rising or
    falling pulse edge is resolved rather than frozen at the step start.
    """
    phi, t, dt, detuning, checked = _validated_step(
        phi, t, dt, detuning, coupling_of_t)
    return _step_rk4(phi, t, dt, detuning, checked)[()]


#: The integrator dispatches to the *unvalidated* kernels: it has already
#: checked everything at its own boundary, and re-checking once per step would
#: dominate the run.  The public wrappers above exist for callers outside it.
STEPPERS = {"euler": _step_euler, "rk4": _step_rk4}
