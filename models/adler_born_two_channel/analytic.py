#!/usr/bin/env python3
"""
analytic.py — the flat-spectrum rate-weighted tongue prediction, and nothing else.

This module is the *companion control*, not part of the experiment.  It works
out, in closed form, what the rate-weighted tongue picture predicts, so that a
later stochastic simulation can be compared against a number that was computed
without ever looking at the simulation.

The prediction
--------------
A clock with detuning ``Delta`` sitting inside a stationary tongue of coupling
``K`` relaxes towards its stable phase at the local rate

    lambda(K, Delta) = sqrt(K^2 - Delta^2)

(the slope of ``Delta - K sin(phi)`` at the stable fixed point, ``K cos(phi*)``).
Summing that rate over a flat, evenly spaced detuning grid of spacing ``h`` and
multiplying by ``h`` is a midpoint quadrature of

    integral over |Delta| < K of sqrt(K^2 - Delta^2) dDelta  =  (pi/2) * K^2

— the semicircle.  So the flat-spectrum flux is a universal constant times the
*square* of an amplitude-linear coupling.  That quadratic is the entire reason
the rate-weighted tongue is a Born candidate at all.

Two independent effects produce it, and they must not be conflated: the tongue
gets **wider** in proportion to ``K`` (more eligible clocks), and its interior
gets **faster** in proportion to ``K`` (each eligible clock relaxes sooner).
Width alone would give an amplitude-*linear* law.

Isolation
---------
Nothing in ``model.py`` or ``simulate.py`` may import or call this module.
``verify.py`` enforces that by parsing their source.  This module returns the
favour: it imports no dynamics, and it re-derives the eligibility predicate
``K > |Delta|`` locally rather than importing it, so that ``verify.py`` can
cross-check the two implementations against each other on a random battery and
catch any silent drift.

The one shared import is ``validation.py``, which holds the argument checks.
That module contains no physics — only type and finiteness tests over anonymous
values — so it cannot carry information between the prediction and the
simulation, and ``verify.py`` asserts that emptiness rather than trusting it.

What this module does NOT show
------------------------------
That the quadratic is the right *shape* is arithmetic about the Adler equation.
It says nothing about whether a noisy race actually commits at a rate
proportional to this flux, and therefore nothing about Born's rule.  That is
what the later stochastic stages are for, and they are allowed to answer "no".
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence, Tuple

import numpy as np

from .validation import (require_finite, require_finite_array,
                         require_instance, require_nonnegative,
                         require_nonnegative_array, require_positive_array,
                         require_positive_int)

__all__ = [
    "SEMICIRCLE_CONSTANT",
    "eligible_mask",
    "local_relaxation_rate",
    "efolds_available",
    "FlatGrid",
    "flat_detuning_grid",
    "tongue_rate_sum",
    "tongue_rate_flux",
    "fit_log_log_exponent",
    "safe_coupling_range",
    "refinement_study",
]

#: The continuum answer: ``flux(K) = SEMICIRCLE_CONSTANT * K**2``.
SEMICIRCLE_CONSTANT = 0.5 * np.pi


def _eligible_mask(coupling, detunings):
    """Unvalidated kernel behind :func:`eligible_mask`."""
    return np.asarray(coupling, dtype=float) > np.abs(
        np.asarray(detunings, dtype=float))


def _local_relaxation_rate(coupling, detunings):
    """Unvalidated kernel behind :func:`local_relaxation_rate`."""
    K = np.asarray(coupling, dtype=float)
    D = np.asarray(detunings, dtype=float)
    elig = _eligible_mask(K, D)
    inside = np.where(elig, K * K - D * D, 0.0)
    return np.where(elig, np.sqrt(inside), np.nan)


def eligible_mask(coupling, detunings):
    """``K > |Delta|``, re-derived locally to keep this module standalone.

    Deliberately *not* imported from ``model.py``; see the module docstring.
    Both arguments must be finite and the coupling non-negative: a NaN or
    negative coupling would make this predicate false everywhere and so imitate
    a physical no-admission result.  Detunings are signed; couplings are not.
    """
    return _eligible_mask(require_nonnegative_array(coupling, "coupling"),
                          require_finite_array(detunings, "detunings"))[()]


def local_relaxation_rate(coupling, detunings):
    """``sqrt(K^2 - Delta^2)`` inside the tongue, NaN outside.

    NaN rather than zero: a clock outside the tongue has no stable fixed point,
    so it has no relaxation rate at all.  Reporting zero would quietly turn "not
    defined" into "defined and slow".  The summation helpers below mask
    explicitly instead of leaning on NaN arithmetic.

    That NaN is an *output* meaning "no rate exists here".  A NaN *input* is a
    different thing and is rejected, so the two cannot be confused.  A negative
    coupling is rejected too: this expression squares it, so it would return a
    positive, entirely physical-looking rate for a configuration the model does
    not define.
    """
    return _local_relaxation_rate(
        require_nonnegative_array(coupling, "coupling"),
        require_finite_array(detunings, "detunings"))[()]


def efolds_available(times, couplings, detuning) -> float:
    """Integral of the local relaxation rate over the eligible part of a pulse.

    This is the number of e-folds a clock is *offered* while its tongue is open:

        A = integral over the eligible window of sqrt(K(t)^2 - Delta^2) dt

    It is the natural finite-time budget.  A clock that starts an order-one
    distance from its stable phase needs several e-folds to arrive, so ``A``
    of order one or less means the pulse closes before the clock can lock —
    the finite-time failure mode, quantified.

    It lives in ``analytic.py`` because it is built from the analytic relaxation
    rate.  The dynamics never see it; ``verify.py`` computes it alongside the
    measured trajectories and compares.  ``times`` must be evenly spaced enough
    for the trapezoid rule to be meaningful.
    """
    t = np.atleast_1d(require_finite_array(times, "times"))
    # Non-negative, not merely finite: this integrand squares the coupling, so
    # a negative one integrates to a positive, physical-looking e-fold budget.
    K = np.atleast_1d(require_nonnegative_array(couplings, "couplings"))
    D = require_finite(detuning, "detuning")
    if t.shape != K.shape:
        raise ValueError(f"times and couplings must have the same shape, got "
                         f"{t.shape} and {K.shape}")
    lam = np.sqrt(np.clip(K * K - D * D, 0.0, None))   # zero where ineligible
    return float(np.trapezoid(lam, t))


@dataclass(frozen=True)
class FlatGrid:
    """A fixed midpoint grid of detunings, shared by every coupling.

    ``half_width`` is the spectral support; the grid runs over
    ``[-half_width, +half_width]`` with ``n_points`` equally spaced midpoints
    and spacing ``h = 2*half_width/n_points``.

    Midpoints, not endpoints: an endpoint grid would place a sample exactly at
    ``Delta = 0`` and exactly at the support edges, and — worse — would make the
    result depend on whether the tongue edge happens to land on a node.  The
    grid is held fixed while ``K`` is swept, so that the coupling dependence
    being measured is the tongue's, not the sampler's.
    """

    half_width: float
    n_points: int

    def __post_init__(self):
        half_width = require_finite(self.half_width, "half_width")
        n_points = require_positive_int(self.n_points, "n_points", minimum=2)
        if not half_width > 0.0:
            raise ValueError(f"half_width must be positive, got {half_width!r}")
        # Store the normalized values, so a frozen grid can never hold
        # something that merely converted once inside the constructor.
        object.__setattr__(self, "half_width", half_width)
        object.__setattr__(self, "n_points", n_points)

    @property
    def spacing(self) -> float:
        return 2.0 * self.half_width / self.n_points

    @property
    def detunings(self) -> np.ndarray:
        h = self.spacing
        return -self.half_width + (np.arange(self.n_points) + 0.5) * h


def flat_detuning_grid(half_width: float, n_points: int) -> Tuple[np.ndarray, float]:
    """Convenience wrapper: ``(detunings, spacing)`` for a :class:`FlatGrid`."""
    grid = FlatGrid(half_width, n_points)
    return grid.detunings, grid.spacing


def tongue_rate_sum(coupling, grid: FlatGrid) -> float:
    """Bare sum of local relaxation rates over the eligible clocks of ``grid``.

    This is the finite-population quantity — "how much total locking speed does
    this coupling command?" — and it scales as ``K**2 / h``.
    """
    coupling = require_nonnegative(coupling, "coupling")
    require_instance(grid, FlatGrid, "grid")
    D = grid.detunings
    elig = _eligible_mask(coupling, D)
    if not np.any(elig):
        return 0.0
    return float(np.sum(_local_relaxation_rate(coupling, D[elig])))


def tongue_rate_flux(coupling, grid: FlatGrid) -> float:
    """Spacing-weighted sum: the midpoint quadrature of the semicircle.

    Converges to ``SEMICIRCLE_CONSTANT * coupling**2`` as ``h -> 0``, provided
    the tongue fits strictly inside the support.
    """
    require_instance(grid, FlatGrid, "grid")
    return grid.spacing * tongue_rate_sum(coupling, grid)


def fit_log_log_exponent(couplings, values) -> float:
    """Unconstrained least-squares slope of ``log(value)`` against ``log(K)``.

    Unconstrained is the point.  The fit is allowed to come back with 1.0, or
    1.7, or anything else; nothing here forces it towards 2.

    Both inputs must be finite and strictly positive.  Non-positive entries are
    **rejected, not dropped**.  Dropping them would be silent post-selection:
    ``fit_log_log_exponent([1, -2, 3], [1, 4, 9])`` would quietly discard the
    middle point and report an exponent of about 2.0 for a sweep the caller
    believed had three points.  A fit that answers a question about a subset
    must not present itself as answering the whole one.
    """
    K = np.atleast_1d(require_positive_array(couplings, "couplings"))
    v = np.atleast_1d(require_positive_array(values, "values"))
    if K.shape != v.shape:
        raise ValueError(f"couplings and values must have the same shape, got "
                         f"{K.shape} and {v.shape}")
    if K.size < 2:
        raise ValueError(f"a log-log fit needs at least two points, got "
                         f"{K.size}")
    slope, _ = np.polyfit(np.log(K), np.log(v), 1)
    return float(slope)


#: Declared safe range for the exponent fit, as multiples of the grid spacing
#: and fractions of the spectral support.  Below ``MIN_TONGUE_CELLS`` the tongue
#: spans too few grid cells for the quadrature to be meaningful; above
#: ``MAX_SUPPORT_FRACTION`` the tongue starts to feel the finite support.
MIN_TONGUE_CELLS = 200.0
MAX_SUPPORT_FRACTION = 0.25


def safe_coupling_range(grid: FlatGrid) -> Tuple[float, float]:
    """``(K_min, K_max)`` over which the flat-grid flux is trustworthy.

    Declared *before* looking at any result, and reported alongside every fit,
    so that a good exponent cannot be manufactured by quietly cropping the
    coupling sweep.
    """
    require_instance(grid, FlatGrid, "grid")
    return (MIN_TONGUE_CELLS * grid.spacing,
            MAX_SUPPORT_FRACTION * grid.half_width)


def refinement_study(couplings: Sequence[float],
                     half_widths: Sequence[float],
                     n_points: Sequence[int]):
    """Refine the grid (and support) and watch the prediction converge.

    Returns a list of dicts, one per refinement level, each carrying the grid
    spacing, the worst relative deviation of ``flux / K**2`` from ``pi/2`` over
    the supplied couplings, and the unconstrained fitted exponent.
    """
    # Strictly positive: these couplings are divided by and log-fitted below.
    K = np.atleast_1d(require_positive_array(couplings, "couplings"))
    half_widths = list(half_widths)
    n_points = list(n_points)
    if len(half_widths) != len(n_points):
        raise ValueError(f"half_widths and n_points must have the same length, "
                         f"got {len(half_widths)} and {len(n_points)}")
    rows = []
    for half_width, n in zip(half_widths, n_points):
        grid = FlatGrid(require_finite(half_width, "half_width"),
                        require_positive_int(n, "n_points", minimum=2))
        flux = np.array([tongue_rate_flux(k, grid) for k in K])
        ratio = flux / (K * K)
        rows.append({
            "half_width": grid.half_width,
            "n_points": grid.n_points,
            "spacing": grid.spacing,
            "safe_range": safe_coupling_range(grid),
            "flux": flux,
            "worst_relative_error": float(
                np.max(np.abs(ratio / SEMICIRCLE_CONSTANT - 1.0))),
            "exponent": fit_log_log_exponent(K, flux),
        })
    return rows
