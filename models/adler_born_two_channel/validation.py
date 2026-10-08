#!/usr/bin/env python3
"""
validation.py — argument checks for the public boundary, and nothing else.

Why this is its own module
--------------------------
``model.py`` and ``analytic.py`` are deliberately kept structurally isolated
from each other: the dynamics must not be able to see the analytic prediction
they are meant to be tested against, and the analytic prediction must not be a
restatement of the simulation.  ``verify.py`` enforces that by parsing both.

Both nevertheless need the same argument checks.  Putting them here lets both
import them without either importing the other.  This module contains **no
physics** — no coupling, no detuning, no phase, no tongue, no rate — only type
and finiteness checks over anonymous values, so sharing it cannot carry
information between the two sides.  ``verify.py`` asserts that emptiness rather
than trusting it.

Why the checks are strict
-------------------------
A NaN is not a small number and infinity is not a large one.  Unchecked they do
not crash; they propagate quietly and then *look like physics*.  A NaN coupling
makes ``K > |Delta|`` false everywhere, so an invalid run reports that no clock
was ever admitted to the tongue — a clean, plausible, entirely fictitious
no-admission result.  The whole point of these functions is that an invalid
input can never impersonate a physical conclusion.

Booleans and strings are rejected as deliberately as NaN.  ``bool`` is a
subclass of ``int`` and ``float("1")`` succeeds, so both would otherwise slip
through a plain ``float()`` conversion and become silent physical parameters.

Cost, and the private-kernel convention
---------------------------------------
These checks scan their input, which is unaffordable four times per Runge-Kutta
step.  So each validated public function in this package is a thin boundary over
an unvalidated private kernel, and the integrator calls the kernel.  The rule is
that a kernel may only be called with values some public boundary has already
checked.
"""
from __future__ import annotations

import numpy as np

__all__ = [
    "require_finite",
    "require_finite_array",
    "require_nonnegative",
    "require_nonnegative_array",
    "require_positive_array",
    "require_boolean_array",
    "require_real_allowing_nan",
    "require_positive_int",
    "require_instance",
    "require_exact_type",
    "require_callable",
]

#: Types that must never be silently coerced into a physical parameter.
_NON_REAL_SCALARS = (bool, np.bool_, str, bytes, bytearray,
                     complex, np.complexfloating)


def require_finite(value, name: str) -> float:
    """Return ``value`` as a float, or raise naming ``name`` as the culprit.

    ``TypeError`` if it is not a real number at all — a bool, a string, a
    complex, an array — and ``ValueError`` if it is a real number that is NaN or
    infinite.  The two are different mistakes and get different exceptions.
    """
    if isinstance(value, _NON_REAL_SCALARS):
        raise TypeError(f"{name} must be a finite real number, got "
                        f"{type(value).__name__} {value!r}")
    if isinstance(value, np.ndarray):
        raise TypeError(f"{name} must be a real scalar, got an array of shape "
                        f"{value.shape}")
    try:
        out = float(value)
    except (TypeError, ValueError):
        raise TypeError(f"{name} must be a finite real number, got "
                        f"{type(value).__name__} {value!r}") from None
    if not np.isfinite(out):
        raise ValueError(f"{name} must be finite, got {value!r}")
    return out


def require_finite_array(values, name: str) -> np.ndarray:
    """Return ``values`` as a finite float array, or raise naming ``name``.

    Accepts scalars and array-likes, so ordinary vectorized use is unaffected;
    rejects boolean, string, object and complex arrays outright rather than
    letting NumPy coerce them into plausible-looking floats.

    The input's **shape is preserved**, including 0-d for a scalar.  That is
    deliberate: the public functions end with ``[()]``, which turns a 0-d result
    back into a NumPy scalar, so ``stable_phase(1.0, 0.5)`` keeps returning a
    scalar rather than a one-element array.  Callers that genuinely need at
    least one dimension say so themselves.
    """
    if isinstance(values, (str, bytes, bytearray)):
        raise TypeError(f"{name} must be real numbers, got "
                        f"{type(values).__name__} {values!r}")
    raw = np.asarray(values)
    if raw.dtype.kind == "b":
        raise TypeError(f"{name} must be real numbers, got booleans")
    if raw.dtype.kind == "c":
        raise TypeError(f"{name} must be real numbers, got complex values")
    if raw.dtype.kind in "USO":
        raise TypeError(f"{name} must be real numbers, got dtype {raw.dtype}")
    try:
        out = raw.astype(float)
    except (TypeError, ValueError):
        raise TypeError(f"{name} must be real numbers, got dtype "
                        f"{raw.dtype}") from None
    finite = np.isfinite(out)
    if not bool(np.all(finite)):
        bad = np.atleast_1d(out)[~np.atleast_1d(finite)]
        raise ValueError(f"{name} must be finite everywhere; {bad.size} of "
                         f"{out.size} entries are not, e.g. {bad.flat[0]!r}")
    return out


def require_positive_int(value, name: str, minimum: int = 1) -> int:
    """Return ``value`` as an int >= ``minimum``, or raise naming ``name``.

    ``bool`` is rejected even though it is a subclass of ``int``: passing
    ``True`` here is a mistake, not a request for one.  Floats are rejected
    rather than truncated, because a silently truncated cadence is a wrong
    answer that looks like a right one.
    """
    if isinstance(value, (bool, np.bool_)):
        raise TypeError(f"{name} must be an integer, got "
                        f"{type(value).__name__} {value!r}")
    if not isinstance(value, (int, np.integer)):
        raise TypeError(f"{name} must be an integer, got "
                        f"{type(value).__name__} {value!r}")
    out = int(value)
    if out < minimum:
        raise ValueError(f"{name} must be at least {minimum}, got {out!r}")
    return out


def require_nonnegative(value, name: str) -> float:
    """Return ``value`` as a float >= 0, or raise naming ``name``.

    For scalars the package treats as a *coupling*.  A negative coupling is not
    a weak coupling: it flips the sign of the restoring term, which moves the
    stable Adler fixed point by pi.  Silently accepting one would return a
    plausible number for a configuration the model does not define.
    """
    out = require_finite(value, name)
    if out < 0.0:
        raise ValueError(f"{name} must be non-negative, got {out!r}")
    return out


def require_nonnegative_array(values, name: str) -> np.ndarray:
    """Return ``values`` as a finite float array with every entry >= 0.

    Used for couplings.  The failure this prevents is quiet rather than loud:
    a negative coupling propagates through ``K > |Delta|`` as a plain ``False``,
    through ``sqrt(K^2 - Delta^2)`` as a NaN, and through anything that squares
    it as a *positive, physical-looking* result.
    """
    out = require_finite_array(values, name)
    negative = out < 0.0
    if bool(np.any(negative)):
        bad = np.atleast_1d(out)[np.atleast_1d(negative)]
        raise ValueError(f"{name} must be non-negative everywhere; {bad.size} "
                         f"of {out.size} entries are negative, e.g. "
                         f"{bad.flat[0]!r}")
    return out


def require_positive_array(values, name: str) -> np.ndarray:
    """Return ``values`` as a finite float array with every entry > 0.

    For inputs to a log-log fit.  Dropping non-positive entries instead would
    be silent post-selection: the fit would quietly answer a question about a
    subset the caller never asked about, and report it as though it were the
    whole sweep.
    """
    out = require_finite_array(values, name)
    bad_mask = out <= 0.0
    if bool(np.any(bad_mask)):
        bad = np.atleast_1d(out)[np.atleast_1d(bad_mask)]
        raise ValueError(f"{name} must be strictly positive everywhere (a "
                         f"log-log fit has no other domain); {bad.size} of "
                         f"{out.size} entries are not, e.g. {bad.flat[0]!r}")
    return out


def require_boolean_array(values, name: str) -> np.ndarray:
    """Return ``values`` as a genuine boolean array, or raise naming ``name``.

    Coercion is refused rather than performed.  ``np.asarray(['False','True'])``
    has string dtype and would coerce to *all true*, and ``[0, 2, 0]`` would
    coerce to ``[False, True, False]`` — both plausible, neither what the caller
    meant.  A predicate flag has to arrive as a predicate.
    """
    if isinstance(values, (str, bytes, bytearray)):
        raise TypeError(f"{name} must be a boolean array, got "
                        f"{type(values).__name__} {values!r}")
    raw = np.asarray(values)
    if raw.dtype.kind != "b":
        raise TypeError(f"{name} must be a boolean array, got dtype "
                        f"{raw.dtype}; numeric and string flags are refused "
                        f"rather than coerced")
    return raw


def require_real_allowing_nan(values, name: str) -> np.ndarray:
    """Return ``values`` as a float array in which **NaN alone** may appear.

    The one place the package accepts a non-finite input, and it accepts
    exactly one kind for exactly one reason: NaN is the documented sentinel for
    "no target existed at this sample", which
    ``Trajectory.tracking_error()`` emits wherever a clock was ineligible.
    That is missing data, and a fit over the eligible band must skip it.

    Infinity is not missing data — it is a numerical blow-up, and silently
    filtering it would hide the very failure it signals.  Booleans, strings and
    complex arrays are wrong-kind input.  All are rejected here, before the
    finite-band mask is applied, so the carve-out stays a carve-out.
    """
    if isinstance(values, (str, bytes, bytearray)):
        raise TypeError(f"{name} must be real numbers, got "
                        f"{type(values).__name__} {values!r}")
    raw = np.asarray(values)
    if raw.dtype.kind == "b":
        raise TypeError(f"{name} must be real numbers, got booleans")
    if raw.dtype.kind == "c":
        raise TypeError(f"{name} must be real numbers, got complex values")
    if raw.dtype.kind in "USO":
        raise TypeError(f"{name} must be real numbers, got dtype {raw.dtype}")
    try:
        out = raw.astype(float)
    except (TypeError, ValueError):
        raise TypeError(f"{name} must be real numbers, got dtype "
                        f"{raw.dtype}") from None
    infinite = np.isinf(out)
    if bool(np.any(infinite)):
        n_inf = int(np.count_nonzero(infinite))
        raise ValueError(f"{name} may contain NaN as the missing-target "
                         f"sentinel, but not infinities: {n_inf} of {out.size} "
                         f"entries are infinite, which indicates a numerical "
                         f"blow-up rather than a missing sample")
    return out


def require_callable(value, name: str):
    """Return ``value`` if it is callable, else raise naming ``name``."""
    if not callable(value):
        raise TypeError(f"{name} must be callable, got "
                        f"{type(value).__name__}")
    return value


def require_instance(value, expected, name: str):
    """Return ``value`` if it is an instance of ``expected``, else raise."""
    if not isinstance(value, expected):
        wanted = getattr(expected, "__name__", str(expected))
        raise TypeError(f"{name} must be a {wanted}, got "
                        f"{type(value).__name__}")
    return value


def require_exact_type(value, expected, name: str):
    """Return ``value`` if its type is *exactly* ``expected``, else raise.

    Stricter than :func:`require_instance`, and used where a subclass would be
    a hole rather than a convenience.

    An ``isinstance`` boundary accepts a subclass, and a subclass may declare
    fields this package never declared and therefore never validates.  That is
    the same hole as an unrecognised configuration key, reached by a different
    route: a frozen subclass adding ``predictor: object`` satisfies every
    ``isinstance`` check, survives ``dataclasses.replace`` — which preserves the
    dynamic class — and arrives at the code that was supposed to be unable to
    see a predictor, carrying one.

    So at boundaries where the whole guarantee is "this object holds these
    fields and nothing else", the type has to be the type.  A subclass gets its
    own message, because "got RawEventConfig" would otherwise be a baffling
    thing to read on a ``RawEventConfig`` subclass.
    """
    wanted = getattr(expected, "__name__", str(expected))
    if type(value) is not expected:
        actual = type(value).__name__
        if isinstance(value, expected):
            raise TypeError(
                f"{name} must be exactly {wanted}, got the subclass {actual}; "
                f"a subclass can carry fields this package never declared and "
                f"so never checked, which is the hole an unrecognised key "
                f"would have been")
        raise TypeError(f"{name} must be exactly {wanted}, got {actual}")
    return value
