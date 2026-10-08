#!/usr/bin/env python3
"""
stochastic.py — the keyed physical white-noise source, and the nested Brownian
tree that splits one of its increments at an exact off-grid crossing.

What this module is
-------------------
One thing only: a *reproducible supply of Gaussian phase kicks*.  It is the
effective bath of the single-channel experiment, and it is deliberately the
dullest possible one — independent, zero-mean, white, with a variance that
depends on nothing but the declared phase-diffusion strength and the duration
of the interval.

    kick over an interval of duration h  ~  N(0, 2 * D_phase * h)

**White noise is a minimal effective bath and nothing more.**  No microscopic
origin is claimed for it, no quantum vacuum, no absorber level structure, no
measurement backaction.  It is the smallest assumption that makes an Adler
clock stochastic at all, and the reason it is written this way — with the noise
law incapable of seeing coupling, mismatch, amplitude, lock tolerance or dwell
— is that a bath which could see any of those could manufacture the coupling
dependence the experiment exists to measure.  Because the noise is additive,
the Ito-versus-Stratonovich distinction does not change this model; that is a
statement about the equation, not an additional physical conclusion.

White noise is also an infinite-bandwidth idealization.  The finest timestep
supplies the numerical cutoff, which is why timestep convergence is a
scientific requirement of the experiment and not merely a software check.

What this module is NOT
----------------------
There is no drift here, no Adler equation, no eligibility, no lock band, no
dwell clock, no commitment, no outcome and no statistics.  Nothing in this file
produces an event, so nothing in it can produce an event *frequency*, so nothing
in it can agree or disagree with Born's rule.  It supplies numbers to a stepper
that does not exist yet.

The two guarantees
------------------
**1. A value belongs to a key, not to a run.**  Every elementary kick is
addressed by an immutable key — dataset namespace, finest mesh, trial, clock,
finest-step index — and is a pure function of that key.  Batch size, block
window, iteration order, how many clocks are asked for, how many trials stopped
early, and which of the two paired models (full Adler or the width-only
control) is consuming the stream are all *absent from the key*, so none of them
can change a value.  That is what makes a paired refinement comparison a
comparison of timestep error rather than of two different noise realizations.

**2. Coarse increments are sums, never separate draws.**  There is exactly one
physical tree.  A coarser refinement level adds up the same finest leaves; an
exact tongue crossing inside a finest step subdivides that step's already
assigned kick.  Both directions refer to the same realization, which is the only
way "the same noise at two timesteps" and "the same noise for two models" can
mean anything.

The conditional Brownian split
------------------------------
An exact deterministic tongue entry or exit falls inside an assigned uniform
timestep.  The realized kick for that step has already been drawn, and the two
sub-intervals need kicks that (a) have the right variances for their own
durations and (b) still add up to the drawn parent.  Dividing the parent in
proportion to the two durations satisfies (b) and fails (a): it would give the
left child variance ``alpha^2 * 2 D h`` where the Brownian path requires
``alpha * 2 D h``.

The correct construction conditions on the parent.  With ``sigma^2 = 2 D``,
parent kick ``X ~ N(0, sigma^2 h)``, crossing fraction ``alpha`` and an
independent standard normal ``Z``:

    left  = alpha * X + sigma * sqrt(h * alpha * (1 - alpha)) * Z
    right = X - left

so that ``Var(left) = alpha sigma^2 h``, ``Var(right) = (1-alpha) sigma^2 h``,
``Cov(left, right) = 0`` — and, since the pair is jointly Gaussian, the two
children are unconditionally independent — while ``left + right = X`` by
construction.  The second child is the *residual*; sampling it independently
would give the right marginals and destroy the parent, and the verifier's
mutation battery says so by requiring that version to fail.

``Z`` is additional information: it is what reveals where the path went inside
the parent interval, and it is physical-noise randomness.  It is not a
bridge-audit uniform and not a control-specific correction, so its key lives in
the same physical namespace as the leaves and excludes the model label.

Floating point, stated honestly
-------------------------------
``left + right = X`` is exact in real arithmetic and *not* bit-exact in binary
floating point: ``fl(left + fl(X - left))`` can differ from ``X`` in the last
place or two.  So the conservation contract is scale-aware rather than
bit-exact — see :data:`CONSERVATION_ULPS` — with two exceptions where bit-exact
is the right claim and is asserted as one: a zero-duration child and a
zero-diffusion tree are exactly, identically zero.

How a key becomes a number
--------------------------
The mapping is fixed here so that it can be recorded in a manifest and
reproduced:

    canonical ASCII key string
      -> SHA-256 digest (32 bytes)
      -> four little-endian 64-bit words
      -> numpy SeedSequence(entropy=words)
      -> Philox counter-based bit generator
      -> Generator.standard_normal(count)

Leaves are generated a frozen :data:`LEAF_BLOCK` at a time because generating
them one at a time would be a hundred times slower; the block size is a
*constant of the schema*, not the caller's batch size, so a caller who asks for
three values or three million gets the same values either way.

Reproducibility of that last arrow is environment-scoped and this is the one
place to say so plainly: NumPy guarantees its bit generators' streams, and does
*not* guarantee that ``standard_normal`` will turn a stream into the same
normals in a future major version.  Identical keys therefore give identical
kicks **in a recorded (Python, NumPy) environment**, which is what the run
manifest exists to record, and identical *structure* — independence, batching
invariance, parent conservation — in any environment at all.
"""
from __future__ import annotations

import hashlib
import math
from dataclasses import dataclass
from typing import Iterator, Tuple

import numpy as np

from .validation import (require_exact_type, require_finite,
                         require_nonnegative, require_positive_int)

__all__ = [
    "CONSERVATION_ULPS",
    "EDGE_LABELS",
    "KEY_SCHEMA",
    "LEAF_BLOCK",
    "Crossing",
    "FinestMesh",
    "PhaseNoiseStream",
    "TreeLeaf",
    "conservation_tolerance",
]

#: The root of every physical-noise key.  Versioned, because the derivation
#: below is part of the experiment's frozen manifest: a change to the key
#: grammar or to the digest is a change of dataset, and it has to be visible as
#: one rather than silently producing different numbers under the same names.
#:
#: ``physical`` is the namespace segment that separates this tree from the
#: diagnostic moving-band audit, which draws its own uniforms and must never be
#: able to reach these values.  That layer is a later ticket; the segment is
#: reserved now so the separation exists before the code that needs it.
#:
#: **v3.**  The duration law moved: a step's kick is scaled by the *realized*
#: width of its canonical interval, ``step_bounds(i)[1] - step_bounds(i)[0]``,
#: rather than by the separately declared nominal ``step``.  On
#: ``FinestMesh(1e15, 0.2)`` those differ — the realized width is ``0.25`` —
#: and the tree was internally inconsistent about which one it meant: the
#: ``TreeLeaf.duration``, the split fractions and the coarse aggregation all
#: said ``0.25`` while the kick variance said ``0.2``.  The Brownian scale is
#: also computed by a stable factorization now, which changes the last bit even
#: where the two durations agree.  Both are changes to the key-to-kick mapping,
#: so the version moves and **every keyed value changes with it**.  No run has
#: been recorded under v1 or v2; nothing needs migrating.
#:
#: **v2 was the address bump.**  v1 published a "leaf key" that ended in
#: ``block=<index // 1024>``, which is the seed of the *block* the value is
#: generated in and not the address of the value.  All 1024 leaves of a block
#: shared it, so a manifest storing a v1 leaf key could not say which kick it
#: meant.  v2 separates the two: :meth:`PhaseNoiseStream.leaf_key` is the public
#: per-leaf address ending in ``step=<index>``, and the block seed is private.
#: The comment above says a change of key grammar must be visible as a change of
#: dataset rather than silently producing different numbers under the same
#: names, so the version moves and **every keyed value changes with it**.  No
#: run has been recorded under v1; nothing needs migrating.
KEY_SCHEMA = "dk-phase-noise/v3/physical"

#: How many finest leaves one keyed block holds.  A constant of the key
#: schema, deliberately not a tunable: the *value* assigned to a step is fixed
#: by its key, and the block is only how the generator gets there efficiently.
#: A caller streaming in windows of 7 and a caller asking for a million at once
#: receive identical arrays.
LEAF_BLOCK = 1024

#: The scale-aware bound on parent/child conservation, in units in the last
#: place.  ``right = fl(parent - left)`` and ``fl(left + right)`` differ from
#: ``parent`` by at most ``u * (|parent| + |right|)`` to first order, with
#: ``u = 2^-53``; at the largest magnitude involved that is at most two ulp, and
#: four is the margin actually asserted.  It is a deterministic bound and never
#: a statistical tolerance: a split that conserves the parent only on average
#: is a broken split.
CONSERVATION_ULPS = 4

#: Which side of the pulse a crossing is on.  Part of the crossing identity
#: because a rising entry and a falling exit at the same clock are different
#: physical events and must not share a split normal.
EDGE_LABELS = ("rising", "falling")

#: The one identifier grammar, shared by every integer that reaches a key:
#: trial, clock and finest-step index, whether they arrive one at a time or in
#: an array.  ``[0, 2**64 - 1]`` is the exact range of an unsigned 64-bit
#: identifier, which is the widest thing a caller can hand this module without
#: a lossy conversion somewhere.
#:
#: It exists because the scalar and vector boundaries used to disagree: the
#: scalar path accepted any non-negative Python integer, while the array path
#: checked for negativity and *then* cast to ``int64``, so ``np.uint64(2**63)``
#: passed the check and arrived at the key as ``-9223372036854775808``.  A
#: clock's noise stream was therefore addressable under an identifier the
#: scalar boundary would have refused.  One bound, both paths, no cast.
_MAX_IDENTIFIER = (1 << 64) - 1

#: Characters a dataset namespace may contain.  A namespace is serialized into
#: an ASCII key, so it has to be a plain token that round-trips exactly.
_NAMESPACE_CHARS = frozenset(
    "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-")

#: Longest permitted namespace.
_MAX_NAMESPACE = 64

#: Characters a Brownian-tree node path may contain: left and right ancestry.
_NODE_CHARS = frozenset("LR")

#: Deepest permitted node path.  A finest step containing more than this many
#: nested crossings is a configuration error rather than a deep tree.
_MAX_NODE_DEPTH = 32

#: Characters a pulse-manifest digest may contain, and its permitted lengths.
_DIGEST_CHARS = frozenset("0123456789abcdef")
_MIN_DIGEST, _MAX_DIGEST = 8, 64

#: How many finest steps :meth:`PhaseNoiseStream._require_window` checks at a
#: time.  The whole window is validated before any randomness is drawn, but a
#: production request can name more steps than fit in memory, so the scan is
#: chunked: eager in coverage, bounded in space.
_VALIDATION_CHUNK = 1 << 16

#: How far :meth:`FinestMesh.step_containing` may correct its initial estimate
#: before declaring the mesh unroutable.  The estimate comes from a division
#: that is inexact by at most an ulp or two, so one correction is the realistic
#: worst case; four is the margin, and a bound is what stops an unresolvable
#: mesh from spinning instead of failing.
_ROUTING_STEPS = 4

#: The name of the empty node path in a key.  An empty string would serialize
#: to ``node=/`` and make two different keys visually similar; the root of the
#: tree is written out.
_ROOT_NODE = "root"


# ---------------------------------------------------------------------------
# Private validation helpers.  Physics-free checks live in validation.py; these
# are specific to this module's own vocabulary and stay here.
# ---------------------------------------------------------------------------
def _require_token(value, name: str) -> str:
    """Return ``value`` as a key-safe namespace token, or raise naming it."""
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, str):
        raise TypeError(f"{name} must be a string, got "
                        f"{type(value).__name__} {value!r}")
    if not value:
        raise ValueError(f"{name} must not be empty")
    if len(value) > _MAX_NAMESPACE:
        raise ValueError(f"{name} must be at most {_MAX_NAMESPACE} characters, "
                         f"got {len(value)}")
    bad = sorted(set(value) - _NAMESPACE_CHARS)
    if bad:
        raise ValueError(f"{name} may contain only letters, digits, underscore "
                         f"and hyphen, so that it round-trips through an ASCII "
                         f"key; got {value!r} containing {bad!r}")
    return value


def _is_canonical_integer(text: str) -> bool:
    """Is ``text`` the one spelling this package emits for an integer?

    Decimal digits only: no sign, no leading zero on a multi-digit number, no
    underscore separator, no surrounding whitespace.  ``str(int)`` produces
    exactly this and nothing else, which is what makes a key round trip.
    """
    if not text or not text.isascii() or not text.isdigit():
        return False
    return text == "0" or not text.startswith("0")


def _require_index(value, name: str) -> int:
    """Return ``value`` as an identifier in ``[0, _MAX_IDENTIFIER]``, or raise.

    The scalar half of the one grammar.  :func:`_require_index_array` is the
    vector half and enforces the same bound on the same terms, because a
    boundary that accepts an identifier one way and refuses it the other is a
    boundary that has two contracts.
    """
    out = require_positive_int(value, name, minimum=0)
    if out > _MAX_IDENTIFIER:
        raise ValueError(f"{name} must be at most {_MAX_IDENTIFIER} (the "
                         f"largest exactly representable unsigned 64-bit "
                         f"identifier), got {out}")
    return out


def _require_index_array(values, name: str):
    """Return ``values`` as a list of exact non-negative integer identifiers.

    Coercion is refused rather than performed.  ``[0.0, 1.0]`` and
    ``[False, True]`` both look like clock identifiers to ``np.asarray`` and
    neither is one; a float identifier that silently truncates would address a
    different clock's noise stream, which is a wrong answer that looks like a
    right one.

    **No lossy cast, in either direction.**  An earlier version checked
    ``any(values < 0)`` and then called ``astype(np.int64)``.  For an unsigned
    array that order is the bug: ``np.uint64(2**63)`` is not negative, so it
    passed, and the cast then wrapped it to ``-2**63``, which was serialized
    into the key and generated a real stream under an identifier the scalar
    boundary refuses.  The range check now happens on the *unconverted* values,
    against the shared bound, and what comes back is a list of Python integers,
    which are exact for every ``int64`` and every ``uint64``.
    """
    if isinstance(values, (str, bytes, bytearray)):
        raise TypeError(f"{name} must be integer identifiers, got "
                        f"{type(values).__name__} {values!r}")
    try:
        raw = np.asarray(values)
    except OverflowError:
        # A Python list holding an integer too large for any NumPy integer
        # dtype.  It is out of the declared grammar, so it fails the same way a
        # representable but out-of-range identifier does.
        raise ValueError(f"{name} must lie in [0, {_MAX_IDENTIFIER}]; one "
                         f"entry is too large for any integer type") from None
    if raw.size == 0:
        # Before the dtype rules: an empty list has float dtype and reporting
        # that would answer a question the caller did not ask.  Naming no clock
        # at all is a wrong value, not a wrong kind of thing.
        raise ValueError(f"{name} must name at least one clock")
    if raw.dtype.kind == "b":
        raise TypeError(f"{name} must be integer identifiers, got booleans")
    if raw.dtype.kind == "O":
        # `np.asarray([1, 2**64])` gives an object array rather than raising.
        # Its entries are Python integers, so they can be range-checked below;
        # anything else in an object array is the wrong kind of thing.
        for entry in np.atleast_1d(raw).ravel():
            if isinstance(entry, bool) or not isinstance(entry, int):
                raise TypeError(f"{name} must be integer identifiers, got "
                                f"{type(entry).__name__} {entry!r}")
    elif raw.dtype.kind not in "iu":
        raise TypeError(f"{name} must be integer identifiers, got dtype "
                        f"{raw.dtype}; floats and strings are refused rather "
                        f"than coerced, because a truncated identifier "
                        f"addresses another clock's stream")
    out = np.atleast_1d(raw)
    if out.ndim != 1:
        raise ValueError(f"{name} must be one-dimensional, got shape "
                         f"{raw.shape}")
    identifiers = [int(entry) for entry in out]
    for entry in identifiers:
        if entry < 0:
            raise ValueError(f"{name} must be non-negative everywhere, got "
                             f"{entry}")
        if entry > _MAX_IDENTIFIER:
            raise ValueError(f"{name} must be at most {_MAX_IDENTIFIER} (the "
                             f"largest exactly representable unsigned 64-bit "
                             f"identifier), got {entry}")
    return identifiers


def _require_node_path(value, name: str) -> str:
    """Return ``value`` as a left/right ancestry path, or raise naming it."""
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, str):
        raise TypeError(f"{name} must be a string of 'L' and 'R', got "
                        f"{type(value).__name__} {value!r}")
    bad = sorted(set(value) - _NODE_CHARS)
    if bad:
        raise ValueError(f"{name} records left/right ancestry and may contain "
                         f"only 'L' and 'R'; got {value!r} containing {bad!r}")
    if len(value) > _MAX_NODE_DEPTH:
        raise ValueError(f"{name} is {len(value)} levels deep; more than "
                         f"{_MAX_NODE_DEPTH} nested crossings inside one "
                         f"finest step is a configuration error, not a tree")
    return value


def _require_same_clock(clock: int, crossing) -> None:
    """Refuse a crossing identity that names a different clock.

    The clock appears twice in a split key — once as the stream it belongs to
    and once inside the canonical crossing identity — and the plan specifies
    both.  Two clocks' tongues open at different times, so a crossing carried
    into the wrong clock's tree is a configuration error rather than a
    redundancy to tolerate: it would key one clock's split to another clock's
    event, quietly and reproducibly.
    """
    if crossing.clock != clock:
        raise ValueError(f"the crossing names clock {crossing.clock} and it is "
                         f"being applied to clock {clock}; a canonical "
                         f"crossing identity belongs to exactly one clock")


def _require_leaf_node(value, name: str) -> str:
    """Return ``value`` as a finished leaf's node label, or raise naming it.

    A leaf carries either the written-out root or the left/right ancestry that
    produced it; the two spellings are the two things a key can hold, and a
    third would name a node no split normal belongs to.
    """
    if value == _ROOT_NODE:
        return value
    path = _require_node_path(value, name)
    if not path:
        raise ValueError(f"{name} must be {_ROOT_NODE!r} for an unsplit step "
                         f"rather than an empty path, so that a leaf's node "
                         f"is written the same way it appears in a key")
    return path


def _require_digest(value, name: str) -> str:
    """Return ``value`` as a lowercase hexadecimal digest, or raise."""
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, str):
        raise TypeError(f"{name} must be a hexadecimal string, got "
                        f"{type(value).__name__} {value!r}")
    if not _MIN_DIGEST <= len(value) <= _MAX_DIGEST:
        raise ValueError(f"{name} must be between {_MIN_DIGEST} and "
                         f"{_MAX_DIGEST} hexadecimal characters, got "
                         f"{len(value)}")
    bad = sorted(set(value) - _DIGEST_CHARS)
    if bad:
        raise ValueError(f"{name} must be lowercase hexadecimal; got {value!r} "
                         f"containing {bad!r}")
    return value


def _require_label(value, name: str) -> str:
    """Return ``value`` as one of :data:`EDGE_LABELS`, or raise naming it."""
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, str):
        raise TypeError(f"{name} must be a string, got "
                        f"{type(value).__name__} {value!r}")
    if value not in EDGE_LABELS:
        raise ValueError(f"{name} must be one of {EDGE_LABELS!r}, got "
                         f"{value!r}; a rising entry and a falling exit are "
                         f"different physical events and cannot share a key")
    return value


def _require_unit_fraction(value, name: str) -> float:
    """Return ``value`` as a float in the closed interval [0, 1], or raise."""
    out = require_finite(value, name)
    if not 0.0 <= out <= 1.0:
        raise ValueError(f"{name} is the position of a crossing inside its "
                         f"parent interval and must lie in [0, 1], got "
                         f"{out!r}")
    return out


# ---------------------------------------------------------------------------
# Key -> value
# ---------------------------------------------------------------------------
def _key_words(key: str):
    """Return the four 64-bit words a canonical key hashes to.

    ASCII is required rather than encoded around: every field that reaches a
    key is a validated token, an integer, or a hexadecimal float, so a
    non-ASCII key means something upstream accepted a value it should not have.
    """
    digest = hashlib.sha256(key.encode("ascii")).digest()
    return [int.from_bytes(digest[i:i + 8], "little") for i in range(0, 32, 8)]


def _normals_from_key(key: str, count: int) -> np.ndarray:
    """Return ``count`` standard normals belonging to ``key``.

    The single place a key becomes a number.  Everything else in this module
    builds key strings; this turns one into values, and it is a pure function
    of the string, which is what makes the whole scheme reproducible and
    order-independent.
    """
    seed = np.random.SeedSequence(_key_words(key))
    generator = np.random.Generator(np.random.Philox(seed))
    return generator.standard_normal(count)


#: ``sqrt(2)``, once, so the factored scale below has one constant rather than
#: a repeated call.
_ROOT_TWO = math.sqrt(2.0)

#: The accepted Brownian-scale domain, ``[_MIN_NORMAL_SCALE, _MAX_SCALE]``.
#:
#: The lower end is the smallest *normal* double: a subnormal scale has fewer
#: than 53 significant bits and therefore cannot carry a Gaussian.  See
#: :func:`_require_representable_scale`.
#:
#: The upper end is ``sqrt(DBL_MAX)``, **not** ``DBL_MAX``, and the difference
#: is the difference between a finite kick and an infinite one.  A kick is
#: ``scale * z`` for a keyed standard normal ``z``; admitting scales up to
#: ``DBL_MAX`` admits ``D = h = DBL_MAX/2``, whose scale is ``1.27e308``, and
#: any ``|z| > 1.41`` then overflows — about a sixteenth of all draws, silently,
#: as ``inf``.  Capping at ``sqrt(DBL_MAX)`` makes the product finite even for
#: ``|z|`` as large as the scale itself, which is ``1.34e154``; real standard
#: normals do not exceed about 9, so the margin is 153 orders of magnitude.
#:
#: Taken from ``numpy.finfo`` rather than written as literals, so the bounds are
#: the platform's actual ones — and from numpy rather than ``sys.float_info``
#: because ``sys`` is not on the raw-module import allowlist and widening that
#: allowlist to spell a constant would trade an isolation contract for a
#: convenience.
_MIN_NORMAL_SCALE = float(np.finfo(np.float64).smallest_normal)
_MAX_SCALE = math.sqrt(float(np.finfo(np.float64).max))


def _exact_significand(value: float):
    """Return ``(m, e)`` with ``value == m * 2**e`` exactly, ``m`` an integer.

    ``math.frexp`` gives a fraction in ``[0.5, 1)`` and an exponent; scaling
    that fraction by ``2**53`` is exact for every binary64, subnormals
    included, because a double has at most 53 significant bits.  Nothing here
    rounds, so the pair is the number rather than an approximation of it.
    """
    fraction, exponent = math.frexp(value)
    return int(math.ldexp(fraction, 53)), exponent - 53


def _sqrt_below(factors, bound: float) -> bool:
    """Is ``sqrt(product of factors)`` strictly less than ``bound``?  Exactly.

    Both sides are non-negative, so the question is an integer one.  Each float
    is an exact ``m * 2**e``; multiplying the significands and adding the
    exponents gives the product exactly, as ``S * 2**p``, with **the unstable
    binary64 product never formed**.  Making ``p`` even lets the square root
    halve it, and the residual comparison ``sqrt(S) < B`` for non-negative
    integers is settled by ``math.isqrt``: ``floor(sqrt(S)) < B`` if and only
    if ``sqrt(S) < B``, because ``B`` is an integer.

    Note what is *not* here: no value is multiplied by itself.  Squaring the
    bound would be the obvious way to remove the square root, and it is also
    the one thing a raw module may not write — the isolation scan forbids it
    so that a squared coupling cannot be formed in the event path, and a
    constant is not worth an exception to that rule.  ``isqrt`` answers the
    same question from the other side.

    This exists because classifying by the *computed* scale is a different
    question.  At ``diffusion = ulp(0)`` and ``duration = nextafter(2**-971,
    0)`` the mathematical target is below the smallest normal double and the
    factored computation rounds *up* to exactly it, so a policy reading its own
    arithmetic would admit a target the mathematics rejects.
    """
    significand, exponent = 1, 0
    for factor in factors:
        factor_significand, factor_exponent = _exact_significand(factor)
        significand *= factor_significand
        exponent += factor_exponent
    if exponent % 2:
        significand <<= 1
        exponent -= 1
    bound_significand, bound_exponent = _exact_significand(bound)
    shift = exponent // 2 - bound_exponent
    if shift >= 0:
        significand <<= 2 * shift
    else:
        bound_significand <<= -shift
    return math.isqrt(significand) < bound_significand


def _scale_verdict(factors) -> str:
    """Classify ``sqrt(product of factors)`` against the accepted domain.

    ``"normal"``, ``"subnormal"`` or ``"overflow"``, decided by exact integer
    comparison against both ends.

    The accepted interval for the exact target is
    ``[MIN_NORMAL, sqrt(DBL_MAX))``, and the two ends are treated
    **asymmetrically, on purpose**.  The lower bound is closed: a target
    exactly equal to the smallest normal double *is* that double, and admitting
    it costs nothing.  The upper bound is open: a target exactly at
    ``sqrt(DBL_MAX)`` is refused along with everything above it.  That is deliberately conservative — it is the single point where the
    mathematical product ``2 D h`` equals ``DBL_MAX`` exactly, so refusing it
    gives up no physically distinguishable configuration and buys a strict
    inequality that is cheap to establish without squaring anything.
    """
    if _sqrt_below(factors, _MIN_NORMAL_SCALE):
        return "subnormal"
    if not _sqrt_below(factors, _MAX_SCALE):
        return "overflow"
    return "normal"


def _brownian_scale(diffusion: float, duration: float) -> float:
    """Return ``sqrt(2 * diffusion * duration)`` without forming the product.

    The product is what breaks.  ``2 * D * h`` underflows to exactly zero long
    before its square root does: at ``D = 0.2`` and ``h = math.ulp(0)`` the
    product is ``0.0`` while the scale ``1.4e-162`` is an ordinary double, and
    at ``D = h = 1e-200`` the product is again ``0.0`` while the scale
    ``1.4e-200`` is not remotely small.  A zero scale is not a small scale: it
    is the deterministic limit, so the old spelling silently turned a physical
    bath into no bath at all.  In the other direction the product overflows at
    large ``D * h`` where the scale itself is perfectly finite.

    Factoring first fixes both, because the square root of each operand is in
    range whenever the operand is: ``sqrt(x)`` of the smallest subnormal is
    ``2.2e-162`` and of the largest finite double is ``1.3e154``.  So

        sqrt(2 * D * h) = sqrt(2) * sqrt(D) * sqrt(h)

    is computed as written, and it can only overflow when the true scale is
    genuinely unrepresentable.  The cost is one extra rounding — three instead
    of two — so the result can differ from the product spelling in the last
    place.  That is a change to the key-to-kick mapping and is why the schema
    version moved.

    Returns ``0.0`` only when an operand is exactly zero.  A *nonzero*
    mathematical scale that underflows to zero is not returned as zero: the
    caller refuses the configuration instead, because "too small to represent"
    and "exactly the deterministic limit" are different physical statements.
    """
    if diffusion == 0.0 or duration == 0.0:
        return 0.0
    return _ROOT_TWO * math.sqrt(diffusion) * math.sqrt(duration)


def _require_representable_scale(scale: float, what: str, **context) -> float:
    """Return ``scale`` if it lies in the accepted domain, else raise.

    **This is the *computed-scale* half of the policy**, and it is closed at
    both ends: the binary64 value that will multiply the keyed normal must lie
    in ``[2.2250738585072014e-308, 1.3407807929942596e+154]``.  Plus the exact
    ``0.0`` of the deterministic limit, which reaches here only when an operand
    is exactly zero and never as the result of rounding.

    The *exact mathematical target* is held to a different interval —
    ``[MIN_NORMAL, sqrt(DBL_MAX))``, **open at the top** — by
    :func:`_scale_verdict`, and :func:`_admit_scale` requires both.  The two
    are not the same question and are deliberately not the same interval.

    Subnormals are excluded, and that exclusion is the physics rather than
    fastidiousness.  A subnormal has fewer than 53 significant bits — at the
    bottom of the range it has *one* — so multiplying a standard normal by one
    does not produce a Gaussian, it produces a handful of quantized magnitudes.
    At ``diffusion = duration = math.ulp(0)`` the required scale is
    ``6.987e-324`` and the nearest double is ``4.94e-324``: **29% low**, a
    variance halved, and a kick that takes exactly two magnitudes over eight
    draws — ``0.0`` and one subnormal.  Worse, the conditional split's bridge
    term ``scale * sqrt(a) * sqrt(1-a)`` rounds to **exactly zero** at
    ``a = 0.25``, so the split silently degenerates into the proportional
    subdivision the plan's mutation contract exists to reject.  None of that is
    a small numerical error; it is a different law.

    A subnormal *input* is fine when the resulting scale is normal:
    ``diffusion = math.ulp(0)`` with ``duration = 1`` gives ``3.14e-162``,
    which agrees with a 60-digit reference to the last bit.  The policy is
    about the output scale, not the operands.

    Called before any key is derived, everywhere a scale is formed —
    constructor, root window, tree node and split — so a configuration outside
    the domain is refused rather than quietly generating from a collapsed,
    quantized or infinite one.
    """
    if scale == 0.0:
        return scale
    if not math.isfinite(scale) or not _MIN_NORMAL_SCALE <= scale <= _MAX_SCALE:
        detail = ", ".join(f"{name}={value!r}" for name, value in context.items())
        raise ValueError(
            f"the {what} {scale!r} is outside the accepted domain "
            f"[{_MIN_NORMAL_SCALE!r}, {_MAX_SCALE!r}] at {detail}; a subnormal "
            f"scale has too few significant bits to carry a Gaussian — it "
            f"quantizes the kick to a few magnitudes and can round the "
            f"conditional split's bridge term to zero, turning it into "
            f"proportional subdivision — and a scale above sqrt(DBL_MAX) "
            f"overflows to infinity when multiplied by an ordinary keyed "
            f"normal")
    return scale


def _admit_scale(factors, scale: float, what: str, **context) -> float:
    """Return ``scale`` if both it and its exact target are admissible.

    **Two tests, and passing one is not enough.**  The *exact* target — the
    mathematical ``sqrt(product of factors)``, classified by integer
    arithmetic — must lie in the accepted domain, and so must the *computed*
    scale that will actually multiply the keyed normal.  Either alone can be
    wrong about the other: at ``diffusion = ulp(0)`` and
    ``duration = nextafter(2**-971, 0)`` the target is below the smallest
    normal double while the computed scale rounds up to exactly it, so reading
    the computation admits a scale the mathematics rejects.

    Requiring both is deliberately conservative at the boundaries.  A target
    inside the domain whose computed scale rounds outside it is refused, even
    though the mathematics allows it, because the computed scale is what the
    run would use; and a target outside whose computation rounds inside is
    refused, because the mathematics is what the run means.

    The two intervals differ at the top and the difference is not cosmetic.
    The exact target must lie in ``[MIN_NORMAL, sqrt(DBL_MAX))`` — lower end
    closed, **upper end open** — while the computed scale must lie in the
    closed representable interval ``[MIN_NORMAL, MAX_SCALE]``.  A target
    exactly at ``sqrt(DBL_MAX)`` is therefore refused, and a target exactly at
    ``MIN_NORMAL`` is admitted.
    """
    if scale == 0.0:
        return scale
    verdict = _scale_verdict(factors)
    if verdict != "normal":
        detail = ", ".join(f"{name}={value!r}"
                           for name, value in context.items())
        raise ValueError(
            f"the exact mathematical {what} is {verdict} at {detail}: it lies "
            f"outside the accepted domain [{_MIN_NORMAL_SCALE!r}, "
            f"{_MAX_SCALE!r}], whatever the factored computation "
            f"({scale!r}) happens to round to")
    return _require_representable_scale(scale, what, **context)


def _scale_for(diffusion: float, duration: float, what: str,
               **context) -> float:
    """Return the accepted Brownian scale for one interval, or raise.

    The one place a diffusion and a duration become a scale: computed by the
    stable factorization, classified exactly, and admitted only if both agree
    it belongs to the accepted domain.  Constructor, root windows, tree nodes
    and splits all come through here, so "accepted here" and "accepted there"
    cannot drift apart.
    """
    if diffusion == 0.0 or duration == 0.0:
        return 0.0
    return _admit_scale((2.0, diffusion, duration),
                        _brownian_scale(diffusion, duration), what,
                        phase_diffusion=diffusion, duration=duration,
                        **context)


def _sum_left_to_right(grouped: np.ndarray) -> np.ndarray:
    """Sum the trailing axis strictly in order, vectorized over the others.

    ``np.sum`` uses pairwise summation, which is a different — and for more
    than eight terms a *differently rounded* — order from walking the group in
    time.  The irregular coarse path accumulates leaf by leaf because that is
    the order a fine consumer would accumulate them in, and this is the same
    order for the uniform path, so the two are bit-identical on identical leaf
    sequences rather than merely close.
    """
    total = np.array(grouped[..., 0], dtype=float)
    for position in range(1, grouped.shape[-1]):
        total += grouped[..., position]
    return total


def conservation_tolerance(*values) -> float:
    """Return the scale-aware conservation bound for the given magnitudes.

    :data:`CONSERVATION_ULPS` units in the last place at the largest magnitude
    involved.  Exposed rather than hidden because "to machine precision" is a
    claim a reader should be able to evaluate: this function is the claim.

    Zero in, zero out — which is what makes the zero-duration and
    zero-diffusion cases *bit-exact* assertions rather than tolerant ones.

    The unit in the last place comes from :func:`math.ulp` rather than from
    ``np.spacing``.  At the top of the finite range ``np.spacing`` computes the
    *forward* gap and overflows: ``np.spacing(sys.float_info.max)`` is ``inf``
    (with a RuntimeWarning), so the returned bound was infinite and every
    conservation error passed it — a tolerance that cannot fail is not a
    tolerance.  ``math.ulp`` is finite for every finite argument, and the result
    is required to be finite before it is returned, so this function can never
    hand back a vacuous bound.
    """
    scale = 0.0
    for index, value in enumerate(values):
        scale = max(scale, abs(require_finite(value, f"values[{index}]")))
    if scale == 0.0:
        return 0.0
    tolerance = CONSERVATION_ULPS * math.ulp(scale)
    if not math.isfinite(tolerance):
        raise ValueError(f"the conservation bound at scale {scale!r} is not "
                         f"finite, so it would accept any error at all")
    return tolerance


# ---------------------------------------------------------------------------
# The frozen mesh, and the identity of a crossing
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class FinestMesh:
    """The finest uniform time grid, which every key is relative to.

    ``origin`` is the time of index zero and ``step`` the finest timestep.  The
    elementary physical mesh of the experiment is this grid *union* the exact
    deterministic tongue crossings; the crossings subdivide leaves of this one
    rather than living on a mesh of their own.

    The mesh is part of the key.  That is deliberate and it is the reason a
    refinement ladder means anything: coarser levels are sums of these leaves,
    so they share the mesh and therefore share the realization, while a run
    with a *different* finest step is a different dataset and gets different
    keys rather than silently reusing values whose variance no longer matches
    their duration.

    ``identity`` writes both numbers as hexadecimal floats, which round-trip
    exactly.  Decimal formatting would let two distinguishable meshes serialize
    to one key.
    """

    origin: float
    step: float

    def __post_init__(self):
        origin = require_finite(self.origin, "origin")
        step = require_finite(self.step, "step")
        if not step > 0.0:
            raise ValueError(f"step must be strictly positive, got {step!r}")
        object.__setattr__(self, "origin", origin)
        object.__setattr__(self, "step", step)

    @property
    def identity(self) -> str:
        """The mesh as it appears inside a key, exactly and reversibly."""
        return f"{self.origin.hex()}:{self.step.hex()}"

    def time_at(self, index) -> float:
        """Return the time of finest-grid index ``index``.

        ``origin + index * step``, evaluated once and the same way everywhere.
        This is the *only* expression in the package that turns an integer step
        index into a time; see :meth:`step_bounds` for why that matters.
        """
        return self.origin + _require_index(index, "index") * self.step

    def step_bounds(self, index) -> Tuple[float, float]:
        """Return the canonical half-open interval of finest step ``index``.

        ``(time_at(index), time_at(index + 1))`` — the *later* boundary is the
        next step's own start, never ``start + step`` computed separately.

        This is the whole of the mesh-boundary rule, and it exists because the
        two spellings are not the same number.  On ``FinestMesh(1.0, 0.1)``,
        ``time_at(1) + step`` is ``1.2000000000000002`` while ``time_at(2)`` is
        ``1.2``: step 1 claimed to end after step 2 began.  A crossing at
        exactly ``1.2`` then satisfied the containment test of *both* steps and
        was split twice, producing a spurious child of width ``2.2e-16`` and a
        seventh leaf where six were due.  Over two hundred steps of that mesh
        there were 31 such overlaps and 36 gaps in the other direction, where a
        crossing would have belonged to no step at all.

        Taking both boundaries from ``time_at`` makes overlaps and gaps
        impossible rather than rare: adjacent steps share one float, exactly.
        The cost is that a step's realized width can differ from the declared
        ``step``, in the last place on an ordinary mesh and by a quarter of the
        nominal step on ``FinestMesh(1e15, 0.2)``.  **The realized width is the
        physical duration**: it scales the root kick, it sets every child's
        variance and split fraction, it is the ``TreeLeaf.duration``, and it is
        what the coarse aggregation sums.  The nominal ``step`` generates the
        grid times and appears in the mesh identity; it is not the noise law.
        Those were once two different answers — the geometry used the realized
        width while the variance used the nominal step — and the tree was
        internally inconsistent about how long its own intervals were.

        A mesh whose steps are not resolvable at this index — where the two
        boundaries collide because ``origin`` swamps ``step`` — is a
        configuration error rather than a zero-width step to route into.
        """
        index = _require_index(index, "index")
        start = self.time_at(index)
        end = self.time_at(index + 1)
        if not end > start:
            raise ValueError(
                f"finest step {index} of this mesh spans [{start!r}, {end!r}], "
                f"which has no width: origin={self.origin!r} is too large for "
                f"step={self.step!r} to be resolved there, so no crossing "
                f"could be routed into it")
        return start, end

    def step_containing(self, time) -> int:
        """Return the unique finest step whose canonical interval holds ``time``.

        The routing half of :meth:`step_bounds`, and the only way a time
        becomes a step index.  Half-open: a time exactly on a boundary belongs
        to the step that *starts* there, so an instant shared by two adjacent
        steps is routed once and only once.

        The index is estimated by division and then *verified* against
        :meth:`step_bounds`, because the division is inexact and can land one
        step out at either end.  Verification is what makes the answer
        canonical; the estimate only makes it quick.
        """
        moment = require_finite(time, "time")
        if moment < self.origin:
            raise ValueError(f"time {moment!r} lies before the mesh origin "
                             f"{self.origin!r}")
        estimate = int(math.floor((moment - self.origin) / self.step))
        index = max(estimate, 0)
        for _ in range(_ROUTING_STEPS):
            start, end = self.step_bounds(index)
            if moment < start:
                if index == 0:
                    raise ValueError(f"time {moment!r} lies before the mesh "
                                     f"origin {self.origin!r}")
                index -= 1
                continue
            if moment >= end:
                index += 1
                continue
            return index
        raise ValueError(f"time {moment!r} could not be routed to a finest "
                         f"step of this mesh within {_ROUTING_STEPS} "
                         f"corrections; the mesh is not resolvable there")


@dataclass(frozen=True)
class Crossing:
    """One exact deterministic tongue crossing, and its canonical identity.

    Four fields, and each is in the key for a reason:

    ``clock``
        Two clocks cross at different times and must not share a split normal.
    ``edge``
        ``"rising"`` or ``"falling"``.  Entry and exit are different events
        even when a symmetric pulse puts them at symmetric times.
    ``pulse_digest``
        The hash of the pulse manifest the crossing time was computed from.  A
        different pulse gives a different tongue geometry, and a key that
        omitted this would let two different experiments collide.
    ``time``
        The crossing time itself, carried in the key as a hexadecimal float so
        that the identity is exact rather than rounded.

    The crossing is supplied by the calibrated deterministic layer.  Nothing
    here computes one, and nothing here knows what a tongue is.
    """

    clock: int
    edge: str
    pulse_digest: str
    time: float

    def __post_init__(self):
        clock = _require_index(self.clock, "clock")
        edge = _require_label(self.edge, "edge")
        pulse_digest = _require_digest(self.pulse_digest, "pulse_digest")
        time = require_finite(self.time, "time")
        object.__setattr__(self, "clock", clock)
        object.__setattr__(self, "edge", edge)
        object.__setattr__(self, "pulse_digest", pulse_digest)
        object.__setattr__(self, "time", time)

    @property
    def identity(self) -> str:
        """The canonical crossing identity, as it appears inside a key."""
        return (f"{self.clock}:{self.edge}:{self.pulse_digest}:"
                f"{self.time.hex()}")


@dataclass(frozen=True)
class TreeLeaf:
    """One elementary interval of the physical mesh, after splitting.

    ``node_path`` is the left/right ancestry inside the finest step that
    produced it, ``"root"`` when the step was never split.  It is carried
    because it is half of the split normal's key, so a reader can check that a
    repeated split never reused a normal.

    A leaf whose ``duration`` is zero is a real member of the list: an exact
    crossing that lands on an interval endpoint produces one, its ``kick`` is
    exactly zero, and the zero-duration coordinate handoff the state machine
    performs there consumes no noise.

    ``step_index`` is the finest uniform step the leaf subdivides.  It is
    carried so that a leaf is self-describing: a coarse consumer aggregating
    the irregular elementary mesh needs to know which uniform step each leaf
    belongs to, and recovering that by comparing times against a grid would be
    a floating-point comparison standing in for an integer that was known
    exactly at construction.
    """

    start_time: float
    duration: float
    kick: float
    node_path: str
    step_index: int

    def __post_init__(self):
        start_time = require_finite(self.start_time, "start_time")
        duration = require_nonnegative(self.duration, "duration")
        kick = require_finite(self.kick, "kick")
        node_path = _require_leaf_node(self.node_path, "node_path")
        step_index = _require_index(self.step_index, "step_index")
        object.__setattr__(self, "start_time", start_time)
        object.__setattr__(self, "duration", duration)
        object.__setattr__(self, "kick", kick)
        object.__setattr__(self, "node_path", node_path)
        object.__setattr__(self, "step_index", step_index)


# ---------------------------------------------------------------------------
# The stream
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class PhaseNoiseStream:
    """The keyed physical white-noise source: keys in, phase kicks out.

    ``namespace``
        The dataset namespace.  Pilot and production must not share one, and
        neither must two coupling sweeps intended to be statistically
        independent.  It is *not* a model label: the full Adler run and the
        width-only control share a namespace on purpose, because sharing the
        realization is what makes them paired.
    ``mesh``
        The finest uniform grid.  See :class:`FinestMesh`.
    ``phase_diffusion``
        ``D_phase``.  A kick over duration ``h`` has variance ``2 * D * h``.
        Zero is permitted and is exact: it produces identically zero kicks and
        consumes no key at all, which is the deterministic limit the Euler
        comparison needs.

    What is *not* here is the point of the class.  There is no coupling, no
    mismatch, no amplitude, no lock tolerance, no dwell, no model label and no
    analytic rate, so none of them can enter the noise law — not through a key,
    not through a scale factor, not through the number of values drawn.

    **The declared scale must lie in the accepted domain.**  That domain is the
    normal positive binary64 values plus the exact zero of the deterministic
    limit, and it is enforced by the same helper the generator uses, on the
    same stable factorization — see :func:`_require_representable_scale`.

    It used to be enforced by requiring ``2 * D * h`` to be finite, which is a
    different and wrong question in both directions.  Left-associated, that
    product overflows at ``D = DBL_MAX`` before the small duration is ever
    applied, so ``DBL_MAX`` with a nominal step of ``0.25`` was refused even
    though its scale, ``9.48e153``, is an ordinary double; and a product that
    merely *underflows* passed, delivering a subnormal or zero scale that no
    longer carries a Gaussian at all.  Asking about the scale directly answers
    both.

    Accepting a scale bounds everything downstream, and the **upper end is
    ``sqrt(DBL_MAX)``, not ``DBL_MAX``**.  That distinction is the difference
    between a finite kick and an infinite one: a kick is ``scale * z``, so
    admitting scales up to ``DBL_MAX`` admitted ``D = h = DBL_MAX/2``, whose
    scale is ``1.27e308``, and any ``|z| > 1.41`` overflowed it to infinity.
    Capped at ``1.3407807929942596e154``, the product is finite even for a
    normal as large as the scale itself; real standard normals do not exceed
    about 8, so the margin is 153 orders of magnitude.

    This checks the *nominal* declaration once, cheaply, at construction.  It
    is not a substitute for the per-window check: what a run actually uses is
    the realized width of each requested step, and those are validated
    separately, every time, before any randomness.
    """

    namespace: str
    mesh: FinestMesh
    phase_diffusion: float

    def __post_init__(self):
        namespace = _require_token(self.namespace, "namespace")
        mesh = require_exact_type(self.mesh, FinestMesh, "mesh")
        phase_diffusion = require_nonnegative(self.phase_diffusion,
                                              "phase_diffusion")
        _scale_for(phase_diffusion, mesh.step, "declared nominal kick scale")
        object.__setattr__(self, "namespace", namespace)
        object.__setattr__(self, "mesh", mesh)
        object.__setattr__(self, "phase_diffusion", phase_diffusion)

    # -- keys ------------------------------------------------------------
    @property
    def stream_prefix(self) -> str:
        """Everything a key holds before the trial: schema, dataset, mesh."""
        return f"{KEY_SCHEMA}/{self.namespace}/mesh={self.mesh.identity}"

    def leaf_key(self, trial, clock, index) -> str:
        """Return the canonical address of **one** finest-grid kick.

        A function of (dataset, mesh, trial, clock, step), and — unlike the v1
        key it replaces — injective in the step: two different steps have two
        different addresses, including the 1024 steps inside one generated
        block.  :meth:`leaf_from_key` turns this string back into the value it
        names, which is the property a manifest or audit record needs and the
        v1 key did not have.

        Blocking is how the value is *computed*, and :meth:`_block_key` is that
        private seed.  Keeping the two apart is the whole fix: the optimization
        stays, and the public vocabulary stops depending on it.
        """
        trial = _require_index(trial, "trial")
        clock = _require_index(clock, "clock")
        index = _require_index(index, "index")
        return (f"{self.stream_prefix}/trial={trial}/clock={clock}/"
                f"step={index}")

    @property
    def ensemble_prefix(self) -> str:
        """Everything an *initial-state* key holds before the trial.

        Deliberately **not** :attr:`stream_prefix`, and the difference is one
        field: the mesh is absent.

        An initial phase is a property of the prepared ensemble, not of the
        integration that follows it, and the plan requires a refinement ladder
        to "reuse the same initial phases and physical clock grid".  Keying the
        preparation to the mesh would break exactly that: halving the timestep
        would hand the same clock a different starting phase, and the paired
        fine/coarse comparison would then be measuring two different ensembles
        as well as two different timesteps.

        The literal ``/ensemble/`` segment is what keeps this disjoint from
        every kick address.  A leaf key is
        ``<schema>/<namespace>/mesh=.../trial=...``; an initial-phase key is
        ``<schema>/<namespace>/ensemble/trial=...``.  No mesh identity begins
        with ``ensemble``, because a mesh identity is two hexadecimal floats
        joined by a colon, so the two families cannot collide.
        """
        return f"{KEY_SCHEMA}/{self.namespace}/ensemble"

    def initial_phase_key(self, trial, clock) -> str:
        """Return the canonical address of one clock's prepared phase.

        A function of (dataset, trial, clock) and of nothing else.  Not the
        coupling, not the mismatch, not the diffusion, not the timestep, not
        the batch, not the position of the clock in whatever order a run walks
        its population — so the prepared ensemble a trial starts from is the
        same ensemble whatever the run does around it.
        """
        trial = _require_index(trial, "trial")
        clock = _require_index(clock, "clock")
        return f"{self.ensemble_prefix}/trial={trial}/clock={clock}/phase"

    def initial_phases(self, trial, clocks) -> np.ndarray:
        """Return the prepared uniform phases of ``clocks`` on one trial.

        Independently uniform on ``[-pi, pi)``, one draw per clock, each from
        its own :meth:`initial_phase_key`.  Half-open at ``+pi`` because that
        is the same half-open convention ``dynamics`` uses to resolve the
        half-turn tie at a coordinate handoff; a phase of exactly ``+pi`` and
        one of exactly ``-pi`` are the same angle, and drawing both would give
        the seam twice the weight of anywhere else.

        Drawn one key at a time rather than one key per trial, so that a run
        walking a subset of the population — a clock block, a single clock, a
        rerun of one trial for a diagnostic — prepares exactly the phases it
        would have prepared inside the whole population.  A per-trial key would
        have made the draw depend on how many clocks the caller asked for.
        """
        trial = _require_index(trial, "trial")
        identifiers = _require_index_array(clocks, "clocks")
        phases = np.empty(len(identifiers), dtype=float)
        for position, clock in enumerate(identifiers):
            key = self.initial_phase_key(trial, clock)
            seed = np.random.SeedSequence(_key_words(key))
            generator = np.random.Generator(np.random.Philox(seed))
            phases[position] = generator.uniform(-math.pi, math.pi)
        return phases

    def _block_key(self, trial: int, clock: int, block: int) -> str:
        """Return the private seed of one generated block of leaves.

        Not a public address and deliberately not reachable as one: it names
        1024 values at once, which is exactly why publishing it as "the key of
        one finest-grid kick" was wrong.
        """
        return (f"{self.stream_prefix}/trial={trial}/clock={clock}/"
                f"block={block}")

    def leaf_from_key(self, key) -> TreeLeaf:
        """Return the elementary leaf a canonical leaf address names.

        The reconstruction half of :meth:`leaf_key`: given the published
        address of a kick and nothing else, this returns that kick — the same
        object, field for field, that :meth:`elementary_leaves` hands a fine
        consumer for an unsplit step.  A key that identified 1024 values could
        not support this method, which is the cleanest statement of what was
        wrong with it.

        The address must belong to *this* stream.  A key from another dataset,
        another mesh or another schema version names a value this stream cannot
        produce, and guessing which one was meant is exactly the silent
        cross-run mixing the namespace exists to prevent.

        **Only the exact canonical spelling is accepted**, and the rule that
        makes that true is the last line rather than the parser: whatever comes
        out is spelled back through :meth:`leaf_key` and required to equal the
        input byte for byte.  ``int()`` is far more generous than a canonical
        grammar — it accepted ``step=007``, ``step=+7``, ``trial=01`` and
        ``step= 7`` as aliases of one address, and it read ``step=7_0`` as
        **seventy**, silently returning a different leaf than the one the
        caller named.  A round trip cannot be fooled by an alias, because an
        alias is by definition not what ``leaf_key`` emits.
        """
        if isinstance(key, (bool, np.bool_)) or not isinstance(key, str):
            raise TypeError(f"key must be a string, got "
                            f"{type(key).__name__} {key!r}")
        prefix = f"{self.stream_prefix}/trial="
        if not key.startswith(prefix):
            raise ValueError(f"key {key!r} does not belong to this stream, "
                             f"whose addresses begin {prefix!r}")
        fields = key[len(prefix):].split("/")
        if len(fields) != 3 or not fields[1].startswith("clock=") \
                or not fields[2].startswith("step="):
            raise ValueError(f"key {key!r} is not a canonical leaf address of "
                             f"the form <stream>/trial=<t>/clock=<c>/step=<j>")
        digits = (fields[0], fields[1][len("clock="):], fields[2][len("step="):])
        for text in digits:
            if not _is_canonical_integer(text):
                raise ValueError(
                    f"key {key!r} spells an integer as {text!r}; a canonical "
                    f"address uses decimal digits with no sign, no leading "
                    f"zero, no underscore and no whitespace")
        trial, clock, index = (int(text) for text in digits)
        if self.leaf_key(trial, clock, index) != key:
            raise ValueError(
                f"key {key!r} is not the canonical spelling of the address it "
                f"parses to, which is "
                f"{self.leaf_key(trial, clock, index)!r}")
        return self.elementary_leaves(trial, clock, index, 1, ())[0]

    def _split_key(self, trial, clock, index, node_path, crossing) -> str:
        """Return the immutable key of one split normal.

        **Private, because a caller who can spell a split key can forge one.**
        The node path is half of this key and the node's geometry is not in it
        at all, so a public builder let one ancestry be paired with two
        different intervals: an independent review split step 2 at
        ``node_path='R'`` over ``[t2, t2+h]`` and again over ``[t2, t2+0.8h]``,
        got the same key, and drew the same normal for two different
        conditional laws.  The tree is now the only thing that produces a node
        path, and it produces one per node it actually walked.

        Physical namespace, dataset, mesh, trial, clock, *finest parent step*,
        Brownian-tree node path, canonical crossing identity.  The finest
        parent step rather than the child interval, because the child interval
        has no index of its own; the node path is what distinguishes repeated
        splits of the same step, so no normal is ever reused.

        There is no model label in it, and there is no way to put one in.
        """
        trial = _require_index(trial, "trial")
        clock = _require_index(clock, "clock")
        index = _require_index(index, "index")
        node_path = _require_node_path(node_path, "node_path")
        require_exact_type(crossing, Crossing, "crossing")
        _require_same_clock(clock, crossing)
        node = node_path if node_path else _ROOT_NODE
        return (f"{self.stream_prefix}/trial={trial}/clock={clock}/"
                f"step={index}/node={node}/crossing={crossing.identity}/split")

    # -- leaves ----------------------------------------------------------
    def step_scale(self, index) -> float:
        """The standard deviation of the kick on finest step ``index``.

        ``sqrt(2 * D * h)`` where ``h`` is the step's **realized** width —
        ``step_bounds(index)[1] - step_bounds(index)[0]`` — and not the
        separately declared nominal ``step``.

        The two are the same number on an ordinary mesh and are not the same
        number in general.  On ``FinestMesh(1e15, 0.2)`` the first step's
        realized width is ``0.25``: the grid times are so far from the origin
        that consecutive ones are five ulp apart rather than the nominal four.
        Scaling by the nominal step there gave a kick of variance ``0.2`` to an
        interval the very same tree recorded as lasting ``0.25`` — the
        ``TreeLeaf.duration``, the split fractions and the coarse aggregation
        all used the realized width, and only the variance did not.  A tree
        cannot be in two minds about how long its own intervals are, and the
        realized width is the one the geometry already used.

        This is a per-step quantity now, so it is a method rather than the
        property it used to be: a mesh does not have *one* leaf scale.
        """
        start, end = self.mesh.step_bounds(index)
        return _scale_for(self.phase_diffusion, end - start,
                          "Brownian scale", step=index)

    def standard_normals(self, trial, clock, start_index, count) -> np.ndarray:
        """Return the keyed standard normals of one clock's finest steps.

        The unscaled layer, exposed because it is what makes the noise
        *shareable*: the same normals scaled by two different diffusion
        strengths are the same realization at two noise levels, and the same
        normals are what the full model and the width-only control both
        receive.
        """
        trial = _require_index(trial, "trial")
        clock = _require_index(clock, "clock")
        start_index = _require_index(start_index, "start_index")
        count = require_positive_int(count, "count", minimum=1)
        self._require_window(start_index, count)
        return self._normals(trial, clock, start_index, count)

    def _window_times(self, start: int, count: int) -> np.ndarray:
        """Return the ``count + 1`` canonical grid times bounding a window.

        Vectorized ``origin + index * step`` — the same expression
        :meth:`FinestMesh.time_at` evaluates, in the same order, so the two
        agree bit for bit.  The suite asserts that agreement rather than
        assuming it, because a second spelling of the mesh is exactly the
        defect the canonical-bounds rule exists to prevent.
        """
        indices = np.arange(start, start + count + 1, dtype=np.float64)
        return self.mesh.origin + indices * self.mesh.step

    def _require_window(self, start: int, count: int) -> None:
        """Refuse a window that is not physically resolvable, before any RNG.

        Every requested step must have a strictly positive realized width.  A
        collapsed step is not a small step: on ``FinestMesh(1e16, 0.1)``
        consecutive grid times are the *same double*, so the window describes
        no time at all — and the generator, which addresses steps by integer
        index and never looked at their widths, was perfectly happy to return
        plausible kicks for it.  A caller reading those numbers has a physical
        result for a physically empty interval.

        Checked eagerly and in full, so a window whose *fifth* step collapses
        fails at the request rather than after four steps of randomness have
        already been drawn and possibly recorded.  The scan runs in bounded
        chunks: a production request can name more steps than fit in memory,
        and validating one is not a reason to materialize it.
        """
        remaining, at = count, start
        while remaining > 0:
            span = min(remaining, _VALIDATION_CHUNK)
            times = self._window_times(at, span)
            if not bool(np.all(np.isfinite(times))):
                bad = int(np.argmax(~np.isfinite(times)))
                raise ValueError(
                    f"finest step {at + bad} of this mesh has a non-finite "
                    f"boundary {times[bad]!r}")
            widths = np.diff(times)
            if not bool(np.all(widths > 0.0)):
                bad = int(np.argmax(~(widths > 0.0)))
                raise ValueError(
                    f"finest step {at + bad} of this mesh spans "
                    f"[{times[bad]!r}, {times[bad + 1]!r}], which has no "
                    f"width: origin={self.mesh.origin!r} is too large for "
                    f"step={self.mesh.step!r} to be resolved there, so the "
                    f"requested window describes no physical time")
            self._require_window_scales(at, widths)
            remaining -= span
            at += span

    def _require_window_scales(self, at: int, widths: np.ndarray) -> None:
        """Refuse a chunk whose root scales leave the accepted domain.

        The scale is monotone increasing in the width at fixed diffusion, so
        the narrowest and widest steps of the chunk bracket every scale in it:
        classifying those two exactly settles all of them.  That matters
        because the exact classifier is integer arithmetic on big numbers and
        a production window can name a billion steps — two comparisons per
        chunk rather than a billion, and no weaker for it.
        """
        if self.phase_diffusion == 0.0:
            return
        narrowest = int(np.argmin(widths))
        widest = int(np.argmax(widths))
        for offset in {narrowest, widest}:
            _scale_for(self.phase_diffusion, float(widths[offset]),
                       "Brownian scale", step=at + offset)

    def _normals(self, trial: int, clock: int, start: int,
                 count: int) -> np.ndarray:
        """Assemble ``count`` normals from the frozen blocks holding them."""
        out = np.empty(count, dtype=float)
        written, index = 0, start
        while written < count:
            block, offset = divmod(index, LEAF_BLOCK)
            take = min(LEAF_BLOCK - offset, count - written)
            values = _normals_from_key(self._block_key(trial, clock, block),
                                       LEAF_BLOCK)
            out[written:written + take] = values[offset:offset + take]
            written += take
            index += take
        return out

    def leaf_kicks(self, trial, clocks, start_index, count) -> np.ndarray:
        """Return finest-grid phase kicks, shaped ``(len(clocks), count)``.

        Each row belongs to the clock identifier in the corresponding position
        of ``clocks``, so asking for a subset of the population — or for the
        survivors of a trial in which other clocks have already stopped — gives
        exactly the values those clocks would have received in the full run.

        Zero diffusion returns an array of exact zeros without deriving a
        single key.  That is not an optimization: "the deterministic limit
        consumes no randomness" is a property the verifier asserts by counting
        key derivations.
        """
        trial = _require_index(trial, "trial")
        identifiers = _require_index_array(clocks, "clocks")
        start_index = _require_index(start_index, "start_index")
        count = require_positive_int(count, "count", minimum=1)

        self._require_window(start_index, count)

        if self.phase_diffusion == 0.0:
            return np.zeros((len(identifiers), count), dtype=float)

        scales = self._window_scales(start_index, count)
        out = np.empty((len(identifiers), count), dtype=float)
        for row, clock in enumerate(identifiers):
            out[row] = self._normals(trial, clock, start_index, count)
        return out * scales

    def _window_scales(self, start: int, count: int) -> np.ndarray:
        """Return each step's realized Brownian scale across a window.

        The vectorized form of :meth:`step_scale`, factored the same way — the
        square roots are taken before anything is multiplied, so a scale that
        the product spelling would underflow to zero survives here.  Every
        entry must be a usable positive number; the window validation above has
        already established that the widths are positive, so a non-finite or
        zero scale at this point means the *scale* is unrepresentable rather
        than the interval, and that is refused before a key is derived.
        """
        widths = np.diff(self._window_times(start, count))
        scales = (_ROOT_TWO * math.sqrt(self.phase_diffusion)
                  * np.sqrt(widths))
        accepted = (np.isfinite(scales) & (scales >= _MIN_NORMAL_SCALE)
                    & (scales <= _MAX_SCALE))
        if not bool(np.all(accepted)):
            bad = int(np.argmax(~accepted))
            # Reported through the scalar policy so the message, the domain
            # and the exception type are the same ones every other path uses.
            # Reaching here means the window preflight and this computation
            # disagree, which the suite asserts cannot happen.
            _scale_for(self.phase_diffusion, float(widths[bad]),
                       "Brownian scale", step=start + bad)
            raise ValueError(f"the Brownian scale of finest step "
                             f"{start + bad} is outside the accepted domain")
        return scales

    def _preflight_schedule(self, start: int, count: int, buckets) -> None:
        """Validate every scale the tree will need, before it draws anything.

        The window check covers the *root* scales.  A split needs two more per
        crossing — the node's own scale and the conditional bridge — and those
        depend on the sub-interval the crossing carves out, which the window
        knows nothing about.  Without this, a schedule whose bridge scale is
        inadmissible was discovered *after* the parent kick had been drawn: one
        key consumed, and on the streaming path not until the first ``next()``,
        long after the call that should have refused it.

        The tree's geometry is a function of the crossing times alone — the
        kick values never affect where it branches — so the whole recursion can
        be walked here with no randomness at all.  Only nodes that actually
        split are checked, because only they draw.
        """
        for index, crossings in buckets.items():
            if not start <= index < start + count:
                continue
            # Depth first, and before anything recurses: the scale walk below
            # descends once per crossing, so a schedule with thousands of them
            # in one step would exhaust the interpreter stack on the way to
            # discovering it was too deep.
            self._require_depth(index, crossings)
            if self.phase_diffusion == 0.0:
                continue
            node_start, node_end = self.mesh.step_bounds(index)
            self._preflight_node(index, node_start, node_end, crossings)

    def _require_depth(self, index: int, crossings) -> None:
        """Refuse a step carrying more crossings than the tree is deep.

        The rule is not new — :meth:`_split_step` has always enforced it — but
        it used to fire *during* generation, so a schedule that was too deep at
        step 2 was discovered only after steps 0 and 1 had drawn their root
        kicks.  Two keys consumed for a request that was always going to be
        refused, and on the streaming path the refusal did not arrive until the
        first ``next()``, long after the call that should have made it.

        Enforced here it is settled with the scales, before any randomness, on
        every schedule-taking path.  Scoped to the requested window, like every
        other preflight rule: a crossing outside the window is not an error,
        because a caller streaming a window at a time holds the whole clock's
        schedule and is not asked to pre-partition it.
        """
        if len(crossings) > _MAX_NODE_DEPTH:
            raise ValueError(
                f"{len(crossings)} crossings inside finest step {index} "
                f"exceeds the {_MAX_NODE_DEPTH}-level tree depth; that many "
                f"exact crossings in one timestep is a configuration error, "
                f"not a deep tree")

    def _preflight_node(self, index: int, start: float, end: float,
                        crossings) -> None:
        """Walk one step's tree geometrically, validating each split's scales."""
        if not crossings:
            return
        crossing = crossings[0]
        duration = end - start
        fraction = (0.0 if duration == 0.0
                    else (crossing.time - start) / duration)
        if duration > 0.0 and 0.0 < fraction < 1.0:
            node_scale = _scale_for(self.phase_diffusion, duration,
                                    "node Brownian scale", step=index)
            _admit_scale(
                (2.0, self.phase_diffusion, duration, fraction,
                 1.0 - fraction),
                node_scale * math.sqrt(fraction) * math.sqrt(1.0 - fraction),
                "conditional split bridge scale", step=index,
                duration=duration, fraction=fraction,
                phase_diffusion=self.phase_diffusion)
        # Only the right child can carry further crossings: the earliest one
        # splits first, so the left child is finished.
        self._preflight_node(index, crossing.time, end, crossings[1:])

    def _coarse_kicks_uniform(self, trial, clocks, start_index, count,
                              stride) -> np.ndarray:
        """Return coarse increments as sums of the **unsplit root** kicks.

        **Private, and only the verifier calls it.**  It cannot be told about
        crossings, so on a window that has any it sums roots rather than the
        elementary leaves the fine path consumes — and root sums and leaf sums
        agree numerically, frequently bit for bit.  A public method with that
        property is a trap: nothing about the call site would look wrong, no
        numeric check would fire, and the two refinement levels would quietly
        stop reading one realization.  It survives as the control that
        :meth:`coarse_kicks` is measured against.

        ``stride`` finest leaves per coarse step, ``count`` finest leaves in
        total, and ``count`` must be a whole number of coarse steps: a partial
        trailing group would be a coarse increment of the wrong duration, and
        silently dropping or keeping it is exactly the kind of choice that
        makes two refinement levels disagree for a reason that is not physics.

        The summation order is fixed — strictly left to right in time within
        each group — so a coarse increment is reproducible bit for bit, and a
        fine consumer and a coarse consumer are demonstrably reading one
        realization rather than two.  It is the same order
        :meth:`coarse_kicks` accumulates in, so the two agree
        *bit for bit* when a step carries no crossing, and the suite asserts
        that rather than asserting they agree to within a tolerance.

        It knows nothing about crossings, so on a window that has any it sums
        the wrong things; :meth:`coarse_kicks` is the public path and takes the
        frozen schedule.
        """
        stride = require_positive_int(stride, "stride", minimum=1)
        count = require_positive_int(count, "count", minimum=1)
        if count % stride:
            raise ValueError(f"count must be a whole number of coarse steps: "
                             f"count={count} is not divisible by "
                             f"stride={stride}")
        fine = self.leaf_kicks(trial, clocks, start_index, count)
        grouped = fine.reshape(fine.shape[0], count // stride, stride)
        return _sum_left_to_right(grouped)

    def stream_leaf_blocks(self, trial, clocks, start_index, count,
                           window) -> Iterator[Tuple[int, np.ndarray]]:
        """Iterate ``(start index, kicks)`` in bounded windows.

        The whole reason this exists: a production sweep is trials by clocks by
        timesteps, and materializing that cube is infeasible.  Nothing here
        holds more than ``len(clocks) * window`` values at once, and — because
        the values belong to keys rather than to a running generator state —
        the windowing is invisible in the numbers.  A caller may use any window
        it likes, including one that divides no block boundary evenly.

        Validation happens now rather than on first iteration.  A generator
        function would defer every check to the first ``next()``, which is
        precisely where an invalid configuration stops being reported at the
        boundary that was supposed to report it.
        """
        trial = _require_index(trial, "trial")
        identifiers = _require_index_array(clocks, "clocks")
        start_index = _require_index(start_index, "start_index")
        count = require_positive_int(count, "count", minimum=1)
        window = require_positive_int(window, "window", minimum=1)
        self._require_window(start_index, count)
        return self._stream(trial, identifiers, start_index, count, window)

    def _stream(self, trial: int, identifiers: np.ndarray, start: int,
                count: int, window: int) -> Iterator[Tuple[int, np.ndarray]]:
        """The bounded iteration itself; arguments are already validated."""
        emitted = 0
        while emitted < count:
            take = min(window, count - emitted)
            index = start + emitted
            yield index, self.leaf_kicks(trial, identifiers, index, take)
            emitted += take

    # -- the tree --------------------------------------------------------
    def _split_node(self, trial, clock, index, node_path, node_start, node_end,
                    parent_kick, crossing) -> Tuple[float, float]:
        """Split a node's assigned kick at the crossing that lies inside it.

        **Private: only the tree walk calls this.**  It validates its geometry
        as strictly as it ever did — and that validation is still what stops an
        out-of-parent crossing reaching a key — but validation cannot make a
        *freely supplied* ancestry consistent, because nothing in the arguments
        says which node ``"R"`` is meant to be.  Ownership is the fix:
        :meth:`elementary_leaves` walks the tree, and every ``node_path`` that
        reaches this method was produced by that walk one recursion earlier.

        The fraction is
        *derived* from ``crossing.time`` and the node's exact interval; it is
        not something a caller supplies alongside a crossing that says
        something else.  The arithmetic primitive underneath is private for
        exactly that reason.

        What that closes: an earlier public primitive took a fraction and a
        crossing as independent arguments and never compared them, so both of
        these succeeded and consumed a split normal keyed to an identity
        inconsistent with the geometry it was used for —

        * a crossing at ``t = 99.0`` used to split the step spanning
          ``[0.04, 0.05]``; and
        * a crossing a quarter of the way through the interval, split at a
          supplied fraction of three quarters.

        Both are now refused, and refused **before** any key is derived: the
        order of the checks is part of the contract, because a rejected
        configuration that has already consumed randomness has already changed
        the run it was rejected from.

        ``node_start`` and ``node_end`` are the node's exact endpoints, which
        for a nested child are the crossing times that produced it rather than
        grid points.  They must lie inside the finest step ``index`` — a node
        of the tree is part of the step it subdivides — and the crossing must
        lie inside them.

        Returns ``(left, right)`` with ``left`` correctly distributed for the
        duration before the crossing, ``right`` for the remainder, zero
        covariance between them, and ``left + right`` equal to ``parent_kick``
        to within :func:`conservation_tolerance`.
        """
        trial = _require_index(trial, "trial")
        clock = _require_index(clock, "clock")
        index = _require_index(index, "index")
        node_path = _require_node_path(node_path, "node_path")
        node_start = require_finite(node_start, "node_start")
        node_end = require_finite(node_end, "node_end")
        parent_kick = require_finite(parent_kick, "parent_kick")
        require_exact_type(crossing, Crossing, "crossing")
        _require_same_clock(clock, crossing)

        if not node_end >= node_start:
            raise ValueError(f"the node interval [{node_start!r}, "
                             f"{node_end!r}] runs backwards")
        step_start, step_end = self.mesh.step_bounds(index)
        if not (step_start <= node_start and node_end <= step_end):
            raise ValueError(
                f"the node interval [{node_start!r}, {node_end!r}] is not "
                f"inside finest step {index}, which spans [{step_start!r}, "
                f"{step_end!r}]; a tree node is part of the step it subdivides")
        if not node_start <= crossing.time <= node_end:
            raise ValueError(
                f"the crossing at t={crossing.time!r} lies outside the node "
                f"interval [{node_start!r}, {node_end!r}] it is being used to "
                f"split")

        duration = node_end - node_start
        fraction = (0.0 if duration == 0.0
                    else (crossing.time - node_start) / duration)
        return self._split_at_fraction(parent_kick, duration, fraction, trial,
                                       clock, index, node_path, crossing)

    def _split_at_fraction(self, parent_kick, duration, fraction, trial, clock,
                           index, node_path, crossing) -> Tuple[float, float]:
        """Split one assigned kick conditionally at ``fraction`` of ``duration``.

        The arithmetic primitive, private because it cannot check that its
        fraction agrees with its crossing: it has no node interval to compare
        them against.  :meth:`_split_node` is the operation that supplies one
        and derives the fraction from it, and the tree walk is the only caller
        of either.

        Returns ``(left, right)`` with ``left`` correctly distributed for
        duration ``fraction * duration``, ``right`` for the remainder, zero
        covariance between them, and ``left + right`` equal to ``parent_kick``
        to within :func:`conservation_tolerance`.

        Four exact branches come first, and none of them draws a normal:

        * **zero diffusion** — every kick in the tree is exactly zero, so both
          children are ``0.0``.  A nonzero ``parent_kick`` at zero diffusion is
          not a rounding artifact but a caller who has mixed two trees, and it
          is refused.
        * **zero duration** — the same, for the same reason: an interval of no
          duration carries no variance, so a nonzero parent is a contradiction
          rather than a value to divide.
        * **fraction 0** — the crossing is at the start.  Left is the exact
          zero of a zero-duration interval; right is the whole parent.
        * **fraction 1** — the mirror image.

        Only a crossing strictly inside the interval reveals anything new, and
        only that case consumes a split normal.
        """
        parent_kick = require_finite(parent_kick, "parent_kick")
        duration = require_nonnegative(duration, "duration")
        fraction = _require_unit_fraction(fraction, "fraction")
        trial = _require_index(trial, "trial")
        clock = _require_index(clock, "clock")
        index = _require_index(index, "index")
        node_path = _require_node_path(node_path, "node_path")
        require_exact_type(crossing, Crossing, "crossing")
        _require_same_clock(clock, crossing)

        if self.phase_diffusion == 0.0:
            if parent_kick != 0.0:
                raise ValueError(
                    f"at zero diffusion every kick in the tree is exactly "
                    f"zero, so a parent kick of {parent_kick!r} cannot belong "
                    f"to it; the children of a zero-diffusion parent are "
                    f"exactly zero and no split normal is drawn")
            return 0.0, 0.0
        if duration == 0.0:
            if parent_kick != 0.0:
                raise ValueError(
                    f"an interval of zero duration carries zero variance, so "
                    f"a parent kick of {parent_kick!r} cannot belong to it")
            return 0.0, 0.0
        if fraction == 0.0:
            return 0.0, parent_kick
        if fraction == 1.0:
            return parent_kick, 0.0

        # Factored, and checked, before a key is derived.  The same underflow
        # that turned a root scale into the deterministic limit turns a split
        # scale into proportional subdivision, which is the one construction
        # the plan's mutation contract exists to reject.
        node_scale = _scale_for(self.phase_diffusion, duration,
                                "node Brownian scale")
        bridge = _admit_scale(
            (2.0, self.phase_diffusion, duration, fraction, 1.0 - fraction),
            node_scale * math.sqrt(fraction) * math.sqrt(1.0 - fraction),
            "conditional split bridge scale", duration=duration,
            fraction=fraction, phase_diffusion=self.phase_diffusion)
        key = self._split_key(trial, clock, index, node_path, crossing)
        split_normal = float(_normals_from_key(key, 1)[0])
        left = fraction * parent_kick + bridge * split_normal
        return left, parent_kick - left

    def _split_step(self, trial, clock, index, crossings) -> Tuple[TreeLeaf, ...]:
        """Return the elementary leaves of one finest step, in time order.

        **Private, because its crossing list is a *subset*.**  A public
        per-step entry point lets a caller choose which crossings a step sees,
        and two different subsets give the same node paths different
        geometries.  :meth:`elementary_leaves` takes the clock's whole frozen
        schedule and routes it canonically, so the subset a step receives is a
        function of the schedule and the mesh rather than of the caller.

        With no crossings that is the step itself.  With crossings it is the
        nested conditional split: the *earliest* crossing divides the step, the
        left child is then complete, and every remaining crossing recurses into
        the right child.  The node path records that ancestry, so a step split
        three times uses three different split normals and never reuses one.

        The crossings are the exact deterministic tongue entry and exit times
        supplied by the calibrated pulse.  They are validated hard rather than
        repaired:

        * a time outside the parent interval is a configuration error, not a
          crossing to clamp;
        * two crossings at the same instant leave the chronological order
          undefined, so they are refused rather than ordered arbitrarily; and
        * a repeated canonical identity would key two different splits to one
          normal, and is refused for that reason.

        The returned leaves sum to the step's assigned kick to within
        :func:`conservation_tolerance` of the largest magnitude in the tree,
        and their durations sum to the step's duration on the same terms.
        """
        trial = _require_index(trial, "trial")
        clock = _require_index(clock, "clock")
        index = _require_index(index, "index")

        if isinstance(crossings, Crossing):
            raise TypeError("crossings must be a sequence of Crossing "
                            "objects, not a single Crossing")
        try:
            ordered = list(crossings)
        except TypeError:
            raise TypeError(f"crossings must be a sequence of Crossing "
                            f"objects, got {type(crossings).__name__}") from None
        if len(ordered) > _MAX_NODE_DEPTH:
            raise ValueError(f"{len(ordered)} crossings inside one finest step "
                             f"exceeds the {_MAX_NODE_DEPTH}-level tree depth; "
                             f"that many exact crossings in one timestep is a "
                             f"configuration error, not a deep tree")
        for position, item in enumerate(ordered):
            require_exact_type(item, Crossing, f"crossings[{position}]")
            _require_same_clock(clock, item)

        start, end = self.mesh.step_bounds(index)
        for item in ordered:
            if not start <= item.time <= end:
                raise ValueError(
                    f"crossing at t={item.time!r} lies outside its parent "
                    f"interval [{start!r}, {end!r}]; an exact crossing belongs "
                    f"to the step that contains it")
        ordered.sort(key=lambda item: item.time)
        for first, second in zip(ordered, ordered[1:]):
            if first.time == second.time:
                raise ValueError(
                    f"two crossings share the instant t={first.time!r}, which "
                    f"leaves their chronological order undefined")
        identities = [item.identity for item in ordered]
        if len(set(identities)) != len(identities):
            raise ValueError("two crossings share one canonical identity, "
                             "which would key two different splits to the "
                             "same split normal")

        parent = float(self.leaf_kicks(trial, [clock], index, 1)[0, 0])
        leaves: list = []
        self._grow(leaves, trial, clock, index, "", start, end, parent, ordered)
        return tuple(leaves)

    def _grow(self, leaves, trial, clock, index, node_path, start, end,
              kick, crossings) -> None:
        """Recurse chronologically, appending finished leaves in time order."""
        if not crossings:
            leaves.append(TreeLeaf(start, end - start, kick,
                                   node_path if node_path else _ROOT_NODE,
                                   index))
            return
        crossing = crossings[0]
        # The node's exact endpoints go to the public geometry-validating
        # operation, which derives the fraction from them rather than being
        # handed one.  Nested children pass their own crossing times as
        # endpoints, so the containment check applies at every level of the
        # tree and not only at the root.
        left, right = self._split_node(trial, clock, index, node_path, start,
                                      end, kick, crossing)
        self._grow(leaves, trial, clock, index, node_path + "L", start,
                   crossing.time, left, [])
        self._grow(leaves, trial, clock, index, node_path + "R", crossing.time,
                   end, right, crossings[1:])

    # -- the irregular elementary mesh ------------------------------------
    def elementary_leaves(self, trial, clock, start_index, count,
                          crossings) -> Tuple[TreeLeaf, ...]:
        """Return the exact elementary leaves of a window of finest steps.

        **The elementary physical mesh is the finest uniform grid union the
        exact deterministic crossings**, and this is that mesh: every finest
        step in ``[start_index, start_index + count)`` in order, each split by
        whichever crossings fall inside it, concatenated chronologically.  It
        is the sequence a fine consumer walks, and — via
        :meth:`coarse_kicks` — the sequence a coarse consumer
        sums.  There is one tree, and this is how both sides read it.

        ``crossings`` is the frozen crossing set for this clock, in any order.
        Each is routed to the step containing it; a crossing outside the
        requested window is not an error, because a caller streaming a window
        at a time holds the whole clock's crossing list and should not have to
        pre-partition it.  A crossing that lands exactly on a shared step
        boundary belongs to the *later* step, so that no instant is split
        twice and the leaves of adjacent steps still abut exactly.
        """
        trial = _require_index(trial, "trial")
        clock = _require_index(clock, "clock")
        start_index = _require_index(start_index, "start_index")
        count = require_positive_int(count, "count", minimum=1)
        self._require_window(start_index, count)
        buckets = self._route(self._checked_crossings(clock, crossings))
        self._preflight_schedule(start_index, count, buckets)

        leaves: list = []
        for index in range(start_index, start_index + count):
            leaves.extend(self._split_step(
                trial, clock, index,
                self._crossings_in_step(buckets, index)))
        return tuple(leaves)

    def stream_elementary_leaves(self, trial, clock, start_index, count,
                                 crossings,
                                 window) -> Iterator[Tuple[TreeLeaf, ...]]:
        """Iterate the elementary leaves in bounded windows of finest steps.

        The same sequence :meth:`elementary_leaves` returns, delivered a window
        at a time so that a long run never holds more than one window of leaves.
        Concatenating the windows reproduces the whole sequence exactly, in any
        window size, because a leaf belongs to a key rather than to the walk
        that reached it.

        Validated eagerly, then iterated, for the reason
        :meth:`stream_leaf_blocks` gives.
        """
        trial = _require_index(trial, "trial")
        clock = _require_index(clock, "clock")
        start_index = _require_index(start_index, "start_index")
        count = require_positive_int(count, "count", minimum=1)
        window = require_positive_int(window, "window", minimum=1)
        self._require_window(start_index, count)
        buckets = self._route(self._checked_crossings(clock, crossings))
        self._preflight_schedule(start_index, count, buckets)
        return self._stream_leaves(trial, clock, start_index, count, buckets,
                                   window)

    def _stream_leaves(self, trial, clock, start, count, buckets, window):
        """The bounded leaf iteration; arguments are already validated."""
        emitted = 0
        while emitted < count:
            take = min(window, count - emitted)
            leaves: list = []
            for index in range(start + emitted, start + emitted + take):
                leaves.extend(self._split_step(
                    trial, clock, index,
                    self._crossings_in_step(buckets, index)))
            yield tuple(leaves)
            emitted += take

    def coarse_kicks(self, trial, clock, start_index, count,
                     stride, crossings) -> np.ndarray:
        """Return coarse increments as sums of the elementary leaves.

        **The one public refinement path.**  The frozen crossing schedule is a
        required argument, not an option: an empty schedule is the uniform
        case, and there is no signature through which a crossing-bearing window
        can be aggregated from unsplit roots.

        The general refinement path, and the one the plan actually specifies:
        a coarse increment spanning ``stride`` finest steps is the sum of every
        elementary leaf inside it — including the pieces an exact tongue
        crossing carved out of a step — accumulated strictly in time order.

        The private root-summing control cannot do this: it has no way to be
        told about crossings, so it sums unsplit roots.  Root sums and leaf
        sums agree to within rounding, which is precisely what makes the
        distinction easy to lose: a later refinement study could sum roots while
        the fine dynamics consumed crossing leaves, and the two levels would no
        longer be reading one realization even though the numbers looked right.
        With no crossings in a group this method reproduces that control **bit
        for bit**, because both accumulate the same leaves in the same order.

        Bounded: one coarse group of leaves is held at a time, never the window
        and never the cube.  Per clock rather than per population, because
        crossing times are a property of one clock's tongue.
        """
        trial = _require_index(trial, "trial")
        clock = _require_index(clock, "clock")
        start_index = _require_index(start_index, "start_index")
        count = require_positive_int(count, "count", minimum=1)
        stride = require_positive_int(stride, "stride", minimum=1)
        if count % stride:
            raise ValueError(f"count must be a whole number of coarse steps: "
                             f"count={count} is not divisible by "
                             f"stride={stride}")
        self._require_window(start_index, count)
        buckets = self._route(self._checked_crossings(clock, crossings))
        self._preflight_schedule(start_index, count, buckets)

        coarse = np.empty(count // stride, dtype=float)
        for group in range(count // stride):
            first = start_index + group * stride
            total = 0.0
            started = False
            for leaves in self._stream_leaves(trial, clock, first, stride,
                                              buckets, 1):
                for leaf in leaves:
                    # Seeded from the first leaf rather than from 0.0, so the
                    # accumulation is the same sequence of additions a fine
                    # consumer performs and the root-summing control can be
                    # matched bit for bit.
                    total = leaf.kick if not started else total + leaf.kick
                    started = True
            coarse[group] = total
        return coarse

    def _checked_crossings(self, clock: int, crossings):
        """Validate a whole clock's frozen crossing set; return it sorted."""
        if isinstance(crossings, Crossing):
            raise TypeError("crossings must be a sequence of Crossing "
                            "objects, not a single Crossing")
        try:
            ordered = list(crossings)
        except TypeError:
            raise TypeError(f"crossings must be a sequence of Crossing "
                            f"objects, got "
                            f"{type(crossings).__name__}") from None
        for position, item in enumerate(ordered):
            require_exact_type(item, Crossing, f"crossings[{position}]")
            _require_same_clock(clock, item)
        ordered.sort(key=lambda item: item.time)
        for first, second in zip(ordered, ordered[1:]):
            if first.time == second.time:
                raise ValueError(
                    f"two crossings share the instant t={first.time!r}, which "
                    f"leaves their chronological order undefined")
        identities = [item.identity for item in ordered]
        if len(set(identities)) != len(identities):
            raise ValueError("two crossings share one canonical identity, "
                             "which would key two different splits to the "
                             "same split normal")
        return ordered

    def _crossings_in_step(self, buckets, index: int):
        """Return the crossings the canonical rule assigns to step ``index``.

        ``buckets`` is the schedule already grouped by step index by
        :meth:`_route`, so this is a dictionary lookup rather than a second
        opinion.  That is the point twice over: a crossing's step is decided
        once, by :meth:`FinestMesh.step_containing`, and every path that needs
        it reads the same answer — and it is decided once *in time* as well.
        The previous form scanned the whole schedule for every step, so a
        window of ``S`` steps carrying ``C`` crossings cost ``S * C``
        comparisons and re-derived the same routing on every pass.
        """
        return buckets.get(index, ())

    def _route(self, ordered):
        """Group the frozen schedule into ``{step index: crossings}``.

        Half-open by construction, because :meth:`FinestMesh.step_containing`
        is: a crossing exactly on a boundary belongs to the step that *starts*
        there.  Comparing against independently computed step ends was what let
        one crossing be owned by two steps at once.

        ``ordered`` arrives sorted by time and ``step_containing`` is monotone
        in time, so each bucket comes out chronological without a second sort,
        and the buckets themselves are in the order their first crossing
        occurs.  Built once per request and handed to the preflight and the
        generation walk alike, so both read one routing.
        """
        buckets: dict = {}
        for crossing in ordered:
            buckets.setdefault(self.mesh.step_containing(crossing.time),
                               []).append(crossing)
        return buckets
