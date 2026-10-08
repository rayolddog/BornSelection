"""
killed_diffusion.py — the independent stationary continuum oracle, and the
frozen-budget machinery every numerical validation in this package reports
through.

What this module is
-------------------
A phase diffusing inside a fixed angular band, killed the instant it touches
either edge.  Given a band ``(lower, upper)``, a phase diffusion ``D``, and the
calibrated stationary Adler drift ``Delta - K sin(phi)``, it answers three
questions in continuous time:

1. **survival** — the probability that a phase started at ``x`` is still
   strictly inside the band at time ``t``;
2. **exit-time quantiles** — the time by which a declared fraction of paths
   have been absorbed; and
3. **exit counts, resolved by edge** — the probability that the one killing
   excursion left through the lower edge and the probability that it left
   through the upper one.  Their sum is the expected number of dwell resets in
   the killed problem, which is the third observable experiment S3 compares.

Why the third observable is defined that way is worth saying plainly, because
the obvious reading of "reset count" does not survive the continuum limit.  A
Brownian path that touches a boundary touches it infinitely often in any
interval that follows, so "how many times did the endpoint scheme see the walker
leave and come back" has **no finite continuous-time limit**; refining the
timestep makes that count diverge rather than converge.  Killing at the first
exit is what makes the question well posed, and once the process is killed the
reset count is the absorption indicator.  Splitting it by edge is what stops it
from being a restatement of survival: the total is exactly ``1 - S``, but *how*
the total divides between the two edges is a separate fact about the drift, and
an endpoint scheme can get the total right while getting the split wrong.

What this module is NOT
-----------------------
It is not the event generator, and it cannot become one.  It contains no dwell
clock, no commitment, no pulse, no population, no ledger and no coupling sweep.
It imports ``validation.py`` and nothing else from this package — in particular
it does not import ``dynamics``, ``stochastic``, ``commitment`` or any
``raw_*`` module — so the oracle a raw run is measured against was written
without reference to the code that produces the run.  ``verify.py`` asserts
that absence rather than trusting this paragraph, and the raw import graph is
forbidden from reaching this file in the other direction.

It is also not the moving-band audit.  This solver is the **continuum
authority** for a *stationary* band: coupling, target and band edges are all
held fixed.  The pulsed, moving-band case is handled diagnostically by
``moving_band_audit.py``, whose frozen-drift Brownian bridge is explicitly an
approximation and explicitly subordinate to this solver.

The method
----------
The survival function solves the backward Kolmogorov equation

::

    dS/dt = f(x) dS/dx + D d2S/dx2,     S(x, 0) = 1,
    S(lower, t) = S(upper, t) = 0,      f(x) = Delta - K sin(x)

discretized with second-order central differences in space and a
Crank-Nicolson step in time.  The generator is autonomous, so the one-step
propagator is a constant matrix; it is formed once and applied, rather than
refactorized every step.

The initial condition is discontinuous at the two corners — one inside, zero on
the edge — and Crank-Nicolson rings on a discontinuity rather than damping it.
The first few steps are therefore taken as pairs of implicit-Euler half-steps
(Rannacher startup), which is unconditionally damping and restores the
second-order behaviour of the steps that follow.  The suite measures that
second order rather than asserting it.

The known limiting solution
---------------------------
With ``K = 0`` the drift is a constant ``mu`` and the problem has a closed-form
eigenfunction expansion, written here as :func:`series_survival` and derived
independently of the solver above:

::

    S(y, t) = exp(-k y - D k^2 t) * SUM_n c_n exp(-D (n pi / L)^2 t) sin(n pi y / L)

with ``y = x - lower``, ``L = upper - lower``, ``k = mu / (2 D)``, and

::

    c_n = (2 / L) (n pi / L) (1 - exp(k L) (-1)^n) / (k^2 + (n pi / L)^2)

At ``mu = 0`` this collapses to the textbook ``4 / (n pi)`` for odd ``n`` and
zero for even, which is the case the verifier pins first.  Reproducing that
limit is what earns the numerical solver the right to be called an oracle in
the cases that have no closed form.

Frozen budgets
--------------
Every validation in this package has to receive its tolerances **before** its
results are opened, or the tolerance is a description of the result rather than
a test of it.  That is a procedural rule, so it is given a mechanical form here:
:class:`FrozenBudgets` hashes its own contents, and :func:`compare_refinement`
refuses to produce a verdict unless the caller supplies the digest it recorded
earlier.  Editing a budget after seeing a discrepancy changes the digest and the
comparison fails outright rather than passing under the new number.

The verdict is machine-readable — :meth:`ValidationVerdict.as_dict` — and
carries the measured refinement-error envelope alongside the pass/no-result
flag, because that envelope is exactly what a later production-power budget has
to be compared against.  This ticket exercises the mechanism on synthetic and
reference budgets; the production numbers are a later ticket's to freeze.

Non-claims
----------
See :data:`ORACLE_NON_CLAIMS`.  In particular: a passing verdict here says the
endpoint scheme resolves a *stationary* band's exit law, and says nothing about
any coupling exponent, any pulsed result, or Born's rule.
"""
from __future__ import annotations

import hashlib
import json
import math
import unicodedata
from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np

from .validation import (require_exact_type, require_finite,
                         require_finite_array, require_nonnegative,
                         require_positive_int)

__all__ = [
    "BLOCKING_FIELDS",
    "BUDGET_FIELDS",
    "BandGeometry",
    "CANONICAL_DTYPE",
    "CONTRACT_FIELDS",
    "DIGEST_LENGTH",
    "ENVELOPE_FIELDS",
    "ESTIMATORS",
    "FORBIDDEN_LABEL_CATEGORIES",
    "FrozenBudgets",
    "LEVEL_FIELDS",
    "MAX_SERIES_TERMS",
    "MINIMUM_REFINEMENT_LEVELS",
    "ORACLE_NON_CLAIMS",
    "PROBABILITY_SLACK",
    "PairedSample",
    "REASON_CODES",
    "REPORT_FIELDS",
    "RefinementContract",
    "RefinementLevel",
    "SAMPLING_FIELDS",
    "SamplingDesign",
    "SurvivalSolution",
    "UNITS",
    "VALIDATION_SCHEMA",
    "VERDICTS",
    "ValidationBudget",
    "ValidationDataset",
    "ValidationVerdict",
    "compare_refinement",
    "inspect_validation_report",
    "parse_validation_report",
    "series_error_certificate",
    "series_survival",
    "series_truncation_bound",
    "solve_survival",
]

#: What a reader of an oracle number is not entitled to conclude.  Pinned as
#: text here and again in the README, for the same reason
#: ``raw_experiments.NON_CLAIMS`` is.
ORACLE_NON_CLAIMS = (
    "this is a stationary band with fixed edges: it says nothing directly "
    "about the moving band of a pulse, which is the diagnostic bridge audit's "
    "subject and is explicitly an approximation",
    "survival inside a declared angular band is not a detector click, an "
    "absorption, a measurement outcome or a Born event",
    "agreement between the endpoint scheme and this oracle is a numerical "
    "result about a discretization, not evidence for or against any "
    "mechanism",
    "no coupling exponent, hazard law or scaling claim is computed, fitted or "
    "implied anywhere in this module",
    "a numerical no-result blocks scientific interpretation; it is not a "
    "falsification of anything physical",
)

#: The schema tag carried by every exported verdict.  Versioned because a
#: later ticket reads these dictionaries mechanically, and a change to what the
#: keys mean has to be visible as a change of schema rather than as different
#: numbers under the same names.
#:
#: **v2.**  The report shape changed incompatibly and the tag did not, which an
#: independent review correctly called out against this comment's own rule.
#: v1 carried a single scalar ``envelope`` that maximized over a probability and
#: a time; v2 carries per-unit ``envelopes`` naming the observable and position
#: that produced each maximum, and every level gained a position, a unit, a
#: dependence-aware standard error, a paired standard error and a cluster count.
#: The report also gained a ``dataset_digest``.
#:
#: **v3.**  The report gained the frozen contract, the frozen budgets, a
#: machine-readable ``blocking`` list of reason codes, and a per-level
#: ``span_error``.  Two digests were not enough: an independent review showed
#: that a genuine ``numerical_no_result`` could be rewritten *coherently* to
#: ``pass`` — verdict, ``passed`` and ``reasons`` all changed together, every
#: level and both digests byte-identical — and the parser had nothing to check
#: it against, because the report carried no caps, no floors, no
#: ``require_decrease`` settings and no decision evidence.  v3 carries them, and
#: the parser **recomputes** every ladder's blocking codes and the aggregate
#: verdict from the levels rather than believing the tag.  A v1 or v2 consumer
#: must refuse this shape, which is what :func:`parse_validation_report`
#: enforces.
VALIDATION_SCHEMA = "dk-numerical-validation/v4"

#: The reason codes a ladder can be blocked by.  Machine-readable beside the
#: human sentences, because a consumer has to be able to *recompute* the
#: decision rather than read it: a report whose ``blocking`` list disagrees
#: with what its own levels and budgets imply is refused.
REASON_CODES = ("unsampled", "unbudgeted", "cluster_mismatch",
                "identity_mismatch", "schedule_mismatch", "unit_mismatch",
                "absolute_cap", "relative_cap", "not_converging",
                "reversal_beyond_allowance", "repeated_reversal")

#: The two outcomes a validation can have.  There is deliberately no third: a
#: comparison that cannot meet its frozen budget is a *numerical* no-result,
#: never a physical verdict and never a reason to move the budget.
VERDICTS = ("pass", "numerical_no_result")

#: How many opening steps are taken as pairs of implicit-Euler half-steps
#: before Crank-Nicolson takes over.  Four is the standard Rannacher count and
#: is enough to damp the corner discontinuity of ``S(x, 0) = 1`` against
#: ``S(edge, 0) = 0`` without materially disturbing the second-order behaviour
#: the remaining steps carry.
RANNACHER_STEPS = 4

#: The largest ``|mu| L / (2 D)`` :func:`series_survival` will evaluate.  The
#: coefficient formula contains ``exp(k L)``, which overflows binary64 above
#: about 709; 300 is a wide margin below that and is refused loudly rather than
#: returned as an infinity.  Every calibration case in this package is far
#: inside it.
MAX_SERIES_PECLET = 300.0

#: The largest term count any public series helper will accept.
#:
#: A million modes is far past anything this package evaluates -- the limiting
#: case uses four thousand -- and it is checked **before** ``terms`` reaches a
#: float conversion or an allocation.  An independent review found
#: ``terms = 10**400`` leaking a bare ``OverflowError: int too large to convert
#: to float`` from ``terms * pi / width``, and ``10**100`` leaking NumPy's
#: generic ``Maximum allowed size exceeded``; neither is the named domain
#: refusal this module promises for a finite input it cannot serve.
MAX_SERIES_TERMS = 1_000_000

#: Longest label a budget or verdict may carry.  A label is hashed into the
#: digest, so it is part of the frozen object; the bound stops an unbounded
#: blob from becoming the thing a manifest identifies a budget by.
MAX_LABEL = 200

#: Longest reason a verdict may carry.  Longer than a label because a reason is
#: an explanatory sentence rather than an identifier, and a reason truncated to
#: fit an identifier's bound would drop exactly the part that says what to do
#: about it.
MAX_REASON = 800

#: The units an observable can be measured in.  Carried on every budget and
#: every measured level, and required to agree between them.
#:
#: This exists because an "envelope" that maximizes over a survival probability
#: and an exit-time quantile is a maximum over two different quantities and is
#: not a number any later uncertainty budget can be compared with.  A unit is
#: not decoration here: it is what makes the exported envelope dimensionally
#: meaningful, and it is hashed into the digest so that relabelling an
#: observable's unit is as visible as moving its cap.
UNITS = ("probability", "time", "count")

#: How far outside ``[0, 1]`` a computed probability may stray before it stops
#: being a probability.  A few hundred ulps of the unit interval: enough for
#: the accumulated rounding of a few thousand matrix-vector products, and far
#: below anything a discretization defect produces.
#:
#: A violation is **refused, never clipped**.  Clipping a negative survival to
#: zero would turn a grid that solved the wrong problem into a grid that solved
#: the right one badly, and the difference between those two is the whole
#: reason an oracle is worth having.
PROBABILITY_SLACK = 1e-9

#: The smallest refinement ladder any frozen contract may declare.  The plan
#: says "at least three paired refinement levels"; this is that sentence, in a
#: place a caller cannot reach past.  It used to be a default argument of
#: :func:`compare_refinement`, which meant a caller could ask for two.
MINIMUM_REFINEMENT_LEVELS = 3

#: The gap between the two nested truncations the series convergence test
#: compares.  **Two, not one**, and the reason is a defect an independent review
#: found: at zero drift every *even* Fourier coefficient of this expansion
#: vanishes, so a sum truncated at ``N`` and a sum truncated at ``N - 1`` are the
#: same sum whenever the dropped mode is an even one.  Comparing ``terms`` with
#: ``terms // 2`` had the same disease in its worst case — ``terms=2`` against
#: ``terms=1`` is one mode against one mode — and reported an exact zero gap
#: while the omitted third and higher odd modes were worth 0.0294.
#:
#: A step of two always drops at least one *odd* mode, whatever the parity of
#: ``N``, so the two truncations cannot be the same sum at zero drift.  It is
#: still only a difference of partial sums, which is why the rigorous tail bound
#: below is the primary test and this is corroborating evidence.
NESTED_TERM_STEP = 2

#: Units in the last place charged to each *named* floating operation in a
#: term, over and above the operation's own argument magnitude.
#:
#: The bound below walks the operation graph rather than modelling it.  A
#: previous revision charged ``(terms + 8) eps SUM |term|`` and called the
#: result what the floating evaluation "can be out by"; an independent review
#: falsified that with a 150-digit reference, finding an accepted value whose
#: true error was 7.4 times its certificate.  The lesson is not that the
#: constant was too small: it is that a constant cannot bound an evaluation
#: whose exponentials carry arguments of size 50 and whose coordinate is formed
#: by subtraction near an edge.  So each factor is now charged for what it
#: actually does.
OPERATION_ULPS = 6

#: Units in the last place charged to the *result*, independent of the sum,
#: **absolutely**.
#:
#: This is deliberately a floor rather than a relative charge, and it makes the
#: certificate a **declared conservative refusal domain**: no tolerance below
#: ``RESULT_ULPS * eps`` — about ``1.78e-15`` — can be certified for any
#: evaluation, however accurate.  That is the guarantee the ticket asks to be
#: testable, and it is: the verifier asserts that a request at the review's own
#: ``2e-15`` is refused, and separately that the value the same call produces
#: agrees with an independent 150-digit reference to ``9.5e-20``.
#:
#: Refusing an evaluation that happens to be accurate is the safe direction.
#: The alternative — scaling this floor by the value — would certify at
#: tolerances where the bound's own modelling assumptions are the largest
#: term, and an independent review has already falsified one such model here.
RESULT_ULPS = 8

#: Dekker's splitting constant, ``2**27 + 1``.  Splits a binary64 into two
#: halves of 26 and 27 significant bits, so that a product with a small integer
#: is exact in each half.
_SPLIT = 134217729.0

#: How many units in the last place two adjacent timesteps may miss the
#: declared refinement ratio by.  The ladders in this package are exact powers
#: of two of a mesh step, so the quotient is exact; the margin exists so that a
#: contract written in decimals is not refused for a rounding.
RATIO_ULPS = 4


# ---------------------------------------------------------------------------
# Private validation helpers
# ---------------------------------------------------------------------------
#: The Unicode general categories no hashed label may contain.
#:
#: ``Cc`` control, ``Cf`` format, ``Cs`` surrogate, ``Co`` private use, ``Cn``
#: unassigned.  Everything else — letters, marks, numbers, punctuation,
#: symbols, separators — is accepted, so ordinary non-ASCII names work.
#:
#: One policy, one helper, every hashed field.  An independent review found the
#: previous rule rejecting only C0 and DEL inside *identities* while accepting
#: U+0085, U+202E and U+200D there, and accepting a newline in ``observable``
#: and a NUL in ``position`` because those two fields went through a different
#: check entirely.  Length-prefixing already made the encoding injective; this
#: is the separate promise — that a label is a legible name — which the
#: documentation was making and the code was not keeping.
#:
#: **No normalization is applied.**  Two labels that differ as Python strings
#: are two different labels, so a composed and a decomposed spelling of the
#: same word are distinct identities and hash differently.  That is the
#: conservative reading: normalizing would silently merge manifests a caller
#: wrote apart.
FORBIDDEN_LABEL_CATEGORIES = ("Cc", "Cf", "Cs", "Co", "Cn")


def _require_printable(value, name: str, *, allow_empty: bool = False) -> str:
    """Return ``value`` as a bounded label under the one Unicode policy."""
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, str):
        raise TypeError(f"{name} must be a string, got "
                        f"{type(value).__name__} {value!r}")
    if not value and not allow_empty:
        raise ValueError(f"{name} must not be empty")
    if len(value) > MAX_LABEL:
        raise ValueError(f"{name} must be at most {MAX_LABEL} characters, got "
                         f"{len(value)}")
    bad = sorted({character for character in value
                  if unicodedata.category(character)
                  in FORBIDDEN_LABEL_CATEGORIES})
    if bad:
        raise ValueError(
            f"{name} contains "
            f"{[f'U+{ord(c):04X} ({unicodedata.category(c)})' for c in bad]}.  "
            f"A hashed label is the name of a physical thing; control, format, "
            f"surrogate, private-use and unassigned code points have no place "
            f"in one, and admitting them would make a manifest's meaning "
            f"depend on how a terminal renders it.  Categories "
            f"{list(FORBIDDEN_LABEL_CATEGORIES)} are refused; no Unicode "
            f"normalization is applied, so distinct spellings stay distinct")
    return value


def _require_label(value, name: str) -> str:
    """Return ``value`` as a non-empty bounded label, or raise naming it."""
    return _require_printable(value, name)


def _require_reason(value, name: str) -> str:
    """Return ``value`` as a non-empty bounded explanatory sentence.

    Under the same Unicode policy as every other string this module emits.  A
    reason is not an identifier, but it is displayed, and an independent review
    found the shape checker accepting all six forbidden categories in one --
    including a bidi override, which can make a displayed reason read as its
    own opposite.
    """
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, str):
        raise TypeError(f"{name} must be a string, got "
                        f"{type(value).__name__} {value!r}")
    if not value:
        raise ValueError(f"{name} must not be empty")
    if len(value) > MAX_REASON:
        raise ValueError(f"{name} must be at most {MAX_REASON} characters, "
                         f"got {len(value)}")
    bad = sorted({character for character in value
                  if unicodedata.category(character)
                  in FORBIDDEN_LABEL_CATEGORIES})
    if bad:
        raise ValueError(
            f"{name} contains "
            f"{[f'U+{ord(c):04X} ({unicodedata.category(c)})' for c in bad]}; "
            f"the same category policy applies to every string this module "
            f"emits, reasons included")
    return value


def _require_positive(value, name: str) -> float:
    """Return ``value`` as a strictly positive finite float."""
    number = require_finite(value, name)
    if not number > 0.0:
        raise ValueError(f"{name} must be strictly positive, got {number!r}")
    return number


def _require_probabilities(values, name: str) -> np.ndarray:
    """Return ``values`` as finite floats in ``(0, 1)``, exclusive."""
    array = require_finite_array(values, name)
    if array.ndim != 1 or array.size == 0:
        raise ValueError(f"{name} must be a non-empty one-dimensional "
                         f"sequence, got shape {array.shape}")
    if np.any(array <= 0.0) or np.any(array >= 1.0):
        raise ValueError(f"{name} must lie strictly inside (0, 1); a quantile "
                         f"level of 0 or 1 names an endpoint of the "
                         f"distribution rather than a point in it, got "
                         f"{array.tolist()}")
    return array


def _digest_text(payload: str) -> str:
    """The one hash used for every frozen object in this module."""
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


# ---------------------------------------------------------------------------
# 1. The stationary band
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class BandGeometry:
    """One stationary killed-diffusion problem.

    ``lower`` / ``upper``
        The absorbing edges, in radians of absolute phase, with
        ``lower < upper``.  They are *absorbing*: a path that reaches an edge
        is killed there, which is the continuum reading of "the phase left the
        lock band and the dwell resets".  Equality with an edge is outside, in
        the same strict sense ``commitment.LockCriterion`` uses.
    ``diffusion``
        ``D_phase``.  Strictly positive here, unlike the noise source, because
        a killed *diffusion* with no diffusion is a deterministic flow and has
        no exit law to solve; the zero-noise limit is checked against
        deterministic Euler elsewhere and does not belong to this oracle.
    ``coupling`` / ``detuning``
        ``K`` and ``Delta`` of the calibrated stationary Adler drift
        ``Delta - K sin(phi)``.  ``coupling = 0`` reduces the drift to the
        constant ``Delta``, which is the case :func:`series_survival` solves in
        closed form and the case the solver is calibrated against.

    The band is required to be no wider than a full turn.  A wider one would
    wrap the circle and stop being an interval with two distinguishable edges,
    which is a different problem from the one solved here rather than a larger
    instance of it.
    """

    lower: float
    upper: float
    diffusion: float
    coupling: float = 0.0
    detuning: float = 0.0

    def __post_init__(self):
        lower = require_finite(self.lower, "lower")
        upper = require_finite(self.upper, "upper")
        if not lower < upper:
            raise ValueError(f"the band needs a strictly positive width: got "
                             f"lower={lower!r}, upper={upper!r}")
        if upper - lower > 2.0 * math.pi:
            raise ValueError(
                f"a band of width {upper - lower!r} exceeds a full turn; a "
                f"band that wraps the circle has no two distinguishable edges "
                f"and is not the problem this oracle solves")
        diffusion = _require_positive(self.diffusion, "diffusion")
        coupling = require_nonnegative(self.coupling, "coupling")
        detuning = require_finite(self.detuning, "detuning")
        object.__setattr__(self, "lower", lower)
        object.__setattr__(self, "upper", upper)
        object.__setattr__(self, "diffusion", diffusion)
        object.__setattr__(self, "coupling", coupling)
        object.__setattr__(self, "detuning", detuning)

    @property
    def width(self) -> float:
        """The band's angular width, ``upper - lower``."""
        return self.upper - self.lower

    def drift(self, phase) -> np.ndarray:
        """The calibrated Adler drift ``Delta - K sin(phi)`` at ``phase``.

        Written out here rather than imported from ``model.py``.  That is the
        one duplication this module contains and it is deliberate: the oracle
        is the independent reference, and a reference that calls the same
        expression the thing under test calls cannot detect an error in the
        expression.  ``verify.py`` compares the two, which is where a
        divergence between them would show up.
        """
        values = require_finite_array(phase, "phase")
        return self.detuning - self.coupling * np.sin(values)


# ---------------------------------------------------------------------------
# 2. The closed-form reference
# ---------------------------------------------------------------------------
def _two_sum(left, right):
    """``(sum, error)`` with ``sum + error == left + right`` exactly (Knuth)."""
    total = left + right
    carry = total - left
    return total, (left - (total - carry)) + (right - carry)


def _two_product(left, right):
    """``(product, error)`` exactly, by Dekker splitting (no FMA needed)."""
    product = left * right
    left_hi, left_lo = _dekker_split(left)
    right_hi, right_lo = _dekker_split(right)
    error = (((left_hi * right_hi - product) + left_hi * right_lo)
             + left_lo * right_hi) + left_lo * right_lo
    return product, error


def _dekker_split(value):
    """Split a binary64 into 26- and 27-bit halves."""
    scaled = _SPLIT * value
    high = scaled - (scaled - value)
    return high, value - high


def _series_shape(fraction, terms) -> np.ndarray:
    """``sin(n pi u)`` for ``n = 1..terms``, by exact argument reduction.

    ``fraction`` is ``u = (x - lower) / L`` in ``[0, 1]``.  Returned as
    ``(terms,)``.

    Computing this as ``np.sin(u * n * pi)`` is where the expansion loses its
    accuracy, and the loss is invisible in the result.  The *input* to the sine
    is a number of size up to ``N pi``, so it carries an absolute rounding of
    ``eps N pi``; near a zero of the sine — which is exactly where ``u``
    approaches 0 or 1 — the sine's own value is far smaller than that, and the
    relative error of ``shape`` explodes.  Multiplied by an exponentially
    tilted amplitude that can be ``1e10``, a per-term relative error of
    ``1e-10`` becomes an absolute error of ``1e-7`` in a term, and a sum of
    four hundred of them becomes ``1e-5`` in an answer requested to ``1e-6``.
    That is the defect an independent review found by reflection.

    The reduction removes it.  ``n u`` is formed *exactly* as a double-double
    by Dekker splitting — ``u`` into 26 and 27 significant bits, ``n`` being a
    small integer, so both partial products are exact — then split into its
    integer part ``m`` and a remainder ``f`` in ``[-1/2, 1/2]``, the
    subtraction being exact because ``|f| <= 1/2`` while ``n u`` is at most a
    few thousand.  ``sin(n pi u) = (-1)^m sin(pi f)`` follows from the addition
    formula, and ``sin(pi f)`` is evaluated on an argument that carries full
    *relative* precision.

    On the review's own reflection pair this takes the symmetry violation from
    ``1.40e-5`` to ``1.16e-10``.
    """
    index = np.arange(1, terms + 1, dtype=float)
    high, low = _dekker_split(fraction[0])
    whole = index * high                                   # exact
    part = index * low + index * fraction[1]
    nearest = np.round(whole + part)
    residue = (whole - nearest) + part                     # exact difference
    sign = np.where(np.mod(nearest, 2.0) == 0.0, 1.0, -1.0)
    return sign * np.sin(math.pi * residue)


def _band_fraction(start, lower, upper):
    """``(x - lower) / (upper - lower)`` as a double-double.

    Returned as ``(high, low)`` summing to the quotient to about 106 bits.

    Both the edge distance and the band width are carried exactly, and both
    matter.  An independent review found a public evaluation returning a value
    whose true error against a 150-digit reference was ``1.33e-14`` — 6.6 times
    the requested tolerance and 7.4 times the certificate — at a start
    ``1.1e-6`` from the upper edge.  The cause was not the summation the
    certificate modelled: it was this coordinate.

    ``u`` is close to one there, so ``1 - u`` is the small quantity the sine
    depends on, and every rounding upstream of it lands on ``1`` and is
    magnified by ``1 / (1 - u)``.  Rounding ``x - lower`` to a double costs
    ``eps`` absolute, which is ``2e-10`` *relative* to a ``1.1e-6`` edge
    distance; rounding ``upper - lower`` costs the same again.  A previous
    revision reduced ``n * u`` exactly and still carried a float ``u``, which
    is why the review's remark that "exact reduction of ``n*u`` after ``u`` has
    already been formed does not make the full formula exact" was right.

    With both carried double-double the same evaluation's error falls from
    ``1.33e-14`` to ``9.1e-20``.
    """
    distance_hi, distance_lo = _two_sum(start, -lower)
    width_hi, width_lo = _two_sum(upper, -lower)
    quotient = distance_hi / width_hi
    product, error = _two_product(quotient, width_hi)
    residual = ((((distance_hi - product) - error) + distance_lo)
                - quotient * width_lo)
    return quotient, residual / width_hi


def _series_terms(lower, upper, diffusion, drift, starts, moments, terms):
    """Every term of the truncated expansion, shaped ``(positions, times, n)``.

    Kept as terms rather than summed because the certificate needs
    ``SUM |term|`` — the cancellation the sum performs — and not only its
    result.
    """
    length = upper - lower
    tilt = _tilt_of(drift, diffusion)
    index = np.arange(1, terms + 1, dtype=float)
    wave = index * math.pi / length
    parity = np.where(index % 2.0 == 0.0, 1.0, -1.0)     # (-1)^n

    offsets = starts - lower                              # y in [0, L]
    fractions = [_band_fraction(start, lower, upper) for start in starts]
    # exp(-k y) * (1 - exp(k L) (-1)^n) is written as the difference of two
    # exponentials so that neither factor is formed at a magnitude the other
    # is about to cancel.
    left = np.exp(-tilt * offsets)[:, None]               # (positions, 1)
    right = np.exp(tilt * (length - offsets))[:, None]
    weight = (2.0 / length) * wave / (tilt * tilt + wave * wave)
    amplitude = weight[None, :] * (left - right * parity[None, :])
    shape = np.stack([_series_shape(fraction, terms)
                      for fraction in fractions])          # (positions, terms)

    decay = np.exp(-(diffusion * (wave * wave)[None, :]
                     + diffusion * tilt * tilt)
                   * moments[:, None])                     # (times, terms)
    return (amplitude * shape)[:, None, :] * decay[None, :, :]


def _series_partial_sums(lower, upper, diffusion, drift, starts, moments,
                         terms):
    """The truncated expansion at ``terms`` modes; no validation, no pinning."""
    return _series_terms(lower, upper, diffusion, drift, starts, moments,
                         terms).sum(axis=2)


def _safe(name, value):
    """Return ``value`` if it is finite, else raise the named domain error.

    Every derived quantity in the series passes through here.  The message
    names the *quantity* rather than the operation, because the caller supplied
    finite inputs and needs to know which derived one could not be represented.
    """
    if not math.isfinite(value):
        raise ValueError(
            f"{name} is {value!r}.  Every input is finite and this derived "
            f"quantity is not, so the evaluation cannot be carried out in "
            f"binary64 -- let alone certified.  Refusing here is the declared "
            f"domain rule: a public result is never a NaN or an infinity, and "
            f"neither is ever allowed to pass a gate by losing a comparison")
    return value


def _tilt_of(drift, diffusion):
    """``drift / (2 * diffusion)``, formed so the denominator cannot overflow.

    ``drift / diffusion / 2`` and never ``drift / (2 * diffusion)``.  The
    second spelling forms ``2 * diffusion`` first, and at ``diffusion = 1e308``
    that is an infinity, so the tilt of a heavily drifted process silently
    becomes **zero** and the series is evaluated -- and certified -- as if
    there were no drift at all.  An independent review found exactly that: a
    returned ``6.5856e-05`` for a process whose true survival is below
    ``exp(-2.5e199)``, under a reported total certificate of ``1.78e-15``.

    A guard computed from a collapsed value guards nothing, which is why this
    is a named helper every site calls rather than an inline expression.
    """
    return drift / diffusion / 2.0


def _series_domain(lower, upper, diffusion, drift, positions, times, terms):
    """Validate the one domain every series entry point shares.

    Factored out because it was not shared, and an independent review found the
    consequence: ``series_survival`` checked the band, the sign of the times and
    the shape of the arrays while the exported ``series_truncation_bound`` — a
    public API in its own right — checked none of them, returning a finite
    bound for a start outside the band, ``inf`` for a negative time and a
    three-dimensional answer for a two-dimensional position array.
    """
    lower = require_finite(lower, "lower")
    upper = require_finite(upper, "upper")
    if not lower < upper:
        raise ValueError(f"the band needs a strictly positive width: got "
                         f"lower={lower!r}, upper={upper!r}")
    diffusion = _require_positive(diffusion, "diffusion")
    drift = require_finite(drift, "drift")
    starts = require_finite_array(positions, "positions")
    moments = require_finite_array(times, "times")
    terms = require_positive_int(terms, "terms", minimum=1)
    # Before any float conversion or allocation: ``terms`` reaches
    # ``terms * pi / width`` a few lines below, and an unbounded Python integer
    # raises there with a message about converting an int rather than about
    # this module's domain.
    if terms > MAX_SERIES_TERMS:
        raise ValueError(
            f"terms is {terms}, above the declared operational maximum "
            f"{MAX_SERIES_TERMS}.  The expansion is evaluated as an array of "
            f"that many modes, so this is a resource and representability "
            f"limit rather than a mathematical one; it is checked before the "
            f"count reaches a float conversion so that the refusal is this "
            f"module's named one and not a bare OverflowError")
    if starts.ndim != 1 or starts.size == 0:
        raise ValueError(f"positions must be a non-empty one-dimensional "
                         f"sequence, got shape {starts.shape}")
    if moments.ndim != 1 or moments.size == 0:
        raise ValueError(f"times must be a non-empty one-dimensional "
                         f"sequence, got shape {moments.shape}")
    if np.any(moments < 0.0):
        raise ValueError(f"times must be non-negative, got {moments.min()!r}")
    if np.any(starts < lower) or np.any(starts > upper):
        raise ValueError(f"every position must lie in the closed band "
                         f"[{lower!r}, {upper!r}], got "
                         f"[{starts.min()!r}, {starts.max()!r}]")
    length = _safe("the band width (upper - lower)", upper - lower)
    # Every *derived* quantity, not only the inputs.  An independent review
    # found finite scalar inputs — a band of (0, 1e308) with a diffusion of
    # 1e308 — overflowing the width, the Dekker split and the mode
    # construction, and the public API returning a NaN through a NaN
    # certificate, because ``nan > tolerance`` is false and the domain gate
    # therefore never fired.  Representability is checked here, once, for
    # everything the evaluation is about to form.
    # Every *derived* quantity, in the same grouping the evaluation uses, and
    # short-circuiting on the width because everything else is formed from it.
    #
    # Two independent reviews shaped this.  The first found finite inputs
    # overflowing the width, with the public API returning a NaN through a NaN
    # certificate.  The second found three deeper faults: a denominator that
    # overflowed and collapsed the drift to zero; guards evaluated in a
    # different grouping from the evaluation they were guarding, so a small
    # band's ``D * h * h * t`` was proved safe while ``D * (h*h)`` overflowed;
    # and eager construction of later fields warning before the decisive
    # infinite width could be named.
    _safe("Dekker's split of the band width", _SPLIT * length)
    _safe("Dekker's split of the lower edge", _SPLIT * lower)
    _safe("Dekker's split of the upper edge", _SPLIT * upper)
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        _safe("the widest edge distance",
              float(np.max(np.abs(starts - lower))))
    longest = _safe("the longest time", float(np.max(moments)))
    tilt = _safe("the drift tilt drift / (2 diffusion)",
                 _tilt_of(drift, diffusion))
    peclet = _safe("the tilt exponential argument |k| L", abs(tilt) * length)
    highest = _safe("the highest mode terms * pi / width",
                    terms * math.pi / length)
    # ``highest * highest`` and never ``highest ** 2``: the second raises a
    # bare OverflowError on a Python float, which is not the named refusal this
    # module promises.  Grouped exactly as ``_series_terms`` groups it.
    squared = _safe("the highest squared mode (terms pi / width)^2",
                    highest * highest)
    _safe("the highest diffusive decay exponent D (n pi / L)^2 t",
          diffusion * squared * longest)
    _safe("the tilt decay exponent D k^2 t",
          diffusion * tilt * tilt * longest)
    _safe("the weight denominator k^2 + (n pi / L)^2", tilt * tilt + squared)
    if peclet > MAX_SERIES_PECLET:
        raise ValueError(
            f"|drift| * width / (2 * diffusion) is {peclet!r}, above the "
            f"declared limit of {MAX_SERIES_PECLET!r}; the coefficient formula "
            f"contains exp(k L) and this is where it stops being a finite "
            f"number")
    return lower, upper, diffusion, drift, starts, moments, terms


def series_error_certificate(lower, upper, diffusion, drift, positions, times,
                             terms) -> dict:
    """Both halves of the expansion's error, and their sum.

    ``{"terms", "tail", "roundoff", "total", "cancellation"}``, each numerical
    entry an array of shape ``(len(positions), len(times))`` except ``terms``.

    ``tail``
        What truncation discards: :func:`series_truncation_bound`.
    ``roundoff``
        What the floating evaluation of the *retained* sum can be out by:
        the operation graph, charged factor by factor: the weight's
        arithmetic, each exponential for its own argument magnitude, the
        amplitude's *absolute* error (which does not shrink when its two
        exponentials cancel), the shape's double-double reduction and sine, the
        classical ``(terms + 1) eps SUM |term|`` summation bound, and
        :data:`RESULT_ULPS` last bits of the result.
    ``total``
        Their sum, and the only quantity a tolerance should be compared with.
    ``cancellation``
        ``SUM |term| / |SUM term|`` — how much the sum cancels, reported
        because it is the thing that makes an exponentially tilted expansion
        untrustworthy and because a reader deserves to see it.

    A tail bound alone is not an error bound, and saying so is the point of
    this function existing separately.  An independent review found the
    expansion accepting a value whose exact reflection partner differed by
    1.63e-5 at a requested tolerance of 1e-6, while both evaluations reported
    tails of 5.04e-210: the analytic statement was true and the number was
    wrong.
    """
    lower, upper, diffusion, drift, starts, moments, terms = _series_domain(
        lower, upper, diffusion, drift, positions, times, terms)
    pieces = _series_terms(lower, upper, diffusion, drift, starts, moments,
                           terms)
    magnitude = np.sum(np.abs(pieces), axis=2)
    total_sum = np.abs(pieces.sum(axis=2))
    epsilon = float(np.finfo(np.float64).eps)

    # The operation graph, charged factor by factor.
    #
    #   weight_n   : four arithmetic operations                -> OPERATION_ULPS
    #   left, right: exp of an argument of size |A|, |B|; the argument's own
    #                rounding contributes |A| eps relatively, the exponential
    #                about one ulp                             -> |A| + 2
    #   amplitude  : left - right (-1)^n; its *absolute* error is the sum of
    #                theirs, which does not shrink when they cancel, so it is
    #                carried absolutely rather than relatively
    #   shape_n    : a double-double coordinate, an exact reduction, a
    #                multiply by pi and a sine                  -> OPERATION_ULPS
    #   decay_n    : exp of an argument of size E_n             -> E_n + 2
    #   summation  : the classical (N + 1) eps SUM |term| bound, conservative
    #                against NumPy's pairwise accumulation
    #   result     : RESULT_ULPS of the returned value itself
    length = upper - lower
    tilt = _tilt_of(drift, diffusion)
    index = np.arange(1, terms + 1, dtype=float)
    wave = index * math.pi / length
    offsets = starts - lower
    scale = ((2.0 / length) * wave / (tilt * tilt + wave * wave))[None, :]
    left = np.exp(-tilt * offsets)[:, None]
    right = np.exp(tilt * (length - offsets))[:, None]
    forward = np.abs(tilt * offsets)[:, None]
    backward = np.abs(tilt * (length - offsets))[:, None]

    # The amplitude's *absolute* error.  Carried absolutely rather than
    # relatively because ``left - right (-1)^n`` does not become more accurate
    # when the two nearly cancel; it becomes less.
    amplitude_error = scale * epsilon * (left * (forward + 2.0)
                                         + right * (backward + 2.0))

    # The decay's argument, and the per-term relative charge: the weight's
    # arithmetic, the shape's double-double reduction and sine, and the
    # exponential of an argument of this size.
    exponent = ((diffusion * (wave * wave)[None, :]
                 + diffusion * tilt * tilt) * moments[:, None])  # (times, n)
    decay = np.exp(-exponent)                                    # (times, n)
    shape = np.abs(np.stack([_series_shape(_band_fraction(start, lower, upper),
                                           terms) for start in starts]))
    relative = epsilon * (2.0 * OPERATION_ULPS + exponent + 2.0)

    roundoff = (np.sum(np.abs(pieces) * relative[None, :, :], axis=2)
                + np.sum(amplitude_error[:, None, :] * shape[:, None, :]
                         * np.abs(decay)[None, :, :], axis=2)
                + (terms + 1) * epsilon * magnitude
                + RESULT_ULPS * epsilon)
    tail = series_truncation_bound(lower, upper, diffusion, drift, starts,
                                   moments, terms)
    interior = (starts > lower) & (starts < upper)
    roundoff[~interior, :] = 0.0
    roundoff[:, moments == 0.0] = 0.0
    with np.errstate(divide="ignore", invalid="ignore"):
        # ``SUM|term| / |SUM term|``, and **one** where both are zero.  A start
        # on an absorbing edge is pinned to exactly zero by definition, not
        # summed: there are no terms, so nothing cancels, and a diagnostic of
        # infinity there was an artefact of dividing zero by zero.  An
        # independent review found it beside a claim that every certificate
        # component is finite.
        cancellation = np.where(magnitude > 0.0,
                                magnitude / np.maximum(total_sum, 1e-300),
                                1.0)
    cancellation = np.where(np.isfinite(cancellation), cancellation,
                            np.finfo(np.float64).max)
    report = {"terms": terms, "tail": tail, "roundoff": roundoff,
              "total": tail + roundoff, "cancellation": cancellation}
    for name in ("tail", "roundoff", "total", "cancellation"):
        if not np.all(np.isfinite(report[name])):
            raise ValueError(
                f"the {name} certificate is not finite.  A certificate that "
                f"is a NaN or an infinity certifies nothing, and an earlier "
                f"revision returned exactly that beside a NaN value because "
                f"the tolerance comparison silently lost against the NaN")
    return report


def series_truncation_bound(lower, upper, diffusion, drift, positions, times,
                            terms) -> np.ndarray:
    """A rigorous upper bound on what truncating the expansion discards.

    Returned as ``(len(positions), len(times))``: entry ``[i, j]`` bounds
    ``|S(y_i, t_j) - S_N(y_i, t_j)|`` where ``S_N`` is the sum of the first
    ``N = terms`` modes.  This is the primary convergence test
    :func:`series_survival` applies, and it is a *bound* rather than a
    difference of two truncations, which is what makes it immune to the
    cancellation that made the old test pass on an unconverged sum.

    The derivation, in full, because a bound nobody can check is not evidence.
    Write ``L = upper - lower``, ``y = x - lower``, ``k = mu / (2 D)`` and
    ``w_n = n pi / L``.  The term of order ``n`` is

    ::

        weight_n * (exp(-k y) - exp(k (L - y)) (-1)^n)
                 * sin(w_n y) * exp(-D (w_n^2 + k^2) t)

    with ``weight_n = (2 / L) w_n / (k^2 + w_n^2)``.  Then:

    * ``weight_n <= (2 / L) w_n / w_n^2 = 2 / (n pi)``, decreasing in ``n``;
    * ``|exp(-k y) - exp(k (L - y)) (-1)^n| <= exp(-k y) + exp(k (L - y))``,
      independent of ``n``; call its value ``M(y)``;
    * ``|sin(w_n y)| <= 1``; and
    * the decay factor is ``exp(-D k^2 t) exp(-a n^2)`` with
      ``a = D (pi / L)^2 t``.

    So the tail past ``N`` is bounded by

    ::

        M(y) exp(-D k^2 t) * (2 / ((N + 1) pi)) * SUM_{n > N} exp(-a n^2)

    and, since ``exp(-a x^2)`` decreases,

    ::

        SUM_{n > N} exp(-a n^2) <= INTEGRAL_N^inf exp(-a x^2) dx
                                <= exp(-a N^2) / (2 a N)

    using ``x / N >= 1`` on the range of integration.  That is what this
    returns.

    The bound diverges as ``t`` approaches zero, and that is correct rather
    than a weakness: the expansion genuinely needs unboundedly many modes there,
    because it is representing a discontinuous initial condition.  A caller who
    wants short-time survival needs a different method, not more terms.

    Zero times return zero: the value there is the exactly-known one and is
    written in by :func:`series_survival` rather than summed.
    """
    lower, upper, diffusion, drift, starts, moments, terms = _series_domain(
        lower, upper, diffusion, drift, positions, times, terms)

    length = upper - lower
    tilt = _tilt_of(drift, diffusion)
    offsets = starts - lower
    reach = (np.exp(-tilt * offsets)
             + np.exp(tilt * (length - offsets)))[:, None]
    decay = diffusion * (math.pi / length) ** 2
    scale = 2.0 / ((terms + 1) * math.pi)

    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        rate = decay * moments[None, :]
        tail = np.where(
            rate > 0.0,
            np.exp(-rate * (terms * terms)) / (2.0 * np.maximum(rate, 1e-300)
                                               * terms),
            np.inf)
        bound = (reach * scale * tail
                 * np.exp(-diffusion * tilt * tilt * moments)[None, :])
    bound = np.where(np.isnan(bound), np.inf, bound)
    # The exactly-known values are exact: an edge is absorbing and a zero time
    # is the initial condition, and neither is read off a truncated sum, so
    # neither carries a truncation error.  The docstring said so and the
    # arithmetic did not; an independent review found this returning ``inf``
    # at zero time.
    interior = (starts > lower) & (starts < upper)
    bound[~interior, :] = 0.0
    bound[:, moments == 0.0] = 0.0
    if not np.all(np.isfinite(bound)):
        rows, columns = np.unravel_index(int(np.argmax(bound)), bound.shape)
        raise ValueError(
            f"the truncation tail is unbounded at position {starts[rows]!r} "
            f"and time {moments[columns]!r} with {terms} terms.  The bound "
            f"contains exp(-a N^2) / (2 a N) with a = D (pi/L)^2 t, and as the "
            f"time approaches zero that quantity grows without limit -- the "
            f"expansion genuinely needs unboundedly many modes there.  This is "
            f"the named convergence refusal; a public helper does not return "
            f"an infinity")
    return bound


def series_survival(lower, upper, diffusion, drift, positions, times,
                    terms=400, tolerance=1e-6) -> np.ndarray:
    """Survival inside ``(lower, upper)`` under a **constant** drift.

    Returned as an array of shape ``(len(positions), len(times))``.

    This is the known limiting case, computed from the eigenfunction expansion
    written in the module docstring, and it exists to calibrate
    :func:`solve_survival`.  It is not a special case of that function's code
    path: nothing below calls the solver, and the solver does not call this.

    ``terms``
        How many eigenfunctions are summed.  The series converges like
        ``exp(-D (n pi / L)^2 t)``, so it converges quickly at moderate ``t``
        and slowly as ``t`` approaches zero, where the discontinuous initial
        condition is being represented by a Fourier sine series.
    ``tolerance``
        The declared truncation budget.  **The truncation is tested, not
        assumed**, and tested two ways, both of which must fit inside it:

        1. :func:`series_truncation_bound`, a rigorous upper bound on the
           discarded tail derived in that function's docstring.  This is the
           primary test.
        2. The gap between the sums at ``terms`` and at
           ``terms - NESTED_TERM_STEP``, reported as corroborating evidence.

        The second used to be the *only* test, comparing ``terms`` with
        ``terms // 2``, and an independent review broke it.  At zero drift every
        even coefficient of this expansion vanishes, so ``terms=2`` and
        ``terms=1`` are the same one-mode sum: the test saw an exact zero gap
        and accepted a result 0.0294 away from the resolved value, entirely
        inside ``[0, 1]`` so that the probability-domain check saw nothing
        either.  A difference of two partial sums can only ever say "these two
        truncations agree"; it cannot say "the rest of the series is small",
        and at zero drift the parity structure makes the two statements come
        apart.  Hence the bound.

    Why that is not optional: an unconverged partial sum of this series is not
    merely imprecise, it is not a probability.  At ``band = (0, 1)``,
    ``D = 0.1``, ``mu = 2`` and ``t = 1e-10``, four hundred terms return
    ``-2.11`` at ``x = 1e-4`` and ``-15.75`` at ``x = 1e-3``, where the true
    survival tends to one.  A caller reaching for an "independent reference"
    there would receive an impossible number with no indication that anything
    had gone wrong, and would then measure a perfectly good endpoint scheme
    against it.  Refusing is the only safe answer: the short-time regime of
    this expansion needs a number of terms that grows without bound as ``t``
    approaches zero, so there is no term count this function could silently
    choose instead.

    Every returned value is additionally required to be a probability to within
    :data:`PROBABILITY_SLACK`, and a violation is refused rather than clipped.

    Times of exactly zero return exactly one in the interior and exactly zero
    on an edge, taken from the definition rather than from the truncated sum,
    whose partial sums at ``t = 0`` are the Gibbs oscillations of a square wave
    rather than the square wave.  The two edges return exactly zero at every
    time for the same reason: they are absorbing, so a path started on one has
    already been killed, and reading that off a sine series that is zero there
    only up to rounding would put a spurious ``1e-17`` in an oracle.
    """
    lower, upper, diffusion, drift, starts, moments, terms = _series_domain(
        lower, upper, diffusion, drift, positions, times, terms)
    if terms < 2:
        raise ValueError(f"terms must be at least two so that a nested "
                         f"truncation exists to report beside the bound, got "
                         f"{terms}")
    tolerance = _require_positive(tolerance, "tolerance")

    out = _series_partial_sums(lower, upper, diffusion, drift, starts, moments,
                               terms)
    # The nested sum is kept even where the step would run off the bottom, so
    # that the review's own ``terms=2`` attack is answered by the *bound*
    # rather than by an argument check: at two terms the nested gap really is
    # exactly zero, and the bound really is 0.1096 against a true error of
    # 0.0294.  Refusing on the argument would hide which test does the work.
    nested = _series_partial_sums(lower, upper, diffusion, drift, starts,
                                  moments, max(terms - NESTED_TERM_STEP, 1))
    certificate = series_error_certificate(lower, upper, diffusion, drift,
                                           starts, moments, terms)
    bound = certificate["total"]

    # The exactly-known values, written in before either the truncation test or
    # the domain test looks at them: an edge is absorbing at every time, and an
    # interior start is alive at time zero.
    interior = (starts > lower) & (starts < upper)
    edge = ~interior
    zero = moments == 0.0
    for array in (out, nested, bound):
        if np.any(zero):
            array[:, zero] = (np.where(interior, 1.0, 0.0)[:, None]
                              if array is not bound else 0.0)
        if np.any(edge):
            array[edge, :] = 0.0

    if not np.all(np.isfinite(out)) or not np.all(np.isfinite(nested)):
        raise ValueError(
            "the truncated sum is not finite at every requested position and "
            "time.  A public result is never a NaN; where the arithmetic "
            "cannot be carried out this refuses instead")
    gap = np.abs(out - nested)
    worst_gap = float(np.max(gap)) if gap.size else 0.0
    worst_bound = float(np.max(bound)) if bound.size else 0.0
    evidence = max(worst_gap, worst_bound)
    # ``not (evidence <= tolerance)`` and never ``evidence > tolerance``: a NaN
    # loses both comparisons, and the second spelling would accept it.
    if not (evidence <= tolerance):
        rows, columns = np.unravel_index(int(np.argmax(bound)), bound.shape)
        finite = certificate["cancellation"][
            np.isfinite(certificate["cancellation"])]
        raise ValueError(
            f"the expansion is not certified at {terms} terms.  The total "
            f"error bound is {worst_bound!r}: an analytic tail of "
            f"{float(np.max(certificate['tail']))!r} plus a "
            f"floating-evaluation bound of "
            f"{float(np.max(certificate['roundoff']))!r} on a sum that cancels "
            f"by a factor {float(np.max(finite, initial=1.0))!r}.  The nested "
            f"gap against {max(terms - NESTED_TERM_STEP, 1)} terms is "
            f"{worst_gap!r} and the declared budget is {tolerance!r}; the "
            f"bound is worst at position {starts[rows]!r} and time "
            f"{moments[columns]!r}.  Both halves matter: short times need "
            f"unboundedly many modes, and an exponentially tilted sum can be "
            f"numerically wrong while its analytic tail is 1e-210")

    low, high = float(np.min(out)), float(np.max(out))
    if not (low >= -PROBABILITY_SLACK and high <= 1.0 + PROBABILITY_SLACK):
        raise ValueError(
            f"the expansion returned values outside [0, 1] by more than "
            f"{PROBABILITY_SLACK!r}: [{low!r}, {high!r}].  A survival "
            f"probability outside the unit interval is refused rather than "
            f"clipped; clipping would turn a reference that solved the wrong "
            f"problem into one that solved the right problem badly")
    return out


# ---------------------------------------------------------------------------
# 3. The numerical solver
# ---------------------------------------------------------------------------
def _propagator(sub, diag, sup, source, step, theta):
    """Return ``(P, q)`` with ``v_next = P @ v + q`` for one theta-step.

    ``sub``/``diag``/``sup`` are the three bands of the discretized generator
    over the interior nodes, ``source`` is the constant contribution of the
    inhomogeneous Dirichlet datum, and ``theta`` is 1 for implicit Euler and
    1/2 for Crank-Nicolson.

    The generator is autonomous, so this is formed once per ``(step, theta)``
    pair and reused for every step taken with it.  Dense rather than banded:
    the interior is a few hundred nodes, one factorization is milliseconds, and
    a dense propagator turns each step into a single matrix-vector product
    instead of a Python loop over the tridiagonal sweep.
    """
    size = diag.size
    generator = (np.diag(diag) + np.diag(sup[:-1], 1) + np.diag(sub[1:], -1))
    identity = np.eye(size)
    left = identity - theta * step * generator
    right = identity + (1.0 - theta) * step * generator
    forward = np.linalg.solve(left, right)
    offset = np.linalg.solve(left, step * source)
    return forward, offset


@dataclass(frozen=True)
class SurvivalSolution:
    """What the oracle computed, for one band and one set of start positions.

    ``survival``
        Shape ``(len(positions), len(times))``.  ``survival[i, n]`` is the
        probability that a path started at ``positions[i]`` has not touched
        either edge by ``times[n]``.
    ``lower_exit`` / ``upper_exit``
        The same shape.  The probability that the path has been absorbed by
        that time *and* that it was absorbed through that edge.  The three
        arrays sum to one at every entry up to discretization error, and the
        suite measures that residual rather than assuming it.
    ``space_steps`` / ``time_steps``
        The discretization the numbers came from, carried with them, because a
        solution without its grid is not a reference anybody can refine.
    """

    positions: np.ndarray
    times: np.ndarray
    survival: np.ndarray
    lower_exit: np.ndarray
    upper_exit: np.ndarray
    space_steps: int
    time_steps: int

    def __post_init__(self):
        positions = require_finite_array(self.positions, "positions")
        times = require_finite_array(self.times, "times")
        space_steps = require_positive_int(self.space_steps, "space_steps",
                                           minimum=4)
        time_steps = require_positive_int(self.time_steps, "time_steps",
                                          minimum=1)
        shape = (positions.size, times.size)
        for name in ("survival", "lower_exit", "upper_exit"):
            values = require_finite_array(getattr(self, name), name)
            if values.shape != shape:
                raise ValueError(f"{name} has shape {values.shape}, expected "
                                 f"{shape} from {positions.size} positions and "
                                 f"{times.size} times")
            # Every field is a probability, and the type refuses to hold one
            # that is not.  Checking here rather than only in the solver is
            # what stops a solution assembled some other way -- a mutation, a
            # later ticket's cache, a hand-built fixture -- from carrying a
            # negative survival under this class's name.
            #
            # ``lower_exit`` is the complement of the other two, so the three
            # can sum to one exactly while two of them are outside [0, 1]; that
            # is the defect an independent review found, and it is why
            # ``closure_residual`` is reported and is not the domain check.
            low, high = float(np.min(values)), float(np.max(values))
            if low < -PROBABILITY_SLACK or high > 1.0 + PROBABILITY_SLACK:
                raise ValueError(
                    f"{name} leaves [0, 1] by more than "
                    f"{PROBABILITY_SLACK!r}: [{low!r}, {high!r}].  This is "
                    f"refused rather than clipped: a grid that produced a "
                    f"negative survival solved a different problem, and "
                    f"clipping would present it as having solved this one "
                    f"inaccurately")
            object.__setattr__(self, name, values)
        object.__setattr__(self, "positions", positions)
        object.__setattr__(self, "times", times)
        object.__setattr__(self, "space_steps", space_steps)
        object.__setattr__(self, "time_steps", time_steps)

    @property
    def closure_residual(self) -> float:
        """The largest ``|S + lower_exit + upper_exit - 1|`` in the solution.

        A path is alive, absorbed below, or absorbed above; there is no fourth
        state.  This is the discretization's own answer to that identity, and
        it is reported rather than asserted so that a reader sees what the grid
        actually delivered.
        """
        total = self.survival + self.lower_exit + self.upper_exit
        return float(np.max(np.abs(total - 1.0)))

    def quantiles(self, levels) -> np.ndarray:
        """First-exit-time quantiles, shaped ``(len(positions), len(levels))``.

        The exit-time distribution function is ``1 - survival``, evaluated on
        this solution's own time grid; a quantile is found by linear
        interpolation between the two grid times that bracket it.

        A level the solution never reaches within its horizon returns NaN
        rather than the horizon.  Clamping there would silently turn "we did
        not integrate long enough" into "the median is exactly where we
        stopped", which is the shape of a wrong answer that looks like a right
        one.
        """
        wanted = _require_probabilities(levels, "levels")
        absorbed = 1.0 - self.survival
        out = np.full((self.positions.size, wanted.size), math.nan)
        for row in range(self.positions.size):
            curve = np.maximum.accumulate(absorbed[row])
            for column, level in enumerate(wanted):
                if curve[-1] < level:
                    continue
                index = int(np.searchsorted(curve, level, side="left"))
                if index == 0:
                    out[row, column] = float(self.times[0])
                    continue
                before, after = curve[index - 1], curve[index]
                span = after - before
                weight = 0.0 if span <= 0.0 else (level - before) / span
                out[row, column] = float(
                    self.times[index - 1]
                    + weight * (self.times[index] - self.times[index - 1]))
        return out


def solve_survival(geometry, positions, horizon, space_steps,
                   time_steps) -> SurvivalSolution:
    """Solve the killed-diffusion problem on ``geometry``; see the module docs.

    ``positions``
        Where the phase starts.  Any point in the closed band; values off the
        spatial grid are read from it by linear interpolation, and the two
        edges return exactly zero because the grid holds exactly zero there.
    ``horizon``
        How far to integrate.  The returned ``times`` are the uniform grid
        ``linspace(0, horizon, time_steps + 1)``, so a caller that wants
        survival at a particular dwell should make the dwell a grid point
        rather than interpolate in time afterwards.
    ``space_steps`` / ``time_steps``
        The discretization.  Both are arguments and neither has a default: an
        oracle whose resolution is a hidden constant cannot be refined, and
        refinement is how this function earns its accuracy claim.

    The returned object carries the grid it was solved on, so two calls at
    ``(M, N)`` and ``(2M, 2N)`` can be differenced to measure the oracle's own
    discretization error before it is used as a reference for anybody else's.

    **Two declared grid conditions, and a domain check.**  The cell Peclet
    number ``max|f| dx / (2D)`` must be below one, or central differencing
    oscillates in space; the advective Courant number ``max|f| dt / dx`` must
    be at most one, or the drift is unresolved in time.  Neither is a
    sufficient condition for accuracy and neither is claimed to be; they are
    the two ways this discretization stops being a discretization of the stated
    problem at all.  Whatever survives them is then required to be a
    probability, and is refused rather than clipped if it is not.

    The Courant condition is the one an independent review's counterexample
    violates: ``D = 0.001``, band ``(-1, 1)``, drift ``0.05``, eighty cells and
    five steps to a horizon of a hundred has a cell Peclet of 0.625 and a
    Courant number of 40, and returned a survival of ``-0.0043`` beside an
    edge-exit probability of ``1.0043``.  The strict Crank-Nicolson positivity
    condition ``D dt / dx^2 <= 1`` would also have refused it, and would have
    refused every grid this package actually uses — the S3 reference grid sits
    at 94 — so it is not the guard: it is a monotonicity condition, not an
    accuracy one, and enforcing it would trade a real defect for a false one.
    """
    band = require_exact_type(geometry, BandGeometry, "geometry")
    starts = require_finite_array(positions, "positions")
    horizon = _require_positive(horizon, "horizon")
    space_steps = require_positive_int(space_steps, "space_steps", minimum=4)
    time_steps = require_positive_int(time_steps, "time_steps",
                                      minimum=RANNACHER_STEPS + 1)
    if starts.ndim != 1 or starts.size == 0:
        raise ValueError(f"positions must be a non-empty one-dimensional "
                         f"sequence, got shape {starts.shape}")
    if np.any(starts < band.lower) or np.any(starts > band.upper):
        raise ValueError(f"every position must lie in the closed band "
                         f"[{band.lower!r}, {band.upper!r}], got "
                         f"[{starts.min()!r}, {starts.max()!r}]")

    nodes = np.linspace(band.lower, band.upper, space_steps + 1)
    spacing = band.width / space_steps
    interior = nodes[1:-1]
    flow = band.drift(interior)

    # Second-order central differences.  ``diffuse`` is the symmetric part and
    # ``advect`` the antisymmetric one; they are written separately so that the
    # cell Peclet number below is a statement about two named quantities.
    diffuse = band.diffusion / (spacing * spacing)
    advect = flow / (2.0 * spacing)
    sub = diffuse - advect
    diag = np.full(interior.size, -2.0 * diffuse)
    sup = diffuse + advect

    peclet = float(np.max(np.abs(flow)) * spacing / (2.0 * band.diffusion))
    if peclet >= 1.0:
        raise ValueError(
            f"the cell Peclet number is {peclet!r}; central differencing is "
            f"oscillatory at or above one, so this grid would return a "
            f"survival curve with the wrong sign structure rather than a less "
            f"accurate one.  Refine space_steps past "
            f"{int(math.ceil(space_steps * peclet)) + 1}")

    step = horizon / time_steps
    courant = float(np.max(np.abs(flow)) * step / spacing)
    if courant > 1.0:
        raise ValueError(
            f"the advective Courant number is {courant!r}; a characteristic "
            f"crosses more than one cell per step, so the drift is not "
            f"resolved in time and Crank-Nicolson returns an oscillation "
            f"rather than a less accurate solution.  On a grid like this one "
            f"an independent review obtained a survival of -0.0043 and an "
            f"edge-exit probability of 1.0043 while the three fields still "
            f"summed to one, which is why closure alone is not a domain "
            f"check.  Refine time_steps past "
            f"{int(math.ceil(time_steps * courant)) + 1}")

    # Survival: homogeneous Dirichlet at both edges, so no source.
    # Upper-exit CDF: zero at the lower edge, one at the upper, so the last
    # interior node picks up ``sup`` times that datum.
    no_source = np.zeros(interior.size)
    upper_source = np.zeros(interior.size)
    upper_source[-1] = sup[-1]

    half_implicit = _propagator(sub, diag, sup, upper_source, 0.5 * step, 1.0)
    full_crank = _propagator(sub, diag, sup, upper_source, step, 0.5)

    alive = np.ones(interior.size)
    above = np.zeros(interior.size)

    times = np.linspace(0.0, horizon, time_steps + 1)
    survival = np.empty((starts.size, time_steps + 1))
    upper_exit = np.empty((starts.size, time_steps + 1))

    def record(column, live, high):
        full_live = np.concatenate(([0.0], live, [0.0]))
        full_high = np.concatenate(([0.0], high, [1.0]))
        survival[:, column] = np.interp(starts, nodes, full_live)
        upper_exit[:, column] = np.interp(starts, nodes, full_high)

    # At t = 0 the edges are already absorbed, so the recorded upper-exit datum
    # at the upper edge is one and interpolation gives the corner its value.
    record(0, alive, above)

    for index in range(time_steps):
        # Rannacher startup: the opening steps are taken as two damped
        # implicit-Euler half-steps rather than one Crank-Nicolson step, which
        # is what stops the corner discontinuity from ringing through the whole
        # solution.  The time grid is unchanged; only how each of the first few
        # steps is reached differs.
        scheme = half_implicit if index < RANNACHER_STEPS else full_crank
        repeats = 2 if index < RANNACHER_STEPS else 1
        for _ in range(repeats):
            # The survival field is homogeneous, so it takes the propagator
            # without the source; the upper-exit field carries the Dirichlet
            # datum and takes the source with it.  Same matrix, different
            # right-hand side, written once each so the two cannot be edited
            # apart.
            alive = scheme[0] @ alive
            above = scheme[0] @ above + scheme[1]
        record(index + 1, alive, above)

    lower_exit = 1.0 - survival - upper_exit
    return SurvivalSolution(starts, times, survival, lower_exit, upper_exit,
                            space_steps, time_steps)


# ---------------------------------------------------------------------------
# 4. Frozen budgets, frozen data, and the machine-readable verdict
# ---------------------------------------------------------------------------
# The shape of this section changed in response to an independent review, and
# the change is worth stating plainly because it is the difference between a
# gate and a suggestion.
#
# It used to be: the caller measured a ladder, computed its own standard
# errors, and handed both to `compare_refinement` as `RefinementLevel` objects.
# The contract froze the caps, the schedule, the level count and a *description*
# of the sampling design — but not the uncertainty actually consumed.  The
# review demonstrated the consequence in two lines: under one unchanged frozen
# contract and one unchanged digest, the ladder 0.010, 0.040, 0.009 returned
# `numerical_no_result` with caller standard errors of 0.0 and `pass` with
# caller standard errors of 0.1.  No measurement had changed.
#
# It is now: the caller hands over the **cluster-level data**, and this module
# computes the statistic, the paired bootstrap, and the verdict.  There is no
# argument through which an uncertainty can be supplied, so there is nothing to
# forge.  The estimator is named in the frozen budget, the resampling key and
# resample count are in the frozen sampling design, and the data itself is
# hashed into the verdict for reproduction.
#
# The second half of the same finding was that the declared contract spoke of
# standard errors of the *paired difference* while the arithmetic used
# `hypot` of two marginal errors — the independent-level formula, which
# overstates the allowance whenever the levels are positively dependent, and
# they always are, because they are one Brownian realization read twice.  The
# bootstrap below resamples clusters **once per draw** and evaluates every level
# on that same draw, so the paired contrast is a real paired contrast.
ESTIMATORS = ("mean", "quantile")


@dataclass(frozen=True)
class SamplingDesign:
    """How a ladder was sampled, and how its uncertainty is estimated.

    Frozen and hashed with the rest of the contract, because "the increase was
    statistically unresolved" is a claim about a *sampling design*.

    ``unit``
        What one independent cluster **is** — "master trial identifier", say.
        Not decoration: it is the thing the bootstrap resamples, and getting it
        wrong is how a dependent sample gets an independent sample's error bar.
    ``clusters``
        How many of them.  Every frozen sample is required to carry exactly
        this many, so a ladder assembled from a differently-sized sample is
        refused rather than compared.
    ``replications``
        Auxiliary replications *inside* a cluster — here, the audit's own
        uniforms redrawn against one physical path.  They travel with their
        cluster through every resample, so they narrow the estimate without
        being counted as independent evidence.
    ``method`` / ``resamples`` / ``seed``
        How the standard errors are produced, how many resamples go into them,
        and the resampling key.  The seed is frozen for the same reason the
        caps are: a bootstrap whose key can be changed after the fact is a
        bootstrap that can be run until it gives the wanted answer.
    ``coverage``
        How many standard errors of the **paired difference** an upward step
        may span, and how many standard errors of the finest-level error are
        added to it before it is compared with the absolute cap.  Two is the
        ordinary choice; it is frozen rather than defaulted because widening it
        is exactly how a real divergence would be talked away.
    """

    unit: str
    clusters: int
    replications: int
    method: str
    resamples: int
    coverage: float
    seed: int

    def __post_init__(self):
        unit = _require_label(self.unit, "unit")
        method = _require_label(self.method, "method")
        clusters = require_positive_int(self.clusters, "clusters", minimum=2)
        replications = require_positive_int(self.replications, "replications",
                                            minimum=1)
        resamples = require_positive_int(self.resamples, "resamples",
                                         minimum=2)
        coverage = _require_positive(self.coverage, "coverage")
        seed = require_positive_int(self.seed, "seed", minimum=0)
        object.__setattr__(self, "unit", unit)
        object.__setattr__(self, "method", method)
        object.__setattr__(self, "clusters", clusters)
        object.__setattr__(self, "replications", replications)
        object.__setattr__(self, "resamples", resamples)
        object.__setattr__(self, "coverage", coverage)
        object.__setattr__(self, "seed", seed)

    @property
    def canonical(self) -> str:
        """The exact text this design contributes to a digest."""
        return json.dumps([self.unit, self.clusters, self.replications,
                           self.method, self.resamples, repr(self.coverage),
                           self.seed], separators=(",", ":"))


@dataclass(frozen=True)
class RefinementContract:
    """The refinement schedule and sampling design a ladder must have run on.

    ``timesteps``
        The **exact** timesteps, coarsest first.  Not a minimum count and not a
        rule about them: the actual numbers, hashed.
    ``ratio``
        The refinement ratio each adjacent pair must realize — two, for the
        paired halving this package uses.  Enforced *here*, at construction, so
        a non-halving schedule cannot be frozen in the first place.
    ``sampling``
        The :class:`SamplingDesign` above.
    """

    timesteps: Tuple[float, ...]
    ratio: float
    sampling: SamplingDesign

    def __post_init__(self):
        require_exact_type(self.sampling, SamplingDesign, "sampling")
        ratio = _require_positive(self.ratio, "ratio")
        if ratio <= 1.0:
            raise ValueError(f"the refinement ratio must exceed one, got "
                             f"{ratio!r}; a ladder whose steps do not shrink "
                             f"is not a refinement")
        try:
            listed = [require_finite(value, f"timesteps[{position}]")
                      for position, value in enumerate(self.timesteps)]
        except TypeError:
            raise TypeError(f"timesteps must be a sequence, got "
                            f"{type(self.timesteps).__name__}") from None
        if len(listed) < MINIMUM_REFINEMENT_LEVELS:
            raise ValueError(
                f"a refinement contract needs at least "
                f"{MINIMUM_REFINEMENT_LEVELS} levels, got {len(listed)}.  The "
                f"plan says 'at least three paired refinement levels', and "
                f"this is that sentence in the one place a caller cannot reach "
                f"past")
        for value in listed:
            if value <= 0.0:
                raise ValueError(f"every timestep must be strictly positive, "
                                 f"got {value!r}")
        for coarse, fine in zip(listed, listed[1:]):
            if fine >= coarse:
                raise ValueError(
                    f"timesteps must be given coarsest first and must strictly "
                    f"shrink: {coarse!r} is followed by {fine!r}")
            realized = coarse / fine
            slack = RATIO_ULPS * abs(np.spacing(ratio))
            if abs(realized - ratio) > slack:
                raise ValueError(
                    f"the step from {coarse!r} to {fine!r} refines by "
                    f"{realized!r}, not the declared ratio {ratio!r}.  A "
                    f"schedule of 0.040, 0.039, 0.038 is three levels and is "
                    f"not a refinement ladder")
        object.__setattr__(self, "timesteps", tuple(listed))
        object.__setattr__(self, "ratio", ratio)

    @property
    def levels(self) -> int:
        """How many refinement levels this contract declares."""
        return len(self.timesteps)

    @property
    def canonical(self) -> str:
        """The exact text this contract contributes to a digest."""
        return json.dumps([[repr(value) for value in self.timesteps],
                           repr(self.ratio), self.sampling.canonical],
                          separators=(",", ":"))


@dataclass(frozen=True)
class ValidationBudget:
    """The caps and the estimator one observable, at one position, is held to.

    ``observable`` / ``position``
        The two halves of the identity.  A position is a separate frozen ladder
        rather than a row inside one: several starts, several regimes, and
        every one of them has to pass, so collapsing them to a per-timestep
        worst case would compare a different path at each level.
    ``unit``
        One of :data:`UNITS`.  Hashed, and required to match nothing else —
        it is what makes an exported envelope a maximum over one quantity.
    ``estimator`` / ``quantile``
        **How the measured value is computed from the frozen data.**  One of
        :data:`ESTIMATORS`, and it lives in the budget rather than beside the
        data so that a caller cannot switch a mean for a quantile after seeing
        which one passes.  ``quantile`` is read only by the ``"quantile"``
        estimator.
    ``noise_floor``
        A **deterministic** resolution floor, in the observable's own unit: the
        granularity below which two levels' errors are not distinguishable for
        reasons that are not sampling.  It only ever excuses a single reversal,
        and the frozen design's paired bootstrap decides whether that reversal
        is unresolved.
    ``require_decrease``
        Whether this observable's error must fall under refinement.  False is a
        narrow declared case: an observable whose *endpoint* value diverges as
        the timestep shrinks cannot have a converging error.  The caps still
        bind, it is judged on its own, and it can never rescue a convergent
        observable.
    """

    observable: str
    unit: str
    absolute: float
    relative: float
    noise_floor: float = 0.0
    require_decrease: bool = True
    position: str = ""
    estimator: str = "mean"
    quantile: float = 0.5

    def __post_init__(self):
        observable = _require_label(self.observable, "observable")
        if self.unit not in UNITS:
            raise ValueError(f"unit must be one of {list(UNITS)}, got "
                             f"{self.unit!r}")
        if self.estimator not in ESTIMATORS:
            raise ValueError(f"estimator must be one of {list(ESTIMATORS)}, "
                             f"got {self.estimator!r}")
        _require_printable(self.position, "position", allow_empty=True)
        absolute = _require_positive(self.absolute, "absolute")
        relative = _require_positive(self.relative, "relative")
        floor = require_nonnegative(self.noise_floor, "noise_floor")
        quantile = require_finite(self.quantile, "quantile")
        if not 0.0 < quantile < 1.0:
            raise ValueError(f"quantile must lie strictly inside (0, 1), got "
                             f"{quantile!r}")
        if not isinstance(self.require_decrease, bool):
            raise TypeError(f"require_decrease must be a bool, got "
                            f"{type(self.require_decrease).__name__}")
        if floor > absolute:
            raise ValueError(
                f"noise_floor {floor!r} exceeds the absolute cap {absolute!r}; "
                f"a floor above the cap would excuse every difference the cap "
                f"exists to refuse")
        if floor > 0.0 and not self.require_decrease:
            raise ValueError(
                f"observable {observable!r} both waives the convergence "
                f"requirement and declares a resolution floor of {floor!r}.  "
                f"The floor only ever excuses an unresolved reversal; with the "
                f"requirement already lifted it excuses nothing")
        object.__setattr__(self, "observable", observable)
        object.__setattr__(self, "absolute", absolute)
        object.__setattr__(self, "relative", relative)
        object.__setattr__(self, "noise_floor", floor)
        object.__setattr__(self, "quantile", quantile)

    @property
    def identity(self) -> Tuple[str, str]:
        """``(observable, position)`` — what a budget and a sample pair on."""
        return (self.observable, self.position)

    @property
    def canonical(self) -> str:
        """The exact text this budget contributes to a digest."""
        return json.dumps([self.observable, self.position, self.unit,
                           repr(self.absolute), repr(self.relative),
                           repr(self.noise_floor), repr(self.require_decrease),
                           self.estimator, repr(self.quantile)],
                          separators=(",", ":"))


@dataclass(frozen=True)
class FrozenBudgets:
    """A whole validation's contract: schedule, sampling design, and caps.

    The digest is the mechanism.  A caller records ``FrozenBudgets(...).digest``
    before any result is computed and passes that string back to
    :func:`compare_refinement`; anything widened after seeing a discrepancy
    hashes differently and the comparison refuses to produce a verdict at all.

    Read the guarantee as "a moved threshold cannot be silent", never as "a
    threshold cannot be moved".

    One budget per ``(observable, position)`` pair, and duplicates are refused.
    """

    label: str
    contract: RefinementContract
    budgets: Tuple[ValidationBudget, ...]

    def __post_init__(self):
        label = _require_label(self.label, "label")
        require_exact_type(self.contract, RefinementContract, "contract")
        if isinstance(self.budgets, ValidationBudget):
            raise TypeError("budgets must be a sequence of ValidationBudget "
                            "objects, not a single ValidationBudget")
        try:
            listed = list(self.budgets)
        except TypeError:
            raise TypeError(f"budgets must be a sequence, got "
                            f"{type(self.budgets).__name__}") from None
        if not listed:
            raise ValueError("budgets must not be empty: a validation with no "
                             "tolerance is not a validation")
        rebuilt = []
        for position, entry in enumerate(listed):
            require_exact_type(entry, ValidationBudget, f"budgets[{position}]")
            rebuilt.append(ValidationBudget(
                entry.observable, entry.unit, entry.absolute, entry.relative,
                entry.noise_floor, entry.require_decrease, entry.position,
                entry.estimator, entry.quantile))
        contract = RefinementContract(
            self.contract.timesteps, self.contract.ratio,
            SamplingDesign(self.contract.sampling.unit,
                           self.contract.sampling.clusters,
                           self.contract.sampling.replications,
                           self.contract.sampling.method,
                           self.contract.sampling.resamples,
                           self.contract.sampling.coverage,
                           self.contract.sampling.seed))
        names = [entry.identity for entry in rebuilt]
        duplicates = sorted({n for n in names if names.count(n) > 1})
        if duplicates:
            raise ValueError(
                f"two budgets name the same observable and position "
                f"{duplicates}; a comparison offered a choice of tolerances is "
                f"not held to either of them")
        object.__setattr__(self, "label", label)
        object.__setattr__(self, "contract", contract)
        object.__setattr__(self, "budgets", tuple(rebuilt))

    @property
    def canonical(self) -> str:
        """The exact text the digest is taken over."""
        body = ",".join(entry.canonical
                        for entry in sorted(self.budgets,
                                            key=lambda e: e.identity))
        return (f"{VALIDATION_SCHEMA}|{self.label}|"
                f"{self.contract.canonical}|[{body}]")

    @property
    def digest(self) -> str:
        """SHA-256 of :attr:`canonical`, as 64 lowercase hexadecimal digits."""
        return _digest_text(self.canonical)

    def budget_for(self, observable, position="") -> ValidationBudget:
        """The budget for one identity, or a ``KeyError`` naming what exists."""
        _require_printable(position, "position", allow_empty=True)
        wanted = (_require_label(observable, "observable"), position)
        for entry in self.budgets:
            if entry.identity == wanted:
                return entry
        raise KeyError(f"no budget for {wanted!r}; this validation froze "
                       f"{sorted(e.identity for e in self.budgets)}")


#: The canonical dtype every stored observation is hashed and reconstructed
#: under: little-endian binary64.  Fixed rather than native, so a digest means
#: the same thing on a big-endian machine as on this one.
CANONICAL_DTYPE = "<f8"


def _sealed(array):
    """Return a read-only ndarray whose authority is an **immutable** object.

    The backing is a ``bytes``.  ``np.frombuffer`` over it yields a read-only
    array whose ``base`` is that ``bytes``; nothing in the chain owns writable
    memory, so no traversal reaches any.

    NumPy write flags alone are not enough, and an independent review showed
    exactly why.  The previous version returned a read-only *view* of a
    read-only owned array, and the owner was publicly reachable::

        sample.values.base.setflags(write=True)
        sample.values.base[...] = replacement

    NumPy permits that, because the owner owns its memory.  Zeros written that
    way turned a genuine ``numerical_no_result`` into a ``pass`` and changed the
    dataset digest, and a mapping subclass whose ``items()`` performed the write
    during the parser's payload serialization got the forged pass accepted.

    ``bytes`` closes the route at the root rather than at each flag: there is
    no writable owner anywhere in the object graph to find.
    """
    canonical = np.ascontiguousarray(array, dtype=CANONICAL_DTYPE)
    payload = canonical.tobytes()
    flat = np.frombuffer(payload, dtype=CANONICAL_DTYPE)
    return flat.reshape(canonical.shape)


def _require_names(values, count, name) -> Tuple[str, ...]:
    """Return ``values`` as exactly ``count`` unique non-empty labels.

    Cluster and member identities are what turn an array of numbers into a
    record of *which physical things* were measured.  Without them a sample of
    twenty-four identical rows is indistinguishable from one trial copied
    twenty-four times — an independent review's example — because identical
    rows are perfectly legitimate for a binary observable, so the observations
    cannot be their own identity.

    Uniqueness is required, not merely encouraged: a duplicated master-trial
    identifier is either a bug or a sample claiming more independent evidence
    than it holds, and a cluster bootstrap over it reports an error bar for a
    sample size that does not exist.
    """
    if isinstance(values, str):
        raise TypeError(f"{name} must be a sequence of strings, not a single "
                        f"string")
    try:
        listed = list(values)
    except TypeError:
        raise TypeError(f"{name} must be a sequence, got "
                        f"{type(values).__name__}") from None
    if len(listed) != count:
        raise ValueError(
            f"{name} holds {len(listed)} labels for {count} rows.  Every "
            f"cluster and every member carries an explicit identity; a sample "
            f"that cannot say which physical trial each row came from cannot "
            f"reject a substituted or duplicated one")
    for position, entry in enumerate(listed):
        _require_printable(entry, f"{name}[{position}]")
    duplicates = sorted({e for e in listed if listed.count(e) > 1})
    if duplicates:
        raise ValueError(
            f"{name} repeats {duplicates}; a duplicated identity is either a "
            f"bug or a sample claiming more independent evidence than it "
            f"holds, and the bootstrap would report an error bar for a sample "
            f"size that does not exist")
    return tuple(listed)


def _encode_labels(field, labels) -> bytes:
    """A canonical injective encoding of one named label sequence.

    ``field name | count | (length, utf-8 bytes) per label``, every integer
    fixed-width big-endian.  Injective: two different sequences cannot produce
    the same bytes, because every boundary is carried by an explicit length
    rather than by a delimiter that could occur inside a label.

    The previous encoding joined labels with NUL.  An independent review found
    the consequence — ``("a\x00b", "c")`` and ``("a", "b\x00c")`` serialize to
    the same bytes and therefore to the same sample **and** dataset digest,
    with identical observations.  SHA-256 had not collided; two different
    physical manifests had been flattened to one byte string before hashing.

    Encoding is UTF-8 with no normalization: two labels that differ as Python
    strings are two different labels here, which is the conservative reading
    and the one the identity checks already assume.
    """
    parts = [field.encode("utf-8")]
    out = bytearray()
    out += len(parts[0]).to_bytes(4, "big")
    out += parts[0]
    out += len(labels).to_bytes(4, "big")
    for label in labels:
        encoded = label.encode("utf-8")
        out += len(encoded).to_bytes(4, "big")
        out += encoded
    return bytes(out)


def _estimate(kind, quantile, block) -> float:
    """The declared estimator applied to one block of frozen observations.

    ``block`` is ``(clusters, members)`` and is pooled before the statistic is
    taken, so a cluster with more members contributes more observations —
    which is what a pooled proportion or a pooled quantile means.  The
    *resampling* is over clusters, so the dependence inside one is respected
    even though the statistic does not weight by it.
    """
    if kind not in ESTIMATORS:
        raise ValueError(f"kind must be one of {list(ESTIMATORS)}, got "
                         f"{kind!r}; an unnamed estimator would fall through to "
                         f"whichever branch happened to be last, which is how a "
                         f"mean becomes a quantile without anybody choosing")
    flat = np.asarray(block, dtype=float).reshape(-1)
    if kind == "mean":
        return float(np.mean(flat))
    return float(np.quantile(flat, quantile))


@dataclass(frozen=True)
class PairedSample:
    """The cluster-level observations one ladder's numbers are computed from.

    ``values``
        ``(clusters, members, levels)``.  One row per independent cluster, one
        column per observation inside it, one page per refinement level — the
        *same* clusters at every level, which is what makes the ladder paired
        and what lets a bootstrap resample of clusters be evaluated at every
        level at once.
    ``baseline``
        The same shape in clusters and levels, possibly with a different member
        count, or ``None``.  When present the measured value is
        ``estimator(values) - estimator(baseline)``, which is how an
        audited-minus-primary shift is expressed; when absent the measured
        value is ``estimator(values)`` and ``reference`` is what it is compared
        against.
    ``reference``
        The oracle value, or zero for a paired shift.

    **This is the frozen input, and it is the only input.**  There is no
    standard-error argument anywhere in this module's comparison path, because
    an independent review showed what one costs: the same ladder returned
    `numerical_no_result` and `pass` depending on a scalar the caller chose.
    The digest below binds the actual observations, so a report can be
    reproduced rather than believed.
    """

    observable: str
    position: str
    values: np.ndarray
    baseline: Optional[np.ndarray] = None
    reference: float = 0.0
    identities: Tuple[str, ...] = ()
    members: Tuple[str, ...] = ()
    baseline_members: Tuple[str, ...] = ()

    def __post_init__(self):
        observable = _require_label(self.observable, "observable")
        _require_printable(self.position, "position", allow_empty=True)
        # A defensive, owned, C-contiguous copy -- never ``np.asarray``, which
        # keeps the caller's array when the dtype already matches.  An
        # independent review turned that aliasing into a verdict bypass: the
        # digest is a property over the current bytes, so mutating the caller's
        # array between ``parse_validation_report``'s digest check and its
        # recomputation turned a genuine no-result into a parser-approved pass.
        # ``frozen=True`` freezes field assignment, not the evidence behind it.
        values = np.array(self.values, dtype=float, order="C", copy=True)
        if values.ndim != 3:
            raise ValueError(f"values must be (clusters, members, levels), "
                             f"got shape {values.shape}")
        if values.shape[0] < 2 or values.shape[1] < 1 \
                or values.shape[2] < MINIMUM_REFINEMENT_LEVELS:
            raise ValueError(
                f"values has shape {values.shape}; a paired sample needs at "
                f"least two clusters, one member and "
                f"{MINIMUM_REFINEMENT_LEVELS} levels")
        if not np.all(np.isfinite(values)):
            raise ValueError("values must be finite; an unresolved run is "
                             "censored at a declared horizon by the caller, "
                             "not carried here as an infinity")
        baseline = self.baseline
        if baseline is not None:
            baseline = np.array(baseline, dtype=float, order="C", copy=True)
            if baseline.ndim != 3:
                raise ValueError(f"baseline must be (clusters, members, "
                                 f"levels), got shape {baseline.shape}")
            if (baseline.shape[0], baseline.shape[2]) != (values.shape[0],
                                                          values.shape[2]):
                raise ValueError(
                    f"baseline {baseline.shape} and values {values.shape} must "
                    f"share their cluster and level axes; the two arms are "
                    f"paired on the cluster, and a bootstrap that resampled "
                    f"different clusters in each arm would not be paired at "
                    f"all")
            if not np.all(np.isfinite(baseline)):
                raise ValueError("baseline must be finite")
        identities = _require_names(self.identities, values.shape[0],
                                    "identities")
        members = _require_names(self.members, values.shape[1], "members")
        if baseline is None:
            if self.baseline_members:
                raise ValueError("baseline_members were given without a "
                                 "baseline arm")
            baseline_members = ()
        else:
            baseline_members = _require_names(self.baseline_members,
                                              baseline.shape[1],
                                              "baseline_members")
        # Read-only, so a caller holding the returned object cannot write
        # through it either.  Both directions of the alias are closed: the
        # copy above stops the source from reaching in, this stops the sample
        # from reaching back out.
        values = _sealed(values)
        if baseline is not None:
            baseline = _sealed(baseline)
        object.__setattr__(self, "observable", observable)
        object.__setattr__(self, "values", values)
        object.__setattr__(self, "baseline", baseline)
        object.__setattr__(self, "identities", identities)
        object.__setattr__(self, "members", members)
        object.__setattr__(self, "baseline_members", baseline_members)
        object.__setattr__(self, "reference",
                           require_finite(self.reference, "reference"))

    @property
    def clusters(self) -> int:
        """How many independent clusters this sample holds."""
        return int(self.values.shape[0])

    @property
    def levels(self) -> int:
        """How many refinement levels this sample holds."""
        return int(self.values.shape[2])

    @property
    def digest(self) -> str:
        """SHA-256 over the identity, the shapes and the observations.

        Byte-exact over the data, so two reports carrying the same digest were
        computed from the same numbers in the same cluster order.  Reordering
        the clusters changes it, which is deliberate: cluster membership and
        order are part of what a bootstrap is reproducible against.
        """
        parts = [_encode_labels("observable", (self.observable,)),
                 _encode_labels("position", (self.position,)),
                 _encode_labels("reference", (repr(self.reference),)),
                 _encode_labels("identities", self.identities),
                 _encode_labels("members", self.members),
                 _encode_labels("baseline_members", self.baseline_members),
                 _encode_labels("values.dtype", (CANONICAL_DTYPE,)),
                 _encode_labels("values.role", ("measured arm",)),
                 _encode_labels("values.shape",
                                tuple(str(n) for n in self.values.shape)),
                 np.ascontiguousarray(self.values,
                                      dtype=CANONICAL_DTYPE).tobytes()]
        if self.baseline is None:
            parts.append(_encode_labels("baseline", ("absent",)))
        else:
            parts.append(_encode_labels("baseline.dtype",
                                        (CANONICAL_DTYPE,)))
            parts.append(_encode_labels("baseline.role", ("reference arm",)))
            parts.append(_encode_labels("baseline.shape",
                                        tuple(str(n)
                                              for n in self.baseline.shape)))
            parts.append(np.ascontiguousarray(
                self.baseline, dtype=CANONICAL_DTYPE).tobytes())
        digest = hashlib.sha256()
        for part in parts:
            digest.update(len(part).to_bytes(8, "big"))
            digest.update(part)
        return digest.hexdigest()

    def measured(self, kind, quantile, rows=None) -> np.ndarray:
        """The **signed** measured value at each level, on one resample.

        ``rows`` is an array of cluster indices drawn with replacement, or
        ``None`` for the whole sample.  The *same* rows are used at every
        level, which is what makes the level-to-level contrast paired.

        Signed because that is what the estimator produces; the *bootstrap*
        does not resample this.  :func:`compare_refinement` resamples
        :meth:`errors`, the folded absolute error, because that is the
        statistic the gate consumes.

        An earlier revision did resample this signed quantity and argued from
        ``| |a| - |b| | <= |a - b|`` that the result was conservative.  An
        independent review showed the argument is invalid — a pointwise
        inequality between two random variables says nothing about their
        standard deviations — with a two-cluster sample whose signed adjacent
        contrast has standard error exactly zero while the folded contrast has
        0.699.  The argument is gone with the practice.
        """
        take = slice(None) if rows is None else rows
        values = self.values[take]
        baseline = None if self.baseline is None else self.baseline[take]
        out = np.empty(self.levels, dtype=float)
        for level in range(self.levels):
            value = _estimate(kind, quantile, values[:, :, level])
            if baseline is not None:
                value -= _estimate(kind, quantile, baseline[:, :, level])
            out[level] = value
        return out

    def errors(self, kind, quantile, rows=None) -> np.ndarray:
        """``|measured - reference|`` at each level, on one resample."""
        return np.abs(self.measured(kind, quantile, rows) - self.reference)


@dataclass(frozen=True)
class ValidationDataset:
    """Every frozen sample a comparison is computed from, plus its contract.

    Hashed, and the digest travels into the report, so that a later reader can
    ask whether two reports describe the same measurements rather than
    inferring it from the numbers.
    """

    label: str
    contract: RefinementContract
    samples: Tuple[PairedSample, ...]

    def __post_init__(self):
        label = _require_label(self.label, "label")
        require_exact_type(self.contract, RefinementContract, "contract")
        if isinstance(self.samples, PairedSample):
            raise TypeError("samples must be a sequence of PairedSample "
                            "objects, not a single PairedSample")
        try:
            listed = list(self.samples)
        except TypeError:
            raise TypeError(f"samples must be a sequence, got "
                            f"{type(self.samples).__name__}") from None
        if not listed:
            raise ValueError("samples must not be empty")
        rebuilt = []
        for position, entry in enumerate(listed):
            require_exact_type(entry, PairedSample, f"samples[{position}]")
            # Rebuilt, not retained: reconstruction runs the defensive copy
            # above, so a dataset never holds a sample somebody else can still
            # write to, however it was assembled.
            rebuilt.append(PairedSample(entry.observable, entry.position,
                                        entry.values, entry.baseline,
                                        entry.reference, entry.identities,
                                        entry.members, entry.baseline_members))
        listed = rebuilt
        names = [entry.observable + "\x00" + entry.position
                 for entry in listed]
        duplicates = sorted({n for n in names if names.count(n) > 1})
        if duplicates:
            raise ValueError(f"two samples name one identity "
                             f"{[n.replace(chr(0), '/') for n in duplicates]}")
        # Every sample under one contract measures the *same ordered physical
        # clusters*.  Without this, two observables could be bootstrapped over
        # different master trials while claiming one sampling design, and a
        # cross-observable substitution would be invisible.  Order is
        # significant and is not canonicalized: a reordered sample is a
        # different dataset, which is the conservative reading and the one the
        # digest already implements.
        expected = listed[0].identities
        for entry in listed[1:]:
            if entry.identities != expected:
                raise ValueError(
                    f"sample {(entry.observable, entry.position)!r} is "
                    f"measured over cluster identities that differ from "
                    f"{(listed[0].observable, listed[0].position)!r}'s.  Every "
                    f"sample governed by one sampling design has to be the "
                    f"same ordered physical clusters, or a bootstrap that "
                    f"resamples an index is resampling two different things")
        object.__setattr__(self, "label", label)
        object.__setattr__(self, "samples", tuple(listed))

    @property
    def digest(self) -> str:
        """SHA-256 over the contract and every sample digest, in order."""
        body = "|".join(entry.digest
                        for entry in sorted(self.samples,
                                            key=lambda e: (e.observable,
                                                           e.position)))
        clusters = _encode_labels("clusters",
                                  self.samples[0].identities).hex()
        return _digest_text(f"{VALIDATION_SCHEMA}|{self.label}|"
                            f"{self.contract.canonical}|{clusters}|{body}")

    def sample_for(self, observable, position="") -> PairedSample:
        """The sample for one identity, or a ``KeyError`` naming what exists."""
        _require_printable(position, "position", allow_empty=True)
        wanted = (_require_label(observable, "observable"), position)
        for entry in self.samples:
            if (entry.observable, entry.position) == wanted:
                return entry
        raise KeyError(f"no sample for {wanted!r}; this dataset holds "
                       f"{sorted((e.observable, e.position) for e in self.samples)}")


@dataclass(frozen=True)
class RefinementLevel:
    """One observable, at one position, at one timestep — **as computed here**.

    This is an *output* record.  It used to be the caller's input, carrying a
    caller-supplied standard error that decided verdicts; it is now built by
    :func:`compare_refinement` from the frozen data, and
    ``compare_refinement`` refuses a sequence of these outright so that the old
    path cannot be reached by accident.

    ``standard_error``
        The cluster-bootstrap standard error of *this level's error*.
    ``paired_error``
        The cluster-bootstrap standard error of the error's change from the
        next coarser level, computed on the same resamples so that the
        dependence between the two levels is carried rather than assumed away.
        Zero at the coarsest level, which has no predecessor.
    """

    observable: str
    unit: str
    timestep: float
    measured: float
    reference: float
    standard_error: float = 0.0
    clusters: int = 2
    position: str = ""
    paired_error: float = 0.0
    span_error: float = 0.0

    def __post_init__(self):
        observable = _require_label(self.observable, "observable")
        if self.unit not in UNITS:
            raise ValueError(f"unit must be one of {list(UNITS)}, got "
                             f"{self.unit!r}")
        _require_printable(self.position, "position", allow_empty=True)
        timestep = _require_positive(self.timestep, "timestep")
        measured = require_finite(self.measured, "measured")
        reference = require_finite(self.reference, "reference")
        error = require_nonnegative(self.standard_error, "standard_error")
        paired = require_nonnegative(self.paired_error, "paired_error")
        span = require_nonnegative(self.span_error, "span_error")
        clusters = require_positive_int(self.clusters, "clusters", minimum=2)
        object.__setattr__(self, "observable", observable)
        object.__setattr__(self, "timestep", timestep)
        object.__setattr__(self, "measured", measured)
        object.__setattr__(self, "reference", reference)
        object.__setattr__(self, "standard_error", error)
        object.__setattr__(self, "paired_error", paired)
        object.__setattr__(self, "span_error", span)
        object.__setattr__(self, "clusters", clusters)

    @property
    def identity(self) -> Tuple[str, str]:
        """``(observable, position)``."""
        return (self.observable, self.position)

    @property
    def absolute_error(self) -> float:
        """``|measured - reference|``, in this level's own unit."""
        return abs(self.measured - self.reference)

    @property
    def relative_error(self) -> float:
        """The absolute error over ``|reference|``, or the absolute error.

        A zero reference has no scale to be relative to, so the absolute
        difference is reported unchanged.  The absolute cap is the binding one
        there, which is why both caps exist.
        """
        scale = abs(self.reference)
        return self.absolute_error if scale == 0.0 else self.absolute_error / scale

    def bounded_error(self, coverage) -> float:
        """The uncertainty-aware upper bound on this level's error.

        ``|error| + coverage * standard_error``.  This, and not the point
        estimate, is what the absolute cap is applied to at the finest level.

        The distinction is a finding.  A point estimate of 0.0112 against a cap
        of 0.10 reads like a comfortable pass; with a bootstrap standard error
        of 0.0829 on the same number, two standard errors reach 0.177 and the
        honest statement is that the sample does not resolve the question.  A
        cap applied to a point estimate says where the estimate landed, not
        where the quantity is.
        """
        return self.absolute_error + coverage * self.standard_error


@dataclass(frozen=True)
class ValidationVerdict:
    """The machine-readable outcome of one refinement comparison.

    ``verdict``
        ``"pass"`` or ``"numerical_no_result"``.  Nothing else, and in
        particular nothing that could be read as a physical finding.
    ``reasons``
        Empty on a pass.  On a no-result, one sentence per failed requirement.
    ``levels``
        The derived per-level records.  Non-empty on a pass; a no-result may
        carry none, because a sample whose cluster count does not match the
        frozen design is refused *before* any statistic is taken and there is
        nothing to report about it.
    ``dataset_digest``
        The digest of the frozen observations the numbers came from, beside the
        digest of the contract they were judged against.  Two digests, because
        they answer two different questions: whether the thresholds were frozen
        before the results, and whether two reports describe the same
        measurements.
    """

    label: str
    digest: str
    verdict: str
    levels: Tuple[RefinementLevel, ...]
    reasons: Tuple[str, ...]
    dataset_digest: str = ""
    blocking: Tuple[Tuple[str, str, str], ...] = ()
    budgets: Optional[FrozenBudgets] = None

    def __post_init__(self):
        label = _require_label(self.label, "label")
        digest = _require_label(self.digest, "digest")
        if not isinstance(self.dataset_digest, str):
            raise TypeError(f"dataset_digest must be a string, got "
                            f"{type(self.dataset_digest).__name__}")
        if self.verdict not in VERDICTS:
            raise ValueError(f"verdict must be one of {list(VERDICTS)}, got "
                             f"{self.verdict!r}")
        listed = tuple(self.levels)
        if not listed and self.verdict == "pass":
            raise ValueError(
                "a passing verdict needs the levels it was reached from; a "
                "comparison that measured nothing has not passed anything")
        for position, entry in enumerate(listed):
            require_exact_type(entry, RefinementLevel, f"levels[{position}]")
        reasons = tuple(self.reasons)
        for position, reason in enumerate(reasons):
            _require_reason(reason, f"reasons[{position}]")
        blocking = tuple(tuple(entry) for entry in self.blocking)
        for entry in blocking:
            if len(entry) != 3 or entry[2] not in REASON_CODES:
                raise ValueError(f"blocking entries are (observable, "
                                 f"position, code) with a code from "
                                 f"{list(REASON_CODES)}, got {entry!r}")
        if len(blocking) != len(reasons):
            raise ValueError(
                f"{len(blocking)} blocking codes beside {len(reasons)} "
                f"reasons; they are one list written twice, machine-readably "
                f"and in words")
        object.__setattr__(self, "blocking", blocking)
        if self.budgets is not None:
            require_exact_type(self.budgets, FrozenBudgets, "budgets")
        if (self.verdict == "pass") != (not reasons):
            raise ValueError(
                f"verdict {self.verdict!r} and {len(reasons)} reasons "
                f"disagree: a pass carries no reason and a no-result must say "
                f"what blocked it")
        object.__setattr__(self, "label", label)
        object.__setattr__(self, "digest", digest)
        object.__setattr__(self, "levels", listed)
        object.__setattr__(self, "reasons", reasons)

    @property
    def passed(self) -> bool:
        """Did the comparison meet every frozen budget?"""
        return self.verdict == "pass"

    @property
    def finest(self) -> Tuple[RefinementLevel, ...]:
        """Each ``(observable, position)`` ladder's finest-timestep level."""
        best = {}
        for level in self.levels:
            current = best.get(level.identity)
            if current is None or level.timestep < current.timestep:
                best[level.identity] = level
        return tuple(best[key] for key in sorted(best))

    def envelope_for(self, unit) -> Optional[dict]:
        """The measured refinement-error envelope **within one unit**.

        ``None`` if this verdict measured nothing in that unit.  The identity
        is part of the answer, not a footnote: a later production budget needs
        to know that the largest error was, say, the exit-time quantile at the
        band's centre.  There is deliberately no single scalar across units.
        """
        if unit not in UNITS:
            raise ValueError(f"unit must be one of {list(UNITS)}, got "
                             f"{unit!r}")
        rows = [level for level in self.finest if level.unit == unit]
        if not rows:
            return None
        worst_absolute = max(rows, key=lambda level: level.absolute_error)
        worst_relative = max(rows, key=lambda level: level.relative_error)
        return {"unit": unit,
                "absolute": worst_absolute.absolute_error,
                "absolute_observable": worst_absolute.observable,
                "absolute_position": worst_absolute.position,
                "standard_error": worst_absolute.standard_error,
                "relative": worst_relative.relative_error,
                "relative_observable": worst_relative.observable,
                "relative_position": worst_relative.position,
                "timestep": worst_absolute.timestep}

    @property
    def envelopes(self) -> dict:
        """:meth:`envelope_for` over every unit this verdict measured."""
        return {unit: self.envelope_for(unit) for unit in UNITS
                if self.envelope_for(unit) is not None}

    def as_dict(self) -> dict:
        """A JSON-serializable report, for a later ticket to read.

        Keys are stable under :data:`VALIDATION_SCHEMA`, which is **v4**.  A
        change to what any of them means is a change of schema version rather
        than a change of numbers under the same names, and the tag has moved
        three times for exactly that reason -- see :data:`VALIDATION_SCHEMA`
        for what each step changed.  A consumer written against any earlier
        tag must refuse this shape, and :func:`parse_validation_report` is the
        parser that shows what refusing looks like.
        """
        if self.budgets is None:
            raise ValueError(
                "this verdict carries no frozen budgets, so a reader could "
                "not recompute it; only a verdict built by "
                "compare_refinement can be serialized")
        contract = self.budgets.contract
        sampling = contract.sampling
        return {
            "schema": VALIDATION_SCHEMA,
            "label": self.label,
            "budget_digest": self.digest,
            "dataset_digest": self.dataset_digest,
            "verdict": self.verdict,
            "passed": self.passed,
            "reasons": list(self.reasons),
            "blocking": [{"observable": observable, "position": position,
                          "code": code}
                         for observable, position, code in self.blocking],
            "contract": {
                "timesteps": list(contract.timesteps),
                "ratio": contract.ratio,
                "sampling": {"unit": sampling.unit,
                             "clusters": sampling.clusters,
                             "replications": sampling.replications,
                             "method": sampling.method,
                             "resamples": sampling.resamples,
                             "coverage": sampling.coverage,
                             "seed": sampling.seed}},
            "budgets": [{"observable": entry.observable,
                         "position": entry.position,
                         "unit": entry.unit,
                         "absolute": entry.absolute,
                         "relative": entry.relative,
                         "noise_floor": entry.noise_floor,
                         "require_decrease": entry.require_decrease,
                         "estimator": entry.estimator,
                         "quantile": entry.quantile}
                        for entry in self.budgets.budgets],
            "envelopes": self.envelopes,
            "levels": [
                {"observable": level.observable,
                 "position": level.position,
                 "unit": level.unit,
                 "timestep": level.timestep,
                 "measured": level.measured,
                 "reference": level.reference,
                 "absolute_error": level.absolute_error,
                 "relative_error": level.relative_error,
                 "standard_error": level.standard_error,
                 "paired_error": level.paired_error,
                 "span_error": level.span_error,
                 "clusters": level.clusters}
                for level in self.levels],
        }


#: The exact top-level keys a report under the live :data:`VALIDATION_SCHEMA`
#: carries.  Closed: an unknown key is a report from some other producer, and
#: reading it as this one is how a downstream consumer ends up trusting a field
#: nobody defined.
REPORT_FIELDS = ("schema", "label", "budget_digest", "dataset_digest",
                 "verdict", "passed", "reasons", "blocking", "contract",
                 "budgets", "envelopes", "levels")

#: The exact keys every level row carries.
LEVEL_FIELDS = ("observable", "position", "unit", "timestep", "measured",
                "reference", "absolute_error", "relative_error",
                "standard_error", "paired_error", "span_error", "clusters")

#: The exact keys every frozen budget row carries in a report.
BUDGET_FIELDS = ("observable", "position", "unit", "absolute", "relative",
                 "noise_floor", "require_decrease", "estimator", "quantile")

#: The exact keys the frozen contract block carries.
CONTRACT_FIELDS = ("timesteps", "ratio", "sampling")
SAMPLING_FIELDS = ("unit", "clusters", "replications", "method", "resamples",
                   "coverage", "seed")

#: The exact keys every blocking entry carries.
BLOCKING_FIELDS = ("observable", "position", "code")

#: The exact keys every per-unit envelope carries.
ENVELOPE_FIELDS = ("unit", "absolute", "absolute_observable",
                   "absolute_position", "standard_error", "relative",
                   "relative_observable", "relative_position", "timestep")

#: How long a digest is, in hexadecimal characters.
DIGEST_LENGTH = 64
_HEX = frozenset("0123456789abcdef")


def _require_digest_field(value, name):
    """A digest is 64 lowercase hexadecimal characters, or it is not one."""
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, str):
        raise ValueError(f"{name} must be a string, got "
                         f"{type(value).__name__} {value!r}")
    if len(value) != DIGEST_LENGTH or set(value) - _HEX:
        raise ValueError(
            f"{name} is {value!r}, not {DIGEST_LENGTH} lowercase hexadecimal "
            f"characters.  An empty or malformed digest is a report that "
            f"cannot be traced to the contract or the data it claims to "
            f"describe, and an independent review found the empty string "
            f"being accepted")
    return value


def _require_number(value, name, *, finite=True, minimum=None):
    """A JSON number, of the kind this field is declared to hold."""
    if isinstance(value, (bool, np.bool_)) or not isinstance(value,
                                                             (int, float)):
        raise ValueError(f"{name} must be a real number, got "
                         f"{type(value).__name__} {value!r}")
    number = float(value)
    if finite and not math.isfinite(number):
        raise ValueError(f"{name} must be finite, got {number!r}")
    if minimum is not None and number < minimum:
        raise ValueError(f"{name} must be at least {minimum}, got {number!r}")
    return number


def inspect_validation_report(payload) -> dict:
    """Read a validation report, refusing anything this schema does not mean.

    Returned unchanged on success.  **This is not the trust boundary and does
    not establish the verdict.**  It checks shape and local consistency: the
    schema tag, the closed key sets, the digests' form, the units, the finite
    and non-negative fields, the derived per-level arithmetic, the one-unit-per
    -ladder rule, the shared schedule, the coarsest paired error, the rebuilt
    envelopes, and that the declared blocking codes agree with what the
    embedded levels and budgets imply under :func:`_ladder_codes`.

    What it cannot do is check that the *levels themselves* describe the data.
    An independent review made that concrete: zeroing every level's standard,
    paired and span error, flipping the verdict to a pass, and emptying the
    reasons produced a payload this function accepts, with both digests
    byte-for-byte unchanged, because the uncertainty summaries are not bound to
    anything.  :func:`parse_validation_report` closes that by recomputing them
    from the frozen dataset, and a caller who wants a verdict they can act on
    must use it.

    An independent review found the previous version checking only the tag,
    the presence of top-level keys, verdict membership and the presence of
    level keys — and accepting, among others, a report with
    ``verdict = "numerical_no_result"`` beside ``passed = True``.  That is the
    concrete downstream bypass: a later ticket calls the advertised parser,
    trusts ``passed``, and advances past a scientifically blocked result.

    What is enforced now:

    * the schema tag is exactly :data:`VALIDATION_SCHEMA`, and the top-level
      key set is exactly :data:`REPORT_FIELDS` — a closed shape, so an unknown
      key is refused rather than ignored;
    * **``passed`` is exactly ``verdict == "pass"``**, and a pass carries no
      reasons while a no-result carries at least one;
    * both digests are 64 lowercase hexadecimal characters, and they are
      separate fields;
    * every level row has exactly :data:`LEVEL_FIELDS`, a string observable and
      position, a declared unit, a strictly positive timestep, finite
      ``measured`` and ``reference``, non-negative errors, at least two
      clusters, ``absolute_error`` equal to ``|measured - reference|`` and
      ``relative_error`` consistent with it;
    * no two level rows share an observable, position and timestep;
    * every ladder present has the same set of timesteps as every other, so a
      report cannot carry a truncated ladder for one identity; and
    * every envelope has exactly :data:`ENVELOPE_FIELDS`, a unit matching its
      key, and names an observable and position that the levels contain.

    The residual tolerance on the derived fields is a few ulps, because they
    are the producer's own floats round-tripped through JSON, not
    recomputations.
    """
    if not isinstance(payload, dict):
        raise TypeError(f"a validation report is a mapping, got "
                        f"{type(payload).__name__}")
    tag = payload.get("schema")
    if tag != VALIDATION_SCHEMA:
        raise ValueError(
            f"this report declares schema {tag!r}; this parser reads "
            f"{VALIDATION_SCHEMA!r} only.  The report shape changed "
            f"incompatibly between versions -- v1's single mixed-unit "
            f"'envelope' became v2's per-unit 'envelopes', and levels gained "
            f"positions, standard errors, paired errors and cluster counts -- "
            f"so a reader that accepted both tags would be interpreting one "
            f"shape under the other's meaning")
    extra = sorted(set(payload) - set(REPORT_FIELDS))
    missing = sorted(set(REPORT_FIELDS) - set(payload))
    if missing:
        raise ValueError(f"the report is missing {missing} and cannot be "
                         f"{VALIDATION_SCHEMA!r}")
    if extra:
        raise ValueError(
            f"the report carries {extra}, which {VALIDATION_SCHEMA!r} does not "
            f"define.  The shape is closed: 'envelope' is v1's mixed-unit "
            f"scalar, and any other unknown key is a report from a different "
            f"producer being read as this one")

    _require_label(payload["label"], "label")
    _require_digest_field(payload["budget_digest"], "budget_digest")
    _require_digest_field(payload["dataset_digest"], "dataset_digest")
    if payload["budget_digest"] == payload["dataset_digest"]:
        raise ValueError(
            "the budget digest and the dataset digest are the same string; "
            "they hash different objects -- the contract frozen before the "
            "results and the observations taken after -- and a report where "
            "they coincide is not carrying both")

    verdict = payload["verdict"]
    if verdict not in VERDICTS:
        raise ValueError(f"verdict {verdict!r} is not one of {list(VERDICTS)}")
    passed = payload["passed"]
    if not isinstance(passed, bool):
        raise ValueError(f"passed must be a bool, got "
                         f"{type(passed).__name__} {passed!r}")
    if passed != (verdict == "pass"):
        raise ValueError(
            f"the report says verdict={verdict!r} with passed={passed!r}.  "
            f"These are one fact written twice, and a consumer that trusts the "
            f"second would advance past a result the first blocks.  This is "
            f"the downstream bypass an independent review found: a "
            f"numerical_no_result carrying passed=True")
    reasons = payload["reasons"]
    if not isinstance(reasons, list):
        raise ValueError("reasons must be a list")
    for position, reason in enumerate(reasons):
        _require_reason(reason, f"reasons[{position}]")
    if bool(reasons) == passed:
        raise ValueError(
            f"the report is a {verdict!r} carrying {len(reasons)} reason(s); a "
            f"pass carries none and a no-result must say what blocked it")

    levels = payload["levels"]
    if not isinstance(levels, list):
        raise ValueError("levels must be a list")
    if not levels and passed:
        raise ValueError("a passing report carries no levels; a comparison "
                         "that measured nothing has not passed anything")
    seen = set()
    ladders = {}
    ladder_units = {}
    rows_by_identity = {}
    for position, level in enumerate(levels):
        where = f"levels[{position}]"
        if not isinstance(level, dict):
            raise ValueError(f"{where} must be a mapping")
        absent = sorted(set(LEVEL_FIELDS) - set(level))
        unknown = sorted(set(level) - set(LEVEL_FIELDS))
        if absent:
            raise ValueError(f"{where} is missing {absent}")
        if unknown:
            raise ValueError(f"{where} carries undefined {unknown}")
        _require_label(level["observable"], f"{where}.observable")
        _require_printable(level["position"], f"{where}.position",
                           allow_empty=True)
        if level["unit"] not in UNITS:
            raise ValueError(f"{where}.unit is {level['unit']!r}, not one of "
                             f"{list(UNITS)}")
        timestep = _require_number(level["timestep"], f"{where}.timestep")
        if timestep <= 0.0:
            raise ValueError(f"{where}.timestep must be strictly positive, "
                             f"got {timestep!r}")
        measured = _require_number(level["measured"], f"{where}.measured")
        reference = _require_number(level["reference"], f"{where}.reference")
        absolute = _require_number(level["absolute_error"],
                                   f"{where}.absolute_error", minimum=0.0)
        relative = _require_number(level["relative_error"],
                                   f"{where}.relative_error", minimum=0.0)
        _require_number(level["standard_error"], f"{where}.standard_error",
                        minimum=0.0)
        _require_number(level["paired_error"], f"{where}.paired_error",
                        minimum=0.0)
        _require_number(level["span_error"], f"{where}.span_error",
                        minimum=0.0)
        clusters = level["clusters"]
        if isinstance(clusters, bool) or not isinstance(clusters, int) \
                or clusters < 2:
            raise ValueError(f"{where}.clusters must be an integer of at least "
                             f"two, got {clusters!r}")
        expected = abs(measured - reference)
        if abs(absolute - expected) > 8.0 * np.spacing(max(expected, 1.0)):
            raise ValueError(
                f"{where}.absolute_error is {absolute!r}, not "
                f"|measured - reference| = {expected!r}; a report whose "
                f"derived fields disagree with its own measurements is not "
                f"describing one comparison")
        scale = abs(reference)
        wanted = expected if scale == 0.0 else expected / scale
        if abs(relative - wanted) > 8.0 * np.spacing(max(wanted, 1.0)):
            raise ValueError(f"{where}.relative_error is {relative!r}, not "
                             f"{wanted!r}")
        key = (level["observable"], level["position"], timestep)
        if key in seen:
            raise ValueError(f"{where} repeats the observable, position and "
                             f"timestep of an earlier level")
        seen.add(key)
        identity = (level["observable"], level["position"])
        ladders.setdefault(identity, set()).add(timestep)
        ladder_units.setdefault(identity, set()).add(level["unit"])
        rows_by_identity.setdefault(identity, []).append(level)
    schedule_seen = ()
    if not ladders:
        if passed:
            raise ValueError("a passing report carries no ladders")
    else:
        schedules = {tuple(sorted(steps, reverse=True))
                     for steps in ladders.values()}
        schedule_seen = sorted(next(iter(schedules)), reverse=True)
        if len(schedules) != 1:
            raise ValueError(
                f"the report carries ladders on {len(schedules)} different "
                f"schedules {sorted(schedules)}; one comparison runs one "
                f"frozen schedule, and a truncated ladder for one identity is "
                f"how a failing observable disappears from a report")

    # ---- the frozen contract, rebuilt --------------------------------------
    contract = payload["contract"]
    if not isinstance(contract, dict) \
            or set(contract) != set(CONTRACT_FIELDS):
        raise ValueError(f"contract must carry exactly {list(CONTRACT_FIELDS)}")
    sampling = contract["sampling"]
    if not isinstance(sampling, dict) \
            or set(sampling) != set(SAMPLING_FIELDS):
        raise ValueError(f"contract.sampling must carry exactly "
                         f"{list(SAMPLING_FIELDS)}")
    try:
        design = SamplingDesign(sampling["unit"], sampling["clusters"],
                                sampling["replications"], sampling["method"],
                                sampling["resamples"], sampling["coverage"],
                                sampling["seed"])
        frozen = RefinementContract(tuple(contract["timesteps"]),
                                    contract["ratio"], design)
    except (TypeError, ValueError) as error:
        raise ValueError(f"the report's frozen contract is not a contract: "
                         f"{error}") from None
    if schedule_seen \
            and tuple(sorted(frozen.timesteps, reverse=True)) \
            != tuple(schedule_seen):
        raise ValueError(
            f"the levels run at timesteps {list(schedule_seen)} while the "
            f"report's own frozen contract declares "
            f"{sorted(frozen.timesteps, reverse=True)}")

    # ---- the frozen budgets, rebuilt ---------------------------------------
    if not isinstance(payload["budgets"], list) or not payload["budgets"]:
        raise ValueError("budgets must be a non-empty list")
    rebuilt = []
    for position, entry in enumerate(payload["budgets"]):
        if not isinstance(entry, dict) or set(entry) != set(BUDGET_FIELDS):
            raise ValueError(f"budgets[{position}] must carry exactly "
                             f"{list(BUDGET_FIELDS)}")
        _require_label(entry["observable"], f"budgets[{position}].observable")
        _require_printable(entry["position"], f"budgets[{position}].position",
                           allow_empty=True)
        try:
            rebuilt.append(ValidationBudget(
                entry["observable"], entry["unit"], entry["absolute"],
                entry["relative"], entry["noise_floor"],
                entry["require_decrease"], entry["position"],
                entry["estimator"], entry["quantile"]))
        except (TypeError, ValueError) as error:
            raise ValueError(f"budgets[{position}] is not a budget: "
                             f"{error}") from None
    try:
        budgets = FrozenBudgets(payload["label"], frozen, tuple(rebuilt))
    except (TypeError, ValueError) as error:
        raise ValueError(f"the report's frozen budgets are not a budget set: "
                         f"{error}") from None
    if budgets.digest != payload["budget_digest"]:
        raise ValueError(
            f"the report's own budgets hash to {budgets.digest!r} while it "
            f"declares {payload['budget_digest']!r}; the contract a verdict "
            f"cites has to be the contract it carries")

    # ---- one unit per ladder, and every ladder budgeted --------------------
    for identity, units in ladder_units.items():
        if len(units) != 1:
            raise ValueError(
                f"ladder {identity!r} carries units {sorted(units)}; a ladder "
                f"is one observable in one unit, and a mixed-unit ladder "
                f"cannot be compared with a cap in either")
    for identity in ladders:
        try:
            budget = budgets.budget_for(*identity)
        except KeyError:
            continue                 # an unbudgeted identity; the codes carry it
        if budget.unit != ladder_units[identity].copy().pop():
            raise ValueError(f"ladder {identity!r} is measured in a unit its "
                             f"own budget does not declare")

    # ---- recompute the verdict --------------------------------------------
    coverage = frozen.sampling.coverage
    # Only the *ladder* codes are locally recomputable.  The machinery codes --
    # a sample whose cluster count or identities do not match the design, an
    # identity that was sampled without a budget or budgeted without a sample --
    # are decided before any level exists, so a report carrying one carries no
    # rows to recompute it from.  They are accepted here and recomputed by
    # :func:`parse_validation_report`, which has the dataset.  An earlier
    # revision tried to recompute them here and crashed on the first legitimate
    # ``cluster_mismatch`` report, which carries no levels at all.
    machinery = {"unsampled", "unbudgeted", "cluster_mismatch",
                 "identity_mismatch", "schedule_mismatch", "unit_mismatch"}
    recomputed = []
    for identity in sorted(ladders):
        try:
            budget = budgets.budget_for(*identity)
        except KeyError:
            continue
        ladder = sorted(rows_by_identity[identity],
                        key=lambda row: row["timestep"], reverse=True)
        if len(ladder) != len(frozen.timesteps):
            raise ValueError(f"ladder {identity!r} has {len(ladder)} levels "
                             f"against the contract's {len(frozen.timesteps)}")
        if ladder[0]["paired_error"] != 0.0:
            raise ValueError(
                f"ladder {identity!r} reports a paired error of "
                f"{ladder[0]['paired_error']!r} at its coarsest level; the "
                f"paired error is the change from the *next coarser* level "
                f"and the coarsest has none, so a nonzero value there is "
                f"evidence nobody produced")
        spans = {row["span_error"] for row in ladder}
        if len(spans) != 1:
            raise ValueError(f"ladder {identity!r} reports {len(spans)} "
                             f"different span errors; the span is one "
                             f"end-to-end quantity per ladder")
        objects = [RefinementLevel(row["observable"], row["unit"],
                                   row["timestep"], row["measured"],
                                   row["reference"], row["standard_error"],
                                   row["clusters"], row["position"],
                                   row["paired_error"], row["span_error"])
                   for row in ladder]
        for code, _ in _ladder_codes(budget, objects, coverage):
            recomputed.append((identity[0], identity[1], code))

    declared_blocking = []
    if not isinstance(payload["blocking"], list):
        raise ValueError("blocking must be a list")
    for position, entry in enumerate(payload["blocking"]):
        if not isinstance(entry, dict) or set(entry) != set(BLOCKING_FIELDS):
            raise ValueError(f"blocking[{position}] must carry exactly "
                             f"{list(BLOCKING_FIELDS)}")
        # The shared Unicode policy, *before* the code is classified.  An
        # independent review found a machinery-only blocking row — one whose
        # code is excluded from local recomputation — carrying an observable
        # and a position with forbidden categories in them, accepted here
        # while every public constructor refuses them.  A shape checker that
        # advertises one policy has to apply it to every label it sees, not
        # only to those a later reconstruction happens to revisit.
        _require_label(entry["observable"], f"blocking[{position}].observable")
        _require_printable(entry["position"], f"blocking[{position}].position",
                           allow_empty=True)
        if entry["code"] not in REASON_CODES:
            raise ValueError(f"blocking[{position}] code {entry['code']!r} is "
                             f"not one of {list(REASON_CODES)}")
        declared_blocking.append((entry["observable"], entry["position"],
                                  entry["code"]))
    local = sorted(entry for entry in declared_blocking
                   if entry[2] not in machinery)
    if local != sorted(recomputed):
        raise ValueError(
            f"the report declares {len(local)} locally recomputable blocking "
            f"code(s) and its own levels, budgets and contract imply "
            f"{len(recomputed)}: {sorted(set(recomputed) ^ set(local))}."
            f"  The verdict is recomputed here rather than believed, which is "
            f"what stops a numerical_no_result from being rewritten to a pass "
            f"with every level and both digests left byte-identical")
    if len(declared_blocking) != len(reasons):
        raise ValueError(f"{len(declared_blocking)} blocking codes beside "
                         f"{len(reasons)} reasons")
    if passed != (not declared_blocking):
        raise ValueError(
            f"the report claims passed={passed!r} beside "
            f"{len(declared_blocking)} blocking code(s).  This is the "
            f"coherent-rewrite bypass: verdict, passed and reasons can be "
            f"changed together and stay internally consistent, so the parser "
            f"recomputes the decision from the levels instead")

    # ---- envelopes rebuilt field for field ---------------------------------
    envelopes = payload["envelopes"]
    if not isinstance(envelopes, dict):
        raise ValueError("envelopes must be a mapping from unit to envelope")
    finest = {}
    for identity in ladders:
        rows = rows_by_identity[identity]
        best = min(rows, key=lambda row: row["timestep"])
        finest[identity] = best
    expected = {}
    for unit in UNITS:
        rows = [row for row in finest.values() if row["unit"] == unit]
        if not rows:
            continue
        worst_absolute = max(rows, key=lambda row: row["absolute_error"])
        worst_relative = max(rows, key=lambda row: row["relative_error"])
        expected[unit] = {
            "unit": unit,
            "absolute": worst_absolute["absolute_error"],
            "absolute_observable": worst_absolute["observable"],
            "absolute_position": worst_absolute["position"],
            "standard_error": worst_absolute["standard_error"],
            "relative": worst_relative["relative_error"],
            "relative_observable": worst_relative["observable"],
            "relative_position": worst_relative["position"],
            "timestep": worst_absolute["timestep"]}
    if set(envelopes) != set(expected):
        raise ValueError(
            f"the report carries envelopes for {sorted(envelopes)} while its "
            f"levels imply exactly {sorted(expected)}; a missing envelope "
            f"hides a unit and an extra one describes rows that are not there")
    for unit, envelope in envelopes.items():
        if not isinstance(envelope, dict) \
                or set(envelope) != set(ENVELOPE_FIELDS):
            raise ValueError(f"envelopes[{unit!r}] must carry exactly "
                             f"{list(ENVELOPE_FIELDS)}")
        for field, wanted in expected[unit].items():
            if envelope[field] != wanted:
                raise ValueError(
                    f"envelopes[{unit!r}].{field} is {envelope[field]!r}; the "
                    f"report's own levels give {wanted!r}.  An envelope is "
                    f"derived from the levels, so a value the levels do not "
                    f"attain is a forgery rather than a summary")
    return payload


def _ladder_codes(budget, ladder, coverage):
    """The blocking reason codes one ladder earns, as ``(code, detail)`` pairs.

    **The single source of the decision.**  :func:`compare_refinement` calls it
    to reach a verdict and :func:`parse_validation_report` calls it to
    recompute one from a serialized report, so the two cannot disagree: a
    report whose tag says ``pass`` while its own levels and budgets imply a
    blocking code is refused at the trust boundary.

    ``ladder`` is that identity's levels, coarsest first, already checked for
    unit, schedule and cluster agreement.
    """
    codes = []
    finest, coarsest = ladder[-1], ladder[0]
    if finest.bounded_error(coverage) > budget.absolute:
        codes.append(("absolute_cap",
                      f"{finest.absolute_error:.6g} {budget.unit} with "
                      f"standard error {finest.standard_error:.6g} gives a "
                      f"{coverage:.6g}-sigma bound "
                      f"{finest.bounded_error(coverage):.6g} past the cap "
                      f"{budget.absolute:.6g}"))
    if finest.relative_error > budget.relative:
        codes.append(("relative_cap",
                      f"relative error {finest.relative_error:.6g} past the "
                      f"cap {budget.relative:.6g}"))
    if not budget.require_decrease:
        return codes
    span = budget.noise_floor + coverage * finest.span_error
    if finest.absolute_error > coarsest.absolute_error + span:
        codes.append(("not_converging",
                      f"{coarsest.absolute_error:.6g} at "
                      f"{coarsest.timestep:.6g} against "
                      f"{finest.absolute_error:.6g} at "
                      f"{finest.timestep:.6g}, past the single whole-ladder "
                      f"allowance {span:.6g}"))
    reversals = 0
    for coarse, fine in zip(ladder, ladder[1:]):
        growth = fine.absolute_error - coarse.absolute_error
        if growth <= 0.0:
            continue
        reversals += 1
        allowance = budget.noise_floor + coverage * fine.paired_error
        if growth > allowance:
            codes.append(("reversal_beyond_allowance",
                          f"a rise of {growth:.6g} at timestep "
                          f"{fine.timestep:.6g} past the resolution "
                          f"{budget.noise_floor:.6g} plus {coverage:.6g} "
                          f"paired standard errors "
                          f"{fine.paired_error:.6g}"))
        elif reversals > 1:
            codes.append(("repeated_reversal",
                          f"{reversals} adjacent reversals; one upward step "
                          f"can be statistically unresolved, a ladder that "
                          f"turns around more than once is not converging"))
    return codes




def parse_validation_report(payload, dataset) -> dict:
    """Read a validation report **and recompute it from the frozen data**.

    Returned unchanged on success.  This is the trust boundary, and it is
    rooted: ``dataset`` is a :class:`ValidationDataset` the caller obtained
    independently, its digest is required to be the one the report declares,
    and the whole comparison is then **re-run** from it and compared with the
    payload field for field.

    ``dataset`` is required and has no default.  A report cannot be validated
    without the evidence it summarizes, and an independent review showed why: a
    genuine ``numerical_no_result`` whose two-cluster bootstrap standard error
    was 0.4999 became a parser-approved pass when only the standard, paired and
    span errors were zeroed, with both digests unchanged.  The parser
    faithfully recomputed a pass from invented summaries because nothing tied
    those summaries to the data.  Every one of them is recomputed here.

    Recomputing the entire report, rather than checking a list of invariants
    over it, is deliberate: it is the only formulation that cannot be
    incomplete.  Any field a forger could reach — a level's measured value, a
    cluster count, an envelope, a reason's wording, a blocking code, the
    verdict — is produced again from the dataset and the embedded frozen
    contract, and any difference at all is a refusal.  It also makes every
    verdict :func:`compare_refinement` can emit round-trip by construction,
    including the machinery no-results (``cluster_mismatch``, ``unbudgeted``)
    that carry no level rows and that an earlier parser crashed on.

    **Why the evidence is not embedded in the report.**  The alternative is to
    carry the cluster observations inside the payload.  For the audit ladders
    that would be a few thousand numbers, but the stationary S3 dataset holds
    twelve samples of six thousand clusters at three levels — 216,000
    observations, several megabytes of JSON for a report whose purpose is to be
    read — so the evidence travels as an object and the report carries its
    digest.  The guarantee is the same: the digest is checked against the
    dataset, and everything else is recomputed from it.
    """
    require_exact_type(dataset, ValidationDataset, "dataset")
    # An owned plain-``dict`` snapshot of the payload, taken before anything is
    # read twice.  A mapping subclass whose ``__getitem__`` mutates on access
    # is how an independent review made its time-of-check/time-of-use race
    # deterministic; a snapshot removes the second read entirely, and the
    # snapshot is what is validated and returned.
    # The trusted snapshot **first**, before any untrusted payload method runs.
    # With bytes-backed evidence either order is safe, but an independent
    # review's attack ran during ``json.dumps(payload)`` -- a mapping subclass
    # whose ``items()`` mutated the dataset the parser had not yet copied -- and
    # taking the snapshot first removes the window rather than relying on the
    # backing to make it harmless.
    dataset = ValidationDataset(dataset.label, dataset.contract,
                                dataset.samples)
    trusted = dataset.digest
    payload = json.loads(json.dumps(payload))
    inspect_validation_report(payload)

    contract = payload["contract"]
    sampling = contract["sampling"]
    design = SamplingDesign(sampling["unit"], sampling["clusters"],
                            sampling["replications"], sampling["method"],
                            sampling["resamples"], sampling["coverage"],
                            sampling["seed"])
    frozen = RefinementContract(tuple(contract["timesteps"]),
                                contract["ratio"], design)
    budgets = FrozenBudgets(
        payload["label"], frozen,
        tuple(ValidationBudget(entry["observable"], entry["unit"],
                               entry["absolute"], entry["relative"],
                               entry["noise_floor"], entry["require_decrease"],
                               entry["position"], entry["estimator"],
                               entry["quantile"])
              for entry in payload["budgets"]))
    if trusted != payload["dataset_digest"]:
        raise ValueError(
            f"the supplied dataset hashes to {trusted!r} and the report "
            f"declares {payload['dataset_digest']!r}.  A report is a summary "
            f"of one set of observations; validating it against a different "
            f"set would establish nothing about either")
    if dataset.contract.canonical != frozen.canonical:
        raise ValueError("the dataset was collected under a different "
                         "refinement contract than the report embeds")

    recomputed = compare_refinement(dataset, budgets,
                                    budgets.digest).as_dict()
    # The digest again, after the recomputation.  With an owned read-only
    # snapshot this cannot fail; it is asserted anyway, because the guarantee
    # being made is "the bytes that produced this verdict are the bytes whose
    # digest was checked", and an assertion is cheaper than trusting that no
    # future edit reintroduces an alias.
    if dataset.digest != trusted:
        raise ValueError(
            f"the dataset's digest changed from {trusted!r} to "
            f"{dataset.digest!r} during recomputation; the evidence a verdict "
            f"rests on moved underneath it")
    if recomputed != payload:
        differing = sorted(key for key in set(recomputed) | set(payload)
                           if recomputed.get(key) != payload.get(key))
        raise ValueError(
            f"the report does not match what its own frozen data and contract "
            f"produce; {differing} differ.  Every number here is recomputed "
            f"from the dataset rather than believed, which is what stops a "
            f"blocked result from being rewritten into a pass by editing the "
            f"summaries nothing was checking")
    return payload


def compare_refinement(dataset, budgets, declared_digest) -> ValidationVerdict:
    """Hold a paired refinement ladder to its frozen contract.

    ``dataset`` is a :class:`ValidationDataset` of frozen cluster-level
    observations.  ``declared_digest`` is the digest the caller recorded
    *before* opening any result; it must equal ``budgets.digest``, and a
    mismatch raises rather than returning a no-result, because a contract that
    changed between freezing and comparison is a broken procedure.

    **Everything numerical is computed here.** The statistic is taken with the
    estimator the frozen budget names, and the uncertainty is a cluster
    bootstrap over the frozen sampling design's own clusters, seeded from its
    frozen key, with every refinement level evaluated on the *same* resample so
    that level-to-level contrasts are paired rather than assumed independent.
    There is no argument through which an uncertainty can be supplied.

    Every ``(observable, position)`` identity is judged **on its own**, and all
    of them must pass:

    1. it is both budgeted and sampled, and its sample has the contract's
       cluster count and level count;
    2. at the finest timestep the **uncertainty-aware** absolute error —
       ``|error| + coverage * standard_error`` — is within the absolute cap,
       and the relative point estimate is within the relative cap;
    3. the error at the finest timestep is no larger than at the coarsest,
       within **one** allowance for the whole ladder: the declared resolution
       plus the coverage factor times the *paired* bootstrap standard error of
       that end-to-end change; and
    4. at most one adjacent pair may go the wrong way, and that one must be
       inside the declared resolution plus the coverage factor times the paired
       bootstrap standard error of that step.

    Clauses 3 and 4 are lifted, and only those two, for an observable whose
    budget declares ``require_decrease=False``.

    Failure of any clause yields ``"numerical_no_result"`` with the reasons
    attached.  That is a legitimate outcome and not a defect: a reference
    sample that cannot resolve a cap should say so.
    """
    require_exact_type(budgets, FrozenBudgets, "budgets")
    if isinstance(dataset, (list, tuple)):
        raise TypeError(
            "compare_refinement takes a ValidationDataset of frozen "
            "cluster-level observations, not a sequence of RefinementLevel "
            "objects.  That signature is gone deliberately: it let a caller "
            "supply the standard error that decided the verdict, and the same "
            "ladder returned a no-result or a pass depending on a scalar "
            "nobody had measured")
    require_exact_type(dataset, ValidationDataset, "dataset")
    # Rebuilt from their own fields, so an object edited past its constructor
    # is refused here rather than used as found.
    frozen = FrozenBudgets(budgets.label, budgets.contract, budgets.budgets)
    data = ValidationDataset(dataset.label, dataset.contract, dataset.samples)
    declared = _require_label(declared_digest, "declared_digest")

    if declared != frozen.digest:
        raise ValueError(
            f"the declared budget digest {declared!r} is not the digest of the "
            f"budget supplied, {frozen.digest!r}.  The contract a validation "
            f"is held to -- its caps, its refinement schedule, its level "
            f"count, its estimators and its sampling design -- has to be "
            f"frozen before its results are opened")
    if data.contract.canonical != frozen.contract.canonical:
        raise ValueError(
            f"the dataset was collected under a different refinement contract "
            f"than the budget declares.  The two have to be the same contract, "
            f"or the schedule and the sampling design a verdict cites are not "
            f"the ones its numbers came from")

    contract = frozen.contract
    schedule = tuple(sorted(contract.timesteps, reverse=True))
    coverage = contract.sampling.coverage
    reasons = []
    blocking = []
    records = []

    def block(identity, code, detail):
        blocking.append((identity[0], identity[1], code))
        reasons.append(f"{identity!r}: [{code}] {detail}")

    sampled = sorted((entry.observable, entry.position)
                     for entry in data.samples)
    budgeted = sorted(entry.identity for entry in frozen.budgets)
    for name in sampled:
        if name not in budgeted:
            block(name, "unbudgeted",
                  "sampled but not budgeted, so no frozen tolerance governs it")
    for name in budgeted:
        if name not in sampled:
            block(name, "unsampled", "budgeted but never sampled")

    for name in budgeted:
        if name not in sampled:
            continue
        budget = frozen.budget_for(*name)
        sample = data.sample_for(*name)

        if sample.identities != data.samples[0].identities:
            block(name, "identity_mismatch",
                  "measured over different physical clusters from the rest of "
                  "the dataset")
            continue
        if sample.clusters != contract.sampling.clusters:
            block(name, "cluster_mismatch",
                  f"{sample.clusters} clusters against the frozen "
                  f"{contract.sampling.clusters}")
            continue
        if sample.levels != len(schedule):
            block(name, "schedule_mismatch",
                  f"{sample.levels} levels against the frozen schedule's "
                  f"{len(schedule)}")
            continue

        errors = sample.errors(budget.estimator, budget.quantile)
        measured = sample.measured(budget.estimator, budget.quantile)

        # One set of cluster resamples, every level evaluated on each of them,
        # and the **folded** absolute error -- the statistic the gate consumes.
        rng = np.random.default_rng(contract.sampling.seed)
        draws = np.empty((contract.sampling.resamples, sample.levels),
                         dtype=float)
        for index in range(contract.sampling.resamples):
            rows = rng.integers(0, sample.clusters, size=sample.clusters)
            draws[index] = sample.errors(budget.estimator, budget.quantile,
                                         rows)
        marginal = np.std(draws, axis=0, ddof=1)
        paired = np.zeros(sample.levels, dtype=float)
        for level in range(1, sample.levels):
            paired[level] = float(np.std(draws[:, level] - draws[:, level - 1],
                                         ddof=1))
        span_error = float(np.std(draws[:, -1] - draws[:, 0], ddof=1))

        ladder = []
        for level, timestep in enumerate(schedule):
            ladder.append(RefinementLevel(
                budget.observable, budget.unit, timestep,
                float(measured[level]), sample.reference,
                float(marginal[level]), sample.clusters, budget.position,
                float(paired[level]), span_error))
        records.extend(ladder)
        for code, detail in _ladder_codes(budget, ladder, coverage):
            block(name, code, detail)

    verdict = "pass" if not reasons else "numerical_no_result"
    return ValidationVerdict(frozen.label, frozen.digest, verdict,
                             tuple(records), tuple(reasons), data.digest,
                             tuple(blocking), frozen)
