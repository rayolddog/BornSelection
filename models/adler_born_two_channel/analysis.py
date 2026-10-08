#!/usr/bin/env python3
"""
analysis.py — the predeclared binomial complementary-log-log exponent, and the
explicit refusal to produce one.

The plan's primary scaling estimand, in words:

    the logarithm of the endpoint cumulative hazard is a constant plus an
    unconstrained exponent times the logarithm of the peak coupling.

In symbols, with ``q_K`` the commitment probability at peak coupling ``K``:

    log(-log(1 - q_K)) = amplitude + exponent * log(K)

which is a binomial generalized linear model with the complementary-log-log
link, fitted **directly to the committed and unresolved counts** of every
frozen coupling cell.

Why the counts and not the proportions
--------------------------------------
The obvious thing to do is compute ``q_K``, take ``-log(1 - q_K)``, take its
logarithm and fit a straight line.  That procedure cannot represent a cell in
which nothing committed (``q = 0``, first logarithm ``-inf``) or one in which
everything did (``q = 1``, second logarithm ``+inf``), so in practice those
cells get dropped — and the cells at the ends of a coupling sweep are exactly
the ones that carry the most information about its slope.  Dropping them is
silent post-selection: the estimate then answers a question about an interior
subset while presenting itself as answering the whole sweep.

This module therefore never transforms a proportion.  ``y_i`` committed out of
``n_i`` trials enters the likelihood as a binomial count, a zero-event cell
contributes ``-n_i H_i`` and an all-event cell contributes ``n_i log q_i``, and
both are as ordinary as any other cell.  ``analytic.fit_log_log_exponent``, the
unweighted least-squares fitter this package already has, is deliberately *not*
used for binomial outcomes; it exists for the analytic prediction, whose points
are exact and positive.

There is no cell filter anywhere in this file.  Every report states how many
cells the design declared and how many entered the likelihood, and the two are
always equal; a design whose data cannot support a fit produces **no valid
exponent**, never a fit over a friendlier subset.

What "no valid exponent" means
------------------------------
It is a result, not an error.  :class:`ExponentReport` always comes back, always
carries the per-cell counts, and states one of two verdicts with the reasons
attached.  The reasons are a closed vocabulary (:data:`NO_EXPONENT_REASONS`)
so that a downstream consumer can branch on them rather than parse prose:

``saturated``
    every cell has zero events, or every cell has all events.  There is no
    interior information about a slope and the maximum-likelihood estimate is
    at infinity.
``separation``
    the linear predictor runs away: the data are perfectly ordered into
    all-zero and all-one cells with nothing in between, so the likelihood is
    monotone and the estimate is again at infinity.
``flat``
    the expected information matrix is numerically singular.  The likelihood
    has no curvature to locate a maximum, so nothing is estimable.

    This reason is about a flat *likelihood*, not a flat *response*.  A sweep
    whose commitment probability genuinely does not vary with coupling, and
    which is adequately informed, returns a **valid exponent near zero** — and
    it should: the declared model contains ``p = 0``, and refusing it would
    erase a real negative scaling result rather than report one.  An
    independent review reached the same disposition, and the case where a flat
    response carries no information is caught by ``uninformative`` instead.
``ill_conditioned``
    the design's own log-coupling matrix is rank-deficient or nearly so — for
    example a "sweep" whose couplings are all equal.
``non_convergence``
    Fisher scoring did not reach the declared gradient tolerance.
``lack_of_fit``
    the residual deviance is too large for its degrees of freedom at the frozen
    level.  The power law does not describe these counts.
``curvature``
    the frozen curvature extension improves the fit significantly.  The
    response is curved in log-coupling, so a single exponent is a local summary
    and is refused as a global one.
``uninformative``
    the fit converged and passed the goodness tests, but the 95 percent
    interval on the exponent is wider than the frozen cap.  An interval that
    contains both one and two answers no question this experiment is asking.

Every threshold that appears in that list is a module constant, declared here
and never taken from a caller, so that a failing fit cannot be turned into a
passing one by widening a budget at the call site.

Common random numbers
---------------------
A coupling sweep run on one keyed noise hierarchy shares master trial
identifiers across cells: trial 7 at ``K = 0.4`` and trial 7 at ``K = 0.8`` are
the *same* Brownian path driven at two couplings.  Those cells are positively
correlated, and a bootstrap that resamples within each cell independently
treats them as though they were not — which understates the uncertainty on the
exponent, sometimes badly.

So the resampling unit is a declared property of the design
(:class:`ScalingDesign`), it is checked rather than assumed — a
``master_trial`` design must have *identical* trial identifiers in every cell,
an ``independent_namespace`` design must have *disjoint* ones — and the wrong
procedure is available under a name that says what it is,
``cell_independent_control``, so that the difference between the two can be
measured rather than argued about.

Non-claims
----------
A statistically adequate exponent is not a Born-rule result.  It is not a
detector click, an outcome, a unique actuality, a two-channel result, or
evidence about the origin of the noise.  It is the coupling dependence of a
stopping rule this study imposed on a model bath.  An exponent near two is
not, by itself, support for the width-times-rate mechanism either; that
attribution requires the negative controls in ``compare.py`` to fail to
reproduce it.  See :data:`ANALYSIS_NON_CLAIMS`.
"""
from __future__ import annotations

import hashlib
import math
from dataclasses import dataclass

import numpy as np

from .validation import (require_finite, require_instance,
                         require_positive_int)

__all__ = [
    "ANALYSIS_NON_CLAIMS",
    "CONFIDENCE_LEVEL",
    "CONTROL_RESAMPLING_UNIT",
    "CURVATURE_ALPHA",
    "DESIGN_RESAMPLING_UNITS",
    "LACK_OF_FIT_ALPHA",
    "LINK_FORMULA",
    "MAX_EXPONENT_HALF_WIDTH",
    "MIN_CELLS",
    "ESTIMABILITY_REASONS",
    "GOODNESS_REASONS",
    "LACK_OF_FIT_REPLICATES",
    "MIN_VALID_RESAMPLE_FRACTION",
    "MIN_RESAMPLE_REPLICATES",
    "NO_EXPONENT_REASONS",
    "RESAMPLE_REPLICATES",
    "RESAMPLING_UNITS",
    "VERDICTS",
    "CalibratedExponent",
    "CouplingCell",
    "ExponentReport",
    "ClusterScoreTest",
    "MasterTrial",
    "ResamplingReport",
    "ScalingDesign",
    "calibrate_exponent",
    "chi_square_upper_tail",
    "cluster_curvature",
    "cluster_lack_of_fit",
    "curvature_column",
    "fit_exponent",
    "master_trials",
    "resample_exponent",
]

#: The token every *derived* record in this layer demands from its factory.
#:
#: An independent review assembled a physical permission out of public
#: constructors: a ``CalibratedExponent`` whose published verdict and interval
#: agreed with neither record it carried, and a ``ScalingVerdict`` whose
#: fingerprints were the string ``"not-a-fingerprint"``.  Those are ordinary
#: supported-API correctness failures, not an adversary model, and the fix is
#: the one this package already uses for the prediction boundary: make the
#: unsupported route *impossible to take by accident and visible when taken*.
#:
#: A derived record therefore cannot be built without this object, and this
#: object cannot be named without writing a leading underscore.  Every factory
#: below passes it; ``verify.py`` imports it deliberately, in one place, to
#: build the evidence whose *shape* it needs to test.  It is not a
#: cryptographic seal and is not described as one — a determined caller can
#: import it, exactly as a determined caller can ``importlib`` the prediction.
#: What it removes is the accident.
_SEAL = object()


def _require_seal(seal, what: str) -> None:
    """Raise unless ``seal`` is the module's own factory token."""
    if seal is not _SEAL:
        raise TypeError(
            f"{what} is a derived record and cannot be constructed directly; "
            f"it is produced by this module's factories, which compute every "
            f"field they publish.  A hand-built one could disagree with the "
            f"records it carries, which is exactly how a forged interval "
            f"became a physical permission")

#: The settled model, written out so a reader never has to infer it from code.
LINK_FORMULA = ("log(-log(1 - q_K)) = amplitude + exponent * log(K), fitted as "
                "a binomial complementary-log-log likelihood directly to the "
                "committed and unresolved counts of every frozen coupling cell")

#: The frozen curvature extension.  Centred on the design's own mean log
#: coupling, which is a property of the frozen grid and not of any outcome, so
#: the extension cannot be shifted after seeing the counts.
CURVATURE_FORMULA = ("the same model plus curvature * (log(K) - mean log(K))^2, "
                     "compared by likelihood ratio on one degree of freedom")

#: Two-sided level of every reported interval.
CONFIDENCE_LEVEL = 0.95

#: Normal quantile for :data:`CONFIDENCE_LEVEL`, written to sixteen digits
#: rather than computed, so the interval a report quotes cannot move.
CONFIDENCE_Z = 1.959963984540054

#: Significance level at which the residual deviance refuses the power law.
LACK_OF_FIT_ALPHA = 0.05

#: Significance level at which the curvature extension refuses the power law.
CURVATURE_ALPHA = 0.05

#: Smallest design this module will fit.  Four cells give the curvature model
#: one residual degree of freedom, so both frozen tests can actually run; a
#: three-cell design would make the curvature model saturated and its
#: lack-of-fit test vacuous.
MIN_CELLS = 4

#: Largest 95-percent half-width on the exponent that still answers a question.
#:
#: An interval of +-1 around a fitted exponent of 1.5 contains both the
#: tongue-width-only expectation and the width-times-rate one, so it
#: discriminates nothing.  A fit that wide is reported as ``uninformative``
#: rather than as a number somebody will quote.
MAX_EXPONENT_HALF_WIDTH = 1.0

#: Linear-predictor magnitude beyond which the estimate is treated as having
#: run away.  ``exp(30)`` is a cumulative hazard of order ``1e13``; no finite
#: data set distinguishes that from infinity.
MAX_ABS_LINEAR_PREDICTOR = 30.0

#: Reciprocal condition number below which an information or design matrix is
#: called singular.
MIN_RECIPROCAL_CONDITION = 1e-10

#: Fisher-scoring stopping rule, applied to the *Newton decrement*
#: ``gradient . information^-1 . gradient`` rather than to the raw gradient.
#:
#: The raw gradient has the units of a count, so an absolute threshold on it
#: means something different for four thousand trials per cell than for forty,
#: and a fit at machine precision in the large-sample case would be reported as
#: unconverged.  The decrement is the quadratic model's remaining improvement
#: in log-likelihood, which is dimensionless and comparable across designs.
GRADIENT_TOLERANCE = 1e-10
MAX_ITERATIONS = 200

#: The two verdicts, and nothing else.
VERDICTS = ("valid_exponent", "no_valid_exponent")

#: The closed vocabulary of refusals.  See the module docstring for what each
#: one means; they are listed here in the order a report emits them.
NO_EXPONENT_REASONS = ("saturated", "separation", "flat", "ill_conditioned",
                       "non_convergence", "lack_of_fit", "curvature",
                       "curvature_unidentified", "uninformative",
                       "inadequate_resampling", "uncalibrated_fit",
                       "uncalibrated_curvature")

#: Multiplier-bootstrap replicates for the cluster-calibrated lack-of-fit test.
#:
#: Frozen here, not taken from a caller, and deliberately not accompanied by a
#: seed argument either: :func:`cluster_lack_of_fit` derives its seed from the
#: design's own digest, so the test has *nothing* a caller could tune after
#: seeing a result.  499 gives a p-value on a 1/500 grid, which resolves the
#: frozen five-percent level with room to spare.
LACK_OF_FIT_REPLICATES = 499

#: Smallest number of independent clusters the cluster-robust covariance is
#: believed at.  Below this the sandwich's finite-cluster bias dominates and
#: the test is reported as *uncalibrated* rather than as a p-value.
MIN_LACK_OF_FIT_CLUSTERS = 30

#: The two decisions the cluster score test makes.
CLUSTER_TEST_KINDS = ("lack_of_fit", "curvature")

#: Bootstrap replicates every **authoritative** exponent record uses.
#:
#: Frozen, and not a caller argument on any decision path.  An independent
#: review showed the previous floor of 100 was adequate as a *floor* and not as
#: an authority: on interior exact-null data the curvature percentile interval
#: excluded zero 7.33 percent of the time at 100 replicates and 4.33 percent at
#: 499.  A record may not be called adequate before its Monte Carlo resolution
#: is adequate for what it is being asked to decide.
RESAMPLE_REPLICATES = 499

#: Smallest resampling record a calibrated exponent will accept.
#:
#: Both numbers are frozen here rather than taken from a caller, for the reason
#: every other threshold in this file is: a resampling record that was allowed
#: to be as small as its caller liked could always be made adequate by asking
#: for less.
MIN_RESAMPLE_REPLICATES = 100

#: Fraction of drawn replicates whose exponent must be **estimable** before the
#: resampling record is believed.
#:
#: An independent review found the previous implementation counting any
#: replicate with a finite exponent, which turned a design whose every replicate
#: was refused into a report saying ``50 / 50`` valid with a percentile interval
#: beside it.  Refusals are now carried by reason and a record below this
#: fraction is *inadequate* — which blocks, rather than producing a narrower
#: interval from the few replicates that happened to survive.
MIN_VALID_RESAMPLE_FRACTION = 0.9

#: The refusal reasons that mean a replicate has **no usable exponent at all**.
#:
#: The distinction matters and is the one place this module departs from
#: counting ``ExponentReport.is_valid`` directly, so it is written down rather
#: than buried.
#:
#: These six say the estimate does not exist: the likelihood is monotone, the
#: information is singular, the design is rank-deficient, the optimizer did not
#: converge, or the required curvature diagnostic is unidentifiable.  A
#: replicate refused for one of them contributes nothing and is counted as a
#: failure.
#:
#: :data:`GOODNESS_REASONS` are different.  ``lack_of_fit`` and ``curvature``
#: are *model-adequacy* verdicts from the marginal deviance and the marginal
#: likelihood ratio — and those are exactly the two statistics whose null
#: distribution common random numbers invalidate.  A paired bootstrap resamples
#: whole master trials, so its replicates are clustered and both tests reject
#: far more often than their nominal level: at two thousand trials a cell, a
#: perfectly specified power law has sixty percent of its replicates rejected
#: for lack of fit.  Counting those as failures would make every genuine
#: common-random-number design "inadequate", which is a blanket refusal rather
#: than a diagnostic.
#:
#: So model adequacy under pairing is decided by the **bootstrap distribution
#: of the curvature coefficient** — :attr:`ResamplingReport.curvature_interval`,
#: which is cluster-robust by construction — and by the marginal tests applied
#: once to the original fit, where they act as a conservative refusal.  The
#: per-replicate goodness rejections are still counted and reported in
#: :attr:`ResamplingReport.refusals`, because their *rate* is informative even
#: though it is not a gate.
ESTIMABILITY_REASONS = ("saturated", "separation", "flat", "ill_conditioned",
                        "non_convergence", "curvature_unidentified")

#: The two model-adequacy refusals.  See :data:`ESTIMABILITY_REASONS`.
GOODNESS_REASONS = ("lack_of_fit", "curvature")

#: Resampling units a *design* may declare.  Both are legitimate answers to
#: common random numbers: share every master trial across cells and resample
#: the shared identifiers, or give each cell its own namespace and resample
#: within it because the cells genuinely are independent.
DESIGN_RESAMPLING_UNITS = ("master_trial", "independent_namespace")

#: The deliberately wrong procedure, kept as a named control.
#:
#: Resampling within each cell independently when the cells share master trials
#: is the mistake the plan calls out: it treats common-random-number cells as
#: independent regression observations.  It is implemented here so that the
#: difference between it and the paired procedure can be *measured* — a control
#: that does not exist cannot discriminate against anything — and every report
#: it produces is flagged ``is_control`` and carries the reason it is wrong.
CONTROL_RESAMPLING_UNIT = "cell_independent_control"

#: Everything :func:`resample_exponent` accepts.
RESAMPLING_UNITS = DESIGN_RESAMPLING_UNITS + (CONTROL_RESAMPLING_UNIT,)

#: What an exponent from this module does not license anybody to say.
ANALYSIS_NON_CLAIMS = (
    "a statistically adequate exponent is not a Born-rule result, a detector "
    "click, an outcome, a unique actuality, or a two-channel result",
    "it is the coupling dependence of a STOPPING RULE this study imposed on a "
    "declared model bath, over one frozen population and one frozen pulse",
    "an exponent near two is not by itself support for the width-times-rate "
    "mechanism: that attribution requires the central-clock and width-only "
    "controls to FAIL to reproduce it, which is compare.py's decision and not "
    "this module's",
    "no cell is ever dropped here, so a 'no valid exponent' verdict is a "
    "result about the whole declared sweep and never about a subset",
    "intervals are Wald intervals from the expected information at the "
    "maximum; they are asymptotic, and a design with few events in every cell "
    "should be read through the resampling report instead",
)


# ---------------------------------------------------------------------------
# The chi-square upper tail, written here rather than imported
# ---------------------------------------------------------------------------
# This package depends on NumPy and nothing else, and both frozen tests below
# need a p-value.  The regularized upper incomplete gamma function is the whole
# requirement, and it is short enough to write correctly and check against two
# closed forms: one degree of freedom is erfc(sqrt(x/2)) and two degrees of
# freedom is exp(-x/2), and ``verify.py`` holds this function to both.
_GAMMA_SERIES_ITERATIONS = 1000
_GAMMA_EPSILON = 1e-16
_GAMMA_TINY = 1e-300


def _gamma_p_series(a: float, x: float) -> float:
    """Regularized lower incomplete gamma by its series.  ``x < a + 1``."""
    term = 1.0 / a
    total = term
    n = 0
    while n < _GAMMA_SERIES_ITERATIONS:
        n += 1
        term *= x / (a + n)
        total += term
        if abs(term) < abs(total) * _GAMMA_EPSILON:
            break
    return total * math.exp(-x + a * math.log(x) - math.lgamma(a))


def _gamma_q_continued_fraction(a: float, x: float) -> float:
    """Regularized upper incomplete gamma by Lentz's method.  ``x >= a + 1``."""
    b = x + 1.0 - a
    c = 1.0 / _GAMMA_TINY
    d = 1.0 / b
    h = d
    for i in range(1, _GAMMA_SERIES_ITERATIONS + 1):
        an = -i * (i - a)
        b += 2.0
        d = an * d + b
        if abs(d) < _GAMMA_TINY:
            d = _GAMMA_TINY
        c = b + an / c
        if abs(c) < _GAMMA_TINY:
            c = _GAMMA_TINY
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < _GAMMA_EPSILON:
            break
    return h * math.exp(-x + a * math.log(x) - math.lgamma(a))


def chi_square_upper_tail(statistic, df) -> float:
    """Return ``P(X > statistic)`` for a chi-square with ``df`` degrees.

    The p-value both frozen tests in this module are decided by.  Written here
    because this package's only dependency is NumPy; it is checked in
    ``verify.py`` against the two closed forms it has —
    ``erfc(sqrt(x/2))`` at one degree of freedom and ``exp(-x/2)`` at two — so
    the implementation is held to something other than itself.

    A statistic of exactly zero returns one; a negative one is refused rather
    than clamped, because a negative deviance or a negative likelihood ratio
    means the fit that produced it went wrong and silently reporting ``p = 1``
    would hide that.
    """
    x = require_finite(statistic, "statistic")
    k = require_positive_int(df, "df", minimum=1)
    if x < 0.0:
        raise ValueError(
            f"statistic must be non-negative, got {x!r}; a negative deviance "
            f"or likelihood ratio means the fit that produced it failed, and "
            f"returning a p-value of one would hide that")
    if x == 0.0:
        return 1.0
    a, half = 0.5 * k, 0.5 * x
    if half < a + 1.0:
        return float(min(1.0, max(0.0, 1.0 - _gamma_p_series(a, half))))
    return float(min(1.0, max(0.0, _gamma_q_continued_fraction(a, half))))


# ---------------------------------------------------------------------------
# Master-trial identity
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class MasterTrial:
    """Which physical path a trial is: a namespace and a number, together.

    An earlier revision identified a master trial by a non-negative integer
    alone.  An independent review pointed out what that cannot say: two runs
    from physically unrelated datasets both number their rows ``0 .. n-1``, so
    ``paired_exponent_contrast`` accepted them as a paired common-random-number
    experiment and returned a contrast.  There was no argument in which the two
    namespaces could have been supplied, so the substitution was not merely
    unchecked — it was unrepresentable.

    ``namespace``
        A canonical token identifying *which stream hierarchy* the trial number
        indexes.  ``compare.crn_identity`` derives it from the frozen manifest
        fields that decide which Brownian path trial ``i`` is — the key schema,
        the stream and ensemble prefixes, the mesh identity, the population and
        the physical parameters — and deliberately **excludes** the peak
        coupling and the model label, because a coupling sweep varies the
        first and a paired full/width-only comparison varies the second while
        sharing one tree.
    ``trial``
        The master trial number inside that namespace.

    :attr:`key` is an injective encoding: the namespace is length-prefixed, so
    no namespace containing the separator can impersonate another
    ``(namespace, trial)`` pair.
    """

    namespace: str
    trial: int

    def __post_init__(self):
        namespace = self.namespace
        if isinstance(namespace, bool) or not isinstance(namespace, str):
            raise TypeError(f"namespace must be a string, got "
                            f"{type(namespace).__name__} {namespace!r}")
        if not namespace:
            raise ValueError(
                "namespace must not be empty; an empty provenance token would "
                "make every dataset the same dataset, which is the defect this "
                "type exists to close")
        if any(character in namespace for character in "\r\n"):
            raise ValueError(f"namespace must be one line, got {namespace!r}")
        trial = self.trial
        if isinstance(trial, (bool, np.bool_)):
            raise TypeError(f"trial must be an integer, got bool {trial!r}")
        if not isinstance(trial, (int, np.integer)):
            raise TypeError(f"trial must be an integer, got "
                            f"{type(trial).__name__} {trial!r}")
        if int(trial) < 0:
            raise ValueError(f"trial must be non-negative, got {int(trial)!r}")
        object.__setattr__(self, "namespace", namespace)
        object.__setattr__(self, "trial", int(trial))

    @property
    def key(self) -> str:
        """The injective canonical spelling of this identity.

        Length-prefixed, so ``("a/1", 2)`` and ``("a", "1/2")`` cannot collide:
        the first is ``3:a/1/2`` and the second is not constructible, but a
        namespace that *did* contain a slash still gets its own unambiguous
        encoding.
        """
        return f"{len(self.namespace)}:{self.namespace}/{self.trial}"


def master_trials(namespace, count):
    """Return ``count`` consecutive :class:`MasterTrial` identities.

    The ordinary way to name the trials of one sweep: every cell of a paired
    design holds this exact tuple, in this exact order.
    """
    if isinstance(namespace, bool) or not isinstance(namespace, str):
        raise TypeError(f"namespace must be a string, got "
                        f"{type(namespace).__name__} {namespace!r}")
    total = require_positive_int(count, "count", minimum=1)
    return tuple(MasterTrial(namespace, index) for index in range(total))


# ---------------------------------------------------------------------------
# The cells and the design
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class CouplingCell:
    """One frozen coupling value and the per-trial outcomes recorded at it.

    ``coupling``
        The peak coupling.  Strictly positive: the model is linear in
        ``log(K)``, and a zero-coupling control — which must produce no
        commitment at any positive dwell — has no logarithm and belongs beside
        the fit as a separate statement, not inside it.
    ``trial_ids``
        The master trial identities, in the order the run recorded them, each
        appearing once.  Every entry is a :class:`MasterTrial`, so a cell knows
        *which dataset* its trial numbers index; they are carried rather than
        summarized because they are what makes paired resampling possible: the
        same identity in two cells is the same Brownian path at two couplings.

        One cell holds one namespace.  A cell assembled from two datasets is
        not a coupling cell, it is two of them, and mixing them silently is how
        an unpaired comparison acquires the appearance of pairing.
    ``outcomes``
        One boolean per identity: did the channel commit in that trial.

    ``trials`` and ``committed`` are derived, never stored, for the reason the
    closed-ledger gate derives its aggregates: a stored count can disagree with
    the rows it summarizes, and then one of them is wrong and nothing says
    which.
    """

    coupling: float
    trial_ids: tuple
    outcomes: tuple

    def __post_init__(self):
        coupling = require_finite(self.coupling, "coupling")
        if not coupling > 0.0:
            raise ValueError(
                f"coupling must be strictly positive, got {coupling!r}; this "
                f"model is linear in log(K), so a zero-coupling control is "
                f"reported beside the fit rather than inside it")
        ids = []
        for position, value in enumerate(tuple(self.trial_ids)):
            if type(value) is not MasterTrial:
                raise TypeError(
                    f"trial_ids[{position}] must be a MasterTrial, got "
                    f"{type(value).__name__} {value!r}; a bare integer cannot "
                    f"say which dataset it numbers, and two unrelated runs "
                    f"both numbered 0..n-1 would pair")
            ids.append(value)
        namespaces = {identity.namespace for identity in ids}
        if len(namespaces) > 1:
            raise ValueError(
                f"this cell mixes {len(namespaces)} namespaces "
                f"{sorted(namespaces)}; one coupling cell is one dataset, and "
                f"a mixed cell is two of them wearing one label")
        outcomes = []
        for position, value in enumerate(tuple(self.outcomes)):
            if not isinstance(value, (bool, np.bool_)):
                raise TypeError(
                    f"outcomes[{position}] must be a bool, got "
                    f"{type(value).__name__} {value!r}; a numeric flag is "
                    f"refused rather than coerced, because 2 would become "
                    f"True and a count would become an outcome")
            outcomes.append(bool(value))
        if len(ids) != len(outcomes):
            raise ValueError(f"this cell has {len(ids)} trial identifiers and "
                             f"{len(outcomes)} outcomes; every trial carries "
                             f"exactly one")
        if not ids:
            raise ValueError("a coupling cell holds at least one trial; an "
                             "empty cell contributes nothing and would make "
                             "the sweep's cell count a fiction")
        if len({identity.key for identity in ids}) != len(ids):
            raise ValueError("trial identities repeat within this cell; a "
                             "master trial appears once per coupling, and a "
                             "repeat would weight one Brownian path twice")
        object.__setattr__(self, "coupling", coupling)
        object.__setattr__(self, "trial_ids", tuple(ids))
        object.__setattr__(self, "outcomes", tuple(outcomes))

    @property
    def namespace(self) -> str:
        """The one dataset this cell's trials come from."""
        return self.trial_ids[0].namespace

    @property
    def keys(self) -> tuple:
        """The canonical identity strings, in recorded order."""
        return tuple(identity.key for identity in self.trial_ids)

    @property
    def trials(self) -> int:
        """How many trials this cell holds."""
        return len(self.outcomes)

    @property
    def committed(self) -> int:
        """How many of them committed."""
        return sum(1 for value in self.outcomes if value)

    @property
    def is_zero_event(self) -> bool:
        """True when nothing committed.  Retained, never dropped."""
        return self.committed == 0

    @property
    def is_all_event(self) -> bool:
        """True when everything committed.  Retained, never dropped."""
        return self.committed == self.trials


@dataclass(frozen=True)
class ScalingDesign:
    """The frozen set of coupling cells, and the resampling unit they imply.

    ``cells``
        At least :data:`MIN_CELLS` of them, with strictly increasing and
        distinct couplings.  Order is fixed here so that two calls on the same
        run cannot produce two different tables.
    ``resampling_unit``
        One of :data:`DESIGN_RESAMPLING_UNITS`, and **checked against the
        identifiers**:

        * ``master_trial`` requires every cell to hold the *identical* tuple of
          trial identifiers.  That is what a common-random-number sweep looks
          like, and it is what makes resampling the shared identifiers a
          correct bootstrap;
        * ``independent_namespace`` requires the cells' identifier sets to be
          pairwise *disjoint*.  Then, and only then, resampling within each
          cell is correct.

        A design that claims one and looks like the other is refused.  The
        declaration is not a hint.
    """

    cells: tuple
    resampling_unit: str

    def __post_init__(self):
        cells = tuple(self.cells)
        for position, cell in enumerate(cells):
            require_instance(cell, CouplingCell, f"cells[{position}]")
        if len(cells) < MIN_CELLS:
            raise ValueError(
                f"a scaling design needs at least {MIN_CELLS} coupling cells, "
                f"got {len(cells)}; with fewer, the frozen curvature model has "
                f"no residual degree of freedom and its lack-of-fit test is "
                f"vacuous")
        couplings = [cell.coupling for cell in cells]
        if any(b <= a for a, b in zip(couplings, couplings[1:])):
            raise ValueError(
                f"couplings must be strictly increasing and distinct, got "
                f"{couplings}; a repeated coupling is one cell recorded twice "
                f"and a shuffled sweep is two different tables")
        unit = self.resampling_unit
        if isinstance(unit, bool) or not isinstance(unit, str):
            raise TypeError(f"resampling_unit must be a string, got "
                            f"{type(unit).__name__} {unit!r}")
        if unit not in DESIGN_RESAMPLING_UNITS:
            raise ValueError(
                f"resampling_unit must be one of "
                f"{list(DESIGN_RESAMPLING_UNITS)}, got {unit!r}; "
                f"{CONTROL_RESAMPLING_UNIT!r} is a control a caller may ask "
                f"resample_exponent for, never a property a design may claim")
        if unit == "master_trial":
            # Compared as identity objects rather than as their encodings:
            # ``MasterTrial`` is a frozen dataclass, so tuple equality *is*
            # identity equality, and building a string per trial per cell would
            # dominate the bootstrap's cost for no additional guarantee.
            first = cells[0].trial_ids
            for position, cell in enumerate(cells[1:], start=1):
                if cell.trial_ids != first:
                    raise ValueError(
                        f"cells[0] and cells[{position}] hold different master "
                        f"trial identities, so this sweep is not paired and "
                        f"{unit!r} would resample identities that do not "
                        f"appear in every cell.  Identity is namespace *and* "
                        f"number: a substituted dataset with the same local "
                        f"numbering is refused here, and used not to be")
        else:
            seen, namespaces = set(), set()
            for position, cell in enumerate(cells):
                keys = set(cell.trial_ids)
                shared = seen & keys
                if shared:
                    raise ValueError(
                        f"cells[{position}] shares {len(shared)} trial "
                        f"identities with an earlier cell, so these cells are "
                        f"not independent namespaces and resampling within "
                        f"each would treat common random numbers as "
                        f"independent")
                if cell.namespace in namespaces:
                    raise ValueError(
                        f"cells[{position}] reuses the namespace "
                        f"{cell.namespace!r}; independent coupling namespaces "
                        f"means one namespace per cell, and two cells drawn "
                        f"from one stream hierarchy are not independent "
                        f"however their trials happen to be numbered")
                seen |= keys
                namespaces.add(cell.namespace)
        object.__setattr__(self, "cells", cells)
        object.__setattr__(self, "resampling_unit", unit)

    @property
    def couplings(self) -> tuple:
        """The frozen coupling grid, low to high."""
        return tuple(cell.coupling for cell in self.cells)

    @property
    def namespaces(self) -> tuple:
        """The dataset each cell came from, in design order."""
        return tuple(cell.namespace for cell in self.cells)

    @property
    def digest(self) -> str:
        """A digest of the exact observations this design is, cell by cell.

        Ordered couplings, canonical master-trial identities *and every
        outcome* — not the aggregate counts.  An independent review pointed out
        that comparing ``(committed, trials)`` per cell cannot distinguish two
        different master-trial outcome matrices with the same totals, and those
        have different paired bootstraps.  So the digest is over the matrix.

        It is what binds a fit, a resampling record, a contrast and a control
        to the data they were computed from, everywhere downstream.
        """
        parts = []
        for cell in self.cells:
            outcomes = "".join("1" if value else "0" for value in cell.outcomes)
            keys = ",".join(identity.key for identity in cell.trial_ids)
            parts.append(f"{cell.coupling!r}|{len(keys)}:{keys}|{outcomes}")
        payload = f"{self.resampling_unit}#" + "#".join(parts)
        return hashlib.sha256(payload.encode("utf-8")).hexdigest()

    @property
    def paired_keys(self) -> tuple:
        """The shared master-trial identities of a ``master_trial`` design.

        Raises for an ``independent_namespace`` design, which by construction
        has no shared identities.  Used by the paired bootstrap and by
        ``compare.paired_exponent_contrast``, both of which need to prove that
        two designs are the *same* trials before differencing them.
        """
        if self.resampling_unit != "master_trial":
            raise ValueError(
                f"this design declares {self.resampling_unit!r}; only a "
                f"master_trial design has identities shared by every cell")
        return self.cells[0].trial_ids

    @property
    def log_couplings(self) -> tuple:
        """``log(K)`` for each cell — the model's only covariate."""
        return tuple(math.log(cell.coupling) for cell in self.cells)

    @property
    def counts(self) -> tuple:
        """``(committed, trials)`` per cell, in design order."""
        return tuple((cell.committed, cell.trials) for cell in self.cells)


# ---------------------------------------------------------------------------
# The likelihood
# ---------------------------------------------------------------------------
def _hazard_terms(eta):
    """Return ``(H, g, log_q)`` for a linear predictor, with the limits right.

    ``H = exp(eta)`` is the cumulative hazard, ``q = 1 - exp(-H)`` the
    commitment probability, and ``g(H) = H exp(-H) / (1 - exp(-H))`` the factor
    that appears in both the score and the expected information.

    Two limits have to be handled rather than left to floating point:
    ``g(0) = 1``, which the ratio spells ``0/0``; and a large ``H``, where
    ``exp(-H)`` underflows and both ``g`` and ``log(q)`` are exactly right at
    ``0`` — no special case needed there, only the reassurance that none is.
    Small ``H`` is safe as written because ``-expm1(-H)`` is accurate to full
    precision where ``1 - exp(-H)`` would not be.
    """
    H = np.exp(eta)
    minus_expm1 = -np.expm1(-H)           # == 1 - exp(-H), accurate at small H
    with np.errstate(divide="ignore", invalid="ignore"):
        g = np.where(H > 0.0, H * np.exp(-H) / minus_expm1, 1.0)
        log_q = np.where(minus_expm1 > 0.0, np.log(minus_expm1), -np.inf)
    return H, g, log_q


def _log_likelihood(eta, committed, trials) -> float:
    """Binomial cloglog log-likelihood, up to the constant binomial factor."""
    H, _, log_q = _hazard_terms(eta)
    unresolved = trials - committed
    terms = np.where(committed > 0, committed * log_q, 0.0) - unresolved * H
    if not np.all(np.isfinite(terms)):
        return -np.inf
    return float(np.sum(terms))


def _saturated_log_likelihood(committed, trials) -> float:
    """The best any model could do: one free probability per cell."""
    total = 0.0
    for y, n in zip(committed, trials):
        if 0 < y < n:
            total += y * math.log(y / n) + (n - y) * math.log((n - y) / n)
    return total


def _fisher_scoring(design_matrix, committed, trials):
    """Fit by Fisher scoring with step halving; return a status dictionary.

    Fisher scoring rather than Newton because the expected information of a
    binomial GLM is positive semi-definite by construction, so the step is
    always an ascent direction where the information is non-singular.  Step
    halving on the log-likelihood then makes divergence a *reported* outcome
    rather than an oscillation that eventually returns a number.
    """
    X = np.asarray(design_matrix, dtype=float)
    y = np.asarray(committed, dtype=float)
    n = np.asarray(trials, dtype=float)
    beta = np.zeros(X.shape[1], dtype=float)
    ll = _log_likelihood(X @ beta, y, n)
    converged, iterations, runaway = False, 0, False

    for iterations in range(1, MAX_ITERATIONS + 1):
        eta = X @ beta
        H, g, _ = _hazard_terms(eta)
        score = y * g - (n - y) * H
        gradient = X.T @ score
        weight = n * H * g
        information = X.T @ (weight[:, None] * X)
        if not np.all(np.isfinite(information)) or not np.all(
                np.isfinite(gradient)):
            break
        try:
            step = np.linalg.solve(information, gradient)
        except np.linalg.LinAlgError:
            break
        if float(gradient @ step) < GRADIENT_TOLERANCE:
            converged = True
            break
        scale, improved = 1.0, False
        for _ in range(40):
            candidate = beta + scale * step
            candidate_ll = _log_likelihood(X @ candidate, y, n)
            if candidate_ll > ll:
                beta, ll, improved = candidate, candidate_ll, True
                break
            scale *= 0.5
        if not improved:
            converged = True          # no ascent direction left: at the top
            break
        if float(np.max(np.abs(X @ beta))) > MAX_ABS_LINEAR_PREDICTOR:
            runaway = True
            break

    eta = X @ beta
    H, g, _ = _hazard_terms(eta)
    weight = n * H * g
    information = X.T @ (weight[:, None] * X)
    return {"beta": beta, "log_likelihood": ll, "information": information,
            "converged": converged, "iterations": iterations,
            "runaway": runaway or float(np.max(np.abs(eta)))
            > MAX_ABS_LINEAR_PREDICTOR,
            "fitted_hazard": H}


def _reciprocal_condition(matrix) -> float:
    """``1 / cond2(matrix)``, and zero for a matrix that is not finite."""
    values = np.asarray(matrix, dtype=float)
    if not np.all(np.isfinite(values)):
        return 0.0
    singular = np.linalg.svd(values, compute_uv=False)
    if singular.size == 0 or singular[0] <= 0.0:
        return 0.0
    return float(singular[-1] / singular[0])


@dataclass(frozen=True)
class ClusterScoreTest:
    """The cluster-calibrated general lack-of-fit test, and what it was run on.

    Why this exists
    ---------------
    The residual deviance of a binomial GLM is compared against a chi-square
    whose degrees of freedom assume the cells are independent observations.  A
    ``master_trial`` sweep violates that by construction, and the distortion is
    not small.  An independent review simulated five hundred exact power-law
    datasets at four hundred trials over six frozen cells:

    ======================  ==========================  ==================
    design                  marginal lack-of-fit        nominal
    ======================  ==========================  ==================
    shared master trials    0.6 %                       5 %
    independent namespaces  4.4 %                       5 %
    ======================  ==========================  ==================

    Five times too conservative under exactly the dependence this experiment
    has.  A test that under-rejects is not a safe extra refusal: it *misses*
    non-power responses, and the review found a ``bend = 0.5`` curved fixture
    passing because of it.

    What this is
    ------------
    A **cluster-robust generalized score test** on the lack-of-fit contrasts,
    with a **Rademacher multiplier bootstrap** p-value.

    Write ``s_t`` for master trial ``t``'s score contribution to the *saturated*
    model — one entry per cell, ``y * g(H) - (1 - y) * H`` at the fitted null.
    Summing over trials gives ``S``, and the maximum-likelihood conditions make
    ``S`` exactly orthogonal to the two columns of the power model's design
    matrix.  All of ``S`` therefore lies in the ``cells - 2`` dimensional
    lack-of-fit space; let ``C`` be an orthonormal basis of it.  The statistic
    is

    ``T = (C' S)' (C' V C)^-1 (C' S)``,  ``V = sum_t s_t s_t'``

    where the sum over ``t`` is over **clusters**, so ``V`` carries the
    master-trial dependence exactly rather than assuming it away.  The
    multiplier bootstrap re-signs whole clusters, ``S* = sum_t w_t s_t`` with
    ``w`` uniform on ``{-1, +1}``, which is centred under the fitted null by
    construction.

    A design with independent namespaces has one observation per cluster, so
    the same formula reduces to the ordinary independent score test — one test
    covers both units.

    Measured behaviour, over four hundred datasets each
    --------------------------------------------------
    ==========================  =============  ====================
    ensemble                    rejection      what it replaces
    ==========================  =============  ====================
    null, shared master trials  3.3 %          0.6 %
    null, independent           4.8 %          4.4 %
    ``bend = 0.5`` curved       82 %           passed
    ``bend = 1.2`` curved       100 %          rejected
    ==========================  =============  ====================

    Still mildly conservative under the strongest dependence tested, and five
    times closer to the declared level than the statistic it replaces, with
    real power against the alternative that statistic missed.

    Nothing here is tunable.  :data:`LACK_OF_FIT_REPLICATES` is frozen, the
    level is :data:`LACK_OF_FIT_ALPHA`, and the seed is derived from
    ``design_digest`` — so there is no argument through which a caller could
    re-roll a p-value it did not like.
    """

    kind: str
    statistic: float
    df: int
    sandwich_p: float
    bootstrap_p: float
    replicates: int
    seed: int
    clusters: int
    unit: str
    design_digest: str
    adequate: bool
    inadequate_reason: str = ""
    seal: object = None

    def __post_init__(self):
        _require_seal(self.seal, "ClusterScoreTest")
        if self.kind not in CLUSTER_TEST_KINDS:
            raise ValueError(f"kind must be one of "
                             f"{list(CLUSTER_TEST_KINDS)}, got {self.kind!r}")
        if not isinstance(self.inadequate_reason, str):
            raise TypeError("inadequate_reason must be a string")
        if bool(self.inadequate_reason) == bool(self.adequate):
            raise ValueError(
                f"adequate={self.adequate} with reason "
                f"{self.inadequate_reason!r}: an adequate test states no "
                f"reason and an inadequate one states exactly why, so that a "
                f"calibration refusal is never mistaken for a pass")
        if self.unit not in DESIGN_RESAMPLING_UNITS:
            raise ValueError(f"unit must be one of "
                             f"{list(DESIGN_RESAMPLING_UNITS)}, got "
                             f"{self.unit!r}")
        require_positive_int(self.df, "df", minimum=1)
        require_positive_int(self.replicates, "replicates", minimum=1)
        require_positive_int(self.clusters, "clusters", minimum=1)
        if not isinstance(self.adequate, bool):
            raise TypeError("adequate must be a bool")
        object.__setattr__(self, "seal", None)

    @property
    def rejects(self) -> bool:
        """Does the authoritative p-value fall below the frozen level?

        ``False`` when the test is inadequate: an inadequate test has not
        failed to reject, it has not run, and ``calibrate_exponent`` refuses
        under ``uncalibrated_fit`` rather than reading silence as a pass.
        """
        return bool(self.adequate
                    and math.isfinite(self.bootstrap_p)
                    and self.bootstrap_p < LACK_OF_FIT_ALPHA)


def _cluster_scores(design, amplitude, exponent):
    """Return the per-cluster saturated-model score matrix at a fitted null.

    Rows are clusters, columns are cells.  For a ``master_trial`` design one
    cluster is one master trial and its row has an entry in every cell —
    which is precisely the dependence the marginal deviance ignores.  For
    independent namespaces one cluster is one trial in one cell, so its row is
    a one-hot, and the sandwich reduces to the independent case.
    """
    cells = design.cells
    x = np.array(design.log_couplings, dtype=float)
    eta = amplitude + exponent * x
    H, g, _ = _hazard_terms(eta)
    if design.resampling_unit == "master_trial":
        y = np.array([[1.0 if value else 0.0 for value in cell.outcomes]
                      for cell in cells], dtype=float).T
        return y * g[None, :] - (1.0 - y) * H[None, :]
    blocks = []
    for index, cell in enumerate(cells):
        block = np.zeros((cell.trials, len(cells)), dtype=float)
        y = np.array([1.0 if value else 0.0 for value in cell.outcomes],
                     dtype=float)
        block[:, index] = y * g[index] - (1.0 - y) * H[index]
        blocks.append(block)
    return np.vstack(blocks)


def _nuisance_orthogonal_basis(design, report, columns=None):
    """Return the contrast basis the *estimated* nuisance parameters leave.

    The correction an independent review asked for, in full.

    Write ``X = [1, log K]`` for the power model's design matrix, ``eta`` for
    the fitted linear predictor, and ``W = diag(n_i H_i g(H_i))`` for the
    expected information in ``eta`` — the same weight :func:`_fisher_scoring`
    uses.  Let ``C`` be any orthonormal basis of the Euclidean complement of
    ``X``, which is the ``cells - 2`` dimensional lack-of-fit space.

    The naive statistic projects each cluster's saturated-model score into
    ``C``.  That is wrong once the two nuisance parameters have been
    *estimated*.  The maximum-likelihood conditions make the **total** score
    orthogonal to ``X``, so the numerator is unaffected — but an individual
    cluster's contribution is not orthogonal to ``X``, and it is the clusters
    that form the sandwich meat and that the multiplier bootstrap re-signs.
    Re-signing raw ``C' s_t`` therefore does not reproduce the score's
    distribution after fitting.

    The efficient score for the contrasts, adjusted for the estimated
    nuisance, replaces ``s_t`` by

    ``e_t = C' s_t - C' W X (X' W X)^-1 X' s_t``

    since the cross-information between the saturated parameters and ``beta``
    is ``W X`` and the nuisance information is ``X' W X``.  Written as a single
    basis so that ``e_t^T = s_t^T B``:

    ``B = C - X (X' W X)^-1 X' W C``

    which satisfies ``X' W B = 0`` exactly — the basis is **W-orthogonal to the
    nuisance directions**, which is the property the naive ``C`` lacks and the
    property ``verify.py`` checks to machine precision.

    Two consequences worth stating because they are what make this safe to
    change:

    * the **numerator is unchanged**.  ``sum_t e_t = C' S - C' W X (X'WX)^-1
      X' S = C' S``, because ``X' S = 0`` at the maximum.  Only the meat and
      the bootstrap draws move.
    * measured over frozen null ensembles, the size moves from mildly
      conservative to nominal and the power *rises*: see the tables in
      ``verify.check_cluster_lack_of_fit``.

    Returns ``None`` when the nuisance information is not invertible, which is
    a calibration refusal rather than a p-value.
    """
    cells = design.cells
    x = np.array(design.log_couplings, dtype=float)
    matrix = np.column_stack([np.ones_like(x), x])
    if columns is None:
        complement = np.linalg.svd(matrix,
                                   full_matrices=True)[0][:, matrix.shape[1]:]
    else:
        complement = np.asarray(columns, dtype=float)
    eta = report.amplitude + report.exponent * x
    hazard, factor, _ = _hazard_terms(eta)
    trials = np.array([cell.trials for cell in cells], dtype=float)
    weights = trials * hazard * factor
    if not np.all(np.isfinite(weights)):
        return None
    weighted = matrix * weights[:, None]             # W X
    information = matrix.T @ weighted                # X' W X
    if _reciprocal_condition(information) <= MIN_RECIPROCAL_CONDITION:
        return None
    adjustment = matrix @ np.linalg.solve(information,
                                          weighted.T @ complement)
    basis = complement - adjustment
    if not np.all(np.isfinite(basis)):
        return None
    return basis


def curvature_column(design):
    """The frozen curvature extension's covariate, centred on its own grid.

    ``(log K - mean log K)^2``.  Centring makes the column near-orthogonal to
    the intercept before any adjustment, and the mean is a property of the
    frozen coupling grid rather than of any outcome, so the extension cannot be
    shifted after seeing the counts.
    """
    require_instance(design, ScalingDesign, "design")
    x = np.array(design.log_couplings, dtype=float)
    centred = x - float(np.mean(x))
    return centred * centred


def cluster_curvature(design, report) -> ClusterScoreTest:
    """The paired curvature decision: a one-dimensional score test.

    The same construction as :func:`cluster_lack_of_fit` with a single contrast
    column — the curvature covariate, W-residualized against ``[1, log K]``.

    It replaces a percentile interval.  An independent review held one dataset
    fixed and changed only the public bootstrap seed, and the curvature
    interval excluded zero for **58 of 100 seeds**; on a thousand interior
    exact-null datasets the same interval excluded zero 6.7 percent of the time
    at a nominal 5.  A percentile interval of bootstrap point estimates is not
    a test, it is rerollable, and it was carrying an authority it could not
    support.

    This is a test: one statistic, one frozen replicate count, one seed derived
    from the design digest, and the same degeneracy refusal.  The percentile
    interval survives on :class:`ResamplingReport` as a **descriptive** figure.
    """
    return _multiplier_score_test(design, report, "curvature",
                                  curvature_column(design)[:, None])


def cluster_lack_of_fit(design, report) -> ClusterScoreTest:
    """Run the cluster-calibrated lack-of-fit test on one fitted design.

    Takes no replicate count and no seed: both are frozen, the seed derived
    from the design's own digest, so the test cannot be re-rolled.  See
    :class:`ClusterScoreTest` for what it computes and how it behaves.
    """
    return _multiplier_score_test(design, report, "lack_of_fit", None)


def _multiplier_score_test(design, report, kind, columns) -> ClusterScoreTest:
    """The nuisance-adjusted cluster score test, over any contrast columns.

    ``columns is None`` means the full ``cells - 2`` lack-of-fit complement;
    a single column means a one-dimensional test of that contrast.  Everything
    else — the weight, the nuisance residualization, the cluster meat, the
    degeneracy refusals, the frozen replicate count and the digest-derived
    seed — is shared, so the two decisions cannot drift apart.
    """
    require_instance(design, ScalingDesign, "design")
    require_instance(report, ExponentReport, "report")
    digest = design.digest
    df = 1 if columns is not None else len(design.cells) - 2
    # One seed per (design, kind): derived, so neither test can be re-rolled
    # and neither reuses the other's signs.
    seed = int(digest[:15], 16) ^ (0 if kind == "lack_of_fit" else 0x5CA1AB1E)

    def refused(reason_clusters, why):
        return ClusterScoreTest(
            kind=kind, statistic=float("nan"), df=max(df, 1),
            sandwich_p=float("nan"), bootstrap_p=float("nan"),
            replicates=LACK_OF_FIT_REPLICATES, seed=seed,
            clusters=max(reason_clusters, 1), unit=design.resampling_unit,
            design_digest=digest, adequate=False, inadequate_reason=why,
            seal=_SEAL)

    if df < 1 or not (math.isfinite(report.amplitude)
                      and math.isfinite(report.exponent)):
        return refused(1, "the fit has no finite parameters, or the design "
                          "has no residual degrees of freedom")

    scores = _cluster_scores(design, report.amplitude, report.exponent)
    clusters = scores.shape[0]
    if clusters < MIN_LACK_OF_FIT_CLUSTERS or not np.all(np.isfinite(scores)):
        return refused(clusters, "too few clusters for a cluster-robust "
                                 "covariance")
    # A cell in which every trial did the same thing contributes a *constant*
    # per-cluster score: it carries a systematic offset and no sampling
    # variation at all.  The sandwich cannot represent that -- it divides a
    # deterministic offset by a variance the multiplier bootstrap invents --
    # and an independent review found the consequence: an exactly specified
    # power law with three saturated cells was rejected 43 times in 120.
    #
    # The cells stay in the likelihood, in the report and in every count.  What
    # is refused is the *calibration*, which is a different thing from deleting
    # data.
    degenerate = [index for index, cell in enumerate(design.cells)
                  if cell.is_zero_event or cell.is_all_event]
    if degenerate:
        return refused(clusters,
                       f"cells {degenerate} have no within-cell variation, so "
                       f"their cluster scores are deterministic and the "
                       f"bootstrap has no variance to calibrate against")

    basis = _nuisance_orthogonal_basis(design, report, columns)
    if basis is None:
        return refused(clusters, "the nuisance-orthogonal contrast basis is "
                                 "not computable at this fit")
    projected = scores @ basis                       # clusters x (cells - 2)
    total = projected.sum(axis=0)
    meat = projected.T @ projected * (clusters / (clusters - 1.0))
    if _reciprocal_condition(meat) <= MIN_RECIPROCAL_CONDITION:
        return refused(clusters, "the cluster-robust covariance of the "
                                 "lack-of-fit contrasts is singular")

    inverse = np.linalg.inv(meat)
    statistic = float(total @ inverse @ total)
    if not math.isfinite(statistic) or statistic < 0.0:
        return refused(clusters)
    sandwich_p = chi_square_upper_tail(statistic, df)

    rng = np.random.default_rng(seed)
    multipliers = rng.integers(0, 2, size=(LACK_OF_FIT_REPLICATES,
                                           clusters)) * 2.0 - 1.0
    resampled = multipliers @ projected
    draws = np.einsum("bi,ij,bj->b", resampled, inverse, resampled)
    bootstrap_p = float((1.0 + np.count_nonzero(draws >= statistic))
                        / (1.0 + LACK_OF_FIT_REPLICATES))
    return ClusterScoreTest(
        kind=kind, statistic=statistic, df=df, sandwich_p=float(sandwich_p),
        bootstrap_p=bootstrap_p, replicates=LACK_OF_FIT_REPLICATES, seed=seed,
        clusters=clusters, unit=design.resampling_unit, design_digest=digest,
        adequate=True, inadequate_reason="", seal=_SEAL)


@dataclass(frozen=True)
class ExponentReport:
    """The structured outcome of one exponent fit — valid or not.

    ``verdict``
        One of :data:`VERDICTS`.  ``no_valid_exponent`` is a result and comes
        back with everything else populated; it is never an exception and never
        a silently narrowed sweep.
    ``reasons``
        A tuple drawn from :data:`NO_EXPONENT_REASONS`, empty exactly when the
        verdict is ``valid_exponent``.
    ``exponent``, ``exponent_stderr``, ``exponent_interval``
        The unconstrained coupling exponent, its standard error from the
        expected information, and the Wald interval at
        :data:`CONFIDENCE_LEVEL`.  Reported whatever the verdict, and whatever
        the value: nothing here pulls an estimate towards two.
    ``amplitude``, ``amplitude_stderr``
        The intercept of the same fit, on the ``log(-log(1-q))`` scale.
    ``deviance``, ``deviance_df``, ``deviance_p``
        Residual deviance against the saturated model, its degrees of freedom
        and the frozen upper-tail p-value.
    ``curvature``, ``curvature_stderr``, ``curvature_lr``, ``curvature_p``
        The frozen curvature extension's coefficient and its likelihood-ratio
        test on one degree of freedom.
    ``cells_declared``, ``cells_fitted``
        Always equal.  They are both reported because "the fit used every cell"
        is the claim most worth being able to check mechanically.
    ``zero_event_cells``, ``all_event_cells``
        How many cells sat at the boundary and were kept anyway.
    ``rows``
        Per cell: coupling, counts, observed probability, observed and fitted
        cumulative hazard.  The raw counts travel with the fit.
    """

    verdict: str
    reasons: tuple
    exponent: float
    exponent_stderr: float
    exponent_interval: tuple
    amplitude: float
    amplitude_stderr: float
    deviance: float
    deviance_df: int
    deviance_p: float
    curvature: float
    curvature_stderr: float
    curvature_lr: float
    curvature_p: float
    cells_declared: int
    cells_fitted: int
    zero_event_cells: int
    all_event_cells: int
    converged: bool
    reciprocal_condition: float
    rows: tuple
    design_digest: str = ""

    def __post_init__(self):
        if not isinstance(self.design_digest, str):
            raise TypeError(f"design_digest must be a string, got "
                            f"{type(self.design_digest).__name__}")
        if self.verdict not in VERDICTS:
            raise ValueError(f"verdict must be one of {list(VERDICTS)}, got "
                             f"{self.verdict!r}")
        reasons = tuple(self.reasons)
        unknown = [r for r in reasons if r not in NO_EXPONENT_REASONS]
        if unknown:
            raise ValueError(
                f"reasons {unknown} are outside the closed vocabulary "
                f"{list(NO_EXPONENT_REASONS)}; a consumer branches on these, "
                f"so a new one is a deliberate edit here rather than a string "
                f"invented at a call site")
        if (self.verdict == "valid_exponent") != (not reasons):
            raise ValueError(
                f"verdict {self.verdict!r} with reasons {list(reasons)}: a "
                f"valid exponent has no refusal reasons and a refusal has at "
                f"least one")
        if self.cells_declared != self.cells_fitted:
            raise ValueError(
                f"{self.cells_declared} cells were declared and "
                f"{self.cells_fitted} entered the likelihood; this module "
                f"never drops a cell, so a difference is a defect and not a "
                f"narrowed sweep")
        object.__setattr__(self, "reasons", reasons)
        object.__setattr__(self, "rows", tuple(self.rows))

    @property
    def is_valid(self) -> bool:
        """``True`` exactly when the verdict is ``valid_exponent``."""
        return self.verdict == "valid_exponent"

    def interval_contains(self, value) -> bool:
        """Is ``value`` inside the reported exponent interval?

        The one question the causal decision rule in ``compare.py`` asks of a
        control's fit — "is this still compatible with quadratic scaling?" — so
        it is answered here, once, rather than by each caller re-deriving the
        endpoints.  A non-finite interval contains everything, which is the
        honest reading of an estimate that did not converge.
        """
        target = require_finite(value, "value")
        low, high = self.exponent_interval
        if not (math.isfinite(low) and math.isfinite(high)):
            return True
        return low <= target <= high


def fit_exponent(design) -> ExponentReport:
    """Fit :data:`LINK_FORMULA` to every cell of ``design`` and report.

    The whole procedure, in order:

    1. every declared cell enters the likelihood — including zero-event and
       all-event cells, with no transformation, no continuity correction, no
       deletion and no window selection;
    2. the two-parameter power model is fitted by Fisher scoring;
    3. the frozen curvature extension is fitted the same way and compared by
       likelihood ratio on one degree of freedom;
    4. the residual deviance is compared against its own degrees of freedom;
       and
    5. the verdict is ``valid_exponent`` only if the estimate is identified,
       converged, informative at the frozen cap, and refused by neither test.

    Failure at any of those points returns ``no_valid_exponent`` with the
    reason attached.  It never returns a fit over a subset, and there is no
    argument to this function that could ask it to.
    """
    require_instance(design, ScalingDesign, "design")
    cells = design.cells
    x = np.array(design.log_couplings, dtype=float)
    y = np.array([cell.committed for cell in cells], dtype=float)
    n = np.array([cell.trials for cell in cells], dtype=float)

    zero_event = sum(1 for cell in cells if cell.is_zero_event)
    all_event = sum(1 for cell in cells if cell.is_all_event)

    power_matrix = np.column_stack([np.ones_like(x), x])
    centred = x - float(np.mean(x))
    curve_matrix = np.column_stack([np.ones_like(x), x, centred * centred])

    design_condition = _reciprocal_condition(power_matrix)
    power = _fisher_scoring(power_matrix, y, n)
    curve = _fisher_scoring(curve_matrix, y, n)

    information = power["information"]
    info_condition = _reciprocal_condition(information)
    if info_condition > MIN_RECIPROCAL_CONDITION:
        covariance = np.linalg.inv(information)
        errors = np.sqrt(np.clip(np.diag(covariance), 0.0, None))
    else:
        errors = np.full(power_matrix.shape[1], float("inf"))

    amplitude, exponent = float(power["beta"][0]), float(power["beta"][1])
    amplitude_err, exponent_err = float(errors[0]), float(errors[1])
    half = CONFIDENCE_Z * exponent_err
    interval = (exponent - half, exponent + half)

    saturated = _saturated_log_likelihood(y, n)
    deviance = 2.0 * (saturated - power["log_likelihood"])
    deviance_df = len(cells) - power_matrix.shape[1]
    if math.isfinite(deviance) and deviance >= 0.0 and deviance_df >= 1:
        deviance_p = chi_square_upper_tail(deviance, deviance_df)
    else:
        deviance_p = float("nan")

    curve_information = curve["information"]
    curve_condition = _reciprocal_condition(curve_information)
    if curve_condition > MIN_RECIPROCAL_CONDITION:
        curve_errors = np.sqrt(np.clip(
            np.diag(np.linalg.inv(curve_information)), 0.0, None))
        curvature_err = float(curve_errors[2])
    else:
        curvature_err = float("inf")
    # The curvature extension is a *required* diagnostic, so it has to be
    # estimable before the power model it tests can be called valid.  An
    # independent review built a six-point grid spanning 1.000 to 1.010 with two
    # hundred thousand trials in every cell and exactly half events: the power
    # exponent was adequately informed, the curvature information was singular,
    # the likelihood ratio was zero, its p-value was one — and the report said
    # `valid_exponent`.  A test that cannot fail is not a test that passed.
    curvature = float(curve["beta"][2])
    curvature_identified = (
        curve_condition > MIN_RECIPROCAL_CONDITION
        and bool(curve["converged"]) and not curve["runaway"]
        and math.isfinite(curvature_err)
        and math.isfinite(curvature)
        and CONFIDENCE_Z * curvature_err <= MAX_EXPONENT_HALF_WIDTH)
    curvature_lr = 2.0 * (curve["log_likelihood"] - power["log_likelihood"])
    if math.isfinite(curvature_lr) and curvature_lr >= 0.0:
        curvature_p = chi_square_upper_tail(curvature_lr, 1)
    else:
        curvature_lr = float("nan")
        curvature_p = float("nan")

    fitted = power["fitted_hazard"]
    rows = []
    for index, cell in enumerate(cells):
        observed_q = cell.committed / cell.trials
        rows.append({
            "coupling": cell.coupling,
            "log_coupling": float(x[index]),
            "trials": cell.trials,
            "committed": cell.committed,
            "unresolved": cell.trials - cell.committed,
            "commitment_probability": observed_q,
            "observed_cumulative_hazard": (float("inf") if observed_q >= 1.0
                                           else -math.log1p(-observed_q)),
            "fitted_cumulative_hazard": float(fitted[index]),
            "zero_event": cell.is_zero_event,
            "all_event": cell.is_all_event,
        })

    reasons = []
    if zero_event == len(cells) or all_event == len(cells):
        reasons.append("saturated")
    if power["runaway"]:
        reasons.append("separation")
    if info_condition <= MIN_RECIPROCAL_CONDITION:
        reasons.append("flat")
    if design_condition <= MIN_RECIPROCAL_CONDITION:
        reasons.append("ill_conditioned")
    if not power["converged"]:
        reasons.append("non_convergence")
    if math.isfinite(deviance_p) and deviance_p < LACK_OF_FIT_ALPHA:
        reasons.append("lack_of_fit")
    if math.isfinite(curvature_p) and curvature_p < CURVATURE_ALPHA:
        reasons.append("curvature")
    if not curvature_identified:
        reasons.append("curvature_unidentified")
    if not (math.isfinite(half) and half <= MAX_EXPONENT_HALF_WIDTH):
        reasons.append("uninformative")
    # One canonical order, so two runs of the same fit report the same tuple.
    ordered = tuple(r for r in NO_EXPONENT_REASONS if r in set(reasons))

    return ExponentReport(
        verdict="no_valid_exponent" if ordered else "valid_exponent",
        reasons=ordered, exponent=exponent, exponent_stderr=exponent_err,
        exponent_interval=interval, amplitude=amplitude,
        amplitude_stderr=amplitude_err, deviance=float(deviance),
        deviance_df=deviance_df, deviance_p=float(deviance_p),
        curvature=curvature, curvature_stderr=curvature_err,
        curvature_lr=float(curvature_lr), curvature_p=float(curvature_p),
        cells_declared=len(cells), cells_fitted=len(cells),
        zero_event_cells=zero_event, all_event_cells=all_event,
        converged=bool(power["converged"]),
        reciprocal_condition=info_condition, rows=tuple(rows),
        design_digest=design.digest)


@dataclass(frozen=True)
class ResamplingReport:
    """Uncertainty on the exponent from resampling, and what was resampled.

    ``unit``
        One of :data:`RESAMPLING_UNITS`.
    ``is_control``
        ``True`` exactly for :data:`CONTROL_RESAMPLING_UNIT`.  A report from
        the control is not an uncertainty estimate; it is the measurement of
        what the wrong procedure would have said.
    ``wrong_because``
        Empty for a legitimate unit, and for the control the reason it is
        wrong, carried with the numbers so they cannot be quoted without it.
    ``replicates``, ``estimable_replicates``
        How many were drawn, and how many produced an **estimable** exponent --
        one refused for none of :data:`ESTIMABILITY_REASONS`.  Named for what
        it counts: an independent review accepted the rule and asked that the
        name stop saying "valid", since a replicate carrying a marginal
        goodness flag is estimable and is deliberately kept.

        Not "produced a finite number".  An independent review pointed a
        four-trial-per-cell design at the earlier implementation: the original
        fit was ``no_valid_exponent``, every one of fifty replicates was
        likewise refused, and the report said ``50 / 50`` valid and printed a
        percentile interval.  Finiteness is not estimability.

        Nor is it full validity, and that difference is deliberate:
        :data:`ESTIMABILITY_REASONS` explains at length why the two
        model-adequacy refusals are counted and reported but do not disqualify
        a replicate.
    ``refusals``
        ``((reason, count), ...)`` over the drawn replicates, so the *shape* of
        the failure is visible rather than only its size.
    ``curvature_interval``
        The percentile interval of the curvature coefficient over the valid
        replicates.  This is the cluster-robust curvature diagnostic: under
        common random numbers the marginal likelihood-ratio test's null
        distribution is not the chi-square it is compared against, and this
        interval is computed under exactly the pairing the design declares.
    ``exponent_mean``, ``exponent_stderr``
        Mean and standard deviation of the valid replicate exponents.
    ``interval``
        The percentile interval at :data:`CONFIDENCE_LEVEL`.
    """

    unit: str
    is_control: bool
    wrong_because: str
    replicates: int
    estimable_replicates: int
    refusals: tuple
    exponent_mean: float
    exponent_stderr: float
    interval: tuple
    curvature_interval: tuple
    seed: int = 0
    design_digest: str = ""

    def __post_init__(self):
        if not isinstance(self.design_digest, str):
            raise TypeError(f"design_digest must be a string, got "
                            f"{type(self.design_digest).__name__}")
        require_positive_int(self.seed, "seed", minimum=0)
        if self.unit not in RESAMPLING_UNITS:
            raise ValueError(f"unit must be one of {list(RESAMPLING_UNITS)}, "
                             f"got {self.unit!r}")
        if self.is_control != (self.unit == CONTROL_RESAMPLING_UNIT):
            raise ValueError(
                f"unit {self.unit!r} with is_control={self.is_control}; the "
                f"control flag is not a label a caller sets, it is what the "
                f"unit means")
        replicates = require_positive_int(self.replicates, "replicates",
                                          minimum=1)
        valid = require_positive_int(self.estimable_replicates,
                                     "estimable_replicates", minimum=0)
        if valid > replicates:
            raise ValueError(f"{valid} estimable of {replicates} drawn")
        refusals = tuple((str(reason), require_positive_int(count, "count",
                                                            minimum=0))
                         for reason, count in self.refusals)
        unknown = [r for r, _ in refusals if r not in NO_EXPONENT_REASONS]
        if unknown:
            raise ValueError(f"refusal reasons {unknown} are outside "
                             f"{list(NO_EXPONENT_REASONS)}")
        object.__setattr__(self, "replicates", replicates)
        object.__setattr__(self, "estimable_replicates", valid)
        object.__setattr__(self, "refusals", refusals)
        object.__setattr__(self, "interval", tuple(self.interval))
        object.__setattr__(self, "curvature_interval",
                           tuple(self.curvature_interval))

    @property
    def estimable_fraction(self) -> float:
        """Share of drawn replicates whose exponent was estimable."""
        return self.estimable_replicates / self.replicates

    @property
    def is_adequate(self) -> bool:
        """Is this record big enough and clean enough to be believed?

        Requires at least :data:`MIN_RESAMPLE_REPLICATES` draws, a valid
        fraction of at least :data:`MIN_VALID_RESAMPLE_FRACTION`, and a finite
        interval.  Both thresholds are frozen in this module, so an inadequate
        record cannot be made adequate by asking for fewer replicates.

        A control record is never adequate: it measures the size of a mistake,
        and a mistake is not an uncertainty estimate.
        """
        if self.is_control:
            return False
        return (self.replicates >= MIN_RESAMPLE_REPLICATES
                and self.estimable_fraction >= MIN_VALID_RESAMPLE_FRACTION
                and all(math.isfinite(value) for value in self.interval))


#: Why the control procedure is wrong, carried in every control report.
_CONTROL_REASON = (
    "resampled each coupling cell independently while the cells share master "
    "trial identifiers.  Common random numbers make those cells positively "
    "correlated, so this understates the uncertainty on the exponent.  It is "
    "a control that measures the size of that mistake, never an estimate.")


def resample_exponent(design, replicates, seed, unit=None) -> ResamplingReport:
    """Bootstrap the exponent, resampling the unit the design declares.

    ``master_trial``
        Draw master trial identifiers with replacement from the shared set and
        rebuild **every** cell from the same draw.  A path that is resampled
        twice is resampled twice at every coupling, which is what preserves the
        correlation common random numbers introduce.
    ``independent_namespace``
        Draw within each cell.  Correct here, and only here, because the design
        gate has already established that the cells' identifiers are disjoint.
    ``cell_independent_control``
        The same within-cell draw applied to a paired design.  Available on
        purpose, flagged on purpose, and never the default: without it there is
        no way to show that the paired procedure differs from the wrong one.

    ``seed`` drives an ordinary NumPy generator.  This is *analysis*
    randomness, over outcomes that are already recorded and hashed; it is not
    physical noise and shares nothing with the counter-keyed Brownian
    hierarchy, which no module in this layer can reach.
    """
    require_instance(design, ScalingDesign, "design")
    count = require_positive_int(replicates, "replicates", minimum=1)
    seed_value = require_positive_int(seed, "seed", minimum=0)
    if unit is None:
        unit = design.resampling_unit
    if isinstance(unit, bool) or not isinstance(unit, str):
        raise TypeError(f"unit must be a string, got {type(unit).__name__} "
                        f"{unit!r}")
    if unit not in RESAMPLING_UNITS:
        raise ValueError(f"unit must be one of {list(RESAMPLING_UNITS)}, got "
                         f"{unit!r}")
    if unit in DESIGN_RESAMPLING_UNITS and unit != design.resampling_unit:
        raise ValueError(
            f"this design declares {design.resampling_unit!r} and the call "
            f"asks for {unit!r}.  The two legitimate units are not "
            f"interchangeable — one resamples identifiers shared by every "
            f"cell, the other resamples inside cells the gate has proved "
            f"disjoint — so applying the wrong one to a design is a mistake "
            f"rather than a choice.  The named way to run the wrong procedure "
            f"deliberately is {CONTROL_RESAMPLING_UNIT!r}")

    rng = np.random.default_rng(seed_value)
    cells = design.cells
    outcomes = [np.array(cell.outcomes, dtype=bool) for cell in cells]
    shared = cells[0].trials

    estimates, curvatures, refusals = [], [], {}
    for index in range(count):
        if unit == "master_trial":
            draw = rng.integers(0, shared, size=shared)
            resampled = [values[draw] for values in outcomes]
        else:
            resampled = [values[rng.integers(0, values.size, size=values.size)]
                         for values in outcomes]
        replicate = _replicate_design(cells, resampled, unit, index)
        if replicate is None:
            refusals["ill_conditioned"] = refusals.get("ill_conditioned",
                                                       0) + 1
            continue
        report = fit_exponent(replicate)
        # Every refusal is counted by reason, never silently dropped.
        for reason in report.reasons:
            refusals[reason] = refusals.get(reason, 0) + 1
        if not _is_estimable(report):
            # No usable exponent at all: the likelihood was monotone, the
            # information singular, the optimizer lost, or the required
            # curvature diagnostic unidentifiable.  Folding such a replicate's
            # number into a percentile interval is what turned a design whose
            # every replicate was refused into a report claiming fifty valid
            # ones.  See ESTIMABILITY_REASONS for why a *goodness* refusal is
            # not treated the same way.
            continue
        estimates.append(report.exponent)
        curvatures.append(report.curvature)

    values = np.array(estimates, dtype=float)
    tail = 0.5 * (1.0 - CONFIDENCE_LEVEL) * 100.0
    if values.size:
        mean = float(np.mean(values))
        spread = float(np.std(values, ddof=1)) if values.size > 1 else 0.0
        interval = (float(np.percentile(values, tail)),
                    float(np.percentile(values, 100.0 - tail)))
        bends = np.array(curvatures, dtype=float)
        curvature_interval = (float(np.percentile(bends, tail)),
                              float(np.percentile(bends, 100.0 - tail)))
    else:
        mean = spread = float("nan")
        interval = (float("nan"), float("nan"))
        curvature_interval = (float("nan"), float("nan"))

    is_control = unit == CONTROL_RESAMPLING_UNIT
    return ResamplingReport(unit, is_control,
                            _CONTROL_REASON if is_control else "",
                            count, int(values.size),
                            tuple(sorted(refusals.items())), mean, spread,
                            interval, curvature_interval, seed_value,
                            design.digest)


def _is_estimable(report) -> bool:
    """Does this fit have a usable exponent, whatever it says about the model?

    ``True`` when the exponent is finite and the fit was refused for none of
    :data:`ESTIMABILITY_REASONS`.  A ``lack_of_fit`` or ``curvature`` refusal
    leaves the estimate perfectly well defined -- it says the power law is the
    wrong model, not that this number does not exist -- and under pairing those
    two tests are the miscalibrated ones, so they do not disqualify a bootstrap
    replicate.
    """
    return (math.isfinite(report.exponent)
            and not set(report.reasons) & set(ESTIMABILITY_REASONS))


def _replicate_design(cells, resampled, unit, index):
    """Rebuild one bootstrap replicate as a design, or ``None`` if it cannot.

    The replicate's identities are synthesized rather than copied, because a
    bootstrap draws the same master trial more than once and a cell may not
    hold a repeated identity.  They are synthesized *under the same rule the
    original design was checked by*: one shared tuple for a paired unit, and
    one distinct namespace per cell for independent namespaces.
    """
    rebuilt, shared_ids = [], None
    for position, (cell, values) in enumerate(zip(cells, resampled)):
        if unit == "independent_namespace":
            namespace = f"resample:{index}:{cell.namespace}:{position}"
            ids = tuple(MasterTrial(namespace, slot)
                        for slot in range(values.size))
        else:
            # One identity tuple, built once and shared by every cell: that is
            # what a paired design *is*, and rebuilding it per cell would make
            # the bootstrap's cost quadratic in the sweep for no guarantee.
            if shared_ids is None or len(shared_ids) != values.size:
                shared_ids = tuple(
                    MasterTrial(f"resample:{index}:{cells[0].namespace}", slot)
                    for slot in range(values.size))
            ids = shared_ids
        rebuilt.append(CouplingCell(cell.coupling, ids,
                                    tuple(bool(v) for v in values)))
    try:
        return ScalingDesign(tuple(rebuilt),
                             unit if unit in DESIGN_RESAMPLING_UNITS
                             else "master_trial")
    except ValueError:
        return None


# ---------------------------------------------------------------------------
# The authoritative, common-random-number-aware record
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class CalibratedExponent:
    """A marginal fit, its paired records, and the verdict all three imply.

    Why this type exists
    --------------------
    :func:`fit_exponent` is a correct *marginal* estimator: with ``y_i``
    committed out of ``n_i``, the binomial likelihood is right whatever the
    cells' joint distribution.  Its expected-information Wald interval, its
    deviance p-value and its curvature likelihood ratio are not: all three
    assume the coupling cells are independent observations, and a
    ``master_trial`` sweep violates that — trial 7 is the same Brownian path at
    every coupling.

    An independent review measured the interval gap on a 600-trial
    common-random-number sweep::

        expected-information Wald stderr      0.06112
        paired master-trial bootstrap         0.07932
        wrong independent-cell bootstrap      0.06219

    and the goodness gap over five hundred exact power-law datasets: the
    marginal deviance rejected 0.6 percent of them against a nominal five.

    What is authoritative here
    --------------------------
    For a ``master_trial`` design, **all three** decisions are made under the
    pairing:

    * the exponent, its standard error and its interval come from the paired
      bootstrap (:class:`ResamplingReport`);
    * curvature is decided by the paired bootstrap interval on the quadratic
      coefficient; and
    * general lack of fit is decided by the cluster-calibrated score test
      (:class:`ClusterScoreTest`).

    The marginal Wald interval, deviance p-value and likelihood ratio remain on
    :attr:`report` and are **display only**.

    For an ``independent_namespace`` design the cells genuinely are independent,
    so the marginal tests are kept as authoritative — and the cluster test
    reduces to the same independent score test anyway, so it is run and
    reported for both.

    Construction
    ------------
    Factory-only, via :func:`calibrate_exponent`, and every published field is
    recomputed from the records it carries and required to agree.  A review
    built one whose verdict was ``valid_exponent`` with no reasons over records
    that said ``curvature_unidentified, uninformative, inadequate_resampling``,
    and handed it to the physical gate; that object is now unconstructible.
    """

    design_unit: str
    report: ExponentReport
    resampling: ResamplingReport
    lack_of_fit: ClusterScoreTest
    curvature_test: ClusterScoreTest
    verdict: str
    reasons: tuple
    exponent: float
    stderr: float
    interval: tuple
    curvature_interval: tuple
    design_digest: str
    seal: object = None

    def __post_init__(self):
        _require_seal(self.seal, "CalibratedExponent")
        require_instance(self.report, ExponentReport, "report")
        require_instance(self.resampling, ResamplingReport, "resampling")
        require_instance(self.lack_of_fit, ClusterScoreTest, "lack_of_fit")
        require_instance(self.curvature_test, ClusterScoreTest,
                         "curvature_test")
        for name, expected in (("lack_of_fit", "lack_of_fit"),
                               ("curvature_test", "curvature")):
            if getattr(self, name).kind != expected:
                raise ValueError(f"{name} carries a {getattr(self, name).kind!r} "
                                 f"test, not a {expected!r} one")
        if self.design_unit not in DESIGN_RESAMPLING_UNITS:
            raise ValueError(f"design_unit must be one of "
                             f"{list(DESIGN_RESAMPLING_UNITS)}, got "
                             f"{self.design_unit!r}")
        if self.resampling.is_control:
            raise ValueError(
                "a calibrated exponent cannot be built on the "
                "cell_independent_control record; that record measures the "
                "size of a known mistake and is never an uncertainty estimate")
        # Every record must be about the same observations.
        digests = {self.design_digest, self.report.design_digest,
                   self.resampling.design_digest,
                   self.lack_of_fit.design_digest,
                   self.curvature_test.design_digest}
        if len(digests) != 1:
            raise ValueError(
                f"this calibration carries records computed from "
                f"{len(digests)} different designs; the fit, the resampling "
                f"record and the lack-of-fit test are bound to one digest over "
                f"the ordered couplings, the master-trial identities and every "
                f"outcome, and a mismatch means they do not describe one "
                f"experiment")
        for name in ("lack_of_fit", "curvature_test"):
            if getattr(self, name).unit != self.design_unit:
                raise ValueError(f"the {name} test was run under a different "
                                 f"resampling unit")
        # Every *published* value is recomputed, so it cannot be asserted.
        verdict, reasons = _calibrated_verdict(self.design_unit, self.report,
                                               self.resampling,
                                               self.lack_of_fit,
                                               self.curvature_test)
        if self.verdict != verdict or tuple(self.reasons) != reasons:
            raise ValueError(
                f"this calibration publishes {self.verdict!r} "
                f"{list(self.reasons)}; its own records imply {verdict!r} "
                f"{list(reasons)}.  A verdict that disagrees with the evidence "
                f"it carries is the forgery this type exists to refuse")
        for name, value in (("exponent", self.resampling.exponent_mean),
                            ("stderr", self.resampling.exponent_stderr)):
            published = getattr(self, name)
            if not (published == value
                    or (math.isnan(published) and math.isnan(value))):
                raise ValueError(f"the published {name} is {published!r} and "
                                 f"the paired record's is {value!r}")
        for name, value in (("interval", self.resampling.interval),
                            ("curvature_interval",
                             self.resampling.curvature_interval)):
            published = tuple(getattr(self, name))
            if not all(a == b or (math.isnan(a) and math.isnan(b))
                       for a, b in zip(published, value)) or \
                    len(published) != len(value):
                raise ValueError(f"the published {name} is {published} and the "
                                 f"paired record's is {value}")
        object.__setattr__(self, "reasons", tuple(self.reasons))
        object.__setattr__(self, "interval", tuple(self.interval))
        object.__setattr__(self, "curvature_interval",
                           tuple(self.curvature_interval))
        object.__setattr__(self, "seal", None)

    @property
    def is_valid(self) -> bool:
        """``True`` exactly when the paired evidence supports an exponent."""
        return self.verdict == "valid_exponent"

    def interval_contains(self, value) -> bool:
        """Is ``value`` inside the **paired** exponent interval?

        The question ``compare.causal_decision`` asks of a control: is this
        still compatible with quadratic scaling?  A non-finite interval
        contains everything, which is the honest reading of a record that could
        not be estimated — but a control is now blocked by its *verdict* before
        this is ever consulted, so an inadequate control cannot be rescued by a
        narrow provisional interval either.
        """
        target = require_finite(value, "value")
        low, high = self.interval
        if not (math.isfinite(low) and math.isfinite(high)):
            return True
        return low <= target <= high


def _calibrated_verdict(unit, report, resampling, lack_of_fit,
                        curvature_test):
    """Return ``(verdict, reasons)`` from a fit and its paired records.

    One function, called by both :func:`calibrate_exponent` and
    :meth:`CalibratedExponent.__post_init__`, so that what a calibration
    *publishes* and what its records *imply* cannot be two different things.

    Which marginal reasons survive depends on the unit, and that is the whole
    point of the type.  Under ``master_trial`` the two model-adequacy verdicts
    are computed from statistics whose null distribution the pairing
    invalidates, so both are dropped and replaced: ``lack_of_fit`` by the
    cluster-calibrated score test, ``curvature`` by the paired bootstrap
    interval on the quadratic coefficient.  Under ``independent_namespace`` the
    marginal tests are correct and are kept.
    """
    reasons = []
    paired = unit == "master_trial"
    for reason in report.reasons:
        if paired and reason in GOODNESS_REASONS:
            continue                       # replaced below, under the pairing
        reasons.append(reason)

    if not resampling.is_adequate:
        reasons.append("inadequate_resampling")

    low, high = resampling.interval
    half = (0.5 * (high - low) if (math.isfinite(low) and math.isfinite(high))
            else float("inf"))
    if not (math.isfinite(half) and half <= MAX_EXPONENT_HALF_WIDTH):
        reasons.append("uninformative")

    if paired:
        # Curvature is decided by its own one-dimensional cluster score test,
        # not by the percentile interval of bootstrap point estimates.  An
        # independent review rerolled that interval across analysis seeds and
        # got "excludes zero" 58 times in 100 on one fixed dataset; it also
        # excluded zero on 6.7 percent of interior exact nulls.  A design whose
        # curvature calibration is inadequate makes **no** positive curvature
        # claim: it reports uncalibrated_curvature and the interval stays
        # descriptive.
        if not curvature_test.adequate:
            reasons.append("uncalibrated_curvature")
        elif curvature_test.rejects:
            reasons.append("curvature")
        if not lack_of_fit.adequate:
            reasons.append("uncalibrated_fit")
        elif lack_of_fit.rejects:
            reasons.append("lack_of_fit")

    ordered = tuple(r for r in NO_EXPONENT_REASONS if r in set(reasons))
    return ("no_valid_exponent" if ordered else "valid_exponent"), ordered


def _authoritative_seed(digest: str) -> int:
    """The one seed an authoritative resampling record may use.

    Derived from the design's own digest, so a decision path has no seed
    argument to reroll and two runs over the same observations give the same
    interval.  ``compare.py``'s factories take neither a seed nor a replicate
    count for exactly this reason.
    """
    return int(digest[15:30], 16)


def calibrate_exponent(design) -> CalibratedExponent:
    """Fit ``design`` and calibrate every decision under the pairing it declares.

    The door every decision in ``compare.py`` goes through.  It runs the
    marginal fit, the paired resampling under the design's own declared unit,
    and the cluster-calibrated lack-of-fit test, and returns a
    :class:`CalibratedExponent` bound to the design's digest.

    **It takes no replicate count and no seed.**  Both are frozen:
    :data:`RESAMPLE_REPLICATES` draws, seeded from the design's own digest.  An
    independent review pointed out that a decision path accepting an analysis
    seed is a decision that can be rerolled, and measured the cost — one fixed
    dataset gave "curvature excludes zero" for 58 of 100 seeds.  Two
    calibrations of the same observations now give the same answer, always.

    An inadequate resampling record — fewer than
    :data:`MIN_RESAMPLE_REPLICATES` draws, or an estimable fraction below
    :data:`MIN_VALID_RESAMPLE_FRACTION` — is ``no_valid_exponent`` under
    ``inadequate_resampling``.  An inadequate lack-of-fit calibration is
    ``uncalibrated_fit`` and an inadequate curvature calibration is
    ``uncalibrated_curvature``.  None of the three is a narrower interval
    computed from whatever survived.
    """
    require_instance(design, ScalingDesign, "design")
    digest = design.digest
    report = fit_exponent(design)
    resampling = resample_exponent(design, RESAMPLE_REPLICATES,
                                   _authoritative_seed(digest),
                                   unit=design.resampling_unit)
    lack_of_fit = cluster_lack_of_fit(design, report)
    curvature_test = cluster_curvature(design, report)
    verdict, reasons = _calibrated_verdict(design.resampling_unit, report,
                                           resampling, lack_of_fit,
                                           curvature_test)
    return CalibratedExponent(
        design_unit=design.resampling_unit, report=report,
        resampling=resampling, lack_of_fit=lack_of_fit,
        curvature_test=curvature_test, verdict=verdict,
        reasons=reasons, exponent=resampling.exponent_mean,
        stderr=resampling.exponent_stderr, interval=resampling.interval,
        curvature_interval=resampling.curvature_interval,
        design_digest=digest, seal=_SEAL)
