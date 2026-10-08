#!/usr/bin/env python3
"""
raw_config.py — everything a raw event run is allowed to be told, and nothing
else.

Why this module exists
----------------------
``verify.py`` parses the raw modules and fails if any of them can reach the
analytic prediction.  That closes the route through the *code*.  It does not
close the route through the *caller*: a runner whose configuration object
accepts arbitrary attributes can be handed the prediction by hand, and the
resulting "measured" commitment law is then a restatement of the input.

So a raw run accepts numbers and labels, and refuses everything else:

* **No callables.**  A callable is an arbitrary computation whose source this
  package cannot parse.  ``predictor=lambda K: ...`` is precisely how a
  prescribed commitment law enters an event generator that is supposed to
  produce one.
* **No opaque objects.**  An object with attributes is a callable with extra
  steps.  The one exception is :class:`RawClockGrid`, which is defined here and
  holds two validated numbers.
* **No arrays.**  An amplitude squared, an intensity weight, a per-clock
  relaxation rate and a prescribed hazard all arrive as arrays of numbers, and
  all of them are the prediction spelled differently.  The clock population is
  named by *its rule* — a support half-width and a count — not by a vector of
  values that could have been computed from anything.
* **No unrecognised keys**, so a forbidden input cannot be smuggled in under a
  name nobody thought to ban.  :data:`FORBIDDEN_RAW_INPUTS` exists on top of
  that, to make the common mistakes fail by name rather than as a generic
  "unexpected keyword".

What is NOT here
----------------
No noise, no stepping, no dwell, no commitment, no ledger, no statistics.  This
is the boundary those stages will be handed; the stages are later tickets, and
nothing in this module produces an event or a number that could be mistaken for
one.

The duplication with ``analytic.FlatGrid``
------------------------------------------
:class:`RawClockGrid` re-derives the symmetric midpoint grid rather than
importing ``analytic.FlatGrid``, for the same reason ``analytic.py``
re-derives ``K > |Delta|`` rather than importing ``model.is_eligible``: a raw
module that imports the analytic module is exactly what this ticket exists to
prevent.  Deliberate duplication is safe only when something compares the two
copies, so ``verify.py`` — which is allowed to see both sides — asserts that
the two grids agree exactly, node for node, with a tolerance of zero.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .validation import (require_exact_type, require_finite, require_instance,
                         require_nonnegative, require_positive_int)

__all__ = [
    "FORBIDDEN_RAW_INPUTS",
    "RAW_INPUT_NAMES",
    "RUN_LABELS",
    "RawClockGrid",
    "RawEventConfig",
    "raw_event_config",
]

#: Input names that are refused by name, with an explanation, rather than as a
#: generic unexpected keyword.  Every one of them is the analytic prediction,
#: some transformation of it, or a callback that could compute it.  The point
#: is the error message: a caller who reaches for one of these is making a
#: category mistake about what a raw run is, and a bare ``TypeError`` naming an
#: unexpected keyword would not tell them so.
FORBIDDEN_RAW_INPUTS = frozenset({
    # Callbacks that could compute anything at all.
    "callback", "predict", "predictor", "predictor_callback", "rate_callback",
    "commitment_rule", "propensity",
    # Prescribed commitment laws.
    "commitment_hazard", "hazard", "hazard_rate", "prescribed_hazard",
    "survival_curve", "expected_probability", "expected_commitments",
    # The prediction, and the ways it is usually written.
    "analytic", "analytic_curve", "reference_curve", "comparator",
    "comparison", "oracle", "killed_diffusion",
    "local_relaxation_rate", "relaxation_rate", "rate_sum", "rate_weights",
    "tongue_rate_flux", "tongue_rate_sum", "flux",
    "efolds", "efolds_available", "semicircle", "semicircle_constant",
    # Amplitude-squared inputs.  The dynamics receive a coupling, never a
    # coupling squared; the same rule applies to their configuration.
    "amplitude2", "amplitude_squared", "intensity", "intensity_weight",
    "peak_intensity", "born_weight", "probability_weight", "weights",
    # Answers to the question being asked.
    "expected_exponent", "target_exponent",
})

#: The complete set of names a raw configuration mapping may contain.  Anything
#: outside it is rejected, so this tuple *is* the input contract.
RAW_INPUT_NAMES = (
    "grid",
    "peak_coupling",
    "pulse_duration",
    "phase_diffusion",
    "lock_tolerance",
    "dwell_time",
    "timestep",
    "trials",
    "stream_namespace",
    "run_label",
    "pulse_centre",
    "shadow_trials",
)

#: A run is a pilot or it is production, and the two may never share a stream
#: namespace or an analysis.  Carrying the label in the configuration means the
#: distinction survives into every manifest written from it.
RUN_LABELS = ("pilot", "production")

#: Characters a stream namespace may contain.  A namespace ends up inside an
#: immutable counter key, so it has to be a plain identifier-like token rather
#: than an arbitrary string that might not round-trip.
_NAMESPACE_CHARS = frozenset(
    "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-")

#: Longest permitted stream namespace.
_MAX_NAMESPACE = 64

#: Scalar types a raw input value may have.  ``bool`` is a subclass of ``int``
#: and so passes this test; the per-field validators reject it, which is where
#: that mistake belongs.
_SCALAR_TYPES = (int, float, str, np.integer, np.floating)


def _require_label(value, name: str, allowed) -> str:
    """Return ``value`` as one of ``allowed``, or raise naming ``name``.

    ``TypeError`` if it is not a string at all, ``ValueError`` if it is a
    string that names nothing.  Both are refusals to guess: a run labelled
    ``"prod"`` is not obviously production, and treating it as production is
    how a pilot result reaches a production manifest.
    """
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, str):
        raise TypeError(f"{name} must be a string, got "
                        f"{type(value).__name__} {value!r}")
    if value not in allowed:
        raise ValueError(f"{name} must be one of {tuple(allowed)!r}, got "
                         f"{value!r}")
    return value


def _require_namespace(value, name: str) -> str:
    """Return ``value`` as a usable stream-namespace token, or raise."""
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
                         f"and hyphen; got {value!r} containing {bad!r}")
    return value


def _require_positive(value, name: str) -> float:
    """Return ``value`` as a float strictly greater than zero, or raise."""
    out = require_finite(value, name)
    if not out > 0.0:
        raise ValueError(f"{name} must be strictly positive, got {out!r}")
    return out


@dataclass(frozen=True)
class RawClockGrid:
    """The finite clock population, named by its rule rather than its values.

    ``half_width`` is the spectral support and ``n_points`` the number of
    clocks; the detunings are ``n_points`` equally spaced midpoints of
    ``[-half_width, +half_width]``, spacing ``h = 2*half_width/n_points``.

    Midpoints, not endpoints, and for a reason that outlives the arithmetic: an
    endpoint grid places a clock exactly at ``Delta = 0`` and exactly at the
    support edges, which makes the result depend on whether a tongue edge lands
    on a node.  The grid is held fixed while the coupling is swept, so that the
    coupling dependence being measured belongs to the tongue and not to the
    sampler.

    Every clock on this grid is a real competitor, not a quadrature node.
    Adding clocks adds commitment opportunities, so ``n_points`` is a physical
    parameter and a change to it is a *population sensitivity*, never a
    convergence study.  The two numbers are stored, and the detunings derived,
    precisely so that no caller can hand a raw run a vector of detunings that
    was computed from something else.

    That distinction is why a single clock is allowed here and is refused by
    ``analytic.FlatGrid``, which otherwise computes the identical grid.  One
    point is not a quadrature, so the analytic side rejects it; one clock is a
    population, and specifically it is the mandatory central-clock negative
    control — the run at exactly zero detuning, in which tongue widening is
    absent and which must *not* look quadratic if the width-times-rate
    explanation is to survive.  ``RawClockGrid(W, 1)`` puts that one clock at
    detuning zero exactly.  ``verify.py`` asserts both the agreement and this
    one deliberate divergence.
    """

    half_width: float
    n_points: int

    def __post_init__(self):
        half_width = _require_positive(self.half_width, "half_width")
        n_points = require_positive_int(self.n_points, "n_points", minimum=1)
        # Store the normalized values, so a frozen grid can never hold
        # something that merely converted once inside the constructor.
        object.__setattr__(self, "half_width", half_width)
        object.__setattr__(self, "n_points", n_points)

    @property
    def spacing(self) -> float:
        """Distance between adjacent clocks, ``2*half_width/n_points``."""
        return 2.0 * self.half_width / self.n_points

    @property
    def detunings(self) -> np.ndarray:
        """The clock detunings, low to high."""
        h = self.spacing
        return -self.half_width + (np.arange(self.n_points) + 0.5) * h

    @property
    def parity(self) -> str:
        """``"even"`` or ``"odd"`` — load-bearing, so it has a name.

        With an even count the grid straddles zero and has no central clock;
        with an odd count it has one exactly at ``Delta = 0``.  At weak
        coupling that difference is the difference between no eligible
        competitor and one, which moves the observed onset.  It is reported
        rather than left to be inferred from ``n_points``.
        """
        return "even" if self.n_points % 2 == 0 else "odd"


@dataclass(frozen=True)
class RawEventConfig:
    """The complete input to a raw event run: physical and numerical, no more.

    Every field is a number, a label, or the clock grid.  There is nowhere in
    here to put a predicted rate, a commitment hazard, or a callback, and that
    is the only reason this type exists rather than a plain dictionary.

    The fields, and what each one fixes:

    ``grid``
        The finite clock population.  See :class:`RawClockGrid`.
    ``peak_coupling``
        Peak of the photon-pulse envelope.  Non-negative rather than positive:
        zero coupling is the deterministic control that must produce no
        commitment at all, and refusing it would remove a required negative
        control from the sweep.
    ``pulse_duration``, ``pulse_centre``
        Support and centre of the raised-cosine envelope, in the one reference
        time system used for the whole sweep.
    ``phase_diffusion``
        Strength of the independent white phase noise, the same for every
        clock, every detuning and every coupling.  Zero is permitted and is the
        deterministic limit: Euler-Maruyama must reduce to Euler exactly.
    ``lock_tolerance``
        Fixed angular half-width of the lock band, in radians, **strictly less
        than half a turn** — the same domain ``commitment.LockCriterion``
        enforces, so that a configuration this boundary accepts can always
        build its own state machine.  Fixed for the whole sweep; a tolerance
        that scaled with coupling would make the selection criterion
        amplitude-dependent, which is the thing being tested.
    ``dwell_time``
        The fixed physical commitment dwell.  Strictly positive: a zero dwell
        makes "remained locked for an uninterrupted interval" vacuous, so it is
        rejected rather than allowed to mean "committed on sight".
    ``timestep``
        The Euler-Maruyama step.  No relation to ``dwell_time`` is imposed —
        dwell shorter than, equal to, and non-integral relative to the step are
        all required to obey the same timestamp rule.
    ``trials``
        Number of independent trials.
    ``stream_namespace``
        Root of the counter-keyed random-number hierarchy.  Pilot and
        production must not share one.
    ``run_label``
        ``"pilot"`` or ``"production"``.  A pilot may choose a coupling range
        under frozen count rules and may never enter a production estimate, so
        which one a run is has to be recorded at the input, not inferred later
        from where its output was filed.
    ``shadow_trials``
        How many of this run's trials are *also* walked with the first-winner
        stop switched off, as the plan's predeclared no-stop shadow sample.
        At least one, and never more than ``trials``.

        It is a configuration field rather than a runtime argument because
        "predeclared" is the whole content of the requirement.  Per-clock
        statistics in a first-passage race are informatively censored — once
        one clock stops the channel, every other clock's later opportunities
        are unobserved — and the shadow sample is what separates a mismatch
        that never commits from a mismatch that never got the chance.  Choosing
        its size after seeing which mismatches looked quiet would make it
        useless for exactly that purpose, so the size is frozen at the input
        boundary and written into the manifest before the run starts.

        The shadow trials are the first ``shadow_trials`` trial identifiers of
        the run.  A *rule* rather than a list, for the same reason the clock
        grid is a rule rather than a vector of detunings: a list could have
        been chosen from the answer.

        One is the minimum rather than zero because a run with no shadow sample
        would report per-clock counts with no way to tell censoring from
        absence, and would look exactly like a run that had checked.

    Two cross-field rules.  The spectral support must be strictly wider than
    the peak coupling: if it is not, the tongue reaches past the edge of the
    population, every clock at the edge is eligible at once, and the eligible
    count stops being a property of the tongue.  And the shadow sample cannot
    be larger than the run it is a sample of.
    """

    grid: RawClockGrid
    peak_coupling: float
    pulse_duration: float
    phase_diffusion: float
    lock_tolerance: float
    dwell_time: float
    timestep: float
    trials: int
    stream_namespace: str
    run_label: str
    pulse_centre: float = 0.0
    shadow_trials: int = 1

    def __post_init__(self):
        # Exactly a RawClockGrid, not merely an instance of one.  A grid
        # subclass could override `detunings` to return anything at all —
        # including a vector computed from the prediction — while satisfying
        # every isinstance check between here and the event generator.
        grid = require_exact_type(self.grid, RawClockGrid, "grid")
        peak_coupling = require_nonnegative(self.peak_coupling, "peak_coupling")
        pulse_duration = _require_positive(self.pulse_duration,
                                           "pulse_duration")
        phase_diffusion = require_nonnegative(self.phase_diffusion,
                                              "phase_diffusion")
        lock_tolerance = _require_positive(self.lock_tolerance,
                                           "lock_tolerance")
        dwell_time = _require_positive(self.dwell_time, "dwell_time")
        timestep = _require_positive(self.timestep, "timestep")
        trials = require_positive_int(self.trials, "trials", minimum=1)
        stream_namespace = _require_namespace(self.stream_namespace,
                                              "stream_namespace")
        run_label = _require_label(self.run_label, "run_label", RUN_LABELS)
        pulse_centre = require_finite(self.pulse_centre, "pulse_centre")
        shadow_trials = require_positive_int(self.shadow_trials,
                                             "shadow_trials", minimum=1)

        # A band of half a turn or more is the whole circle read through a
        # shortest-distance comparison, so "inside the band" would be true
        # everywhere and dwell would measure nothing.
        #
        # **Strictly** less than pi, and the strictness is not a nicety.  This
        # bound used to be `> np.pi`, which admitted `lock_tolerance == pi`;
        # `validate_raw_config` and `raw_noise_stream` then accepted the
        # configuration and `raw_lock_criterion` raised, because
        # `commitment.LockCriterion` has always required a strict inequality.
        # The declared single input boundary of the raw process was therefore
        # accepting a run configuration that could not instantiate its own
        # state machine.  One domain, stated once, checked at every factory.
        if not lock_tolerance < np.pi:
            raise ValueError(f"lock_tolerance must be strictly less than pi "
                             f"radians (a circular distance never exceeds pi, "
                             f"so a band of pi or more admits every phase and "
                             f"is not a band), got {lock_tolerance!r}")
        if not grid.half_width > peak_coupling:
            raise ValueError(
                f"the spectral support must be strictly wider than the peak "
                f"coupling, or the tongue reaches past the edge of the clock "
                f"population and the eligible count stops measuring the "
                f"tongue; got half_width={grid.half_width!r} and "
                f"peak_coupling={peak_coupling!r}")
        if shadow_trials > trials:
            raise ValueError(
                f"shadow_trials={shadow_trials} exceeds trials={trials}; the "
                f"no-stop shadow sample is a subset of the run's own trials — "
                f"it re-walks them with the first-winner stop switched off — "
                f"so it cannot contain trials the run never performs")

        for name, value in (("grid", grid),
                            ("peak_coupling", peak_coupling),
                            ("pulse_duration", pulse_duration),
                            ("phase_diffusion", phase_diffusion),
                            ("lock_tolerance", lock_tolerance),
                            ("dwell_time", dwell_time),
                            ("timestep", timestep),
                            ("trials", trials),
                            ("stream_namespace", stream_namespace),
                            ("run_label", run_label),
                            ("pulse_centre", pulse_centre),
                            ("shadow_trials", shadow_trials)):
            object.__setattr__(self, name, value)


def raw_event_config(inputs) -> RawEventConfig:
    """Build a :class:`RawEventConfig` from a mapping, refusing by name.

    ``RawEventConfig(**inputs)`` already rejects an unknown keyword, but it
    rejects it as an unexpected keyword argument, which reads like a typo.  It
    is not a typo.  Passing ``hazard=`` or ``predictor=`` to an event generator
    is a category mistake about what the run is for, and this function is here
    so that the mistake fails saying so.

    Four refusals, in order, so the most informative one wins:

    1. a key in :data:`FORBIDDEN_RAW_INPUTS`, which names the prediction or a
       callback that could compute it;
    2. a key outside :data:`RAW_INPUT_NAMES`;
    3. a callable value, whatever it is called; and
    4. a value that is neither a number, a string, nor the clock grid — an
       array of weights and an object with a ``__call__`` both land here.

    ``grid`` may be given as a nested mapping of ``half_width`` and
    ``n_points``, which is what a JSON manifest will contain.

    Every refusal is a ``TypeError`` or a ``ValueError``, matching the rest of
    the package: ``TypeError`` for "not that kind of thing", ``ValueError`` for
    "right kind, impossible value".
    """
    require_instance(inputs, dict, "inputs")

    for key in inputs:
        if not isinstance(key, str):
            raise TypeError(f"raw configuration keys must be strings, got "
                            f"{type(key).__name__} {key!r}")

    forbidden = sorted(set(inputs) & FORBIDDEN_RAW_INPUTS)
    if forbidden:
        raise TypeError(
            f"a raw event run may not be given {forbidden!r}: these name the "
            f"analytic prediction, a transformation of it, or a callback that "
            f"could compute one.  The raw process exists to produce events "
            f"without seeing the prediction; a separate comparison process "
            f"reads the closed ledger afterwards")

    unknown = sorted(set(inputs) - set(RAW_INPUT_NAMES))
    if unknown:
        raise TypeError(f"unknown raw configuration input(s) {unknown!r}; the "
                        f"accepted inputs are {sorted(RAW_INPUT_NAMES)!r}")

    prepared = dict(inputs)
    grid = prepared.get("grid")
    if isinstance(grid, dict):
        grid_unknown = sorted(set(grid) - {"half_width", "n_points"})
        if grid_unknown:
            raise TypeError(f"unknown clock-grid input(s) {grid_unknown!r}; a "
                            f"grid is named by half_width and n_points only")
        prepared["grid"] = RawClockGrid(**grid)

    for key, value in prepared.items():
        if callable(value):
            raise TypeError(
                f"{key} must be a number, not a callable: a raw event run "
                f"takes physical and numerical inputs, and a callback is an "
                f"arbitrary computation this package cannot inspect")
        if key != "grid" and not isinstance(value, _SCALAR_TYPES):
            raise TypeError(
                f"{key} must be a single number or label, got "
                f"{type(value).__name__}; arrays and objects are refused "
                f"because an intensity, an amplitude squared, a rate weight "
                f"and a prescribed hazard all arrive as one")

    missing = sorted(set(RAW_INPUT_NAMES) - set(prepared)
                     - {"pulse_centre", "shadow_trials"})
    if missing:
        raise TypeError(f"missing raw configuration input(s) {missing!r}")

    return RawEventConfig(**prepared)
