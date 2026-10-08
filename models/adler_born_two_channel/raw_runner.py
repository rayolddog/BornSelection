#!/usr/bin/env python3
"""
raw_runner.py — the entry point of the raw event process, and the root of the
import graph that must never reach the analytic prediction.

What "raw" means here
---------------------
The experiment asks whether noisy Adler clocks commit in a pattern matching an
independently computed prediction.  A "raw" run is the half of that experiment
which must not know the answer: it evolves phases, applies the dwell rule,
records what happened, closes its ledger, and stops.  A *separate* comparison
process then opens the closed ledger and — only there — imports ``analytic``.

This module is the root of the first half.  ``verify.py`` walks its imports
transitively and fails if any module reachable from here can name the
prediction, the killed-diffusion oracle, the moving-band audit, the analysis
layer, the comparison layer, or the reporting layer.  It also runs a subprocess
that imports this module and asserts that none of those are in ``sys.modules``
afterwards, because a structural scan cannot see an import that happens at run
time.

What is implemented, and what is not
------------------------------------
**This boundary now runs a population.**  ``dynamics`` steps a phase with
Euler-Maruyama over the exact elementary mesh, under either a pulsed or a
stationary drive, with exact tongue entry, exit and re-entry handoffs and the
fixed-contraction width-only control; ``commitment`` runs the
contraction-guarded lock predicate and the timestamp dwell clock over the
endpoints; ``raw_experiments`` turns the two into an in-memory one-clock
record; ``raw_race`` puts ``N`` of them on one immutable grid and stops the
channel at the first dwell completion; and this module freezes that run into a
manifest, three CSV tables and a close marker on disk.
:data:`IMPLEMENTED_RAW_STAGES` is the list, and this paragraph is checked
against it.

**A first clock is a model stopping event.**  It is not a detector click, an
absorption, a measurement outcome, a unique actuality, a microscopic bath, one
side of a two-channel competition, or a Born-rule result.  The trial ends
because this one-channel first-passage study defines it to end there, and no
physical exclusivity is claimed.  ``raw_race.RACE_NON_CLAIMS`` carries the
full list and :data:`RAW_NON_CLAIMS` re-exports it.

**What is written here is diagnostic.**  Nothing produced through this module
is production-cleared.  The stationary killed-diffusion gate passes at its
reference tolerances; the moving-band audit returns a frozen
``numerical_no_result``, pooled and in fourteen of fifteen per-regime cells;
and the production numerical budget ticket 07 froze is **not met** — its
disposition is a numerical no-result and no matrix was feasible among the
options it searched, so nothing written here has passed a production gate and
nothing can until that block is lifted.  Every manifest written here says so in
its own fields — see
:data:`NUMERICAL_GATE` — so a ledger cannot be read later as though it had
passed a gate that was never applied to it.  Everything still absent is in
:data:`UNIMPLEMENTED_RAW_STAGES`, and it is absent rather than stubbed.

:func:`validate_raw_config` remains the boundary function — the single place a
raw run accepts its inputs.  :func:`raw_noise_stream`,
:func:`raw_lock_criterion`, :func:`raw_one_clock_path`,
:func:`raw_run_window`, :func:`raw_clock_paths` and
:func:`raw_clock_population` are the only routes from a validated
configuration to the objects a run needs, and none of them has an
argument through which a predicted rate, a weight, a prescribed commitment law,
or a per-run contraction rate could arrive.

``stochastic.py`` still supplies numbers rather than behaviour: it can say what
phase kick a given clock receives on a given finest step of a given trial, and
what phase that clock was prepared in, and has no idea what a phase is.

Why the entry point is a module rather than a script
-----------------------------------------------------
The raw process is launched separately from anything that imports the
prediction, which is what makes the isolation real at run time rather than
merely at import time.  :func:`write_raw_run` is the launcher: a caller hands
it a validated configuration and a **run name**, and it returns having closed
the run under the package's own ignored results directory.  It takes a name
rather than a directory because generated output belongs in exactly one place;
an earlier revision accepted an arbitrary directory and an independent review
correctly called that a contradiction of the contract the scoped path helper
beside it was declaring.
"""
from __future__ import annotations

import hashlib
import math
import os
import platform
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from . import raw_experiments
from .commitment import LockCriterion
from .dynamics import (DRIVE_LABELS, FIXED_CONTRACTION_RATE, MODEL_LABELS,
                       ClockPath, PulseTrain, StationaryDrive)
from .model import RaisedCosinePulse
from .raw_config import RawClockGrid, RawEventConfig
from .raw_ledger import (CLOCK_COLUMNS, DIGEST_ALGORITHM, LEDGER_COLUMNS,
                         LEDGER_SCHEMA_VERSION, MANIFEST_SCHEMA_VERSION,
                         SHADOW_COLUMNS, TABLE_NAMES, TRIAL_CATEGORIES,
                         CloseMarker, RawManifest,
                         canonical_manifest_bytes, canonical_marker_bytes,
                         digest_bytes, encode_row, header_row, parse_marker,
                         require_closed_ledger)
from .raw_race import (DEFAULT_CLOCK_BLOCK, DEFAULT_STEP_WINDOW,
                       RACE_NON_CLAIMS, ClockPopulation, PopulationIdentity,
                       no_stop_exposure, race_one_channel)
from .stochastic import KEY_SCHEMA, FinestMesh, PhaseNoiseStream
from .validation import (require_exact_type, require_finite,
                         require_positive_int)

__all__ = [
    "EXPERIMENT_PILOT_NAMESPACE_PREFIX",
    "IMPLEMENTED_RAW_STAGES",
    "MOVING_BAND_STATUS",
    "NUMERICAL_GATE",
    "RAW_NON_CLAIMS",
    "RESULTS_DIRNAME",
    "RUN_FILENAMES",
    "SHADOW_RULE",
    "STATIONARY_ORACLE_STATUS",
    "UNIMPLEMENTED_RAW_STAGES",
    "VALIDATION_NOTE",
    "VALIDATION_NOTE_REVISION",
    "RawRunReport",
    "open_raw_run",
    "raw_clock_paths",
    "raw_clock_population",
    "raw_lock_criterion",
    "raw_manifest",
    "raw_noise_stream",
    "raw_one_clock_path",
    "raw_results_directory",
    "raw_run_window",
    "validate_raw_config",
    "write_raw_run",
]

#: What the raw side of this package actually does today.
IMPLEMENTED_RAW_STAGES = (
    "validated raw configuration boundary (raw_config.py)",
    "closed-ledger manifest, table schema, marker and digest protocol "
    "(raw_ledger.py)",
    "counter-keyed streamed phase kicks, the nested physical Brownian tree, "
    "and the keyed prepared initial ensemble (stochastic.py)",
    "Euler-Maruyama one-clock phase evolution over the exact elementary mesh, "
    "under a pulsed or a stationary drive, with exact tongue "
    "entry/exit/re-entry handoffs (dynamics.py)",
    "the fixed-contraction width-only control, sharing one Brownian tree with "
    "the full Adler model (dynamics.py)",
    "the contraction-guarded lock predicate and the timestamp dwell clock "
    "(commitment.py)",
    "the S1 deterministic-limit, S2 synthetic dwell and S4 one-clock "
    "experiments, as in-memory records (raw_experiments.py)",
    "the finite-population one-channel race over a declared population and "
    "grid identity, first-completion stopping, retained same-endpoint "
    "co-completers, censored per-clock exposure and the predeclared no-stop "
    "shadow sample (raw_race.py)",
    "the immutable hashed run on disk — manifest, trial ledger, one "
    "authoritative row per (trial, clock), the no-stop shadow sample and the "
    "close marker — and the gate a consumer must pass before reading one, "
    "which derives every channel-level number from those rows rather than "
    "believing it (raw_runner.py, raw_ledger.py)",
    "transitive structural and runtime isolation of this import graph",
)

#: What a raw run needs and does not have.  This tuple is the honest answer to
#: "can I interpret the experiment yet", and the answer is still no.
UNIMPLEMENTED_RAW_STAGES = (
    # These two exist.  They are listed here because they are absent *from
    # this import graph*: the entry point cannot reach either, no verdict from
    # either reaches a record produced here, and a raw run therefore still does
    # not have them.  An oracle a run could import would not be an oracle.
    "the stationary killed-diffusion oracle and the endpoint convergence "
    "gate, which live outside this import graph and are unreachable from it",
    "the diagnostic moving-band Brownian-bridge audit, likewise outside this "
    "import graph and unreachable from it",
    "survival, hazard and the binomial complementary-log-log estimator, which "
    "live in the read-only analysis layer outside this import graph and are "
    "unreachable from it; a run written here carries no fitted anything",
    "every isolated comparison curve — the eligible-clock count, the bare "
    "relaxation-rate sum, the continuum control and the constant-rate "
    "control — likewise outside this import graph, computed only after a run "
    "has closed and hashed",
    "the throughput and memory benchmark, the power calculation and the "
    "finite production run matrix, which live in the feasibility layer "
    "outside this import graph and are unreachable from it; what exists there "
    "is a matrix *proposal* awaiting a human checkpoint, not a schedule",
    "the pilot firewall and its count-only range rule, likewise outside this "
    "import graph; no pilot has been run, and neither a signed production "
    "manifest nor a production coupling sweep exists",
    "the production numerical budget, which ticket 07 froze and which is "
    "**not met**; ticket 08 is not implemented in this package.  Nothing "
    "written here has passed one, and every manifest written here says so",
)

#: What a record produced through this entry point does not license anybody to
#: say.  Both lists, joined rather than restated, so the one-clock records and
#: the channel outcomes cannot drift into two different sets of disclaimers.
RAW_NON_CLAIMS = raw_experiments.NON_CLAIMS + RACE_NON_CLAIMS

#: The stream-namespace prefix reserved for **experiment pilot runs**.
#:
#: Raw-owned on purpose.  ``raw_runner`` must be able to recognize a pilot
#: identity without importing ``experiments`` — the import graph is one-way and
#: stays that way — so the authority is this constant rather than prose or a
#: label.  ``run_label = "pilot"`` is deliberately *not* the authority: ordinary
#: historical pilot-labelled fixtures predate ticket 07 and must keep working.
#:
#: Two consequences, both enforced below:
#:
#: * :func:`write_raw_run` **refuses** a configuration on this namespace, so a
#:   quarantined pilot can never be created under the production results root,
#:   not even for an instant; and
#: * :func:`open_raw_run` **refuses** a closed run whose manifest carries it, so
#:   a run manually placed in the production root cannot be read through the
#:   production door either.
#:
#: The only writer that accepts it is the private scoped one, which is handed an
#: already-opened parent descriptor by the caller that owns the quarantine.
EXPERIMENT_PILOT_NAMESPACE_PREFIX = "xpilot-"

#: Where a run is written: a package-scoped directory that ``.gitignore``
#: already excludes.  Generated output never lands beside source.
RESULTS_DIRNAME = "results"

#: The five files one closed run consists of, by role.  The close marker is
#: written last and its *presence* is half of what it means: a run that was
#: interrupted has no marker at all, so a reader knows the run is unfinished
#: without having to hash anything.
RUN_FILENAMES = {"manifest": "manifest.json",
                 "ledger": "ledger.csv",
                 "clocks": "clocks.csv",
                 "shadow": "shadow.csv",
                 "marker": "CLOSED.json"}

#: The predeclared rule that selects the no-stop shadow sample, written into
#: every manifest.  A rule rather than a list, for the same reason the clock
#: grid is a rule rather than a vector of detunings: a list could have been
#: chosen after seeing which mismatches looked quiet.
SHADOW_RULE = ("the first shadow_trials trial identifiers of this run "
               "re-walked on the same keys with the first-winner stop "
               "switched off; declared before the run and diagnostic only")

#: What the stationary continuum gate concluded, as of ticket 04's closure.
STATIONARY_ORACLE_STATUS = "passed at the declared reference tolerances"

#: What the moving-band audit concluded, as of ticket 04's closure.  Recorded
#: verbatim rather than summarised as "checked", because the audit's own
#: verdict is that this sample cannot bound the discrepancy — which is neither
#: a pass nor evidence that the discrepancy is small.
MOVING_BAND_STATUS = ("numerical_no_result pooled and in fourteen of fifteen "
                      "per-regime cells at the frozen reference sample")

#: The one field a later reader is most likely to want and least likely to
#: ask for.  A run written by this module is diagnostic: it exercises the race
#: and the ledger, and it has not passed a production numerical budget.  One
#: exists -- ticket 07 froze it -- and it is **not met**: the disposition is a
#: numerical no-result and no matrix was feasible among the options searched.
NUMERICAL_GATE = "diagnostic_only"

#: The content revision of :data:`VALIDATION_NOTE`.  Bumped when the sentence
#: changes meaning rather than wording, and named in the note itself so that a
#: manifest read on its own says which revision it carries.
#:
#: Revision 2 corrects a statement that had become false.  Revision 1 asserted
#: that no such budget existed; ticket 07 froze one, tied
#: it to planned statistical uncertainty, evaluated ticket 04's measured
#: envelopes against it, and concluded that it is not met.  An independent
#: review judged that a manifest statement known to be false is worse than a
#: changed digest, and that pre-checkpoint — with no durable pilot or
#: production run in existence — is the least disruptive moment to fix it.
#: Every manifest written from here therefore carries new bytes and a new
#: digest.  **The schema is unchanged**: 45 keys, version 3, same layout.
VALIDATION_NOTE_REVISION = 2

#: The long form of the three status fields above, so that a manifest read on
#: its own says what it means rather than requiring this file beside it.
VALIDATION_NOTE = (
    "Ticket 04's implementation is closed and its machinery checks pass, but "
    "its scientific comparison is a frozen numerical_no_result: the pooled "
    "moving-band ladder and fourteen of its fifteen per-regime cells cannot "
    "bound the discrepancy below their frozen caps at the reference sample. "
    "The correct reading is that this sample cannot bound the moving-band "
    "discrepancy, not that the discrepancy is small and not that it is large. "
    "This run is therefore diagnostic. It exercises the race and the ledger; "
    "it has not passed a production numerical budget. Ticket 07 froze such a "
    "budget, tied it to the planned production statistical uncertainty, and "
    "evaluated ticket 04's measured envelopes against it; the disposition is "
    "numerical_no_result and no matrix was feasible among the options it "
    "searched, so no run has passed a production budget and none can until "
    "that block is lifted. Validation note revision 2. No survival curve, no "
    "commitment law, no coupling exponent and no scaling claim may be drawn "
    "from these tables.")


def validate_raw_config(config) -> RawEventConfig:
    """Return a freshly built, exact :class:`~.raw_config.RawEventConfig`.

    The single input boundary of the raw process: a raw run accepts its inputs
    here and nowhere else.  Two separate things are guaranteed, and the second
    was missing until an independent review found it.

    **Nothing stale.**  ``RawEventConfig`` is frozen against ordinary
    assignment and not against ``object.__setattr__``.  A configuration
    validated once at construction and then edited in place would carry, say, a
    NaN peak coupling into a run — where it makes ``K > |Delta|`` false
    everywhere and produces a clean, physical-looking, entirely fictitious "no
    clock was ever admitted" result.  So this rebuilds, re-running every field
    check and every cross-field rule at the moment the run begins, which is the
    moment that matters.  The grid is rebuilt for the same reason: re-running
    the outer constructor alone would check only that the grid is a grid, not
    that it is still a valid one.

    **Nothing extra.**  The rebuild used to be ``dataclasses.replace``, which
    preserves the *dynamic* class.  A frozen subclass adding
    ``predictor: object`` therefore passed the boundary, was rebuilt as the
    same subclass, and arrived at the raw process still holding its predictor —
    an arbitrary opaque object, reaching the code that exists specifically to
    not have one.  Every other route for that object is closed: unknown keys,
    callables, arrays, non-string keys and about forty forbidden names are all
    refused by ``raw_config.raw_event_config``.  Subclassing was the way round
    all of it.

    So the type must be *exactly* ``RawEventConfig``, and the rebuild names its
    twelve fields one at a time rather than copying whatever the object
    happens to hold.  Both halves matter: the exact-type check stops the
    subclass, and the explicit field list means that even if some future edit
    loosened that check, an undeclared field would still not be copied.  A
    field added to the real dataclass and forgotten here is caught by
    ``verify.py``, which requires the rebuilt configuration to equal the
    original.

    ``TypeError`` if this is not exactly a ``RawEventConfig`` — a mapping, a
    namespace object, or a subclass of any kind.  ``ValueError`` if it is one
    but no longer holds valid values.
    """
    require_exact_type(config, RawEventConfig, "config")
    require_exact_type(config.grid, RawClockGrid, "config.grid")
    grid = RawClockGrid(config.grid.half_width, config.grid.n_points)
    return RawEventConfig(
        grid=grid,
        peak_coupling=config.peak_coupling,
        pulse_duration=config.pulse_duration,
        phase_diffusion=config.phase_diffusion,
        lock_tolerance=config.lock_tolerance,
        dwell_time=config.dwell_time,
        timestep=config.timestep,
        trials=config.trials,
        stream_namespace=config.stream_namespace,
        run_label=config.run_label,
        pulse_centre=config.pulse_centre,
        shadow_trials=config.shadow_trials,
    )


def raw_noise_stream(config, origin) -> PhaseNoiseStream:
    """Return the physical noise stream a validated raw configuration names.

    The one supported way to get from a run's inputs to its randomness, and it
    exists so that the connection is short enough to read: the dataset
    namespace, the finest timestep and the phase-diffusion strength come from
    the configuration and from nowhere else, and there is no argument through
    which a coupling, a mismatch, an amplitude or a predicted rate could reach
    the noise law.  The configuration is re-validated here rather than trusted,
    for the reason :func:`validate_raw_config` gives.

    ``origin`` is the time of finest-grid index zero.  It is explicit because
    the finest mesh is part of every key: two runs that agree on their
    timestep and disagree on where step zero begins are different datasets and
    must not silently share a realization.

    Note what is *not* passed: the run's coupling.  Reusing this stream across
    a coupling sweep is the intended use, and it is what makes a paired sweep
    paired.
    """
    checked = validate_raw_config(config)
    mesh = FinestMesh(require_finite(origin, "origin"), checked.timestep)
    return PhaseNoiseStream(checked.stream_namespace, mesh,
                            checked.phase_diffusion)


def raw_lock_criterion(config) -> LockCriterion:
    """Return the lock band and dwell a validated raw configuration names.

    The companion of :func:`raw_noise_stream` on the other side of the run:
    that one is the only route from a configuration to its randomness, and this
    is the only route from a configuration to its commitment rule.

    Both numbers come from the configuration and from nowhere else, and the
    configuration cannot hold a coupling-dependent one — it holds a single
    tolerance and a single dwell for the whole sweep.  That is the structural
    half of "dwell and tolerance never scale with coupling"; the verifier's
    mutation battery is the behavioural half.
    """
    checked = validate_raw_config(config)
    return LockCriterion(checked.lock_tolerance, checked.dwell_time)


def raw_one_clock_path(config, detuning, model: str = "full",
                       drive: str = "pulsed") -> ClockPath:
    """Return the one-clock data-generating process a configuration names.

    The drive is rebuilt from the validated configuration, so the coupling
    history a run integrates is the one its manifest declares: a pulsed run
    takes the peak, duration and centre and gets the calibrated raised-cosine
    tongue geometry; a stationary run takes the peak alone and is eligible for
    all time or none of it.  ``detuning`` is separate because it belongs to a
    clock rather than to the run: the configuration carries the whole frozen
    grid, and a one-clock experiment picks one mismatch out of it — or, for the
    central-clock control, uses exactly zero.

    **There is no ``fixed_rate`` argument, and its absence is the contract.**
    This factory used to take one, so the same configuration could be handed
    ``0.1``, ``1.0`` and ``10.0`` and would return three different width-only
    controls — which is precisely the per-run tuning the control exists to
    exclude.  The rate is now derived in one place,
    :data:`dynamics.FIXED_CONTRACTION_RATE`, from the frozen rule in
    :func:`dynamics.fixed_contraction_rate`: the central clock's Adler
    relaxation rate at the geometric midpoint of the declared reference
    peak-coupling range.  It is common to every clock, coupling, mismatch,
    trial, dwell, tolerance and noise strength, and there is no argument here
    through which a caller can vary it.

    No sensitivity factory is provided.  Lower and upper fixed-rate choices are
    a later ticket's *labelled* sensitivity against a frozen production
    manifest; an unlabelled one now would be the bypass this change removes,
    wearing a different name.

    Note what else cannot be passed: there is no argument for a predicted rate,
    a weight, a hazard, or a callback of any kind, and the configuration
    boundary refuses all of those before this function is ever reached.
    """
    checked = validate_raw_config(config)
    if model not in MODEL_LABELS:
        raise ValueError(f"model must be one of {list(MODEL_LABELS)}, got "
                         f"{model!r}")
    if drive not in DRIVE_LABELS:
        raise ValueError(f"drive must be one of {list(DRIVE_LABELS)}, got "
                         f"{drive!r}")
    if drive == "stationary":
        source = StationaryDrive(checked.peak_coupling)
    else:
        source = PulseTrain((RaisedCosinePulse(checked.peak_coupling,
                                               checked.pulse_duration,
                                               checked.pulse_centre),))
    rate = FIXED_CONTRACTION_RATE if model == "width_only" else 0.0
    return ClockPath(require_finite(detuning, "detuning"), source, model, rate)


def raw_results_directory(run_name):
    """Return the package-scoped path a run called ``name`` belongs at.

    ``<package>/results/<name>``, and the directory is not created here: this
    names a place, and :func:`write_raw_run` is what fills one.

    The plan says generated output goes only beneath a package-scoped ignored
    results directory, and this is the only function that knows where that is.
    ``.gitignore`` in this package already excludes ``results/`` along with
    every CSV and JSON below it, so a run written through here cannot be
    committed beside the source by accident.

    ``run_name`` is one path segment.  A separator, a parent reference, an absolute
    path or a leading dot is refused rather than resolved, because every one of
    them is a way for "a run called X" to land somewhere other than under this
    directory — and the whole value of a scoped output location is that nothing
    can quietly leave it.

    :func:`write_raw_run` takes a run name and calls this; there is no public
    writer that accepts a directory.  An earlier revision had one, and an
    independent review correctly called it a contradiction of the contract this
    helper was declaring — the scoped path existed beside a writer that did not
    use it, and every verification left complete run fixtures outside the
    package.  A caller who genuinely needs to write elsewhere is doing
    something the plan does not describe and should say so in its own module.
    """
    # The parameter is ``run_name`` rather than ``name`` because ``name`` is
    # this package's convention for the *label* a validator quotes in its own
    # error message, and a real argument wearing that name would be exempted
    # from the exported-surface probe as though it were one.
    if isinstance(run_name, bool) or not isinstance(run_name, str):
        raise TypeError(f"run_name must be a string, got "
                        f"{type(run_name).__name__} {run_name!r}")
    if not run_name:
        raise ValueError("run_name must not be empty")
    bad = sorted(set(run_name) & set("/\\:"))
    if bad or run_name.startswith(".") or Path(run_name).is_absolute() \
            or len(Path(run_name).parts) != 1:
        raise ValueError(
            f"a run name is one path segment, got {run_name!r}; a separator, "
            f"a "
            f"parent reference, a leading dot or an absolute path would put "
            f"the run somewhere other than under the package's own ignored "
            f"results directory, which is the only thing that makes the "
            f"location worth having")
    return Path(__file__).resolve().parent / RESULTS_DIRNAME / run_name


def raw_run_window(config):
    """Return ``(mesh, start_index, steps)``: the exact grid a run walks.

    Derived from the configuration and from nothing else, so that the window is
    a consequence of the frozen pulse rather than a separate number a caller
    could have chosen.  The mesh origin is the opening of the pulse support,
    step zero begins there, and the run covers the support in whole finest
    steps — rounded *up*, so that the last step never falls short of the pulse
    end.  Overshooting costs nothing: the tongue has closed, no clock is
    eligible, and no dwell can start.  Falling short would silently truncate
    the experiment's own interaction.

    Returning the mesh rather than only the origin is deliberate.  A caller
    that rebuilt the mesh from the origin would have a second spelling of it,
    and every noise key contains the mesh identity: two spellings that rounded
    differently would be two datasets under one name.
    """
    checked = validate_raw_config(config)
    origin = checked.pulse_centre - 0.5 * checked.pulse_duration
    steps = int(math.ceil(checked.pulse_duration / checked.timestep))
    return FinestMesh(origin, checked.timestep), 0, max(steps, 1)


def raw_clock_paths(config, model: str = "full"):
    """Return the finite population a configuration names, in clock order.

    One :class:`~.dynamics.ClockPath` per clock, at the exact detunings the
    frozen grid rule produces, all sharing one pulse train.  The clock at index
    ``i`` has detuning ``grid.detunings[i]``, and that correspondence is the
    only thing anywhere that says which clock is which: a clock identifier is
    half of every key the clock consumes, so a population built in a different
    order would be a different experiment wearing the same manifest.

    Every clock is a real competitor.  Changing ``n_points`` changes the number
    of physical commitment opportunities and is a **population sensitivity**,
    never a convergence study, and the manifest records the grid rule so that a
    later reader cannot mistake one for the other.

    There is no argument for a per-clock weight, a rate, or a contraction rate.
    The width-only control's rate comes from :data:`dynamics.FIXED_CONTRACTION_RATE`
    and from nowhere else, exactly as it does for one clock.
    """
    checked = validate_raw_config(config)
    if model not in MODEL_LABELS:
        raise ValueError(f"model must be one of {list(MODEL_LABELS)}, got "
                         f"{model!r}")
    train = PulseTrain((RaisedCosinePulse(checked.peak_coupling,
                                          checked.pulse_duration,
                                          checked.pulse_centre),))
    rate = FIXED_CONTRACTION_RATE if model == "width_only" else 0.0
    return tuple(ClockPath(float(detuning), train, model, rate)
                 for detuning in checked.grid.detunings)


def raw_clock_population(config, model: str = "full") -> ClockPopulation:
    """Return the declared population a configuration names, bound to its grid.

    :func:`raw_clock_paths` returns the paths; this returns the paths *together
    with the declaration they are supposed to satisfy* — the population size,
    the clock order, and every clock's mismatch, all taken from the frozen grid
    rule rather than from the paths themselves.

    That distinction is the whole of an independent review's third finding.
    Handed a bare tuple, the race could be given the first path of a declared
    two-clock grid and would run, return an ordinary outcome, and claim a
    category for a population it did not have — because a trailing omission
    leaves a complete range of clock identifiers, just of a smaller population.
    A declaration cannot be truncated by dropping a path, so this is the object
    the race takes.

    The identity also carries the drive's digest, so two populations with
    identical mismatches under different pulses are different experiments
    rather than the same one.
    """
    checked = validate_raw_config(config)
    paths = raw_clock_paths(checked, model)
    identity = PopulationIdentity(
        model=model,
        pulse_digest=paths[0].train.digest,
        detunings=tuple(float(value) for value in checked.grid.detunings))
    return ClockPopulation(paths, identity)


def raw_manifest(config, model: str = "full") -> RawManifest:
    """Return the frozen description a run of this configuration writes.

    Everything a later comparison needs in order to know what it is comparing,
    and nothing a run could have chosen after seeing its own results.  The
    grid rule and its realized detunings, the pulse, the noise, the lock band,
    the dwell, the mesh and the exact step window, the key namespaces, the
    column layout of all three tables, the software versions, and the three
    numerical-status fields that say this run is diagnostic.

    Two absences are deliberate.  There is no timestamp, and there are no
    execution parameters — no clock-block size, no step-window size.  A
    timestamp would make two runs of one frozen configuration differ in their
    digests while being the same run; execution parameters bound the live
    working set and cannot change a recorded number, so freezing them into the
    immutable description would imply that they were part of it.  Both
    absences are what make "the same configuration writes byte-identical
    tables" a checkable statement rather than an aspiration.
    """
    checked = validate_raw_config(config)
    mesh, start_index, steps = raw_run_window(checked)
    stream = PhaseNoiseStream(checked.stream_namespace, mesh,
                              checked.phase_diffusion)
    grid = checked.grid
    half = 0.5 * checked.pulse_duration
    return RawManifest({
        "schema_version": MANIFEST_SCHEMA_VERSION,
        "ledger_schema_version": LEDGER_SCHEMA_VERSION,
        "run_label": checked.run_label,
        "model": model if model in MODEL_LABELS else str(model),
        "drive": "pulsed",
        "trials": checked.trials,
        "shadow_trials": checked.shadow_trials,
        "shadow_rule": SHADOW_RULE,
        "clocks": grid.n_points,
        "grid_half_width": grid.half_width,
        "grid_spacing": grid.spacing,
        "grid_parity": grid.parity,
        "grid_origin": 0.0,
        "grid_support": [-grid.half_width, grid.half_width],
        "detunings": [float(value) for value in grid.detunings],
        "peak_coupling": checked.peak_coupling,
        "pulse_duration": checked.pulse_duration,
        "pulse_centre": checked.pulse_centre,
        "pulse_support": [checked.pulse_centre - half,
                          checked.pulse_centre + half],
        "pulse_digest": PulseTrain((RaisedCosinePulse(
            checked.peak_coupling, checked.pulse_duration,
            checked.pulse_centre),)).digest,
        "phase_diffusion": checked.phase_diffusion,
        "lock_tolerance": checked.lock_tolerance,
        "dwell_time": checked.dwell_time,
        "fixed_contraction_rate": (FIXED_CONTRACTION_RATE
                                   if model == "width_only" else 0.0),
        "timestep": checked.timestep,
        "mesh_origin": mesh.origin,
        "mesh_identity": mesh.identity,
        "start_index": start_index,
        "steps": steps,
        "window_open": mesh.time_at(start_index),
        "window_close": mesh.time_at(start_index + steps),
        "key_schema": KEY_SCHEMA,
        "stream_namespace": checked.stream_namespace,
        "stream_prefix": stream.stream_prefix,
        "ensemble_prefix": stream.ensemble_prefix,
        "ledger_columns": [name for name, _ in LEDGER_COLUMNS],
        "clock_columns": [name for name, _ in CLOCK_COLUMNS],
        "shadow_columns": [name for name, _ in SHADOW_COLUMNS],
        "endpoint_status": STATIONARY_ORACLE_STATUS,
        "moving_band_status": MOVING_BAND_STATUS,
        "numerical_gate": NUMERICAL_GATE,
        "validation_note": VALIDATION_NOTE,
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
        "non_claims": list(RAW_NON_CLAIMS),
    })


@dataclass(frozen=True)
class RawRunReport:
    """What :func:`write_raw_run` did, for the caller and for nobody else.

    Deliberately *not* part of the immutable record.  The execution parameters
    and the category tally are here because a caller wants them; they are not
    in the manifest because the block sizes cannot change a result and the
    tally is derivable from the ledger the manifest already binds.  A number
    that lives in two places is a number two readers can disagree about.

    ``directory``
        Where the five files were written.
    ``marker``
        The close marker, as written.
    ``categories``
        How many trials landed in each of the four channel categories.  They
        sum to the trial count, always: no trial is conditioned away.
    ``clock_block`` / ``step_window``
        The live-block sizes this execution used.
    """

    directory: object
    marker: CloseMarker
    categories: tuple
    clock_block: int
    step_window: int

    def __post_init__(self):
        require_exact_type(self.marker, CloseMarker, "marker")
        if not isinstance(self.directory, (str, Path)):
            raise TypeError(f"directory must be a path or a string, got "
                            f"{type(self.directory).__name__}")
        tally = []
        for position, entry in enumerate(self.categories):
            pair = tuple(entry)
            if len(pair) != 2:
                raise ValueError(f"categories[{position}] must be a "
                                 f"(category, count) pair, got {pair!r}")
            if pair[0] not in TRIAL_CATEGORIES:
                raise ValueError(f"categories[{position}] names {pair[0]!r}, "
                                 f"which is not one of "
                                 f"{list(TRIAL_CATEGORIES)}")
            tally.append((pair[0],
                          require_positive_int(pair[1],
                                               f"categories[{position}] count",
                                               minimum=0)))
        if sorted(name for name, _ in tally) != sorted(TRIAL_CATEGORIES):
            raise ValueError(
                f"the tally names {[name for name, _ in tally]}; every one of "
                f"{list(TRIAL_CATEGORIES)} appears in a run's tally, because a "
                f"category with no trials in it is a zero rather than an "
                f"absence and the four together are what makes the tally "
                f"exhaustive")
        object.__setattr__(self, "categories", tuple(tally))
        for name in ("clock_block", "step_window"):
            object.__setattr__(
                self, name,
                require_positive_int(getattr(self, name), name, minimum=1))

    @property
    def trials(self) -> int:
        """How many trials the run wrote.  The tally's own sum."""
        return sum(count for _, count in self.categories)


class _Table:
    """A CSV table being written, hashed as it is written.

    The digest is folded in one line at a time rather than computed from the
    finished file, so nothing ever holds the whole table in memory and the
    marker's digest is over exactly the bytes that reached the disk.

    It is created **relative to a pinned run-directory descriptor**, and
    exclusively: nothing may already stand at that name, including a symbolic
    link, including a broken one.
    """

    def __init__(self, run_fd, filename, columns):
        self.run_fd = run_fd
        self.filename = filename
        self.columns = columns
        self.rows = 0
        self._digest = None
        self._hash = None
        self._handle = None

    def __enter__(self):
        self._hash = hashlib.new(DIGEST_ALGORITHM)
        self._handle = _create_run_file(self.run_fd, self.filename,
                                        f"the {self.filename} table")
        self._write(header_row(self.columns))
        return self

    def _write(self, line: str) -> None:
        payload = (line + "\n").encode("utf-8")
        self._handle.write(payload)
        self._hash.update(payload)

    def add(self, values) -> None:
        """Append one row, in the declared column order."""
        self._write(encode_row(self.columns, values))
        self.rows += 1

    def __exit__(self, kind, value, trace):
        self._handle.close()
        self._handle = None
        if kind is None:
            self._digest = self._hash.hexdigest()
        return False

    @property
    def digest(self) -> str:
        """The digest of the bytes written, available once closed."""
        if self._digest is None:
            raise ValueError(f"{self.filename} has not been closed, so its "
                             f"digest is not the digest of a finished table")
        return self._digest


def _trial_values(outcome):
    """Return the trial ledger cells one channel outcome spells.

    Every multi-valued cell is indexed by the same winner list, in clock
    order, so that a reader zipping ``winners`` against ``winner_detunings``
    is zipping one clock's fields together.  Nothing is dropped when there is
    more than one winner and nothing is invented when there are none.

    ``winner_last_reset`` is here because an independent review found it
    missing: the plan's committed record asks for the winner's tongue-entry,
    band-entry *and* last-reset times, and a run that dropped the third lost
    permanently the one number that says whether the winning dwell was the
    clock's first attempt or its fifth.
    """
    winners = [outcome.clocks[clock] for clock in outcome.winners]
    return {
        "trial": outcome.trial,
        "category": outcome.category,
        "stop_time": outcome.stop_time,
        "commit_time": outcome.commit_time,
        "censored": outcome.censored,
        "co_completion": outcome.co_completion,
        "winners": [w.clock for w in winners],
        "winner_detunings": [w.detuning for w in winners],
        "winner_phases": [w.final_phase for w in winners],
        "winner_band_entry": [w.first_inside_at for w in winners],
        "winner_first_eligible": [w.first_eligible_at for w in winners],
        "winner_last_reset": [w.last_reset_at for w in winners],
        "clocks_ever_eligible": outcome.clocks_ever_eligible,
        "clocks_ever_inside": outcome.clocks_ever_inside,
        "band_entries": outcome.band_entries,
        "dwell_resets": outcome.dwell_resets,
        "eligible_time": outcome.eligible_time,
    }


def _clock_values(trial, exposure, co_completion):
    """Return the authoritative per-``(trial, clock)`` cells for one exposure.

    The durable record of what this clock was doing when the channel stopped.
    An earlier revision wrote an aggregate instead — one row per clock, totals
    across the run — and an independent review pointed out that it could not
    say which trial contributed an entry, a reset, an exposure or a win, so a
    closed run could never afterwards support the censoring analysis it exists
    for.

    ``co_completion`` is the trial's winner count, and is passed in rather than
    read off the exposure because whether a clock *shared* its win is a fact
    about the channel rather than about the clock.
    """
    return {
        "trial": trial,
        "clock": exposure.clock,
        "detuning": exposure.detuning,
        "category": exposure.category,
        "committed_at": exposure.committed_at,
        "first_eligible_at": exposure.first_eligible_at,
        "first_inside_at": exposure.first_inside_at,
        "last_reset_at": exposure.last_reset_at,
        "eligibility_entries": exposure.eligibility_entries,
        "band_entries": exposure.band_entries,
        "dwell_resets": exposure.dwell_resets,
        "endpoints": exposure.endpoints,
        "eligible_endpoints": exposure.eligible_endpoints,
        "inside_endpoints": exposure.inside_endpoints,
        "exposure_time": exposure.exposure,
        "eligible_time": exposure.eligible_time,
        "final_phase": exposure.final_phase,
        "won": exposure.committed,
        "co_completed": exposure.committed and co_completion > 1,
    }


def _shadow_values(trial, exposure):
    """Return the no-stop cells for one uncensored clock exposure.

    The same fields as the censored row minus the two that only a race can
    have, so that the two tables subtract cleanly.  A commitment recorded here
    is one that in the actual race did not happen, because the channel had
    already stopped.
    """
    values = _clock_values(trial, exposure, 0)
    return {name: values[name] for name, _ in SHADOW_COLUMNS}


#: The bits ``os.open`` needs to open a directory, and to refuse a symbolic
#: link at the final component.
_DIRECTORY_FLAGS = os.O_RDONLY | getattr(os, "O_DIRECTORY", 0)
_NOFOLLOW = getattr(os, "O_NOFOLLOW", 0)

#: ``st_mode`` masks for "is this a regular file", spelled out rather than
#: imported.  ``stat`` would be a third module on the raw allowlist to answer
#: one question, and these two constants have not moved since POSIX.
_FILE_TYPE_MASK = 0o170000
_REGULAR_FILE = 0o100000

#: Permissions a written run file gets.  Owner-writable, world-readable: a
#: result is meant to be read.
_FILE_MODE = 0o644
_DIRECTORY_MODE = 0o755


def _require_scoped_io():
    """Raise unless this platform can pin a directory and refuse a link.

    Fail closed.  The whole point of the descriptor chain below is that a name
    is resolved once and everything afterwards happens relative to the thing
    that resolution opened; on a platform without ``O_NOFOLLOW``,
    ``O_DIRECTORY`` or ``dir_fd`` support, that guarantee simply does not
    exist, and running anyway would mean claiming a property the code cannot
    deliver.
    """
    missing = [name for name in ("O_DIRECTORY", "O_NOFOLLOW", "O_CREAT",
                                 "O_EXCL") if not getattr(os, name, 0)]
    if os.open not in os.supports_dir_fd:
        missing.append("os.open(dir_fd=...)")
    if os.mkdir not in os.supports_dir_fd:
        missing.append("os.mkdir(dir_fd=...)")
    if missing:
        raise ValueError(
            f"this platform is missing {missing}, so a run directory cannot be "
            f"pinned and a symbolic link at a final component cannot be "
            f"refused.  Raw runs are not written or read without that, because "
            f"the alternative is to validate a name and then reopen it — which "
            f"is exactly the gap this refuses to have")


def _open_results_root():
    """Return a descriptor for the package's results directory.

    Opened *relative to the package directory*, which comes from ``__file__``
    and is therefore fixed at import.  The only name resolved by string is the
    package's own location; ``results`` and everything below it are reached
    through descriptors, so nothing that happens afterwards — a rename, a
    replacement, a symbolic link dropped in place — can redirect the
    operations that follow.

    ``O_NOFOLLOW`` on the results directory itself, so a link standing in for
    it is refused rather than followed.
    """
    _require_scoped_io()
    package_fd = os.open(str(Path(__file__).resolve().parent),
                         _DIRECTORY_FLAGS)
    try:
        try:
            os.mkdir(RESULTS_DIRNAME, _DIRECTORY_MODE, dir_fd=package_fd)
        except FileExistsError:
            pass
        try:
            return os.open(RESULTS_DIRNAME, _DIRECTORY_FLAGS | _NOFOLLOW,
                           dir_fd=package_fd)
        except OSError as exc:
            raise ValueError(
                f"the package results directory could not be opened as a real "
                f"directory ({exc.strerror}); a symbolic link standing in for "
                f"it is refused rather than followed") from None
    finally:
        os.close(package_fd)


def _open_run_directory(results_fd, run_name, create: bool) -> int:
    """Return a descriptor for one run directory beneath ``results_fd``.

    ``O_NOFOLLOW | O_DIRECTORY``, so a symbolic link at ``results/<run>`` is
    refused whatever it points at, and a plain file under that name is refused
    too.  An independent review replaced a run directory *after* a lexical
    check had passed and redirected both the reader and the writer; there is
    nothing to replace here, because the descriptor is what the rest of the
    operation uses.
    """
    if create:
        try:
            os.mkdir(run_name, _DIRECTORY_MODE, dir_fd=results_fd)
        except FileExistsError:
            pass                    # every file is created exclusively below
    try:
        return os.open(run_name, _DIRECTORY_FLAGS | _NOFOLLOW,
                       dir_fd=results_fd)
    except OSError as exc:
        raise ValueError(
            f"the run directory {run_name!r} could not be opened as a real "
            f"directory beneath the package results root "
            f"({exc.strerror}).  A symbolic link, a plain file, or a missing "
            f"directory all land here, and a link is refused rather than "
            f"followed") from None


def _read_run_file(run_fd, filename, label):
    """Return the bytes of one run file, opened relative to ``run_fd``.

    ``O_NOFOLLOW``, so an individual table name that is a symbolic link to an
    external file is refused — the reproduction an independent review used to
    read twelve rows from outside the package through five linked filenames.
    The opened object must also be a regular file, checked through ``fstat`` on
    the descriptor rather than on the name.
    """
    try:
        handle = os.open(filename, os.O_RDONLY | _NOFOLLOW, dir_fd=run_fd)
    except OSError as exc:
        raise ValueError(
            f"{label} could not be read from the run directory "
            f"({exc.strerror}); a symbolic link at that name is refused "
            f"rather than followed, and a missing file means the run is not "
            f"the five files a closed run has") from None
    with os.fdopen(handle, "rb") as stream:
        info = os.fstat(stream.fileno())
        if info.st_mode & _FILE_TYPE_MASK != _REGULAR_FILE:
            raise ValueError(f"{label} is not a regular file")
        return stream.read()


def _create_run_file(run_fd, filename, label):
    """Return a binary stream for one new run file, relative to ``run_fd``.

    ``O_CREAT | O_EXCL | O_NOFOLLOW``.  ``O_EXCL`` refuses anything already at
    that name — including a *broken* symbolic link, which ``Path.exists()``
    reports as absent and which an independent review used to make the writer
    create all five files outside the package.
    """
    try:
        handle = os.open(filename,
                         os.O_WRONLY | os.O_CREAT | os.O_EXCL | _NOFOLLOW,
                         _FILE_MODE, dir_fd=run_fd)
    except FileExistsError:
        raise ValueError(
            f"{label} already exists in this run directory.  A raw run writes "
            f"a fresh one rather than overwriting, and this refusal covers a "
            f"symbolic link left at that name — including a broken one, which "
            f"an ordinary existence test reports as absent") from None
    except OSError as exc:
        raise ValueError(f"{label} could not be created ({exc.strerror})") \
            from None
    return os.fdopen(handle, "wb")


def _require_inside_results(target):
    """Return the resolved run directory, or raise if it is not scoped.

    The plan says generated output goes only beneath a package-scoped ignored
    results directory.  Version 2's writer accepted an arbitrary directory and
    said so in its own documentation, which an independent review correctly
    called a contradiction of the contract: the scoped path helper existed
    beside a writer that did not use it.

    The check is on the **resolved** path, so a symbolic link that points out
    of the results root is refused rather than followed.  ``Path.resolve``
    handles the ``..`` and link cases together, and the comparison is
    ``is_relative_to`` against the resolved root — a prefix test on the
    unresolved strings would accept ``results/../../elsewhere``.
    """
    root = Path(__file__).resolve().parent / RESULTS_DIRNAME
    resolved = Path(target).resolve()
    if resolved == root or not resolved.is_relative_to(root):
        raise ValueError(
            f"{resolved} is not a run directory beneath {root}.  Generated "
            f"output goes only under the package's own ignored results "
            f"directory: an absolute path, a parent traversal and a symbolic "
            f"link that leaves the root all resolve outside it and are "
            f"refused rather than followed")
    return resolved


def write_raw_run(config, run_name, model: str = "full",
                  clock_block: int = DEFAULT_CLOCK_BLOCK,
                  step_window: int = DEFAULT_STEP_WINDOW) -> RawRunReport:
    """Run the finite one-channel race and close it into a named run directory.

    The launcher, and the only supported way to write a run.  ``run_name`` is
    one path segment; the run lands at
    ``<package>/results/<run_name>``, which :func:`raw_results_directory`
    names and this function creates.  There is deliberately **no argument for
    an arbitrary directory**: an independent review found version 2's writer
    accepting one, which contradicted the plan's package-scoped output
    contract and left complete run fixtures outside the package after every
    verification.

    It writes five files — the manifest, the trial ledger, the authoritative
    per-``(trial, clock)`` ledger, the no-stop shadow sample, and last of all
    the close marker — and returns a :class:`RawRunReport` describing what it
    did.

    **The marker is written last, and that ordering is the protocol.**  A run
    that was interrupted, killed, or is still going has no marker on disk at
    all, so a reader knows the run is unfinished before it hashes anything, and
    :func:`open_raw_run` refuses it for that reason rather than reporting a
    digest mismatch for a file that was simply not finished.  Everything the
    marker asserts is computed from the bytes as they were written: each table
    is hashed a line at a time while it is being written, so the digest is over
    exactly what reached the disk and never over a reconstruction of it.

    **Nothing is held.**  Trials are raced one at a time and their rows are
    written and dropped; the leaves of a clock's walk arrive in windows and are
    dropped as they are consumed.  Version 2 kept ``N`` running counters
    between trials; version 3 does not keep even those, because the per-clock
    aggregate is now derived from the durable rows at read time.  There is no
    array indexed by trials and none indexed by both clocks and steps.

    ``clock_block`` and ``step_window`` bound the live working set and cannot
    change a single recorded number.  They are arguments rather than
    configuration fields for exactly that reason, and they are absent from the
    manifest so that the same configuration written at any block size produces
    byte-identical files and an identical marker.

    A run name that already has files under it is refused.  Overwriting one in
    place is how a rerun's rows end up mixed with a previous run's, and the
    refusal is loud rather than a silent truncate.
    """
    checked = validate_raw_config(config)
    population = raw_clock_population(checked, model)
    manifest = raw_manifest(checked, model)
    mesh, start_index, steps = raw_run_window(checked)
    stream_source = PhaseNoiseStream(checked.stream_namespace, mesh,
                                     checked.phase_diffusion)
    criterion = raw_lock_criterion(checked)
    clock_block = require_positive_int(clock_block, "clock_block", minimum=1)
    step_window = require_positive_int(step_window, "step_window", minimum=1)

    # Name first, then descriptors.  The name check is cheap and gives a good
    # error; the *authority* is the chain below, in which the package
    # directory is the only thing resolved by string and everything under it
    # is reached through a pinned descriptor.
    if checked.stream_namespace.startswith(EXPERIMENT_PILOT_NAMESPACE_PREFIX):
        raise ValueError(
            f"stream_namespace {checked.stream_namespace!r} is on the reserved "
            f"experiment-pilot prefix {EXPERIMENT_PILOT_NAMESPACE_PREFIX!r}; a "
            f"pilot run is written directly into its quarantine by the scoped "
            f"writer and must never exist under the production results root, "
            f"not even for an instant")
    root = _require_inside_results(raw_results_directory(run_name))
    results_fd = _open_results_root()
    try:
        return _write_run_beneath(checked, run_name, results_fd, root, model,
                                  clock_block, step_window)
    finally:
        os.close(results_fd)


def _write_run_beneath(checked, run_name, parent_fd, directory, model,
                       clock_block, step_window) -> RawRunReport:
    """Create the run directory beneath ``parent_fd`` and write the five files.

    The shared body of both writers.  ``parent_fd`` stays owned by the caller;
    everything below it is opened relative to a descriptor this function pins,
    with ``O_EXCL | O_NOFOLLOW`` on creation and a regular-file check on read,
    and the close marker is written last.  Nothing is moved afterwards, which
    is what lets a quarantined pilot exist only where it belongs.
    """
    population = raw_clock_population(checked, model)
    manifest = raw_manifest(checked, model)
    mesh, start_index, steps = raw_run_window(checked)
    stream_source = PhaseNoiseStream(checked.stream_namespace, mesh,
                                     checked.phase_diffusion)
    criterion = raw_lock_criterion(checked)
    run_fd = _open_run_directory(parent_fd, run_name, create=True)
    try:
        manifest_bytes = canonical_manifest_bytes(manifest)
        with _create_run_file(run_fd, RUN_FILENAMES["manifest"],
                              "the manifest") as stream:
            stream.write(manifest_bytes)

        counted = {name: 0 for name in TRIAL_CATEGORIES}
        ledger = _Table(run_fd, RUN_FILENAMES["ledger"], LEDGER_COLUMNS)
        clocks = _Table(run_fd, RUN_FILENAMES["clocks"], CLOCK_COLUMNS)
        shadow = _Table(run_fd, RUN_FILENAMES["shadow"], SHADOW_COLUMNS)
        with ledger, clocks, shadow:
            for trial in range(checked.trials):
                outcome = race_one_channel(population, criterion, stream_source,
                                           trial, start_index, steps,
                                           clock_block, step_window)
                ledger.add(_trial_values(outcome))
                for exposure in outcome.clocks:
                    clocks.add(_clock_values(trial, exposure,
                                             outcome.co_completion))
                counted[outcome.category] += 1
                if trial < checked.shadow_trials:
                    for exposure in no_stop_exposure(population, criterion,
                                                     stream_source, trial,
                                                     start_index, steps,
                                                     clock_block, step_window):
                        shadow.add(_shadow_values(trial, exposure))

        tables = {"ledger": ledger, "clocks": clocks, "shadow": shadow}
        marker = CloseMarker(
            LEDGER_SCHEMA_VERSION,
            tuple((name, tables[name].rows, tables[name].digest)
                  for name in TABLE_NAMES),
            digest_bytes(manifest_bytes))
        # Last, and through the same pinned descriptor.
        with _create_run_file(run_fd, RUN_FILENAMES["marker"],
                              "the close marker") as stream:
            stream.write(canonical_marker_bytes(marker))
    finally:
        os.close(run_fd)

    return RawRunReport(directory=directory, marker=marker,
                        categories=tuple((name, counted[name])
                                         for name in sorted(counted)),
                        clock_block=clock_block, step_window=step_window)


def _write_scoped_run(config, run_name, parent_fd, directory, model="full",
                      clock_block=DEFAULT_CLOCK_BLOCK,
                      step_window=DEFAULT_STEP_WINDOW) -> RawRunReport:
    """Write one closed run under an **already-vetted parent descriptor**.

    Private, and the only door that accepts the reserved experiment-pilot
    namespace.  ``parent_fd`` is a directory descriptor the caller opened and
    owns; this function creates the run directory beneath it with the same
    ``O_EXCL | O_NOFOLLOW`` semantics, writes the same five files in the same
    order with the marker last, and returns the same report.  Nothing is moved
    afterwards, so no entry ever appears in the production results root.

    ``directory`` is the path the report should carry; it is used for the
    record and never resolved for I/O.
    """
    checked = validate_raw_config(config)
    if not checked.stream_namespace.startswith(
            EXPERIMENT_PILOT_NAMESPACE_PREFIX):
        raise ValueError(
            f"the scoped writer writes experiment pilots only; "
            f"stream_namespace {checked.stream_namespace!r} does not carry the "
            f"reserved prefix {EXPERIMENT_PILOT_NAMESPACE_PREFIX!r}")
    clock_block = require_positive_int(clock_block, "clock_block", minimum=1)
    step_window = require_positive_int(step_window, "step_window", minimum=1)
    return _write_run_beneath(checked, run_name, parent_fd, directory, model,
                              clock_block, step_window)


def _refuse_reserved_pilot(manifest):
    """Refuse a closed run whose manifest carries the reserved pilot identity."""
    namespace = manifest.value("stream_namespace")
    if isinstance(namespace, str) and namespace.startswith(
            EXPERIMENT_PILOT_NAMESPACE_PREFIX):
        raise ValueError(
            f"this run is on the reserved experiment-pilot namespace "
            f"{namespace!r}; the production reader does not open experiment "
            f"pilots, whichever directory somebody has placed one in")
    return manifest


def open_raw_run(run_name):
    """Return the verified :class:`~.raw_ledger.ClosedLedger` of a named run.

    The consumer's door, and the only supported way through it.  It reads the
    five files once, hands the bytes to
    :func:`~.raw_ledger.require_closed_ledger`, and returns what that gate
    returns — which is the parsed rows, not the bytes, so that a caller cannot
    reparse them with a looser parser than the one that was checked.

    Scope is a *capability*, not a check.  ``run_name`` must be one path
    segment; the results root and the run directory are then opened as pinned
    directory descriptors with ``O_NOFOLLOW``, and all five files are opened
    relative to the run descriptor, also with ``O_NOFOLLOW``, and required to
    be regular files.  A symbolic link standing in for the run directory, or
    for any individual table, is refused rather than followed — and nothing
    that happens after the directory is opened can redirect the reads, because
    there is no name left to replace.

    A missing close marker fails next, and fails plainly: it is the ordinary
    state of a run that is still going, and a reader that hashed the tables
    first would report a digest mismatch for a file that was simply not
    finished yet.

    Scoped the same way the writer is: ``run_name`` is one path segment under
    the package's results root.  A consumer holding bytes from somewhere else
    calls :func:`~.raw_ledger.require_closed_ledger` directly, which is the
    general door; this is the convenience for the runs this package writes.

    This lives beside the writer rather than in ``raw_ledger`` because it
    touches the filesystem and ``raw_ledger`` is a pure description of a
    format.  It does not make the raw process a reader of its own results: no
    raw run calls this, and the comparison layer that will is a later ticket
    outside this import graph.
    """
    # Name first, then descriptors — and the descriptors are the authority.
    #
    # An earlier revision resolved the path and then reopened it by string.  An
    # independent review replaced the run directory after that resolution and
    # read three trial rows from an external run; it also built an in-root
    # directory whose five filenames were symbolic links to external files and
    # read all twelve rows through them.  Resolving a name to another name pins
    # nothing.
    #
    # So: the package directory is the only thing resolved by string, the
    # results root and the run directory are opened beneath it with
    # ``O_NOFOLLOW``, and every one of the five files is opened relative to the
    # pinned run descriptor, also with ``O_NOFOLLOW``, and required to be a
    # regular file.  There is no name left for anything to replace.
    _require_inside_results(raw_results_directory(run_name))
    results_fd = _open_results_root()
    try:
        run_fd = _open_run_directory(results_fd, run_name, create=False)
    except BaseException:
        os.close(results_fd)
        raise
    try:
        names = {role: RUN_FILENAMES[role]
                 for role in ("marker", "manifest", "ledger", "clocks",
                              "shadow")}
        try:
            marker_bytes = _read_run_file(run_fd, names["marker"],
                                          "the close marker")
        except ValueError:
            raise ValueError(
                f"{run_name!r} has no readable close marker, so this run was "
                f"never closed.  A raw run writes its marker last: no marker "
                f"means the run is still going, was interrupted, or was "
                f"killed, and its tables are not a result") from None
        payloads = {role: _read_run_file(run_fd, filename,
                                         f"the {role} table")
                    for role, filename in names.items() if role != "marker"}
    finally:
        os.close(run_fd)
        os.close(results_fd)

    closed = require_closed_ledger(
        parse_marker(marker_bytes),
        payloads["manifest"], payloads["ledger"], payloads["clocks"],
        payloads["shadow"])
    _refuse_reserved_pilot(closed.manifest)
    return closed
