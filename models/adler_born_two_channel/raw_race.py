#!/usr/bin/env python3
"""
raw_race.py — the finite population of actual competing clocks, and the first
dwell completion that stops the channel.

What changes here
-----------------
Everything before this file walked **one** clock.  ``dynamics`` moved a phase,
``commitment`` decided whether that phase had stayed close enough for long
enough, and ``raw_experiments`` put the two together and returned a record.
There was nothing for the clock to win against, so nothing won.

This file puts ``N`` of them on one exact immutable grid and lets them run
against each other.  The channel stops at the earliest endpoint on which *any*
clock completes its dwell, and that endpoint is the channel's event time.

Read the next paragraph before reading anything else.

**A first clock is a model stopping event.**  It is not a detector click, not
an absorption, not a measurement outcome, not a unique actuality, not a
microscopic bath, not one side of a two-channel competition, and not a
Born-rule result.  The stop is an experimental endpoint for a one-channel
first-passage study: the run is over because the study defined it to be over,
and no physical exclusivity is claimed or implied.  :data:`RACE_NON_CLAIMS`
says the same thing in a form a summary cannot drop silently.

The population is declared, not counted
---------------------------------------
A race takes a :class:`ClockPopulation`: the paths, *and* the declaration they
are supposed to satisfy — how many clocks there are, in what order, with what
mismatches, under which drive, from which data-generating process.

An independent review is why.  Handed a bare tuple, the race could be given the
first path of a declared two-clock grid and would run, return an ordinary
outcome, and claim a category for a population it did not have.  A clock missing
from the *end* leaves a complete range of identifiers — of a smaller population
— so ``range(len(clocks))`` could never have caught it, and the missing clock
might have been the winner.  Completeness is therefore a declared fact the
clocks are checked against, at the race door and again in
:class:`ChannelOutcome`.

No hidden winner
----------------
When several clocks complete on the same endpoint, **all of them are
retained**, with their mismatches and their phases.  Nothing anywhere in this
file picks one of them, and there is no argument through which a caller could
ask it to.  A same-step tie is a fact about the timestep, not a coin to be
tossed: the plan requires co-commit multiplicity to shrink or reach a stable
negligible level under refinement *before* winner-mismatch statistics mean
anything, and a race that silently broke its own ties would have destroyed the
evidence needed to check that.

Why the race is walked twice
----------------------------
Once to find out when the channel stops, and once to record what every clock
was doing up to that moment.

The stop time is not known until every clock has been asked, and the per-clock
diagnostics have to be censored *at* the stop — that is the whole point of
recording them.  A single pass would have to keep every clock's time-resolved
history in order to truncate it afterwards, which is the trial-by-clock-by-time
cube the plan forbids.  So:

* **Pass one** (:func:`first_completion`) walks each clock and keeps one number
  from it, that clock's own completion time, discarding everything else.  Once
  some clock has completed at time ``t``, later clocks are walked only as far
  as ``t`` — an endpoint *at* ``t`` is still processed, so a tie is still seen.
* **Pass two** (:func:`censored_exposure`) walks each clock again, stopping at
  the channel stop, and keeps the diagnostics.

The second pass is where the reproducibility contract earns its keep.  A leaf
belongs to its key — dataset, trial, clock, finest step — and not to the walk
that reached it, so pass two re-derives bit-for-bit the values pass one saw,
and a clock that was cut short in pass one is cut short at exactly the same
place.  If early stopping could move a stream, these two passes would disagree
and the ledger would be built from a different realization than the one that
decided the stop.

The cost is real and is named rather than hidden: the race does up to twice the
clock-steps of a single-pass implementation that did not censor honestly.

Censoring, and the shadow sample
--------------------------------
Per-clock statistics in a first-passage race are *informatively* censored.
Once one clock stops the channel, every other clock's later opportunities are
unobserved — so a mismatch that produced no commitments and a mismatch that
never got the chance look identical in a raw count, and the difference between
them is exactly what a winner-mismatch distribution is supposed to measure.

Two things are recorded against that.  Each clock carries its **exposure**, the
observed time it actually had, and its **eligible time**, the part of that
exposure it spent strictly inside its tongue — the opportunity, rather than the
wall clock.  And :func:`no_stop_exposure` re-walks a small predeclared sample
of trials with the stop switched off, on *the same physical keys*, so the same
clocks on the same realizations can be compared with and without truncation.

The shadow sample is diagnostic.  It never contributes a channel outcome, and
nothing in it enters a category count.

Memory
------
Bounded by the live block, not by the run.  Leaves arrive in windows of finest
steps through ``PhaseNoiseStream.stream_elementary_leaves``, one clock at a
time, and are discarded as they are consumed; clocks are walked in blocks; a
trial holds ``N`` small diagnostic records and a channel row, which the writer
turns into rows and drops.  Nothing here allocates an array indexed by trials,
and nothing allocates an array indexed by both clocks and steps.  Doubling the
trial count or the step count does not change the peak.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional, Tuple

from .commitment import (TRIAL_CATEGORIES, DwellState, _advance_trusted,
                         validate_lock_criterion)
from .dynamics import validate_clock_path
from .raw_experiments import elementary_segments
from .stochastic import PhaseNoiseStream
from .validation import (require_exact_type, require_finite,
                         require_positive_int)

__all__ = [
    "CHANNEL_CATEGORIES",
    "DEFAULT_CLOCK_BLOCK",
    "DEFAULT_STEP_WINDOW",
    "RACE_NON_CLAIMS",
    "ChannelOutcome",
    "ClockExposure",
    "ClockPopulation",
    "PopulationIdentity",
    "censored_exposure",
    "exposure_of",
    "first_completion",
    "no_stop_exposure",
    "race_one_channel",
]

#: The channel-level outcome vocabulary, which is deliberately the *same* four
#: names one clock's exposure already used.
#:
#: ``committed`` — at least one clock completed an uninterrupted dwell before
#: the window ended.  ``never_eligible`` — no clock ever entered its tongue.
#: ``lock_failed`` — at least one clock was eligible and none ever reached the
#: lock band.  ``dwell_failed`` — at least one clock reached the band and none
#: completed the dwell.
#:
#: Mutually exclusive and exhaustive, so every trial lands in exactly one and
#: no trial is conditioned away.  Imported from ``commitment`` rather than
#: restated, because a channel category and a clock category that drifted into
#: two different vocabularies would make the per-clock reconciliation in the
#: ledger meaningless.
CHANNEL_CATEGORIES = TRIAL_CATEGORIES

#: How many clocks a race walks before releasing their working state, and how
#: many finest steps of leaves it holds at once.
#:
#: Both are *execution* parameters and neither is physical.  They bound the
#: live working set and they may not change a single recorded number: the same
#: run at any block size and any window size produces the same ledger, because
#: a leaf belongs to its key rather than to the walk that reached it.  That is
#: asserted rather than asserted-to-be-obvious — see the verifier's
#: batch-invariance check — and it is why these are module constants a caller
#: may override per call, rather than fields of the frozen configuration: a
#: number that cannot affect a result has no business being frozen into the
#: manifest as though it could.
DEFAULT_CLOCK_BLOCK = 16
DEFAULT_STEP_WINDOW = 64

#: What a reader of a channel outcome is not entitled to conclude.  Pinned as
#: text, and pinned again in the README, because these are the sentences a
#: summary drops first and they are the ones carrying the whole meaning of the
#: word "commitment" in this package.
RACE_NON_CLAIMS = (
    "a first clock is a model stopping event: the trial ends because this "
    "one-channel first-passage study defines it to end there",
    "it is not a detector click, an absorption, a measurement outcome or a "
    "Born event",
    "it is not a unique actuality and not global one-world exclusivity; no "
    "physical rule forbids the other clocks from having continued, and the "
    "no-stop shadow sample walks trials in which they do",
    "there is one channel here, so there is no second outcome and no "
    "two-channel competition",
    "the white noise is a minimal effective bath with no microscopic origin",
    "same-endpoint co-completion is retained, never resolved by an internal "
    "random winner; co-commit multiplicity must shrink or reach a stable "
    "negligible level under timestep refinement before winner-mismatch "
    "statistics may be interpreted",
    "per-clock counts are informatively censored by the first winner; read "
    "them against exposure and eligible time, never as raw evidence of "
    "impossibility",
    "no survival curve, no time-dependent commitment law, no fitted coupling "
    "exponent and no scaling claim of any kind is produced here",
)


@dataclass(frozen=True)
class PopulationIdentity:
    """Which population a channel is, stated rather than inferred.

    An independent review took a declared two-clock grid, handed the race the
    first path only, and got back a perfectly ordinary ``ChannelOutcome``: one
    clock, identifier ``[0]``, category ``never_eligible``.  Nothing was wrong
    with any field.  The outcome checked its clock identifiers against
    ``range(len(clocks))``, which catches a clock dropped from the front or the
    middle and *cannot* catch one dropped from the end, because what is left is
    still a complete range — of a smaller population.

    A trailing omission is not a cosmetic error.  The missing clock might have
    been the winner; without it the channel category, the stop time, the
    opportunity and the whole censoring record are those of a different
    experiment.  So completeness stops being inferred from the clocks present
    and becomes a declared fact the clocks are checked against.

    ``model``
        ``"full"`` or ``"width_only"``.  One channel races one
        data-generating process; the width-only control is a separate paired
        run on the same keys, never a competitor inside the same channel.
    ``pulse_digest``
        The drive's own digest.  Two populations with identical detunings under
        different pulses are different experiments, and this is the field that
        says so.
    ``detunings``
        Every clock's fixed mismatch, in clock order.  The clock at index ``i``
        has ``detunings[i]``, and that correspondence is the whole of what a
        clock identifier means: it is half of every key the clock consumes.

    The identity carries the *rule's* detunings — the ones the frozen grid
    produces — not the ones the supplied paths happen to hold.  Comparing the
    two is the point.
    """

    model: str
    pulse_digest: str
    detunings: Tuple[float, ...]

    def __post_init__(self):
        if isinstance(self.model, bool) or not isinstance(self.model, str):
            raise TypeError(f"model must be a string, got "
                            f"{type(self.model).__name__} {self.model!r}")
        if isinstance(self.pulse_digest, bool) \
                or not isinstance(self.pulse_digest, str):
            raise TypeError(f"pulse_digest must be a string, got "
                            f"{type(self.pulse_digest).__name__}")
        if not self.pulse_digest:
            raise ValueError("pulse_digest must not be empty")
        detunings = tuple(require_finite(value, f"detunings[{index}]")
                          for index, value in enumerate(self.detunings))
        if not detunings:
            raise ValueError(
                "a population has at least one clock; an empty declared "
                "population is a run that did not happen, not a race nobody "
                "won")
        object.__setattr__(self, "detunings", detunings)

    @property
    def size(self) -> int:
        """How many clocks this population declares."""
        return len(self.detunings)

    @property
    def clocks(self) -> Tuple[int, ...]:
        """The declared clock identifiers, in order."""
        return tuple(range(len(self.detunings)))

    def require_complete(self, clocks, name: str) -> None:
        """Raise unless ``clocks`` is exactly this declared population.

        ``clocks`` is any sequence of objects carrying ``clock`` and
        ``detuning``.  Length, identifiers, order and mismatches are all
        checked, and the mismatches are compared **exactly**: a same-length
        substitution that moved one clock's detuning would otherwise produce a
        valid-looking outcome for a population that was never declared.
        """
        present = tuple(getattr(item, "clock") for item in clocks)
        if present != self.clocks:
            raise ValueError(
                f"{name} names clocks {list(present)}; the declared "
                f"population is {list(self.clocks)}.  A clock missing from "
                f"the end leaves a complete range of a smaller population, so "
                f"completeness is checked against the declaration rather than "
                f"inferred from what is present")
        for index, item in enumerate(clocks):
            declared = self.detunings[index]
            supplied = getattr(item, "detuning")
            if supplied != declared:
                raise ValueError(
                    f"{name} gives clock {index} the mismatch {supplied!r}; "
                    f"the declared grid gives it {declared!r}.  A clock "
                    f"identifier means the clock at that node of the frozen "
                    f"grid, and a substitution at the same length would "
                    f"otherwise race a population nobody declared")


@dataclass(frozen=True)
class ClockPopulation:
    """The declared population, and the paths that are supposed to be it.

    Two things that a bare tuple of paths cannot be: a *declaration* of how
    many clocks there are and what their mismatches are, and a *checked* claim
    that the supplied paths are those clocks.

    Built from the frozen grid rule — see ``raw_runner.raw_clock_population``
    — so the declaration comes from the same place the manifest's detunings do.
    Every path is rebuilt through ``validate_clock_path`` here, at the door,
    before any key is derived, for the reason ``raw_experiments.run_one_clock``
    gives: a run refused after consuming keys has already changed the
    realization it was refused from.
    """

    paths: Tuple[object, ...]
    identity: PopulationIdentity

    def __post_init__(self):
        require_exact_type(self.identity, PopulationIdentity, "identity")
        rebuilt = tuple(validate_clock_path(path) for path in self.paths)
        if len(rebuilt) != self.identity.size:
            raise ValueError(
                f"the population declares {self.identity.size} clocks and "
                f"{len(rebuilt)} paths were supplied; a race between fewer "
                f"clocks than were declared is a different experiment, "
                f"whichever end the missing ones came from")
        models = sorted({path.model for path in rebuilt})
        if models != [self.identity.model]:
            raise ValueError(
                f"the population declares the {self.identity.model!r} model "
                f"and the supplied paths are {models}; a race is between "
                f"clocks of one data-generating process, and the width-only "
                f"control is a separate paired run on the same keys rather "
                f"than a competitor inside the same channel")
        digests = sorted({path.train.digest for path in rebuilt})
        if digests != [self.identity.pulse_digest]:
            raise ValueError(
                f"the population declares the drive {self.identity.pulse_digest!r} "
                f"and the supplied paths carry {digests}; two populations with "
                f"identical mismatches under different drives are different "
                f"experiments")
        for index, path in enumerate(rebuilt):
            declared = self.identity.detunings[index]
            if path.detuning != declared:
                raise ValueError(
                    f"path {index} has the mismatch {path.detuning!r} and the "
                    f"declared grid gives clock {index} {declared!r}")
        object.__setattr__(self, "paths", rebuilt)

    @property
    def size(self) -> int:
        """How many clocks race in this channel."""
        return len(self.paths)


@dataclass(frozen=True)
class ClockExposure:
    """One clock's censored diagnostics for one trial.

    Everything is measured *as observed*: the walk that produced this stopped
    at the channel stop, so every count here is a count over the opportunity
    the clock actually had.  A zero is therefore never by itself evidence that
    something could not happen.

    ``clock`` / ``detuning``
        Which competitor this is, and its fixed frequency mismatch.  Both are
        carried so a per-clock row is readable without the grid beside it.
    ``category``
        One of :data:`CHANNEL_CATEGORIES`, for this clock alone, as of the
        stop.  A clock is ``committed`` here only if *it* completed a dwell,
        which for every clock but the winners is false in a committed trial.
    ``committed_at``
        The endpoint at which this clock's dwell completed, or ``None``.
    ``first_eligible_at``
        The first observed endpoint at which this clock was strictly inside its
        tongue, or ``None``.  The tongue-entry time of the plan's committed
        record.
    ``first_inside_at``
        The timestamp the first dwell started at, or ``None``.  The band-entry
        time.
    ``last_reset_at``
        The endpoint at which the most recent dwell was cleared, or ``None``.
    ``eligibility_entries``
        How many times the clock *entered* its tongue — transitions from
        ineligible to strictly eligible.  Distinct from
        ``eligible_endpoints``, which counts samples rather than entries: a
        clock that enters once and stays has one entry and many endpoints, and
        a clock the falling edge keeps pushing in and out has many of both.
        Under a pulse a clock normally enters once; more than one is the
        signature of a clock sitting on the tongue boundary.
    ``band_entries`` / ``dwell_resets``
        Fresh dwells started, and started dwells cleared before completing.
        Read together they separate noise-assisted band entry from
        noise-broken dwell, which are different physical statements about the
        same unresolved clock.
    ``endpoints`` / ``eligible_endpoints`` / ``inside_endpoints``
        How many endpoints were observed, how many of them were strictly
        eligible, and how many of them the clock was inside the lock band at.
        The third is the one an exposure-normalized inside-fraction needs, and
        it is not recoverable from the entry and reset counts: those say how
        many dwells there were, not how long they lasted.
    ``exposure``
        Observed time: the channel stop minus the start of the window, or this
        clock's own commitment minus that start, whichever came first.
    ``eligible_time``
        The part of ``exposure`` spent strictly inside the tongue, summed over
        elementary intervals by their endpoint state.  This is the
        *opportunity*, and it is the denominator an exposure-normalized
        per-clock statistic needs.  An endpoint convention, stated because it
        is one: an interval counts as eligible when its own endpoint is, which
        is the same convention every other observable in this package uses.
    ``final_phase``
        The authoritative unwrapped phase at the last observed endpoint.
    """

    clock: int
    detuning: float
    category: str
    committed_at: Optional[float]
    first_eligible_at: Optional[float]
    first_inside_at: Optional[float]
    last_reset_at: Optional[float]
    eligibility_entries: int
    band_entries: int
    dwell_resets: int
    endpoints: int
    eligible_endpoints: int
    inside_endpoints: int
    exposure: float
    eligible_time: float
    final_phase: float

    def __post_init__(self):
        if self.category not in CHANNEL_CATEGORIES:
            raise ValueError(f"category must be one of "
                             f"{list(CHANNEL_CATEGORIES)}, got "
                             f"{self.category!r}")
        for name in ("clock", "eligibility_entries", "band_entries",
                     "dwell_resets", "endpoints", "eligible_endpoints",
                     "inside_endpoints"):
            object.__setattr__(
                self, name,
                require_positive_int(getattr(self, name), name, minimum=0))
        for name in ("detuning", "final_phase"):
            object.__setattr__(self, name,
                               require_finite(getattr(self, name), name))
        for name in ("exposure", "eligible_time"):
            value = require_finite(getattr(self, name), name)
            if value < 0.0:
                raise ValueError(f"{name} must be non-negative, got {value!r}")
            object.__setattr__(self, name, value)
        for name in ("committed_at", "first_eligible_at", "first_inside_at",
                     "last_reset_at"):
            value = getattr(self, name)
            if value is not None:
                object.__setattr__(self, name, require_finite(value, name))

        # Cross-field rules, for the reason ``OneClockRecord`` has them: this
        # object is what a ledger row is built from, and an internally
        # impossible row is worse than a missing one because somebody will
        # later reconcile it.
        if (self.committed_at is None) != (self.category != "committed"):
            raise ValueError(
                f"category {self.category!r} and committed_at "
                f"{self.committed_at!r} disagree about whether this clock "
                f"completed a dwell")
        if self.eligible_time > self.exposure:
            raise ValueError(
                f"eligible_time={self.eligible_time!r} exceeds exposure="
                f"{self.exposure!r}; the time a clock spent inside its tongue "
                f"is part of the time it was observed at all")
        if self.eligible_endpoints > self.endpoints:
            raise ValueError(
                f"eligible_endpoints={self.eligible_endpoints} exceeds the "
                f"{self.endpoints} endpoints it is a subtotal of")
        if self.inside_endpoints > self.eligible_endpoints:
            raise ValueError(
                f"inside_endpoints={self.inside_endpoints} exceeds "
                f"eligible_endpoints={self.eligible_endpoints}; a clock is "
                f"inside the band only at endpoints where it is eligible")
        if self.eligibility_entries > self.eligible_endpoints:
            raise ValueError(
                f"eligibility_entries={self.eligibility_entries} exceeds "
                f"eligible_endpoints={self.eligible_endpoints}; every entry "
                f"lands on an eligible endpoint")
        if (self.eligibility_entries == 0) != (self.eligible_endpoints == 0):
            raise ValueError(
                f"eligibility_entries={self.eligibility_entries} and "
                f"eligible_endpoints={self.eligible_endpoints} disagree about "
                f"whether the clock was ever inside its tongue")
        if self.band_entries > self.inside_endpoints:
            raise ValueError(
                f"band_entries={self.band_entries} exceeds "
                f"inside_endpoints={self.inside_endpoints}; every dwell that "
                f"started was inside the band at the endpoint it started on")
        if self.band_entries > self.eligible_endpoints:
            raise ValueError(
                f"band_entries={self.band_entries} exceeds eligible_endpoints="
                f"{self.eligible_endpoints}; a dwell can only start at a "
                f"strictly eligible endpoint")
        live = self.band_entries - self.dwell_resets
        if live not in (0, 1):
            raise ValueError(
                f"band_entries={self.band_entries} minus dwell_resets="
                f"{self.dwell_resets} is {live}; one clock has at most one "
                f"live dwell")
        if self.category == "committed" and live != 1:
            raise ValueError(
                f"a committed clock ends on the live dwell that completed, so "
                f"band_entries must be exactly one more than dwell_resets; "
                f"got {self.band_entries} and {self.dwell_resets}")
        if self.category == "never_eligible" and self.eligible_endpoints:
            raise ValueError(
                f"category 'never_eligible' with {self.eligible_endpoints} "
                f"strictly eligible endpoints")
        if self.category != "never_eligible" and self.eligible_endpoints < 1:
            raise ValueError(
                f"category {self.category!r} means the clock was eligible "
                f"somewhere, but eligible_endpoints is "
                f"{self.eligible_endpoints}")
        if self.category == "lock_failed" and self.band_entries:
            raise ValueError(
                f"category 'lock_failed' with {self.band_entries} band "
                f"entries; a clock that reached the band and lost it is "
                f"'dwell_failed'")
        if self.category in ("committed", "dwell_failed") \
                and self.band_entries < 1:
            raise ValueError(
                f"category {self.category!r} needs at least one band entry")
        if (self.first_eligible_at is None) != (self.eligible_endpoints == 0):
            raise ValueError(
                f"first_eligible_at={self.first_eligible_at!r} disagrees with "
                f"eligible_endpoints={self.eligible_endpoints}")
        if (self.first_inside_at is None) != (self.band_entries == 0):
            raise ValueError(
                f"first_inside_at={self.first_inside_at!r} disagrees with "
                f"band_entries={self.band_entries}")
        if (self.last_reset_at is None) != (self.dwell_resets == 0):
            raise ValueError(
                f"last_reset_at={self.last_reset_at!r} disagrees with "
                f"dwell_resets={self.dwell_resets}")

    @property
    def committed(self) -> bool:
        """Did this clock complete a dwell within its observed exposure?"""
        return self.committed_at is not None


@dataclass(frozen=True)
class ChannelOutcome:
    """One trial of the one-channel race, complete.

    ``trial``
        The master trial identifier, which is also half of every key the trial
        consumed.
    ``category``
        One of :data:`CHANNEL_CATEGORIES`.  Exactly one, always.
    ``stop_time``
        When observation of this trial ended: the first completion if there was
        one, the end of the pulse window otherwise.
    ``commit_time``
        The first completion, or ``None``.  Equal to ``stop_time`` whenever it
        is not ``None``, and carried separately so that a reader never has to
        infer which kind of stop this was from a comparison of two floats.
    ``censored``
        True when the window was cut short by a commitment.  The channel-level
        statement of the thing every per-clock count in this trial is subject
        to.
    ``winners``
        Every clock that completed its dwell on ``commit_time``, in clock
        order.  Empty for an unresolved trial, and never trimmed to one: see
        this module's docstring.
    ``identity``
        The declared :class:`PopulationIdentity` this channel raced.  Present
        so that completeness is a checked claim rather than an inference from
        the clocks that happen to be here — see that class for the trailing
        omission an independent review got accepted without it.
    ``clocks``
        One :class:`ClockExposure` per clock in the declared population, in
        clock order.
    """

    trial: int
    category: str
    stop_time: float
    commit_time: Optional[float]
    censored: bool
    winners: Tuple[int, ...]
    identity: PopulationIdentity
    clocks: Tuple[ClockExposure, ...]

    def __post_init__(self):
        trial = require_positive_int(self.trial, "trial", minimum=0)
        require_exact_type(self.identity, PopulationIdentity, "identity")
        if self.category not in CHANNEL_CATEGORIES:
            raise ValueError(f"category must be one of "
                             f"{list(CHANNEL_CATEGORIES)}, got "
                             f"{self.category!r}")
        if not isinstance(self.censored, bool):
            raise TypeError(f"censored must be a bool, got "
                            f"{type(self.censored).__name__}")
        stop_time = require_finite(self.stop_time, "stop_time")
        commit_time = self.commit_time
        if commit_time is not None:
            commit_time = require_finite(commit_time, "commit_time")
        clocks = tuple(self.clocks)
        for position, exposure in enumerate(clocks):
            require_exact_type(exposure, ClockExposure, f"clocks[{position}]")
        if not clocks:
            raise ValueError("a channel has at least one clock; an empty "
                             "population is not a race with no winner, it is "
                             "a run that did not happen")
        self.identity.require_complete(clocks, "this channel's clocks")
        winners = tuple(require_positive_int(w, "winners", minimum=0)
                        for w in self.winners)

        if (commit_time is None) != (self.category != "committed"):
            raise ValueError(
                f"category {self.category!r} and commit_time {commit_time!r} "
                f"disagree about whether this channel committed")
        if commit_time is not None and commit_time != stop_time:
            raise ValueError(
                f"commit_time={commit_time!r} is not stop_time={stop_time!r}; "
                f"a committed trial stops at its commitment and nowhere else")
        if self.censored != (commit_time is not None):
            raise ValueError(
                f"censored={self.censored} disagrees with commit_time="
                f"{commit_time!r}; observation is cut short exactly when a "
                f"clock completes")

        committed = tuple(e.clock for e in clocks if e.committed)
        if winners != committed:
            raise ValueError(
                f"the declared winners {winners} are not the clocks whose own "
                f"exposure says they completed, {committed}; a winner list "
                f"that is not derived from the per-clock ledger is a hidden "
                f"selection")
        for clock in committed:
            if clocks[clock].committed_at != commit_time:
                raise ValueError(
                    f"clock {clock} completed at "
                    f"{clocks[clock].committed_at!r} but the channel stopped "
                    f"at {commit_time!r}; every retained co-completer "
                    f"completed on the stopping endpoint itself")
        if winners and self.category != "committed":
            raise ValueError(f"category {self.category!r} with winners "
                             f"{winners}")

        derived = _channel_category(clocks)
        if derived != self.category:
            raise ValueError(
                f"category {self.category!r} contradicts the per-clock "
                f"ledger, which says {derived!r}; the channel category is "
                f"derived from the clocks and cannot be set beside them")

        object.__setattr__(self, "trial", trial)
        object.__setattr__(self, "stop_time", stop_time)
        object.__setattr__(self, "commit_time", commit_time)
        object.__setattr__(self, "winners", winners)
        object.__setattr__(self, "clocks", clocks)

    @property
    def committed(self) -> bool:
        """Did any clock complete a dwell before the window ended?"""
        return self.commit_time is not None

    @property
    def co_completion(self) -> int:
        """How many clocks completed on the stopping endpoint.

        Zero for an unresolved trial and one for the ordinary committed case.
        Two or more is a same-endpoint tie: a real, retained, reportable state
        of this discretization, and the quantity whose behaviour under timestep
        refinement gates any winner-mismatch reading.
        """
        return len(self.winners)

    @property
    def clocks_ever_eligible(self) -> int:
        """How many clocks entered their tongue within the observed window."""
        return sum(1 for e in self.clocks if e.eligible_endpoints > 0)

    @property
    def clocks_ever_inside(self) -> int:
        """How many clocks reached the lock band within the observed window."""
        return sum(1 for e in self.clocks if e.band_entries > 0)

    @property
    def band_entries(self) -> int:
        """Fresh dwells started by any clock, over the observed window."""
        return sum(e.band_entries for e in self.clocks)

    @property
    def dwell_resets(self) -> int:
        """Started dwells cleared by any clock, over the observed window."""
        return sum(e.dwell_resets for e in self.clocks)

    @property
    def eligible_time(self) -> float:
        """Total opportunity: eligible time summed over the population."""
        return math.fsum(e.eligible_time for e in self.clocks)


def _channel_category(clocks) -> str:
    """The one category this population's exposures put the channel in.

    Written once, used by the constructor and by the race, so that the value
    recorded and the value validated cannot come from two different rules.
    """
    if any(e.committed for e in clocks):
        return "committed"
    if not any(e.eligible_endpoints for e in clocks):
        return "never_eligible"
    if not any(e.band_entries for e in clocks):
        return "lock_failed"
    return "dwell_failed"


def _windows(start_index: int, steps: int, window: int):
    """Yield ``(first_step, count)`` for each live block of finest steps.

    The same partition ``PhaseNoiseStream._stream_leaves`` uses, so the leaves
    of one window and the elementary intervals of that window line up without
    either side knowing about the other.
    """
    emitted = 0
    while emitted < steps:
        take = min(window, steps - emitted)
        yield start_index + emitted, take
        emitted += take


@dataclass(frozen=True)
class _Walk:
    """What walking one clock produced.  Private: it is not a ledger row."""

    dwell: DwellState
    first_eligible_at: Optional[float]
    first_inside_at: Optional[float]
    last_reset_at: Optional[float]
    eligibility_entries: int
    endpoints: int
    eligible_endpoints: int
    inside_endpoints: int
    eligible_time: float
    final_phase: float
    observed_until: float


def _walk_clock(path, criterion, stream, trial, clock, phase, start_index,
                steps, window, stop_time) -> _Walk:
    """Walk one clock over the elementary mesh; stop at ``stop_time``.

    ``stop_time`` is ``None`` for an unbounded walk — to the end of the window
    or to this clock's own commitment, whichever comes first — and otherwise a
    time the walk will not observe past.  An endpoint exactly *at* ``stop_time``
    is observed, which is what makes a same-endpoint tie visible: a clock
    completing on the stopping endpoint is a co-completer, not a latecomer.

    That comparison is exact rather than tolerant, and can be: ``stop_time`` is
    always an endpoint some other clock's walk produced on this same mesh, and
    the elementary mesh is derived from ``FinestMesh.step_bounds`` and the
    frozen crossing times identically for every clock.  Two clocks therefore
    have bit-identical endpoint sequences wherever their crossing schedules
    agree, and the stopping endpoint is a member of both.

    Everything is streamed.  Leaves arrive one window of finest steps at a
    time, the matching elementary intervals are derived for that same window,
    and both are dropped when the window is done.
    """
    state = path._initial_state_trusted(stream.mesh.time_at(start_index), phase)
    crossings = path._cached_schedule.crossings(clock, path.train.digest)

    dwell = DwellState()
    first_eligible_at = None
    first_inside_at = None
    last_reset_at = None
    eligibility_entries = 0
    endpoints = 0
    eligible_endpoints = 0
    inside_endpoints = 0
    eligible_time = 0.0
    observed_until = state.time
    was_eligible = False

    def finished():
        """The walk so far, as a record.  One construction, three exits."""
        return _Walk(dwell, first_eligible_at, first_inside_at, last_reset_at,
                     eligibility_entries, endpoints, eligible_endpoints,
                     inside_endpoints, eligible_time, state.phase,
                     observed_until)

    blocks = stream.stream_elementary_leaves(trial, clock, start_index, steps,
                                             crossings, window)
    for (first_step, take), leaves in zip(_windows(start_index, steps, window),
                                          blocks):
        segments = elementary_segments(stream.mesh, first_step, take, crossings)
        if len(leaves) != len(segments):
            raise ValueError(
                f"the noise tree produced {len(leaves)} elementary leaves for "
                f"finest steps [{first_step}, {first_step + take}) and the "
                f"schedule {len(segments)} elementary intervals; the phase and "
                f"the noise are not being walked on one mesh")
        for leaf, (opening, closing) in zip(leaves, segments):
            if stop_time is not None and closing > stop_time:
                return finished()
            state = path._advance_trusted(state, leaf, closing)
            if closing == opening:
                continue                  # zero-duration handoff: no endpoint
            seen = path._observe_trusted(state)
            endpoints += 1
            observed_until = closing
            if seen.eligible:
                eligible_endpoints += 1
                eligible_time += closing - opening
                if not was_eligible:
                    eligibility_entries += 1
                if first_eligible_at is None:
                    first_eligible_at = closing
            was_eligible = seen.eligible
            before = dwell
            dwell = _advance_trusted(criterion, dwell, seen.time, seen.eligible,
                                     seen.lock_error, seen.contraction)
            if dwell.inside_since is not None:
                inside_endpoints += 1
            if dwell.band_entries > before.band_entries:
                if first_inside_at is None:
                    first_inside_at = dwell.inside_since
            if dwell.dwell_resets > before.dwell_resets:
                last_reset_at = closing
            if dwell.committed:
                return finished()
    return finished()


def _validated_race(population, criterion, stream, trial, start_index, steps,
                    clock_block, step_window):
    """Rebuild and re-check every frozen input; return the trusted tuple.

    The same trusted-snapshot rule ``raw_experiments.run_one_clock`` follows,
    at the population's door: the population and the criterion are rebuilt from
    their own fields rather than accepted on the strength of their type, and
    the rebuilt copies are what the race uses.  It happens here, before
    anything else, so that a refused run has not yet derived a key and cannot
    have changed the realization it was refused from.

    The first argument is a :class:`ClockPopulation` rather than a tuple of
    paths, and the change is the whole of an independent review's third
    finding.  A tuple cannot say how many clocks there were *supposed* to be,
    so a race handed the first path of a declared two-clock grid ran, returned
    a valid outcome, and claimed a category for a population it did not have.
    A population declares its size, its clock order and its mismatches, and
    rebuilding it here is what checks the paths against that declaration.
    """
    require_exact_type(population, ClockPopulation, "population")
    rebuilt = ClockPopulation(population.paths, population.identity)
    criterion = validate_lock_criterion(criterion)
    require_exact_type(stream, PhaseNoiseStream, "stream")
    trial = require_positive_int(trial, "trial", minimum=0)
    start_index = require_positive_int(start_index, "start_index", minimum=0)
    steps = require_positive_int(steps, "steps", minimum=1)
    clock_block = require_positive_int(clock_block, "clock_block", minimum=1)
    step_window = require_positive_int(step_window, "step_window", minimum=1)
    return rebuilt, criterion, stream, trial, start_index, steps, \
        clock_block, step_window


def first_completion(population, criterion, stream, trial, start_index, steps,
                     clock_block=DEFAULT_CLOCK_BLOCK,
                     step_window=DEFAULT_STEP_WINDOW) -> Optional[float]:
    """Return the earliest endpoint on which any clock completes its dwell.

    ``None`` when no clock completes within the window, which is the
    unresolved trial and is a result rather than a missing one.

    Pass one of the race.  One number comes out and everything else is
    discarded, because everything else has to be measured against the answer
    this pass produces and cannot be recorded before it exists.

    The population is the declared one: :func:`first_completion` takes a
    :class:`ClockPopulation` and rebuilds it, so a truncated or substituted
    set of paths is refused here rather than producing a stop for a channel
    nobody declared.

    Clocks are walked in blocks and, once some clock has completed at ``t``,
    every later clock is walked only as far as ``t``.  That truncation is an
    optimization and is required to be nothing more: the returned time does not
    depend on the block size, on the window size, or on the order the
    population is walked in, because a clock's realization is a function of its
    key alone.  The verifier asserts that directly rather than leaving it to
    this paragraph.
    """
    return _first_completion_trusted(*_validated_race(
        population, criterion, stream, trial, start_index, steps, clock_block,
        step_window))


def _first_completion_trusted(population, criterion, stream, trial,
                              start_index, steps, clock_block, step_window):
    """Unvalidated kernel behind :func:`first_completion`."""
    paths = population.paths
    earliest = None
    for block in range(0, len(paths), clock_block):
        for clock in range(block, min(block + clock_block, len(paths))):
            phase = float(stream.initial_phases(trial, [clock])[0])
            walk = _walk_clock(paths[clock], criterion, stream, trial, clock,
                               phase, start_index, steps, step_window,
                               earliest)
            moment = walk.dwell.committed_at
            if moment is not None and (earliest is None or moment < earliest):
                earliest = moment
    return earliest


def censored_exposure(population, criterion, stream, trial, start_index, steps,
                      stop_time, clock_block=DEFAULT_CLOCK_BLOCK,
                      step_window=DEFAULT_STEP_WINDOW
                      ) -> Tuple[ClockExposure, ...]:
    """Return every clock's diagnostics, observed up to ``stop_time``.

    Pass two of the race.  ``stop_time`` is the channel stop that pass one
    established — a commitment, or the end of the window — and no clock is
    observed past it.  A clock whose own dwell completes exactly there is
    recorded as having completed; that is the co-completion case, and there is
    no tie-break anywhere in this function.
    """
    checked = _validated_race(population, criterion, stream, trial,
                              start_index, steps, clock_block, step_window)
    stop_time = require_finite(stop_time, "stop_time")
    if stop_time < stream.mesh.time_at(start_index):
        raise ValueError(
            f"stop_time={stop_time!r} precedes the start of the window "
            f"{stream.mesh.time_at(start_index)!r}; a channel cannot stop "
            f"before it opens")
    return _exposure_trusted(*checked, stop_time=stop_time, no_stop=False)


def no_stop_exposure(population, criterion, stream, trial, start_index, steps,
                     clock_block=DEFAULT_CLOCK_BLOCK,
                     step_window=DEFAULT_STEP_WINDOW
                     ) -> Tuple[ClockExposure, ...]:
    """Return every clock's diagnostics with the first-winner stop switched off.

    The plan's predeclared **no-stop shadow sample**, for one trial.  Every
    clock is walked to the end of the pulse window — or to its own commitment,
    which is where its own dwell machine stops — regardless of what any other
    clock did.

    On the same keys.  Same dataset namespace, same trial, same clock
    identifiers, same mesh, therefore the same prepared phases and the same
    Brownian leaves as the censored run of the same trial.  That is what makes
    the two comparable at all: a shadow sample drawn from a fresh realization
    would differ from the censored run for two reasons at once and could not
    separate them.

    Diagnostic only.  Nothing here is a channel outcome, nothing here enters a
    category count, and a commitment recorded here is a commitment that in the
    actual race *did not happen*, because the channel had already stopped.
    """
    checked = _validated_race(population, criterion, stream, trial,
                              start_index, steps, clock_block, step_window)
    closed = stream.mesh.time_at(start_index + steps)
    return _exposure_trusted(*checked, stop_time=closed, no_stop=True)


def _exposure_trusted(population, criterion, stream, trial, start_index, steps,
                      clock_block, step_window, stop_time, no_stop):
    """Unvalidated kernel behind the two exposure walks.

    One body, two doors, and the only difference between them is whether the
    walk is allowed to run past ``stop_time``.  Written once so that the
    censored run and the shadow run cannot drift into two different definitions
    of a band entry, a reset, or an eligible interval — which would make the
    comparison between them measure the implementations rather than the
    censoring.
    """
    paths = population.paths
    opened = stream.mesh.time_at(start_index)
    bound = None if no_stop else stop_time
    exposures = []
    for block in range(0, len(paths), clock_block):
        for clock in range(block, min(block + clock_block, len(paths))):
            phase = float(stream.initial_phases(trial, [clock])[0])
            walk = _walk_clock(paths[clock], criterion, stream, trial, clock,
                               phase, start_index, steps, step_window, bound)
            exposures.append(_exposure_of(paths[clock], clock, walk, opened,
                                          stop_time))
    return tuple(exposures)


def _exposure_of(path, clock, walk, opened, horizon) -> ClockExposure:
    """Turn one walk into one clock's ledger row.

    ``horizon`` is when observation of this clock ended for reasons outside
    itself — the channel stop, or the end of the window in a no-stop walk.  A
    clock that committed stopped earlier, at its own commitment, and its
    exposure is that shorter interval: it was not observed through a time it
    had already ended.
    """
    dwell = walk.dwell
    ended = dwell.committed_at if dwell.committed_at is not None else horizon
    return ClockExposure(
        clock=clock,
        detuning=path.detuning,
        category=dwell.category,
        committed_at=dwell.committed_at,
        first_eligible_at=walk.first_eligible_at,
        first_inside_at=walk.first_inside_at,
        last_reset_at=walk.last_reset_at,
        eligibility_entries=walk.eligibility_entries,
        band_entries=dwell.band_entries,
        dwell_resets=dwell.dwell_resets,
        endpoints=walk.endpoints,
        eligible_endpoints=walk.eligible_endpoints,
        inside_endpoints=walk.inside_endpoints,
        exposure=exposure_of(opened, ended),
        eligible_time=walk.eligible_time,
        final_phase=walk.final_phase)


def exposure_of(window_open, ended) -> float:
    """Return the observed time between ``window_open`` and ``ended``.

    Three lines with a name, because a ledger consumer has to be able to
    recompute this exact number from the manifest window and the channel stop
    and get the identical binary64 back.  An independent review forged a run
    reporting 999 units of exposure inside a 0.4-unit window and it was
    accepted, because nothing on the reading side derived the quantity — it
    only checked that it was a non-negative float.  One expression, used by the
    writer and by the gate, is what makes that forgery a mismatch rather than a
    plausible number.
    """
    opened = require_finite(window_open, "window_open")
    closed = require_finite(ended, "ended")
    return max(closed - opened, 0.0)


def race_one_channel(population, criterion, stream, trial, start_index, steps,
                     clock_block=DEFAULT_CLOCK_BLOCK,
                     step_window=DEFAULT_STEP_WINDOW) -> ChannelOutcome:
    """Run one trial of the finite one-channel race.

    Both passes, and the outcome they agree on.  The channel stops at the
    earliest dwell completion across the whole population, every clock that
    completed on that endpoint is retained, and every clock — winner, loser and
    clock whose tongue never opened — carries the diagnostics of the exposure
    it actually had.

    The trial lands in exactly one of :data:`CHANNEL_CATEGORIES`, derived from
    the per-clock ledger rather than recorded beside it, and no trial is
    conditioned away.
    """
    checked = _validated_race(population, criterion, stream, trial,
                              start_index, steps, clock_block, step_window)
    population, _, stream, trial, start_index, steps, _, _ = checked
    commit_time = _first_completion_trusted(*checked)
    closed = stream.mesh.time_at(start_index + steps)
    stop_time = closed if commit_time is None else commit_time
    clocks = _exposure_trusted(*checked, stop_time=stop_time, no_stop=False)
    winners = tuple(e.clock for e in clocks if e.committed)
    return ChannelOutcome(trial=trial,
                          category=_channel_category(clocks),
                          stop_time=stop_time,
                          commit_time=commit_time,
                          censored=commit_time is not None,
                          winners=winners,
                          identity=population.identity,
                          clocks=clocks)
