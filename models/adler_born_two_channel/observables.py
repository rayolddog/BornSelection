#!/usr/bin/env python3
"""
observables.py — survival, risk sets and hazard, derived from a closed ledger.

This is the first module of the **read-only analysis layer**.  It sits outside
the raw event graph on purpose: ``verify.py`` names ``observables``,
``analysis`` and ``compare`` in the raw graph's forbidden list, walks the import
graph from ``raw_runner.py``, and starts a fresh interpreter to confirm that
importing the raw entry point loads none of the three.  A module that the
generator can reach is a module the generator can be influenced by, and the
whole point of a survival table is that it is computed *after* the events, from
bytes that are already closed and hashed.

What this module does
---------------------
It turns one :class:`~.raw_ledger.ClosedLedger` into exact counts:

* how many trials are at risk at the start of each interval;
* how many committed inside it and how many were censored out of it;
* how wide the interval is, in the run's own time units;
* the discrete conditional failure probability, ``events / at risk``;
* the time-rate hazard, which is that probability **divided by the width**;
* the product-limit survival and the Nelson-Aalen cumulative hazard;
* the unresolved fraction with its three separate causes; and
* which clocks supplied the first commitment, including every co-completer.

Nothing here fits anything.  There is no exponent, no coupling sweep and no
predictor in this file; those are ``analysis.py`` and ``compare.py``.

Three conventions, frozen
-------------------------
An observable is only comparable across runs if the rules that produced it are
the same rules, so the three that decide what a number means are constants in
this file rather than choices made at each call site:

``EVENT_CONVENTION``
    Bins are half-open on the left: an event at ``t`` belongs to the one
    interval with ``left < t <= right``.  The first edge is the window's
    opening, at which every trial is at risk and nothing has happened yet.

``CENSOR_CONVENTION``
    A trial that never commits is *administratively right-censored* at the
    manifest's window close.  It is never dropped, never conditioned away, and
    never counted as an event.

``TIE_CONVENTION``
    Several clocks completing on one endpoint is **one** channel event.  All of
    them are retained as co-completers, none is chosen as *the* winner, and the
    winner-mismatch counts therefore sum to the number of *co-completers*, not
    to the number of committed trials.  Both totals are reported side by side
    so the difference cannot be read as an error.

The word "hazard"
-----------------
``events / at risk`` is a **discrete conditional failure probability**.  It is
dimensionless, it depends on how wide the interval happens to be, and it is not
a rate.  Dividing it by the interval width gives something with units of
inverse time, and only that may be compared against a quantity that also has
those units — a relaxation rate, a flux, a constant hazard.

This distinction is not pedantry: the same events binned twice as coarsely
halve the rate and leave the conditional probability roughly doubled, so a
comparison that used the wrong one would be reporting the binning.  Every row
this module produces therefore carries the risk count, the event count, the
censoring count *and* the width, and :data:`DISPLAY_COLUMNS` is the complete
list, so a table cannot be printed with the width left out.

One thing the ledger's own ``censored`` column does not mean
------------------------------------------------------------
The trial ledger has a ``censored`` boolean, and it is ``True`` exactly when
the trial **committed**.  That is a statement about the *clocks*: a first
commitment stops the channel, so every other clock's observation is cut short
there.  It is not the survival-analysis censoring indicator, which is the
opposite — the trials censored *for the channel-level event* are the ones that
never committed.  This module never reads that column, and computes the
censoring it needs from the category.

Non-claims
----------
A survival curve here is the survival of a *model stopping event*.  It is not a
detector click, an absorption, a measurement outcome, a unique actuality, one
side of a two-channel competition, or a Born-rule result.  A hazard that rises
and falls with a pulse is a fact about this dwell machine under this noise law,
not evidence about microscopic physics.  See :data:`OBSERVABLE_NON_CLAIMS`.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .raw_ledger import (ClosedLedger, TRIAL_CATEGORIES,
                         require_unmutated_ledger)
from .validation import (require_finite, require_finite_array,
                         require_instance, require_positive_int)

__all__ = [
    "CENSOR_CONVENTION",
    "DISPLAY_COLUMNS",
    "EVENT_CONVENTION",
    "OBSERVABLE_NON_CLAIMS",
    "TIE_CONVENTION",
    "UNRESOLVED_CATEGORIES",
    "CommitmentCounts",
    "RiskInterval",
    "RiskTable",
    "WinnerMismatch",
    "commitment_counts",
    "endpoint_cumulative_hazard",
    "risk_table",
    "uniform_edges",
    "winner_mismatch",
]

#: How an instant is assigned to an interval.  Stated as data so that a reader
#: of a table does not have to reconstruct it from the code that made one.
EVENT_CONVENTION = (
    "intervals are half-open on the left: an event at t belongs to the one "
    "interval with left < t <= right.  The first edge is the manifest's "
    "window_open, at which every trial is at risk and nothing has happened.")

#: What happens to a trial that never commits.
CENSOR_CONVENTION = (
    "a trial in any of the three unresolved categories is administratively "
    "right-censored at the manifest's window_close.  It is never discarded, "
    "never conditioned away, and never counted as an event.")

#: What happens when several clocks complete on one endpoint.
TIE_CONVENTION = (
    "several clocks completing on one endpoint is ONE channel event.  Every "
    "co-completer is retained and none is chosen as the winner, so "
    "winner-mismatch counts sum to the number of co-completers rather than to "
    "the number of committed trials.  Both totals are reported.")

#: The three unresolved causes, in the order the ledger declares them.
UNRESOLVED_CATEGORIES = tuple(c for c in TRIAL_CATEGORIES if c != "committed")

#: Every field a binned display row carries, in order.
#:
#: This is the anti-mislabelling rule in one place.  ``conditional_failure`` is
#: ``events / at_risk`` and is dimensionless; ``hazard_rate`` is that divided by
#: ``width`` and has units of inverse time.  A row that carried one of them and
#: not the risk count and the width would let the two be confused, so the list
#: is complete and :meth:`RiskTable.rows` emits all of it.
DISPLAY_COLUMNS = ("left", "right", "width", "at_risk", "events", "censored",
                   "exits", "conditional_failure", "hazard_rate",
                   "survival", "cumulative_hazard")

#: What a curve produced here does not license anybody to say.
OBSERVABLE_NON_CLAIMS = (
    "a survival curve here is the survival of a MODEL STOPPING EVENT: the "
    "trial ends because this one-channel first-passage study defines it to "
    "end at the first completed dwell",
    "it is not a detector click, an absorption, a measurement outcome, a "
    "unique actuality, one side of a two-channel competition, or a Born-rule "
    "result",
    "events divided by risk is a discrete conditional failure probability, "
    "not a rate; only the width-divided hazard has units of inverse time",
    "a hazard that rises and falls with the pulse is a fact about this dwell "
    "machine under this noise law, not evidence about microscopic physics",
    "no coupling exponent, no scaling law and no comparison with any "
    "prediction is computed in this module",
)


def _require_int_count(value, name: str) -> int:
    """Return ``value`` as a non-negative plain ``int``, or raise."""
    if isinstance(value, (bool, np.bool_)):
        raise TypeError(f"{name} must be an integer count, got "
                        f"{type(value).__name__} {value!r}")
    if not isinstance(value, (int, np.integer)):
        raise TypeError(f"{name} must be an integer count, got "
                        f"{type(value).__name__} {value!r}")
    out = int(value)
    if out < 0:
        raise ValueError(f"{name} must be non-negative, got {out!r}")
    return out


@dataclass(frozen=True)
class RiskInterval:
    """One row of a risk table: what happened between two edges.

    ``left``, ``right``
        The interval's endpoints, in the run's own time units.  ``right`` must
        strictly follow ``left``; a zero-width interval would make the time
        rate below infinite for any event at all.
    ``at_risk``
        Trials still unresolved at ``left``.  Not "trials in the interval": a
        trial censored in an earlier interval is not at risk here.
    ``events``
        Channel commitments with ``left < t <= right``.  A same-endpoint tie is
        one event; see :data:`TIE_CONVENTION`.
    ``censored``
        Trials whose observation ended in this interval without a commitment.

    ``events + censored <= at_risk`` is enforced: more exits than entrants is
    not a small error, it is a different table.

    The two derived quantities are deliberately separate attributes with
    different units, and both are always present:

    * :attr:`conditional_failure` is ``events / at_risk`` — dimensionless, and
      the honest thing to report for a discrete-time table;
    * :attr:`hazard_rate` is that divided by :attr:`width` — inverse time, and
      the only one of the two that may be compared against a relaxation rate,
      a flux or a constant hazard.

    Both are NaN when nothing is at risk, because zero events out of zero
    trials is not a probability of zero.
    """

    left: float
    right: float
    at_risk: int
    events: int
    censored: int

    def __post_init__(self):
        left = require_finite(self.left, "left")
        right = require_finite(self.right, "right")
        if not right > left:
            raise ValueError(
                f"right={right!r} must strictly follow left={left!r}; a "
                f"zero-width or reversed interval has no time rate, and a "
                f"width of zero would make one infinite for any event at all")
        at_risk = _require_int_count(self.at_risk, "at_risk")
        events = _require_int_count(self.events, "events")
        censored = _require_int_count(self.censored, "censored")
        if events + censored > at_risk:
            raise ValueError(
                f"{events} events and {censored} censorings leave the "
                f"interval, but only {at_risk} trials were at risk entering "
                f"it; more exits than entrants is a different table, not a "
                f"rounding error")
        object.__setattr__(self, "left", left)
        object.__setattr__(self, "right", right)
        object.__setattr__(self, "at_risk", at_risk)
        object.__setattr__(self, "events", events)
        object.__setattr__(self, "censored", censored)

    @property
    def width(self) -> float:
        """``right - left``, the interval's duration."""
        return self.right - self.left

    @property
    def exits(self) -> int:
        """Trials that left the risk set here, for any reason."""
        return self.events + self.censored

    @property
    def conditional_failure(self) -> float:
        """``events / at_risk`` — a probability, **not** a rate.

        NaN when nothing was at risk: zero events out of zero trials is not a
        probability of zero, and reporting one would put a fabricated point on
        every plot of a curve that had run out of survivors.
        """
        if self.at_risk == 0:
            return float("nan")
        return self.events / self.at_risk

    @property
    def hazard_rate(self) -> float:
        """The conditional failure probability **divided by the width**.

        This is the only quantity in this class with units of inverse time, and
        the division is written here rather than left to a caller precisely so
        that it cannot be forgotten.  Comparing ``events / at_risk`` against a
        relaxation rate compares a probability with a rate, and the answer
        changes if the binning does.
        """
        return self.conditional_failure / self.width


@dataclass(frozen=True)
class RiskTable:
    """A complete risk table over one closed run's observation window.

    ``intervals``
        The rows, in time order and edge-to-edge: each interval's ``left`` is
        the previous one's ``right``, with no gap and no overlap.  A table with
        a gap would be a table over a *chosen* sub-window, which is the
        post-hoc selection this layer exists not to do.
    ``trials``
        The run's trial count.  Every one of them enters at ``entry_time``.
    ``entry_time``
        The window's opening: the first interval's left edge.

    Four invariants are checked at construction, and each of them is a way a
    plausible-looking table can be wrong:

    1. the intervals tile the window contiguously;
    2. the risk set is exactly ``trials`` at entry and decreases by each
       interval's exits — so a risk count cannot be asserted beside the counts
       that determine it;
    3. every count is non-negative and no interval loses more trials than it
       had; and
    4. the totals conserve: ``events + censored`` over the whole table equals
       ``trials``, so nothing has been dropped and nothing invented.

    Survival is the product-limit estimator and is non-increasing by
    construction; :meth:`survival` returns it with the risk counts beside it,
    for the same reason every row carries its width.
    """

    intervals: tuple
    trials: int
    entry_time: float

    def __post_init__(self):
        intervals = tuple(self.intervals)
        if not intervals:
            raise ValueError("a risk table has at least one interval; an "
                             "empty one describes no window at all")
        for position, interval in enumerate(intervals):
            require_instance(interval, RiskInterval, f"intervals[{position}]")
        trials = require_positive_int(self.trials, "trials", minimum=1)
        entry = require_finite(self.entry_time, "entry_time")

        if intervals[0].left != entry:
            raise ValueError(
                f"the first interval opens at {intervals[0].left!r} and the "
                f"table declares entry at {entry!r}; a table that starts after "
                f"its own entry is a table over a chosen sub-window")
        for position in range(1, len(intervals)):
            previous, current = intervals[position - 1], intervals[position]
            if current.left != previous.right:
                raise ValueError(
                    f"interval {position} opens at {current.left!r} and "
                    f"interval {position - 1} closed at {previous.right!r}; "
                    f"the rows must tile the window with no gap and no "
                    f"overlap, or the table is over a sub-window")

        remaining = trials
        for position, interval in enumerate(intervals):
            if interval.at_risk != remaining:
                raise ValueError(
                    f"interval {position} declares {interval.at_risk} at risk; "
                    f"{remaining} trials had not yet left.  The risk set is "
                    f"derived from the exits before it and cannot be asserted "
                    f"beside them")
            remaining -= interval.exits
        if remaining != 0:
            raise ValueError(
                f"{remaining} of {trials} trials never leave this table; every "
                f"trial either commits or is censored at the window close, so "
                f"a survivor at the end means an event or a censoring was "
                f"dropped")
        object.__setattr__(self, "intervals", intervals)
        object.__setattr__(self, "trials", trials)
        object.__setattr__(self, "entry_time", entry)

    @property
    def events(self) -> int:
        """Total channel commitments in the table."""
        return sum(interval.events for interval in self.intervals)

    @property
    def censored(self) -> int:
        """Total administratively censored trials."""
        return sum(interval.censored for interval in self.intervals)

    def survival(self):
        """Return the product-limit survival, one entry per interval edge.

        A tuple of ``(time, at_risk, survival)`` triples: the entry time with
        the full risk set and survival one, then one entry per interval right
        edge.  The risk count travels with the estimate because a survival of
        ``0.2`` computed from five survivors and from five thousand are
        different statements, and a plot of the value alone cannot tell them
        apart.

        Non-increasing by construction — each factor is ``1 - events/at_risk``,
        which lies in ``[0, 1]`` — and the constructor has already established
        that the risk counts are the ones the exits imply.
        """
        out = [(self.entry_time, self.trials, 1.0)]
        value = 1.0
        for interval in self.intervals:
            if interval.at_risk:
                value *= 1.0 - interval.events / interval.at_risk
            out.append((interval.right, interval.at_risk - interval.exits,
                        value))
        return tuple(out)

    def cumulative_hazard(self):
        """Return the Nelson-Aalen cumulative hazard at each interval edge.

        A tuple of ``(time, cumulative_hazard)`` pairs beginning at zero.  This
        is the *dimensionless* accumulation ``sum events/at_risk``, not an
        integral of the time rate; the two agree only in the limit of equal
        narrow bins, and the difference is exactly the binning artefact this
        module keeps visible.
        """
        out = [(self.entry_time, 0.0)]
        total = 0.0
        for interval in self.intervals:
            if interval.at_risk:
                total += interval.events / interval.at_risk
            out.append((interval.right, total))
        return tuple(out)

    def rows(self):
        """Return one dictionary per interval, carrying every display field.

        The keys are exactly :data:`DISPLAY_COLUMNS`, always, so a table cannot
        be printed with the risk count or the width left out — which is the one
        omission that lets ``events / at_risk`` be read as a rate.
        """
        survival = self.survival()[1:]
        hazard = self.cumulative_hazard()[1:]
        out = []
        for interval, (_, _, s), (_, h) in zip(self.intervals, survival,
                                               hazard):
            out.append({
                "left": interval.left,
                "right": interval.right,
                "width": interval.width,
                "at_risk": interval.at_risk,
                "events": interval.events,
                "censored": interval.censored,
                "exits": interval.exits,
                "conditional_failure": interval.conditional_failure,
                "hazard_rate": interval.hazard_rate,
                "survival": s,
                "cumulative_hazard": h,
            })
        return tuple(out)


@dataclass(frozen=True)
class CommitmentCounts:
    """The channel-level tally of one closed run, by category.

    The three unresolved causes are kept apart because they mean different
    things: ``never_eligible`` is a tongue that never opened wide enough,
    ``lock_failed`` is a tongue entered and a band never reached, and
    ``dwell_failed`` is a band reached and a dwell never completed.  Pooling
    them into "did not commit" throws away the only diagnostic that separates a
    detector that never saw anything from one that saw and did not latch.
    """

    trials: int
    committed: int
    never_eligible: int
    lock_failed: int
    dwell_failed: int

    def __post_init__(self):
        trials = require_positive_int(self.trials, "trials", minimum=1)
        counts = {}
        for name in ("committed",) + UNRESOLVED_CATEGORIES:
            counts[name] = _require_int_count(getattr(self, name), name)
        total = sum(counts.values())
        if total != trials:
            raise ValueError(
                f"the four categories hold {total} trials and the run has "
                f"{trials}; they are mutually exclusive and exhaustive, so a "
                f"disagreement means a trial was counted twice or not at all")
        object.__setattr__(self, "trials", trials)
        for name, value in counts.items():
            object.__setattr__(self, name, value)

    @property
    def unresolved(self) -> int:
        """Trials in any of the three unresolved categories."""
        return self.never_eligible + self.lock_failed + self.dwell_failed

    @property
    def commitment_probability(self) -> float:
        """``committed / trials`` — the raw endpoint proportion.

        Reported, and never fed to a logarithm.  ``analysis.py`` fits the
        counts themselves through a complementary-log-log likelihood precisely
        so that a zero or a one here stays in the estimate instead of being
        transformed into an infinity and quietly dropped.
        """
        return self.committed / self.trials


@dataclass(frozen=True)
class WinnerMismatch:
    """Which clocks supplied first commitments, and how often they shared.

    ``clocks``, ``detunings``
        The full declared population, in population order — every clock,
        including those that never won.  A table over only the winners would
        make an untouched part of the tongue invisible.
    ``wins``
        Per clock, the number of trials in which it completed the dwell that
        stopped the channel.
    ``committed_trials``
        Trials with at least one winner.
    ``co_completers``
        Total winners summed over trials, which is ``sum(wins)``.  It exceeds
        ``committed_trials`` exactly when ties occurred, and the two are
        reported together so that the difference reads as the tie count it is
        rather than as an inconsistency.
    ``multiplicity``
        ``{k: trials whose winner count was k}`` for every ``k >= 1``, which is
        the distribution the plan's refinement gate on co-commit multiplicity
        is about.
    """

    clocks: tuple
    detunings: tuple
    wins: tuple
    committed_trials: int
    multiplicity: tuple

    def __post_init__(self):
        clocks = tuple(_require_int_count(c, f"clocks[{i}]")
                       for i, c in enumerate(self.clocks))
        detunings = tuple(require_finite(d, f"detunings[{i}]")
                          for i, d in enumerate(self.detunings))
        wins = tuple(_require_int_count(w, f"wins[{i}]")
                     for i, w in enumerate(self.wins))
        if not (len(clocks) == len(detunings) == len(wins)):
            raise ValueError(
                f"clocks, detunings and wins describe {len(clocks)}, "
                f"{len(detunings)} and {len(wins)} clocks; every clock in the "
                f"declared population carries all three, including the ones "
                f"that never won")
        committed = _require_int_count(self.committed_trials,
                                       "committed_trials")
        multiplicity = tuple(
            (_require_int_count(k, "multiplicity key"),
             _require_int_count(n, "multiplicity count"))
            for k, n in self.multiplicity)
        if any(k < 1 for k, _ in multiplicity):
            raise ValueError("a committed trial has at least one winner, so "
                             "multiplicity is keyed from one upwards")
        if sum(n for _, n in multiplicity) != committed:
            raise ValueError(
                f"the multiplicity histogram covers "
                f"{sum(n for _, n in multiplicity)} trials and "
                f"{committed} committed; every committed trial has exactly one "
                f"winner count")
        if sum(k * n for k, n in multiplicity) != sum(wins):
            raise ValueError(
                f"the multiplicity histogram implies "
                f"{sum(k * n for k, n in multiplicity)} co-completers and the "
                f"per-clock wins sum to {sum(wins)}; a tie is retained in both "
                f"places or in neither")
        object.__setattr__(self, "clocks", clocks)
        object.__setattr__(self, "detunings", detunings)
        object.__setattr__(self, "wins", wins)
        object.__setattr__(self, "committed_trials", committed)
        object.__setattr__(self, "multiplicity", multiplicity)

    @property
    def co_completers(self) -> int:
        """Total winners over the run: ``sum(wins)``, ties included."""
        return sum(self.wins)

    @property
    def tied_trials(self) -> int:
        """Committed trials in which more than one clock completed."""
        return sum(n for k, n in self.multiplicity if k > 1)


def commitment_counts(ledger) -> CommitmentCounts:
    """Return the four-category tally of a verified closed run.

    Reads the trial ledger's ``category`` column, which
    :func:`~.raw_ledger.require_closed_ledger` has already re-derived from the
    authoritative per-clock rows — so this is a count of categories that were
    checked, not of categories that were asserted.
    """
    # Bound, not merely called.  The gate returns a *freshly reconciled*
    # snapshot and everything below reads that one, so this door cannot
    # validate one object and go on computing from another.
    ledger = require_unmutated_ledger(require_instance(ledger, ClosedLedger,
                                                       "ledger"))
    tally = {name: 0 for name in TRIAL_CATEGORIES}
    for row in ledger.rows("ledger"):
        tally[row["category"]] += 1
    return CommitmentCounts(ledger.manifest.trials, tally["committed"],
                            tally["never_eligible"], tally["lock_failed"],
                            tally["dwell_failed"])


def endpoint_cumulative_hazard(counts) -> float:
    """``-log(1 - committed/trials)`` for one closed run, as a display value.

    Returned as ``+inf`` when every trial committed, which is the honest answer
    and the reason this number is a *display* quantity rather than a regression
    input.  ``analysis.py`` never calls it: fitting the logarithm of an
    empirical proportion is what loses the zero-event and all-event cells, and
    the complementary-log-log likelihood there works from the counts.
    """
    require_instance(counts, CommitmentCounts, "counts")
    q = counts.commitment_probability
    if q >= 1.0:
        return float("inf")
    return -math.log1p(-q)


def uniform_edges(window_open, window_close, bins):
    """Return ``bins + 1`` equally spaced edges spanning the whole window.

    The endpoints are exact — ``edges[0]`` is ``window_open`` and ``edges[-1]``
    is ``window_close`` — because :func:`risk_table` requires a table to cover
    the run's entire observation window.  That requirement is what stops a
    binning from becoming a window choice: there is no argument here that could
    crop the sweep, and a caller who wants different resolution changes the
    number of bins, never the span.
    """
    opened = require_finite(window_open, "window_open")
    closed = require_finite(window_close, "window_close")
    count = require_positive_int(bins, "bins", minimum=1)
    if not closed > opened:
        raise ValueError(f"window_close={closed!r} must follow "
                         f"window_open={opened!r}")
    edges = np.linspace(opened, closed, count + 1)
    # linspace is exact at both endpoints; restated rather than assumed,
    # because the table gate below compares them for equality.
    edges[0], edges[-1] = opened, closed
    return tuple(float(edge) for edge in edges)


def risk_table(ledger, edges) -> RiskTable:
    """Return the exact risk table of a closed run over ``edges``.

    ``edges`` must be strictly increasing, must begin at the manifest's
    ``window_open`` and must end at its ``window_close``.  Both requirements are
    deliberate and both are refusals rather than adjustments:

    * a table that began late would silently condition on surviving to its
      first edge, and a table that ended early would silently discard the
      trials still at risk; and
    * either one is a *window choice made after seeing the data*, which is the
      selection the plan forbids and which no amount of care at the call site
      would make visible in the resulting table.

    Events are the committed trials' ``commit_time``.  Censorings are the
    unresolved trials, all of them at ``window_close`` — this run's
    administrative right-censoring — which therefore all land in the final
    interval.  The two are counted separately and never pooled.
    """
    # Bound, not merely called.  The gate returns a *freshly reconciled*
    # snapshot and everything below reads that one, so this door cannot
    # validate one object and go on computing from another.
    ledger = require_unmutated_ledger(require_instance(ledger, ClosedLedger,
                                                       "ledger"))
    raw = np.atleast_1d(require_finite_array(edges, "edges"))
    if raw.ndim != 1:
        raise ValueError(f"edges must be one-dimensional, got shape "
                         f"{raw.shape}")
    if raw.size < 2:
        raise ValueError(f"a risk table needs at least two edges, got "
                         f"{raw.size}")
    if not bool(np.all(np.diff(raw) > 0.0)):
        raise ValueError("edges must be strictly increasing; a repeated edge "
                         "is a zero-width interval and a decreasing one is "
                         "not a time axis")

    manifest = ledger.manifest
    opened = manifest.value("window_open")
    closed = manifest.value("window_close")
    if float(raw[0]) != opened or float(raw[-1]) != closed:
        raise ValueError(
            f"edges span [{float(raw[0])!r}, {float(raw[-1])!r}] and this run's "
            f"observation window is [{opened!r}, {closed!r}].  A risk table "
            f"covers the whole window: a table starting late conditions on "
            f"survival to its first edge and one ending early discards the "
            f"trials still at risk, and either is a window chosen after the "
            f"fact")

    events, censorings = [], []
    for row in ledger.rows("ledger"):
        if row["category"] == "committed":
            events.append(row["commit_time"])
        else:
            censorings.append(closed)

    bounds = np.asarray(raw, dtype=float)
    # ``searchsorted`` with side="left" over the interior edges puts an instant
    # t into bucket i when edges[i] < t <= edges[i+1], which is exactly
    # EVENT_CONVENTION.  window_close itself lands in the final bucket.
    interior = bounds[1:-1]
    n_bins = bounds.size - 1

    def bucket(values):
        counts = np.zeros(n_bins, dtype=np.int64)
        if not values:
            return counts
        found = np.searchsorted(interior, np.asarray(values, dtype=float),
                                side="left")
        for index in found:
            counts[int(index)] += 1
        return counts

    event_counts = bucket(events)
    censor_counts = bucket(censorings)

    intervals, at_risk = [], manifest.trials
    for index in range(n_bins):
        interval = RiskInterval(float(bounds[index]), float(bounds[index + 1]),
                                at_risk, int(event_counts[index]),
                                int(censor_counts[index]))
        intervals.append(interval)
        at_risk -= interval.exits
    return RiskTable(tuple(intervals), manifest.trials, opened)


def winner_mismatch(ledger) -> WinnerMismatch:
    """Return which clocks supplied first commitments in a closed run.

    Every clock in the declared population appears, whether it ever won or not,
    and every co-completer of a tie is counted.  The two totals that a reader
    is most likely to confuse — committed trials and co-completers — are
    carried side by side by :class:`WinnerMismatch`, together with the full
    multiplicity histogram, so a tie shows up as a tie rather than as a count
    that does not add up.
    """
    # Bound, not merely called.  The gate returns a *freshly reconciled*
    # snapshot and everything below reads that one, so this door cannot
    # validate one object and go on computing from another.
    ledger = require_unmutated_ledger(require_instance(ledger, ClosedLedger,
                                                       "ledger"))
    manifest = ledger.manifest
    detunings = manifest.value("detunings")
    wins = [0] * manifest.clocks
    multiplicity = {}
    committed = 0
    for row in ledger.rows("ledger"):
        winners = row["winners"]
        if not winners:
            continue
        committed += 1
        multiplicity[len(winners)] = multiplicity.get(len(winners), 0) + 1
        for clock in winners:
            wins[clock] += 1
    return WinnerMismatch(tuple(range(manifest.clocks)), tuple(detunings),
                          tuple(wins), committed,
                          tuple(sorted(multiplicity.items())))
