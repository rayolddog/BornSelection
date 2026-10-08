#!/usr/bin/env python3
"""
raw_ledger.py — what a raw run writes, and what a reader must prove before it
may believe a single row of it.

The problem this solves
-----------------------
Structural isolation stops the raw modules from importing the prediction.  It
says nothing about *time*.  A comparison process that opens a ledger while the
raw run is still writing it can, entirely by accident, become part of the
experiment: trials that have not happened yet are missing rather than
unresolved, an interrupted sweep looks like a completed one with fewer
commitments, and a rerun that appends to an existing file mixes two random
streams into one apparent result.  None of those announce themselves.  Each
produces a smaller, cleaner-looking, entirely fictitious ledger.

So the contract is one sentence: **a raw ledger may be read only after the
process that wrote it has closed it, and only if its bytes, its manifest, its
row counts and the content of its rows all still agree with what the close
marker says.**  The comparison process — the one allowed to import the
prediction — calls :func:`require_closed_ledger` before it looks at a single
row, and gets an exception rather than a plausible number if any table is
truncated, still being written, appended to, edited by hand, paired with
another run's manifest, missing a trial, holding a trial twice, interleaved
with a second run's rows, or accompanied by a marker whose own fields were
changed after it was built.

What changed at schema 3
------------------------
An independent review found that version 2 kept the *aggregate* per-clock table
and threw the per-trial histories away.  Totals are not a substitute: different
trial-level histories produce the same totals and different censoring
conclusions, and a winner's last-reset time was gone permanently.  The same
review forged three records that passed every version-2 check after their
digests were recomputed — an unresolved category contradicting its own clock
counts, a manifest detuning bound to different clock detunings, and 999 units
of exposure inside a 0.4-unit window.

So version 3 makes two changes and they are the same change.  ``clocks.csv``
becomes **one authoritative row per (trial, clock)** — the censored history of
every clock in every trial, with its entry, reset, endpoint, exposure and
opportunity fields and its own category — and the shadow table carries the same
fields uncensored.  There is then no unbound summary anywhere in a closed run:
:meth:`ClosedLedger.clock_totals` derives the aggregate *after* reading, from
rows the gate has reconciled, rather than reading a total nothing could check.

And because every channel-level number is now derivable from those rows,
:func:`_require_consistent_rows` derives it and compares.  Categories, winners,
co-completers, detunings, counts, exposure and opportunity are recomputed from
the authoritative table and required to match what the trial view claims.  A
forger who recomputes every digest now has to produce a *consistent* run, which
is a different and much larger problem than producing a well-typed one.

What changed at schema 2
------------------------
Ticket 01 built this module as an *interface*: a close marker, a three-field
manifest, and a gate that verified the ledger bytes against the marker.  It
could not do more, because no ledger existed and a row parser written before
the rows would have been a guess.

The rows exist now, so three things move.

**The manifest describes the run.**  Three fields — layout version, pilot/
production label, trial count — were the right amount to freeze when there were
no run parameters to freeze.  A manifest that cannot state the clock grid, the
pulse, the noise, the dwell, the mesh, the key namespaces or the numerical
status it was produced under cannot support the comparison it exists for.
:data:`MANIFEST_FIELDS` is the whole declared schema, and
:func:`parse_manifest` refuses any key set that is not exactly it.

**The trial count is checked against rows, not against a claim.**  Ticket 01's
gate could only establish *the marker's count agrees with the frozen manifest
bound to that marker*, and said so in its own docstring: "counting parsed rows
is the stronger form of the same check and should replace this one when the
parser exists, rather than sit beside it."  The parser exists.  The gate now
parses every table and requires the manifest's count, the marker's count and
the number of rows actually present to be one number — and requires the *row
identities* to be exactly the ones the run declared, which is what catches a
ledger that has the right number of rows because two runs' rows were
interleaved.

**There are three tables, not one.**  A trial ledger, a per-clock ledger and
the no-stop shadow sample.  The marker carries a count and a digest for each,
because a run whose trial ledger closed and whose shadow sample did not is not
a closed run.

At schema 2 both versions moved to 2, and a version-1 marker or manifest was
refused rather than partially read.  Nothing had been recorded under version 1.

What the digests do and do not establish
----------------------------------------
They establish that the bytes a reader has are the bytes the writer finished
with.  They establish nothing about whether those bytes are *correct* physics,
and they are not signatures: anyone able to rewrite a ledger can recompute
them.  This is protection against truncation, interleaving and accidental edit,
which is what actually goes wrong, and it is not protection against a
determined author of the file.
"""
from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass

from .raw_config import RUN_LABELS
from .validation import (require_exact_type, require_finite,
                         require_positive_int)

__all__ = [
    "CLOCK_COLUMNS",
    "COLUMN_KINDS",
    "DIGEST_ALGORITHM",
    "DIGEST_HEX_LENGTH",
    "DRIVE_LABELS",
    "LEDGER_COLUMNS",
    "LEDGER_SCHEMA_VERSION",
    "LIST_SEPARATOR",
    "MANIFEST_FIELDS",
    "MANIFEST_KEYS",
    "MANIFEST_SCHEMA_VERSION",
    "MISSING_ENTRY",
    "MODEL_LABELS",
    "SHADOW_COLUMNS",
    "TABLE_NAMES",
    "TRIAL_CATEGORIES",
    "CloseMarker",
    "ClosedLedger",
    "RawManifest",
    "canonical_manifest_bytes",
    "canonical_marker_bytes",
    "digest_bytes",
    "encode_row",
    "encode_table",
    "header_row",
    "parse_manifest",
    "parse_marker",
    "parse_table",
    "require_closed_ledger",
    "require_unmutated_ledger",
]

#: The one hash this package uses.  Named so that a future change to it is a
#: visible edit here rather than a silent difference between writer and reader.
DIGEST_ALGORITHM = "sha256"

#: Length of a :data:`DIGEST_ALGORITHM` digest written as lowercase hex.
DIGEST_HEX_LENGTH = 64

#: Version of the ledger layout a marker describes.  A reader that understands
#: version 3 must refuse a version 1, 2 or 4 ledger rather than read the
#: columns it recognises and quietly ignore the rest.
LEDGER_SCHEMA_VERSION = 3

#: Version of the *manifest* layout, which is a different thing from the ledger
#: layout above and versioned separately because the two will not grow at the
#: same rate.
MANIFEST_SCHEMA_VERSION = 3

#: The three tables a closed raw run consists of, in the order a marker lists
#: them.  Naming them here is what lets the gate say "the shadow sample is
#: missing" rather than silently verifying two files out of three.
TABLE_NAMES = ("ledger", "clocks", "shadow")

#: How a multi-valued cell is written inside one CSV field.  A pipe rather than
#: a comma for the obvious reason, and rather than a semicolon so that a reader
#: opening the file in a spreadsheet sees one column holding a list instead of
#: a row whose width depends on how many clocks happened to tie.
LIST_SEPARATOR = "|"

#: The four mutually exclusive, exhaustive channel categories.
#:
#: Duplicated from ``commitment.TRIAL_CATEGORIES`` rather than imported, for
#: the reason ``raw_config.RawClockGrid`` re-derives the analytic grid: this
#: module is the *format*, and a reader validating a ledger should not have to
#: load the state machine that produced it.  ``verify.py`` is allowed to see
#: both and asserts the two tuples are identical, so the duplication is checked
#: rather than trusted.
TRIAL_CATEGORIES = ("committed", "never_eligible", "lock_failed",
                    "dwell_failed")

#: The data-generating processes a manifest may name, duplicated from
#: ``dynamics.MODEL_LABELS`` under the same rule and cross-checked the same way.
MODEL_LABELS = ("full", "width_only")

#: The drives a manifest may name, duplicated from ``dynamics.DRIVE_LABELS``.
DRIVE_LABELS = ("pulsed", "stationary")

#: Every value kind a manifest field or a ledger column may hold.
#:
#: ``int``      any integer, ``int0`` non-negative, ``int+`` at least one
#: ``float``    finite, ``float0`` finite and non-negative, ``float+`` positive
#: ``time``     a finite instant, ``time?`` one or absent
#: ``bool``     a genuine boolean, written ``true``/``false``
#: ``text``     a non-empty string
#: ``texts``    a list of non-empty strings
#: ``floats``   a list of finite numbers
#: ``times``    a list of instants, any of which may be absent
#: ``ints``     a list of non-negative integers, possibly empty
#: ``category`` one of :data:`TRIAL_CATEGORIES`
#: ``run_label``/``model``/``drive`` one of the declared vocabularies
#: ``digest``   a lowercase hex digest of :data:`DIGEST_HEX_LENGTH` characters
COLUMN_KINDS = ("int", "int0", "int+", "float", "float0", "float+", "time",
                "time?", "bool", "text", "texts", "floats", "times", "ints",
                "category", "run_label", "model", "drive", "digest")

#: How an absent entry is spelled *inside* a list cell.
#:
#: An empty field already means the empty list, so a one-element list holding
#: nothing could not be written as one -- ``[None]`` and ``[]`` would be the
#: same bytes.  A winner that completed its dwell on the first attempt has no
#: last-reset time and is an ordinary case, so the two have to be tellable
#: apart.
MISSING_ENTRY = "-"

#: The complete version-3 manifest schema: every key, in order, with its kind.
#:
#: This is the whole freeze.  A comparison process computes its prediction from
#: these numbers, so a manifest that could omit one, hold one of the wrong
#: kind, or carry an extra key this version does not understand would let a
#: ledger be compared against parameters that are not its own.
#:
#: Deliberately **absent**: the clock-block and step-window sizes the run was
#: executed with.  Those bound the live working set and cannot change a
#: recorded number — the same run at any block size produces the same three
#: tables and the same digests — so freezing them into the immutable
#: description of the run would imply they were part of it.  They are reported
#: by the runner as execution facts instead.
MANIFEST_FIELDS = (
    ("schema_version", "int+"),
    ("ledger_schema_version", "int+"),
    ("run_label", "run_label"),
    ("model", "model"),
    ("drive", "drive"),
    ("trials", "int+"),
    ("shadow_trials", "int+"),
    ("shadow_rule", "text"),
    ("clocks", "int+"),
    ("grid_half_width", "float+"),
    ("grid_spacing", "float+"),
    ("grid_parity", "text"),
    ("grid_origin", "float"),
    ("grid_support", "floats"),
    ("detunings", "floats"),
    ("peak_coupling", "float0"),
    ("pulse_duration", "float+"),
    ("pulse_centre", "float"),
    ("pulse_support", "floats"),
    ("pulse_digest", "text"),
    ("phase_diffusion", "float0"),
    ("lock_tolerance", "float+"),
    ("dwell_time", "float+"),
    ("fixed_contraction_rate", "float0"),
    ("timestep", "float+"),
    ("mesh_origin", "float"),
    ("mesh_identity", "text"),
    ("start_index", "int0"),
    ("steps", "int+"),
    ("window_open", "float"),
    ("window_close", "float"),
    ("key_schema", "text"),
    ("stream_namespace", "text"),
    ("stream_prefix", "text"),
    ("ensemble_prefix", "text"),
    ("ledger_columns", "texts"),
    ("clock_columns", "texts"),
    ("shadow_columns", "texts"),
    ("endpoint_status", "text"),
    ("moving_band_status", "text"),
    ("numerical_gate", "text"),
    ("validation_note", "text"),
    ("python_version", "text"),
    ("numpy_version", "text"),
    ("non_claims", "texts"),
)

#: Just the names, for the key-set comparison :func:`parse_manifest` makes.
MANIFEST_KEYS = tuple(name for name, _ in MANIFEST_FIELDS)

#: One row per trial: the channel-level ledger.
#:
#: ``winners`` is a *list*, and it is a list on purpose.  Several clocks can
#: complete on one endpoint, all of them are retained, and no random internal
#: winner is ever chosen — so the column that would have held "the winner" holds
#: every co-completer instead, and ``co_completion`` states how many there were.
#: A schema with room for exactly one winner would have made the tie
#: unrecordable and the plan's refinement gate on co-commit multiplicity
#: unanswerable.
LEDGER_COLUMNS = (
    ("trial", "int0"),
    ("category", "category"),
    ("stop_time", "time"),
    ("commit_time", "time?"),
    ("censored", "bool"),
    ("co_completion", "int0"),
    ("winners", "ints"),
    ("winner_detunings", "floats"),
    ("winner_phases", "floats"),
    ("winner_band_entry", "times"),
    ("winner_first_eligible", "times"),
    ("winner_last_reset", "times"),
    ("clocks_ever_eligible", "int0"),
    ("clocks_ever_inside", "int0"),
    ("band_entries", "int0"),
    ("dwell_resets", "int0"),
    ("eligible_time", "float0"),
)

#: One authoritative row per ``(trial, clock)``: what every clock in the
#: population was doing when the channel stopped.
#:
#: This is the durable record the plan's censoring analysis needs, and version 2
#: did not have it.  It held an aggregate — one row per clock, totals across
#: every trial — which cannot say which trial contributed an entry, a reset, an
#: exposure or a win.  Two runs with identical totals and different per-trial
#: histories support different conclusions about first-winner censoring, and a
#: closed run that kept only the totals could never be asked again.
#:
#: ``exposure_time`` and ``eligible_time`` are why the table exists.  A
#: first-passage race censors every clock but the winners, so every count here
#: is a count over the opportunity that clock actually had, and dividing by
#: anything else — trials, wall time, the pulse duration — quietly treats a
#: stopped clock as fully observed.  Both are *derivable*: exposure is the
#: channel stop (or this clock's own commitment, whichever came first) minus the
#: window's opening, and the gate recomputes it rather than believing it.
CLOCK_COLUMNS = (
    ("trial", "int0"),
    ("clock", "int0"),
    ("detuning", "float"),
    ("category", "category"),
    ("committed_at", "time?"),
    ("first_eligible_at", "time?"),
    ("first_inside_at", "time?"),
    ("last_reset_at", "time?"),
    ("eligibility_entries", "int0"),
    ("band_entries", "int0"),
    ("dwell_resets", "int0"),
    ("endpoints", "int0"),
    ("eligible_endpoints", "int0"),
    ("inside_endpoints", "int0"),
    ("exposure_time", "float0"),
    ("eligible_time", "float0"),
    ("final_phase", "float"),
    ("won", "bool"),
    ("co_completed", "bool"),
)

#: One row per ``(shadow trial, clock)``: the predeclared no-stop sample, with
#: the same fields as the censored table above and without the two that only a
#: race can have.
#:
#: The same physical keys as the censored run of the same trial, with the
#: first-winner stop switched off.  A commitment recorded here is one that in
#: the actual race did *not* happen, because the channel had already stopped;
#: that difference is the measurement, and holding the two tables in the same
#: shape is what makes it a subtraction rather than an interpretation.
#:
#: ``won`` and ``co_completed`` are absent on purpose.  There is no channel stop
#: in a no-stop walk, so there is nothing to win.
SHADOW_COLUMNS = tuple(
    column for column in CLOCK_COLUMNS
    if column[0] not in ("won", "co_completed"))

#: The columns of each declared table, by name.
_TABLE_COLUMNS = {"ledger": LEDGER_COLUMNS, "clocks": CLOCK_COLUMNS,
                  "shadow": SHADOW_COLUMNS}

#: Lowercase hex digits.  A digest written in uppercase, or with whitespace, is
#: rejected rather than normalised: two spellings of one digest is how a
#: comparison that should have failed ends up passing.
_HEX_DIGITS = frozenset("0123456789abcdef")

#: The two spellings a boolean cell may have.  Not ``True``/``False``, so that
#: a Python ``repr`` accidentally reaching a CSV writer is a visible defect
#: rather than a second accepted spelling.
_BOOL_TEXT = {True: "true", False: "false"}
_TEXT_BOOL = {"true": True, "false": False}


def _require_digest(value, name: str) -> str:
    """Return ``value`` as a canonical lowercase hex digest, or raise."""
    if isinstance(value, bool) or not isinstance(value, str):
        raise TypeError(f"{name} must be a hex string, got "
                        f"{type(value).__name__} {value!r}")
    if len(value) != DIGEST_HEX_LENGTH:
        raise ValueError(f"{name} must be {DIGEST_HEX_LENGTH} hex characters "
                         f"({DIGEST_ALGORITHM}), got {len(value)}")
    bad = sorted(set(value) - _HEX_DIGITS)
    if bad:
        raise ValueError(f"{name} must be lowercase hex, got {value!r} "
                         f"containing {bad!r}")
    return value


def digest_bytes(payload) -> str:
    """Return the canonical digest of ``payload`` as lowercase hex.

    Bytes only.  A ``str`` is refused rather than encoded, because the encoding
    would be this function's guess and the writer's guess might differ, which
    turns a matching ledger into a mismatched one for a reason that has nothing
    to do with the ledger.
    """
    if not isinstance(payload, (bytes, bytearray, memoryview)):
        raise TypeError(f"payload must be bytes, got {type(payload).__name__}; "
                        f"text is refused rather than encoded, so that writer "
                        f"and reader cannot disagree about the encoding")
    return hashlib.new(DIGEST_ALGORITHM, bytes(payload)).hexdigest()


# ---------------------------------------------------------------------------
# Values: one kind system for manifest fields and ledger cells alike
# ---------------------------------------------------------------------------
def _require_int(value, name: str, minimum) -> int:
    """Return ``value`` as an exact ``int`` at or above ``minimum``."""
    if isinstance(value, bool):
        raise TypeError(f"{name} must be an integer, got bool {value!r}")
    if not isinstance(value, int):
        raise TypeError(f"{name} must be an integer, got "
                        f"{type(value).__name__} {value!r}")
    if value < minimum:
        raise ValueError(f"{name} must be at least {minimum}, got {value!r}")
    return value


def _require_float(value, name: str, kind: str) -> float:
    """Return ``value`` as a finite float obeying ``kind``'s sign rule."""
    out = require_finite(value, name)
    if kind == "float0" and out < 0.0:
        raise ValueError(f"{name} must be non-negative, got {out!r}")
    if kind == "float+" and not out > 0.0:
        raise ValueError(f"{name} must be strictly positive, got {out!r}")
    return out


def _require_text(value, name: str) -> str:
    """Return ``value`` as a non-empty single-line string with no separator.

    A newline would end a row and the list separator would end a field, so
    either one silently changes how many values a reader finds.  Both are
    refused rather than quoted or escaped, so that one row is always one line
    and one list entry is always one field.

    A comma is *not* refused here, because manifest values are JSON and prose
    with commas in it belongs there.  What a comma must never do is reach a CSV
    cell, and :func:`encode_row` checks that directly rather than banning the
    character everywhere it could not do harm.
    """
    if isinstance(value, bool) or not isinstance(value, str):
        raise TypeError(f"{name} must be a string, got "
                        f"{type(value).__name__} {value!r}")
    if not value:
        raise ValueError(f"{name} must not be empty")
    bad = sorted(set(value) & set("\r\n" + LIST_SEPARATOR))
    if bad:
        raise ValueError(
            f"{name} contains {bad!r}, which would change how many rows or how "
            f"many list entries a reader finds; the value is refused rather "
            f"than quoted or escaped")
    return value


def _require_member(value, name: str, allowed) -> str:
    """Return ``value`` as one of ``allowed``, or raise naming ``name``."""
    if isinstance(value, bool) or not isinstance(value, str):
        raise TypeError(f"{name} must be a string, got "
                        f"{type(value).__name__} {value!r}")
    if value not in allowed:
        raise ValueError(f"{name} must be one of {tuple(allowed)!r}, got "
                         f"{value!r}")
    return value


def _require_sequence(value, name: str):
    """Return ``value`` as a list, refusing the things that look like one."""
    if isinstance(value, (str, bytes, bytearray, dict)):
        raise TypeError(f"{name} must be a list of values, got "
                        f"{type(value).__name__} {value!r}")
    try:
        return list(value)
    except TypeError:
        raise TypeError(f"{name} must be a list of values, got "
                        f"{type(value).__name__} {value!r}") from None


def _normalize(value, name: str, kind: str):
    """Return ``value`` normalized to ``kind``, or raise naming ``name``.

    The single place a declared kind becomes a rule.  Every manifest field and
    every ledger cell passes through here, on the way in and on the way back
    out, so a value that survives a write/read round trip is a value both sides
    validated by the same code.
    """
    if kind == "int":
        return _require_int(value, name, -(1 << 62))
    if kind == "int0":
        return _require_int(value, name, 0)
    if kind == "int+":
        return _require_int(value, name, 1)
    if kind in ("float", "float0", "float+", "time"):
        return _require_float(value, name, kind)
    if kind == "time?":
        return None if value is None else _require_float(value, name, "float")
    if kind == "bool":
        if not isinstance(value, bool):
            raise TypeError(f"{name} must be a bool, got "
                            f"{type(value).__name__} {value!r}")
        return value
    if kind == "text":
        return _require_text(value, name)
    if kind == "texts":
        return tuple(_require_text(item, f"{name}[{i}]")
                     for i, item in enumerate(_require_sequence(value, name)))
    if kind == "floats":
        return tuple(_require_float(item, f"{name}[{i}]", "float")
                     for i, item in enumerate(_require_sequence(value, name)))
    if kind == "times":
        return tuple(None if item is None
                     else _require_float(item, f"{name}[{i}]", "float")
                     for i, item in enumerate(_require_sequence(value, name)))
    if kind == "ints":
        return tuple(_require_int(item, f"{name}[{i}]", 0)
                     for i, item in enumerate(_require_sequence(value, name)))
    if kind == "category":
        return _require_member(value, name, TRIAL_CATEGORIES)
    if kind == "run_label":
        return _require_member(value, name, RUN_LABELS)
    if kind == "model":
        return _require_member(value, name, MODEL_LABELS)
    if kind == "drive":
        return _require_member(value, name, DRIVE_LABELS)
    if kind == "digest":
        return _require_digest(value, name)
    raise ValueError(f"unknown value kind {kind!r} for {name}; the declared "
                     f"kinds are {list(COLUMN_KINDS)}")


def _cell_text(value, kind: str) -> str:
    """Return the one spelling of ``value`` a ledger cell may carry.

    Floats go through ``repr``, which is the shortest decimal that reads back
    as the identical binary64.  A rounded cell would make the ledger a lossy
    record of its own run, and the digest would then be certifying bytes that
    no longer say what the run said.
    """
    if kind == "bool":
        return _BOOL_TEXT[value]
    if kind == "time?":
        return "" if value is None else repr(value)
    if kind in ("float", "float0", "float+", "time"):
        return repr(value)
    if kind in ("texts", "floats", "ints", "times"):
        inner = {"texts": "text", "floats": "float", "ints": "int0",
                 "times": "float"}[kind]
        return LIST_SEPARATOR.join(
            MISSING_ENTRY if item is None else _cell_text(item, inner)
            for item in value)
    return str(value)


def _cell_value(text: str, name: str, kind: str):
    """Return the value a ledger cell spells, or raise naming its column."""
    if kind in ("int", "int0", "int+"):
        try:
            parsed = int(text)
        except ValueError:
            raise ValueError(f"{name} is not an integer: {text!r}") from None
        if str(parsed) != text:
            raise ValueError(
                f"{name} spells an integer as {text!r}; a canonical cell has "
                f"no sign on a non-negative value, no leading zero, no "
                f"underscore and no whitespace")
        return _normalize(parsed, name, kind)
    if kind == "bool":
        if text not in _TEXT_BOOL:
            raise ValueError(f"{name} must be 'true' or 'false', got {text!r}")
        return _TEXT_BOOL[text]
    if kind == "time?" and text == "":
        return None
    if kind in ("float", "float0", "float+", "time", "time?"):
        try:
            parsed = float(text)
        except ValueError:
            raise ValueError(f"{name} is not a number: {text!r}") from None
        value = _normalize(parsed, name, "float" if kind == "time?" else kind)
        if repr(value) != text:
            raise ValueError(
                f"{name} spells {value!r} as {text!r}; a canonical cell is the "
                f"shortest decimal that reads back as the identical binary64, "
                f"so an alias would let two files disagree about one number "
                f"while both parsing")
        return value
    if kind in ("texts", "floats", "ints", "times"):
        inner = {"texts": "text", "floats": "float", "ints": "int0",
                 "times": "float"}[kind]
        if text == "":
            return ()
        return tuple(None if part == MISSING_ENTRY and kind == "times"
                     else _cell_value(part, f"{name}[{i}]", inner)
                     for i, part in enumerate(text.split(LIST_SEPARATOR)))
    return _normalize(text, name, kind)


def _require_columns(columns, name: str):
    """Return ``columns`` as a validated ``((name, kind), ...)`` tuple."""
    prepared = []
    seen = set()
    for position, column in enumerate(_require_sequence(columns, name)):
        pair = tuple(_require_sequence(column, f"{name}[{position}]"))
        if len(pair) != 2:
            raise ValueError(f"{name}[{position}] must be a (name, kind) pair, "
                             f"got {pair!r}")
        label = _require_text(pair[0], f"{name}[{position}] name")
        kind = _require_member(pair[1], f"{name}[{position}] kind",
                               COLUMN_KINDS)
        if label in seen:
            raise ValueError(f"{name} names the column {label!r} twice; a "
                             f"repeated column makes a row's meaning depend on "
                             f"which occurrence a reader keeps")
        seen.add(label)
        prepared.append((label, kind))
    if not prepared:
        raise ValueError(f"{name} must declare at least one column")
    return tuple(prepared)


def header_row(columns) -> str:
    """Return the header line of a table with these columns.

    Written into every table and checked on every read.  It is what makes a
    file self-describing, and — because :func:`parse_table` requires it to be
    exactly the declared column names in the declared order — what stops a
    reader silently mapping a run's columns onto a different schema's.
    """
    return ",".join(name for name, _ in _require_columns(columns, "columns"))


def encode_row(columns, values) -> str:
    """Return the one canonical line these column values spell.

    ``values`` is a mapping from column name to value, and its key set must be
    *exactly* the declared columns: a missing key is a hole and an extra key is
    a value the schema has no place for, and guessing which of the two happened
    is how a table ends up with rows of two different shapes.
    """
    columns = _require_columns(columns, "columns")
    if not isinstance(values, dict):
        raise TypeError(f"values must be a mapping of column name to value, "
                        f"got {type(values).__name__}")
    declared = {name for name, _ in columns}
    present = set(values)
    missing = sorted(declared - present)
    unknown = sorted(present - declared)
    if missing or unknown:
        raise ValueError(f"a row has exactly the columns {sorted(declared)}; "
                         f"missing {missing}, unexpected {unknown}")
    cells = [_cell_text(_normalize(values[name], name, kind), kind)
             for name, kind in columns]
    for (name, _), cell in zip(columns, cells):
        if "," in cell or "\n" in cell or "\r" in cell:
            raise ValueError(
                f"the {name!r} cell spells {cell!r}, which contains a "
                f"separator or a line break and would change how many fields "
                f"a reader finds in this row")
    return ",".join(cells)


def parse_table(payload, columns):
    """Return the rows these bytes hold, typed, or raise.

    A tuple of dictionaries, one per data row, with every cell converted to the
    value its declared kind names — so a caller receives ``None`` for an
    absent commitment rather than an empty string it has to remember to check,
    and a category that is not one of the four declared ones fails here rather
    than downstream.

    Strict about everything a table could be wrong about: the bytes must be
    UTF-8, the file must end in exactly one newline per row, the header must be
    the declared columns in the declared order, every row must have exactly as
    many fields as there are columns, and every cell must be the canonical
    spelling of its value.  A ledger with a row too few, a row too many, a
    column renamed, a number rewritten by a spreadsheet, or a stray blank line
    is refused rather than partially read.
    """
    columns = _require_columns(columns, "columns")
    if not isinstance(payload, (bytes, bytearray, memoryview)):
        raise TypeError(f"a table must be bytes, got "
                        f"{type(payload).__name__}; text is refused rather "
                        f"than encoded, so that writer and reader cannot "
                        f"disagree about the encoding")
    try:
        text = bytes(payload).decode("utf-8")
    except UnicodeDecodeError as exc:
        raise ValueError(f"the table is not valid UTF-8: {exc}") from None
    if not text:
        raise ValueError("the table is empty; a closed raw run wrote at least "
                         "a header, so empty means the run was interrupted or "
                         "the wrong file was opened")
    if not text.endswith("\n"):
        raise ValueError(
            "the table does not end in a newline, which is what a truncated "
            "file looks like: a writer that was still running had not written "
            "the last line's terminator yet")
    if "\r" in text:
        raise ValueError("the table contains a carriage return; it has been "
                         "rewritten by a tool that changed its line endings, "
                         "so its bytes are no longer the writer's bytes")
    lines = text[:-1].split("\n")
    expected = header_row(columns)
    if lines[0] != expected:
        raise ValueError(
            f"the header is {lines[0]!r}; a version-{LEDGER_SCHEMA_VERSION} "
            f"table of this kind has exactly {expected!r}.  A table with a "
            f"different column set or a different column order belongs to a "
            f"different schema and must not be read as this one")

    rows = []
    for number, line in enumerate(lines[1:], start=2):
        fields = line.split(",")
        if len(fields) != len(columns):
            raise ValueError(
                f"line {number} has {len(fields)} fields, not the "
                f"{len(columns)} the header declares: {line!r}")
        rows.append({name: _cell_value(field, f"line {number} column {name}",
                                       kind)
                     for (name, kind), field in zip(columns, fields)})
    return tuple(rows)


# ---------------------------------------------------------------------------
# The manifest
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class RawManifest:
    """The frozen description of one raw run, in the only shape there is.

    Built from a mapping and stored as an ordered tuple of validated pairs, so
    that a manifest is immutable in the way the rest of this package means the
    word: the fields cannot be reassigned, the values are normalized copies,
    and every list is a tuple.

    One constructor argument rather than forty-four fields, because forty-four
    fields would be forty-four places for a caller to pass the right value to
    the wrong name.  A mapping is checked against :data:`MANIFEST_FIELDS` as a
    whole: exactly those keys, exactly those kinds, no more and no fewer.
    """

    fields: tuple

    def __post_init__(self):
        supplied = self.fields
        if isinstance(supplied, dict):
            given = dict(supplied)
        else:
            pairs = _require_sequence(supplied, "fields")
            given = {}
            for position, pair in enumerate(pairs):
                item = tuple(_require_sequence(pair, f"fields[{position}]"))
                if len(item) != 2:
                    raise ValueError(f"fields[{position}] must be a "
                                     f"(name, value) pair, got {item!r}")
                if item[0] in given:
                    raise ValueError(f"fields names {item[0]!r} twice")
                given[item[0]] = item[1]
        for key in given:
            if not isinstance(key, str):
                raise TypeError(f"manifest field names must be strings, got "
                                f"{type(key).__name__} {key!r}")
        missing = sorted(set(MANIFEST_KEYS) - set(given))
        unknown = sorted(set(given) - set(MANIFEST_KEYS))
        if missing or unknown:
            raise ValueError(
                f"a version-{MANIFEST_SCHEMA_VERSION} manifest has exactly the "
                f"{len(MANIFEST_KEYS)} keys {list(MANIFEST_KEYS)}; missing "
                f"{missing}, unexpected {unknown}.  A manifest with a "
                f"different key set belongs to a different schema and must not "
                f"be read as this one")

        normalized = tuple((name, _normalize(given[name], name, kind))
                           for name, kind in MANIFEST_FIELDS)
        settled = dict(normalized)
        if settled["schema_version"] != MANIFEST_SCHEMA_VERSION:
            raise ValueError(
                f"schema_version must be {MANIFEST_SCHEMA_VERSION}, got "
                f"{settled['schema_version']}; a manifest written to a layout "
                f"this code does not implement must be refused, not partially "
                f"read")
        if settled["ledger_schema_version"] != LEDGER_SCHEMA_VERSION:
            raise ValueError(
                f"ledger_schema_version must be {LEDGER_SCHEMA_VERSION}, got "
                f"{settled['ledger_schema_version']}")
        if settled["shadow_trials"] > settled["trials"]:
            raise ValueError(
                f"shadow_trials={settled['shadow_trials']} exceeds "
                f"trials={settled['trials']}; the no-stop sample is a subset "
                f"of the run's own trials")
        if len(settled["detunings"]) != settled["clocks"]:
            raise ValueError(
                f"the manifest lists {len(settled['detunings'])} detunings for "
                f"{settled['clocks']} clocks; the population and the grid it "
                f"is drawn from cannot disagree")
        for name, declared in (("ledger_columns", LEDGER_COLUMNS),
                               ("clock_columns", CLOCK_COLUMNS),
                               ("shadow_columns", SHADOW_COLUMNS)):
            expected = tuple(label for label, _ in declared)
            if settled[name] != expected:
                raise ValueError(
                    f"{name} is {list(settled[name])}; this code writes and "
                    f"reads {list(expected)}.  The manifest states the layout "
                    f"of its own tables, so a mismatch means the files beside "
                    f"it were written by different code than is reading them")
        if settled["window_close"] <= settled["window_open"]:
            raise ValueError(
                f"window_close={settled['window_close']!r} does not follow "
                f"window_open={settled['window_open']!r}")
        object.__setattr__(self, "fields", normalized)

    def value(self, name):
        """Return the frozen value of one manifest field, or raise.

        Access by name rather than by attribute, so that a typo is a
        ``ValueError`` naming the field it looked for rather than an
        ``AttributeError`` a caller might have caught for another reason.
        """
        if isinstance(name, bool) or not isinstance(name, str):
            raise TypeError(f"name must be a string, got "
                            f"{type(name).__name__} {name!r}")
        for key, value in self.fields:
            if key == name:
                return value
        raise ValueError(f"no manifest field {name!r}; the declared fields are "
                         f"{list(MANIFEST_KEYS)}")

    def as_dict(self):
        """Return the manifest as an ordinary dictionary, for reading."""
        return dict(self.fields)

    @property
    def run_label(self) -> str:
        """``"pilot"`` or ``"production"`` — the firewall, in one word."""
        return self.value("run_label")

    @property
    def trials(self) -> int:
        """The authoritative trial count of this run."""
        return self.value("trials")

    @property
    def shadow_trials(self) -> int:
        """How many trials the predeclared no-stop sample covers."""
        return self.value("shadow_trials")

    @property
    def clocks(self) -> int:
        """The physical size of the competing population."""
        return self.value("clocks")


def canonical_manifest_bytes(manifest) -> bytes:
    """Return the bytes a raw run should write for ``manifest``.

    Human-readable on purpose: two-space indentation and one key per line, in
    the declared schema order rather than sorted, so that the file reads as a
    description of a run rather than as a serialization of a dictionary.  A
    person is expected to open this.

    The writer digests exactly these bytes into its close marker; a reader
    digests the bytes it finds and compares.  Reading does **not** require this
    exact spelling — :func:`parse_manifest` accepts any valid JSON object with
    the right keys, because the digest in the marker is over the bytes as they
    were actually written and reformatting them would break it.
    """
    require_exact_type(manifest, RawManifest, "manifest")
    body = {name: (list(value) if isinstance(value, tuple) else value)
            for name, value in manifest.fields}
    return (json.dumps(body, indent=2, sort_keys=False, ensure_ascii=True)
            + "\n").encode("utf-8")


class _DuplicateMember(ValueError):
    """A JSON object names the same member twice.

    Its own class so that the parsers can let it past their
    ``JSONDecodeError`` handlers unchanged.  It is a ``ValueError``, which is
    the convention for "right kind of thing, impossible value", and callers
    have no reason to catch anything narrower.
    """


def _reject_duplicate_members(pairs):
    """``object_pairs_hook``: build the mapping, refusing a repeated name.

    ``json.loads`` hands every JSON object to this hook as the *ordered list of
    pairs it actually parsed*, before any dictionary exists.  That is the only
    place a duplicate is still visible: ``dict(pairs)`` silently keeps the last
    value for a repeated key, so by the time the parser returns,
    ``{"trials": 2, "trials": 999}`` and ``{"trials": 999}`` are the same
    object — and the count the gate enforces and the count a person reads in
    the same file could differ, with nothing to say so.

    Whether the repeated values *agree* is not consulted.  A file that names a
    member twice is malformed however it is read, and accepting the
    harmless-looking case would make the rule about the values rather than
    about the file.
    """
    seen = set()
    for name, _ in pairs:
        if name in seen:
            raise _DuplicateMember(
                f"the document names {name!r} more than once.  A repeated "
                f"member makes its meaning depend on which occurrence a reader "
                f"keeps")
        seen.add(name)
    return dict(pairs)


def _load_object(payload, what: str):
    """Return the JSON object these bytes hold, or raise naming ``what``."""
    if not isinstance(payload, (bytes, bytearray, memoryview)):
        raise TypeError(f"the {what} must be bytes, got "
                        f"{type(payload).__name__}; text is refused rather "
                        f"than encoded, so that writer and reader cannot "
                        f"disagree about the encoding")
    try:
        text = bytes(payload).decode("utf-8")
    except UnicodeDecodeError as exc:
        raise ValueError(f"the {what} is not valid UTF-8: {exc}") from None
    try:
        loaded = json.loads(text, object_pairs_hook=_reject_duplicate_members)
    except _DuplicateMember:
        raise
    except json.JSONDecodeError as exc:
        raise ValueError(f"the {what} is not valid JSON: {exc}") from None
    if not isinstance(loaded, dict):
        raise ValueError(f"a {what} is a JSON object, got "
                         f"{type(loaded).__name__}")
    return loaded


def parse_manifest(payload) -> RawManifest:
    """Return the :class:`RawManifest` these bytes describe, or raise.

    Strict about the schema and relaxed about the formatting: key order,
    indentation and spacing are the writer's business, but the key set must be
    exactly :data:`MANIFEST_KEYS`, every value must satisfy its declared kind,
    and every cross-field rule must hold.

    Refusing an *unknown* key matters as much as requiring a known one.  A
    manifest carrying a field this version does not understand is a manifest
    from a different schema, and reading the keys it recognises while ignoring
    the rest is how a run gets compared against parameters that are not its
    own.

    A *repeated* key is refused too, during parsing rather than after it — see
    :func:`_reject_duplicate_members` for why the distinction is not pedantic.
    """
    return RawManifest(_load_object(payload, "manifest"))


# ---------------------------------------------------------------------------
# The close marker
# ---------------------------------------------------------------------------
#: Every key a close-marker document has, and no others.
CLOSE_MARKER_FIELDS = ("schema_version", "tables", "manifest_digest")

#: The keys of one table's entry inside a marker.
_TABLE_FIELDS = ("name", "rows", "digest")


@dataclass(frozen=True)
class CloseMarker:
    """What a finished raw run states about the tables it wrote.

    ``schema_version``
        Layout of the ledger.  Must equal :data:`LEDGER_SCHEMA_VERSION`; a
        marker for a layout this code does not implement cannot be built, so it
        cannot be half-understood.
    ``tables``
        One ``(name, rows, digest)`` triple for each of :data:`TABLE_NAMES`, in
        that order.  All three, always: a run whose trial ledger closed and
        whose shadow sample did not is not a closed run, and a marker with room
        for only the table somebody remembered would make that unsayable.

        ``rows`` counts *data* rows, not lines, so a table with a header and no
        rows is zero.  The trial ledger and the clock ledger must both have at
        least one; the shadow table may not, only because a shadow sample
        covering zero clocks is arithmetically expressible and the manifest is
        where that is refused.
    ``manifest_digest``
        Digest of the frozen run manifest.  Carried here so a reader can tell
        that the tables in front of it and the manifest in front of it came
        from the same run, rather than assuming it because they sit in the same
        directory.

    The marker is a statement, not a proof.  :func:`require_closed_ledger` is
    what turns it into one, by recomputing every digest from the bytes the
    reader actually has and counting the rows those bytes actually contain.
    """

    schema_version: int
    tables: tuple
    manifest_digest: str

    def __post_init__(self):
        schema_version = require_positive_int(self.schema_version,
                                              "schema_version", minimum=1)
        if schema_version != LEDGER_SCHEMA_VERSION:
            raise ValueError(
                f"schema_version must be {LEDGER_SCHEMA_VERSION}, got "
                f"{schema_version}; a ledger written to a layout this code "
                f"does not implement must be refused, not partially read")
        entries = []
        for position, entry in enumerate(_require_sequence(self.tables,
                                                           "tables")):
            item = tuple(_require_sequence(entry, f"tables[{position}]"))
            if len(item) != 3:
                raise ValueError(f"tables[{position}] must be a "
                                 f"(name, rows, digest) triple, got {item!r}")
            entries.append((_require_text(item[0], f"tables[{position}] name"),
                            _require_int(item[1], f"tables[{position}] rows",
                                         0),
                            _require_digest(item[2],
                                            f"tables[{position}] digest")))
        named = tuple(name for name, _, _ in entries)
        if named != TABLE_NAMES:
            raise ValueError(
                f"a close marker names exactly {list(TABLE_NAMES)} in that "
                f"order, got {list(named)}; a run is closed when all three of "
                f"its tables are, and a marker missing one would certify a "
                f"partial run")
        counts = dict((name, rows) for name, rows, _ in entries)
        for name in ("ledger", "clocks"):
            if counts[name] < 1:
                raise ValueError(
                    f"the {name} table declares {counts[name]} rows; a closed "
                    f"raw run wrote at least one, so zero means the run did "
                    f"not happen rather than that it found nothing")
        manifest_digest = _require_digest(self.manifest_digest,
                                          "manifest_digest")
        object.__setattr__(self, "schema_version", schema_version)
        object.__setattr__(self, "tables", tuple(entries))
        object.__setattr__(self, "manifest_digest", manifest_digest)

    def rows(self, name) -> int:
        """How many data rows this marker declares for one table."""
        return self._entry(name)[1]

    def digest(self, name) -> str:
        """The digest this marker declares for one table's bytes."""
        return self._entry(name)[2]

    def _entry(self, name):
        if isinstance(name, bool) or not isinstance(name, str):
            raise TypeError(f"name must be a string, got "
                            f"{type(name).__name__} {name!r}")
        for entry in self.tables:
            if entry[0] == name:
                return entry
        raise ValueError(f"no table {name!r}; a closed run has "
                         f"{list(TABLE_NAMES)}")


def canonical_marker_bytes(marker) -> bytes:
    """Return the bytes a raw run should write for its close marker.

    The close marker is a *file*, and its presence is half of what it means.
    A run that was interrupted has no marker at all, so a reader that cannot
    find one knows the run is unfinished without having to hash anything.
    """
    require_exact_type(marker, CloseMarker, "marker")
    body = {"schema_version": marker.schema_version,
            "tables": [{"name": name, "rows": rows, "digest": digest}
                       for name, rows, digest in marker.tables],
            "manifest_digest": marker.manifest_digest}
    return (json.dumps(body, indent=2, sort_keys=False, ensure_ascii=True)
            + "\n").encode("utf-8")


def parse_marker(payload) -> CloseMarker:
    """Return the :class:`CloseMarker` these bytes describe, or raise."""
    loaded = _load_object(payload, "close marker")
    present = set(loaded)
    missing = sorted(set(CLOSE_MARKER_FIELDS) - present)
    unknown = sorted(present - set(CLOSE_MARKER_FIELDS))
    if missing or unknown:
        raise ValueError(
            f"a close marker has exactly the keys "
            f"{sorted(CLOSE_MARKER_FIELDS)}; missing {missing}, unexpected "
            f"{unknown}")
    entries = []
    for position, entry in enumerate(_require_sequence(loaded["tables"],
                                                       "tables")):
        if not isinstance(entry, dict):
            raise ValueError(f"tables[{position}] must be a JSON object, got "
                             f"{type(entry).__name__}")
        missing = sorted(set(_TABLE_FIELDS) - set(entry))
        unknown = sorted(set(entry) - set(_TABLE_FIELDS))
        if missing or unknown:
            raise ValueError(f"tables[{position}] has exactly the keys "
                             f"{sorted(_TABLE_FIELDS)}; missing {missing}, "
                             f"unexpected {unknown}")
        entries.append((entry["name"], entry["rows"], entry["digest"]))
    return CloseMarker(loaded["schema_version"], tuple(entries),
                       loaded["manifest_digest"])


# ---------------------------------------------------------------------------
# The gate
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class ClosedLedger:
    """A raw run that has been proved closed, with its rows already parsed.

    What :func:`require_closed_ledger` hands back, and the only object a
    comparison process should ever compute anything from.  Holding the parsed
    rows rather than the bytes is deliberate: a caller that received bytes
    could reparse them with a looser parser, and the point of the gate is that
    there is exactly one parse and it is the one that was checked.
    """

    marker: CloseMarker
    manifest: RawManifest
    ledger: tuple
    clocks: tuple
    shadow: tuple

    def __post_init__(self):
        require_exact_type(self.marker, CloseMarker, "marker")
        require_exact_type(self.manifest, RawManifest, "manifest")
        for name in TABLE_NAMES:
            rows = getattr(self, name)
            if not isinstance(rows, tuple):
                raise TypeError(f"{name} must be a tuple of parsed rows, got "
                                f"{type(rows).__name__}")
            for row in rows:
                if not isinstance(row, dict):
                    raise TypeError(f"every {name} row is a parsed mapping, "
                                    f"got {type(row).__name__}")

    def clock_totals(self):
        """Return the per-clock aggregate, derived from the rows just verified.

        A tuple of dictionaries, one per clock in population order, holding the
        totals a summary wants: how many trials the clock was eligible in and
        reached the band in, how many it won and how many of those it shared,
        its entry and reset counts, and its total exposure and opportunity.

        **Derived, never read.**  Version 2 wrote these totals to disk, and an
        independent review pointed out that nothing could check them: a
        fabricated aggregate passed every byte, schema and count rule the gate
        had.  Version 3 keeps the per-trial rows instead and computes the
        summary here, from rows :func:`require_closed_ledger` has already
        reconciled against the manifest and against each other.  There is then
        one place each number lives, and it is the place that was checked.
        """
        totals = [{"clock": clock, "detuning": None, "trials_eligible": 0,
                   "trials_inside": 0, "wins": 0, "co_completions": 0,
                   "eligibility_entries": 0, "band_entries": 0,
                   "dwell_resets": 0, "exposure_time": [], "eligible_time": []}
                  for clock in range(self.manifest.clocks)]
        for row in self.clocks:
            entry = totals[row["clock"]]
            entry["detuning"] = row["detuning"]
            entry["trials_eligible"] += 1 if row["eligible_endpoints"] else 0
            entry["trials_inside"] += 1 if row["band_entries"] else 0
            entry["wins"] += 1 if row["won"] else 0
            entry["co_completions"] += 1 if row["co_completed"] else 0
            for field in ("eligibility_entries", "band_entries",
                          "dwell_resets"):
                entry[field] += row[field]
            entry["exposure_time"].append(row["exposure_time"])
            entry["eligible_time"].append(row["eligible_time"])
        for entry in totals:
            for field in ("exposure_time", "eligible_time"):
                entry[field] = math.fsum(entry[field])
        return tuple(totals)

    def rows(self, name):
        """Return one table's parsed rows by name."""
        if isinstance(name, bool) or not isinstance(name, str):
            raise TypeError(f"name must be a string, got "
                            f"{type(name).__name__} {name!r}")
        if name not in TABLE_NAMES:
            raise ValueError(f"no table {name!r}; a closed run has "
                             f"{list(TABLE_NAMES)}")
        return getattr(self, name)


def encode_table(columns, rows) -> bytes:
    """Return the exact bytes a table with these rows was written as.

    The inverse of :func:`parse_table`, byte for byte: the header line, then
    one canonical line per row, each terminated by a single newline — which is
    precisely what ``raw_runner._Table`` folded into the digest as it wrote.

    It exists so that a *parsed* ledger can be re-derived into bytes and
    checked against the marker that came with it.  See
    :func:`require_unmutated_ledger`.
    """
    lines = [header_row(columns)]
    lines += [encode_row(columns, row) for row in rows]
    return "".join(line + "\n" for line in lines).encode("utf-8")


def require_unmutated_ledger(ledger):
    """Return the **freshly reconciled** snapshot of ``ledger``, or raise.

    Not a predicate.  The value it returns is the value a caller must go on to
    use, and every door in the analysis layer binds it.

    Why it is not a predicate
    -------------------------
    A :class:`ClosedLedger` is a frozen dataclass, but "frozen" freezes the
    *references*: its rows are ordinary dictionaries.  An independent review
    first mutated one in place, and an earlier version of this function caught
    that by re-encoding the rows and comparing the digests with the marker the
    object carried.

    The review then pointed out what that could not catch.  If a caller edits
    the rows, re-encodes all three tables with :func:`encode_table`, rebuilds a
    :class:`CloseMarker` over the new digests and constructs a new
    ``ClosedLedger``, then the snapshot agrees with *its own* marker perfectly.
    Four all-public variants were accepted::

        trial 0 rewritten dwell_failed -> committed, rehashed
        marker declaring 7 ledger rows beside the genuine 6
        manifest detunings shifted +0.125 beside unchanged clock rows
        ledger rows reversed and rehashed

    The three observables then answered inconsistently — a tally reporting five
    commitments beside a winner table holding four — while
    :func:`require_closed_ledger` refused all four outright.  Canonical
    agreement with a caller-supplied marker is simply weaker than the
    authoritative gate's row-count, identity, ordering and cross-table
    reconciliation.

    What it does now
    ----------------
    Reconstructs the canonical manifest and table bytes from the snapshot, and
    puts those bytes and the snapshot's marker through the whole of
    :func:`require_closed_ledger` — the same eight checks a run read from disk
    passes.  The object handed back is the one *that* gate parsed, so there is
    no check-then-use gap: a caller cannot validate one snapshot and keep
    computing from another.

    There is no recursion here.  :func:`require_closed_ledger` reads bytes and
    knows nothing about this function; this one builds bytes and calls it once.

    This is not a substitute for re-reading the run.  The decision factories in
    ``compare.py`` reopen the named run through the descriptor-scoped reader
    before they build anything, because a *file* can change under an object
    that is internally perfect.  This is for the doors that hold a ledger and
    no run name — the observables — and for the design builder's own belt and
    braces.
    """
    require_exact_type(ledger, ClosedLedger, "ledger")
    payloads = {name: encode_table(_TABLE_COLUMNS[name], ledger.rows(name))
                for name in TABLE_NAMES}
    return require_closed_ledger(ledger.marker,
                                 canonical_manifest_bytes(ledger.manifest),
                                 payloads["ledger"], payloads["clocks"],
                                 payloads["shadow"])


def require_closed_ledger(marker, manifest, ledger, clocks,
                          shadow) -> ClosedLedger:
    """Return a fully verified :class:`ClosedLedger`, or raise.

    The gate a comparison process passes before it reads a row.  It takes
    *bytes* rather than paths because a path can be reopened between the check
    and the read; the caller reads once, verifies the bytes it read, and
    receives the parse of those same bytes.

    Eight things are checked, in this order, and each one is here because
    skipping it lets a specific wrong number look right.

    1. **The marker is exactly a** :class:`CloseMarker`, not an instance of
       one: a subclass could override anything.
    2. **The marker is rebuilt here**, which re-runs every field check.  A
       marker is frozen against ordinary assignment and not against
       ``object.__setattr__``, and a marker edited after construction would
       otherwise sail through.  The value returned holds the rebuilt marker,
       not the one passed in, so a caller cannot keep a reference to something
       the gate never validated.
    3. **The manifest bytes are verified** against ``manifest_digest``, and are
       verified *before* they are parsed, so everything read out of them is
       read out of bytes this marker has already vouched for.  Without this,
       nothing binds the two together and a ledger can be compared against a
       prediction computed from a different run's parameters.
    4. **Every table's bytes are verified** against its own declared digest.
       Truncation, extension, a hand edit and a single flipped byte all fail
       here.
    5. **Every table is parsed**, header and all, by :func:`parse_table`.  A
       renamed column, a reordered column, a row of the wrong width, a number
       rewritten by a spreadsheet or a category outside the declared four fails
       here.
    6. **The counts are reconciled three ways**: what the manifest froze, what
       the marker declares, and how many rows are actually present.  This is
       what replaces ticket 01's manifest-level-only check, which could
       establish only that a marker agreed with a manifest about a number
       neither had counted.
    7. **The row identities are exactly the declared ones.**  Trials
       ``0 .. trials-1``, each once; the authoritative clock table exactly the
       full trial-by-clock grid, each pair once and in order; shadow rows
       exactly the predeclared shadow trial-by-clock grid.  A ledger
       holding one trial twice and another not at all has the right row count
       and is not the run's ledger, and this is the check that says so.
    8. **Every row is internally possible, and the tables agree with each
       other.**  A cell can be the right kind and still be a lie: one word
       changed from ``dwell_failed`` to ``committed`` passes the digests — a
       forger recomputes them — the schema and every count, and turns an
       unresolved trial into a commitment.  :func:`_require_consistent_rows`
       restates the identities the writer already enforced, over the parsed
       table, so a reader establishes them for the bytes in front of it.

    ``TypeError`` for "not that kind of thing" — no marker, a marker subclass,
    text where bytes belong.  ``ValueError`` for "right kind, impossible or
    disagreeing value".

    None of this is authentication.  The digests establish that the bytes a
    reader holds are the bytes a writer finished with; anyone able to rewrite a
    table can recompute them.  It is protection against truncation,
    interleaving, mismatched pairs and accidental edit, which is what actually
    goes wrong.  It is not a signature.
    """
    require_exact_type(marker, CloseMarker, "marker")

    # (2) rebuilding is the validation.
    canonical = CloseMarker(marker.schema_version, marker.tables,
                            marker.manifest_digest)

    # (3) the manifest, verified before it is parsed.
    manifest_actual = digest_bytes(manifest)
    if manifest_actual != canonical.manifest_digest:
        raise ValueError(
            f"manifest digest mismatch: the close marker states "
            f"{canonical.manifest_digest}, the manifest supplied hashes to "
            f"{manifest_actual}.  These tables and this manifest are not from "
            f"the same run, so the parameters a prediction would be computed "
            f"from are not the parameters the events were generated under")
    frozen = parse_manifest(manifest)

    payloads = {"ledger": ledger, "clocks": clocks, "shadow": shadow}
    parsed = {}
    for name in TABLE_NAMES:
        payload = payloads[name]
        # (4) the bytes.
        actual = digest_bytes(payload)
        if actual != canonical.digest(name):
            raise ValueError(
                f"{name} digest mismatch: the close marker states "
                f"{canonical.digest(name)}, the bytes read hash to {actual}.  "
                f"This table has been truncated, appended to, or edited since "
                f"the raw run closed it, and must not be read")
        # (5) the rows.
        parsed[name] = parse_table(payload, _TABLE_COLUMNS[name])

    # (6) the three counts.
    expected = {"ledger": frozen.trials,
                "clocks": frozen.trials * frozen.clocks,
                "shadow": frozen.shadow_trials * frozen.clocks}
    for name in TABLE_NAMES:
        declared, present, promised = (canonical.rows(name),
                                       len(parsed[name]), expected[name])
        if not declared == present == promised:
            raise ValueError(
                f"the {name} table has {present} rows; the close marker "
                f"declares {declared} and the frozen manifest implies "
                f"{promised}.  A marker that states its own row count and is "
                f"believed proves nothing, so all three must agree — the "
                f"table is from a different run, a partial run, or a rerun")

    # (7) the row identities.
    _require_identities(tuple(row["trial"] for row in parsed["ledger"]),
                        tuple(range(frozen.trials)), "ledger", "trial")
    _require_identities(
        tuple((row["trial"], row["clock"]) for row in parsed["clocks"]),
        tuple((trial, clock)
              for trial in range(frozen.trials)
              for clock in range(frozen.clocks)),
        "clocks", "(trial, clock)")
    _require_identities(
        tuple((row["trial"], row["clock"]) for row in parsed["shadow"]),
        tuple((trial, clock)
              for trial in range(frozen.shadow_trials)
              for clock in range(frozen.clocks)),
        "shadow", "(trial, clock)")

    # (8) the rows have to be internally possible, not merely well typed.
    _require_consistent_rows(parsed, frozen)

    return ClosedLedger(canonical, frozen, parsed["ledger"], parsed["clocks"],
                        parsed["shadow"])


def _require_consistent_rows(parsed, frozen):
    """Derive every channel-level number from the authoritative rows, and
    compare.

    A cell can be the right *kind* and still be a lie.  ``category`` is one of
    the declared four and ``commit_time`` is a number or absent, and a row
    saying ``committed`` with no commitment time passes the digests — a forger
    recomputes them — the schema, and every count.  An independent review went
    further, twice.

    **Round one** got three records accepted, each with every hash and marker
    field recomputed after the edit: ``never_eligible`` rewritten as
    ``lock_failed`` while zero clocks were ever eligible or inside; a manifest
    detuning of ``0.0`` bound to clock and shadow detunings of ``0.75`` and
    ``-0.75``; and one eligible trial with ``999`` units of exposure in a
    two-trial run whose whole observation window is ``0.4``.  Those close by
    deriving the category from the clock rows, binding every indexed detuning
    to the manifest, and recomputing censored exposure.

    **Round two** got four more accepted, and they are all about *time*:

    1. a no-stop shadow row's exposure changed from ``4.0`` to ``2.0`` — the
       row spans the whole window and nothing derived what that costs;
    2. a censored loser's ``first_eligible_at`` moved to the manifest close,
       long after the trial it belongs to had stopped;
    3. a same-key shadow row's ``first_eligible_at`` moved away from the value
       the censored walk had already observed on the same physical path; and
    4. a shadow row's ``first_inside_at`` moved before its own
       ``first_eligible_at`` — a clock inside a band it had not yet entered.

    All four passed because the gate bounded times against the *whole* manifest
    window rather than against each row's own observation interval, never
    derived shadow exposure at all, and compared only counts across the
    shadow prefix.

    So each row now carries a derived **observation end**: its own commitment
    if it made one, and otherwise the channel stop for a censored row or the
    window close for a no-stop row.  Exposure is that end minus the window's
    opening, exactly, through the same expression the writer used.  Every
    recorded instant must fall inside that row's own interval.  Event times
    must be in causal order.  And where the censored and no-stop walks of one
    ``(trial, clock)`` cover the same interval, the two rows must be identical
    — because they are the same walk.

    Nothing here re-derives physics.  Each rule is an identity the writer
    already enforced at the object that produced the row, restated over the
    parsed table so that a reader establishes it for the bytes in front of it
    rather than inheriting it from a process it never saw.
    """
    clocks = frozen.clocks
    detunings = frozen.value("detunings")
    opened = frozen.value("window_open")
    closed = frozen.value("window_close")

    # The stop each trial claims.  Used to give every censored row its own
    # observation end, and independently re-derived from the clock rows below,
    # so a forger who moves a stop to make its clock rows fit is caught by the
    # category and winner derivation rather than believed here.
    stops = {row["trial"]: row["stop_time"] for row in parsed["ledger"]}

    by_trial = {}
    for row in parsed["clocks"]:
        trial = row["trial"]
        if trial not in stops:
            raise ValueError(f"clocks trial {trial} has no trial row")
        ended = (row["committed_at"] if row["committed_at"] is not None
                 else stops[trial])
        _require_possible_clock_row(row, detunings, opened, closed, ended,
                                    "clocks", won_known=True)
        by_trial.setdefault(trial, []).append(row)
    shadow_rows = {}
    for row in parsed["shadow"]:
        # A no-stop walk runs to the window's close unless the clock's own
        # dwell completed first.  Deriving that is what refuses the review's
        # first round-two variant.
        ended = (row["committed_at"] if row["committed_at"] is not None
                 else closed)
        _require_possible_clock_row(row, detunings, opened, closed, ended,
                                    "shadow", won_known=False)
        shadow_rows[(row["trial"], row["clock"])] = (row, ended)

    for row in parsed["ledger"]:
        trial = row["trial"]
        where = f"ledger trial {trial}"
        history = by_trial.get(trial)
        if history is None or len(history) != clocks:
            raise ValueError(
                f"{where} has {0 if history is None else len(history)} "
                f"authoritative clock rows for a population of {clocks}")

        winners = tuple(item["clock"] for item in history if item["won"])
        committed = bool(winners)
        derived = _derived_category(history)
        if row["category"] != derived:
            raise ValueError(
                f"{where} is recorded as {row['category']!r}; its own "
                f"{clocks} clock rows say {derived!r}.  The channel category "
                f"is derived from the clocks and cannot be asserted beside "
                f"them")
        if committed != (row["commit_time"] is not None):
            raise ValueError(f"{where} is {row['category']!r} and its commit "
                             f"time is {row['commit_time']!r}")
        if row["censored"] != committed:
            raise ValueError(f"{where} is {row['category']!r} and censored="
                             f"{row['censored']}; observation is cut short "
                             f"exactly when a clock completes")
        if committed and row["commit_time"] != row["stop_time"]:
            raise ValueError(f"{where} committed at {row['commit_time']!r} and "
                             f"stopped at {row['stop_time']!r}")
        expected_stop = row["commit_time"] if committed else closed
        if row["stop_time"] != expected_stop:
            raise ValueError(
                f"{where} stops at {row['stop_time']!r}; a committed trial "
                f"stops at its own commitment and an unresolved one at the "
                f"manifest's window close {closed!r}")
        if not opened < row["stop_time"] <= closed:
            raise ValueError(
                f"{where} stops at {row['stop_time']!r}, outside the "
                f"manifest's observation window ({opened!r}, {closed!r}]")
        if tuple(row["winners"]) != winners:
            raise ValueError(
                f"{where} names winners {list(row['winners'])}; the clocks "
                f"whose own authoritative rows say they completed are "
                f"{list(winners)}")
        if row["co_completion"] != len(winners):
            raise ValueError(f"{where} declares {row['co_completion']} "
                             f"co-completers and {len(winners)} clocks won")
        for item in history:
            if item["won"] and item["committed_at"] != row["commit_time"]:
                raise ValueError(
                    f"{where} clock {item['clock']} won but completed at "
                    f"{item['committed_at']!r} rather than at the stop "
                    f"{row['commit_time']!r}")
            if item["co_completed"] != (item["won"] and len(winners) > 1):
                raise ValueError(
                    f"{where} clock {item['clock']} is marked "
                    f"co_completed={item['co_completed']} with "
                    f"{len(winners)} winners in the trial")

        won_rows = [item for item in history if item["won"]]
        for field, source in (("winner_detunings", "detuning"),
                              ("winner_phases", "final_phase"),
                              ("winner_band_entry", "first_inside_at"),
                              ("winner_first_eligible", "first_eligible_at"),
                              ("winner_last_reset", "last_reset_at")):
            declared = tuple(row[field])
            if len(declared) != len(won_rows):
                raise ValueError(
                    f"{where} names {len(won_rows)} winners and "
                    f"{len(declared)} values in {field}; every per-winner "
                    f"column is indexed by the same winner list")
            for item, value in zip(won_rows, declared):
                # ``last_reset_at`` is legitimately absent for a winner that
                # completed on its first attempt, so absence is compared
                # rather than refused.
                if value != item[source]:
                    raise ValueError(
                        f"{where} gives clock {item['clock']} the {field} "
                        f"value {value!r}; its own row says {item[source]!r}")

        for field, derived_value in (
                ("clocks_ever_eligible",
                 sum(1 for i in history if i["eligible_endpoints"])),
                ("clocks_ever_inside",
                 sum(1 for i in history if i["band_entries"])),
                ("band_entries", sum(i["band_entries"] for i in history)),
                ("dwell_resets", sum(i["dwell_resets"] for i in history))):
            if row[field] != derived_value:
                raise ValueError(
                    f"{where} reports {field}={row[field]}; its "
                    f"authoritative clock rows give {derived_value}")
        opportunity = math.fsum(i["eligible_time"] for i in history)
        if row["eligible_time"] != opportunity:
            raise ValueError(
                f"{where} reports {row['eligible_time']!r} of total "
                f"opportunity; its clock rows sum to {opportunity!r}")

    # The shadow sample is the same clocks on the same keys without the stop.
    for (trial, clock), (free, shadow_end) in shadow_rows.items():
        observed = next((item for item in by_trial.get(trial, ())
                         if item["clock"] == clock), None)
        if observed is None:
            raise ValueError(f"the shadow sample covers trial {trial} clock "
                             f"{clock}, which has no censored row")
        censored_end = (observed["committed_at"]
                        if observed["committed_at"] is not None
                        else stops[trial])
        _require_shadow_prefix(observed, free, censored_end, shadow_end,
                               trial, clock)


def _require_shadow_prefix(observed, free, censored_end, shadow_end, trial,
                           clock):
    """Require the no-stop row to be the same physical path, continued.

    The shadow walk is the censored walk with the stop switched off, on the
    same dataset, trial, clock and mesh — therefore the same prepared phase and
    the same Brownian leaves.  It is not a second realization and it must not
    read like one.

    Two regimes, decided by the observation ends:

    * **The same end.**  Either the clock's own dwell completed, in which case
      both walks stopped there, or the trial ran to the window's close, in
      which case neither walk stopped early.  The two walks are then the same
      walk, and every shared field must be **identical**.
    * **A shorter censored end.**  The shadow saw everything the censored walk
      saw and then continued.  Counts may only grow; the first eligibility and
      the first band entry, once observed, may not move; a reset already
      observed may not move, and any additional one must fall after the
      censored walk stopped looking.

    An independent review changed a same-key shadow ``first_eligible_at`` away
    from the value the censored row had already recorded, and the previous gate
    — which compared only counts, opportunity and detuning — accepted it.
    """
    where = f"trial {trial} clock {clock}"
    if free["detuning"] != observed["detuning"]:
        raise ValueError(f"{where} has two mismatches, one in each table")

    # A commitment already observed under the stop is the same event on the
    # same path, so it recurs at the identical instant without the stop.  This
    # is checked *before* the ends are compared, and that ordering is the
    # point: an independent review moved a winner's shadow commitment from
    # -0.5 to -0.4 and derived the new exposure to match, which moved the
    # shadow's own observation end and stepped out of the equal-ends branch
    # where the rule used to live.
    if observed["committed_at"] is not None \
            and free["committed_at"] != observed["committed_at"]:
        raise ValueError(
            f"{where} completed its dwell at {observed['committed_at']!r} "
            f"under the stop and at {free['committed_at']!r} without it.  Both "
            f"walks consume the same keys on the same path and both stop at "
            f"that clock's own completion, so the instant cannot move")

    if censored_end > shadow_end:
        raise ValueError(
            f"{where} was observed until {censored_end!r} under the stop and "
            f"only until {shadow_end!r} without it; switching the stop off "
            f"cannot shorten a walk")

    shared = [name for name, _ in SHADOW_COLUMNS]
    if censored_end == shadow_end:
        differing = [name for name in shared if observed[name] != free[name]]
        if differing:
            raise ValueError(
                f"{where} was observed over the same interval "
                f"({censored_end!r}) with and without the stop, so the two "
                f"rows are the same walk, but they differ in {differing}")
        return

    for field in ("eligibility_entries", "band_entries", "dwell_resets",
                  "endpoints", "eligible_endpoints", "inside_endpoints"):
        if observed[field] > free[field]:
            raise ValueError(
                f"{where}: the censored {field} ({observed[field]}) exceeds "
                f"the uncensored one ({free[field]}).  The stop can only "
                f"remove observation, so the censored walk is a prefix of the "
                f"free one")
    for field in ("eligible_time", "exposure_time"):
        if observed[field] > free[field]:
            raise ValueError(f"{where}: censored {field} exceeds uncensored "
                             f"{field}")

    # A strictly longer observation is a longer walk over the same mesh, so it
    # must have seen at least one more endpoint.  Counts and times are then
    # *jointly* monotone rather than independently so: opportunity can only
    # grow by adding eligible endpoints, and adding an eligible endpoint adds
    # the strictly positive duration of the interval it closes.  An independent
    # review moved a shadow endpoint count down to its censored value, changed
    # opportunity while the eligible count stayed equal, and added seven
    # eligible endpoints worth nothing -- all three passed a gate whose counts
    # and times were monotone but unrelated.
    if not free["endpoints"] > observed["endpoints"]:
        raise ValueError(
            f"{where} was observed until {censored_end!r} under the stop and "
            f"until {shadow_end!r} without it, yet reports {free['endpoints']} "
            f"endpoints against {observed['endpoints']}.  A longer walk over "
            f"the same mesh sees at least one more endpoint")
    if free["eligible_endpoints"] == observed["eligible_endpoints"]:
        if free["eligible_time"] != observed["eligible_time"]:
            raise ValueError(
                f"{where} was eligible at the same "
                f"{free['eligible_endpoints']} endpoints in both walks and "
                f"reports {observed['eligible_time']!r} of opportunity under "
                f"the stop against {free['eligible_time']!r} without it.  The "
                f"same endpoints close the same intervals")
    elif not free["eligible_time"] > observed["eligible_time"]:
        raise ValueError(
            f"{where} gains "
            f"{free['eligible_endpoints'] - observed['eligible_endpoints']} "
            f"eligible endpoints without gaining opportunity: "
            f"{observed['eligible_time']!r} against "
            f"{free['eligible_time']!r}.  Every added eligible endpoint closes "
            f"an interval of strictly positive duration")
    if free["band_entries"] == observed["band_entries"] \
            and free["inside_endpoints"] < observed["inside_endpoints"]:
        raise ValueError(
            f"{where} started the same {free['band_entries']} dwells in both "
            f"walks and was inside the band at fewer endpoints without the "
            f"stop")

    for field in ("first_eligible_at", "first_inside_at"):
        seen, later = observed[field], free[field]
        if seen is not None:
            if later != seen:
                raise ValueError(
                    f"{where}: the censored walk recorded {field}={seen!r} and "
                    f"the no-stop walk of the same keys says {later!r}.  A "
                    f"first event already observed on a path cannot move when "
                    f"the path is walked further")
        elif later is not None and not later > censored_end:
            raise ValueError(
                f"{where}: the no-stop walk reports {field}={later!r}, at or "
                f"before the censored walk stopped looking at "
                f"{censored_end!r}, yet the censored row records none")

    if free["dwell_resets"] == observed["dwell_resets"]:
        if free["last_reset_at"] != observed["last_reset_at"]:
            raise ValueError(
                f"{where}: both walks count {free['dwell_resets']} resets and "
                f"they disagree about when the last one was")
    else:
        if free["last_reset_at"] is None:
            raise ValueError(f"{where}: the no-stop walk counts "
                             f"{free['dwell_resets']} resets and names none")
        if not free["last_reset_at"] > censored_end:
            raise ValueError(
                f"{where}: the no-stop walk has more resets than the censored "
                f"one but its last is at {free['last_reset_at']!r}, at or "
                f"before the censored walk stopped looking at "
                f"{censored_end!r}")

    if free["committed_at"] is not None \
            and observed["committed_at"] is None \
            and not free["committed_at"] > censored_end:  # noqa: E501
        raise ValueError(
            f"{where}: the no-stop walk commits at {free['committed_at']!r}, "
            f"at or before the censored walk stopped looking at "
            f"{censored_end!r}, yet the censored row records no commitment")


def _derived_category(history):
    """The one category a trial's authoritative clock rows put it in."""
    if any(item["won"] for item in history):
        return "committed"
    if not any(item["eligible_endpoints"] for item in history):
        return "never_eligible"
    if not any(item["band_entries"] for item in history):
        return "lock_failed"
    return "dwell_failed"


def _require_possible_clock_row(row, detunings, opened, closed, ended,
                                table: str, won_known: bool):
    """Require one authoritative clock row to be a state a walk can reach.

    Every rule here is one the in-memory ``ClockExposure`` already enforces,
    restated over a parsed row, plus the temporal rules that only a *reader*
    can state — because only a reader has the manifest window and the trial's
    stop beside the row.

    ``ended`` is this row's own observation end: its commitment if it made one,
    the channel stop if it is a censored row, and the window close if it is a
    no-stop row.  Exposure is derived from it and every recorded instant is
    bounded by it.  An independent review moved a censored loser's first
    eligibility to the manifest close — long after the trial it belongs to had
    stopped — and a gate that bounded times by the whole window accepted it.
    """
    clock = row["clock"]
    where = f"{table} trial {row['trial']} clock {clock}"
    if clock >= len(detunings):
        raise ValueError(f"{where} is outside the declared population of "
                         f"{len(detunings)}")
    if row["detuning"] != detunings[clock]:
        raise ValueError(
            f"{where} has the mismatch {row['detuning']!r}; the manifest's "
            f"frozen grid gives clock {clock} the mismatch "
            f"{detunings[clock]!r}.  A clock identifier means the clock at "
            f"that node of the declared grid, and a row naming a different "
            f"mismatch describes a population this manifest never froze")

    committed = row["category"] == "committed"
    if committed != (row["committed_at"] is not None):
        raise ValueError(f"{where} is {row['category']!r} and its commitment "
                         f"time is {row['committed_at']!r}")
    if won_known and row["won"] != committed:
        raise ValueError(f"{where} is {row['category']!r} and won="
                         f"{row['won']}; a clock wins exactly when its own "
                         f"dwell completed within its exposure")
    if not won_known and "won" in row:
        raise ValueError(f"{where} carries a win flag; a no-stop walk has no "
                         f"channel stop and therefore nothing to win")

    # This row's own observation interval, and every instant inside it.
    if not opened < ended <= closed:
        raise ValueError(
            f"{where} was observed until {ended!r}, outside the manifest's "
            f"window ({opened!r}, {closed!r}]")
    expected_exposure = max(ended - opened, 0.0)
    if row["exposure_time"] != expected_exposure:
        raise ValueError(
            f"{where} claims {row['exposure_time']!r} of exposure; this row "
            f"was observed from {opened!r} until {ended!r}, which is "
            f"{expected_exposure!r}.  Exposure is derived from the window and "
            f"this row's own observation end, never read")
    for field in ("committed_at", "first_eligible_at", "first_inside_at",
                  "last_reset_at"):
        moment = row[field]
        if moment is None:
            continue
        if not opened < moment <= ended:
            raise ValueError(
                f"{where} reports {field}={moment!r}, outside its own "
                f"observation interval ({opened!r}, {ended!r}].  A row's "
                f"events are bounded by when that row stopped being watched, "
                f"not by when the run ended")
    if row["eligible_time"] > row["exposure_time"]:
        raise ValueError(f"{where} spent {row['eligible_time']!r} inside its "
                         f"tongue out of {row['exposure_time']!r} observed; "
                         f"opportunity is part of exposure")
    # Opportunity is the summed duration of the intervals whose endpoints were
    # eligible, and every counted interval has strictly positive duration -- a
    # zero-duration coordinate handoff is not an endpoint.  So the count and
    # the time are zero together and positive together, and an independent
    # review got both halves of that accepted: 0.1 units of opportunity with no
    # eligible endpoint, and two eligible endpoints with none.
    if (row["eligible_endpoints"] == 0) != (row["eligible_time"] == 0.0):
        raise ValueError(
            f"{where} reports {row['eligible_endpoints']} eligible endpoints "
            f"and {row['eligible_time']!r} of opportunity.  Every counted "
            f"endpoint closes an interval of strictly positive duration, so "
            f"the count and the time are zero together and positive together")

    ends = row["endpoints"]
    eligible = row["eligible_endpoints"]
    inside = row["inside_endpoints"]
    entries = row["eligibility_entries"]
    band = row["band_entries"]
    resets = row["dwell_resets"]
    if eligible > ends:
        raise ValueError(f"{where} has {eligible} eligible endpoints out of "
                         f"{ends}")
    if inside > eligible:
        raise ValueError(f"{where} was inside the band at {inside} endpoints "
                         f"and eligible at {eligible}")
    if entries > eligible:
        raise ValueError(f"{where} entered its tongue {entries} times with "
                         f"{eligible} eligible endpoints")
    if (entries == 0) != (eligible == 0):
        raise ValueError(f"{where} has {entries} tongue entries and "
                         f"{eligible} eligible endpoints")
    if band > inside:
        raise ValueError(f"{where} started {band} dwells and was inside the "
                         f"band at {inside} endpoints")
    if (inside == 0) != (band == 0):
        raise ValueError(f"{where} was inside the band at {inside} endpoints "
                         f"and started {band} dwells; a clock is inside "
                         f"exactly when some dwell of its own is running")
    if resets > band:
        raise ValueError(f"{where} counts {resets} resets against {band} "
                         f"entries")
    live = band - resets
    if live not in (0, 1):
        raise ValueError(f"{where} has {live} live dwells; one clock has at "
                         f"most one")
    if committed and live != 1:
        raise ValueError(f"{where} committed without ending on a live dwell")
    if row["category"] == "never_eligible" and eligible:
        raise ValueError(f"{where} is 'never_eligible' with {eligible} "
                         f"eligible endpoints")
    if row["category"] != "never_eligible" and eligible < 1:
        raise ValueError(f"{where} is {row['category']!r} with no eligible "
                         f"endpoint")
    if row["category"] == "lock_failed" and band:
        raise ValueError(f"{where} is 'lock_failed' with {band} band entries")
    if row["category"] in ("committed", "dwell_failed") and band < 1:
        raise ValueError(f"{where} is {row['category']!r} with no band entry")
    if (row["first_eligible_at"] is None) != (eligible == 0):
        raise ValueError(f"{where} disagrees with itself about whether it was "
                         f"ever eligible")
    if (row["first_inside_at"] is None) != (band == 0):
        raise ValueError(f"{where} disagrees with itself about whether it ever "
                         f"reached the band")
    if (row["last_reset_at"] is None) != (resets == 0):
        raise ValueError(f"{where} disagrees with itself about whether a dwell "
                         f"was ever reset")

    # Causal order.  A clock is inside a band it has entered, resets a dwell it
    # has started, and completes a dwell it has not just reset.
    first_eligible = row["first_eligible_at"]
    first_inside = row["first_inside_at"]
    last_reset = row["last_reset_at"]
    moment = row["committed_at"]
    if first_inside is not None:
        if first_eligible is None:
            raise ValueError(f"{where} reached the band without ever being "
                             f"eligible")
        if first_eligible > first_inside:
            raise ValueError(
                f"{where} was first inside the band at {first_inside!r} and "
                f"first eligible only at {first_eligible!r}; a clock is inside "
                f"a band it has entered")
    if last_reset is not None:
        if first_inside is None:
            raise ValueError(f"{where} reset a dwell it never started")
        if not last_reset > first_inside:
            raise ValueError(
                f"{where} reset a dwell at {last_reset!r}, at or before the "
                f"first one started at {first_inside!r}; a reset ends a dwell "
                f"at a later endpoint than the one that began it")
    if moment is not None:
        if first_inside is None:
            raise ValueError(f"{where} completed a dwell it never started")
        if not moment > first_inside:
            raise ValueError(
                f"{where} committed at {moment!r}, at or before its first band "
                f"entry at {first_inside!r}; a fixed positive dwell takes at "
                f"least one endpoint to elapse")
        if last_reset is not None and not moment > last_reset:
            raise ValueError(
                f"{where} committed at {moment!r}, at or before its last reset "
                f"at {last_reset!r}; the dwell that completed began after that "
                f"reset cleared the one before it")


def _require_identities(present, expected, table: str, what: str):
    """Require a table's row identities to be exactly the declared ones.

    Order included.  A ledger whose trials are complete but shuffled was
    written by something other than the run that declared them, and the cheap
    way for two runs' output to end up in one file is for their rows to
    interleave — which produces exactly a complete-but-shuffled table when the
    two runs cover the same trials.
    """
    if present == expected:
        return
    missing = sorted(set(expected) - set(present), key=repr)
    unexpected = sorted(set(present) - set(expected), key=repr)
    repeated = sorted({item for item in present
                       if present.count(item) > 1}, key=repr)
    if missing or unexpected or repeated:
        raise ValueError(
            f"the {table} table does not hold exactly the declared {what} "
            f"identities: missing {missing[:8]}, unexpected {unexpected[:8]}, "
            f"repeated {repeated[:8]}.  A table with the right number of rows "
            f"and the wrong rows in it is what two runs' output interleaved "
            f"into one file looks like")
    raise ValueError(
        f"the {table} table holds every declared {what} exactly once but not "
        f"in the declared order; rows are written in identity order, so a "
        f"reordered table was not written by the run that declared it")
