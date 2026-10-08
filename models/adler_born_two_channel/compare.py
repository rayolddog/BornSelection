#!/usr/bin/env python3
"""
compare.py — the separate comparison process: read a closed run, never write
one, and refuse to turn a fitted number into a physical claim.

This is the only module in the package allowed to hold a closed raw ledger and
the analytic prediction at the same time.  That is the whole reason the raw
side was built the way it was: ``verify.py`` walks the import graph from
``raw_runner.py`` and requires that it reach neither ``analytic`` nor this
file, and starts a fresh interpreter to confirm that importing the raw entry
point loads none of ``observables``, ``analysis``, ``compare``,
``killed_diffusion`` or ``moving_band_audit``.  A generator that could see the
prediction could not test it.

Read-only, and checked both ways
--------------------------------
Everything here goes through :func:`~.raw_runner.open_raw_run`, which reads the
five files once behind pinned directory descriptors, hands the bytes to
:func:`~.raw_ledger.require_closed_ledger`, and returns the parse those exact
bytes produced.  Nothing in this module opens a file for writing, and there is
no path here that could: the closed run is a tuple of parsed rows by the time
this file sees it.

That is an argument, and an argument is not evidence, so
:class:`ReadOnlyRun` also carries the run's :class:`LedgerFingerprint` — the
manifest digest and all three table digests, as the close marker declares them
and as the gate has already recomputed them from the bytes.
:meth:`ReadOnlyRun.reverify` re-opens the run from disk and requires an
identical fingerprint, so a comparison can be bracketed by a before-and-after
hash rather than by a promise.

The four comparators are four different quantities
--------------------------------------------------
They are separately named, separately typed and carry their own units, because
three of the four are numerically similar and mean entirely different things:

======================  ==============================  ==================
kind                    what it is                      finite population?
======================  ==============================  ==================
``eligible_count``      how many of the exact production yes
                        clocks are inside the tongue
``bare_rate_sum``       the sum of their local           yes
                        relaxation rates
``continuum_flux``      that sum times the grid          **no**
                        spacing: a midpoint quadrature
``constant_hazard``     a declared constant rate         n/a — a control
======================  ==============================  ==================

``bare_rate_sum`` is the **primary** comparator, because the experiment has a
finite population of actual competitors and that is what the sum is over.
``continuum_flux`` is the same numbers multiplied by a spacing; it converges to
``(pi/2) K^2`` as the grid refines and is a *continuum control*, never the
finite-population predictor.  They differ by a population-density factor and
answer different questions the moment the clock count changes, so
:class:`ComparatorCurve` refuses to be built with the wrong unit or the wrong
finite-population flag for its kind, and :class:`ComparatorSet` refuses to hold
two curves of one kind.  The conflation the plan warns about is a construction
error here rather than a reading mistake later.

Causal attribution is a decision rule, and it defaults to blocked
----------------------------------------------------------------
A quadratic population result supports the width-times-rate mechanism only if
the two negative controls do **not** also look quadratic.  :func:`causal_decision`
implements exactly the plan's rule and returns ``blocked`` with reasons whenever:

* the full model has no valid exponent at all;
* the central single-clock control's interval still contains two;
* the width-only control's interval still contains two;
* the **paired** full-minus-control exponent contrast fails to exceed the
  frozen minimum mechanistic difference; or
* the control matches or outperforms the full model on the frozen
  out-of-sample rising/falling survival and winner-mismatch comparison.

``not_blocked`` is not support.  It is the absence of one specific defeater,
and :data:`ATTRIBUTION_VERDICTS` is deliberately not spelled "supported".

Nothing here can produce a physical exponent
--------------------------------------------
Ticket 04 closed with the stationary killed-diffusion endpoint validation
*passed* at its reference tolerances and the moving-band audit returning a
``numerical_no_result`` pooled and in fourteen of its fifteen per-regime cells.
Every run this package writes therefore records ``numerical_gate =
"diagnostic_only"``.  Ticket 07 has since **frozen** a production numerical
budget and evaluated it: the disposition is ``numerical_no_result`` and no
matrix was feasible among the options searched, so the budget exists and is
**not met**, which is why no run is production-cleared.

:func:`scaling_verdict` is the gate that keeps that fact load-bearing: it
returns ``machinery_only`` unless *every* contributing run is
``production_cleared``, whatever the fit says.  The machinery in this layer may
be — and is — exercised on synthetic data with known answers.  A physical
exponent verdict cannot be produced by any argument to any function in this
file.  See :data:`COMPARISON_NON_CLAIMS`.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from . import analytic as _analytic
from . import model as _model
from . import observables as _observables
from .analysis import (CalibratedExponent, CouplingCell, MasterTrial,
                       MIN_RESAMPLE_REPLICATES, MIN_VALID_RESAMPLE_FRACTION,
                       RESAMPLE_REPLICATES, ScalingDesign, calibrate_exponent,
                       fit_exponent)
from .analysis import _authoritative_seed
from .analysis import _is_estimable as _estimable
from .analysis import _SEAL, _require_seal
from .raw_ledger import (ClosedLedger, MANIFEST_KEYS, MANIFEST_SCHEMA_VERSION,
                         RawManifest, TABLE_NAMES, digest_bytes,
                         require_unmutated_ledger)
from .raw_runner import open_raw_run
from .validation import (require_exact_type, require_finite,
                         require_finite_array, require_instance,
                         require_positive_int)

__all__ = [
    "ATTRIBUTION_VERDICTS",
    "ATTRIBUTION_BLOCKERS",
    "COMPARATOR_FINITE_POPULATION",
    "COMPARATOR_KINDS",
    "COMPARATOR_UNITS",
    "COMPARISON_NON_CLAIMS",
    "CONTROL_ROLES",
    "EXPERIMENT_FAMILY_EXCLUSIONS",
    "EXPERIMENT_FAMILY_FIELDS",
    "EXPERIMENT_VARYING_FIELDS",
    "PRODUCTION_STATUS_SCHEMA",
    "STATUS_FIELDS",
    "CRN_IDENTITY_FIELDS",
    "FROZEN_SEGMENTS",
    "PRIMARY_COMPARATOR",
    "PRODUCTION_GATE",
    "QUADRATIC_EXPONENT",
    "SCALING_VERDICTS",
    "SWEEP_VARYING_FIELDS",
    "CausalDecision",
    "ComparatorCurve",
    "ComparatorSet",
    "ComparisonScore",
    "ContrastReport",
    "ControlEvidence",
    "MachineryReport",
    "HoldoutEvidence",
    "LedgerFingerprint",
    "ReadOnlyRun",
    "RisingFallingSplit",
    "RunSweep",
    "ScalingVerdict",
    "causal_decision",
    "comparator_curves",
    "control_defeats",
    "control_evidence",
    "POPULATION_FIELDS",
    "crn_identity",
    "design_from_sweep",
    "experiment_family",
    "discrepancy_score",
    "frozen_edges",
    "holdout_evidence",
    "missing_holdout",
    "open_for_comparison",
    "paired_exponent_contrast",
    "rising_falling_split",
    "scaling_verdict",
    "sweep_master_trials",
    "synthetic_machinery_verdict",
]

#: The four comparison curves, separately named for their whole lives.
COMPARATOR_KINDS = ("eligible_count", "bare_rate_sum", "continuum_flux",
                    "constant_hazard")

#: The one the plan selects as the finite-population predictor.
PRIMARY_COMPARATOR = "bare_rate_sum"

#: What each curve is measured in.  Carried on the curve so that a plot axis
#: cannot be labelled from memory.
COMPARATOR_UNITS = {
    "eligible_count": "clocks (dimensionless count on the exact grid)",
    "bare_rate_sum": "inverse time (sum of local relaxation rates over the "
                     "exact finite production clocks)",
    "continuum_flux": "inverse time times detuning (spacing-weighted midpoint "
                      "quadrature; a CONTINUUM CONTROL, not a "
                      "finite-population predictor)",
    "constant_hazard": "inverse time (a declared constant, retained as a "
                       "falsifiable control)",
}

#: Whether each curve is a statement about the finite population that actually
#: ran.  ``continuum_flux`` is the one that is not, and this flag is what makes
#: that unmissable at every use site.
COMPARATOR_FINITE_POPULATION = {"eligible_count": True, "bare_rate_sum": True,
                                "continuum_flux": False,
                                "constant_hazard": False}

#: The frozen split of the pulse for the out-of-sample time-resolved
#: comparison.  The boundary is the manifest's own pulse centre, so there is no
#: argument anywhere in this module that could move it after seeing a result.
FROZEN_SEGMENTS = ("rising", "falling")

#: The exponent the mechanism story predicts, and the value the two negative
#: controls must *not* remain compatible with.
QUADRATIC_EXPONENT = 2.0

#: The two attribution outcomes.  ``not_blocked`` is the absence of a specific
#: defeater and is deliberately not spelled "supported".
ATTRIBUTION_VERDICTS = ("blocked", "not_blocked")

#: The closed vocabulary of reasons an attribution is blocked.
#:
#: The three ``no_valid_exponent`` blockers and ``holdout_evidence_missing``
#: were added after an independent review demonstrated ``not_blocked`` being
#: returned while *both* negative controls reported ``no_valid_exponent`` and
#: while the holdout was a bare boolean nobody had computed.  The premise that
#: an inadequate control has an unbounded interval and therefore blocks was
#: simply false for the ``lack_of_fit`` and ``curvature`` refusals, which retain
#: a finite and often narrow interval.
ATTRIBUTION_BLOCKERS = ("full_model_no_valid_exponent",
                        "central_control_no_valid_exponent",
                        "width_control_no_valid_exponent",
                        "central_control_compatible_with_quadratic",
                        "width_control_compatible_with_quadratic",
                        "contrast_inadequate",
                        "contrast_underpowered",
                        "holdout_evidence_missing",
                        "control_matches_out_of_sample",
                        "numerical_status_not_cleared")

#: The manifest value a run must carry before any physical scaling verdict is
#: permitted.  No run this package writes carries it.  Ticket 07 froze a
#: production numerical budget and evaluated it: the disposition is
#: ``numerical_no_result``, the budget is **not met**, and ticket 08 owns what
#: would follow a met one.
PRODUCTION_GATE = "production_cleared"

#: The manifest schema version that would carry a production-budget status.
#:
#: **It does not exist, and that is the point.**  An independent review took
#: this module's positive control — a real closed run whose manifest was
#: re-emitted with ``numerical_gate = "production_cleared"`` and re-hashed —
#: and observed that the *same hashed manifest* still said
#:
#:     endpoint_status     passed at the declared reference tolerances
#:     moving_band_status  numerical_no_result pooled and in fourteen of
#:                         fifteen per-regime cells
#:     validation_note     This run is therefore diagnostic ... no scaling
#:                         claim may be drawn from these tables.
#:
#: A gate field is free text.  Rewriting one word of it cannot create a
#: production budget, and a schema whose only coherent status tuple is
#: "diagnostic, moving band unbounded" must not be able to express a permission
#: at all.  So schema v3 **fails closed**: :class:`ScalingVerdict` refuses to
#: hold ``physical_verdict_permitted`` under it, and a future permission needs a
#: schema revision carrying a closed-vocabulary production-budget status rather
#: than a rewritten sentence.
PRODUCTION_STATUS_SCHEMA = 4

#: The manifest fields that decide whether a run's numbers may be read as
#: physics.  All four, jointly: an independent review found the gate checking
#: only the first while the other three contradicted it inside the same hash.
STATUS_FIELDS = ("numerical_gate", "endpoint_status", "moving_band_status",
                 "validation_note")

#: Substrings whose presence in a status field means "not production cleared".
#: Matched case-insensitively, so a reworded note that still says
#: ``numerical_no_result`` or ``diagnostic`` still blocks.
_BLOCKING_STATUS_MARKERS = ("numerical_no_result", "diagnostic",
                            "no scaling claim", "not been production",
                            "has not passed")

#: The two scaling outcomes this module can reach.
SCALING_VERDICTS = ("machinery_only", "physical_verdict_permitted")


#: Manifest fields the runs of one coupling sweep are allowed to differ in.
#:
#: A sweep varies the peak coupling, and the pulse digest is a function of it.
#: Nothing else -- and in particular not the model label or the fixed
#: contraction rate.  An independent review assembled a ``RunSweep`` out of
#: alternating full and width-only runs and fitted it; the ``model`` property
#: raised only if somebody happened to read it, and the same structure let the
#: one-rate width-only control carry a different rate in every cell.
#: :class:`RunSweep` now checks every field outside this pair, so one sweep is
#: one data-generating process at one contraction rate.
SWEEP_VARYING_FIELDS = ("peak_coupling", "pulse_digest")

#: Manifest fields two runs of one **experiment** may differ in.
#:
#: Wider than :data:`SWEEP_VARYING_FIELDS` by exactly the two a paired
#: full/width-only comparison is *about*.  Used where runs from the two sweeps
#: are compared with each other -- the holdout, and the cross-sweep
#: comparability checks -- and never to decide what one sweep is.
EXPERIMENT_VARYING_FIELDS = SWEEP_VARYING_FIELDS + ("model",
                                                    "fixed_contraction_rate")

#: Manifest fields that decide *which physical path* master trial ``i`` is.
#:
#: This is the provenance token behind ``analysis.MasterTrial``.  It contains
#: the counter-key hierarchy, the mesh, the population and the physical and
#: numerical parameters — everything that fixes the Brownian path and the
#: prepared initial phase of a given trial number.
#:
#: It deliberately excludes ``peak_coupling`` (a sweep varies it), ``model`` and
#: ``fixed_contraction_rate`` (full and width-only share one tree by
#: construction — the plan excludes the model label from the physical noise
#: key), and the status/environment fields, which describe the *run* rather than
#: the path.  Two runs agreeing here are two runs whose trial 7 is the same
#: path; two runs disagreeing here are not paired data however their rows happen
#: to be numbered.
CRN_IDENTITY_FIELDS = (
    "key_schema", "stream_namespace", "stream_prefix", "ensemble_prefix",
    "mesh_identity", "mesh_origin", "timestep", "start_index", "steps",
    "window_open", "window_close", "clocks", "grid_half_width",
    "grid_spacing", "grid_parity", "grid_origin", "grid_support", "detunings",
    "pulse_duration", "pulse_centre", "pulse_support", "phase_diffusion",
    "lock_tolerance", "dwell_time", "trials", "shadow_trials", "shadow_rule",
    "run_label",
)

#: Manifest fields that describe the clock *population* rather than the
#: experiment around it.  The central-clock control is the one experiment that
#: legitimately differs in all of them: it is one clock at detuning zero where
#: the full model has many.
POPULATION_FIELDS = ("clocks", "grid_half_width", "grid_spacing",
                     "grid_parity", "grid_origin", "grid_support", "detunings")

#: The manifest fields two roles of one causal decision may legitimately
#: differ in, and **the only ones**.
#:
#: Four sweep-level differences and the seven population fields:
#:
#: ``peak_coupling``, ``pulse_digest``
#:     a sweep varies the coupling and the digest is a function of it.
#: ``model``, ``fixed_contraction_rate``
#:     the width-only control *is* the changed model at the frozen rate.
#: :data:`POPULATION_FIELDS`
#:     the central control *is* the one zero-detuning clock.
#:
#: Everything else must match, and :data:`EXPERIMENT_FAMILY_FIELDS` is built by
#: subtraction from the whole manifest rather than by listing what to include.
#: That direction is the fix.  An earlier revision derived the family from
#: :data:`CRN_IDENTITY_FIELDS` — an allowlist assembled for a different purpose
#: — and an independent review found ``drive`` missing from it: a genuine
#: one-clock sweep whose manifest said ``stationary`` joined a ``pulsed``
#: causal experiment with no provenance refusal at all.  Subtracting a reviewed
#: exclusion list from ``MANIFEST_KEYS`` cannot omit a field, because a field
#: nobody thought about is *included* by default.
EXPERIMENT_FAMILY_EXCLUSIONS = (("peak_coupling", "pulse_digest", "model",
                                 "fixed_contraction_rate")
                                + POPULATION_FIELDS)

#: The fields that must agree across **every** role of one causal decision:
#: the whole frozen manifest minus :data:`EXPERIMENT_FAMILY_EXCLUSIONS`.
#:
#: It therefore includes the drive, the bath, the pulse geometry, the window
#: and mesh, the dwell, the lock tolerance, the timestep, the trial and shadow
#: counts, the stream hierarchy, the pilot/production label, the schema and
#: column layout, the validation statuses, the recorded environment and the
#: non-claims.  ``verify.py`` mutates **every** manifest field in turn and
#: requires each one either to change the family or to be on the exclusion list
#: above, so a future schema field is covered the day it is added.
EXPERIMENT_FAMILY_FIELDS = tuple(name for name in MANIFEST_KEYS
                                 if name not in EXPERIMENT_FAMILY_EXCLUSIONS)

#: The three roles a fitted sweep may play in the causal decision.
CONTROL_ROLES = ("full", "central", "width_only")

#: What a comparison produced here does not license anybody to say.
COMPARISON_NON_CLAIMS = (
    "nothing here demonstrates Born's rule, a detector click, an absorption, "
    "a measurement outcome, a unique actuality, a two-channel outcome, or a "
    "microscopic origin for the noise",
    "the bare relaxation-rate sum is the finite-population comparator and the "
    "spacing-weighted flux is a continuum control; they are different "
    "quantities with different units and neither predicts the other",
    "agreement between a hazard and a shifted instantaneous predictor is a "
    "heuristic: a commitment at t depends on the whole path over the "
    "preceding dwell, and no free time shift is fitted anywhere here",
    "'not_blocked' means one named defeater did not fire.  It is not causal "
    "support for the width-times-rate mechanism",
    "every run this package writes records numerical_gate='diagnostic_only', "
    "so scaling_verdict returns 'machinery_only' and no physical exponent "
    "verdict can be produced by any argument to any function in this file",
)


# ---------------------------------------------------------------------------
# Read-only access, with the hash carried
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class LedgerFingerprint:
    """The digests that identify one closed run's bytes.

    Taken from the close marker the gate has already validated, which means
    every one of these has been recomputed from the bytes actually read rather
    than believed from the file that declares them.

    Comparing two fingerprints is how this layer demonstrates that it did not
    modify what it read.  It is not authentication — anybody who can rewrite a
    table can recompute a digest — and it is not claimed to be; it is the
    check that catches truncation, an accidental edit, a half-written rerun and
    a comparison that opened a different run than it thought it had.
    """

    manifest_digest: str
    table_digests: tuple

    def __post_init__(self):
        digests = tuple(self.table_digests)
        if len(digests) != len(TABLE_NAMES):
            raise ValueError(
                f"a fingerprint carries one digest per table "
                f"{list(TABLE_NAMES)}, got {len(digests)}")
        for name, entry in zip(TABLE_NAMES, digests):
            pair = tuple(entry)
            if len(pair) != 2 or pair[0] != name:
                raise ValueError(f"table_digests must be "
                                 f"(name, digest) pairs in the order "
                                 f"{list(TABLE_NAMES)}, got {entry!r}")
            if not isinstance(pair[1], str) or not pair[1]:
                raise TypeError(f"the digest for {name!r} must be a non-empty "
                                f"string, got {pair[1]!r}")
        if not isinstance(self.manifest_digest, str) or not \
                self.manifest_digest:
            raise TypeError(f"manifest_digest must be a non-empty string, got "
                            f"{self.manifest_digest!r}")
        object.__setattr__(self, "table_digests",
                           tuple((name, str(digest))
                                 for name, digest in digests))


def _fingerprint_of(ledger) -> LedgerFingerprint:
    """The fingerprint a verified closed ledger's own marker declares."""
    marker = ledger.marker
    return LedgerFingerprint(marker.manifest_digest,
                             tuple((name, marker.digest(name))
                                   for name in TABLE_NAMES))


@dataclass(frozen=True)
class ReadOnlyRun:
    """One closed run, opened for comparison, with its fingerprint.

    There is no method here that writes anything, and no field that holds a
    file handle or a path that could be reopened for writing: the run is a
    tuple of parsed rows and a set of digests by the time it reaches this
    class.  :meth:`reverify` is the evidence for that claim rather than the
    mechanism enforcing it — it re-reads the run through the same gate and
    requires the same digests.
    """

    run_name: str
    ledger: ClosedLedger
    fingerprint: LedgerFingerprint
    seal: object = None

    def __post_init__(self):
        _require_seal(self.seal, "ReadOnlyRun")
        if isinstance(self.run_name, bool) or not isinstance(self.run_name,
                                                             str):
            raise TypeError(f"run_name must be a string, got "
                            f"{type(self.run_name).__name__} "
                            f"{self.run_name!r}")
        if not self.run_name:
            raise ValueError("run_name must not be empty")
        # Exact, both of them, for the reason `RunSweep` requires an exact
        # `ReadOnlyRun`.  A `ClosedLedger` subclass could override ``rows`` or
        # ``manifest`` and be read before `require_unmutated_ledger`'s own
        # exact-type gate ever sees it; a `LedgerFingerprint` subclass could
        # override ``__eq__`` to return `True`, which is every reopen
        # comparison in this module answered by the object being compared.
        require_exact_type(self.ledger, ClosedLedger, "ledger")
        require_exact_type(self.fingerprint, LedgerFingerprint, "fingerprint")
        declared = _fingerprint_of(self.ledger)
        if declared != self.fingerprint:
            raise ValueError(
                "the fingerprint supplied is not the one this ledger's own "
                "close marker declares; a run and a fingerprint from "
                "different reads prove nothing about either")
        object.__setattr__(self, "seal", None)

    @property
    def manifest(self) -> RawManifest:
        """The frozen manifest of this run."""
        return self.ledger.manifest

    def reopen(self):
        """Re-read this run from disk and return the **fresh** object.

        This is the operation every authoritative path uses, and the word
        *fresh* is the whole of it.  An earlier revision had ``reverify``
        re-open the run, compare fingerprints and then return only the
        fingerprint — after which the caller went on computing from the
        in-memory ledger it already had.  An independent review took a genuine
        closed sweep, built public ``ClosedLedger`` objects carrying the
        genuine close markers with every trial category rewritten to
        ``committed`` and all four status fields rewritten to
        ``production_cleared``, and observed:

            forged_sweep.reverify()   True
            genuine disk counts       (2/6, 4/6, 5/6, 6/6)
            fitted counts             (6/6, 6/6, 6/6, 6/6)
            verdict gates             production_cleared x 4

        The disk was never touched.  The fingerprint came from the marker,
        which was genuine; everything that was *used* came from the forged
        rows beside it.

        So the supplied handle is now a **locator**: its name says which run,
        its fingerprint is compared, and then its data are discarded.  What
        comes back is parsed from the bytes on disk by
        :func:`~.raw_runner.open_raw_run` — five files opened under pinned
        descriptors — and that is the only object the caller may compute from.
        """
        fresh = open_for_comparison(self.run_name)
        if fresh.fingerprint != self.fingerprint:
            raise ValueError(
                f"the run {self.run_name!r} does not hash to what it hashed "
                f"when it was opened.  Something between those two reads "
                f"modified the closed ledger, and nothing in this comparison "
                f"layer opens a file for writing")
        return fresh

    def reverify(self) -> LedgerFingerprint:
        """Re-read this run and return the fresh fingerprint.

        A thin reading of :meth:`reopen`, kept because "the bytes did not
        change" is a question worth asking on its own — the read-only proof
        brackets a whole comparison with it.  **It is not the authoritative
        read**: a path that will compute something calls :meth:`reopen` and
        uses what it returns.
        """
        return ReadOnlyRun.reopen(self).fingerprint


def _canonical_fields(manifest, names) -> str:
    """A canonical, injective spelling of some of a manifest's frozen fields.

    Length-prefixed on both the name and the value, so no combination of field
    values can be re-read as a different combination.
    """
    parts = []
    for name in names:
        value = manifest.value(name)
        if isinstance(value, tuple):
            text = "[" + ",".join(repr(item) for item in value) + "]"
        else:
            text = repr(value)
        parts.append(f"{len(name)}:{name}={len(text)}:{text}")
    return "|".join(parts)


def _status_blockers(manifest) -> tuple:
    """Every reason one manifest's own status forbids reading it as physics.

    Four fields, jointly, plus the schema itself.  An independent review found
    the earlier gate checking ``numerical_gate`` alone, and demonstrated a
    hash-valid manifest in which that one word said ``production_cleared``
    while the three beside it still said the run was diagnostic and that no
    scaling claim could be drawn from its tables.

    The schema check is the fail-closed one.  Version 3 has no
    closed-vocabulary production-budget status; its status fields are free
    text, and free text rewritten is not a budget passed.  So a permission
    needs :data:`PRODUCTION_STATUS_SCHEMA`, which does not exist, and this
    function therefore never returns empty for any run this package can read.
    """
    require_instance(manifest, RawManifest, "manifest")
    blockers = []
    version = manifest.value("schema_version")
    if version < PRODUCTION_STATUS_SCHEMA:
        blockers.append(
            f"manifest schema v{version} has no closed-vocabulary "
            f"production-budget status; a permission needs schema "
            f"v{PRODUCTION_STATUS_SCHEMA}, and rewriting a free-text status "
            f"field does not create one")
    gate = manifest.value("numerical_gate")
    if gate != PRODUCTION_GATE:
        blockers.append(f"numerical_gate is {gate!r}, not {PRODUCTION_GATE!r}")
    for field in STATUS_FIELDS:
        value = str(manifest.value(field)).lower()
        for marker in _BLOCKING_STATUS_MARKERS:
            if marker in value:
                blockers.append(f"{field} states {marker!r}")
    return tuple(sorted(set(blockers)))


def experiment_family(manifest) -> str:
    """Return the token every role of one causal decision must share.

    The digest of :data:`EXPERIMENT_FAMILY_FIELDS` — everything about the
    experiment except the clock population, which the central control
    legitimately differs in.  See :data:`EXPERIMENT_FAMILY_FIELDS`.
    """
    require_instance(manifest, RawManifest, "manifest")
    return digest_bytes(_canonical_fields(
        manifest, EXPERIMENT_FAMILY_FIELDS).encode("utf-8"))


def crn_identity(manifest) -> str:
    """Return the common-random-number provenance token of one run.

    The digest of :data:`CRN_IDENTITY_FIELDS` as this manifest froze them.  It
    is what ``analysis.MasterTrial`` carries as its namespace, so a trial number
    can no longer be paired across two datasets that merely happen to number
    their rows the same way.
    """
    require_instance(manifest, RawManifest, "manifest")
    return digest_bytes(_canonical_fields(manifest,
                                          CRN_IDENTITY_FIELDS).encode("utf-8"))


@dataclass(frozen=True)
class RunSweep:
    """The closed runs one coupling sweep consists of, verified as a set.

    This is the object that makes the numerical gate *manifest-bound*.  An
    earlier revision let :func:`scaling_verdict` take a tuple of strings, and an
    independent review pointed out the consequence: a good synthetic fit and the
    literal string ``"production_cleared"`` returned
    ``physical_verdict_permitted`` without a ledger being opened at all, which
    is exactly the permission every real ``diagnostic_only`` run must be denied.

    A sweep can only be built from :class:`ReadOnlyRun` objects, which can only
    be built by :func:`open_for_comparison`, which reads five files through the
    verified closed-ledger gate.  The construction then requires:

    * at least one run, each an exact ``ReadOnlyRun``;
    * strictly increasing, distinct peak couplings — one run per cell, in the
      order the fit will see them;
    * distinct fingerprints, so one run cannot be counted twice; and
    * **every** frozen manifest field outside :data:`SWEEP_VARYING_FIELDS`
      identical across the runs.  Not a chosen handful: the whole schema minus
      the four a sweep is allowed to vary.

    :meth:`reverify` re-opens all of them and requires unchanged bytes, and
    :func:`scaling_verdict` calls it before it will say anything.
    """

    runs: tuple

    def __post_init__(self):
        runs = tuple(self.runs)
        if not runs:
            raise ValueError(
                "a sweep holds at least one closed run; an empty sweep would "
                "let a fit with no run behind it be cleared, which is the "
                "defect this type exists to close")
        for position, run in enumerate(runs):
            if type(run) is not ReadOnlyRun:
                raise TypeError(
                    f"runs[{position}] must be exactly a ReadOnlyRun opened "
                    f"through open_for_comparison, got "
                    f"{type(run).__name__}; a gate that accepted a description "
                    f"of a run rather than a run is not a gate")
        couplings = [run.manifest.value("peak_coupling") for run in runs]
        if any(b <= a for a, b in zip(couplings, couplings[1:])):
            raise ValueError(
                f"the runs' peak couplings are {couplings}; a sweep lists its "
                f"runs once each in strictly increasing coupling order, so a "
                f"repeat or a shuffle is refused rather than sorted")
        prints = [run.fingerprint for run in runs]
        if len(set(prints)) != len(prints):
            raise ValueError("two runs in this sweep have the same "
                             "fingerprint, so one run is being counted twice")
        models = {run.manifest.value("model") for run in runs}
        if len(models) != 1:
            raise ValueError(
                f"this sweep mixes the models {sorted(models)}.  The full "
                f"model and its width-only control are two paired sweeps, not "
                f"one, and an alternating mixture fitted as a single coupling "
                f"sweep is four cells of two different experiments")
        rates = {run.manifest.value("fixed_contraction_rate") for run in runs}
        if len(rates) != 1:
            raise ValueError(
                f"this sweep carries {len(rates)} fixed contraction rates "
                f"{sorted(rates)}; the plan's width-only control uses ONE "
                f"reference rate common to every clock and every coupling, and "
                f"a rate that varied across the sweep would be a different "
                f"control at every cell")
        rate = rates.pop()
        if models == {"width_only"} and not rate > 0.0:
            raise ValueError(
                f"a width-only sweep records a fixed contraction rate of "
                f"{rate!r}; the control is defined by a single positive "
                f"manifest-frozen rate")
        if models == {"full"} and rate != 0.0:
            raise ValueError(
                f"a full-model sweep records a fixed contraction rate of "
                f"{rate!r}; the full Adler model has no substituted rate")
        shared = tuple(name for name in MANIFEST_KEYS
                       if name not in SWEEP_VARYING_FIELDS)
        reference = _canonical_fields(runs[0].manifest, shared)
        for position, run in enumerate(runs[1:], start=1):
            if _canonical_fields(run.manifest, shared) != reference:
                differing = sorted(
                    name for name in shared
                    if run.manifest.value(name) != runs[0].manifest.value(name))
                raise ValueError(
                    f"runs[0] and runs[{position}] disagree about the frozen "
                    f"fields {differing}; a sweep varies "
                    f"{list(SWEEP_VARYING_FIELDS)} and nothing else, so these "
                    f"are runs of two different experiments")
        object.__setattr__(self, "runs", runs)

    @property
    def couplings(self) -> tuple:
        """The peak couplings, low to high."""
        return tuple(run.manifest.value("peak_coupling") for run in self.runs)

    @property
    def numerical_gates(self) -> tuple:
        """The gate each run's own hashed manifest records."""
        return tuple(run.manifest.value("numerical_gate") for run in self.runs)

    @property
    def crn_identity(self) -> str:
        """The provenance token shared by every run in this sweep."""
        return crn_identity(self.runs[0].manifest)

    @property
    def experiment_family(self) -> str:
        """The token this sweep shares with the other roles of its decision."""
        return experiment_family(self.runs[0].manifest)

    @property
    def model(self) -> str:
        """The one data-generating process these runs share.

        Settled at construction rather than raised here: an independent review
        found that a mixed sweep was accepted and fitted, and that this
        property raised only if somebody happened to read it.
        """
        return self.runs[0].manifest.value("model")

    @property
    def fixed_contraction_rate(self) -> float:
        """The one contraction rate these runs share."""
        return self.runs[0].manifest.value("fixed_contraction_rate")

    @property
    def status_blockers(self) -> tuple:
        """Why this sweep's own manifests forbid a physical reading.

        Empty only for a sweep every one of whose runs is production cleared
        under a schema that can express one.  Currently never empty; see
        :data:`PRODUCTION_STATUS_SCHEMA`.
        """
        blockers = []
        for run in self.runs:
            blockers.extend(_status_blockers(run.manifest))
        return tuple(sorted(set(blockers)))

    def reopen(self):
        """Re-read every run from disk and return the **fresh** sweep.

        The authoritative read.  Each run is re-opened by name through
        :func:`~.raw_runner.open_raw_run`, its fingerprint is compared against
        the one this handle carried, and the handle's own parsed rows and
        manifest are then discarded.  The returned :class:`RunSweep` is rebuilt
        from the fresh runs, so every field this class checks — one model, one
        contraction rate, one experiment — is re-established on the bytes that
        are actually on disk rather than on whatever was in memory.

        Every factory that computes something calls this immediately before
        building a design.  See :meth:`ReadOnlyRun.reopen` for the replay that
        made it necessary.

        The per-run call is written ``ReadOnlyRun.reopen(run)`` rather than
        ``run.reopen()``.  ``__post_init__`` already requires each run to be an
        exact :class:`ReadOnlyRun`, so the two are the same call today; writing
        the unbound form means the authoritative read stays this class's read
        even if the exactness above is ever relaxed, instead of becoming
        whatever the object in the tuple decided ``reopen`` should mean.
        """
        return RunSweep(tuple(ReadOnlyRun.reopen(run) for run in self.runs))

    def reverify(self) -> tuple:
        """Re-open every run and return the fresh fingerprints.

        The "did the bytes change" question, asked on its own.  A path that
        computes something calls :meth:`reopen` and uses what it returns.
        """
        return tuple(run.fingerprint for run in RunSweep.reopen(self).runs)


def _fresh_sweep(sweep, name="sweep"):
    """Re-read ``sweep`` from disk, refusing anything but an exact ``RunSweep``.

    The single door every sweep-taking factory goes through, and it closes two
    holes an independent review opened with seven public lines::

        class StaleSweep(compare.RunSweep):
            def reopen(self):
                return self

    ``RunSweep`` is public and ordinary, so that subclass is public and
    ordinary too.  Against an ``isinstance`` boundary followed by
    ``sweep.reopen()`` it turned the authoritative disk read into an assertion:
    with one byte appended to a contributing ``manifest.json``, the genuine
    sweep refused with ``manifest digest mismatch`` while ``StaleSweep`` was
    accepted by ``design_from_sweep``, ``scaling_verdict``,
    ``control_evidence`` and ``paired_exponent_contrast`` alike.  The rows it
    carried were internally valid, so no re-encode guard could see that the
    named files had changed; only the skipped reopen could.

    Both halves matter.  The exact type is what refuses the subclass, and
    ``RunSweep.reopen(sweep)`` — the unbound, trusted implementation — is what
    makes the refusal the *only* thing standing between a caller and this
    package's own read, rather than one of two.  A boundary whose correctness
    depends on a method the caller supplied is not a boundary.
    """
    require_exact_type(sweep, RunSweep, name)
    return RunSweep.reopen(sweep)


def _reconciled_ledgers(sweep, name="sweep"):
    """Reconcile every run of a sweep and re-establish the sweep's invariants.

    Returns the tuple of freshly parsed :class:`~.raw_ledger.ClosedLedger`
    objects :func:`~.raw_ledger.require_unmutated_ledger` hands back — one per
    run, in coupling order — and **nothing** of the caller's.  Everything a
    design, a count or an identity is built from is then read out of these.

    The round-six contract was that every value come from the returned object,
    and an independent review found two places where it did not.  ``_sweep_counts``
    took its numerator from the reconciled ledger and its denominator from the
    caller's ``run.manifest.trials``: an interception that returned the genuine
    fresh ledger while changing the caller's alias from six trials to seven
    moved the helper from ``((2, 6), (4, 6), (5, 6), (6, 6))`` to
    ``((2, 7), (4, 7), (5, 7), (6, 7))``.  ``_design_from_fresh_sweep`` built
    its master-trial identities from the caller's first manifest *before*
    reconciling anything, so a design could carry one namespace's identities
    beside another snapshot's cells.

    The invariants ``RunSweep.__post_init__`` established over the caller's
    manifests are therefore re-established here over the reconciled ones.  They
    are the same rules, deliberately: a rule checked on one set of objects and
    relied on for another is not checked at all.
    """
    require_exact_type(sweep, RunSweep, name)
    ledgers = tuple(require_unmutated_ledger(run.ledger) for run in sweep.runs)
    manifests = tuple(ledger.manifest for ledger in ledgers)

    couplings = [manifest.value("peak_coupling") for manifest in manifests]
    if any(b <= a for a, b in zip(couplings, couplings[1:])):
        raise ValueError(
            f"the reconciled ledgers' peak couplings are {couplings}; a sweep "
            f"lists its runs once each in strictly increasing coupling order, "
            f"and the order the fit will see is the order on disk")
    models = {manifest.value("model") for manifest in manifests}
    if len(models) != 1:
        raise ValueError(
            f"the reconciled ledgers mix the models {sorted(models)}; the full "
            f"model and its width-only control are two paired sweeps, not one")
    rates = {manifest.value("fixed_contraction_rate") for manifest in manifests}
    if len(rates) != 1:
        raise ValueError(
            f"the reconciled ledgers carry {len(rates)} fixed contraction "
            f"rates {sorted(rates)}; the width-only control uses ONE reference "
            f"rate common to every clock and every coupling")
    shared = tuple(field for field in MANIFEST_KEYS
                   if field not in SWEEP_VARYING_FIELDS)
    reference = _canonical_fields(manifests[0], shared)
    for position, manifest in enumerate(manifests[1:], start=1):
        if _canonical_fields(manifest, shared) != reference:
            differing = sorted(field for field in shared
                               if manifest.value(field)
                               != manifests[0].value(field))
            raise ValueError(
                f"the reconciled ledgers of runs[0] and runs[{position}] "
                f"disagree about the frozen fields {differing}; a sweep varies "
                f"{list(SWEEP_VARYING_FIELDS)} and nothing else, so these are "
                f"runs of two different experiments")
    return ledgers


def _fresh_master_trials(ledgers):
    """The master-trial identities the *reconciled* ledgers declare.

    Both the count and the namespace come from ``ledgers[0].manifest``, which
    :func:`_reconciled_ledgers` has just re-derived from the bytes and required
    every other run to agree with outside :data:`SWEEP_VARYING_FIELDS`.
    """
    manifest = ledgers[0].manifest
    namespace = crn_identity(manifest)
    return tuple(MasterTrial(namespace, index)
                 for index in range(manifest.trials))


def sweep_master_trials(sweep):
    """Return the master-trial identities of one sweep's shared trial count.

    The bridge between a verified sweep and a :class:`~.analysis.ScalingDesign`:
    the namespace is the sweep's own manifest-derived
    :func:`crn_identity`, so a design built from these identities cannot be
    paired against a design built from another dataset's.

    Public, and therefore reached the same way every other public sweep door is
    reached: the argument must be an exact :class:`RunSweep`, it is re-read from
    disk, and the count and the namespace are taken from the reconciled
    manifests rather than from the handle that was passed in.
    """
    return _fresh_master_trials(
        _reconciled_ledgers(_fresh_sweep(sweep), "sweep"))


def design_from_sweep(sweep) -> ScalingDesign:
    """Build the scaling design of a sweep, **re-read from disk**.

    The public door, and it reopens.  An independent review mutated one parsed
    row of a verified ledger — a plain dictionary, in a frozen dataclass that
    freezes only the reference — and watched this function produce a design
    with a different digest, which then calibrated::

        row = sweep.runs[0].ledger.rows("ledger")[0]
        row["category"] = "committed"
        design_from_sweep(sweep).digest != disk_digest       # True

    No file changed and no private name was used.  So the sweep handed in is a
    locator here as it is everywhere else: it is reopened, and the design is
    built from the bytes on disk.

    One cell per run, in coupling order; one outcome per trial, read from the
    trial ledger's own ``category`` column, which
    :func:`~.raw_ledger.require_closed_ledger` has already re-derived from the
    authoritative per-clock rows.  The identities are
    :func:`sweep_master_trials`, so the design's namespace *is* the sweep's
    manifest-derived provenance.

    The four decision factories reopen once themselves and then call
    :func:`_design_from_fresh_sweep`, so a verdict costs one read rather than
    two.
    """
    return _design_from_fresh_sweep(_fresh_sweep(sweep))


def _design_from_fresh_sweep(sweep) -> ScalingDesign:
    """Build the design of a sweep that has *just* been reopened.

    The private half of :func:`design_from_sweep`.  Every value below — the
    trial count, the common-random-number namespace, each cell's coupling and
    each cell's outcomes — is read out of the tuple
    :func:`_reconciled_ledgers` returns, and none out of the sweep handed in.

    An earlier revision took the identities from ``sweep_master_trials(sweep)``
    *before* reconciling any run, which left a check/use interval an
    interception could open: the design came back carrying one manifest's
    master namespace beside another snapshot's cells.  Reconciling first and
    deriving afterwards is what closes it, and it is why this helper's own
    guarantee no longer depends on its caller having reopened.
    """
    ledgers = _reconciled_ledgers(sweep)
    identities = _fresh_master_trials(ledgers)
    cells = []
    for run, ledger in zip(sweep.runs, ledgers):
        outcomes = tuple(row["category"] == "committed"
                         for row in ledger.rows("ledger"))
        if len(outcomes) != len(identities):
            raise ValueError(
                f"run {run.run_name!r} holds {len(outcomes)} trials and the "
                f"sweep's shared identities cover {len(identities)}")
        cells.append(CouplingCell(ledger.manifest.value("peak_coupling"),
                                  identities, outcomes))
    return ScalingDesign(tuple(cells), "master_trial")


def _sweep_counts(sweep):
    """``(committed, trials)`` per run, derived from the closed ledgers.

    Both numbers, from the same reconciled object.  The denominator used to be
    ``run.manifest.trials`` — the caller's alias — which is how an independent
    review moved this helper's answer from ``(2, 6), (4, 6), (5, 6), (6, 6)``
    to ``(2, 7), (4, 7), (5, 7), (6, 7)`` while the numerators stayed fresh.
    """
    counts = []
    for ledger in _reconciled_ledgers(sweep):
        committed = sum(1 for row in ledger.rows("ledger")
                        if row["category"] == "committed")
        counts.append((committed, ledger.manifest.trials))
    return tuple(counts)


def open_for_comparison(run_name) -> ReadOnlyRun:
    """Open a closed run for comparison, recording its fingerprint.

    The only door into this module.  It refuses everything
    :func:`~.raw_runner.open_raw_run` refuses — an open run with no close
    marker, a manifest that does not hash to what the marker declares, a table
    that has been truncated or edited, a row count that disagrees three ways, a
    trial recorded twice, a semantically impossible row — before any comparison
    exists to be performed on it.
    """
    if isinstance(run_name, bool) or not isinstance(run_name, str):
        raise TypeError(f"run_name must be a string, got "
                        f"{type(run_name).__name__} {run_name!r}")
    ledger = open_raw_run(run_name)
    return ReadOnlyRun(run_name, ledger, _fingerprint_of(ledger), seal=_SEAL)


# ---------------------------------------------------------------------------
# The four comparators
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class ComparatorCurve:
    """One named comparison curve, with its units and its population status.

    The three fields after ``values`` are not decoration.  A curve that could
    be built with any unit and any finite-population flag would let the
    spacing-weighted continuum quadrature be handed to a caller expecting the
    finite-population predictor, which is precisely the substitution the plan
    forbids — and it would look identical at the call site.  So ``kind``
    determines both, they are checked against :data:`COMPARATOR_UNITS` and
    :data:`COMPARATOR_FINITE_POPULATION`, and a mismatch is a construction
    error.
    """

    kind: str
    times: tuple
    values: tuple
    unit: str
    finite_population: bool

    def __post_init__(self):
        if self.kind not in COMPARATOR_KINDS:
            raise ValueError(f"kind must be one of {list(COMPARATOR_KINDS)}, "
                             f"got {self.kind!r}")
        times = tuple(require_finite(t, f"times[{i}]")
                      for i, t in enumerate(self.times))
        values = tuple(require_finite(v, f"values[{i}]")
                       for i, v in enumerate(self.values))
        if len(times) != len(values):
            raise ValueError(f"this curve has {len(times)} times and "
                             f"{len(values)} values")
        if not times:
            raise ValueError("a comparator curve has at least one sample")
        if any(b <= a for a, b in zip(times, times[1:])):
            raise ValueError("times must be strictly increasing")
        if self.unit != COMPARATOR_UNITS[self.kind]:
            raise ValueError(
                f"a {self.kind!r} curve is measured in "
                f"{COMPARATOR_UNITS[self.kind]!r}, and this one claims "
                f"{self.unit!r}.  The unit is what the kind means, not a "
                f"label a caller chooses")
        if not isinstance(self.finite_population, bool):
            raise TypeError(f"finite_population must be a bool, got "
                            f"{type(self.finite_population).__name__}")
        declared = COMPARATOR_FINITE_POPULATION[self.kind]
        if self.finite_population != declared:
            is_or_is_not = "is" if declared else "is not"
            raise ValueError(
                f"a {self.kind!r} curve {is_or_is_not} a statement about the "
                f"finite population that ran, and this one claims otherwise.  "
                f"Continuum flux presented as a finite-population prediction "
                f"is the exact conflation this type exists to refuse")
        object.__setattr__(self, "times", times)
        object.__setattr__(self, "values", values)

    @property
    def is_primary(self) -> bool:
        """``True`` only for :data:`PRIMARY_COMPARATOR`."""
        return self.kind == PRIMARY_COMPARATOR


@dataclass(frozen=True)
class ComparatorSet:
    """All four curves of one run, on one time axis, one of each kind.

    Holding them together rather than returning a bare list is what makes the
    completeness checkable: the constructor requires exactly the four declared
    kinds, exactly once each, on identical times.  A set missing the continuum
    control would let the bare sum be reported as though nothing else had been
    computed; a set holding two curves of one kind would let a plot legend say
    one thing and a table another.
    """

    curves: tuple

    def __post_init__(self):
        curves = tuple(self.curves)
        for position, curve in enumerate(curves):
            require_instance(curve, ComparatorCurve, f"curves[{position}]")
        kinds = tuple(curve.kind for curve in curves)
        if sorted(kinds) != sorted(COMPARATOR_KINDS):
            raise ValueError(
                f"a comparator set holds exactly one of each of "
                f"{list(COMPARATOR_KINDS)}, got {list(kinds)}")
        times = curves[0].times
        for curve in curves[1:]:
            if curve.times != times:
                raise ValueError(
                    f"the {curve.kind!r} curve is sampled at different times "
                    f"than the {curves[0].kind!r} curve; curves compared "
                    f"against one hazard share one time axis")
        object.__setattr__(self, "curves", curves)

    def curve(self, kind) -> ComparatorCurve:
        """Return one curve by kind, or raise naming the four that exist."""
        if isinstance(kind, bool) or not isinstance(kind, str):
            raise TypeError(f"kind must be a string, got "
                            f"{type(kind).__name__} {kind!r}")
        for curve in self.curves:
            if curve.kind == kind:
                return curve
        raise ValueError(f"no comparator {kind!r}; this set holds "
                         f"{list(COMPARATOR_KINDS)}")


def comparator_curves(manifest, times, constant_hazard,
                      dwell_shifted=True) -> ComparatorSet:
    """Build all four comparison curves for one run's frozen parameters.

    Every quantity is computed on the run's **exact production detunings**, read
    from the manifest — not on a refined grid, not on a continuum
    approximation, and not on a grid rebuilt from a half-width and a count that
    might round differently.  The staircase the finite grid produces at weak
    coupling is a physical feature of a finite population of competitors and is
    reported rather than smoothed away.

    ``dwell_shifted`` selects between the two declared time conventions and
    nothing else.  The fixed dwell introduces a real delay, so a predictor may
    be evaluated at ``t - dwell_time``; that shift is the manifest's own dwell
    and the only alternative is no shift at all.  **No free time offset is
    fitted anywhere in this module**, because a fitted shift would convert a
    heuristic into a two-parameter family that can follow almost any hazard.
    """
    require_instance(manifest, RawManifest, "manifest")
    sample_times = np.atleast_1d(require_finite_array(times, "times"))
    if sample_times.ndim != 1 or sample_times.size < 1:
        raise ValueError(f"times must be a non-empty one-dimensional array, "
                         f"got shape {sample_times.shape}")
    if sample_times.size > 1 and not bool(np.all(np.diff(sample_times) > 0.0)):
        raise ValueError("times must be strictly increasing")
    rate = require_finite(constant_hazard, "constant_hazard")
    if rate < 0.0:
        raise ValueError(f"constant_hazard must be non-negative, got {rate!r}; "
                         f"a negative rate is not a slow one")
    if not isinstance(dwell_shifted, bool):
        raise TypeError(f"dwell_shifted must be a bool, got "
                        f"{type(dwell_shifted).__name__}; it selects between "
                        f"the two declared conventions and is not a free "
                        f"offset")

    pulse = _model.RaisedCosinePulse(manifest.value("peak_coupling"),
                                     manifest.value("pulse_duration"),
                                     manifest.value("pulse_centre"))
    detunings = np.array(manifest.value("detunings"), dtype=float)
    spacing = manifest.value("grid_spacing")
    shift = manifest.value("dwell_time") if dwell_shifted else 0.0

    evaluated = sample_times - shift
    couplings = np.asarray(pulse.coupling(evaluated), dtype=float)

    counts, sums = [], []
    for coupling in np.atleast_1d(couplings):
        eligible = _analytic.eligible_mask(float(coupling), detunings)
        mask = np.atleast_1d(np.asarray(eligible, dtype=bool))
        counts.append(float(np.count_nonzero(mask)))
        if bool(np.any(mask)):
            rates = _analytic.local_relaxation_rate(float(coupling),
                                                    detunings[mask])
            sums.append(float(np.sum(np.atleast_1d(rates))))
        else:
            sums.append(0.0)

    axis = tuple(float(t) for t in sample_times)
    bare = tuple(sums)
    return ComparatorSet((
        ComparatorCurve("eligible_count", axis, tuple(counts),
                        COMPARATOR_UNITS["eligible_count"], True),
        ComparatorCurve("bare_rate_sum", axis, bare,
                        COMPARATOR_UNITS["bare_rate_sum"], True),
        ComparatorCurve("continuum_flux", axis,
                        tuple(spacing * value for value in bare),
                        COMPARATOR_UNITS["continuum_flux"], False),
        ComparatorCurve("constant_hazard", axis,
                        tuple(rate for _ in axis),
                        COMPARATOR_UNITS["constant_hazard"], False),
    ))


# ---------------------------------------------------------------------------
# The frozen out-of-sample comparison
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class RisingFallingSplit:
    """One run's survival, partitioned at the pulse centre.

    The partition is applied to a risk table that already covers the *whole*
    observation window, so this is a labelling of rows rather than a choice of
    window.  There is no argument here that could move the boundary: it is the
    manifest's ``pulse_centre``, frozen before the run was written.

    ``rising`` and ``falling`` hold the display rows of
    :meth:`~.observables.RiskTable.rows`, so each one still carries its risk
    count, event count, censoring count and width.
    """

    centre: float
    rising: tuple
    falling: tuple

    def __post_init__(self):
        centre = require_finite(self.centre, "centre")
        rising, falling = tuple(self.rising), tuple(self.falling)
        if not rising or not falling:
            raise ValueError(
                f"the pulse centre {centre!r} leaves "
                f"{len(rising)} rising and {len(falling)} falling rows; a "
                f"comparison of rising against falling behaviour needs both, "
                f"and an empty side means the table does not straddle the "
                f"centre")
        for name, rows in (("rising", rising), ("falling", falling)):
            for position, row in enumerate(rows):
                missing = [c for c in _observables.DISPLAY_COLUMNS
                           if c not in row]
                if missing:
                    raise ValueError(
                        f"{name}[{position}] is missing {missing}; every row "
                        f"in a split carries the full display record, so a "
                        f"conditional failure probability can never be read "
                        f"without its risk count and its width")
        object.__setattr__(self, "centre", centre)
        object.__setattr__(self, "rising", rising)
        object.__setattr__(self, "falling", falling)

    def survival(self, segment) -> tuple:
        """Return the survival values of one frozen segment."""
        return tuple(row["survival"] for row in self._rows(segment))

    def bins(self, segment) -> tuple:
        """Return ``(left, right, width)`` for each row of one segment.

        Carried so that :func:`discrepancy_score` can require two splits to be
        binned *identically* rather than merely to the same length.  An
        independent review found the earlier comparison checking only that the
        two survival sequences had the same number of entries, which let two
        differently binned tables be differenced as though they described the
        same instants.
        """
        return tuple((row["left"], row["right"], row["width"])
                     for row in self._rows(segment))

    def _rows(self, segment):
        if segment not in FROZEN_SEGMENTS:
            raise ValueError(f"segment must be one of "
                             f"{list(FROZEN_SEGMENTS)}, got {segment!r}")
        return self.rising if segment == "rising" else self.falling


def frozen_edges(manifest, bins_per_side):
    """Return whole-window edges with the pulse centre exactly on one of them.

    The frozen binning the out-of-sample comparison uses.  Uniform edges over
    the whole window need not put the pulse centre on a boundary, and an
    interval that straddles the centre has to be assigned wholly to one side —
    which makes the rising/falling score depend on the bin count after the
    fact.  So the two halves are binned separately, each uniformly, and the
    centre is an edge by construction.

    ``bins_per_side`` is the only freedom, it applies to both sides equally,
    and the span is always the run's whole observation window.
    """
    require_instance(manifest, RawManifest, "manifest")
    count = require_positive_int(bins_per_side, "bins_per_side", minimum=1)
    opened = manifest.value("window_open")
    closed = manifest.value("window_close")
    centre = manifest.value("pulse_centre")
    if not opened < centre < closed:
        raise ValueError(
            f"the pulse centre {centre!r} is not strictly inside the "
            f"observation window ({opened!r}, {closed!r}), so this run has no "
            f"rising and falling halves to compare")
    left = _obs_uniform(opened, centre, count)
    right = _obs_uniform(centre, closed, count)
    return left + right[1:]


def _obs_uniform(start, stop, count):
    """``observables.uniform_edges`` under a name this module owns."""
    return _observables.uniform_edges(start, stop, count)


def rising_falling_split(manifest, table) -> RisingFallingSplit:
    """Partition a whole-window risk table at the manifest's pulse centre.

    Two requirements, both refusals:

    * the table must be the run's own — opening at ``window_open`` and closing
      at ``window_close``.  A table over a sub-window would make the split a
      second, hidden selection on top of the first; and
    * **the pulse centre must be one of the table's edges.**  An interval that
      straddles the centre belongs wholly to neither side, and the earlier
      implementation assigned it to the falling half, so an arbitrary
      whole-window bin count could move rows across the boundary and change the
      frozen score after the fact.  :func:`frozen_edges` produces edges that
      satisfy this by construction.
    """
    require_instance(manifest, RawManifest, "manifest")
    require_instance(table, _observables.RiskTable, "table")
    opened = manifest.value("window_open")
    closed = manifest.value("window_close")
    if table.entry_time != opened or table.intervals[-1].right != closed:
        raise ValueError(
            f"this table spans [{table.entry_time!r}, "
            f"{table.intervals[-1].right!r}] and the run's window is "
            f"[{opened!r}, {closed!r}]; the rising/falling split is a "
            f"labelling of a whole-window table, never a second window choice")
    centre = manifest.value("pulse_centre")
    edges = ((table.entry_time,)
             + tuple(interval.right for interval in table.intervals))
    if centre not in edges:
        raise ValueError(
            f"the pulse centre {centre!r} is not one of this table's edges "
            f"{list(edges)}.  An interval straddling the centre belongs wholly "
            f"to neither half, and assigning it to one of them would let the "
            f"bin count move rows across the frozen boundary; use "
            f"frozen_edges(manifest, bins_per_side)")
    rows = table.rows()
    rising = tuple(row for row in rows if row["right"] <= centre)
    falling = tuple(row for row in rows if row["right"] > centre)
    return RisingFallingSplit(centre, rising, falling)


@dataclass(frozen=True)
class ComparisonScore:
    """How far one run's time-resolved behaviour sits from a reference.

    Three non-negative discrepancies with no free weights between them:
    mean absolute survival difference on the rising side, the same on the
    falling side, and the total-variation distance between winner-mismatch
    distributions.  :attr:`total` is their plain sum, which is the frozen
    scoring rule — a weighted sum would be a knob, and a knob is what turns a
    losing control into a winning one.
    """

    rising: float
    falling: float
    winner: float

    def __post_init__(self):
        for name in ("rising", "falling", "winner"):
            value = require_finite(getattr(self, name), name)
            if value < 0.0:
                raise ValueError(f"{name} must be non-negative, got {value!r}")
            object.__setattr__(self, name, value)

    @property
    def total(self) -> float:
        """The unweighted sum of the three discrepancies."""
        return self.rising + self.falling + self.winner


def _mean_absolute(a, b) -> float:
    if len(a) != len(b):
        raise ValueError(f"survival sequences have {len(a)} and {len(b)} "
                         f"entries; a discrepancy is only defined between "
                         f"tables binned identically")
    if not a:
        return 0.0
    return float(np.mean(np.abs(np.array(a, dtype=float)
                                - np.array(b, dtype=float))))


def _require_same_bins(reference, candidate, segment):
    """Require two splits to be binned identically on one segment."""
    left, right = reference.bins(segment), candidate.bins(segment)
    if left != right:
        raise ValueError(
            f"the {segment} halves are binned differently — {len(left)} and "
            f"{len(right)} intervals, first difference at "
            f"{next((i for i, (a, b) in enumerate(zip(left, right)) if a != b), len(left))}"
            f".  A survival discrepancy is only defined between tables whose "
            f"intervals have the same left edge, right edge and width; equal "
            f"lengths are not equal binnings")


def discrepancy_score(reference_split, candidate_split, reference_winners,
                      candidate_winners) -> ComparisonScore:
    """Score a candidate run against a reference, on the frozen quantities.

    The comparison is *out of sample* in the sense the plan means: the
    quantities scored here — rising and falling survival, and the distribution
    of which part of the tongue supplied the first commitment — are not the
    endpoint counts the exponent was fitted to.  An endpoint exponent that
    matches while these do not is coincidence or incompleteness, not
    mechanistic support.

    Winner distributions are compared as total-variation distance between the
    normalized win counts, which is in ``[0, 1]`` and is defined even when one
    of the runs never committed at all — that case scores zero shared mass and
    a distance of one, rather than raising.
    """
    require_instance(reference_split, RisingFallingSplit, "reference_split")
    require_instance(candidate_split, RisingFallingSplit, "candidate_split")
    require_instance(reference_winners, _observables.WinnerMismatch,
                     "reference_winners")
    require_instance(candidate_winners, _observables.WinnerMismatch,
                     "candidate_winners")
    if reference_split.centre != candidate_split.centre:
        raise ValueError(
            f"the two splits are cut at {reference_split.centre!r} and "
            f"{candidate_split.centre!r}; a paired comparison uses one frozen "
            f"boundary")
    if reference_winners.detunings != candidate_winners.detunings:
        raise ValueError(
            "the two runs have different clock populations, so their "
            "winner-mismatch distributions are over different supports and "
            "the distance between them is not defined")

    for segment in FROZEN_SEGMENTS:
        _require_same_bins(reference_split, candidate_split, segment)
    rising = _mean_absolute(reference_split.survival("rising"),
                            candidate_split.survival("rising"))
    falling = _mean_absolute(reference_split.survival("falling"),
                             candidate_split.survival("falling"))

    a_total, b_total = (reference_winners.co_completers,
                        candidate_winners.co_completers)
    if a_total == 0 and b_total == 0:
        winner = 0.0
    elif a_total == 0 or b_total == 0:
        winner = 1.0
    else:
        a = np.array(reference_winners.wins, dtype=float) / a_total
        b = np.array(candidate_winners.wins, dtype=float) / b_total
        winner = 0.5 * float(np.sum(np.abs(a - b)))
    return ComparisonScore(rising, falling, winner)


def control_defeats(full_score, control_score) -> bool:
    """Does the control match or outperform the full model out of sample?

    ``True`` when the control's total discrepancy is no larger than the full
    model's.  The comparison is ``<=`` rather than ``<`` on purpose: a control
    that does exactly as well as the model it is supposed to lack the mechanism
    of has not been distinguished from it, and calling that a pass would let a
    tie certify the mechanism.
    """
    require_instance(full_score, ComparisonScore, "full_score")
    require_instance(control_score, ComparisonScore, "control_score")
    return control_score.total <= full_score.total


# ---------------------------------------------------------------------------
# The paired contrast and the causal decision
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class ContrastReport:
    """The paired full-minus-control exponent difference, bound to its designs.

    Paired, because the two designs share master trial identities by
    construction: the width-only control is run on the same initial phases and
    the same Brownian tree as the full model, so the two exponents are
    correlated and their difference is estimated far more precisely than
    ``sqrt(se_full^2 + se_control^2)`` would suggest.

    ``estimable_replicates`` counts a replicate only when **both** rebuilt
    exponents are estimable.

    ``full_digest``, ``control_digest``, ``couplings``, ``crn_identity`` and
    the two fingerprint tuples are why this record is trustworthy at all.  An
    independent review observed that the earlier one "carries no coupling grid,
    CRN identity, design digest, or fingerprints", so "a caller-authored
    adequate interval is indistinguishable from the output of
    ``paired_exponent_contrast``" — and used one to manufacture a
    ``not_blocked``.  :func:`causal_decision` now requires every one of these
    to match the evidence it is deciding about.
    """

    replicates: int
    estimable_replicates: int
    mean: float
    stderr: float
    interval: tuple
    full_digest: str
    control_digest: str
    couplings: tuple
    crn_identity: str
    full_fingerprints: tuple
    control_fingerprints: tuple
    seal: object = None

    def __post_init__(self):
        _require_seal(self.seal, "ContrastReport")
        require_positive_int(self.replicates, "replicates", minimum=1)
        if not 0 <= self.estimable_replicates <= self.replicates:
            raise ValueError(f"{self.estimable_replicates} estimable of "
                             f"{self.replicates} replicates")
        for name in ("full_digest", "control_digest", "crn_identity"):
            value = getattr(self, name)
            if not isinstance(value, str) or not value:
                raise TypeError(f"{name} must be a non-empty string")
        if self.full_digest == self.control_digest:
            raise ValueError(
                "the full model and the control were fitted to the same "
                "observations, so their difference is identically zero and "
                "there is no contrast here")
        for name in ("full_fingerprints", "control_fingerprints"):
            for position, entry in enumerate(tuple(getattr(self, name))):
                require_instance(entry, LedgerFingerprint,
                                 f"{name}[{position}]")
        object.__setattr__(self, "interval", tuple(self.interval))
        object.__setattr__(self, "couplings", tuple(self.couplings))
        object.__setattr__(self, "full_fingerprints",
                           tuple(self.full_fingerprints))
        object.__setattr__(self, "control_fingerprints",
                           tuple(self.control_fingerprints))
        object.__setattr__(self, "seal", None)

    @property
    def lower(self) -> float:
        """The interval's lower endpoint — the underpowered test's subject."""
        return self.interval[0]

    @property
    def estimable_fraction(self) -> float:
        """Share of drawn replicates in which both exponents were estimable."""
        return self.estimable_replicates / self.replicates

    @property
    def is_adequate(self) -> bool:
        """Is this contrast big enough and clean enough to decide anything?

        The same frozen floors the exponent's own resampling record is held to:
        :data:`~.analysis.MIN_RESAMPLE_REPLICATES` draws and an estimable
        fraction of at least
        :data:`~.analysis.MIN_VALID_RESAMPLE_FRACTION`, with a finite interval.
        """
        return (self.replicates >= MIN_RESAMPLE_REPLICATES
                and self.estimable_fraction >= MIN_VALID_RESAMPLE_FRACTION
                and all(math.isfinite(value) for value in self.interval))


def paired_exponent_contrast(full_sweep, control_sweep) -> ContrastReport:
    """Bootstrap ``exponent(full) - exponent(control)`` on shared trial IDs.

    Takes **sweeps**, not designs.  Both are re-opened and re-hashed, both
    designs are built from their own ledgers, and the resulting record carries
    the two design digests, the coupling grid, the shared
    common-random-number identity and both fingerprint tuples — so a decision
    can require that this contrast is about the runs it is deciding on.

    The two sweeps must share their CRN identity, which is what a paired
    full/width-only pair has by construction, and must be the two different
    models.  One draw of identities is applied to both designs in every
    replicate, so the pairing survives the resampling instead of being averaged
    away by it, and a replicate counts only when both exponents are estimable.
    """
    # Locators in, disk objects out — and each must be an exact RunSweep, so
    # that neither the read nor the type is the caller's to choose.
    full_sweep = _fresh_sweep(full_sweep, "full_sweep")
    control_sweep = _fresh_sweep(control_sweep, "control_sweep")

    if full_sweep.model != "full":
        raise ValueError(f"the full sweep records "
                         f"model={full_sweep.model!r}")
    if control_sweep.model != "width_only":
        raise ValueError(f"the control sweep records "
                         f"model={control_sweep.model!r}")
    if full_sweep.crn_identity != control_sweep.crn_identity:
        raise ValueError(
            "the two sweeps come from different common-random-number "
            "hierarchies, so they are not the paired run this contrast is "
            "defined for")
    if full_sweep.couplings != control_sweep.couplings:
        raise ValueError(
            f"the full model was swept over {list(full_sweep.couplings)} and "
            f"the control over {list(control_sweep.couplings)}; a paired "
            f"contrast compares the same coupling grid")

    full_design = _design_from_fresh_sweep(full_sweep)
    control_design = _design_from_fresh_sweep(control_sweep)
    keys = full_design.paired_keys
    if control_design.paired_keys != keys:
        raise ValueError(
            "the two designs hold different master trial identities, so they "
            "are not the paired run this contrast is defined for")

    # Frozen count, digest-derived seed: a contrast that could be rerolled is
    # a contrast that can be selected.
    count = RESAMPLE_REPLICATES
    rng = np.random.default_rng(_authoritative_seed(full_design.digest))
    full_outcomes = [np.array(c.outcomes, dtype=bool)
                     for c in full_design.cells]
    control_outcomes = [np.array(c.outcomes, dtype=bool)
                        for c in control_design.cells]
    size = len(keys)
    namespace = full_design.cells[0].namespace

    differences = []
    for index in range(count):
        draw = rng.integers(0, size, size=size)
        replicate_ids = tuple(
            MasterTrial(f"contrast:{index}:{namespace}", slot)
            for slot in range(size))
        reports = []
        for design, outcomes in ((full_design, full_outcomes),
                                 (control_design, control_outcomes)):
            cells = tuple(
                CouplingCell(cell.coupling, replicate_ids,
                             tuple(bool(v) for v in values[draw]))
                for cell, values in zip(design.cells, outcomes))
            reports.append(fit_exponent(ScalingDesign(cells, "master_trial")))
        if all(_estimable(report) for report in reports):
            differences.append(reports[0].exponent - reports[1].exponent)

    values = np.array(differences, dtype=float)
    if values.size:
        mean = float(np.mean(values))
        spread = float(np.std(values, ddof=1)) if values.size > 1 else 0.0
        interval = (float(np.percentile(values, 2.5)),
                    float(np.percentile(values, 97.5)))
    else:
        mean = spread = float("nan")
        interval = (float("nan"), float("nan"))
    return ContrastReport(
        count, int(values.size), mean, spread, interval,
        full_design.digest, control_design.digest, full_sweep.couplings,
        full_sweep.crn_identity,
        tuple(run.fingerprint for run in full_sweep.runs),
        tuple(run.fingerprint for run in control_sweep.runs), seal=_SEAL)


@dataclass(frozen=True)
class ControlEvidence:
    """One calibration, bound to the sweep, role, data and status it came from.

    An independent review found the earlier version checking couplings and role
    but neither the ledger counts nor the outcome matrix: it accepted "an
    honest synthetic calibration with counts ``[(54,120), (66,120), (76,120),
    (87,120)]`` as evidence for a closed sweep holding ``[(2,6), (4,6), (5,6),
    (6,6)]``".  It also carried no numerical gate, "so a causal decision cannot
    block a diagnostic-only experiment".

    Both are closed by construction: :func:`control_evidence` *computes* the
    calibration from the sweep's own ledgers, and the record carries the design
    digest — over identities and outcomes, not totals — together with the
    sweep's status blockers.
    """

    role: str
    calibrated: CalibratedExponent
    crn_identity: str
    experiment_family: str
    model: str
    clocks: int
    detunings: tuple
    couplings: tuple
    fingerprints: tuple
    design_digest: str
    status_blockers: tuple
    seal: object = None

    def __post_init__(self):
        _require_seal(self.seal, "ControlEvidence")
        if self.role not in CONTROL_ROLES:
            raise ValueError(f"role must be one of {list(CONTROL_ROLES)}, got "
                             f"{self.role!r}")
        require_instance(self.calibrated, CalibratedExponent, "calibrated")
        for name in ("crn_identity", "experiment_family", "model",
                     "design_digest"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, str) or not \
                    value:
                raise TypeError(f"{name} must be a non-empty string, got "
                                f"{value!r}")
        if self.calibrated.design_digest != self.design_digest:
            raise ValueError(
                "this evidence's calibration was computed from different "
                "observations than the design it names")
        clocks = require_positive_int(self.clocks, "clocks", minimum=1)
        detunings = tuple(require_finite(d, f"detunings[{i}]")
                          for i, d in enumerate(self.detunings))
        if len(detunings) != clocks:
            raise ValueError(f"{len(detunings)} detunings for {clocks} clocks")
        couplings = tuple(require_finite(k, f"couplings[{i}]")
                          for i, k in enumerate(self.couplings))
        if not couplings:
            raise ValueError("a control covers at least one coupling")
        for position, entry in enumerate(tuple(self.fingerprints)):
            require_instance(entry, LedgerFingerprint,
                             f"fingerprints[{position}]")
        if self.role == "central":
            if clocks != 1 or detunings != (0.0,):
                raise ValueError(
                    f"the central-clock control is the ONE zero-detuning "
                    f"clock, in which tongue widening is absent; this sweep "
                    f"has {clocks} clock(s) at {list(detunings)}")
            if self.model != "full":
                raise ValueError(f"the central control runs the full Adler "
                                 f"model, not {self.model!r}")
        if self.role == "width_only" and self.model != "width_only":
            raise ValueError(
                f"the width-only control must be the fixed-contraction model; "
                f"this sweep records model={self.model!r}")
        if self.role == "full":
            if self.model != "full":
                raise ValueError(f"the full model role requires model='full', "
                                 f"got {self.model!r}")
            if clocks < 2:
                raise ValueError(
                    f"the full model role is the finite competing population; "
                    f"{clocks} clock is the central control, which has its own "
                    f"role")
        object.__setattr__(self, "detunings", detunings)
        object.__setattr__(self, "couplings", couplings)
        object.__setattr__(self, "fingerprints", tuple(self.fingerprints))
        object.__setattr__(self, "status_blockers",
                           tuple(self.status_blockers))
        object.__setattr__(self, "seal", None)


def control_evidence(role, sweep) -> ControlEvidence:
    """Compute a control's calibration from its own verified sweep.

    **The calibration is not an argument.**  The sweep is re-opened and
    re-hashed, its design is built from its own ledgers — binding the exact
    master-trial identities and outcomes — and the calibration is computed
    here.  The role is then required to be true of that sweep: ``central`` is
    the ONE zero-detuning clock under the full model, ``width_only`` is the
    fixed-contraction model at one positive rate, ``full`` is the multi-clock
    population under the full model.
    """
    require_exact_type(sweep, RunSweep, "sweep")
    if isinstance(role, bool) or not isinstance(role, str):
        raise TypeError(f"role must be a string, got {type(role).__name__} "
                        f"{role!r}")
    sweep = _fresh_sweep(sweep)     # locator in, fresh disk object out
    design = _design_from_fresh_sweep(sweep)
    calibrated = calibrate_exponent(design)
    manifest = sweep.runs[0].manifest
    return ControlEvidence(
        role=role, calibrated=calibrated, crn_identity=sweep.crn_identity,
        experiment_family=sweep.experiment_family,
        model=sweep.model, clocks=manifest.clocks,
        detunings=tuple(manifest.value("detunings")),
        couplings=sweep.couplings,
        fingerprints=tuple(run.fingerprint for run in sweep.runs),
        design_digest=design.digest, status_blockers=sweep.status_blockers,
        seal=_SEAL)


@dataclass(frozen=True)
class HoldoutEvidence:
    """The frozen out-of-sample comparison, computed rather than asserted.

    ``source`` is ``"computed"`` or ``"missing"`` and there is no caller-
    writable completeness bit.  An independent review found ``complete`` and
    both scores settable, so "missing holdout evidence is indistinguishable
    from evidence that the control was worse" — and used a hand-built one, with
    fingerprints belonging to no control in the decision, to manufacture a
    ``not_blocked``.

    :func:`holdout_evidence` is the only producer of a computed record and
    :func:`missing_holdout` the only producer of a missing one, and
    :func:`causal_decision` requires the fingerprints to belong to the
    evidence it is deciding about.
    """

    source: str
    reference_fingerprint: object
    full_fingerprint: object
    control_fingerprint: object
    segments: tuple
    bins_per_side: int
    full_score: object
    control_score: object
    reason: str = ""
    seal: object = None

    def __post_init__(self):
        _require_seal(self.seal, "HoldoutEvidence")
        if self.source not in ("computed", "missing"):
            raise ValueError(f"source must be 'computed' or 'missing', got "
                             f"{self.source!r}")
        if self.source == "missing":
            if not self.reason:
                raise ValueError("a missing holdout states why it is missing")
        else:
            for name in ("reference_fingerprint", "full_fingerprint",
                         "control_fingerprint"):
                require_instance(getattr(self, name), LedgerFingerprint, name)
            for name in ("full_score", "control_score"):
                require_instance(getattr(self, name), ComparisonScore, name)
            segments = tuple(self.segments)
            if segments != FROZEN_SEGMENTS:
                raise ValueError(f"holdout segments are the frozen "
                                 f"{list(FROZEN_SEGMENTS)}, got "
                                 f"{list(segments)}")
            require_positive_int(self.bins_per_side, "bins_per_side",
                                 minimum=1)
            if self.full_fingerprint == self.control_fingerprint:
                raise ValueError("the holdout scored one run against itself "
                                 "as both the full model and the control")
        object.__setattr__(self, "segments", tuple(self.segments))
        object.__setattr__(self, "seal", None)

    @property
    def complete(self) -> bool:
        """Was this holdout computed, or is it a declared absence?"""
        return self.source == "computed"

    @property
    def control_not_worse(self) -> bool:
        """Did the control match or outperform the full model?

        ``True`` for a missing record: an absence is not evidence that the
        control was worse, and reading it as one is the defect this type
        exists to close.
        """
        if not self.complete:
            return True
        return control_defeats(self.full_score, self.control_score)


def missing_holdout(reason) -> HoldoutEvidence:
    """A declared absence of holdout evidence, which blocks."""
    if isinstance(reason, bool) or not isinstance(reason, str) or not reason:
        raise TypeError("reason must be a non-empty string saying why the "
                        "holdout is missing")
    return HoldoutEvidence("missing", None, None, None, FROZEN_SEGMENTS, 1,
                           None, None, reason, seal=_SEAL)


def holdout_evidence(reference, full, control, bins_per_side):
    """Compute the frozen rising/falling and winner-mismatch comparison.

    All three arguments are :class:`ReadOnlyRun` objects.  They are re-hashed,
    required to agree on every frozen field outside
    :data:`EXPERIMENT_VARYING_FIELDS`, and required to play their declared
    roles: ``full`` runs the full Adler model and ``control`` the width-only
    one.

    The binning is :func:`frozen_edges` with the pulse centre on an edge, so
    the rising/falling boundary cannot move; both candidates are scored against
    the *same* reference run with the *same* edges.
    """
    fresh = {}
    for name, run in (("reference", reference), ("full", full),
                      ("control", control)):
        if type(run) is not ReadOnlyRun:
            raise TypeError(f"{name} must be exactly a ReadOnlyRun opened "
                            f"through open_for_comparison, got "
                            f"{type(run).__name__}")
        # Unbound and trusted: the exactness above is what refuses a subclass,
        # and this is what keeps the read this package's read regardless.
        fresh[name] = ReadOnlyRun.reopen(run)   # locator in, disk object out
    reference, full, control = fresh["reference"], fresh["full"], \
        fresh["control"]
    count = require_positive_int(bins_per_side, "bins_per_side", minimum=1)

    shared = tuple(name for name in MANIFEST_KEYS
                   if name not in EXPERIMENT_VARYING_FIELDS)
    reference_fields = _canonical_fields(reference.manifest, shared)
    for name, run in (("full", full), ("control", control)):
        if _canonical_fields(run.manifest, shared) != reference_fields:
            raise ValueError(
                f"the {name} run and the reference disagree about the frozen "
                f"fields outside {list(EXPERIMENT_VARYING_FIELDS)}; a holdout "
                f"compares one experiment's runs, not two experiments'")
    if full.manifest.value("model") != "full":
        raise ValueError(f"the full run records "
                         f"model={full.manifest.value('model')!r}")
    if control.manifest.value("model") != "width_only":
        raise ValueError(f"the width-only control records "
                         f"model={control.manifest.value('model')!r}")

    edges = frozen_edges(reference.manifest, count)
    splits, winners = {}, {}
    for name, run in (("reference", reference), ("full", full),
                      ("control", control)):
        table = _observables.risk_table(run.ledger, edges)
        splits[name] = rising_falling_split(run.manifest, table)
        winners[name] = _observables.winner_mismatch(run.ledger)

    full_score = discrepancy_score(splits["reference"], splits["full"],
                                   winners["reference"], winners["full"])
    control_score = discrepancy_score(splits["reference"], splits["control"],
                                      winners["reference"], winners["control"])
    return HoldoutEvidence(
        "computed", reference.fingerprint, full.fingerprint,
        control.fingerprint, FROZEN_SEGMENTS, count, full_score,
        control_score, "", seal=_SEAL)


@dataclass(frozen=True)
class CausalDecision:
    """Whether the width-times-rate attribution survives its negative controls.

    Sealed, and carrying the digests of the evidence it was made from.  An
    independent review constructed ``CausalDecision("not_blocked", (), 1.0,
    0.25)`` directly; that is now impossible, and a decision that says
    ``not_blocked`` names the exact designs and holdout it cleared.
    """

    attribution: str
    reasons: tuple
    contrast_lower: float
    minimum_difference: float
    full_digest: str
    central_digest: str
    control_digest: str
    holdout_source: str
    seal: object = None

    def __post_init__(self):
        _require_seal(self.seal, "CausalDecision")
        if self.attribution not in ATTRIBUTION_VERDICTS:
            raise ValueError(f"attribution must be one of "
                             f"{list(ATTRIBUTION_VERDICTS)}, got "
                             f"{self.attribution!r}")
        reasons = tuple(self.reasons)
        unknown = [r for r in reasons if r not in ATTRIBUTION_BLOCKERS]
        if unknown:
            raise ValueError(f"reasons {unknown} are outside "
                             f"{list(ATTRIBUTION_BLOCKERS)}")
        if (self.attribution == "blocked") != bool(reasons):
            raise ValueError(f"attribution {self.attribution!r} with reasons "
                             f"{list(reasons)}: a block names its blockers and "
                             f"an unblocked decision names none")
        object.__setattr__(self, "reasons", reasons)
        object.__setattr__(self, "seal", None)


def causal_decision(full, central, width_only, contrast, minimum_difference,
                    holdout) -> CausalDecision:
    """Apply the plan's frozen defeat rule to one full/control comparison.

    Every argument is rooted, and the rule **joins them** before it decides.
    An independent review combined public constructors into three unrelated
    controls, a caller-authored contrast and a holdout whose fingerprints
    belonged to neither, and received ``not_blocked ()`` over underlying
    calibrations that all said ``no_valid_exponent``.  So:

    * the contrast's ``full_digest`` and ``control_digest`` must be the full
      and width-only evidence's own design digests, and its coupling grid and
      CRN identity must be theirs;
    * the holdout's full and control fingerprints must **belong to** the full
      and width-only evidence's fingerprint tuples; and
    * every control must come from the same common-random-number hierarchy and
      the same coupling grid.

    Then the rule, in this order:

    1. a full model with **no valid exponent** cannot attribute anything;
    2. either control with **no valid exponent** blocks — a ``lack_of_fit`` or
       ``curvature`` refusal keeps a finite, often narrow interval, so the
       interval test alone cannot catch it;
    3. a **central single-clock** control whose paired interval still contains
       two shows that a quadratic crossover does not require tongue widening;
    4. a **width-only** control whose paired interval still contains two shows
       that the eligibility geometry, the thresholding and the censoring can
       manufacture the exponent without the second factor;
    5. a paired contrast that is inadequate, or whose lower endpoint does not
       exceed the frozen minimum mechanistic difference, is a
       causal-identification failure;
    6. missing holdout evidence blocks, and a control that matches or
       outperforms the full model on the frozen out-of-sample comparison
       defeats the attribution outright; and
    7. **any control whose own hashed manifests forbid a physical reading
       blocks.**  A diagnostic-only experiment cannot support a causal
       attribution, and the earlier evidence objects carried no status at all.
    """
    roles = {"full": full, "central": central, "width_only": width_only}
    for name, evidence in roles.items():
        require_instance(evidence, ControlEvidence, name)
        if evidence.role != name:
            raise ValueError(
                f"the {name!r} argument carries a control whose declared role "
                f"is {evidence.role!r}; roles are bound to the sweeps they "
                f"were fitted from and cannot be relabelled at the call site")
    require_instance(contrast, ContrastReport, "contrast")
    require_instance(holdout, HoldoutEvidence, "holdout")
    minimum = require_finite(minimum_difference, "minimum_difference")
    if minimum <= 0.0:
        raise ValueError(
            f"minimum_difference must be positive, got {minimum!r}; a "
            f"non-positive minimum mechanistic difference makes rule 5 "
            f"unfalsifiable")
    if full.crn_identity != width_only.crn_identity:
        raise ValueError(
            "the full model and the width-only control come from different "
            "common-random-number hierarchies, so they are not the paired run "
            "this contrast is defined for")
    # Every role must be the same *experiment*, population aside.  The central
    # control legitimately has a different population and therefore a different
    # CRN token, so equality of that token is the wrong check for it; the
    # experiment-family digest is the right one.  Without this a genuine
    # one-clock sweep whose manifest-frozen stream namespace had been replaced
    # by an unrelated hierarchy was accepted as the central control.
    family = full.experiment_family
    for name in ("central", "width_only"):
        if roles[name].experiment_family != family:
            raise ValueError(
                f"the {name} control belongs to a different experiment family "
                f"than the full model: they disagree about at least one of "
                f"{list(EXPERIMENT_FAMILY_FIELDS)}.  A control from another "
                f"bath, pulse, lock band, dwell, mesh or stream hierarchy is "
                f"not a control for this experiment, whatever its couplings "
                f"and counts look like")
        if roles[name].couplings != full.couplings:
            raise ValueError(
                f"the {name} control was swept over "
                f"{list(roles[name].couplings)} and the full model over "
                f"{list(full.couplings)}; a defeat rule compares one coupling "
                f"grid")
    # The contrast has to be about these two controls -- the same observations
    # *and* the same runs.  Deterministic rows can give two different runs the
    # same design digest, so a review built a contrast from a gate-rewritten
    # copy of the full sweep and it joined on digest alone.
    if contrast.full_digest != full.design_digest or \
            contrast.control_digest != width_only.design_digest:
        raise ValueError(
            "this contrast was computed from different observations than the "
            "full model and width-only control it is being used to compare; a "
            "caller-authored interval is not a paired contrast")
    if contrast.full_fingerprints != full.fingerprints:
        raise ValueError(
            f"the contrast was computed over {len(contrast.full_fingerprints)} "
            f"full-model runs that are not this evidence's "
            f"{len(full.fingerprints)}.  Two runs can hold identical rows and "
            f"still be different runs -- a copy with one manifest field "
            f"rewritten has the same design digest and a different "
            f"fingerprint -- so the ordered fingerprint tuples must match "
            f"exactly, not merely the digests")
    if contrast.control_fingerprints != width_only.fingerprints:
        raise ValueError(
            "the contrast's width-only runs are not this evidence's runs")
    if contrast.couplings != full.couplings:
        raise ValueError("the contrast covers a different coupling grid")
    if contrast.crn_identity != full.crn_identity:
        raise ValueError("the contrast comes from a different "
                         "common-random-number hierarchy")
    # The holdout has to be about these two sweeps.
    if holdout.complete:
        if holdout.full_fingerprint not in full.fingerprints:
            raise ValueError(
                "the holdout's full run is not one of the full model's own "
                "runs, so it is evidence about a different experiment")
        if holdout.control_fingerprint not in width_only.fingerprints:
            raise ValueError(
                "the holdout's control run is not one of the width-only "
                "control's own runs")
        if holdout.reference_fingerprint not in (full.fingerprints
                                                 + width_only.fingerprints):
            raise ValueError(
                "the holdout was scored against a reference run belonging to "
                "neither sweep, so both of its scores are about a comparison "
                "this decision cannot see")

    reasons = []
    if not full.calibrated.is_valid:
        reasons.append("full_model_no_valid_exponent")
    if not central.calibrated.is_valid:
        reasons.append("central_control_no_valid_exponent")
    if not width_only.calibrated.is_valid:
        reasons.append("width_control_no_valid_exponent")
    if central.calibrated.interval_contains(QUADRATIC_EXPONENT):
        reasons.append("central_control_compatible_with_quadratic")
    if width_only.calibrated.interval_contains(QUADRATIC_EXPONENT):
        reasons.append("width_control_compatible_with_quadratic")
    if not contrast.is_adequate:
        reasons.append("contrast_inadequate")
    lower = contrast.lower
    if not (math.isfinite(lower) and lower > minimum):
        reasons.append("contrast_underpowered")
    if not holdout.complete:
        reasons.append("holdout_evidence_missing")
    elif holdout.control_not_worse:
        reasons.append("control_matches_out_of_sample")
    if any(evidence.status_blockers for evidence in roles.values()):
        reasons.append("numerical_status_not_cleared")

    ordered = tuple(r for r in ATTRIBUTION_BLOCKERS if r in set(reasons))
    return CausalDecision(
        "blocked" if ordered else "not_blocked", ordered, float(lower),
        minimum, full.design_digest, central.design_digest,
        width_only.design_digest, holdout.source, seal=_SEAL)


@dataclass(frozen=True)
class ScalingVerdict:
    """What a fitted exponent may be read as, given the runs behind it.

    Produced only by :func:`scaling_verdict`, which recomputes the calibration
    from the re-opened sweep.  An independent review constructed one directly::

        ScalingVerdict("physical_verdict_permitted", "closed_runs", (),
                       ("production_cleared",), ("not-a-fingerprint",))
        # accepted

    That is now impossible twice over: the record is sealed, and under manifest
    schema v3 ``physical_verdict_permitted`` cannot be held at all.  See
    :data:`PRODUCTION_STATUS_SCHEMA` for why a rewritten free-text gate is not
    a production budget.
    """

    verdict: str
    reasons: tuple
    gates: tuple
    fingerprints: tuple
    design_digest: str
    status_blockers: tuple
    seal: object = None

    def __post_init__(self):
        _require_seal(self.seal, "ScalingVerdict")
        if self.verdict not in SCALING_VERDICTS:
            raise ValueError(f"verdict must be one of "
                             f"{list(SCALING_VERDICTS)}, got {self.verdict!r}")
        if self.verdict == "physical_verdict_permitted":
            raise ValueError(
                f"manifest schema v{MANIFEST_SCHEMA_VERSION} cannot express "
                f"a production-budget status, so no run it describes can carry "
                f"a physical permission.  A permission needs schema "
                f"v{PRODUCTION_STATUS_SCHEMA}; the vocabulary retains the word "
                f"so that the boundary is nameable, and this refusal is what "
                f"makes it fail closed")
        prints = tuple(self.fingerprints)
        for position, entry in enumerate(prints):
            require_instance(entry, LedgerFingerprint,
                             f"fingerprints[{position}]")
        gates = tuple(self.gates)
        if not gates or not prints or len(gates) != len(prints):
            raise ValueError(
                f"a scaling verdict carries one gate and one fingerprint per "
                f"contributing run, got {len(gates)} and {len(prints)}")
        if not isinstance(self.design_digest, str) or not self.design_digest:
            raise TypeError("design_digest must be a non-empty string")
        object.__setattr__(self, "reasons", tuple(self.reasons))
        object.__setattr__(self, "gates", gates)
        object.__setattr__(self, "fingerprints", prints)
        object.__setattr__(self, "status_blockers",
                           tuple(self.status_blockers))
        object.__setattr__(self, "seal", None)


@dataclass(frozen=True)
class MachineryReport:
    """The result of exercising the estimator on data that never came from a run.

    A *separate type*, deliberately.  Fitting declared counts is how this
    package shows the estimator works, and it is the right thing to do — but if
    it came back through the same door a closed run does, "the machinery
    recovered exponent two" and "this experiment measured exponent two" would
    be the same sentence with the same type on it.

    There is no verdict field.  There is a calibration, a description of what
    was exercised, and the statement that no run stands behind it.
    """

    calibrated: CalibratedExponent
    note: str
    seal: object = None

    def __post_init__(self):
        _require_seal(self.seal, "MachineryReport")
        require_instance(self.calibrated, CalibratedExponent, "calibrated")
        if isinstance(self.note, bool) or not isinstance(self.note, str) or \
                not self.note:
            raise TypeError("note must be a non-empty string describing what "
                            "was exercised")
        object.__setattr__(self, "seal", None)

    @property
    def is_physical(self) -> bool:
        """Always ``False``.  Present so the question has an answer."""
        return False

    @property
    def summary(self) -> str:
        """One line a reader cannot quote without the disclaimer."""
        return (f"synthetic machinery exercise ({self.note}): the calibrated "
                f"fit reports {self.calibrated.verdict}"
                + (f" {list(self.calibrated.reasons)}"
                   if self.calibrated.reasons else "")
                + "; no closed run stands behind this number and it is not a "
                  "physical coupling exponent")


def scaling_verdict(sweep) -> ScalingVerdict:
    """Recompute a sweep's calibration from its own ledgers, and rule on it.

    **The calibration is not an argument.**  An earlier signature took one, and
    an independent review built a ``CalibratedExponent`` whose published
    verdict and interval agreed with neither record it carried and handed it
    here, receiving ``physical_verdict_permitted`` over a sweep whose own
    records said ``curvature_unidentified, uninformative,
    inadequate_resampling``.  So this function now:

    1. re-opens and re-hashes every run in ``sweep``;
    2. builds the design from those ledgers with
       :func:`design_from_sweep`, binding the exact master-trial identities and
       the exact per-trial outcomes — not the aggregate counts, which two
       different outcome matrices can share while having different paired
       bootstraps;
    3. calibrates it with :func:`~.analysis.calibrate_exponent`; and
    4. reads every load-bearing status out of the hashed manifests.

    It returns ``machinery_only`` — always, at present, and not because the fit
    happens to fail.  Schema v3 cannot express a production-budget status, so
    :data:`PRODUCTION_STATUS_SCHEMA` blocks every run it describes.  The
    *reasons* are what discriminate: rewriting one manifest's ``numerical_gate``
    removes exactly the gate reason and leaves the schema and the contradictory
    moving-band and validation-note reasons standing.

    Exercising the estimator on data that never came from a run is a different
    activity with a different type: :func:`synthetic_machinery_verdict`.
    """
    # The supplied sweep is a locator, and an exact RunSweep: a subclass whose
    # ``reopen`` returned itself would make the read below an assertion.
    # Everything after it is computed from the object the read returns, which
    # was parsed from the bytes on disk moments ago.  There is no replicate
    # count and no seed: an authoritative decision that accepted either could
    # be rerolled until it said something else.
    sweep = _fresh_sweep(sweep)
    design = _design_from_fresh_sweep(sweep)
    calibrated = calibrate_exponent(design)

    reasons = []
    blockers = sweep.status_blockers
    if blockers:
        reasons.append(
            f"{len(blockers)} status condition(s) in the contributing runs' "
            f"own hashed manifests forbid a physical reading: {list(blockers)}")
    if not calibrated.is_valid:
        reasons.append(f"the calibrated fit reports no valid exponent: "
                       f"{list(calibrated.reasons)}")
    if not reasons:                                    # unreachable under v3
        reasons.append("no status condition blocked, but manifest schema v3 "
                       "cannot express a production budget")
    return ScalingVerdict(
        "machinery_only", tuple(reasons), sweep.numerical_gates,
        tuple(run.fingerprint for run in sweep.runs), design.digest, blockers,
        seal=_SEAL)


def synthetic_machinery_verdict(calibrated, note) -> MachineryReport:
    """The separate, clearly non-physical door for synthetic exercises."""
    require_instance(calibrated, CalibratedExponent, "calibrated")
    if isinstance(note, bool) or not isinstance(note, str):
        raise TypeError(f"note must be a string describing what was "
                        f"exercised, got {type(note).__name__} {note!r}")
    if not note:
        raise ValueError("note must not be empty; a synthetic result carries "
                         "the description of what it exercised")
    return MachineryReport(calibrated, note, seal=_SEAL)
