"""
experiments.py — ticket 07: measured feasibility, frozen power, the finite
production matrix proposal, and the pilot firewall.

This module is **orchestration and arithmetic**.  It measures what the real
streamed race costs, works out how many trials the frozen precision target
needs, decides whether the numerical evidence ticket 04 produced is small
enough against that precision, proposes one finite production matrix, and
builds the firewall a range-only pilot would have to run behind.

It sits **outside the raw import graph**, deliberately and in both directions:

* nothing under ``raw_runner.py`` imports this file, so the event generator
  cannot see a power calculation, a cost model, or a production proposal; and
* this file imports the raw runner, the ledger gate, the analysis layer and the
  numerical-validation vocabulary, because orchestrating them is what it is
  for.

What this module does **not** do
--------------------------------
It does not run the pilot.  It does not open pilot exponent-related output.  It
cannot construct a signed or production-cleared manifest: the only statuses a
proposal may carry are in :data:`PROPOSAL_STATUSES`, ``production_cleared`` is
not among them, and :data:`compare.PRODUCTION_STATUS_SCHEMA` is untouched — the
schema boundary that makes ``physical_verdict_permitted`` unreachable under
manifest schema v3 is ticket 06's, and nothing here widens it.

The proposal this module produces is exactly that: a **proposal**, at status
``checkpoint_blocked``, waiting on a named human approval that this file has no
way to forge.

Three vocabularies, and why each is closed
------------------------------------------
:data:`FEASIBILITY_VERDICTS`, :data:`BUDGET_VERDICTS` and
:data:`PROPOSAL_STATUSES` are closed tuples for the reason ``analysis.py``'s
refusal reasons are: a consumer branches on them, so a new one is a deliberate
edit here rather than a string invented at a call site.  In particular there is
no verdict spelled ``feasible_and_approved``.  Approval is not a computation.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
import shutil
import time
import tracemalloc
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from . import analysis as _analysis
from . import compare as _compare
from . import dynamics as _dynamics
from . import killed_diffusion as _oracle
from . import raw_config as _raw_config
from . import raw_ledger as _raw_ledger
from . import raw_runner as _raw_runner
from . import stochastic as _stochastic
from .validation import (require_exact_type, require_finite,
                         require_nonnegative, require_positive_int)

__all__ = [
    # Identity and vocabularies.
    "EXPERIMENTS_SCHEMA",
    "CONTENT_REVISION",
    "EXPERIMENT_NON_CLAIMS",
    "PROPOSAL_STATUSES",
    "FEASIBILITY_VERDICTS",
    "TICKET_08_DERIVATION_STATUSES",
    "TICKET_08_AUTHORIZATION_STATUSES",
    "CHECKPOINT_DECISION_CLAUSES",
    "CHECKPOINT_REQUESTED_DECISION",
    "BUDGET_VERDICTS",
    "BUDGET_BLOCKERS",
    "POWER_STAGES",
    "POWER_VERDICTS",
    "POST_RANGE_STATUSES",
    "POST_RANGE_MATRIX_NOT_DERIVED",
    "CONFIGURATION_LABELS",
    "PERMITTED_PILOT_FIELDS",
    "FORBIDDEN_PILOT_FIELDS",
    "PROPOSED_MANIFEST_KEYS",
    "PROPOSAL_FIELD_SOURCES",
    "SEAL_IS_NOT_THE_AUTHORITY",
    "DESIGN_DIGEST_KEYS",
    "SELECTION_REFUSALS",
    "ALTERNATIVE_LABELS",
    # Frozen numbers.
    "CONFIDENCE_LEVEL",
    "CONFIDENCE_Z",
    "COVERAGE_SIGMA",
    "ALLOWED_NUMERICAL_FRACTION",
    "RACE_PASSES",
    "CONTRAST_COVARIANCE_FLOOR",
    "MIN_PILOT_EVENTS",
    "MIN_PILOT_SURVIVORS",
    "MIN_SELECTED_CELLS",
    "MIN_ELIGIBLE_CLOCKS",
    "PILOT_NAMESPACE_PREFIX",
    "PRODUCTION_NAMESPACE_PREFIX",
    "QUARANTINE_DIRNAME",
    "REFERENCE_RATE_RULE",
    "SOURCE_FILES",
    # Hashing.
    "canonical_bytes",
    # Provenance.
    "SourceFingerprint",
    "source_fingerprint",
    # Feasibility.
    "BenchmarkFixture",
    "BenchmarkResult",
    "run_benchmark",
    # Power and precision.
    "ArmEnvelope",
    "SamplingTarget",
    "PowerEstimate",
    "power_estimate",
    "post_range_power",
    "cloglog_slope_variance",
    "RecoveryReport",
    "recover_exponent",
    # Numerical budget and the finite feasibility search.
    "NumericalEvidence",
    "NumericalBudget",
    "BudgetDisposition",
    "numerical_disposition",
    "FeasibilityOption",
    "SearchReport",
    "search_feasible_options",
    # The finite matrix and what it costs.
    "PrimaryMatrix",
    "CostLine",
    "ResourceEstimate",
    "matrix_resource_estimate",
    "pilot_cost_line",
    # The pilot firewall and the frozen pilot plan.
    "row_width_bound",
    "PilotFirewall",
    "PILOT_ROLES",
    "PILOT_SELECTION_RULE",
    "PILOT_PURPOSES",
    "FIXTURE_PREFIX_MARKER",
    "PRODUCTION_PILOT_SPEC",
    "SHARED_PHYSICAL_FIELDS",
    "shared_physical_identity",
    "require_shared_physical_identity",
    "PilotPlan",
    "pilot_plan",
    "canonical_pilot_plan",
    "require_canonical_pilot_plan",
    "require_plan_firewall",
    "require_authoritative_plan",
    "QuarantinedPilotRun",
    "write_pilot_run",
    "PilotCellCounts",
    "pilot_counts",
    "RangeSelection",
    "PostRangePlan",
    "post_range_plan",
    "PreviewRangeSelection",
    "preview_range",
    "select_range",
    "reference_rate",
    # The proposal and the handoff.
    "ProposedProductionManifest",
    "proposed_production_manifest",
    "require_authoritative_proposal",
    "require_unchanged_proposal",
    "CheckpointAlternative",
    "SUPPORT_CLASSES",
    "ALTERNATIVE_CONTRACT",
    "intended_configuration_option",
    "narrowed_claim_option",
    "checkpoint_alternatives",
    "pilot_cost_from_plan",
    "Checkpoint1Handoff",
    "checkpoint_one_handoff",
]

#: The identity of everything this module hashes.  A record built under one
#: schema and read under another is a different record, and the tag is inside
#: every digest so that it cannot be read as the other one.
EXPERIMENTS_SCHEMA = "dk-experiments/v2"

#: The content revision of this layer's frozen statements.  Bumped whenever a
#: recorded sentence changes meaning rather than wording, so that a reader who
#: has an older record can tell that the two disagree on purpose.  Revision 2
#: is the independent review's closure round: the benchmark work model changed
#: from nominal to walked clock-steps, the coverage sigma became authoritative,
#: the feasibility conclusion was narrowed to the options actually searched,
#: and the pilot ledger moved behind a quarantine.
CONTENT_REVISION = 2

#: Repeated here rather than imported, because a reader who opens this file to
#: find out what a production proposal would mean should meet them before the
#: first number.
EXPERIMENT_NON_CLAIMS = (
    "no Born-rule result",
    "no detector click, absorption or measurement outcome",
    "no unique actuality or global one-world exclusivity",
    "no microscopic origin for the phase noise",
    "no two-channel competition and no second outcome channel",
    "no production-cleared coupling exponent",
    "no approval: a proposal is not a decision, and this module cannot make "
    "one",
)

#: What a proposed production manifest may say about itself.  There is no
#: ``production_cleared`` and no ``signed``: the first belongs to the manifest
#: schema ``compare.py`` gates on and is unreachable under schema v3, and the
#: second is a human act.
#: There is exactly one, and that is the finding.  An earlier revision offered
#: ``checkpoint_one_approved`` beside it, and an independent review pointed out
#: that a public type able to *hold* an approved status is a public type able to
#: manufacture one — the refusal that was supposed to prevent it lived in a
#: factory a caller could go around.  Approval is a later ticket's schema and a
#: human act, and neither of those is a string this module can spell.
PROPOSAL_STATUSES = ("checkpoint_blocked",)

#: The outcome of asking whether a finite production matrix exists at all.
#: Note what is *not* here: ``no_feasible_matrix``, unqualified.  An
#: independent review was right that the evidence supports only a statement
#: about the options that were actually examined.  A claim that no timestep,
#: sample size or configuration could ever work needs a search or a bound, and
#: :func:`search_feasible_options` performs a *finite declared* search — so the
#: verdict it can support says exactly that and no more.
FEASIBILITY_VERDICTS = ("feasible_matrix", "no_feasible_matrix_among_searched",
                        "unresolved")

#: Where a power estimate sits in the pilot sequence.  A pre-pilot number is
#: **provisional**: it assumes a coupling range and a response curve that the
#: count-only pilot has not yet chosen, and it assumes every arm carries the
#: same information, which the two negative controls need not.
POWER_STAGES = ("provisional_pre_pilot", "post_range_join")

#: What a power calculation can conclude.  ``no_powered_matrix`` is a result:
#: a selected range can be information-bearing under the count rule and still
#: fail to support the causal contrast the plan requires.
POWER_VERDICTS = ("powered", "no_powered_matrix")

#: Whether a piece of numerical evidence, or a searched option, refers to the
#: configuration the matrix proposes or to the reference cell ticket 04
#: actually ran.  Kept as a closed label rather than a boolean because
#: "unknown" is a third state and silently reading it as either is how a
#: reference cell becomes a production measurement.
CONFIGURATION_LABELS = ("intended", "non_intended", "unmeasured")

#: The outcome of comparing ticket 04's measured numerical evidence with the
#: statistical uncertainty this ticket plans for.  ``numerical_no_result`` is
#: ticket 04's own word and is carried through unchanged; ``unresolved`` is the
#: weaker statement that the evidence exists but was never measured at the
#: configuration being proposed.
BUDGET_VERDICTS = ("satisfied", "numerical_no_result", "unresolved")

#: Why a budget disposition blocks.  Closed for the same reason the analysis
#: layer's refusal reasons are closed.
BUDGET_BLOCKERS = (
    "moving_band_numerical_no_result",
    "endpoint_envelope_exceeds_allowance",
    "audit_envelope_exceeds_allowance",
    "no_evidence_at_intended_configuration",
    "probability_window_empty",
    "time_bounds_fail_at_every_trial_count",
    "evidence_missing",
)

#: The only fields a range-only pilot is allowed to expose before the
#: production manifest is signed.  Counts, and the coupling they were counted
#: at.  Nothing else.
PERMITTED_PILOT_FIELDS = ("coupling", "trials", "committed", "survivors")

#: Everything a pilot may not be asked for.  This is not a denylist standing in
#: for a missing allowlist — :data:`PERMITTED_PILOT_FIELDS` is the allowlist and
#: is what :class:`PilotCellCounts` structurally has.  This tuple exists so that
#: a *mutation test* has something to iterate over, and so that a reader can see
#: the specific inspections the plan forbids named rather than implied.
FORBIDDEN_PILOT_FIELDS = (
    "exponent", "exponent_interval", "exponent_stderr", "curvature",
    "amplitude", "deviance", "lack_of_fit", "likelihood", "survival",
    "survival_shape", "hazard", "cumulative_hazard", "commit_time",
    "commit_time_quantile", "dwell", "dwell_resets", "band_entries",
    "noise", "phase_diffusion", "pulse", "population", "state_machine",
    "winner_mismatch", "exclusions", "co_completion", "detunings",
    "final_phase", "eligible_time", "exposure_time",
)

#: Why a count-only range selection can refuse.  A refusal is a result.
SELECTION_REFUSALS = ("no_eligible_cells", "run_too_short",
                      "not_contiguous_enough")

#: What a post-range join can conclude.  Derived from the selection and the
#: re-sized power estimate, never declared: a refused range, a selected range
#: that cannot be powered, and a selected range that can.
POST_RANGE_STATUSES = ("post_range_refused", "post_range_no_powered_matrix",
                       "post_range_powered")

#: The named refusal a post-range proposal gets until its target, matrix and
#: resource estimate have been re-derived from the post-range state.  Deriving
#: them is ticket-08 work; serializing the pre-pilot ones beside a post-range
#: power estimate would put two conflicting designs in one record, which is
#: what an independent review found this factory doing.
POST_RANGE_MATRIX_NOT_DERIVED = "post_range_matrix_not_derived"

#: What a checkpoint may say about ticket 08's matrix derivation.  Closed, and
#: **derived** rather than declared: this ticket does not derive a post-range
#: target, matrix or resource estimate, so the only value it can produce is
#: ``not_started``.  It says nothing about ticket 08 being authorized -- it
#: says the work has not begun, which is what a reader of a checkpoint
#: otherwise has to infer from an absence.
TICKET_08_DERIVATION_STATUSES = ("not_started",)

#: What a checkpoint may say about ticket 08's **authorization**.  Closed at one
#: member, derived, and with no factory parameter.
#:
#: An independent review pointed out that the summary carried
#: ``ticket_08_derivation_status`` and a null post-range proposal but nothing a
#: reader could *observe* about authorization: the non-authorization statement
#: lived only inside the presented sentence.  "Not started" and "not
#: authorized" are different facts -- work can be unauthorized and also not
#: begun, and neither implies the other -- so the second is now a field of its
#: own rather than a clause someone has to parse out of prose.
TICKET_08_AUTHORIZATION_STATUSES = ("not_authorized",)

#: The three things a checkpoint-1 presentation must say, each on its own, so
#: that each can be removed and the removal caught.
#:
#: An independent review handed :func:`checkpoint_one_handoff` coherent
#: authoritative components and the caller text ``"approve the presented
#: choice"``.  The handoff still carried
#: ``ticket_08_derivation_status = "not_started"`` and no post-range proposal
#: -- the *fields* were right -- and the sentence a reader would actually be
#: shown said none of it.  A derived field nobody reads is not a presentation.
CHECKPOINT_DECISION_CLAUSES = (
    "There is no post-range production proposal.",
    "Ticket 08 matrix derivation has not begun.",
    "Ticket 08 is not authorized.",
)

#: The **only** decision text a checkpoint-1 handoff can present.  Derived from
#: one constant, required by :class:`Checkpoint1Handoff` to be exactly this,
#: and with no parameter through which a caller could supply another.
CHECKPOINT_REQUESTED_DECISION = (
    "Choose between gathering intended-configuration numerical evidence and "
    "narrowing the claim.  No pilot has been run, no pilot exponent-related "
    "output has been opened, and no production manifest has been signed.  "
    + "  ".join(CHECKPOINT_DECISION_CLAUSES)
    + "  The last two state what exists and what has not been done; neither "
      "is a request to authorize ticket 08.")

#: The two things John can be asked to choose between at checkpoint 1: one
#: **blocked and explicitly unpriced** option and one **machinery-only**
#: diagnostic.  They are not two costed alternatives; Option A carries no cost,
#: no range and no bound.
#: Neither is "run the pilot": with the numerical disposition blocking, a range
#: chosen now feeds a sweep that cannot be interpreted whichever range it picks.
ALTERNATIVE_LABELS = ("intended_configuration_validation",
                      "narrowed_claim")

#: What an alternative may say it supports, as a closed enum.  There is **no**
#: member that licenses a physical, causal or scaling claim, so an independent
#: review's forged ``supports = "a physical scaling claim"`` is not a value
#: this vocabulary can spell rather than one a guard happens to catch.
SUPPORT_CLASSES = ("unpriced_no_claim", "machinery_and_diagnostic_only")

#: ``label -> (required support class, must it be priced)``.  Both alternatives
#: are fully determined by their label: the intended-configuration outline is
#: unpriced and supports nothing, and the narrowed diagnostic is priced and
#: supports a machinery statement.  A record disagreeing with its own label is
#: refused at construction.
ALTERNATIVE_CONTRACT = {
    "intended_configuration_validation": ("unpriced_no_claim", False),
    "narrowed_claim": ("machinery_and_diagnostic_only", True),
}

#: The complete key set of a proposed production manifest.  Closed: a proposal
#: missing one of these is refused, and a proposal carrying anything else is
#: refused, so "every required frozen field is enumerated" is a property of the
#: type rather than of the caller who filled it in.
PROPOSED_MANIFEST_KEYS = (
    # Identity and status.
    "schema", "content_revision", "status", "status_reason", "proposal_label",
    # Provenance.
    "source_digest", "python_version", "numpy_version", "platform",
    "ledger_schema_version", "manifest_schema_version",
    "production_status_schema",
    # Physical inputs.
    "peak_couplings", "pulse_duration", "pulse_centre", "phase_diffusion",
    "lock_tolerance", "dwell_time", "clocks", "arm_clocks", "grid_half_width",
    "grid_spacing", "grid_parity", "grid_origin", "grid_support",
    "minimum_eligible_clocks", "eligible_at_weakest",
    "discrete_staircase_fallback",
    # Numerical inputs.
    "timestep", "timestep_ladder", "refinement_ratio",
    "refinement_couplings", "refinement_arm",
    "numerical_gate", "endpoint_status", "moving_band_status",
    "numerical_budget_digest", "coverage_sigma", "allowed_fraction",
    "numerical_verdict", "numerical_blockers", "evidence_digest",
    "search_digest", "search_verdict",
    # Statistical inputs.
    "trials_per_cell", "arm_trials", "arms", "shadow_fraction",
    "resampling_unit",
    "fit_model", "confidence_level", "exponent_half_width_target",
    "minimum_causal_contrast", "contrast_covariance_floor",
    "measured_arm_latent_correlation", "recovery_digest",
    "curvature_test", "lack_of_fit_test", "inclusion_rule",
    "sampling_target_digest", "power_digest", "power_stage", "power_verdict",
    # Streams, folders and the firewall.
    "stream_namespace", "pilot_stream_namespace", "run_directory_prefix",
    "pilot_run_directory_prefix", "pilot_quarantine_directory",
    "pilot_trials_eligible", "firewall_digest", "pilot_plan_digest",
    "pilot_candidate_couplings", "pilot_trials_per_cell",
    "reference_rate_rule",
    "fixed_contraction_rate", "selection_digest", "post_range_digest",
    "post_range_status", "pilot_ledger_fingerprints",
    # What it costs, the design freeze, and what it is not.
    "benchmark_digest", "matrix_digest", "resource_digest", "design_digest",
    "stop_rule", "failure_rule", "non_claims",
)

#: The subset of :data:`PROPOSED_MANIFEST_KEYS` that is the **scientific
#: freeze**, and contains no measurement.
#:
#: An independent review found that the proposal digest moved when the same
#: source was re-run, because benchmark wall time sat inside it.  That is right
#: for a *record of what this machine did* and wrong for *the thing a reviewer
#: freezes*.  So there are two digests now: ``design_digest`` over these keys
#: only, which is reproducible on any machine from the same declared inputs,
#: and the whole-proposal digest, which additionally binds the measurement and
#: therefore does move.
DESIGN_DIGEST_KEYS = tuple(
    key for key in PROPOSED_MANIFEST_KEYS
    if key not in ("design_digest", "status_reason", "proposal_label",
                   "benchmark_digest", "resource_digest"))

#: The five keys deliberately outside the reproducible design freeze, and why.
#:
#: ``design_digest`` is the freeze itself and cannot contain its own value.
#: ``status_reason`` and ``proposal_label`` are free-text descriptions of the
#: record rather than statements about the experiment.
#: ``benchmark_digest`` and ``resource_digest`` carry **wall time**: they are
#: measurements of a machine, so binding them would make the scientific freeze
#: move whenever the same source was re-run on the same inputs.
#:
#: That is a two-level distinction and it is deliberate:
#:
#: * :attr:`ProposedProductionManifest.design_digest` is the **reproducible
#:   scientific freeze** — source, environment, firewall, selection, reference
#:   rate, numerical budget and evidence, power, matrix, and every
#:   measurement-free resource assumption.  The same source and the same
#:   declared inputs reproduce it.
#: * :attr:`ProposedProductionManifest.digest` is the **measured proposal
#:   envelope**: the design freeze *and* the benchmark and resource records
#:   this machine produced.  It moves with the machine, which is correct for a
#:   record of what a machine did.
#:
#: Semantic authority is not sacrificed to make the first one stable.  A
#: substituted measurement cannot enter a proposal at all, because the factory
#: recomputes the construction-time joins — the resource estimate must cite the
#: matrix and the benchmark it was actually scaled from, and the power estimate
#: must cite the target it was sized from — and refuses otherwise.
DESIGN_DIGEST_EXCLUSIONS = ("design_digest", "status_reason",
                            "proposal_label", "benchmark_digest",
                            "resource_digest")

#: Where every serialized proposal field comes from, as a durable table.
#:
#: ``(authoritative component, classification)``.  ``classification`` is
#: ``derived`` when :func:`require_authoritative_proposal` rebuilds the value
#: from that component and compares it, and ``free_text`` when the value is a
#: human description the rebuild takes from the record itself.  Only two fields
#: are free text, and neither makes a statement about the experiment.
#:
#: Membership of the design freeze is not repeated here: it is
#: :data:`DESIGN_DIGEST_KEYS`, and the suite checks the two agree rather than
#: letting a second list drift.  Every field is inside the whole-proposal
#: digest without exception.
PROPOSAL_FIELD_SOURCES = {
    "schema": ("experiments module", "derived"),
    "content_revision": ("experiments module", "derived"),
    "status": ("caller, from PROPOSAL_STATUSES", "derived"),
    "status_reason": ("human description", "free_text"),
    "proposal_label": ("human description", "free_text"),
    "source_digest": ("SourceFingerprint", "derived"),
    "python_version": ("SourceFingerprint", "derived"),
    "numpy_version": ("SourceFingerprint", "derived"),
    "platform": ("SourceFingerprint", "derived"),
    "ledger_schema_version": ("raw_ledger", "derived"),
    "manifest_schema_version": ("raw_ledger", "derived"),
    "production_status_schema": ("compare", "derived"),
    "peak_couplings": ("PrimaryMatrix", "derived"),
    "pulse_duration": ("PrimaryMatrix", "derived"),
    "pulse_centre": ("PrimaryMatrix", "derived"),
    "phase_diffusion": ("PrimaryMatrix", "derived"),
    "lock_tolerance": ("PrimaryMatrix", "derived"),
    "dwell_time": ("PrimaryMatrix", "derived"),
    "clocks": ("PrimaryMatrix", "derived"),
    "arm_clocks": ("PrimaryMatrix", "derived"),
    "grid_half_width": ("PrimaryMatrix", "derived"),
    "grid_spacing": ("PrimaryMatrix", "derived"),
    "grid_parity": ("PrimaryMatrix", "derived"),
    "grid_origin": ("PrimaryMatrix", "derived"),
    "grid_support": ("PrimaryMatrix", "derived"),
    "minimum_eligible_clocks": ("PrimaryMatrix", "derived"),
    "eligible_at_weakest": ("PrimaryMatrix", "derived"),
    "discrete_staircase_fallback": ("PrimaryMatrix", "derived"),
    "timestep": ("PrimaryMatrix", "derived"),
    "timestep_ladder": ("PrimaryMatrix", "derived"),
    "refinement_ratio": ("PrimaryMatrix", "derived"),
    "refinement_couplings": ("PrimaryMatrix", "derived"),
    "refinement_arm": ("PrimaryMatrix", "derived"),
    "numerical_gate": ("raw_runner", "derived"),
    "endpoint_status": ("raw_runner", "derived"),
    "moving_band_status": ("raw_runner", "derived"),
    "numerical_budget_digest": ("BudgetDisposition", "derived"),
    "coverage_sigma": ("experiments module", "derived"),
    "allowed_fraction": ("experiments module", "derived"),
    "numerical_verdict": ("BudgetDisposition", "derived"),
    "numerical_blockers": ("BudgetDisposition", "derived"),
    "evidence_digest": ("BudgetDisposition", "derived"),
    "search_digest": ("SearchReport", "derived"),
    "search_verdict": ("SearchReport", "derived"),
    "trials_per_cell": ("PrimaryMatrix", "derived"),
    "arm_trials": ("PrimaryMatrix", "derived"),
    "arms": ("PrimaryMatrix", "derived"),
    "shadow_fraction": ("PrimaryMatrix", "derived"),
    "resampling_unit": ("SamplingTarget", "derived"),
    "fit_model": ("analysis", "derived"),
    "confidence_level": ("SamplingTarget", "derived"),
    "exponent_half_width_target": ("SamplingTarget", "derived"),
    "minimum_causal_contrast": ("SamplingTarget", "derived"),
    "contrast_covariance_floor": ("SamplingTarget", "derived"),
    "measured_arm_latent_correlation": ("RecoveryReport", "derived"),
    "recovery_digest": ("RecoveryReport", "derived"),
    "curvature_test": ("analysis", "derived"),
    "lack_of_fit_test": ("analysis", "derived"),
    "inclusion_rule": ("analysis", "derived"),
    "sampling_target_digest": ("SamplingTarget", "derived"),
    "power_digest": ("PowerEstimate", "derived"),
    "power_stage": ("PowerEstimate", "derived"),
    "power_verdict": ("PowerEstimate", "derived"),
    "stream_namespace": ("PilotFirewall", "derived"),
    "pilot_stream_namespace": ("PilotFirewall", "derived"),
    "run_directory_prefix": ("PilotFirewall", "derived"),
    "pilot_run_directory_prefix": ("PilotFirewall", "derived"),
    "pilot_quarantine_directory": ("experiments module", "derived"),
    "pilot_trials_eligible": ("PilotFirewall", "derived"),
    "firewall_digest": ("PilotFirewall", "derived"),
    "pilot_plan_digest": ("PilotPlan", "derived"),
    "pilot_candidate_couplings": ("PilotPlan", "derived"),
    "pilot_trials_per_cell": ("PilotPlan", "derived"),
    "reference_rate_rule": ("experiments module", "derived"),
    "fixed_contraction_rate": ("RangeSelection via reference_rate", "derived"),
    "selection_digest": ("RangeSelection", "derived"),
    "post_range_digest": ("PostRangePlan", "derived"),
    "post_range_status": ("PostRangePlan", "derived"),
    "pilot_ledger_fingerprints": ("RangeSelection", "derived"),
    "benchmark_digest": ("BenchmarkResult", "derived"),
    "matrix_digest": ("PrimaryMatrix", "derived"),
    "resource_digest": ("ResourceEstimate", "derived"),
    "design_digest": ("recomputed over DESIGN_DIGEST_KEYS", "derived"),
    "stop_rule": ("PrimaryMatrix", "derived"),
    "failure_rule": ("PrimaryMatrix", "derived"),
    "non_claims": ("experiments module", "derived"),
}

#: The seal is a module-private attribute and **nothing depends on its
#: secrecy**.  A caller who reaches ``experiments._SEAL`` can build a
#: self-consistent proposal; what they cannot do is make it survive
#: :func:`require_authoritative_proposal`, which rebuilds every derived field
#: from the component records.  The same holds downstream: caller-built counts
#: are display-only through :func:`preview_range`, and the authoritative
#: selection and post-range power doors re-open and re-count the quarantined
#: runs themselves.  Secrecy would be a weak property; a rebuild is a checkable
#: one, and the suite tests exactly that.
SEAL_IS_NOT_THE_AUTHORITY = (
    "the module-private seal permits construction, never authority; every "
    "derived field is rebuilt from its component record and compared, "
    "caller-built counts are display-only, and the authoritative pilot doors "
    "reopen and recount")

#: The confidence level every interval in this module is stated at, and the
#: normal quantile that goes with it.  The level is ``analysis.CONFIDENCE_LEVEL``
#: rather than a second opinion about it: two layers disagreeing about what 95
#: per cent means is a defect nobody would find.
CONFIDENCE_LEVEL = _analysis.CONFIDENCE_LEVEL
CONFIDENCE_Z = 1.959963984540054

#: How many bootstrap standard errors a measured numerical shift is inflated by
#: before it is compared with an allowance.  Two, matching the coverage the
#: ticket-04 validation contract already uses: a shift that is small only
#: because nobody measured it precisely is not small.
COVERAGE_SIGMA = 2.0

#: The fraction of the planned production 95-per-cent confidence half-width
#: that numerical error is allowed to occupy.  One quarter, which is the closed
#: plan's number and is not a parameter of this module's decisions — it is
#: carried in the frozen budget so that changing it changes a digest.
ALLOWED_NUMERICAL_FRACTION = 0.25

#: The frozen count thresholds a range-only pilot may use, and nothing else.
#: Declared here, before any pilot exists, because that is the whole content of
#: "predeclared".
MIN_PILOT_EVENTS = 20
MIN_PILOT_SURVIVORS = 20
MIN_SELECTED_CELLS = 4

#: The minimum number of clocks that must be inside the tongue at the weakest
#: coupling of a proposed matrix.  Below this the eligible count is a visible
#: staircase rather than a population, and the claim narrows to that staircase
#: instead of quietly invoking a continuum.
MIN_ELIGIBLE_CLOCKS = 8

#: Pilot and production stream namespaces and run directories are prefixed, and
#: the prefixes are disjoint by construction rather than by convention.  A
#: :class:`PilotFirewall` refuses to exist if they overlap.
#: Taken from ``raw_runner`` rather than declared here.  The raw side owns the
#: reserved identity because it must recognize a pilot without importing this
#: module; a second constant would be a second opinion about the same fact.
PILOT_NAMESPACE_PREFIX = _raw_runner.EXPERIMENT_PILOT_NAMESPACE_PREFIX
PRODUCTION_NAMESPACE_PREFIX = "prod-"

#: Where a pilot run is **written**, so that the ordinary production reader
#: cannot address it.  A pilot run is created directly beneath this directory
#: through a pinned descriptor; it is not written under the production results
#: root and moved, and no transient production-root entry exists at any point.
#:
#: ``raw_runner.open_raw_run`` takes **one path segment** and opens
#: ``<package>/results/<segment>`` through pinned directory descriptors.  A run
#: parked at ``<package>/results/<QUARANTINE_DIRNAME>/<name>`` is therefore two
#: segments deep and is unreachable through that door: the two-segment spelling
#: is refused as a name, and the quarantine directory itself holds no run files.
#: An independent review found that a typed narrow return does not prevent a
#: second supported import path from opening the same five files, and it was
#: right — the fix is to move the files out of that path's reach rather than to
#: add a rule to it.
QUARANTINE_DIRNAME = "pilot_quarantine"

#: How many times the race walks every clock.  Ticket 05 walks each clock once
#: to find the channel stop and once again to censor the per-clock diagnostics
#: at it, because a single pass would have to retain the trial-by-clock-by-time
#: history the plan forbids.  The cost model charges both, and charges the
#: shadow rewalk separately on top.
RACE_PASSES = 2

#: The most adverse full-minus-control correlation the contrast is sized
#: against.  An earlier revision assumed the two arms' estimates were
#: independent and described that as conservative because "pairing can only
#: reduce" the contrast variance.  An independent review pointed out that this
#: is false: two arms sharing one Brownian tree can be *negatively* correlated
#: in their fitted exponents, and then pairing inflates the contrast rather
#: than shrinking it.  ``Var(delta) = 2 Var(p) (1 - rho)``, so a floor of
#: ``-0.25`` sizes the contrast at ``2.5 Var(p)`` instead of ``2 Var(p)``.
CONTRAST_COVARIANCE_FLOOR = -0.25

#: The rule that fixes the width-only control's one contraction rate, frozen
#: *before* the pilot and computed mechanically *after* it.  It is a sentence
#: rather than a number because the number does not exist until a range does.
REFERENCE_RATE_RULE = (
    "fixed_contraction_rate = dynamics.fixed_contraction_rate(low, high), the "
    "full-Adler central-clock relaxation rate at the geometric midpoint of the "
    "peak-coupling range the pilot selects under the frozen count-only rule; "
    "common to every clock, coupling, pulse time and trial, and never refitted")

#: The package sources whose bytes a proposal is bound to.  Every module the
#: run or its analysis passes through, so that "the same source produced this"
#: is checkable rather than asserted.  ``verify.py`` is deliberately absent: it
#: is the checker, not the experiment, and pinning it would make every new
#: check invalidate every proposal.
SOURCE_FILES = (
    "__init__.py", "analysis.py", "analytic.py", "commitment.py", "compare.py",
    "dynamics.py", "experiments.py", "killed_diffusion.py", "model.py",
    "moving_band_audit.py", "observables.py", "raw_config.py",
    "raw_experiments.py", "raw_ledger.py", "raw_race.py", "raw_runner.py",
    "simulate.py", "stochastic.py", "validation.py",
)

_PACKAGE_DIR = Path(__file__).resolve().parent

#: How many bytes one float64 occupies.  Named because the noise-cube argument
#: is entirely about this number times three loop bounds.
_FLOAT_BYTES = 8


# ---------------------------------------------------------------------------
# Canonical bytes
# ---------------------------------------------------------------------------
def canonical_bytes(payload) -> bytes:
    """Return the one byte string that stands for ``payload``.

    Sorted keys, tight separators, ASCII only, and **no non-finite numbers**.
    Two records that say the same thing hash the same; a record that grew a key
    hashes differently; and a NaN cannot enter a digest at all, because a NaN
    that reaches a frozen record is a number nobody can compare and the digest
    would make it look decided.

    Deliberately *not* the human-readable spelling
    ``raw_ledger.canonical_manifest_bytes`` uses.  That file is opened by
    people; these digests are compared by machines, and the two jobs want
    opposite formatting.  The proposal's own human-readable rendering is
    :meth:`ProposedProductionManifest.as_dict`, which a caller may print with
    whatever indentation it likes without changing the digest.
    """
    if not isinstance(payload, dict):
        raise TypeError(f"payload must be a dict, got "
                        f"{type(payload).__name__}")
    for key in payload:
        if not isinstance(key, str):
            raise TypeError(f"every key must be a string, got "
                            f"{type(key).__name__} {key!r}")
    try:
        text = json.dumps(payload, sort_keys=True, ensure_ascii=True,
                          allow_nan=False, separators=(",", ":"),
                          default=_reject_value)
    except ValueError as exc:                       # allow_nan=False raises it
        raise ValueError(f"a non-finite number cannot enter a canonical "
                         f"record: {exc}") from exc
    return text.encode("utf-8")


def _reject_value(value):
    """``json.dumps`` default hook: refuse rather than repr an unknown type."""
    raise TypeError(f"a canonical record holds only strings, finite numbers, "
                    f"booleans, null and lists of those; got "
                    f"{type(value).__name__} {value!r}")


def _digest(payload) -> str:
    """The hex digest of :func:`canonical_bytes`, under the schema tag."""
    tagged = dict(payload)
    tagged["__schema__"] = EXPERIMENTS_SCHEMA
    return hashlib.new(_raw_ledger.DIGEST_ALGORITHM,
                       canonical_bytes(tagged)).hexdigest()


def _require_fraction(value, name: str) -> float:
    """A number in ``[0, 1]``."""
    out = require_finite(value, name)
    if not 0.0 <= out <= 1.0:
        raise ValueError(f"{name} must lie in [0, 1], got {out!r}")
    return out


def _require_positive(value, name: str) -> float:
    out = require_finite(value, name)
    if not out > 0.0:
        raise ValueError(f"{name} must be strictly positive, got {out!r}")
    return out


def _require_hex_digest(value, name: str) -> str:
    """A real lowercase 64-character hex digest, or a refusal.

    An independent review found every provenance type accepting arbitrary
    64-character strings — ``"z" * 64`` passed as a ledger digest — so a
    fabricated record was indistinguishable from a measured one by shape alone.
    Length was never the property that mattered.
    """
    text = _require_text(value, name, maximum=_raw_ledger.DIGEST_HEX_LENGTH)
    if len(text) != _raw_ledger.DIGEST_HEX_LENGTH:
        raise ValueError(f"{name} must be {_raw_ledger.DIGEST_HEX_LENGTH} hex "
                         f"characters, got {len(text)}")
    if any(character not in "0123456789abcdef" for character in text):
        raise ValueError(
            f"{name}={text!r} is not lowercase hexadecimal; a digest is "
            f"the "
            f"output of a hash function, and a string of the right length is "
            f"not one")
    return text


def _require_text(value, name: str, *, maximum: int = 400) -> str:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, str):
        raise TypeError(f"{name} must be a string, got "
                        f"{type(value).__name__} {value!r}")
    if not value:
        raise ValueError(f"{name} must not be empty")
    if len(value) > maximum:
        raise ValueError(f"{name} must be at most {maximum} characters, got "
                         f"{len(value)}")
    if not value.isprintable():
        raise ValueError(f"{name} contains a non-printable character")
    return value


def _require_member(value, allowed, name: str) -> str:
    value = _require_text(value, name)
    if value not in allowed:
        raise ValueError(f"{name} must be one of {list(allowed)}, got "
                         f"{value!r}")
    return value


def _require_couplings(values, name: str) -> tuple:
    """A strictly increasing tuple of at least ``MIN_CELLS`` positive nodes."""
    out = []
    for position, value in enumerate(tuple(values)):
        out.append(_require_positive(value, f"{name}[{position}]"))
    if len(out) < _analysis.MIN_CELLS:
        raise ValueError(
            f"{name} holds {len(out)} nodes; the frozen cloglog model needs at "
            f"least {_analysis.MIN_CELLS} coupling cells before an exponent "
            f"and its curvature extension are identified")
    if any(later <= earlier for earlier, later in zip(out, out[1:])):
        raise ValueError(f"{name} must be strictly increasing, got {out!r}; a "
                         f"repeated node is one cell wearing two labels")
    return tuple(out)


# ---------------------------------------------------------------------------
# Provenance
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class SourceFingerprint:
    """The bytes and the interpreter a measurement or a proposal came from.

    ``files``
        ``((name, sha256), ...)`` over :data:`SOURCE_FILES`, in that order.
    ``python_version``, ``numpy_version``, ``platform``
        The recorded environment.  A run is reproducible *in a recorded
        environment*; a claim of reproducibility that does not say which one is
        not checkable.

    Built by :func:`source_fingerprint`.  The constructor is public because a
    reviewer reconstructing a proposal from an archive has to be able to state
    the fingerprint the archive carried rather than the one this machine has.
    """

    files: tuple
    python_version: str
    numpy_version: str
    platform: str

    def __post_init__(self):
        entries = []
        for position, entry in enumerate(tuple(self.files)):
            pair = tuple(entry)
            if len(pair) != 2:
                raise ValueError(f"files[{position}] must be a (name, digest) "
                                 f"pair, got {pair!r}")
            name = _require_text(pair[0], f"files[{position}] name",
                                 maximum=120)
            entries.append((name, _require_hex_digest(
                pair[1], f"files[{position}] digest")))
        if not entries:
            raise ValueError("a source fingerprint over no files fingerprints "
                             "nothing")
        if len({name for name, _ in entries}) != len(entries):
            raise ValueError("a file appears twice in the fingerprint")
        object.__setattr__(self, "files", tuple(entries))
        for field in ("python_version", "numpy_version", "platform"):
            object.__setattr__(self, field,
                               _require_text(getattr(self, field), field))

    @property
    def digest(self) -> str:
        """The one hash a proposal cites for "these bytes, this interpreter"."""
        return _digest({"files": [list(pair) for pair in self.files],
                        "python_version": self.python_version,
                        "numpy_version": self.numpy_version,
                        "platform": self.platform})


def source_fingerprint() -> SourceFingerprint:
    """Fingerprint the package as it stands on disk right now.

    Every file in :data:`SOURCE_FILES` must exist.  A missing one is a refusal
    rather than a shorter fingerprint, because a fingerprint that silently
    covers less is worse than none: it still looks like provenance.
    """
    import platform as _platform                    # noqa: PLC0415 (local)

    entries = []
    for name in SOURCE_FILES:
        path = _PACKAGE_DIR / name
        if not path.is_file():
            raise ValueError(f"{name} is named in SOURCE_FILES but is not a "
                             f"regular file under {_PACKAGE_DIR}; a "
                             f"fingerprint that quietly covers fewer files "
                             f"than it claims is not provenance")
        entries.append((name, hashlib.new(_raw_ledger.DIGEST_ALGORITHM,
                                          path.read_bytes()).hexdigest()))
    return SourceFingerprint(tuple(entries), _platform.python_version(),
                             np.__version__, _platform.platform())


# ---------------------------------------------------------------------------
# Feasibility: what the real streamed race actually costs
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class BenchmarkFixture:
    """One bounded, reproducible benchmark run of the **real** race.

    ``config``
        Exactly a :class:`~.raw_config.RawEventConfig`.  Not a surrogate, not a
        reduced kernel, not a vectorized reimplementation: the thing being
        measured is ``raw_runner.write_raw_run``, because that is the thing
        production would run, and ticket 05 chose the scalar path deliberately.
    ``clock_block``, ``step_window``
        The live-block sizes the run executes with.  They cannot change a
        recorded number — that is ticket 05's contract and this file relies on
        it — so they belong to the measurement rather than to the run.
    ``warmup``, ``repeats``
        How many runs are discarded before timing starts, and how many are
        timed.  Both recorded, and both are now **required to be at least
        one**: an independent review found the scaling fixtures running with no
        warmup and a single timed repeat, which cannot establish machine
        variation or support a conservative slow-repeat rate.

    ``label`` doubles as the run-directory stem, so it is restricted to the
    characters a run name may contain.
    """

    label: str
    config: object
    model: str
    clock_block: int
    step_window: int
    warmup: int
    repeats: int

    def __post_init__(self):
        label = _require_text(self.label, "label", maximum=48)
        if not all(character.isalnum() or character in "-_"
                   for character in label):
            raise ValueError(f"label must be letters, digits, hyphen or "
                             f"underscore, got {label!r}; it becomes a run "
                             f"directory name")
        config = require_exact_type(self.config, _raw_config.RawEventConfig,
                                    "config")
        model = _require_member(self.model, _raw_ledger.MODEL_LABELS, "model")
        object.__setattr__(self, "label", label)
        object.__setattr__(self, "config", config)
        object.__setattr__(self, "model", model)
        for field in ("clock_block", "step_window"):
            object.__setattr__(self, field,
                               require_positive_int(getattr(self, field),
                                                    field, minimum=1))
        object.__setattr__(self, "warmup",
                           require_positive_int(self.warmup, "warmup",
                                                minimum=1))
        object.__setattr__(self, "repeats",
                           require_positive_int(self.repeats, "repeats",
                                                minimum=2))

    @property
    def steps(self) -> int:
        """Finest steps in this fixture's run window."""
        _, _, steps = _raw_runner.raw_run_window(self.config)
        return int(steps)

    @property
    def nominal_clock_steps(self) -> int:
        """``trials * clocks * steps``.

        The *shape* of the workload, and deliberately **not** the work.  It is
        kept because the trial-by-clock-by-time cube this design refuses to
        materialize is exactly this many float64 values, so it is the number
        the no-cube argument is about.
        """
        return (int(self.config.trials) * int(self.config.grid.n_points)
                * self.steps)

    @property
    def shadow_clock_steps(self) -> int:
        """What the no-stop shadow rewalk adds, on its own.

        Charged separately because an independent review found the cost model
        multiplying a rate whose benchmark *already contained* a shadow
        fraction by another ``1.05``.  A shadow sample is
        ``shadow_trials * clocks * steps`` of extra walking; the fixture's own
        shadow fraction and the matrix's need not be the same number, so the
        two are normalized here rather than assumed equal.
        """
        return (int(self.config.shadow_trials) * int(self.config.grid.n_points)
                * self.steps)

    @property
    def walked_clock_steps(self) -> int:
        """The work model: ``passes * nominal + shadow``.

        ``RACE_PASSES`` is two, because ticket 05 walks every clock once to
        find the channel stop and once more to censor at it.  This is an
        **upper bound**: the first pass truncates each clock at the best time
        found so far, so a run with commitments does less than this.

        That direction is deliberate.  A fixture in which nothing commits does
        exactly this much work, so its measured rate is the true rate; a
        fixture with commitments appears *faster* under this model.  Taking the
        slowest such rate across fixtures therefore prices a matrix
        conservatively: any real matrix that truncates costs less than the
        estimate, never more.
        """
        return RACE_PASSES * self.nominal_clock_steps + self.shadow_clock_steps

    @property
    def digest(self) -> str:
        """Everything that decides what this fixture measures."""
        grid = self.config.grid
        return _digest({
            "label": self.label, "model": self.model,
            "clock_block": self.clock_block, "step_window": self.step_window,
            "warmup": self.warmup, "repeats": self.repeats,
            "passes": RACE_PASSES,
            "half_width": grid.half_width, "n_points": grid.n_points,
            "peak_coupling": self.config.peak_coupling,
            "pulse_duration": self.config.pulse_duration,
            "pulse_centre": self.config.pulse_centre,
            "phase_diffusion": self.config.phase_diffusion,
            "lock_tolerance": self.config.lock_tolerance,
            "dwell_time": self.config.dwell_time,
            "timestep": self.config.timestep,
            "trials": self.config.trials,
            "shadow_trials": self.config.shadow_trials,
            "stream_namespace": self.config.stream_namespace,
            "run_label": self.config.run_label,
        })


#: What ``run_benchmark`` did, in one sentence, recorded inside every result so
#: that a reader does not have to reconstruct the protocol from field names.
BENCHMARK_PROTOCOL = (
    "raw_runner.write_raw_run at the fixture's own clock_block and "
    "step_window; warmup runs written and removed before timing; every timed "
    "repeat run under tracemalloc and its wall time recorded individually; "
    "throughput quoted per walked clock-step (RACE_PASSES * trials * clocks * "
    "steps + shadow_trials * clocks * steps) at the slowest timed repeat; "
    "every run directory removed on every exit path")


@dataclass(frozen=True)
class BenchmarkResult:
    """What one fixture measured, and what it proves about memory.

    ``repeat_seconds``
        **Every** timed repeat, in order.  An independent review asked for the
        raw repeats rather than a summary, and it is right: a single number
        cannot say whether the machine was quiet, and a mean cannot say whether
        one repeat was three times another.
    ``wall_seconds``
        The slowest of them.  A cost model built on the fastest repeat is a
        cost model that has already assumed the machine will stay quiet.
    ``peak_bytes``
        ``tracemalloc``'s peak over the timed runs.  NumPy's own buffers are
        inside that number — the allocator is traced under
        ``numpy.lib.tracemalloc_domain`` — so this is a measurement of the
        noise blocks and the row objects together, not of Python overhead
        beside an untracked array.
    ``storage_bytes``, ``rows``, ``fixed_bytes``
        The five written files, the ledger + per-clock + shadow rows they hold,
        and the part of the file bytes that is **not** proportional to rows —
        the manifest, the close marker and three header lines.  Separated
        because a per-row figure that silently contains a five-kilobyte
        manifest overstates storage on a big run and understates it on a small
        one.
    ``traced``
        Always true, and recorded anyway.  Every timed repeat runs with
        ``tracemalloc`` on, which slows it down; the throughput this record
        reports is therefore conservative by construction, and a reader is
        entitled to know that rather than infer it.
    """

    fixture_label: str
    fixture_digest: str
    source_digest: str
    protocol: str
    nominal_clock_steps: int
    shadow_clock_steps: int
    walked_clock_steps: int
    passes: int
    trials: int
    shadow_trials: int
    clocks: int
    steps: int
    repeat_seconds: tuple
    peak_bytes: int
    storage_bytes: int
    fixed_bytes: int
    rows: int
    committed_trials: int
    traced: bool = True

    def __post_init__(self):
        object.__setattr__(self, "fixture_label",
                           _require_text(self.fixture_label, "fixture_label",
                                         maximum=48))
        for field in ("fixture_digest", "source_digest"):
            object.__setattr__(self, field,
                               _require_hex_digest(getattr(self, field),
                                                   field))
        object.__setattr__(self, "protocol",
                           _require_text(self.protocol, "protocol",
                                         maximum=800))
        for field in ("nominal_clock_steps", "shadow_clock_steps",
                      "walked_clock_steps", "passes", "trials",
                      "shadow_trials", "clocks", "steps", "peak_bytes",
                      "storage_bytes", "fixed_bytes", "rows"):
            object.__setattr__(self, field,
                               require_positive_int(getattr(self, field),
                                                    field, minimum=1))
        object.__setattr__(self, "committed_trials",
                           require_positive_int(self.committed_trials,
                                                "committed_trials", minimum=0))
        if self.committed_trials > self.trials:
            raise ValueError(f"committed_trials={self.committed_trials} "
                             f"exceeds trials={self.trials}")
        if self.shadow_trials > self.trials:
            raise ValueError(f"shadow_trials={self.shadow_trials} exceeds "
                             f"trials={self.trials}")
        if self.fixed_bytes >= self.storage_bytes:
            raise ValueError(
                f"fixed_bytes={self.fixed_bytes} is not below "
                f"storage_bytes={self.storage_bytes}; the fixed part is the "
                f"manifest, the marker and the three header lines, and a run "
                f"whose rows contributed nothing wrote no rows")
        expected = (self.passes * self.nominal_clock_steps
                    + self.shadow_clock_steps)
        if self.walked_clock_steps != expected:
            raise ValueError(
                f"walked_clock_steps={self.walked_clock_steps} is not "
                f"{self.passes} * {self.nominal_clock_steps} + "
                f"{self.shadow_clock_steps} = {expected}; the work model is "
                f"the whole content of a throughput figure and it is not a "
                f"free field")
        repeats = tuple(_require_positive(value, f"repeat_seconds[{position}]")
                        for position, value in
                        enumerate(tuple(self.repeat_seconds)))
        if len(repeats) < 2:
            raise ValueError(
                f"repeat_seconds holds {len(repeats)} timed repeat(s); a "
                f"single repeat cannot establish machine variation, and an "
                f"independent review found exactly that being reported as a "
                f"throughput measurement")
        object.__setattr__(self, "repeat_seconds", repeats)
        if not isinstance(self.traced, (bool, np.bool_)):
            raise TypeError(f"traced must be a bool, got "
                            f"{type(self.traced).__name__}")
        object.__setattr__(self, "traced", bool(self.traced))

    @property
    def wall_seconds(self) -> float:
        """The slowest timed repeat."""
        return max(self.repeat_seconds)

    @property
    def fastest_seconds(self) -> float:
        """The fastest timed repeat."""
        return min(self.repeat_seconds)

    @property
    def repeat_spread(self) -> float:
        """Slowest over fastest.  How quiet the machine was, as one number."""
        return self.wall_seconds / self.fastest_seconds

    @property
    def clock_steps_per_second(self) -> float:
        """**Walked** clock-steps per second at the slowest timed repeat.

        This is the planning rate, and the unit matters: a matrix priced with
        it must compute its own walked clock-steps the same way — both passes
        and its own shadow fraction — rather than multiplying a
        shadow-bearing rate by a shadow factor again.
        """
        return self.walked_clock_steps / self.wall_seconds

    @property
    def nominal_clock_steps_per_second(self) -> float:
        """Nominal clock-steps per second.  Reported, never used for costing."""
        return self.nominal_clock_steps / self.wall_seconds

    @property
    def truncating(self) -> bool:
        """Did any trial commit, so that the first pass truncated?

        A fixture with no commitments does exactly ``walked_clock_steps`` of
        work, so its rate is the honest one.  A fixture with commitments does
        less and therefore looks faster.  Costing takes the slowest rate across
        fixtures, and this flag is what lets a reader see which kind each one
        was.
        """
        return self.committed_trials > 0

    @property
    def cube_bytes(self) -> int:
        """What a trial-by-clock-by-time float64 noise cube would have cost."""
        return self.nominal_clock_steps * _FLOAT_BYTES

    @property
    def cube_margin(self) -> float:
        """How many times smaller the peak was than that cube.

        Not a proof on its own — a single ratio can be made large by choosing a
        big fixture — which is why the suite pairs it with the *flatness* test:
        quadrupling trials, or quadrupling steps, quadruples this cube and must
        leave the peak nearly where it was.
        """
        return self.cube_bytes / self.peak_bytes

    @property
    def bytes_per_row(self) -> float:
        """Row-proportional bytes per durable row, with the fixed part removed."""
        return (self.storage_bytes - self.fixed_bytes) / self.rows

    @property
    def digest(self) -> str:
        return _digest({
            "fixture_digest": self.fixture_digest,
            "source_digest": self.source_digest, "protocol": self.protocol,
            "nominal_clock_steps": self.nominal_clock_steps,
            "shadow_clock_steps": self.shadow_clock_steps,
            "walked_clock_steps": self.walked_clock_steps,
            "passes": self.passes, "trials": self.trials,
            "shadow_trials": self.shadow_trials, "clocks": self.clocks,
            "steps": self.steps,
            "repeat_seconds": list(self.repeat_seconds),
            "peak_bytes": self.peak_bytes,
            "storage_bytes": self.storage_bytes,
            "fixed_bytes": self.fixed_bytes, "rows": self.rows,
            "committed_trials": self.committed_trials, "traced": self.traced,
        })


def _remove_run(run_name: str) -> None:
    """Remove a benchmark run directory.  Idempotent, and never a surprise."""
    try:
        target = _raw_runner.raw_results_directory(run_name)
    except (TypeError, ValueError):
        return
    shutil.rmtree(target, ignore_errors=True)


def run_benchmark(fixture, fingerprint) -> BenchmarkResult:
    """Measure ``fixture`` against the real streamed race and leave nothing.

    The measurement, in order:

    1. ``warmup`` runs, written and removed, timed and discarded.  The first
       call into the race pays for module import, the pulse schedule and the
       first NumPy allocations, and charging production for that would be wrong
       in the optimistic direction.  At least one, always.
    2. ``repeats`` timed runs — at least two — each under ``tracemalloc``, each
       written to its own run directory and removed before the next begins.
       Every wall time is retained, not just the extremes.
    3. The record: all repeats, the largest traced peak, the file sizes split
       into their fixed and row-proportional parts, and the committed-trial
       count that says whether the first pass truncated.

    **Every run directory is removed on every exit path**, including an
    exception, so a benchmark leaves the package results root exactly as it
    found it.  That is not tidiness: ``verify.py`` asserts the results root is
    empty when the suite ends, and a benchmark that littered would fail the
    suite somewhere else entirely.
    """
    require_exact_type(fixture, BenchmarkFixture, "fixture")
    require_exact_type(fingerprint, SourceFingerprint, "fingerprint")

    config = fixture.config
    names = []
    timings = []
    peak = 0
    storage = 0
    fixed = 0
    rows = 0
    committed = 0
    try:
        for index in range(fixture.warmup + fixture.repeats):
            timed = index >= fixture.warmup
            name = f"{fixture.label}-{index:03d}"
            _remove_run(name)
            names.append(name)
            if timed:
                tracemalloc.start()
            start = time.perf_counter()
            report = _raw_runner.write_raw_run(config, name, fixture.model,
                                               fixture.clock_block,
                                               fixture.step_window)
            elapsed = time.perf_counter() - start
            if timed:
                peak = max(peak, tracemalloc.get_traced_memory()[1])
                tracemalloc.stop()
                timings.append(elapsed)
                directory = Path(report.directory)
                sizes = {entry.name: entry.stat().st_size
                         for entry in directory.iterdir() if entry.is_file()}
                storage = sum(sizes.values())
                # The fixed part: the manifest, the close marker, and the three
                # CSV header lines.  Everything else scales with rows.
                fixed = (sizes.get(_raw_runner.RUN_FILENAMES["manifest"], 0)
                         + sizes.get(_raw_runner.RUN_FILENAMES["marker"], 0))
                for table in ("ledger", "clocks", "shadow"):
                    columns = {"ledger": _raw_ledger.LEDGER_COLUMNS,
                               "clocks": _raw_ledger.CLOCK_COLUMNS,
                               "shadow": _raw_ledger.SHADOW_COLUMNS}[table]
                    fixed += len(_raw_ledger.header_row(columns)) + 1
                rows = sum(count for _, count, _ in report.marker.tables)
                committed = dict(report.categories)["committed"]
            _remove_run(name)
    finally:
        if tracemalloc.is_tracing():
            tracemalloc.stop()
        for name in names:
            _remove_run(name)

    return BenchmarkResult(
        fixture_label=fixture.label, fixture_digest=fixture.digest,
        source_digest=fingerprint.digest, protocol=BENCHMARK_PROTOCOL,
        nominal_clock_steps=fixture.nominal_clock_steps,
        shadow_clock_steps=fixture.shadow_clock_steps,
        walked_clock_steps=fixture.walked_clock_steps, passes=RACE_PASSES,
        trials=int(config.trials), shadow_trials=int(config.shadow_trials),
        clocks=int(config.grid.n_points), steps=fixture.steps,
        repeat_seconds=tuple(timings), peak_bytes=int(peak),
        storage_bytes=int(storage), fixed_bytes=int(fixed), rows=int(rows),
        committed_trials=int(committed), traced=True)


# ---------------------------------------------------------------------------
# Power and precision for the binomial cloglog exponent
# ---------------------------------------------------------------------------

def _cloglog_weight(cumulative_hazard):
    """Per-trial Fisher information for the linear predictor.

    With ``eta = log H``, ``q = 1 - exp(-H)`` and one binomial trial, the
    information carried about ``eta`` is

    ``(dq/deta)^2 / (q (1 - q)) = H^2 exp(-H) / (1 - exp(-H)) = H^2/(e^H - 1)``.

    It vanishes at both ends and peaks near ``H = 1.6``, which is the whole
    reason a coupling sweep has to *straddle* the crossover: a cell with no
    events and a cell with nothing but events both contribute almost nothing to
    the slope, however many trials are spent in them.  They are still kept in
    the likelihood — ``analysis.fit_exponent`` never drops a cell — but they do
    not buy precision, and a power calculation that ignored that would size the
    experiment from cells that carry no information.
    """
    hazard = np.asarray(cumulative_hazard, dtype=float)
    return np.where(hazard > 0.0, hazard * hazard / np.expm1(hazard), 0.0)


def _weight_at_probability(probability):
    """The same weight, written in terms of a commitment probability."""
    probability = np.asarray(probability, dtype=float)
    hazard = -np.log1p(-probability)
    return _cloglog_weight(hazard)


def cloglog_slope_variance(couplings, trials, amplitude, exponent) -> float:
    """Return ``Var(exponent-hat)`` for one trial allocation, independent cells.

    The design-based variance of the frozen model
    ``log(-log(1 - q_K)) = amplitude + exponent * log(K)`` under independent
    binomial sampling: ``S0 / (S0 S2 - S1^2)`` with
    ``Sk = sum_i n_i w_i (log K_i)^k`` and ``w_i`` the weight above.

    **This is the analytic independent-cell number and nothing else.**  Three
    quantities are deliberately kept apart everywhere in this module, because
    an independent review found them running together:

    * *this* — exact, checkable, and about a design nobody is proposing;
    * the **measured** common-random-number inflation, which
      :func:`recover_exponent` estimates and which has its own standard error;
      and
    * the **declared conservative factor**, which is a judgement written down
      before the measurement and required to exceed it.

    Reporting one number would hide which of the three moved.
    """
    nodes = _require_couplings(couplings, "couplings")
    counts = []
    for position, value in enumerate(tuple(trials)):
        counts.append(require_positive_int(value, f"trials[{position}]",
                                           minimum=1))
    if len(counts) != len(nodes):
        raise ValueError(f"{len(nodes)} couplings and {len(counts)} trial "
                         f"counts; every cell carries exactly one")
    amplitude = require_finite(amplitude, "amplitude")
    exponent = require_finite(exponent, "exponent")

    logk = np.log(np.asarray(nodes, dtype=float))
    hazard = np.exp(amplitude + exponent * logk)
    if not np.all(np.isfinite(hazard)):
        raise ValueError(
            f"amplitude={amplitude!r} and exponent={exponent!r} put the "
            f"cumulative hazard outside the representable range at one of "
            f"{nodes!r}; a power calculation on an overflowed design is a "
            f"number about nothing")
    weight = _cloglog_weight(hazard) * np.asarray(counts, dtype=float)
    return _slope_variance_from_weights(weight, logk)


def _slope_variance_from_weights(weight, logk):
    """``S0 / (S0 S2 - S1^2)``, refusing a singular information matrix."""
    weight = np.asarray(weight, dtype=float)
    logk = np.asarray(logk, dtype=float)
    s0 = float(weight.sum())
    s1 = float((weight * logk).sum())
    s2 = float((weight * logk * logk).sum())
    determinant = s0 * s2 - s1 * s1
    if not determinant > 0.0:
        raise ValueError(
            f"the expected information matrix is singular (determinant "
            f"{determinant!r}); with these couplings and this hazard every "
            f"informative cell has collapsed to one point and no slope is "
            f"identified")
    return s0 / determinant


@dataclass(frozen=True)
class ArmEnvelope:
    """The predeclared conservative information rule for **one** model arm.

    An independent review's finding, restated: sizing the whole matrix from one
    assumed response curve prices the full-Adler arm and silently assumes the
    two mandatory negative controls are equally informative.  They need not be.
    A one-clock central control commits far less often than a sixty-four-clock
    population; a width-only control with one fixed contraction rate has a
    different curve again.  An arm that is less informative needs *more*
    trials, and a design sized from the most informative arm is under-powered
    for the contrast it exists to support.

    So every arm declares an envelope, **before** any pilot:

    ``probability_low``, ``probability_high``
        The interval the arm's per-cell commitment probability is assumed to
        lie in.  The sizing uses the **least informative** probability in that
        interval — the one minimizing ``H^2/(e^H - 1)`` — so a wider envelope
        is a more conservative one and cannot be narrowed after seeing counts.
    ``information_factor``
        A multiplier on that weight, at most one, for anything the envelope
        does not capture.  One means "the envelope is the whole story".
    ``invalid_fit_probability``
        The declared probability that a replicate of this arm returns
        ``no_valid_exponent`` rather than an estimate.  Trials are inflated by
        ``1/(1 - p)`` to pay for it, because a design that budgets only for
        successful fits has not budgeted for the experiment.
    ``clocks``
        How many clocks the arm runs, which is what makes its *cost* different
        from its *precision*.
    ``pairs_with_full``
        Whether this arm shares master trials with the full-Adler arm.  The
        width-only control does, by construction — it is the same population on
        the same Brownian tree with one drift term changed — so the two must
        carry the **same** trial count or the pairing that makes their contrast
        cheap does not exist.  The central-clock control does not: it is a
        one-clock population, a different experiment family, and tying its
        count to the full arm's would either over-run it or under-power it for
        no reason.  Sizing each arm from its own envelope and then tying only
        the genuinely paired ones is what stops one conservative envelope
        inflating every arm at once.
    """

    name: str
    clocks: int
    probability_low: float
    probability_high: float
    information_factor: float
    invalid_fit_probability: float
    pairs_with_full: bool = True

    def __post_init__(self):
        object.__setattr__(self, "name",
                           _require_text(self.name, "name", maximum=40))
        object.__setattr__(self, "clocks",
                           require_positive_int(self.clocks, "clocks",
                                                minimum=1))
        low = _require_fraction(self.probability_low, "probability_low")
        high = _require_fraction(self.probability_high, "probability_high")
        if not 0.0 < low < high < 1.0:
            raise ValueError(
                f"the probability envelope must satisfy "
                f"0 < low < high < 1, got low={low!r}, high={high!r}; a "
                f"degenerate or reversed envelope is not a conservative rule, "
                f"and an endpoint at exactly zero or one is a cell that "
                f"carries no slope information at all")
        object.__setattr__(self, "probability_low", low)
        object.__setattr__(self, "probability_high", high)
        factor = _require_positive(self.information_factor,
                                   "information_factor")
        if factor > 1.0:
            raise ValueError(
                f"information_factor={factor!r} exceeds one, which claims the "
                f"arm carries more information than its own envelope implies; "
                f"this multiplier exists to be conservative and cannot be "
                f"used to make an arm cheaper")
        object.__setattr__(self, "information_factor", factor)
        invalid = _require_fraction(self.invalid_fit_probability,
                                    "invalid_fit_probability")
        if invalid >= 1.0:
            raise ValueError(
                f"invalid_fit_probability={invalid!r} says every replicate "
                f"fails, so no trial count is enough and the arm is not "
                f"sizable")
        object.__setattr__(self, "invalid_fit_probability", invalid)
        if not isinstance(self.pairs_with_full, (bool, np.bool_)):
            raise TypeError(f"pairs_with_full must be a bool, got "
                            f"{type(self.pairs_with_full).__name__}")
        object.__setattr__(self, "pairs_with_full", bool(self.pairs_with_full))

    @property
    def weight(self) -> float:
        """The least information one trial of this arm is assumed to carry."""
        endpoints = np.array([self.probability_low, self.probability_high])
        return float(_weight_at_probability(endpoints).min()
                     * self.information_factor)

    @property
    def failure_inflation(self) -> float:
        """``1 / (1 - invalid_fit_probability)`` — the price of refused fits."""
        return 1.0 / (1.0 - self.invalid_fit_probability)

    @property
    def digest(self) -> str:
        return _digest({"name": self.name, "clocks": self.clocks,
                        "probability_low": self.probability_low,
                        "probability_high": self.probability_high,
                        "information_factor": self.information_factor,
                        "invalid_fit_probability":
                            self.invalid_fit_probability,
                        "pairs_with_full": self.pairs_with_full})


@dataclass(frozen=True)
class SamplingTarget:
    """The frozen precision the production sweep is being sized for.

    Everything here is chosen **before** any production data exists and is
    hashed into the proposal, which is what makes "the trial count follows the
    frozen precision target" a checkable statement rather than a description of
    whatever count was affordable.

    ``stage``
        ``provisional_pre_pilot`` or ``post_range_join``.  A pre-pilot target
        assumes a coupling range the count-only pilot has not chosen, so the
        number it produces is **provisional** and is labelled as such
        everywhere it appears.  Only :func:`post_range_power` produces the
        second stage, and it takes nothing but a selected range and the count
        summaries the firewall permits.
    ``exponent_half_width``
        The 95-per-cent confidence half-width the fitted coupling exponent must
        achieve.  Capped by ``analysis.MAX_EXPONENT_HALF_WIDTH``, because a
        target looser than the analysis layer's own validity cap would size an
        experiment whose result that layer would refuse.
    ``minimum_causal_contrast``
        The smallest full-model-minus-control exponent difference the sweep
        must resolve.  The plan's causal defeat rule turns on it.
    ``contrast_covariance_floor``
        The most adverse full-minus-control correlation the contrast is sized
        against.  See :data:`CONTRAST_COVARIANCE_FLOOR`.
    ``pairing_inflation``
        The **declared** factor by which common-random-number pairing
        multiplies the exponent's sampling variance.  It is required to be at
        least the measured one, never derived from it.
    ``conservative_factor``
        A separate declared multiplier on the final trial count, for everything
        neither the envelope nor the inflation captures.
    ``arm_envelopes``
        One :class:`ArmEnvelope` per arm, and the arms are exactly its names.
    ``maximum_trials_per_cell``
        The affordability ceiling.  A design that needs more is
        ``no_powered_matrix`` — a result, not an exception.
    """

    exponent_half_width: float
    minimum_causal_contrast: float
    coupling_low: float
    coupling_high: float
    cells: int
    pairing_inflation: float
    conservative_factor: float
    resampling_unit: str
    arm_envelopes: tuple
    shadow_fraction: float
    maximum_trials_per_cell: int
    stage: str = "provisional_pre_pilot"
    contrast_covariance_floor: float = CONTRAST_COVARIANCE_FLOOR
    confidence_level: float = CONFIDENCE_LEVEL

    def __post_init__(self):
        half = _require_positive(self.exponent_half_width,
                                 "exponent_half_width")
        if half > _analysis.MAX_EXPONENT_HALF_WIDTH:
            raise ValueError(
                f"exponent_half_width={half!r} is looser than "
                f"analysis.MAX_EXPONENT_HALF_WIDTH="
                f"{_analysis.MAX_EXPONENT_HALF_WIDTH!r}; sizing an experiment "
                f"to a precision the analysis layer would refuse produces a "
                f"run whose own estimator reports no valid exponent")
        object.__setattr__(self, "exponent_half_width", half)
        object.__setattr__(self, "minimum_causal_contrast",
                           _require_positive(self.minimum_causal_contrast,
                                             "minimum_causal_contrast"))
        low = _require_positive(self.coupling_low, "coupling_low")
        high = _require_positive(self.coupling_high, "coupling_high")
        if not high > low:
            raise ValueError(f"coupling_high={high!r} must exceed "
                             f"coupling_low={low!r}; a range of one point has "
                             f"no lever arm and no exponent")
        object.__setattr__(self, "coupling_low", low)
        object.__setattr__(self, "coupling_high", high)
        object.__setattr__(self, "cells",
                           require_positive_int(self.cells, "cells",
                                                minimum=_analysis.MIN_CELLS))
        inflation = _require_positive(self.pairing_inflation,
                                      "pairing_inflation")
        if inflation < 1.0:
            raise ValueError(
                f"pairing_inflation={inflation!r} is below one, which claims "
                f"that reusing one Brownian path across couplings makes the "
                f"slope cheaper than independent sampling; measure it before "
                f"asserting it, and if a measurement really says so, record "
                f"the measurement and set this to one")
        object.__setattr__(self, "pairing_inflation", inflation)
        factor = _require_positive(self.conservative_factor,
                                   "conservative_factor")
        if factor < 1.0:
            raise ValueError(f"conservative_factor={factor!r} is below one, "
                             f"which is not conservative")
        object.__setattr__(self, "conservative_factor", factor)
        object.__setattr__(self, "resampling_unit",
                           _require_member(self.resampling_unit,
                                           _analysis.DESIGN_RESAMPLING_UNITS,
                                           "resampling_unit"))
        envelopes = tuple(
            require_exact_type(value, ArmEnvelope, f"arm_envelopes[{position}]")
            for position, value in enumerate(tuple(self.arm_envelopes)))
        names = [envelope.name for envelope in envelopes]
        if len(set(names)) != len(names):
            raise ValueError(f"arm envelopes repeat a name: {names}")
        for required in ("full", "central_control", "width_only_control"):
            if required not in names:
                raise ValueError(
                    f"arm_envelopes must contain {required!r}; the plan makes "
                    f"the central-clock and fixed-contraction width-only "
                    f"controls mandatory, and a matrix sized without them is a "
                    f"matrix that cannot pay for its own causal rule")
        object.__setattr__(self, "arm_envelopes", envelopes)
        object.__setattr__(self, "shadow_fraction",
                           _require_fraction(self.shadow_fraction,
                                             "shadow_fraction"))
        object.__setattr__(self, "maximum_trials_per_cell",
                           require_positive_int(self.maximum_trials_per_cell,
                                                "maximum_trials_per_cell",
                                                minimum=1))
        object.__setattr__(self, "stage",
                           _require_member(self.stage, POWER_STAGES, "stage"))
        floor = require_finite(self.contrast_covariance_floor,
                               "contrast_covariance_floor")
        if not -1.0 <= floor <= 1.0:
            raise ValueError(f"contrast_covariance_floor={floor!r} is not a "
                             f"correlation")
        if floor > 0.0:
            raise ValueError(
                f"contrast_covariance_floor={floor!r} is positive, which "
                f"assumes the two arms help each other.  An independent review "
                f"found exactly that assumption stated as 'pairing can only "
                f"reduce' the contrast variance, and it is false: two arms on "
                f"one Brownian tree can be negatively correlated in their "
                f"fitted exponents, and then pairing inflates the contrast")
        object.__setattr__(self, "contrast_covariance_floor", floor)
        level = _require_fraction(self.confidence_level, "confidence_level")
        if level != CONFIDENCE_LEVEL:
            raise ValueError(
                f"confidence_level={level!r} differs from the package's "
                f"{CONFIDENCE_LEVEL!r}; two layers disagreeing about what a "
                f"confidence level means is a defect nobody would find")
        object.__setattr__(self, "confidence_level", level)

    @property
    def arms(self) -> tuple:
        """The arm names, in declared order."""
        return tuple(envelope.name for envelope in self.arm_envelopes)

    @property
    def couplings(self) -> tuple:
        """The frozen coupling nodes: log-spaced, endpoints included."""
        return tuple(float(value) for value in np.exp(
            np.linspace(math.log(self.coupling_low),
                        math.log(self.coupling_high), self.cells)))

    @property
    def lever_arm(self) -> float:
        """``sum (log K_i - mean log K)^2`` — the sweep's whole leverage."""
        logk = np.log(np.asarray(self.couplings, dtype=float))
        return float(((logk - logk.mean()) ** 2).sum())

    @property
    def digest(self) -> str:
        return _digest({
            "stage": self.stage,
            "exponent_half_width": self.exponent_half_width,
            "minimum_causal_contrast": self.minimum_causal_contrast,
            "contrast_covariance_floor": self.contrast_covariance_floor,
            "coupling_low": self.coupling_low,
            "coupling_high": self.coupling_high, "cells": self.cells,
            "pairing_inflation": self.pairing_inflation,
            "conservative_factor": self.conservative_factor,
            "resampling_unit": self.resampling_unit,
            "arm_envelopes": [envelope.digest
                              for envelope in self.arm_envelopes],
            "shadow_fraction": self.shadow_fraction,
            "maximum_trials_per_cell": self.maximum_trials_per_cell,
            "confidence_level": self.confidence_level,
        })


@dataclass(frozen=True)
class PowerEstimate:
    """How many trials the frozen target needs, arm by arm, and what it buys.

    ``stage``
        Carried from the target.  A ``provisional_pre_pilot`` estimate is a
        planning number about an assumed range; it is **not** a validated
        all-arm sample size and every consumer of this record is required to
        say so.
    ``verdict``
        ``powered`` or ``no_powered_matrix``.
    ``arm_rows``
        ``(name, weight, trials_per_cell)`` per arm — **each arm's own count**,
        so that a conservative envelope on one control does not silently
        inflate every arm.  ``trials_per_cell`` below is the count the *paired*
        arms carry, because the full and width-only arms must share master
        trials for their contrast to be paired at all.
    ``analytic_independent_trials``
        What the exact independent-cell design variance alone would need, using
        the **full** arm's envelope.  Kept beside the sized number so the three
        contributions — analytic design, measured pairing inflation, declared
        conservative factor — can be read apart.
    ``physical_trial_executions``, ``shadow_rewalks``
        Counted separately, because a shadow rewalk is validation work over an
        existing master trial and not a new independent trial.  An independent
        review found the two added together and called "total trials".
    ``selection_digest``
        **The range this estimate was sized from**, or empty before a pilot.
        An independent review pointed out that a post-range power record could
        not be joined back to the selection that sized it, so a forged range
        and an exact-typed power record labelled ``post_range_join`` could be
        presented together with nothing to contradict them.  A post-range
        estimate now carries the digest of the selection it came from, and the
        relation is checked here: a provisional estimate has none and a joined
        one must have one.
    """

    target_digest: str
    stage: str
    verdict: str
    reasons: tuple
    arm_rows: tuple
    analytic_independent_trials: int
    pairing_inflation: float
    conservative_factor: float
    trials_per_cell: int
    achieved_half_width: float
    contrast_half_width: float
    contrast_covariance_floor: float
    cells: int
    shadow_fraction: float
    physical_trial_executions: int
    shadow_rewalks: int
    selection_digest: str = ""

    def __post_init__(self):
        object.__setattr__(self, "target_digest",
                           _require_hex_digest(self.target_digest,
                                               "target_digest"))
        object.__setattr__(self, "stage",
                           _require_member(self.stage, POWER_STAGES, "stage"))
        object.__setattr__(self, "verdict",
                           _require_member(self.verdict, POWER_VERDICTS,
                                           "verdict"))
        reasons = tuple(_require_text(value, f"reasons[{position}]",
                                      maximum=200)
                        for position, value in enumerate(tuple(self.reasons)))
        if (self.verdict == "powered") != (not reasons):
            raise ValueError(
                f"verdict {self.verdict!r} with reasons {list(reasons)}: a "
                f"powered design has no refusal reasons and a refusal has at "
                f"least one")
        object.__setattr__(self, "reasons", reasons)
        rows = []
        for position, entry in enumerate(tuple(self.arm_rows)):
            row = tuple(entry)
            if len(row) != 3:
                raise ValueError(f"arm_rows[{position}] must be a "
                                 f"(name, weight, trials) triple, got {row!r}")
            rows.append((_require_text(row[0], f"arm_rows[{position}] name",
                                       maximum=40),
                         _require_positive(row[1],
                                           f"arm_rows[{position}] weight"),
                         require_positive_int(row[2],
                                              f"arm_rows[{position}] trials",
                                              minimum=1)))
        if not rows:
            raise ValueError("a power estimate sizes at least one arm")
        object.__setattr__(self, "arm_rows", tuple(rows))
        for field in ("analytic_independent_trials", "trials_per_cell",
                      "cells", "physical_trial_executions"):
            object.__setattr__(self, field,
                               require_positive_int(getattr(self, field),
                                                    field, minimum=1))
        object.__setattr__(self, "shadow_rewalks",
                           require_positive_int(self.shadow_rewalks,
                                                "shadow_rewalks", minimum=0))
        for field in ("pairing_inflation", "conservative_factor",
                      "achieved_half_width", "contrast_half_width"):
            object.__setattr__(self, field,
                               _require_positive(getattr(self, field), field))
        object.__setattr__(self, "contrast_covariance_floor",
                           require_finite(self.contrast_covariance_floor,
                                          "contrast_covariance_floor"))
        object.__setattr__(self, "shadow_fraction",
                           _require_fraction(self.shadow_fraction,
                                             "shadow_fraction"))
        if self.stage == "post_range_join":
            object.__setattr__(self, "selection_digest",
                               _require_hex_digest(self.selection_digest,
                                                   "selection_digest"))
        elif self.selection_digest != "":
            raise ValueError(
                f"a {self.stage!r} estimate carries selection digest "
                f"{self.selection_digest!r}; before a pilot there is no "
                f"selected range, and a provisional estimate that names one "
                f"is claiming a join that did not happen")
        if self.trials_per_cell < self.analytic_independent_trials:
            raise ValueError(
                f"trials_per_cell={self.trials_per_cell} is below the "
                f"analytic independent requirement "
                f"{self.analytic_independent_trials}; pairing, arm envelopes "
                f"and the conservative factor all inflate the count and none "
                f"of them can reduce it")

    @property
    def is_provisional(self) -> bool:
        """``True`` until a count-only range selection has joined to it."""
        return self.stage == "provisional_pre_pilot"

    @property
    def driving_arm(self) -> str:
        """The arm whose envelope needs the most trials."""
        return max(self.arm_rows, key=lambda row: row[2])[0]

    @property
    def arm_trials(self) -> dict:
        """``{arm name: trials per cell}`` — what each arm actually runs."""
        return {name: trials for name, _, trials in self.arm_rows}

    @property
    def digest(self) -> str:
        return _digest({
            "target_digest": self.target_digest, "stage": self.stage,
            "verdict": self.verdict, "reasons": list(self.reasons),
            "arm_rows": [list(row) for row in self.arm_rows],
            "analytic_independent_trials": self.analytic_independent_trials,
            "pairing_inflation": self.pairing_inflation,
            "conservative_factor": self.conservative_factor,
            "trials_per_cell": self.trials_per_cell,
            "achieved_half_width": self.achieved_half_width,
            "contrast_half_width": self.contrast_half_width,
            "contrast_covariance_floor": self.contrast_covariance_floor,
            "cells": self.cells, "shadow_fraction": self.shadow_fraction,
            "physical_trial_executions": self.physical_trial_executions,
            "shadow_rewalks": self.shadow_rewalks,
            "selection_digest": self.selection_digest,
        })


def _size_from_target(target, couplings, lever_arm):
    """Shared sizing arithmetic for both power stages.

    With a constant per-trial weight ``w`` for an arm, the information matrix
    is ``n w [[m, sum x], [sum x, sum x^2]]`` and the slope variance collapses
    to ``1 / (n w S_xx)`` with ``S_xx`` the lever arm.  So

    ``n = z^2 / (w S_xx half_width^2)``,

    inflated by the declared pairing factor, the arm's own failure inflation,
    and the declared conservative factor.  Every one of those four is a
    separate number and every one is reported.
    """
    rows = []
    for envelope in target.arm_envelopes:
        weight = envelope.weight
        base = (CONFIDENCE_Z ** 2) / (
            weight * lever_arm * target.exponent_half_width ** 2)
        needed = (base * target.pairing_inflation * envelope.failure_inflation
                  * target.conservative_factor)
        rows.append((envelope.name, weight, max(1, int(math.ceil(needed)))))
    return tuple(rows)


def _finish_power(target, couplings, lever_arm, stage,
                  selection_digest=""):
    """Assemble a :class:`PowerEstimate` from a target and a coupling grid.

    Each arm is sized from its own envelope.  The arms that declare
    ``pairs_with_full`` are then raised to their common maximum, because a
    paired contrast between two arms with different trial counts is not a
    paired contrast; an arm that does not pair keeps its own count.
    """
    sized = dict((name, trials) for name, _, trials
                 in _size_from_target(target, couplings, lever_arm))
    paired_names = [envelope.name for envelope in target.arm_envelopes
                    if envelope.pairs_with_full]
    trials = max(sized[name] for name in paired_names) if paired_names else \
        max(sized.values())
    for name in paired_names:
        sized[name] = trials
    rows = tuple((envelope.name, envelope.weight, sized[envelope.name])
                 for envelope in target.arm_envelopes)

    full = next(envelope for envelope in target.arm_envelopes
                if envelope.name == "full")
    analytic = max(1, int(math.ceil(
        (CONFIDENCE_Z ** 2) / (full.weight * lever_arm
                               * target.exponent_half_width ** 2))))

    reasons = []
    # **Every arm, not the paired count.**  An independent review pointed out
    # that this compared only the paired arms' shared count against the
    # affordability ceiling, so a post-range join that sized the unpaired
    # central control at 22 211 trials per cell against a declared maximum of
    # 20 000 was reported ``powered``.  A matrix is affordable when *no* arm
    # exceeds the ceiling: the central control is a run that has to happen too.
    over = sorted((name, count) for name, count in sized.items()
                  if count > target.maximum_trials_per_cell)
    if over:
        reasons.append(
            f"{', '.join(f'{name} at {count}' for name, count in over)} "
            f"exceed the affordability ceiling "
            f"{target.maximum_trials_per_cell}; it binds every arm")

    # The achieved half-width, recomputed at the rounded count under the
    # weakest *paired* arm's own weight -- the contrast is between those two,
    # so the arm that is not paired does not set the contrast's precision.
    paired_weight = min((envelope.weight for envelope in target.arm_envelopes
                         if envelope.pairs_with_full),
                        default=min(row[1] for row in rows))
    achieved = CONFIDENCE_Z * math.sqrt(
        target.pairing_inflation / (trials * paired_weight * lever_arm))
    contrast = achieved * math.sqrt(
        2.0 * (1.0 - target.contrast_covariance_floor))
    if not contrast < target.minimum_causal_contrast:
        reasons.append(
            f"the full-minus-control contrast half-width {contrast:.4f} at the "
            f"declared adverse covariance {target.contrast_covariance_floor} "
            f"does not clear the frozen minimum causal contrast "
            f"{target.minimum_causal_contrast}")

    # Counted per arm, and the shadow is counted per *run* and then summed --
    # an independent review found a 5% shadow taken on the grand total, which
    # is a different and smaller number than 5% of each run rounded up.
    physical = sum(len(couplings) * sized[envelope.name]
                   for envelope in target.arm_envelopes)
    shadow = sum(len(couplings)
                 * int(math.ceil(sized[envelope.name]
                                 * target.shadow_fraction))
                 for envelope in target.arm_envelopes)
    return PowerEstimate(
        target_digest=target.digest, stage=stage,
        verdict="powered" if not reasons else "no_powered_matrix",
        reasons=tuple(reasons), arm_rows=rows,
        analytic_independent_trials=analytic,
        pairing_inflation=target.pairing_inflation,
        conservative_factor=target.conservative_factor,
        trials_per_cell=trials, achieved_half_width=achieved,
        contrast_half_width=contrast,
        contrast_covariance_floor=target.contrast_covariance_floor,
        cells=len(couplings), shadow_fraction=target.shadow_fraction,
        physical_trial_executions=physical, shadow_rewalks=shadow,
        selection_digest=selection_digest)


def power_estimate(target) -> PowerEstimate:
    """Size the sweep from ``target`` alone — the **provisional** stage.

    It assumes the coupling range in the target, which the count-only pilot has
    not chosen yet.  The record it returns says ``provisional_pre_pilot`` and
    every consumer carries that word.
    """
    require_exact_type(target, SamplingTarget, "target")
    if target.stage != "provisional_pre_pilot":
        raise ValueError(
            f"power_estimate produces the provisional stage; this target is "
            f"{target.stage!r} and belongs to post_range_power")
    return _finish_power(target, target.couplings, target.lever_arm,
                         "provisional_pre_pilot")


def post_range_power(target, plan, firewall, run_names) -> PowerEstimate:
    """Re-size the sweep from a count-only range selection.

    **The frozen mechanical post-range join.**  It takes the frozen
    :class:`PilotPlan`, the :class:`PilotFirewall` it was frozen against and
    the ordered quarantined run names, and derives everything from **one**
    authoritative reopen-and-recount snapshot: the counts, the selection the
    frozen rule makes from them, and the re-sized estimate.  There is no
    argument through which a caller's counts or a caller's selection could
    arrive, and none through which an exponent, a curvature, a survival shape,
    a commit time or any other forbidden observable could — the recount reads
    four numbers per cell and the closure digest, and nothing else.

    Two things change relative to the provisional estimate:

    1. the coupling nodes become the selected range's own eligible cells, so
       the lever arm is the one the experiment will actually have; and
    2. nothing else.  The arm envelopes, the inflation, the conservative
       factor, the contrast floor and the half-width target are all still the
       predeclared ones.  The pilot narrows the *range*; it does not get to
       narrow the rules.

    Zero-event and all-event cells inside the selected range are **kept**: they
    enter the returned cell list and therefore the eventual likelihood, exactly
    as ``analysis.fit_exponent`` requires.  They simply carry little weight,
    which the envelope already accounts for.

    The returned estimate carries the digest of the selection that sized it, so
    it can be joined back to the range rather than merely labelled with a
    stage.
    """
    return _post_range_snapshot(target, plan, firewall, run_names)[1]


def _power_from_snapshot(target, selection, pairs) -> PowerEstimate:
    """Re-size from one already-taken authoritative snapshot."""
    require_exact_type(target, SamplingTarget, "target")
    cells = [item for item, _, _ in pairs]
    if not selection.selected:
        # A refusal is a result, not an exception: the count rule declined to
        # name an information-bearing range, so there is no powered matrix.
        return PowerEstimate(
            target_digest=target.digest, stage="post_range_join",
            verdict="no_powered_matrix",
            reasons=(f"the count-only rule refused with "
                     f"{selection.refusal!r}, so there is no selected range to "
                     f"re-size from",),
            arm_rows=_size_from_target(target, (1.0,), 1.0),
            analytic_independent_trials=1,
            pairing_inflation=target.pairing_inflation,
            conservative_factor=target.conservative_factor,
            trials_per_cell=target.maximum_trials_per_cell,
            achieved_half_width=target.exponent_half_width,
            contrast_half_width=target.minimum_causal_contrast,
            contrast_covariance_floor=target.contrast_covariance_floor,
            cells=1, shadow_fraction=target.shadow_fraction,
            physical_trial_executions=1, shadow_rewalks=0,
            selection_digest=selection.digest)

    inside = tuple(sorted(item.coupling for item in cells
                          if selection.low <= item.coupling <= selection.high))
    if len(inside) < _analysis.MIN_CELLS:
        return PowerEstimate(
            target_digest=target.digest, stage="post_range_join",
            verdict="no_powered_matrix",
            reasons=(f"the selected range holds {len(inside)} cells, below the "
                     f"{_analysis.MIN_CELLS} the frozen model needs before an "
                     f"exponent and its curvature extension are identified",),
            arm_rows=_size_from_target(target, inside or (1.0,), 1.0),
            analytic_independent_trials=1,
            pairing_inflation=target.pairing_inflation,
            conservative_factor=target.conservative_factor,
            trials_per_cell=target.maximum_trials_per_cell,
            achieved_half_width=target.exponent_half_width,
            contrast_half_width=target.minimum_causal_contrast,
            contrast_covariance_floor=target.contrast_covariance_floor,
            cells=max(1, len(inside)), shadow_fraction=target.shadow_fraction,
            physical_trial_executions=1, shadow_rewalks=0,
            selection_digest=selection.digest)

    logk = np.log(np.asarray(inside, dtype=float))
    lever = float(((logk - logk.mean()) ** 2).sum())
    joined = SamplingTarget(
        exponent_half_width=target.exponent_half_width,
        minimum_causal_contrast=target.minimum_causal_contrast,
        coupling_low=inside[0], coupling_high=inside[-1], cells=len(inside),
        pairing_inflation=target.pairing_inflation,
        conservative_factor=target.conservative_factor,
        resampling_unit=target.resampling_unit,
        arm_envelopes=target.arm_envelopes,
        shadow_fraction=target.shadow_fraction,
        maximum_trials_per_cell=target.maximum_trials_per_cell,
        stage="post_range_join",
        contrast_covariance_floor=target.contrast_covariance_floor)
    return _finish_power(joined, inside, lever, "post_range_join",
                         selection.digest)


@dataclass(frozen=True)
class RecoveryReport:
    """A frozen synthetic exponent, recovered under a declared dependence.

    Three ensembles over the *same* declared truth:

    ``independent_stdev``
        Every cell draws its own randomness, in its own namespace.
    ``paired_stdev``
        Cells share one master-trial latent drawn from the **keyed physical
        noise stream**, combined with an independent per-cell component through
        a Gaussian copula at ``correlation``.  At ``correlation = 1`` this is
        the perfectly monotone nested-uniform case an earlier revision used
        alone; below one it produces genuinely non-monotone realizations across
        coupling, which is what an independent review pointed out the earlier
        generator could not.
    ``contrast_stdev``
        The spread of a full-minus-control exponent *difference* under the
        declared between-arm correlation, so that adverse covariance is
        measured rather than assumed away.

    Every statistic is computed over **valid fits only**.  An earlier revision
    included refused fits in its mean and standard deviation, which reports the
    fitter's failure mode as though it were the design's variance.
    """

    exponent_under_test: float
    trials_per_cell: int
    replicates: int
    correlation: float
    arm_correlation: float
    paired_mean: float
    paired_stdev: float
    independent_mean: float
    independent_stdev: float
    contrast_mean: float
    contrast_stdev: float
    predicted_stdev: float
    paired_valid: int
    independent_valid: int
    contrast_valid: int

    def __post_init__(self):
        object.__setattr__(self, "exponent_under_test",
                           require_finite(self.exponent_under_test,
                                          "exponent_under_test"))
        for field in ("trials_per_cell", "replicates"):
            object.__setattr__(self, field,
                               require_positive_int(getattr(self, field),
                                                    field, minimum=2))
        for field in ("correlation", "arm_correlation"):
            value = require_finite(getattr(self, field), field)
            if not -1.0 <= value <= 1.0:
                raise ValueError(f"{field}={value!r} is not a correlation")
            object.__setattr__(self, field, value)
        for field in ("paired_mean", "independent_mean", "contrast_mean"):
            object.__setattr__(self, field,
                               require_finite(getattr(self, field), field))
        for field in ("paired_stdev", "independent_stdev", "contrast_stdev",
                      "predicted_stdev"):
            object.__setattr__(self, field,
                               _require_positive(getattr(self, field), field))
        for field in ("paired_valid", "independent_valid", "contrast_valid"):
            value = require_positive_int(getattr(self, field), field,
                                         minimum=2)
            if value > self.replicates:
                raise ValueError(f"{field}={value} exceeds "
                                 f"replicates={self.replicates}")
            object.__setattr__(self, field, value)

    @property
    def inflation(self) -> float:
        """``Var(paired) / Var(independent)`` — the resource price of CRN."""
        return (self.paired_stdev / self.independent_stdev) ** 2

    @property
    def contrast_inflation(self) -> float:
        """``Var(contrast) / (2 Var(paired))``.

        One when the two arms are independent, above one when they are
        adversely correlated.  This is the number the declared
        :data:`CONTRAST_COVARIANCE_FLOOR` has to cover.
        """
        return (self.contrast_stdev ** 2) / (2.0 * self.paired_stdev ** 2)

    @property
    def invalid_fraction(self) -> float:
        """The worst refused-fit fraction across the three ensembles."""
        worst = min(self.paired_valid, self.independent_valid,
                    self.contrast_valid)
        return 1.0 - worst / self.replicates

    @property
    def digest(self) -> str:
        return _digest({
            "exponent_under_test": self.exponent_under_test,
            "trials_per_cell": self.trials_per_cell,
            "replicates": self.replicates, "correlation": self.correlation,
            "arm_correlation": self.arm_correlation,
            "paired_mean": self.paired_mean, "paired_stdev": self.paired_stdev,
            "independent_mean": self.independent_mean,
            "independent_stdev": self.independent_stdev,
            "contrast_mean": self.contrast_mean,
            "contrast_stdev": self.contrast_stdev,
            "predicted_stdev": self.predicted_stdev,
            "paired_valid": self.paired_valid,
            "independent_valid": self.independent_valid,
            "contrast_valid": self.contrast_valid,
        })


def _normal_cdf(values):
    """The standard normal distribution function, vectorized."""
    return 0.5 * (1.0 + np.vectorize(math.erf)(
        np.asarray(values, dtype=float) / math.sqrt(2.0)))


def recover_exponent(couplings, amplitude, exponent, trials_per_cell,
                     replicates, seed, correlation=1.0,
                     arm_correlation=CONTRAST_COVARIANCE_FLOOR,
                     control_exponent=1.0) -> RecoveryReport:
    """Recover a declared exponent under a declared dependence structure.

    The data-generating process is the frozen model itself: cell ``i`` commits
    with probability ``1 - exp(-exp(amplitude + exponent log K_i))``.  What
    changed after an independent review is *how the randomness is shared*.

    * The **shared master-trial latent** comes from
      ``stochastic.PhaseNoiseStream.standard_normals`` — the package's own
      counter-keyed physical source, keyed by trial and clock exactly as a real
      sweep keys it — rather than from an ad-hoc uniform.  One draw per master
      trial, reused at every coupling, which is what a shared Brownian path
      does.
    * Each cell adds its **own independent** standard normal, and the two are
      combined as ``rho * shared + sqrt(1 - rho^2) * own``.  At ``rho = 1`` the
      outcomes are perfectly monotone in coupling, which is the nested-uniform
      case; below one they are not, which is the case the earlier generator
      could not produce and which an independent review specifically asked for.
    * A **control arm** is generated beside the full arm at
      ``control_exponent``, sharing the master-trial latent with a declared
      ``arm_correlation``.  Its exponent is fitted independently and the
      full-minus-control difference is retained, so the contrast's spread under
      *adverse* covariance is measured rather than assumed.

    Every ensemble is fitted with ``analysis.fit_exponent`` — the production
    estimator — and only replicates whose verdict is ``valid_exponent``
    contribute to the reported means and standard deviations.  Zero-event and
    all-event cells are retained in every design, as the analysis layer
    requires; they are never dropped to make a fit succeed.
    """
    nodes = _require_couplings(couplings, "couplings")
    amplitude = require_finite(amplitude, "amplitude")
    exponent = require_finite(exponent, "exponent")
    control_exponent = require_finite(control_exponent, "control_exponent")
    trials_per_cell = require_positive_int(trials_per_cell, "trials_per_cell",
                                           minimum=2)
    replicates = require_positive_int(replicates, "replicates", minimum=2)
    seed = require_positive_int(seed, "seed", minimum=0)
    correlation = require_finite(correlation, "correlation")
    arm_correlation = require_finite(arm_correlation, "arm_correlation")
    for name, value in (("correlation", correlation),
                        ("arm_correlation", arm_correlation)):
        if not -1.0 <= value <= 1.0:
            raise ValueError(f"{name}={value!r} is not a correlation")

    logk = np.log(np.asarray(nodes, dtype=float))
    for label, slope in (("exponent", exponent),
                         ("control_exponent", control_exponent)):
        # Checked on the *linear predictor* rather than after exponentiating,
        # so an overflowing design is refused without raising a warning on the
        # way -- the suite runs one path under ``-W error``.
        predictor = amplitude + slope * logk
        if not np.all(np.isfinite(predictor)) or float(predictor.max()) > 700.0:
            raise ValueError(f"amplitude={amplitude!r} and {label}={slope!r} "
                             f"overflow the cumulative hazard at one of "
                             f"{nodes!r}")
    thresholds = _normal_inverse(-np.expm1(-np.exp(amplitude
                                                   + exponent * logk)))
    control_thresholds = _normal_inverse(
        -np.expm1(-np.exp(amplitude + control_exponent * logk)))

    mesh = _stochastic.FinestMesh(0.0, 1.0)
    stream = _stochastic.PhaseNoiseStream(f"t07-recovery-{seed}", mesh, 0.5)
    rng = np.random.default_rng(seed)
    shared_ids = _analysis.master_trials("experiments-crn", trials_per_cell)
    control_ids = _analysis.master_trials("experiments-crn-control",
                                          trials_per_cell)
    per_cell = tuple(_analysis.master_trials(f"experiments-indep-{index}",
                                             trials_per_cell)
                     for index in range(len(nodes)))

    def fitted(cells, unit):
        report = _analysis.fit_exponent(
            _analysis.ScalingDesign(tuple(cells), unit))
        return (report.exponent, report.is_valid)

    root = math.sqrt(max(0.0, 1.0 - correlation * correlation))
    arm_root = math.sqrt(max(0.0, 1.0 - arm_correlation * arm_correlation))
    paired, independent, contrast = [], [], []
    for replicate in range(replicates):
        # The shared master-trial latent, from the package's keyed source.
        shared = stream.standard_normals(replicate, 0, 0, trials_per_cell)
        control_shared = (arm_correlation * shared
                          + arm_root * stream.standard_normals(
                              replicate, 1, 0, trials_per_cell))

        full_cells, control_cells = [], []
        for index, node in enumerate(nodes):
            own = rng.standard_normal(trials_per_cell)
            latent = correlation * shared + root * own
            full_cells.append(_analysis.CouplingCell(
                node, shared_ids,
                tuple(bool(value) for value in latent < thresholds[index])))
            own = rng.standard_normal(trials_per_cell)
            latent = correlation * control_shared + root * own
            control_cells.append(_analysis.CouplingCell(
                node, control_ids,
                tuple(bool(value)
                      for value in latent < control_thresholds[index])))
        full_value, full_ok = fitted(full_cells, "master_trial")
        control_value, control_ok = fitted(control_cells, "master_trial")
        if full_ok:
            paired.append(full_value)
        if full_ok and control_ok:
            contrast.append(full_value - control_value)

        cells = []
        for index, node in enumerate(nodes):
            draw = rng.standard_normal(trials_per_cell)
            cells.append(_analysis.CouplingCell(
                node, per_cell[index],
                tuple(bool(value) for value in draw < thresholds[index])))
        value, ok = fitted(cells, "independent_namespace")
        if ok:
            independent.append(value)

    arrays = {}
    for label, values in (("paired", paired), ("independent", independent),
                          ("contrast", contrast)):
        block = np.asarray(values, dtype=float)
        if block.size < 2 or not np.all(np.isfinite(block)):
            raise ValueError(
                f"the {label} ensemble produced {block.size} valid finite "
                f"fits out of {replicates}; a recovery study with nothing left "
                f"to summarize measures the fitter's failure mode rather than "
                f"the design's variance")
        arrays[label] = block

    predicted = math.sqrt(cloglog_slope_variance(
        nodes, (trials_per_cell,) * len(nodes), amplitude, exponent))
    return RecoveryReport(
        exponent_under_test=exponent, trials_per_cell=trials_per_cell,
        replicates=replicates, correlation=correlation,
        arm_correlation=arm_correlation,
        paired_mean=float(arrays["paired"].mean()),
        paired_stdev=float(arrays["paired"].std(ddof=1)),
        independent_mean=float(arrays["independent"].mean()),
        independent_stdev=float(arrays["independent"].std(ddof=1)),
        contrast_mean=float(arrays["contrast"].mean()),
        contrast_stdev=float(arrays["contrast"].std(ddof=1)),
        predicted_stdev=predicted, paired_valid=int(arrays["paired"].size),
        independent_valid=int(arrays["independent"].size),
        contrast_valid=int(arrays["contrast"].size))


def _normal_inverse(probability):
    """The standard normal quantile, by bisection on :func:`_normal_cdf`.

    Written out rather than imported because this module may not depend on
    SciPy and ``statistics.NormalDist`` would be a second opinion about the
    same function.  Bisection to 1e-12 is far below anything a Monte Carlo
    threshold notices, and it is monotone, which is the property that matters:
    a higher commitment probability must give a higher threshold.
    """
    values = np.asarray(probability, dtype=float)
    low = np.full(values.shape, -12.0)
    high = np.full(values.shape, 12.0)
    for _ in range(80):
        middle = 0.5 * (low + high)
        below = _normal_cdf(middle) < values
        low = np.where(below, middle, low)
        high = np.where(below, high, middle)
    return 0.5 * (low + high)


# ---------------------------------------------------------------------------
# The numerical budget, tied to the planned statistical uncertainty
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class NumericalEvidence:
    """One measured numerical discrepancy, as ticket 04 reported it.

    ``measured``, ``standard_error``
        The discrepancy at ``timestep`` and its bootstrap standard error, kept
        **apart**.  An independent review found the two being scaled together
        under a hypothetical refinement, which is wrong in the optimistic
        direction: refining the timestep moves the *bias*, and enlarging the
        sample moves the *standard error*, and no single factor does both.
        :meth:`projected_bound` scales each by its own lever.
    ``sample_clusters``
        How many independent clusters the standard error was computed from, so
        that a projected enlargement has something to divide.
    ``verdict``
        The ticket-04 verdict the evidence came with.  ``numerical_no_result``
        is carried through unchanged and is blocking on its own.
    ``configuration``
        One of :data:`CONFIGURATION_LABELS`.  ``non_intended`` is not an error;
        it is a blocker, because transferring a stationary reference cell's
        envelope to a pulsed production grid is an assumption and not a
        measurement.
    """

    source: str
    observable: str
    unit: str
    position: str
    timestep: float
    measured: float
    standard_error: float
    sample_clusters: int
    verdict: str
    configuration: str
    description: str

    SOURCES = ("stationary_oracle", "moving_band_audit")

    def __post_init__(self):
        object.__setattr__(self, "source",
                           _require_member(self.source,
                                           NumericalEvidence.SOURCES,
                                           "source"))
        object.__setattr__(self, "observable",
                           _require_text(self.observable, "observable",
                                         maximum=80))
        object.__setattr__(self, "unit",
                           _require_member(self.unit, _oracle.UNITS, "unit"))
        object.__setattr__(self, "position",
                           _require_text(self.position, "position",
                                         maximum=80))
        object.__setattr__(self, "timestep",
                           _require_positive(self.timestep, "timestep"))
        object.__setattr__(self, "measured",
                           require_finite(self.measured, "measured"))
        object.__setattr__(self, "standard_error",
                           require_nonnegative(self.standard_error,
                                               "standard_error"))
        object.__setattr__(self, "sample_clusters",
                           require_positive_int(self.sample_clusters,
                                                "sample_clusters", minimum=1))
        object.__setattr__(self, "verdict",
                           _require_member(self.verdict, _oracle.VERDICTS,
                                           "verdict"))
        object.__setattr__(self, "configuration",
                           _require_member(self.configuration,
                                           CONFIGURATION_LABELS,
                                           "configuration"))
        object.__setattr__(self, "description",
                           _require_text(self.description, "description",
                                         maximum=300))

    @property
    def matches_intended(self) -> bool:
        """Was this measured at the configuration the matrix proposes?"""
        return self.configuration == "intended"

    def bound_at(self, coverage_sigma) -> float:
        """``|measured| + coverage_sigma * standard_error``.

        The uncertainty-aware magnitude, and **the sigma is an argument**.  An
        independent review found this reading a module constant while the
        budget stored and hashed its own value, so two budgets differing only
        in coverage produced identical rows and identical dispositions.  The
        frozen number now has to be supplied by whoever froze it.
        """
        sigma = _require_positive(coverage_sigma, "coverage_sigma")
        return abs(self.measured) + sigma * self.standard_error

    def projected_bound(self, coverage_sigma, timestep_factor=1.0,
                        sample_factor=1.0) -> float:
        """The bound this row would have at a finer step and a larger sample.

        Two separate levers, because they act on two separate terms: the
        **bias** is assumed to fall like ``sqrt(timestep)`` — the rate the
        ticket-04 endpoint ladder measured and asserted — and the **standard
        error** falls like ``1/sqrt(sample_factor)`` and is untouched by
        refinement.  An independent review showed that scaling the whole bound
        by ``sqrt(dt)`` makes a refinement look like an escape when it is not.

        Both are planning projections and neither is evidence.
        """
        sigma = _require_positive(coverage_sigma, "coverage_sigma")
        timestep_factor = _require_positive(timestep_factor, "timestep_factor")
        sample_factor = _require_positive(sample_factor, "sample_factor")
        bias = abs(self.measured) * math.sqrt(timestep_factor)
        error = self.standard_error / math.sqrt(sample_factor)
        return bias + sigma * error

    @property
    def digest(self) -> str:
        return _digest({
            "source": self.source, "observable": self.observable,
            "unit": self.unit, "position": self.position,
            "timestep": self.timestep, "measured": self.measured,
            "standard_error": self.standard_error,
            "sample_clusters": self.sample_clusters, "verdict": self.verdict,
            "configuration": self.configuration,
            "description": self.description,
        })


@dataclass(frozen=True)
class NumericalBudget:
    """The production-specific allowance, derived from planned uncertainty.

    Frozen **before** any evidence is opened, and hashed, so that widening it
    is a visible act.

    ``probability_half_width``
        The planned 95-per-cent half-width of a per-cell commitment
        probability at the proposed trial count.  It shrinks like
        ``1/sqrt(n)``, so tolerating a larger numerical error means planning a
        *less* precise experiment.
    ``time_half_width``
        The planned half-width of a commit-time quantile.  It does **not**
        depend on the trial count, so a time-unit failure fails at every trial
        count — which is why collapsing it into a scalar "admissible n" makes
        the failure invisible, as an independent review pointed out.
    """

    label: str
    allowed_fraction: float
    coverage_sigma: float
    probability_half_width: float
    time_half_width: float
    trials_per_cell: int
    reference_probability: float
    intended_timestep: float

    def __post_init__(self):
        object.__setattr__(self, "label",
                           _require_text(self.label, "label", maximum=120))
        fraction = _require_fraction(self.allowed_fraction, "allowed_fraction")
        if not 0.0 < fraction <= 1.0:
            raise ValueError(f"allowed_fraction={fraction!r} must lie in "
                             f"(0, 1]; a zero allowance can never be met and "
                             f"is not a budget")
        object.__setattr__(self, "allowed_fraction", fraction)
        object.__setattr__(self, "coverage_sigma",
                           _require_positive(self.coverage_sigma,
                                             "coverage_sigma"))
        for field in ("probability_half_width", "time_half_width",
                      "intended_timestep"):
            object.__setattr__(self, field,
                               _require_positive(getattr(self, field), field))
        object.__setattr__(self, "trials_per_cell",
                           require_positive_int(self.trials_per_cell,
                                                "trials_per_cell", minimum=1))
        object.__setattr__(self, "reference_probability",
                           _require_fraction(self.reference_probability,
                                             "reference_probability"))

    def allowance(self, unit: str) -> float:
        """The allowance for one of :data:`killed_diffusion.UNITS`.

        ``count`` maps to the probability allowance: the ticket-04 count
        observables are edge-resolved *probabilities* of at least one exit.
        ``time`` additionally receives one intended timestep, the closed plan's
        rule for commit-time quantiles.
        """
        unit = _require_member(unit, _oracle.UNITS, "unit")
        if unit == "time":
            return (self.allowed_fraction * self.time_half_width
                    + self.intended_timestep)
        return self.allowed_fraction * self.probability_half_width

    def probability_allowance_at(self, trials_per_cell) -> float:
        """The probability allowance a design of ``trials_per_cell`` would have."""
        trials = require_positive_int(trials_per_cell, "trials_per_cell",
                                      minimum=1)
        spread = self.reference_probability * (1.0
                                               - self.reference_probability)
        return self.allowed_fraction * CONFIDENCE_Z * math.sqrt(spread / trials)

    @property
    def digest(self) -> str:
        return _digest({
            "label": self.label, "allowed_fraction": self.allowed_fraction,
            "coverage_sigma": self.coverage_sigma,
            "probability_half_width": self.probability_half_width,
            "time_half_width": self.time_half_width,
            "trials_per_cell": self.trials_per_cell,
            "reference_probability": self.reference_probability,
            "intended_timestep": self.intended_timestep,
        })


@dataclass(frozen=True)
class BudgetDisposition:
    """Whether the numerical evidence fits the planned statistical precision.

    ``rows``
        **Every** authoritative evidence row.  An independent review found the
        composition rule selecting one maximum per unit by *absolute point
        error*, which can drop a row with a smaller point estimate and a larger
        uncertainty-aware bound.  Nothing is selected now.
    ``probability_admissible_trials``
        The largest per-cell trial count at which the worst **probability**
        bound would still fit.  ``0`` means none.
    ``time_admissible``
        Whether every **time** row fits.  Time allowances do not depend on the
        trial count, so a time failure fails at every ``n``.
    """

    verdict: str
    blockers: tuple
    rows: tuple
    budget_digest: str
    evidence_digest: str
    probability_admissible_trials: int
    time_admissible: bool

    def __post_init__(self):
        object.__setattr__(self, "verdict",
                           _require_member(self.verdict, BUDGET_VERDICTS,
                                           "verdict"))
        blockers = tuple(self.blockers)
        unknown = [value for value in blockers if value not in BUDGET_BLOCKERS]
        if unknown:
            raise ValueError(
                f"blockers {unknown} are outside the closed vocabulary "
                f"{list(BUDGET_BLOCKERS)}; a consumer branches on these, so a "
                f"new one is a deliberate edit here rather than a string "
                f"invented at a call site")
        if (self.verdict == "satisfied") != (not blockers):
            raise ValueError(
                f"verdict {self.verdict!r} with blockers {list(blockers)}: a "
                f"satisfied budget has no blockers and a blocked one has at "
                f"least one")
        object.__setattr__(self, "blockers", blockers)
        object.__setattr__(self, "rows", tuple(tuple(row) for row in self.rows))
        for field in ("budget_digest", "evidence_digest"):
            object.__setattr__(self, field,
                               _require_hex_digest(getattr(self, field),
                                                   field))
        object.__setattr__(self, "probability_admissible_trials",
                           require_positive_int(
                               self.probability_admissible_trials,
                               "probability_admissible_trials", minimum=0))
        if not isinstance(self.time_admissible, (bool, np.bool_)):
            raise TypeError(f"time_admissible must be a bool, got "
                            f"{type(self.time_admissible).__name__}")
        object.__setattr__(self, "time_admissible", bool(self.time_admissible))

    @property
    def is_blocking(self) -> bool:
        """``True`` unless the verdict is ``satisfied``."""
        return self.verdict != "satisfied"

    @property
    def admissible_trials(self) -> int:
        """The largest admissible trial count **overall**.

        Zero whenever any time row fails, because a time failure is
        ``n``-independent.  Reading the probability scalar as an overall
        maximum is exactly the misreading an independent review flagged, so the
        overall number is derived here and the probability-only one keeps its
        own explicit name.
        """
        return (self.probability_admissible_trials if self.time_admissible
                else 0)

    def limiting(self, unit: str):
        """The row with the largest uncertainty-aware bound in one unit."""
        unit = _require_member(unit, _oracle.UNITS, "unit")
        candidates = [row for row in self.rows if row[3] == unit]
        if not candidates:
            return None
        return max(candidates, key=lambda row: row[4])

    @property
    def digest(self) -> str:
        return _digest({
            "verdict": self.verdict, "blockers": list(self.blockers),
            "rows": [list(row) for row in self.rows],
            "budget_digest": self.budget_digest,
            "evidence_digest": self.evidence_digest,
            "probability_admissible_trials":
                self.probability_admissible_trials,
            "time_admissible": self.time_admissible,
        })


def _evidence_digest(items):
    return _digest({"evidence": [item.digest for item in items]})


def numerical_disposition(budget, evidence) -> BudgetDisposition:
    """Compare ticket 04's measured evidence with this ticket's allowance.

    Every authoritative row is carried and judged; nothing is selected away.
    The bound uses **the budget's own frozen coverage sigma**.  Five ways to
    block, and none can be argued away by a point estimate:

    1. a row arrives with ticket 04's own ``numerical_no_result`` verdict —
       carried through, never re-decided;
    2. its uncertainty-aware bound exceeds the allowance;
    3. it was not measured at the intended configuration;
    4. every positive trial count fails the probability rule; or
    5. a time row fails, which fails at every trial count at once.
    """
    require_exact_type(budget, NumericalBudget, "budget")
    items = tuple(require_exact_type(item, NumericalEvidence,
                                     f"evidence[{position}]")
                  for position, item in enumerate(tuple(evidence)))
    if not items:
        return BudgetDisposition("unresolved", ("evidence_missing",), (),
                                 budget.digest, _evidence_digest(()), 0, False)

    rows = []
    blockers = []
    worst_probability = 0.0
    time_admissible = True
    for item in items:
        bound = item.bound_at(budget.coverage_sigma)
        allowance = budget.allowance(item.unit)
        fits = bound <= allowance
        rows.append((item.source, item.observable, item.position, item.unit,
                     bound, allowance, fits))
        if item.verdict == "numerical_no_result":
            blockers.append("moving_band_numerical_no_result"
                            if item.source == "moving_band_audit"
                            else "endpoint_envelope_exceeds_allowance")
        if not fits:
            blockers.append("endpoint_envelope_exceeds_allowance"
                            if item.source == "stationary_oracle"
                            else "audit_envelope_exceeds_allowance")
        if not item.matches_intended:
            blockers.append("no_evidence_at_intended_configuration")
        if item.unit == "time":
            time_admissible = time_admissible and fits
        else:
            worst_probability = max(worst_probability, bound)

    admissible = 0
    if worst_probability > 0.0:
        spread = budget.reference_probability * (
            1.0 - budget.reference_probability)
        ratio = budget.allowed_fraction * CONFIDENCE_Z / worst_probability
        admissible = int(math.floor(spread * ratio * ratio))
    if admissible < budget.trials_per_cell:
        blockers.append("probability_window_empty")
    if not time_admissible:
        blockers.append("time_bounds_fail_at_every_trial_count")

    order = {name: index for index, name in enumerate(BUDGET_BLOCKERS)}
    blockers = tuple(sorted(set(blockers), key=order.__getitem__))
    verdict = "satisfied"
    if any(item.verdict == "numerical_no_result" for item in items):
        verdict = "numerical_no_result"
    elif blockers:
        verdict = "unresolved"
    return BudgetDisposition(verdict, blockers, tuple(rows), budget.digest,
                             _evidence_digest(items), max(0, admissible),
                             time_admissible)


@dataclass(frozen=True)
class FeasibilityOption:
    """One declared point in the finite search over ways out of the block.

    ``timestep_factor``
        The proposed step as a fraction of the one the evidence was measured
        at.  ``0.125`` is "eight times finer".
    ``sample_factor``
        How many times the validation sample would be enlarged.
    ``trials_per_cell``
        The production trial count this option plans for.  Fewer trials means a
        wider planned half-width and a larger allowance, which is the trade the
        search exists to make visible.
    ``configuration``
        Whether the option's evidence would be at the intended configuration.
        An option that keeps ``non_intended`` evidence cannot be *feasible*
        however small its bounds get.
    """

    label: str
    timestep_factor: float
    sample_factor: float
    trials_per_cell: int
    configuration: str

    def __post_init__(self):
        object.__setattr__(self, "label",
                           _require_text(self.label, "label", maximum=120))
        timestep_factor = _require_positive(self.timestep_factor,
                                            "timestep_factor")
        if timestep_factor > 1.0:
            raise ValueError(
                f"timestep_factor={timestep_factor!r} is above one, which "
                f"proposes a *coarser* step than the evidence was measured "
                f"at; that is not a way out of a numerical block and the "
                f"search refuses to enumerate it")
        object.__setattr__(self, "timestep_factor", timestep_factor)
        sample_factor = _require_positive(self.sample_factor, "sample_factor")
        if sample_factor < 1.0:
            raise ValueError(
                f"sample_factor={sample_factor!r} is below one, which proposes "
                f"throwing away validation data")
        object.__setattr__(self, "sample_factor", sample_factor)
        object.__setattr__(self, "trials_per_cell",
                           require_positive_int(self.trials_per_cell,
                                                "trials_per_cell", minimum=1))
        object.__setattr__(self, "configuration",
                           _require_member(self.configuration,
                                           CONFIGURATION_LABELS,
                                           "configuration"))

    @property
    def digest(self) -> str:
        return _digest({"label": self.label,
                        "timestep_factor": self.timestep_factor,
                        "sample_factor": self.sample_factor,
                        "trials_per_cell": self.trials_per_cell,
                        "configuration": self.configuration})


@dataclass(frozen=True)
class SearchReport:
    """The result of a **finite, declared** search over feasibility options.

    ``verdict``
        ``feasible_matrix`` if some searched option fits every row at the
        intended configuration; ``unresolved`` if an option would fit but only
        on transferred evidence; otherwise
        ``no_feasible_matrix_among_searched``.

    The last one is the whole point of this record.  An independent review was
    right that the earlier implementation asserted no feasible
    timestep/sample/configuration exists without performing a search or proving
    a bound.  What can honestly be said is a statement about the options that
    were enumerated — and they are enumerated here, hashed, and reported row by
    row.
    """

    verdict: str
    rows: tuple
    budget_digest: str
    evidence_digest: str
    options_digest: str

    def __post_init__(self):
        object.__setattr__(self, "verdict",
                           _require_member(self.verdict, FEASIBILITY_VERDICTS,
                                           "verdict"))
        rows = tuple(tuple(row) for row in self.rows)
        if not rows:
            raise ValueError("a search over no options has searched nothing, "
                             "and reporting a verdict from it would be the "
                             "unqualified claim this record exists to replace")
        object.__setattr__(self, "rows", rows)
        for field in ("budget_digest", "evidence_digest", "options_digest"):
            object.__setattr__(self, field,
                               _require_hex_digest(getattr(self, field),
                                                   field))

    @property
    def searched(self) -> int:
        """How many declared options were examined."""
        return len(self.rows)

    def row(self, label: str):
        """One searched option's row, by label."""
        label = _require_text(label, "label", maximum=120)
        for entry in self.rows:
            if entry[0] == label:
                return entry
        raise KeyError(f"{label!r} was not among the searched options")

    @property
    def digest(self) -> str:
        return _digest({"verdict": self.verdict,
                        "rows": [list(row) for row in self.rows],
                        "budget_digest": self.budget_digest,
                        "evidence_digest": self.evidence_digest,
                        "options_digest": self.options_digest})


def search_feasible_options(budget, evidence, options) -> SearchReport:
    """Enumerate a finite declared set of escapes and report what each does.

    For every option, every evidence row is re-projected: the bias by
    ``sqrt(timestep_factor)``, the standard error by ``1/sqrt(sample_factor)``
    — separately, because they have different levers — and compared against the
    allowance the option's own trial count implies.

    Nothing here is evidence.  A row that "fits" under an option is a statement
    that *if* the projection held, the option would satisfy the budget; the
    projection is a planning assumption, and an option is labelled
    ``non_intended`` unless it also proposes measuring at the configuration the
    matrix runs.  That is why an option can fit every row and still leave the
    verdict ``unresolved``.
    """
    require_exact_type(budget, NumericalBudget, "budget")
    items = tuple(require_exact_type(item, NumericalEvidence,
                                     f"evidence[{position}]")
                  for position, item in enumerate(tuple(evidence)))
    declared = tuple(require_exact_type(item, FeasibilityOption,
                                        f"options[{position}]")
                     for position, item in enumerate(tuple(options)))
    if not items:
        raise ValueError("a search needs evidence to project; with none, every "
                         "option trivially fits and the report would be a "
                         "statement about nothing")
    if not declared:
        raise ValueError("a search over no options has searched nothing")

    rows = []
    feasible_intended = False
    feasible_transferred = False
    for option in declared:
        worst = 0.0
        fits = True
        for item in items:
            bound = item.projected_bound(budget.coverage_sigma,
                                         option.timestep_factor,
                                         option.sample_factor)
            allowance = (budget.allowance("time") if item.unit == "time"
                         else budget.probability_allowance_at(
                             option.trials_per_cell))
            fits = fits and bound <= allowance
            worst = max(worst, bound / allowance)
        rows.append((option.label, option.timestep_factor,
                     option.sample_factor, option.trials_per_cell,
                     option.configuration, worst, fits))
        if fits and option.configuration == "intended":
            feasible_intended = True
        elif fits:
            feasible_transferred = True

    if feasible_intended:
        verdict = "feasible_matrix"
    elif feasible_transferred:
        verdict = "unresolved"
    else:
        verdict = "no_feasible_matrix_among_searched"
    return SearchReport(verdict, tuple(rows), budget.digest,
                        _evidence_digest(items),
                        _digest({"options": [item.digest
                                             for item in declared]}))


# ---------------------------------------------------------------------------
# The finite primary matrix, and what it costs
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class PrimaryMatrix:
    """One finite production matrix proposal.  No factorial, no continuum.

    Everything the plan requires to be explicit is a field here, and an
    independent review's two findings about *identity* are closed by three of
    them:

    ``arm_clocks``
        ``((name, clocks), ...)`` — the **arm-specific grid**.  The central
        control is a one-clock population containing zero detuning; the full
        and width-only arms are the sixty-four-clock midpoint grid.  An earlier
        revision carried only the global grid, so the manifest did not describe
        the population one of its three arms would actually run.
    ``arm_trials``
        ``((name, trials_per_cell), ...)``.  Arms sized from different
        envelopes run different counts, and the paired arms run the same one.
    ``refinement_couplings``, ``refinement_arm``, ``refinement_trials_per_cell``
        **Exactly which** cells and **which** arm the timestep ladder refines,
        and at what count.  "2 of 6 cells, 1 of 3 arms" is not an execution
        graph: it leaves the choice to be made after data are opened, which is
        the one thing a frozen matrix exists to prevent.  The refinement count
        is separate because a convergence study does not need production
        precision, and pricing it at production precision was making the ladder
        as expensive as the sweep.
    """

    label: str
    couplings: tuple
    clocks: int
    arm_clocks: tuple
    arm_trials: tuple
    grid_half_width: float
    pulse_duration: float
    pulse_centre: float
    phase_diffusion: float
    lock_tolerance: float
    dwell_time: float
    timestep: float
    timestep_ladder: tuple
    refinement_couplings: tuple
    refinement_arm: str
    refinement_trials_per_cell: int
    trials_per_cell: int
    shadow_fraction: float
    minimum_eligible_clocks: int
    discrete_staircase_fallback: str
    stop_rule: str
    failure_rule: str

    def __post_init__(self):
        object.__setattr__(self, "label",
                           _require_text(self.label, "label", maximum=120))
        object.__setattr__(self, "couplings",
                           _require_couplings(self.couplings, "couplings"))
        object.__setattr__(self, "clocks",
                           require_positive_int(self.clocks, "clocks",
                                                minimum=1))
        for field in ("grid_half_width", "pulse_duration", "lock_tolerance",
                      "dwell_time", "timestep"):
            object.__setattr__(self, field,
                               _require_positive(getattr(self, field), field))
        object.__setattr__(self, "pulse_centre",
                           require_finite(self.pulse_centre, "pulse_centre"))
        object.__setattr__(self, "phase_diffusion",
                           require_nonnegative(self.phase_diffusion,
                                               "phase_diffusion"))
        ladder = tuple(_require_positive(value, f"timestep_ladder[{position}]")
                       for position, value in
                       enumerate(tuple(self.timestep_ladder)))
        if len(ladder) < _oracle.MINIMUM_REFINEMENT_LEVELS:
            raise ValueError(
                f"timestep_ladder holds {len(ladder)} levels; the numerical "
                f"contract needs at least "
                f"{_oracle.MINIMUM_REFINEMENT_LEVELS} paired levels before a "
                f"trend is a trend")
        if any(later >= earlier for earlier, later in zip(ladder, ladder[1:])):
            raise ValueError(f"timestep_ladder must be strictly decreasing, "
                             f"coarsest first, got {ladder!r}")
        if ladder[0] != self.timestep:
            raise ValueError(
                f"timestep_ladder starts at {ladder[0]!r} but the production "
                f"timestep is {self.timestep!r}; the ladder refines the step "
                f"the matrix runs at, and a ladder that starts somewhere else "
                f"measures a different experiment")
        object.__setattr__(self, "timestep_ladder", ladder)

        clocks_by_arm = {}
        order = []
        for position, entry in enumerate(tuple(self.arm_clocks)):
            pair = tuple(entry)
            if len(pair) != 2:
                raise ValueError(f"arm_clocks[{position}] must be a "
                                 f"(name, clocks) pair, got {pair!r}")
            name = _require_text(pair[0], f"arm_clocks[{position}] name",
                                 maximum=40)
            if name in clocks_by_arm:
                raise ValueError(f"arm {name!r} appears twice in arm_clocks")
            clocks_by_arm[name] = require_positive_int(
                pair[1], f"arm_clocks[{position}] clocks", minimum=1)
            order.append(name)
        if not order:
            raise ValueError("a matrix runs at least one arm")
        for required in ("full", "central_control", "width_only_control"):
            if required not in clocks_by_arm:
                raise ValueError(
                    f"arm_clocks must name {required!r}; the plan makes the "
                    f"central-clock and fixed-contraction width-only controls "
                    f"mandatory")
        if clocks_by_arm["central_control"] != 1:
            raise ValueError(
                f"the central-clock control runs "
                f"{clocks_by_arm['central_control']} clocks; it is *the one "
                f"clock at zero detuning*, and a multi-clock version of it is "
                f"a different control that tests something else")
        if clocks_by_arm["full"] != self.clocks:
            raise ValueError(
                f"the full arm runs {clocks_by_arm['full']} clocks while the "
                f"matrix declares a {self.clocks}-clock population")
        object.__setattr__(self, "arm_clocks",
                           tuple((name, clocks_by_arm[name]) for name in order))

        trials_by_arm = {}
        for position, entry in enumerate(tuple(self.arm_trials)):
            pair = tuple(entry)
            if len(pair) != 2:
                raise ValueError(f"arm_trials[{position}] must be a "
                                 f"(name, trials) pair, got {pair!r}")
            name = _require_text(pair[0], f"arm_trials[{position}] name",
                                 maximum=40)
            trials_by_arm[name] = require_positive_int(
                pair[1], f"arm_trials[{position}] trials", minimum=1)
        if sorted(trials_by_arm) != sorted(clocks_by_arm):
            raise ValueError(
                f"arm_trials names {sorted(trials_by_arm)} and arm_clocks "
                f"names {sorted(clocks_by_arm)}; every arm has exactly one "
                f"grid and exactly one trial count")
        object.__setattr__(self, "arm_trials",
                           tuple((name, trials_by_arm[name])
                                 for name in order))

        for field in ("refinement_trials_per_cell", "trials_per_cell",
                      "minimum_eligible_clocks"):
            object.__setattr__(self, field,
                               require_positive_int(getattr(self, field),
                                                    field, minimum=1))
        refinement = tuple(
            _require_positive(value, f"refinement_couplings[{position}]")
            for position, value in enumerate(tuple(self.refinement_couplings)))
        if len(refinement) < 2:
            raise ValueError(
                f"refinement_couplings holds {len(refinement)} node(s); a "
                f"convergence study at a single coupling cannot show that the "
                f"trend is the same across the sweep")
        if len(set(refinement)) != len(refinement):
            raise ValueError("refinement_couplings repeats a node")
        outside = [value for value in refinement if value not in self.couplings]
        if outside:
            raise ValueError(
                f"refinement_couplings {outside} are not matrix couplings; the "
                f"ladder refines cells the sweep runs, and a node outside it "
                f"is a different experiment")
        if len(refinement) >= len(self.couplings):
            raise ValueError(
                f"the ladder refines all {len(refinement)} couplings, which is "
                f"the automatic full factorial the plan forbids")
        object.__setattr__(self, "refinement_couplings", tuple(sorted(refinement)))
        arm = _require_text(self.refinement_arm, "refinement_arm", maximum=40)
        if arm not in clocks_by_arm:
            raise ValueError(f"refinement_arm={arm!r} is not one of "
                             f"{sorted(clocks_by_arm)}")
        object.__setattr__(self, "refinement_arm", arm)
        if trials_by_arm.get("full") != self.trials_per_cell:
            raise ValueError(
                f"trials_per_cell={self.trials_per_cell} is not the full arm's "
                f"own count {trials_by_arm.get('full')}; the scalar is the "
                f"paired-arm count and cannot disagree with the table")
        object.__setattr__(self, "shadow_fraction",
                           _require_fraction(self.shadow_fraction,
                                             "shadow_fraction"))
        for field in ("discrete_staircase_fallback", "stop_rule",
                      "failure_rule"):
            object.__setattr__(self, field,
                               _require_text(getattr(self, field), field,
                                             maximum=400))
        if not self.grid_half_width > self.couplings[-1]:
            raise ValueError(
                f"grid_half_width={self.grid_half_width!r} does not exceed the "
                f"strongest coupling {self.couplings[-1]!r}; the tongue would "
                f"reach past the edge of the population and the eligible count "
                f"would stop measuring the tongue — the same rule "
                f"raw_config.RawEventConfig enforces")

    @property
    def arms(self) -> tuple:
        """The arm names, in declared order."""
        return tuple(name for name, _ in self.arm_clocks)

    def clocks_of(self, arm: str) -> int:
        """How many clocks one named arm runs."""
        arm = _require_text(arm, "arm", maximum=40)
        for name, clocks in self.arm_clocks:
            if name == arm:
                return clocks
        raise KeyError(f"{arm!r} is not one of {list(self.arms)}")

    def trials_of(self, arm: str) -> int:
        """How many trials per cell one named arm runs."""
        arm = _require_text(arm, "arm", maximum=40)
        for name, trials in self.arm_trials:
            if name == arm:
                return trials
        raise KeyError(f"{arm!r} is not one of {list(self.arms)}")

    @property
    def grid_spacing(self) -> float:
        """``2 * half_width / clocks`` — the midpoint grid's spacing."""
        return 2.0 * self.grid_half_width / self.clocks

    @property
    def grid_parity(self) -> str:
        """``even`` or ``odd``.  Load-bearing at low coupling, so it is named."""
        return "even" if self.clocks % 2 == 0 else "odd"

    @property
    def eligible_at_weakest(self) -> int:
        """How many clocks are inside the tongue at the weakest coupling."""
        detunings = (np.arange(self.clocks) + 0.5) * self.grid_spacing \
            - self.grid_half_width
        return int(np.count_nonzero(np.abs(detunings) < self.couplings[0]))

    @property
    def meets_minimum_cells(self) -> bool:
        """Whether the weakest coupling admits enough clocks to be a population."""
        return self.eligible_at_weakest >= self.minimum_eligible_clocks

    def steps_at(self, timestep) -> int:
        """Finest steps the run window holds at one timestep."""
        timestep = _require_positive(timestep, "timestep")
        return max(1, int(math.ceil(self.pulse_duration / timestep)))

    @property
    def digest(self) -> str:
        return _digest({
            "label": self.label, "couplings": list(self.couplings),
            "clocks": self.clocks,
            "arm_clocks": [list(pair) for pair in self.arm_clocks],
            "arm_trials": [list(pair) for pair in self.arm_trials],
            "grid_half_width": self.grid_half_width,
            "pulse_duration": self.pulse_duration,
            "pulse_centre": self.pulse_centre,
            "phase_diffusion": self.phase_diffusion,
            "lock_tolerance": self.lock_tolerance,
            "dwell_time": self.dwell_time, "timestep": self.timestep,
            "timestep_ladder": list(self.timestep_ladder),
            "refinement_couplings": list(self.refinement_couplings),
            "refinement_arm": self.refinement_arm,
            "refinement_trials_per_cell": self.refinement_trials_per_cell,
            "trials_per_cell": self.trials_per_cell,
            "shadow_fraction": self.shadow_fraction,
            "minimum_eligible_clocks": self.minimum_eligible_clocks,
            "discrete_staircase_fallback": self.discrete_staircase_fallback,
            "stop_rule": self.stop_rule, "failure_rule": self.failure_rule,
        })


@dataclass(frozen=True)
class CostLine:
    """One named, separately costed piece of work.

    An independent review found the resource ledger reporting one number that
    silently contained a mispriced refinement arm, a shadow multiplier applied
    to a rate that already carried a shadow, and no line at all for the pilot
    or the numerical validation that would have to precede it.  A total is only
    honest when the things inside it are visible, so nothing is summed here
    that is not first named.
    """

    name: str
    runs: int
    physical_trials: int
    shadow_rewalks: int
    clock_steps: float
    durable_rows: float
    storage_bytes: float
    wall_seconds: float

    def __post_init__(self):
        object.__setattr__(self, "name",
                           _require_text(self.name, "name", maximum=120))
        for field in ("runs", "physical_trials"):
            object.__setattr__(self, field,
                               require_positive_int(getattr(self, field),
                                                    field, minimum=1))
        object.__setattr__(self, "shadow_rewalks",
                           require_positive_int(self.shadow_rewalks,
                                                "shadow_rewalks", minimum=0))
        for field in ("clock_steps", "durable_rows", "storage_bytes",
                      "wall_seconds"):
            object.__setattr__(self, field,
                               _require_positive(getattr(self, field), field))

    @property
    def digest(self) -> str:
        return _digest({"name": self.name, "runs": self.runs,
                        "physical_trials": self.physical_trials,
                        "shadow_rewalks": self.shadow_rewalks,
                        "clock_steps": self.clock_steps,
                        "durable_rows": self.durable_rows,
                        "storage_bytes": self.storage_bytes,
                        "wall_seconds": self.wall_seconds})


@dataclass(frozen=True)
class ResourceEstimate:
    """What one matrix would cost, scaled from a measurement, line by line.

    ``lines``
        The named :class:`CostLine` items.  The totals below are their sums and
        nothing else, so a reader can check the addition.
    ``safety_factor``
        Multiplies the wall-time estimate.  Not a fudge: the benchmark is a
        bounded configuration on one quiet machine.  Recorded so it can be
        divided out.
    """

    label: str
    matrix_digest: str
    benchmark_digest: str
    lines: tuple
    safety_factor: float
    clock_steps_per_second: float
    bytes_per_row: float
    fixed_bytes_per_run: int
    peak_bytes: int

    def __post_init__(self):
        object.__setattr__(self, "label",
                           _require_text(self.label, "label", maximum=120))
        for field in ("matrix_digest", "benchmark_digest"):
            object.__setattr__(self, field,
                               _require_hex_digest(getattr(self, field),
                                                   field))
        lines = tuple(require_exact_type(item, CostLine, f"lines[{position}]")
                      for position, item in enumerate(tuple(self.lines)))
        if not lines:
            raise ValueError("a resource estimate with no line items is a "
                             "total nobody can check")
        names = [item.name for item in lines]
        if len(set(names)) != len(names):
            raise ValueError(f"cost lines repeat a name: {names}")
        object.__setattr__(self, "lines", lines)
        for field in ("safety_factor", "clock_steps_per_second",
                      "bytes_per_row"):
            object.__setattr__(self, field,
                               _require_positive(getattr(self, field), field))
        for field in ("fixed_bytes_per_run", "peak_bytes"):
            object.__setattr__(self, field,
                               require_positive_int(getattr(self, field),
                                                    field, minimum=1))

    def _total(self, field):
        return sum(getattr(item, field) for item in self.lines)

    @property
    def clock_steps(self) -> float:
        """Walked clock-steps, summed over the line items."""
        return self._total("clock_steps")

    @property
    def wall_seconds(self) -> float:
        return self._total("wall_seconds")

    @property
    def wall_hours(self) -> float:
        """The wall-time estimate in hours, which is how a person reads it."""
        return self.wall_seconds / 3600.0

    @property
    def storage_bytes(self) -> float:
        return self._total("storage_bytes")

    @property
    def rows(self) -> float:
        return self._total("durable_rows")

    @property
    def physical_trials(self) -> int:
        return int(self._total("physical_trials"))

    @property
    def shadow_rewalks(self) -> int:
        return int(self._total("shadow_rewalks"))

    @property
    def runs(self) -> int:
        return int(self._total("runs"))

    def line(self, name: str):
        """One named line item."""
        name = _require_text(name, "name", maximum=120)
        for item in self.lines:
            if item.name == name:
                return item
        raise KeyError(f"{name!r} is not one of "
                       f"{[item.name for item in self.lines]}")

    @property
    def digest(self) -> str:
        return _digest({
            "label": self.label, "matrix_digest": self.matrix_digest,
            "benchmark_digest": self.benchmark_digest,
            "lines": [item.digest for item in self.lines],
            "safety_factor": self.safety_factor,
            "clock_steps_per_second": self.clock_steps_per_second,
            "bytes_per_row": self.bytes_per_row,
            "fixed_bytes_per_run": self.fixed_bytes_per_run,
            "peak_bytes": self.peak_bytes,
        })


def _cost_line(name, runs, trials, clocks, steps, shadow_fraction, rate,
               bytes_per_row, fixed_bytes, safety_factor):
    """Price one family of identically shaped runs, exactly.

    One run of ``trials`` trials over ``clocks`` clocks and ``steps`` steps:

    * physical trial executions: ``trials``;
    * shadow rewalks: ``ceil(trials * shadow_fraction)`` — **per run**, rounded
      up per run, because that is what ``raw_runner`` writes and because
      rounding once on a grand total is a different, smaller number;
    * walked clock-steps: ``RACE_PASSES * trials * clocks * steps`` for the two
      censoring passes, plus ``shadow * clocks * steps`` for the rewalk.  The
      shadow is charged here, once, against a rate whose own benchmark charged
      its shadow the same way — never as a multiplier on top of it;
    * durable rows: ``trials`` ledger rows, ``trials * clocks`` authoritative
      per-clock rows, and ``shadow * clocks`` shadow rows; and
    * storage: those rows at the measured row cost, plus the measured fixed
      per-run overhead (manifest, marker, three headers).
    """
    shadow = int(math.ceil(trials * shadow_fraction))
    clock_steps = runs * (RACE_PASSES * trials * clocks * steps
                          + shadow * clocks * steps)
    rows = runs * (trials + trials * clocks + shadow * clocks)
    return CostLine(name=name, runs=runs, physical_trials=runs * trials,
                    shadow_rewalks=runs * shadow, clock_steps=float(clock_steps),
                    durable_rows=float(rows),
                    storage_bytes=rows * bytes_per_row + runs * fixed_bytes,
                    wall_seconds=clock_steps / rate * safety_factor)


def matrix_resource_estimate(matrix, benchmark, safety_factor=1.0,
                             extra_lines=(), bytes_per_row=None,
                             fixed_bytes=None, peak_bytes=None
                             ) -> ResourceEstimate:
    """Scale a measured benchmark to a proposed matrix, arm by arm.

    One line per arm at the primary timestep, one line per refinement level,
    and whatever ``extra_lines`` the caller adds — the range pilot, an
    intended-configuration validation ladder — so that a total never silently
    omits work that has to happen.

    Units are normalized once, explicitly: the benchmark's rate is per
    **walked** clock-step and every line computes its own walked clock-steps
    the same way, including its own shadow.  An independent review found a
    shadow-bearing rate being multiplied by a shadow factor again, and this is
    the shape that cannot do that.
    """
    require_exact_type(matrix, PrimaryMatrix, "matrix")
    require_exact_type(benchmark, BenchmarkResult, "benchmark")
    safety_factor = _require_positive(safety_factor, "safety_factor")
    extra = tuple(require_exact_type(item, CostLine, f"extra_lines[{position}]")
                  for position, item in enumerate(tuple(extra_lines)))

    rate = benchmark.clock_steps_per_second
    # The *timing* fixture is the slowest one; its serialization is not
    # necessarily the widest.  An independent review found the storage estimate
    # taking 91.16 bytes/row from the slow-timing fixture while another fixture
    # in the same family measured 106.22, and the handoff summary quoting the
    # larger number beside a total built on the smaller one.  These three
    # overrides exist so a caller can hand in the **suite maxima** and have the
    # storage and memory figures be conservative in their own right rather than
    # inheriting whichever fixture happened to be slowest.
    per_row = (benchmark.bytes_per_row if bytes_per_row is None
               else _require_positive(bytes_per_row, "bytes_per_row"))
    fixed = (benchmark.fixed_bytes if fixed_bytes is None
             else require_positive_int(fixed_bytes, "fixed_bytes", minimum=1))
    peak = (benchmark.peak_bytes if peak_bytes is None
            else require_positive_int(peak_bytes, "peak_bytes", minimum=1))
    if per_row < benchmark.bytes_per_row:
        raise ValueError(
            f"bytes_per_row={per_row!r} is below the benchmark's own measured "
            f"{benchmark.bytes_per_row!r}; this override exists to be "
            f"conservative and cannot be used to make a run look smaller")
    primary_steps = matrix.steps_at(matrix.timestep)

    lines = []
    for name in matrix.arms:
        lines.append(_cost_line(
            f"primary/{name}", len(matrix.couplings), matrix.trials_of(name),
            matrix.clocks_of(name), primary_steps, matrix.shadow_fraction,
            rate, per_row, fixed, safety_factor))
    for level, timestep in enumerate(matrix.timestep_ladder[1:], start=1):
        lines.append(_cost_line(
            f"refinement/level{level}/{matrix.refinement_arm}",
            len(matrix.refinement_couplings),
            matrix.refinement_trials_per_cell,
            matrix.clocks_of(matrix.refinement_arm), matrix.steps_at(timestep),
            matrix.shadow_fraction, rate, per_row, fixed, safety_factor))
    lines.extend(extra)

    return ResourceEstimate(
        label=f"resource estimate for {matrix.label}",
        matrix_digest=matrix.digest, benchmark_digest=benchmark.digest,
        lines=tuple(lines), safety_factor=safety_factor,
        clock_steps_per_second=rate, bytes_per_row=per_row,
        fixed_bytes_per_run=fixed, peak_bytes=peak)


def pilot_cost_line(matrix, benchmark, trials_per_cell, cells,
                    safety_factor=1.0, bytes_per_row=None) -> CostLine:
    """Price the range-only pilot as its own line item.

    The pilot runs the **full** arm's grid over ``cells`` candidate couplings
    at ``trials_per_cell``.  It is separate work, it is not part of the
    production estimate, and its trials are permanently ineligible for
    production analysis — so it is a line a reader can see and subtract, never
    a term inside the matrix's own total.
    """
    require_exact_type(matrix, PrimaryMatrix, "matrix")
    require_exact_type(benchmark, BenchmarkResult, "benchmark")
    return _cost_line(
        "pilot/range-only", require_positive_int(cells, "cells", minimum=2),
        require_positive_int(trials_per_cell, "trials_per_cell", minimum=1),
        matrix.clocks_of("full"), matrix.steps_at(matrix.timestep),
        matrix.shadow_fraction, benchmark.clock_steps_per_second,
        max(benchmark.bytes_per_row,
            _require_positive(bytes_per_row, "bytes_per_row"))
        if bytes_per_row is not None else benchmark.bytes_per_row,
        benchmark.fixed_bytes,
        _require_positive(safety_factor, "safety_factor"))


# ---------------------------------------------------------------------------
# The validation kernels, the staged campaign, and the row-width bound
# ---------------------------------------------------------------------------
#: ``KernelBenchmark``, ``benchmark_kernel``, ``ValidationStage`` and
#: ``ValidationCampaign`` are deliberately **not** in ``__all__``.  They are
#: planning records for work that has not been authorized: they describe what a
#: validation campaign would cost and what its stop rules would be, they are
#: built and exercised by ``verify.py``, and nothing in the experiment's own
#: data path consumes them.  Keeping them off the exported surface says that a
#: campaign is a proposal rather than an API, and it keeps the public surface
#: to the records a decision is actually taken over.
#:
#: The worst-case identifiers a proposed run could serialize.  Used to build a
#: *schema-derived* row-width bound rather than trusting the small fixtures.
#: The **widest finite** spellings the encoder produces: 24 characters signed,
#: 23 unsigned.  An independent review found the earlier probe values
#: (``-1.2345678901234567e-16`` at 23, ``1234.5678901234567`` at 18) were not
#: the widest the same schema accepts, so the "bound" was not one.
_WIDTH_FLOAT = -1.7976931348623157e+308
_WIDTH_TIME = 1.7976931348623157e+308


@dataclass(frozen=True)
class KernelBenchmark:
    """A timing of one **actual** ticket-04 validation kernel.

    Not the raw writer.  An independent review found the validation options
    priced with raw-writer-equivalent work, and the killed-diffusion solve and
    the moving-band replay/bootstrap have different work and storage shapes
    entirely.  These records time the real public entry points.

    ``unit`` names what ``operations`` counts, because the two kernels do not
    share a unit: the oracle is priced per **space-time cell** of its PDE grid,
    the audit per **replayed interval** and its bootstrap per **resample**.
    Nothing in this module divides one by the other.
    """

    label: str
    kernel: str
    unit: str
    operations: int
    repeat_seconds: tuple
    peak_bytes: int
    warmup: int
    problem: str
    source_digest: str

    def __post_init__(self):
        for field in ("label", "kernel", "unit", "problem"):
            object.__setattr__(self, field,
                               _require_text(getattr(self, field), field,
                                             maximum=200))
        object.__setattr__(self, "source_digest",
                           _require_hex_digest(self.source_digest,
                                               "source_digest"))
        for field in ("operations", "peak_bytes"):
            object.__setattr__(self, field,
                               require_positive_int(getattr(self, field),
                                                    field, minimum=1))
        object.__setattr__(self, "warmup",
                           require_positive_int(self.warmup, "warmup",
                                                minimum=1))
        repeats = tuple(_require_positive(value, f"repeat_seconds[{position}]")
                        for position, value in
                        enumerate(tuple(self.repeat_seconds)))
        if len(repeats) < 2:
            raise ValueError(f"repeat_seconds holds {len(repeats)} timed "
                             f"repeat(s); one cannot establish machine "
                             f"variation")
        object.__setattr__(self, "repeat_seconds", repeats)

    @property
    def wall_seconds(self) -> float:
        """The slowest timed repeat."""
        return max(self.repeat_seconds)

    @property
    def operations_per_second(self) -> float:
        """Operations of this kernel's **own unit** per second, at the slowest."""
        return self.operations / self.wall_seconds

    @property
    def digest(self) -> str:
        return _digest({"label": self.label, "kernel": self.kernel,
                        "unit": self.unit, "operations": self.operations,
                        "repeat_seconds": list(self.repeat_seconds),
                        "peak_bytes": self.peak_bytes, "warmup": self.warmup,
                        "problem": self.problem,
                        "source_digest": self.source_digest})


def benchmark_kernel(label, kernel, unit, operations, run, problem,
                     fingerprint, warmup=1, repeats=2) -> KernelBenchmark:
    """Time ``run`` under ``tracemalloc``, discarding ``warmup`` calls first.

    ``run`` is a zero-argument callable that performs one whole invocation of
    the kernel being priced.  ``operations`` is how many of ``unit`` that
    invocation performs, supplied by the caller because only the caller knows
    the kernel's own shape.
    """
    require_exact_type(fingerprint, SourceFingerprint, "fingerprint")
    warmup = require_positive_int(warmup, "warmup", minimum=1)
    repeats = require_positive_int(repeats, "repeats", minimum=2)
    if not callable(run):
        raise TypeError(f"run must be callable, got {type(run).__name__}")
    timings = []
    peak = 0
    for index in range(warmup + repeats):
        timed = index >= warmup
        if timed:
            tracemalloc.start()
        start = time.perf_counter()
        run()
        elapsed = time.perf_counter() - start
        if timed:
            peak = max(peak, tracemalloc.get_traced_memory()[1])
            tracemalloc.stop()
            timings.append(elapsed)
    return KernelBenchmark(label=label, kernel=kernel, unit=unit,
                           operations=operations,
                           repeat_seconds=tuple(timings),
                           peak_bytes=int(peak), warmup=warmup,
                           problem=problem, source_digest=fingerprint.digest)


@dataclass(frozen=True)
class ValidationStage:
    """One finite, predeclared stage of the numerical-validation campaign.

    ``unit``
        Which budget unit this stage attacks — ``probability`` or ``time``.
        The two are attacked separately because a time-unit failure is
        trial-count independent and a probability failure is not.
    ``attacks``
        ``timestep_bias`` or ``sampling_error``.  Refining the step moves the
        bias; enlarging the sample moves the standard error; and a stage that
        claimed to do both at once would be the conflation an independent
        review found in the earlier projection.
    ``candidate``
        A short description of the resolution or sample size the stage would
        use.  Qualitative on purpose: an independent review showed that the
        only kernels timed so far omit the moving replay and the stationary
        endpoint generation entirely, so **there is no numeric price path in
        this record at all** — no seconds, no cap, no aggregation.  A number
        that is not a bound is worse than no number.

    ``success_rule``, ``no_result_rule``, ``stop_rule``
        Written down **before** the stage runs.  ``no_result_rule`` always
        preserves ``numerical_no_result``: a stage that misses its target is a
        numerical no-result, never a reason to widen the budget.
    ``sufficiency``
        ``not_promised`` for every stage in this campaign.  There is no field
        that can say a stage is guaranteed to resolve the block, because
        nothing has run.
    """

    label: str
    unit: str
    attacks: str
    kernel: str
    operations: int
    timestep_factor: float
    sample_factor: float
    depends_on: str
    candidate: str
    success_rule: str
    no_result_rule: str
    stop_rule: str
    sufficiency: str = "not_promised"

    ATTACKS = ("timestep_bias", "sampling_error")
    SUFFICIENCY = ("not_promised",)

    def __post_init__(self):
        for field in ("label", "kernel", "candidate", "success_rule",
                      "no_result_rule", "stop_rule"):
            object.__setattr__(self, field,
                               _require_text(getattr(self, field), field,
                                             maximum=400))
        object.__setattr__(self, "unit",
                           _require_member(self.unit, _oracle.UNITS, "unit"))
        object.__setattr__(self, "attacks",
                           _require_member(self.attacks,
                                           ValidationStage.ATTACKS, "attacks"))
        object.__setattr__(self, "sufficiency",
                           _require_member(self.sufficiency,
                                           ValidationStage.SUFFICIENCY,
                                           "sufficiency"))
        object.__setattr__(self, "operations",
                           require_positive_int(self.operations, "operations",
                                                minimum=1))
        for field in ("timestep_factor", "sample_factor"):
            object.__setattr__(self, field,
                               _require_positive(getattr(self, field), field))
        if self.timestep_factor > 1.0:
            raise ValueError("a validation stage does not propose a coarser "
                             "step than the evidence was measured at")
        if self.sample_factor < 1.0:
            raise ValueError("a validation stage does not propose a smaller "
                             "sample than the evidence was measured with")
        if self.attacks == "timestep_bias" and self.timestep_factor == 1.0:
            raise ValueError("a stage that attacks the timestep bias must "
                             "refine the timestep")
        if self.attacks == "sampling_error" and self.sample_factor == 1.0:
            raise ValueError("a stage that attacks the sampling error must "
                             "enlarge the sample")
        if "numerical_no_result" not in self.no_result_rule:
            raise ValueError(
                "no_result_rule must name numerical_no_result literally; a "
                "stage that misses its target is a numerical no-result, and "
                "the plan does not permit that to be re-decided")
        object.__setattr__(self, "depends_on",
                           _require_text(self.depends_on, "depends_on",
                                         maximum=200)
                           if self.depends_on else "")

    @property
    def digest(self) -> str:
        return _digest({
            "label": self.label, "unit": self.unit, "attacks": self.attacks,
            "kernel": self.kernel, "operations": self.operations,
            "timestep_factor": self.timestep_factor,
            "sample_factor": self.sample_factor,
            "depends_on": self.depends_on, "success_rule": self.success_rule,
            "no_result_rule": self.no_result_rule, "stop_rule": self.stop_rule,
            "candidate": self.candidate,
            "sufficiency": self.sufficiency})


@dataclass(frozen=True)
class ValidationCampaign:
    """A finite staged plan for gathering the evidence the budget lacks.

    It promises nothing, and it **carries no cost at all**.  ``sufficiency``
    is ``not_promised`` on every stage, :attr:`is_priced` is ``False``
    unconditionally, and there is no code path in this module that turns a
    stage's relative operation shape into seconds.  The campaign's own
    docstring-level claim is that running it *might* resolve the block — which
    is the only claim the present evidence supports.
    """

    label: str
    stages: tuple
    benchmarks: tuple
    scaling_assumption: str
    note: str

    def __post_init__(self):
        for field in ("label", "scaling_assumption", "note"):
            object.__setattr__(self, field,
                               _require_text(getattr(self, field), field,
                                             maximum=800))
        stages = tuple(require_exact_type(item, ValidationStage,
                                          f"stages[{position}]")
                       for position, item in enumerate(tuple(self.stages)))
        if not stages:
            raise ValueError("a campaign with no stages plans nothing")
        names = [item.label for item in stages]
        if len(set(names)) != len(names):
            raise ValueError(f"stages repeat a label: {names}")
        for item in stages:
            if item.depends_on and item.depends_on not in names:
                raise ValueError(f"stage {item.label!r} depends on "
                                 f"{item.depends_on!r}, which is not a stage")
        for unit in ("probability", "time"):
            if not any(item.unit == unit for item in stages):
                raise ValueError(
                    f"no stage attacks the {unit} unit; both fail in the "
                    f"present disposition and a campaign that addresses one "
                    f"cannot resolve the block")
        for attack in ValidationStage.ATTACKS:
            if not any(item.attacks == attack for item in stages):
                raise ValueError(f"no stage attacks {attack}")
        object.__setattr__(self, "stages", stages)
        benchmarks = tuple(
            require_exact_type(item, KernelBenchmark,
                               f"benchmarks[{position}]")
            for position, item in enumerate(tuple(self.benchmarks)))
        if not benchmarks:
            raise ValueError(
                "a campaign costed from no kernel measurement is the "
                "raw-writer-equivalent proxy this record exists to replace")
        object.__setattr__(self, "benchmarks", benchmarks)

    @property
    def is_priced(self) -> bool:
        """Always ``False``.  There is no code path that prices this."""
        return False

    @property
    def digest(self) -> str:
        return _digest({"label": self.label,
                        "stages": [item.digest for item in self.stages],
                        "benchmarks": [item.digest for item in self.benchmarks],
                        "scaling_assumption": self.scaling_assumption,
                        "note": self.note})


def row_width_bound(trials, clocks, steps) -> dict:
    """A **schema-derived** worst-case byte cost per durable row.

    **An intended-configuration bound, not a schema-wide one.**  There is no
    schema-wide worst row: ``raw_ledger``'s integer normalization has a lower
    bound and no upper bound, so a "widest possible" integer does not exist.
    What does exist is a bound *for the configuration the matrix proposes*, and
    every counter here is derived from it — trial and clock identifiers from
    ``trials`` and ``clocks``, and every endpoint, entry and reset counter from
    ``steps``, because a clock cannot record more endpoints than the run window
    has and cannot enter or reset more often than it has endpoints.  Changing
    the trial count, the clock count or the step count therefore forces a
    recomputation rather than silently reusing a stale coefficient.

    An independent review found the storage total using the row coefficient of
    whichever small fixture happened to be slowest.  Even the suite maximum is
    a measurement of *small* runs: a production ledger carries four-digit trial
    identifiers, two-digit clock identifiers, full seventeen-significant-digit
    float spellings and list cells as long as the co-completion list can get.

    So a worst-case valid schema-v3 row is constructed for each of the three
    tables at the proposed maxima and encoded through the writer's own
    ``raw_ledger.encode_row``.  The returned mapping gives the encoded width of
    each table's worst row, the three header lines, and the widest row overall
    — which a caller takes the maximum of against the empirical figure.

    **What is included**: the trial ledger, the authoritative per-``(trial,
    clock)`` table and the no-stop shadow sample, which are the three tables a
    closed run writes.  **What is excluded**: the manifest and the close
    marker, which are per-run fixed overhead and are returned separately; and
    any derived product — binned survival tables, plots, analysis output —
    which no run writes and which this ticket does not cost.
    """
    trials = require_positive_int(trials, "trials", minimum=1)
    clocks = require_positive_int(clocks, "clocks", minimum=1)
    steps = require_positive_int(steps, "steps", minimum=1)
    worst_trial = trials - 1
    worst_clock = clocks - 1
    winners = tuple(range(clocks))
    # Every counter's own intended-configuration maximum.  A clock records at
    # most one endpoint per step, and cannot enter a band or reset a dwell more
    # often than it has endpoints; the channel's totals are those summed over
    # the population.
    per_clock = steps
    per_channel = steps * clocks

    def widest(columns, values):
        row = dict(zip((name for name, _ in columns), values))
        return len(_raw_ledger.encode_row(columns, row).encode("utf-8")) + 1

    ledger = widest(_raw_ledger.LEDGER_COLUMNS, (
        worst_trial, "never_eligible", _WIDTH_TIME, _WIDTH_TIME, True, clocks,
        winners, (_WIDTH_FLOAT,) * clocks, (_WIDTH_FLOAT,) * clocks,
        (_WIDTH_TIME,) * clocks, (_WIDTH_TIME,) * clocks,
        (_WIDTH_TIME,) * clocks, clocks, clocks, per_channel,
        per_channel, _WIDTH_TIME))
    clock_row = widest(_raw_ledger.CLOCK_COLUMNS, (
        worst_trial, worst_clock, _WIDTH_FLOAT, "never_eligible", _WIDTH_TIME,
        _WIDTH_TIME, _WIDTH_TIME, _WIDTH_TIME, per_clock, per_clock,
        per_clock, per_clock, per_clock, per_clock, _WIDTH_TIME, _WIDTH_TIME,
        _WIDTH_FLOAT, True, True))
    shadow = widest(_raw_ledger.SHADOW_COLUMNS, (
        worst_trial, worst_clock, _WIDTH_FLOAT, "never_eligible", _WIDTH_TIME,
        _WIDTH_TIME, _WIDTH_TIME, _WIDTH_TIME, per_clock, per_clock,
        per_clock, per_clock, per_clock, per_clock, _WIDTH_TIME, _WIDTH_TIME,
        _WIDTH_FLOAT))
    headers = sum(
        len(_raw_ledger.header_row(columns).encode("utf-8")) + 1
        for columns in (_raw_ledger.LEDGER_COLUMNS, _raw_ledger.CLOCK_COLUMNS,
                        _raw_ledger.SHADOW_COLUMNS))
    # The *mix* is what a run actually writes: one ledger row and ``clocks``
    # authoritative rows per trial.  Quoting the widest single row would price
    # every row as a 64-winner ledger line, which no clocks row ever is; the
    # blended figure is the conservative bound for the mix the matrix runs.
    blended = (ledger + clocks * clock_row) / (1.0 + clocks)
    central = (ledger + clock_row) / 2.0
    return {"ledger": float(ledger), "clocks": float(clock_row),
            "shadow": float(shadow),
            "widest": float(max(ledger, clock_row, shadow)),
            "blended": float(blended), "central_blended": float(central),
            "header_bytes": float(headers)}


# ---------------------------------------------------------------------------
# The pilot firewall
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class PilotFirewall:
    """The wall between a range-only pilot and everything downstream.

    Five separations, and the fifth is the one an independent review found
    missing:

    1. **Streams.**  ``pilot_namespace`` and ``production_namespace`` are
       distinct and carry disjoint prefixes, so no counter key can be shared.
    2. **Folders.**  ``pilot_prefix`` and ``production_prefix`` name disjoint
       run-directory families, and :meth:`require_pilot_run` now checks the
       **run name** against them as well as the manifest.  A pilot-labelled,
       pilot-keyed run sitting in a production-named directory used to be
       accepted; it is not a naming quibble, because the directory family is
       what a later reader uses to decide what a run *is*.
    3. **Quarantine.**  A pilot run is **written directly** beneath
       :data:`QUARANTINE_DIRNAME`, which ``raw_runner.open_raw_run`` cannot
       address, so the unrestricted ledger reader is not an alternate route to
       forbidden pilot data.  A typed narrow return never prevented that; only
       the files' location does, and there is no write-then-move window in
       which the production root holds a pilot run.
    4. **Eligibility.**  ``pilot_trials_eligible`` is ``False`` and cannot be
       set otherwise.
    5. **Information.**  :meth:`permits` answers, for one field name, whether a
       pilot may expose it — by allowlist, so a field nobody thought about is
       forbidden by default.
    """

    pilot_namespace: str
    production_namespace: str
    pilot_prefix: str
    production_prefix: str
    minimum_events: int
    minimum_survivors: int
    minimum_cells: int
    pilot_trials_eligible: bool = False

    def __post_init__(self):
        for field, prefix in (("pilot_namespace", PILOT_NAMESPACE_PREFIX),
                              ("production_namespace",
                               PRODUCTION_NAMESPACE_PREFIX),
                              ("pilot_prefix", PILOT_NAMESPACE_PREFIX),
                              ("production_prefix",
                               PRODUCTION_NAMESPACE_PREFIX)):
            value = _require_text(getattr(self, field), field, maximum=64)
            if not all(character.isalnum() or character in "-_"
                       for character in value):
                raise ValueError(f"{field} must be letters, digits, hyphen or "
                                 f"underscore, got {value!r}")
            if not value.startswith(prefix):
                raise ValueError(
                    f"{field}={value!r} must start with {prefix!r}; the pilot "
                    f"and production families are told apart by their "
                    f"prefixes, and a name that does not carry one can be "
                    f"filed on either side")
            object.__setattr__(self, field, value)
        if self.pilot_namespace == self.production_namespace:
            raise ValueError(
                f"the pilot and production stream namespaces are both "
                f"{self.pilot_namespace!r}; sharing one means the pilot and "
                f"the production run walk the same Brownian paths, and the "
                f"pilot's trials would be inside the production estimate "
                f"whatever any label said")
        if self.pilot_prefix == self.production_prefix:
            raise ValueError(f"the pilot and production run-directory prefixes "
                             f"are both {self.pilot_prefix!r}")
        if (self.pilot_prefix.startswith(self.production_prefix)
                or self.production_prefix.startswith(self.pilot_prefix)):
            raise ValueError(
                f"run-directory prefixes {self.pilot_prefix!r} and "
                f"{self.production_prefix!r} are nested; a name under one "
                f"would also be under the other")
        for field in ("minimum_events", "minimum_survivors", "minimum_cells"):
            object.__setattr__(self, field,
                               require_positive_int(getattr(self, field),
                                                    field, minimum=1))
        if self.pilot_trials_eligible is not False:
            raise ValueError(
                "pilot_trials_eligible must be False; a pilot's trials are "
                "permanently ineligible for production analysis, and the field "
                "is recorded so that the manifest says so rather than so that "
                "a caller can choose")

    def permits(self, field: str) -> bool:
        """May a pilot expose ``field``?  Allowlist, not denylist."""
        return _require_text(field, "field", maximum=80) in \
            PERMITTED_PILOT_FIELDS

    def require_pilot_name(self, run_name: str) -> str:
        """Return ``run_name`` if it belongs to the pilot directory family.

        One path segment, carrying the pilot prefix and **not** the production
        one.  An independent review wrote a pilot-labelled, pilot-keyed
        four-trial run into a directory called ``prod-t07-review-bypass`` and
        the count door accepted it, because nothing looked at the name.  The
        directory family is how a later reader decides what a run is, so it is
        checked here with the same force as the stream.
        """
        name = _require_text(run_name, "run_name", maximum=96)
        if "/" in name or "\\" in name or name in (".", ".."):
            raise ValueError(f"run_name must be one path segment, got "
                             f"{name!r}")
        if name.startswith(self.production_prefix):
            raise ValueError(
                f"run name {name!r} is in the production directory family "
                f"{self.production_prefix!r}; a run filed under production "
                f"cannot be read through the pilot door whatever its manifest "
                f"says, because the directory family is what a later reader "
                f"uses to decide what it is")
        if not name.startswith(self.pilot_prefix):
            raise ValueError(
                f"run name {name!r} does not carry the pilot directory prefix "
                f"{self.pilot_prefix!r}; a run in neither family is a third "
                f"dataset nobody declared")
        return name

    def require_pilot_run(self, manifest, run_name):
        """Return ``manifest`` if the run belongs to this firewall's pilot.

        Three things about the manifest and one about the name: the run must be
        *labelled* a pilot, it must be on the pilot stream namespace, it must
        not be on the production one, and its directory must be in the pilot
        family.  The third is not implied by the second — a namespace that is
        neither is a third dataset nobody declared — and the fourth is not
        implied by any of them.
        """
        require_exact_type(manifest, _raw_ledger.RawManifest, "manifest")
        self.require_pilot_name(run_name)
        label = manifest.value("run_label")
        if label != "pilot":
            raise ValueError(
                f"this run is labelled {label!r}; the count-only door reads "
                f"pilot runs, and reading a production run through it would "
                f"expose production counts before the manifest is signed")
        namespace = manifest.value("stream_namespace")
        if namespace == self.production_namespace:
            raise ValueError(
                f"this run is on the production stream namespace "
                f"{namespace!r} while calling itself a pilot; the label and "
                f"the stream disagree and the stream is what decides which "
                f"Brownian paths were walked")
        if namespace != self.pilot_namespace:
            raise ValueError(
                f"this run is on stream namespace {namespace!r}, which is "
                f"neither the declared pilot namespace "
                f"{self.pilot_namespace!r} nor anything this firewall knows")
        return manifest

    @property
    def digest(self) -> str:
        return _digest({
            "pilot_namespace": self.pilot_namespace,
            "production_namespace": self.production_namespace,
            "pilot_prefix": self.pilot_prefix,
            "production_prefix": self.production_prefix,
            "quarantine": QUARANTINE_DIRNAME,
            "minimum_events": self.minimum_events,
            "minimum_survivors": self.minimum_survivors,
            "minimum_cells": self.minimum_cells,
            "pilot_trials_eligible": self.pilot_trials_eligible,
            "permitted_fields": list(PERMITTED_PILOT_FIELDS),
        })


# ---------------------------------------------------------------------------
# The frozen pilot plan
# ---------------------------------------------------------------------------
#: The only role a pilot may be given.  A range-only pilot chooses a coupling
#: range from counts and does nothing else; there is no second member for a
#: caller to move to, so "which pilot was approved" cannot be answered with a
#: different kind of run.
PILOT_ROLES = ("range_only",)

#: What a pilot plan is **for**.  An independent review pointed out that the
#: verifier's deliberately cheap eight-clock fixture plan could reach the
#: production proposal and the checkpoint, which are the two places a plan
#: becomes a claim about the experiment John would authorize.  The purpose is a
#: declared, digest-bound field and the production doors accept exactly one
#: value, so a fixture plan is refused *before* any decision rather than
#: distinguished by a reader noticing its clock count.
PILOT_PURPOSES = ("approved_production_pilot", "verifier_fixture")

#: The reserved marker a **verifier fixture** plan's run-name family must
#: carry, and an approved production pilot's must not.
#:
#: The purpose would otherwise be a declaration a caller could simply restate:
#: reaching the construction token and relabelling a fixture
#: ``approved_production_pilot`` would have passed the rebuild, because the
#: rebuild would have been handed the same relabelled purpose.  Tying it to the
#: run-name family makes the two statements contradict each other -- a plan
#: claiming to be the approved pilot while writing its cells under
#: ``...-fixture-...`` names does not build -- so the separation survives a
#: caller who has the token.
FIXTURE_PREFIX_MARKER = "-fixture-"

#: **The** approved production pilot, frozen here as declared constants.
#:
#: An independent review pointed out that ``purpose`` plus a run-name marker
#: only separates a plan that *calls itself* a fixture: an eight-clock,
#: ``dt = 0.05``, four-cell plan declaring ``approved_production_pilot``
#: satisfied every self-consistency rule and drove the production proposal.
#: Self-consistency is not identity.  There is exactly one pilot ticket 07
#: would put to a checkpoint, and these are its numbers; every production door
#: requires field-for-field and digest equality with the plan
#: :func:`canonical_pilot_plan` builds from them, the trusted live source
#: fingerprint and the live schema versions.
PRODUCTION_PILOT_SPEC = {
    "label": "T07 range-only pilot, frozen before checkpoint 1",
    "couplings": (0.5, 0.6597539553864471, 0.8705505632961241,
                  1.1486983549970349, 1.515716566510398, 2.0),
    "clocks": 64,
    "grid_half_width": 3.0,
    "timestep": 0.001953125,
    "pulse_duration": 4.0,
    "pulse_centre": 0.0,
    "phase_diffusion": 0.08,
    "lock_tolerance": 0.35,
    "dwell_time": 0.5,
    "trials_per_cell": 200,
    "shadow_fraction": 0.05,
    "run_name_prefix": "xpilot-t07-cell",
    "model": "full",
}

#: The pilot plan's construction token.  A separate object from the proposal's
#: seal, so that reaching one does not reach the other, and — like that one — a
#: permission to *construct* and never an authority: a reachable-token plan
#: still has to agree with the manifest of every run it authorizes and with the
#: proposal that binds its digest.
_PLAN_SEAL = object()

#: The count-only selection rule, in words, frozen inside the plan digest.
#: Carried as text as well as as three numbers, because what a checkpoint
#: approves is the rule and not only its thresholds.
PILOT_SELECTION_RULE = (
    "a cell is eligible when its committed count is at least minimum_events "
    "and its survivor count is at least minimum_survivors; the selection is "
    "the longest run of consecutive eligible cells, ties broken to the "
    "weakest such run; and it must hold at least minimum_cells cells or the "
    "pilot refuses with a named reason.  Nothing in the rule reads an "
    "exponent, a curvature, a survival shape, a commit time, a dwell or reset "
    "count, a noise strength, a band, a pulse, a population, a state machine, "
    "a likelihood or an exclusion")


@dataclass(frozen=True)
class PilotPlan:
    """The exact pilot a checkpoint would approve, frozen before it runs.

    An independent review found that the pilot writer authenticated **where** a
    run was written and not **which pilot John approved**: it accepted an
    arbitrary ``RawEventConfig`` and either raw model, so a width-only,
    two-clock, ``dt=1``, two-trial diagnostic run could be written to the right
    place and then projected into authoritative range evidence.  The firewall
    checked namespaces, prefixes and thresholds; nothing checked the physics.

    This record is that missing half.  It names the role, the model, the run
    label, the ordered candidate couplings, the whole clock grid, the intended
    timestep and the window it implies, the pulse and noise and band and dwell,
    the trial and shadow counts, the execution blocking, the reserved stream
    namespace and run-name family, the count-only selection rule and its three
    thresholds, and the source, environment and schema versions the plan was
    frozen against.  :meth:`config_for` derives the one configuration each
    coupling may be run at, and :meth:`require_manifest` rebuilds the whole
    45-key raw manifest from that configuration and requires **field-for-field
    equality**, so a run whose physics is not this plan's physics is refused at
    the count door as well as at the writer.

    It is **factory-only**: build it with :func:`pilot_plan`, which derives
    every value that can be derived and checks every relationship between the
    ones that cannot.  The seal is a module attribute and nothing depends on
    its secrecy — a reachable-token plan still has to agree with the manifest
    of every run it is asked to authorize, and with the proposal that binds its
    digest.
    """

    label: str
    role: str
    model: str
    run_label: str
    couplings: tuple
    clocks: int
    grid_half_width: float
    grid_spacing: float
    grid_parity: str
    grid_origin: float
    grid_support: tuple
    timestep: float
    steps: int
    window_open: float
    window_close: float
    pulse_duration: float
    pulse_centre: float
    phase_diffusion: float
    lock_tolerance: float
    dwell_time: float
    trials_per_cell: int
    shadow_trials: int
    shadow_fraction: float
    clock_block: int
    step_window: int
    stream_namespace: str
    namespace_prefix: str
    run_name_prefix: str
    purpose: str
    minimum_events: int
    minimum_survivors: int
    minimum_cells: int
    selection_rule: str
    firewall_digest: str
    source_digest: str
    python_version: str
    numpy_version: str
    platform: str
    ledger_schema_version: int
    manifest_schema_version: int
    config_digest: str
    seal: object = None

    def __post_init__(self):
        if self.seal is not _PLAN_SEAL:
            raise TypeError(
                "PilotPlan is factory-only: build it with pilot_plan, which "
                "derives the grid, the window and the configuration digest "
                "from the declared inputs rather than taking them on trust.  A "
                "plan whose grid spacing disagrees with its own clock count "
                "would authorize a pilot nobody approved")
        object.__setattr__(self, "seal", None)

    @property
    def cells(self) -> int:
        """How many candidate couplings the pilot would run."""
        return len(self.couplings)

    @property
    def run_names(self) -> tuple:
        """The **exact ordered** run names this plan's cells are written under.

        Derived from the run-name family and the position of each declared
        coupling, so the name/coupling mapping is frozen inside the plan digest
        rather than left to a caller.  An independent review handed the
        authoritative doors the same four runs *in reverse*, and got a
        selection back; and omitted one, and got a scientific
        ``run_too_short`` refusal rather than an authority refusal.  Both are
        questions about which runs the approved pilot wrote, and the answer is
        here.
        """
        return tuple(f"{self.run_name_prefix}-{index}"
                     for index in range(len(self.couplings)))

    def coupling_for(self, run_name):
        """The one declared coupling this run name carries."""
        names = self.run_names
        name = _require_text(run_name, "run_name", maximum=96)
        if name not in names:
            raise ValueError(
                f"{name!r} is not one of this plan's declared run names "
                f"{list(names)}; a run the plan does not name is not a cell of "
                f"the approved pilot")
        return self.couplings[names.index(name)]

    def run_name_for(self, coupling):
        """The one run name this plan's cell at ``coupling`` is written under."""
        value = _require_positive(coupling, "coupling")
        for index, declared in enumerate(self.couplings):
            if value == declared:
                return self.run_names[index]
        raise ValueError(
            f"coupling {value!r} is not one of this plan's declared "
            f"candidates {list(self.couplings)}")

    def require_run_names(self, run_names) -> tuple:
        """Return ``run_names`` only if they are **exactly** this plan's, in order.

        Set, length, order, uniqueness and the one-to-one coupling mapping, all
        at once.  A missing, extra, duplicated, reordered or foreign run is an
        **authority** refusal — the caller did not name the approved pilot —
        and never a scientific no-result, because "these are not that pilot's
        runs" is not a statement about coupling ranges.
        """
        supplied = tuple(_require_text(value, f"run_names[{position}]",
                                       maximum=96)
                         for position, value in enumerate(tuple(run_names)))
        expected = self.run_names
        if supplied == expected:
            return supplied
        missing = [name for name in expected if name not in supplied]
        extra = [name for name in supplied if name not in expected]
        duplicated = sorted({name for name in supplied
                             if supplied.count(name) > 1})
        raise ValueError(
            f"these are not the approved pilot's runs.  The plan writes "
            f"{list(expected)}, in that order, one per declared coupling; the "
            f"names supplied are {list(supplied)}"
            + (f", missing {missing}" if missing else "")
            + (f", with the undeclared {extra}" if extra else "")
            + (f", repeating {duplicated}" if duplicated else "")
            + (", in a different order" if not missing and not extra
               and not duplicated else "")
            + ".  A range chosen from a different set of runs, or from these "
              "runs read in a different order, is not the range this pilot "
              "measured")

    def config_for(self, coupling):
        """The one raw configuration this plan permits at ``coupling``.

        Derived, never supplied.  Every field of the returned configuration
        comes from the plan; the only freedom is *which declared coupling*, and
        an undeclared one is refused rather than run.
        """
        value = _require_positive(coupling, "coupling")
        if not any(value == declared for declared in self.couplings):
            raise ValueError(
                f"coupling {value!r} is not one of this plan's declared "
                f"candidates {list(self.couplings)}; a pilot cell at a "
                f"coupling the plan does not name is a cell nobody approved")
        return _raw_config.RawEventConfig(
            grid=_raw_config.RawClockGrid(self.grid_half_width, self.clocks),
            peak_coupling=value, pulse_duration=self.pulse_duration,
            pulse_centre=self.pulse_centre,
            phase_diffusion=self.phase_diffusion,
            lock_tolerance=self.lock_tolerance, dwell_time=self.dwell_time,
            timestep=self.timestep, trials=self.trials_per_cell,
            stream_namespace=self.stream_namespace,
            run_label=self.run_label, shadow_trials=self.shadow_trials)

    def require_run_name(self, run_name: str) -> str:
        """Return ``run_name`` if it is one of this plan's **declared** names."""
        name = _require_text(run_name, "run_name", maximum=96)
        if "/" in name or "\\" in name or name in (".", ".."):
            raise ValueError(f"run_name must be one path segment, got "
                             f"{name!r}")
        if name not in self.run_names:
            raise ValueError(
                f"run name {name!r} is not one of this plan's declared run "
                f"names {list(self.run_names)}; the plan freezes which cell is "
                f"written under which name, so a run outside that mapping is a "
                f"run the plan does not describe")
        return name

    def require_manifest(self, manifest, coupling):
        """Return ``manifest`` only if it is **this plan's** manifest.

        The whole 45-key raw manifest is rebuilt from :meth:`config_for` and
        compared field for field, so every physical, numerical, grid, window,
        stream and schema value has to agree.  A run at the right place, on the
        right stream, with the right label and the wrong model, clock count,
        timestep, trial count or pulse is refused here — which is the finding
        this record exists to close.
        """
        require_exact_type(manifest, _raw_ledger.RawManifest, "manifest")
        expected = _raw_runner.raw_manifest(self.config_for(coupling),
                                            self.model)
        disagree = [key for key in _raw_ledger.MANIFEST_KEYS
                    if manifest.value(key) != expected.value(key)]
        if disagree:
            raise ValueError(
                f"this run's manifest disagrees with the approved pilot plan "
                f"on {disagree}; the plan names the exact model, grid, "
                f"timestep, window, pulse, noise, band, dwell, trial and "
                f"shadow counts a pilot may be run at, and a run that is not "
                f"that run cannot become range evidence for it")
        return manifest

    @property
    def digest(self) -> str:
        return _digest({
            "label": self.label, "role": self.role, "model": self.model,
            "run_label": self.run_label, "couplings": list(self.couplings),
            "clocks": self.clocks, "grid_half_width": self.grid_half_width,
            "grid_spacing": self.grid_spacing, "grid_parity": self.grid_parity,
            "grid_origin": self.grid_origin,
            "grid_support": list(self.grid_support),
            "timestep": self.timestep, "steps": self.steps,
            "window_open": self.window_open,
            "window_close": self.window_close,
            "pulse_duration": self.pulse_duration,
            "pulse_centre": self.pulse_centre,
            "phase_diffusion": self.phase_diffusion,
            "lock_tolerance": self.lock_tolerance,
            "dwell_time": self.dwell_time,
            "trials_per_cell": self.trials_per_cell,
            "shadow_trials": self.shadow_trials,
            "shadow_fraction": self.shadow_fraction,
            "clock_block": self.clock_block,
            "step_window": self.step_window,
            "stream_namespace": self.stream_namespace,
            "namespace_prefix": self.namespace_prefix,
            "run_name_prefix": self.run_name_prefix,
            "run_names": list(self.run_names),
            "purpose": self.purpose,
            "minimum_events": self.minimum_events,
            "minimum_survivors": self.minimum_survivors,
            "minimum_cells": self.minimum_cells,
            "selection_rule": self.selection_rule,
            "firewall_digest": self.firewall_digest,
            "source_digest": self.source_digest,
            "python_version": self.python_version,
            "numpy_version": self.numpy_version,
            "platform": self.platform,
            "ledger_schema_version": self.ledger_schema_version,
            "manifest_schema_version": self.manifest_schema_version,
            "config_digest": self.config_digest,
        })


def pilot_plan(label, firewall, fingerprint, couplings, clocks,
               grid_half_width, timestep, pulse_duration, pulse_centre,
               phase_diffusion, lock_tolerance, dwell_time, trials_per_cell,
               shadow_fraction, run_name_prefix, clock_block=None,
               step_window=None, model="full",
               purpose=PILOT_PURPOSES[0]) -> PilotPlan:
    """Freeze the exact pilot a checkpoint would approve.

    Everything derivable is derived: the grid spacing, parity and support come
    from ``raw_config.RawClockGrid``; the step count and the observation window
    come from ``raw_runner.raw_run_window``; the shadow count is the declared
    fraction of the trial count, rounded up **per run** exactly as the writer
    rounds it; and the configuration digest is taken over a real
    ``RawEventConfig`` built at the weakest declared coupling, so the plan
    cannot describe a configuration that does not construct.

    The plan is bound to a :class:`PilotFirewall` and to a
    :class:`SourceFingerprint` at construction, so the approved pilot is
    specific to a stream, a run-name family, a set of count thresholds, and a
    set of source bytes.
    """
    require_exact_type(firewall, PilotFirewall, "firewall")
    require_exact_type(fingerprint, SourceFingerprint, "fingerprint")
    label = _require_text(label, "label", maximum=120)
    purpose = _require_member(purpose, PILOT_PURPOSES, "purpose")
    model = _require_member(model, _raw_ledger.MODEL_LABELS, "model")
    if model != "full":
        raise ValueError(
            f"a range-only pilot runs the {'full'!r} arm; {model!r} is a "
            f"control arm, and choosing a range from a control's counts would "
            f"choose it from a different experiment")
    values = _require_couplings(couplings, "couplings")
    if len(values) < firewall.minimum_cells:
        raise ValueError(
            f"the plan declares {len(values)} candidate couplings and the "
            f"count rule needs at least {firewall.minimum_cells} consecutive "
            f"eligible cells; a pilot that cannot satisfy its own selection "
            f"rule is a pilot that refuses whatever it measures")
    prefix = _require_text(run_name_prefix, "run_name_prefix", maximum=64)
    if not prefix.startswith(firewall.pilot_prefix):
        raise ValueError(
            f"the run-name family {prefix!r} is not inside the firewall's "
            f"pilot directory family {firewall.pilot_prefix!r}")
    # The purpose and the run-name family have to agree, so that relabelling
    # one without the other does not build.
    marked = FIXTURE_PREFIX_MARKER in prefix
    if purpose == "verifier_fixture" and not marked:
        raise ValueError(
            f"a {purpose!r} plan writes its cells under a run-name family "
            f"carrying {FIXTURE_PREFIX_MARKER!r}; {prefix!r} does not, and a "
            f"fixture whose runs are indistinguishable by name from the "
            f"approved pilot's is not separated from it")
    if purpose == "approved_production_pilot" and marked:
        raise ValueError(
            f"the run-name family {prefix!r} carries the reserved fixture "
            f"marker {FIXTURE_PREFIX_MARKER!r} while the plan calls itself "
            f"{purpose!r}; the purpose and the names it writes under must "
            f"agree, or the purpose is a label a caller can simply restate")
    fraction = _require_fraction(shadow_fraction, "shadow_fraction")
    trials = require_positive_int(trials_per_cell, "trials_per_cell",
                                  minimum=1)
    shadow = max(1, int(math.ceil(trials * fraction)))
    grid = _raw_config.RawClockGrid(
        _require_positive(grid_half_width, "grid_half_width"),
        require_positive_int(clocks, "clocks", minimum=1))
    block = (_raw_runner.DEFAULT_CLOCK_BLOCK if clock_block is None
             else require_positive_int(clock_block, "clock_block", minimum=1))
    window = (_raw_runner.DEFAULT_STEP_WINDOW if step_window is None
              else require_positive_int(step_window, "step_window", minimum=1))
    config = _raw_config.RawEventConfig(
        grid=grid, peak_coupling=values[0], pulse_duration=pulse_duration,
        pulse_centre=pulse_centre, phase_diffusion=phase_diffusion,
        lock_tolerance=lock_tolerance, dwell_time=dwell_time,
        timestep=timestep, trials=trials,
        stream_namespace=firewall.pilot_namespace, run_label="pilot",
        shadow_trials=shadow)
    mesh, start_index, steps = _raw_runner.raw_run_window(config)
    return PilotPlan(
        label=label, role=PILOT_ROLES[0], model=model, run_label="pilot",
        couplings=values, clocks=grid.n_points,
        grid_half_width=grid.half_width, grid_spacing=grid.spacing,
        grid_parity=grid.parity, grid_origin=0.0,
        grid_support=(-grid.half_width, grid.half_width),
        timestep=config.timestep, steps=steps,
        window_open=mesh.time_at(start_index),
        window_close=mesh.time_at(start_index + steps),
        pulse_duration=config.pulse_duration,
        pulse_centre=config.pulse_centre,
        phase_diffusion=config.phase_diffusion,
        lock_tolerance=config.lock_tolerance, dwell_time=config.dwell_time,
        trials_per_cell=trials, shadow_trials=shadow, shadow_fraction=fraction,
        clock_block=block, step_window=window,
        stream_namespace=firewall.pilot_namespace,
        namespace_prefix=PILOT_NAMESPACE_PREFIX, run_name_prefix=prefix,
        purpose=purpose,
        minimum_events=firewall.minimum_events,
        minimum_survivors=firewall.minimum_survivors,
        minimum_cells=firewall.minimum_cells,
        selection_rule=PILOT_SELECTION_RULE,
        firewall_digest=firewall.digest, source_digest=fingerprint.digest,
        python_version=fingerprint.python_version,
        numpy_version=fingerprint.numpy_version, platform=fingerprint.platform,
        ledger_schema_version=_raw_ledger.LEDGER_SCHEMA_VERSION,
        manifest_schema_version=_raw_ledger.MANIFEST_SCHEMA_VERSION,
        config_digest=_digest({
            "grid": [grid.half_width, grid.n_points],
            "pulse": [config.pulse_duration, config.pulse_centre],
            "noise": config.phase_diffusion,
            "band": config.lock_tolerance, "dwell": config.dwell_time,
            "timestep": config.timestep, "trials": config.trials,
            "shadow_trials": config.shadow_trials,
            "stream_namespace": config.stream_namespace,
            "run_label": config.run_label, "model": model,
            "clock_block": block, "step_window": window}),
        seal=_PLAN_SEAL)


def canonical_pilot_plan(firewall, fingerprint) -> PilotPlan:
    """Build **the** approved production pilot, from frozen constants.

    :data:`PRODUCTION_PILOT_SPEC` plus the trusted live source fingerprint and
    the live schema versions.  There is one of these, it is what a checkpoint
    would approve, and :func:`require_canonical_pilot_plan` is what every
    production door compares against.
    """
    return pilot_plan(firewall=firewall, fingerprint=fingerprint,
                      purpose="approved_production_pilot",
                      **PRODUCTION_PILOT_SPEC)


def require_canonical_pilot_plan(plan, firewall, fingerprint) -> PilotPlan:
    """Return ``plan`` only if it **is** the approved production pilot.

    Field for field and digest for digest against
    :func:`canonical_pilot_plan`.  An independent review's eight-clock,
    ``dt = 0.05``, four-cell plan declared the production purpose, carried a
    run-name family with no fixture marker, rebuilt from its own inputs
    consistently, and drove the proposal — every rule it met was a rule about
    self-consistency, and none of them asked *which* pilot it was.

    Generic, fixture and probe plans keep their use: they exercise the count
    and power mechanics at the lower doors.  They reach no production decision.
    """
    require_authoritative_plan(plan, firewall, fingerprint,
                               purpose="approved_production_pilot")
    canonical = canonical_pilot_plan(firewall, fingerprint)
    disagree = [name for name in PilotPlan.__dataclass_fields__
                if name != "seal"
                and getattr(plan, name) != getattr(canonical, name)]
    if disagree or plan.digest != canonical.digest:
        raise ValueError(
            f"this is not the approved production pilot.  It disagrees with "
            f"the one PRODUCTION_PILOT_SPEC declares on "
            f"{disagree or ['digest']}: the approved pilot runs "
            f"{canonical.clocks} clocks at timestep {canonical.timestep!r} "
            f"over {canonical.cells} declared couplings at "
            f"{canonical.trials_per_cell} trials per cell, and a plan that is "
            f"internally consistent is not thereby that plan")
    return plan


def require_authoritative_plan(plan, firewall, fingerprint,
                               purpose=PILOT_PURPOSES[0]) -> PilotPlan:
    """Return ``plan`` only if the plan factory **rebuilds it exactly**.

    The seal is a construction token and nothing more.  An independent review
    reached ``_PLAN_SEAL``, built an exact-typed :class:`PilotPlan` carrying a
    false Python version, a false NumPy version, a false platform and a false
    configuration digest, and passed it to the production proposal — which
    accepted it and then recorded the *live* environment beside it, so the
    record disagreed with the plan it claimed to be frozen against and nothing
    said so.

    So the plan is **reconstructed**, exactly as the proposal is: its declared
    design inputs are fed back through :func:`pilot_plan` together with the
    **trusted live** source fingerprint and the live schema versions, and every
    field is compared.  A plan whose environment, source digest, schema
    versions, derived grid, derived window, derived shadow count, derived run
    names or configuration digest is not what those inputs imply is refused
    before any decision is taken from it.

    ``purpose`` is checked here too, and it is the reason a verifier fixture
    plan cannot reach a production decision: the production doors ask for
    ``approved_production_pilot`` and a fixture declares itself otherwise.
    """
    require_exact_type(plan, PilotPlan, "plan")
    require_exact_type(fingerprint, SourceFingerprint, "fingerprint")
    expected_purpose = _require_member(purpose, PILOT_PURPOSES, "purpose")
    if plan.purpose != expected_purpose:
        raise ValueError(
            f"this pilot plan declares purpose {plan.purpose!r} and this door "
            f"takes {expected_purpose!r}; a plan built to exercise the count "
            f"mechanics is not the pilot a checkpoint would approve, and the "
            f"difference is declared rather than inferred from its clock count")
    require_plan_firewall(plan, firewall)
    rebuilt = pilot_plan(
        label=plan.label, firewall=firewall, fingerprint=fingerprint,
        couplings=plan.couplings, clocks=plan.clocks,
        grid_half_width=plan.grid_half_width, timestep=plan.timestep,
        pulse_duration=plan.pulse_duration, pulse_centre=plan.pulse_centre,
        phase_diffusion=plan.phase_diffusion,
        lock_tolerance=plan.lock_tolerance, dwell_time=plan.dwell_time,
        trials_per_cell=plan.trials_per_cell,
        shadow_fraction=plan.shadow_fraction,
        run_name_prefix=plan.run_name_prefix, clock_block=plan.clock_block,
        step_window=plan.step_window, model=plan.model, purpose=plan.purpose)
    disagree = [name for name in PilotPlan.__dataclass_fields__
                if name != "seal"
                and getattr(plan, name) != getattr(rebuilt, name)]
    if disagree:
        raise ValueError(
            f"this pilot plan disagrees with the plan its own declared inputs "
            f"build, on {disagree}; the seal permits construction and proves "
            f"nothing about provenance, so the environment, the source digest, "
            f"the schema versions and every derived value are reconstructed "
            f"from the live package and compared")
    if plan.digest != rebuilt.digest:
        raise ValueError("the pilot plan's digest is not the digest its own "
                         "fields imply")
    return plan


def require_plan_firewall(plan, firewall):
    """Return ``plan`` if it was frozen against **this** firewall.

    The two records describe one boundary from two sides — the firewall owns
    the streams, the directory families and the count thresholds, the plan owns
    the physics — and every authoritative door takes both, so the join has to
    be checked rather than assumed.
    """
    require_exact_type(plan, PilotPlan, "plan")
    require_exact_type(firewall, PilotFirewall, "firewall")
    if plan.firewall_digest != firewall.digest:
        raise ValueError(
            f"this pilot plan was frozen against firewall "
            f"{plan.firewall_digest} and the firewall supplied is "
            f"{firewall.digest}; the approved pilot is specific to a stream, a "
            f"run-name family and a set of count thresholds")
    if plan.stream_namespace != firewall.pilot_namespace:
        raise ValueError("the plan and the firewall name different pilot "
                         "stream namespaces")
    return plan


def _quarantine_root():
    """The quarantine directory, resolved from the package directory only."""
    return (_PACKAGE_DIR / _raw_runner.RESULTS_DIRNAME / QUARANTINE_DIRNAME)


def _counts_digest(pairs):
    """The counts digest, over recounted values, run names **and** closures.

    One ordered record per run, name and counts and closure together, so the
    digest cannot be assembled by pairing one run's name with another run's
    closure.
    """
    return _digest({"counts": [
        {"run_name": name, "coupling": item.coupling, "trials": item.trials,
         "committed": item.committed, "closure": closure}
        for item, closure, name in pairs]})


@dataclass(frozen=True)
class QuarantinedPilotRun:
    """A closed pilot run addressed **by name** under the quarantine root.

    There is deliberately no directory field.  An independent review built one
    of these pointing at an ordinary, unquarantined run in the production
    reader's own root and the count door accepted it, because ``_read_quarantined``
    trusted the path it was handed.  The path is now derived from the run name
    and resolved through pinned directory descriptors with ``O_NOFOLLOW``, so
    constructing this record is harmless: it names a run, and naming a run that
    is not in the quarantine simply fails at the door.

    ``closure_digest`` binds the **whole closure** — the marker's schema
    version, all three ``(table, rows, digest)`` triples and the manifest
    digest — not one table.  The clock and shadow tables the gate checked are
    now durably bound to any later selection.
    """

    run_name: str
    closure_digest: str
    trials: int
    coupling: float

    def __post_init__(self):
        object.__setattr__(self, "run_name",
                           _require_text(self.run_name, "run_name",
                                         maximum=96))
        object.__setattr__(self, "closure_digest",
                           _require_hex_digest(self.closure_digest,
                                               "closure_digest"))
        object.__setattr__(self, "trials",
                           require_positive_int(self.trials, "trials",
                                                minimum=1))
        object.__setattr__(self, "coupling",
                           _require_positive(self.coupling, "coupling"))

    @property
    def digest(self) -> str:
        return _digest({"run_name": self.run_name,
                        "closure_digest": self.closure_digest,
                        "trials": self.trials, "coupling": self.coupling})


_DIR_FLAGS = os.O_RDONLY | getattr(os, "O_DIRECTORY", 0)
_NOFOLLOW = getattr(os, "O_NOFOLLOW", 0)


def _open_quarantined(run_name):
    """Yield the five payloads of a quarantined run, by name and nothing else.

    ``results/`` and the quarantine directory and the run directory are opened
    as pinned descriptors with ``O_NOFOLLOW | O_DIRECTORY``, and every file is
    opened relative to the run descriptor, also with ``O_NOFOLLOW``, and
    required to be a regular file.  A symbolic link standing in for the run
    directory or for any table is refused rather than followed, and nothing
    that happens after the directory is opened can redirect the reads.

    The only name resolved by string is the package directory, exactly as
    ``raw_runner`` does it.
    """
    name = _require_text(run_name, "run_name", maximum=96)
    if "/" in name or "\\" in name or name in (".", ".."):
        raise ValueError(f"run_name must be one path segment, got {name!r}")
    if not (_NOFOLLOW and getattr(os, "O_DIRECTORY", 0)):
        raise ValueError(
            "this platform cannot open a directory without following a "
            "symlink, so the quarantine cannot deliver the property it claims")
    try:
        results = os.open(str(_PACKAGE_DIR / _raw_runner.RESULTS_DIRNAME),
                          _DIR_FLAGS | _NOFOLLOW)
    except OSError as exc:
        raise ValueError(f"there is no quarantine to read {name!r} from") \
            from exc
    try:
        root = os.open(QUARANTINE_DIRNAME, _DIR_FLAGS | _NOFOLLOW,
                       dir_fd=results)
    except OSError as exc:
        raise ValueError(f"there is no quarantine directory, so {name!r} is "
                         f"not a quarantined run") from exc
    finally:
        os.close(results)
    try:
        run = os.open(name, _DIR_FLAGS | _NOFOLLOW, dir_fd=root)
    except OSError as exc:
        raise ValueError(
            f"{name!r} is not a run inside the quarantine; a run that has not "
            f"been quarantined is not readable through this door, whatever "
            f"else it may be readable through") from exc
    finally:
        os.close(root)
    payloads = {}
    try:
        for role, filename in _raw_runner.RUN_FILENAMES.items():
            try:
                handle = os.open(filename, os.O_RDONLY | _NOFOLLOW,
                                 dir_fd=run)
            except OSError as exc:
                raise ValueError(f"the quarantined run {name!r} is missing "
                                 f"{filename}") from exc
            with os.fdopen(handle, "rb") as stream:
                if not os.path.stat.S_ISREG(os.fstat(stream.fileno()).st_mode):
                    raise ValueError(f"{filename} is not a regular file")
                payloads[role] = stream.read()
    finally:
        os.close(run)
    return payloads


def _closure_digest(marker, manifest_bytes):
    """One digest over the whole five-file closure, not one table."""
    return _digest({
        "schema_version": marker.schema_version,
        "tables": [[name, rows, digest] for name, rows, digest in
                   marker.tables],
        "manifest": _raw_ledger.digest_bytes(manifest_bytes),
    })


def _read_quarantined(run_name):
    """Read a quarantined run by name and run the whole closure gate on it."""
    payloads = _open_quarantined(run_name)
    marker = _raw_ledger.parse_marker(payloads["marker"])
    ledger = _raw_ledger.require_closed_ledger(
        marker, payloads["manifest"], payloads["ledger"], payloads["clocks"],
        payloads["shadow"])
    return ledger, marker, _closure_digest(marker, payloads["manifest"])


def write_pilot_run(plan, firewall, coupling, run_name) -> QuarantinedPilotRun:
    """Write **the approved pilot** at one declared coupling, and return a receipt.

    Two things are authenticated, and an independent review found only the
    first of them present.

    *Where.*  No write-then-move.  The quarantine directory is opened as a
    pinned descriptor and ``raw_runner``'s private scoped writer creates the
    run beneath it, so a pilot run never exists under the production results
    root — not for an instant.  The production writer refuses the reserved
    namespace outright and the production reader refuses a manifest carrying
    it.

    *Which pilot.*  There is **no configuration parameter and no model
    parameter**.  The configuration is derived from the frozen
    :class:`PilotPlan` at the declared coupling, so the width-only, two-clock,
    ``dt=1``, two-trial run the review wrote through this door cannot be
    expressed: the plan names the model, the grid, the timestep, the pulse, the
    noise, the band, the dwell, the trial count and the shadow count, and a
    coupling it does not declare is refused.

    **Every check happens before any directory is created.**  The plan/firewall
    join, the run-name family on both records, and the configuration itself —
    which is where an impossible grid, timestep, trial count or pulse raises —
    are all resolved first, so a refusal leaves no quarantined directory behind
    to block a correct retry.  The review found a wrong ``run_label`` refused
    *after* the run had been written and closed; a label is no longer
    expressible, and nothing else can be refused late either.

    What comes back is a **receipt**, never a ``ClosedLedger``: the run name,
    the whole-closure digest, the trial count and the coupling.  The rows the
    closure gate read do not leave this function.
    """
    require_exact_type(plan, PilotPlan, "plan")
    # The plan is rebuilt from its own declared inputs and the live package
    # here too, so a false environment, source digest or schema version is
    # refused at the door as well as at the proposal.  The *purpose* is the
    # plan's own: a verifier fixture is allowed to exercise the count
    # mechanics, and is refused only where a decision is taken from it.
    require_authoritative_plan(plan, firewall, source_fingerprint(),
                               purpose=plan.purpose)
    # Validation first, in full: the plan rebuild above, the name against both
    # records *and* against the plan's frozen name/coupling mapping, and the
    # configuration derived from the plan.  Nothing below this line can raise
    # for a reason a caller controls, so nothing below it can leave residue.
    name = firewall.require_pilot_name(plan.require_run_name(run_name))
    if name != plan.run_name_for(coupling):
        raise ValueError(
            f"the plan writes coupling {coupling!r} under "
            f"{plan.run_name_for(coupling)!r} and this call names {name!r}; "
            f"which cell is written under which name is frozen in the plan")
    config = plan.config_for(coupling)
    root = _quarantine_root()
    root.mkdir(parents=True, exist_ok=True)
    results = os.open(str(_PACKAGE_DIR / _raw_runner.RESULTS_DIRNAME),
                      _DIR_FLAGS | _NOFOLLOW)
    try:
        parent = os.open(QUARANTINE_DIRNAME, _DIR_FLAGS | _NOFOLLOW,
                         dir_fd=results)
    finally:
        os.close(results)
    try:
        _raw_runner._write_scoped_run(
            config, name, parent, str(root / name), plan.model,
            plan.clock_block, plan.step_window)
    finally:
        os.close(parent)
    ledger, _, closure = _read_quarantined(name)
    firewall.require_pilot_run(ledger.manifest, name)
    plan.require_manifest(ledger.manifest, coupling)
    return QuarantinedPilotRun(
        run_name=name, closure_digest=closure,
        trials=int(ledger.manifest.value("trials")),
        coupling=float(ledger.manifest.value("peak_coupling")))


def _quarantine_pilot_run(firewall, run_name) -> QuarantinedPilotRun:
    """**Deprecated.** Move an already-written run into the quarantine.

    Private and deprecated, and it **refuses every current reserved identity**:
    a run on the ``xpilot-`` namespace cannot be created under the production
    root at all, so a mover for one would be a mover for something that cannot
    exist.  It is retained only for runs written before the direct writer, and
    it cannot become an alternate current workflow.

    New pilots go through :func:`write_pilot_run`, which never creates a
    production-root entry.

    ``raw_runner.open_raw_run`` takes one path segment under ``results/``.
    After this call the run lives at ``results/<QUARANTINE_DIRNAME>/<name>``,
    which that door cannot name.  The run's closure and hashes are untouched;
    what changes is who can open it.
    """
    require_exact_type(firewall, PilotFirewall, "firewall")
    name = firewall.require_pilot_name(run_name)
    source = Path(_raw_runner.raw_results_directory(name))
    if name.startswith(_raw_runner.EXPERIMENT_PILOT_NAMESPACE_PREFIX) or \
            firewall.pilot_namespace.startswith(
                _raw_runner.EXPERIMENT_PILOT_NAMESPACE_PREFIX):
        raise ValueError(
            f"{name!r} is a current reserved experiment pilot; those are "
            f"written directly into the quarantine by write_pilot_run and "
            f"never exist under the production root, so there is nothing here "
            f"to move and this deprecated path must not become an alternate "
            f"current workflow")
    if not source.is_dir():
        raise ValueError(f"there is no run at {name!r} to quarantine")
    root = _quarantine_root()
    root.mkdir(parents=True, exist_ok=True)
    target = root / name
    if target.exists():
        raise ValueError(f"a quarantined run already stands at {name!r}; "
                         f"overwriting one in place is how a rerun's rows end "
                         f"up mixed with a previous run's")
    shutil.move(str(source), str(target))

    ledger, _, closure = _read_quarantined(name)
    firewall.require_pilot_run(ledger.manifest, name)
    return QuarantinedPilotRun(
        run_name=name, closure_digest=closure,
        trials=int(ledger.manifest.value("trials")),
        coupling=float(ledger.manifest.value("peak_coupling")))


@dataclass(frozen=True)
class PilotCellCounts:
    """The **entire** result a pilot is allowed to produce for one coupling.

    **Exactly four values**: coupling, trials, committed, and the survivors
    derived from them.  No exponent, no curvature, no survival curve, no commit
    time, no dwell or reset count, no winner mismatch, no population detail —
    and, since an independent review counted it as a fifth visible field, no
    provenance token either.

    That means this record proves nothing about where it came from, and it is
    not asked to.  Provenance lives where it can be *checked*:
    :func:`select_range` and :func:`post_range_power` take quarantined **run
    names** and re-open and re-count every one of them internally, so a
    caller-built or ``dataclasses.replace``-edited record is a display value
    that the authoritative path never consults.
    """

    coupling: float
    trials: int
    committed: int

    def __post_init__(self):
        object.__setattr__(self, "coupling",
                           _require_positive(self.coupling, "coupling"))
        object.__setattr__(self, "trials",
                           require_positive_int(self.trials, "trials",
                                                minimum=1))
        object.__setattr__(self, "committed",
                           require_positive_int(self.committed, "committed",
                                                minimum=0))
        if self.committed > self.trials:
            raise ValueError(f"committed={self.committed} exceeds "
                             f"trials={self.trials}")

    @property
    def survivors(self) -> int:
        """Trials that did not commit — the unresolved fraction, as a count."""
        return self.trials - self.committed


def pilot_counts(firewall, quarantined) -> PilotCellCounts:
    """Project one quarantined pilot run down to its counts.

    The only supported path to a quarantined run in this package.  It resolves
    the run **by name** under the quarantine root through pinned descriptors,
    re-reads the five files, puts them through
    ``raw_ledger.require_closed_ledger``, confirms the run against the firewall
    by both manifest and directory family, checks the whole-closure digest
    against the handle, and then reduces the reconciled ledger to a category
    tally and throws the rest away.
    """
    require_exact_type(firewall, PilotFirewall, "firewall")
    require_exact_type(quarantined, QuarantinedPilotRun, "quarantined")
    ledger, _, closure = _read_quarantined(quarantined.run_name)
    firewall.require_pilot_run(ledger.manifest, quarantined.run_name)
    if closure != quarantined.closure_digest:
        raise ValueError(
            f"the quarantined run closes to {closure} and the handle carries "
            f"{quarantined.closure_digest}; the files changed after they were "
            f"quarantined")
    committed = 0
    for row in ledger.rows("ledger"):
        if row["category"] == "committed":
            committed += 1
    return PilotCellCounts(
        coupling=float(ledger.manifest.value("peak_coupling")),
        trials=int(ledger.manifest.value("trials")),
        committed=committed)


def _authoritative_cells(plan, firewall, run_names):
    """Re-open every named quarantined run and recount it, in coupling order.

    The authoritative join.  ``select_range``, ``reference_rate``,
    ``post_range_power`` and the proposal factory all go through here, so none
    of them can be driven by a count record a caller made up: the values they
    use are the ones this function read off disk through the closure gate a
    moment earlier.

    **The plan is rechecked here, not only at the writer.**  Every run's
    manifest is rebuilt from the frozen :class:`PilotPlan` at that run's own
    coupling and compared field for field, so a directory placed in the
    quarantine by any other route — a hand-written fixture, a legacy move, a
    run written before the plan was frozen — is refused at the count door as
    well.  Authenticating the writer alone would leave the recount believing
    whatever was on disk.
    """
    require_exact_type(plan, PilotPlan, "plan")
    # The plan is rebuilt from its own declared inputs and the live package
    # here too, so a false environment, source digest or schema version is
    # refused at the door as well as at the proposal.  The *purpose* is the
    # plan's own: a verifier fixture is allowed to exercise the count
    # mechanics, and is refused only where a decision is taken from it.
    require_authoritative_plan(plan, firewall, source_fingerprint(),
                               purpose=plan.purpose)
    # **Exactly the approved pilot's runs, in the approved order.**  Set,
    # length, order, uniqueness and the one-to-one coupling mapping are all
    # checked here, and a miss is an *authority* refusal rather than a
    # scientific one: an independent review read these four runs in reverse and
    # got a selection back, and omitted one and got ``run_too_short`` -- a
    # statement about coupling ranges standing in for "these are not that
    # pilot's runs".
    names = list(plan.require_run_names(run_names))
    pairs = []
    for name in names:
        ledger, _, closure = _read_quarantined(name)
        firewall.require_pilot_run(ledger.manifest, name)
        coupling = float(ledger.manifest.value("peak_coupling"))
        # The name/coupling mapping is frozen, so a run filed under another
        # cell's name is refused here rather than silently re-sorted into it.
        if coupling != plan.coupling_for(name):
            raise ValueError(
                f"the run {name!r} carries peak coupling {coupling!r} and the "
                f"plan writes {plan.coupling_for(name)!r} under that name; "
                f"which cell is written where is frozen in the plan")
        plan.require_manifest(ledger.manifest, coupling)
        committed = sum(1 for row in ledger.rows("ledger")
                        if row["category"] == "committed")
        # **The name travels with its own counts and its own closure.**  An
        # independent review found the closure digest built by zipping the
        # caller's original name order against a coupling-sorted closure list,
        # so a reordering silently paired a name with another run's closure.
        # There is one ordered record per run now, and nothing is zipped.
        pairs.append((PilotCellCounts(
            coupling=coupling,
            trials=int(ledger.manifest.value("trials")),
            committed=committed), closure, name))
    pairs.sort(key=lambda entry: entry[0].coupling)
    if len({item.coupling for item, _, _ in pairs}) != len(pairs):
        raise ValueError("two quarantined runs share one coupling; a pilot "
                         "cell is one coupling")
    if len({closure for _, closure, _ in pairs}) != len(pairs):
        raise ValueError("two quarantined runs close to the same digest")
    return tuple(pairs)


@dataclass(frozen=True)
class RangeSelection:
    """The one thing a pilot decides: which coupling range carries information.

    ``closure_digests``
        The whole-closure digest of every pilot run the counts came from, in
        coupling order — the marker, all three tables and the manifest, not one
        table.
    """

    low: float
    high: float
    cells: int
    refusal: str
    rule_digest: str
    counts_digest: str
    closure_digests: tuple

    def __post_init__(self):
        refusal = self.refusal
        if refusal:
            refusal = _require_member(refusal, SELECTION_REFUSALS, "refusal")
        elif not isinstance(refusal, str):
            raise TypeError(f"refusal must be a string, got "
                            f"{type(refusal).__name__}")
        object.__setattr__(self, "refusal", refusal)
        if refusal:
            for field in ("low", "high"):
                if require_finite(getattr(self, field), field) != 0.0:
                    raise ValueError(f"a refused selection carries {field}=0.0")
                object.__setattr__(self, field, 0.0)
            object.__setattr__(self, "cells",
                               require_positive_int(self.cells, "cells",
                                                    minimum=0))
            if self.cells != 0:
                raise ValueError("a refused selection selected zero cells")
        else:
            low = _require_positive(self.low, "low")
            high = _require_positive(self.high, "high")
            if not high > low:
                raise ValueError(f"high={high!r} must exceed low={low!r}")
            object.__setattr__(self, "low", low)
            object.__setattr__(self, "high", high)
            object.__setattr__(self, "cells",
                               require_positive_int(self.cells, "cells",
                                                    minimum=2))
        for field in ("rule_digest", "counts_digest"):
            object.__setattr__(self, field,
                               _require_hex_digest(getattr(self, field),
                                                   field))
        digests = [_require_hex_digest(value, f"closure_digests[{position}]")
                   for position, value in
                   enumerate(tuple(self.closure_digests))]
        if len(set(digests)) != len(digests):
            raise ValueError("two pilot cells carry the same closure digest; "
                             "one closed run is one coupling")
        if not digests:
            raise ValueError(
                "a range selection with no pilot closure digests behind it is "
                "a pair of numbers anybody could have typed")
        object.__setattr__(self, "closure_digests", tuple(digests))

    @property
    def selected(self) -> bool:
        """``True`` when a range was chosen."""
        return not self.refusal

    @property
    def digest(self) -> str:
        return _digest({"low": self.low, "high": self.high,
                        "cells": self.cells, "refusal": self.refusal,
                        "rule_digest": self.rule_digest,
                        "counts_digest": self.counts_digest,
                        "closure_digests": list(self.closure_digests)})


def _range_from_pairs(firewall, pairs):
    """The frozen count rule, over ``((counts, closure), ...)`` in order."""
    cells = [item for item, _, _ in pairs]
    rule = _digest({"minimum_events": firewall.minimum_events,
                    "minimum_survivors": firewall.minimum_survivors,
                    "minimum_cells": firewall.minimum_cells,
                    "permitted_fields": list(PERMITTED_PILOT_FIELDS)})
    fingerprint = _counts_digest(pairs)
    digests = tuple(closure for _, closure, _ in pairs)
    if len(cells) < firewall.minimum_cells:
        return RangeSelection(0.0, 0.0, 0, "run_too_short", rule, fingerprint,
                              digests)
    eligible = [item.committed >= firewall.minimum_events
                and item.survivors >= firewall.minimum_survivors
                for item in cells]
    if not any(eligible):
        return RangeSelection(0.0, 0.0, 0, "no_eligible_cells", rule,
                              fingerprint, digests)
    best = (0, 0)
    start = None
    for index, flag in enumerate(list(eligible) + [False]):
        if flag and start is None:
            start = index
        elif not flag and start is not None:
            if index - start > best[0]:
                best = (index - start, start)
            start = None
    length, begin = best
    if length < firewall.minimum_cells:
        return RangeSelection(0.0, 0.0, 0, "not_contiguous_enough", rule,
                              fingerprint, digests)
    return RangeSelection(cells[begin].coupling,
                          cells[begin + length - 1].coupling, length, "",
                          rule, fingerprint, digests)


@dataclass(frozen=True)
class PreviewRangeSelection:
    """The count rule applied to **display** counts.  A different type.

    An independent review asked for this to be distinguishable by type rather
    than by a convention about digests, and it is: no authoritative door
    accepts it, because every one of them demands exactly a
    :class:`RangeSelection`.  It carries the same fields so a reader can see
    what a pattern would mean, and it can enter nothing.
    """

    low: float
    high: float
    cells: int
    refusal: str

    def __post_init__(self):
        object.__setattr__(self, "low", require_nonnegative(self.low, "low"))
        object.__setattr__(self, "high", require_nonnegative(self.high, "high"))
        object.__setattr__(self, "cells",
                           require_positive_int(self.cells, "cells",
                                                minimum=0))
        refusal = self.refusal
        if refusal:
            _require_member(refusal, SELECTION_REFUSALS, "refusal")
        elif not isinstance(refusal, str):
            raise TypeError(f"refusal must be a string, got "
                            f"{type(refusal).__name__}")

    @property
    def selected(self) -> bool:
        """``True`` when the rule would have chosen a range."""
        return not self.refusal


def preview_range(firewall, counts) -> PreviewRangeSelection:
    """Apply the count rule to **display** counts.  Not authoritative.

    The ruling is explicit: caller-built or ``dataclasses.replace``-edited
    records may be display inputs only, and must mismatch an authoritative
    recomputation.  This door exists so that the rule's *arithmetic* can be
    exercised on declared patterns without writing runs to disk — and every
    selection it returns carries a closure digest of all zeros, which no real
    quarantined run can produce, so it can never be mistaken for one.
    """
    require_exact_type(firewall, PilotFirewall, "firewall")
    cells = [require_exact_type(item, PilotCellCounts, f"counts[{position}]")
             for position, item in enumerate(tuple(counts))]
    if not cells:
        raise ValueError("a preview needs at least one count record")
    if len({item.coupling for item in cells}) != len(cells):
        raise ValueError("two count records share one coupling; a pilot cell "
                         "is one coupling, and a duplicate would weight it "
                         "twice in a rule that counts consecutive cells")
    # Ordered ``(counts, closure, run name)`` records, the same shape the
    # authoritative recount builds -- with an all-zero-ish closure no real run
    # can produce and a display-only name, so a preview can never be mistaken
    # for a snapshot.
    pairs = tuple((item, f"{index:064x}", f"preview-cell-{index}")
                  for index, item in enumerate(
                      sorted(cells, key=lambda entry: entry.coupling)))
    outcome = _range_from_pairs(firewall, pairs)
    return PreviewRangeSelection(outcome.low, outcome.high, outcome.cells,
                                 outcome.refusal)


def select_range(plan, firewall, run_names) -> RangeSelection:
    """Choose a coupling range from counts alone, under the frozen rule.

    The rule, entire:

    1. a cell is **eligible** when its committed count is at least
       ``firewall.minimum_events`` *and* its survivor count is at least
       ``firewall.minimum_survivors`` — one threshold for events and one for
       survivors, because a cell that saturated carries as little slope
       information as one that never fired;
    2. the selection is the **longest run of consecutive eligible cells**, ties
       broken by taking the weakest such run, so the rule is deterministic; and
    3. it must contain at least ``firewall.minimum_cells`` cells, or the pilot
       refuses.

    Nothing in it reads an exponent, a curvature, a survival shape, a commit
    time, a dwell or reset count, a noise strength, a band, a pulse, a
    population, a state machine, a likelihood or an exclusion — because
    :class:`PilotCellCounts`, the only input type this function accepts, does
    not carry any of them.
    """
    return _range_from_pairs(firewall,
                             _authoritative_cells(plan, firewall, run_names))


def _rate_of(selection):
    """:data:`REFERENCE_RATE_RULE` in code, over an authoritative selection.

    One expression, called by every door that needs the number, so a manifest
    and a checkpoint cannot quote two different rates for one range.  It is a
    call into ``dynamics.fixed_contraction_rate`` rather than a re-derivation,
    so the rate written into a manifest is the same expression the width-only
    control's own module uses.
    """
    require_exact_type(selection, RangeSelection, "selection")
    if not selection.selected:
        raise ValueError(
            f"the count-only rule refused with {selection.refusal!r}, so there "
            f"is no range and no geometric midpoint; the reference rate is "
            f"mechanically computable from a selected range and from nothing "
            f"else")
    return _dynamics.fixed_contraction_rate(selection.low, selection.high)


#: The post-range plan's construction token.  A third object, for the same
#: reason as the other two: reaching one is not reaching the others, and none
#: of them is an authority.
_POST_RANGE_SEAL = object()


@dataclass(frozen=True)
class PostRangePlan:
    """The whole post-pilot state, under one digest, from one snapshot.

    An independent review found the post-range join spread across records that
    could not be joined back to one another: a ``RangeSelection`` that
    authenticated itself, a ``PowerEstimate`` with no way to name the range
    that sized it, and a proposal whose source table said the rate came from a
    count-only rule while the component handing it over was self-authenticating.

    This record is the join.  It is **factory-only** and every field in it is
    derived, in **one** authoritative reopen-and-recount snapshot, from the
    frozen :class:`PilotPlan`, the :class:`PilotFirewall` it was frozen
    against, and the ordered quarantined run names:

    ``pilot_plan_digest``, ``firewall_digest``
        Which pilot was approved, and on which boundary.
    ``run_names``, ``closure_digest``
        The ordered runs the snapshot read, and one digest over all of their
        whole-five-file closures, so a different set of runs — or the same runs
        after an edit — is a different plan.
    ``selection_digest``, ``counts_digest``
        The range the frozen count rule chose from that snapshot, and the
        recounted numbers it chose from.
    ``target_digest``, ``power_digest``, ``power_verdict``
        The sampling target the sweep was re-sized against and the estimate
        that came out.
    ``status``
        **Derived**, and re-derived in ``__post_init__`` from ``refusal`` and
        ``power_verdict``, so a record whose status disagrees with its own
        fields is refused rather than believed.
    """

    pilot_plan_digest: str
    firewall_digest: str
    run_names: tuple
    closure_digest: str
    counts_digest: str
    selection_digest: str
    low: float
    high: float
    cells: int
    refusal: str
    reference_rate: object
    target_digest: str
    power_digest: str
    power_verdict: str
    status: str
    seal: object = None

    def __post_init__(self):
        if self.seal is not _POST_RANGE_SEAL:
            raise TypeError(
                "PostRangePlan is factory-only: build it with post_range_plan, "
                "which takes the frozen plan, the firewall and the ordered "
                "quarantined run names and derives every field from one "
                "authoritative recount.  A caller-built post-range record is "
                "exactly the self-authenticating component this type replaces")
        for field in ("pilot_plan_digest", "firewall_digest",
                      "closure_digest", "counts_digest", "target_digest",
                      "power_digest"):
            object.__setattr__(self, field,
                               _require_hex_digest(getattr(self, field),
                                                   field))
        names = tuple(_require_text(value, f"run_names[{position}]",
                                    maximum=96)
                      for position, value in enumerate(tuple(self.run_names)))
        if not names:
            raise ValueError("a post-range plan reads at least one "
                             "quarantined pilot run")
        object.__setattr__(self, "run_names", names)
        refusal = self.refusal
        if refusal:
            refusal = _require_member(refusal, SELECTION_REFUSALS, "refusal")
        elif not isinstance(refusal, str):
            raise TypeError(f"refusal must be a string, got "
                            f"{type(refusal).__name__}")
        object.__setattr__(self, "refusal", refusal)
        object.__setattr__(self, "power_verdict",
                           _require_member(self.power_verdict, POWER_VERDICTS,
                                           "power_verdict"))
        if refusal:
            object.__setattr__(self, "selection_digest",
                               _require_hex_digest(self.selection_digest,
                                                   "selection_digest"))
            if self.reference_rate is not None:
                raise ValueError(
                    "a refused range has no geometric midpoint and therefore "
                    "no reference rate")
        else:
            object.__setattr__(self, "selection_digest",
                               _require_hex_digest(self.selection_digest,
                                                   "selection_digest"))
            object.__setattr__(self, "reference_rate",
                               _require_positive(self.reference_rate,
                                                 "reference_rate"))
        derived = ("post_range_refused" if refusal else
                   "post_range_powered" if self.power_verdict == "powered"
                   else "post_range_no_powered_matrix")
        if self.status != derived:
            raise ValueError(
                f"this post-range plan declares {self.status!r} while its own "
                f"refusal {refusal!r} and power verdict "
                f"{self.power_verdict!r} give {derived!r}; the status is "
                f"derived from the join and is not a field to be chosen")
        object.__setattr__(self, "status",
                           _require_member(derived, POST_RANGE_STATUSES,
                                           "status"))

    @property
    def selected(self) -> bool:
        """``True`` when the count rule named an information-bearing range."""
        return not self.refusal

    @property
    def digest(self) -> str:
        return _digest({
            "pilot_plan_digest": self.pilot_plan_digest,
            "firewall_digest": self.firewall_digest,
            "run_names": list(self.run_names),
            "closure_digest": self.closure_digest,
            "counts_digest": self.counts_digest,
            "selection_digest": self.selection_digest,
            "low": self.low, "high": self.high, "cells": self.cells,
            "refusal": self.refusal, "reference_rate": self.reference_rate,
            "target_digest": self.target_digest,
            "power_digest": self.power_digest,
            "power_verdict": self.power_verdict, "status": self.status,
        })


def _post_range_snapshot(target, plan, firewall, run_names):
    """One reopen-and-recount, and everything the post-range state derives.

    Returns ``(selection, power, plan_record)``.  Every caller in this module
    goes through here, so the selection a proposal writes, the rate it writes,
    the estimate that sized it and the record that binds them all come from
    **one** read of the quarantine rather than from four.
    """
    require_exact_type(target, SamplingTarget, "target")
    pairs = _authoritative_cells(plan, firewall, run_names)
    selection = _range_from_pairs(firewall, pairs)
    power = _power_from_snapshot(target, selection, pairs)
    rate = _rate_of(selection) if selection.selected else None
    # The names, the couplings and the closures come out of the **same ordered
    # records** the recount built.  Nothing is zipped across two orderings: an
    # independent review found the caller's original name order zipped against
    # a coupling-sorted closure list, so a reordered call paired each name with
    # another run's closure and the digest still looked well formed.
    record = PostRangePlan(
        pilot_plan_digest=plan.digest, firewall_digest=firewall.digest,
        run_names=tuple(name for _, _, name in pairs),
        closure_digest=_digest({"closures": [
            [name, item.coupling, closure] for item, closure, name in pairs]}),
        counts_digest=selection.counts_digest,
        selection_digest=selection.digest, low=selection.low,
        high=selection.high, cells=selection.cells, refusal=selection.refusal,
        reference_rate=rate, target_digest=target.digest,
        power_digest=power.digest, power_verdict=power.verdict,
        status=("post_range_refused" if selection.refusal else
                "post_range_powered" if power.verdict == "powered"
                else "post_range_no_powered_matrix"),
        seal=_POST_RANGE_SEAL)
    return selection, power, record


#: Every physical field that more than one of the sampling target, the primary
#: matrix and the frozen pilot plan carries, named once.
#:
#: An independent review multiplied every coupling of an otherwise valid,
#: honestly resource-priced six-node matrix by 1.1 and put it through the
#: proposal, the authoritative rebuild and the checkpoint.  All three accepted
#: it, and the record it produced said ``peak_couplings[0] = 0.55`` beside
#: ``pilot_candidate_couplings[0] = 0.5`` -- a matrix and a pilot describing two
#: different experiments, each internally consistent, each with its own honest
#: digest.  Every join in this module compared *counts and digests*; none
#: compared the physics the records share.
#:
#: This table is that comparison, written once so a field cannot be forgotten:
#: the couplings themselves at exact ordered binary64 equality, the grid, the
#: timing, the pulse, the noise, the band, the dwell, the arm names and their
#: clock counts, the arm the pilot runs, and the arm whose contraction rate the
#: reference-rate rule is about.
SHARED_PHYSICAL_FIELDS = (
    "peak_couplings", "cells", "clocks", "grid_half_width", "grid_support",
    "timestep", "pulse_duration", "pulse_centre", "phase_diffusion",
    "lock_tolerance", "dwell_time", "shadow_fraction", "arms", "arm_clocks",
    "pilot_arm", "reference_rate_arm",
)


def shared_physical_identity(target, matrix, plan) -> dict:
    """``{field: {record: value}}`` over every shared physical field.

    Only the records that actually carry a field appear under it, so the
    comparison is over genuine overlaps rather than over invented defaults.
    """
    require_exact_type(target, SamplingTarget, "target")
    require_exact_type(matrix, PrimaryMatrix, "matrix")
    require_exact_type(plan, PilotPlan, "plan")
    envelopes = tuple(target.arm_envelopes)
    table = {
        "peak_couplings": {
            "SamplingTarget": tuple(target.couplings),
            "PrimaryMatrix": tuple(matrix.couplings),
            "PilotPlan": tuple(plan.couplings)},
        "cells": {
            "SamplingTarget": int(target.cells),
            "PrimaryMatrix": len(matrix.couplings),
            "PilotPlan": plan.cells},
        "clocks": {
            "PrimaryMatrix": int(matrix.clocks),
            "PilotPlan": int(plan.clocks)},
        "grid_half_width": {
            "PrimaryMatrix": float(matrix.grid_half_width),
            "PilotPlan": float(plan.grid_half_width)},
        "grid_support": {
            "PrimaryMatrix": (-float(matrix.grid_half_width),
                              float(matrix.grid_half_width)),
            "PilotPlan": tuple(float(value) for value in plan.grid_support)},
        "timestep": {
            "PrimaryMatrix": float(matrix.timestep),
            "PilotPlan": float(plan.timestep)},
        "pulse_duration": {
            "PrimaryMatrix": float(matrix.pulse_duration),
            "PilotPlan": float(plan.pulse_duration)},
        "pulse_centre": {
            "PrimaryMatrix": float(matrix.pulse_centre),
            "PilotPlan": float(plan.pulse_centre)},
        "phase_diffusion": {
            "PrimaryMatrix": float(matrix.phase_diffusion),
            "PilotPlan": float(plan.phase_diffusion)},
        "lock_tolerance": {
            "PrimaryMatrix": float(matrix.lock_tolerance),
            "PilotPlan": float(plan.lock_tolerance)},
        "dwell_time": {
            "PrimaryMatrix": float(matrix.dwell_time),
            "PilotPlan": float(plan.dwell_time)},
        "shadow_fraction": {
            "SamplingTarget": float(target.shadow_fraction),
            "PrimaryMatrix": float(matrix.shadow_fraction),
            "PilotPlan": float(plan.shadow_fraction)},
        "arms": {
            "SamplingTarget": tuple(item.name for item in envelopes),
            "PrimaryMatrix": tuple(matrix.arms)},
        "arm_clocks": {
            "SamplingTarget": tuple((item.name, int(item.clocks))
                                    for item in envelopes),
            "PrimaryMatrix": tuple((name, int(clocks))
                                   for name, clocks in matrix.arm_clocks)},
        # The arm the pilot actually runs, and the arm the reference-rate rule
        # is *about*.  A pilot on a control arm, or a matrix with no
        # width-only control for the contraction rate to describe, is a
        # different experiment from the one the rule was frozen for.
        "pilot_arm": {
            "PrimaryMatrix": str(matrix.refinement_arm),
            "PilotPlan": str(plan.model)},
        "reference_rate_arm": {
            "SamplingTarget": tuple(sorted(
                item.name for item in envelopes
                if not item.pairs_with_full or item.name.endswith("control"))),
            "PrimaryMatrix": tuple(sorted(
                name for name in matrix.arms if name.endswith("control")))},
    }
    return {field: table[field] for field in SHARED_PHYSICAL_FIELDS}


def require_shared_physical_identity(target, matrix, plan):
    """Refuse three records that do not describe **one** experiment.

    Exact equality, field by field, including exact ordered binary64 equality
    of the coupling grids.  The refusal names the disagreeing field and shows
    what each record said, because "these describe different experiments" is
    useless without saying where.
    """
    disagree = []
    for field, values in shared_physical_identity(target, matrix,
                                                  plan).items():
        distinct = list(values.values())
        if any(value != distinct[0] for value in distinct[1:]):
            disagree.append(
                f"{field}: " + ", ".join(f"{record} says {value!r}"
                                         for record, value in values.items()))
    if disagree:
        raise ValueError(
            "the sampling target, the primary matrix and the frozen pilot plan "
            "do not describe one experiment.  They disagree on "
            + "; ".join(disagree)
            + ".  Each record is internally consistent and carries its own "
              "honest digest; agreeing digests are not agreeing physics, and a "
              "proposal that serializes two coupling grids is a proposal for "
              "two experiments")
    return target


def _require_post_range_design(power, target, matrix):
    """Refuse a post-range estimate beside a target and matrix it disagrees with.

    An independent review found the proposal factory serializing two
    conflicting designs at once: a four-cell post-range estimate calling for
    8 419/22 211/8 419 trials per cell beside the supplied pre-pilot six-cell
    matrix carrying 2 406/6 346/2 406, each with its own digest, in one record.
    A proposal is a single design or it is not a proposal.

    Deriving the coherent post-range target, matrix and resource estimate is
    ticket-08 work and this ticket does not do it, so this refuses **by name**
    until those records exist, rather than skipping the joins because run names
    happen to be present.  It is factored out so the branch can be exercised
    directly: the outer canonical-plan gate means no fixture plan reaches the
    factory, and a rule that cannot be reached cannot be trusted.
    """
    require_exact_type(power, PowerEstimate, "power")
    require_exact_type(target, SamplingTarget, "target")
    require_exact_type(matrix, PrimaryMatrix, "matrix")
    if power.target_digest != target.digest or \
            dict(matrix.arm_trials) != power.arm_trials:
        raise ValueError(
            f"{POST_RANGE_MATRIX_NOT_DERIVED}: the count-only rule selected a "
            f"range and the sweep re-sizes to {power.arm_trials} over "
            f"{power.cells} cells, while the matrix supplied is the pre-pilot "
            f"one carrying {dict(matrix.arm_trials)} over "
            f"{len(matrix.couplings)} cells.  A post-range proposal needs a "
            f"target, a matrix and a resource estimate re-derived from the "
            f"post-range state and joined exactly; until those exist this "
            f"refuses rather than serializing two designs in one record")
    return power


def _require_snapshot_power(power, joined):
    """Refuse a power estimate that is not the authoritative snapshot's.

    Factored out so the join can be exercised directly: an independent review
    handed the proposal factory a ``dataclasses.replace``-edited estimate --
    six cells and the provisional counts, relabelled ``post_range_join`` with
    an invented selection digest -- beside four real quarantined runs, and the
    factory recorded the forgery.  After a pilot the sizing is whatever
    re-reading the quarantined runs says it is, and a record handed in beside
    them is not evidence about them.
    """
    require_exact_type(power, PowerEstimate, "power")
    require_exact_type(joined, PowerEstimate, "joined")
    if power.digest != joined.digest:
        raise ValueError(
            f"this proposal carries a power estimate that is not the one the "
            f"authoritative post-range recount produces: it hashes to "
            f"{power.digest} and the snapshot's own estimate hashes to "
            f"{joined.digest}")
    return joined


def post_range_plan(target, plan, firewall, run_names) -> PostRangePlan:
    """The authoritative post-range state, from one snapshot.  Factory-only."""
    return _post_range_snapshot(target, plan, firewall, run_names)[2]


def reference_rate(plan, firewall, run_names) -> float:
    """The width-only control's one contraction rate, recomputed from disk.

    **An authoritative door, not a helper.**  An independent review built a
    :class:`RangeSelection` out of six invented digests, handed it to the
    previous signature, and got a contraction rate back; the rate then went
    into a proposal whose source table said it came from a count-only rule.
    Nothing about that record was true and nothing in this module could tell.

    So there is no selection parameter.  This takes the frozen
    :class:`PilotPlan`, the :class:`PilotFirewall` it was frozen against, and
    the ordered names of the quarantined pilot runs; it re-opens and re-counts
    every one of them through the five-file closure gate and the plan's own
    manifest rebuild in **one** snapshot, applies the frozen count rule to that
    snapshot, and computes the rate from the range the rule chose.  A
    caller-built selection cannot reach it, because a selection is no longer
    something a caller can hand it.
    """
    pairs = _authoritative_cells(plan, firewall, run_names)
    return _rate_of(_range_from_pairs(firewall, pairs))


# ---------------------------------------------------------------------------
# The proposed production manifest
# ---------------------------------------------------------------------------

#: The seal a proposal factory holds and a caller does not.
#:
#: The same line this package already draws at the analysis layer's derived
#: records: the unsupported route requires writing a leading underscore, which
#: makes it impossible by accident and visible when taken.  It is not a
#: cryptographic boundary and is not described as one.  An independent review
#: built a complete proposal in which every non-special field was the string
#: ``"wrong-but-json-valid"`` and the public type accepted and hashed it — a
#: complete wrong record being more dangerous than an omitted field, because it
#: looks frozen.  Two things close that: this token, and the per-field rules
#: below, which refuse a wrong value even when the token is present.
_SEAL = object()

#: What every declared proposal field must *be*, checked by value and not only
#: by presence.
#:
#: ``kind`` is one of ``digest``, ``text``, ``enum``, ``int``, ``nonneg_int``,
#: ``probability``, ``positive``, ``finite``, ``bool_false``, ``list_positive``,
#: ``list_text``, ``list_number``, ``pair``, ``optional_positive``.  The tuple
#: after it, where present, is the enum's members or the exact required value.
_FIELD_RULES = {
    "schema": ("enum", (EXPERIMENTS_SCHEMA,)),
    "content_revision": ("enum", (CONTENT_REVISION,)),
    "status": ("enum", PROPOSAL_STATUSES),
    "status_reason": ("text", None),
    "proposal_label": ("text", None),
    "source_digest": ("digest", None),
    "python_version": ("text", None),
    "numpy_version": ("text", None),
    "platform": ("text", None),
    "ledger_schema_version": ("enum",
                              (_raw_ledger.LEDGER_SCHEMA_VERSION,)),
    "manifest_schema_version": ("enum",
                                (_raw_ledger.MANIFEST_SCHEMA_VERSION,)),
    "production_status_schema": ("enum",
                                 (_compare.PRODUCTION_STATUS_SCHEMA,)),
    "peak_couplings": ("list_positive", None),
    "pulse_duration": ("positive", None),
    "pulse_centre": ("finite", None),
    "phase_diffusion": ("positive", None),
    "lock_tolerance": ("positive", None),
    "dwell_time": ("positive", None),
    "clocks": ("int", None),
    "arm_clocks": ("pair", None),
    "grid_half_width": ("positive", None),
    "grid_spacing": ("positive", None),
    "grid_parity": ("enum", ("even", "odd")),
    "grid_origin": ("finite", None),
    "grid_support": ("list_number", None),
    "minimum_eligible_clocks": ("int", None),
    "eligible_at_weakest": ("int", None),
    "discrete_staircase_fallback": ("text", None),
    "timestep": ("positive", None),
    "timestep_ladder": ("list_positive", None),
    "refinement_ratio": ("positive", None),
    "refinement_couplings": ("list_positive", None),
    "refinement_arm": ("text", None),
    "numerical_gate": ("enum", (_raw_runner.NUMERICAL_GATE,)),
    "endpoint_status": ("enum", (_raw_runner.STATIONARY_ORACLE_STATUS,)),
    "moving_band_status": ("enum", (_raw_runner.MOVING_BAND_STATUS,)),
    "numerical_budget_digest": ("digest", None),
    "coverage_sigma": ("positive", None),
    "allowed_fraction": ("probability", None),
    "numerical_verdict": ("enum", BUDGET_VERDICTS),
    "numerical_blockers": ("list_blockers", None),
    "evidence_digest": ("digest", None),
    "search_digest": ("digest", None),
    "search_verdict": ("enum", FEASIBILITY_VERDICTS),
    "trials_per_cell": ("int", None),
    "arm_trials": ("pair", None),
    "arms": ("list_text", None),
    "shadow_fraction": ("probability", None),
    "resampling_unit": ("enum", _analysis.DESIGN_RESAMPLING_UNITS),
    "fit_model": ("enum", (_analysis.LINK_FORMULA,)),
    "confidence_level": ("enum", (CONFIDENCE_LEVEL,)),
    "exponent_half_width_target": ("positive", None),
    "minimum_causal_contrast": ("positive", None),
    "contrast_covariance_floor": ("finite", None),
    "measured_arm_latent_correlation": ("finite", None),
    "recovery_digest": ("digest", None),
    "curvature_test": ("text", None),
    "lack_of_fit_test": ("text", None),
    "inclusion_rule": ("text", None),
    "sampling_target_digest": ("digest", None),
    "power_digest": ("digest", None),
    "power_stage": ("enum", POWER_STAGES),
    "power_verdict": ("enum", POWER_VERDICTS),
    "stream_namespace": ("text", None),
    "pilot_stream_namespace": ("text", None),
    "run_directory_prefix": ("text", None),
    "pilot_run_directory_prefix": ("text", None),
    "pilot_quarantine_directory": ("enum", (QUARANTINE_DIRNAME,)),
    "pilot_trials_eligible": ("bool_false", None),
    "firewall_digest": ("digest", None),
    "pilot_plan_digest": ("digest", None),
    "pilot_candidate_couplings": ("list_positive", None),
    "pilot_trials_per_cell": ("int", None),
    "reference_rate_rule": ("enum", (REFERENCE_RATE_RULE,)),
    "fixed_contraction_rate": ("optional_positive", None),
    "selection_digest": ("optional_digest", None),
    "post_range_digest": ("optional_digest", None),
    "post_range_status": ("optional_enum", POST_RANGE_STATUSES),
    "pilot_ledger_fingerprints": ("list_text", None),
    "benchmark_digest": ("digest", None),
    "matrix_digest": ("digest", None),
    "resource_digest": ("digest", None),
    "design_digest": ("digest", None),
    "stop_rule": ("text", None),
    "failure_rule": ("text", None),
    "non_claims": ("list_text", None),
}


def _check_field(key, value, kind, allowed):
    """Refuse one proposal field that is the wrong *kind of thing*."""
    if kind == "digest":
        text = _require_text(value, key, maximum=_raw_ledger.DIGEST_HEX_LENGTH)
        if len(text) != _raw_ledger.DIGEST_HEX_LENGTH or \
                any(character not in "0123456789abcdef" for character in text):
            raise ValueError(f"{key} must be a "
                             f"{_raw_ledger.DIGEST_HEX_LENGTH}-character hex "
                             f"digest, got {value!r}")
    elif kind == "optional_digest":
        if value is not None:
            _check_field(key, value, "digest", None)
    elif kind == "optional_enum":
        if value is not None:
            _check_field(key, value, "enum", allowed)
    elif kind == "text":
        _require_text(value, key, maximum=800)
    elif kind == "enum":
        if value not in allowed:
            raise ValueError(f"{key} must be one of {list(allowed)}, got "
                             f"{value!r}")
    elif kind == "int":
        require_positive_int(value, key, minimum=1)
    elif kind == "nonneg_int":
        require_positive_int(value, key, minimum=0)
    elif kind == "probability":
        _require_fraction(value, key)
    elif kind == "positive":
        _require_positive(value, key)
    elif kind == "optional_positive":
        if value is not None:
            _require_positive(value, key)
    elif kind == "finite":
        require_finite(value, key)
    elif kind == "bool_false":
        if value is not False:
            raise ValueError(f"{key} must be False, got {value!r}")
    elif kind == "list_positive":
        entries = list(value) if isinstance(value, (list, tuple)) else None
        if not entries:
            raise ValueError(f"{key} must be a non-empty list, got {value!r}")
        for position, entry in enumerate(entries):
            _require_positive(entry, f"{key}[{position}]")
    elif kind == "list_number":
        entries = list(value) if isinstance(value, (list, tuple)) else None
        if not entries:
            raise ValueError(f"{key} must be a non-empty list, got {value!r}")
        for position, entry in enumerate(entries):
            require_finite(entry, f"{key}[{position}]")
    elif kind == "list_blockers":
        entries = list(value) if isinstance(value, (list, tuple)) else None
        if entries is None:
            raise ValueError(f"{key} must be a list, got {value!r}")
        unknown = [entry for entry in entries
                   if entry not in BUDGET_BLOCKERS]
        if unknown:
            raise ValueError(
                f"{key} names {unknown}, which are outside the closed "
                f"vocabulary {list(BUDGET_BLOCKERS)}; an independent review "
                f"put ['looks_good'] into a sealed proposal and it was hashed, "
                f"because this field was checked for being a list of strings "
                f"and not for what the strings meant")
    elif kind == "list_text":
        entries = list(value) if isinstance(value, (list, tuple)) else None
        if entries is None:
            raise ValueError(f"{key} must be a list, got {value!r}")
        for position, entry in enumerate(entries):
            _require_text(entry, f"{key}[{position}]", maximum=800)
    elif kind == "pair":
        entries = list(value) if isinstance(value, (list, tuple)) else None
        if not entries:
            raise ValueError(f"{key} must be a non-empty list of pairs, got "
                             f"{value!r}")
        for position, entry in enumerate(entries):
            pair = list(entry) if isinstance(entry, (list, tuple)) else None
            if pair is None or len(pair) != 2:
                raise ValueError(f"{key}[{position}] must be a pair, got "
                                 f"{entry!r}")
            _require_text(pair[0], f"{key}[{position}][0]", maximum=80)
            require_positive_int(pair[1], f"{key}[{position}][1]", minimum=1)
    else:                                            # pragma: no cover
        raise ValueError(f"unknown field rule {kind!r} for {key!r}")


@dataclass(frozen=True)
class ProposedProductionManifest:
    """A production manifest **proposal**, complete and deliberately unsigned.

    Its key set is exactly :data:`PROPOSED_MANIFEST_KEYS`: one missing key or
    one extra key is a refusal.  Every value is additionally checked against
    :data:`_FIELD_RULES` — exact enums, digests that are digests, positive
    numbers that are positive, lists whose entries are the right kind of thing
    — so a record that is JSON-valid and scientifically false is refused rather
    than hashed.

    It is **factory-only**.  ``proposed_production_manifest`` derives every
    field from component records that are themselves exact types, so the values
    are not merely *checked*, they are *reconstructed*; the public constructor
    demands a token it does not export.

    Three things it structurally cannot be:

    * it cannot say ``production_cleared``.  The gate is an enum of exactly the
      diagnostic value ticket 05 writes;
    * it cannot claim ticket 06's production status schema.  It records
      ``compare.PRODUCTION_STATUS_SCHEMA`` as a *number it does not meet*; and
    * it cannot be approved or signed.  There is no signature field and
      :data:`PROPOSAL_STATUSES` has exactly one member.
    """

    fields: tuple
    seal: object = None

    def __post_init__(self):
        if self.seal is not _SEAL:
            raise TypeError(
                "ProposedProductionManifest is factory-only: build it with "
                "proposed_production_manifest, which derives every field from "
                "component records rather than taking values on trust.  An "
                "independent review assembled a complete proposal out of the "
                "string 'wrong-but-json-valid' and this type hashed it")
        entries = {}
        order = []
        for position, entry in enumerate(tuple(self.fields)):
            pair = tuple(entry)
            if len(pair) != 2:
                raise ValueError(f"fields[{position}] must be a (key, value) "
                                 f"pair, got {pair!r}")
            key = _require_text(pair[0], f"fields[{position}] key", maximum=80)
            if key in entries:
                raise ValueError(f"{key!r} appears twice; a manifest that "
                                 f"names one field twice has two answers and "
                                 f"a reader has no way to tell which is meant")
            entries[key] = pair[1]
            order.append(key)
        missing = [key for key in PROPOSED_MANIFEST_KEYS if key not in entries]
        extra = [key for key in order if key not in PROPOSED_MANIFEST_KEYS]
        if missing:
            raise ValueError(
                f"the proposal omits {missing}; every field in "
                f"PROPOSED_MANIFEST_KEYS is required, because a proposal that "
                f"can be built without one is a proposal that can be signed "
                f"without one")
        if extra:
            raise ValueError(
                f"the proposal carries {extra}, which are outside "
                f"PROPOSED_MANIFEST_KEYS; a field nobody declared is a field "
                f"nobody reviewed")
        for key in PROPOSED_MANIFEST_KEYS:
            kind, allowed = _FIELD_RULES[key]
            _check_field(key, entries[key], kind, allowed)

        # Cross-field relationships.  Each of these is a way two individually
        # legal values can disagree about the same fact.
        if entries["numerical_gate"] == _compare.PRODUCTION_GATE:
            raise ValueError("a proposal cannot clear itself")
        couplings = list(entries["peak_couplings"])
        if any(later <= earlier
               for earlier, later in zip(couplings, couplings[1:])):
            raise ValueError(f"peak_couplings must be strictly increasing, got "
                             f"{couplings}")
        if entries["grid_half_width"] <= couplings[-1]:
            raise ValueError("grid_half_width does not exceed the strongest "
                             "coupling")
        spacing = 2.0 * entries["grid_half_width"] / entries["clocks"]
        if abs(entries["grid_spacing"] - spacing) > 1e-12 * spacing:
            raise ValueError(f"grid_spacing={entries['grid_spacing']!r} is not "
                             f"2 * half_width / clocks = {spacing!r}")
        parity = "even" if entries["clocks"] % 2 == 0 else "odd"
        if entries["grid_parity"] != parity:
            raise ValueError(f"grid_parity says {entries['grid_parity']!r} for "
                             f"{entries['clocks']} clocks")
        if list(entries["grid_support"]) != [-entries["grid_half_width"],
                                             entries["grid_half_width"]]:
            raise ValueError("grid_support is not the declared half-width")
        ladder = list(entries["timestep_ladder"])
        if ladder[0] != entries["timestep"]:
            raise ValueError("the ladder does not start at the production step")
        if any(later >= earlier for earlier, later in zip(ladder, ladder[1:])):
            raise ValueError("the ladder is not strictly refining")
        ratio = ladder[0] / ladder[1]
        if abs(entries["refinement_ratio"] - ratio) > 1e-12 * ratio:
            raise ValueError("refinement_ratio is not the ladder's own ratio")
        if not set(entries["refinement_couplings"]).issubset(set(couplings)):
            raise ValueError("refinement_couplings are not matrix couplings")
        if len(entries["refinement_couplings"]) >= len(couplings):
            raise ValueError("the ladder refines every cell, which is the "
                             "automatic full factorial the plan forbids")
        arm_names = [pair[0] for pair in entries["arm_clocks"]]
        if sorted(arm_names) != sorted(entries["arms"]):
            raise ValueError(f"arm_clocks names {sorted(arm_names)} and arms "
                             f"says {sorted(entries['arms'])}")
        if entries["refinement_arm"] not in arm_names:
            raise ValueError("refinement_arm is not one of the arms")
        for required in ("full", "central_control", "width_only_control"):
            if required not in arm_names:
                raise ValueError(f"the arms omit {required!r}")
        if dict(pair for pair in entries["arm_clocks"])["full"] != \
                entries["clocks"]:
            raise ValueError("the full arm's grid is not the declared "
                             "population")
        if dict(pair for pair in entries["arm_clocks"])["central_control"] != 1:
            raise ValueError("the central control is not one clock")
        if entries["eligible_at_weakest"] > entries["clocks"]:
            raise ValueError("more clocks are eligible than exist")
        if entries["numerical_verdict"] != "satisfied" and \
                not entries["numerical_blockers"]:
            raise ValueError("a blocked numerical verdict with no blockers")
        if entries["production_status_schema"] <= \
                entries["manifest_schema_version"]:
            raise ValueError(
                "production_status_schema must sit above the live manifest "
                "schema, which is what keeps physical_verdict_permitted "
                "unreachable")
        if entries["power_stage"] == "post_range_join" and \
                entries["selection_digest"] is None:
            raise ValueError("a post-range power stage with no selection")
        if entries["power_stage"] == "provisional_pre_pilot" and \
                entries["fixed_contraction_rate"] is not None:
            raise ValueError(
                "a provisional proposal names a contraction rate; the rule is "
                "frozen and the number is not computable until a range exists")

        # The serialized design digest must be the digest of the design.  An
        # earlier revision wrote a placeholder of sixty-four zeros into the
        # record and computed the real value in a property beside it, so the
        # persisted manifest carried a design digest that was false on
        # construction.  It is recomputed here and compared.
        expected = _digest({key: (list(entries[key])
                                  if isinstance(entries[key], tuple)
                                  else entries[key])
                            for key in DESIGN_DIGEST_KEYS})
        if entries["design_digest"] != expected:
            raise ValueError(
                f"the serialized design_digest is {entries['design_digest']} "
                f"and the design hashes to {expected}; a proposal cannot carry "
                f"a freeze that is not its own")

        # Derived relationships the seal must not be able to carry past.
        if (entries["numerical_verdict"] == "satisfied") != (
                not entries["numerical_blockers"]):
            raise ValueError(
                f"numerical_verdict={entries['numerical_verdict']!r} with "
                f"blockers {entries['numerical_blockers']}: a satisfied budget "
                f"has no blockers and a blocked one has at least one")
        if entries["search_verdict"] == "feasible_matrix" and \
                entries["numerical_verdict"] != "satisfied":
            raise ValueError("a feasible search verdict over a blocked budget")
        if entries["power_verdict"] == "no_powered_matrix" and \
                entries["status"] != "checkpoint_blocked":
            raise ValueError("an unpowered design that is not blocked")
        declared = dict(pair for pair in entries["arm_trials"])
        if sorted(declared) != sorted(entries["arms"]):
            raise ValueError(f"arm_trials names {sorted(declared)} and arms "
                             f"says {sorted(entries['arms'])}")
        if declared.get("full") != entries["trials_per_cell"]:
            raise ValueError(
                f"trials_per_cell={entries['trials_per_cell']} is not the full "
                f"arm's own count {declared.get('full')}")
        floor = entries["contrast_covariance_floor"]
        measured = entries["measured_arm_latent_correlation"]
        if not -1.0 <= measured <= 1.0:
            raise ValueError("the measured arm latent correlation is not a "
                             "correlation")
        if measured < floor:
            raise ValueError(
                f"the measured arm latent correlation {measured} is below the "
                f"frozen conservative floor {floor}; the floor exists to be "
                f"the adverse bound the ensemble does not reach, and a "
                f"measurement past it means the design was sized against a "
                f"covariance the experiment can beat")

        object.__setattr__(self, "fields",
                           tuple((key, entries[key])
                                 for key in PROPOSED_MANIFEST_KEYS))
        object.__setattr__(self, "seal", None)

    def value(self, key: str):
        """The value of one declared field, or ``KeyError`` naming it."""
        key = _require_text(key, "key", maximum=80)
        for name, value in self.fields:
            if name == key:
                return value
        raise KeyError(f"{key!r} is not one of the "
                       f"{len(PROPOSED_MANIFEST_KEYS)} declared proposal "
                       f"fields")

    def as_dict(self) -> dict:
        """A plain mapping, in declared order, for printing."""
        return {key: (list(value) if isinstance(value, tuple) else value)
                for key, value in self.fields}

    @property
    def status(self) -> str:
        """``checkpoint_blocked``.  There is no other member."""
        return self.value("status")

    @property
    def design_digest(self) -> str:
        """The **measurement-free** freeze, over :data:`DESIGN_DIGEST_KEYS`.

        Reproducible on any machine from the same declared inputs, because no
        wall time and no traced peak is inside it.  This is the hash a reviewer
        freezes; :attr:`digest` additionally binds the measurement and
        therefore moves when the benchmark is re-run, which is correct for a
        record of what a machine did and wrong for a scientific freeze.
        """
        return self.value("design_digest")

    @property
    def digest(self) -> str:
        """The whole-proposal hash: design **and** the measurement it carries."""
        return _digest(self.as_dict())


def proposed_production_manifest(label, fingerprint, matrix, target, power,
                                 disposition, search, benchmark, resources,
                                 firewall, status_reason, recovery, plan,
                                 pilot_run_names=(),
                                 status="checkpoint_blocked"):
    """Assemble the dry-run proposal from the records that justify it.

    Every argument is an *exact* type built by this module's own factories, and
    every field is **derived** from them — no value is taken on trust and none
    is passed in.  ``fixed_contraction_rate`` is not an argument at all: it is
    computed here from ``selection`` through :func:`reference_rate`, so the
    manifest cannot carry a rate that is not the one the count-only rule
    implies.  Before a pilot there is no selection, so the rate is ``None`` and
    the rule alone is recorded.
    """
    require_exact_type(fingerprint, SourceFingerprint, "fingerprint")
    require_exact_type(matrix, PrimaryMatrix, "matrix")
    require_exact_type(target, SamplingTarget, "target")
    require_exact_type(power, PowerEstimate, "power")
    require_exact_type(disposition, BudgetDisposition, "disposition")
    require_exact_type(search, SearchReport, "search")
    require_exact_type(benchmark, BenchmarkResult, "benchmark")
    require_exact_type(resources, ResourceEstimate, "resources")
    require_exact_type(firewall, PilotFirewall, "firewall")
    require_exact_type(recovery, RecoveryReport, "recovery")
    # The approved pilot is part of the proposal, not a footnote beside it.
    # An independent review found the proposal's whole pilot vocabulary to be
    # two namespaces, two prefixes and an eligibility flag -- nothing that says
    # *which* pilot, at which couplings, at which configuration, for how many
    # trials.  The plan is bound here, by digest, inside the design freeze.
    # The plan is **rebuilt**, not merely type-checked, and the production
    # doors take exactly the approved-production purpose -- so a verifier
    # fixture plan, or one carrying a false environment, source digest, schema
    # version or configuration digest, is refused before any field is derived
    # from it.
    require_canonical_pilot_plan(plan, firewall, fingerprint)
    # **One experiment, or a named refusal.**  Every other join in this factory
    # compares counts and digests; this one compares the physics the three
    # records share, at exact ordered binary64 equality on the coupling grids.
    require_shared_physical_identity(target, matrix, plan)
    if plan.source_digest != fingerprint.digest:
        raise ValueError(
            f"the pilot plan was frozen against source {plan.source_digest} "
            f"and the fingerprint is {fingerprint.digest}; a plan frozen "
            f"against different bytes describes a different package")
    label = _require_text(label, "label", maximum=120)
    # **Derived, never supplied.**  An independent review rehashed a proposal
    # carrying a measured correlation of 0.9 and the guard accepted it, because
    # the same false scalar could be handed to the rebuild.  The Gaussian
    # copula realizes ``r^2 * rho`` in the generated variables, so the number
    # is computed here from the recovery record and from nothing else.
    measured_arm_correlation = recovery.correlation ** 2 \
        * recovery.arm_correlation
    # The benchmark must be a measurement of *these* bytes.  A record whose
    # source digest is not the current fingerprint is a measurement of another
    # package, however well-formed it looks.
    if benchmark.source_digest != fingerprint.digest:
        raise ValueError(
            f"the benchmark records source {benchmark.source_digest} and the "
            f"fingerprint is {fingerprint.digest}; a timing of different bytes "
            f"is not a timing of this proposal")
    status_reason = _require_text(status_reason, "status_reason", maximum=400)
    status = _require_member(status, PROPOSAL_STATUSES, "status")
    names = tuple(pilot_run_names)

    # **One snapshot, taken before anything is joined.**  After a pilot the
    # authoritative estimate is whatever re-reading the quarantined runs says,
    # so it is taken first and every join below is made against it.
    if names:
        selection, joined, post_range = _post_range_snapshot(
            target, plan, firewall, names)
        _require_snapshot_power(power, joined)
        power = joined
        rate = post_range.reference_rate
        # **Fail closed, by name.**  An independent review found this factory
        # serializing two conflicting designs at once: a four-cell post-range
        # power estimate calling for 8 419/22 211/8 419 trials per cell beside
        # the supplied pre-pilot six-cell matrix carrying 2 406/6 346/2 406,
        # each with its own digest, in one record.  A proposal is a single
        # design or it is not a proposal.
        #
        # Deriving the coherent post-range target, matrix and resource
        # estimate is ticket-08 work and this ticket does not do it, so the
        # honest behaviour is to refuse until those records exist rather than
        # to skip the joins because run names happen to be present.
        _require_post_range_design(power, target, matrix)
    else:
        selection = None
        post_range = None
        rate = None
        if power.stage != "provisional_pre_pilot":
            raise ValueError(
                f"this proposal names no quarantined pilot runs and carries a "
                f"{power.stage!r} power estimate; a post-range join without a "
                f"pilot to join to is a claim about a range nobody measured")

    if power.target_digest != target.digest:
        raise ValueError("the power estimate was computed from a different "
                         "sampling target than the one supplied")
    if resources.matrix_digest != matrix.digest:
        raise ValueError("the resource estimate was computed from a different "
                         "matrix than the one supplied")
    if resources.benchmark_digest != benchmark.digest:
        raise ValueError("the resource estimate was scaled from a different "
                         "benchmark than the one supplied")
    if search.budget_digest != disposition.budget_digest:
        raise ValueError("the feasibility search used a different numerical "
                         "budget than the disposition")
    if search.evidence_digest != disposition.evidence_digest:
        raise ValueError("the feasibility search used different evidence than "
                         "the disposition")
    # The matrix and the power estimate must agree arm by arm, at either
    # stage.  After a pilot both have to be the *post-range* records, which is
    # what the named refusal above requires.
    if dict(matrix.arm_trials) != power.arm_trials:
        raise ValueError(
            f"the matrix runs {dict(matrix.arm_trials)} and the power estimate "
            f"calls for {power.arm_trials}; the proposed trial counts follow "
            f"the frozen precision target arm by arm, and a matrix that does "
            f"not carry those numbers is not the matrix the target sized")
    # One snapshot, one selection, one rate — all recomputed here from the
    # named quarantined runs, never accepted from a caller.
    body = ((
        ("schema", EXPERIMENTS_SCHEMA),
        ("content_revision", CONTENT_REVISION),
        ("status", status),
        ("status_reason", status_reason),
        ("proposal_label", label),
        ("source_digest", fingerprint.digest),
        ("python_version", fingerprint.python_version),
        ("numpy_version", fingerprint.numpy_version),
        ("platform", fingerprint.platform),
        ("ledger_schema_version", _raw_ledger.LEDGER_SCHEMA_VERSION),
        ("manifest_schema_version", _raw_ledger.MANIFEST_SCHEMA_VERSION),
        ("production_status_schema", _compare.PRODUCTION_STATUS_SCHEMA),
        ("peak_couplings", list(matrix.couplings)),
        ("pulse_duration", matrix.pulse_duration),
        ("pulse_centre", matrix.pulse_centre),
        ("phase_diffusion", matrix.phase_diffusion),
        ("lock_tolerance", matrix.lock_tolerance),
        ("dwell_time", matrix.dwell_time),
        ("clocks", matrix.clocks),
        ("arm_clocks", [list(pair) for pair in matrix.arm_clocks]),
        ("grid_half_width", matrix.grid_half_width),
        ("grid_spacing", matrix.grid_spacing),
        ("grid_parity", matrix.grid_parity),
        ("grid_origin", 0.0),
        ("grid_support", [-matrix.grid_half_width, matrix.grid_half_width]),
        ("minimum_eligible_clocks", matrix.minimum_eligible_clocks),
        ("eligible_at_weakest", matrix.eligible_at_weakest),
        ("discrete_staircase_fallback", matrix.discrete_staircase_fallback),
        ("timestep", matrix.timestep),
        ("timestep_ladder", list(matrix.timestep_ladder)),
        ("refinement_ratio", matrix.timestep_ladder[0]
         / matrix.timestep_ladder[1]),
        ("refinement_couplings", list(matrix.refinement_couplings)),
        ("refinement_arm", matrix.refinement_arm),
        ("numerical_gate", _raw_runner.NUMERICAL_GATE),
        ("endpoint_status", _raw_runner.STATIONARY_ORACLE_STATUS),
        ("moving_band_status", _raw_runner.MOVING_BAND_STATUS),
        ("numerical_budget_digest", disposition.budget_digest),
        ("coverage_sigma", COVERAGE_SIGMA),
        ("allowed_fraction", ALLOWED_NUMERICAL_FRACTION),
        ("numerical_verdict", disposition.verdict),
        ("numerical_blockers", list(disposition.blockers)),
        ("evidence_digest", disposition.evidence_digest),
        ("search_digest", search.digest),
        ("search_verdict", search.verdict),
        ("trials_per_cell", matrix.trials_per_cell),
        ("arm_trials", [list(pair) for pair in matrix.arm_trials]),
        ("arms", list(matrix.arms)),
        ("shadow_fraction", matrix.shadow_fraction),
        ("resampling_unit", target.resampling_unit),
        ("fit_model", _analysis.LINK_FORMULA),
        ("confidence_level", target.confidence_level),
        ("exponent_half_width_target", target.exponent_half_width),
        ("minimum_causal_contrast", target.minimum_causal_contrast),
        ("contrast_covariance_floor", target.contrast_covariance_floor),
        # Two different numbers, named apart and never conflated.  The floor is
        # the frozen adverse covariance the contrast is *sized* against; the
        # measurement is what the recovery ensemble's generator actually
        # realizes, which for a Gaussian copula at per-cell correlation r and a
        # requested between-arm correlation rho is r^2 * rho.  An independent
        # review pointed out that the ensemble therefore does not exercise the
        # declared floor, and it is right -- so both travel, and the record
        # refuses a measurement past the floor.
        ("measured_arm_latent_correlation", measured_arm_correlation),
        ("recovery_digest", recovery.digest),
        ("curvature_test", "analysis.cluster_curvature, a nuisance-adjusted "
                           "cluster-robust score test on the W-residualized "
                           "quadratic column"),
        ("lack_of_fit_test", "analysis.cluster_lack_of_fit, the same test on "
                             "the cells-minus-two lack-of-fit contrasts"),
        ("inclusion_rule", "every declared coupling cell enters the "
                           "likelihood, including zero-event and all-event "
                           "cells; no cell is dropped, transformed or "
                           "window-selected"),
        ("sampling_target_digest", target.digest),
        ("power_digest", power.digest),
        ("power_stage", power.stage),
        ("power_verdict", power.verdict),
        ("stream_namespace", firewall.production_namespace),
        ("pilot_stream_namespace", firewall.pilot_namespace),
        ("run_directory_prefix", firewall.production_prefix),
        ("pilot_run_directory_prefix", firewall.pilot_prefix),
        ("pilot_quarantine_directory", QUARANTINE_DIRNAME),
        ("pilot_trials_eligible", False),
        ("firewall_digest", firewall.digest),
        ("pilot_plan_digest", plan.digest),
        ("pilot_candidate_couplings", list(plan.couplings)),
        ("pilot_trials_per_cell", plan.trials_per_cell),
        ("reference_rate_rule", REFERENCE_RATE_RULE),
        ("fixed_contraction_rate", rate),
        ("selection_digest",
         selection.digest if selection is not None else None),
        ("post_range_digest",
         post_range.digest if post_range is not None else None),
        ("post_range_status",
         post_range.status if post_range is not None else None),
        ("pilot_ledger_fingerprints",
         list(selection.closure_digests) if selection is not None else []),
        ("benchmark_digest", benchmark.digest),
        ("matrix_digest", matrix.digest),
        ("resource_digest", resources.digest),
        ("stop_rule", matrix.stop_rule),
        ("failure_rule", matrix.failure_rule),
        ("non_claims", list(EXPERIMENT_NON_CLAIMS)),
    ))
    held = dict(body)
    design = _digest({key: (list(held[key]) if isinstance(held[key], tuple)
                            else held[key])
                      for key in DESIGN_DIGEST_KEYS})
    return ProposedProductionManifest(body + (("design_digest", design),),
                                      _SEAL)


def require_authoritative_proposal(proposal, label, fingerprint, matrix,
                                   target, power, disposition, search,
                                   benchmark, resources, firewall,
                                   status_reason, recovery, plan,
                                   pilot_run_names=(),
                                   status="checkpoint_blocked"):
    """Return ``proposal`` only if the components **rebuild it exactly**.

    The authority.  An independent review took a valid factory record, changed
    same-typed fields — ``numerical_verdict="satisfied"`` with an empty blocker
    list and ``search_verdict="feasible_matrix"``; arm trials 7/8/7; a measured
    arm correlation of 0.9 — recomputed the declared design digest exactly as
    the class does, supplied the module's own reachable seal, and had all three
    accepted.  Every one of those records was internally consistent.  That is
    precisely the point: **a self-digest proves byte consistency, not
    provenance**, and neither does a token.

    So a serialized proposal is accepted here only when the authoritative
    component records — the numerical disposition and the feasibility search,
    the power plan and the matrix's own arm and cell identities, the recovery
    measurement, the firewall, the selection and the resource estimate —
    rebuild it **field for field**.  Nothing is trusted about the record; it is
    reconstructed and compared.

    A reviewer may freeze a design digest only after this join.
    """
    require_exact_type(proposal, ProposedProductionManifest, "proposal")
    rebuilt = proposed_production_manifest(
        label=label, fingerprint=fingerprint, matrix=matrix, target=target,
        power=power, disposition=disposition, search=search,
        benchmark=benchmark, resources=resources, firewall=firewall,
        status_reason=status_reason, recovery=recovery, plan=plan,
        pilot_run_names=pilot_run_names, status=status)
    disagree = [key for key, value in proposal.fields
                if rebuilt.value(key) != value]
    if disagree:
        raise ValueError(
            f"this proposal disagrees with the components that should have "
            f"produced it, on {disagree}; a record that hashes to itself is "
            f"internally consistent and says nothing about where its values "
            f"came from, which is how three same-typed false records passed "
            f"an earlier revision of this guard")
    return proposal


def require_unchanged_proposal(proposal, digest):
    """Return ``proposal`` if its **design** still hashes to ``digest``.

    The post-freeze guard, and it compares the measurement-free
    :attr:`~ProposedProductionManifest.design_digest` rather than the whole
    record.  That is the fix for an independent review's finding that
    re-running the same source changed the putative freeze: what a reviewer
    freezes is the declared design, and the benchmark beside it is evidence
    about a machine.  A field edited after the digest was quoted still fails
    here, because every declared design field is inside that hash.
    """
    require_exact_type(proposal, ProposedProductionManifest, "proposal")
    digest = _require_hex_digest(digest, "digest")
    actual = proposal.design_digest
    if actual != digest:
        raise ValueError(
            f"this proposal's design hashes to {actual} and the frozen digest "
            f"is {digest}; something changed after the freeze, and a proposal "
            f"that changed after it was quoted is a different proposal")
    return proposal


# ---------------------------------------------------------------------------
# The checkpoint-1 handoff: one blocked option, one machinery-only diagnostic
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class CheckpointAlternative:
    """One thing John can be asked to choose, described and -- for one of the
    two -- priced.

    They are **not two costed alternatives**: the intended-configuration
    outline carries no cost, no range and no bound, and the vocabulary below
    makes that a property of the label rather than of the text.

    Neither alternative is "run the pilot".  With the numerical disposition
    blocking, a range chosen now feeds a sweep that cannot be interpreted
    whichever range it picks, so the choice on the table is between *gathering
    the missing evidence* and *narrowing the claim to what the present evidence
    can carry*.

    ``supports_class``
        **What the alternative would let the experiment say, as a closed
        enum.**  An independent review presented this record with
        ``supports = "a physical scaling claim"`` beside a costed Option A and
        the checkpoint accepted it, because the field was free text and free
        text can say anything.  The decision-bearing value is now a member of
        :data:`SUPPORT_CLASSES`, and there is no member that licenses a
        physical claim -- so the forged record is not a record this type can
        express, rather than one it happens to reject.
    ``supports``
        The same statement in plain words, subordinate to the enum and
        required to agree with it.
    """

    label: str
    summary: str
    supports: str
    cost: object
    prerequisites: tuple
    supports_class: str = SUPPORT_CLASSES[0]

    @property
    def is_priced(self) -> bool:
        """Whether a defensible cost exists.  ``False`` is a statement."""
        return self.cost is not None

    def __post_init__(self):
        object.__setattr__(self, "label",
                           _require_member(self.label, ALTERNATIVE_LABELS,
                                           "label"))
        for field in ("summary", "supports"):
            object.__setattr__(self, field,
                               _require_text(getattr(self, field), field,
                                             maximum=800))
        # ``None`` means **explicitly unpriced**.  An independent review
        # showed that Option A's quoted range timed the PDE oracle and the
        # post-generation bootstrap only, omitting the moving replay and the
        # stationary observation generation entirely, so it was neither an
        # estimate nor an upper bound.  A number that is not a bound is worse
        # than no number, and this field can now be absent.
        if self.cost is not None:
            require_exact_type(self.cost, ResourceEstimate, "cost")
        object.__setattr__(self, "supports_class",
                           _require_member(self.supports_class,
                                           SUPPORT_CLASSES, "supports_class"))
        # The two cross-rules that make the review's forged record
        # inexpressible: the unpriced option may not carry a cost, and the
        # priced one may not be presented as unpriced.
        required_class, must_be_priced = ALTERNATIVE_CONTRACT[self.label]
        if self.supports_class != required_class:
            raise ValueError(
                f"alternative {self.label!r} declares support class "
                f"{self.supports_class!r}; it is {required_class!r} by "
                f"construction, and an option that says it supports something "
                f"else is a different decision wearing this label")
        if self.is_priced != must_be_priced:
            raise ValueError(
                f"alternative {self.label!r} is "
                f"{'priced' if self.is_priced else 'unpriced'} and must be "
                f"{'priced' if must_be_priced else 'unpriced'}; the "
                f"intended-configuration outline has no defensible cost -- the "
                f"two timed kernels omit the moving replay and the stationary "
                f"endpoint generation entirely -- and the narrowed diagnostic "
                f"has an exact one")
        prerequisites = tuple(
            _require_text(value, f"prerequisites[{position}]", maximum=400)
            for position, value in enumerate(tuple(self.prerequisites)))
        if not prerequisites:
            raise ValueError("an alternative with no prerequisites is a "
                             "decision with nothing to decide")
        object.__setattr__(self, "prerequisites", prerequisites)

    @property
    def digest(self) -> str:
        return _digest({"label": self.label, "summary": self.summary,
                        "supports": self.supports,
                        "supports_class": self.supports_class,
                        "cost": (self.cost.digest
                                 if self.cost is not None else None),
                        "prerequisites": list(self.prerequisites)})


def intended_configuration_option() -> CheckpointAlternative:
    """Option A, built here rather than handed in: an **unpriced** outline.

    **It takes nothing.**  An alignment review pointed out that this record was
    still built from a caller-supplied :class:`ValidationCampaign` while the
    checkpoint's docstring claimed every derivation input was bound to the
    frozen proposal.  It was not: the 86-field proposal has no campaign digest,
    so there was no authoritative join to make, and a foreign but internally
    valid campaign reached the presented text through the reachable
    construction token.  Rather than claim a join that does not exist, the
    campaign is gone from this record and from the checkpoint: what Option A
    says is a statement about what was *not* measured, and that statement does
    not depend on how many stages someone's campaign happens to declare.

    An independent review handed ``checkpoint_one_handoff`` a *costed* Option A
    whose ``supports`` string was ``"a physical scaling claim"``, using genuine
    components everywhere else, and the checkpoint accepted it — the same
    self-consistency error the proposal rebuild was added to prevent, one layer
    later.  The fix is that there is nowhere to hand one in.  This record
    takes **no input at all**: it is the one fixed unpriced outline, ``cost``
    is ``None`` by construction and the support class is the one member that
    licenses nothing.

    The omission is the closure.  The two kernels that were timed — the PDE
    oracle solve and the paired bootstrap over already-generated observations —
    omit the moving replay that generates those observations and the stationary
    endpoint sampling entirely, so no range built on them is an estimate or an
    upper bound.  The four kernels that would have to be measured end to end
    before any number exists are named in the prerequisites.
    """
    return CheckpointAlternative(
        label="intended_configuration_validation",
        summary=(
            f"an **unpriced** intended-configuration validation outline: "
            f"predeclared stages, each attacking one "
            f"budget unit by one lever, with stationary refinement to dt/256 "
            f"and space refined separately from time, and moving stages "
            f"separating timestep bias from sampling error.  It carries **no "
            f"cost and no range**: the two kernels timed here -- the PDE "
            f"oracle solve and the paired bootstrap over already-generated "
            f"observations -- omit the moving replay that generates those "
            f"observations and the stationary endpoint sampling entirely, so "
            f"any number built on them would be neither an estimate nor an "
            f"upper bound.  Pricing it needs end-to-end benchmarks of four "
            f"kernels: the oracle solve, stationary endpoint generation, "
            f"moving replay/paired-leaf generation, and the bootstrap "
            f"comparison"),
        supports=(
            "nothing, and it is not costed.  No stage is promised sufficient, "
            "every miss remains a numerical_no_result, and no physical, "
            "causal or scaling claim is licensed under any outcome"),
        supports_class="unpriced_no_claim",
        cost=None,
        prerequisites=(
            "this option is explicitly unpriced and must not be presented as a "
            "costed alternative, or as one carrying a range",
            "end-to-end benchmarks of the stationary endpoint generation and "
            "the moving replay/paired-leaf generation are prerequisites to any "
            "number at all",
            "each stage's success, no-result and stop rules are frozen before "
            "it runs, and a miss remains a numerical_no_result",
            "the kernel units are never divided by one another"))


def narrowed_claim_option(matrix, benchmark, disposition,
                          safety_factor) -> CheckpointAlternative:
    """Option B, built here: the exact primary-only diagnostic sweep.

    Every number is derived.  The trial count is the numerical disposition's
    own ``probability_admissible_trials`` — the largest count at which the
    present worst probability bound still fits its allowance — not a constant
    beside it.  The run count is the matrix's arms times its coupling nodes.
    The rows are priced arm by arm through :func:`row_width_bound` at *this*
    option's trial count, with each distinct clock count taking its own blend,
    because a run whose rows are half ledger lines is not priced by a 64-clock
    blend.

    **Excluded, in the summary and in the prerequisites**: the refinement runs,
    the range pilot and every numerical-validation run.  What it supports is a
    machinery and diagnostic statement, and the support class says so.
    """
    require_exact_type(matrix, PrimaryMatrix, "matrix")
    require_exact_type(benchmark, BenchmarkResult, "benchmark")
    require_exact_type(disposition, BudgetDisposition, "disposition")
    safety = _require_positive(safety_factor, "safety_factor")
    trials = require_positive_int(disposition.probability_admissible_trials,
                                  "probability_admissible_trials", minimum=1)
    steps = matrix.steps_at(matrix.timestep)
    bound = row_width_bound(trials, matrix.clocks, steps)
    fixed_bytes = max(benchmark.fixed_bytes,
                      int(math.ceil(bound["header_bytes"])))
    rate = benchmark.clock_steps_per_second
    cells = len(matrix.couplings)
    # One line per distinct clock count, so the one-clock central control is
    # priced by its own blend rather than by the 64-clock arms'.
    grouped = {}
    for name, clocks in matrix.arm_clocks:
        grouped.setdefault(int(clocks), []).append(name)
    lines = []
    for clocks in sorted(grouped, reverse=True):
        arms = sorted(grouped[clocks])
        per_row = (bound["blended"] if clocks == matrix.clocks
                   else bound["central_blended"])
        lines.append(_cost_line(
            f"narrowed/{'-and-'.join(arms)}", len(arms) * cells, trials,
            clocks, steps, matrix.shadow_fraction, rate, per_row, fixed_bytes,
            safety))
    return CheckpointAlternative(
        label="narrowed_claim",
        summary=(
            f"a {sum(line.runs for line in lines)}-run primary-only "
            f"diagnostic sweep: the {cells}-cell {len(matrix.arms)}-arm matrix "
            f"at {trials} trials per cell -- the largest count at which the "
            f"present worst probability bound still fits its allowance -- and "
            f"nothing else.  It **excludes** the refinement runs, the range "
            f"pilot and any numerical validation, and it accepts an exponent "
            f"half-width far wider than the frozen target"),
        supports=(
            f"a machinery and diagnostic statement only: at {trials} trials "
            f"per cell the exponent interval is too wide for the frozen causal "
            f"contrast, and both ticket-04 time-unit bounds still fail at "
            f"every trial count, so this option does not restore a scaling "
            f"claim -- it delimits precisely what the current evidence can "
            f"carry, which is the behaviour of the pipeline and not a law"),
        supports_class="machinery_and_diagnostic_only",
        cost=ResourceEstimate(
            label="cost of the narrowed diagnostic matrix",
            matrix_digest=matrix.digest, benchmark_digest=benchmark.digest,
            lines=tuple(lines), safety_factor=safety,
            clock_steps_per_second=rate, bytes_per_row=bound["blended"],
            fixed_bytes_per_run=fixed_bytes, peak_bytes=benchmark.peak_bytes),
        prerequisites=(
            "the narrowed claim is written down before the run, not after",
            "the refinement runs, the range pilot and all numerical "
            "validation are explicitly outside this option and outside its "
            "cost",
            "the two time-unit bounds still fail, so no commit-time quantile "
            "or survival-shape statement is licensed under this option either",
            "the moving-band numerical_no_result still blocks scientific "
            "interpretation and is not waived by a wider allowance"))


def checkpoint_alternatives(matrix, benchmark, disposition,
                            safety_factor) -> tuple:
    """Both alternatives, in declared order, derived from the components.

    There is no ``campaign`` parameter: Option A is a statement about what was
    not measured and takes nothing, and Option B is priced from the matrix,
    the benchmark and the numerical disposition, each of which the checkpoint
    binds to the frozen proposal.
    """
    return (intended_configuration_option(),
            narrowed_claim_option(matrix, benchmark, disposition,
                                  safety_factor))


def pilot_cost_from_plan(plan, benchmark, safety_factor) -> CostLine:
    """Price the **approved** pilot, from the plan and the benchmark alone.

    An independent review handed the checkpoint a pilot ``CostLine`` of one
    second and it was accepted, because the boundary checked only that the
    argument was a ``CostLine``.  There is no argument now: the line is built
    here from the frozen :class:`PilotPlan`'s own cell count, trial count,
    clock count, step count and shadow fraction, and from the measured rate and
    the row bound at *that* configuration.  A pilot that costs a second is not
    a value this function can return for a plan that does not.
    """
    require_exact_type(plan, PilotPlan, "plan")
    require_exact_type(benchmark, BenchmarkResult, "benchmark")
    safety = _require_positive(safety_factor, "safety_factor")
    bound = row_width_bound(plan.trials_per_cell, plan.clocks, plan.steps)
    per_row = (bound["blended"] if plan.clocks > 1
               else bound["central_blended"])
    return _cost_line(
        "pilot/range-only", plan.cells, plan.trials_per_cell, plan.clocks,
        plan.steps, plan.shadow_fraction, benchmark.clock_steps_per_second,
        max(benchmark.bytes_per_row, per_row),
        max(benchmark.fixed_bytes, int(math.ceil(bound["header_bytes"]))),
        safety)


#: The checkpoint's construction token.  A separate object again, so that
#: reaching the plan's seal or the proposal's does not reach this one -- and,
#: like both of those, a permission to *construct* and never an authority: a
#: sealed handoff still has to rebuild, component for component, through
#: :func:`_require_authoritative_checkpoint`.
_CHECKPOINT_SEAL = object()


# ``slots=True`` is load-bearing and is the whole of the closure below: a
# frozen dataclass without it still has an instance dictionary, and a plain
# method is a non-data descriptor, so ``object.__setattr__(handoff, "summary",
# ...)`` shadows the class implementation outright.  With slots there is no
# dictionary to write into and no attribute outside the declared fields can
# exist at all.
@dataclass(frozen=True, slots=True)
class Checkpoint1Handoff:
    """Everything checkpoint 1 is asked to decide about, under one digest.

    ``feasibility`` is derived, never passed, and it is
    ``no_feasible_matrix_among_searched`` rather than an unqualified
    impossibility: what has been established is a statement about the finite
    declared option set the search enumerated.

    **Nothing a caller can substitute.**  An independent review constructed
    this record directly with every canonical component and one substitution --
    Option B's comparison ``CostLine`` in place of the pilot line -- and
    :meth:`summary` reported the approved pilot at about eleven hours and
    28 836 rows instead of its own thirty-two hours and 81 840.  The factory
    had already removed the ``alternatives`` and ``pilot_line`` *parameters*;
    the record still stored them.  A first repair sealed the record, and an
    alignment review pointed out that a seal is ordinary module state: a
    caller who reaches it builds the same forgery, calls the public
    :meth:`summary` immediately, and never presents it to any door.

    So the three values a caller could once choose are **not fields any more**.
    ``feasibility``, ``alternatives`` and ``pilot_line`` are properties derived
    from the components this record carries, so there is nothing to substitute
    and no path -- sealed, replaced, or otherwise -- by which
    :meth:`summary` or :attr:`digest` can report a caller's number.

    A second alignment review found two derivation *inputs* still open where
    the prose claimed otherwise.  The slowest benchmark -- which prices the
    pilot line, Option B and the presented rate -- was never compared with the
    proposal's own ``benchmark_digest``, and it is now.  And Option A was
    still built from a caller-supplied ``ValidationCampaign`` that the
    86-field proposal does not freeze at all; rather than claim a join that
    does not exist, the campaign is gone from this record and Option A takes
    nothing.

    A third alignment review then showed that deriving the values was still
    not enough, because the *record* could be poisoned: this class was frozen
    but not slotted, so ``object.__setattr__(handoff, "summary", ...)``
    installed an instance attribute that shadowed the method -- a plain
    method is a non-data descriptor and the instance dictionary wins -- and
    the authority door, which compared fields and digests, returned the
    poisoned object unchanged.  The record is **slotted** now: there is no
    instance dictionary, no undeclared attribute can exist, and no method or
    property on it can be shadowed.  A subclass could still carry one, and is
    refused by exact type at the door.

    What the
    record does carry, it **binds**: the plan, matrix, resources, firewall,
    power, recovery, search, target and benchmark it derives from must each be
    the one the canonical proposal froze, checked here at construction.

    The seal remains, and :func:`_require_authoritative_checkpoint` remains,
    as defence in depth rather than as the defence.
    """

    label: str
    ticket_08_derivation_status: str
    ticket_08_authorization_status: str
    benchmarks: tuple
    target: object
    power: object
    recovery: object
    disposition: object
    search: object
    matrix: object
    resources: object
    firewall: object
    proposal: object
    plan: object
    requested_decision: str
    seal: object = None

    def __post_init__(self):
        object.__setattr__(self, "label",
                           _require_text(self.label, "label", maximum=120))
        # Derived, never passed: the factory computes it and the vocabulary has
        # one member, so a checkpoint cannot imply that the ticket-08 matrix
        # derivation has started, or that ticket 08 is authorized.
        object.__setattr__(self, "ticket_08_derivation_status",
                           _require_member(self.ticket_08_derivation_status,
                                           TICKET_08_DERIVATION_STATUSES,
                                           "ticket_08_derivation_status"))
        # Derived and closed at one member, for the same reason: this ticket
        # cannot authorize ticket 08, so ``authorized`` is not a value the
        # vocabulary can spell rather than one a guard happens to reject.
        object.__setattr__(self, "ticket_08_authorization_status",
                           _require_member(
                               self.ticket_08_authorization_status,
                               TICKET_08_AUTHORIZATION_STATUSES,
                               "ticket_08_authorization_status"))
        benchmarks = tuple(
            require_exact_type(item, BenchmarkResult, f"benchmarks[{position}]")
            for position, item in enumerate(tuple(self.benchmarks)))
        if not benchmarks:
            raise ValueError("a feasibility handoff carries at least one "
                             "measured benchmark; an unmeasured feasibility "
                             "claim is the thing this ticket exists to replace")
        object.__setattr__(self, "benchmarks", benchmarks)
        for field, expected in (("target", SamplingTarget),
                                ("power", PowerEstimate),
                                ("recovery", RecoveryReport),
                                ("disposition", BudgetDisposition),
                                ("search", SearchReport),
                                ("matrix", PrimaryMatrix),
                                ("resources", ResourceEstimate),
                                ("firewall", PilotFirewall),
                                ("proposal", ProposedProductionManifest),
                                ("plan", PilotPlan)):
            require_exact_type(getattr(self, field), expected, field)
        # **The derivation's inputs are bound to the frozen proposal.**  The
        # three presented values below are computed from these components, so
        # a component that is not the one the canonical proposal froze would
        # compute a different answer under the same name.  Every digest the
        # proposal carries for a component this record derives from is
        # compared here, at construction.
        for field, key in (("plan", "pilot_plan_digest"),
                           ("matrix", "matrix_digest"),
                           ("resources", "resource_digest"),
                           ("firewall", "firewall_digest"),
                           ("power", "power_digest"),
                           ("recovery", "recovery_digest"),
                           ("search", "search_digest"),
                           ("target", "sampling_target_digest")):
            carried = getattr(self, field).digest
            frozen = self.proposal.value(key)
            if carried != frozen:
                raise ValueError(
                    f"this handoff's {field} hashes to {carried} and the "
                    f"proposal it presents froze {key} as {frozen}; the "
                    f"feasibility verdict, both alternatives and the pilot "
                    f"cost line are derived from these components, so a "
                    f"component the proposal did not freeze would present a "
                    f"different experiment under the same record")
        # **The one benchmark every derivation reads.**  An alignment review
        # found this binding claimed in prose and missing from the loop: the
        # pilot line, Option B and the presented rate all come from
        # ``_slowest``, and nothing compared it with the benchmark the
        # proposal froze, so a caller could swap the tuple for another valid
        # measurement and move every derived number.  The slowest supplied
        # benchmark **is** the frozen one, or this is not that checkpoint.
        slowest = sorted(benchmarks,
                         key=lambda item: item.clock_steps_per_second)[0]
        frozen_benchmark = self.proposal.value("benchmark_digest")
        if slowest.digest != frozen_benchmark:
            raise ValueError(
                f"the slowest supplied benchmark hashes to {slowest.digest} "
                f"and the proposal it presents froze benchmark_digest as "
                f"{frozen_benchmark}; the pilot cost line, Option B and the "
                f"presented clock-step rate are all derived from that one "
                f"measurement, so a different one is a different checkpoint")
        # The numerical disposition has no digest field of its own in the
        # proposal; what the proposal freezes is the verdict and the ordered
        # blockers, and those are exactly what the derived verdict reads.
        for label, carried, frozen in (
                ("numerical verdict", self.disposition.verdict,
                 self.proposal.value("numerical_verdict")),
                ("numerical blockers", list(self.disposition.blockers),
                 list(self.proposal.value("numerical_blockers")))):
            if carried != frozen:
                raise ValueError(
                    f"this handoff's {label} is {carried!r} and the proposal "
                    f"it presents froze {frozen!r}; the derived feasibility "
                    f"verdict reads the disposition, so the two cannot be "
                    f"different records")
        require_plan_firewall(self.plan, self.firewall)
        # **Not a caller's sentence.**  There is one decision text, it is
        # derived from :data:`CHECKPOINT_DECISION_CLAUSES`, and a handoff
        # presenting anything else -- including a shorter one that merely omits
        # a clause -- is refused here rather than shown.
        presented = _require_text(self.requested_decision,
                                  "requested_decision", maximum=800)
        if presented != CHECKPOINT_REQUESTED_DECISION:
            missing = [clause for clause in CHECKPOINT_DECISION_CLAUSES
                       if clause not in presented]
            raise ValueError(
                f"this handoff presents a decision text that is not the "
                f"canonical one"
                + (f", omitting {missing}" if missing else "")
                + ".  The presented sentence is what a reader is actually "
                  "shown, so it is derived from CHECKPOINT_DECISION_CLAUSES "
                  "and compared exactly; carrying the right fields beside a "
                  "sentence that says none of them is the defect this "
                  "refuses")
        object.__setattr__(self, "requested_decision", presented)
        if self.proposal.status != "checkpoint_blocked":
            raise ValueError("the proposal is not blocked")
        # Derived, so this can no longer be a contradiction a caller writes;
        # it is asserted because the *rule* is what a reader is being shown.
        if self.disposition.is_blocking and \
                self.feasibility == "feasible_matrix":
            raise AssertionError(
                f"the numerical disposition is {self.disposition.verdict!r} "
                f"with blockers {list(self.disposition.blockers)} while the "
                f"derived verdict is {self.feasibility!r}; a blocked "
                f"numerical budget is a blocked matrix")
        # **Last, deliberately.**  Every field rule above is a statement about
        # what a checkpoint may say, and each has to stay reachable and
        # separately testable; the seal is a further statement that even a
        # record saying only permitted things is not one a caller may mint.
        # It is **not** the defence: the three values a caller once chose are
        # properties now, so reaching this token buys nothing.
        if self.seal is not _CHECKPOINT_SEAL:
            raise TypeError(
                "Checkpoint1Handoff is factory-only: build it with "
                "checkpoint_one_handoff.  The feasibility verdict, both "
                "alternatives and the pilot cost line are derived from the "
                "components this record carries and cannot be supplied, so "
                "reaching this token gains a caller nothing; the token exists "
                "so that a record presented as a checkpoint is one the "
                "factory made")
        object.__setattr__(self, "seal", None)

    @property
    def _slowest(self):
        """The slowest measured benchmark, as every derivation here uses it."""
        return sorted(self.benchmarks,
                      key=lambda item: item.clock_steps_per_second)[0]

    @property
    def feasibility(self) -> str:
        """The verdict, **derived**: the search's, downgraded when anything blocks.

        Not a field.  A caller cannot present a friendlier verdict than the
        numerical disposition, the power estimate and the matrix's own minimum
        eligible-cell count allow, because there is nowhere to present one.
        """
        blocked = (self.disposition.is_blocking
                   or self.power.verdict != "powered"
                   or not self.matrix.meets_minimum_cells)
        if blocked and self.search.verdict == "feasible_matrix":
            return "unresolved"
        return self.search.verdict

    @property
    def alternatives(self) -> tuple:
        """Both alternatives, **derived** from the components in declared order.

        An independent review substituted a costed Option A here, and a later
        one reordered the pair in a directly built record.  Neither is a value
        this record accepts any more.  Option A is the fixed, no-input,
        unpriced outline :func:`intended_configuration_option` returns; Option
        B is computed from the matrix, the authoritative benchmark, the
        numerical disposition and the resource estimate's own safety factor,
        every one of which is bound above to what the canonical proposal
        froze.  There is no campaign here and none is claimed: the proposal
        has no campaign digest, so there would be nothing to bind one to.
        """
        return checkpoint_alternatives(self.matrix, self._slowest,
                                       self.disposition,
                                       self.resources.safety_factor)

    @property
    def pilot_line(self) -> CostLine:
        """The approved pilot's cost, **derived** from the frozen plan.

        The escape this closes: a caller built this record with Option B's
        comparison line here and :meth:`summary` reported the pilot at eleven
        hours instead of thirty-two.  There is no field to substitute now --
        the line is :func:`pilot_cost_from_plan` over the bound plan and the
        slowest measured benchmark, computed on every read.
        """
        return pilot_cost_from_plan(self.plan, self._slowest,
                                    self.resources.safety_factor)

    def alternative(self, label: str):
        """One named alternative."""
        label = _require_member(label, ALTERNATIVE_LABELS, "label")
        for item in self.alternatives:
            if item.label == label:
                return item
        raise KeyError(label)

    @property
    def digest(self) -> str:
        return _digest({
            "label": self.label, "feasibility": self.feasibility,
            "ticket_08_derivation_status": self.ticket_08_derivation_status,
            "ticket_08_authorization_status":
                self.ticket_08_authorization_status,
            "benchmarks": [item.digest for item in self.benchmarks],
            "target": self.target.digest, "power": self.power.digest,
            "recovery": self.recovery.digest,
            "disposition": self.disposition.digest,
            "search": self.search.digest, "matrix": self.matrix.digest,
            "resources": self.resources.digest,
            "firewall": self.firewall.digest, "proposal": self.proposal.digest,
            "plan": self.plan.digest,
            # Derived, so these hash what the components imply rather than
            # what a caller supplied.
            "alternatives": [item.digest for item in self.alternatives],
            "pilot_line": self.pilot_line.digest,
            "requested_decision": self.requested_decision,
        })

    def summary(self) -> dict:
        """The plain-language numbers, as a mapping, for a report to print."""
        slowest = min(item.clock_steps_per_second for item in self.benchmarks)
        peak = max(item.peak_bytes for item in self.benchmarks)
        return {
            # **Provenance, so the supplied set is not read as frozen.**  The
            # slowest benchmark is bound above to the proposal's own
            # ``benchmark_digest`` and is what every derived number reads.
            # The rest of the tuple is supplied measurement context, and the
            # two aggregates it feeds are maxima: an extra measurement can
            # only raise the peak and the empirical row width, which shrinks
            # the cube margin and widens the row bound.  A supplied benchmark
            # can therefore make this evidence more conservative and never
            # more favourable, and its digest is listed here either way.
            "authoritative_benchmark_digest": self._slowest.digest,
            "supplied_benchmark_digests": [item.digest
                                           for item in self.benchmarks],
            "walked_clock_steps_per_second": slowest,
            "peak_bytes": peak,
            "projected_cube_margin": (self.resources.clock_steps * 8.0
                                      / RACE_PASSES / peak),
            # The resource estimate's own coefficient, not the benchmark
            # maximum beside it.  An independent review found the summary
            # quoting one number while the total was built on another.
            "bytes_per_row": self.resources.bytes_per_row,
            "empirical_bytes_per_row": max(item.bytes_per_row
                                           for item in self.benchmarks),
            "power_stage": self.power.stage,
            "power_verdict": self.power.verdict,
            "arm_trials": dict(self.power.arm_trials),
            "paired_trials_per_cell": self.power.trials_per_cell,
            "physical_trial_executions": self.resources.physical_trials,
            "shadow_rewalks": self.resources.shadow_rewalks,
            "durable_rows": self.resources.rows,
            "exponent_half_width": self.power.achieved_half_width,
            "contrast_half_width": self.power.contrast_half_width,
            "pairing_inflation_measured": self.recovery.inflation,
            "contrast_inflation_measured": self.recovery.contrast_inflation,
            "numerical_verdict": self.disposition.verdict,
            "numerical_blockers": list(self.disposition.blockers),
            "probability_admissible_trials":
                self.disposition.probability_admissible_trials,
            "time_admissible": self.disposition.time_admissible,
            "search_verdict": self.search.verdict,
            "options_searched": self.search.searched,
            "matrix_wall_hours": self.resources.wall_hours,
            "matrix_storage_bytes": self.resources.storage_bytes,
            # The range pilot is deliberately outside the production total and
            # was therefore invisible at the checkpoint.  Keeping it out of the
            # total is right; leaving its number out of the decision is not.
            "pilot_wall_hours": self.pilot_line.wall_seconds / 3600.0,
            "pilot_storage_bytes": self.pilot_line.storage_bytes,
            "pilot_durable_rows": self.pilot_line.durable_rows,
            "feasibility": self.feasibility,
            # Stated rather than left to be inferred from an absence: no
            # post-range production proposal exists, because the target,
            # matrix and resource estimate a post-range record needs have not
            # been derived.  That derivation is ticket 08's, and it has not
            # begun.  This is not a statement that ticket 08 is authorized.
            "ticket_08_derivation_status": self.ticket_08_derivation_status,
            # A separate, observable fact.  Not-begun and not-authorized are
            # different statements and a reader should not have to infer the
            # second from a sentence.
            "ticket_08_authorization_status":
                self.ticket_08_authorization_status,
            "post_range_production_proposal": None,
            "proposal_status": self.proposal.status,
            "proposal_design_digest": self.proposal.design_digest,
            "proposal_digest": self.proposal.digest,
        }


def checkpoint_one_handoff(label, benchmarks, target, power, recovery,
                           disposition, search, matrix, resources, firewall,
                           proposal, plan,
                           pilot_run_names=()) -> Checkpoint1Handoff:
    """Bundle the checkpoint-1 evidence, deriving **everything** it presents.

    The verdict is the **search's** verdict, downgraded whenever the numerical
    disposition blocks, the power estimate is not powered, or the matrix cannot
    meet its own minimum eligible-clock count.  A caller cannot supply a
    friendlier one, and the friendliest one available is still a statement
    about the options that were searched.

    **Neither can a caller supply the decision.**  An independent review passed
    this boundary genuine proposal, recovery, numerical, power, matrix,
    resource and firewall components and substituted a *costed* Option A whose
    ``supports`` string was ``"a physical scaling claim"``, plus a pilot cost
    line of one second, and it was accepted -- the boundary checked only the
    types.  There are no ``alternatives`` and no ``pilot_line`` parameters any
    more.  Both alternatives are rebuilt through
    :func:`checkpoint_alternatives`: Option A is the fixed no-input unpriced
    outline, and Option B comes from the matrix, the benchmark and the
    numerical disposition.  The pilot line comes through
    :func:`pilot_cost_from_plan` from the frozen plan and the same benchmark.
    The safety factor is taken from the resource estimate rather than chosen,
    so the pilot and the matrix are priced on one convention.
    """
    require_exact_type(disposition, BudgetDisposition, "disposition")
    require_exact_type(search, SearchReport, "search")
    require_exact_type(power, PowerEstimate, "power")
    require_exact_type(matrix, PrimaryMatrix, "matrix")
    require_exact_type(recovery, RecoveryReport, "recovery")
    require_exact_type(proposal, ProposedProductionManifest, "proposal")
    measured = tuple(benchmarks)
    if not measured:
        raise ValueError("a feasibility handoff carries at least one measured "
                         "benchmark; an unmeasured feasibility claim is the "
                         "thing this ticket exists to replace")
    # **The consumption boundary.**  A proposal reaching a checkpoint is
    # rebuilt from the very components this handoff carries and compared field
    # for field; a bare self-consistent record is refused here rather than
    # presented.  The two free-text fields are taken from the record, because
    # they describe it rather than the experiment.
    require_exact_type(resources, ResourceEstimate, "resources")
    require_plan_firewall(plan, firewall)
    require_authoritative_proposal(
        proposal, label=proposal.value("proposal_label"),
        fingerprint=source_fingerprint(), matrix=matrix, target=target,
        power=power, disposition=disposition, search=search,
        benchmark=sorted(measured,
                         key=lambda item: item.clock_steps_per_second)[0],
        resources=resources, firewall=firewall,
        status_reason=proposal.value("status_reason"), recovery=recovery,
        plan=plan, pilot_run_names=pilot_run_names, status=proposal.status)
    # The feasibility verdict, both alternatives and the pilot cost line are
    # not built here any more: they are properties of the record, computed
    # from the components below on every read.  An alignment review pointed
    # out that deriving them *in the factory* still left them stored, and a
    # caller who reached the construction token could substitute one and call
    # the public ``summary()`` without presenting the record to any door.
    return Checkpoint1Handoff(
        label=label,
        # Derived: this ticket builds no post-range target, matrix or resource
        # estimate, and the proposal factory refuses by name until they exist.
        ticket_08_derivation_status=TICKET_08_DERIVATION_STATUSES[0],
        ticket_08_authorization_status=TICKET_08_AUTHORIZATION_STATUSES[0],
        benchmarks=tuple(benchmarks),
        target=target, power=power, recovery=recovery,
        disposition=disposition, search=search, matrix=matrix,
        resources=resources, firewall=firewall, proposal=proposal,
        plan=plan,
        # Derived from the one canonical constant; there is no parameter.
        requested_decision=CHECKPOINT_REQUESTED_DECISION,
        seal=_CHECKPOINT_SEAL)


def _require_authoritative_checkpoint(handoff, label, benchmarks, target,
                                      power, recovery, disposition, search,
                                      matrix, resources, firewall, proposal,
                                      plan, pilot_run_names=()):
    """Return ``handoff`` only if the factory **rebuilds it exactly**.

    The seal permits construction and proves nothing about provenance, so the
    handoff is reconstructed the way a proposal is: its own components are fed
    back through :func:`checkpoint_one_handoff` and every field is compared.

    This is what a reachable seal runs into.  An independent review's escape
    was a directly constructed record carrying Option B's comparison cost line
    where the approved pilot's belonged; the same substitution made through
    the seal, or through ``dataclasses.replace`` on a sealed record, or by
    reordering the two alternatives, or by presenting anything but the
    canonical decision text, disagrees with what these components imply and is
    refused here.  A subclass is refused by type, because a subclass reinstates
    the instance dictionary the exact type does not have and could therefore
    answer ``summary`` or ``digest`` with whatever it liked.

    **What comes back is the rebuild, never the supplied object.**  An
    alignment review poisoned a genuine handoff with an instance-attribute
    ``summary`` and this door returned it: the comparison was over fields and
    digests, and the caller's behaviour travelled straight through.  The exact
    type is slotted now, so that particular poisoning cannot be applied; the
    door returns the freshly rebuilt record anyway, so a boundary that
    authorizes something never forwards the object it was handed.
    """
    require_exact_type(handoff, Checkpoint1Handoff, "handoff")
    _require_no_undeclared_state(handoff)
    rebuilt = checkpoint_one_handoff(
        label=label, benchmarks=benchmarks, target=target, power=power,
        recovery=recovery, disposition=disposition, search=search,
        matrix=matrix, resources=resources, firewall=firewall,
        proposal=proposal, plan=plan,
        pilot_run_names=pilot_run_names)
    disagree = [name for name in Checkpoint1Handoff.__dataclass_fields__
                if name != "seal"
                and getattr(handoff, name) != getattr(rebuilt, name)]
    if disagree or handoff.digest != rebuilt.digest:
        raise ValueError(
            f"this checkpoint handoff disagrees with the one its own "
            f"components build, on {disagree or ['digest']}.  The feasibility "
            f"verdict, both alternatives, the pilot cost line, the requested "
            f"decision and the two ticket-08 states are derived by "
            f"checkpoint_one_handoff from the benchmark, the matrix, the "
            f"numerical disposition and the frozen pilot plan; "
            f"a record carrying different ones is presenting a decision "
            f"nobody's evidence supports")
    _require_no_undeclared_state(rebuilt)
    return rebuilt


def _require_no_undeclared_state(handoff):
    """Return ``handoff`` only if it carries nothing outside its declared slots.

    Structurally guaranteed by ``slots=True`` and asserted anyway, so that a
    future edit removing the slots is a test failure rather than a silent
    reopening of the shadowing route.
    """
    require_exact_type(handoff, Checkpoint1Handoff, "handoff")
    if hasattr(handoff, "__dict__"):
        raise AssertionError(
            f"a Checkpoint1Handoff carries an instance dictionary holding "
            f"{sorted(vars(handoff))}; the record is slotted precisely so "
            f"that no attribute outside its declared fields can exist, "
            f"because such an attribute shadows the method of the same name "
            f"and summary() would then be a caller's function")
    return handoff


def _require_unchanged_checkpoint(handoff, frozen_digest):
    """Return ``handoff`` only if it still hashes to ``frozen_digest``.

    There is nothing to rebuild from here -- a digest is not a set of
    components -- so this door returns the record it was given, and the
    no-undeclared-state rule above is what makes that safe.
    """
    require_exact_type(handoff, Checkpoint1Handoff, "handoff")
    _require_no_undeclared_state(handoff)
    expected = _require_hex_digest(frozen_digest, "frozen_digest")
    if handoff.digest != expected:
        raise ValueError(
            f"this checkpoint handoff hashes to {handoff.digest} and the "
            f"frozen record hashes to {expected}; a checkpoint that changed "
            f"after it was frozen is a different checkpoint")
    return handoff
