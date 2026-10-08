# Scope, selection decisions, and provenance

## Inclusion rule

Material was included when it directly addresses at least one of the following:

- complex-amplitude propagation and interference at a detector;
- microscopic absorption, reversible capture, or first irreversible registration;
- dynamical outcome selection, exclusivity, or energy routing;
- attempts to obtain Born-weighted records without drawing outcomes from the Born distribution;
- negative results that identify where Born weights were implicitly assumed or where a proposed mechanism fails.

Files were copied from `/Users/john-bramble/Projects/Physics/DiracKuramotoFramework/`; the source repository was not modified. The pre-existing contents of this `BornSelection` folder were retained.

## Included articles and notes

### Core Born-selection paper

- `papers/born-selection/Born_Selection.pdf` — latest named compiled paper in the source project.
- `papers/born-selection/PAPER1_DRAFT_born_selection.md` — editable working manuscript and its explicit revision ledger.
- `papers/born-selection/PAPER1_OUTLINE_born_selection.md`, `PRIOR_WORK_born_selection.md`, `BORN_COUPLING_SUBSECTION_DRAFT.md`, and the two walkthrough/technical notes — supporting structure and derivations.
- `papers/born-selection/NEGATIVE_RESULT.md` — essential negative result documenting that synchronization dynamics did not by themselves reproduce Born basin measures.
- `papers/born-selection/single_photon_derivation.md` — detector-specific attempt connecting photon capture to avalanche registration.

### Closely coupled measurement paper

- `papers/measurement-cut/DIRAC_HEISENBERG_CUT_PAPER.md` and `Heisenberg_Cut_v0.4.pdf` — included because the reversible-to-irreversible detector boundary, registration, and exclusivity are inseparable from the microscopic Born-selection question.
- `papers/measurement-cut/measurement_paper_outline.md` — retained as the conceptual map for that companion work.

The measurement-cut material is a companion track, not evidence that the Born-selection mechanism has been derived.

## Included simulations

- `models/born_selection_sims/` — early outcome-selection tests and controls.
- `models/adler_two_channel_exploratory/` — the two-channel locking race, including its recorded negative and diagnostic results.
- `models/adler_born_two_channel/` — the later, carefully gated stochastic/analytic package. Bulky generated `results/` and bytecode caches were not copied.
- `models/first_mark_two_absorber/` — microscopic two-absorber unitary baseline. It models reversible capture and explicitly contains no collapse or outcome rule.
- `models/heisenberg_cut_recoverability/` — tests of the reversible/irreversible registration boundary.
- `models/legacy_code/` — compact Born-selection scripts and plots used by the manuscript.
- `models/resolution/dispute-001-*` — the specific resolution loop concerning whether basin measures reproduce Born weights. Unrelated interpretive disputes were excluded.

## Included sessions and reviews

- `sessions/2026-07-20-born-gamblers-ruin-outcome-selection.md` — directly focused on Born selection.
- `sessions/2026-07-16-t5-stage4-and-selection-background-reconciliation.md` — included only because it explicitly re-homes an effect to the selection background; its sidereal content is contextual rather than a new Born-selection result.
- `reviews/born_selection_AI_review_2026-07/` — the preserved multi-model review record, including the internal Claude/Fable review with its stated down-weighting.

The existing `PROJECT_CONTEXT.txt` records three additional ChatGPT/Codex conversation identifiers and explains what was recovered from them. The current repository contains their durable outputs and summaries, not a fabricated transcript export.

## Deliberately excluded from BornSelection

The following source-project article families were reviewed and excluded because their principal question is not microscopic detector selection:

- `PAPER_UNIFIED.md` / the broad Many Clocks interpretation;
- `DISCRETIZATION_AS_SYNC_PAPER.md`;
- `EMERGENT_FIELDS_PAPER.md`;
- `COSMIC_EXPANSION_PAPER.md`;
- `AB_VISIBILITY_PAPER.md`;
- `SIDEREAL_DECOHERENCE_PAPER.md` and `LIGO_SIDEREAL_TEST_T5.md`;
- general submission, cover-letter, authorship, and repository-administration files;
- the `ZBW_synchronization`, `limiting_curvature`, hologram, cosmology, and general reference collections;
- large generated caches, logs, and raw result trees that can be regenerated or remain available in the source repository.

## Recommended separate repository: complexified 4D and 3D + spacetime

The following files form a coherent separate project and were intentionally not mixed into BornSelection:

- `4D-3DandSpacetime.md` and its PDF;
- `4D-3D+spacetime-parametors-and-experiments.md`;
- dated/numbered backups of the 4D–3D manuscript, retained only as history if desired;
- any future focused extraction of light-cone discussion from broad framework sessions.

Suggested working title: `Complex4D-Spacetime`.

That future repository should ask how an unoriented or complexified four-dimensional substrate relates to emergent 3+1 spacetime, causal ordering, and the speed-limited light cone. BornSelection should import a result from that project only when it becomes a concrete constraint on detector interaction or exclusivity.

## Provenance snapshot

- Source repository: `/Users/john-bramble/Projects/Physics/DiracKuramotoFramework/`
- Source Git remote at collection time: `https://github.com/rayolddog/DiracKuramotoFramework.git`
- Collection date: 2026-10-08
- Existing BornSelection context assembled: 2026-10-03

No source files were moved or deleted.
