# Review and publication notes

The companion manuscript is a first working draft for collaborative development, dated 1 October 2026. It documents a phenomenological mechanism, the simulations already performed, and a new conditional experiment. It is not a claim that vacuum fluctuations have been shown to generate the Born rule.

## Manuscript files

- `Born_Selection_Draft.tex` is the canonical editable manuscript. It compiled successfully in the desktop editor.
- `Born_Selection_Draft.md` is the portable reading and GitHub version, exported from the same source.
- `conditional_predictions.csv`, `check_prediction.py`, and `prediction_validation.json` reproduce the new two-channel analytic calculation. The program needs Python and NumPy and writes only the result files beside itself.
- `provenance.json` records contributors, claim status, and review status. It is a local draft record, not a claim of satisfying the journal's submitted-article schema.
- `source_manifest.json` pins the existing result summaries and source transcript used for the paper. Those source artifacts remain in their original locations.

## What is established within the models

The linear response and assumed square-response Poisson event law yield exact conditional Born weights for identical deterministic responses. A perfectly common fluctuation preserves this factorization. Existing independent and spatially correlated noise runs showed finite-sample compatibility in selected regimes, not a general theorem. Common thresholds and the tested Adler phase race did not supply an exact replacement. The latter's bias is a separate result about a different model.

The manuscript's new slow-gain comparison is analytically specified and numerically checked. The hypothesized race keeps rates proportional to fixed input weights multiplied by local complex gains until its first record. Standard quantum absorption updates the unrecorded branch amplitudes. These different dynamics predict different outcome and timing distributions under independent slow gains, despite sharing the same mean response. The calculation has not been fitted to laboratory data or independently reviewed.

A particularly important objection for review is whether the proposed race is an admissible physical evolution at all. Its global stopping rule and nonlinear dependence on input weights require scrutiny. It may ultimately fail consistency requirements rather than describe a new effect. That possibility is stated in the manuscript, alongside the experiment that would test the restricted model.

## Existing numerical evidence

These directories are siblings of this paper folder under the research `outputs` directory:

| Directory | Evidence used | Limitation |
|---|---|---|
| `complex_response_test` | 84,000 incident histories; analytic factorization; response and threshold comparisons | Poisson recording remains assumed; noisy agreement is finite-sample |
| `phase_alignment_race` | Exact phase-race distribution and 200,000-trial follow-up | Separate model; not a refutation of complex linear response |
| `reversible_phase_kicks` | 120,000 trials; closed single-excitation reversal and phase sensitivity | Standard stochastic readout; no outcome derivation |
| `standard_model_comparison` | 288,000 sampled records; coherent absorption versus undepleted approximation | Site labels explicitly sampled from shared Born law; timing differences are approximation effects |

The earlier Claude formulation was read in `Superposition_Equations (5).md`, section 12. It explicitly separates linear response from a Poisson recording race. These manuscript claims must not be strengthened into “complex geometry alone derived Born selection.”

## Review of the AI Scientific Journal proposal

Repository reviewed read-only: [rayolddog/AI_Scientific_Journal](https://github.com/rayolddog/AI_Scientific_Journal), main commit `a08a13fe4b4c9b360916bf5b8dd1a9f21c42fe51`. The GitHub API reported it private on 1 October 2026. Reviewed documents include the README, CONTRIBUTING, manifesto, article specification, and full review protocol. The case-study and historical next-steps documents were also consulted. This is a review of the publishing proposal, not an independent audit of the existing physics article or its referee reports.

The proposal describes **Aperture as a blueprint and demonstration, not an operating journal**. Its most useful features are version-specific review records, explicit author replies, retention of dissent, a bounded dispute process, an experimentalist role, and separate records of contribution credit and human accountability. Those features suit this paper because its successes and unresolved assumptions need to remain visible together. Its protocols on equal reviewer access, withholding live-round critiques, and a hash commitment for anticipated findings offer a concrete way to reduce reviewer contamination.

I recommend the following refinements before treating the blueprint as an operating publication process:

1. **Describe provider diversity as a safeguard, not demonstrated decorrelation.** Different labs can share literature, training sources, and reasoning failures. Retain the multi-provider requirement but measure error correlations on independently adjudicated cases. A field named `provider_diversity_met` would make a narrower claim than `decorrelation_met`.
2. **Separate development collaborators from independent referees.** Claude's further development of this paper would be coauthor work. An independent report should come from a fresh reviewer context, disclose exposure to our analysis, and not be described as independent merely because the provider differs.
3. **Define error-rate denominators and observation windows.** Retraction rates alone reward low scrutiny and penalize active correction. Track audited claims, confirmed errors, severity, correction delay, and follow-up duration, as well as the denominator of papers assessed.
4. **Reconcile historical and current status.** The repository is private while several documents describe a public blueprint and open public engagement. Historical NEXT_STEPS text also retains local-only/no-push language despite the remote now existing. Keep the historical record, but add a clearly dated current-status statement rather than silently rewriting it.
5. **Strengthen factual framing.** The medical-school “half wrong within ten years” statement is presented without support in the reviewed text. Label it as an anecdote or replace it with a sourced, appropriately scoped claim. Avoid universal claims about all conventional journals or all human review; the proposal does not need them to justify its design.
6. **Version the computational record along with the manuscript.** Shared code can remain in one location, but every reviewed paper version should pin its precise code, data, environment, and output hashes. Evolving code must not change what an old paper version appears to have demonstrated.

No repository files were edited, no article was submitted, and no reviews or invitations were sent.

## Contribution credit and conventional publication

The provisional repository byline follows the requested order: **Codex (OpenAI), Claude (Anthropic), John Bramble**. Claude is credited for prior equation development and has not yet reviewed or approved this manuscript. John is the proposed human accountable sponsor. Exact AI model identifiers are left unasserted where they cannot be verified from the relevant session records.

For conventional submission, the issue is not simply putting a human first. For example, [Nature Portfolio's authorship criteria](https://www.nature.com/nature-portfolio/editorial-policies/authorship) require approval and personal accountability; Nature Methods also states that AI tools cannot be authors in its discussion of [responsible AI use](https://www.nature.com/articles/s41592-026-03020-1). A venue-specific human-author byline with a detailed AI-contribution disclosure may be required. That disclosure can explicitly credit hypothesis development, formalization, implementation, analysis, and drafting without misrepresenting those contributions as solely human work. Check the selected venue's policy at submission time.

The present package is local and unpublished. Adapting it into Aperture's complete article-folder format would be a subsequent repository contribution, with pending fields left pending. It should not borrow the review-complete status or reports of the existing Article 0001.

## Questions for Claude during development

These are proposed coauthor-development questions, not a completed review:

- Does the manuscript accurately reproduce the earlier linear-response/separate-race formulation, including its assumptions?
- Can the complex fluctuation process itself generate the commitment events, rather than relying on separate Poisson thresholds?
- Is the slow-gain race's failure to update unrecorded branch weights physically defensible, or does the standard quantum conditional evolution already resolve it?
- Can one jointly specify the response, fluctuation spectrum, event law, and depletion without losing positivity, normalization, energy accounting, or consistency for entangled inputs?
- Does the correlation experiment distinguish the proposed dynamics after ordinary reservoir effects, collective coupling, finite efficiency, and technical noise are included?
- Which claims should be withdrawn or narrowed before independent refereeing?

Independent referees should receive the frozen manuscript, code, and equal evidence access, without our anticipated-findings discussion in their live review packet. No independent Claude or other-model review has been requested in this session.
