# BornSelection

Research repository for a microscopic account of how a single quantum excitation interacts with a detector and produces one recorded outcome.

## Guiding constraint

The project does **not** sample repeated microscopic outcomes from the Born distribution. Complex amplitudes remain complex throughout propagation and interference; squared magnitude is used where the geometry of complex amplitudes requires it, especially in interference and conserved intensity. The research question is whether detector dynamics, microscopic ready-state variation, and a physically justified commitment mechanism can generate individual records whose aggregate frequencies agree with Born weights without inserting those frequencies as a sampling rule.

This distinction is the repository's central honesty test. Several preserved models are negative controls or conditional constructions: they reproduce Born weights only after assuming a linear hazard, martingale, threshold distribution, or other structure equivalent to the desired measure. Such results are useful, but they are not microscopic derivations.

## Repository map

- `papers/born-selection/` — the Born-selection manuscript, working source, and supporting notes.
- `papers/measurement-cut/` — closely related detector-registration and Heisenberg-cut work.
- `models/` — simulation families and their documentation. Generated caches and bulky run directories are intentionally excluded.
- `studies/` — current BornSelection studies already assembled in this project.
- `sessions/` — preserved research-session records relevant to outcome selection.
- `reviews/` — AI review record for the Born-selection manuscript.
- `reference/` — earlier Claude notes, experiments, and manuscript material already collected here.
- `docs/SCOPE_AND_PROVENANCE.md` — inclusion decisions, exclusions, provenance, and the boundary with the proposed spacetime repository.

## Current scientific position

The present reference model is complex linear detector response followed by a stochastic first-event race. It reproduces Born proportions conditional on an event when the event hazard is proportional to local response intensity. That hazard law and physical exclusivity are assumptions still requiring a detector-level derivation. Equal-threshold and several deterministic locking races are preserved because they fail in informative ways.

The immediate target is therefore not another fit to the Born curve. It is a semiconductor- or other detector-specific microscopic account of:

1. coherent capture and interference;
2. the detector ready-state ensemble;
3. the first irreversible commitment event;
4. one-record exclusivity and energy routing; and
5. why the commitment law has the needed dependence on complex amplitude without presupposing an outcome distribution.

## Running the current studies

The existing study scripts use Python with NumPy, SciPy, and Matplotlib. See `PROJECT_CONTEXT.txt` before extending them. Individual imported model directories contain their own status and run instructions; their caveats are part of the scientific record.

## Related but separate topic

The emergence or interpretation of a light cone from complexified four-dimensional space versus three-dimensional space plus spacetime is not included here as a primary research track. It may constrain a future microscopic model, but it deserves a separate repository. See `docs/SCOPE_AND_PROVENANCE.md` for the files identified for that future project.
