# Growth on a Finite Surface: Land Competition, Subsistence Requirements, and Regenerative Natural Capital

## Reproducibility materials

This repository accompanies the theoretical manuscript **“Growth on a Finite Surface: Land Competition, Subsistence Requirements, and Regenerative Natural Capital”** (prepared for submission to *Ecological Economics*).

The model allocates a fixed land endowment among agricultural production, housing, non-agricultural market production, and conservation. Functional natural capital is a stock distinct from conservation land. The numerical exercises illustrate the reduced constant-saving dynamical system and an explicitly exogenous demographic extension; they **do not** solve the infinite-horizon planner's optimization problem.

**Author:** Guamán-Lazo, Luis - Universidad Técnica de Manabí.

## Repository contents

- `code/Reproducible_model_code.R`: portable Base R script. This is the original submitted computational logic with only the machine-specific `setwd()` removed and the run instruction corrected.
- `code/Reproducible_model_code_original.R`: unmodified original script supplied by the author (contains an absolute Windows path).
- `figures/`: five author-supplied PNG images corresponding to Figures 1–5 of the manuscript. They were not regenerated or altered in preparing this repository.
- `documentation/FIGURES.md`: figure-to-code mapping and interpretation.
- `documentation/REPRODUCIBILITY.md`: execution and numerical scope.
- `finite_surface_outputs/`: created when the script is run; intentionally not versioned.

## Quick start (RStudio Desktop / Windows)

1. Clone or download the repository and open its top-level folder as an RStudio project or set the RStudio working directory to that folder.
2. Run the following command in the R console:

```r
source("code/Reproducible_model_code.R")
```

The script uses **Base R only**, with no additional package installation required. It creates `finite_surface_outputs/` containing PNG, TIFF (600 dpi), PDF graphics and CSV numerical checks. The five supplied images in `figures/` are the author-verified manuscript copies; rerunning the script creates output files with matching descriptive stems in `finite_surface_outputs/`.

## Baseline parameters

The numerical experiments use normalized illustrative values, not empirical estimates. Core values include total land `Tbar=1`, agricultural ecological elasticity `xi=0.50`, production elasticities `alpha=0.35`, `beta=0.25`, `gamma=0.20`, saving rate `s=0.22`, depreciation `delta=0.05`, ecological regeneration parameter `rhoN=0.08`, ecological carrying capacity `Nbar=1`, and ecological floor `Ncrit=0.20`. The full parameter list is declared at the start of the R script.

## Scope and limitations

- Existence, uniqueness, and local stability calculations concern the **reduced constant-saving system** with fixed population and conservation allocation.
- The demographic exercise allows exogenous population growth; the fixed-population local-stability theorem does not automatically apply to that nonautonomous system.
- RK4 trajectories use a fixed time step `dt=0.05`; stopping times are the last feasible grid points, not exact analytically located boundary-hitting times.
- The numerical Dulac grid is a diagnostic and not a substitute for the analytic Bendixson–Dulac argument; exclusion of periodic orbits does not prove global convergence.
- Figure files supplied by the author have been preserved exactly. Their equivalence to fresh outputs has not been independently tested in this preparation environment.

## Citation and versioning

For a submitted or published article, cite a **tagged release** (for example `v1.0.0`) rather than the moving `main` branch. An archival DOI can subsequently be minted by linking a GitHub release to Zenodo. The manuscript's permanent repository URL and DOI should be inserted only after publication of the release.

## License

The code is distributed under the MIT License; see `LICENSE`. The manuscript figures remain copyright of their respective rights holder(s) and are provided for scholarly inspection and reproduction of this study; see `figures/README.md`. Reuse of figures may require separate permission.
