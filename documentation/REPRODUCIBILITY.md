# Reproducibility notes

**Environment:** Base R. No external R packages are called by the supplied script. R version was not specified by the author; record `sessionInfo()` after a successful run.

**Working directory:** Start RStudio in the repository root. The portable version removes the original absolute Windows `setwd()` command and writes all outputs to a relative folder.

**Checks performed by script:** finite-difference derivative checks, interior equilibrium conditions, Jacobian trace/determinant/eigenvalues, reduced-system trajectories and status logs, parameter-grid sensitivity, direct versus indirect natural-capital production channel, numerical Dulac divergence sign check, and exogenous demographic scenarios.

**Important distinction:** The analytic existence/stability statements apply to the autonomous constant-saving reduced system with fixed population. Exogenous demographic paths do not establish an autonomous stability theorem.

**Numerical event detection:** RK4 with `dt=0.05` checks stages and endpoints; an infeasibility stop indicates the first detected problematic grid step, not the exact crossing time.

**Verification status:** The uploaded script has been inspected, and the five author-supplied images have been copied byte-for-byte. R execution and pixel-by-pixel comparison with regenerated figures have not been performed in this environment.
