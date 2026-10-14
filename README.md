# Graph-Fused Multinomial (GFM)

This project estimates category probabilities from multinomial counts observed
on a graph. It includes the GFM estimator, graph/simulation utilities,
within-node count-splitting cross-validation, illustration scripts, Chicago
crime analysis, and Census of India mother-tongue preparation. The supplied
real-data caches, summary tables, and historical figures are retained.

At node `i`, `Y[i, ]` contains raw counts and `n_i = rowSums(Y)[i]`. The
estimator uses `D - 1` logits relative to the final category and minimizes the
**total** multinomial negative log-likelihood plus
`lambda * sum(abs(D_G %*% beta))`. The likelihood is not normalized by trial
count. The same lambda can behave differently when counts, graph degrees,
or category counts change.

See [data/README.md](data/README.md) for data sources, schemas, and verified inventories.

## Folder structure

| Location | Purpose |
| --- | --- |
| `R/` | Reusable model, graph, evaluation, plotting, and preparation functions |
| `simulations/` | Grid, kNN, SBM, and oracle-versus-CV illustrations |
| `applications/` | Chicago comparison and India preparation entry points |
| `scripts/` | Dependency installation, grid signal figure, and summary plots |
| `data/` | Original real-data caches and their descriptions |
| `results/` | Supplied summary CSVs, unchanged |
| `results/generated/<workflow>/` | New analysis objects, metrics, and session information |
| `figures/Comparison/`, `figures/Real Data Analysis/` | Historical figures, unchanged |
| `figures/generated/<workflow>/` | Newly produced figures |
| `tests/` | Numerical, preparation, and integration checks; synthetic fixtures |
| `archive/` | Historical code, pre-review reference scripts, and unrelated coursework |

## Installation

```bash
Rscript scripts/install_dependencies.R core plots data
# Optional MCAR comparisons:
Rscript scripts/install_dependencies.R mcar
# Optional oracle-table LaTeX export:
Rscript scripts/install_dependencies.R latex
```

| Group | Packages |
| --- | --- |
| `core` | Matrix, igraph, FNN |
| `plots` | ggplot2, patchwork |
| `data` | sf, readxl |
| `mcar` | CARBayes |
| `latex` | xtable |

The installer installs only missing packages. Active analysis scripts check
packages but never install or attach them. Unused dependencies on reshape2,
doParallel, foreach, dplyr, tidyr, viridis, ggh4x, and httr2 were removed from
active workflows. `sf` may need system libraries when installed from source;
follow its installation instructions for your R platform.

Run the examples below from the project root. Entry scripts also accept an
absolute script path from another working directory; no personal `setwd()`
is needed. R 4.3.3 was used for verification.

## Quick start

```r
source("R/load_all.R")
set.seed(1)
g <- build_grid_graph(5, 5)
D_G <- build_incidence_matrix(g)
truth <- generate_true_beta_X(g, D = 4, nrow_grid = 5, ncol_grid = 5)
Y <- simulate_data(truth$P_true, n_i = rep(100, igraph::vcount(g)))
cv <- cv_graph_fused_multinomial(
  Y, D_G, lambda_grid = c(0.5, 1, 2.5, 5), K = 5, seed = 1,
  T_max = 8000, verbose = FALSE
)
fit <- fit_graph_fused_multinomial(
  Y, D_G, lambda = cv$lambda_min, T_max = 8000, verbose = FALSE
)
fit$P_hat
fit$converged
fit$residual
```

`load_all.R` defines core and preparation functions. It requires only Matrix,
igraph, and FNN up front; optional data/plot packages are checked when needed.
Package namespaces are explicit.

## Execution order and workflows

Start with the validation suites and small operational profiles:

```bash
Rscript tests/run_all.R
Rscript simulations/01_grid_x_pattern.R --quick
Rscript simulations/02_knn_graph.R --quick
Rscript simulations/03_sbm_graph.R --quick
Rscript simulations/05_oracle_vs_cv.R --quick
Rscript scripts/grid_signal.R --quick --no-mcar
Rscript scripts/plot_summary_curves.R
```

Quick simulation profiles use smaller graphs, `tol = 1e-4`, and fewer lambda
candidates. They check execution, not reproduce the research experiments.
Full profiles retain the supplied default graph sizes and `tol = 1e-6`.

| Script | Full default | Interpretation |
| --- | --- | --- |
| `simulations/01_grid_x_pattern.R` | 50 by 50 grid, D = 4 | Rook CV and oracle sweep; rook/king fits at the same selected lambda |
| `simulations/02_knn_graph.R` | 2,000 locations, k = 5 and 15, D = 4 | Lambda selected on smaller-k graph; both graphs fitted to identical counts |
| `simulations/03_sbm_graph.R` | Five blocks: 300, 250, 400, 350, 300; D = 4 | Noisy block signal; CV fit versus lambda = 0 |
| `simulations/05_oracle_vs_cv.R` | Target n = 100, 500, 1,000, 2,000; D = 4 | One dataset per target; truth-selected versus CV-selected lambda |
| `scripts/grid_signal.R` | 20 by 20 grid, D = 4, lambda = 2.5 | Truth/MLE/GFM/MCAR panels and calculated relative errors |

Run full versions by omitting `--quick`. The signal script requires CARBayes
unless `--no-mcar` is supplied. These demonstrations are not the complete
20-replication paper study; that study's full generators and replicate logs
are not present. The kNN demonstration does not jointly select k and lambda.
The SBM truth includes Gaussian logit perturbations rather than exact block
constancy. The oracle comparison is an illustration, not a CV consistency proof.

Sourcing a workflow defines its function without starting a fit:

```r
source("simulations/02_knn_graph.R")
result <- run_knn_demo(n_pts = 100, k_small = 5, k_large = 10, seed = 7)
```

This replaces the former immediate execution on sourcing. Use `Rscript` or
call the named function explicitly. Outputs are written to workflow-specific
`results/generated/` and `figures/generated/` folders.

### Chicago

```bash
Rscript applications/chicago_comparison.R --no-mcar
Rscript applications/chicago_comparison.R --largest-component --no-mcar
# Full MCAR comparison once CARBayes is installed:
Rscript applications/chicago_comparison.R --largest-component
```

The cached full dataset has 652 cells, two components, and no isolated nodes.
The smaller component has two cells and four training records, with several
absent categories. A finite unconstrained logit minimizer need not exist on
this domain. Reviewed full-graph GFM runs at 5,000 and 20,000 iterations did
not meet `tol = 1e-6`. The workflow reports this; metrics with `converged = FALSE`
describe the final iterate rather than a certified optimum.

`--largest-component` explicitly retains 650 cells and excludes five full-data
records. `chicago_largest_component()` keeps the existing holdout assignments
for retained cells and records exclusions. It does not add neighbors or
replace the original cache. Use the same selected object for both fits.
Component selection does not by itself guarantee convergence at a chosen
iteration budget; inspect `converged`, `residual`, and `n_iter`.

```r
source("applications/chicago_comparison.R")
result <- run_chicago_comparison(
  include_mcar = FALSE, largest_component = TRUE,
  lambda = 1, T_max = 20000, tol = 1e-6
)
result$fit$converged
result$fit$residual
```

Figure labels use scores calculated from their own counts and predictions.
Historical MSE and timing constants are never inserted into computed outputs.

### India

```bash
Rscript applications/india_preparation.R
Rscript applications/india_preparation.R --largest-component
```

Preparation checks the 36 cached workbooks and district boundaries, constructs
counts and a graph, and saves `prepared_india.rds` under `results/generated/`.
It does not fit either statistical model.

```r
source("R/load_all.R")
india_full <- prepare_india_mother_tongue(
  cache_dir = "data/india_mother_tongue_cache"
)
india <- india_largest_component(india_full)  # Optional common MCAR domain.
# Choose lambda for this count scale before model fitting:
# fit <- fit_graph_fused_multinomial(india$Y_train, india$D_G, lambda = ...)
```

The full graph has nine isolated districts. Its largest component retains 625
districts, with all 15 exclusions recorded. The MCAR wrapper checks positive
neighbor totals, count/graph order, and the original intercept-only multinomial
specification. See the [CARBayes reference manual](https://cran.r-project.org/web/packages/CARBayes/refman/CARBayes.html)
for fitted-value and MCMC conventions.

## Model, score, and CV conventions

- `Y`, `Y_train`, `Y_test` are node-by-category raw counts. OTHER is last for
  real data. Pseudocount `a` affects GFM initialization only.
- `D_G` is sparse, edge-by-node, with one row per undirected edge. Orientation
  does not change the absolute-difference penalty. `W` follows the same nodes.
- Simulation MSE is `mean((P_hat - P_true)^2)`; relative error is
  `sqrt(mean((P_hat - P_true)^2) / mean(P_true^2))`.
- Chicago test MSE averages squared differences from raw test proportions
  over nonempty test cells and categories. Empty test cells are grey in all
  map panels and omitted from this score.
- Individual log loss and multicategory Brier scores weight nodes by test
  count. Brier is not MSE between node proportions.
- Named probabilities are aligned explicitly; unnamed matrices must already
  follow the count order. MCAR fitted counts are normalized by their row totals.

CV partitions trials within nodes and retains the full graph. It averages
fold negative log-likelihood per trial equally over nonempty folds. Its SE is
fold-rate SD divided by the square root of active-fold count. The one-SE rule
is a heuristic; overlapping training folds do not yield a confidence interval.
Lambda is not automatically rescaled at the full-data refit. Small iterate
changes are a stopping diagnostic, not an optimality-gap certificate.

## Imported summaries and original figures

Summary plots read CSV axes, means, and `_sig` fields without changing them.
They produce distinct n, D, joint-growth, and timing figures. Error axes use
"Reported error" because complete metric/provenance metadata are missing.
The default interval multiplier is one; missing replicate logs prevent the
code from establishing whether `_sig` represents SD or SE.

The old plot script overwrote x values, transformed uncertainty as though it
were a mean, divided widths by 40, multiplied others by six, and overwrote an
error plot with a timing plot. These operations were removed. Timing is now
shown in seconds on a logarithmic axis. Original CSVs and historical PDFs are
preserved and are not presented as newly reproduced results.

## Reproducibility, tests, and review scope

Workflows save settings, analysis objects, metrics, measured elapsed time, and
`sessionInfo()`. Seeds for splitting, simulations, and MCAR are local to those
operations. Grid and oracle illustrations now have explicit seeds; the supplied
versions were unseeded. Synthetic test workbooks are labeled and are not Census
observations.

`tests/run_all.R` combines 13 numerical, 15 preparation, and seven project
integration groups. They cover analytic solver cases, extreme logits, graphs,
RNG/split conservation, workbook/shapefile fixtures, independent score
calculations, cached Chicago alignment, unchanged summaries, and project paths.

The review also parsed all 35 state workbooks and the national workbook,
constructed the full India graph, ran small simulation workflows, rendered
figures, and exercised GFM on cached Chicago counts. No new live-data download,
full MCAR MCMC run, or full historical Monte Carlo benchmark was performed.
MCMC diagnostics and reproduction of paper tables require a separate run.

Verification versions: R 4.3.3; Matrix 1.6.5; igraph 1.6.0; FNN 1.1.4;
sf 1.0.15; readxl 1.4.3; ggplot2 3.4.4; patchwork 1.2.0. This is a record of
tested versions, not a lockfile for all transitive dependencies.

Existing public estimator/preparation signatures and the statistical objective
are retained. This review added reusable path/score/plot helpers, callable
workflows, explicit output locations, optional MCAR checks, and a holdout-
preserving Chicago selection helper. Preparation files moved into `R/`, and
analysis scripts moved out of data/figure output folders. Historical scripts
and unrelated coursework remain in `archive/` outside the active execution path.
