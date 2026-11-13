# kronecker-graphical-model-selection

Code for *Model Selection in Kronecker Graphical Models* by Kyung Han Lee and Min Jin Ha.

- `kgmselect/`: the installable R package.
- `replication/simulations/`: simulation drivers, settings, results and figure scripts.
- `replication/adni/`: ADNI analysis and figures for authorised inputs.

```text
replication/
├── simulations/
│   ├── scripts/
│   ├── R/
│   ├── python/
│   ├── settings/
│   ├── summaries/
│   ├── figures/
│   ├── results/
│   ├── tests/
│   └── requirements.txt
├── adni/
│   ├── scripts/
│   ├── R/
│   ├── python/
│   ├── figures/
│   ├── tests/
│   └── README.md
├── io.py
├── renv.lock
└── run_simulations.sh
```

Fresh simulation outputs are written under `replication/outputs/`.
Stored reference results remain under `replication/simulations/results/`.

Simulation replication counts every nonzero off-diagonal fitted entry as an
edge and passes `tol = 0` to the package. The separate ADNI scripts retain
their `1e-6` fitted-edge threshold. This is distinct from the partial-correlation
threshold used only to display the ADNI networks.

## Quick review

To inspect the reported simulation results without running code, open
[the table CSVs](replication/simulations/results/tables) and
[Figures 1–3](replication/simulations/results/figures).
[The table and figure mapping](#table-and-figure-mapping) identifies the
corresponding items in the main article and supplement.

For an execution check, install the package, run the small example, and then
use the quick pipeline below. Run shell commands from the repository root.

## Installation

The R dependencies for the simulations, package tests and vignette are
recorded in [renv.lock](replication/renv.lock).
Install R and the system build tools separately.

### Restore the R packages

Run this in R from the repository root:

~~~r
install.packages("renv", repos = "https://cloud.r-project.org")
replication_dir <- normalizePath("replication")
renv::restore(
    project = replication_dir,
    lockfile = file.path(replication_dir, "renv.lock"),
    library = file.path(replication_dir, ".R-library"),
    prompt = FALSE
)
~~~

The restored packages are installed in `replication/.R-library`.

In the same shell used for building and running the code:

~~~bash
export R_LIBS_USER="$PWD/replication/.R-library"
export R_LIBS_SITE="$R_LIBS_USER"
unset R_LIBS KGMSELECT_SRC
export R_ENVIRON=/dev/null R_ENVIRON_USER=/dev/null
export R_PROFILE=/dev/null R_PROFILE_USER=/dev/null
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1

R CMD build kgmselect
R CMD INSTALL --library="$R_LIBS_USER" kgmselect_0.1.4.tar.gz
Rscript -e 'print(.libPaths()); print(sapply(c("huge", "igraph", "Rcpp"), function(p) as.character(packageVersion(p))))'
~~~

Keep these shell settings for the examples and simulations so that R uses
the restored library.

Source installation on Linux needs C, C++ and Fortran compilers, make and the
system libraries required by the dependency packages.
Building the vignette also needs Pandoc. BLAS and LAPACK are supplied by the
R installation and are not installed by renv.

Simulation seed rules are defined in [generate.R](replication/simulations/R/data_generation.R).
StARS subsampling seeds are listed in [stars_seeds.csv](replication/simulations/settings/stars_seeds.csv).
Keep the settings and seed rules when comparing a rerun with stored results.

### Numerical environment

Use R 4.3.3 with huge 1.6, glasso 1.11, expm 1.0-0, rTensor 1.4.9
and igraph 2.2.3. The lockfile includes Rcpp 1.1.1.
Use one BLAS thread per worker.

The manuscript figures use Latin Modern Roman from TeX Live and Computer
Modern mathematics. ADNI dependencies are listed in
[the ADNI environment instructions](replication/adni/README.md#environment).

### Python and ADNI

Install the simulation figure dependencies separately:

~~~bash
python3 -m pip install -r replication/simulations/requirements.txt
~~~

This requirements file pins NumPy, SciPy, pandas and matplotlib.
It is not a complete operating-system or Python dependency lock.
The ADNI scripts have a separate environment and require authorised inputs.
See [their instructions](replication/adni/README.md).

The package can also be installed from GitHub with
`remotes::install_github("Kyung-Han-Lee/kronecker-graphical-model-selection", subdir = "kgmselect")`.

## Small example

Run this in R after installation:

~~~r
library(kgmselect)
truth <- list(kgm_simulate_precision(10, "er", degree = 2, seed = 1),
              kgm_simulate_precision(8, "er", degree = 2, seed = 2))
Y <- kgm_sample(60, lapply(truth, solve), seed = 3)
path <- kgm_path(Y, kgm_eta_grid(15, c(1e-2, 1e2)))
sel_ebic <- kgm_select(path, Y, criterion = "ebic", gamma = 0.5)
sel_relax <- kgm_select(path, Y, criterion = "relaxation", gamma = 0.5)
c(ebic = sel_ebic$index, relaxation = sel_relax$index)
unlist(kgm_metrics(truth, path$fits[[sel_relax$index]])[c("mcc", "fdp", "tpr")])
~~~

The expected selected indices are 15 and 8.
The Relaxation fit has MCC 0.7644869, FDP 0 and TPR 0.6470588.

An optional package test is:

~~~bash
Rscript -e 'res <- tinytest::test_package("kgmselect", verbose = 0); print(res); stopifnot(all(as.logical(res)))'
~~~

## Quick simulation pipeline

~~~bash
bash replication/run_simulations.sh --quick --workers=1 --skip=refit_pilot,refit_comparison
~~~

This runs the small matrix, tensor, cross-validation, StARS and unit-rescaling
settings, then aggregates their results and draws three figures.
It leaves out the Refit comparisons, whose time-limited fits can take longer.

Expected outputs under `replication/outputs/quick/results/`:

- Seven table CSVs: Table1, TableB1–B4 and TableB6–B7. TableB5 is skipped because
  the Refit pilot was excluded.
- Four figure-value CSVs and Figure1.pdf, Figure2.pdf and Figure3.pdf.
- manifest.json with the output hashes.

Each executed simulation step should report zero failed datasets.
A fresh run reports zero reused datasets. Repeating the same command reuses
compatible saved outputs.
The quick settings exercise the workflow; their numerical values are not the
full-study results stored under `replication/simulations/results/`.

To include the Refit stages in the small run:

~~~bash
bash replication/run_simulations.sh --quick --workers=1
~~~

## Full reproduction

Running `bash replication/run_simulations.sh --workers=1` without `--quick`
uses the full simulation settings and can take substantially longer.
Its outputs go to `replication/outputs/`.
ADNI analysis is separate and requires authorised inputs; see
[the ADNI instructions](replication/adni/README.md).

## Stored results

The reference simulation results use `tol = 0`.
The files in `replication/simulations/results/per_dataset/` contain stored selections
and summaries for regenerating the tables and simulation figures.

| File | Contents |
|---|---|
| `selected.csv.gz` | 304,900 selection rows with grid points, metrics and support counts |
| `paths.csv.gz` | 16,400 paths with counts of true-support candidates and adjacent pairs |
| `pilot_candidates.csv.gz` | 10,000 Refit pilot candidates and their statuses |
| `figure1_scores.csv.gz` | 12,000 criterion values for the candidates in Figure 1 |

The stored tensor EBIC selections cover 400 three-way BA datasets.
Stored unit-rescaling results include Scale and Relaxation at constant 1 only.
The simulation drivers compute the additional datasets and constants.
The result manifest lists the output files and their hashes.

## Table and figure mapping

CSV filenames follow the table numbers in the main manuscript and supplement.

| Result file | Manuscript item | Simulation driver |
|---|---|---|
| `tables/Table1.csv` | Main Table 1 | tensor, cv |
| `tables/TableB1.csv` | Supplement Table B.1 | matrix |
| `tables/TableB2.csv` | Supplement Table B.2 | units |
| `tables/TableB3.csv` | Supplement Table B.3 | matrix, tensor |
| `tables/TableB4.csv` | Supplement Table B.4 | matrix, cv, stars |
| `tables/TableB5.csv` | Supplement Table B.5 | refit_pilot |
| `tables/TableB6.csv` | Supplement Table B.6 | matrix, cv |
| `tables/TableB7.csv` | Supplement Table B.7 | matrix |
| `figures/Figure1.pdf` | Main Figure 1 | matrix |
| `figures/Figure2.pdf` | Main Figure 2 | analytic example drawn by figure2.py |
| `figures/Figure3.pdf` | Main Figure 3 | matrix, cv, stars |

These paths are relative to `replication/simulations/results/`.
The table CSVs contain the values used in the tables embedded in the manuscript
and supplement, retaining unrounded values and sample standard deviations.
The FDP/TPR table and Figure 3 include both StARS instability thresholds, 0.05 and 0.1.

ADNI commands for Main Figure 4, Supplement Figure C.1 and Table C.1 are
documented in [the ADNI instructions](replication/adni/README.md).
ADNI participant data and results are not distributed.

### Regenerate tables and figures from stored results

After [installation](#installation), run these commands from the repository
root to regenerate the eight table CSVs and three simulation figures:

~~~bash
mkdir -p replication/outputs/review-from-stored/per_dataset
cp replication/simulations/results/per_dataset/* replication/outputs/review-from-stored/per_dataset/
Rscript replication/simulations/scripts/summarize_results.R --from=per_dataset --results=replication/outputs/review-from-stored
python3 -m replication.simulations.figures.draw_figures --results=replication/outputs/review-from-stored
~~~

The generated files go into `replication/outputs/review-from-stored/`.
Repeating the commands replaces the generated review files and preserves
the distributed results. This aggregation uses the stored selections and
does not fit models or optimise likelihoods.

### Run an individual simulation

An experiment's driver is `replication/simulations/scripts/run_<driver>.R`.
For example:

~~~bash
Rscript replication/simulations/scripts/run_matrix.R --out=replication/outputs/matrix --workers=1
Rscript replication/simulations/scripts/summarize_results.R --outputs=replication/outputs
python3 -m replication.simulations.figures.draw_figures --results=replication/outputs/results
~~~

Aggregation reports skipped outputs when required experiment results are
absent. Driver names in the mapping also identify the corresponding steps
of `replication/run_simulations.sh`.

## Individual empirical calculations

Each of the eight calculations has its own command. The R drivers optionally
retain candidate fits and partial-covariance spectra with `--save_inputs=1`.
The default is 0. Saved input settings and seeds regenerate observations when
needed; individual simulated observations are not copied into this input bundle.

From the repository root, prepare inputs in a fresh output location:

~~~bash
Rscript replication/simulations/scripts/run_matrix.R --out=replication/outputs/with-inputs/matrix --save_inputs=1 --workers=1
Rscript replication/simulations/scripts/run_tensor.R --out=replication/outputs/with-inputs/tensor --save_inputs=1 --workers=1
Rscript replication/simulations/scripts/run_refit_comparison.R --out=replication/outputs/with-inputs/refit_comparison --save_inputs=1 --workers=1
Rscript replication/simulations/scripts/prepare_inputs.R --outputs=replication/outputs/with-inputs --out=replication/outputs/with-inputs/inputs
~~~

Alternatively, `bash replication/run_simulations.sh --save-inputs --workers=1`
retains these inputs while running the full simulation pipeline.
Use a new output location when changing `--save_inputs`.
Preparing inputs reads the retained fits and evaluates their profiles. It does
not refit the paths. The stored selected-result CSVs in
`simulations/results/per_dataset/` are sufficient for the tables and Figures
1–3; they do not contain the full candidate inputs for these eight calculations.

The prepared folder contains:

| File | Contents |
|---|---|
| `paths.csv` | Dataset IDs, experiment, dimensions, graph family, sample size, replicate and cache path |
| `scores.csv` | Candidate indices, criterion, gamma, log likelihood, EBIC and recovery metrics |
| `supports.csv` | Exact-support hashes and edge counts |
| `intervals.csv` | Candidate relaxation dimensions |
| `scales.csv` | Profiled scales, true scales and certificate bounds |
| `covariances.csv` | Partial-covariance eigenvalues by sweep and mode |
| `conditions.csv` | Fitted and relaxed condition numbers |
| `cases.json`, `profiles/` | Fixed-profile inputs for the optimisation calculation |
| `inputs.json` | Schema version and input counts |

The default `--profiles=paper` prepares the 400 two-mode 116×102 ER
(degree 2) and BA paths at sample sizes 100 and 3200.
`--profiles=all` and `--profiles=none` are also available.
The scale, covariance and condition-number commands select the paper's
equal-dimension ER settings at dimensions 40, 60 and 80 and sample sizes
100, 200, 800 and 3200. Their scope can be changed with `--dimensions`
and `--sample-sizes`.

Run the commands independently. Each output must have a new name:

~~~bash
INPUT=replication/outputs/with-inputs/inputs

python3 -m replication.simulations.scripts.check_relaxation_intervals --input "$INPUT" --output replication/outputs/intervals.json
python3 -m replication.simulations.scripts.calculate_unit_multiplier --input "$INPUT" --output replication/outputs/unit_multiplier
python3 -m replication.simulations.scripts.check_scale_certificate --input "$INPUT" --output replication/outputs/scale_certificate.json
python3 -m replication.simulations.scripts.check_partial_covariances --input "$INPUT" --output replication/outputs/partial_covariances.json
python3 -m replication.simulations.scripts.check_condition_numbers --input "$INPUT" --output replication/outputs/condition_numbers.json
python3 -m replication.simulations.scripts.check_profile_optimisation --input "$INPUT" --output replication/outputs/profile_optimisation --workers=1 --fraction=0.01 --max-boxes=20000
python3 -m replication.simulations.scripts.compare_refit_likelihoods --input "$INPUT" --output replication/outputs/refit_likelihoods
~~~

The profile calculation refines the saved fixed candidates and reports
unresolved bounds as well as selection changes. It does not fit new paths.

The Refit comparison writes `refit_jobs.csv`. Prepare and inspect the stricter
Refit convergence results separately:

~~~bash
Rscript replication/simulations/scripts/prepare_refit_convergence.R --outputs=replication/outputs/with-inputs --inputs="$INPUT" --jobs=replication/outputs/refit_likelihoods/refit_jobs.csv --out=replication/outputs/refit_diagnostics.csv
python3 -m replication.simulations.scripts.check_refit_convergence --input replication/outputs/refit_diagnostics.csv --output replication/outputs/refit_convergence.json
~~~

This step refits only the listed candidates with a maximum of 2000 iterations
and tolerance `1e-8`.

## Tests

After installation, run from the repository root:

~~~bash
Rscript replication/simulations/scripts/run_tests.R
python3 -m unittest discover -s replication/simulations/tests/python -q
python3 -m unittest discover -s replication/simulations/figures/tests -q
~~~

These tests use small synthetic inputs and cover input preparation.
ADNI tests are documented [separately](replication/adni/README.md).

## Documentation

- [Worked example](kgmselect/vignettes/kgmselect.Rmd).
- Function help: `?kgm_fit`, `?kgm_score` and the other exported function names.
