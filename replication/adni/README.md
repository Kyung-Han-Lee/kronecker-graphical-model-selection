# ADNI analysis

Analysis and figures for *Model Selection in Kronecker Graphical Models*.
Run all commands below from the repository root.

## Files

- `scripts/`: commands for analysis, bootstrap, projection, figures and table export.
- `R/`: fitting, likelihood, selection, relaxation, Refit and bootstrap functions.
- `python/`: passband projection, preprocessing settings and bootstrap table data.
- `figures/`: criterion curves, connectome rendering and figure styles.
- `tests/`: synthetic R and Python tests.

The ADNI numerical helpers retain their fitted-edge tolerance `1e-6`.
The network display threshold is a separate absolute partial correlation of 0.2.

## Inputs

The starting inputs are harmonised ROI tensors and covariates after the image
preprocessing and ComBat steps described in Appendix C. Participant data and
individual analysis results are not distributed.

Obtain ADNI data through the [LONI Image and Data Archive](https://ida.loni.usc.edu)
under the applicable access approval and Data Use Agreement.

Set `ADNI_INPUT` to a local aggregated-input folder with this structure:

~~~text
aggregated/C/
├── combat/
│   ├── Y_CN_combat.npy
│   ├── Y_MCI_combat.npy
│   ├── Y_Dementia_combat.npy
│   ├── covars_CN.csv
│   ├── covars_MCI.csv
│   └── covars_Dementia.csv
└── combat_pb102/
    ├── Y_CN_combat.npy
    ├── Y_MCI_combat.npy
    ├── Y_Dementia_combat.npy
    └── resample_manifest.json
~~~

The `combat/` arrays have shape (participants, 116, 187). Covariate rows must
follow the same participant order and contain Age, Sex and APOE4 without missing
values. The manuscript uses 289 CN, 158 MCI and 73 Dementia participants.

The projection command below creates `combat_pb102/`, whose arrays have shape
(participants, 116, 102). An existing projection can be used with its matching
`resample_manifest.json`. The R analysis checks this manifest before fitting.

## Environment

Restore the simulation R dependencies as described in the
[root README](../../README.md#installation), then install reticulate 1.46.0
into the same R library. The analysis also uses huge 1.6, glasso 1.11,
expm 1.0-0, rTensor 1.4.9, RhpcBLASctl 0.23-42, jsonlite 2.0.0 and digest.
The numerical environment uses R 4.3.3 and one BLAS thread per worker.

~~~bash
python3 -m pip install -r replication/adni/requirements.txt
export RETICULATE_PYTHON="$(command -v python3)"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
~~~

The Python environment uses Python 3.12, NumPy 2.3.5, SciPy 1.15.3,
pandas 2.3.3, matplotlib 3.10.8, nibabel 5.4.2, PyVista 0.48.4 and
Pillow 12.1.1. Renderer checks use a working headless VTK implementation
(VTK OSMesa 9.3.1 in the reference environment).
Use Latin Modern Roman from TeX Live for the manuscript figure typography.

## Analysis

Set `ADNI_INPUT` to the absolute input folder and choose a fresh output folder:

~~~bash
ADNI_INPUT=/absolute/path/to/aggregated/C
ADNI_OUTPUT=replication/outputs/adni

python3 -m replication.adni.scripts.project_passband --input "$ADNI_INPUT/combat" --output "$ADNI_INPUT/combat_pb102"

Rscript replication/adni/scripts/run_analysis.R --aggregated "$ADNI_INPUT" --groups CN,MCI,Dementia --arms projected --n-eta 40 --out "$ADNI_OUTPUT/joint_profile"

Rscript replication/adni/scripts/save_selected_fits.R --aggregated "$ADNI_INPUT" --groups CN,MCI,Dementia --arm projected --gamma 0.5 --verdict "$ADNI_OUTPUT/joint_profile/probe_verdict.csv" --out "$ADNI_OUTPUT/selected_fits"

Rscript replication/adni/scripts/run_bootstrap.R --aggregated "$ADNI_INPUT" --groups CN,MCI,Dementia --arm projected --B 500 --seed 2026 --workers 1 --pairs reported --refit 0 --out "$ADNI_OUTPUT/bootstrap"

Rscript replication/adni/scripts/check_covariance.R --aggregated "$ADNI_INPUT" --groups CN,MCI,Dementia --arms projected --n-eta 40 --out-csv "$ADNI_OUTPUT/diagnostics/assumption_check.csv"
~~~

Each command runs separately. The bootstrap repeats fitting and selection.
The all-homotopic-pairs calculation uses the same bootstrap command with
`--pairs homotopic --atlas-labels /absolute/path/to/ROI_MNI_V4.txt` and a
separate `--out` folder. It writes `bootstrap_family_summary.csv`.
The reported-pair output above supplies the bootstrap table.

## Frequency-band diagnostics

The following command reproduces the out-of-band variance fraction in the
main article and the comparisons between the lower, middle and upper thirds
of the pass band in Supplementary Appendix C.4:

~~~bash
python3 -m replication.adni.scripts.check_bands --aggregated "$ADNI_INPUT" --output "$ADNI_OUTPUT/band_diagnostics" --seed 20261101 --splits 100
~~~

It reads the same `combat/` tensors, `combat_pb102/` tensors and
`combat/covars_<group>.csv` files used above. The covariate rows must match the
tensor participant order. It removes an intercept and Age, Sex and APOE4 from
the projected coefficients before computing ROI correlations and partial
correlations. The variance fraction is computed before covariate adjustment.

The groups run in CN, MCI, Dementia order using one random-number stream.
Each random split partitions the 102 components into three sets of 34.
Keep the default group order, seed and 100 splits to reproduce the reported
random-split means.

The fresh output directory contains `summary.csv`, `bands.csv`,
`random_splits.csv` and `metadata.json`. The files contain group summaries
and comparison results, without participant records. Existing output
folders are rejected. The calculation uses inverse sample covariances and
does not refit graphical models or repeat the participant bootstrap.

## Figures and table

The connectome additionally requires external atlas and surface files:

| Asset | SHA-256 |
|---|---|
| ROI_MNI_V4.txt | fc3aa3001131895e74898b1034a915a93d5dbc644167f02110e895c913f3860e |
| ROI_MNI_V4.nii | 91e8ec62a293c2a28c8e7343f37abb978ab941c5e26fbc496945c66d073a2cd8 |
| BrainMesh_ICBM152.nv | fdbc3770838952193b0bcdabe38ce092afe718607b4edd9dfccab98a6647d4a9 |

Use the AAL V4 image and label order used in preprocessing.
The surface is `Data/SurfTemplate/BrainMesh_ICBM152.nv` from
[BrainNet Viewer commit 5163519](https://github.com/mingruixia/BrainNet-Viewer/tree/5163519457729c909029fe9df6c259f652dc7cdd).
Retain the upstream licence when obtaining the surface.
These assets are not redistributed.

~~~bash
python3 -m replication.adni.scripts.plot_connectome --fits "$ADNI_OUTPUT/selected_fits" --atlas-image /absolute/path/to/ROI_MNI_V4.nii --atlas-labels /absolute/path/to/ROI_MNI_V4.txt --surface /absolute/path/to/BrainMesh_ICBM152.nv --output "$ADNI_OUTPUT/figures/connectome"

python3 -m replication.adni.scripts.plot_criterion_curves --paths "$ADNI_OUTPUT/joint_profile/probe_path.csv" --verdicts "$ADNI_OUTPUT/joint_profile/probe_verdict.csv" --output "$ADNI_OUTPUT/figures/criteria"

python3 -m replication.adni.scripts.export_bootstrap_table --summary "$ADNI_OUTPUT/bootstrap/bootstrap_summary.csv" --metadata "$ADNI_OUTPUT/bootstrap/RUN_bootstrap.json" --output "$ADNI_OUTPUT/tables/TableC1.csv"
~~~

| Manuscript item | Output |
|---|---|
| Main Figure 4 | figures/connectome/Figure4.pdf |
| Supplement Figure C.1 | figures/criteria/FigureC1.pdf |
| Supplement Table C.1 | tables/TableC1.csv |

The curve figure plots the stored `ebic_plug` and
`ebic_joint_uncounted` columns for EBIC and Relaxation.

Connectome edges present in all three groups are dark grey.
Other edges are blue for CN, gold for MCI and red for Dementia.
The coloured edges are wider at the same absolute partial correlation.
BrainNet node and edge files are written alongside the connectome figure
under `brainnet/`.

## Tests

These commands use synthetic inputs and require no participant data:

~~~bash
Rscript replication/adni/scripts/run_tests.R
python3 -m pytest replication/adni/tests/python -q
~~~
