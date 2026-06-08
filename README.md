# R Code for "From Pest Traps to Management Maps..."
This repository contains the R scripts to reproduce the analyses in the manuscript:

**Bang et al. (In progress).** "From Pest Traps to Management Maps: Predicting the Abundance and Phenology of Japanese Pine Bast Scale to Guide National Forest Adaptation and Timely Control".

The code was written by **Seunguk Kim** (adrenaline@snu.ac.kr).

--- 
## How to Reproduce the Analysis
Follow these three steps to set up the environment and reproduce the results.

### 1. Download Data
All required input data (occurrence .csv files and predictor .tif rasters) are permanently archived on Zenodo. Download the complete dataset from the link below and unzip it.
Zenodo Archive: https://doi.org/10.5281/zenodo.17403184

### 2. Set Up the R Environment
**2.1.** Clone this GitHub repository to your local machine.
**2.2.** Place the downloaded data from Zenodo into the data directory. Your project folder should now look like this:
```
/pest-spatiotemporal-modeling-JPBS
|-- codes/
|   |-- 01_process_phenology_data.R
|   |-- 02_run_xgboost_pipeline.R
|   |-- 02b_tidymodels_reproduction.R   (optional; reproduces script 02's modeling on a current xgboost)
|   |-- 03_run_statistical_analysis.R
|   |-- 04_evaluate_spatial_temporal_performance.R
|-- data/
|   |-- occurrence/
|   |-- predictors/1km/
|-- ... etc.
```
**2.3.** Install the required R packages. Run the following command in your R console:
```r
install.packages(c(
  "here", "tidyverse", "SuppDists", "sf", "terra", "xgboost", "caret",	
  "doParallel", "doRNG", "openxlsx", "ggbeeswarm", "patchwork",
  "RColorBrewer", "matrixStats", "ggcorrplot", "viridis", "scales",
  "ggtext", "glue",
  "parsnip", "tune", "rsample", "yardstick", "dials", "workflows"
))
```
**2.4.** *(Recommended.)* The repository includes an `renv.lock` file that pins the package versions of the environment used to run this repository. To restore it, run:
```r
install.packages("renv")
renv::restore()
```
This environment uses a current version of `xgboost`. Under it, scripts 01, 04, and 02b run directly, while scripts 02 and 03 require an older `xgboost` (< 2.0); see the note at the end of Section 3. Script 02b reproduces script 02's model building, and the reported performance, under the current (locked) environment using `tidymodels`.

### 3. Run the Analysis Scripts
Run the scripts from the project root directory (pest-spatiotemporal-modeling-JPBS). The here() package ensures paths work correctly even though the scripts are in /codes/. Which scripts you run depends on your `xgboost` version (see the note at the end of this section).

**With an older `xgboost` (< 2.0)** — the environment in which the published results were produced — run the full pipeline in order:
```r
source("codes/01_process_phenology_data.R")
source("codes/02_run_xgboost_pipeline.R")
source("codes/03_run_statistical_analysis.R")
source("codes/04_evaluate_spatial_temporal_performance.R")
```

**With a current `xgboost` (>= 2.0)** — for example, the pinned `renv` environment — scripts 02 and 03 are not compatible (see the note below). Run instead:
```r
source("codes/01_process_phenology_data.R")
source("codes/02b_tidymodels_reproduction.R")   # model building, in place of script 02
source("codes/04_evaluate_spatial_temporal_performance.R")
```
The small trained artifacts that script 04 needs (`master_data.rds`, the tuned hyperparameters, and the feature sets) are included in `output/`, so script 04 runs without re-running script 02.
After running all scripts, the output/ directory (created at the project root) will contain all reproducible artifacts, including model files (.rds), prediction rasters (.tif), and statistical tables (.csv). Script 04 writes the supplementary evaluation outputs to output/evaluation/, including the spatial-versus-temporal cross-validation table (Table S6), the observed-versus-predicted figure (Figure S4), and the supporting diagnostics reported in Text S4 (year-as-predictor, paired-site interannual shift, out-of-year rank correlation, and the extreme-value/outbreak metrics).

### A note on `xgboost` versions (and optional script 02b)
The published models were built with `xgboost` < 2.0. Two scripts depend on that environment:
- **Script 02** tunes the models with `caret::train(method = "xgbTree")`, which is not compatible with `xgboost` >= 2.0 (the tuning step errors). Script 02 still builds and saves `master_data.rds` before that step.
- **Script 03** loads the fitted `xgboost` model objects that script 02 saved (`definitive_model_*.rds`); model objects serialized under `xgboost` < 2.0 cannot be reloaded under `xgboost` >= 2.0.

Scripts 01, 04, and the optional script 02b run on a current `xgboost`. Script 02b reproduces script 02's model building (hyperparameter tuning, the 500-repeat 85/15 evaluation, and the final model fit) using the actively maintained `tidymodels`/`parsnip` framework. It reproduces the reported predictive performance to within about 0.01 R². Because gradient boosting has a flat optimum, the exact selected hyperparameters and metrics vary slightly with the framework, the cross-validation randomization, and the `xgboost` version; this is expected and does not change the conclusions.
```r
source("codes/02b_tidymodels_reproduction.R")
```

## Citation
If you use this code or data in your research, please cite both the dataset and the associated manuscript.

### Dataset
Bang et al. (2025). Data for: "From Pest Traps to Management Maps: Predicting the Abundance and Phenology of Japanese Pine Bast Scale to Guide National Forest Adaptation and Timely Control" [Dataset]. Zenodo. https://doi.org/10.5281/zenodo.17403184
### Manuscript
Bang et al. (In progress). "From Pest Traps to Management Maps: Predicting the Abundance and Phenology of Japanese Pine Bast Scale to Guide National Forest Adaptation and Timely Control".

## Acknowledgments
Google Gemini, OpenAI ChatGPT, and Anthropic Claude Code were used to assist with code refinement and debugging.

## License
This project is licensed under the MIT License. See the LICENSE file for details.