# R Code for "From Pest Traps to Management Maps..."
This repository contains the R scripts to reproduce the analyses in the manuscript:

**Bang et al. (In progress).** "From Pest Traps to Management Maps: A Predictive Framework for the Temporal Abundance Change of Japanese Pine Bast Scale to Guide National Forest Adaptation and Timely Control".

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
|   |-- 03_run_statistical_analysis.R
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
  "ggtext", "glue"
))
```

### 3. Run the Analysis Scripts
The scripts must be executed from the project root directory (pest-spatiotemporal-modeling-JPBS) in the following numerical order. The here() package ensures paths work correctly even though the scripts are in /codes/.
```r
source("codes/01_process_phenology_data.R")
source("codes/02_run_xgboost_pipeline.R")
source("codes/03_run_statistical_analysis.R")
```
After running all scripts, the output/ directory (created at the project root) will contain all reproducible artifacts, including model files (.rds), prediction rasters (.tif), and statistical tables (.csv).

## Citation
If you use this code or data in your research, please cite both the dataset and the associated manuscript.

### Dataset
Bang et al. (2025). Data for: "From Pest Traps to Management Maps: A Predictive Framework for the Temporal Abundance Change of Japanese Pine Bast Scale to Guide National Forest Adaptation and Timely Control" [Dataset]. Zenodo. https://doi.org/10.5281/zenodo.17403184
### Manuscript
Bang et al. (In progress). "From Pest Traps to Management Maps: A Predictive Framework for the Temporal Abundance Change of Japanese Pine Bast Scale to Guide National Forest Adaptation and Timely Control".

## Acknowledgments
Google Gemini and OpenAI ChatGPT were used to assist with code refinement and debugging.

## License
This project is licensed under the MIT License. See the LICENSE file for details.