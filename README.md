# R Code for "From Pest Traps to Management Maps..."
## This repository contains the R scripts to reproduce the analyses in the manuscript:
## Bang, S., Kim, S., et al. "From Pest Traps to Management Maps: A Predictive Framework for the Temporal Abundance Change of Japanese Pine Bast Scale to Guide National Forest Adaptation and Timely Control"
## The code was written by Seunguk Kim (adrenaline@snu.ac.kr).

How to Reproduce the Analysis
Follow these three steps to set up the environment and reproduce the results.

1. Download Data
All required input data (occurrence .csv files and predictor .tif rasters) are permanently archived on Zenodo. Download the complete dataset from the link below and unzip it.
Zenodo Archive: https://doi.org/10.5281/zenodo.17403184

2. Set Up the R Environment
Clone this GitHub repository to your local machine.
Place the downloaded data from Zenodo into the data directory. Your project folder should now look like this:
/yourprojectfolder
|-- data/
|   |-- occurrence/
|   |-- predictors/
|-- ... etc.
Install the required R packages. Run the following command in your R console:

3. Run the Analysis Scripts
The scripts must be executed in the following numerical order:
01_process_phenology_data.R
02_run_xgboost_pipeline.R
03_run_statistical_analysis.R
After running all scripts, the output/ directory will contain all reproducible artifacts, including model files (.rds), prediction rasters (.tif), and statistical tables (.csv).

Citation
If you use this code or data in your research, please cite both the dataset and the associated manuscript.

Dataset:
Seunguk Kim. (2025). Data for: "From Pest Traps to Management Maps: A Predictive Framework for the Temporal Abundance Change of Japanese Pine Bast Scale" [Data set]. Zenodo. https://doi.org/10.5281/zenodo.17403184
Manuscript:
Bang, S., Kim, S., et al. (Year). "From Pest Traps to Management Maps...". Earth's Future. [Add final paper citation/DOI when available]

License
This project is licensed under the MIT License. See the LICENSE file for details.