# -------------------------------------------------------------------------- #
# SCRIPT 2: XGBoost Modeling and Prediction Pipeline
# -------------------------------------------------------------------------- #
#
# Author: Seunguk Kim
# Contact: adrenaline@snu.ac.kr
# Last Modified: 2025-10-21
#
# Description:
# This script reproduces the core analysis for the manuscript, including
# model tuning, robust evaluation, final model training, spatial
# prediction (.tif raster generation), and SHAP analysis (Figure 6).
#
# Associated Manuscript:
# Bang, S., Kim, S., et al. "From Pest Traps to Management Maps..."
#
# Usage Notes:
# This is the MAIN script. It requires two sets of inputs:
# 1. The output from '01_process_phenology_data.R' (phenology_metrics_clean.csv).
# 2. The provided predictor rasters (Host maps, Bioclim, LST, etc.)
#    placed in the 'data/predictors/1km/' directory.
#
# -------------------------------------------------------------------------- #

# RStudio Document Outline (####) Pipeline Flow:
#
# 0. Setup: Load libraries, define global variables, and guide package management.
# 1. Data Preparation: Load and merge occurrence/phenology data with predictors.
#
# --- Part 1: Model Training & Evaluation (Executes if Excel report is missing) ---
# 2. Tuning & Robust Evaluation: Hyperparameter tuning and 500-repeat validation.
# 3. Final Model Training: Train models on all data and save them.
# 4. Aggregate & Save Results: Collate all performance metrics into an Excel file.
#
# --- Part 2: Spatial Prediction & Visualization ---
# 5. Prediction & Uncertainty Assessment: Generate prediction/uncertainty rasters (.tif)
#    using a 50-bootstrap ensemble.
#

#### 0. Setup ####
print("--- 0. Setup ---")
print(paste("Pipeline started at:", Sys.time()))

#---- 0.1. Load Libraries ----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(terra)
  library(xgboost)
  library(caret)
  library(doParallel)
  library(doRNG)
  library(openxlsx)
  library(ggbeeswarm)
  library(patchwork)
  library(sf)
  library(RColorBrewer)
  library(matrixStats)
  library(ggcorrplot)
  library(viridis)
  library(scales)
  library(ggtext)
  library(glue)
})

#---- 0.2. Global Settings ----
target_resolution <- "1km"
response_vars <- c("Total_Abundance", "Onset_Week", "Peak_Week")
N_REPEATS <- 500      # Number of repeats for robust performance evaluation
N_BOOTSTRAPS <- 50    # Number of bootstrap models for uncertainty assessment

n_cores <- parallel::detectCores() - 1
if (n_cores < 1)
  n_cores <- 1

final_results_list <- list()

#---- 0.3. Package Version Management (renv) ----
# This project uses 'renv' to ensure reproducibility.
# If you are running this for the first time, run:
# > renv::restore()
# This will install all required packages at their correct versions
# from the 'renv.lock' file.
#
# If you modify the code (e.g., add new packages):
# 1. install.packages("new_package_name")
# 2. renv::snapshot()
# 3. Commit the updated 'renv.lock' file to your repository.


#### 1. Data Preparation ####
print("--- 1. Data Preparation ---")

#---- 1.1. Paths and Directories ----
# Base output folder
dir_output_base <- here("output")
dir_data_prep <- file.path(dir_output_base, "prepared_data")
dir_models <- file.path(dir_output_base, "models")
dir_evaluation <- file.path(dir_output_base, "evaluation")
dir_predictions <- file.path(dir_output_base, "predictions")

# Create directories if they don't exist
dir.create(dir_data_prep, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_models, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_evaluation,
           showWarnings = FALSE,
           recursive = TRUE)
dir.create(dir_predictions,
           showWarnings = FALSE,
           recursive = TRUE)

# Input data path
dir_raw_predictors <- here("data/predictors", target_resolution)

#---- 1.2. Variable Name Map ----
# Map for renaming Bioclim variables for plots/analysis
bio_name_map <- c(
  "bio1" = "Annual Mean Temp",
  "bio2" = "Mean Diurnal Range",
  "bio3" = "Isothermality",
  "bio4" = "Temp Seasonality",
  "bio5" = "Max Temp of Warmest Month",
  "bio6" = "Min Temp of Coldest Month",
  "bio7" = "Temp Annual Range",
  "bio8" = "Mean Temp of Wettest Quarter",
  "bio9" = "Mean Temp of Driest Quarter",
  "bio10" = "Mean Temp of Warmest Quarter",
  "bio11" = "Mean Temp of Coldest Quarter",
  "bio12" = "Annual Precipitation",
  "bio13" = "Precip of Wettest Month",
  "bio14" = "Precip of Driest Month",
  "bio15" = "Precip Seasonality",
  "bio16" = "Precip of Wettest Quarter",
  "bio17" = "Precip of Driest Quarter",
  "bio18" = "Precip of Warmest Quarter",
  "bio19" = "Precip of Coldest Quarter"
)
host_abundance_map <- c("Ab_AG" = "Host Abundance (Aggregate)",
                        "Ab_PD" = "Host Abundance (Pinus densiflora)",
                        "Ab_PT" = "Host Abundance (Pinus thunbergii)")
full_name_map <- c(bio_name_map, host_abundance_map)

#---- 1.3. Load or Create Master Data Frame ----
master_data_path <- file.path(dir_data_prep, "master_data.rds")

if (file.exists(master_data_path)) {
  print("Master data file (master_data.rds) already exists. Loading from disk.")
  master_data <- readRDS(master_data_path)
} else {
  print("Creating master data file from scratch.")
  
  # NOTE: Ensure all input CSV files (occurrence*.csv, phenology_metrics_clean.csv)
  # are saved with UTF-8 encoding to avoid read errors.
  # The 'fileEncoding = "EUC-KR"' argument has been removed.
  
  dir_raw_occurrence <- here("data/occurrence")
  occ2022_path <-
    file.path(dir_raw_occurrence, "occurrence2022.csv")
  occ2023_path <-
    file.path(dir_raw_occurrence, "occurrence2023.csv")
  
  abundance_data_full <- bind_rows(
    read.csv(occ2022_path) %>% dplyr::select(Area, Region, trap_no, lon, lat, Total_Abundance = total) %>% mutate(Year = 2022),
    read.csv(occ2023_path) %>% dplyr::select(Area, Region, trap_no, lon, lat, Total_Abundance = total) %>% mutate(Year = 2023)
  ) %>%
    mutate(point_id = paste(Year, Area, Region, trap_no, lon, lat, sep = "_"))
  
  phenology_metrics_path <-
    file.path(dir_raw_occurrence, "phenology_metrics_clean.csv")
  phenology_data_full <-
    read.csv(phenology_metrics_path) %>%
    mutate(point_id = paste(Year, Area, Region, trap_no, lon, lat, sep = "_"))
  
  #--- Preliminary check for missing data ---
  predictor_files_pre <- list.files(dir_raw_predictors,
                                    pattern = "\\.tif$",
                                    full.names = TRUE)
  if (length(predictor_files_pre) == 0)
    stop("Predictor .tif files not found for ", target_resolution)
  
  predictor_stack_pre <- terra::rast(predictor_files_pre)
  
  points_sf_ab <- sf::st_as_sf(
    abundance_data_full,
    coords = c("lon", "lat"),
    crs = 4326,
    remove = FALSE
  )
  extracted_vals_ab <-
    terra::extract(predictor_stack_pre, terra::vect(points_sf_ab), ID = FALSE)
  valid_points_ab <-
    abundance_data_full %>% bind_cols(extracted_vals_ab) %>% drop_na()
  
  points_sf_ph <- sf::st_as_sf(
    phenology_data_full,
    coords = c("lon", "lat"),
    crs = 4326,
    remove = FALSE
  )
  extracted_vals_ph <-
    terra::extract(predictor_stack_pre, terra::vect(points_sf_ph), ID = FALSE)
  valid_points_ph <-
    phenology_data_full %>% bind_cols(extracted_vals_ph) %>% drop_na()
  
  # Filter data to common valid points
  common_ids_abundance <- valid_points_ab$point_id
  common_ids_phenology <- valid_points_ph$point_id
  
  abundance_data_common <-
    abundance_data_full %>% filter(point_id %in% common_ids_abundance)
  phenology_data_common <-
    phenology_data_full %>% filter(point_id %in% common_ids_phenology)
  
  combined_data_common <- abundance_data_common %>%
    left_join(
      phenology_data_common %>% dplyr::select(point_id, Onset_Week = First_Occurrence_5pct, Peak_Week = Peak_Occurrence_Mode),
      by = "point_id"
    )
  
  #--- Load and stack predictors by year ---
  predictor_files <- list.files(dir_raw_predictors,
                                pattern = "\\.tif$",
                                full.names = TRUE)
  static_files <-
    str_subset(predictor_files, "_2022|_2023", negate = TRUE)
  dynamic_files_2022 <- str_subset(predictor_files, "_2022")
  dynamic_files_2023 <- str_subset(predictor_files, "_2023")
  
  static_stack <- terra::rast(static_files)
  dynamic_stack_2022 <- terra::rast(dynamic_files_2022)
  dynamic_stack_2023 <- terra::rast(dynamic_files_2023)
  
  names(static_stack) <-
    tools::file_path_sans_ext(basename(static_files))
  names(dynamic_stack_2022) <-
    str_remove(tools::file_path_sans_ext(basename(dynamic_files_2022)), "_2022")
  names(dynamic_stack_2023) <-
    str_remove(tools::file_path_sans_ext(basename(dynamic_files_2023)), "_2023")
  
  predictor_stack_2022 <- c(static_stack, dynamic_stack_2022)
  predictor_stack_2023 <- c(static_stack, dynamic_stack_2023)
  
  #--- Extract year-specific values ---
  points_sf_2022 <-
    combined_data_common %>% filter(Year == 2022) %>% sf::st_as_sf(
      coords = c("lon", "lat"),
      crs = 4326,
      remove = FALSE
    )
  points_sf_2023 <-
    combined_data_common %>% filter(Year == 2023) %>% sf::st_as_sf(
      coords = c("lon", "lat"),
      crs = 4326,
      remove = FALSE
    )
  
  extracted_2022 <-
    terra::extract(predictor_stack_2022, terra::vect(points_sf_2022), bind = TRUE) %>% as.data.frame()
  extracted_2023 <-
    terra::extract(predictor_stack_2023, terra::vect(points_sf_2023), bind = TRUE) %>% as.data.frame()
  
  master_data_raw <- bind_rows(extracted_2022, extracted_2023)
  
  #--- Final data processing ---
  master_data <- master_data_raw %>%
    dplyr::rename_with( ~ ifelse(. %in% names(full_name_map), full_name_map[.], .)) %>%
    dplyr::mutate(Total_Abundance = base::log1p(Total_Abundance)) %>% # Log-transform abundance
    dplyr::select(-point_id)
  
  saveRDS(master_data, file = master_data_path)
  master_csv_path <- file.path(dir_data_prep, "master_data.csv")
  write.csv(master_data, file = master_csv_path, row.names = FALSE)
  print("Master data file created and saved.")
}


#### 2. Tuning and Robust Evaluation ####
print(paste("\n--- PART 1: Model Tuning & Robust Evaluation ---"))
print(paste("Part 1 started at:", Sys.time()))

output_excel_path <- file.path(dir_output_base,
                               paste0("XGBoost_Analysis_Robust_", target_resolution, ".xlsx"))
if (file.exists(output_excel_path)) {
  message("Final Excel report already exists. Skipping Part 1 (Tuning, Evaluation, Final Training).")
} else {
  cl <- parallel::makeCluster(n_cores)
  registerDoParallel(cl)
  registerDoRNG(123) # for reproducible parallel loops
  print(paste(n_cores, "cores will be used for parallel processing."))
  
  for (response_var in response_vars) {
    print(paste("\n--- Processing for:", response_var, "---"))
    
    # Use all data for Abundance, filter NAs for phenology
    current_data <- if (response_var == "Total_Abundance")
      master_data
    else
      master_data %>% filter(!is.na(.data[[response_var]]))
    
    if (nrow(current_data) < 20) {
      print("Not enough data. Skipping.")
      next
    }
    
    predictors_df <- current_data %>% dplyr::select(-dplyr::any_of(
      c(
        "Year",
        "Area",
        "Region",
        "trap_no",
        "lon",
        "lat",
        response_vars
      )
    ))
    
    #---- 2.1. Hyperparameter Tuning ----
    print("Hyperparameter Tuning with all features...")
    tune_grid <- expand.grid(
      nrounds = c(100, 200),
      max_depth = c(2, 3, 5),
      eta = c(0.01, 0.05, 0.1),
      gamma = c(0.1, 1),
      colsample_bytree = 0.8,
      subsample = 0.8,
      min_child_weight = c(1, 5)
    )
    
    train_control <- caret::trainControl(
      method = "repeatedcv",
      number = 5,
      repeats = 2,
      allowParallel = FALSE # Parallelism is handled by foreach
    )
    
    tuned_model <- caret::train(
      x = predictors_df,
      y = current_data[[response_var]],
      method = "xgbTree",
      trControl = train_control,
      tuneGrid = tune_grid,
      verbose = FALSE,
      verbosity = 0,
      nthread = n_cores
    )
    
    definitive_best_params <- tuned_model$bestTune
    final_results_list[[response_var]]$best_params <-
      definitive_best_params
    saveRDS(
      definitive_best_params,
      file.path(
        dir_models,
        paste0("definitive_best_params_", response_var, ".rds")
      )
    )
    cat("Optimal hyperparameters determined and saved.\n")
    
    #---- 2.2. Robust Performance Evaluation ----
    print(paste(
      "Starting robust performance evaluation with",
      N_REPEATS,
      "repeats..."
    ))
    
    results_chunk <- foreach(
      i = 1:N_REPEATS,
      .packages = c("caret", "xgboost", "dplyr"),
      .combine = 'c'
    ) %dopar% {
      set.seed(i)
      train_val_index <-
        createDataPartition(current_data[[response_var]], p = 0.85, list = FALSE)
      train_data <- current_data[train_val_index, ]
      test_data <- current_data[-train_val_index, ]
      
      train_predictors <-
        train_data %>% dplyr::select(names(predictors_df))
      dtrain <-
        xgb.DMatrix(data = data.matrix(train_predictors), label = train_data[[response_var]])
      
      params <- as.list(definitive_best_params)
      params$objective <- "reg:squarederror"
      params$nthread <- 1 # Each parallel worker uses 1 thread
      
      model_iter <- xgb.train(
        params = params,
        data = dtrain,
        nrounds = definitive_best_params$nrounds,
        verbose = 0
      )
      
      predictions <- predict(model_iter, data.matrix(test_data %>% dplyr::select(names(
        predictors_df
      ))))
      
      list(
        performance = postResample(pred = predictions, obs = test_data[[response_var]]),
        importance = xgb.importance(model = model_iter)
      )
    }
    
    performance_results <-
      results_chunk[names(results_chunk) == "performance"]
    importance_results <-
      results_chunk[names(results_chunk) == "importance"]
    
    perf_df <- do.call(rbind, performance_results) %>% as.data.frame()
    final_results_list[[response_var]]$robust_performance <- perf_df
    
    imp_df <- bind_rows(importance_results, .id = "repetition")
    agg_importance <- imp_df %>%
      group_by(Feature) %>%
      summarise(
        Mean_Gain = mean(Gain),
        SD_Gain = sd(Gain),
        .groups = 'drop'
      ) %>%
      arrange(desc(Mean_Gain))
    final_results_list[[response_var]]$aggregated_importance <-
      agg_importance
    
    cat("\nRobust performance evaluation finished.\n")
  }
  stopCluster(cl)
  print("Parallel processing cluster stopped for Tuning and Evaluation.")
  
  
  #### 3. Final Model Training and Interpretation ####
  print("\n--- 3. Final Model Training & Interpretation ---")
  for (response_var in response_vars) {
    cat(paste("\n--- Finalizing model for:", response_var, "---\n"))
    if (is.null(final_results_list[[response_var]])) {
      next
    }
    
    current_data <- if (response_var == "Total_Abundance")
      master_data
    else
      master_data %>% filter(!is.na(.data[[response_var]]))
    
    final_predictors_df <-
      current_data %>% dplyr::select(-dplyr::any_of(
        c(
          "Year",
          "Area",
          "Region",
          "trap_no",
          "lon",
          "lat",
          response_vars
        )
      ))
    final_target <- current_data[[response_var]]
    
    #---- 3.1. Train Definitive Model ----
    print("Training the definitive model on the full dataset...")
    dtrain_definitive <-
      xgb.DMatrix(data = data.matrix(final_predictors_df), label = final_target)
    definitive_best_params <-
      final_results_list[[response_var]]$best_params
    
    final_params <- as.list(definitive_best_params)
    final_params$objective <- "reg:squarederror"
    final_params$nthread <- detectCores()
    
    definitive_model <- xgb.train(
      params = final_params,
      data = dtrain_definitive,
      nrounds = definitive_best_params$nrounds,
      verbose = 0
    )
    
    all_features <- names(final_predictors_df)
    saveRDS(
      definitive_model,
      file.path(dir_models, paste0(
        "definitive_model_", response_var, ".rds"
      ))
    )
    saveRDS(
      all_features,
      file.path(dir_models, paste0(
        "definitive_features_", response_var, ".rds"
      ))
    )
    final_results_list[[response_var]]$definitive_model <-
      definitive_model
    cat("Final model and feature list saved.\n")
    
    #---- 3.2. SHAP Analysis (for saving) ----
    cat("Calculating and saving SHAP values...\n")
    shap_values <- predict(definitive_model,
                           data.matrix(final_predictors_df),
                           predcontrib = TRUE)
    
    mean_abs_shap <- as.data.frame(shap_values) %>%
      select(-BIAS) %>% abs() %>% colMeans() %>% sort(decreasing = TRUE) %>%
      tibble::enframe(name = "Feature", value = "Mean_Abs_SHAP")
    
    final_results_list[[response_var]]$shap_importance <-
      mean_abs_shap
    
    cat("Generating individual SHAP plots (for diagnostics)...\n")
    shap_long <-
      as.data.frame(shap_values) %>% select(-BIAS) %>% mutate(row_id = row_number()) %>%
      pivot_longer(-row_id, names_to = "feature", values_to = "shap_value")
    
    feature_values_long <-
      as.data.frame(final_predictors_df) %>% mutate(row_id = row_number()) %>%
      pivot_longer(-row_id, names_to = "feature", values_to = "feature_value")
    
    shap_plot_data <-
      left_join(shap_long, feature_values_long, by = c("row_id", "feature"))
    
    shap_plot_data_scaled <- shap_plot_data %>% group_by(feature) %>%
      mutate(feature_value_scaled = as.numeric(scale(feature_value))) %>% ungroup() %>%
      mutate(feature_value_scaled = if_else(is.nan(feature_value_scaled), 0, feature_value_scaled))
    
    top_20_features <- mean_abs_shap %>% head(20) %>% pull(Feature)
    
    shap_plot_data_scaled_top20 <- shap_plot_data_scaled %>%
      filter(feature %in% top_20_features) %>%
      mutate(feature = factor(feature, levels = rev(top_20_features)))
    
    shap_plot <- ggplot(shap_plot_data_scaled_top20,
                        aes(x = shap_value, y = feature)) +
      geom_vline(
        xintercept = 0,
        color = "gray",
        linetype = "dashed"
      ) +
      ggbeeswarm::geom_quasirandom(
        aes(color = feature_value_scaled),
        groupOnX = FALSE,
        alpha = 0.7,
        orientation = "y"
      ) +
      scale_color_gradient2(
        low = "blue",
        mid = "gray",
        high = "red",
        midpoint = 0,
        name = "Feature\nValue\n(Scaled)"
      ) +
      labs(
        title = paste("SHAP Summary Plot for", response_var, "(Top 20 Features)"),
        subtitle = "Definitive Model (trained on all data)",
        x = "SHAP Value (Impact on model output)",
        y = "Feature"
      ) +
      theme_bw(base_size = 12)
    
    ggsave(
      file.path(
        dir_evaluation,
        paste0("definitive_shap_plot_", response_var, ".png")
      ),
      plot = shap_plot,
      width = 10,
      height = 8
    )
    cat("Diagnostic SHAP plots saved.\n")
  }
  
  
  #### 4. Aggregate and Save Results ####
  print("\n--- 4. Aggregate and Save Results ---")
  tryCatch({
    #---- 4.1. Summarize Results ----
    robust_perf_summary <-
      map_dfr(final_results_list, "robust_performance", .id = "Response_Variable") %>%
      group_by(Response_Variable) %>%
      summarise(across(everything(), list(mean = mean, sd = sd)), .groups = 'drop')
    agg_importance_df <-
      map_dfr(final_results_list, "aggregated_importance", .id = "Response_Variable")
    shap_importance_df <-
      map_dfr(final_results_list, "shap_importance", .id = "Response_Variable")
    
    #---- 4.2. Descriptive Statistics ----
    renamed_vars <- unname(full_name_map)
    data_summary <- master_data %>%
      select(any_of(
        c(
          response_vars,
          renamed_vars,
          "Elevation",
          "EVI",
          "LST_Day",
          "LST_Night",
          "NDVI",
          "SWIR"
        )
      )) %>%
      pivot_longer(everything(), names_to = "Variable", values_to = "Value") %>%
      group_by(Variable) %>%
      summarise(
        N = sum(!is.na(Value)),
        Mean = mean(Value, na.rm = TRUE),
        SD = sd(Value, na.rm = TRUE),
        Min = min(Value, na.rm = TRUE),
        Median = median(Value, na.rm = TRUE),
        Max = max(Value, na.rm = TRUE),
        .groups = 'drop'
      )
    
    #---- 4.3. Save to Excel ----
    wb <- createWorkbook()
    addWorksheet(wb, "Robust_Performance")
    writeData(wb, "Robust_Performance", robust_perf_summary, borders = "all")
    addWorksheet(wb, "SHAP_Importance")
    writeData(wb, "SHAP_Importance", shap_importance_df, borders = "all")
    addWorksheet(wb, "Gain_Importance")
    writeData(wb, "Gain_Importance", agg_importance_df, borders = "all")
    addWorksheet(wb, "Data_Summary")
    writeData(wb, "Data_Summary", data_summary, borders = "all")
    for (sheet_name in names(wb)) {
      setColWidths(wb,
                   sheet = sheet_name,
                   cols = 1:20,
                   widths = "auto")
    }
    saveWorkbook(wb, output_excel_path, overwrite = TRUE)
    print(paste("All results aggregated and saved to:", output_excel_path))
  }, error = function(e) {
    print(paste(
      "An error occurred during final result aggregation:",
      e$message
    ))
  })
}
print(paste("Part 1 finished at:", Sys.time()))


#### 5. Prediction and Uncertainty Assessment ####
print(paste(
  "\n--- PART 2: Spatial Prediction & Uncertainty Assessment (",
  N_BOOTSTRAPS,
  " repeats) ---"
))
print(paste("Part 2 started at:", Sys.time()))

cl <- parallel::makeCluster(n_cores)
registerDoParallel(cl)
registerDoRNG(456)
message(paste(n_cores, "cores will be used for parallel processing."))

#---- 5.1. Prediction Wrapper Function ----
# Wrapper for terra::predict to work with xgb.DMatrix
predict_wrapper <- function(model, ...) {
  newdata <- list(...)[[1]]
  predict(model, data.matrix(newdata))
}

tryCatch({
  for (response_var in response_vars) {
    cat(paste("\n--- Processing for:", response_var, "---\n"))
    
    # Check if final output files already exist
    expected_files <- c(
      file.path(
        dir_predictions,
        paste0("prediction_mean_", response_var, "_2022.tif")
      ),
      file.path(
        dir_predictions,
        paste0("prediction_mean_", response_var, "_2023.tif")
      ),
      file.path(
        dir_predictions,
        paste0("uncertainty_width_", response_var, "_2022.tif")
      ),
      file.path(
        dir_predictions,
        paste0("uncertainty_width_", response_var, "_2023.tif")
      )
    )
    if (all(sapply(expected_files, file.exists))) {
      message(
        paste(
          "  - All prediction and uncertainty rasters for",
          response_var,
          "already exist. Skipping."
        )
      )
      next
    }
    
    features_path <- file.path(dir_models,
                               paste0("definitive_features_", response_var, ".rds"))
    params_path <- file.path(dir_models,
                             paste0("definitive_best_params_", response_var, ".rds"))
    model_path <- file.path(dir_models,
                            paste0("definitive_model_", response_var, ".rds"))
    if (!file.exists(features_path) ||
        !file.exists(params_path) || !file.exists(model_path)) {
      message("  - Final model files not found. Skipping prediction.")
      next
    }
    
    #---- 5.2. Load Model and Train Bootstrap Ensemble ----
    final_features <- readRDS(features_path)
    best_params <- readRDS(params_path)
    definitive_model <- readRDS(model_path)
    
    current_data <- if (response_var == "Total_Abundance")
      master_data
    else
      master_data %>% filter(!is.na(.data[[response_var]]))
    training_df <-
      current_data %>% select(all_of(c(response_var, final_features)))
    
    message(paste(
      "  - Training bootstrap ensemble (",
      N_BOOTSTRAPS,
      "models)..."
    ))
    tuned_params <- as.list(best_params)
    tuned_params$nrounds <- NULL # nrounds is specified in xgb.train
    tuned_params$objective <- "reg:squarederror"
    tuned_params$nthread <- 1 # Each parallel worker uses 1 thread
    
    bootstrap_models <- foreach(b = 1:N_BOOTSTRAPS,
                                .packages = c('xgboost', 'dplyr')) %dopar% {
                                  bootstrap_sample_df <-
                                    training_df %>% slice_sample(n = nrow(training_df), replace = TRUE)
                                  dtrain_boot <-
                                    xgb.DMatrix(data = data.matrix(bootstrap_sample_df[, final_features]),
                                                label = bootstrap_sample_df[[response_var]])
                                  xgb.train(
                                    params = tuned_params,
                                    data = dtrain_boot,
                                    nrounds = best_params$nrounds,
                                    verbose = 0
                                  )
                                }
    
    #---- 5.3. Generate Prediction and Uncertainty Maps by Year ----
    predictor_files <- list.files(dir_raw_predictors,
                                  pattern = "\\.tif$",
                                  full.names = TRUE)
    static_files <-
      str_subset(predictor_files, "_2022|_2023", negate = TRUE)
    static_stack <- terra::rast(static_files)
    names(static_stack) <-
      dplyr::recode(tools::file_path_sans_ext(basename(static_files)),!!!full_name_map)
    
    for (year in c(2022, 2023)) {
      message(paste("  - Processing maps for year", year, "..."))
      
      dynamic_files <- str_subset(predictor_files, as.character(year))
      dynamic_stack <- terra::rast(dynamic_files)
      names(dynamic_stack) <-
        str_remove(tools::file_path_sans_ext(basename(dynamic_files)),
                   paste0("_", year))
      names(dynamic_stack) <-
        dplyr::recode(names(dynamic_stack),!!!full_name_map)
      
      predictor_stack_year <- c(static_stack, dynamic_stack)
      # Ensure order matches training
      predictor_stack_final <-
        predictor_stack_year[[final_features]]
      
      message("    - Generating final prediction map (mean)...")
      pred_map <- terra::predict(
        predictor_stack_final,
        definitive_model,
        fun = predict_wrapper,
        na.rm = TRUE
      )
      if (response_var == "Total_Abundance") {
        pred_map <- expm1(pred_map) # Back-transform from log1p
      }
      pred_fname <- file.path(
        dir_predictions,
        paste0("prediction_mean_", response_var, "_", year, ".tif")
      )
      terra::writeRaster(
        pred_map,
        pred_fname,
        overwrite = TRUE,
        gdal = c("COMPRESS=LZW")
      )
      
      message("    - Generating uncertainty maps (processing by block)...")
      terra::readStart(predictor_stack_final)
      tryCatch({
        out_lower <- terra::rast(predictor_stack_final, nlyrs = 1)
        out_upper <- terra::rast(predictor_stack_final, nlyrs = 1)
        out_width <- terra::rast(predictor_stack_final, nlyrs = 1)
        
        b <- terra::blocks(predictor_stack_final)
        
        lower_fname <- file.path(
          dir_predictions,
          paste0("uncertainty_lower_", response_var, "_", year, ".tif")
        )
        upper_fname <- file.path(
          dir_predictions,
          paste0("uncertainty_upper_", response_var, "_", year, ".tif")
        )
        width_fname <- file.path(
          dir_predictions,
          paste0("uncertainty_width_", response_var, "_", year, ".tif")
        )
        
        # Clean up old files if they exist
        if (file.exists(lower_fname))
          unlink(lower_fname, force = TRUE)
        if (file.exists(upper_fname))
          unlink(upper_fname, force = TRUE)
        if (file.exists(width_fname))
          unlink(width_fname, force = TRUE)
        
        terra::writeStart(
          out_lower,
          filename = lower_fname,
          overwrite = TRUE,
          gdal = c("COMPRESS=LZW")
        )
        terra::writeStart(
          out_upper,
          filename = upper_fname,
          overwrite = TRUE,
          gdal = c("COMPRESS=LZW")
        )
        terra::writeStart(
          out_width,
          filename = width_fname,
          overwrite = TRUE,
          gdal = c("COMPRESS=LZW")
        )
        
        for (i in 1:b$n) {
          cat(paste("\r      - Processing block", i, "of", b$n, "..."))
          block_vals <- terra::readValues(
            predictor_stack_final,
            row = b$row[i],
            nrows = b$nrows[i],
            col = 1,
            ncols = ncol(predictor_stack_final),
            mat = TRUE
          )
          
          valid_rows <- complete.cases(block_vals)
          if (sum(valid_rows) > 0) {
            preds_matrix <- sapply(bootstrap_models, function(m)
              predict(m, block_vals[valid_rows, , drop = FALSE]))
            
            # Calculate 5th and 95th percentiles
            quantiles <- matrixStats::rowQuantiles(
              preds_matrix,
              probs = c(0.05, 0.95),
              na.rm = TRUE,
              useNames = FALSE
            )
            
            if (response_var == "Total_Abundance") {
              lower_vals_orig <- expm1(quantiles[, 1])
              upper_vals_orig <- expm1(quantiles[, 2])
            } else {
              lower_vals_orig <- quantiles[, 1]
              upper_vals_orig <- quantiles[, 2]
            }
            
            width_vals_orig <- upper_vals_orig - lower_vals_orig
            
            block_lower <- rep(NA_real_, nrow(block_vals))
            block_upper <- rep(NA_real_, nrow(block_vals))
            block_width <- rep(NA_real_, nrow(block_vals))
            block_lower[valid_rows] <- lower_vals_orig
            block_upper[valid_rows] <- upper_vals_orig
            block_width[valid_rows] <- width_vals_orig
          } else {
            block_lower <- rep(NA_real_, nrow(block_vals))
            block_upper <- rep(NA_real_, nrow(block_vals))
            block_width <- rep(NA_real_, nrow(block_vals))
          }
          
          terra::writeValues(out_lower, block_lower, b$row[i], b$nrows[i])
          terra::writeValues(out_upper, block_upper, b$row[i], b$nrows[i])
          terra::writeValues(out_width, block_width, b$row[i], b$nrows[i])
        }
        terra::writeStop(out_lower)
        terra::writeStop(out_upper)
        terra::writeStop(out_width)
        
        cat("\n")
        message(paste(
          "    - Uncertainty rasters for",
          year,
          "saved successfully."
        ))
      }, finally = {
        terra::readStop(predictor_stack_final)
      })
    }
    rm(bootstrap_models) # Free up memory
    gc()
  }
}, finally = {
  stopCluster(cl)
  message("Parallel processing cluster stopped for Prediction.")
})


#### 6. Visualization and Finalization ####
print("\n--- 6. Visualization and Finalization ---")

#---- 6.1. Summarize Prediction Rasters ----
# This summarizes the .tif files generated in Part 5
all_prediction_files <-
  list.files(dir_predictions, pattern = "\\.tif$", full.names = TRUE)

if (length(all_prediction_files) > 0) {
  prediction_summary <- map_dfr(all_prediction_files, function(f) {
    r <- terra::rast(f)
    vals <- values(r, na.rm = TRUE)
    tibble(
      File = basename(f),
      Mean = mean(vals),
      SD = sd(vals),
      Min = min(vals),
      Median = median(vals),
      Max = max(vals)
    )
  })
  
  output_excel_path <- file.path(dir_output_base,
                                 paste0("XGBoost_Analysis_Robust_", target_resolution, ".xlsx"))
  
  if (file.exists(output_excel_path)) {
    wb <- loadWorkbook(output_excel_path)
    
    if (!("Prediction_Summary" %in% names(wb))) {
      message("Adding Prediction_Summary sheet to the existing Excel file.")
      addWorksheet(wb, "Prediction_Summary")
      writeData(wb, "Prediction_Summary", prediction_summary, borders = "all")
      setColWidths(
        wb,
        sheet = "Prediction_Summary",
        cols = 1:ncol(prediction_summary),
        widths = "auto"
      )
      saveWorkbook(wb, output_excel_path, overwrite = TRUE)
    } else {
      message("Prediction_Summary sheet already exists. Skipping update.")
    }
  }
} else {
  message("No prediction maps found to summarize.")
}

#---- 6.2. Create Correlation Plot ----
corr_plot_path <-
  file.path(dir_evaluation, "predictor_correlation_plot.png")
if (!file.exists(corr_plot_path)) {
  message("Creating predictor correlation plot...")
  renamed_vars <- unname(full_name_map)
  
  corr_matrix <- master_data %>%
    select(any_of(
      c(
        renamed_vars,
        "Elevation",
        "EVI",
        "LST_Day",
        "LST_Night",
        "NDVI",
        "SWIR"
      )
    )) %>%
    select(where(is.numeric)) %>%
    cor(use = "pairwise.complete.obs")
  
  corr_plot <- ggcorrplot::ggcorrplot(
    corr_matrix,
    hc.order = TRUE,
    type = "lower",
    lab = TRUE,
    lab_size = 2,
    method = "circle",
    colors = c("#6D9EC1", "white", "#E46726")
  ) +
    labs(title = "Predictor Correlation Matrix") +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  ggsave(
    corr_plot_path,
    plot = corr_plot,
    width = 9,
    height = 8,
    dpi = 600
  )
  message(paste("Correlation plot saved to:", corr_plot_path))
} else {
  message("Correlation plot already exists. Skipping.")
}

print(paste("Part 2 finished at:", Sys.time()))


#### 7. SHAP Analysis Re-generation and Visualization ####
print("\n--- 7. SHAP Analysis Re-generation and Visualization (for Figure 6) ---")
# This section can be run independently as long as the final model
# files (.rds) from Part 1 exist.

all_shap_plot_data <- list()
all_shap_importance <- list()
shap_data_available <- FALSE

#---- 7.1. Load or Generate SHAP Data ----
# This loop checks for existing SHAP data (.rds).
# If not found, it generates it from the saved definitive models.
for (response_var in response_vars) {
  cat(paste("\n--- Checking SHAP data for:", response_var, "---\n"))
  
  model_path <-
    file.path(dir_models,
              paste0("definitive_model_", response_var, ".rds"))
  features_path <-
    file.path(dir_models,
              paste0("definitive_features_", response_var, ".rds"))
  shap_data_path <-
    file.path(dir_evaluation,
              paste0("shap_plot_data_", response_var, ".rds"))
  shap_imp_path <-
    file.path(dir_evaluation,
              paste0("shap_importance_", response_var, ".rds"))
  
  if (!file.exists(model_path) || !file.exists(features_path)) {
    cat(
      "  - Definitive model or feature list not found. Skipping SHAP analysis for this variable.\n"
    )
    next
  }
  
  if (file.exists(shap_data_path) && file.exists(shap_imp_path)) {
    cat("  - Loading pre-existing SHAP data from disk...\n")
    all_shap_plot_data[[response_var]] <- readRDS(shap_data_path)
    all_shap_importance[[response_var]] <- readRDS(shap_imp_path)
    shap_data_available <- TRUE
  } else {
    cat("  - SHAP data file not found. Generating from the definitive model...\n")
    
    # Load model and data
    definitive_model <- readRDS(model_path)
    final_features <- readRDS(features_path)
    current_data <-
      if (response_var == "Total_Abundance")
        master_data
    else
      master_data %>% filter(!is.na(.data[[response_var]]))
    final_predictors_df <-
      current_data %>% select(all_of(final_features))
    
    # Calculate SHAP values (can be slow)
    shap_values <- predict(definitive_model,
                           data.matrix(final_predictors_df),
                           predcontrib = TRUE)
    
    # Calculate mean absolute SHAP (importance)
    mean_abs_shap <- as.data.frame(shap_values) %>%
      select(-BIAS) %>%
      abs() %>%
      colMeans() %>%
      sort(decreasing = TRUE) %>%
      tibble::enframe(name = "Feature", value = "Mean_Abs_SHAP")
    
    # Prepare data for beeswarm plot
    shap_long <- as.data.frame(shap_values) %>%
      select(-BIAS) %>%
      mutate(row_id = row_number()) %>%
      pivot_longer(-row_id, names_to = "feature", values_to = "shap_value")
    feature_values_long <-
      as.data.frame(final_predictors_df) %>%
      mutate(row_id = row_number()) %>%
      pivot_longer(-row_id, names_to = "feature", values_to = "feature_value")
    shap_plot_data <-
      left_join(shap_long, feature_values_long, by = c("row_id", "feature"))
    shap_plot_data_scaled <- shap_plot_data %>%
      group_by(feature) %>%
      mutate(feature_value_scaled = as.numeric(scale(feature_value))) %>%
      ungroup() %>%
      mutate(feature_value_scaled = if_else(is.nan(feature_value_scaled), 0, feature_value_scaled))
    
    # Save generated data to avoid re-calculation
    saveRDS(shap_plot_data_scaled, shap_data_path)
    saveRDS(mean_abs_shap, shap_imp_path)
    cat("  - New SHAP data generated and saved to disk.\n")
    
    all_shap_plot_data[[response_var]] <- shap_plot_data_scaled
    all_shap_importance[[response_var]] <- mean_abs_shap
    shap_data_available <- TRUE
  }
}

#---- 7.2. Create Combined SHAP Plots (for Figure 6) ----
if (shap_data_available) {
  print(
    "\n--- 7.2. Creating combined SHAP plot (Individually Ranked, custom labels) ---"
  )
  
  #--- Define full names for plot labels ---
  plot_label_map <- c(
    "LST_Day"   = "Land Surface Temperature (Day)",
    "LST_Night" = "Land Surface Temperature (Night)",
    "NDVI"      = "Normalized Difference Vegetation Index",
    "EVI"       = "Enhanced Vegetation Index",
    "SWIR"      = "Short-Wave Infrared Reflectance"
    # Other variables will use their names from 'full_name_map'
  )
  
  #--- Helper function for y-axis label colors ---
  feature_colors <- function(features) {
    # Define variable groups and colors
    host_vars <-
      c(
        "Host Abundance (Aggregate)",
        "Host Abundance (Pinus densiflora)",
        "Host Abundance (Pinus thunbergii)"
      )
    lst_vars <- c("LST_Day", "LST_Night")
    spectral_vars <- c("NDVI", "EVI", "SWIR")
    elevation_var <- "Elevation"
    
    sapply(features, function(feature) {
      if (feature %in% host_vars)
        return("darkgreen") # Host
      if (feature %in% lst_vars)
        return("orange")  # Land Surface Temp
      if (feature %in% spectral_vars)
        return("#00008B") # Dark Blue (Vegetation Vigor)
      if (feature == elevation_var)
        return("brown")     # Topography
      return("black")    # Default (Bioclimatic)
    }, USE.NAMES = FALSE)
  }
  
  # Calculate global range for a consistent color legend
  all_shap_dfs <- list()
  for (rv in response_vars) {
    if (!is.null(all_shap_plot_data[[rv]])) {
      all_shap_dfs[[rv]] <- all_shap_plot_data[[rv]]
    }
  }
  combined_shap_data <- bind_rows(all_shap_dfs)
  global_scaled_range <-
    range(combined_shap_data$feature_value_scaled, na.rm = TRUE)
  
  tryCatch({
    # Generate plots for Top 10, 15, and 20 features
    for (n_top in c(10, 15, 20)) {
      cat(paste0(
        "  - Creating combined SHAP plot for Top ",
        n_top,
        " features...\n"
      ))
      
      plot_list <- list() # To store individual plots (a, b, c)
      plot_idx <- 1      # Index for plot labels (a), (b), (c)
      
      # Create one plot per response variable
      for (response_var in response_vars) {
        shap_plot_data_scaled <- all_shap_plot_data[[response_var]]
        shap_importance <- all_shap_importance[[response_var]]
        
        if (is.null(shap_plot_data_scaled) ||
            is.null(shap_importance)) {
          next
        }
        
        # --- 1. Prepare plot data ---
        top_features_local <-
          shap_importance %>% head(n_top) %>% pull(Feature)
        feature_order <- rev(top_features_local)
        
        plot_data_top_local <- shap_plot_data_scaled %>%
          filter(feature %in% top_features_local) %>%
          mutate(feature = factor(feature, levels = feature_order))
        
        bar_plot_data <- shap_importance %>%
          filter(Feature %in% top_features_local) %>%
          mutate(Feature = factor(Feature, levels = feature_order))
        
        #--- Create rich text y-axis labels (colored) ---
        y_axis_labels_markdown <-
          glue::glue(
            "<span style='color:{feature_colors(feature_order)};'>{recode(feature_order, !!!plot_label_map)}</span>"
          )
        
        # --- 2. Create Bar Plot (mean |SHAP value|) ---
        plot_title <- paste0(
          "(",
          letters[plot_idx],
          ") ",
          recode(
            response_var,
            "Total_Abundance" = "Abundance",
            "Onset_Week" = "Onset Week",
            "Peak_Week" = "Peak Week"
          )
        )
        
        p_bar <-
          ggplot(bar_plot_data, aes(x = Mean_Abs_SHAP, y = Feature)) +
          geom_col(fill = "grey50") +
          scale_y_discrete(labels = y_axis_labels_markdown) +
          labs(x = "mean(|SHAP value|)",
               y = NULL,
               title = plot_title) +
          theme_bw(base_size = 14) +
          theme(
            plot.title = element_text(
              hjust = 0,
              size = rel(1.1),
              face = "bold"
            ),
            axis.text.y = ggtext::element_markdown(size = 11),
            panel.grid.minor = element_blank()
          )
        
        # --- 3. Create Beeswarm Plot (SHAP values) ---
        p_beeswarm <-
          ggplot(
            plot_data_top_local,
            aes(
              x = shap_value,
              y = feature,
              color = feature_value_scaled
            )
          ) +
          ggbeeswarm::geom_quasirandom(
            groupOnX = FALSE,
            alpha = 0.7,
            orientation = "y",
            size = 1
          ) +
          geom_vline(
            xintercept = 0,
            color = "black",
            linetype = "dashed"
          ) +
          scale_color_gradient2(
            low = "blue",
            mid = "gray",
            high = "red",
            midpoint = 0,
            name = "Feature\nvalue",
            limits = global_scaled_range,
            breaks = global_scaled_range,
            labels = c("Low", "High"),
            oob = scales::squish
          ) +
          labs(x = "SHAP value", y = NULL) +
          theme_bw(base_size = 14) +
          theme(
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            panel.grid.minor = element_blank()
          )
        
        # --- 4. Combine bar and beeswarm plots ---
        plot_list[[response_var]] <-
          p_bar + p_beeswarm + plot_layout(widths = c(0.5, 1))
        
        plot_idx <- plot_idx + 1
      }
      
      if (length(plot_list) > 0) {
        # --- 5. Combine all plots (a, b, c) and save ---
        final_plot <- wrap_plots(plot_list, ncol = 1, guides = "collect") &
          theme(
            legend.position = "right",
            legend.key.height = unit(1.5, "cm"),
            legend.key.width = unit(0.5, "cm")
          )
        
        ggsave(
          file.path(
            dir_evaluation,
            paste0(
              "definitive_shap_plot_customized_top",
              n_top,
              "_combined.png"
            )
          ),
          plot = final_plot,
          width = 9,
          height = 2 + (n_top * 3 * 0.15),
          dpi = 600
        )
      }
    }
    print("Customized combined SHAP plots generated and saved.")
  }, error = function(e) {
    print(paste(
      "An error occurred during combined SHAP plot generation:",
      e$message
    ))
  })
} else {
  print("\nNo SHAP data was available or generated, skipping combined plots.")
}

print("\n--- All pipelines finished. ---")