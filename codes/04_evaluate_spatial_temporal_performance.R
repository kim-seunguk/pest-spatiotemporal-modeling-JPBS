# -------------------------------------------------------------------------- #
# SCRIPT 4: Spatial vs. Temporal Predictive Performance
# -------------------------------------------------------------------------- #
#
# Author: Seunguk Kim
# Contact: adrenaline@snu.ac.kr
# Last Modified: 2026-06-08
#
# Description:
# This script decomposes the reported predictive accuracy into its spatial
# and temporal components, and visualizes observed vs. predicted values.
# It produces the data for Supplementary Table S6 and Supplementary
# Figure S4, added in response to the reviewers.
#
# Associated Manuscript:
# Bang, S., Kim, S., et al. "From Pest Traps to Management Maps..."
#
# Usage Notes:
# This script must be run AFTER '02_run_xgboost_pipeline.R'. It re-uses the
# tuned hyperparameters and predictor sets saved by that pipeline and refits
# models on data subsets; the final (definitive) models are not modified.
#
# -------------------------------------------------------------------------- #

# The evaluations are:
# 1. Pooled cross-validation: repeated random 85/15 splits on the pooled
#    2022-2023 data (reproduces the accuracy reported in the main text).
# 2. Within-year cross-validation: repeated random 85/15 splits within each
#    year, isolating the spatial component (abundance only, as the 2023
#    phenology sample is too small for repeated splitting).
# 3. Leave-one-year-out (LOYO): trains on one year and predicts the other,
#    assessing across-year transfer.
# 4. Observed vs. predicted: out-of-fold k-fold predictions for all three
#    response variables (Figure S4).
#

#### 0. Setup ####
print("--- 0. Setup ---")
print(paste("Pipeline started at:", Sys.time()))

#---- 0.1. Load Libraries ----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(xgboost)
  library(caret)      # For createDataPartition, createFolds, postResample
  library(patchwork)  # For combining the scatter-plot panels
})

#---- 0.2. Define File Paths ----
dir_output_base <- here("output")
dir_data_prep <- file.path(dir_output_base, "prepared_data")
dir_models <- file.path(dir_output_base, "models")
dir_evaluation <- file.path(dir_output_base, "evaluation")

if (!dir.exists(dir_evaluation)) {
  dir.create(dir_evaluation, recursive = TRUE)
}

#---- 0.3. Package Version Management (renv) ----
# This project uses 'renv' to ensure reproducibility.
# If you are running this for the first time, run:
# > renv::restore()

#---- 0.4. Analysis Settings ----
response_vars <- c("Total_Abundance", "Onset_Week", "Peak_Week")
N_REPEATS <- 500   # Repeats for the pooled and within-year random-split evaluation
K_FOLDS <- 10      # Folds for the out-of-fold observed-vs-predicted plot


#### 1. Data Loading ####
#
# Description:
# Loads the modeling dataset and the tuned hyperparameters / predictor sets
# saved by the main pipeline ('02_run_xgboost_pipeline.R').
#

print("--- 1. Loading Data and Tuned Models from Main Pipeline ---")

#---- 1.1. Load Modeling Dataset ----
master_data_path <- file.path(dir_data_prep, "master_data.rds")
if (!file.exists(master_data_path)) {
  stop(
    "Master data file (master_data.rds) not found. Please run '02_run_xgboost_pipeline.R' first."
  )
}
master_data <- readRDS(master_data_path)

#---- 1.2. Load Tuned Hyperparameters and Predictor Sets ----
best_params_list <- list()
features_list <- list()
for (response_var in response_vars) {
  params_path <-
    file.path(dir_models, paste0("definitive_best_params_", response_var, ".rds"))
  features_path <-
    file.path(dir_models, paste0("definitive_features_", response_var, ".rds"))
  if (!file.exists(params_path) || !file.exists(features_path)) {
    stop(
      paste0(
        "Tuned parameter/feature files for '",
        response_var,
        "' not found. Please run '02_run_xgboost_pipeline.R' first."
      )
    )
  }
  best_params_list[[response_var]] <- readRDS(params_path)
  features_list[[response_var]] <- readRDS(features_path)
}

#---- 1.3. Define Helpers ----
# Subset the data for a given response variable, following the same rule as
# the main pipeline (phenology rows with missing values are dropped).
get_response_data <- function(data, response_var) {
  if (response_var == "Total_Abundance") {
    data
  } else {
    data %>% filter(!is.na(.data[[response_var]]))
  }
}

# Fit a single XGBoost model, replicating the training call used in
# '02_run_xgboost_pipeline.R' (reg:squarederror; nrounds from the tuned set).
fit_xgb_model <- function(train_data, features, response_var, best_params) {
  predictor_matrix <-
    data.matrix(train_data %>% dplyr::select(all_of(features)))
  dtrain <- xgb.DMatrix(data = predictor_matrix,
                        label = train_data[[response_var]])
  params <- as.list(best_params)
  params$nrounds <- NULL # nrounds is specified in xgb.train
  params$objective <- "reg:squarederror"
  params$nthread <- 1
  xgb.train(
    params = params,
    data = dtrain,
    nrounds = best_params$nrounds,
    verbose = 0
  )
}

# Predict on new data using the model's feature set.
predict_xgb_model <- function(model, new_data, features) {
  predict(model, data.matrix(new_data %>% dplyr::select(all_of(features))))
}


#### 2. Spatial vs. Temporal Cross-Validation (Table S6) ####
#
# Description:
# Compares pooled, within-year, and leave-one-year-out cross-validation to
# separate the spatial and temporal components of predictive performance.
#

print("\n--- 2. Decomposing Predictive Performance (for Table S6) ---")

cv_results_list <- list()

for (response_var in response_vars) {
  cat(paste("\n\n===== Evaluating:", response_var, "=====\n"))
  current_data <- get_response_data(master_data, response_var)
  features <- features_list[[response_var]]
  best_params <- best_params_list[[response_var]]

  #---- 2.1. Pooled Repeated 85/15 Split (reproduces main-text accuracy) ----
  cat("-> Pooled cross-validation...\n")
  pooled_metrics <- map_dfr(seq_len(N_REPEATS), function(i) {
    set.seed(i)
    train_index <-
      createDataPartition(current_data[[response_var]], p = 0.85, list = FALSE)
    train_data <- current_data[train_index, ]
    test_data <- current_data[-train_index, ]
    model <- fit_xgb_model(train_data, features, response_var, best_params)
    predictions <- predict_xgb_model(model, test_data, features)
    as.data.frame(t(postResample(predictions, test_data[[response_var]])))
  })
  cv_results_list[[length(cv_results_list) + 1]] <- tibble(
    Response_Var = response_var,
    Evaluation = "Pooled",
    N = nrow(current_data),
    R2 = mean(pooled_metrics$Rsquared, na.rm = TRUE),
    RMSE = mean(pooled_metrics$RMSE, na.rm = TRUE)
  )

  #---- 2.2. Within-Year Repeated 85/15 Split (spatial component) ----
  # Restricted to abundance; the 2023 phenology sample is too small.
  if (response_var == "Total_Abundance") {
    for (yr in c(2022, 2023)) {
      cat(paste("-> Within-year cross-validation (", yr, ")...\n"))
      year_data <- current_data %>% filter(Year == yr)
      within_metrics <- map_dfr(seq_len(N_REPEATS), function(i) {
        set.seed(i)
        train_index <-
          createDataPartition(year_data[[response_var]], p = 0.85, list = FALSE)
        train_data <- year_data[train_index, ]
        test_data <- year_data[-train_index, ]
        model <- fit_xgb_model(train_data, features, response_var, best_params)
        predictions <- predict_xgb_model(model, test_data, features)
        as.data.frame(t(postResample(predictions, test_data[[response_var]])))
      })
      cv_results_list[[length(cv_results_list) + 1]] <- tibble(
        Response_Var = response_var,
        Evaluation = paste0("Within-year ", yr),
        N = nrow(year_data),
        R2 = mean(within_metrics$Rsquared, na.rm = TRUE),
        RMSE = mean(within_metrics$RMSE, na.rm = TRUE)
      )
    }
  }

  #---- 2.3. Leave-One-Year-Out (across-year transfer) ----
  for (direction in list(c(2022, 2023), c(2023, 2022))) {
    train_year <- direction[1]
    test_year <- direction[2]
    cat(paste("-> Leave-one-year-out (", train_year, "->", test_year, ")...\n"))
    train_data <- current_data %>% filter(Year == train_year)
    test_data <- current_data %>% filter(Year == test_year)
    model <- fit_xgb_model(train_data, features, response_var, best_params)
    predictions <- predict_xgb_model(model, test_data, features)
    loyo_metrics <- postResample(predictions, test_data[[response_var]])
    cv_results_list[[length(cv_results_list) + 1]] <- tibble(
      Response_Var = response_var,
      Evaluation = paste0("LOYO ", train_year, "->", test_year),
      N = nrow(test_data),
      R2 = unname(loyo_metrics["Rsquared"]),
      RMSE = unname(loyo_metrics["RMSE"])
    )
  }
}

#---- 2.4. Save and Display Decomposition Table ----
cv_decomposition <- bind_rows(cv_results_list)
output_table_path <-
  file.path(dir_evaluation, "cross_validation_decomposition.csv")
write.csv(cv_decomposition, output_table_path, row.names = FALSE)
print(paste("--- Cross-validation table saved to:", output_table_path, "---"))
print(cv_decomposition %>% mutate(across(where(is.numeric), ~ round(., 3))), n = 50)


#### 3. Observed vs. Predicted Scatter Plots (Figure S4) ####
#
# Description:
# Generates out-of-fold (k-fold) predictions and plots them against the
# observed values for each response variable.
#

print("\n--- 3. Generating Observed vs. Predicted Plots (for Figure S4) ---")

#---- 3.1. Generate Out-of-Fold Predictions ----
oof_list <- list()
for (response_var in response_vars) {
  cat(paste("-> Out-of-fold predictions for:", response_var, "\n"))
  current_data <- get_response_data(master_data, response_var)
  features <- features_list[[response_var]]
  best_params <- best_params_list[[response_var]]

  set.seed(42)
  folds <- createFolds(current_data[[response_var]], k = K_FOLDS)
  oof_predictions <- rep(NA_real_, nrow(current_data))
  for (fold_index in folds) {
    train_data <- current_data[-fold_index, ]
    test_data <- current_data[fold_index, ]
    model <- fit_xgb_model(train_data, features, response_var, best_params)
    oof_predictions[fold_index] <-
      predict_xgb_model(model, test_data, features)
  }

  oof_list[[response_var]] <- tibble(
    Response_Var = response_var,
    Year = factor(current_data$Year),
    Observed = current_data[[response_var]],
    Predicted = oof_predictions
  )
}
oof_data <- bind_rows(oof_list)

#---- 3.2. Create a Plotting Function ----
# Standardized observed-vs-predicted panel with a shared (square) axis range,
# a 1:1 reference line, a linear fit, and accuracy annotations.
create_obs_pred_panel <-
  function(data, plot_title, axis_limits, axis_breaks) {
    performance <- postResample(data$Predicted, data$Observed)
    stats_label <- sprintf("R² = %.2f\nRMSE = %.2f",
                           performance["Rsquared"], performance["RMSE"])

    ggplot(data, aes(x = Observed, y = Predicted)) +
      geom_abline(
        slope = 1,
        intercept = 0,
        linetype = "dashed",
        color = "grey40"
      ) +
      geom_point(aes(color = Year), alpha = 0.55, size = 1.2) +
      geom_smooth(
        method = "lm",
        se = FALSE,
        color = "black",
        linewidth = 0.6
      ) +
      annotate(
        "text",
        x = axis_limits[1],
        y = axis_limits[2],
        label = stats_label,
        hjust = 0,
        vjust = 1,
        size = 3
      ) +
      scale_x_continuous(limits = axis_limits, breaks = axis_breaks) +
      scale_y_continuous(limits = axis_limits, breaks = axis_breaks) +
      scale_color_brewer(palette = "Set1", name = "Year") +
      labs(title = plot_title, x = "Observed", y = "Predicted") +
      theme_bw(base_size = 12) +
      theme(
        aspect.ratio = 1,
        plot.title = element_text(face = "bold", hjust = 0, size = 10),
        plot.margin = margin(2, 2, 2, 2),
        panel.grid.minor = element_blank()
      )
  }

#---- 3.3. Build Panels with Matched Axes ----
# Onset and peak share a common week range; abundance uses its own log1p range.
abundance_data <- oof_data %>% filter(Response_Var == "Total_Abundance")
phenology_data <-
  oof_data %>% filter(Response_Var %in% c("Onset_Week", "Peak_Week"))

week_breaks <-
  pretty(range(c(phenology_data$Observed, phenology_data$Predicted), na.rm = TRUE), n = 5)
week_limits <-
  range(c(phenology_data$Observed, phenology_data$Predicted, week_breaks), na.rm = TRUE)
abundance_breaks <-
  pretty(range(c(abundance_data$Observed, abundance_data$Predicted), na.rm = TRUE), n = 5)
abundance_limits <-
  range(c(abundance_data$Observed, abundance_data$Predicted, abundance_breaks), na.rm = TRUE)

panel_abundance <- create_obs_pred_panel(
  abundance_data,
  "(a) Abundance (log1p)",
  abundance_limits,
  abundance_breaks
)
panel_onset <- create_obs_pred_panel(
  phenology_data %>% filter(Response_Var == "Onset_Week"),
  "(b) Onset week",
  week_limits,
  week_breaks
)
panel_peak <- create_obs_pred_panel(
  phenology_data %>% filter(Response_Var == "Peak_Week"),
  "(c) Peak week",
  week_limits,
  week_breaks
)

#---- 3.4. Combine and Save Final Figure ----
figure_s4 <- (panel_abundance | panel_onset | panel_peak) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

output_figure_path <-
  file.path(dir_evaluation, "figure_S4_observed_vs_predicted.png")
ggsave(
  filename = output_figure_path,
  plot = figure_s4,
  width = 8.5,
  height = 3.9,
  dpi = 600
)
print(paste("--- Observed vs. predicted figure saved to:", output_figure_path, "---"))

print("\n--- All supplementary validation analyses complete. ---")
