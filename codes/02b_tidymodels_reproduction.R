# -------------------------------------------------------------------------- #
# SCRIPT 2b (OPTIONAL): Reproduction of the Model Building with tidymodels
# -------------------------------------------------------------------------- #
#
# Author: Seunguk Kim
# Contact: adrenaline@snu.ac.kr
#
# Description:
# This is a standalone, optional companion to the main pipeline. The main
# pipeline tunes the XGBoost models with caret (`caret::train(method =
# "xgbTree")`, script 02), which was run with xgboost < 2.0. caret's xgbTree
# method is not compatible with xgboost >= 2.0, so script 02 cannot be re-run
# on a current xgboost installation.
#
# This script reproduces the same model-building step (hyperparameter tuning,
# repeated-split performance evaluation, and final model fit) using the
# actively maintained tidymodels / parsnip framework, which runs on current
# xgboost. It reproduces the reported predictive performance to within about
# 0.01 R-squared. Because gradient boosting has a flat optimum, the exact
# selected hyperparameters and metrics vary slightly with the framework, the
# cross-validation fold randomization, and the xgboost version; this is
# expected and does not change the conclusions.
#
# It reuses the same prepared data (output/prepared_data/master_data.rds,
# created by script 02) and the same predictor set, grid, cross-validation
# scheme (5-fold x 2 repeats), and 500-repeat 85/15 evaluation as script 02.
#
# Usage Notes:
# Run from the project root, after script 02 has created master_data.rds.
# Requires the tidymodels packages (install.packages("tidymodels")).
#
# -------------------------------------------------------------------------- #

#### 0. Setup ####
print("--- 0. Setup ---")
print(paste("Pipeline started at:", Sys.time()))

#---- 0.1. Load Libraries ----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(parsnip)
  library(tune)
  library(rsample)
  library(yardstick)
  library(dials)
  library(workflows)
  library(xgboost)
  library(caret)      # createDataPartition / postResample only (no xgbTree)
})

#---- 0.2. Define File Paths ----
dir_output_base <- here("output")
dir_data_prep <- file.path(dir_output_base, "prepared_data")
dir_models <- file.path(dir_output_base, "models")
dir_evaluation <- file.path(dir_output_base, "evaluation")
if (!dir.exists(dir_evaluation)) dir.create(dir_evaluation, recursive = TRUE)

#---- 0.3. Analysis Settings ----
response_vars <- c("Total_Abundance", "Onset_Week", "Peak_Week")
N_REPEATS <- 500   # 85/15 repeated-split evaluation, matching script 02

# Non-predictor columns removed before modeling, matching script 02 (line 305).
non_predictor_cols <- c("Year", "Area", "Region", "trap_no", "lon", "lat")

# Published headline performance (manuscript), shown as a reference only.
published_R2   <- c(Total_Abundance = 0.78, Onset_Week = 0.91, Peak_Week = 0.90)
published_RMSE <- c(Total_Abundance = 0.78, Onset_Week = 0.91, Peak_Week = 1.11)

# The tuning grid of script 02, expressed in parsnip parameter names. parsnip
# maps trees -> nrounds, tree_depth -> max_depth, learn_rate -> eta,
# loss_reduction -> gamma, min_n -> min_child_weight.
tuning_grid <- expand.grid(
  trees          = c(100, 200),
  tree_depth     = c(2, 3, 5),
  learn_rate     = c(0.01, 0.05, 0.1),
  loss_reduction = c(0.1, 1),
  min_n          = c(1, 5)
)


#### 1. Data Loading ####
print("--- 1. Loading Prepared Data ---")
master_data_path <- file.path(dir_data_prep, "master_data.rds")
if (!file.exists(master_data_path)) {
  stop("master_data.rds not found. Please run '02_run_xgboost_pipeline.R' first (it writes master_data.rds before the tuning step).")
}
master_data <- readRDS(master_data_path)

get_response_data <- function(data, response_var) {
  if (response_var == "Total_Abundance") data
  else data %>% filter(!is.na(.data[[response_var]]))
}

# Build the modeling frame for a response: drop non-predictor and response
# columns, then rename predictors to syntactic names (same values and order as
# script 02's data.matrix; xgboost is column-name agnostic, and renaming avoids
# formula issues with spaces/parentheses in the original column names).
build_model_frame <- function(data, response_var) {
  predictors <- data %>%
    dplyr::select(-dplyr::any_of(c(non_predictor_cols, response_vars)))
  colnames(predictors) <- sprintf("p%02d", seq_len(ncol(predictors)))
  cbind(data.frame(y = data[[response_var]]), predictors)
}

# A boost_tree spec with colsample_bytree and subsample fixed at 0.8. parsnip
# maps mtry -> colsample_bytree and sample_size -> subsample; counts = FALSE
# makes mtry a proportion, so mtry = 0.8 corresponds to colsample_bytree = 0.8.
make_spec <- function(params = NULL) {
  if (is.null(params)) {
    boost_tree(
      mode = "regression", mtry = 0.8, sample_size = 0.8,
      trees = tune(), tree_depth = tune(), learn_rate = tune(),
      loss_reduction = tune(), min_n = tune()
    ) %>% set_engine("xgboost", counts = FALSE, nthread = 1)
  } else {
    boost_tree(
      mode = "regression", mtry = 0.8, sample_size = 0.8,
      trees = params$trees, tree_depth = params$tree_depth,
      learn_rate = params$learn_rate, loss_reduction = params$loss_reduction,
      min_n = params$min_n
    ) %>% set_engine("xgboost", counts = FALSE, nthread = 1)
  }
}


#### 2. Tuning, Evaluation, and Final Fit ####
performance_list <- list()
best_params_list <- list()

for (response_var in response_vars) {
  cat(paste("\n\n===== Modeling:", response_var, "=====\n"))
  current_data <- get_response_data(master_data, response_var)
  model_df <- build_model_frame(current_data, response_var)

  #---- 2.1. Hyperparameter Tuning (5-fold x 2 repeats, RMSE) ----
  cat("-> Tuning hyperparameters with tidymodels...\n")
  workflow_spec <- workflow() %>%
    add_formula(y ~ .) %>%
    add_model(make_spec())
  set.seed(42)
  folds <- vfold_cv(model_df, v = 5, repeats = 2)
  tune_results <- tune_grid(
    workflow_spec, resamples = folds, grid = tuning_grid,
    metrics = metric_set(rmse, rsq),
    control = control_grid(verbose = FALSE, save_pred = FALSE)
  )
  best_params <- select_best(tune_results, metric = "rmse")
  best_params_list[[response_var]] <- best_params %>%
    mutate(Response_Var = response_var, .before = 1)

  #---- 2.2. Robust Performance (500 repeats of an 85/15 split) ----
  cat("-> Evaluating performance over", N_REPEATS, "repeated splits...\n")
  final_spec <- make_spec(best_params)
  r2_vec <- numeric(N_REPEATS)
  rmse_vec <- numeric(N_REPEATS)
  for (i in seq_len(N_REPEATS)) {
    set.seed(i)
    train_index <- caret::createDataPartition(model_df$y, p = 0.85, list = FALSE)
    fitted <- fit(final_spec, y ~ ., data = model_df[train_index, ])
    preds <- predict(fitted, model_df[-train_index, ])$.pred
    metrics <- caret::postResample(preds, model_df$y[-train_index])
    r2_vec[i] <- metrics["Rsquared"]
    rmse_vec[i] <- metrics["RMSE"]
  }
  performance_list[[response_var]] <- tibble(
    Response_Var = response_var,
    R2_tidymodels = mean(r2_vec, na.rm = TRUE),
    RMSE_tidymodels = mean(rmse_vec, na.rm = TRUE),
    R2_published_reference = published_R2[[response_var]],
    RMSE_published_reference = published_RMSE[[response_var]]
  )

  #---- 2.3. Fit and Save the Final Model ----
  final_model <- fit(final_spec, y ~ ., data = model_df)
  saveRDS(final_model,
          file.path(dir_models, paste0("tidymodels_model_", response_var, ".rds")))
}


#### 3. Save and Display ####
best_params_table <- bind_rows(best_params_list)
performance_table <- bind_rows(performance_list)

write.csv(best_params_table,
          file.path(dir_evaluation, "tidymodels_best_params.csv"), row.names = FALSE)
write.csv(performance_table,
          file.path(dir_evaluation, "tidymodels_reproduction_performance.csv"), row.names = FALSE)

cat("\n=================== Selected hyperparameters (tidymodels) ===================\n")
print(best_params_table)
cat("\n=================== Performance vs published reference ======================\n")
print(performance_table %>% mutate(across(where(is.numeric), ~ round(., 3))))

print("\n--- Reproduction with tidymodels complete. ---")
