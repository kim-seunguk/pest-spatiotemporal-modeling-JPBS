# -------------------------------------------------------------------------- #
# SCRIPT 1: Phenology Data Processing Pipeline
# -------------------------------------------------------------------------- #
#
# Author: Seunguk Kim
# Contact: adrenaline@snu.ac.kr
# Last Modified: 2025-10-21
#
# Description:
# This script processes raw biweekly trap data (occurrence*.csv) to
# calculate the final phenological metrics (Onset_Week, Peak_Week)
# used in the main XGBoost modeling pipeline. It corresponds to the
# "Estimation of phenological metrics" section of the manuscript.
#
# Associated Manuscript:
# Bang, S., Kim, S., et al. "From Pest Traps to Management Maps..."
#
# Usage Notes:
# This is the FIRST script in the analysis pipeline. It must be run
# before '03_run_xgboost_pipeline.R'.
#
# -------------------------------------------------------------------------- #

# Pipeline Flow:
# 0. Setup: Load libraries and define functions.
# 1. Calculate Phenology Metrics: Load raw data, fit Johnson SB
#    distributions, validate models, and save the final clean output.
#

#### 0. Setup ####
#### 0. Setup ####
print("--- 0. Setup ---")
print(paste("Pipeline started at:", Sys.time()))

#---- 0.1. Load Libraries ----
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(SuppDists) # For Johnson SB distribution
  library(sf)        # Used in Part 5 (now removed), but kept if other scripts depend on it
})

#---- 0.2. Global Path and File Configuration ----
# Input data path
dir_raw_data <- here("data", "occurrence")

# Standardized output directory structure
dir_output_base <- here("output", "phenology")
dir_output_data <- file.path(dir_output_base, "data")

# Create all output directories
dir.create(dir_output_data, showWarnings = FALSE, recursive = TRUE)

# Paths for Part 1 outputs
path_phenology_raw_data <- file.path(dir_output_data, "phenology_data.rds")
path_phenology_clean <- file.path(dir_output_data, "phenology_metrics_clean.csv")
path_phenology_all_metrics <- file.path(dir_output_data, "phenology_metrics.csv")
path_phenology_all_attempts <- file.path(dir_output_data, "phenology_metrics_all_attempts.csv")

#---- 0.3. Package Version Management (renv) ----
# This project uses 'renv' to ensure reproducibility.
# If you are running this for the first time, run:
# > renv::restore()
# This will install all required packages at their correct versions
# from the 'renv.lock' file.

#---- 0.4. Common Function Definitions ----

#' Fit Johnson SB Distribution
#'
#' Fits a Johnson SB cumulative distribution function to observed
#' cumulative proportion data using non-linear least squares (nls).
fit_johnson_sb <- function(df) {
  df_for_fitting <- df %>% filter(cumulative_proportion > 0 & cumulative_proportion < 1)
  
  # Require at least 4 data points between 0 and 1 for a robust fit
  if (nrow(df_for_fitting) < 4) {
    return(tibble(gamma = NA, delta = NA, xi = NA, lambda = NA, r_squared = NA))
  }
  
  min_time <- min(df_for_fitting$time_variable, na.rm = TRUE)
  max_time <- max(df_for_fitting$time_variable, na.rm = TRUE)
  range_time <- max_time - min_time
  
  fit <- try(nls(
    cumulative_proportion ~ pnorm(gamma + delta * log((time_variable - xi) / (lambda - (time_variable - xi)))),
    data = df_for_fitting,
    start = list(gamma = 0, delta = 1, xi = min_time - 1, lambda = range_time + 2),
    algorithm = "port",
    lower = list(gamma = -7, delta = 1e-6, xi = -Inf, lambda = range_time + 1e-6),
    upper = list(gamma = 7, delta = 10, xi = min_time - 1e-6, lambda = Inf),
    control = nls.control(maxiter = 100, warnOnly = TRUE)
  ), silent = TRUE)
  
  if (inherits(fit, "nls")) {
    observed <- df_for_fitting$cumulative_proportion
    predicted <- predict(fit)
    rss <- sum((observed - predicted) ^ 2)
    tss <- sum((observed - mean(observed)) ^ 2)
    r_squared <- 1 - (rss / tss)
    params <- as_tibble(as.list(coef(fit)))
    params$r_squared <- r_squared
    return(params)
  } else {
    return(tibble(gamma = NA, delta = NA, xi = NA, lambda = NA, r_squared = NA))
  }
}

#' Get Mode of Johnson SB Distribution
#'
#' Finds the peak (mode) of the Johnson SB probability density function
#' by numerical optimization.
get_mode <- function(gamma, delta, xi, lambda) {
  opt <- try(optimize(
    f = function(x) dJohnson(x, parms = list(gamma = gamma, delta = delta, xi = xi, lambda = lambda, type = "SB")),
    interval = c(xi, xi + lambda),
    maximum = TRUE
  ), silent = TRUE)
  
  if (inherits(opt, "try-error")) return(NA)
  return(opt$maximum)
}

#' Safe Wrapper for pJohnson
#'
#' Handles edge cases for the cumulative distribution function
#' where x is outside the bounds [xi, xi + lambda].
safe_pJohnson <- function(x, parms) {
  sapply(x, function(val) {
    if (is.na(val) || is.na(parms$xi) || is.na(parms$lambda)) return(NA)
    if (val <= parms$xi) return(0)
    if (val >= (parms$xi + parms$lambda)) return(1)
    return(pJohnson(val, parms = parms))
  })
}


cat("--- Phenology Data Processing Pipeline Started ---\n")

#### Part 1: Calculate Phenology Metrics ####
if (file.exists(path_phenology_clean) && file.exists(path_phenology_all_metrics)) {
  
  cat("✅ Part 1: Phenology metric files already exist. Skipping.\n")
  
} else {
  
  cat("🚀 Part 1: Starting phenology analysis...\n")
  
  #---- 1.1. Prepare Data ----
  # NOTE: Ensure 'occurrence2022.csv' and 'occurrence2023.csv'
  # are saved with UTF-8 encoding.
  occ2022csv <- read.csv(file.path(dir_raw_data, "occurrence2022.csv"))
  occ2023csv <- read.csv(file.path(dir_raw_data, "occurrence2023.csv"))
  
  # Convert wide date columns (e.g., "X2022.1.25..2.8") to Julian day ("D025")
  data2022_mod <- occ2022csv %>%
    rename_with(.cols = starts_with("X"),
                .fn = ~ sapply(.x, function(col_name) {
                  date_parts <- str_extract_all(col_name, "\\d+")[[1]]
                  start_date <- as.Date(paste("2022", date_parts[1], date_parts[2], sep = "-"))
                  paste0("D", format(start_date, "%j"))
                }, USE.NAMES = FALSE)) %>%
    dplyr::select(Area, Region, trap_no, lon, lat, starts_with("D"))
  
  data2023_mod <- occ2023csv %>%
    rename_with(.cols = starts_with("X"),
                .fn = ~ sapply(.x, function(col_name) {
                  date_parts <- str_extract_all(col_name, "\\d+")[[1]]
                  start_date <- as.Date(paste("2023", date_parts[1], date_parts[2], sep = "-"))
                  paste0("D", format(start_date, "%j"))
                }, USE.NAMES = FALSE)) %>%
    dplyr::select(Area, Region, trap_no, lon, lat, starts_with("D"))
  
  data_all <- bind_rows("2022" = data2022_mod, "2023" = data2023_mod, .id = "Year")
  
  # Convert to long format
  plot_data_trap_level <- data_all %>%
    pivot_longer(cols = starts_with("D"), names_to = "period", values_to = "captures")
  
  # Calculate cumulative proportion for fitting
  phenology_data <- plot_data_trap_level %>%
    filter(!is.na(captures)) %>%
    mutate(julian_day = as.numeric(str_remove(period, "D"))) %>%
    group_by(Year, Area, Region, trap_no) %>%
    filter(sum(captures) > 0) %>% # Keep only sites with captures
    mutate(
      lon = first(lon),
      lat = first(lat),
      total_abundance = sum(captures),
      cumulative_captures = cumsum(captures),
      cumulative_proportion = cumulative_captures / total_abundance,
      time_variable = julian_day / 7 # Convert Julian day to week of year
    ) %>%
    ungroup()
  
  #---- 1.2. Run Phenology Analysis ----
  # Count number of positive capture intervals per site
  phenology_data_with_counts <- phenology_data %>%
    group_by(Year, Area, Region, trap_no) %>%
    summarise(n_obs = sum(captures > 0), .groups = 'drop')
  
  # [cite_start]Define minimum observations needed for modeling (as per manuscript [cite: 131])
  MIN_OBS_FOR_MODEL <- 5
  sites_for_model <- phenology_data_with_counts %>% filter(n_obs >= MIN_OBS_FOR_MODEL)
  sites_for_direct <- phenology_data_with_counts %>% filter(n_obs < MIN_OBS_FOR_MODEL)
  
  data_for_model <- phenology_data %>% semi_join(sites_for_model, by = c("Year", "Area", "Region", "trap_no"))
  
  phenology_from_model <- tibble()
  phenology_metrics_unfiltered <- tibble()
  
  if (nrow(data_for_model) > 0) {
    # Nest data and apply the fitting function to each site
    fitted_params <- data_for_model %>%
      group_by(Year, Area, Region, trap_no) %>%
      nest() %>%
      mutate(params = map(data, fit_johnson_sb)) %>%
      unnest(params) %>%
      dplyr::select(-data)
    
    # Get metadata (abundance, location) for all sites
    abundance_lookup <- data_for_model %>%
      group_by(Year, Area, Region, trap_no) %>%
      summarise(lon = first(lon), lat = first(lat), Total_Abundance = first(total_abundance), .groups = 'drop')
    
    # Join abundance to ALL fit attempts (for diagnostics)
    all_attempts_data <- fitted_params %>%
      left_join(abundance_lookup, by = c("Year", "Area", "Region", "trap_no"))
    
    # Calculate metrics (Onset/Peak) from all successful parameter fits
    phenology_metrics_unfiltered <- all_attempts_data %>%
      filter(!is.na(gamma)) %>% # Filter for successful nls fits
      rowwise() %>%
      mutate(onset = qJohnson(0.05, parms = list(gamma = gamma, delta = delta, xi = xi, lambda = lambda, type = "SB")),
             peak = get_mode(gamma, delta, xi, lambda)) %>%
      ungroup()
    
    # [cite_start]Apply strict validation protocol (as per manuscript [cite: 137-141])
    phenology_from_model <- phenology_metrics_unfiltered %>%
      filter(lambda >= 2, xi >= 0, delta <= 8, delta >= 0.71, onset >= 0, peak >= 0, peak >= onset) %>%
      mutate(method = "model") %>%
      dplyr::select(Year, Area, Region, trap_no, lon, lat, Total_Abundance, onset, peak, method)
  }
  
  # Identify sites that were modeled but failed validation
  successful_sites <- phenology_from_model %>% distinct(Year, Area, Region, trap_no)
  unfit_sites <- sites_for_model %>% anti_join(successful_sites, by = c("Year", "Area", "Region", "trap_no"))
  
  # Combine sites that failed validation with sites that had too few obs
  sites_for_direct_combined <- bind_rows(sites_for_direct, unfit_sites)
  data_for_direct_combined <- phenology_data %>% semi_join(sites_for_direct_combined, by = c("Year", "Area", "Region", "trap_no"))
  
  # Calculate metrics directly (min day, max capture day) for these sites
  phenology_from_direct <- tibble()
  if (nrow(data_for_direct_combined) > 0) {
    phenology_from_direct <- data_for_direct_combined %>%
      group_by(Year, Area, Region, trap_no) %>%
      summarise(lon = first(lon), lat = first(lat), first_day = min(julian_day[captures > 0]),
                peak_day = julian_day[which.max(captures)][1], Total_Abundance = first(total_abundance), .groups = 'drop') %>%
      mutate(onset = first_day / 7, peak = peak_day / 7, method = "direct") %>%
      dplyr::select(Year, Area, Region, trap_no, lon, lat, Total_Abundance, onset, peak, method)
  }
  
  # Combine all results
  phenology_metrics <- bind_rows(phenology_from_model, phenology_from_direct) %>%
    arrange(Year, Area, Region, trap_no)
  
  #---- 1.3. Save Final Results ----
  # Save the long-format processed data
  saveRDS(phenology_data, file = path_phenology_raw_data)
  
  # Create the FINAL clean file for the XGBoost model
  # This file contains ONLY the validated model results
  phenology_metrics_clean <- phenology_metrics %>%
    filter(method == "model") %>%
    rename(First_Occurrence_5pct = onset, Peak_Occurrence_Mode = peak) %>%
    dplyr::select(Year, Area, Region, trap_no, lon, lat, Total_Abundance, First_Occurrence_5pct, Peak_Occurrence_Mode)
  
  # Save the clean data
  write.csv(phenology_metrics_clean, file = path_phenology_clean, row.names = FALSE)
  
  # Save the diagnostic files (all attempts and all metrics)
  write.csv(phenology_metrics, file = path_phenology_all_metrics, row.names = FALSE)
  
  if(exists("all_attempts_data") && nrow(all_attempts_data) > 0) {
    write.csv(all_attempts_data, file = path_phenology_all_attempts, row.names = FALSE)
  } else {
    # Write an empty file if no models were attempted
    tibble() %>% write.csv(file = path_phenology_all_attempts, row.names = FALSE)
  }
  
  cat(paste("  -> Part 1 Complete: Phenology analysis results saved to:\n     ", dir_output_data, "\n"))
  cat(paste("  -> FINAL OUTPUT (for modeling):", path_phenology_clean, "\n"))
}

cat("--- Phenology Pipeline Finished ---\n")