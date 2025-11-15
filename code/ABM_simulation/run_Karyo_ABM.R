#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(yaml)
  library(Rcpp)
  library(parallel)
})

## 1) Source core ABM utilities (read_config, init_simulation, run_abm_simulations, etc.)
source("/Users/taolee/Documents/GitHub/ABM_Ploidy/code/ABM_simulation/Karyo_ABM.R")  # Ensure path points to the correct Karyo_ABM.R

## 2) Read configuration file
# Modify this to your actual config.yaml path
cfg_path <- "/Users/taolee/Documents/GitHub/ABM_Ploidy/ReqiuredData/config.yaml"
cfg <- read_config(cfg_path)

## 3) Set C++ source path
cpp_path <- "/Users/taolee/Documents/GitHub/ABM_Ploidy/code/ABM_simulation/abm_core.cpp"

## 4) Set number of simulation steps
steps <- 4320L   # Adjustable

## 5) Loop over WGDr values (if provided) and run batch simulations
wgdr_values <- cfg$WGDr

# If WGDr is not specified in config, fall back to a single default run
if (is.null(wgdr_values) || length(wgdr_values) == 0L) {
  wgdr_values <- NA_real_
}

# Prepare result container
res_all <- vector("list", length(wgdr_values))

for (ii in seq_along(wgdr_values)) {
  cfg_i <- cfg

  # Set WGDr and corresponding subdirectory name
  if (is.na(wgdr_values[ii])) {
    # No explicit WGDr in cfg: use whatever is already in cfg and a default label
    wgdr_label <- "WGDr_default"
  } else {
    cfg_i$WGDr <- wgdr_values[ii]
    wgdr_label <- sprintf("WGDr_%g", wgdr_values[ii])
  }

  # Base output directory from config, then append WGDp and WGDr-specific subfolders
  base_outdir <- cfg$OutDir
  if (is.null(base_outdir) || is.na(base_outdir)) {
    stop("cfg$OutDir must be defined in the configuration file.")
  }
  wgp_label <- sprintf("WGDp_%g", cfg_i$WGDp)
  cfg_i$OutDir <- file.path(base_outdir, wgp_label, wgdr_label)
  dir.create(cfg_i$OutDir, recursive = TRUE, showWarnings = FALSE)

  message("Running simulations for ", wgdr_label,
          " (WGDr = ", ifelse(is.na(wgdr_values[ii]), "cfg default", wgdr_values[ii]),
          "), OutDir = ", cfg_i$OutDir)

  # Run the existing batch simulation function (unchanged)
  res_all[[ii]] <- run_abm_simulations(
    cfg      = cfg_i,
    cpp_path = cpp_path,
    steps    = steps
  )
}

# Optionally return the list of results (one entry per WGDr value) when sourced
invisible(res_all)
gc()