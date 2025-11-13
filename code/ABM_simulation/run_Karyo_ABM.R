#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(yaml)
  library(Rcpp)
  library(parallel)
})

## 1) Source core ABM utilities (read_config, init_simulation, run_abm_simulations, etc.)
source("/Users/4482173/Documents/GitHub/ABM_Ploidy/code/ABM_simulation/Karyo_ABM.R")  # Ensure path points to the correct Karyo_ABM.R

## 2) Read configuration file
# Modify this to your actual config.yaml path
cfg_path <- "/Users/4482173/Documents/GitHub/ABM_Ploidy/ReqiuredData/config.yaml"
cfg <- read_config(cfg_path)

## 3) Set C++ source path
cpp_path <- "/Users/4482173/Documents/GitHub/ABM_Ploidy/code/ABM_simulation/abm_core.cpp"

## 4) Set number of simulation steps
steps <- 1440L   # Adjustable

## 5) Run batch simulations
res_list <- run_abm_simulations(
  cfg      = cfg,
  cpp_path = cpp_path,
  steps    = steps
)

