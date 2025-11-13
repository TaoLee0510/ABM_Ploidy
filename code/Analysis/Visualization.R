# ===== Setup =====
# Load packages
suppressPackageStartupMessages({
  library(tidyverse)
  library(fs)
  library(stringr)
  library(viridis)   # for color scale
  library(magick)   # for GIF assembly
})

# ---- User-configurable paths ----
# Top-level directory where your simulation outputs are written.
# E.g., "data" contains subfolders like Coxy_0.123456/Rep_001/Sim_t000001.rds
base_dir <- "/Users/4482173/Documents/IMO_workshop13/ABM_results/WGDp_0"
setwd('/Users/4482173/Documents/IMO_workshop13/ABM_results/WGDp_0')
# Output folders for figures
out_scatter_dir <- file.path("analysis_plots_WGDp_0", "scatter_by_time")
out_violin_dir  <- file.path("analysis_plots_WGDp_0", "violin_P_over_time")
out_cells_dir <- file.path("analysis_plots_WGDp_0", "scatter_cells_by_time")
out_ploidy_hist_dir <- file.path("analysis_plots_WGDp_0", "ploidy_hist_by_time")
dir_create <- function(p) if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
dir_create(out_scatter_dir)
dir_create(out_violin_dir)
dir_create(out_cells_dir)
dir_create(out_ploidy_hist_dir)

# ---- Preload toggle ----
USE_PRELOAD <- TRUE  # set FALSE to use on-demand disk IO

#
# ---- Helpers ----

# Extract integer time t from a filename like Sim_t000123.rds
parse_t_from_filename <- function(f) {
  mt <- stringr::str_match(basename(f), "Sim_t(\\d+)\\.rds$")
  if (!is.na(mt[1,2])) as.integer(mt[1,2]) else NA_integer_
}

# Keep only Sim_t000001 ... Sim_t000720 files
filter_rds_by_t_range <- function(rds_files, tmin = 1L, tmax = 720L) {
  tt <- vapply(rds_files, parse_t_from_filename, integer(1))
  rds_files[which(is.finite(tt) & tt >= tmin & tt <= tmax)]
}

# Safely read state from RDS and return a standardized list:
#   list(cells = data.frame, coxy = numeric|NA, rep = integer|NA, t = integer)
read_state_file <- function(f, coxy = NA_real_, rep_id = NA_integer_) {
  st <- readRDS(f)
  # accept state$cells or state$Cells
  cells <- NULL
  if (is.list(st) && !is.null(st$cells)) cells <- st$cells
  if (is.null(cells) && !is.null(st$Cells)) cells <- st$Cells
  if (is.null(cells)) stop("Neither state$cells nor state$Cells found in: ", f)
  
  # accept state$O2 (oxygen field) if present
  O2 <- NULL
  if (is.list(st) && !is.null(st$O2)) O2 <- st$O2
  
  # ensure required columns exist
  must_have <- c("X","Y","P","Status")
  missing <- setdiff(must_have, names(cells))
  if (length(missing) > 0) {
    stop("Missing columns in cells: ", paste(missing, collapse = ", "), " in file: ", f)
  }
  
  # parse time t from filename: Sim_t000123.rds -> t = 123 (integer hours)
  fn <- basename(f)
  t <- NA_integer_
  mt <- str_match(fn, "Sim_t(\\d+)\\.rds$")
  if (!is.na(mt[1,2])) t <- as.integer(mt[1,2])
  
  list(cells = cells, coxy = coxy, rep = rep_id, t = t, o2 = O2)
}

# List all run leaf directories that contain Sim_t*.rds
# Pattern: base_dir/Coxy_xxxxxx/Rep_###/*.rds
list_run_dirs <- function(base_dir) {
  # depth=3 to reach ".../Rep_###"
  cand <- dir_ls(base_dir, type = "directory", recurse = TRUE)
  # keep dirs that actually contain Sim_t*.rds
  keep <- cand[map_lgl(cand, ~ length(dir_ls(.x, regexp = "Sim_t\\d+\\.rds$")) > 0)]
  # Return data.frame with parsed coxy and rep
  tibble(run_dir = keep) |>
    mutate(
      coxy_str = basename(dirname(run_dir)),     # Coxy_xxxxxx
      rep_str  = basename(run_dir)               # Rep_###
    ) |>
    mutate(
      coxy = suppressWarnings(as.numeric(str_remove(coxy_str, "^Coxy_"))),
      rep  = suppressWarnings(as.integer(str_remove(rep_str,  "^Rep_")))
    )
}

# ---- Preload all runs into memory (optional) ----
coxy_key <- function(x) sprintf("Coxy_%s", format(x, scientific = FALSE))
rep_key  <- function(r) sprintf("Rep_%03d", r)

# Load all Sim_t*.rds into a nested list: data[[Coxy_key]][[Rep_key]] = list of records
# Each record is a list(cells=..., coxy=..., rep=..., t=...)
load_all_runs <- function(base_dir, verbose = TRUE) {
  runs <- list_run_dirs(base_dir)
  if (nrow(runs) == 0) stop("No run directories found under: ", base_dir)
  
  data_env <- new.env(parent = emptyenv())
  index_list <- vector("list", nrow(runs))
  
  for (i in seq_len(nrow(runs))) {
    run_dir <- runs$run_dir[i]
    cx      <- runs$coxy[i]
    rp      <- runs$rep[i]
    cxk     <- coxy_key(cx)
    rpk     <- rep_key(rp)
    
    rds_files <- dir_ls(run_dir, regexp = "Sim_t\\d+\\.rds$", type = "file") |> sort()
    if (length(rds_files) == 0) next
    
    # Prepare nested list containers in environment
    if (!exists(cxk, envir = data_env, inherits = FALSE)) assign(cxk, list(), envir = data_env)
    cx_list <- get(cxk, envir = data_env, inherits = FALSE)
    if (is.null(cx_list[[rpk]])) cx_list[[rpk]] <- list()
    
    recs <- vector("list", length(rds_files))
    for (k in seq_along(rds_files)) {
      rec <- read_state_file(rds_files[k], coxy = cx, rep_id = rp)
      recs[[k]] <- rec
    }
    cx_list[[rpk]] <- recs
    assign(cxk, cx_list, envir = data_env)
    
    index_list[[i]] <- tibble(
      run_dir = run_dir,
      coxy    = cx,
      rep     = rp,
      n_files = length(rds_files)
    )
    if (verbose) message(sprintf("Loaded %d files into %s/%s", length(rds_files), cxk, rpk))
  }
  
  # Build a simple list view from the environment
  coxy_names <- ls(envir = data_env)
  data_list <- lapply(coxy_names, function(nm) get(nm, envir = data_env, inherits = FALSE))
  names(data_list) <- coxy_names
  
  list(
    runs  = list_run_dirs(base_dir),
    data  = data_list,
    index = bind_rows(index_list)
  )
}

# ---- O2 scatter helpers (fixed [0,1] scale, palette distinct from P) ----
make_o2_palette <- function(n = 256L) {
  # Distinct from viridis "plasma" used for P: deep blue -> cyan -> light -> yellow -> orange -> red
  grDevices::colorRampPalette(c("#002147", "#0066CC", "#33CCFF", "#CCFFFF", "#FFE680", "#FF9900", "#CC3300"))(n)
}

# Given an O2 matrix, return a tidy data.frame (x,y,val in [0,1])
o2_to_df <- function(O2) {
  if (is.null(O2)) return(NULL)
  N <- nrow(O2); M <- ncol(O2)
  # clamp to [0,1]
  Z <- pmin(pmax(O2, 0), 1)
  expand.grid(X = seq_len(M), Y = seq_len(N)) |>
    dplyr::mutate(O2 = as.numeric(t(Z)))  # row-major to x-right, y-up
}

# ---- GIF helpers ----
list_pngs_sorted <- function(dir) {
  files <- fs::dir_ls(dir, type = "file", glob = "*.png") |> as.character()
  sort(files)
}

# Extract integer t from PNG filename suffix like *_t000123.png
extract_t_from_png <- function(f) {
  mt <- stringr::str_match(basename(f), "_t(\\d+)\\.png$")
  if (!is.na(mt[1,2])) as.integer(mt[1,2]) else NA_integer_
}

# Keep t = 1 and multiples of 'stride' (e.g., 24)
select_sampled_pngs <- function(pngs, stride = 24L) {
  if (length(pngs) == 0) return(pngs)
  tt <- vapply(pngs, extract_t_from_png, integer(1))
  keep_idx <- which(is.finite(tt) & (tt == 1L | (tt %% stride) == 0L))
  pngs[keep_idx]
}

make_gif_in_dir <- function(in_dir, gif_name = NULL, fps = 12, loop = 0L, stride = 24L) {
  if (!dir.exists(in_dir)) return(invisible(NULL))
  pngs <- list_pngs_sorted(in_dir)
  pngs <- select_sampled_pngs(pngs, stride = stride)
  if (length(pngs) == 0) return(invisible(NULL))
  if (is.null(gif_name)) {
    gif_name <- paste0(basename(in_dir), ".gif")
  }
  gif_path <- file.path(in_dir, gif_name)
  img <- magick::image_read(pngs)
  anim <- magick::image_animate(img, fps = fps, loop = loop)
  magick::image_write(anim, path = gif_path)
  message(sprintf("GIF written: %s  (frames=%d, fps=%d)", gif_path, length(pngs), fps))
  invisible(gif_path)
}

make_gifs_for_parent <- function(parent_dir, fps = 12, loop = 0L, stride = 24L) {
  if (!dir.exists(parent_dir)) return(invisible(NULL))
  subdirs <- fs::dir_ls(parent_dir, type = "directory", recurse = FALSE)
  for (d in subdirs) {
    make_gif_in_dir(d, fps = fps, loop = loop, stride = stride)
  }
  invisible(NULL)
}


# ---- helper: robust cells accessor ----
# Get cells table from either $cells or $Cells; normalize x/y/p/status to X/Y/P/Status
get_cells_tbl <- function(rec_or_state) {
  df <- NULL
  if (!is.null(rec_or_state$cells)) df <- rec_or_state$cells
  else if (!is.null(rec_or_state$Cells)) df <- rec_or_state$Cells
  else return(NULL)
  cn <- names(df)
  if ("x" %in% cn && !"X" %in% cn) df$X <- df$x
  if ("y" %in% cn && !"Y" %in% cn) df$Y <- df$y
  if ("p" %in% cn && !"P" %in% cn) df$P <- df$p
  if ("status" %in% cn && !"Status" %in% cn) df$Status <- df$status
  df
}

# ---- O2 field scatter: white background, no axes ----
plot_O2_scatter_for_run <- function(run_dir, coxy, rep, preloaded = NULL) {
  subdir <- file.path(out_scatter_dir,
                      sprintf("Coxy_%s__Rep_%03d_O2", format(coxy, scientific = FALSE), rep))
  dir_create(subdir)
  
  pal <- make_o2_palette()
  
  plot_one <- function(rec_or_state, tt, o2_mat) {
    if (is.null(o2_mat)) return()
    o2df <- o2_to_df(o2_mat); if (is.null(o2df)) return()
    p <- ggplot(o2df, aes(x = .data$X, y = .data$Y, color = .data$O2)) +
      geom_point(shape = 15, size = 1.0, alpha = 0.95) +
      coord_fixed() +
      scale_color_gradientn(colors = pal, limits = c(0, 1), oob = scales::squish, name = "O2") +
      theme_void(base_size = 12) +
      theme(plot.background = element_rect(fill = "white", colour = NA),
            panel.background = element_rect(fill = "white", colour = NA))
    png_path <- file.path(subdir, sprintf("O2_scatter_t%06d.png", tt))
    ggsave(png_path, p, width = 6.8, height = 6.2, dpi = 150, bg = "white")
  }
  
  if (!is.null(preloaded)) {
    cxk <- coxy_key(coxy); rpk <- rep_key(rep)
    if (!is.null(preloaded$data[[cxk]]) && !is.null(preloaded$data[[cxk]][[rpk]])) {
      for (rec in preloaded$data[[cxk]][[rpk]]) {
        if (!is.finite(rec$t) || rec$t < 1L || rec$t > 720L) next
        plot_one(rec, rec$t, rec$o2)
      }
      return(invisible(NULL))
    }
  }
  
  rds_files <- dir_ls(run_dir, regexp = "Sim_t\\d+\\.rds$", type = "file") |> sort()
  rds_files <- filter_rds_by_t_range(rds_files, 1L, 720L)
  if (length(rds_files) == 0) return(invisible(NULL))
  for (f in rds_files) {
    st <- read_state_file(f, coxy = coxy, rep_id = rep)
    plot_one(st, st$t, st$o2)
  }
  invisible(NULL)
}

# ---- P-colored scatter (plasma), limits [0,10], white background, no axes ----
plot_scatter_for_run <- function(run_dir, coxy, rep, preloaded = NULL) {
  subdir <- file.path(out_scatter_dir,
                      sprintf("Coxy_%s__Rep_%03d", format(coxy, scientific = FALSE), rep))
  dir_create(subdir)
  
  plot_one <- function(tt, cells_df) {
    cells_df <- dplyr::filter(cells_df, .data$Status == 1L)
    if (nrow(cells_df) == 0) return()
    p <- ggplot(cells_df, aes(x = .data$X, y = .data$Y, color = .data$P)) +
      geom_point(size = 0.8, alpha = 0.8) +
      coord_fixed() +
      scale_color_viridis(option = "plasma", limits = c(0, 10), oob = scales::squish, name = "P") +
      theme_void(base_size = 12) +
      theme(plot.background = element_rect(fill = "white", colour = NA),
            panel.background = element_rect(fill = "white", colour = NA))
    png_path <- file.path(subdir, sprintf("ploidy_scatter_t%06d.png", tt))
    ggsave(png_path, p, width = 6.8, height = 6.2, dpi = 150, bg = "white")
  }
  
  if (!is.null(preloaded)) {
    cxk <- coxy_key(coxy); rpk <- rep_key(rep)
    if (!is.null(preloaded$data[[cxk]]) && !is.null(preloaded$data[[cxk]][[rpk]])) {
      for (rec in preloaded$data[[cxk]][[rpk]]) {
        if (!is.finite(rec$t) || rec$t < 1L || rec$t > 720L) next
        cells <- get_cells_tbl(rec); if (is.null(cells)) next
        plot_one(rec$t, cells)
      }
      return(invisible(NULL))
    }
  }
  
  rds_files <- dir_ls(run_dir, regexp = "Sim_t\\d+\\.rds$", type = "file") |> sort()
  rds_files <- filter_rds_by_t_range(rds_files, 1L, 720L)
  if (length(rds_files) == 0) return(invisible(NULL))
  for (f in rds_files) {
    st <- read_state_file(f, coxy = coxy, rep_id = rep)
    cells <- get_cells_tbl(st); if (is.null(cells)) next
    plot_one(st$t, cells)
  }
  invisible(NULL)
}

# ---- Cells-only scatter (orange), white background, no axes ----
plot_cells_scatter_for_run <- function(run_dir, coxy, rep, preloaded = NULL) {
  subdir <- file.path(out_cells_dir,
                      sprintf("Coxy_%s__Rep_%03d_cells", format(coxy, scientific = FALSE), rep))
  dir_create(subdir)
  
  plot_one <- function(tt, cells_df) {
    cells_df <- dplyr::filter(cells_df, .data$Status == 1L)
    if (nrow(cells_df) == 0) return()
    p <- ggplot(cells_df, aes(x = .data$X, y = .data$Y, color = "Cells")) +
      geom_point(shape = 15, size = 0.9, alpha = 0.9, color = "orange") +
      scale_color_manual(values = c("Cells" = "orange"), name = NULL) +
      coord_fixed() +
      theme_void(base_size = 12) +
      theme(plot.background = element_rect(fill = "white", colour = NA),
            panel.background = element_rect(fill = "white", colour = NA))
    png_path <- file.path(subdir, sprintf("cells_scatter_t%06d.png", tt))
    ggsave(png_path, p, width = 6.8, height = 6.2, dpi = 150, bg = "white")
  }
  
  if (!is.null(preloaded)) {
    cxk <- coxy_key(coxy); rpk <- rep_key(rep)
    if (!is.null(preloaded$data[[cxk]]) && !is.null(preloaded$data[[cxk]][[rpk]])) {
      for (rec in preloaded$data[[cxk]][[rpk]]) {
        if (!is.finite(rec$t) || rec$t < 1L || rec$t > 720L) next
        cells <- get_cells_tbl(rec); if (is.null(cells)) next
        plot_one(rec$t, cells)
      }
      return(invisible(NULL))
    }
  }
  
  rds_files <- dir_ls(run_dir, regexp = "Sim_t\\d+\\.rds$", type = "file") |> sort()
  rds_files <- filter_rds_by_t_range(rds_files, 1L, 720L)
  if (length(rds_files) == 0) return(invisible(NULL))
  for (f in rds_files) {
    st <- read_state_file(f, coxy = coxy, rep_id = rep)
    cells <- get_cells_tbl(st); if (is.null(cells)) next
    plot_one(st$t, cells)
  }
  invisible(NULL)
}

# ---- Ploidy distribution per timepoint (x:[0,10], y proportion [0,1]) ----
plot_ploidy_hist_for_run <- function(run_dir, coxy, rep, preloaded = NULL, binwidth = 0.2) {
  subdir <- file.path(out_ploidy_hist_dir,
                      sprintf("Coxy_%s__Rep_%03d_hist", format(coxy, scientific = FALSE), rep))
  dir_create(subdir)

  # A small local plotting helper
  .plot_one <- function(df, tt) {
    if (nrow(df) == 0) return()
    p <- ggplot(df, aes(x = P)) +
      geom_histogram(binwidth = binwidth, boundary = 0, closed = "left",
                     aes(y = after_stat(count / sum(count)),
                         fill = after_stat(x))) +
      scale_x_continuous(limits = c(0, 10), expand = expansion(mult = c(0, 0))) +
      scale_y_continuous(limits = c(0, 0.5), expand = expansion(mult = c(0, 0))) +
      scale_fill_viridis(option = "plasma", limits = c(0, 10), oob = scales::squish, name = "P") +
      labs(x = "Ploidy (P)", y = "Proportion") +
      theme_minimal(base_size = 12) +
      theme(plot.background = element_rect(fill = "white", colour = NA),
            panel.background = element_rect(fill = "white", colour = NA),
            panel.grid = element_blank())
    png_path <- file.path(subdir, sprintf("ploidy_hist_t%06d.png", tt))
    ggsave(png_path, p, width = 6, height = 6, dpi = 150, bg = "white")
  }

  # Preloaded
  if (!is.null(preloaded)) {
    cxk <- coxy_key(coxy); rpk <- rep_key(rep)
    if (!is.null(preloaded$data[[cxk]]) && !is.null(preloaded$data[[cxk]][[rpk]])) {
      recs <- preloaded$data[[cxk]][[rpk]]
      for (rec in recs) {
        if (!is.finite(rec$t) || rec$t < 1L || rec$t > 720L) next
        cells <- rec$cells |> dplyr::filter(.data$Status == 1L, is.finite(.data$P))
        .plot_one(tibble(P = as.numeric(cells$P)), rec$t)
      }
      return(invisible(NULL))
    }
  }

  # Disk
  rds_files <- dir_ls(run_dir, regexp = "Sim_t\\d+\\.rds$", type = "file") |> sort()
  rds_files <- filter_rds_by_t_range(rds_files, 1L, 720L)
  if (length(rds_files) == 0) return(invisible(NULL))
  for (f in rds_files) {
    st <- read_state_file(f, coxy = coxy, rep_id = rep)
    cells <- st$cells |> dplyr::filter(.data$Status == 1L, is.finite(.data$P))
    .plot_one(tibble(P = as.numeric(cells$P)), st$t)
  }
  invisible(NULL)
}

# ---- ② Violin: all Coxy & all reps over time ----
# Build long table: one row per cell per timepoint across all runs
collect_long_P <- function(base_dir, preloaded = NULL) {
  # Path A: build from preloaded object if available
  if (!is.null(preloaded)) {
    out <- vector("list", length(preloaded$data))
    ii <- 0L
    for (cx_name in names(preloaded$data)) {
      cx_val <- suppressWarnings(as.numeric(sub("^Coxy_", "", cx_name)))
      reps_list <- preloaded$data[[cx_name]]
      for (rp_name in names(reps_list)) {
        rp_val <- suppressWarnings(as.integer(sub("^Rep_", "", rp_name)))
        recs <- reps_list[[rp_name]]
        if (length(recs) == 0) next
        df_list <- vector("list", length(recs))
        for (k in seq_along(recs)) {
          rec <- recs[[k]]
          cells <- rec$cells |> dplyr::filter(.data$Status == 1L)
          if (!is.null(cells) && nrow(cells) > 0) {
            df_list[[k]] <- tibble(
              coxy = cx_val,
              rep  = rp_val,
              t    = rec$t,
              P    = as.numeric(cells$P)
            )
          }
        }
        ii <- ii + 1L
        out[[ii]] <- bind_rows(df_list)
      }
    }
    return(
      bind_rows(out) |>
        dplyr::filter(is.finite(.data$P), .data$P >= 0, !is.na(.data$t), .data$t %% 24 == 0) |>
        dplyr::mutate(
          coxy_f = factor(sprintf("%.6f", .data$coxy),
                          levels = sort(unique(sprintf("%.6f", .data$coxy))))
        )
    )
  }
  
  # Path B (fallback): original disk-based approach
  runs <- list_run_dirs(base_dir)
  if (nrow(runs) == 0) {
    stop("No run directories found under: ", base_dir)
  }
  out_list <- vector("list", nrow(runs))
  for (i in seq_len(nrow(runs))) {
    run_dir <- runs$run_dir[i]
    coxy    <- runs$coxy[i]
    rep     <- runs$rep[i]
    rds_files <- dir_ls(run_dir, regexp = "Sim_t\\d+\\.rds$", type = "file") |> sort()
    if (length(rds_files) == 0) next
    
    df_list <- vector("list", length(rds_files))
    for (k in seq_along(rds_files)) {
      rec <- read_state_file(rds_files[k], coxy = coxy, rep_id = rep)
      cells <- rec$cells |> dplyr::filter(.data$Status == 1L)
      if (!is.null(cells) && nrow(cells) > 0) {
        df_list[[k]] <- tibble(
          coxy = rec$coxy,
          rep  = rec$rep,
          t    = rec$t,
          P    = as.numeric(cells$P)
        )
      }
    }
    out_list[[i]] <- bind_rows(df_list)
  }
  bind_rows(out_list) |>
    dplyr::filter(is.finite(.data$P), .data$P >= 0, !is.na(.data$t), .data$t %% 24 == 0) |>
    dplyr::mutate(
      coxy_f = factor(sprintf("%.6f", .data$coxy), levels = sort(unique(sprintf("%.6f", .data$coxy))))
    )
}

# Draw violin plot:
plot_violin_P_over_time <- function(long_df) {
  # We want x = time, y = P, color = Coxy; use violin + overlay median
  # If time has many unique values, consider thinning/breaks. Here we plot all.
  p <- ggplot(long_df, aes(x = factor(t), y = P, fill = coxy_f, color = coxy_f)) +
    geom_violin(trim = TRUE, alpha = 0.35, linewidth = 0.25) +
    stat_summary(fun = median, geom = "point", size = 0.8, position = position_dodge(width = 0.9)) +
    #scale_y_continuous(limits = c(0, 10), expand = expansion(mult = c(0.02, 0.02))) +
    scale_fill_viridis(discrete = TRUE, option = "plasma") +
    scale_color_viridis(discrete = TRUE, option = "plasma") +
    labs(
      title = "Distribution of P over time by Coxy",
      x = "Time (t, simulation step)",
      y = "P",
      fill = "Coxy",
      color = "Coxy"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
    )
  p
}

# ===== Run analysis =====

# 1 Optional preload of all runs
preloaded <- NULL
if (isTRUE(USE_PRELOAD)) {
  message("Preloading all runs into memory...")
  preloaded <- load_all_runs(base_dir)
}

# Derive runs list either from preloaded or from disk
runs <- if (!is.null(preloaded)) preloaded$runs else list_run_dirs(base_dir)
if (nrow(runs) == 0) stop("No runs found under base_dir: ", base_dir)

message("Generating ploidy-colored scatter plots for each time point (t=1..720)...")
for (i in seq_len(nrow(runs))) {
  plot_scatter_for_run(runs$run_dir[i], coxy = runs$coxy[i], rep = runs$rep[i], preloaded = preloaded)
}

message("Generating cells-only scatter plots for each time point (t=1..720)...")
for (i in seq_len(nrow(runs))) {
  plot_cells_scatter_for_run(runs$run_dir[i], coxy = runs$coxy[i], rep = runs$rep[i], preloaded = preloaded)
}

message("Generating O2-field scatter plots for each time point (t=1..720)...")
for (i in seq_len(nrow(runs))) {
  plot_O2_scatter_for_run(runs$run_dir[i], coxy = runs$coxy[i], rep = runs$rep[i], preloaded = preloaded)
}

message("Generating ploidy distribution histograms for each time point (t=1..720)...")
for (i in seq_len(nrow(runs))) {
  plot_ploidy_hist_for_run(runs$run_dir[i], coxy = runs$coxy[i], rep = runs$rep[i], preloaded = preloaded)
}

# 2 Violin across all runs
message("Collecting P across all runs for violin plot...")
long_df <- collect_long_P(base_dir, preloaded = preloaded)

long_df_use <- long_df[which(long_df$t <= 720),]



message("Drawing violin plot...")
p_violin <- plot_violin_P_over_time(long_df_use)
violin_path <- file.path(out_violin_dir, "Violin_P_over_time_by_Coxy.pdf")
ggsave(violin_path, p_violin, width = 20, height = 10)

message("Assembling GIFs from PNG sequences (smooth animation, ~20 fps, stride = 24)...")
make_gifs_for_parent(out_scatter_dir, fps = 2, loop = 0, stride = 24)
make_gifs_for_parent(out_cells_dir, fps = 2, loop = 0, stride = 24)
make_gifs_for_parent(out_ploidy_hist_dir, fps = 2, loop = 0, stride = 24)

message("Done.")
cat(sprintf("\nPloidy-colored scatter PNGs saved under: %s\nCells-only scatter PNGs saved under: %s\nO2 scatter PNGs saved under: %s\nPloidy histogram PNGs saved under: %s\nViolin PDF saved at: %s\n(GIFs created in each corresponding subfolder as *.gif)\n",
            out_scatter_dir,
            out_cells_dir,
            out_scatter_dir,
            out_ploidy_hist_dir,
            violin_path))



create_gif_strip_1x4 <- function(gif_paths, out_path = "strip_1x4.gif",
                                 target_height = 600,   
                                 fps = 2,             
                                 loop = 0) {            # 0=无限循环
  stopifnot(length(gif_paths) == 4)
  
  # 读入并“固化”每一帧（处理 dispose），然后统一高度
  imgs <- lapply(gif_paths, function(p) {
    x <- magick::image_read(p)
    x <- magick::image_coalesce(x)
    if (!is.null(target_height)) {
      x <- magick::image_scale(x, paste0("x", target_height))  # 只给高度，宽度按比例
    }
    x
  })
  
  # 取最短的帧数对齐（也可改成 max 并用最后一帧补齐）
  nframes <- min(length(imgs[[1]]), length(imgs[[2]]), length(imgs[[3]]), length(imgs[[4]]))
  
  # 逐帧横向拼接
  frames <- vector("list", nframes)
  for (i in seq_len(nframes)) {
    row_i <- magick::image_join(imgs[[1]][i], imgs[[2]][i], imgs[[3]][i], imgs[[4]][i])
    frames[[i]] <- magick::image_append(row_i, stack = FALSE)  # FALSE = 水平拼接
  }
  
  # 合成动画并写出
  anim <- magick::image_animate(magick::image_join(frames), fps = fps, loop = loop)
  magick::image_write(anim, path = out_path)
  message(sprintf("1x4 GIF saved: %s (frames=%d, fps=%d)", out_path, nframes, fps))
  invisible(out_path)
}


create_gif_strip_1x4(
  gif_paths = c(
    "/Users/4482173/Documents/IMO_workshop13/Hpoxia/Low/analysis_plots_low/Coxy_1__Rep_001_cells.gif",
    "/Users/4482173/Documents/IMO_workshop13/Hpoxia/Low/analysis_plots_low/Coxy_1__Rep_001.gif",
    "/Users/4482173/Documents/IMO_workshop13/Hpoxia/Low/analysis_plots_low/Coxy_1__Rep_001_O2.gif",
    "/Users/4482173/Documents/IMO_workshop13/Hpoxia/Low/analysis_plots_low/Coxy_1__Rep_001_hist.gif"
  ),
  out_path = "/Users/4482173/Documents/IMO_workshop13/Hpoxia/Low/analysis_plots_low/HP_HO.gif",
  target_height = 600,
  fps = 2,   
  loop = 0
)


create_gif_strip_1x4(
  gif_paths = c(
    "/Users/4482173/Documents/IMO_workshop13/Hpoxia/Low/analysis_plots_low/Coxy_0.5__Rep_001_cells.gif",
    "/Users/4482173/Documents/IMO_workshop13/Hpoxia/Low/analysis_plots_low/Coxy_0.5__Rep_001.gif",
    "/Users/4482173/Documents/IMO_workshop13/Hpoxia/Low/analysis_plots_low/Coxy_0.5__Rep_001_O2.gif",
    "/Users/4482173/Documents/IMO_workshop13/Hpoxia/Low/analysis_plots_low/Coxy_0.5__Rep_001_hist.gif"
  ),
  out_path = "/Users/4482173/Documents/IMO_workshop13/Hpoxia/Low/analysis_plots_low/HP_MO.gif",
  target_height = 600,
  fps = 2,   
  loop = 0
)


create_gif_strip_1x4(
  gif_paths = c(
    "/Users/4482173/Documents/IMO_workshop13/Hpoxia/Low/analysis_plots_low/Coxy_0.05__Rep_001_cells.gif",
    "/Users/4482173/Documents/IMO_workshop13/Hpoxia/Low/analysis_plots_low/Coxy_0.05__Rep_001.gif",
    "/Users/4482173/Documents/IMO_workshop13/Hpoxia/Low/analysis_plots_low/Coxy_0.05__Rep_001_O2.gif",
    "/Users/4482173/Documents/IMO_workshop13/Hpoxia/Low/analysis_plots_low/Coxy_0.05__Rep_001_hist.gif"
  ),
  out_path = "/Users/4482173/Documents/IMO_workshop13/Hpoxia/Low/analysis_plots_low/HP_LO.gif",
  target_height = 600,
  fps = 2,   
  loop = 0
)


create_gif_strip_1x2 <- function(gif_paths, out_path = "strip_1x2.gif",
                                 target_height = 600,   
                                 fps = 2,             
                                 loop = 0) {            # 0=无限循环
  stopifnot(length(gif_paths) == 2)
  
  # 读入并“固化”每一帧（处理 dispose），然后统一高度
  imgs <- lapply(gif_paths, function(p) {
    x <- magick::image_read(p)
    x <- magick::image_coalesce(x)
    if (!is.null(target_height)) {
      x <- magick::image_scale(x, paste0("x", target_height))  # 只给高度，宽度按比例
    }
    x
  })
  
  # 取最短的帧数对齐（也可改成 max 并用最后一帧补齐）
  nframes <- min(length(imgs[[1]]), length(imgs[[2]]))
  
  # 逐帧横向拼接
  frames <- vector("list", nframes)
  for (i in seq_len(nframes)) {
    row_i <- magick::image_join(imgs[[1]][i], imgs[[2]][i])
    frames[[i]] <- magick::image_append(row_i, stack = FALSE)  # FALSE = 水平拼接
  }
  
  # 合成动画并写出
  anim <- magick::image_animate(magick::image_join(frames), fps = fps, loop = loop)
  magick::image_write(anim, path = out_path)
  message(sprintf("1x2 GIF saved: %s (frames=%d, fps=%d)", out_path, nframes, fps))
  invisible(out_path)
}


create_gif_strip_1x2(
  gif_paths = c(
    "/Users/4482173/Documents/IMO_workshop13/Hpoxia/Low/analysis_plots_low/Coxy_1__Rep_001_hist.gif",
    "/Users/4482173/Documents/IMO_workshop13/Hpoxia/High/analysis_plots_high/Coxy_1__Rep_001_hist.gif"
  ),
  out_path = "/Users/4482173/Documents/IMO_workshop13/Hpoxia/Ploidy_hist.gif",
  target_height = 600,
  fps = 2,   
  loop = 0
)


create_gif_strip_2x3 <- function(gif_paths,
                                 out_path = "strip_2x3.gif",
                                 target_height = 600,   # tile height (each of the 6 GIFs after scaling)
                                 fps = 2,               # desired frames per second
                                 loop = 0) {            # 0 = infinite loop
  # ---- Basic checks ----
  stopifnot(length(gif_paths) == 6)
  if (!requireNamespace("magick", quietly = TRUE)) {
    stop("Package 'magick' is required.")
  }
  
  # ---- Read, coalesce (fix disposal), and scale each GIF to the same tile height ----
  imgs <- lapply(gif_paths, function(p) {
    x <- magick::image_read(p)
    x <- magick::image_coalesce(x)
    if (!is.null(target_height)) {
      # Scale by height only; width keeps aspect ratio
      x <- magick::image_scale(x, paste0("x", target_height))
    }
    x
  })
  
  # ---- Make every tile the same width by letterboxing to the max width (transparent padding) ----
  # Using the first frame of each to compute widths (after height scaling)
  tile_width <- max(vapply(imgs, function(x) magick::image_info(x[1])$width, numeric(1)))
  # Extend all frames of each GIF to (tile_width x target_height), centered, transparent background
  imgs <- lapply(imgs, function(x) {
    magick::image_extent(
      x,
      geometry = sprintf("%dx%d", tile_width, target_height),
      gravity = "center",
      color = "none"
    )
  })
  
  # ---- Synchronize by the shortest length ----
  nframes <- min(vapply(imgs, length, integer(1)))
  
  # ---- Build frames: (1,2,3) -> row1 (horizontal); (4,5,6) -> row2; then stack row1 over row2 ----
  frames <- vector("list", nframes)
  for (i in seq_len(nframes)) {
    row1 <- magick::image_append(
      magick::image_join(imgs[[1]][i], imgs[[2]][i], imgs[[3]][i]),
      stack = FALSE # horizontal
    )
    row2 <- magick::image_append(
      magick::image_join(imgs[[4]][i], imgs[[5]][i], imgs[[6]][i]),
      stack = FALSE # horizontal
    )
    # Two rows have identical width thanks to tile_width; safe to append vertically
    frames[[i]] <- magick::image_append(
      magick::image_join(row1, row2),
      stack = TRUE  # vertical
    )
  }
  
  # ---- Animate ----
  # Use 'delay' to avoid the GIF constraint that fps must divide 100 evenly
  delay <- 100 / fps
  anim <- magick::image_animate(magick::image_join(frames), delay = delay, loop = loop)
  
  # ---- Write out ----
  magick::image_write(anim, path = out_path)
  message(sprintf("2x3 GIF saved: %s (frames=%d, fps=%.3f, delay=%.3f/100s)",
                  out_path, nframes, fps, delay))
  invisible(out_path)
}


create_gif_strip_2x3(
  gif_paths = c(
    "/Users/4482173/Documents/IMO_workshop13/Hpoxia/Low/analysis_plots_low/Coxy_1__Rep_001_cells.gif",
    "/Users/4482173/Documents/IMO_workshop13/Hpoxia/Low/analysis_plots_low/Coxy_1__Rep_001.gif",
    "/Users/4482173/Documents/IMO_workshop13/Hpoxia/Low/analysis_plots_low/Coxy_1__Rep_001_O2.gif",
    "/Users/4482173/Documents/IMO_workshop13/Drug/2N/analysis_plots_high/Coxy_1__Rep_001_cells.gif",
    "/Users/4482173/Documents/IMO_workshop13/Drug/2N/analysis_plots_high/Coxy_1__Rep_001.gif",
    "/Users/4482173/Documents/IMO_workshop13/Drug/2N/analysis_plots_high/Coxy_1__Rep_001_O2.gif"
  ),
  out_path = "/Users/4482173/Documents/IMO_workshop13/ABM.gif",
  target_height = 600,
  fps = 2,   
  loop = 0
)





