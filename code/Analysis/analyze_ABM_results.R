#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(viridis)
  library(dplyr)
  library(purrr)
  library(scales)
  library(parallel)
})

# --------- Helper: discover parameter combinations ---------

find_param_combos <- function(base_dir) {
  wgdp_dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)
  wgdp_dirs <- wgdp_dirs[grepl("^WGDp_", basename(wgdp_dirs))]
  
  combos <- list()
  
  for (wgdp in wgdp_dirs) {
    wgdr_dirs <- list.dirs(wgdp, full.names = TRUE, recursive = FALSE)
    wgdr_dirs <- wgdr_dirs[grepl("^WGDr_", basename(wgdr_dirs))]
    if (length(wgdr_dirs) == 0L) next
    
    for (wgdr in wgdr_dirs) {
      coxy_dirs <- list.dirs(wgdr, full.names = TRUE, recursive = FALSE)
      coxy_dirs <- coxy_dirs[grepl("^Coxy_", basename(coxy_dirs))]
      if (length(coxy_dirs) == 0L) next
      
      for (coxy in coxy_dirs) {
        rep_dirs <- list.dirs(coxy, full.names = TRUE, recursive = FALSE)
        rep_dirs <- rep_dirs[grepl("^Rep_", basename(rep_dirs))]
        if (length(rep_dirs) == 0L) next
        
        for (rep_dir in rep_dirs) {
          combos[[length(combos) + 1L]] <- list(
            wgdp = basename(wgdp),
            wgdr = basename(wgdr),
            coxy = basename(coxy),
            rep  = basename(rep_dir),
            path = rep_dir
          )
        }
      }
    }
  }
  
  combos
}

make_param_label <- function(combo) {
  sprintf("%s__%s__%s__%s",
          combo$wgdp, combo$wgdr, combo$coxy, combo$rep)
}

# --------- Helper: load states for one parameter combination ---------

load_states_for_combo <- function(combo) {
  rds_files <- list.files(combo$path, pattern = "\\.rds$", full.names = TRUE)
  if (length(rds_files) == 0L) {
    warning("No .rds files found in ", combo$path)
    return(list())
  }
  rds_files <- sort(rds_files)
  states <- lapply(rds_files, readRDS)
  states
}

# --------- Helper: directory creation ---------

ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  invisible(path)
}

# --------- Plot helpers ---------

plot_cells <- function(state, out_path) {
  df <- state[["cells"]]
  if (is.null(df) || nrow(df) == 0L) return(invisible(NULL))
  
  p <- ggplot(df, aes(x = X, y = Y)) +
    geom_point(color = "#FFA500", size = 0.6) +
    coord_fixed() +
    theme_void() +
    theme(
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA)
    )
  
  ggsave(out_path, p, width = 4, height = 4, dpi = 150)
}

plot_ploidy_cells <- function(state, out_path) {
  df <- state[["cells"]]
  if (is.null(df) || nrow(df) == 0L) return(invisible(NULL))
  if (!("P" %in% names(df))) return(invisible(NULL))
  
  p <- ggplot(df, aes(x = X, y = Y, colour = pmin(pmax(P, 0), 10))) +
    geom_point(size = 0.6) +
    coord_fixed() +
    scale_color_viridis(option = "plasma", limits = c(0, 10), oob = squish) +
    theme_void() +
    theme(
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA),
      legend.position  = "none"
    )
  
  ggsave(out_path, p, width = 4, height = 4, dpi = 150)
}

plot_O2 <- function(state, out_path) {
  O2 <- state[["O2"]]
  if (is.null(O2)) return(invisible(NULL))
  nr <- nrow(O2)
  nc <- ncol(O2)
  if (nr == 0L || nc == 0L) return(invisible(NULL))
  
  df_O2 <- expand.grid(X = seq_len(nr), Y = seq_len(nc))
  df_O2$O2 <- as.vector(O2)
  
  p <- ggplot(df_O2, aes(x = X, y = Y, fill = pmin(pmax(O2, 0), 1))) +
    geom_raster() +
    coord_fixed() +
    scale_fill_viridis(option = "turbo", limits = c(0, 1), oob = squish) +
    theme_void() +
    theme(
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA),
      legend.position  = "none"
    )
  
  ggsave(out_path, p, width = 4, height = 4, dpi = 150)
}

plot_ploidy_hist <- function(state, out_path) {
  df <- state[["cells"]]
  if (is.null(df) || nrow(df) == 0L) return(invisible(NULL))
  if (!("P" %in% names(df))) return(invisible(NULL))
  
  df <- df %>%
    mutate(P_clip = pmin(pmax(P, 0), 10))
  
  p <- ggplot(df, aes(x = P_clip)) +
    geom_histogram(
      aes(y = after_stat(..count.. / sum(..count..)),
          fill = after_stat(x)),
      binwidth = 0.1,
      boundary = 0,
      closed = "right"
    ) +
    scale_fill_viridis(option = "plasma", limits = c(0, 10), oob = squish) +
    xlim(0, 10) +
    scale_y_continuous(limits = c(0, 0.25), expand = expansion(mult = c(0, 0.02))) +
    labs(x = "Ploidy", y = "Proportion per bin") +
    theme_minimal(base_size = 10) +
    theme(
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA),
      legend.position  = "none"
    )
  
  ggsave(out_path, p, width = 4, height = 4, dpi = 150)
}

# --------- Aggregation for trends ---------

compute_dead_ploidy_trend <- function(states) {
  dead_list <- lapply(seq_along(states), function(idx) {
    st <- states[[idx]]
    dl <- st[["dead_log"]]
    if (is.null(dl) || nrow(dl) == 0L) return(NULL)
    step_val <- if (!is.null(st[["step"]])) st[["step"]] else idx
    if (!("P" %in% names(dl))) return(NULL)
    data.frame(step = step_val, P = dl$P)
  })
  dead_df <- bind_rows(dead_list)
  if (nrow(dead_df) == 0L) return(NULL)
  
  dead_df <- dead_df %>%
    mutate(
      P_clip = pmin(pmax(P, 0), 9.9999),
      P_bin  = floor(P_clip),
      P_mid  = P_bin + 0.5
    )
  
  by_step <- dead_df %>%
    group_by(step) %>%
    summarise(total = n(), .groups = "drop")
  
  agg <- dead_df %>%
    group_by(step, P_mid) %>%
    summarise(n_bin = n(), .groups = "drop") %>%
    left_join(by_step, by = "step") %>%
    mutate(prop = n_bin / total) %>%
    filter(is.finite(prop))
  
  agg
}

plot_dead_ploidy_trend <- function(trend_df, out_path) {
  if (is.null(trend_df) || nrow(trend_df) == 0L) return(invisible(NULL))
  
  p <- ggplot(trend_df, aes(x = step, y = prop, colour = P_mid, group = P_mid)) +
    stat_smooth(se = TRUE, method = "loess", span = 0.4) +
    scale_color_viridis(option = "plasma", limits = c(0, 10), oob = squish) +
    ylab("Death proportion") +
    xlab("Step") +
    theme_minimal(base_size = 10) +
    theme(
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA)
    )
  
  ggsave(out_path, p, width = 6, height = 4, dpi = 150)
}

compute_G_ploidy_trend <- function(states) {
  cell_list <- lapply(seq_along(states), function(idx) {
    st <- states[[idx]]
    df <- st[["cells"]]
    if (is.null(df) || nrow(df) == 0L) return(NULL)
    step_val <- if (!is.null(st[["step"]])) st[["step"]] else idx
    if (!("P" %in% names(df)) || !("G" %in% names(df))) return(NULL)
    if ("Status" %in% names(df)) {
      df <- df[df$Status == 1L, , drop = FALSE]
    }
    if (nrow(df) == 0L) return(NULL)
    df$step <- step_val
    df
  })
  
  cells_all <- bind_rows(cell_list)
  if (nrow(cells_all) == 0L) return(NULL)
  
  cells_all <- cells_all %>%
    mutate(
      P_clip = pmin(pmax(P, 0), 9.9999),
      P_bin  = floor(P_clip),
      P_mid  = P_bin + 0.5
    )
  
  agg <- cells_all %>%
    group_by(step, P_mid) %>%
    summarise(
      G_mean = mean(G, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(is.finite(G_mean))
  
  agg
}

plot_G_ploidy_trend <- function(trend_df, out_path) {
  if (is.null(trend_df) || nrow(trend_df) == 0L) return(invisible(NULL))
  
  p <- ggplot(trend_df, aes(x = step, y = G_mean, colour = P_mid, group = P_mid)) +
    stat_smooth(se = TRUE, method = "loess", span = 0.4) +
    scale_color_viridis(option = "plasma", limits = c(0, 10), oob = squish) +
    ylab("G (growth rate)") +
    xlab("Step") +
    theme_minimal(base_size = 10) +
    theme(
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA)
    )
  
  ggsave(out_path, p, width = 6, height = 4, dpi = 150)
}

# --------- Main driver for one parameter combo ---------

analyze_one_combo <- function(combo,
                              base_dir = "/Users/4482173/Documents/IMO_workshop13/ABM_results") {
  analysis_dir <- file.path(base_dir, "Analysis")
  ensure_dir(analysis_dir)
  
  states <- load_states_for_combo(combo)
  if (length(states) == 0L) return(invisible(NULL))
  
  label <- make_param_label(combo)
  
  # Save aggregated list as RDS
  rds_out <- file.path(analysis_dir, paste0("states_", label, ".rds"))
  saveRDS(states, rds_out)
  
  # Create per-combo subfolders
  combo_dir <- file.path(analysis_dir, label)
  ensure_dir(combo_dir)
  cell_dir   <- ensure_dir(file.path(combo_dir, "cell_plot"))
  ploidy_dir <- ensure_dir(file.path(combo_dir, "ploidy_plot"))
  O2_dir     <- ensure_dir(file.path(combo_dir, "O2_plot"))
  pl_stat_dir<- ensure_dir(file.path(combo_dir, "ploidy_stat"))
  
  # Per-step plots
  for (i in seq_along(states)) {
    st <- states[[i]]
    step_val <- if (!is.null(st[["step"]])) st[["step"]] else i
    
    cell_path   <- file.path(cell_dir,   sprintf("cell_step_%04d.png",   step_val))
    ploidy_path <- file.path(ploidy_dir, sprintf("ploidy_step_%04d.png", step_val))
    O2_path     <- file.path(O2_dir,     sprintf("O2_step_%04d.png",     step_val))
    hist_path   <- file.path(pl_stat_dir,sprintf("ploidy_hist_step_%04d.png", step_val))
    
    plot_cells(st, cell_path)
    plot_ploidy_cells(st, ploidy_path)
    plot_O2(st, O2_path)
    plot_ploidy_hist(st, hist_path)
  }
  
  # Dead ploidy trend (saved in Analysis root, parameter-labeled file name)
  dead_trend <- compute_dead_ploidy_trend(states)
  dead_plot_path <- file.path(analysis_dir,
                              paste0("dead_ploidy_trend_", label, ".png"))
  plot_dead_ploidy_trend(dead_trend, dead_plot_path)
  
  # G vs ploidy trend
  G_trend <- compute_G_ploidy_trend(states)
  G_plot_path <- file.path(analysis_dir,
                           paste0("G_ploidy_trend_", label, ".png"))
  plot_G_ploidy_trend(G_trend, G_plot_path)
  
  invisible(list(states = states,
                 dead_trend = dead_trend,
                 G_trend = G_trend))
}

# --------- Analyze all combos ---------

analyze_all_abm_results <- function(
    base_dir = "/Users/4482173/Documents/IMO_workshop13/ABM_results"
) {
  combos <- find_param_combos(base_dir)
  if (length(combos) == 0L) {
    stop("No parameter combinations (WGDp_*/WGDr_*/Coxy_*/Rep_*) found under ", base_dir)
  }
  
  n_cores <- max(1L, parallel::detectCores() - 1L)
  message("Using ", n_cores, " cores for analysis.")
  
  res <- parallel::mclapply(combos, function(cb) {
    message("Analyzing combo: ", make_param_label(cb),
            " (path: ", cb$path, ")")
    analyze_one_combo(cb, base_dir = base_dir)
  }, mc.cores = n_cores)
  
  invisible(res)
}

# If this script is run via Rscript, run everything:
if (identical(environment(), globalenv()) && !interactive()) {
  analyze_all_abm_results()
}
