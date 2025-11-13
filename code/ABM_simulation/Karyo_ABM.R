# =========================================
# ABM for tumor growth under hypoxia
# R layer = thin wrapper:
#   1) Read YAML config
#   2) Initialize state (grid, O2, cells, vessel_mask, cfg, step, dead_log)
#   3) Call C++ core (simulate_step_cpp) in a loop
#   4) Save per‑step RDS snapshots
# All cell dynamics / oxygen update / division / death are handled in abm_core.cpp
# =========================================

suppressPackageStartupMessages({
  library(yaml)
  library(MASS)   # fitdistr
  library(Rcpp)
  library(parallel)
})

# -------------------------------
# Genome weights for ploidy (must exist in R_GlobalEnv for C++)
# -------------------------------

# Human autosome lengths (GRCh38), chromosomes 1..22 (bp)
CHR_LEN <- c(
  248956422, 242193529, 198295559, 190214555, 181538259, 170805979,
  159345973, 145138636, 138394717, 133797422, 135086622, 133275309,
  114364328, 107043718, 101991189, 90338345, 83257441, 80373285,
  58617616, 64444167, 46709983, 50818468
)
# Normalized weights that sum to 1 (used for weighted ploidy)
CHR_W <- CHR_LEN / sum(CHR_LEN)

# -------------------------------
# Config I/O
# -------------------------------

read_config <- function(path) {
  cfg <- yaml::read_yaml(path)

  required <- c("N1","R","m","MSR","alpha","Db","beta","Dv","Coxy",
                "WGDp","WGDr","KaryoPath","Dtime")
  miss <- setdiff(required, names(cfg))
  if (length(miss) > 0) {
    stop("Missing config keys: ", paste(miss, collapse = ", "))
  }

  # basic sanity clamps
  cfg$alpha <- max(-1, min(1, cfg$alpha))
  cfg$beta  <- max(1e-6, min(1, cfg$beta))

  # Coxy may be scalar/sequence; coerce to numeric vector and map to [0,1]
  cx <- cfg$Coxy
  if (is.list(cx)) cx <- unlist(cx, recursive = TRUE, use.names = FALSE)
  cx <- suppressWarnings(as.numeric(cx))
  if (length(cx) == 0L || all(is.na(cx))) {
    stop("Coxy must be numeric (scalar or numeric vector)")
  }
  cx <- cx[is.finite(cx)]
  if (length(cx) == 0L) stop("Coxy contains no finite numeric values")

  if (length(cx) > 1L) {
    rng <- range(cx, finite = TRUE)
    if ((rng[2] - rng[1]) > 0) {
      cfg$Coxy <- (cx - rng[1]) / (rng[2] - rng[1])
    } else {
      cfg$Coxy <- rep(0, length(cx))
    }
  } else {
    cfg$Coxy <- max(0, min(1, cx))
  }

  cfg$Dv    <- max(0, min(1, cfg$Dv))
  cfg$WGDp  <- max(0, min(1, cfg$WGDp))
  
  wd <- cfg$WGDr
  if (is.list(wd)) wd <- unlist(wd, recursive = TRUE, use.names = FALSE)
  wd <- suppressWarnings(as.numeric(wd))
  if (length(wd) == 0L || all(is.na(wd))) {
    stop("WGDr must be numeric (scalar or numeric vector)")
  }
  wd <- wd[is.finite(wd)]
  if (length(wd) == 0L) stop("WGDr contains no finite numeric values")
  cfg$WGDr <- wd
  
  cfg$m     <- max(0, min(100, cfg$m))

  # MSR interpreted as probability per chromosome per mitosis
  # Normalize to probability if given as percentage (>1)
  if (!is.null(cfg$MSR) && is.finite(cfg$MSR) && cfg$MSR > 1) {
    cfg$MSR <- cfg$MSR / 100
  }
  cfg$MSR <- max(0, cfg$MSR)

  # Dtime: maximum allowed time since last division before death (hours)
  if (is.null(cfg$Dtime)) stop("Config must include 'Dtime'.")
  cfg$Dtime <- as.numeric(cfg$Dtime)
  if (!is.finite(cfg$Dtime)) stop("'Dtime' must be numeric.")
  cfg$Dtime <- max(0, cfg$Dtime)

  # SimsPerCoxy: how many independent simulations per Coxy value (default = 1)
  if (is.null(cfg$SimsPerCoxy)) {
    cfg$SimsPerCoxy <- 1L
  } else {
    cfg$SimsPerCoxy <- as.integer(cfg$SimsPerCoxy)
    if (!is.finite(cfg$SimsPerCoxy) || cfg$SimsPerCoxy < 1L) cfg$SimsPerCoxy <- 1L
  }

  # PloidyMax: per‑chromosome upper bound. Use Inf for no cap.
  if (is.null(cfg$PloidyMax)) {
    cfg$PloidyMax <- Inf
  } else if (is.character(cfg$PloidyMax) &&
             tolower(cfg$PloidyMax) %in% c("inf", ".inf", "infinity")) {
    cfg$PloidyMax <- Inf
  } else {
    cfg$PloidyMax <- suppressWarnings(as.numeric(cfg$PloidyMax))
    if (!is.finite(cfg$PloidyMax) || cfg$PloidyMax <= 0) {
      cfg$PloidyMax <- Inf
    } else {
      cfg$PloidyMax <- floor(cfg$PloidyMax)
    }
  }

  # Quiescent death threshold (hours): default 72 if unspecified
  if (is.null(cfg$QuiescentDeathHours)) {
    cfg$QuiescentDeathHours <- 72
  } else {
    cfg$QuiescentDeathHours <- suppressWarnings(as.numeric(cfg$QuiescentDeathHours))
    if (!is.finite(cfg$QuiescentDeathHours) || cfg$QuiescentDeathHours < 0) {
      cfg$QuiescentDeathHours <- 72
    }
  }

  # Oxygen PDE defaults (used by C++ core)
  if (is.null(cfg$O2DiffRate))    cfg$O2DiffRate    <- 0.2
  if (is.null(cfg$O2SupplyRate))  cfg$O2SupplyRate  <- 0.05
  if (is.null(cfg$O2Consume))     cfg$O2Consume     <- 0.02
  if (is.null(cfg$O2JacobiIters)) cfg$O2JacobiIters <- 1L else {
    cfg$O2JacobiIters <- as.integer(cfg$O2JacobiIters)
    if (!is.finite(cfg$O2JacobiIters) || cfg$O2JacobiIters < 1L) cfg$O2JacobiIters <- 1L
  }

  # Boundary vs vessel mode flags (C++ side can choose behavior)
  if (is.null(cfg$O2UseBoundary))   cfg$O2UseBoundary   <- TRUE
  if (is.null(cfg$O2BoundaryMode))  cfg$O2BoundaryMode  <- "dirichlet_edges"
  if (is.null(cfg$O2BoundaryValue)) {
    cfg$O2BoundaryValue <- 1.0
  } else {
    cfg$O2BoundaryValue <- suppressWarnings(as.numeric(cfg$O2BoundaryValue))
    if (!is.finite(cfg$O2BoundaryValue)) cfg$O2BoundaryValue <- 1.0
    cfg$O2BoundaryValue <- max(0, min(1, cfg$O2BoundaryValue))
  }
  if (is.null(cfg$O2VesselSoftSupply)) cfg$O2VesselSoftSupply <- FALSE

  cfg
}

# -------------------------------
# Load K and fit LogNormal to log(K)
# -------------------------------

load_K_column <- function(path) {
  try_read <- function(p) {
    tryCatch(
      read.csv(p, stringsAsFactors = FALSE),
      error = function(e) {
        tryCatch(
          read.delim(p, stringsAsFactors = FALSE),
          error = function(e2) stop("Failed to read KaryoPath: ", e2$message)
        )
      }
    )
  }
  df <- try_read(path)
  if (!("K" %in% colnames(df))) stop("KaryoPath must contain a numeric column named 'K'.")
  K <- suppressWarnings(as.numeric(df$K))
  K <- K[is.finite(K) & K > 0]
  if (length(K) < 10) stop("Not enough positive finite K values to fit a LogNormal.")
  K
}

fit_lognormal_from_K <- function(K) {
  fit <- MASS::fitdistr(log(K), "normal")
  list(mu = fit$estimate[1], sigma = fit$estimate[2])
}

sample_P_from_fit <- function(n, fit, clamp = c(0, 10)) {
  P <- exp(rnorm(n, mean = fit$mu, sd = fit$sigma))
  P <- pmin(pmax(P, clamp[1]), clamp[2])
  P
}

# -------------------------------
# Vessel oxygen sources (geometry only)
# -------------------------------

make_vessel_mask <- function(N, centers = NULL, diameter = 10L) {
  r <- as.integer(ceiling(as.numeric(diameter) / 2))
  if (is.null(centers)) {
    # Default: five vessels (center + four cardinal directions)
    clamp_inside <- function(ci, cj) {
      ci <- as.integer(max(1L + r, min(N - r, ci)))
      cj <- as.integer(max(1L + r, min(N - r, cj)))
      c(ci, cj)
    }
    half <- as.integer(round(N / 2))
    q1   <- as.integer(round(N / 4))
    q3   <- as.integer(round(3 * N / 4))
    centers <- list(
      clamp_inside(half, half),
      clamp_inside(q1,  half),
      clamp_inside(q3,  half),
      clamp_inside(half, q1),
      clamp_inside(half, q3)
    )
  }
  mask <- matrix(0L, nrow = N, ncol = N)
  for (cij in centers) {
    ci <- as.integer(cij[1]); cj <- as.integer(cij[2])
    if (!is.finite(ci) || !is.finite(cj)) {
      stop(sprintf("make_vessel_mask: invalid center with NA/NaN: (%s)", paste(cij, collapse = ",")))
    }
    if (ci < 1L || ci > N || cj < 1L || cj > N) {
      stop(sprintf("make_vessel_mask: center out of bounds after clamping: (%d,%d)", ci, cj))
    }
    imin <- max(1L, ci - r); imax <- min(N, ci + r)
    jmin <- max(1L, cj - r); jmax <- min(N, cj + r)
    for (ii in imin:imax) {
      for (jj in jmin:jmax) {
        if ((ii - ci)^2 + (jj - cj)^2 <= r^2) mask[ii, jj] <- 1L
      }
    }
  }
  mask
}

resolve_vessel_centers <- function(N, centers_cfg, diameter) {
  r <- as.integer(ceiling(as.numeric(diameter) / 2))
  clamp_inside <- function(xi, xj) {
    ci <- as.integer(max(1L + r, min(N - r, xi)))
    cj <- as.integer(max(1L + r, min(N - r, xj)))
    c(ci, cj)
  }
  key_to_pair <- function(key) {
    key <- tolower(trimws(key))
    key <- gsub("[ _-]+", "", key)
    if (key %in% c("center","centre","middle","mid")) {
      return(clamp_inside(round(N/2), round(N/2)))
    } else if (key %in% c("topleft","ul")) {
      return(clamp_inside(1L, 1L))
    } else if (key %in% c("topright","ur")) {
      return(clamp_inside(1L, N))
    } else if (key %in% c("bottomleft","bl")) {
      return(clamp_inside(N, 1L))
    } else if (key %in% c("bottomright","br")) {
      return(clamp_inside(N, N))
    } else {
      stop(sprintf("Unsupported vessel center keyword: %s", key))
    }
  }
  parse_single <- function(tok) {
    tok <- trimws(as.character(tok))
    if (!nzchar(tok)) stop("Empty token in VesselCenters")
    if (grepl("^\\s*[-+]?[0-9]+[ ,;]+[-+]?[0-9]+\\s*$", tok)) {
      nums <- as.integer(strsplit(tok, "[ ,;]+", perl = TRUE)[[1]])
      return(clamp_inside(nums[1], nums[2]))
    }
    key_to_pair(tok)
  }
  as_pair <- function(v) {
    if (is.list(v) && !is.null(names(v))) {
      nms <- tolower(names(v))
      if (all(c("i","j") %in% nms)) {
        return(clamp_inside(
          as.integer(v[[which(nms == "i")[1]]]),
          as.integer(v[[which(nms == "j")[1]]])
        ))
      }
      if (all(c("x","y") %in% nms)) {
        return(clamp_inside(
          as.integer(v[[which(nms == "x")[1]]]),
          as.integer(v[[which(nms == "y")[1]]])
        ))
      }
      v <- unlist(v, use.names = FALSE)
    }
    if (is.character(v)) {
      if (length(v) == 1L) return(parse_single(v))
      stop("Character vector must be normalized before as_pair()")
    }
    if (is.numeric(v)) {
      vals <- suppressWarnings(as.integer(v))
      if (length(vals) < 2L || any(!is.finite(vals[1:2]))) {
        stop("Numeric pair in VesselCenters must have two finite integers")
      }
      return(clamp_inside(vals[1L], vals[2L]))
    }
    stop("VesselCenters entries must be keywords, numeric pairs, or named lists with i/j or x/y")
  }

  if (is.null(centers_cfg)) return(NULL)

  centers_list <- NULL
  if (is.character(centers_cfg) && length(centers_cfg) == 1L) {
    toks <- trimws(unlist(strsplit(centers_cfg, "[ ,;]+", perl = TRUE)))
    toks <- toks[nzchar(toks)]
    centers_list <- lapply(as.list(toks), parse_single)
  } else if (is.character(centers_cfg) && length(centers_cfg) > 1L) {
    centers_list <- lapply(as.list(centers_cfg), parse_single)
  } else if (is.numeric(centers_cfg)) {
    vals <- as.integer(centers_cfg)
    if ((length(vals) %% 2L) != 0L) stop("Numeric VesselCenters length must be even (pairs of i,j)")
    centers_list <- lapply(
      seq(1L, length(vals), by = 2L),
      function(k) clamp_inside(vals[k], vals[k + 1L])
    )
  } else if (is.list(centers_cfg)) {
    flat <- list()
    for (el in centers_cfg) {
      if (is.character(el) && length(el) > 1L) {
        for (tk in el) flat[[length(flat) + 1L]] <- tk
      } else {
        flat[[length(flat) + 1L]] <- el
      }
    }
    centers_list <- lapply(flat, as_pair)
  } else {
    stop(sprintf("Unsupported VesselCenters type: %s", class(centers_cfg)[1]))
  }

  for (idx in seq_along(centers_list)) {
    ci <- centers_list[[idx]][1]; cj <- centers_list[[idx]][2]
    if (!is.finite(ci) || !is.finite(cj)) {
      stop(sprintf("Invalid vessel center at index %d (NA/NaN)", idx))
    }
    if (ci < 1L || ci > N || cj < 1L || cj > N) {
      stop(sprintf("Vessel center at index %d is out of grid after clamping: (%d,%d)", idx, ci, cj))
    }
  }
  centers_list
}

# -------------------------------
# Data containers
# -------------------------------

make_cells_df <- function(n_rows) {
  c_names <- c(
    "X","Y",
    paste0("C", 1:22),
    "R","D","P","G",
    "WGDLabel","WGDCount",
    "DivisionTime","Time","Status",
    "Quiescent","QuiescentTime",
    "DeathReason"
  )
  df <- as.data.frame(matrix(NA_real_, nrow = n_rows, ncol = length(c_names)))
  colnames(df) <- c_names

  int_cols <- c("X","Y", paste0("C", 1:22),
                "WGDLabel","WGDCount","Status","Quiescent")
  for (nm in int_cols) df[[nm]] <- as.integer(df[[nm]])

  df$Quiescent      <- ifelse(is.na(df$Quiescent), 0L, df$Quiescent)
  df$QuiescentTime  <- ifelse(is.na(df$QuiescentTime), 0.0, df$QuiescentTime)
  df$WGDCount       <- ifelse(is.na(df$WGDCount), 0L, df$WGDCount)
  df$DeathReason    <- NA_character_

  df
}

init_dead_log <- function(cells_template) {
  df <- cells_template[0, , drop = FALSE]
  if (!("DeathReason" %in% names(df))) df$DeathReason <- character(0)
  df$Step <- integer(0)
  df
}

# -------------------------------
# Initialization (R -> initial state; C++ handles dynamics)
# -------------------------------

init_simulation <- function(cfg, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  N <- as.integer(cfg$N1)
  grid <- matrix(NA_integer_, nrow = N, ncol = N)

  # Determine scalar oxygen supply for this run
  if (is.null(cfg$Coxy_scalar)) {
    cfg$Coxy_scalar <- if (length(cfg$Coxy) <= 1L) {
      as.numeric(cfg$Coxy)
    } else {
      mean(as.numeric(cfg$Coxy), na.rm = TRUE)
    }
  }

  # Vessel geometry
  diam <- if (!is.null(cfg$VesselDiameter)) as.integer(cfg$VesselDiameter) else 10L
  centers <- NULL
  if (!is.null(cfg$VesselCenters)) {
    centers <- resolve_vessel_centers(N, cfg$VesselCenters, diam)
  }
  vessel_mask <- make_vessel_mask(N, centers = centers, diameter = diam)

  # In vessel mode, default to no boundary source (C++ may still honor cfg if needed)
  cfg$O2UseBoundary   <- FALSE
  cfg$O2BoundaryMode  <- "neumann"
  cfg$O2BoundaryValue <- 0.0

  # Initial O2 field: uniform Coxy_scalar; C++ will evolve it with vessels
  O2 <- matrix(cfg$Coxy_scalar, nrow = N, ncol = N)

  # Initialize karyotypes from K distribution
  K <- load_K_column(cfg$KaryoPath)
  fit <- fit_lognormal_from_K(K)

  # Initial seeding density: m% of non-vessel grid
  n_init <- max(1L, round(N * N * (cfg$m / 100)))
  allowed <- which(as.vector(vessel_mask) == 0L)
  if (length(allowed) < n_init) {
    stop("Not enough non-vessel locations to place initial cells. Reduce m or vessel diameter.")
  }
  pos <- sample(allowed, n_init, replace = FALSE)
  X <- as.integer((pos - 1L) %% N + 1L)
  Y <- as.integer((pos - 1L) %/% N + 1L)

  cells <- make_cells_df(n_init)
  cells$X <- X
  cells$Y <- Y

  # Initialize P from fitted LogNormal, clamp to configured ploidy cap
  pmax_cfg <- if (is.finite(cfg$PloidyMax)) as.numeric(cfg$PloidyMax) else Inf
  cells$P <- sample_P_from_fit(n_init, fit = fit, clamp = c(0, pmax_cfg))

  # Invert P into integer karyotype (C++ helper), with per-chromosome cap
  max_copy_cpp <- if (is.finite(cfg$PloidyMax)) {
    as.integer(cfg$PloidyMax)
  } else {
    .Machine$integer.max
  }
  kt_mat <- karyotype_from_P_batch_weighted_cpp(
    as.numeric(cells$P),
    as.numeric(CHR_W),
    1e-3,
    20000L,
    max_copy_cpp
  )

  # Initial karyotypes: enforce a softer cap at 5 unless PloidyMax is lower
  kt_init_cap <- 5L
  if (is.finite(cfg$PloidyMax)) {
    kt_init_cap <- as.integer(min(kt_init_cap, cfg$PloidyMax))
  }
  kt_mat <- pmin(kt_mat, kt_init_cap)
  storage.mode(kt_mat) <- "integer"

  cells[, paste0("C", 1:22)] <- kt_mat

  # Backfill P with actual integer-weighted mean (C++ helper)
  cells$P <- mean_copy_number_rows_weighted_cpp(kt_mat, as.numeric(CHR_W))

  # Scalar attributes
  cells$R <- cfg$R
  cells$D <- cfg$Db
  cells$WGDLabel <- 0L
  if (is.null(cells$WGDCount)) cells$WGDCount <- integer(nrow(cells))
  cells$DivisionTime <- 0.0
  cells$Time <- 0.0
  cells$Status <- 1L
  if (is.null(cells$DeathReason)) cells$DeathReason <- NA_character_

  # Initial WGD proportion WGDp
  if (cfg$WGDp > 0) {
    k <- max(0L, floor(n_init * cfg$WGDp))
    if (k > 0) {
      idx <- sample(seq_len(n_init), k, replace = FALSE)
      for (i in idx) {
        v <- as.integer(cells[i, paste0("C", 1:22)])
        if (is.finite(cfg$PloidyMax)) {
          v <- pmin(as.integer(cfg$PloidyMax), v * 2L)
        } else {
          v <- v * 2L
        }
        cells[i, paste0("C", 1:22)] <- v
        cells$WGDLabel[i] <- 1L
        cells$P[i] <- mean_copy_number_weighted_cpp(as.integer(v), as.numeric(CHR_W))
        cells$WGDCount[i] <- as.integer(cells$WGDCount[i]) + 1L
      }
    }
  }

  # Place cells on grid (C++ helper)
  grid <- rebuild_grid_from_cells_cpp(grid, as.integer(cells$X), as.integer(cells$Y))

  # Initial G & DivisionTime via C++ batch; randomize Time in [0,24]
  res0 <- compute_G_and_Div_cpp(
    grid,
    O2,
    as.integer(cells$X),
    as.integer(cells$Y),
    as.numeric(cells$R),
    as.numeric(cells$P),
    as.numeric(cfg$beta),
    as.integer(121L),
    as.integer(5L)
  )
  cells$G <- as.numeric(res0$G)
  cells$DivisionTime <- as.numeric(res0$DivisionTime)
  cells$Time <- runif(n_init, min = 0, max = 24)

  list(
    grid        = grid,
    O2          = O2,
    cells       = cells,
    cfg         = cfg,
    step        = 0L,
    vessel_mask = vessel_mask,
    dead_log    = init_dead_log(cells)
  )
}

# -------------------------------
# Snapshot saving helpers
# -------------------------------

ensure_data_dir <- function(dir = "data") {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  dir
}

save_cells_snapshot <- function(state, t, dir = "data") {
  ensure_data_dir(dir)
  fname <- file.path(dir, sprintf("Sim_t%06d.rds", as.integer(t)))
  saveRDS(state, fname)
  invisible(fname)
}

# -------------------------------
# Thin wrapper to C++ core
# -------------------------------

simulate_step <- function(state) {
  simulate_step_cpp(state)
}

run_simulation <- function(state,
                           steps = 240L,
                           verbose_every = 24L,
                           out_dir = "data") {
  for (t in seq_len(steps)) {
    state <- simulate_step(state)
    save_cells_snapshot(state, t, dir = out_dir)

    if (!is.null(verbose_every) &&
        verbose_every > 0L &&
        (t %% verbose_every == 0L)) {
      alive <- sum(state$cells$Status == 1L)
      meanG <- if (alive > 0L) {
        mean(state$cells$G[state$cells$Status == 1L], na.rm = TRUE)
      } else {
        NA_real_
      }
      cat(sprintf(
        "[t=%4dh] alive=%d, meanG=%s\n",
        state$step, alive,
        ifelse(is.na(meanG), "NA", sprintf("%.4f", meanG))
      ))
    }
  }
  state
}


# -------------------------------
# Batch simulation driver
# -------------------------------
# Wrap the full "build jobs -> run in parallel -> print summary" pipeline
# into a reusable function. The caller must provide:
#   - cfg: list from read_config()
#   - cpp_path: path to abm_core.cpp for sourceCpp()
#   - steps: number of simulation steps per run
run_abm_simulations <- function(cfg, cpp_path, steps = 720L) {
  if (!file.exists(cpp_path)) {
    stop("cpp_path not found: ", cpp_path)
  }

  # Compile / load C++ core
  Rcpp::sourceCpp(cpp_path)

  # Output base directory from config (fallback to "data")
  out_base_dir <- if (!is.null(cfg$OutDir) && nzchar(cfg$OutDir)) cfg$OutDir else "data"

  if (!is.numeric(cfg$Coxy)) stop("cfg$Coxy must be numeric after read_config")

  coxy_vec_full <- as.numeric(cfg$Coxy)
  # Keep 0 inside cfg$Coxy for scaling/mapping; skip zero in actual simulation runs
  coxy_vec <- coxy_vec_full[is.finite(coxy_vec_full) & coxy_vec_full > 0]
  if (length(coxy_vec) == 0L) {
    stop("After excluding Coxy==0, no positive Coxy values remain to simulate.")
  }

  n_runs <- length(coxy_vec)
  S      <- as.integer(cfg$SimsPerCoxy)

  # Build layered job list: ri = 1..S, for each ci = 1..n_runs
  jobs <- vector("list", n_runs * S)
  idx  <- 1L
  for (ri in seq_len(S)) {
    for (ci in seq_len(n_runs)) {
      jobs[[idx]] <- list(ci = ci, ri = ri)
      idx <- idx + 1L
    }
  }

  run_job <- function(job) {
    tryCatch({
      ci <- job$ci
      ri <- job$ri

      cfg_run <- cfg
      cfg_run$Coxy_scalar <- coxy_vec[ci]

      coxy_dir <- sprintf("Coxy_%0.6f", cfg_run$Coxy_scalar)
      out_dir  <- file.path(out_base_dir, coxy_dir, sprintf("Rep_%03d", ri))

      state <- init_simulation(cfg_run, seed = NULL)
      state <- run_simulation(state, steps = steps, out_dir = out_dir)

      list(
        status = "ok",
        ci     = ci,
        ri     = ri,
        coxy   = cfg_run$Coxy_scalar,
        out_dir = out_dir,
        alive   = sum(state$cells$Status == 1L)
      )
    }, error = function(e) {
      list(
        status = "error",
        ci     = job$ci,
        ri     = job$ri,
        coxy   = if (job$ci >= 1 && job$ci <= length(coxy_vec))
                   coxy_vec[job$ci] else NA_real_,
        out_dir = NA_character_,
        alive   = NA_integer_,
        error   = conditionMessage(e)
      )
    })
  }

  if (.Platform$OS.type == "windows") {
    res_list <- lapply(jobs, run_job)
  } else {
    ncores <- max(1L, parallel::detectCores(logical = TRUE) - 1L)
    res_list <- parallel::mclapply(
      jobs, run_job,
      mc.cores       = ncores,
      mc.preschedule = FALSE
    )
  }

  invisible(lapply(seq_along(res_list), function(i) {
    res <- res_list[[i]]
    if (is.list(res) && !is.null(res$coxy)) {
      msg_err <- if (!is.null(res$status) && res$status == "error") {
        sprintf(" [ERROR: %s]", res$error)
      } else {
        ""
      }
      cat(sprintf(
        "\n[RUN Coxy=%0.6f rep=%03d] finished. Alive cells: %s%s\n",
        as.numeric(res$coxy),
        as.integer(res$ri),
        ifelse(is.na(res$alive), "NA", as.character(res$alive)),
        msg_err
      ))
    } else {
      cat(sprintf("\n[RUN #%d] unexpected result type: %s\n", i, typeof(res)))
    }
  }))

  invisible(res_list)
}

