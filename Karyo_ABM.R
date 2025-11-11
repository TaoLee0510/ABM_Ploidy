# =========================================
# ABM for tumor growth under hypoxia



suppressPackageStartupMessages({
  library(yaml)
  library(MASS)   # fitdistr
  library(Rcpp)
  library(parallel)
})

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
  required <- c("N1","R","m","MSR","alpha","Db","beta","Dv","Coxy","WGDp","WGDr","KaryoPath","Dtime")
  miss <- setdiff(required, names(cfg))
  if (length(miss) > 0) stop("Missing config keys: ", paste(miss, collapse=", "))
  # basic sanity clamps
  cfg$alpha <- max(-1,min(1,cfg$alpha))
  cfg$beta  <- max(1e-6,min(1,cfg$beta))
  # Coxy may be scalar/sequence; coerce to numeric vector and map to [0,1]
  cx <- cfg$Coxy
  # YAML sequences can arrive as lists; unlist and coerce
  if (is.list(cx)) cx <- unlist(cx, recursive = TRUE, use.names = FALSE)
  cx <- suppressWarnings(as.numeric(cx))
  if (length(cx) == 0L || all(is.na(cx))) stop("Coxy must be numeric (scalar or numeric vector)")
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
  cfg$Dv    <- max(0,min(1,cfg$Dv))
  cfg$WGDp  <- max(0,min(1,cfg$WGDp))
  cfg$WGDr  <- max(0,min(1,cfg$WGDr))
  cfg$m     <- max(0,min(100,cfg$m))
  # MSR now interpreted as % per chromosome per mitosis
  # Normalize to probability if given as percentage (>1)
  if (!is.null(cfg$MSR) && is.finite(cfg$MSR) && cfg$MSR > 1) cfg$MSR <- cfg$MSR / 100
  cfg$MSR <- max(0, cfg$MSR)
  # Dtime: maximum allowed Time since last division before death (in hours)
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
  } else {
    if (is.character(cfg$PloidyMax) && tolower(cfg$PloidyMax) %in% c("inf", ".inf", "infinity")) {
      cfg$PloidyMax <- Inf
    } else {
      cfg$PloidyMax <- suppressWarnings(as.numeric(cfg$PloidyMax))
      if (!is.finite(cfg$PloidyMax) || cfg$PloidyMax <= 0) cfg$PloidyMax <- Inf else cfg$PloidyMax <- floor(cfg$PloidyMax)
    }
  }
  # Quiescent death threshold (hours): default 72 if unspecified
  if (is.null(cfg$QuiescentDeathHours)) {
    cfg$QuiescentDeathHours <- 72
  } else {
    cfg$QuiescentDeathHours <- suppressWarnings(as.numeric(cfg$QuiescentDeathHours))
    if (!is.finite(cfg$QuiescentDeathHours) || cfg$QuiescentDeathHours < 0) cfg$QuiescentDeathHours <- 72
  }

  # Default to boundary-driven inward diffusion (Dirichlet edges)
  if (is.null(cfg$O2UseBoundary)) cfg$O2UseBoundary <- TRUE
  if (is.null(cfg$O2BoundaryMode)) cfg$O2BoundaryMode <- "dirichlet_edges"  # "neumann" or "dirichlet_edges"
  if (is.null(cfg$O2BoundaryValue)) {
    cfg$O2BoundaryValue <- 1.0
  } else {
    cfg$O2BoundaryValue <- suppressWarnings(as.numeric(cfg$O2BoundaryValue))
    if (!is.finite(cfg$O2BoundaryValue)) cfg$O2BoundaryValue <- 1.0
    cfg$O2BoundaryValue <- max(0, min(1, cfg$O2BoundaryValue))
  }
  # Per-step PDE coefficients (tunable via YAML)
  if (is.null(cfg$O2DiffRate))   cfg$O2DiffRate   <- 0.2
  if (is.null(cfg$O2SupplyRate)) cfg$O2SupplyRate <- 0.05
  if (is.null(cfg$O2Consume))    cfg$O2Consume    <- 0.02
  # PDE sub-iterations per hour (speeds up apparent diffusion range)
  if (is.null(cfg$O2JacobiIters)) {
    cfg$O2JacobiIters <- 1L
  } else {
    cfg$O2JacobiIters <- as.integer(cfg$O2JacobiIters)
    if (!is.finite(cfg$O2JacobiIters) || cfg$O2JacobiIters < 1L) cfg$O2JacobiIters <- 1L
  }
  # Optional soft supply term in vessel mode (disabled by default to preserve vessel-only physics)
  if (is.null(cfg$O2VesselSoftSupply)) cfg$O2VesselSoftSupply <- FALSE
  cfg
}

# -------------------------------
# Load K and fit LogNormal to log(K)
# -------------------------------
load_K_column <- function(path) {
  try_read <- function(p) {
    tryCatch(read.csv(p, stringsAsFactors = FALSE),
             error = function(e) tryCatch(read.delim(p, stringsAsFactors = FALSE),
                                          error = function(e2) stop("Failed to read KaryoPath: ", e2$message)))
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

sample_P_from_fit <- function(n, fit, clamp = c(0,10)) {
  P <- exp(rnorm(n, mean = fit$mu, sd = fit$sigma))
  P <- pmin(pmax(P, clamp[1]), clamp[2])
  P
}

# -------------------------------
# P <-> karyotype utilities
# -------------------------------
weighted_P_from_karyotype <- function(kt, w = NULL) {
  if (is.null(w)) w <- CHR_W
  sum(kt * w) / sum(w)
}

karyotype_from_P <- function(P, w = NULL, tol = 1e-3, max_iter = 1e5L) {
  # Build an integer karyotype (22 chromosomes in [0,10]) to approximate P.
  # Start from diploid (all 2s) and greedily nudge +/-1.
  chr_n <- 22L
  kt <- rep(2L, chr_n)
  if (is.null(w)) w <- CHR_W
  currP <- weighted_P_from_karyotype(kt, w)

  iter <- 0L
  best <- list(P = currP, kt = kt, err = abs(currP - P))

  while (iter < max_iter) {
    iter <- iter + 1L
    dir <- ifelse(currP < P, +1L, -1L)
    idx <- sample.int(chr_n, chr_n)
    improved <- FALSE
    for (j in idx) {
      newval <- kt[j] + dir
      if (newval >= 0L && newval <= 10L) {
        kt2 <- kt; kt2[j] <- newval
        P2 <- weighted_P_from_karyotype(kt2, w)
        if (abs(P2 - P) + 1e-12 < abs(currP - P)) {
          kt <- kt2; currP <- P2; improved <- TRUE
          if (abs(currP - P) < best$err) best <- list(P=currP, kt=kt, err=abs(currP-P))
          if (abs(currP - P) <= tol) break
        }
      }
    }
    if (!improved) break
  }
  best$kt
}
Rcpp::cppFunction(code='
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
double mean_copy_number_weighted_cpp(IntegerVector kt, NumericVector w){
  int n = kt.size();
  if (w.size() != n) stop("mean_copy_number_weighted_cpp: weight length mismatch");
  long double num = 0.0L, den = 0.0L;
  for(int i=0;i<n;++i){ num += (long double)kt[i] * (long double)w[i]; den += (long double)w[i]; }
  if (den <= 0.0L) return NA_REAL;
  return (double)(num / den);
}
')
Rcpp::cppFunction(code='
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
NumericVector mean_copy_number_rows_weighted_cpp(IntegerMatrix kt_mat, NumericVector w){
  int n = kt_mat.nrow();
  int m = kt_mat.ncol();
  if (w.size() != m) stop("mean_copy_number_rows_weighted_cpp: weight length mismatch");
  NumericVector out(n);
  long double den = 0.0L;
  for(int j=0;j<m;++j) den += (long double)w[j];
  if (den <= 0.0L) {
    for(int i=0;i<n;++i) out[i] = NA_REAL;
    return out;
  }
  for(int i=0;i<n;++i){
    long double num = 0.0L;
    for(int j=0;j<m;++j){ num += (long double)kt_mat(i,j) * (long double)w[j]; }
    out[i] = (double)(num / den);
  }
  return out;
}
')
Rcpp::cppFunction(code='
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
IntegerVector karyotype_from_P_weighted_cpp(double P, NumericVector w, double tol=1e-3, int max_iter=20000, int max_copy=-1){
  const int chr_n = w.size();
  if (chr_n != 22) Rcpp::stop("karyotype_from_P_weighted_cpp expects 22 weights (autosomes 1..22)");
  IntegerVector kt(chr_n); for(int i=0;i<chr_n;++i) kt[i]=2;
  auto wsum = [&]()->long double{ long double s=0.0L; for(int t=0;t<chr_n;++t) s += (long double)w[t]; return s; };
  auto meanW = [&](const IntegerVector& v)->double{
    long double num=0.0L, den=wsum();
    for(int t=0;t<chr_n;++t) num += (long double)v[t] * (long double)w[t];
    if (den<=0.0L) return NA_REAL;
    return (double)(num/den);
  };
  long double SW = wsum();
  if (SW <= 0.0L) return kt;
  double currP = meanW(kt); double err = std::fabs(currP - P); int iter=0;
  while (iter < max_iter && err > tol){
    int dir = (currP < P) ? +1 : -1;
    int best_j = -1; double best_err = err;
    for(int j=0;j<chr_n;++j){
      int nv = kt[j] + dir;
      if (nv < 0) continue;
      if (max_copy > 0 && nv > max_copy) continue;
      double newP = currP + dir * ((double)w[j] / (double)SW);
      double e2 = std::fabs(newP - P);
      if (e2 + 1e-12 < best_err){ best_err = e2; best_j = j; }
    }
    if (best_j < 0) break;
    kt[best_j] = kt[best_j] + dir;
    currP = currP + dir * ((double)w[best_j] / (double)SW);
    err = std::fabs(currP - P);
    ++iter;
  }
  return kt;
}
')
Rcpp::cppFunction(code='
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
IntegerMatrix karyotype_from_P_batch_weighted_cpp(NumericVector Pvec, NumericVector w, double tol=1e-3, int max_iter=20000, int max_copy=-1){
  const int chr_n = w.size();
  if (chr_n != 22) Rcpp::stop("karyotype_from_P_batch_weighted_cpp expects 22 weights");
  int n = Pvec.size();
  IntegerMatrix out(n, chr_n);
  long double SW = 0.0L; for(int t=0;t<chr_n;++t) SW += (long double)w[t];
  if (SW <= 0.0L) return out;
  for (int i=0;i<n;++i){
    double P = Pvec[i];
    IntegerVector kt(chr_n); for(int c=0;c<chr_n;++c) kt[c]=2;
    double currP = 0.0; { long double num=0.0L; for(int t=0;t<chr_n;++t) num += (long double)kt[t]*(long double)w[t]; currP = (double)(num/SW); }
    double err = std::fabs(currP - P); int iter=0;
    while (iter < max_iter && err > tol){
      int dir = (currP < P) ? +1 : -1; int best_j = -1; double best_err = err;
      for(int j=0;j<chr_n;++j){
        int nv = kt[j] + dir;
        if (nv < 0) continue;
        if (max_copy > 0 && nv > max_copy) continue;
        double newP = currP + dir * ((double)w[j] / (double)SW);
        double e2 = std::fabs(newP - P);
        if (e2 + 1e-12 < best_err){ best_err = e2; best_j = j; }
      }
      if (best_j < 0) break;
      kt[best_j] = kt[best_j] + dir;
      currP = currP + dir * ((double)w[best_j] / (double)SW);
      err = std::fabs(currP - P);
      ++iter;
    }
    for (int c=0;c<chr_n;++c) out(i,c) = kt[c];
  }
  return out;
}
')

# -------------------------------
# Oxygen field (explicit diffusion-consumption; lightweight)
# -------------------------------
init_oxygen_field <- function(N, Coxy) {
  # Initialize oxygen in [0,1] as a uniform field near Coxy
  matrix(Coxy, nrow = N, ncol = N)
}

update_oxygen_field <- function(O2, grid, Coxy, diff_rate = 0.2, supply_rate = 0.05, consume = 0.02, iters = 1L) {
  # One or several Jacobi-like smoothing steps with source and consumption.
  # O2 in [0,1]; 'grid' NA means empty; non-NA means occupied -> consumes oxygen.
  N <- nrow(O2)
  for (k in seq_len(iters)) {
    lap <- matrix(0, N, N)
    # 5-point Laplacian (Neumann-like by clamping edges)
    up    <- rbind(O2[1, , drop=FALSE],  O2[-N, ])
    down  <- rbind(O2[-1,],              O2[N,  , drop=FALSE])
    left  <- cbind(O2[,1, drop=FALSE],   O2[,-ncol(O2)])
    right <- cbind(O2[,-1],              O2[,ncol(O2), drop=FALSE])
    lap <- (up + down + left + right - 4*O2)

    occ <- !is.na(grid)
    O2 <- O2 + diff_rate * lap + supply_rate * (Coxy - O2) - consume * occ
    O2 <- pmin(pmax(O2, 0), 1)  # clamp to [0,1]
  }
  O2
}

oxygen_theta <- function(O2, i, j) {
  # theta in [0,1] at cell center
  O2[i, j]
}

# -------------------------------
# Vessel oxygen sources
# Two circular vessels with given diameter (in grid cells). Oxygen supplied at vessel cells only.
# -------------------------------

make_vessel_mask <- function(N, centers = NULL, diameter = 10L) {
  # centers: list of length-2 integer vectors (i, j). Defaults to 5 vessels at center and cardinal quarters.
  r <- as.integer(ceiling(as.numeric(diameter) / 2))
  if (is.null(centers)) {
    # Place 5 vessels by default:
    #  - One at the center (N/2, N/2)
    #  - Four at quarter positions along up/down/left/right:
    #      (N/4, N/2), (3N/4, N/2), (N/2, N/4), (N/2, 3N/4)
    # Clamp each center so the full circle of radius r stays within [1, N]
    clamp_inside <- function(ci, cj) {
      ci <- as.integer(max(1L + r, min(N - r, ci)))
      cj <- as.integer(max(1L + r, min(N - r, cj)))
      c(ci, cj)
    }
    half <- as.integer(round(N / 2))
    q1   <- as.integer(round(N / 4))
    q3   <- as.integer(round(3 * N / 4))
    centers <- list(
      clamp_inside(half, half),  # center
      clamp_inside(q1,  half),   # up (top)
      clamp_inside(q3,  half),   # down (bottom)
      clamp_inside(half, q1),    # left
      clamp_inside(half, q3)     # right
    )
  }
  mask <- matrix(0L, nrow = N, ncol = N)
  for (cij in centers) {
    ci <- as.integer(cij[1]); cj <- as.integer(cij[2])
    if (!is.finite(ci) || !is.finite(cj)) {
      stop(sprintf("make_vessel_mask: encountered invalid center with NA/NaN: (%s)", paste(cij, collapse=",")))
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

# Resolve vessel centers from config; robustly parse multiple formats
resolve_vessel_centers <- function(N, centers_cfg, diameter) {
  r <- as.integer(ceiling(as.numeric(diameter) / 2))
  clamp_inside <- function(xi, xj) {
    # keep full circle inside grid bounds when possible
    ci <- as.integer(max(1L + r, min(N - r, xi)))
    cj <- as.integer(max(1L + r, min(N - r, xj)))
    c(ci, cj)
  }
  # keyword -> pair
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
  # parse a single token that might be a keyword or a numeric pair like "40,50" or "40 50"
  parse_single <- function(tok) {
    tok <- trimws(as.character(tok))
    if (!nzchar(tok)) stop("Empty token in VesselCenters")
    # try numeric pair pattern
    if (grepl("^\\s*[-+]?[0-9]+[ ,;]+[-+]?[0-9]+\\s*$", tok)) {
      nums <- as.integer(strsplit(tok, "[ ,;]+", perl = TRUE)[[1]])
      return(clamp_inside(nums[1], nums[2]))
    }
    # otherwise treat as keyword
    return(key_to_pair(tok))
  }
  # convert one entry to a pair
  as_pair <- function(v) {
    # named list like list(i=10,j=20) or list(x=..., y=...)
    if (is.list(v) && !is.null(names(v))) {
      nms <- tolower(names(v))
      if (all(c("i","j") %in% nms)) {
        return(clamp_inside(as.integer(v[[which(nms=="i")[1]]]), as.integer(v[[which(nms=="j")[1]]])))
      }
      if (all(c("x","y") %in% nms)) {
        return(clamp_inside(as.integer(v[[which(nms=="x")[1]]]), as.integer(v[[which(nms=="y")[1]]])))
      }
      # fall through: try to unlist and treat as numeric
      v <- unlist(v, use.names = FALSE)
    }
    if (is.character(v)) {
      if (length(v) == 1L) return(parse_single(v))
      # if a character vector of multiple tokens, parse each separately (handled in outer normalizer)
      stop("Character vector must be normalized before as_pair()")
    }
    if (is.numeric(v)) {
      vals <- suppressWarnings(as.integer(v))
      if (length(vals) < 2L || any(!is.finite(vals[1:2]))) stop("Numeric pair in VesselCenters must have two finite integers")
      return(clamp_inside(vals[1L], vals[2L]))
    }
    stop("VesselCenters entries must be keywords, numeric pairs, or named lists with i/j or x/y")
  }
  if (is.null(centers_cfg)) return(NULL)
  centers_list <- NULL
  # 1) single string like "topleft, center" -> split tokens
  if (is.character(centers_cfg) && length(centers_cfg) == 1L) {
    toks <- trimws(unlist(strsplit(centers_cfg, "[ ,;]+", perl = TRUE)))
    toks <- toks[nzchar(toks)]
    centers_list <- lapply(as.list(toks), parse_single)
  # 2) character vector c("topleft","center")
  } else if (is.character(centers_cfg) && length(centers_cfg) > 1L) {
    centers_list <- lapply(as.list(centers_cfg), parse_single)
  # 3) flat numeric vector (length 2 -> one center, length 4 -> two centers, etc.)
  } else if (is.numeric(centers_cfg)) {
    vals <- as.integer(centers_cfg)
    if ((length(vals) %% 2L) != 0L) stop("Numeric VesselCenters length must be even (pairs of i,j)")
    centers_list <- lapply(seq(1L, length(vals), by = 2L), function(k) clamp_inside(vals[k], vals[k+1L]))
  # 4) list of items (each may be keyword string, numeric pair, or named list)
  } else if (is.list(centers_cfg)) {
    # flatten any character vectors inside the list
    flat <- list()
    for (el in centers_cfg) {
      if (is.character(el) && length(el) > 1L) {
        for (tk in el) flat[[length(flat)+1L]] <- tk
      } else {
        flat[[length(flat)+1L]] <- el
      }
    }
    centers_list <- lapply(flat, as_pair)
  } else {
    stop(sprintf("Unsupported VesselCenters type: %s", class(centers_cfg)[1]))
  }
  # basic validation
  for (idx in seq_along(centers_list)) {
    ci <- centers_list[[idx]][1]; cj <- centers_list[[idx]][2]
    if (!is.finite(ci) || !is.finite(cj)) stop(sprintf("Invalid vessel center at index %d (NA/NaN)", idx))
    if (ci < 1L || ci > N || cj < 1L || cj > N) stop(sprintf("Vessel center at index %d is out of grid after clamping: (%d,%d)", idx, ci, cj))
  }
  centers_list
}

# -------------------------------
# Rcpp accelerators (oxygen update + batch G & DivisionTime)
# -------------------------------
Rcpp::cppFunction(code='
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
NumericMatrix update_oxygen_field_cpp(NumericMatrix O2, IntegerMatrix grid,
                                      NumericVector Pv,
                                      double Coxy, double diff_rate,
                                      double supply_rate, double consume,
                                      int iters) {
  int N = O2.nrow();
  int M = O2.ncol();
  int nP = Pv.size();
  for (int it = 0; it < iters; ++it) {
    NumericMatrix nextO2(N, M);
    for (int i = 0; i < N; ++i) {
      for (int j = 0; j < M; ++j) {
        int iu = (i-1<0)?0:i-1;
        int id = (i+1>=N)?N-1:i+1;
        int jl = (j-1<0)?0:j-1;
        int jr = (j+1>=M)?M-1:j+1;
        double up    = O2(iu,j);
        double down  = O2(id,j);
        double left  = O2(i,jl);
        double right = O2(i,jr);
        double center = O2(i,j);
        double lap = (up + down + left + right - 4.0*center);
        // per-cell consumption scaled by P: consume_eff = consume * (P / 2)
        double cons = 0.0;
        int idx = grid(i,j);
        if (idx != NA_INTEGER) {
          int k = idx - 1; // 0-based
          if (k >= 0 && k < nP) {
            double P = Pv[k];
            if (!(R_finite(P)) || P <= 0.0) P = 1e-6;
            cons = consume * (P / 2.0);
          } else {
            cons = consume; // fallback
          }
        }
        double val = center + diff_rate * lap + supply_rate * (Coxy - center) - cons;
        if (val < 0.0) val = 0.0;
        if (val > 1.0) val = 1.0;
        nextO2(i,j) = val;
      }
    }
    O2 = nextO2;
  }
  return O2;
}
')

Rcpp::cppFunction(code='\
#include <Rcpp.h>\
using namespace Rcpp;\
// [[Rcpp::plugins(cpp11)]]\
\
// [[Rcpp::export]]\
NumericMatrix update_oxygen_field_vessels_cpp(NumericMatrix O2, IntegerMatrix grid,\
                                              NumericVector Pv, IntegerMatrix vessel_mask,\
                                              double Coxy, double diff_rate,\
                                              double supply_rate, double consume,\
                                              int iters, int use_soft_supply) {\
  int N = O2.nrow();\
  int M = O2.ncol();\
  int nP = Pv.size();\
  for (int it = 0; it < iters; ++it) {\
    NumericMatrix nextO2(N, M);\
    for (int i = 0; i < N; ++i) {\
      for (int j = 0; j < M; ++j) {\
        int iu = (i-1<0)?0:i-1;\
        int id = (i+1>=N)?N-1:i+1;\
        int jl = (j-1<0)?0:j-1;\
        int jr = (j+1>=M)?M-1:j+1;\
        double up    = O2(iu,j);\
        double down  = O2(id,j);\
        double left  = O2(i,jl);\
        double right = O2(i,jr);\
        double center = O2(i,j);\
        double lap = (up + down + left + right - 4.0*center);\
        // per-cell consumption scaled by P: consume_eff = consume * (P / 2)\
        double cons = 0.0;\
        int idx = grid(i,j);\
        if (idx != NA_INTEGER) {\
          int k = idx - 1;\
          if (k >= 0 && k < nP) {\
            double P = Pv[k];\
            if (!(R_finite(P)) || P <= 0.0) P = 1e-6;\
            cons = consume * (P / 2.0);\
          } else {\
            cons = consume;\
          }\
        }\
        if (vessel_mask(i,j) == 1) {\
          // HARD pin at vessel cells: override diffusion/consumption, keep vessel O2 fixed at Coxy\
          nextO2(i,j) = Coxy;\
        } else {\
          // Non-vessel cells: diffusion from pinned vessels drives gradients;\
          // optionally add soft supply relaxation toward Coxy\
          double val = center + diff_rate * lap - cons;\
          if (use_soft_supply == 1) {\
            val += supply_rate * (Coxy - center);\
          }\
          if (val < 0.0) val = 0.0;\
          if (val > 1.0) val = 1.0;\
          nextO2(i,j) = val;\
        }\
      }\
    }\
    O2 = nextO2;\
  }\
  return O2;\
}\
')

Rcpp::cppFunction(code='
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
NumericMatrix update_oxygen_field_bc_cpp(NumericMatrix O2, IntegerMatrix grid,
                                         NumericVector Pv,
                                         double Coxy, double diff_rate,
                                         double supply_rate, double consume,
                                         int iters, int bc_mode, double bc_value) {
  // bc_mode: 0 = Neumann-like (reflect by clamping), 1 = Dirichlet edges (fixed edge value = bc_value)
  int N = O2.nrow();
  int M = O2.ncol();
  int nP = Pv.size();
  for (int it = 0; it < iters; ++it) {
    NumericMatrix nextO2(N, M);
    for (int i = 0; i < N; ++i) {
      for (int j = 0; j < M; ++j) {
        // Neighbor sampling with boundary condition
        double up    = (i-1 >= 0) ? O2(i-1,j) : (bc_mode==1 ? bc_value : O2(i,j));
        double down  = (i+1 <  N) ? O2(i+1,j) : (bc_mode==1 ? bc_value : O2(i,j));
        double left  = (j-1 >= 0) ? O2(i,j-1) : (bc_mode==1 ? bc_value : O2(i,j));
        double right = (j+1 <  M) ? O2(i,j+1) : (bc_mode==1 ? bc_value : O2(i,j));
        double center = O2(i,j);
        double lap = (up + down + left + right - 4.0*center);
        // per-cell consumption scaled by P: consume_eff = consume * (P / 2)
        double cons = 0.0;
        int idx = grid(i,j);
        if (idx != NA_INTEGER) {
          int k = idx - 1;
          if (k >= 0 && k < nP) {
            double P = Pv[k];
            if (!(R_finite(P)) || P <= 0.0) P = 1e-6;
            cons = consume * (P / 2.0);
          } else {
            cons = consume;
          }
        }
        double val = center + diff_rate * lap + supply_rate * (Coxy - center) - cons;
        if (val < 0.0) val = 0.0;
        if (val > 1.0) val = 1.0;
        // Enforce Dirichlet boundary value at the physical edge cells if requested
        if (bc_mode == 1 && (i==0 || j==0 || i==N-1 || j==M-1)) {
          nextO2(i,j) = bc_value;
        } else {
          nextO2(i,j) = val;
        }
      }
    }
    O2 = nextO2;
  }
  return O2;
}
')

Rcpp::cppFunction(code='
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
Rcpp::List compute_G_and_Div_cpp(Rcpp::IntegerMatrix grid, Rcpp::NumericMatrix O2,
                                 Rcpp::IntegerVector X, Rcpp::IntegerVector Y,
                                 Rcpp::NumericVector Rv, Rcpp::NumericVector Pv,
                                 double beta, int Kcap, int radius) {
  int N = grid.nrow();
  int M = grid.ncol();
  int n = X.size();
  Rcpp::NumericVector G(n);
  Rcpp::NumericVector Div(n);
  for (int idx = 0; idx < n; ++idx) {
    int xi = X[idx] - 1;
    int yj = Y[idx] - 1;
    if (xi < 0 || xi >= N || yj < 0 || yj >= M) {
      G[idx] = 0.0;
      Div[idx] = R_PosInf;
      continue;
    }
    int rmin = std::max(0, xi - radius);
    int rmax = std::min(N - 1, xi + radius);
    int cmin = std::max(0, yj - radius);
    int cmax = std::min(M - 1, yj + radius);
    int n_occ = 0;
    for (int ii = rmin; ii <= rmax; ++ii) {
      for (int jj = cmin; jj <= cmax; ++jj) {
        if (grid(ii,jj) != NA_INTEGER) ++n_occ;
      }
    }
    double Rcell = Rv[idx];
    double th = 1.0 - O2(xi, yj);
    double Pval = Pv[idx];
    if (Pval < 0.0) Pval = 0.0;
    double base = Rcell * (1.0 - th * (beta * Pval));
    if (base < 0.0) base = 0.0;
    double dens = 1.0 - ( (double)n_occ / (double)Kcap );
    if (dens < 0.0) dens = 0.0;
    double g = base * dens;
    G[idx] = g;
    if (g > 0.0) Div[idx] = (std::log(2.0) / g) * 24.0; else Div[idx] = R_PosInf;
  }
  return Rcpp::List::create(Rcpp::_["G"] = G, Rcpp::_["DivisionTime"] = Div);
}
')

# ---- New C++ function: compute_G_Div_quiescence_cpp ----
Rcpp::cppFunction(code='
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

// Fast recompute using integral image of occupancy and quiescence bookkeeping
// X,Y are 1-based indices for ALIVE cells only. Vectors must be same length n.
// Returns: list(G, DivisionTime, Quiescent, QuiescentTime, AliveStatus)
// - Quiescent/QuiescentTime are updated in-place logic
// - AliveStatus: 1 = alive, 0 = quiescence-dead (> thresh)
// Kcap should be (2*radius+1)^2 (e.g., 121 for radius=5)
// q_death_thresh_hours is the death threshold for quiescent time (e.g., 72.0)
// [[Rcpp::export]]
Rcpp::List compute_G_Div_quiescence_cpp(int N, Rcpp::NumericMatrix O2,
                                        Rcpp::IntegerVector X, Rcpp::IntegerVector Y,
                                        Rcpp::NumericVector Rv, Rcpp::NumericVector Pv,
                                        double beta, int Kcap, int radius,
                                        Rcpp::IntegerVector Quiescent, Rcpp::NumericVector QuiescentTime,
                                        double q_death_thresh_hours) {
  int n = X.size();
  if (Y.size()!=n || Rv.size()!=n || Pv.size()!=n || Quiescent.size()!=n || QuiescentTime.size()!=n) {
    Rcpp::stop("Input vector lengths must match");
  }
  // Build occupancy from alive cells only
  Rcpp::IntegerMatrix occ(N,N);
  for (int i=0;i<n;++i){
    int xi = X[i]-1; int yj = Y[i]-1;
    if (xi>=0 && xi<N && yj>=0 && yj<N) occ(xi,yj) = 1;
  }
  // Integral image (prefix sum) of occ, size (N+1)x(N+1)
  Rcpp::NumericMatrix pref(N+1, N+1);
  for (int i=1;i<=N;++i){
    double rowsum = 0.0;
    for (int j=1;j<=N;++j){
      rowsum += (double)occ(i-1,j-1);
      pref(i,j) = pref(i-1,j) + rowsum;
    }
  }
  auto rect_sum = [&](int r1,int c1,int r2,int c2)->int{
    // all 1-based inclusive; clamp to [1,N]
    if (r1<1) r1=1; if (c1<1) c1=1; if (r2>N) r2=N; if (c2>N) c2=N;
    if (r1>r2 || c1>c2) return 0;
    double s = pref(r2,c2) - pref(r1-1,c2) - pref(r2,c1-1) + pref(r1-1,c1-1);
    if (s < 0) s = 0; // numerical safety
    return (int)std::llround(s);
  };

  Rcpp::NumericVector G(n); Rcpp::NumericVector Div(n);
  Rcpp::IntegerVector Qout = clone(Quiescent);
  Rcpp::NumericVector QTout = clone(QuiescentTime);
  Rcpp::IntegerVector Alive(n, 1);

  for (int idx=0; idx<n; ++idx){
    int xi = X[idx]; int yj = Y[idx];
    int rmin = xi - radius; int rmax = xi + radius;
    int cmin = yj - radius; int cmax = yj + radius;
    int n_occ = rect_sum(rmin, cmin, rmax, cmax);

    // growth components
    double Rcell = Rv[idx];
    double th = 1.0 - O2(xi-1, yj-1);
    double Pval = Pv[idx];
    if (!(R_finite(Pval)) || Pval < 0.0) Pval = 0.0;
    // Align growth with compute_local_G and compute_G_and_Div_cpp:
    // base = R * (1 - theta * (beta * P))
    double base = Rcell * (1.0 - th * (beta * Pval));
    if (base < 0.0) base = 0.0;
    double dens = 1.0 - ((double)n_occ / (double)Kcap);
    if (dens < 0.0) dens = 0.0;
    double g = base * dens;

    // quiescence: fully-occupied window => quiescent
    if (n_occ >= Kcap) {
      Qout[idx] = 1;
      g = 0.0; // ensure G=0 under full occupancy
    }

    G[idx] = g;
    Div[idx] = (g > 0.0) ? ((std::log(2.0) / g) * 24.0) : R_PosInf;

    // recovery: if previously quiescent and now g>0, clear flags
    if (Qout[idx] == 1 && g > 0.0) {
      Qout[idx] = 0;
      QTout[idx] = 0.0;
    }
    // accumulate quiescent time for those still quiescent
    if (Qout[idx] == 1) {
      QTout[idx] = QTout[idx] + 1.0; // +1 hour per step
      if (QTout[idx] > q_death_thresh_hours) {
        Alive[idx] = 0; // mark for death
      }
    }
  }
  return Rcpp::List::create(Rcpp::_["G"]=G,
                            Rcpp::_["DivisionTime"]=Div,
                            Rcpp::_["Quiescent"]=Qout,
                            Rcpp::_["QuiescentTime"]=QTout,
                            Rcpp::_["AliveStatus"]=Alive);
}
')

# ==== Additional C++ helpers ====

Rcpp::cppFunction(code='
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

// Return a logical mask marking rows that contain any NA or any chromosome copy < 1
// kt_mat: IntegerMatrix with 22 chromosome columns
// [[Rcpp::export]]
Rcpp::LogicalVector lethal_zero_copy_mask_cpp(Rcpp::IntegerMatrix kt_mat){
  int n = kt_mat.nrow();
  int m = kt_mat.ncol();
  Rcpp::LogicalVector out(n, false);
  for (int i=0; i<n; ++i){
    for (int j=0; j<m; ++j){
      int v = kt_mat(i,j);
      if (v == NA_INTEGER || v < 1){
        out[i] = true;
        break;
      }
    }
  }
  return out;
}
')

Rcpp::cppFunction(code='
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

// Return a logical mask for stochastic baseline death with per-cell probability D[i]
// [[Rcpp::export]]
Rcpp::LogicalVector stochastic_death_mask_cpp(Rcpp::NumericVector D){
  int n = D.size();
  Rcpp::LogicalVector out(n);
  for (int i=0; i<n; ++i){
    double p = D[i];
    if (!R_finite(p) || p <= 0.0) { out[i] = false; continue; }
    if (p >= 1.0) { out[i] = true;  continue; }
    out[i] = (R::runif(0.0, 1.0) < p);
  }
  return out;
}
')

Rcpp::cppFunction(code='
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

// Return a logical mask where P[i] > threshold
// [[Rcpp::export]]
Rcpp::LogicalVector ploidy_over_cap_mask_cpp(Rcpp::NumericVector P, double thr){
  int n = P.size();
  Rcpp::LogicalVector out(n, false);
  for (int i=0; i<n; ++i){
    double v = P[i];
    if (R_finite(v) && v > thr) out[i] = true;
  }
  return out;
}
')
Rcpp::cppFunction(code='
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
IntegerMatrix empty_neighbors_cpp(IntegerMatrix grid, int i, int j) {
  int N = grid.nrow();
  int M = grid.ncol();
  static const int di[8] = {-1,-1,-1, 0,0, 1,1,1};
  static const int dj[8] = {-1, 0, 1,-1,1,-1,0,1};
  std::vector<int> ris, cjs;
  for (int k=0;k<8;++k){
    int ii = i + di[k];
    int jj = j + dj[k];
    if (ii>=1 && ii<=N && jj>=1 && jj<=M){
      if (IntegerVector::is_na(grid(ii-1,jj-1))) {
        ris.push_back(ii);
        cjs.push_back(jj);
      }
    }
  }
  int n = ris.size();
  IntegerMatrix out(n,2);
  for (int t=0;t<n;++t){ out(t,0)=ris[t]; out(t,1)=cjs[t]; }
  return out;
}
')

Rcpp::cppFunction(code='
#include <Rcpp.h>
#include <limits>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
List apply_MS_to_two_cpp(IntegerVector parent_vec, int k_events, int max_copy){
  int n = parent_vec.size();
  IntegerVector a = clone(parent_vec);
  IntegerVector b = clone(parent_vec);
  if (k_events <= 0) return List::create(_["a"]=a, _["b"]=b);
  int upper = (max_copy > 0) ? max_copy : std::numeric_limits<int>::max();
  std::vector<int> valid; valid.reserve(n);
  for(int j=0;j<n;++j){ int v = parent_vec[j]; if (v >= 1 && v <= upper - 1) valid.push_back(j); }
  if (valid.empty()) return List::create(_["a"]=a, _["b"]=b);
  int k = std::min((int)valid.size(), k_events);
  for (int t = (int)valid.size()-1; t>0; --t){ int u = (int)std::floor(R::runif(0.0, t+1)); std::swap(valid[t], valid[u]); }
  for (int t=0; t<k; ++t){
    int j = valid[t];
    if (R::runif(0.0,1.0) < 0.5){
      a[j] = a[j] + 1; if (a[j] > upper) a[j] = upper; b[j] = std::max(0,  b[j]-1);
    } else {
      a[j] = std::max(0,  a[j]-1); b[j] = b[j] + 1; if (b[j] > upper) b[j] = upper;
    }
  }
  return List::create(_["a"]=a, _["b"]=b);
}
')

Rcpp::cppFunction(code='
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
double mean_copy_number_cpp(IntegerVector kt){
  int n = kt.size(); double s = 0.0; for(int i=0;i<n;++i) s += kt[i]; return s / (double)n;
}
')

Rcpp::cppFunction(code='
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
NumericVector mean_copy_number_rows_cpp(IntegerMatrix kt_mat){
  int n = kt_mat.nrow(); int m = kt_mat.ncol(); NumericVector out(n);
  for(int i=0;i<n;++i){ double s=0.0; for(int j=0;j<m;++j) s += kt_mat(i,j); out[i] = s / (double)m; }
  return out;
}
')

Rcpp::cppFunction(code='
#include <Rcpp.h>
#include <limits>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
IntegerVector apply_WGD_to_row_vec_cpp(IntegerVector kt, int max_copy){
  int n = kt.size(); IntegerVector out = clone(kt);
  int upper = (max_copy > 0) ? max_copy : std::numeric_limits<int>::max();
  for(int i=0;i<n;++i){ long v = (long)out[i] * 2L; if (v > upper) v = upper; out[i] = (int)v; }
  return out;
}
')

Rcpp::cppFunction(code='
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
IntegerVector karyotype_from_P_cpp(double P, double tol=1e-3, int max_iter=20000, int max_copy=-1){
  const int chr_n = 22; IntegerVector kt(chr_n); for(int i=0;i<chr_n;++i) kt[i]=2;
  auto mean22 = [&](const IntegerVector& v)->double{ double s=0.0; for(int i=0;i<chr_n;++i) s += v[i]; return s/22.0; };
  double currP = mean22(kt); double err = std::fabs(currP - P); int iter=0;
  while (iter < max_iter && err > tol){
    int dir = (currP < P) ? +1 : -1; int best_j = -1; double best_err = err;
    for(int j=0;j<chr_n;++j){
      int nv = kt[j] + dir;
      if (nv < 0) continue;
      if (max_copy > 0 && nv > max_copy) continue;
      double newP = currP + (nv - kt[j]) / 22.0; double e2 = std::fabs(newP - P);
      if (e2 + 1e-12 < best_err){ best_err = e2; best_j = j; }
    }
    if (best_j < 0) break; kt[best_j] = kt[best_j] + dir; currP = currP + dir / 22.0; err = best_err; ++iter;
  }
  return kt;
}
')

Rcpp::cppFunction(code='
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
IntegerMatrix karyotype_from_P_batch_cpp(NumericVector Pvec, double tol=1e-3, int max_iter=20000, int max_copy=-1){
  const int chr_n = 22; int n = Pvec.size(); IntegerMatrix out(n,22);
  for (int i=0;i<n;++i){
    double P = Pvec[i];
    IntegerVector kt(chr_n); for(int c=0;c<chr_n;++c) kt[c]=2;
    auto mean22 = [&](const IntegerVector& v)->double{ double s=0.0; for(int t=0;t<chr_n;++t) s += v[t]; return s/22.0; };
    double currP = mean22(kt); double err = std::fabs(currP - P); int iter=0;
    while (iter < max_iter && err > tol){
      int dir = (currP < P) ? +1 : -1; int best_j = -1; double best_err = err;
      for(int j=0;j<chr_n;++j){
        int nv = kt[j] + dir;
        if (nv < 0) continue;
        if (max_copy > 0 && nv > max_copy) continue;
        double newP = currP + (nv - kt[j]) / 22.0; double e2 = std::fabs(newP - P);
        if (e2 + 1e-12 < best_err){ best_err = e2; best_j = j; }
      }
      if (best_j < 0) break; kt[best_j] = kt[best_j] + dir; currP = currP + dir / 22.0; err = best_err; ++iter;
    }
    for (int c=0;c<chr_n;++c) out(i,c) = kt[c];
  }
  return out;
}
')

Rcpp::cppFunction(code='
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
IntegerMatrix rebuild_grid_from_cells_cpp(IntegerMatrix grid, IntegerVector X, IntegerVector Y){
  int N = grid.nrow(); int M = grid.ncol();
  for (int i=0;i<N;++i) for (int j=0;j<M;++j) grid(i,j) = NA_INTEGER;
  int n = X.size();
  for (int r=0;r<n;++r){ int xi = X[r]-1; int yj = Y[r]-1; if (xi>=0 && xi<N && yj>=0 && yj<M){ grid(xi,yj) = r+1; } }
  return grid;
}
')

# R wrapper for batch recompute using quiescence/death-aware C++ function
recompute_all_G_and_divtime_cpp <- function(state, radius = 5L, Kcap = 121L) {
  cells <- state$cells; O2 <- state$O2; cfg <- state$cfg
  alive_idx <- which(cells$Status == 1L)
  if (length(alive_idx) == 0L) return(state)
  res <- compute_G_Div_quiescence_cpp(
    as.integer(nrow(state$grid)), O2,
    as.integer(cells$X[alive_idx]), as.integer(cells$Y[alive_idx]),
    as.numeric(cells$R[alive_idx]), as.numeric(cells$P[alive_idx]),
    as.numeric(cfg$beta), as.integer(Kcap), as.integer(radius),
    as.integer(cells$Quiescent[alive_idx]), as.numeric(cells$QuiescentTime[alive_idx]),
    as.numeric(cfg$QuiescentDeathHours)
  )
  cells$G[alive_idx]            <- as.numeric(res$G)
  cells$DivisionTime[alive_idx] <- as.numeric(res$DivisionTime)
  cells$Quiescent[alive_idx]    <- as.integer(res$Quiescent)
  cells$QuiescentTime[alive_idx]<- as.numeric(res$QuiescentTime)
  alive_flag <- as.integer(res$AliveStatus)
  if (any(alive_flag == 0L)) {
    die_idx <- alive_idx[alive_flag == 0L]
    cells$Status[die_idx] <- 0L
    cells$DeathReason[die_idx] <- ifelse(is.na(cells$DeathReason[die_idx]) | cells$DeathReason[die_idx]=="",
                                         "quiescence_timeout", cells$DeathReason[die_idx])
  }
  state$cells <- cells
  state
}

# -------------------------------
# Neighborhood helpers
# -------------------------------
get_window_indices <- function(i, j, N, radius = 5L) {
  rmin <- max(1L, i - radius); rmax <- min(N, i + radius)
  cmin <- max(1L, j - radius); cmax <- min(N, j + radius)
  list(r = rmin:rmax, c = cmin:cmax)
}

moore_neighbors8 <- function(i, j, N) {
  di <- c(-1,-1,-1, 0, 0, 1, 1, 1)
  dj <- c(-1, 0, 1,-1, 1,-1, 0, 1)
  ii <- i + di; jj <- j + dj
  keep <- which(ii >= 1 & ii <= N & jj >= 1 & jj <= N)
  cbind(ii[keep], jj[keep], deparse.level = 0)
}

# -------------------------------
# Growth G in 11x11 window (K=121) with hypoxia–ploidy coupling
# G = max( R * (1 - theta * (beta * P)), 0 ) * max(1 - n_occ / 121, 0)
# where theta = 1 - O2 (hypoxia) in [0,1]; beta multiplies ploidy P.
# -------------------------------
compute_local_G <- function(grid, O2, cells, idx, cfg, Kcap = 121L, radius = 5L) {
  i <- cells$X[idx]; j <- cells$Y[idx]
  N <- nrow(grid)
  rc <- get_window_indices(i, j, N, radius)
  n_occ <- sum(!is.na(grid[rc$r, rc$c]))

  R_cell <- cells$R[idx]
  # theta as hypoxia from oxygen field: theta = 1 - O2
  th <- 1 - oxygen_theta(O2, i, j)
  Pval <- max(cells$P[idx], 0)                 # non-negative ploidy
  base <- R_cell * (1 - th * (cfg$beta * Pval))
  base <- max(base, 0)                         # clamp if over-suppressed
  dens <- max(1 - n_occ / Kcap, 0)             # logistic density term

  base * dens
}

# -------------------------------
# Data containers (32 columns)
# -------------------------------
make_cells_df <- function(n_rows) {
  c_names <- c("X","Y",
               paste0("C",1:22),
               "R","D","P","G","WGDLabel","WGDCount","DivisionTime","Time","Status",
               "Quiescent","QuiescentTime",
               "DeathReason")
  df <- as.data.frame(matrix(NA_real_, nrow = n_rows, ncol = length(c_names)))
  colnames(df) <- c_names
  int_cols <- c("X","Y", paste0("C",1:22), "WGDLabel","WGDCount","Status","Quiescent")
  for (nm in int_cols) df[[nm]] <- as.integer(df[[nm]])
  # Initialize the new columns
  df$Quiescent      <- ifelse(is.na(df$Quiescent), 0L, df$Quiescent)
  df$QuiescentTime  <- ifelse(is.na(df$QuiescentTime), 0.0, df$QuiescentTime)
  df$WGDCount       <- ifelse(is.na(df$WGDCount), 0L, df$WGDCount)
  df$DeathReason    <- NA_character_
  df
}

recompute_P_from_row <- function(cells, idx) {
  kt <- as.integer(cells[idx, paste0("C",1:22)])
  cells$P[idx] <- weighted_P_from_karyotype(kt)
  cells
}

recompute_G_and_divtime <- function(grid, O2, cells, idx, cfg) {
  g <- compute_local_G(grid, O2, cells, idx, cfg = cfg)
  cells$G[idx] <- g
  # Spec: DivisionTime = (log(2) / G) * 24  (use Inf when G <= 0 to avoid division by zero)
  cells$DivisionTime[idx] <- ifelse(g > 0, (log(2) / g) * 24.0, Inf)
  cells
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

init_dead_log <- function(cells_template) {
  df <- cells_template[0, , drop = FALSE]
  if (!("DeathReason" %in% names(df))) df$DeathReason <- character(0)
  df$Step <- integer(0)
  df
}

# -------------------------------
# Initialization
# -------------------------------
init_simulation <- function(cfg, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  N <- as.integer(cfg$N1)
  grid <- matrix(NA_integer_, nrow = N, ncol = N)
  # Determine scalar oxygen supply for this run
  if (is.null(cfg$Coxy_scalar)) {
    cfg$Coxy_scalar <- if (length(cfg$Coxy) <= 1L) as.numeric(cfg$Coxy) else mean(as.numeric(cfg$Coxy), na.rm = TRUE)
  }
  # Build two-vessel mask (diameter = 10 by default, configurable via cfg$VesselDiameter and cfg$VesselCenters)
  diam <- if (!is.null(cfg$VesselDiameter)) as.integer(cfg$VesselDiameter) else 10L
  centers <- NULL
  if (!is.null(cfg$VesselCenters)) {
    centers <- resolve_vessel_centers(N, cfg$VesselCenters, diam)
  }
  vessel_mask <- make_vessel_mask(N, centers = centers, diameter = diam)
  state_vessel_mask <- vessel_mask  # keep for state
  # No boundary source in vessel mode
  cfg$O2UseBoundary <- FALSE
  cfg$O2BoundaryMode <- "neumann"
  cfg$O2BoundaryValue <- 0.0
  # Initialize O2: uniform Coxy everywhere at t=0; vessels will pin at Coxy each update
  O2 <- matrix(cfg$Coxy_scalar, nrow = N, ncol = N)

  K <- load_K_column(cfg$KaryoPath)
  fit <- fit_lognormal_from_K(K)

  n_init <- max(1L, round(N*N * (cfg$m/100)))
  allowed <- which(as.vector(vessel_mask) == 0L)
  if (length(allowed) < n_init) stop("Not enough non-vessel locations to place initial cells. Reduce m or vessel diameter.")
  pos <- sample(allowed, n_init, replace = FALSE)
  X <- as.integer((pos - 1L) %% N + 1L)
  Y <- as.integer((pos - 1L) %/% N + 1L)

  cells <- make_cells_df(n_init)
  cells$X <- X; cells$Y <- Y

  # Initialize P from fitted LogNormal, clamp to configured ploidy cap if present
  pmax_cfg <- if (is.finite(cfg$PloidyMax)) as.numeric(cfg$PloidyMax) else Inf
  cells$P <- sample_P_from_fit(n_init, fit = fit, clamp = c(0, pmax_cfg))

  # Invert P into integer karyotype close to P (weighted C++ batch), pass cap to C++
  max_copy_cpp <- if (is.finite(cfg$PloidyMax)) as.integer(cfg$PloidyMax) else .Machine$integer.max
  kt_mat <- karyotype_from_P_batch_weighted_cpp(as.numeric(cells$P), as.numeric(CHR_W), 1e-3, 20000L, max_copy_cpp)
  # Clamp initial karyotypes so each chromosome copy number ≤ 5 or PloidyMax (initialization constraint)
  kt_init_cap <- 5L
  if (is.finite(cfg$PloidyMax)) kt_init_cap <- as.integer(min(kt_init_cap, cfg$PloidyMax))
  kt_mat <- pmin(kt_mat, kt_init_cap)
  storage.mode(kt_mat) <- "integer"
  cells[, paste0("C",1:22)] <- kt_mat
  # Backfill P with actual integer-mean (weighted C++)
  cells$P <- mean_copy_number_rows_weighted_cpp(kt_mat, as.numeric(CHR_W))

  # Set scalar attributes
  cells$R <- cfg$R
  cells$D <- cfg$Db
  cells$WGDLabel <- 0L
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
        v <- as.integer(cells[i, paste0("C",1:22)])
        if (is.finite(cfg$PloidyMax)) {
          v <- pmin(as.integer(cfg$PloidyMax), v * 2L)
        } else {
          v <- v * 2L
        }
        cells[i, paste0("C",1:22)] <- v
        cells$WGDLabel[i] <- 1L
        cells$P[i] <- mean_copy_number_weighted_cpp(as.integer(v), as.numeric(CHR_W))
        # Count this as one WGD event at t=0
        if (is.null(cells$WGDCount)) cells$WGDCount <- integer(nrow(cells))
        cells$WGDCount[i] <- as.integer(cells$WGDCount[i]) + 1L
      }
    }
  }

  # Place on grid (C++)
  grid <- rebuild_grid_from_cells_cpp(grid, as.integer(cells$X), as.integer(cells$Y))

  # Initial G & DivisionTime via C++ batch (then randomize Time in [0,24])
  res0 <- compute_G_and_Div_cpp(grid, O2,
                                as.integer(cells$X), as.integer(cells$Y),
                                as.numeric(cells$R), as.numeric(cells$P),
                                as.numeric(cfg$beta), as.integer(121L), as.integer(5L))
  cells$G <- as.numeric(res0$G)
  cells$DivisionTime <- as.numeric(res0$DivisionTime)
  cells$Time <- runif(n_init, min = 0, max = 24)

  list(grid = grid, O2 = O2, cells = cells, cfg = cfg, step = 0L, vessel_mask = state_vessel_mask,
       dead_log = init_dead_log(cells))
}

# -------------------------------
# WGD & MS helpers
# -------------------------------


apply_MS_to_two <- function(parent_vec, k_events) {
  # Perform k_events mis-segregations on distinct chromosomes (if available).
  # Each event picks a chromosome in [1,9] copies and assigns +1 to one child and -1 to the other.
  childA <- parent_vec
  childB <- parent_vec
  k_events <- as.integer(k_events)
  if (!is.finite(k_events) || k_events <= 0L) return(list(a = childA, b = childB))

  valid <- which(parent_vec >= 1L & parent_vec <= 9L)
  if (length(valid) == 0L) return(list(a = childA, b = childB))

  k <- min(k_events, length(valid))
  js <- sample(valid, k, replace = FALSE)
  for (j in js) {
    if (runif(1) < 0.5) {
      # Child A gains, child B loses
      childA[j] <- min(10L, childA[j] + 1L)
      childB[j] <- max(0L, childB[j] - 1L)
    } else {
      # Child B gains, child A loses
      childA[j] <- max(0L, childA[j] - 1L)
      childB[j] <- min(10L, childB[j] + 1L)
    }
  }
  list(a = childA, b = childB)
}


# Effective MSR with environment:
# spec mentions alpha in [-1,1] and also "(1 + theta)"
msr_effective <- function(MSR, theta, alpha) {
  # Effective per-chromosome mis-segregation rate under hypoxia coupling
  # Clamp to [0,1] and handle non-finite inputs defensively
  p <- MSR * (1 + theta * alpha)
  if (!is.finite(p)) return(0)
  max(0, min(1, p))
}

# -------------------------------
# Division (strict per spec, plus Dv-branch when no space)
# -------------------------------
attempt_division <- function(state, irow) {
  grid <- state$grid
  O2   <- state$O2
  cells <- state$cells
  cfg <- state$cfg
  # (Optional safety) Ensure WGDCount exists before use
  if (is.null(cells$WGDCount)) {
    cells$WGDCount <- integer(nrow(cells))
  }

  if (cells$Status[irow] != 1L) return(state)

  # Quiescence check: fully occupied 11x11 window -> set G=0 and Quiescent=1
  {
    N <- nrow(grid)
    i <- cells$X[irow]; j <- cells$Y[irow]
    rc <- get_window_indices(i, j, N, radius = 5L)
    n_occ <- sum(!is.na(grid[rc$r, rc$c]))
    if (n_occ >= 121L) {
      cells$G[irow] <- 0.0
      if (is.null(cells$Quiescent)) cells$Quiescent <- integer(nrow(cells))
      if (is.null(cells$QuiescentTime)) cells$QuiescentTime <- numeric(nrow(cells))
      cells$Quiescent[irow] <- 1L
      # do not attempt division this step
      state$cells <- cells
      return(state)
    }
  }

  # Divide when G>0 and Time <= DivisionTime (one attempt per time step)
  if (!(cells$G[irow] > 0 && cells$Time[irow] <= cells$DivisionTime[irow])) return(state)

  N <- nrow(grid)
  i <- cells$X[irow]; j <- cells$Y[irow]
  empty <- empty_neighbors_cpp(grid, as.integer(i), as.integer(j))

  # Hard forbid using vessel cells as targets: filter them out unconditionally
  if (!is.null(state$vessel_mask) && nrow(empty) > 0L) {
    # keep only neighbors where vessel_mask == 0 (paired indices)
    idx <- cbind(empty[,1], empty[,2])  # paired (row,col) indices
    keep <- state$vessel_mask[idx] == 0L
    empty <- empty[which(keep), , drop = FALSE]
  }

  if (nrow(empty) == 0L) {
    # No neighbor space: try division-death; else enter quiescence
    if (runif(1) < cfg$Dv) {
      cells$Status[irow] <- 0L
      cells$DeathReason[irow] <- ifelse(is.na(cells$DeathReason[irow]) || cells$DeathReason[irow]=="",
                                        "division_death_Dv", cells$DeathReason[irow])
    } else {
      if (is.null(cells$Quiescent)) cells$Quiescent <- integer(nrow(cells))
      if (is.null(cells$QuiescentTime)) cells$QuiescentTime <- numeric(nrow(cells))
      cells$Quiescent[irow] <- 1L
    }
    state$cells <- cells
    return(state)
  }
  choice <- if (nrow(empty) == 1L) 1L else sample(seq_len(nrow(empty)), 1L)
  xi <- empty[choice, 1]; yj <- empty[choice, 2]

  # WGD branch (coordinate unchanged) decided by WGDr; no per-cell cap
  if (runif(1) < cfg$WGDr) {
    v <- as.integer(cells[irow, paste0("C",1:22)])
    v <- apply_WGD_to_row_vec_cpp(v, if (is.finite(cfg$PloidyMax)) as.integer(cfg$PloidyMax) else .Machine$integer.max)
    cells[irow, paste0("C",1:22)] <- v
    cells$WGDLabel[irow] <- 1L
    cells$P[irow] <- mean_copy_number_weighted_cpp(v, as.numeric(CHR_W))
    cells$Time[irow] <- 0.0
    # Increment per-cell WGD count
    cells$WGDCount[irow] <- as.integer(cells$WGDCount[irow]) + 1L
    # Defer G/DivisionTime recomputation to next step-wide refresh
    cells$G[irow] <- NA_real_
    cells$DivisionTime[irow] <- Inf
    state$cells <- cells
    return(state)
  }

  # --- Mitosis branch ---
  temp <- cells[irow, , drop = FALSE]
  temp$DeathReason <- NA_character_
  temp$X <- xi; temp$Y <- yj
  # Reset new columns for newborn
  temp$Quiescent <- 0L
  temp$QuiescentTime <- 0.0
  if (!is.null(temp$WGDCount)) temp$WGDCount <- as.integer(cells$WGDCount[irow])

  # Poisson-distributed number of mis-segregating chromosomes per division
  # Use environment-adjusted effective MSR under hypoxia
  Pval <- max(cells$P[irow], 0)
  val_o2 <- oxygen_theta(O2, i, j)
  if (!is.finite(val_o2)) val_o2 <- 0
  theta_local <- 1 - val_o2  # hypoxia level at the cell (in [0,1])
  eff_msr <- msr_effective(cfg$MSR, theta_local, cfg$alpha)
  if (!is.finite(eff_msr) || eff_msr < 0) eff_msr <- 0
  lambda_ms <- eff_msr * 22.0 * Pval
  if (!is.finite(lambda_ms) || lambda_ms < 0) lambda_ms <- 0
  k_ms <- rpois(1, lambda_ms)

  parent_kt <- as.integer(cells[irow, paste0("C",1:22)])
  two <- apply_MS_to_two_cpp(
    parent_kt,
    as.integer(k_ms),
    if (is.finite(cfg$PloidyMax)) as.integer(cfg$PloidyMax) else .Machine$integer.max
  )

  a_vec <- as.integer(two$a)
  b_vec <- as.integer(two$b)

  # Lethality: any chromosome < 1 (i.e., 0) or NA => lethal
  lethal_a <- any(!is.finite(a_vec)) || any(a_vec < 1L)
  lethal_b <- any(!is.finite(b_vec)) || any(b_vec < 1L)

  # Case 1: both lethal -> kill parent, do not spawn daughter
  if (lethal_a && lethal_b) {
    cells$Status[irow] <- 0L
    cells$DeathReason[irow] <- ifelse(is.na(cells$DeathReason[irow]) || cells$DeathReason[irow]=="",
                                      "division_death_Dv", cells$DeathReason[irow])
    state$cells <- cells
    state$grid  <- grid
    return(state)
  }

  # Case 2: parent lethal, daughter viable -> parent dies; only place daughter
  if (lethal_a && !lethal_b) {
    temp[, paste0("C",1:22)] <- b_vec
    temp$P    <- mean_copy_number_weighted_cpp(as.integer(b_vec), as.numeric(CHR_W))
    temp$Time <- 0.0

    new_row <- nrow(cells) + 1L
    cells <- rbind(cells, temp)
    grid[xi, yj] <- new_row

    cells$Status[irow] <- 0L
    cells$DeathReason[irow] <- ifelse(is.na(cells$DeathReason[irow]) || cells$DeathReason[irow]=="",
                                      "division_death_Dv", cells$DeathReason[irow])
    # Defer newborn G/DivisionTime recomputation to end-of-step
    cells$G[new_row]            <- NA_real_
    cells$DivisionTime[new_row] <- Inf

    state$grid  <- grid
    state$cells <- cells
    return(state)
  }

  # Case 3: parent viable, daughter lethal -> keep parent only; do not place daughter
  if (!lethal_a && lethal_b) {
    cells[irow, paste0("C",1:22)] <- a_vec
    cells$P[irow]    <- mean_copy_number_weighted_cpp(as.integer(a_vec), as.numeric(CHR_W))
    cells$Time[irow] <- 0.0
    # Defer parent's G/DivisionTime recomputation to end-of-step
    cells$G[irow]            <- NA_real_
    cells$DivisionTime[irow] <- Inf

    state$cells <- cells
    return(state)
  }

  # Case 4: both viable -> proceed with normal mitosis
  cells[irow, paste0("C",1:22)] <- a_vec
  temp[,     paste0("C",1:22)] <- b_vec

  # Reset times
  cells$Time[irow] <- 0.0
  temp$Time        <- 0.0

  # Recompute P for both (weighted C++)
  cells$P[irow] <- mean_copy_number_weighted_cpp(as.integer(a_vec), as.numeric(CHR_W))
  temp$P        <- mean_copy_number_weighted_cpp(as.integer(b_vec), as.numeric(CHR_W))

  # Insert daughter and update grid
  new_row <- nrow(cells) + 1L
  cells <- rbind(cells, temp)
  grid[xi, yj] <- new_row

  # Defer G/DivisionTime recomputation for both mother and daughter
  cells$G[c(irow, new_row)]            <- NA_real_
  cells$DivisionTime[c(irow, new_row)] <- Inf

  state$grid  <- grid
  state$cells <- cells
  return(state)
}

# -------------------------------
# One simulation step = 1 hour
# -------------------------------
simulate_step <- function(state) {
  grid <- state$grid; cells <- state$cells; cfg <- state$cfg; O2 <- state$O2

  # Safety: ensure vessel cells are never occupied
  if (!is.null(state$vessel_mask)) {
    # any(!is.na(grid) & state$vessel_mask == 1) should always be FALSE
    bad <- grid[state$vessel_mask == 1L]
    bad <- bad[!is.na(bad)]
    if (length(bad) > 0L) {
      # Mark those cells dead and clear grid at vessel locations
      cells$Status[bad] <- 0L
      cells$DeathReason[bad] <- ifelse(is.na(cells$DeathReason[bad]) | cells$DeathReason[bad]=="",
                                       "vessel_mask_violation", cells$DeathReason[bad])
      grid[state$vessel_mask == 1L] <- NA_integer_
    }
  }

  state$grid <- grid
  state$cells <- cells

  # 0) Update oxygen field (C++ accelerated)
  Coxy_run <- if (!is.null(cfg$Coxy_scalar)) cfg$Coxy_scalar else if (length(cfg$Coxy) <= 1L) as.numeric(cfg$Coxy) else mean(as.numeric(cfg$Coxy), na.rm = TRUE)
  # Tunable PDE coefficients from cfg
  dr <- as.numeric(if (!is.null(cfg$O2DiffRate)) cfg$O2DiffRate else 0.2)
  sr <- as.numeric(if (!is.null(cfg$O2SupplyRate)) cfg$O2SupplyRate else 0.05)
  cs <- as.numeric(if (!is.null(cfg$O2Consume)) cfg$O2Consume else 0.02)
  if (!is.null(state$vessel_mask)) {
    # Vessel-driven oxygen supply: source only at vessel cells
    iters <- if (!is.null(cfg$O2JacobiIters)) as.integer(cfg$O2JacobiIters) else 1L
    if (!is.finite(iters) || iters < 1L) iters <- 1L
    use_soft <- if (!is.null(cfg$O2VesselSoftSupply) && isTRUE(cfg$O2VesselSoftSupply)) 1L else 0L
    O2 <- update_oxygen_field_vessels_cpp(O2, grid, as.numeric(cells$P),
                                          as.matrix(state$vessel_mask),
                                          Coxy_run, dr, sr, cs, iters, use_soft)
  } else if (isTRUE(cfg$O2UseBoundary)) {
    # Boundary-driven inward diffusion (legacy)
    bc_mode <- if (!is.null(cfg$O2BoundaryMode) && tolower(cfg$O2BoundaryMode) %in% c("dirichlet", "dirichlet_edges")) 1L else 0L
    bc_val  <- Coxy_run
    cfg$O2BoundaryValue <- bc_val
    state$cfg <- cfg
    O2 <- update_oxygen_field_bc_cpp(O2, grid, as.numeric(cells$P),
                                     Coxy_run, dr, sr, cs, 1L, bc_mode, bc_val)
  } else {
    # Legacy isotropic diffusion with Neumann-like boundary and uniform supply
    O2 <- update_oxygen_field_cpp(O2, grid, as.numeric(cells$P),
                                  Coxy_run, dr, sr, cs, 1L)
  }
  state$O2 <- O2

  cells <- state$cells
  alive_idx <- which(cells$Status == 1L)
  if (length(alive_idx) == 0L) {
    state$step <- state$step + 1L
    return(state)
  }

  # 1) Time elapses
  cells$Time[alive_idx] <- cells$Time[alive_idx] + 1.0
  state$cells <- cells

  # 2) Attempt division (snapshot indices so newborns won't re-run this step)
  snapshot <- alive_idx
  for (irow in snapshot) {
    if (irow <= nrow(state$cells) && state$cells$Status[irow] == 1L) {
      state <- attempt_division(state, irow)
    }
  }

  # 3) Death process: (a) G<=0 cull; (b) Time-threshold death; (c) baseline stochastic death Db
  cells <- state$cells
  alive_idx <- which(cells$Status == 1L)
  # (new) Death A: kill if any chromosome copy number < 1 (0 or NA treated as lethal)
  if (length(alive_idx) > 0L) {
    kt_cols <- paste0("C", 1:22)
    kt_mat <- as.matrix(cells[alive_idx, kt_cols, drop = FALSE])
    storage.mode(kt_mat) <- "integer"
    lethal_mask <- lethal_zero_copy_mask_cpp(kt_mat)
    if (any(lethal_mask)) {
      zero_chr <- alive_idx[which(lethal_mask)]
      cells$Status[zero_chr] <- 0L
      cells$DeathReason[zero_chr] <- ifelse(is.na(cells$DeathReason[zero_chr]) | cells$DeathReason[zero_chr]=="",
                                            "lethal_zero_copy", cells$DeathReason[zero_chr])
    }
  }

  # (new) Death B: kill if ploidy exceeds threshold (use cfg$PloidyMax; if Inf, treat as 10)
  alive_idx <- which(cells$Status == 1L)
  if (length(alive_idx) > 0L) {
    thr <- if (is.finite(state$cfg$PloidyMax)) as.numeric(state$cfg$PloidyMax) else 10
    over_mask <- ploidy_over_cap_mask_cpp(cells$P[alive_idx], thr)
    if (any(over_mask)) {
      over_ploidy <- alive_idx[which(over_mask)]
      cells$Status[over_ploidy] <- 0L
      cells$DeathReason[over_ploidy] <- ifelse(is.na(cells$DeathReason[over_ploidy]) | cells$DeathReason[over_ploidy]=="",
                                               "ploidy_over_cap", cells$DeathReason[over_ploidy])
    }
  }
  
  # Refresh alive index after the two new death rules
  alive_idx <- which(cells$Status == 1L)
  if (length(alive_idx) > 0L) {
    # (a) G-based death: if a cell's growth intensity G <= 0, mark it dead
    gvals <- cells$G[alive_idx]
    nonprolif <- alive_idx[is.finite(gvals) & gvals <= 0]
    if (length(nonprolif) > 0L) {
      cells$Status[nonprolif] <- 0L
      cells$DeathReason[nonprolif] <- ifelse(is.na(cells$DeathReason[nonprolif]) | cells$DeathReason[nonprolif]=="",
                                             "nonproliferative_G_le_0", cells$DeathReason[nonprolif])
    }

    # refresh alive index after G-based culling
    alive_idx <- which(cells$Status == 1L)

    # (b) Time-based death: if a cell's Time since last division exceeds cfg$Dtime, mark it dead
    if (length(alive_idx) > 0L && !is.null(cfg$Dtime) && is.finite(cfg$Dtime)) {
      over_age <- alive_idx[cells$Time[alive_idx] > cfg$Dtime]
      if (length(over_age) > 0L) {
        cells$Status[over_age] <- 0L
        cells$DeathReason[over_age] <- ifelse(is.na(cells$DeathReason[over_age]) | cells$DeathReason[over_age]=="",
                                              "over_Dtime", cells$DeathReason[over_age])
      }
    }

    # refresh alive index after time-threshold culling
    alive_idx <- which(cells$Status == 1L)

    # (c) Baseline stochastic death with per-cell D (Db baseline)
    if (length(alive_idx) > 0L) {
      kill_mask <- stochastic_death_mask_cpp(cells$D[alive_idx])
      if (any(kill_mask)) {
        deaths <- alive_idx[which(kill_mask)]
        cells$Status[deaths] <- 0L
        cells$DeathReason[deaths] <- ifelse(is.na(cells$DeathReason[deaths]) | cells$DeathReason[deaths]=="",
                                            "baseline_death_Db", cells$DeathReason[deaths])
      }
    }
  }

  # 4) Defer grid rebuild to end-of-step (after quiescence bookkeeping)
  #    We intentionally skip rebuilding here to avoid double work.

  # === End-of-step recompute & quiescence bookkeeping ===
  # Recompute G and DivisionTime once per step (after all updates)
  state$cells <- cells
  state$grid  <- grid
  state <- recompute_all_G_and_divtime_cpp(state)
  cells <- state$cells; grid <- state$grid

  # If a quiescent cell now has G > 0, clear its quiescent flags
  if (!is.null(cells$Quiescent)) {
    alive_idx <- which(cells$Status == 1L)
    if (length(alive_idx) > 0L) {
      recovering <- alive_idx[ (cells$Quiescent[alive_idx] == 1L) & is.finite(cells$G[alive_idx]) & (cells$G[alive_idx] > 0) ]
      if (length(recovering) > 0L) {
        cells$Quiescent[recovering] <- 0L
        cells$QuiescentTime[recovering] <- 0.0
      }
      # Accumulate quiescent time (+1 hour per step) for those still quiescent
      still_q <- alive_idx[ cells$Quiescent[alive_idx] == 1L ]
      if (length(still_q) > 0L) cells$QuiescentTime[still_q] <- cells$QuiescentTime[still_q] + 1.0
      # Death rule: quiescent time > cfg$QuiescentDeathHours -> cell dies
      to_die_q <- alive_idx[ cells$QuiescentTime[alive_idx] > as.numeric(cfg$QuiescentDeathHours) ]
      if (length(to_die_q) > 0L) {
        cells$Status[to_die_q] <- 0L
      }
    }
  }

  # If any cells died (any reason), log them and rebuild grid once
  if (any(cells$Status == 0L)) {
    dead_idx <- which(cells$Status == 0L)
    # Fill missing death reasons
    miss <- is.na(cells$DeathReason[dead_idx]) | cells$DeathReason[dead_idx] == ""
    if (any(miss)) cells$DeathReason[dead_idx[miss]] <- "unspecified"
    # Record step and append to state$dead_log
    step_tag <- state$step + 1L
    tolog <- cells[dead_idx, , drop = FALSE]
    tolog$Step <- step_tag
    if (is.null(state$dead_log)) {
      state$dead_log <- init_dead_log(cells)
    }
    state$dead_log <- rbind(state$dead_log, tolog)
    # Remove from the pool of living cells and rebuild grid
    cells <- cells[cells$Status == 1L, , drop = FALSE]
    grid  <- rebuild_grid_from_cells_cpp(grid, as.integer(cells$X), as.integer(cells$Y))
  }

  state$cells <- cells
  state$grid <- grid
  state$step <- state$step + 1L
  state
}

# -------------------------------
# Driver
# -------------------------------
run_simulation <- function(state, steps = 240L, verbose_every = 24L, out_dir = "data") {
  for (t in seq_len(steps)) {
    state <- simulate_step(state)
    save_cells_snapshot(state, t, dir = out_dir)
    if (!is.null(verbose_every) && verbose_every > 0L && (t %% verbose_every == 0L)) {
      alive <- sum(state$cells$Status == 1L)
      meanG <- if (alive > 0L) mean(state$cells$G[state$cells$Status==1L], na.rm=TRUE) else NA_real_
      cat(sprintf("[t=%4dh] alive=%d, meanG=%s\n",
                  state$step, alive,
                  ifelse(is.na(meanG),"NA", sprintf("%.4f", meanG))))
    }
  }
  state
}












# -------------------------------
#  usage
# -------------------------------
cfg <- read_config('/Users/4482173/Documents/IMO_workshop13/config.yaml')

# Parallel runs over mapped Coxy values with dynamic resource allocation:
# - Priority breadth across Coxy (ri layered)
# - Spare cores are used for additional reps (ri>1)
if (!is.numeric(cfg$Coxy)) stop("cfg$Coxy must be numeric after read_config")
coxy_vec_full <- as.numeric(cfg$Coxy)
# Keep 0 inside cfg$Coxy for scaling/mapping; skip zero in actual simulation runs
coxy_vec <- coxy_vec_full[is.finite(coxy_vec_full) & coxy_vec_full > 0]
if (length(coxy_vec) == 0L) {
  stop("After excluding Coxy==0, no positive Coxy values remain to simulate. Keep 0 in cfg$Coxy for scaling, but include at least one > 0 value.")
}
n_runs <- length(coxy_vec)
S <- as.integer(cfg$SimsPerCoxy)

# Build layered job list: ri = 1..S, for each ci = 1..n_runs
jobs <- vector("list", n_runs * S)
idx <- 1L
for (ri in seq_len(S)) {
  for (ci in seq_len(n_runs)) {
    jobs[[idx]] <- list(ci = ci, ri = ri)
    idx <- idx + 1L
  }
}

run_job <- function(job) {
  tryCatch({
    ci <- job$ci; ri <- job$ri
    cfg_run <- cfg
    cfg_run$Coxy_scalar <- coxy_vec[ci]
    coxy_dir <- sprintf("Coxy_%0.6f", cfg_run$Coxy_scalar)
    out_dir <- file.path("data", coxy_dir, sprintf("Rep_%03d", ri))
    # No reproducibility required: do not set seeds; allow default RNG
    state <- init_simulation(cfg_run, seed = NULL)
    state <- run_simulation(state, steps = 720, out_dir = out_dir)
    list(status = "ok", ci = ci, ri = ri, coxy = cfg_run$Coxy_scalar,
         out_dir = out_dir,
         alive = sum(state$cells$Status == 1L))
  }, error = function(e) {
    list(status = "error", ci = job$ci, ri = job$ri,
         coxy = if (job$ci >= 1 && job$ci <= length(coxy_vec)) coxy_vec[job$ci] else NA_real_,
         out_dir = NA_character_, alive = NA_integer_,
         error = conditionMessage(e))
  })
}

if (.Platform$OS.type == "windows") {
  # Fallback: sequential over jobs on Windows
  res_list <- lapply(jobs, run_job)
} else {
  ncores <- max(1L, parallel::detectCores(logical = TRUE) - 1L)
  # Dynamic scheduling lets idle workers pull next job immediately
  res_list <- parallel::mclapply(jobs, run_job, mc.cores = ncores, mc.preschedule = FALSE)
}

invisible(lapply(seq_along(res_list), function(i) {
  res <- res_list[[i]]
  if (is.list(res) && !is.null(res$coxy)) {
    msg_err <- if (!is.null(res$status) && res$status == "error") sprintf(" [ERROR: %s]", res$error) else ""
    cat(sprintf("\n[RUN Coxy=%0.6f rep=%03d] finished. Alive cells: %s%s\n",
                as.numeric(res$coxy), as.integer(res$ri),
                ifelse(is.na(res$alive), "NA", as.character(res$alive)),
                msg_err))
  } else {
    cat(sprintf("\n[RUN #%d] unexpected result type: %s\n", i, typeof(res)))
  }
}))

