# beta_estimation.R
# Input: a CSV with columns at least including:
# passage_id, Ploidy, g, correctedCount, sublabel_value, K, oxygen_pct
# Output: console summaries + CSVs under <input_dir>/beta_out
# Method: control-based r(P); oxygen stratification; robust filters; beta from G = r*(1 - theta*beta*P)*dens with dens ≈ (1 - N/K)

suppressPackageStartupMessages(library(dplyr))

# ---- I/O ----
infile <- '/Users/4482173/Documents/IMO_workshop13/meta_data.csv'  # <-- adjust if needed
stopifnot(file.exists(infile))

df <- utils::read.csv(infile, stringsAsFactors = FALSE)

# Ensure required columns exist
req_cols <- c('passage_id','Ploidy','g','correctedCount','sublabel_value','K','oxygen_pct')
missing_cols <- setdiff(req_cols, names(df))
if (length(missing_cols) > 0) stop('Missing required columns: ', paste(missing_cols, collapse=', '))

# Coerce numerics conservatively
num_cols <- c('Ploidy','g','correctedCount','K')
for (cc in num_cols) df[[cc]] <- suppressWarnings(as.numeric(df[[cc]]))

# ---- Oxygen scaling: map oxygen_pct to fraction in [0,1] ----
scale_O2_frac <- function(x) {
  # Accept numbers or strings like "80%"; clamp to [0,1]
  x_chr <- as.character(x)
  x_chr <- trimws(x_chr)
  x_chr <- gsub('%','', x_chr, fixed = TRUE)
  x_num <- suppressWarnings(as.numeric(x_chr))
  prop_gt1 <- mean(x_num > 1, na.rm = TRUE)
  if (!is.finite(prop_gt1)) prop_gt1 <- 0
  if (prop_gt1 > 0.3) x_num <- x_num / 100
  x_num[x_num < 0] <- NA_real_
  x_num[x_num > 1] <- 1
  x_num
}

# Add oxygen fraction and hypoxia level
# O2_frac \in [0,1], theta = 1 - O2_frac

O2_frac <- scale_O2_frac(df$oxygen_pct)
df <- df %>% mutate(O2_frac = O2_frac,
                    theta = 1 - O2_frac)

# ---- Step 1: Estimate r(P) from CONTROL only ----
ctrl <- df %>%
  mutate(sublabel_value = tolower(sublabel_value)) %>%
  filter(sublabel_value == 'control') %>%
  mutate(logi = 1 - correctedCount / K)

r_rep <- ctrl %>%
  filter(is.finite(g), is.finite(K), is.finite(correctedCount),
         K > 0, correctedCount > 0,
         is.finite(logi), logi > 0.05, logi < 1.2) %>%
  group_by(Ploidy, passage_id) %>%
  summarise(r_rep = stats::median(g / logi, na.rm = TRUE), .groups = 'drop')

r_tab <- r_rep %>%
  group_by(Ploidy) %>%
  summarise(r_hat = stats::median(r_rep, na.rm = TRUE), .groups = 'drop') %>%
  arrange(Ploidy)

cat('\nPloidy-specific r_hat (control-based):\n')
print(r_tab)

# Optional: replicate-level linear trend r_rep ~ Ploidy (for reporting only)
if (nrow(r_rep) >= 2 && length(unique(r_rep$Ploidy)) >= 2) {
  r_lm <- try(stats::lm(r_rep ~ Ploidy, data = r_rep), silent = TRUE)
  if (!inherits(r_lm, 'try-error')) {
    cat('\nEstimated linear relationship between Ploidy and r (replicate-level):\n')
    print(summary(r_lm))
  }
}

# ---- Step 2: Row-wise beta using CONTROL + STARVATION ----
starv <- df %>%
  mutate(sublabel_value = tolower(sublabel_value)) %>%
  filter(sublabel_value %in% c('starvation','control')) %>%
  mutate(logi = 1 - correctedCount / K) %>%
  # density proxy for ABM crowding term: dens := max(0, min(1, 1 - N/K))
  mutate(dens = pmin(pmax(logi, 0), 1)) %>%
  left_join(r_tab, by = 'Ploidy') %>%
  # new denominator from G = r*(1 - theta*beta*P)*dens
  mutate(denom = 1 - g / (r_hat * dens))

starv <- starv %>%
  mutate(beta_hat = ifelse(
    # guards: finite inputs
    is.finite(g) & is.finite(r_hat) & is.finite(dens) & is.finite(theta) &
    is.finite(Ploidy) & is.finite(K) & is.finite(correctedCount),
    {
      # validity: positive capacities, sensible density proxy and theta
      ok <- (K > 0) & (correctedCount > 0) &
            (dens > 1e-6) & (dens <= 1.0 + 1e-8) &
            (theta >= 0) & (theta <= 1) &
            (r_hat > 0) & (Ploidy > 0)
      val <- (1 - g / (r_hat * dens)) / (theta * Ploidy)
      # only keep when denom>0 and theta>0 to avoid negative/inf beta
      ifelse(ok & (denom > 0) & (theta > 0), val, NA_real_)
    },
    NA_real_
  ))

# ---- Step 2.5: Oxygen stratification (bin by identical O2 up to 1e-3) ----
starv <- starv %>% mutate(O2_level = round(O2_frac, 3))

beta_by_O2 <- starv %>%
  filter(is.finite(beta_hat), is.finite(O2_level)) %>%
  group_by(O2_level) %>%
  summarise(beta_median = stats::median(beta_hat, na.rm = TRUE),
            beta_mean   = mean(beta_hat, na.rm = TRUE),
            beta_sd     = stats::sd(beta_hat,   na.rm = TRUE),
            n           = dplyr::n(),
            .groups = 'drop') %>%
  arrange(O2_level)

cat('\nEstimated beta per oxygen level (O2_level rounded to 1e-3):\n')
print(beta_by_O2)

# ---- Step 3: Aggregate across oxygen levels ----
beta_agg <- beta_by_O2 %>%
  summarise(beta_median_over_levels = stats::median(beta_median, na.rm = TRUE),
            beta_mean_over_levels   = mean(beta_median,   na.rm = TRUE),
            beta_weighted_mean      = ifelse(sum(n) > 0, sum(beta_median * n) / sum(n), NA_real_),
            levels                  = dplyr::n())

cat('\nAggregated beta across oxygen levels (oxygen influence removed):\n')
print(beta_agg)

cat('\nNOTE: beta computed from G = r*(1 - theta*beta*P)*dens with dens ≈ clamp(1 - N/K, 0, 1).\n')

# ---- Write outputs ----
outdir <- file.path(dirname(infile), 'beta_out')
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

utils::write.csv(r_tab,       file.path(outdir, 'r_tab.csv'),       row.names = FALSE)
utils::write.csv(beta_by_O2,  file.path(outdir, 'beta_by_O2.csv'),  row.names = FALSE)
utils::write.csv(beta_agg,    file.path(outdir, 'beta_agg.csv'),    row.names = FALSE)
utils::write.csv(starv %>%
  dplyr::select(passage_id, Ploidy, g, correctedCount, K, sublabel_value,
                O2_frac, theta, logi, dens, r_hat, denom, beta_hat, O2_level),
  file.path(outdir, 'beta_rowwise.csv'), row.names = FALSE)

cat('\nOutputs written to: ', outdir, '\n')
