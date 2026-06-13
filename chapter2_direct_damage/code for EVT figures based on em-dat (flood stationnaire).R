# ================================================
# EVT on EM-DAT (Annual Max) — Hardened Version
# (with correct L-moment fallback using lmomco::lmoms/pargev)
# ================================================

options(stringsAsFactors = FALSE)

# ---- Packages (install if missing) ----
needs <- c("tidyverse","lubridate","extRemes","lmomco","ggplot2","scales","readr")
to_get <- needs[!needs %in% rownames(installed.packages())]
if (length(to_get)) install.packages(to_get, dependencies = TRUE)
invisible(lapply(needs, library, character.only = TRUE))

# ---- Quick sanity on data presence ----
if (!exists("em_dat")) stop("Object `em_dat` not found in your R session.")

# ---- Column mapping (matches your schema) ----
COL_ID   <- "DisNo."
COL_TYPE <- "Disaster Type"
COL_YR   <- "Start Year"
COL_MO   <- "Start Month"
COL_DY   <- "Start Day"
LOSS_ADJ <- "Total Damage, Adjusted ('000 US$)"
LOSS_RAW <- "Total Damage ('000 US$)"
COL_LOSS <- if (LOSS_ADJ %in% names(em_dat)) LOSS_ADJ else LOSS_RAW

# ---- User controls ----
HAZARD_KEEP <- c("Flood")
T_vec       <- c(2,5,10,20,50,100)

# ---- Helper: fail early with context ----
need_cols <- c(COL_ID, COL_TYPE, COL_YR, COL_LOSS)
miss <- setdiff(need_cols, names(em_dat))
if (length(miss)) stop("Missing expected columns in `em_dat`: ", paste(miss, collapse=", "))

# ---- 1) Tidy events ----
raw <- em_dat
yr <- suppressWarnings(as.integer(raw[[COL_YR]]))
mo <- if (COL_MO %in% names(raw)) suppressWarnings(as.integer(raw[[COL_MO]])) else 1L
dy <- if (COL_DY %in% names(raw)) suppressWarnings(as.integer(raw[[COL_DY]])) else 1L

date_vec <- lubridate::make_date(
  year  = yr,
  month = ifelse(is.na(mo) | mo < 1, 1L, mo),
  day   = ifelse(is.na(dy) | dy < 1, 1L, dy)
)

loss_k   <- readr::parse_number(as.character(raw[[COL_LOSS]]))  # '000 US$
loss_usd <- 1000 * loss_k

df_event0 <- tibble::tibble(
  event_id = raw[[COL_ID]],
  type     = trimws(as.character(raw[[COL_TYPE]])),
  year     = yr,
  date     = date_vec,
  loss_usd = as.numeric(loss_usd)
) |>
  dplyr::filter(!is.na(year), year > 0, !is.na(loss_usd), loss_usd > 0) |>
  dplyr::filter(type %in% HAZARD_KEEP)

if (nrow(df_event0) == 0) {
  cat("\nTop 10 types with losses > 0 (to help choose HAZARD_KEEP):\n")
  print(
    tibble(type = trimws(as.character(raw[[COL_TYPE]])),
           loss = as.numeric(1000*readr::parse_number(as.character(raw[[COL_LOSS]])))) |>
      mutate(has_loss = !is.na(loss) & loss > 0) |>
      count(type, has_loss) |>
      filter(has_loss) |>
      arrange(desc(n)) |>
      head(10)
  )
  stop("No rows after filtering. Check HAZARD_KEEP or the loss column mapping.")
}

# Declustering (1 row per event id)
df_event <- df_event0 |>
  dplyr::group_by(event_id) |>
  dplyr::slice_max(order_by = loss_usd, n = 1, with_ties = FALSE) |>
  dplyr::ungroup()

# ---- 2) Annual maxima ----
annual_max <- df_event |>
  dplyr::group_by(year) |>
  dplyr::summarise(x = max(loss_usd, na.rm = TRUE), .groups = "drop") |>
  dplyr::arrange(year) |>
  dplyr::filter(is.finite(x), x > 0)

if (nrow(annual_max) < 10)
  warning("Only ", nrow(annual_max), " years available; GEV may be unstable.")

# ---- First plot (keep as-is) ----
png("annual_max_timeseries.png", width = 1600, height = 1000, res = 180)
plot(annual_max$year, annual_max$x, type = "b", pch = 19,
     xlab = "Year", ylab = "Annual maximum loss (USD)",
     main = "Annual Maximum Losses — Time Series")
dev.off()
cat("Saved: annual_max_timeseries.png\n")

# =========================================================
# 3) Robust GEV fit (MLE with L-moments starts; fallback)
# =========================================================
fit_gev_lmom <- try(extRemes::fevd(annual_max$x, type = "GEV", method = "Lmoments"),
                    silent = TRUE)

get_num <- function(co, nm) {
  if (is.null(co)) return(NA_real_)
  if (nm %in% names(co)) return(as.numeric(co[[nm]]))
  m <- grep(paste0("^", nm), names(co))
  if (length(m)) return(as.numeric(co[[m[1]]]))
  NA_real_
}

init_list <- NULL
if (!inherits(fit_gev_lmom, "try-error")) {
  co0 <- coef(fit_gev_lmom)
  init_list <- list(
    location = get_num(co0, "location"),
    scale    = get_num(co0, "scale"),
    shape    = get_num(co0, "shape")
  )
}

fit_gev <- try(
  if (!is.null(init_list) && all(is.finite(unlist(init_list)))) {
    extRemes::fevd(annual_max$x, type = "GEV", method = "MLE",
                   initial = init_list, control = list(maxit = 10000))
  } else {
    extRemes::fevd(annual_max$x, type = "GEV", method = "MLE",
                   control = list(maxit = 10000))
  },
  silent = TRUE
)

used_method <- "MLE"
if (inherits(fit_gev, "try-error")) {
  fit_gev <- fit_gev_lmom
  used_method <- "Lmoments"
}

# ---- Safely extract parameters; if missing, use lmomco fallback ----
co <- try(coef(fit_gev), silent = TRUE)

if (inherits(co, "try-error") || is.null(co)) {
  # L-moments (lmomco): sample L-moments + GEV from L-moments
  LMs  <- lmomco::lmoms(annual_max$x)        # <— correct sampler in lmomco
  para <- lmomco::pargev(LMs)                # estimates GEV params
  # lmomco::pargev returns para$para = c(xi, alpha, kappa) = (shape, scale, location)
  xi_hat    <- as.numeric(para$para[1])
  sigma_hat <- as.numeric(para$para[2])
  mu_hat    <- as.numeric(para$para[3])
  used_method <- "Lmoments"
} else {
  mu_hat    <- as.numeric(co[grep("^loc|^location", names(co))][1])
  sigma_hat <- as.numeric(co[grep("^sca|^scale",   names(co))][1])
  xi_hat    <- as.numeric(co[grep("^sha|^shape",   names(co))][1])
}

if (!all(is.finite(c(mu_hat, sigma_hat, xi_hat))) || sigma_hat <= 0)
  stop("GEV parameters not finite or invalid. Investigate data or filtering.")

cat(sprintf("\nGEV fit (%s): mu = %s, sigma = %s, xi = %s\n",
            used_method, signif(mu_hat,4), signif(sigma_hat,4), signif(xi_hat,4)))

# Wrap parameters for stable CDF/PDF/Quantiles
para_lm  <- lmomco::vec2par(c(xi_hat, sigma_hat, mu_hat), type = "gev")  # (shape, scale, location)
pgev_lm  <- function(q) lmomco::cdfgev(q, para = para_lm)
qgev_lm  <- function(p) lmomco::quagev(p, para = para_lm)
dgev_lm  <- function(x) lmomco::pdfgev(x, para = para_lm)

# =========================
# 4) Diagnostics 2×2
# =========================
x   <- sort(as.numeric(annual_max$x))
n   <- length(x)
if (n < 5) stop("Too few annual maxima (n<5) to build diagnostics.")

pi  <- (1:n)/(n+1)           # plotting positions
q_m <- qgev_lm(pi)           # model quantiles
p_m <- pgev_lm(x)            # model probs at observed data

pp_L <- qbeta(0.025, 1:n, n - (1:n) + 1)
pp_U <- qbeta(0.975, 1:n, n - (1:n) + 1)

# densities
xr  <- seq(min(x), max(x), length.out = 400)
den_data  <- density(x, n = 512, bw = "nrd0")
den_model <- data.frame(x = xr, y = dgev_lm(xr))

# RL curve and (optional) CI
T_vec <- sort(unique(T_vec))
rl_model <- qgev_lm(1 - 1/T_vec)
has_ci <- FALSE; rl_L <- rl_U <- rep(NA_real_, length(T_vec))
if (used_method == "MLE") {
  for (i in seq_along(T_vec)) {
    tmp <- try(extRemes::ci(fit_gev, return.period = T_vec[i],
                            type = "return.level", method = "delta"),
               silent = TRUE)
    if (!inherits(tmp, "try-error")) {
      vals <- suppressWarnings(as.numeric(tmp))
      names(vals) <- tolower(names(tmp))
      rl_L[i] <- vals[grep("lower", names(vals))[1]]
      rl_U[i] <- vals[grep("upper", names(vals))[1]]
      has_ci <- TRUE
    }
  }
}
rl_emp <- data.frame(T = 1/(1 - pi), RL = x)

# Save panel
png("gev_diagnostic_panel.png", width = 1800, height = 1350, res = 220)
op <- par(mfrow = c(2,2), mar = c(4.5,4.8,2.5,1.5), mgp = c(2.2,0.8,0))

# (1) QQ
plot(q_m, x, xlab = "Model Quantiles (GEV)", ylab = "Empirical Quantiles",
     pch = 1, cex = 0.9)
abline(0, 1, lwd = 1.2)
title("Q–Q Plot")

# (2) PP with 95% beta bands
plot(pi, p_m,
     xlab = "Empirical Probabilities (plotting positions)",
     ylab = "Model CDF at data (GEV)",
     pch = 19, cex = 0.7)
abline(0, 1, col = "orange", lty = 2)                     # 1–1 line
abline(lm(p_m ~ 0 + pi), col = "grey40", lwd = 1.2)       # regression
lines(pi, pp_L, col = "grey70", lty = 3)                  # 95% band
lines(pi, pp_U, col = "grey70", lty = 3)
title("P–P Plot (95% beta bands)")

# (3) Density comparison
plot(den_data, main = "Density: Data (solid) vs GEV (dashed)",
     xlab = "Annual maxima", ylab = "Density")
lines(den_model$x, den_model$y, lty = 2, col = "blue")

# (4) Return level (log-x), add CI if available
plot(T_vec, rl_model, type = "l", log = "x",
     xlab = "Return Period (years, log scale)",
     ylab = "Return Level", lwd = 1.4)
points(rl_emp$T, rl_emp$RL, pch = 1, cex = 0.9)
if (has_ci) {
  lines(T_vec, rl_L, lty = 2, col = "grey50")
  lines(T_vec, rl_U, lty = 2, col = "grey50")
}
title(sprintf("Return Levels — GEV (%s)", used_method))
mtext("GEV 2×2 diagnostics (QQ, PP, density, RL)", outer = TRUE, line = -2, cex = 1.0)

par(op)
dev.off()
cat("Saved: gev_diagnostic_panel.png\n")

# ---- 5) Print indicators ----
cat("\nGEV parameters (fitted):\n")
print(data.frame(mu = mu_hat, sigma = sigma_hat, xi = xi_hat))

rl_tbl <- data.frame(
  `T (years)`    = T_vec,
  `Return level` = rl_model,
  `Lower 95%`    = if (has_ci) rl_L else NA_real_,
  `Upper 95%`    = if (has_ci) rl_U else NA_real_
)
cat("\nReturn levels (model) and 95% CI (if available):\n")
print(rl_tbl, row.names = FALSE)


set.seed(123)
T_vec <- c(2,5,10,20,50,100)
B <- 1000
boot_RL <- matrix(NA_real_, nrow=B, ncol=length(T_vec))

for(b in 1:B){
  xs <- sample(annual_max$x, replace=TRUE)           # resample years
  para_b <- lmomco::pargev(lmomco::lmoms(xs))        # refit by L-moments
  boot_RL[b, ] <- lmomco::quagev(1 - 1/T_vec, para_b)
}

ci_L <- apply(boot_RL, 2, quantile, probs=0.025, na.rm=TRUE)
ci_U <- apply(boot_RL, 2, quantile, probs=0.975, na.rm=TRUE)

out <- data.frame(`T (years)`=T_vec,
                  `Return level`=lmomco::quagev(1 - 1/T_vec,
                                                lmomco::vec2par(c(xi_hat, sigma_hat, mu_hat),
                                                                type="gev")),
                  `Lower 95%`=ci_L, `Upper 95%`=ci_U)
print(out, row.names=FALSE)





# ============================
# VaR & CVaR from POT–GPD
# for a catalog already truncated at u0 = min(losses)
# ============================
suppressPackageStartupMessages({
  library(dplyr); library(purrr); library(scales)
})

# ---- 0) Inputs ----
x <- df_event$loss_usd                  # event losses (USD), already X > u0
years_span <- diff(range(df_event$year, na.rm = TRUE)) + 1
T_vec <- c(2,5,10,20,50,100)            # return periods to report

stopifnot(is.numeric(x), length(x) > 20, years_span >= 1)

# ---- 1) Base catalog threshold & rate (already truncated) ----
u0 <- min(x, na.rm = TRUE)              # base threshold of the catalog (USD)
lambda_base <- length(x) / years_span   # rate for X > u0  (events/year within catalog)

# ---- 2) Rescale to improve conditioning (billions of USD) ----
sf <- 1e9
xs <- x / sf
u0s <- u0 / sf

# Candidate thresholds (include u0s)
q_grid <- seq(0.90, 0.99, by = 0.01)
u_grid <- sort(unique(c(u0s, as.numeric(quantile(xs, q_grid, na.rm = TRUE)))))

min_exc     <- max(30, ceiling(0.05 * length(xs)))  # at least 5% of points or 30
min_unique  <- 8
max_abs_xi  <- 2

safe_fit_one_u <- function(u) {
  exc <- xs[xs > u]
  n_exc <- length(exc)
  out <- list(u = u, n_exc = n_exc, xi = NA_real_, beta = NA_real_, src = NA_character_)
  if (n_exc < min_exc || length(unique(exc)) < min_unique) return(out)
  
  # Try extRemes::fevd
  fevd_ok <- try(extRemes::fevd(xs, threshold = u, type = "GP", method = "MLE"), silent = TRUE)
  if (!inherits(fevd_ok, "try-error")) {
    co <- try(coef(fevd_ok), silent = TRUE)
    if (!inherits(co, "try-error") && length(co)) {
      beta <- suppressWarnings(as.numeric(co["scale"]))
      xi   <- suppressWarnings(as.numeric(co["shape"]))
      if (is.finite(beta) && beta > 0 && is.finite(xi) && abs(xi) <= max_abs_xi) {
        out$beta <- beta; out$xi <- xi; out$src <- "extRemes"
        return(out)
      }
    }
  }
  
  # Fallback: ismev::gpd.fit if available
  if (requireNamespace("ismev", quietly = TRUE)) {
    fit2 <- try(ismev::gpd.fit(xs, threshold = u, show = FALSE), silent = TRUE)
    if (!inherits(fit2, "try-error") && !is.null(fit2$mle)) {
      xi   <- suppressWarnings(as.numeric(fit2$mle[1]))
      beta <- suppressWarnings(as.numeric(fit2$mle[2]))
      if (is.finite(beta) && beta > 0 && is.finite(xi) && abs(xi) <= max_abs_xi) {
        out$beta <- beta; out$xi <- xi; out$src <- "ismev"
      }
    }
  }
  
  out
}

# Fit across thresholds safely
stab <- bind_rows(
  lapply(u_grid, function(u) {
    res <- try(safe_fit_one_u(u), silent = TRUE)
    if (is.list(res)) res else NULL
  })
) %>%
  mutate(ok = is.finite(xi) & is.finite(beta) & beta > 0 & n_exc >= min_exc & abs(xi) <= max_abs_xi)

if (!any(stab$ok, na.rm = TRUE)) {
  message("No usable GPD fits across thresholds. Reporting empirical VaR/ES only.")
  # Empirical VaR/ES on the truncated catalog
  out_tab <- tibble(
    `T (years)`     = T_vec,
    `VaR z_T (USD)` = sapply(T_vec, function(T) {
      p <- 1 - 1/T
      as.numeric(quantile(x, p, na.rm = TRUE, type = 7))
    }),
    `ES_emp (USD)`  = sapply(T_vec, function(T) {
      p <- 1 - 1/T
      thr <- as.numeric(quantile(x, p, na.rm = TRUE, type = 7))
      mean(x[x > thr], na.rm = TRUE)
    })
  ) %>% mutate(
    `VaR z_T (USD)` = comma(`VaR z_T (USD)`),
    `ES_emp (USD)`  = comma(`ES_emp (USD)`)
  )
  print(out_tab, n = nrow(out_tab))
} else {
  # Choose u* (prefer xi<1 for finite ES; otherwise best ok at highest u)
  cand <- stab %>% filter(ok, xi < 1)
  picked <- if (nrow(cand) > 0) cand %>% arrange(desc(u)) %>% slice(1) else
    stab %>% filter(ok) %>% arrange(desc(u)) %>% slice(1)
  
  u_star    <- picked$u                  # in billions
  beta_hat  <- picked$beta
  xi_hat    <- picked$xi
  n_exc_star<- picked$n_exc
  
  # Rate at u* for a truncated catalog:
  # lambda_u0 = length(x)/years_span ; P(X>u* | X>u0) = mean(x > u*)
  lambda_u <- lambda_base * mean(x > (u_star * sf), na.rm = TRUE)
  
  cat("\n--- Selected threshold & params (on rescaled data) ---\n")
  params <- tibble(
    param = c("u0 (USD)", "u* (USD)", "beta (USD)", "xi", "lambda_u* (per year)", "n_exc(u*)", "years_span"),
    value = c(comma(u0), comma(u_star * sf), comma(beta_hat * sf),
              signif(xi_hat, 4), signif(lambda_u, 4), n_exc_star, years_span)
  )
  print(params, n = nrow(params))
  
  # VaR = z_T (parametric GP) back on USD scale
  if (is.finite(xi_hat) && abs(xi_hat) > 1e-8) {
    zT_bn <- u_star + (beta_hat/xi_hat) * ((lambda_u * T_vec)^xi_hat - 1)
  } else if (is.finite(xi_hat)) {
    zT_bn <- u_star + beta_hat * log(lambda_u * T_vec)
  } else {
    zT_bn <- rep(NA_real_, length(T_vec))
  }
  zT <- zT_bn * sf
  
  # CVaR/ES (parametric if xi<1), else empirical fallback
  if (is.finite(xi_hat) && xi_hat < 1) {
    ES_T_bn <- zT_bn + (beta_hat + xi_hat * (zT_bn - u_star)) / (1 - xi_hat)
    ES_T    <- ES_T_bn * sf
    ES_emp  <- rep(NA_real_, length(T_vec))
    es_note <- "Parametric ES (ξ<1)"
  } else {
    ES_T <- rep(NA_real_, length(T_vec))
    ES_emp <- sapply(T_vec, function(T) {
      p <- 1 - 1/T
      thr <- as.numeric(quantile(x, p, na.rm = TRUE, type = 7))
      mean(x[x > thr], na.rm = TRUE)
    })
    es_note <- "ξ ≥ 1 → ES infinite; showing empirical ES"
  }
  
  # add empirical ES alongside the parametric one (even when xi < 1)
  ES_emp_vec <- sapply(zT, function(th) mean(x[x > th], na.rm = TRUE))
  
  out_tab <- tibble::tibble(
    `T (years)`        = T_vec,
    `VaR z_T (USD)`    = zT,
    `CVaR / ES (USD)`  = ES_T,
    `ES_emp (USD)`     = ES_emp_vec
  ) %>%
    dplyr::mutate(
      `VaR z_T (USD)`   = scales::comma(`VaR z_T (USD)`),
      `CVaR / ES (USD)` = ifelse(is.na(`CVaR / ES (USD)`), NA, scales::comma(`CVaR / ES (USD)`)),
      `ES_emp (USD)`    = ifelse(is.na(`ES_emp (USD)`),   NA, scales::comma(`ES_emp (USD)`))
    )
  
  cat("\n--- VaR & ES (parametric + empirical) ---\n")
  print(out_tab, n = nrow(out_tab))
  


