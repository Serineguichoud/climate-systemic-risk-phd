# ANALYSE COVARIABLES CLIMATIQUES (ERA5/indices) AVEC EMDAT (2000–2023)

library(dplyr)
library(ggplot2)
library(Hmisc)
library(stringr)
library(tidyr)

select <- dplyr::select
filter <- dplyr::filter

# 1) Nettoyage / harmonisation des covariables climatiques ---------------------

# Fonction utilitaire: convertir nombres au format "FR" (virgule décimale, espaces)
clean_num <- function(x) {
  if (is.numeric(x)) return(x)
  x <- as.character(x)
  x <- gsub("\\u00A0", " ", x, useBytes = TRUE)       # NBSP -> espace
  x <- gsub("[[:space:]]", "", x)                     # retirer espaces de milliers
  x <- gsub(",", ".", x, fixed = TRUE)                # virgule -> point
  suppressWarnings(as.numeric(x))
}

# Standardiser quelques noms de colonnes possibles
name_map <- c(
  "year"     = "year",
  "pr"       = "pr",
  "txx"      = "txx",
  "pib"      = "pib",
  "population" = "population",
  "rx5day"   = "rx5day",
  "cwd"      = "cwd",
  "rx1day"   = "rx1day",
  "r20mm"    = "r20mm",
  "cdd"      = "cdd",
  "tas"      = "tas",
  "prpcnt"   = "prpcnt",      # ton en-tête actuel
  "prpercnt" = "prpcnt",      # variante fréquente
  "prperc"   = "prpcnt"
)

# Appliquer la normalisation de noms (sans casser les colonnes non listées)
clim_covars <- covar %>%
  rename_with(.fn = function(nm) ifelse(nm %in% names(name_map), name_map[nm], nm))

# Colonnes attendues
expected_cols <- c("year","pr","txx","rx5day","cwd","rx1day","r20mm","cdd","tas","prpcnt","pib","population")
missing_cols <- setdiff(expected_cols, names(clim_covars))
if (length(missing_cols) > 0) {
  warning("⚠️ Colonnes manquantes dans 'clim_covars': ", paste(missing_cols, collapse=", "))
}

# Conversion numérique (en gardant year en entier)
num_cols <- intersect(c("pr","txx","rx5day","cwd","rx1day","r20mm","cdd","tas","prpcnt","pib","population"),
                      names(clim_covars))

clim_covars_clean <- clim_covars %>%
  mutate(across(all_of(num_cols), clean_num)) %>%
  mutate(year = clean_num(year) %>% as.integer()) %>%
  filter(!is.na(year), year >= 2000, year <= 2023) %>%
  arrange(year)

cat("✅ Covariables climatiques prêtes sur 2000–2023. Années:", clim_covars_clean$year[1], "→", clim_covars_clean$year[nrow(clim_covars_clean)], "\n")

cat("\n=== PRÉPARATION EMDAT 2000–2023 ===\n")

# Helper pour retrouver une colonne par motifs (insensible à la casse)
find_col <- function(df, patterns, stop_if_missing = TRUE) {
  nm <- names(df)
  cand <- unlist(lapply(patterns, function(p) nm[grepl(p, nm, ignore.case = TRUE)]))
  cand <- unique(cand)
  if (length(cand) == 0) {
    if (stop_if_missing) {
      stop("❌ Colonne non trouvée dans em_dat. Motifs essayés: ",
           paste(patterns, collapse = " | "))
    } else {
      return(NA_character_)
    }
  }
  cand[1]
}

# 1) Colonne année (priorité "Start Year", sinon "Year", ou colonnes proches)
year_col <- find_col(
  em_dat,
  patterns = c("^Start\\s*Year$", "^Start\\.?Year$", "\\bStart.*Year\\b", "^Year$")
)
cat("➡️ Colonne année détectée dans EMDAT:", year_col, "\n")

# 2) Colonne dommages ajustés (selon export)
# Motifs courants EM-DAT (anglais) :
# "Total Damage, Adjusted ('000 US$)" ou "Total Damages, Adjusted ('000 US$)"
# On prévoit aussi les variantes sans "Adjusted" au cas où (fallback).
damage_adj_col <- find_col(
  em_dat,
  patterns = c("Total\\s*Damage.*Adjusted.*US\\$", "Total\\s*Damages.*Adjusted.*US\\$")
  , stop_if_missing = FALSE
)

if (is.na(damage_adj_col)) {
  # Fallback non ajusté si ajusté absent
  damage_adj_col <- find_col(
    em_dat,
    patterns = c("Total\\s*Damage.*\\('000\\s*US\\$\\)", "Total\\s*Damages.*\\('000\\s*US\\$\\)")
  )
  cat("⚠️ Colonne dommages *ajustés* introuvable, fallback sur dommages non ajustés:\n   →", damage_adj_col, "\n")
} else {
  cat("➡️ Colonne dommages ajustés détectée:", damage_adj_col, "\n")
}

# 3) Préparation / nettoyage
emdat_prepared <- em_dat %>%
  mutate(
    year = as.integer(.data[[year_col]]),
    damages_thousand_usd = clean_num(.data[[damage_adj_col]]) # en milliers US$
  ) %>%
  filter(!is.na(year), year >= 2000, year <= 2023, !is.na(damages_thousand_usd)) %>%
  mutate(damages_musd = damages_thousand_usd / 1000) %>%      # millions US$
  select(year, damages_musd)

if (nrow(emdat_prepared) == 0) {
  stop("❌ Aucune donnée EMDAT exploitable entre 2000 et 2023 (colonnes année/dommages).")
}

# 4) Agrégation annuelle globale
emdat_annual <- emdat_prepared %>%
  group_by(year) %>%
  summarise(
    total_damage_adj_musd = sum(damages_musd, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(year)

cat("✅ EMDAT agrégé. Années:", emdat_annual$year[1], "→", emdat_annual$year[nrow(emdat_annual)],
    "| N années:", nrow(emdat_annual), "\n")


# 3) Construire sous-indices climatiques + indice composite -------------------

cat("\n=== CONSTRUCTION DES SOUS-INDICES CLIMATIQUES ===\n")

# Helper pour z-score par la période (ici 2000–2023)
zscore <- function(x) {
  if (all(is.na(x))) return(x)
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

# Orientations "risque +":
# - Pluie totale (pr) et prpcnt (hausse %) -> diminuent la sécheresse => on inverse
# - Les autres (txx, rx5day, rx1day, r20mm, cwd, cdd, tas) restent positifs
cc <- clim_covars_clean %>%
  mutate(
    z_pr      = -zscore(pr),
    z_prpcnt  = -zscore(prpcnt),
    z_txx     =  zscore(txx),
    z_tas     =  zscore(tas),
    z_rx5day  =  zscore(rx5day),
    z_rx1day  =  zscore(rx1day),
    z_r20mm   =  zscore(r20mm),
    z_cwd     =  zscore(cwd),
    z_cdd     =  zscore(cdd)
  ) %>%
  # Sous-indices (moyenne robuste simple)
  rowwise() %>%
  mutate(
    # Inondations / Tempêtes (proxy pluie extrême + persistance humide)
    SUB_FLOOD_STORM = mean(c(z_rx5day, z_rx1day, z_r20mm, z_cwd), na.rm = TRUE),
    
    # Chaleur (intensité + fond thermique)
    SUB_HEAT        = mean(c(z_txx, z_tas), na.rm = TRUE),
    
    # Sécheresse / Feux (sécheresse + chaleur) – pr & prpcnt inversés
    SUB_DROUGHT_FIRE = mean(c(z_cdd, z_pr, z_prpcnt, z_tas), na.rm = TRUE),
    
    # Indice composite climatique global (non pondéré)
    CLIMATE_INDEX = mean(c(SUB_FLOOD_STORM, SUB_HEAT, SUB_DROUGHT_FIRE), na.rm = TRUE)
  ) %>%
  ungroup()

cat("✅ Sous-indices construits: SUB_FLOOD_STORM, SUB_HEAT, SUB_DROUGHT_FIRE, CLIMATE_INDEX\n")

# 4) Fusion avec EMDAT ---------------------------------------------------------

analysis_df <- cc %>%
  select(year, SUB_FLOOD_STORM, SUB_HEAT, SUB_DROUGHT_FIRE, CLIMATE_INDEX) %>%
  inner_join(emdat_annual, by = "year") %>%
  arrange(year)

cat("\n=== DONNÉES FUSIONNÉES ===\n")
cat("Années communes:", paste(range(analysis_df$year), collapse=" - "),
    "| N années:", nrow(analysis_df), "\n")

# 5) Analyses de corrélation (Pearson & Spearman) -----------------------------

corr_wrap <- function(x, y) {
  p <- suppressWarnings(cor.test(x, y, method = "pearson"))
  s <- suppressWarnings(cor.test(x, y, method = "spearman"))
  c(Pearson_r = unname(p$estimate), Pearson_p = p$p.value,
    Spearman_rho = unname(s$estimate), Spearman_p = s$p.value)
}

cat("\n=== ANALYSE DE CORRÉLATION AVEC EMDAT (dommages ajustés, M$) ===\n")
cors <- rbind(
  CLIMATE_INDEX   = corr_wrap(analysis_df$CLIMATE_INDEX,   analysis_df$total_damage_adj_musd),
  SUB_FLOOD_STORM = corr_wrap(analysis_df$SUB_FLOOD_STORM, analysis_df$total_damage_adj_musd),
  SUB_HEAT        = corr_wrap(analysis_df$SUB_HEAT,        analysis_df$total_damage_adj_musd),
  SUB_DROUGHT_FIRE= corr_wrap(analysis_df$SUB_DROUGHT_FIRE,analysis_df$total_damage_adj_musd)
) %>% as.data.frame()

print(round(cors, 3))

# 6) (Optionnel) Graphique rapide de la relation composite vs dommages --------

ggplot(analysis_df, aes(CLIMATE_INDEX, total_damage_adj_musd, label = year)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  ggrepel::geom_text_repel(size = 3, max.overlaps = 15) +
  labs(
    x = "Indice climatique composite (z-score)",
    y = "Dommages EMDAT ajustés (M$ 2019 ou série ajustée)",
    title = "Relation indice climatique global vs dommages (2000–2023)"
  ) +
  theme_minimal()





# === COMPARAISON ERA5 vs EM-DAT PAR PÉRIL (2000–2023) =======================

library(dplyr)
library(stringr)
library(tidyr)
library(purrr)
library(broom)

# ---------- Helpers (si pas déjà définis) ------------------------------------
clean_num <- function(x) {
  if (is.numeric(x)) return(x)
  x <- as.character(x)
  x <- gsub("\\u00A0", " ", x, useBytes = TRUE)
  x <- gsub("[[:space:]]", "", x)
  x <- gsub(",", ".", x, fixed = TRUE)
  suppressWarnings(as.numeric(x))
}
zscore <- function(x) if (all(is.na(x))) x else (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)

find_col <- function(df, patterns, stop_if_missing = TRUE) {
  nm <- names(df)
  cand <- unlist(lapply(patterns, function(p) nm[grepl(p, nm, ignore.case = TRUE)]))
  cand <- unique(cand)
  if (length(cand) == 0) {
    if (stop_if_missing) stop("Colonne non trouvée. Motifs: ", paste(patterns, collapse=" | "))
    return(NA_character_)
  }
  cand[1]
}

# ---------- 0) Préparer les covariables climatiques --------------------------
# Accepte 'clim_covars_clean' sinon nettoie 'clim_covars'
if (!exists("clim_covars_clean")) {
  stopifnot(exists("clim_covars"))
  name_map <- c(
    "year"="year","pr"="pr","txx"="txx","rx5day"="rx5day","cwd"="cwd",
    "rx1day"="rx1day","r20mm"="r20mm","cdd"="cdd","tas"="tas",
    "prpcnt"="prpcnt","prpercnt"="prpcnt","prperc"="prpcnt"
  )
  clim_covars_clean <- clim_covars %>%
    rename_with(~ ifelse(.x %in% names(name_map), name_map[.x], .x)) %>%
    mutate(across(c(pr,txx,rx5day,cwd,rx1day,r20mm,cdd,tas,prpcnt), clean_num),
           year = as.integer(clean_num(year))) %>%
    filter(!is.na(year), year>=2000, year<=2023) %>%
    arrange(year)
}

# Z-scores orientés "risque +"
cc <- clim_covars_clean %>%
  mutate(
    z_pr      = -zscore(pr),       # moins de pluie => plus de risque sécheresse
    z_prpcnt  = -zscore(prpcnt),   # baisse % pluie => plus de risque
    z_txx     =  zscore(txx),
    z_tas     =  zscore(tas),
    z_rx5day  =  zscore(rx5day),
    z_rx1day  =  zscore(rx1day),
    z_r20mm   =  zscore(r20mm),
    z_cwd     =  zscore(cwd),
    z_cdd     =  zscore(cdd)
  ) %>%
  rowwise() %>%
  mutate(
    SUB_FLOOD_STORM   = mean(c(z_rx5day, z_rx1day, z_r20mm, z_cwd), na.rm = TRUE),
    SUB_HEAT          = mean(c(z_txx, z_tas), na.rm = TRUE),
    SUB_DROUGHT_FIRE  = mean(c(z_cdd, z_pr, z_prpcnt, z_tas), na.rm = TRUE),
    COLD_INDEX        = mean(c(-z_tas, -z_txx), na.rm = TRUE)  # proxy froid (faute d'indices dédiés)
  ) %>% ungroup()

# ---------- 1) Préparer EM-DAT par péril (via Disaster Subtype) --------------
stopifnot(exists("em_dat"))

year_col <- find_col(em_dat, c("^Start\\s*Year$","^Start\\.?Year$","\\bStart.*Year\\b","^Year$"))
damage_adj_col <- find_col(
  em_dat,
  c("Total\\s*Damage.*Adjusted.*US\\$","Total\\s*Damages.*Adjusted.*US\\$")
  , stop_if_missing = FALSE
)
if (is.na(damage_adj_col)) {
  damage_adj_col <- find_col(em_dat, c("Total\\s*Damage.*\\('000\\s*US\\$\\)","Total\\s*Damages.*\\('000\\s*US\\$\\)"))
}

emdat_perils <- em_dat %>%
  mutate(
    year  = as.integer(.data[[year_col]]),
    dmg_kusd = clean_num(.data[[damage_adj_col]]),
    subtype  = trimws(tolower(`Disaster Subtype`))
  ) %>%
  filter(!is.na(year), year>=2000, year<=2023, !is.na(dmg_kusd))

# Mapping de sous-types -> périls analytiques
map_list <- list(
  flood    = c("flood","flash flood","riverine flood","coastal flood","urban flood","landslide.*flood"),
  storm    = c("storm","tropical cyclone","hurricane","convective storm","extratropical storm","wind storm","severe storm","storm surge"),
  heat     = c("extreme temperature","heat wave","heatwave","hot spell"),
  cold     = c("cold wave","extreme cold","frost","cold spell"),
  drought  = c("drought"),
  wildfire = c("wildfire","forest fire","bush fire","wild fire")
)

classify_peril <- function(x) {
  for (p in names(map_list)) {
    if (any(stringr::str_detect(x, paste0("\\b(", paste(map_list[[p]], collapse="|"), ")\\b")))) return(p)
  }
  return(NA_character_)
}

emdat_perils <- emdat_perils %>%
  mutate(peril = vapply(subtype, classify_peril, character(1))) %>%
  filter(!is.na(peril))

emdat_peril_annual <- emdat_perils %>%
  transmute(year, peril, dmg_musd = dmg_kusd/1000) %>%
  group_by(year, peril) %>%
  summarise(total_damage_musd = sum(dmg_musd, na.rm=TRUE), .groups="drop")

# ---------- 2) Joindre indices climatiques (avec option de lag) --------------
# Choisir un décalage (lag) pour certaines paires (0 = simultané)
lag_setup <- tibble::tribble(
  ~peril,    ~lag_years,
  "flood",      0,
  "storm",      0,
  "heat",       0,
  "cold",       0,
  "drought",    1,  # sécheresse souvent à décalage saison/année
  "wildfire",   0   # tu peux tester 1 aussi si besoin
)

cc_long <- cc %>%
  select(year, SUB_FLOOD_STORM, SUB_HEAT, SUB_DROUGHT_FIRE, COLD_INDEX,
         z_rx5day, z_rx1day, z_r20mm, z_cwd, z_txx, z_tas, z_cdd, z_pr, z_prpcnt) 

# Appliquer les lags par péril
make_lagged <- function(df, peril, lag_years) {
  # on décalera seulement les colonnes climatiques (année + clim)
  df %>% mutate(year = year + lag_years, peril = peril)
}

cc_lagged <- lag_setup %>%
  group_by(peril) %>%
  group_map(~ make_lagged(cc_long, .y$peril, .x$lag_years), .keep=TRUE) %>%
  bind_rows()

analysis_peril <- emdat_peril_annual %>%
  inner_join(cc_lagged, by = c("year","peril"))

# ---------- 3) Définir les "bons" indicateurs par péril ----------------------
# (Sous-indice + composantes pertinentes)
peril_vars <- list(
  flood    = c("SUB_FLOOD_STORM","z_rx5day","z_rx1day","z_r20mm","z_cwd"),
  storm    = c("SUB_FLOOD_STORM","z_rx1day","z_rx5day"),
  heat     = c("SUB_HEAT","z_txx","z_tas"),
  cold     = c("COLD_INDEX"),  # proxy basé sur températures basses (faute de CSDI/FD)
  drought  = c("SUB_DROUGHT_FIRE","z_cdd","z_pr","z_prpcnt"),
  wildfire = c("SUB_DROUGHT_FIRE","z_cdd","z_tas")
)

# ---------- 4) Corrélations (Spearman & Pearson) par péril × variable --------
corr_one <- function(df, yvar, xvar) {
  x <- df[[xvar]]; y <- df[[yvar]]
  tibble::tibble(
    var = xvar,
    n   = sum(complete.cases(x,y)),
    spearman_r = suppressWarnings(cor(x, y, method="spearman", use="pairwise")),
    pearson_r  = suppressWarnings(cor(x, y, method="pearson",  use="pairwise")),
    spearman_p = tryCatch(broom::glance(cor.test(x, y, method="spearman"))$p.value, error=function(e) NA_real_),
    pearson_p  = tryCatch(broom::glance(cor.test(x, y, method="pearson"))$p.value,  error=function(e) NA_real_)
  )
}

results <- analysis_peril %>%
  group_by(peril) %>%
  group_modify(~ {
    vars <- peril_vars[[unique(.y$peril)]]
    bind_rows(lapply(vars, function(v) corr_one(.x, "total_damage_musd", v)))
  }) %>%
  ungroup() %>%
  arrange(peril, desc(abs(spearman_r))) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

print(results)

# ---------- 5) Lecture rapide (optionnel) ------------------------------------
cat("\nNOTE: lags appliqués ->\n")
print(lag_setup)
cat("\nInterprétation attendue :\n",
    "- flood/storm : SUB_FLOOD_STORM et composantes pluie extrême devraient dominer (+, p-val faibles)\n",
    "- heat : lien plus modéré avec SUB_HEAT ; cold : COLD_INDEX (proxy) possiblement +\n",
    "- drought : corrélation améliorée avec lag=1 ; wildfire : SUB_DROUGHT_FIRE & z_cdd\n")






# --- Étape 1 : Préparer les données climatologiques ---------------------------
library(extRemes)  # pour fevd()
library(dplyr)
library(ggplot2)

# Extraire les événements climatologiques (heat + drought + wildfire)
climato_df <- analysis_peril %>%
  filter(peril %in% c("heat", "drought", "wildfire")) %>%
  group_by(year) %>%
  summarise(
    mean_damage_musd = mean(total_damage_musd, na.rm = TRUE),
    rx5day = mean(z_rx5day, na.rm = TRUE)  # covariable ERA5 normalisée
  ) %>%
  filter(!is.na(mean_damage_musd), !is.na(rx5day))

cat("✅ Données climatologiques prêtes pour EVT non stationnaire :",
    nrow(climato_df), "années\n")

# --- Étape 2 : Ajuster un modèle GEV non stationnaire --------------------------
# μ = α0 + α1 * rx5day (position variable avec rx5day)
fit_ns <- fevd(
  x = climato_df$mean_damage_musd,
  data = climato_df,
  location.fun = ~ rx5day,
  type = "GEV"
)

# Modèle stationnaire de comparaison
fit_s <- fevd(
  x = climato_df$mean_damage_musd,
  data = climato_df,
  type = "GEV"
)

# --- Étape 3 : Comparaison statistique (LRT) -----------------------------------
lrt <- anova(fit_s, fit_ns)
print(lrt)

# --- Étape 4 : Affichage des coefficients --------------------------------------
summary(fit_ns)

# --- Étape 5 : Visualisation ---------------------------------------------------
rx_seq <- seq(min(climato_df$rx5day, na.rm = TRUE),
              max(climato_df$rx5day, na.rm = TRUE),
              length.out = 100)

pred_mu <- fit_ns$results$par[1] + fit_ns$results$par[2] * rx_seq
pred_sigma <- fit_ns$results$par[3]

plot_df <- data.frame(rx5day = rx_seq, mu = pred_mu)

ggplot(climato_df, aes(rx5day, mean_damage_musd)) +
  geom_point(size = 3) +
  geom_line(data = plot_df, aes(rx5day, mu), color = "red", linewidth = 1.1) +
  labs(
    x = "ERA5 rx5day (z-score)",
    y = "Mean climatological damage (M$)",
    title = "Non-stationary GEV: damage location parameter vs rx5day"
  ) +
  theme_minimal()





# === Non-stationary EVT on annual BLOCK SUMS (climatological) with rx5day ===
library(extRemes)
library(dplyr)
library(ggplot2)
library(broom)

# 1) Annual block sums for climatological events (heat + drought + wildfire)
climato_events <- analysis_peril %>%
  filter(peril %in% c("heat", "drought", "wildfire"))

blocksum_df <- climato_events %>%
  group_by(year) %>%
  summarise(
    block_sum_musd = sum(total_damage_musd, na.rm = TRUE),
    rx5day = mean(z_rx5day, na.rm = TRUE)   # covariate (annual z-score)
  ) %>%
  filter(is.finite(block_sum_musd), is.finite(rx5day)) %>%
  arrange(year)

cat("Annual block sums (climatological):", nrow(blocksum_df), "years\n")

# 2) Fit GEV to BLOCK SUMS: stationary vs non-stationary (location/scale ~ rx5day)
#    We'll try four specs and compare:
#    S: stationary; L: location~rx5day; Sc: scale~rx5day; LSc: both

fit_S   <- fevd(x = blocksum_df$block_sum_musd, data = blocksum_df, type = "GEV")
fit_L   <- fevd(x = blocksum_df$block_sum_musd, data = blocksum_df, type = "GEV",
                location.fun = ~ rx5day)
fit_Sc  <- fevd(x = blocksum_df$block_sum_musd, data = blocksum_df, type = "GEV",
                scale.fun    = ~ rx5day)
fit_LSc <- fevd(x = blocksum_df$block_sum_musd, data = blocksum_df, type = "GEV",
                location.fun = ~ rx5day, scale.fun = ~ rx5day)

# 3) Likelihood ratio tests vs stationary (S)
lrt_L   <- lr.test(fit_S, fit_L)
lrt_Sc  <- lr.test(fit_S, fit_Sc)
lrt_LSc <- lr.test(fit_S, fit_LSc)

cat("\n--- Likelihood ratio tests (vs stationary) ---\n")
print(lrt_L)
print(lrt_Sc)
print(lrt_LSc)

# 4) Information criteria
aics <- c(S  = AIC(fit_S), L  = AIC(fit_L), Sc = AIC(fit_Sc), LSc = AIC(fit_LSc))
bics <- c(S  = BIC(fit_S), L  = BIC(fit_L), Sc = BIC(fit_Sc), LSc = BIC(fit_LSc))

cat("\nAIC:\n"); print(round(aics, 2))
cat("\nBIC:\n"); print(round(bics, 2))

# 5) Summaries (coefficients) for each non-stationary spec
cat("\n--- Summary: Location ~ rx5day ---\n"); print(summary(fit_L))
cat("\n--- Summary: Scale ~ rx5day ---\n");    print(summary(fit_Sc))
cat("\n--- Summary: Location & Scale ~ rx5day ---\n"); print(summary(fit_LSc))

# 6) Visual: fitted μ(rx5day) for the best of L or LSc (if applicable)
#    (If the best spec modifies only scale, this plot is less relevant.)
best <- names(which.min(aics))  # pick by AIC
if (best %in% c("L", "LSc")) {
  f <- if (best == "L") fit_L else fit_LSc
  # Extract coefficients by name (extRemes names: "location:(Intercept)", "location:rx5day", etc.)
  mu0 <- f$results$par["location:(Intercept)"]
  mu1 <- f$results$par["location:rx5day"]
  rx_seq <- seq(min(blocksum_df$rx5day), max(blocksum_df$rx5day), length.out = 100)
  pred_mu <- as.numeric(mu0 + mu1 * rx_seq)
  ggplot(blocksum_df, aes(rx5day, block_sum_musd)) +
    geom_point(size = 3) +
    geom_line(data = data.frame(rx5day = rx_seq, mu = pred_mu),
              aes(rx5day, mu), linewidth = 1.1) +
    labs(x = "ERA5 rx5day (z-score)",
         y = "Annual block SUM of climatological damages (M$)",
         title = paste0("GEV on block sums (", best, "): location ~ rx5day")) +
    theme_minimal() -> p_mu
  print(p_mu)
}

# 7) Basic diagnostic plots for the chosen spec
if (best == "S") chosen <- fit_S
if (best == "L") chosen <- fit_L
if (best == "Sc") chosen <- fit_Sc
if (best == "LSc") chosen <- fit_LSc

par(mfrow = c(2,2)); plot(chosen); par(mfrow = c(1,1))








# === Non-stationary EVT on annual BLOCK SUMS (STORM) with rx5day ==============
# Prérequis: 'analysis_peril' contient (year, peril, total_damage_musd, z_rx5day)

library(extRemes)
library(dplyr)
library(ggplot2)

# 1) Annual block sums for STORM
storm_events <- analysis_peril %>%
  filter(peril == "storm")

blocksum_storm <- storm_events %>%
  group_by(year) %>%
  summarise(
    block_sum_musd = sum(total_damage_musd, na.rm = TRUE),
    rx5day = mean(z_rx5day, na.rm = TRUE)   # covariate (annual z-score)
  ) %>%
  filter(is.finite(block_sum_musd), is.finite(rx5day)) %>%
  arrange(year)

cat("Annual block sums (STORM):", nrow(blocksum_storm), "years\n")

# 2) Fit GEV to BLOCK SUMS: stationary vs non-stationary (location/scale ~ rx5day)
fit_S   <- fevd(x = blocksum_storm$block_sum_musd, data = blocksum_storm, type = "GEV")
fit_L   <- fevd(x = blocksum_storm$block_sum_musd, data = blocksum_storm,
                type = "GEV", location.fun = ~ rx5day)
fit_Sc  <- fevd(x = blocksum_storm$block_sum_musd, data = blocksum_storm,
                type = "GEV", scale.fun    = ~ rx5day)
fit_LSc <- fevd(x = blocksum_storm$block_sum_musd, data = blocksum_storm,
                type = "GEV", location.fun = ~ rx5day, scale.fun = ~ rx5day)

# 3) Likelihood ratio tests (vs stationary)
lrt_L   <- lr.test(fit_S, fit_L)
lrt_Sc  <- lr.test(fit_S, fit_Sc)
lrt_LSc <- lr.test(fit_S, fit_LSc)

cat("\n--- Likelihood ratio tests (vs stationary) ---\n")
print(lrt_L); print(lrt_Sc); print(lrt_LSc)

# 4) Information criteria  (utiliser les champs $AIC / $BIC des objets fevd)
aics <- c(S  = fit_S$AIC,  L  = fit_L$AIC,  Sc = fit_Sc$AIC,  LSc = fit_LSc$AIC)
bics <- c(S  = fit_S$BIC,  L  = fit_L$BIC,  Sc = fit_Sc$BIC,  LSc = fit_LSc$BIC)

cat("\nAIC:\n"); print(round(aics, 3))
cat("\nBIC:\n"); print(round(bics, 3))

# 5) Résumés des modèles non stationnaires (pour lecture des coefficients)
cat("\n--- Summary: Location ~ rx5day ---\n"); print(summary(fit_L))
cat("\n--- Summary: Scale ~ rx5day ---\n");    print(summary(fit_Sc))
cat("\n--- Summary: Location & Scale ~ rx5day ---\n"); print(summary(fit_LSc))

# 6) Choix du modèle par AIC et visualisation de μ(rx5day) si pertinent
best <- names(which.min(aics))
cat("\nBest by AIC:", best, "\n")

# Fonction utilitaire pour extraire mu0/mu1 quels que soient les noms utilisés
get_mu_coefs <- function(fit) {
  parn <- names(fit$results$par)
  mu0_name <- if ("mu0" %in% parn) "mu0" else "location:(Intercept)"
  mu1_name <- if ("mu1" %in% parn) "mu1" else "location:rx5day"
  c(mu0 = as.numeric(fit$results$par[mu0_name]),
    mu1 = as.numeric(fit$results$par[mu1_name]))
}

if (best %in% c("L", "LSc")) {
  f <- if (best == "L") fit_L else fit_LSc
  mu <- get_mu_coefs(f)
  rx_seq  <- seq(min(blocksum_storm$rx5day), max(blocksum_storm$rx5day), length.out = 100)
  pred_mu <- mu["mu0"] + mu["mu1"] * rx_seq
  
  p_mu <- ggplot(blocksum_storm, aes(rx5day, block_sum_musd)) +
    geom_point(size = 3) +
    geom_line(data = data.frame(rx5day = rx_seq, mu = pred_mu),
              aes(rx5day, mu), linewidth = 1.1) +
    labs(x = "ERA5 rx5day (z-score)",
         y = "Annual block SUM of STORM damages (M$)",
         title = paste0("GEV on block sums (", best, "): location ~ rx5day")) +
    theme_minimal()
  print(p_mu)
}

# 7) Diagnostics pour le modèle retenu
chosen <- switch(best, S = fit_S, L = fit_L, Sc = fit_Sc, LSc = fit_LSc)
par(mfrow = c(2,2)); plot(chosen); par(mfrow = c(1,1))





# --- LRT (déjà ok) ---
lrt_L   <- lr.test(fit_S, fit_L)
lrt_Sc  <- lr.test(fit_S, fit_Sc)
lrt_LSc <- lr.test(fit_S, fit_LSc)
cat("\n--- LRT vs stationary ---\n"); print(lrt_L); print(lrt_Sc); print(lrt_LSc)

# --- AIC/BIC robustes : utiliser $AIC/$BIC si présents, sinon recalculer ---
get_ic <- function(fit, k, n) {
  nllh <- if (!is.null(fit$nllh)) fit$nllh else fit$results$nllh
  AICv <- if (!is.null(fit$AIC)) fit$AIC else 2*nllh + 2*k
  BICv <- if (!is.null(fit$BIC)) fit$BIC else 2*nllh + k*log(n)
  c(AIC = as.numeric(AICv), BIC = as.numeric(BICv))
}

n <- nrow(blocksum_storm)
ic_S   <- get_ic(fit_S,   k = 3, n = n)   # (mu, sigma, xi)
ic_L   <- get_ic(fit_L,   k = 4, n = n)   # (+ mu1)
ic_Sc  <- get_ic(fit_Sc,  k = 4, n = n)   # (+ sigma1)
ic_LSc <- get_ic(fit_LSc, k = 5, n = n)   # (+ mu1 + sigma1)

aics <- c(S = ic_S["AIC"], L = ic_L["AIC"], Sc = ic_Sc["AIC"], LSc = ic_LSc["AIC"])
bics <- c(S = ic_S["BIC"], L = ic_L["BIC"], Sc = ic_Sc["BIC"], LSc = ic_LSc["BIC"])

cat("\nAIC:\n"); print(round(aics, 3))
cat("\nBIC:\n"); print(round(bics, 3))

# --- Sélection safe ---
best <- names(which.min(aics))
cat("\nBest by AIC:", best, "\n")

# --- Visualiser μ(rx5day) seulement si μ non stationnaire ---
if (!is.null(best) && length(best) == 1 && best %in% c("L", "LSc")) {
  get_mu <- function(fit) {
    parn <- names(fit$results$par)
    mu0_name <- if ("mu0" %in% parn) "mu0" else "location:(Intercept)"
    mu1_name <- if ("mu1" %in% parn) "mu1" else "location:rx5day"
    c(mu0 = as.numeric(fit$results$par[mu0_name]),
      mu1 = as.numeric(fit$results$par[mu1_name]))
  }
  f <- if (best == "L") fit_L else fit_LSc
  mu <- get_mu(f)
  rx_seq  <- seq(min(blocksum_storm$rx5day), max(blocksum_storm$rx5day), length.out = 100)
  pred_mu <- mu["mu0"] + mu["mu1"] * rx_seq
  
  library(ggplot2)
  p_mu <- ggplot(blocksum_storm, aes(rx5day, block_sum_musd)) +
    geom_point(size = 3) +
    geom_line(data = data.frame(rx5day = rx_seq, mu = pred_mu),
              aes(rx5day, mu), linewidth = 1.1) +
    labs(x = "ERA5 rx5day (z-score)",
         y = "Annual block SUM of STORM damages (M$)",
         title = paste0("GEV on block sums (", best, "): location ~ rx5day")) +
    theme_minimal()
  print(p_mu)
}

# --- Diagnostics pour le modèle retenu ---
chosen <- switch(best, S = fit_S, L = fit_L, Sc = fit_Sc, LSc = fit_LSc, fit_S)
par(mfrow = c(2,2)); plot(chosen); par(mfrow = c(1,1))


