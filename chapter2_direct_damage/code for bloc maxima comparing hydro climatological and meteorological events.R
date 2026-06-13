# =========================================================
# EM-DAT → Annual block maxima by Disaster Type
# Keep ONLY Climatological / Hydrological / Meteorological
# One plot in R graphics window
# Requires `emdat` or `em_dat` in memory.
# =========================================================
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(lubridate)
})

# ---------- pick the data object ----------
raw <- if (exists("emdat")) emdat else if (exists("em_dat")) em_dat else
  stop("I can't find `emdat` or `em_dat` in your session.")

# ---------- small helpers ----------
pick_col <- function(cands, nm) {
  for (c in cands) if (c %in% nm) return(c)
  norm <- function(x) gsub("\\s+"," ", gsub("\\.+"," ", trimws(tolower(x))))
  nm_n <- norm(nm)
  for (cand in cands) {
    j <- which(nm_n == norm(cand)); if (length(j)) return(nm[j[1]])
  }
  NA_character_
}
`%||%` <- function(a,b) if (is.null(a)) b else a

nm <- names(raw)
col_id    <- pick_col(c("DisNo.","Dis No.","Disaster No.","DisNo","Disaster No"), nm)
col_year  <- pick_col(c("Start Year","Year","StartYear"), nm)
col_month <- pick_col(c("Start Month","Month","StartMonth"), nm)
col_day   <- pick_col(c("Start Day","Day","StartDay"), nm)
col_type  <- pick_col(c("Disaster Type","Type"), nm)  # category we plot
col_name  <- pick_col(c("Disaster Name","Event Name","Name"), nm)
col_loss  <- pick_col(c(
  "Total Damage, Adjusted ('000 US$)","Total Damages, Adjusted ('000 US$)",
  "Total Damage ('000 US$)","Total Damages ('000 US$)",
  "Total Damage, Adjusted (000 US$)","Total Damages, Adjusted (000 US$)",
  "Total Damage (000 US$)","Total Damages (000 US$)"
), nm)

# NEW: Group/Subgroup columns to filter the universe of events
col_group    <- pick_col(c("Disaster Group","Group"), nm)
col_subgroup <- pick_col(c("Disaster Subgroup","Subgroup"), nm)

if (any(is.na(c(col_year, col_type, col_loss))))
  stop("Missing essential columns. Check EM-DAT headers (year/type/damage).")

# ---------- build base table ----------
yr <- suppressWarnings(as.integer(raw[[col_year]]))
mo <- if (!is.na(col_month)) suppressWarnings(as.integer(raw[[col_month]])) else 1L
dy <- if (!is.na(col_day))   suppressWarnings(as.integer(raw[[col_day]]))   else 1L
date_vec <- make_date(
  year  = yr,
  month = ifelse(is.na(mo) | mo < 1, 1L, mo),
  day   = ifelse(is.na(dy) | dy < 1, 1L, dy)
)

# EM-DAT losses are in '000 US$
loss_k   <- readr::parse_number(as.character(raw[[col_loss]]),
                                locale = readr::locale(decimal_mark = "."))
loss_usd <- 1000 * loss_k

df_all <- tibble::tibble(
  event_id = if (!is.na(col_id))   as.character(raw[[col_id]])   else NA_character_,
  year     = yr,
  date     = date_vec,
  loss_usd = as.numeric(loss_usd),
  type     = as.character(raw[[col_type]] %||% NA_character_),
  name     = if (!is.na(col_name)) as.character(raw[[col_name]]) else NA_character_,
  group    = if (!is.na(col_subgroup)) as.character(raw[[col_subgroup]])
  else if (!is.na(col_group)) as.character(raw[[col_group]])
  else NA_character_
) %>%
  mutate(type = trimws(type),
         group = trimws(tolower(group))) %>%
  filter(!is.na(year), year > 0, is.finite(loss_usd), loss_usd > 0, !is.na(type))

# ---------- FILTER to Clim/Hydro/Meteo universe ----------
keep_groups <- c("climatological","hydrological","meteorological")
df <- df_all %>%
  filter(!is.na(group), group %in% keep_groups)

if (nrow(df) == 0) stop("After filtering to Climatological/Hydrological/Meteorological, no rows remain.")

# Decluster by event (keep max per id) if id exists
if (all(!is.na(df$event_id))) {
  df <- df %>%
    group_by(event_id) %>%
    slice_max(loss_usd, n = 1, with_ties = FALSE) %>%
    ungroup()
}

# ---------- overall annual maxima (within filtered universe) ----------
annual_all <- df %>%
  group_by(year) %>%
  slice_max(loss_usd, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(year, loss_usd, type_overall = type, name_overall = name)

# ---------- which types to overlay ----------
types_to_plot <- c("Drought","Wildfire","Flood","Storm","Extreme temperature")
types_present <- intersect(types_to_plot, sort(unique(df$type)))
if (!length(types_present)) {
  types_present <- names(sort(table(df$type), decreasing = TRUE))[1:min(6, dplyr::n_distinct(df$type))]
  message("Requested types not present after filtering; showing top present types: ",
          paste(types_present, collapse=", "))
}

# Annual maxima per selected type
annual_by_type <- lapply(types_present, function(tp) {
  sub <- df %>% filter(tolower(type) == tolower(tp))
  if (!nrow(sub)) return(NULL)
  sub %>% group_by(year) %>% slice_max(loss_usd, n = 1, with_ties = FALSE) %>%
    ungroup() %>% mutate(type = tp) %>% select(year, loss_usd, type, name)
})
names(annual_by_type) <- types_present
annual_by_type <- annual_by_type[!vapply(annual_by_type, is.null, FALSE)]

# Driver years (type max == overall max) within filtered universe
drivers <- lapply(annual_by_type, function(ann) {
  j <- inner_join(ann %>% rename(loss_t = loss_usd),
                  annual_all %>% rename(loss_o = loss_usd),
                  by = "year")
  idx <- is.finite(j$loss_t) & is.finite(j$loss_o) & abs(j$loss_t - j$loss_o) < 1e-6
  list(year = j$year[idx], loss = j$loss_t[idx])
})

# ---------- plotting ----------
xrange <- range(df$year, na.rm = TRUE)
yrange <- range(df$loss_usd, na.rm = TRUE)

base_cols <- c("#1b9e77","#d95f02","#2c7fb8","#7570b3","#e7298a","#66a61e","#e6ab02","#a6761d")
cols <- setNames(base_cols[seq_along(annual_by_type)], names(annual_by_type))
pch_map <- setNames(c(19,17,15,18,0,1,2,5)[seq_along(annual_by_type)], names(annual_by_type))

op <- par(mar = c(4.5, 5, 4, 10), xaxs = "i"); on.exit(par(op), add = TRUE)
plot(NA, xlim = xrange, ylim = yrange,
     xlab = "Year", ylab = "Loss (USD)",
     main = "EM-DAT Annual Block Maxima (Climatological/Hydrological/Meteorological)")

# context points = only selected universe
points(df$year, df$loss_usd, pch = 20, col = rgb(0,0,0,0.12))

# overall annual maxima (within selected universe)
aa <- annual_all[order(annual_all$year), ]
lines(aa$year, aa$loss_usd, lwd = 2, col = "black")
points(aa$year, aa$loss_usd, pch = 19, col = "black")

# each selected type
for (tp in names(annual_by_type)) {
  an <- annual_by_type[[tp]][order(annual_by_type[[tp]]$year), ]
  if (!nrow(an)) next
  if (nrow(an) > 1) {
    for (i in 2:nrow(an)) {
      segments(an$year[i-1], an$loss_usd[i-1], an$year[i], an$loss_usd[i],
               col = adjustcolor(cols[tp], alpha.f = 0.8), lwd = 1.7)
    }
  }
  points(an$year, an$loss_usd, pch = pch_map[tp], cex = 1.3, col = cols[tp])
  d <- drivers[[tp]]
  if (length(d$year)) points(d$year, d$loss, pch = 21, bg = "white", col = cols[tp], lwd = 2, cex = 1.6)
}

legend("topleft",
       legend = c("All selected events (context)", "Overall annual max", names(annual_by_type), "○ = type drives overall max"),
       col    = c(rgb(0,0,0,0.4), "black", cols, "black"),
       pch    = c(20, 19, unname(pch_map), 21),
       pt.bg  = c(NA, NA, rep(NA, length(annual_by_type)), "white"),
       lwd    = c(NA, 2, rep(1.7, length(annual_by_type)), 2),
       bty = "n")

# quick console recap
cat("\n--- Summary (filtered to Clim/Hydro/Meteo) ---\n")
cat("Years: ", min(aa$year), "–", max(aa$year), " (", nrow(aa), " years)\n", sep = "")
cat("# events after declustering: ", nrow(df), "\n", sep = "")
med_all <- format(round(median(aa$loss_usd)), big.mark=",")
cat("Median annual max (selected universe): ", med_all, " USD\n", sep = "")
for (tp in names(annual_by_type)) {
  d <- drivers[[tp]]
  msg <- if (!length(d$year)) "none" else paste(d$year, collapse=", ")
  cat(sprintf("Type %-20s driver years: %s\n", paste0("'",tp,"'"), msg))
}

