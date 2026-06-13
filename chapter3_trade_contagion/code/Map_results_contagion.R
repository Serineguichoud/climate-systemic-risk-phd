# ---------------------------
# Packages
# ---------------------------
pkgs <- c("dplyr","ggplot2","sf","rnaturalearth","rnaturalearthdata","readr","tibble","scales")
to_install <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if(length(to_install) > 0) install.packages(to_install)

library(dplyr)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(scales)

# 1) Carte monde
world <- ne_countries(scale = "medium", returnclass = "sf") %>%
  select(iso_a3, name_long, geometry)

# 2) Préparer tes données (ISO3)
df2 <- df %>%
  transmute(
    iso3 = toupper(trimws(Noeud)),
    impact = as.numeric(Impact_Total)
  )

# 3) Join
map_df <- world %>%
  left_join(df2, by = c("iso_a3" = "iso3"))

# 4) Plot (souvent mieux en log si très dispersé)
p_map <- ggplot(map_df) +
  geom_sf(aes(fill = impact), color = "gray75", linewidth = 0.1) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Impact total by country",
    subtitle = "Source: df (ISO3 join)",
    fill = "Impact",
    caption = "Base map: Natural Earth"
  ) +
  scale_fill_continuous(labels = label_number(big.mark = " ", decimal.mark = ","))

print(p_map)

# 5) Export
ggsave("map_world.png", plot = p_map, width = 12, height = 7, dpi = 300)
ggsave("map_world.pdf", plot = p_map, width = 12, height = 7)






library(dplyr)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(scales)

world <- ne_countries(scale = "medium", returnclass = "sf") %>%
  select(iso_a3, name_long, geometry)

df_pct <- df %>%
  transmute(
    iso3 = toupper(trimws(Noeud)),
    impact = as.numeric(Impact_Total)
  ) %>%
  group_by(iso3) %>%
  summarise(impact = sum(impact), .groups = "drop") %>%
  mutate(share_pct = 100 * impact / sum(impact, na.rm = TRUE))

# Patch ISO côté carte (pour éviter les iso_a3 atypiques)
world2 <- world %>%
  mutate(
    iso_fix = case_when(
      grepl("^France", name_long) ~ "FRA",
      name_long == "Norway" ~ "NOR",
      TRUE ~ iso_a3
    )
  )

map_df <- world2 %>%
  left_join(df_pct, by = c("iso_fix" = "iso3"))

p_pct <- ggplot(map_df) +
  geom_sf(aes(fill = share_pct), color = "gray75", linewidth = 0.1) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Share of total losses by country",
    subtitle = "Indicator = 100 × Impact_Total / sum(Impact_Total)",
    fill = "% of total",
    caption = "Base map: Natural Earth"
  ) +
  scale_fill_continuous(labels = label_number(suffix = "%", accuracy = 0.01))

print(p_pct)
ggsave("map_world_share_pct.png", plot = p_pct, width = 12, height = 7, dpi = 300)

