# Reproduce fig_geo_vs_network_final.pdf
# Run from the repo/ directory:
#   Rscript make_figure5.R
#
# Expects (after Rscript precompute_metro_metrics.R):
#   data/{cuebiq,lodes,uk,spain,latam}/
#     metro_index.csv
#     od/{metro_id}.parquet
#     centroids/{metro_id}.parquet
#     metrics/{metro_id}.parquet

suppressPackageStartupMessages({
  library(dplyr)
  library(arrow)
  library(ggplot2)
  library(patchwork)
  library(sf)
  library(ggrepel)
  library(scales)
  library(grid)
})

sf::sf_use_s2(FALSE)
source("R/od_zone_metrics.R")
dir.create("output", recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Load precomputed metro metrics (entropy + local ρ)
# ---------------------------------------------------------------------------
DATASETS <- tribble(
  ~folder,  ~area_type,
  "cuebiq", "US Mobility",
  "lodes",  "US Commuting",
  "uk",     "UK Commuting",
  "spain",  "Spain Mobility",
  "latam",  "LatAm Mobility"
)

metric_files <- unlist(lapply(DATASETS$folder, function(folder) {
  list.files(
    file.path("data", folder, "metrics"),
    pattern = "\\.parquet$",
    full.names = TRUE
  )
}))
if (length(metric_files) == 0L) {
  stop(
    "No metrics found under data/*/metrics/. ",
    "Run: Rscript precompute_metro_metrics.R",
    call. = FALSE
  )
}

metro <- bind_rows(lapply(metric_files, read_parquet)) |>
  filter(
    n_zones > 40,
    msa_name != "Buenos Aires",
    is.finite(entropy),
    is.finite(mean_rho_local_d5)
  ) |>
  mutate(
    area_type = factor(as.character(area_type), levels = DATASETS$area_type),
    mean_rho = mean_rho_local_d5
  )

message(
  "Loaded metrics for ", nrow(metro), " metros:\n",
  paste(capture.output(count(metro, area_type)), collapse = "\n")
)

# ---------------------------------------------------------------------------
# City example panels: build pairs from OD + centroids
# ---------------------------------------------------------------------------
centroids_from_parquet <- function(path) {
  tab <- read_parquet(path)
  st_as_sf(tab, coords = c("lon", "lat"), crs = 4326, remove = FALSE) |>
    transmute(GEOID = as.character(GEOID), geometry)
}

load_city_pairs <- function(folder, metro_id, msa_name) {
  od <- read_parquet(file.path("data", folder, "od", paste0(metro_id, ".parquet")))
  cents <- centroids_from_parquet(
    file.path("data", folder, "centroids", paste0(metro_id, ".parquet"))
  )
  res <- compute_city_error_od_zones(
    od = od,
    centroids_sf = cents,
    area_name = msa_name,
    cbsa = metro_id,
    min_zones = 40L,
    beta = 0.8,
    K = 5L,
    return_pairs = TRUE
  )
  if (is.null(res)) {
    stop("Failed to compute pairs for ", metro_id, call. = FALSE)
  }
  res
}

message("Computing NYC / Mexico City / London pairs from OD networks...")
paper_nyc <- load_city_pairs("cuebiq", "35620", "New York")
paper_mexico <- load_city_pairs("latam", "ciudad-de-mexico", "Mexico City")
paper_london <- load_city_pairs("uk", "london", "London")

nyc_pairs <- paper_nyc$pairs
mexico_pairs <- paper_mexico$pairs
london_pairs <- paper_london$pairs

# ---------------------------------------------------------------------------
# NYC rank-mismatch inset (focal tract + neighbors within 10 km)
# ---------------------------------------------------------------------------
radius_km <- 10
near_km <- 2
min_neighbors <- 50L
min_neighbors_near <- 25L
nyc_borough_counties <- c("36005", "36047", "36061", "36081", "36085")
exclude_geoids <- c("36005000100")

sym_pairs <- nyc_pairs |>
  filter(is.finite(d_geo_km), d_geo_km > 0, is.finite(d_net), d_net > 0) |>
  bind_rows(nyc_pairs |> rename(geoid_i = geoid_j, geoid_j = geoid_i))

nbr <- sym_pairs |>
  filter(d_geo_km <= radius_km) |>
  group_by(geoid_i) |>
  mutate(
    geo_rank = rank(d_geo_km, ties.method = "average"),
    net_rank = rank(d_net, ties.method = "average")
  ) |>
  ungroup()

focal_id <- "36081047600"

focal_nbrs <- nbr |>
  filter(geoid_i == focal_id) |>
  arrange(geo_rank) |>
  mutate(
    geo_rank_plot = rank(d_geo_km, ties.method = "first"),
    net_rank_plot = rank(d_net, ties.method = "first")
  )

zone_ids <- unique(c(focal_id, focal_nbrs$geoid_j))

# Tract polygons for the inset: prefer shipped subset, refresh if incomplete.
load_nyc_inset_tracts <- function(zone_ids) {
  shipped <- "data/cuebiq/nyc_inset_tracts.geojson"
  zone_ids <- unique(as.character(zone_ids))

  read_or_download_states <- function(states) {
    cache_dir <- file.path(path.expand("~"), "Library", "Caches", "tigris")
    pieces <- lapply(states, function(st) {
      shp <- file.path(cache_dir, paste0("cb_2019_", st, "_tract_500k.shp"))
      if (file.exists(shp)) {
        message("  reading tracts from tigris cache: ", basename(shp))
        return(st_read(shp, quiet = TRUE))
      }
      message("  downloading tracts via tigris for state ", st)
      options(tigris_use_cache = TRUE)
      tigris::tracts(state = st, cb = TRUE, year = 2019, progress = FALSE)
    })
    bind_rows(pieces) |>
      filter(GEOID %in% zone_ids) |>
      st_transform(4326) |>
      st_make_valid() |>
      transmute(GEOID = as.character(GEOID), geometry = geometry)
  }

  tracts <- NULL
  if (file.exists(shipped)) {
    tracts <- st_read(shipped, quiet = TRUE) |>
      mutate(GEOID = as.character(GEOID)) |>
      filter(GEOID %in% zone_ids)
  }
  missing <- setdiff(zone_ids, if (is.null(tracts)) character() else tracts$GEOID)
  if (length(missing) > 0L) {
    message(
      "  inset tracts incomplete (", length(missing),
      " missing); rebuilding for current neighborhood"
    )
    tracts <- read_or_download_states(unique(substr(zone_ids, 1L, 2L)))
    dir.create(dirname(shipped), recursive = TRUE, showWarnings = FALSE)
    st_write(tracts, shipped, delete_dsn = TRUE, quiet = TRUE)
    message("  wrote ", shipped, " (", nrow(tracts), " tracts)")
  }

  tracts
}

tracts <- load_nyc_inset_tracts(zone_ids)


map_df <- tracts |>
  left_join(
    focal_nbrs |> select(GEOID = geoid_j, geo_rank_plot, net_rank_plot, d_geo_km, d_net),
    by = "GEOID"
  ) |>
  mutate(role = if_else(GEOID == focal_id, "focal", "neighbor")) |>
  # Neighbors must have ranks; drop any leftover polygons without them.
  filter(role == "focal" | (is.finite(geo_rank_plot) & is.finite(net_rank_plot)))

rho_entropy_colors <- c(
  "US Mobility" = "#2C7FB8",
  "US Commuting" = "#7B3294",
  "UK Commuting" = "#D95F02",
  "Spain Mobility" = "#1B9E77",
  "LatAm Mobility" = "#E69F00"
)
CITY_FILL <- c(
  nyc = unname(rho_entropy_colors[["US Mobility"]]),
  mexico = unname(rho_entropy_colors[["LatAm Mobility"]]),
  london = unname(rho_entropy_colors[["UK Commuting"]])
)
HIGHLIGHT_FILL <- "#CC79A7"
HIGHLIGHT_COLOR <- "#4A1942"
HIGHLIGHT_DOT_FILL <- "#9E9E9E"
HIGHLIGHT_DOT_COLOR <- "#424242"

# ---------------------------------------------------------------------------
# Entropy panel (right)
# ---------------------------------------------------------------------------
plot_rho <- metro |>
  mutate(
    label_city = case_when(
      area_type == "US Mobility" & grepl("^New York", msa_name, ignore.case = TRUE) ~ "New York",
      area_type == "US Mobility" & grepl("^Los Angeles", msa_name, ignore.case = TRUE) ~ "Los Angeles",
      area_type == "US Mobility" & grepl("^Chicago", msa_name, ignore.case = TRUE) ~ "Chicago",
      area_type == "US Mobility" & grepl("^Boston", msa_name, ignore.case = TRUE) ~ "Boston",
      area_type == "US Mobility" & grepl("^Houston", msa_name, ignore.case = TRUE) ~ "Houston",
      area_type == "US Commuting" & grepl("^New York", msa_name, ignore.case = TRUE) ~ "New York",
      area_type == "US Commuting" & grepl("^Los Angeles", msa_name, ignore.case = TRUE) ~ "Los Angeles",
      area_type == "US Commuting" & grepl("^Chicago", msa_name, ignore.case = TRUE) ~ "Chicago",
      area_type == "US Commuting" & grepl("^Boston", msa_name, ignore.case = TRUE) ~ "Boston",
      area_type == "UK Commuting" & grepl("London", msa_name, ignore.case = TRUE) ~ "London",
      area_type == "UK Commuting" & grepl("Manchester", msa_name, ignore.case = TRUE) ~ "Manchester",
      area_type == "UK Commuting" & grepl("Birmingham", msa_name, ignore.case = TRUE) ~ "Birmingham",
      area_type == "UK Commuting" & grepl("Leeds", msa_name, ignore.case = TRUE) ~ "Leeds",
      area_type == "UK Commuting" & grepl("Liverpool", msa_name, ignore.case = TRUE) ~ "Liverpool",
      area_type == "Spain Mobility" & msa_name == "Madrid" ~ "Madrid",
      area_type == "Spain Mobility" & msa_name == "Barcelona" ~ "Barcelona",
      area_type == "Spain Mobility" & msa_name == "Valencia" ~ "Valencia",
      area_type == "Spain Mobility" & msa_name == "Sevilla" ~ "Sevilla",
      area_type == "LatAm Mobility" & msa_name == "Mexico City" ~ "Mexico City",
      area_type == "LatAm Mobility" & msa_name == "São Paulo" ~ "São Paulo",
      area_type == "LatAm Mobility" & msa_name == "Bogotá" ~ "Bogotá",
      area_type == "LatAm Mobility" & msa_name == "Santiago" ~ "Santiago",
      TRUE ~ NA_character_
    )
  )

label_rho <- plot_rho |>
  filter(!is.na(label_city)) |>
  arrange(desc(population), label_city) |>
  mutate(
    .nudge_x = c(0.018, -0.018, 0.022, -0.015, 0.012, -0.022)[
      ((row_number() - 1L) %% 6L) + 1L
    ],
    .nudge_y = c(0.055, -0.045, 0.040, -0.055, 0.048, -0.035)[
      ((row_number() - 1L) %% 6L) + 1L
    ]
  )
bg_rho <- plot_rho |> filter(is.na(label_city))

p_entropy <- ggplot(plot_rho, aes(x = entropy, y = mean_rho, color = area_type, fill = area_type)) +
  geom_point(data = bg_rho, aes(size = population), shape = 16, alpha = 0.72, show.legend = TRUE) +
  geom_point(
    data = label_rho, aes(size = population),
    shape = 21, color = "grey15", stroke = 0.85, alpha = 1
  ) +
  geom_label_repel(
    data = label_rho,
    aes(label = label_city, color = area_type),
    size = 3.0, fontface = "bold",
    fill = alpha("white", 0.9), label.size = 0,
    label.padding = unit(0.12, "lines"),
    point.size = 12,
    nudge_x = label_rho$.nudge_x,
    nudge_y = label_rho$.nudge_y,
    min.segment.length = 0,
    segment.color = "grey50", segment.size = 0.35,
    box.padding = 0.35, point.padding = 0.5,
    force = 1.5, force_pull = 0,
    max.overlaps = Inf, max.iter = 100000, max.time = 5,
    seed = 42, show.legend = FALSE
  ) +
  scale_color_manual(values = rho_entropy_colors, name = NULL, drop = FALSE) +
  scale_fill_manual(values = rho_entropy_colors, name = NULL, drop = FALSE) +
  scale_size_area(
    name = "Population", max_size = 11,
    breaks = c(1e6, 3e6, 6e6), labels = c("1M", "3M", "6M")
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.06, 0.10))) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0.04, 0.14))) +
  coord_cartesian(clip = "off") +
  labs(x = "Entropy of visit distribution", y = expression(Mean ~ local ~ Spearman ~ rho)) +
  guides(
    color = guide_legend(order = 1, override.aes = list(size = 3.5, alpha = 1, shape = 16)),
    fill = "none",
    size = guide_legend(order = 2, title = "Population")
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.box = "vertical",
    legend.background = element_rect(fill = alpha("white", 0.85), color = NA),
    legend.key = element_rect(fill = NA, color = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey88", linewidth = 0.4),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

# ---------------------------------------------------------------------------
# NYC rank-mismatch inset maps
# ---------------------------------------------------------------------------
n_rank <- nrow(focal_nbrs)
rank_fill <- scale_fill_gradient(
  name = "Rank", low = "#1A1A1A", high = "#F0F0F0",
  breaks = pretty(c(1, n_rank), n = 4),
  limits = c(1, n_rank), oob = squish,
  guide = guide_colorbar(barwidth = 0.4, barheight = 3.5)
)
map_theme <- theme_void(base_size = 9) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    plot.margin = margin(2, 4, 2, 2)
  )

p_geo_map <- ggplot() +
  geom_sf(
    data = filter(map_df, role == "neighbor"),
    aes(fill = geo_rank_plot), color = "white", linewidth = 0.05
  ) +
  geom_sf(
    data = filter(map_df, role == "focal"),
    fill = HIGHLIGHT_FILL, color = HIGHLIGHT_COLOR, linewidth = 0.51
  ) +
  rank_fill + map_theme

p_net_map <- ggplot() +
  geom_sf(
    data = filter(map_df, role == "neighbor"),
    aes(fill = net_rank_plot), color = "white", linewidth = 0.05
  ) +
  geom_sf(
    data = filter(map_df, role == "focal"),
    fill = HIGHLIGHT_FILL, color = HIGHLIGHT_COLOR, linewidth = 0.55
  ) +
  rank_fill + map_theme

p_inset <- (p_geo_map | p_net_map) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

# ---------------------------------------------------------------------------
# Geo vs network hex scatters (NYC, Mexico City, London)
# ---------------------------------------------------------------------------
nyc_df <- nyc_pairs |>
  filter(is.finite(d_geo_km), d_geo_km > 0, is.finite(d_net), d_net > 0)
mexico_df <- mexico_pairs |>
  filter(is.finite(d_geo_km), d_geo_km > 0, is.finite(d_net), d_net > 0)
london_df <- london_pairs |>
  filter(is.finite(d_geo_km), d_geo_km > 0, is.finite(d_net), d_net > 0)
highlight_df <- focal_nbrs |>
  filter(is.finite(d_geo_km), d_geo_km > 0, is.finite(d_net), d_net > 0)

p_a_nyc <- ggplot(nyc_df, aes(x = d_geo_km, y = d_net)) +
  geom_bin2d(bins = 45) +
  geom_point(data = highlight_df, aes(x = d_geo_km, y = d_net), color = HIGHLIGHT_COLOR, shape = 16, size = 0.5, alpha = 0.75) +
  scale_fill_gradient(low = "#DEEBF7", high = CITY_FILL[["nyc"]], trans = "log10", guide = "none") +
  scale_x_log10() + scale_y_log10() +
  labs(x = "Geographic distance (km)", y = "Network distance") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "none")

p_a_mexico <- ggplot(mexico_df, aes(x = d_geo_km, y = d_net)) +
  geom_bin2d(bins = 45) +
  scale_fill_gradient(low = "#FEEDDE", high = CITY_FILL[["mexico"]], trans = "log10", guide = "none") +
  scale_x_log10() + scale_y_log10() +
  labs(x = "Geographic distance (km)", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "none")

legend_grob <- {
  g <- ggplotGrob(
    ggplot(data.frame(x = 1, y = 1, n = c(1, 10, 100, 1000)), aes(x, y, fill = n)) +
      geom_raster(show.legend = TRUE) +
      scale_fill_gradient(
        low = "#F0F0F0", high = "#4A4A4A", trans = "log10", name = NULL,
        breaks = c(1, 10, 100, 1000), labels = c("1", "10", "100", "1k")
      ) +
      theme_void(base_size = 11) +
      theme(
        legend.position = "right",
        legend.key.height = unit(0.42, "cm"),
        legend.key.width = unit(0.28, "cm"),
        legend.text = element_text(size = 8),
        legend.background = element_rect(fill = alpha("white", 0.9), color = NA)
      )
  )
  guide_i <- which(vapply(g$grobs, function(x) {
    identical(x$name, "guide-box") || startsWith(x$name, "guide-box")
  }, logical(1)))
  if (length(guide_i)) g$grobs[[guide_i[[1]]]] else g
}
p_a_mexico <- p_a_mexico +
  inset_element(
    legend_grob,
    left = 0.08, bottom = 0.58, right = 0.36, top = 0.98,
    align_to = "panel", on_top = TRUE
  )

p_a_london <- ggplot(london_df, aes(x = d_geo_km, y = d_net)) +
  geom_bin2d(bins = 45) +
  scale_fill_gradient(low = "#FEEDDE", high = CITY_FILL[["london"]], trans = "log10", guide = "none") +
  scale_x_log10() + scale_y_log10() +
  labs(x = "Geographic distance (km)", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "none")

p_a <- (p_a_nyc | p_a_mexico | p_a_london) + plot_layout(widths = c(1, 1, 1))

# ---------------------------------------------------------------------------
# Assemble and save
# ---------------------------------------------------------------------------
fig <- ((wrap_elements(full = p_inset) / p_a) + plot_layout(heights = c(0.95, 1.15))) | p_entropy

ggsave("output/fig_geo_vs_network_final.png", fig, width = 14, height = 8, dpi = 300, bg = "white")
tryCatch(
  ggsave("output/fig_geo_vs_network_final.pdf", fig, width = 14, height = 8, device = cairo_pdf, bg = "white"),
  error = function(e) {
    ggsave("output/fig_geo_vs_network_final.pdf", fig, width = 14, height = 8, bg = "white")
  }
)
message("Saved: output/fig_geo_vs_network_final.pdf")
