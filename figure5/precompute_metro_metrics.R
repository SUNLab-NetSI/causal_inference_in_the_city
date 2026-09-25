# Precompute metro-level entropy and mean local Spearman ρ (within 5 km)
# from the OD networks in data/{cuebiq,lodes,uk,spain,latam}/.
#
# Run from the repo/ directory:
#   Rscript precompute_metro_metrics.R
#
# Expected layout for every dataset folder:
#   metro_index.csv
#   od/{metro_id}.parquet            # home_zone, visit_zone, count
#   centroids/{metro_id}.parquet     # GEOID, lon, lat
#   metrics/{metro_id}.parquet       # written by this script (resumable)
#
# Also writes a combined table: data/metro_metrics.csv  (+ .parquet)

suppressPackageStartupMessages({
  library(dplyr)
  library(arrow)
  library(sf)
})
sf::sf_use_s2(FALSE)

source("R/od_zone_metrics.R")

BETA <- 0.8
K_DIFF <- 5L
MIN_ZONES <- 40L

DATASETS <- tribble(
  ~folder,  ~area_type,
  "cuebiq", "US Mobility",
  "lodes",  "US Commuting",
  "uk",     "UK Commuting",
  "spain",  "Spain Mobility",
  "latam",  "LatAm Mobility"
)

centroids_from_parquet <- function(path) {
  tab <- read_parquet(path)
  st_as_sf(tab, coords = c("lon", "lat"), crs = 4326, remove = FALSE) |>
    transmute(GEOID = as.character(GEOID), geometry)
}

# Used once per metro across all datasets.
compute_metro_metrics <- function(od, centroids_sf, metro_id, msa_name,
                                  area_type, population) {
  zids <- intersect(
    unique(c(od$home_zone, od$visit_zone)),
    as.character(centroids_sf$GEOID)
  )
  centroids_sf <- centroids_sf |> filter(GEOID %in% zids)
  od <- od |> filter(home_zone %in% zids, visit_zone %in% zids)
  if (length(zids) <= MIN_ZONES) {
    return(NULL)
  }

  res <- compute_city_error_od_zones(
    od = od,
    centroids_sf = centroids_sf,
    population = population,
    area_name = msa_name,
    cbsa = metro_id,
    min_zones = MIN_ZONES,
    beta = BETA,
    K = K_DIFF,
    min_pairs_tract = 30L,
    return_pairs = FALSE
  )
  if (is.null(res)) {
    return(NULL)
  }

  tibble(
    area_type = area_type,
    metro_id = metro_id,
    msa_name = msa_name,
    population = as.numeric(population),
    n_zones = as.integer(res$msa$n_tracts),
    entropy = as.numeric(res$msa$entropy),
    mean_rho_local_d5 = as.numeric(res$msa$mean_rho_local_d5),
    mean_rho_local = as.numeric(res$msa$mean_rho_local),
    mean_rho_local_d10 = as.numeric(res$msa$mean_rho_local_d10),
    mean_dist_nn = as.numeric(res$msa$mean_dist_nn)
  )
}

for (d in seq_len(nrow(DATASETS))) {
  folder <- DATASETS$folder[[d]]
  area_type <- DATASETS$area_type[[d]]
  message("=== ", area_type, " (", folder, ") ===")

  metrics_dir <- file.path("data", folder, "metrics")
  dir.create(metrics_dir, recursive = TRUE, showWarnings = FALSE)

  idx <- read.csv(
    file.path("data", folder, "metro_index.csv"),
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(idx))) {
    mid <- as.character(idx$metro_id[[i]])
    od_path <- file.path("data", folder, "od", paste0(mid, ".parquet"))
    cent_path <- file.path("data", folder, "centroids", paste0(mid, ".parquet"))
    out_path <- file.path(metrics_dir, paste0(mid, ".parquet"))

    if (file.exists(out_path)) {
      if (i %% 20L == 0L || i == nrow(idx)) {
        message(sprintf("  [%d/%d] cached %s", i, nrow(idx), mid))
      }
      next
    }
    if (!file.exists(od_path) || !file.exists(cent_path)) {
      message("  missing files for ", mid)
      next
    }

    message(sprintf("  [%d/%d] compute %s", i, nrow(idx), mid))
    od <- read_parquet(od_path) |>
      transmute(
        home_zone = as.character(home_zone),
        visit_zone = as.character(visit_zone),
        count = as.numeric(count)
      ) |>
      filter(home_zone != visit_zone, is.finite(count), count > 0)
    cents <- centroids_from_parquet(cent_path)

    row <- tryCatch(
      compute_metro_metrics(
        od, cents, mid, idx$msa_name[[i]], area_type, idx$population[[i]]
      ),
      error = function(e) {
        message("    FAIL: ", conditionMessage(e))
        NULL
      }
    )
    if (!is.null(row)) {
      write_parquet(row, out_path)
    }
    invisible(gc())
  }
}

# --- assemble combined table from each dataset's metrics/ ---------------------
metric_files <- unlist(lapply(DATASETS$folder, function(folder) {
  list.files(
    file.path("data", folder, "metrics"),
    pattern = "\\.parquet$",
    full.names = TRUE
  )
}))
metro_metrics <- bind_rows(lapply(metric_files, read_parquet)) |>
  mutate(
    area_type = factor(
      area_type,
      levels = c(
        "US Mobility", "US Commuting", "UK Commuting",
        "Spain Mobility", "LatAm Mobility"
      )
    )
  ) |>
  arrange(area_type, desc(population))

write.csv(metro_metrics, "data/metro_metrics.csv", row.names = FALSE)
write_parquet(metro_metrics, "data/metro_metrics.parquet")

message("\nSaved data/metro_metrics.csv  (", nrow(metro_metrics), " metros)")
print(as.data.frame(count(metro_metrics, area_type)))
print(
  metro_metrics |>
    group_by(area_type) |>
    summarise(
      median_entropy = round(median(entropy, na.rm = TRUE), 3),
      median_rho_d5 = round(median(mean_rho_local_d5, na.rm = TRUE), 3),
      .groups = "drop"
    ) |>
    as.data.frame()
)
