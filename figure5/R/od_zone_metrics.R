# Metro metrics from an OD network + zone centroids.
# Used by precompute_metro_metrics.R (and later make_figure5.R).

#' Mean Spearman ρ(d_geo, d_net) over origins, among neighbors within dist_km.
mean_local_rho_by_distance <- function(sym_pairs, dist_km, min_neighbors = 3L) {
  local <- sym_pairs |>
    dplyr::filter(d_geo_km <= dist_km) |>
    dplyr::group_by(geoid_i) |>
    dplyr::summarise(
      n_local = dplyr::n(),
      rho = suppressWarnings(stats::cor(d_geo_km, d_net, method = "spearman")),
      .groups = "drop"
    ) |>
    dplyr::filter(is.finite(rho), n_local >= min_neighbors)
  if (nrow(local) == 0L) NA_real_ else mean(local$rho)
}

#' Mean Spearman ρ(d_geo, d_net) over origins, among the k nearest neighbors by d_geo.
mean_local_rho_by_k <- function(sym_pairs, k_local, min_neighbors = 3L) {
  local <- sym_pairs |>
    dplyr::group_by(geoid_i) |>
    dplyr::mutate(geo_rank = rank(d_geo_km, ties.method = "average")) |>
    dplyr::filter(geo_rank <= k_local) |>
    dplyr::summarise(
      n_local = dplyr::n(),
      rho = suppressWarnings(stats::cor(d_geo_km, d_net, method = "spearman")),
      .groups = "drop"
    ) |>
    dplyr::filter(is.finite(rho), n_local >= min_neighbors)
  if (nrow(local) == 0L) NA_real_ else mean(local$rho)
}

#' Diffusion-kernel network distance, entropy, and local Spearman ρ for one metro.
#'
#' @param od data frame with home_zone, visit_zone, count
#' @param centroids_sf sf POINTS with column GEOID
#' @param return_pairs if TRUE, also return the undirected pairs table (for scatters)
compute_city_error_od_zones <- function(
    od,
    centroids_sf,
    population = NA_real_,
    area_name = "city",
    cbsa = NULL,
    min_zones = 40L,
    beta = 0.8,
    K = 5L,
    min_pairs_tract = 30L,
    return_pairs = FALSE) {
  population <- as.numeric(population)[1]
  if (length(population) == 0L) population <- NA_real_
  area_name <- as.character(area_name)[1]
  if (is.null(cbsa) || length(cbsa) == 0L || is.na(cbsa)) {
    cbsa <- tolower(gsub("[^a-zA-Z0-9]+", "-", area_name))
  }
  cbsa <- as.character(cbsa)[1]

  geoids <- as.character(centroids_sf$GEOID)
  n <- length(geoids)
  if (n <= min_zones) {
    return(NULL)
  }

  flows_zone <- od |>
    dplyr::filter(home_zone %in% geoids, visit_zone %in% geoids, count > 0) |>
    dplyr::group_by(home_zone, visit_zone) |>
    dplyr::summarise(count = sum(as.numeric(count), na.rm = TRUE), .groups = "drop")

  if (nrow(flows_zone) < min_pairs_tract) {
    return(NULL)
  }

  entropy_val <- flows_zone |>
    dplyr::group_by(visit_zone) |>
    dplyr::summarise(visits = sum(count), .groups = "drop") |>
    dplyr::mutate(p = visits / sum(visits)) |>
    dplyr::summarise(entropy = -sum(p * log(p)) / log(dplyr::n())) |>
    dplyr::pull(entropy)

  W <- matrix(0, n, n, dimnames = list(geoids, geoids))
  W[cbind(
    match(flows_zone$home_zone, geoids),
    match(flows_zone$visit_zone, geoids)
  )] <- flows_zone$count
  rs <- rowSums(W)
  nz <- rs > 0
  W[nz, ] <- W[nz, ] / rs[nz]

  W_diff <- matrix(0, n, n)
  Wk <- diag(n)
  for (k in seq_len(K)) {
    Wk <- Wk %*% W
    W_diff <- W_diff + (beta^k) * Wk
  }

  geo_dist_m <- matrix(
    as.numeric(sf::st_distance(centroids_sf, centroids_sf)),
    n, n,
    dimnames = list(geoids, geoids)
  )
  conn <- (W_diff + t(W_diff)) / 2
  d_net <- ifelse(conn > 0, 1 / conn, NA_real_)
  diag(d_net) <- NA_real_

  pair_idx <- which(upper.tri(geo_dist_m), arr.ind = TRUE)
  pairs <- tibble::tibble(
    geoid_i = geoids[pair_idx[, 1]],
    geoid_j = geoids[pair_idx[, 2]],
    d_geo_km = geo_dist_m[cbind(pair_idx[, 1], pair_idx[, 2])] / 1000,
    d_net = d_net[cbind(pair_idx[, 1], pair_idx[, 2])]
  ) |>
    dplyr::filter(
      is.finite(d_geo_km), d_geo_km > 0,
      is.finite(d_net), d_net > 0
    )

  sym_pairs <- dplyr::bind_rows(
    pairs,
    dplyr::rename(pairs, geoid_i = geoid_j, geoid_j = geoid_i)
  )

  # nearest-neighbor spacing from the geo distance matrix (no nngeo needed)
  diag(geo_dist_m) <- Inf
  mean_dist_nn <- mean(apply(geo_dist_m, 1, min)) / 1000

  msa <- tibble::tibble(
    cbsa = cbsa,
    msa_name = area_name,
    n_tracts = as.integer(n),
    population = population,
    entropy = as.numeric(entropy_val)[1],
    mean_rho_local = mean_local_rho_by_k(sym_pairs, k_local = 15L),
    mean_rho_local_d5 = mean_local_rho_by_distance(sym_pairs, dist_km = 5),
    mean_rho_local_d10 = mean_local_rho_by_distance(sym_pairs, dist_km = 10),
    mean_dist_nn = mean_dist_nn
  )

  out <- list(msa = msa)
  if (return_pairs) {
    out$pairs = pairs
  }
  out
}
