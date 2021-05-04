#' Fit a base model to a species
#' @inheritParams load_relevant
#' @inheritParams base_model
#' @inheritParams get_country_grid
#' @export
#' @importFrom assertthat assert_that is.flag noNA
#' @importFrom dplyr arrange bind_cols bind_rows distinct inner_join mutate
#' select starts_with row_number %>%
#' @importFrom git2rdata read_vc
#' @importFrom magrittr %<>%
#' @importFrom rlang .data !!
#' @importFrom sf as_Spatial st_area st_as_sf st_buffer st_coordinates st_crs
#' st_distance st_drop_geometry st_intersection st_transform st_union
#' @importFrom stats setNames
#' @importFrom tidyr complete extract
#' @importFrom units as_units
#' @importFrom utils head
smooth_model <- function(
  species = "Harm_axyr", min_occurrences = 1000, min_species = 3,
  country = c("BE", "GB", "NL", "all"), path = ".", cellsize = 10e3,
  buffer_distance = 50e3, buffer_locations = 50
) {
  which_country <- match.arg(country)
  crs <- c(BE = 31370, GB = 27700, NL = 28992, all = 3035)
  base_loc <- read_vc("location", system.file(package = "ladybird"))
  if (which_country != "all") {
    base_loc %<>%
      filter(.data$country == which_country)
  }
  base_loc %>%
    st_as_sf(coords = c("long", "lat"), crs = 4326) %>%
    st_transform(crs = crs[which_country]) -> base_loc
  base_loc %>%
    bind_cols(
      st_coordinates(base_loc) %>%
        as.data.frame()
    ) %>%
    inner_join(
      load_relevant(
        min_occurrences = min_occurrences, min_species = min_species,
        country = which_country, buffer_distance = buffer_distance,
        buffer_locations = buffer_locations
      ),
      by = "location"
    ) %>%
    select(
      .data$year, .data$location, .data$X, .data$Y, .data$visits,
      .data$buffer_count, occurrence = !!species
    ) %>%
    filter(.data$year >= min(.data$year[.data$occurrence == 1]) - 1) %>%
    mutate(
      log_visits = log(.data$visits),
      iyear = .data$year - min(.data$year) + 1
    ) -> base_data
  base_data %>%
    filter(.data$buffer_count >= buffer_locations) %>%
    st_buffer(buffer_distance) %>%
    st_union() -> base_buffer
  base_buffer %>%
    st_buffer(buffer_distance) %>%
    st_buffer(-buffer_distance) %>%
    as_Spatial() -> mesh_buffer
  time_mesh <- inla.mesh.1d(loc = pretty(base_data$year, 3))
  inla.mesh.2d(
    boundary = mesh_buffer, max.edge = cellsize, cutoff = cellsize / 2
  ) -> mesh
  spde <- inla.spde2.pcmatern(
    mesh = mesh, prior.range = c(50e3, 0.5), prior.sigma = c(0.1, 0.05)
  )
  base_data %<>%
    st_drop_geometry()
  base_data %>%
    select(.data$X, .data$Y) %>%
    as.matrix() %>%
    inla.spde.make.A(
      mesh = mesh, group = base_data$year, mesh.group = time_mesh
    ) -> a_estimate
  site_index <- inla.spde.make.index(
    name = "site", n.spde = mesh$n, n.group = time_mesh$n
  )
  inla.stack(
    data = select(base_data, .data$occurrence),
    A = list(a_estimate, 1),
    effects = list(
      c(site_index, intercept = 1),
      list(
        base_data %>%
          select(.data$year, .data$iyear, .data$log_visits)
      )
    ),
    tag = "estimate"
  ) -> stack_estimate

  rbind(
    moving_trend(
      n_year = max(base_data$iyear), first_year = min(base_data$year),
      trend_year = 2
    ),
    moving_trend(
      n_year = max(base_data$iyear), first_year = min(base_data$year),
      trend_year = 5
    ),
    moving_trend(
      n_year = max(base_data$iyear), first_year = min(base_data$year),
      trend_year = 10
    )
  ) -> lc_iyear
  inla.make.lincombs(iyear = lc_iyear) %>%
    setNames(rownames(lc_iyear)) -> lc
  model_formula <- occurrence ~ 0 + intercept + log_visits +
    f(
      iyear, model = "rw2",
      hyper = list(theta = list(prior = "pc.prec", param = c(0.05, 0.01)))
    ) +
    f(
      site, model = spde, group = site.group,
      control.group = list(
        model = "ar1",
        hyper = list(theta = list(prior = "pc.cor1", param = c(0.6, 0.7)))
      )
    )
  m0 <- inla(
    model_formula, family = "binomial", Ntrials = 1,
    data = inla.stack.data(stack_estimate),
    control.predictor = list(A = inla.stack.A(stack_estimate), compute = FALSE)
  )
  fixed <- m0$summary.fixed
  hyperpar <- m0$summary.hyperpar

  base_data %>%
    distinct(.data$year, .data$iyear) %>%
    arrange(.data$year) %>%
    mutate(
      intercept = 1,
      log_visits = log(4)
    ) -> trend_prediction
  inla.stack(
    data = data.frame(occurrence = NA),
    A = list(1),
    effects = list(list(trend_prediction)),
    tag = "trend"
  ) -> stack_trend
  inla.stack(stack_estimate, stack_trend) -> stack
  m1 <- inla(
    model_formula, family = "binomial", Ntrials = 1, lincomb = lc,
    data = inla.stack.data(stack), control.compute = list(waic = TRUE),
    control.predictor = list(A = inla.stack.A(stack), compute = TRUE, link = 1),
    control.mode = list(theta = m0$mode$theta, restart = FALSE, fixed = TRUE)
  )
  index_auc <- inla.stack.index(stack, "estimate")$data
  roc_curve <- roc(
    predictor = m1$summary.fitted.values[index_auc, "mean"],
    response = base_data$occurrence,
    levels = c(0, 1), direction = "<", ci = TRUE
  )
  index_trend <- inla.stack.index(stack, "trend")$data
  m1$summary.fitted.values[index_trend, ] %>%
    select(.data$mean, median = 4, lcl = 3, ucl = 5) %>%
    as_tibble() %>%
    bind_cols(trend_prediction) -> trend
  waic <- m1$waic$waic
  m1$summary.lincomb.derived %>%
    rownames_to_column(var = "parameter") %>%
    extract(
      .data$parameter, c("center", "duration"), "trend_(.*)_(.*)",
      convert = TRUE
    ) %>%
    transmute(
      .data$center, .data$duration, median = .data$`0.5quant`,
      lcl = .data$`0.025quant`, ucl = .data$`0.975quant`
    ) -> lincomb
  rm(m1)
  gc()

  get_country_grid(
    path = path, country = country, cellsize = cellsize, what = "centers"
  ) %>%
    st_intersection(base_buffer) %>%
    st_coordinates() %>%
    as.data.frame() %>%
    mutate(
      location = row_number()
    ) -> prediction_grid
  expand.grid(
    year = time_mesh$loc,
    location = prediction_grid$location
  ) %>%
    inner_join(prediction_grid, by = "location") -> field_grid
  field_grid %>%
    select(.data$X, .data$Y) %>%
    as.matrix() %>%
    inla.spde.make.A(
      mesh = mesh, group = field_grid$year, mesh.group = time_mesh
    ) -> a_field
  inla.stack(
    data = data.frame(occurrence = NA),
    A = list(a_field),
    effects = list(
      c(site_index, list(intercept = 0))
    ),
    tag = "field"
  ) -> stack_field
  expand.grid(
    year = sort(unique(base_data$year)),
    location = prediction_grid$location
  ) %>%
    inner_join(prediction_grid, by = "location") %>%
    mutate(
      iyear = .data$year - min(.data$year) + 1,
      log_visits = log(4)
    ) -> prediction_grid
  prediction_grid %>%
    select(.data$X, .data$Y) %>%
    as.matrix() %>%
    inla.spde.make.A(
      mesh = mesh, group = prediction_grid$year, mesh.group = time_mesh
    ) -> a_prediction
  inla.stack(
    data = data.frame(occurrence = NA),
    A = list(a_prediction, 1),
    effects = list(c(site_index, intercept = 1), list(prediction_grid)),
    tag = "prediction"
  ) -> stack_prediction
  inla.stack(stack_estimate, stack_prediction, stack_field) -> stack2
  m2 <- inla(
    model_formula, family = "binomial", Ntrials = 1,
    data = inla.stack.data(stack2),
    control.predictor = list(
      A = inla.stack.A(stack2), compute = TRUE, link = 1
    ),
    control.mode = list(theta = m0$mode$theta, restart = FALSE, fixed = TRUE)
  )
  prediction_field <- inla.stack.index(stack2, "prediction")$data
  m2$summary.fitted.values[prediction_field, ] %>%
    select(median = 4, lcl = 3, ucl = 5) %>%
    as_tibble() %>%
    bind_cols(prediction_grid) -> prediction
  index_field <- inla.stack.index(stack2, "field")$data
  m2$summary.linear.predictor[index_field, ] %>%
    select(median = 4, lcl = 3, ucl = 5) %>%
    as_tibble() %>%
    bind_cols(field_grid) -> field
  rm(m0, m2)
  gc()

  return(
    list(
      fixed = fixed, trend = trend, field = field, hyperpar = hyperpar,
      waic = waic, roc = roc_curve, lincomb = lincomb, species = species,
      country = country, min_occurrences = min_occurrences,
      min_species = min_species, prediction = prediction
    )
  )
}

moving_trend <- function(n_year = 21, trend_year = 5, first_year = 1991) {
  trend_coef <- seq_len(trend_year) - (trend_year + 1) / 2
  trend_coef <- trend_coef / sum(trend_coef ^ 2)
  lc <- vapply(
    seq_len(n_year - trend_year + 1),
    function(i) {
      c(rep(0, i - 1), trend_coef, rep(0, n_year - trend_year - i + 1))
    },
    numeric(n_year)
  )
  colnames(lc) <- sprintf(
    "trend_%.1f_%i",
    seq_len(ncol(lc)) + first_year - 1 + (trend_year - 1) / 2, trend_year
  )
  t(lc)
}

#' Get a grid of points within the country
#'
#' We use the borders from [Natural Earth
#' Data](https://www.naturalearthdata.com/)
#' @inheritParams base_model
#' @param path Location of the Natural Earth Data.
#'   The function will download and unzip the data when not available.
#' @param cellsize Size of the grid cells in meter.
#'    Must be at least 1000 m.
#'    Defaults to 10.000 m.
#' @inheritParams sf::st_make_grid
#' @export
#' @importFrom assertthat assert_that is.number
#' @importFrom curl curl_download
#' @importFrom sf read_sf st_cast st_coordinates st_intersects st_make_grid
#' st_transform st_union
#' @importFrom dplyr filter %>%
#' @importFrom utils unzip
get_country_grid <- function(
  path = ".", country = c("BE", "GB", "NL", "all"), cellsize = 10e3,
  what = "polygons"
) {
  which_country <- match.arg(country)
  path <- normalizePath(path)
  assert_that(is.number(cellsize))
  assert_that(cellsize >= 1e3)
  assert_that(is.string(what))
  if (!file.exists(file.path(path, "ne_10m_admin_0_map_units.shp"))) {
    if (!file.exists(file.path(path, "ne_10m_admin_0_map_untis.zip"))) {
      file.path(
        "https://www.naturalearthdata.com", "http", "",
        "www.naturalearthdata.com", "download", "10m", "cultural",
        "ne_10m_admin_0_map_units.zip", fsep = "/"
      ) %>%
        curl_download(file.path(path, "ne_10m_admin_0_map_units.zip"))
    }
    file.path(path, "ne_10m_admin_0_map_units.zip") %>%
      unzip(exdir = path)
    assert_that(file.exists(file.path(path, "ne_10m_admin_0_map_units.shp")))
  }
  crs <- c(BE = 31370, NL = 28992, GB = 27700, all = 3035)
  country_name <- c(BE = "Belgium", NL = "Netherlands", GB = "United Kingdom")
  file.path(path, "ne_10m_admin_0_map_units.shp") %>%
    read_sf() %>%
    filter(
      .data$ADMIN == country_name[which_country], .data$scalerank == 0,
      .data$GEOUNIT != "Northern Ireland"
    ) %>%
    st_transform(crs = crs[which_country]) %>%
    st_union() %>%
    st_cast("POLYGON") -> borders
  if (what == "borders") {
    return(borders)
  }
  hex_grid <- st_make_grid(
    x = borders, cellsize = cellsize, what = what, square = FALSE,
    flat_topped = TRUE
  )
  hex_grid %>%
    st_intersects(borders, sparse = FALSE) %>%
    apply(1, any) -> relevant
  hex_grid[relevant, ]
}
