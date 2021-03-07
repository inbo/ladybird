#' Fit a model to a species using the predictions for a secundary species
#' @param species Name of the species
#' @inheritParams load_relevant
#' @param secondary the output of `base_model()` for a different species
#' @export
#' @importFrom assertthat assert_that has_name is.string
#' @importFrom dplyr arrange as_tibble bind_cols distinct inner_join mutate
#' select %>%
#' @importFrom git2rdata read_vc
#' @importFrom INLA inla inla.mesh.1d inla.mesh.2d inla.spde2.pcmatern
#' inla.spde.make.A inla.spde.make.index inla.stack inla.stack.A inla.stack.data
#' inla.stack.index
#' @importFrom rlang .data !!
#' @importFrom sf st_as_sf st_coordinates st_drop_geometry st_transform
secondary_model <- function(
  species = "Adal_dece", min_occurrences = 1000, min_species = 3, secondary
) {
  assert_that(is.string(species))
  if (missing(secondary)) {
    secondary <- base_model(
      species = "Harm_axyr", min_occurrences = min_occurrences,
      min_species = min_species
    )
  } else {
    assert_that(is.list(secondary))
    assert_that(has_name(secondary, "species"))
    assert_that(species != secondary$species)
    assert_that(has_name(secondary, "min_occurrences"))
    assert_that(min_occurrences != secondary$min_occurrences)
    assert_that(has_name(secondary, "min_species"))
    assert_that(min_species != secondary$min_species)
    assert_that(has_name(secondary, "predictions"))
  }
  read_vc("location", system.file(package = "ladybird")) %>%
    st_as_sf(coords = c("long", "lat"), crs = 4326) %>%
    st_transform(crs = 31370) %>%
    bind_cols(
      st_coordinates(.) %>%
        `/`(1e3) %>%
        as.data.frame()
    ) %>%
    st_drop_geometry() %>%
    inner_join(
      load_relevant(
        min_occurrences = min_occurrences, min_species = min_species
      ),
      by = "location"
    ) %>%
    select(
      .data$year, .data$location, .data$X, .data$Y, occurrence = !!species
    ) %>%
    mutate(cyear = .data$year - min(.data$year) + 1) %>%
    inner_join(
      secondary$predictions %>%
        select(.data$year, .data$location, secondary = .data$mean),
      by = c("year", "location")
    ) -> base_data
  n_year <- max(base_data$cyear)
  year_interval <- 5
  seq(year_interval / 2, n_year, by = year_interval) %>%
    inla.mesh.1d() -> time_mesh
  inla.mesh.2d(
    loc.domain = distinct(base_data, .data$X, .data$Y),
    max.edge = 10
  ) -> mesh
  spde <- inla.spde2.pcmatern(
    mesh = mesh, prior.range = c(50, 0.5), prior.sigma = c(0.1, 0.05)
  )
  base_data %>%
    select(.data$X, .data$Y) %>%
    as.matrix() %>%
    inla.spde.make.A(
      mesh = mesh, group = base_data$cyear, mesh.group = time_mesh
    ) -> a_estimate
  site_index <- inla.spde.make.index(
    name = "site", n.spde = mesh$n, n.group = time_mesh$n
  )
  inla.stack(
    data = select(base_data, .data$occurrence),
    A = list(a_estimate, 1),
    effects = list(
      c(site_index, list(intercept = 1)),
      list(select(base_data, .data$cyear, .data$secondary))
    ),
    tag = "estimate"
  ) -> stack_estimate
  base_data %>%
    group_by(.data$year) %>%
    summarise(
      median = quantile(.data$secondary, probs = 0.5),
      lcl = quantile(.data$secondary, probs = 0.025),
      ucl = quantile(.data$secondary, probs = 0.975)
    ) %>%
    pivot_longer(-.data$year, names_to = "type", values_to = "secondary") %>%
    mutate(cyear = .data$year - min(.data$year) + 1) %>%
    arrange(.data$cyear, .data$type) -> trend_prediction
  inla.stack(
    data = data.frame(occurrence = NA),
    A = list(1),
    effects = list(list(trend_prediction)),
    tag = "trend"
  ) -> stack_trend
  inla.stack(stack_estimate, stack_trend) -> stack
  model <- inla(
    occurrence ~ 0 + intercept + secondary +
      f(
        cyear, model = "rw1",
        hyper = list(theta = list(prior = "pc.prec", param = c(0.1, 0.05)))
      ) +
      f(
        site, model = spde, group = site.group,
        control.group = list(
          model = "ar1",
          hyper = list(theta = list(prior = "pc.cor1", param = c(0.6, 0.7)))
        )
      ),
    family = "binomial",
    data = inla.stack.data(stack),
    control.compute = list(waic = TRUE),
    control.predictor = list(A = inla.stack.A(stack), compute = TRUE, link = 1)
  )
  index_estimate <- inla.stack.index(stack, "estimate")$data
  model$summary.fitted.values[index_estimate, ] %>%
    select(.data$mean, median = 4, lcl = 3, ucl = 5) %>%
    bind_cols(base_data) %>%
    as_tibble() %>%
    select(
      .data$year, .data$location, .data$mean, .data$median, .data$lcl, .data$ucl
    ) -> predictions
  index_trend <- inla.stack.index(stack, "trend")$data
  model$summary.fitted.values[index_trend, ] %>%
    select(.data$mean, median = 4, lcl = 3, ucl = 5) %>%
    as_tibble() %>%
    bind_cols(trend_prediction) -> trend
  return(
    list(
      species = species, secondary = secondary$species,
      min_occurrences = min_occurrences, min_species = min_species,
      fixed = model$summary.fixed, trend = trend,
      hyperpar = model$summary.hyperpar, predictions = predictions,
      waic = model$waic$waic
    )
  )
}
