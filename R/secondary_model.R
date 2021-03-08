#' Fit a model to a species using the predictions for a secundary species
#' @inheritParams base_model
#' @inheritParams load_relevant
#' @param secondary The output of `base_model()` for a different species.
#' @export
#' @importFrom assertthat assert_that are_equal has_name is.flag is.string noNA
#' @importFrom dplyr arrange as_tibble bind_cols distinct inner_join mutate
#' select summarise %>%
#' @importFrom git2rdata read_vc
#' @importFrom INLA inla inla.mesh.1d inla.mesh.2d inla.spde2.pcmatern
#' inla.spde.make.A inla.spde.make.index inla.stack inla.stack.A inla.stack.data
#' inla.stack.index
#' @importFrom rlang .data !!
#' @importFrom sf st_as_sf st_coordinates st_drop_geometry st_transform
#' @importFrom tidyr pivot_longer
#' @importFrom stats median
secondary_model <- function(
  species = "Adal_dece", min_occurrences = 1000, min_species = 3, secondary,
  first_order = TRUE
) {
  assert_that(is.string(species))
  assert_that(is.flag(first_order), noNA(first_order))
  if (missing(secondary)) {
    secondary <- base_model(
      species = "Harm_axyr", min_occurrences = min_occurrences,
      min_species = min_species, first_order = first_order
    )
  } else {
    assert_that(is.list(secondary))
    assert_that(has_name(secondary, "species"))
    assert_that(!are_equal(species, secondary$species))
    assert_that(has_name(secondary, "min_occurrences"))
    assert_that(are_equal(min_occurrences, secondary$min_occurrences))
    assert_that(has_name(secondary, "min_species"))
    assert_that(are_equal(min_species, secondary$min_species))
    assert_that(has_name(secondary, "first_order"))
    assert_that(are_equal(first_order, secondary$first_order))
    assert_that(has_name(secondary, "predictions"))
  }
  read_vc("location", system.file(package = "ladybird")) %>%
    st_as_sf(coords = c("long", "lat"), crs = 4326) %>%
    st_transform(crs = 31370) %>%
    st_transform(crs = 31370) -> base_data
  base_data %>%
    bind_cols(
      st_coordinates(base_data) %>%
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
    mutate(
      X = .data$X / 1e3, Y = .data$Y / 1e3,
      cyear = .data$year - min(.data$year) + 1
    ) %>%
    inner_join(
      secondary$predictions %>%
        select(.data$year, .data$location, secondary = .data$lp_mean),
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
      median = median(.data$secondary, na.rm = TRUE),
      min = min(.data$secondary, na.rm = TRUE),
      max = max(.data$secondary, na.rm = TRUE), without = 0
    ) %>%
    pivot_longer(-.data$year, names_to = "type", values_to = "secondary") %>%
    mutate(
      intercept = 1,
      cyear = .data$year - min(.data$year) + 1
    ) %>%
    arrange(.data$cyear, .data$type) -> trend_prediction
  inla.stack(
    data = data.frame(occurrence = NA),
    A = list(1),
    effects = list(list(trend_prediction)),
    tag = "trend"
  ) -> stack_trend
  base_data %>%
    distinct(.data$location, .data$year) %>%
    complete(.data$location, .data$year) %>%
    inner_join(
      base_data %>%
        distinct(.data$location, .data$X, .data$Y),
      by = "location"
    ) %>%
    mutate(cyear = .data$year - min(.data$year) + 1) %>%
    arrange(.data$location, .data$year) -> base_prediction
  base_prediction %>%
    select(.data$X, .data$Y) %>%
    as.matrix() %>%
    inla.spde.make.A(
      mesh = mesh, group = base_prediction$cyear, mesh.group = time_mesh
    ) -> a_prediction
  inla.stack(
    data = data.frame(occurrence = NA),
    A = list(a_prediction, 1),
    effects = list(
      c(site_index, list(intercept = 1)),
      list(select(base_prediction, .data$cyear))
    ),
    tag = "prediction"
  ) -> stack_prediction
  inla.stack(stack_estimate, stack_trend) -> stack
  inla.stack(stack_estimate, stack_prediction) -> stack2
  fixed_formula <- "occurrence ~ 0 + intercept + secondary"
  rw_formula <- ifelse(
    first_order,
    "f(
      cyear, model = \"rw1\",
      hyper = list(theta = list(prior = \"pc.prec\", param = c(0.1, 0.05)))
    ) +",
    "f(
      cyear, model = \"rw2\",
      hyper = list(theta = list(prior = \"pc.prec\", param = c(0.02, 0.05)))
    )"
  )
  st_formula <- "f(
      site, model = spde, group = site.group,
      control.group = list(
        model = \"ar1\",
        hyper = list(theta = list(prior = \"pc.cor1\", param = c(0.6, 0.7)))
      )
    )"
  paste(fixed_formula, rw_formula, st_formula, sep = " +\n") %>%
    as.formula() -> model_formula
  m0 <- inla(
    model_formula, family = "binomial", data = inla.stack.data(stack_estimate),
    control.predictor = list(A = inla.stack.A(stack_estimate), compute = FALSE)
  )
  m1 <- inla(
    model_formula, family = "binomial", data = inla.stack.data(stack),
    control.compute = list(waic = TRUE),
    control.predictor = list(A = inla.stack.A(stack), compute = TRUE, link = 1),
    control.mode = list(theta = m0$mode$theta, restart = FALSE, fixed = TRUE)
  )
  m2 <- inla(
    model_formula, family = "binomial", data = inla.stack.data(stack2),
    control.predictor = list(
      A = inla.stack.A(stack2), compute = TRUE, link = 1
    ),
    control.mode = list(theta = m0$mode$theta, restart = FALSE, fixed = TRUE)
  )
  index_estimate <- inla.stack.index(stack2, "prediction")$data
  m2$summary.fitted.values[index_estimate, ] %>%
    select(.data$mean, median = 4, lcl = 3, ucl = 5) %>%
    as_tibble() %>%
    bind_cols(
      m2$summary.linear.predictor[index_estimate, ] %>%
        select(lp_mean = .data$mean, lp_median = 4, lp_lcl = 3, lp_ucl = 5),
      base_prediction %>%
        select(.data$year, .data$location)
    ) -> predictions
  index_trend <- inla.stack.index(stack, "trend")$data
  m1$summary.fitted.values[index_trend, ] %>%
    select(.data$mean, median = 4, lcl = 3, ucl = 5) %>%
    as_tibble() %>%
    bind_cols(trend_prediction) -> trend
  return(
    list(
      species = species, secondary = secondary$species,
      min_occurrences = min_occurrences, min_species = min_species,
      fixed = m0$summary.fixed, trend = trend,
      hyperpar = m0$summary.hyperpar, predictions = predictions,
      waic = m1$waic$waic, first_order = first_order
    )
  )
}
