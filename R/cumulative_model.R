#' Fit a model to a species using the cumulative predictions for a secundary
#' species
#' @inheritParams base_model
#' @inheritParams load_relevant
#' @inheritParams probability_model
#' @export
#' @importFrom assertthat assert_that are_equal has_name is.flag is.string noNA
#' @importFrom dplyr arrange as_tibble bind_cols distinct inner_join mutate
#' pull select summarise %>%
#' @importFrom git2rdata read_vc
#' @importFrom INLA inla inla.mesh.1d inla.mesh.2d inla.spde2.pcmatern
#' inla.spde.make.A inla.spde.make.index inla.stack inla.stack.A inla.stack.data
#' inla.stack.index
#' @importFrom pROC roc
#' @importFrom rlang .data !!
#' @importFrom sf st_as_sf st_coordinates st_drop_geometry st_transform
#' @importFrom tidyr pivot_longer
#' @importFrom stats median
cumulative_model <- function(
  species = "Adal_bipu", min_occurrences = 1000, min_species = 3, secondary,
  knots = c(1990, 2000, 2010, 2020)
) {
  assert_that(is.string(species))
  if (missing(secondary)) {
    secondary <- base_model(
      species = "Harm_axyr", min_occurrences = min_occurrences,
      min_species = min_species, knots = knots
    )
  } else {
    assert_that(is.list(secondary))
    assert_that(has_name(secondary, "species"))
    assert_that(!are_equal(species, secondary$species))
    assert_that(has_name(secondary, "min_occurrences"))
    assert_that(are_equal(min_occurrences, secondary$min_occurrences))
    assert_that(has_name(secondary, "min_species"))
    assert_that(are_equal(min_species, secondary$min_species))
    assert_that(has_name(secondary, "predictions"))
  }
  read_vc("location", system.file(package = "ladybird")) %>%
    st_as_sf(coords = c("long", "lat"), crs = 4326) %>%
    st_transform(crs = 31370) -> base_loc
  base_loc %>%
    bind_cols(
      st_coordinates(base_loc) %>%
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
      .data$year, .data$location, .data$X, .data$Y, .data$visits,
      occurrence = !!species
    ) %>%
    filter(.data$year >= min(.data$year[.data$occurrence == 1]) - 1) %>%
    left_join(
      secondary$predictions %>%
        arrange(.data$year) %>%
        group_by(.data$location) %>%
        transmute(
          .data$year, .data$location, secondary = cumsum(.data$mean)
        ) %>%
        ungroup(),
      by = c("year", "location")
    ) %>%
    ladybird:::add_knots(knots = knots) %>%
    mutate(
      X = .data$X / 1e3, Y = .data$Y / 1e3,
      iyear = .data$year - min(.data$year) + 1,
      secondary = replace_na(.data$secondary, 0),
      visits = pmin(.data$visits, 30)
    ) -> base_data

  base_data %>%
    group_by(.data$year) %>%
    summarise(
      median = median(.data$secondary, na.rm = TRUE),
      min = min(.data$secondary, na.rm = TRUE),
      max = max(.data$secondary, na.rm = TRUE),
      everywhere = 1,
      without = 0
    ) %>%
    pivot_longer(-.data$year, names_to = "type", values_to = "secondary") %>%
    mutate(
      visits = 1
    ) %>%
    arrange(.data$year, .data$type) %>%
    add_knots(knots = knots) -> trend_prediction_base
  trend_prediction_base %>%
    mutate(iyear = .data$year - min(.data$year) + 1) %>%
    bind_rows(trend_prediction_base) -> trend_prediction

  base_data %>%
    distinct(.data$location, .data$year) %>%
    complete(.data$location, .data$year) %>%
    inner_join(
      base_data %>%
        distinct(.data$location, .data$X, .data$Y),
      by = "location"
    ) %>%
    add_knots(knots = knots) %>%
    inner_join(
      secondary$predictions %>%
        select(.data$year, .data$location, secondary = .data$mean),
      by = c("year", "location")
    ) %>%
    mutate(
      iyear = .data$year - min(.data$year) + 1, visits = 1
    ) %>%
    arrange(.data$location, .data$year) -> base_prediction

  expand.grid(
    X = pretty(base_data$X, 20),
    Y = pretty(base_data$Y, 20),
    year = knots,
    visits = 1
  ) -> field_prediction

  results <- fit_model(
    base_data = base_data, trend_prediction = trend_prediction, knots = knots,
    base_prediction = base_prediction, field_prediction = field_prediction
  )

  return(
    c(
      species = species, secondary = secondary$species,
      min_occurrences = min_occurrences, min_species = min_species,
      results, type = "cumulative"
    )
  )
}
