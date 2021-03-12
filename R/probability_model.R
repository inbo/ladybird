#' Fit a model to a species using the predictions for a secondary species
#' @inheritParams base_model
#' @inheritParams load_relevant
#' @param secondary The output of `base_model()` for a different species.
#' @export
#' @importFrom assertthat assert_that are_equal has_name is.flag is.string noNA
#' @importFrom dplyr arrange as_tibble bind_cols distinct inner_join mutate
#' select summarise %>%
#' @importFrom git2rdata read_vc
#' @importFrom rlang .data !!
#' @importFrom sf st_as_sf st_coordinates st_drop_geometry st_transform
#' @importFrom tidyr pivot_longer
#' @importFrom stats median
probability_model <- function(
  species = "Adal_dece", min_occurrences = 1000, min_species = 3, secondary,
  first_order = TRUE, center_year = 2001
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
      iyear = .data$year - min(.data$year) + 1,
      iyear2 = .data$iyear,
      before = pmin(.data$year - center_year, 0),
      after = pmax(.data$year - center_year, 0)
    ) %>%
    inner_join(
      secondary$predictions %>%
        select(.data$year, .data$location, secondary = .data$mean),
      by = c("year", "location")
    ) -> base_data
  base_data %>%
    group_by(.data$year) %>%
    summarise(
      median = median(.data$secondary, na.rm = TRUE),
      min = min(.data$secondary, na.rm = TRUE),
      max = max(.data$secondary, na.rm = TRUE),
      without = 0
    ) %>%
    pivot_longer(-.data$year, names_to = "type", values_to = "secondary") %>%
    mutate(
      intercept = 1,
      iyear = .data$year - min(.data$year) + 1,
      iyear2 = .data$iyear,
      before = pmin(.data$year - center_year, 0),
      after = pmax(.data$year - center_year, 0)
    ) %>%
    arrange(.data$year, .data$type) -> trend_prediction
  base_data %>%
    distinct(.data$location, .data$year) %>%
    complete(.data$location, .data$year) %>%
    inner_join(
      base_data %>%
        distinct(.data$location, .data$X, .data$Y),
      by = "location"
    ) %>%
    inner_join(
      secondary$predictions %>%
        select(.data$year, .data$location, secondary = .data$mean),
      by = c("year", "location")
    ) %>%
    mutate(
      iyear = .data$year - min(.data$year) + 1,
      iyear2 = .data$iyear,
      before = pmin(.data$year - center_year, 0),
      after = pmax(.data$year - center_year, 0)
    ) %>%
    arrange(.data$location, .data$year) -> base_prediction
  results <- fit_model(
    first_order = first_order, base_data = base_data,
    trend_prediction = trend_prediction, base_prediction = base_prediction
  )
  return(
    c(
      species = species, secondary = secondary$species,
      min_occurrences = min_occurrences, min_species = min_species,
      results, type = "proportion"
    )
  )
}
