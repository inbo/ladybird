#' Fit a base model to a species
#' @param species Name of the species.
#' @inheritParams load_relevant
#' @param first_order Use first (`TRUE`) or second (`FALSE`) order random walk
#' for the year component.
#' Defaults to `TRUE`.
#' @param center_year The year to center to.
#' Defaults to `2001`.
#' @export
#' @importFrom assertthat assert_that is.flag noNA
#' @importFrom dplyr arrange bind_cols distinct inner_join mutate select %>%
#' @importFrom git2rdata read_vc
#' @importFrom rlang .data !!
#' @importFrom sf st_as_sf st_coordinates st_drop_geometry st_transform
#' @importFrom tidyr complete
base_model <- function(
  species = "Harm_axyr", min_occurrences = 1000, min_species = 3,
  first_order = TRUE, center_year = 2001
) {
  assert_that(is.flag(first_order), noNA(first_order))
  read_vc("location", system.file(package = "ladybird")) %>%
    st_as_sf(coords = c("long", "lat"), crs = 4326) %>%
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
      iyear = .data$year - min(.data$year) + 1,
      iyear2 = .data$iyear,
      cyear = .data$year - center_year,
      X = .data$X / 1e3, Y = .data$Y / 1e3,
      secondary = NA_real_
    ) -> base_data
  base_data %>%
    distinct(.data$year) %>%
    arrange(.data$year) %>%
    mutate(
      iyear = .data$year - min(.data$year) + 1,
      iyear2 = .data$iyear,
      cyear = .data$year - center_year,
      intercept = 1,
      secondary = NA_real_
    ) -> trend_prediction
  base_data %>%
    distinct(.data$location, .data$year) %>%
    complete(.data$location, .data$year) %>%
    inner_join(
      base_data %>%
        distinct(.data$location, .data$X, .data$Y),
      by = "location"
    ) %>%
    mutate(
      iyear = .data$year - min(.data$year) + 1,
      iyear2 = .data$iyear,
      cyear = .data$year - center_year,
      secondary = NA
    ) %>%
    arrange(.data$location, .data$year) -> base_prediction
  results <- fit_model(
    first_order = first_order, base_data = base_data,
    trend_prediction = trend_prediction, base_prediction = base_prediction
  )
  return(
    c(
      species = species, min_occurrences = min_occurrences,
      min_species = min_species, results, type = "base"
    )
  )
}
