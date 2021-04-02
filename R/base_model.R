#' Fit a base model to a species
#' @param species Name of the species.
#' @param country Data from which country to select.
#' @inheritParams load_relevant
#' @param knots Which years to use a knots for the piecewise linear regression.
#' @export
#' @importFrom assertthat assert_that is.flag noNA
#' @importFrom dplyr arrange bind_cols bind_rows distinct inner_join mutate
#' select starts_with %>%
#' @importFrom git2rdata read_vc
#' @importFrom rlang .data !!
#' @importFrom sf st_as_sf st_coordinates st_drop_geometry st_transform
#' @importFrom stats setNames
#' @importFrom tidyr complete
#' @importFrom utils head
base_model <- function(
  species = "Harm_axyr", min_occurrences = 1000, min_species = 3,
  knots = c(1990, 2000, 2010, 2020), country = c("BE", "NL")
) {
  which_country <- match.arg(country)
  crs <- c(BE = 31370, NL = 28992)
  read_vc("location", system.file(package = "ladybird")) %>%
    filter(.data$country == which_country) %>%
    st_as_sf(coords = c("long", "lat"), crs = 4326) %>%
    st_transform(crs = crs[which_country]) -> base_loc
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
    mutate(
      X = .data$X / 1e3, Y = .data$Y / 1e3, secondary = NA_real_,
      iyear = .data$year - min(.data$year) + 1,
      log_visits = log(.data$visits), secondary = NA_real_
    ) %>%
    add_knots(knots = knots) -> base_data

  base_data %>%
    select(.data$year, .data$iyear, starts_with("knot")) %>%
    distinct() %>%
    bind_rows(
      base_data %>%
        select(.data$year, starts_with("knot")) %>%
        distinct()
    ) %>%
    mutate(secondary = NA_real_) -> trend_prediction

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
      secondary = NA_real_
    ) %>%
    add_knots(knots = knots) %>%
    mutate(secondary = NA_real_, visits = 1) %>%
    arrange(.data$location, .data$year) -> base_prediction

  expand.grid(
    X = pretty(base_data$X, 20),
    Y = pretty(base_data$Y, 20),
    year = knots
  ) -> field_prediction

  results <- fit_model(
    base_data = base_data, trend_prediction = trend_prediction, knots = knots,
    base_prediction = base_prediction, field_prediction = field_prediction
  )

  return(
    c(
      species = species, min_occurrences = min_occurrences,
      min_species = min_species, results, type = "base", country = which_country
    )
  )
}

#' @importFrom assertthat assert_that has_name
#' @importFrom dplyr bind_cols select_if %>%
#' @importFrom INLA inla.mesh.1d inla.mesh.1d.A
add_knots <- function(data, knots =  c(1990, 2000, 2010, 2020)) {
  assert_that(inherits(data, "data.frame"), is.numeric(knots))
  knots <- sort(unique(knots))
  assert_that(length(knots) > 1, has_name(data, "year"))
  if (!is.integer(knots)) {
    assert_that(
      all(abs(as.integer(knots) - knots) < 1e-6),
      msg = "knots must be integer"
    )
    knots <- as.integer(knots)
  }
  mesh <- inla.mesh.1d(loc = knots)
  inla.mesh.1d.A(mesh, data$year) %>%
    as.matrix() %>%
    as.data.frame() %>%
    `colnames<-`(paste0("knot_", knots)) %>%
    select_if(function(x){max(x) > 0}) %>%
    bind_cols(data)
}
