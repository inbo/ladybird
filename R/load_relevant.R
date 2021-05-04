#' Load the relevant occurrence data
#' @param min_occurrences The minimum number of occurrences per species.
#' @param min_species The minimum number of species recorded at the combination
#'   of location and year.
#' @param country Data from which country to select.
#' @param buffer_distance The range in which to look for other locations.
#' @param buffer_locations  The minimal number of locations with suffient data
#' within the buffer range
#' @export
#' @importFrom assertthat assert_that is.count
#' @importFrom dplyr distinct filter group_by inner_join n transmute ungroup %>%
#' @importFrom git2rdata read_vc
#' @importFrom magrittr %<>%
#' @importFrom sf st_as_sf st_buffer st_distance st_intersection st_transform
#' st_union
#' @importFrom tidyr pivot_wider
#' @importFrom units as_units
load_relevant <- function(
  min_occurrences = 1000, min_species = 3, country = c("BE", "GB", "NL", "all"),
  buffer_distance = 50e3, buffer_locations = 50
) {
  assert_that(is.count(min_occurrences), is.count(min_species))
  which_country <- match.arg(country)
  base_data <- read_vc("location", system.file(package = "ladybird"))
  if (which_country != "all") {
    base_data %<>%
      filter(.data$country == which_country)
  }
  base_data %>%
    inner_join(
      x = read_vc("occurrence", system.file(package = "ladybird")),
      by = "location"
    ) %>%
    group_by(.data$taxon_key) %>%
    filter(n() >= min_occurrences) -> relevant
  relevant %>%
    filter(.data$taxon_key != "4989904") %>%
    group_by(.data$location, .data$year) %>%
    filter(n() >= min_species) %>%
    semi_join(x = relevant, by = c("location", "year")) %>%
    group_by(.data$taxon_key) %>%
    filter(n() >= min_occurrences) -> relevant
  relevant %>%
    filter(.data$taxon_key != "4989904") %>%
    group_by(.data$location, .data$year) %>%
    filter(n() >= min_species) %>%
    semi_join(x = relevant, by = c("location", "year")) %>%
    ungroup() -> relevant
  crs <- c(BE = 31370, GB = 27700, NL = 28992, all = 3035)
  relevant %>%
    distinct(.data$location, .data$long, .data$lat) %>%
    st_as_sf(coords = c("long", "lat"), crs = 4326) %>%
    st_transform(crs = crs[which_country]) -> locations
  obs_dist <- st_distance(locations)
  locations %<>%
    mutate(
      buffer_count = rowSums(obs_dist <= as_units(buffer_distance, "m"))
    )
  suppressWarnings(
    {
      locations %>%
        filter(.data$buffer_count >= buffer_locations) %>%
        st_buffer(buffer_distance) %>%
        st_union() %>%
        st_intersection(x = locations) -> locations
    }
  )

  locations %>%
    st_drop_geometry() %>%
    inner_join(x = relevant, by = "location") %>%
    inner_join(
      read_vc("species", system.file(package = "ladybird")),
      by = "taxon_key"
    ) %>%
    transmute(
      .data$year, .data$location, .data$buffer_count, species = .data$code,
      observed = 1L
    ) %>%
    pivot_wider(
      names_from = .data$species, values_from = .data$observed, values_fill = 0L
    ) %>%
    inner_join(
      read_vc("visits", system.file(package = "ladybird")),
      by = c("year", "location")
    )
}
