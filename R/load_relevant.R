#' Load the relevant occurrence data
#' @param min_occurrences The minimum number of occurrences per species.
#' @param min_species The minimum number of species recorded at the combination
#'   of location and year.
#' @param country Data from which country to select.
#' @export
#' @importFrom assertthat assert_that is.count
#' @importFrom dplyr filter group_by inner_join n transmute ungroup %>%
#' @importFrom git2rdata read_vc
#' @importFrom magrittr %<>%
#' @importFrom tidyr pivot_wider
load_relevant <- function(
  min_occurrences = 1000, min_species = 3, country = c("BE", "GB", "NL", "all")
) {
  assert_that(is.count(min_occurrences), is.count(min_species))
  which_country <- match.arg(country)
  base_data <- read_vc("location", system.file(package = "ladybird"))
  if (which_country != "all") {
    base_data %<>%
      filter(.data$country == which_country)
  }
  base_data %>%
    semi_join(
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
    ungroup() %>%
    inner_join(
      read_vc("species", system.file(package = "ladybird")),
      by = "taxon_key"
    ) %>%
    transmute(
      .data$year, .data$location, species = .data$code, observed = 1L
    ) %>%
    pivot_wider(
      names_from = .data$species, values_from = .data$observed, values_fill = 0L
    ) %>%
    inner_join(
      read_vc("visits", system.file(package = "ladybird")),
      by = c("year", "location")
    )
}
