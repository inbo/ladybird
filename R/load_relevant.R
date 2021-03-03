#' Load the relevant occurrence data
#' @param min_occurrences The minimum number of occurrences per species.
#' @param min_species The minimum number of species recorded at the combination
#'   of location and year.
#' @export
#' @importFrom assertthat assert_that is.count
#' @importFrom dplyr filter group_by n ungroup
#' @importFrom git2rdata read_vc
load_relevant <- function(min_occurrences = 1000, min_species = 3) {
  assert_that(is.count(min_occurrences), is.count(min_species))
  read_vc("occurrence", system.file(package = "ladybird")) %>%
    group_by(.data$taxon_key) %>%
    filter(n() >= min_occurrences) %>%
    group_by(.data$location, .data$year) %>%
    filter(n() >= min_species) %>%
    group_by(.data$taxon_key) %>%
    filter(n() >= min_occurrences) %>%
    group_by(.data$location, .data$year) %>%
    filter(n() >= min_species) %>%
    ungroup()
}
