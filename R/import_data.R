#' Import and standardise the raw data
#' @param belgium_occurrence path to the CSV file with the Belgian occurrence
#' data.
#' @param belgium_visits path to the CSV file with the Belgian visit data.
#' @param output path to the root of the data package
#' @inheritParams git2rdata::write_vc
#' @export
#' @importFrom assertthat assert_that has_name is.string
#' @importFrom git2rdata repository write_vc
#' @importFrom dplyr anti_join bind_cols distinct filter mutate semi_join
#' transmute %>%
#' @importFrom readr locale read_delim
#' @importFrom rlang .data
#' @importFrom sf st_as_sf st_coordinates st_transform
import_data <- function(
  belgium_occurrence, belgium_visits, output, strict = TRUE
) {
  assert_that(is.string(belgium_occurrence), is.string(belgium_visits))
  root <- repository(output)
  raw_belgium <- read_delim(belgium_occurrence, delim = ";", locale = locale())
  assert_that(
    all(
      has_name(
        raw_belgium,
        c("taxon_key", "scientific_name", "year", "utm1km", "centroid_x",
          "centroid_y")
      )
    )
  )
  raw_bel_visit <- read_delim(belgium_visits, delim = ";", locale = locale())
  assert_that(all(has_name(raw_bel_visit, c("year", "utm1km", "nDayVisits"))))
  raw_belgium %>%
    filter(.data$year >= 1991) %>%
    distinct(.data$year, .data$utm1km) %>%
    anti_join(raw_bel_visit, by = c("year", "utm1km")) %>%
    nrow() -> missing_visits
  assert_that(
    missing_visits == 0,
    msg = "Some Belgian occurrence data have no matching number of visits."
  )

  st_as_sf(
    raw_belgium, coords = c("centroid_x", "centroid_y"), crs = 31370
  ) %>%
    st_transform(crs = 4326) %>%
    st_coordinates() %>%
    as.data.frame() %>%
    bind_cols(raw_belgium) %>%
    transmute(
      .data$taxon_key, .data$scientific_name, .data$year,
      country = factor("BE", levels = c("BE", "NL", "GB")),
      location = .data$utm1km, long = .data$X, lat = .data$Y
    ) -> raw_belgium
  raw_belgium %>%
    filter(.data$year >= 1991) %>%
    mutate(
      location = factor(.data$location),
      taxon_key = factor(.data$taxon_key)
    ) -> raw_data
  raw_bel_visit %>%
    transmute(
      .data$year,
      location = factor(.data$utm1km, levels = levels(raw_data$location)),
      visits = .data$nDayVisits
    ) %>%
    filter(!is.na(.data$location)) %>%
    write_vc(
      "inst/visits", root = root, sorting = c("year", "location"), stage = TRUE,
      strict = strict
    )
  raw_data %>%
    distinct(.data$taxon_key, .data$scientific_name) %>%
    mutate(
      code = gsub("^(.{4}).* (.{4}).*$", "\\1_\\2", .data$scientific_name)
    ) %>%
    write_vc(
      file = "inst/species", root = root, sorting = "taxon_key", stage = TRUE,
      strict = strict
    )
  raw_data %>%
    distinct(.data$country, .data$location, .data$long, .data$lat) %>%
    write_vc(
      file = "inst/location", root = root, sorting = "location", stage = TRUE,
      strict = strict
    )
  raw_data %>%
    distinct(location = factor(.data$location), .data$year, .data$taxon_key) %>%
    write_vc(
      file = "inst/occurrence", root = root, strict = strict,
      sorting = c("location", "year", "taxon_key"), stage = TRUE
    )
}
