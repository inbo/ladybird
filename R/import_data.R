#' Import and standardise the raw data
#' @param belgium path to the CSV file with the Belgian data.
#' @param output path to the root of the data package
#' @inheritParams git2rdata::write_vc
#' @export
#' @importFrom assertthat assert_that has_name is.string
#' @importFrom git2rdata repository write_vc
#' @importFrom dplyr bind_cols distinct filter mutate transmute %>%
#' @importFrom readr locale read_delim
#' @importFrom rlang .data
#' @importFrom sf st_as_sf st_coordinates st_transform
import_data <- function(belgium, output, strict = TRUE) {
  assert_that(is.string(belgium))
  root <- repository(output)
  raw_belgium <- read_delim(belgium, delim = ";", locale = locale())
  assert_that(
    all(
      has_name(
        raw_belgium,
        c("taxon_key", "scientific_name", "year", "utm1km", "centroid_x",
          "centroid_y")
      )
    )
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
    filter(.data$year > 1990) %>%
    mutate(
      location = factor(.data$location),
      taxon_key = factor(.data$taxon_key)
    ) -> raw_data
  raw_data %>%
    distinct(.data$taxon_key, .data$scientific_name) %>%
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
