#' Prepare the results from the smoother models for the report
#'
#' @param models the path to the individual model results
#' @param maps Location of the Natural Earth Data.
#'  See `get_country_grid()`.
#' @param output the path for the output file
#' @export
#' @importFrom assertthat assert_that is.string
#' @importFrom dplyr arrange bind_cols bind_rows mutate %>%
#' @importFrom effectclass classification coarse_classification
#' @importFrom sf st_as_sf st_centroid st_coordinates st_transform
#' @importFrom stringr str_replace_all
#' @importFrom purrr map map_chr map_dbl map2
prepare_smoother_results <- function(models = ".", maps = ".", output = ".") {
  assert_that(is.string(models))
  models <- normalizePath(models)
  assert_that(is.string(maps))
  maps <- normalizePath(maps)
  assert_that(is.string(output))
  output <- normalizePath(output)
  list.files(models, full.names = TRUE, recursive = TRUE) %>%
    map(readRDS) -> base_models
  species <- map_chr(base_models, "species")
  country <- map_chr(base_models, "country")
  map(base_models, "hyperpar") %>%
    map(rownames_to_column, var = "parameter") %>%
    map2(species, ~mutate(.x, species = .y)) %>%
    map2(country, ~mutate(.x, country = .y)) %>%
    bind_rows() -> hyperpar
  map(base_models, "fixed") %>%
    map(rownames_to_column, var = "parameter") %>%
    map2(species, ~mutate(.x, species = .y)) %>%
    map2(country, ~mutate(.x, country = .y)) %>%
    bind_rows() -> fixed
  map(base_models, "trend") %>%
    map2(species, ~mutate(.x, species = .y)) %>%
    map2(country, ~mutate(.x, country = .y)) %>%
    bind_rows() %>%
    arrange(species, country) %>%
    mutate(
      id = str_replace_all(.data$species, "_", "-")
    ) -> trend
  map(base_models, "roc") %>%
    map("ci") %>%
    map(
      ~data.frame(parameter = c("lcl", "estimate", "ucl"), auc = as.vector(.x))
    ) %>%
    map2(species, ~mutate(.x, species = .y)) %>%
    map2(country, ~mutate(.x, country = .y)) %>%
    bind_rows() -> auc
  thresholds <- log(0.9) * c(-1, 1)
  map(base_models, "lincomb") %>%
    map2(species, ~mutate(.x, species = .y)) %>%
    map2(country, ~mutate(.x, country = .y)) %>%
    bind_rows() %>%
    mutate(
      change = classification(.data$lcl, .data$ucl, thresholds) %>%
        as.character(),
      duration = factor(.data$duration)
    ) -> lincomb
  map(base_models, "prediction") %>%
    map2(species, ~mutate(.x, species = .y)) %>%
    map2(country, ~mutate(.x, country = .y)) %>%
    bind_rows() %>%
    mutate(
      X = round(.data$X),
      Y = round(.data$Y)
    ) -> predictions
  get_country_grid(
    path = maps, country = "BE", cellsize = 10e3, what = "polygons"
  ) %>%
    st_as_sf() -> be_poly
  get_country_grid(
    path = maps, country = "NL", cellsize = 10e3, what = "polygons"
  ) %>%
    st_as_sf() -> nl_poly
  get_country_grid(
    path = maps, country = "GB", cellsize = 15e3, what = "polygons"
  ) %>%
    st_as_sf() -> gb_poly
  bind_cols(
    be_poly,
    st_centroid(be_poly) %>%
      st_coordinates() %>%
      as.data.frame() %>%
      mutate(country = "BE")
    ) %>%
    st_transform(crs = 4326) %>%
    bind_rows(
      bind_cols(
        nl_poly,
        st_centroid(nl_poly) %>%
          st_coordinates() %>%
          as.data.frame() %>%
          mutate(country = "NL")
        ) %>%
        st_transform(crs = 4326),
      bind_cols(
        gb_poly,
        st_centroid(gb_poly) %>%
          st_coordinates() %>%
          as.data.frame() %>%
          mutate(country = "GB")
        ) %>%
        st_transform(crs = 4326)
    ) %>%
    mutate(
      X = round(.data$X),
      Y = round(.data$Y)
    ) -> hex_grid
  hex_grid %>%
    inner_join(predictions, by = c("X", "Y", "country")) %>%
    mutate(
      presence = classification(.data$lcl, .data$ucl, c(0.25, 0.75), 0.5),
      presence2 = coarse_classification(.data$presence) %>%
        factor(
          levels = c("+", "~", "-", "?"),
          labels = c("50% - 100%", "25% - 75%", "0% - 50%", "0% - 100%")
        ),
      presence = factor(
        .data$presence,
        levels = c("++", "+", "+~", "~", "-~", "-", "--", "?+", "?-", "?"),
        labels = c(
          "75% - 100%", "50% - 100%", "50% - 75%", "25% - 75%",  "25% - 50%",
          "0% - 50%", "0% - 25%", "25% - 100%", "0% - 75%", "0% - 100%"
        )
      )
    ) -> predictions
  be_segment <- create_mesh_segments("BE")
  be_border <- get_country_grid(maps, "BE", what = "borders")
  nl_segment <- create_mesh_segments("NL")
  nl_border <- get_country_grid(maps, "NL", what = "borders")
  gb_segment <- create_mesh_segments("GB")
  gb_border <- get_country_grid(maps, "GB", what = "borders")

  save(
    auc, be_border, be_segment, fixed, gb_border, gb_segment, hyperpar, lincomb,
    nl_border, nl_segment, predictions, thresholds, trend,
    file = file.path(output, "smoothers.Rdata")
  )
  return(invisible(NULL))
}

#' @importFrom dplyr filter semi_join %>%
#' @importFrom INLA inla.mesh.2d
#' @importFrom sf as_Spatial st_as_sf st_as_sfc st_buffer st_linestring
#' st_transform st_union
create_mesh_segments <- function(country) {
  crs <- c(BE = 31370, GB = 27700, NL = 28992, all = 3035)
  load_relevant(
    min_occurrences = 1000, min_species = 3, buffer_distance = 50e3,
    buffer_locations = 50, country = country
  ) %>%
    filter(.data$buffer_count >= 50) %>%
    semi_join(
      x = read_vc("location", system.file(package = "ladybird")),
      by = "location"
    ) %>%
    st_as_sf(coords = c("long", "lat"), crs = 4326) %>%
    st_transform(crs = crs[country]) %>%
    st_buffer(50e3) %>%
    st_union() %>%
    st_buffer(50e3) %>%
    st_buffer(-50e3) -> closing
  inla.mesh.2d(
    boundary = as_Spatial(closing), max.edge = 10e3, cutoff = 5e3
  ) -> mesh
  apply(
    mesh$graph$tv[, 1:2], 1, function(x) {
      list(st_linestring(mesh$loc[x, ]))
    }
  ) %>%
    c(
      apply(
        mesh$graph$tv[, 2:3], 1, function(x) {
          list(st_linestring(mesh$loc[x, ]))
        }
      ),
      apply(
        mesh$graph$tv[, c(1, 3)], 1, function(x) {
          list(st_linestring(mesh$loc[x, ]))
        }
      )
    ) %>%
    unlist(recursive = FALSE) %>%
    st_as_sfc(crs = crs[country])
}
