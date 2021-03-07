#' Fit a base model to a species
#' @param species Name of the species
#' @inheritParams load_relevant
#' @export
#' @importFrom dplyr arrange as_tibble bind_cols distinct inner_join mutate
#' select %>%
#' @importFrom git2rdata read_vc
#' @importFrom INLA inla inla.mesh.1d inla.mesh.2d inla.spde2.pcmatern
#' inla.spde.make.A inla.spde.make.index inla.stack inla.stack.A inla.stack.data
#' inla.stack.index
#' @importFrom rlang .data !!
#' @importFrom sf st_as_sf st_coordinates st_drop_geometry st_transform
base_model <- function(
  species = "Hrma", min_occurrences = 1000, min_species = 3
) {
  read_vc("location", system.file(package = "ladybird")) %>%
    st_as_sf(coords = c("long", "lat"), crs = 4326) %>%
    st_transform(crs = 31370) %>%
    bind_cols(
      st_coordinates(.) %>%
        `/`(1e3) %>%
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
    mutate(cyear = .data$year - min(.data$year) + 1) -> base_data
  n_year <- max(base_data$cyear)
  year_interval <- 5
  seq(year_interval / 2, n_year, by = year_interval) %>%
    inla.mesh.1d() -> time_mesh
  inla.mesh.2d(
    loc.domain = distinct(base_data, .data$X, .data$Y),
    max.edge = 10
  ) -> mesh
  spde <- inla.spde2.pcmatern(
    mesh = mesh, prior.range = c(50, 0.5), prior.sigma = c(0.1, 0.05)
  )
  base_data %>%
    select(.data$X, .data$Y) %>%
    as.matrix() %>%
    inla.spde.make.A(
      mesh = mesh, group = base_data$cyear, mesh.group = time_mesh
    ) -> a_estimate
  site_index <- inla.spde.make.index(
    name = "site", n.spde = mesh$n, n.group = time_mesh$n
  )
  inla.stack(
    data = select(base_data, .data$occurrence),
    A = list(a_estimate, 1),
    effects = list(
      c(site_index, list(intercept = 1)),
      list(select(base_data, .data$cyear))
    ),
    tag = "estimate"
  ) -> stack_estimate
  inla.stack(
    data = data.frame(occurrence = NA),
    A = list(1),
    effects = list(
      list(
        base_data %>%
          distinct(.data$cyear) %>%
          arrange(.data$cyear)
      )
    ),
    tag = "trend"
  ) -> stack_trend
  inla.stack(stack_estimate, stack_trend) -> stack
  model <- inla(
    occurrence ~ 0 + intercept +
      f(
        cyear, model = "rw1",
        hyper = list(theta = list(prior = "pc.prec", param = c(0.1, 0.05)))
      ) +
      f(
        site, model = spde, group = site.group,
        control.group = list(
          model = "ar1",
          hyper = list(theta = list(prior = "pc.cor1", param = c(0.6, 0.7)))
        )
      ),
    family = "binomial",
    data = inla.stack.data(stack),
    control.compute = list(waic = TRUE),
    control.predictor = list(A = inla.stack.A(stack), compute = TRUE, link = 1)
  )
  index_estimate <- inla.stack.index(stack, "estimate")$data
  model$summary.fitted.values[index_estimate, ] %>%
    select(.data$mean, median = 4, lcl = 3, ucl = 5) %>%
    bind_cols(base_data) %>%
    as_tibble() %>%
    select(
      .data$year, .data$location, .data$mean, .data$median, .data$lcl, .data$ucl
    ) -> predictions
  index_trend <- inla.stack.index(stack, "trend")$data
  model$summary.fitted.values[index_trend, ] %>%
    select(.data$mean, median = 4, lcl = 3, ucl = 5) %>%
    as_tibble() %>%
    bind_cols(
      base_data %>%
        distinct(.data$year) %>%
        arrange(.data$year)
    ) -> trend
  return(
    list(
      species = species, min_occurrences = min_occurrences,
      min_species = min_species, fixed = model$summary.fixed, trend = trend,
      hyperpar = model$summary.hyperpar, predictions = predictions,
      waic = model$waic$waic
    )
  )
}
