#' Fit a model to a species using the predictions for a secondary species
#' @inheritParams base_model
#' @param base_data A dataframe with the base data.
#' @param trend_prediction A dataframe with the timestamps to predict the trend.
#' @param base_prediction A dataframe with the locations and timestamps to
#' predict.
#' @export
#' @importFrom assertthat assert_that has_name is.flag noNA
#' @importFrom dplyr as_tibble bind_cols distinct select %>%
#' @importFrom INLA inla inla.mesh.1d inla.mesh.2d inla.spde2.pcmatern
#' inla.spde.make.A inla.spde.make.index inla.stack inla.stack.A inla.stack.data
#' inla.stack.index
#' @importFrom pROC roc
#' @importFrom rlang .data
#' @importFrom stats as.formula
fit_model <- function(
  first_order = TRUE, base_data, trend_prediction, base_prediction
) {
  assert_that(is.flag(first_order), noNA(first_order))
  assert_that(inherits(base_data, "data.frame"))
  assert_that(inherits(trend_prediction, "data.frame"))
  assert_that(inherits(base_prediction, "data.frame"))
  assert_that(
    has_name(base_data, "iyear"), has_name(base_data, "iyear2"),
    has_name(base_data, "X"), has_name(base_data, "Y"),
    has_name(base_data, "secondary")
  )
  assert_that(
    has_name(trend_prediction, "iyear"), has_name(trend_prediction, "iyear2"),
    has_name(base_prediction, "secondary")
  )
  assert_that(
    has_name(base_prediction, "iyear"), has_name(base_prediction, "iyear2"),
    has_name(base_prediction, "X"), has_name(base_prediction, "Y"),
    has_name(base_prediction, "secondary")
  )
  is_secondary <- any(!is.na(base_data$secondary))
  if (is_secondary) {
    time_vars <- c("before", "after")
    assert_that(
      has_name(base_data, "before"), has_name(base_data, "after"),
      has_name(trend_prediction, "before"), has_name(trend_prediction, "after"),
      has_name(base_prediction, "before"), has_name(base_prediction, "after")
    )
  } else {
    time_vars <- "cyear"
    assert_that(
      has_name(base_data, "cyear"), has_name(trend_prediction, "cyear"),
      has_name(base_prediction, "cyear")
    )
  }
  n_year <- max(base_data$iyear)
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
      mesh = mesh, group = base_data$iyear, mesh.group = time_mesh
    ) -> a_estimate
  site_index <- inla.spde.make.index(
    name = "site", n.spde = mesh$n, n.group = time_mesh$n
  )
  inla.stack(
    data = select(base_data, .data$occurrence),
    A = list(a_estimate, 1),
    effects = list(
      c(site_index, list(intercept = 1)),
      list(
        base_data %>%
          select(!!time_vars, .data$iyear, .data$iyear2, .data$secondary)
      )
    ),
    tag = "estimate"
  ) -> stack_estimate
  inla.stack(
    data = data.frame(occurrence = NA),
    A = list(1),
    effects = list(list(trend_prediction)),
    tag = "trend"
  ) -> stack_trend
  base_prediction %>%
    select(.data$X, .data$Y) %>%
    as.matrix() %>%
    inla.spde.make.A(
      mesh = mesh, group = base_prediction$iyear, mesh.group = time_mesh
    ) -> a_prediction
  inla.stack(
    data = data.frame(occurrence = NA),
    A = list(a_prediction, 1),
    effects = list(
      c(site_index, list(intercept = 1)),
      list(
        base_prediction %>%
          select(!!time_vars, .data$iyear, .data$iyear2, .data$secondary )
      )
    ),
    tag = "prediction"
  ) -> stack_prediction
  inla.stack(stack_estimate, stack_trend) -> stack
  inla.stack(stack_estimate, stack_prediction) -> stack2
  fixed_formula <- ifelse(
    is_secondary,
    "occurrence ~ 0 + intercept + before + after",
    "occurrence ~ 0 + intercept + cyear"
  )
  rw_formula <- ifelse(
    first_order,
    ifelse(
      is_secondary,
      "f(
        iyear, model = \"rw1\",
        hyper = list(theta = list(prior = \"pc.prec\", param = c(0.1, 0.05)))
      ) +
      f(
        iyear2, secondary, copy = \"iyear\",
      hyper = list(beta = list(fixed = FALSE))
      )",
      "f(
        iyear, model = \"rw1\",
        hyper = list(theta = list(prior = \"pc.prec\", param = c(0.1, 0.05)))
      )"
    ),
    ifelse(
      is_secondary,
      "f(
        iyear, model = \"rw2\",
        hyper = list(theta = list(prior = \"pc.prec\", param = c(0.02, 0.05)))
      ) +
      f(
        iyear2, secondary, copy = \"iyear\",
        hyper = list(beta = list(fixed = FALSE))
      )",
      "f(
        iyear, model = \"rw2\",
        hyper = list(theta = list(prior = \"pc.prec\", param = c(0.02, 0.05)))
      )"
    )
  )
  st_formula <- "f(
      site, model = spde, group = site.group,
      control.group = list(
        model = \"ar1\",
        hyper = list(theta = list(prior = \"pc.cor1\", param = c(0.6, 0.7)))
      )
    )"
  paste(fixed_formula, rw_formula, st_formula, sep = " +\n") %>%
    as.formula() -> model_formula
  m0 <- inla(
    model_formula, family = "binomial", data = inla.stack.data(stack_estimate),
    control.predictor = list(A = inla.stack.A(stack_estimate), compute = FALSE)
  )
  m1 <- inla(
    model_formula, family = "binomial", data = inla.stack.data(stack),
    control.compute = list(waic = TRUE),
    control.predictor = list(A = inla.stack.A(stack), compute = TRUE, link = 1),
    control.mode = list(theta = m0$mode$theta, restart = FALSE, fixed = TRUE)
  )
  m2 <- inla(
    model_formula, family = "binomial", data = inla.stack.data(stack2),
    control.predictor = list(
      A = inla.stack.A(stack2), compute = TRUE, link = 1
    ),
    control.mode = list(theta = m0$mode$theta, restart = FALSE, fixed = TRUE)
  )
  index_auc <- inla.stack.index(stack2, "estimate")$data
  roc_curve <- roc(
    predictor = m2$summary.fitted.values[index_auc, "mean"],
    response = base_data$occurrence,
    levels = c(0, 1), direction = "<", ci = TRUE
  )
  index_estimate <- inla.stack.index(stack2, "prediction")$data
  m2$summary.fitted.values[index_estimate, ] %>%
    select(.data$mean, median = 4, lcl = 3, ucl = 5) %>%
    as_tibble() %>%
    bind_cols(
      m2$summary.linear.predictor[index_estimate, ] %>%
        select(lp_mean = .data$mean, lp_median = 4, lp_lcl = 3, lp_ucl = 5),
      base_prediction %>%
        select(.data$year, .data$location)
    ) -> predictions
  index_trend <- inla.stack.index(stack, "trend")$data
  m1$summary.fitted.values[index_trend, ] %>%
    select(.data$mean, median = 4, lcl = 3, ucl = 5) %>%
    as_tibble() %>%
    bind_cols(trend_prediction) -> trend
  return(
    list(
      fixed = m0$summary.fixed, trend = trend,
      hyperpar = m0$summary.hyperpar, predictions = predictions,
      waic = m1$waic$waic, first_order = first_order, roc = roc_curve
    )
  )
}
