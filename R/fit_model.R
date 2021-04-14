#' Fit a model to a species using the predictions for a secondary species
#' @inheritParams base_model
#' @param base_data A dataframe with the base data.
#' @param trend_prediction A dataframe with the timestamps to predict the trend.
#' @param base_prediction A dataframe with the locations and timestamps to
#' predict.
#' @param field_prediction A dataframe with locations and timestamps to predict
#' the spatiotemporal field.
#' @export
#' @importFrom assertthat assert_that has_name is.flag noNA
#' @importFrom dplyr as_tibble bind_cols distinct select %>%
#' @importFrom INLA inla inla.make.lincombs inla.mesh.1d inla.mesh.2d
#' inla.spde2.pcmatern inla.spde.make.A inla.spde.make.index inla.stack
#' inla.stack.A inla.stack.data inla.stack.index
#' @importFrom pROC roc
#' @importFrom rlang .data
#' @importFrom stats as.formula
#' @importFrom tibble rownames_to_column
#' @importFrom tidyselect ends_with
fit_model <- function(
  base_data, trend_prediction, base_prediction, field_prediction, knots
) {
  assert_that(inherits(base_data, "data.frame"))
  assert_that(inherits(trend_prediction, "data.frame"))
  assert_that(inherits(base_prediction, "data.frame"))
  assert_that(inherits(field_prediction, "data.frame"))
  assert_that(
    has_name(base_data, "iyear"),
    has_name(base_data, "X"), has_name(base_data, "Y"),
    has_name(base_data, "secondary"), has_name(base_data, "log_visits")
  )
  assert_that(
    has_name(trend_prediction, "iyear"), has_name(trend_prediction, "secondary")
  )
  assert_that(
    has_name(base_prediction, "iyear"), has_name(base_prediction, "secondary"),
    has_name(base_prediction, "X"), has_name(base_prediction, "Y")
  )
  is_secondary <- any(!is.na(base_data$secondary))

  time_mesh <- inla.mesh.1d(loc = knots)
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
      mesh = mesh, group = base_data$year, mesh.group = time_mesh
    ) -> a_estimate
  site_index <- inla.spde.make.index(
    name = "site", n.spde = mesh$n, n.group = time_mesh$n
  )
  time_vars <- colnames(base_data)
  time_vars <- time_vars[grepl("knot", time_vars)]
  inla.stack(
    data = select(base_data, .data$occurrence),
    A = list(a_estimate, 1),
    effects = list(
      c(site_index, intercept = 1),
      list(
        base_data %>%
          select(
            !!time_vars, .data$year, .data$iyear, .data$secondary,
            .data$log_visits
          )
      )
    ),
    tag = "estimate"
  ) -> stack_estimate

  fixed_formula <- sprintf(
    "occurrence ~ 0 + log_visits + %s", paste(time_vars, collapse = " + ")
  )
  rw_formula <- "f(
    iyear, model = \"rw1\",
    hyper = list(theta = list(prior = \"pc.prec\", param = c(0.1, 0.05)))
  )"
  st_formula <- "f(
    site, model = spde, group = site.group,
    control.group = list(
      model = \"ar1\",
      hyper = list(theta = list(prior = \"pc.cor1\", param = c(0.6, 0.7)))
    )
  )"
  paste(fixed_formula, rw_formula, st_formula, sep = " +\n") %>%
    as.formula() -> model_formula

  data.frame(
    after = tail(time_vars, -1),
    before = head(time_vars, -1)
  ) %>%
    mutate(
      decade = gsub("^knot_([0-9]{4}).*", "\\1", .data$before),
      period = gsub("^knot_([0-9]{4}).*", "\\1", .data$after) %>%
        as.integer() %>%
        `-`(as.integer(.data$decade))
    ) %>%
    filter(.data$period > 0) %>%
    pivot_longer(c("after", "before")) %>%
    transmute(
      .data$decade, .data$value,
      weight = ifelse(.data$name == "after", 1, -1) / .data$period
    ) %>%
    pivot_wider(names_from = "value", values_from = "weight") %>%
    as.data.frame() -> lc_base
  if (is_secondary) {
    lc_base %>%
      mutate(
        across(
          ends_with("secondary"), function(x){
            NA
          }
        )
      ) %>%
      bind_rows(
        lc_base %>%
          mutate(
            across(
              paste0("knot_", knots), function(x){
                NA
              }
            ),
            decade = paste0(.data$decade, ":0")
          ),
        lc_base %>%
          mutate(decade = paste0(.data$decade, ":1"))
      ) -> lc_base
    lc_base %>%
      select(-decade) %>%
      is.na() %>%
      apply(1, all) -> all_na
    lc_base <- lc_base[!all_na, ]
    lc_base %>%
      select(-decade) %>%
      anyDuplicated() -> duplicates
    lc_base <- lc_base[-duplicates, ]
  }
  lc_base %>%
    select(-.data$decade) %>%
    inla.make.lincombs() %>%
    setNames(lc_base$decade) -> lc

  m0 <- inla(
    model_formula, family = "binomial",
    Ntrials = 1,
    data = inla.stack.data(stack_estimate),
    control.predictor = list(A = inla.stack.A(stack_estimate), compute = FALSE),
    lincomb = lc
  )
  m0$summary.lincomb.derived %>%
    rownames_to_column(var = "decade") %>%
    select(
      "decade", "mean", median = "0.5quant", lcl = "0.025quant",
      ucl = "0.975quant"
    ) -> lincomb
  fixed <- m0$summary.fixed
  hyperpar <- m0$summary.hyperpar

  inla.stack(
    data = data.frame(occurrence = NA),
    A = list(1),
    effects = list(list(trend_prediction)),
    tag = "trend"
  ) -> stack_trend
  inla.stack(stack_estimate, stack_trend) -> stack
  m1 <- inla(
    model_formula, family = "binomial",
    Ntrials = 1,
    data = inla.stack.data(stack), control.compute = list(waic = TRUE),
    control.predictor = list(A = inla.stack.A(stack), compute = TRUE, link = 1),
    control.mode = list(theta = m0$mode$theta, restart = FALSE, fixed = TRUE)
  )
  index_auc <- inla.stack.index(stack, "estimate")$data
  roc_curve <- roc(
    predictor = m1$summary.fitted.values[index_auc, "mean"],
    response = base_data$occurrence,
    levels = c(0, 1), direction = "<", ci = TRUE
  )
  index_trend <- inla.stack.index(stack, "trend")$data
  m1$summary.fitted.values[index_trend, ] %>%
    select(.data$mean, median = 4, lcl = 3, ucl = 5) %>%
    as_tibble() %>%
    bind_cols(trend_prediction) -> trend
  waic <- m1$waic$waic
  rm(m1)
  gc()

  base_prediction %>%
    select(.data$X, .data$Y) %>%
    as.matrix() %>%
    inla.spde.make.A(
      mesh = mesh, group = base_prediction$year, mesh.group = time_mesh
    ) -> a_prediction
  inla.stack(
    data = data.frame(occurrence = NA),
    A = list(a_prediction, 1),
    effects = list(
      c(site_index, list(intercept = 1)),
      list(
        base_prediction %>%
          select(!!time_vars, .data$iyear, .data$secondary)
      )
    ),
    tag = "prediction"
  ) -> stack_prediction
  field_prediction %>%
    select(.data$X, .data$Y) %>%
    as.matrix() %>%
    inla.spde.make.A(
      mesh = mesh, group = field_prediction$year, mesh.group = time_mesh
    ) -> a_field
  inla.stack(
    data = data.frame(occurrence = NA),
    A = list(a_field),
    effects = list(
      c(site_index, list(intercept = 1))
    ),
    tag = "field"
  ) -> stack_field
  inla.stack(stack_estimate, stack_prediction, stack_field) -> stack2
  m2 <- inla(
    model_formula, family = "binomial",
    Ntrials = 1,
    data = inla.stack.data(stack2),
    control.predictor = list(
      A = inla.stack.A(stack2), compute = TRUE, link = 1
    ),
    control.mode = list(theta = m0$mode$theta, restart = FALSE, fixed = TRUE)
  )
  index_estimate <- inla.stack.index(stack2, "prediction")$data
  m2$summary.fitted.values[index_estimate, ] %>%
    select(.data$mean, median = 4, lcl = 3, ucl = 5) %>%
    as_tibble() %>%
    bind_cols(
      m2$summary.linear.predictor[index_estimate, ] %>%
        select(lp_mean = .data$mean, lp_median = 4, lp_lcl = 3, lp_ucl = 5),
      base_prediction %>%
        select(.data$year, .data$location, .data$X, .data$Y)
    ) -> predictions
  index_field <- inla.stack.index(stack2, "field")$data
  m2$summary.linear.predictor[index_field, ] %>%
    select(.data$mean, median = 4, lcl = 3, ucl = 5) %>%
    as_tibble() %>%
    bind_cols(field_prediction) -> field
  rm(m0, m2)
  gc()


  return(
    list(
      fixed = fixed, trend = trend, field = field, hyperpar = hyperpar,
      predictions = predictions, waic = waic, roc = roc_curve, lincomb = lincomb
    )
  )
}
