# fit_single_variable_lm.R -- fit one linear model per predictor on a CLR response
#
# The reported regression is marginal in every predictor: each coefficient comes from a separate response ~ predictor model, so no estimate is adjusted for the other predictors.

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

# ---- fit_gene_centric_lm_effects -------------------------------------------
# Fit one linear model of a response on a set of predictors.
#
# This is the per-model fitting engine; fit_gene_centric_single_variable_lm_effects() calls it once per predictor, which is what makes every reported coefficient a single-predictor estimate.
#
# Parameters:
#   df_model       — data.frame. Model table containing the response and predictor columns.
#   response       — character scalar. Response column to regress.
#   predictor_cols — character vector. Predictor columns entered into the model.
#   ...            — additional arguments passed to stats::lm().
#
# Returns:
#   list with:
#     model           — the stats::lm fit
#     coefficients    — data.frame with columns term, estimate, std_error, statistic, p_value, conf_low, conf_high; the intercept included
#     predictions     — numeric vector of fitted values for the complete-case rows
#     observed        — numeric vector of observed response values for those rows
#     r_squared       — numeric scalar
#     adj_r_squared   — numeric scalar
#     sigma           — numeric scalar, residual standard error
#     predictors      — character vector of predictor columns fitted
#     response        — character scalar
#     data            — the complete-case rows actually fitted
#     input_row_count — integer, rows offered before incomplete rows were dropped
fit_gene_centric_lm_effects <- function(
  df_model,
  response,
  predictor_cols,
  ...
) {
  stopifnot(
    is.data.frame(df_model),
    is.character(response),
    length(response) == 1L,
    is.character(predictor_cols)
  )

  if (!response %in% names(df_model)) {
    stop(sprintf("Response column '%s' not found in model table.", response))
  }
  missing_predictors <- setdiff(predictor_cols, names(df_model))
  if (length(missing_predictors) > 0L) {
    stop(sprintf(
      "LM-effects model is missing predictor columns: %s",
      paste(missing_predictors, collapse = ", ")
    ))
  }
  if (nrow(df_model) < 3L) {
    stop(sprintf(
      "LM-effects fitting for response '%s' requires at least 3 rows after filtering; received %d.",
      response,
      nrow(df_model)
    ))
  }

  input_row_count <- nrow(df_model)
  fit_df <- df_model[
    stats::complete.cases(df_model[,
      c(response, predictor_cols),
      drop = FALSE
    ]),
    c(response, predictor_cols),
    drop = FALSE
  ]
  if (nrow(fit_df) < 3L) {
    stop(sprintf(
      "LM-effects fitting for response '%s' requires at least 3 complete rows after dropping missing response/predictor values; received %d.",
      response,
      nrow(fit_df)
    ))
  }

  model_formula <- stats::as.formula(sprintf(
    "%s ~ %s",
    response,
    paste(predictor_cols, collapse = " + ")
  ))
  lm_fit <- stats::lm(model_formula, data = fit_df, ...)
  summary_fit <- summary(lm_fit)
  coef_mat <- summary_fit$coefficients
  conf_int <- stats::confint(lm_fit)
  conf_low <- conf_int[rownames(coef_mat), 1]
  conf_high <- conf_int[rownames(coef_mat), 2]

  coefficients <- tibble::tibble(
    term = rownames(coef_mat),
    estimate = coef_mat[, "Estimate"],
    std_error = coef_mat[, "Std. Error"],
    statistic = coef_mat[, "t value"],
    p_value = coef_mat[, "Pr(>|t|)"],
    conf_low = conf_low,
    conf_high = conf_high
  )

  list(
    model = lm_fit,
    coefficients = coefficients,
    predictions = stats::predict(lm_fit, newdata = fit_df),
    observed = fit_df[[response]],
    r_squared = unname(summary_fit$r.squared),
    adj_r_squared = unname(summary_fit$adj.r.squared),
    sigma = unname(summary_fit$sigma),
    predictors = predictor_cols,
    response = response,
    data = fit_df,
    input_row_count = input_row_count
  )
}

# ---- fit_gene_centric_single_variable_lm_effects ---------------------------
# Fit one linear model per predictor and report each predictor's own marginal coefficient.
#
# Every returned coefficient row comes from a separate response ~ predictor model.
# Predictors that are constant or entirely missing should be filtered before calling this function.
#
# Parameters:
#   df_model       — data.frame. Model table containing the response and predictor columns.
#   response       — character scalar. Response column to regress on each predictor.
#   predictor_cols — non-empty character vector. Predictor columns fitted one at a time.
#   ...            — additional arguments passed to stats::lm().
#
# Returns:
#   list with:
#     coefficients    — data.frame with columns single_variable_feature, term, estimate, std_error, statistic, p_value, conf_low, conf_high; the intercept of each model dropped
#     diagnostics     — data.frame with columns single_variable_feature, metric (r_squared, adj_r_squared, sigma), value
#     model_summaries — data.frame with columns single_variable_feature, r_squared, adj_r_squared, sigma, predictions (observed, predicted) and refit (the stats::lm object)
#     predictors      — character vector of predictor columns fitted
#     response        — character scalar
#     data            — the model table that was offered, before each model's complete-case filtering
fit_gene_centric_single_variable_lm_effects <- function(
  df_model,
  response,
  predictor_cols,
  ...
) {
  # === Input Validation (MANDATORY) ===
  stopifnot(
    "df_model must be a data.frame" = is.data.frame(df_model),
    "response must be a character scalar" = is.character(response) &&
      length(response) == 1L &&
      nzchar(response),
    "predictor_cols must be a non-empty character vector" = is.character(
      predictor_cols
    ) &&
      length(predictor_cols) > 0L
  )

  if (!response %in% names(df_model)) {
    stop(sprintf(
      "Single-variable LM-effects fitting cannot find response column '%s'. Expected it in df_model; check response transformation before fitting.",
      response
    ))
  }
  missing_predictors <- setdiff(predictor_cols, names(df_model))
  if (length(missing_predictors) > 0L) {
    stop(sprintf(
      "Single-variable LM-effects fitting is missing predictor columns: %s. Check predictor matrix construction.",
      paste(missing_predictors, collapse = ", ")
    ))
  }

  # === Function Logic ===
  fits <- lapply(predictor_cols, function(predictor) {
    fit <- fit_gene_centric_lm_effects(
      df_model = df_model,
      response = response,
      predictor_cols = predictor,
      ...
    )
    fit$single_variable_feature <- predictor
    fit
  })
  names(fits) <- predictor_cols

  coefficients <- dplyr::bind_rows(lapply(fits, function(fit) {
    single_variable_feature <- fit$single_variable_feature
    fit$coefficients %>%
      dplyr::filter(.data$term != "(Intercept)") %>%
      dplyr::mutate(
        single_variable_feature = .env$single_variable_feature,
        .before = 1L
      )
  }))

  diagnostics <- dplyr::bind_rows(lapply(fits, function(fit) {
    tibble::tibble(
      single_variable_feature = fit$single_variable_feature,
      metric = c("r_squared", "adj_r_squared", "sigma"),
      value = c(fit$r_squared, fit$adj_r_squared, fit$sigma)
    )
  }))

  model_summaries <- dplyr::bind_rows(lapply(fits, function(fit) {
    tibble::tibble(
      single_variable_feature = fit$single_variable_feature,
      r_squared = fit$r_squared,
      adj_r_squared = fit$adj_r_squared,
      sigma = fit$sigma,
      predictions = list(tibble::tibble(
        observed = fit$observed,
        predicted = as.numeric(fit$predictions)
      )),
      refit = list(fit$model)
    )
  }))

  # === Return Value ===
  return(list(
    coefficients = coefficients,
    diagnostics = diagnostics,
    model_summaries = model_summaries,
    predictors = predictor_cols,
    response = response,
    data = df_model
  ))
}
