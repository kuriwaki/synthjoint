#' Decompose formula into parts useful for subsequent analyses
#'
#' Outputs a list. You can use `base::list2env` on the output to release
#'  the items of the list on your environment.
#'
#' @param formula A representation of the aggregate imputation or "outcome" model,
#'  of the form `X_{K} ~ X_1 + ... X_{K - 1}`
#'
#' @importFrom Formula as.Formula
#' @keywords internal
#' @examples
#'  synthjoint:::formula_parts(race ~ female + age + edu)
formula_parts <- function(formula) {
  Form         <- as.Formula(formula)
  outcome_var  <- all.vars(formula(Form, lhs = 1, rhs = 0))
  outcome_form <- stats::as.formula(paste0("~ ", outcome_var, "- 1")) # ~ X_{K} - 1

  X_form       <- formula(Form, lhs = 0, rhs = 1) # ~ X1 + .... X_{K - 1}
  X_vars       <- all.vars(X_form)

  list(Form = Form,
       outcome_var = outcome_var,
       outcome_form = outcome_form,
       X_form = X_form,
       X_vars = X_vars)
}


#' Tidy joint probabilities
#'
#'
#' @details Internal function to predict from a emlogit/bmlogit object, then reshape to
#' the population of interest. See `synth_mlogit()` and `synth_bmlogit()`
#'
#' @param fit A model of class bmlogit/emlogit
#' @param outcome_names A character vector of names that correspond to each level
#'  of the outcome of the multinomial. Must be manually set to the same ordering as
#'  the fitted objects.
#'
#' @inheritParams synth_mlogit
#' @importFrom rlang :=
#' @import emlogit
#' @import bmlogit
#' @keywords internal
#'
predict_longer <- function(fit, poptable, microdata, X_form, X_vars, area_var, count_var, outcome_var) {

  # Data for area var {A, X_{1}, ..., X_{K-1}}
  X_pred_df  <- collapse_table(
    poptable,
    area_var = area_var, X_vars = X_vars, count_var = count_var,
    new_name = count_var
  ) %>%
    group_by(!!sym(area_var)) %>%
    mutate(prX = !!sym(count_var) / sum(!!sym(count_var))) %>%
    ungroup()

  X_p_mat <- stats::model.matrix(X_form, X_pred_df)[, -1]

  # predicted values
  pred_X_p <- predict(fit, newdata = X_p_mat)

  out <- as_tibble(pred_X_p) %>%
    bind_cols(X_pred_df) %>%
    pivot_longer(cols = -c(X_vars, area_var, count_var, "prX"),
                 names_to = outcome_var,
                 values_to = "prZ_givenX") %>%
    mutate(prXZ = .data$prX * .data$prZ_givenX,
           !!sym(count_var) := !!sym(count_var)*.data$prZ_givenX) %>%
    dplyr::relocate(!!sym(count_var), .after = dplyr::last_col())

  # if original factor, make it back into a factor
  # (it was deconstructed in stats::model.matrix)
  if (inherits(microdata[[outcome_var]], "factor")) {
    out[[outcome_var]] <- factor(out[[outcome_var]],
                                 levels = levels(microdata[[outcome_var]]))

  }
  out
}
