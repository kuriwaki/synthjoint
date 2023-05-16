#' Synthetic joint estimation with balancing constraint
#'
#' Imputes cells with a balancing constraint, using Yamauchi's algorithm.
#'
#' @param fix_to A dataset with only marginal counts or proportions of the outcome
#'  in question, by each area. Proportions will be corrected so that the margins
#'  of the synthetic joint will match these, with a simple ratio.
#' @param fix_by_area logical, whether to fix to targets area by area. Defaults to
#'  TRUE if `area_var` is a variable in `fix_to`. If `FALSE`, collapses the
#'  input to a single target.
#' @param tol Tolerance for balance
#'
#' @inheritParams synth_mlogit
#' @seealso `synth_mlogit()`
#'
#' @source
#'   Soichiro Yamauchi and Shiro Kuriwaki (2021). bmlogit: Multinomial logit with
#'   balancing constraints. R package version 0.0.3.
#'
#' @importFrom bmlogit bmlogit
#' @importFrom furrr future_map_dfr
#' @importFrom dplyr progress_estimated select filter
#' @importFrom tibble deframe
#' @examples
#'
#' library(dplyr)
#' library(ccesMRPprep)
#'
#' # can take a few minutes if fix_by_area = TRUE (the default)
#' educ_target <- count(acs_educ_NY, cd, educ, wt = count, name = "count")
#'
#' educ_target
#' acs_race_NY
#'
#' pop_syn <- synth_bmlogit(educ ~ race + age + female,
#'                          microdata = cc18_NY,
#'                          fix_to = educ_target,
#'                          poptable = acs_race_NY,
#'                          area_var = "cd")
#' pop_syn
#'
#' @export
synth_bmlogit <- function(formula,
                          microdata,
                          poptable,
                          fix_to,
                          fix_by_area = any(area_var %in% colnames(fix_to)),
                          area_var,
                          count_var = "count",
                          tol = 0.05) {
  # formula setup
  list2env(formula_parts(formula), envir = environment())

  # checks
  stopifnot(all(c(outcome_var, X_vars) %in% colnames(microdata)))
  stopifnot(all(c(area_var, X_vars) %in% colnames(poptable)))

  # Drop NAs
  microdata <- select(microdata, !!!syms(c(outcome_var, X_vars)))
  if (nrow(microdata) > sum(stats::complete.cases(microdata))) {
    warning("NAs in the microdata -- dropping data")
    microdata <- filter(microdata, stats::complete.cases(microdata))
  }

  # microdata ----
  # ys (in microdata)
  y_m_mat <- stats::model.matrix(outcome_form, microdata)

  # binary data
  if (all(microdata[[outcome_var]] %in% c(0, 1))) {
    y_m_mat <- cbind(1 - y_m_mat, y_m_mat)
    colnames(y_m_mat) <- c("0", "1")
  }

  if (is.factor(microdata[[outcome_var]]))
    colnames(y_m_mat) <- levels(microdata[[outcome_var]])

  # Xs setup microdata
  X_m_mat <- stats::model.matrix(X_form, microdata)[, -1]

  # population ----
  # Xs setup population table -- aggregate up to {X_1, ..., X_{K -1 }}
  X_p_df  <- collapse_table(poptable, area_var, X_vars, count_var, new_name = "N_X")
  X_p_mat <- stats::model.matrix(X_form, X_p_df)[, -1]

  # Ns of the Xs
  X_counts_vec <- X_p_df[["N_X"]]



   # when you wait to the aggregate thing

   if (isFALSE(fix_by_area)) {
     outcome_df <- collapse_table(
       fix_to,
       area_var = NULL, # only place that differs from other case
       X_vars = outcome_var,
       count_var = count_var,
       report = "proportions",
       new_name = "pr_outcome_tgt")

     pr_outcome_tgt <- outcome_df %>%
       select(!!sym(outcome_var), pr_outcome_tgt) %>%
       deframe()

     fit <- bmlogit(
       Y = y_m_mat,
       X = X_m_mat,
       target_Y = pr_outcome_tgt, # vector
       pop_X = X_p_mat, # matrix
       count_X = X_counts_vec, # vector
       control = list(tol_pred = tol)
     )

     out <- predict_longer(fit,
                           poptable = poptable,
                           microdata = microdata,
                           X_form = X_form,
                           X_vars = X_vars,
                           area_var = area_var,
                           count_var = count_var,
                           outcome_var = outcome_var)
   }


   # area by area, loop
   if (isTRUE(fix_by_area)) {
     outcome_df <- collapse_table(
       fix_to,
       area_var = area_var,
       X_vars = outcome_var,
       count_var = count_var,
       report = "proportions",
       new_name = "pr_outcome_tgt")

     areas <- unique(outcome_df[[area_var]])
     pb <- progress_estimated(length(areas))

     out <- future_map_dfr(
       .x = areas,

       .f = function(a) {
         outcome_df <- collapse_table(
           fix_to,
           area_var = area_var, # only place that differs from other case
           X_vars = outcome_var,
           count_var = count_var,
           report = "proportions",
           new_name = "pr_outcome_tgt")

         # overwrite
         outcome_A <- filter(outcome_df, !!sym(area_var) == a)
         pr_outcome_tgt <- outcome_A %>%
           select(!!sym(outcome_var), pr_outcome_tgt) %>%
           deframe()

         # overwrite this to area subset
         X_p_df  <- filter(X_p_df, !!sym(area_var) == a)
         X_p_mat <- stats::model.matrix(X_form, X_p_df)[, -1]
         X_counts_vec <- X_p_df[["N_X"]]

         # fit the model
         fit <- bmlogit(
           Y = y_m_mat,
           X = X_m_mat,
           target_Y = pr_outcome_tgt,
           pop_X   = X_p_mat,
           count_X = X_counts_vec,
           control = list(tol_pred = tol)
         )
         pb$tick()$print()

         # give it only the area population
         out <- predict_longer(fit,
                               poptable = filter(poptable, !!sym(area_var) == a),
                               microdata = microdata,
                               X_form = X_form,
                               X_vars = X_vars,
                               area_var = area_var,
                               count_var = count_var,
                               outcome_var = outcome_var)
       }
     ) # end map_dfr
   }

   out
}

