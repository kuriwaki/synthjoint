
#' Model synthetic joint distribution with multinomial logit
#'
#' Imputes the counts of a joint distribution of count variables for small areas
#' based on microdata. `synth_smoothfix()` is basically `synth_mlogit()` which is then corrected by
#' if one knows the marginals.See Details.
#'
#' @details
#'
#' In this setup, the population distribution table (`poptable`)
#' has the joint distribution of `(A, X_{1}, ..., X_{K - 1})` categorical variables where
#' `A` denotes a categorical small area, `X`s denote categorical covariates, and
#' the missing covariate is `X_{K}`.
#'
#' Now, the survey data (`microdata`) has a sample joint distribution of
#' `(X_{1}, .., X_{K - 1}, X_{K})` categorical variables but the sample size is too small
#' for small areas. Therefore, the function models a multinomial outcome model roughly of the
#' form `X_{K} ~ X_{1} + ... X_{K}` and predicts onto `poptable` to estimate the
#' joint distribution of `(A, X_{1}, .., X_{K})`
#'
#' Currently, this function does not support post-stratifiation based on a known aggregate
#' distribution -- that is, further adjusting the probabilities based on
#' a known population distribution (see e.g. Leeman and Wasserfallen AJPS <https://doi.org/10.1111/ajps.12319>)
#'
#' @param formula A representation of the aggregate imputation or "outcome" model,
#'  of the form `X_{K} ~ X_1 + ... X_{K - 1}`
#' @param microdata The survey table that the multinomial model will be built off.
#'  Must contain all variables in the LHS and RHS of `formula`.
#' @param poptable The population table, collapsed in terms of counts. Must contain
#'  all variables in the RHS of `formula`, as well as the variables specified in
#'   `area_var` and `count_var` below.
#' @param area_var A character vector of the area of interest.
#' @param count_var A character variable that specifies which variable in `poptable`
#'  indicates the count
#'
#' @returns A dataframe with a similar format as the `poptarget` table
#'  but with rows expanded to serve as a joint distribution. In general,
#'  if the variable of interest has `L` values, the final dataset will have
#'  `L` times more rows than `poptarget`. The data will have additional variables:
#'
#'  * The outcome variable of interest (`Z`). For example if the LHS of the formula
#'     was `party_id`, then there would be a column called `party_id` containing the
#'     values of that variable in long form.
#'  * `prX`: The known distribution of the covriates (RHS) within the area, i.e.
#'     `Pr(X | A)`.
#'  * `prZ_givenX`: The main estimate from the multinomial logit model. Formally,
#'     `Pr(Z | X, A)`, although this is usually the same value for every `A`
#'     and thus equal to `Pr(Z | X)` unless `A` is on the RHS as well.
#'  * `prXZ`: A new estimate for the joint distribution within the area, i.e.
#'     `Pr(Z, X | A)`. Computed by `prX * prZ_givenX`.
#'  * `count`: A new count variable. Simply the product of `n_aggregate` and `pr_pred`.
#'
#'
#' @export
#'
#' @importFrom emlogit emlogit
#' @importFrom dplyr bind_cols as_tibble count sym syms left_join select mutate matches
#' @importFrom tidyr pivot_longer
#' @importFrom stats complete.cases
#' @importFrom rlang :=
#'
#' @source
#'  Jonathan P. Kastellec, Jeffrey R. Lax, Michael Malecki, and Justin H.
#'   Phillips (2015). Polarizing the electoral connection: Partisan representation in
#'   Supreme Court confirmation politics. _The Journal of Politics_ 77:3 <http://dx.doi.org/10.1086/681261>
#'
#'  Soichiro Yamauchi (2021). emlogit: Implementing the ECM algorithm for multinomial logit
#'   model. R package version 0.1.1. <https://github.com/soichiroy/emlogit>
#'
#'  Yair Ghitza and Mark Steitz (2020). DEEP-MAPS Model of the Labor Force.
#'    Working Paper. <https://github.com/Catalist-LLC/unemployment>
#'
#'
#'
#' @examples
#'  library(dplyr)
#'  library(ccesMRPprep)
#'
#'  # Impute the joint distribution of party ID with race, sex, and age, using
#'  # survey data in NY.
#'
#'  synth_acs <- synth_mlogit(pid3 ~ race + age + female,
#'                            microdata = cc18_NY,
#'                            poptable = acs_race_NY,
#'                            area_var = "cd")
#'
#'  # original (27 districts x 2 sex x 5 age x 6 race categories)
#'  count(acs_race_NY, cd, female, age, race, wt = count)
#'
#'  # new, modeled (original x 5 party categories)
#'  synth_acs
#'
#'  # See the data elec_NY to see if these numbers look reasonable.
#'
#'
#'\dontrun{
#'   # another example -- imputing education
#'   library(ccesMRPrun)
#'   synth_mlogit(educ ~ age + female,
#'               microdata = cces_GA,
#'               poptable = acs_GA,
#'               area_var = "cd")
#'}
#'
synth_mlogit <- function(formula,
                         microdata,
                         poptable,
                         area_var,
                         count_var = "count") {
  # formula setup
  list2env(formula_parts(formula), envir = environment())

  # checks
  stopifnot(all(c(outcome_var, X_vars) %in% colnames(microdata)))
  stopifnot(all(c(area_var, X_vars) %in% colnames(poptable)))

  # Drop NAs
  microdata <- select(microdata, !!!syms(c(outcome_var, X_vars)))
  if (nrow(microdata) > sum(complete.cases(microdata))) {
    warning("NAs in the microdata -- dropping data")
    microdata <- filter(microdata, complete.cases(microdata))
  }


  # ys (in microdata)
  y_m_mat <- stats::model.matrix(outcome_form, microdata)
  colnames(y_m_mat) <- levels(microdata[[outcome_var]])

  # Xs setup microdata
  X_m_mat <- stats::model.matrix(X_form, microdata)[, -1]

  # Xs setup population table -- aggregate up to {A, X_1, ..., X_{K -1 }}
  X_p_df <- collapse_table(poptable, area_var, X_vars, count_var)

  X_p_mat <- stats::model.matrix(X_form, X_p_df)[, -1]


  # fit model
  fit <- emlogit(Y = y_m_mat, X = X_m_mat)

  # predict and get predictions, formats
  predict_longer(fit,
                 poptable = poptable,
                 microdata = microdata,
                 X_form = X_form,
                 X_vars = X_vars,
                 area_var = area_var,
                 count_var = count_var,
                 outcome_var = outcome_var)
}



#' @rdname synth_mlogit
#'
#' @param fix_to A dataset with only marginal counts or proportions of the outcome
#'  in question, by each area. Proportions will be corrected so that the margins
#'  of the synthetic joint will match these, with a simple ratio.
#'
#' @examples
#'
#' # synth_mlogit WITH MARGINS CORRECTION -----
#'
#' library(dplyr)
#' # suppose we want know the distribution of (age x female) and we know the
#' # distribution of (race), by CD, but we don't know the joint of the two.
#'
#' educ_target <- count(acs_educ_NY, cd, educ, wt = count, name = "count")
#' pop_syn <- synth_smoothfix(educ ~ race + age + female,
#'                          microdata = cc18_NY,
#'                          fix_to = educ_target,
#'                          poptable = acs_race_NY,
#'                          area_var = "cd")
#'
#'
#' @export
synth_smoothfix <- function(formula,
                            microdata,
                            poptable,
                            fix_to,
                            area_var,
                            count_var = "count") {
  # formula setup
  list2env(formula_parts(formula), envir = environment())

  # estimate cells
  smooth_tbl <- synth_mlogit(formula, microdata, poptable, area_var, count_var)

  change_agg <- rake_spf(outcome_var = outcome_var,
                         data = smooth_tbl,
                         fix_to = fix_to,
                         area_var, count_var)

  # back to estimates
  smooth_tbl %>%
    left_join(change_agg, by = c(area_var, outcome_var)) %>%
    mutate(!!sym(count_var) := !!sym(count_var)*correction) %>%
    select(-correction) %>%
    select(-matches("prZ"), -matches("prXZ"))

  # fix margins
  # WHY DOES THIS GIVE SAME ANSWER AS THE PROD without survey
  # fixed_tbl <- synth_prod(formula,
  #                         poptable = smooth_agg,
  #                         newtable = fix_to,
  #                         area_var,
  #                         count_var = "n_aggregate")
}


