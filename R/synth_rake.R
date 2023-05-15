#' Model synthetic joint distribution with simple product
#'
#' Estimates joint distribution by simply assuming independence and multiplying
#' proportions.
#'
#' @details That is, we already know `p(X_{1}, ..., X_{K - 1}, A)` from `poptable`
#'  and a marginal `p(X_{K}, A)`from the additional distribution to weight to. Then
#'  `p(X_{1}, .., X_{K - 1}, X_{K}, A) = p(X_{1}, ..., X_{K - 1}, A) x p(X_{K}, A)`.
#'
#' @inheritParams synth_mlogit
#' @param newtable A dataset that contains marginal counts or proportions. Will be
#'  collapsed internally to get simple proportions.
#'
#' @seealso `synth_mlogit()` for a more nuanced model that uses survey data as
#'  the basis of the joint estimation.
#'
#' @importFrom dplyr full_join mutate
#' @export
#'
#' @examples
#'
#' library(dplyr)
#' library(ccesMRPprep)
#'
#' # suppose we want know the distribution of (age x female) and we know the
#' # distribution of (race), by CD, but we don't know the joint of the two.
#'
#' race_target <- count(acs_race_NY, cd, race, wt = count, name = "count")
#'
#' pop_prod <- synth_prod(race ~ age + female,
#'                        poptable = acs_race_NY,
#'                        newtable = race_target,
#'                        area_var = "cd")
#'
#' # In this example, we know the true joint. Does it match?
#' pop_val <- left_join(pop_prod,
#'                      count(acs_race_NY,  cd, age, female, race, wt = count, name = "count"),
#'                      by = c("cd", "age", "female", "race"),
#'                      suffix = c("_est", "_truth"))
#'
#' # AOC's district in the bronx
#' pop_val %>%
#'   filter(cd == "NY-14", age == "35 to 44 years", female == 0) %>%
#'   select(cd, race, count_est, count_truth)
#'
synth_prod <- function(formula,
                       poptable,
                       newtable,
                       area_var,
                       count_var = "count") {
  # formula setup
  list2env(formula_parts(formula), envir = environment())

  # N_Area
  N_area <- collapse_table(newtable, area_var, X_vars = NULL, count_var,
                           report = "counts", new_name = "N_area")

  # margins to start with
  X_p_df <- collapse_table(poptable, area_var, X_vars, count_var,
                           report = "proportions", new_name = "pr_Xs")

  # target proportions
  X_t_df <- collapse_table(newtable, area_var, X_vars = outcome_var, count_var,
                           report = "proportions", new_name = "pr_outcomes")

  # merge and get new counts
  N_area %>% # |A|
    full_join(X_p_df, by = area_var) %>% # |X1| * .... |X{K-1}|
    full_join(X_t_df, by = area_var) %>% # * |XK|
    mutate(!!sym(count_var) := N_area * pr_Xs * pr_outcomes)
}



#' Create long table for rake weighitng
#'
#'
#' @inheritParams synth_prod
#'
#' @importFrom dplyr full_join mutate
#' @importFrom stringr str_c
#'
#' @returns A long table that, for each districts, list the target distribution
#'  of the outcome (`variable = "outcome"`) and the cell of the known Xs (`variable = "Xs`).
#'  This is consistent with the dataframe format for targets in the
#'  `autumn` package <https://github.com/aaronrudkin/autumn>
#'
#' @export
#' @examples
#' library(dplyr)
#' library(ccesMRPprep)
#'
#' # suppose we want know the distribution of (age x female) and we know the
#' # distribution of (race), by CD, but we don't know the joint of the two.
#'
#' race_target <- count(acs_race_NY, cd, race, wt = count, name = "count")
#'
#' rake_target(race ~ age + female,
#'             poptable = acs_race_NY,
#'             newtable = race_target,
#'             area_var = "cd")
#'
rake_target <- function(formula,
                       poptable,
                       newtable,
                       area_var,
                       count_var = "count") {
  # formula setup
  list2env(formula_parts(formula), envir = environment())


  # aggregate this to the estimated X_{K}s
  Xs_agg <- collapse_table(poptable,
                           area_var = area_var,
                           X_vars = X_vars,
                           count_var = count_var,
                           report = "proportions",
                           new_name = "proportion") %>%
    transmute(!!sym(area_var),
              variable = "Xs",
              label = str_c(!!!syms(X_vars), sep = "_"),
              proportion
              )


  outcome_agg <- collapse_table(newtable,
                                area_var = area_var,
                                X_vars = outcome_var,
                                count_var = count_var,
                                report = "proportions",
                                new_name = "proportion") %>%
    transmute(!!sym(area_var),
              variable = "outcome",
              label = !!sym(outcome_var),
              proportion
              )

  # target
  tgt_stacked <- dplyr::bind_rows(Xs_agg, outcome_agg) %>%
    dplyr::arrange(!!sym(area_var), variable)


  # rake
  tgt_stacked
}

#' Single Proportional Fitting (as opposed to Iterative) with microdata
#'
#' With microdata or a table that includes both
#'
#' @param outcome_var A string for the variable to match the proportion to
#' @param data Data of microdata or tables whose margin on `outcom_var` must be fixed
#'
#'
#' @inheritParams synth_mlogit
#' @export
rake_spf <- function(outcome_var,
                     data,
                     fix_to,
                     area_var,
                     count_var = "count") {

  data_agg <- collapse_table(data,
                             area_var = area_var,
                             X_vars = outcome_var,
                             count_var = count_var,
                             report = "proportions",
                             new_name = "pr_outcome_data")

  target_agg <- collapse_table(fix_to,
                               area_var = area_var,
                               X_vars = outcome_var,
                               count_var = count_var,
                               report = "proportions",
                               new_name = "pr_outcome_tgt")

  # correction factor
  left_join(data_agg,
            target_agg,
            by = c(area_var, outcome_var)) %>%
    transmute(
      !!!syms(area_var),
      !!sym(outcome_var),
      correction = pr_outcome_tgt / pr_outcome_data
    )
}
