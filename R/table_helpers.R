#' Reduce the dimensionality of the table
#'
#' @param poptable The population table, collapsed in terms of counts. Must contain
#'  all variables in the RHS of `formula`, as well as the variables specified in
#'   `area_var` and `count_var` below.
#' @param area_var A character vector of the area of interest.
#' @param count_var A character variable that specifies which variable in `poptable`
#'  indicates the count
#' @param X_vars A character vector of the variables to keep. By design, this
#'  should be fewer variables than what is available in `poptable`.
#' @param report Should the output be in simple counts (the default, `"counts"`) or
#'  the proportion that count represents in the area (`"proportions"`)?
#' @param new_name What should the new count (or proportion) variable be called?
#' @param warnings Show some warnings about zero cells in the original data?
#'  Defaults to `FALSE`
#'
#' @returns A dataframe of counts or proportions. By design, it will include the
#'  variables `area_var`, `X_vars`, and whatever is specified in `new_name`.
#'
#' @importFrom dplyr count transmute group_by ungroup rename
#' @importFrom tidyr complete
#' @export
#'
#' @examples
#'
#'  # If you want to estimate education by female and age
#'  collapse_table(acs_race_NY, area_var = "cd", X_vars = c("female", "age"),
#'                 count_var = "count")
#'
#'  # Report proportions
#'  collapse_table(acs_race_NY, area_var = "cd", X_vars = c("female", "age"),
#'                 count_var = "count",
#'                 report = "proportions", new_name = "prop_in_cd")
#'
collapse_table <- function(poptable,
                           area_var,
                           X_vars,
                           count_var,
                           report = "counts",
                           new_name = "n_aggregate",
                           warnings = FALSE) {

  # collapse
  out <- count(poptable,
               !!!syms(c(area_var, X_vars)), #  by
               wt = !!sym(count_var), # sum the counts
               name = "new_count") %>%
    complete(!!!syms(c(area_var, X_vars)), fill = list(new_count = 0))

  if (any(out$new_count == 0) & warnings)
    warning("Some population combinations have zero people")

  # change into proportions
  if (report == "proportions") {
    if (grepl("(n_|count_|_n$|_count)", new_name)) {
      warning("You asked for proporitons but the value of `new_name` looks more appropriate for a count.")
    }

    out <- out %>%
      group_by(!!!syms(area_var)) %>%
      transmute(!!!syms(area_var),
                !!!syms(X_vars),
                !!sym(new_name) := new_count/sum(new_count)) %>% # turn into fraction
      ungroup()
  }

  if (report == "counts")
    out <- rename(out, !!sym(new_name) := "new_count")

  out
}
