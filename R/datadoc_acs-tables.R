#' Three-way ACS Tables for New York
#'
#'
#' Table of counts by CD from the 2018 one-year ACS in New York. We load the
#' variables that roughly partition gender, age, and race. We get another
#' three way table of gender, age, and education. That is, this is the
#' output of `ccesMRPprep::get_acs_cces(acscodes_age_sex_race)` and
#' `ccesMRPprep::get_acs_cces(acscodes_age_sex_educ)`.
#'
#' @rdname acs_NY
#' @source
#' Kyle Walker and Matt Herman (2021). tidycensus: Load US Census Boundary and Attribute Data as
#'  'tidyverse' and 'sf'-Ready Data Frames. R package version 0.11.4.
#'  <https://CRAN.R-project.org/package=tidycensus>
#'
#' @examples
#'  library(tibble)
#'  acs_race_NY
#'  acs_educ_NY
#'
"acs_educ_NY"

#' @rdname acs_NY
"acs_race_NY"
