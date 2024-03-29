% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/synth_rake.R
\name{synth_prod}
\alias{synth_prod}
\title{Model synthetic joint distribution with simple product}
\usage{
synth_prod(formula, poptable, newtable, area_var, count_var = "count")
}
\arguments{
\item{formula}{A representation of the aggregate imputation or "outcome" model,
of the form \verb{X_\{K\} ~ X_1 + ... X_\{K - 1\}}}

\item{poptable}{The population table, collapsed in terms of counts. Must contain
all variables in the RHS of \code{formula}, as well as the variables specified in
\code{area_var} and \code{count_var} below.}

\item{newtable}{A dataset that contains marginal counts or proportions. Will be
collapsed internally to get simple proportions.}

\item{area_var}{A character vector of the area of interest.}

\item{count_var}{A character variable that specifies which variable in \code{poptable}
indicates the count}
}
\description{
Estimates joint distribution by simply assuming independence and multiplying
proportions.
}
\details{
That is, we already know \verb{p(X_\{1\}, ..., X_\{K - 1\}, A)} from \code{poptable}
and a marginal \verb{p(X_\{K\}, A)}from the additional distribution to weight to. Then
\verb{p(X_\{1\}, .., X_\{K - 1\}, X_\{K\}, A) = p(X_\{1\}, ..., X_\{K - 1\}, A) x p(X_\{K\}, A)}.
}
\examples{

library(dplyr)
library(ccesMRPprep)

# suppose we want know the distribution of (age x female) and we know the
# distribution of (race), by CD, but we don't know the joint of the two.

race_target <- count(acs_race_NY, cd, race, wt = count, name = "count")

pop_prod <- synth_prod(race ~ age + female,
                       poptable = acs_race_NY,
                       newtable = race_target,
                       area_var = "cd")

# In this example, we know the true joint. Does it match?
pop_val <- left_join(pop_prod,
                     count(acs_race_NY,  cd, age, female, race, wt = count, name = "count"),
                     by = c("cd", "age", "female", "race"),
                     suffix = c("_est", "_truth"))

# AOC's district in the bronx
pop_val \%>\%
  filter(cd == "NY-14", age == "35 to 44 years", female == 0) \%>\%
  select(cd, race, count_est, count_truth)

}
\seealso{
\code{synth_mlogit()} for a more nuanced model that uses survey data as
the basis of the joint estimation.
}
