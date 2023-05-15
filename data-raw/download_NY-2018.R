library(ccesMRPprep)
library(dplyr)
library(stringr)

acs_race_NY <- get_acs_cces(acscodes_age_sex_race) %>%
  filter(str_detect(cd, "NY")) %>%
  transmute(year, cd, female, race, age, count)

acs_educ_NY <- get_acs_cces(acscodes_age_sex_educ) %>%
  filter(str_detect(cd, "NY")) %>%
  transmute(year, cd, female, educ, age, count)

usethis::use_data(acs_educ_NY, overwrite = TRUE)
usethis::use_data(acs_race_NY, overwrite = TRUE)

