library(tidyverse)
library(lubridate)
library(fs)
library(coda)
library(tidybayes)
library(zoo)
source("code/stemr_functions.R")

old_data_for_calcat <- read_csv("data/data_for_calcat.csv")
oc_data <- read_csv("data/oc_data.csv")

results_folder <-
  tibble(path = dir_ls("code/results")) %>%
  mutate(folder = path_file(path)) %>%
  separate(folder, c("start_date", "end_date"), sep = "_", remove = F) %>%
  mutate(across(ends_with("date"), ymd)) %>%
  tail(1) %>%
  pull(path)

ci_width <- 0.95

first_day <- ymd(str_sub(results_folder, start = -21, end = -12))
last_day <- ymd(str_sub(results_folder, start = -10))

multi_chain_stem_fit <- read_rds(path(results_folder, "original", ext = "rds"))
popsize <- multi_chain_stem_fit$stem_fit_list[[1]]$dynamics$popsize
forecast_obj <- read_rds(path(results_folder, "forecast", ext = "rds"))

epi_curves <- forecast_obj$forecast_results %>%
  select(starts_with("."), natural_paths) %>%
  mutate(natural_paths = map(natural_paths, as_tibble)) %>%
  unnest(natural_paths) %>%
  filter(time != 0)

params <- tidy_draws(extract_stem_parameter_posterior(multi_chain_stem_fit, to_human_scale))

deaths_at_t0 <- sum(oc_data[oc_data$date <= first_day, "deaths_calcat"])

new_cummort <-
  forecast_obj$forecast_results %>%
  select(starts_with("."), datasets) %>%
  mutate(datasets = map(datasets, as_tibble)) %>%
  unnest(datasets) %>%
  group_by(.draw) %>%
  mutate(cumulative_deaths = cumsum(deaths) + deaths_at_t0) %>%
  select(time, cumulative_deaths) %>%
  group_by(time) %>%
  median_qi(.width = ci_width) %>%
  left_join(select(forecast_obj$data, time, date = end_date)) %>%
  select(-time) %>%
  select(date, cummort_mean = cumulative_deaths, cummort_CI95l = .lower, cummort_CI95u = .upper)

new_re <-
  epi_curves %>%
  mutate(prop_S = S / popsize) %>%
  select(starts_with("."), time, prop_S) %>%
  left_join(select(params, starts_with("."), R0)) %>%
  mutate(Reff = prop_S * R0) %>%
  select(time, Reff) %>%
  group_by(time) %>%
  median_qi(.width = ci_width) %>%
  left_join(select(forecast_obj$data, time, date = end_date)) %>%
  select(date, re_mean = Reff, re_CI95l = .lower, re_CI95u = .upper)

new_data_for_calcat <-
  full_join(new_cummort, new_re) %>%
  right_join(., tibble(date = seq(min(.$date), max(.$date), 1))) %>%
  arrange(date) %>%
  pivot_longer(-date) %>%
  arrange(name, date) %>%
  mutate(value = na.approx(value, maxgap = 2)) %>%
  pivot_wider(names_from = name, values_from = value) %>%
  mutate(county = "Orange County",
         fips = 6059) %>%
  select(names(old_data_for_calcat))

updated_data_for_calcat <-
  old_data_for_calcat %>%
  filter(date < min(new_data_for_calcat$date)) %>%
  bind_rows(new_data_for_calcat)

write_csv(updated_data_for_calcat, "data/data_for_calcat.csv")

# new_data_for_calcat %>%
#   pivot_longer(-date) %>%
#   separate(col = name, sep = "_", into = c("stat", "type")) %>%
#   ggplot(aes(date, value, group = type, color = stat, linetype = type)) +
#   facet_wrap(. ~ stat, scales = "free_y") +
#   geom_line()
#
# updated_data_for_calcat %>%
#   select(-county, -fips) %>%
#   pivot_longer(-date) %>%
#   separate(col = name, sep = "_", into = c("stat", "type")) %>%
#   ggplot(aes(date, value, group = type, color = stat, linetype = type)) +
#   facet_wrap(. ~ stat, scales = "free_y") +
#   geom_line()
