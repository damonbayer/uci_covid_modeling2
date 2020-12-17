library(stemr)
library(tidyverse)
library(tidybayes)
library(fs)
library(scales)
library(lubridate)
library(coda)
source("code/stemr_functions.R")

results_folder <- "code/results_5x"

tmp <- tibble(path = dir_ls(results_folder)) %>%
  mutate(folder = path_file(path)) %>%
  separate(col = folder, into = c("first_day", "last_day"), sep = "_") %>%
  mutate(across(ends_with("day"), ymd)) %>%
  mutate(multi_chain_stem_fit = map(path, ~read_rds(fs::path(., "original.rds"))))


popsize <- tmp[[1, "multi_chain_stem_fit"]][[1]]$stem_fit_list[[1]]$dynamics$popsize
ci_width <- c(0.5, 0.8, 0.95)

tmp_epi_curves <- tmp %>%
  mutate(epi_curves = map(multi_chain_stem_fit, extract_epi_curves)) %>%
  mutate(epi_curves = map2(epi_curves, multi_chain_stem_fit, ~left_join(.x, select(.y$data, time, date = end_date))))


# Epi Curves --------------------------------------------------------------
tmp_epi_curves %>%
  select(last_day, epi_curves) %>%
  unnest(epi_curves) %>%
  group_by(date) %>%
  filter(last_day == max(last_day)) %>%
  drop_na() %>%
  select(-starts_with("."), -time, -last_day) %>%
  pivot_longer(-date) %>%
  mutate(name = fct_inorder(name)) %>%
  group_by(date, name) %>%
  median_qi(.width = ci_width) %>%
  ggplot(aes(date, value, ymin = .lower, ymax = .upper)) +
  facet_wrap(. ~ name, scales = "free_y") +
  geom_lineribbon() +
  cowplot::theme_minimal_grid() +
  scale_fill_brewer(guide = NULL) +
  scale_y_continuous(label = comma, name = NULL) +
  scale_x_date(date_labels = "%b\n%d", name = NULL, date_breaks = "1 month") +
  ggtitle(results_folder)


# Prevalence --------------------------------------------------------------
tmp_epi_curves %>%
  select(first_day, last_day, epi_curves) %>%
  mutate(new_last_day = lead(first_day)) %>%
  unnest(epi_curves) %>%
  filter(date <= new_last_day) %>%
  mutate(prevalence = E + Ie + Ip) %>%
  select(prevalence, date) %>%
  group_by(date) %>%
  median_qi(.width = ci_width) %>%
  ggplot(aes(date, prevalence, ymin = .lower, ymax = .upper)) +
  geom_lineribbon() +
  cowplot::theme_minimal_grid() +
  scale_fill_brewer(guide = NULL) +
  scale_y_continuous(label = comma, name = "Prevalence") +
  scale_x_date(date_labels = "%b\n%d", name = NULL, date_breaks = "2 weeks") +
  ggtitle(results_folder)


# Seroprev ----------------------------------------------------------------
tmp_epi_curves %>%
  select(first_day, last_day, epi_curves) %>%
  mutate(new_last_day = lead(first_day)) %>%
  unnest(epi_curves) %>%
  filter(date <= new_last_day) %>%
  mutate(seroprevalence = R + D) %>%
  select(seroprevalence, date) %>%
  group_by(date) %>%
  median_qi(.width = ci_width) %>%
  ggplot(aes(date, seroprevalence, ymin = .lower, ymax = .upper)) +
  geom_lineribbon() +
  cowplot::theme_minimal_grid() +
  scale_fill_brewer(guide = NULL) +
  scale_y_continuous(label = comma, name = "seroprevalence (R + D)") +
  scale_x_date(date_labels = "%b\n%d", name = NULL, date_breaks = "2 weeks") +
  geom_segment(aes(x = ymd("2020-08-16"), xend = ymd("2020-08-16"),  y = 0.105 * popsize, yend = 0.125 * popsize), color = "red") +
  geom_point(aes(x = ymd("2020-08-16"), y = .115*popsize), color = "red") +
  ggtitle(results_folder)

tmp_epi_curves %>%
  select(first_day, last_day, epi_curves) %>%
  mutate(new_last_day = lead(first_day)) %>%
  unnest(epi_curves) %>%
  filter(date <= new_last_day) %>%
  mutate(seroprevalence = popsize - S) %>%
  select(seroprevalence, date) %>%
  group_by(date) %>%
  median_qi(.width = ci_width) %>%
  ggplot(aes(date, seroprevalence, ymin = .lower, ymax = .upper)) +
  geom_lineribbon() +
  cowplot::theme_minimal_grid() +
  scale_fill_brewer(guide = NULL) +
  scale_y_continuous(label = comma, name = "seroprevalence (pop - S)") +
  scale_x_date(date_labels = "%b\n%d", name = NULL, date_breaks = "2 weeks") +
  geom_segment(aes(x = ymd("2020-08-16"), xend = ymd("2020-08-16"),  y = 0.105 * popsize, yend = 0.125 * popsize), color = "red") +
  geom_point(aes(x = ymd("2020-08-16"), y = .115*popsize), color = "red") +
  ggtitle(results_folder)

# R0 ----------------------------------------------------------------------
tmp %>%
  mutate(post_params = map(multi_chain_stem_fit, ~extract_stem_parameter_posterior(.) %>% tidy_draws())) %>%
  unnest(post_params) %>%
  select(ends_with("day"), R0_est) %>%
  mutate(R0 = exp(R0_est)) %>%
  select(date = last_day, R0) %>%
  group_by(date) %>%
  median_qi(.width = ci_width) %>%
  ggplot(aes(date, R0, ymin = .lower, ymax = .upper)) +
  geom_lineribbon() +
  cowplot::theme_minimal_grid() +
  scale_fill_brewer(guide = NULL) +
  scale_x_date(date_labels = "%b\n%d", name = NULL, date_breaks = "2 weeks") +
  ggtitle(results_folder)


# IFR ---------------------------------------------------------------------
tmp %>%
  mutate(post_params = map(multi_chain_stem_fit, ~extract_stem_parameter_posterior(.) %>% tidy_draws())) %>%
  unnest(post_params) %>%
  select(ends_with("day"), ifr_est) %>%
  mutate(ifr = expit(ifr_est)) %>%
  select(date = last_day, ifr) %>%
  group_by(date) %>%
  median_qi(.width = ci_width) %>%
  ggplot(aes(date, ifr, ymin = .lower, ymax = .upper)) +
  geom_lineribbon() +
  cowplot::theme_minimal_grid() +
  scale_fill_brewer(guide = NULL) +
  scale_x_date(date_labels = "%b\n%d", name = NULL, date_breaks = "2 weeks") +
  scale_y_continuous(labels = percent) +
  ggtitle(results_folder)



# PP ----------------------------------------------------------------------
pp_function <- function(multi_chain_stem_fit) {
  imap_dfr(multi_chain_stem_fit$stem_fit_list,  ~{
    .x$results$posterior$latent_paths %>%
      gather_array(value, row, name, .iteration) %>%
      as_tibble() %>%
      mutate(name = dimnames(.x$results$posterior$latent_paths)[[2]][name],
             time = .x$results$posterior$latent_paths[,1,1][row]) %>%
      filter(name != "time") %>%
      select(.iteration, time, name, value, -row) %>%
      pivot_wider(names_from = name, values_from = value) %>%
      mutate(.chain = .y)
  }) %>%
    mutate(.draw = tidybayes:::draw_from_chain_and_iteration_(.chain, .iteration)) %>%
    select(.chain, .iteration, .draw, everything()) %>%
    left_join(tidy_draws(extract_stem_parameter_posterior(multi_chain_stem_fit = multi_chain_stem_fit, transform = "natural"))) %>%
    filter(time != 0) %>%
    left_join(select(multi_chain_stem_fit$dat, time, tests)) %>%
    mutate(tests = tests,
           alpha = kappa * (exp(alpha0) * (Ie2Ip / popsize) ^ alpha1) / (exp(alpha0) * (Ie2Ip / popsize) ^ alpha1 + (1 - Ie2Ip/popsize) ^ alpha1),
           beta = kappa * ((1 - Ie2Ip/popsize) ^ alpha1) / (exp(alpha0) * (Ie2Ip / popsize) ^ alpha1 + (1 - Ie2Ip/popsize) ^ alpha1),
           size = phi_death,
           mu = rho_death * Ip2D) %>%
    mutate(
      cases = pmap_dbl(
        list(tests = as.list(tests),
             alpha = as.list(alpha),
             beta = as.list(beta)),
        function(tests, alpha, beta) rbbinom(n = 1, size = tests, alpha = alpha, beta = beta)),
      deaths = pmap_dbl(
        list(mu = as.list(mu),
             size = as.list(size)),
        function(mu, size) rnbinom(n = 1, size = size, mu = mu))
    ) %>%
    select(starts_with("."), time, tests, cases, deaths) %>%
    mutate(pos = cases / tests) %>%
    select(.chain, time, deaths, pos) %>%
    pivot_longer(-time)
}

tmp_pp <- tmp %>%
  mutate(pp_dat = map(multi_chain_stem_fit, pp_function))

tmp_pp
