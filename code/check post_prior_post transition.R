library(tidyverse)
library(lubridate)
library(tidybayes)
library(stemr)
library(fs)
library(extraDistr)
library(scales)
source("code/stemr_functions.R")

compartments <- c("S", "E", "Ie", "Ip", "R", "D")

tmp <- tibble(path = dir_ls("code/results_10x")) %>%
  mutate(folder = path_file(path)) %>%
  separate(col = folder, into = c("first_day", "last_day"), sep = "_") %>%
  mutate(multi_chain_stem_fit = map(path, ~read_rds(fs::path(., "original.rds")))) %>%
  tail(2)

popsize <- 3175692

samples_to_plot <- bind_rows(
  extract_epi_curves(multi_chain_stem_fit = tmp[[1,"multi_chain_stem_fit"]][[1]]) %>%
    filter(time == max(time)) %>%
    select(!!compartments) %>%
    nest(samples = everything()) %>%
    mutate(source = "old post"),
  rdirmnom(n = 8000, size = tmp[[2,"multi_chain_stem_fit"]][[1]]$stem_fit_list[[1]]$dynamics$popsize, alpha = tmp[[2,"multi_chain_stem_fit"]][[1]]$stem_fit_list[[1]]$dynamics$initdist_priors) %>%
    `colnames<-`(compartments) %>%
    as_tibble() %>%
    nest(samples = everything()) %>%
    mutate(source = "new prior"),
  extract_epi_curves(multi_chain_stem_fit = tmp[[2,"multi_chain_stem_fit"]][[1]]) %>%
    filter(time == min(time)) %>%
    select(!!compartments) %>%
    nest(samples = everything()) %>%
    mutate(source = "new post"))

# vglm(cbind(S, E, Ie, Ip, R, D) ~ 1, dirmultinomial,
#      data = samples_to_plot[[1,1]][[1]] %>%
#        round() %>%
#        mutate(S = popsize - (E + Ie + Ip + R + D)), trace = TRUE)

samples_to_plot %>%
  mutate(source = fct_inorder(source)) %>%
  unnest(samples) %>%
  pivot_longer(-source) %>%
  mutate(name = fct_inorder(name)) %>%
  ggplot(aes(source, value)) +
  facet_wrap(. ~ name, scales = "free_y") +
  stat_eye(normalize = "xy") +
  cowplot::theme_cowplot() +
  scale_y_continuous(labels = comma)


tmp2 <- tibble(path = dir_ls("code/results")) %>%
  mutate(folder = path_file(path)) %>%
  separate(col = folder, into = c("first_day", "last_day"), sep = "_") %>%
  mutate(multi_chain_stem_fit = map(path, ~read_rds(fs::path(., "original.rds")))) %>%
  mutate(prev_multi_chain_stem_fit = lag(multi_chain_stem_fit)) %>%
  slice(-1) %>%
  mutate(old_post = map(prev_multi_chain_stem_fit, ~extract_epi_curves(.) %>%
                               filter(time == max(time)) %>%
                               select(!!compartments))) %>%
  mutate(new_prior = map(multi_chain_stem_fit,
                         ~rdirmnom(n = 8000, size = .$stem_fit_list[[1]]$dynamics$popsize, alpha =.$stem_fit_list[[1]]$dynamics$initdist_priors) %>%
                           `colnames<-`(compartments) %>%
                           as_tibble())) %>%
  mutate(new_post = map(multi_chain_stem_fit, ~extract_epi_curves(.) %>%
                               filter(time == min(time)) %>%
                               select(!!compartments)))



# tmp2 %>%
#   select(first_day, old_post, new_prior, new_post) %>%
#   pivot_longer(-first_day) %>%
#   unnest(value) %>%
#   pivot_longer(cols = !!compartments, names_to = "compartment") %>%
#   ggplot(aes(name, value)) +
#   facet_grid(compartment ~ first_day) +
#   stat_eye(normalize = "xy") +
#   cowplot::theme_cowplot() +
#   scale_y_continuous(labels = comma)

tmp2 %>%
  select(first_day, old_post, new_prior, new_post) %>%
  pivot_longer(-first_day) %>%
  unnest(value) %>%
  pivot_longer(cols = !!compartments, names_to = "compartment") %>%
  mutate(name = str_replace(name, "_", "\n")) %>%
  mutate(first_day = fct_inorder(first_day),
         name = fct_inorder(name),
         compartment = fct_inorder(compartment)) %>%
  filter(compartment == "D") %>%
  ggplot(aes(name, value)) +
  facet_grid(. ~first_day) +
  stat_eye(normalize = "xy") +
  cowplot::theme_cowplot() +
  scale_y_continuous(labels = comma) +
  # coord_cartesian(ylim = c(0, 300000)) +
  xlab(NULL) +
  ylab(NULL)
