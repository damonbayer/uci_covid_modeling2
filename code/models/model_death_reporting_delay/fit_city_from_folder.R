library(tidyverse)
library(lubridate)
library(stemr)
library(extraDistr)
library(fs)
library(foreach)
library(doRNG)
library(doParallel)
registerDoParallel(cores = future::availableCores())
source('code/helper_functions.R')
source("code/stemr_functions.R")
source('code/forecast.R')
death_delay_ecdf <- read_rds("data/death_delay_ecdf.rds")

time_interval_in_days <- 3


first_day <- county_results_folder %>% path_file() %>% str_sub(end = 10)
last_day <- county_results_folder %>% path_file() %>% str_sub(start = 12)

county_multi_chain_stem_fit <- read_rds(path(county_results_folder, "original", ext = "rds"))

city_incid <- read_csv("data/oc_city_incidence.csv") %>%
  filter(city == city_name)

city_data <- read_csv("data/oc_city_data.csv") %>%
  filter(city == city_name)

county_popsize <- county_multi_chain_stem_fit$stem_fit_list[[1]]$dynamics$popsize
county_initdist <- county_multi_chain_stem_fit$stem_fit_list[[1]]$dynamics$initdist_params
county_initprior <- county_multi_chain_stem_fit$stem_fit_list[[1]]$dynamics$initdist_priors
county_C <- county_popsize / sum(county_initprior)

city_popsize <- city_incid$population

city_initidist <-
  c(S_0 = city_popsize - sum(county_initdist[-1] * city_incid$prop_incid),
    county_initdist[-1] * city_incid$prop_incid) %>%
  `names<-`(., str_sub(names(.), end = -3))

target_raw <- rdirmnom(n = 8000, size = county_popsize, alpha = county_initprior)
target <- cbind(city_popsize - rowSums(round(target_raw[,-1] * city_incid$prop_incid)),
                round(target_raw[,-1] * city_incid$prop_incid)) %>%
  `colnames<-`(names(city_initidist))


city_C <- optimize(f = function(C) sum(extraDistr::ddirmnom(x = target, size = city_popsize, alpha = city_initidist / C, log = T)), lower = 1, upper = 100000, maximum = T)$maximum


dat <- city_data %>%
  lump_oc_data(time_interval_in_days,
               first_day,
               last_day) %>%
  mutate(prop_deaths_reported = death_delay_ecdf(as.numeric(max(city_data$date) - end_date)))

init_states <- city_initidist
C <- city_C
popsize <- city_popsize

# -------------------------------------------------------------------------


strata <- NULL
compartments <- c("S", "E", "Ie", "Ip", "R", "D")
obs_times <- dat$time

rates <-
  list(rate(rate = "beta * (Ie + 0.8 * Ip)", # individual level rate (unlumped)
            from = "S",        # source compartment
            to   = "E",        # destination compartment
            incidence = F),    # compute incidence of S2I transitions, required for simulating incidence data
       rate(rate = "gamma",
            from = "E",
            to = "Ie",
            incidence = F),
       rate(rate = "nu_early",
            from = "Ie",
            to = "Ip",
            incidence = T),
       rate(rate = "mu_rec",
            from = "Ip",
            to = "R",
            incidence = F),
       rate(rate = "mu_death",       # individual level rate
            from = "Ip",        # source compartment
            to   = "D",        # destination compartment
            incidence = TRUE)) # compute incidence of I2R transitions (not required for simulating data)

state_initializer <-
  list(stem_initializer(
    init_states = init_states, # must match compartment names
    fixed = F,
    prior = init_states / C,
    dist = "dirmultinom"
  )) # initial state fixed for simulation, we'll change this later


parameters <- numeric(10); names(parameters) <- c("beta", "gamma", "nu_early", "mu_rec", "mu_death", "rho_death", "phi_death", "alpha0", "alpha1", "kappa")
constants <- c(t0 = 0)
tcovar <- data.frame(time = obs_times,
                     tests = dat$tests,
                     prop_deaths_reported = dat$prop_deaths_reported)
tmax <- max(tcovar$time)

emissions <-
  list(emission(meas_var = "cases", # transition or compartment being measured (S->I transitions)
                distribution    = "betabinomial",        # emission distribution
                emission_params =
                  c("tests",
                    "kappa * (exp(alpha0) * (Ie2Ip / popsize) ^ alpha1) / (exp(alpha0) * (Ie2Ip / popsize) ^ alpha1 + (1 - Ie2Ip/popsize) ^ alpha1)",
                    "kappa * ((1 - Ie2Ip/popsize) ^ alpha1) / (exp(alpha0) * (Ie2Ip / popsize) ^ alpha1 + (1 - Ie2Ip/popsize) ^ alpha1)"), # distribution pars, here overdispersion and mean
                incidence       = TRUE,                  # is the data incidence
                obstimes        = obs_times), # vector of observation times
       emission(meas_var = "deaths",
                distribution = "negbinomial",
                emission_params = c("phi_death",
                                    "rho_death * Ip2D * prop_deaths_reported"),
                incidence = T,
                obstimes        = obs_times)) # vector of observation times)
# list of emission distribution lists (analogous to rate specification)

dynamics <-
  stem_dynamics(
    rates = rates,
    parameters = parameters,
    state_initializer = state_initializer,
    compartments = compartments,
    constants = constants,
    tcovar = tcovar,
    tmax = tmax,
    compile_ode = T,   # compile ODE functions
    compile_rates = F, # compile MJP functions for Gillespie simulation
    compile_lna = T,   # compile LNA functions
    messages = F       # don't print messages
  )

measurement_process <-
  stem_measure(emissions = emissions,
               dynamics = dynamics,
               data = dat %>%
                 mutate(t = obs_times) %>%
                 select(t, cases, deaths))

stem_object <- make_stem(dynamics = dynamics, measurement_process = measurement_process)

# Build Priors ------------------------------------------------------------
to_estimation_scale = function(params_nat) {
  c(R0_est = log(params_nat[["beta"]]) + log(popsize) + log(1 / params_nat[["nu_early"]] + 0.8 / (params_nat[["mu_rec"]] + params_nat[["mu_death"]])), # log(R0)
    dur_latent_est = log(params_nat[["gamma"]]), # -log(dur_latent)
    dur_early_est = log(params_nat[["nu_early"]]), # -log(dur_early)
    dur_progress_est = log(params_nat[["mu_rec"]] + params_nat[["mu_death"]]), # -log(dur_progress)
    ifr_est = log(params_nat[["mu_death"]]) - log(params_nat[["mu_rec"]]), # logit(ifr)
    rho_death_est = logit(params_nat[["rho_death"]]), # logit(rho_death)
    phi_death_est = -0.5 * log(params_nat[["phi_death"]]), # -0.5 * log(phi_death)
    alpha0_est = log(params_nat[["alpha0"]]), # log(alpha0)
    alpha1_est = logit(params_nat[["alpha1"]]), # logit(alpha1)
    kappa_est = -0.5 * log(params_nat[["kappa"]])) # -0.5 * log(kappa)
}


from_estimation_scale = function(params_est) {
  c(beta = exp(params_est[["R0_est"]] - log(popsize) - log(exp(-params_est[["dur_early_est"]]) + 0.8 * exp(-params_est[["dur_progress_est"]]))),
    gamma = exp(params_est[["dur_latent_est"]]),
    nu_early = exp(params_est[["dur_early_est"]]),
    mu_rec = exp(params_est[["dur_progress_est"]]) / (1 + exp(params_est[["ifr_est"]])),
    mu_death = exp(params_est[["dur_progress_est"]]) / (1 + exp(-params_est[["ifr_est"]])),
    rho_death = expit(params_est[["rho_death_est"]]),
    phi_death = exp(-2 * params_est[["phi_death_est"]]),
    alpha0 = exp(params_est[["alpha0_est"]]),
    alpha1 = expit(params_est[["alpha1_est"]]),
    kappa = exp(-2 * params_est[["kappa_est"]]))
}


logprior =
  function(params_est) {
    sum(dnorm(params_est["R0_est"], -0.2554128198465173693599, 0.7, log = TRUE), # log(R0)
        dnorm(-params_est["dur_latent_est"], 0, 0.22, log = TRUE), # -log(dur_latent)
        dnorm(-params_est["dur_early_est"], 0, 0.22, log = TRUE), # -log(dur_early)
        dnorm(-params_est["dur_progress_est"], 0, 0.22, log = TRUE), # -log(dur_progress)
        # dbeta(expit(params_est["ifr_est"]), 1.5, 200, log = TRUE) + params_est["ifr_est"] - 2 * log(exp(params_est["ifr_est"]) + 1), # logit(ifr)
        dbeta(expit(params_est["ifr_est"]), 2, 460, log = TRUE) + params_est["ifr_est"] - 2 * log(exp(params_est["ifr_est"]) + 1), # logit(ifr)
        dbeta(expit(params_est["rho_death_est"]), 8, 2, log = TRUE) + params_est["rho_death_est"] - 2 * log(exp(params_est["rho_death_est"]) + 1) , # logit(rho_death)
        dexp(exp(params_est["phi_death_est"]), 1, log = TRUE) +  params_est["phi_death_est"], # -0.5 * log(phi_death)
        dtnorm(exp(params_est["alpha0_est"]), mean = 4, sd = 2, a = 0, log = TRUE) + params_est["alpha0_est"], # log(alpha0)
        dbeta(expit(params_est["alpha1_est"]), 3, 1, log = TRUE) + params_est["alpha1_est"] - 2 * log(exp(params_est["alpha1_est"]) + 1), # logit(alpha1)
        dexp(exp(params_est["kappa_est"]), 1, log = T) +  params_est["kappa_est"]) # -0.5 * log(kappa)
  }

priors <- list(logprior = logprior,
               to_estimation_scale = to_estimation_scale,
               from_estimation_scale = from_estimation_scale)

# Build par_initializer ---------------------------------------------------
true_pars =
  c(R0       = 1.5,    # basic reproduction number
    dur_latent = 1,
    dur_early   = 1,      # infectious period duration = 2 days
    dur_progress = 1,
    ifr = 0.06,
    rho_death = 0.7,
    phi_death = 2.2,
    alpha0   = 4, # beta-binomial intercept
    alpha1   = 0.8,    # beta-binomial slope
    kappa    = 2.2)

parameters =
  c(beta = true_pars[["R0"]] / popsize / (true_pars[["dur_early"]] + true_pars[["dur_progress"]]), # R0 = beta * P / mu
    gamma = 1 / true_pars[["dur_latent"]],
    nu_early = 1 / true_pars[["dur_early"]],
    mu_rec = (1 - true_pars[["ifr"]]) / true_pars[["dur_progress"]],
    mu_death = true_pars[["ifr"]] / true_pars[["dur_progress"]],
    rho_death = true_pars[["rho_death"]],
    phi_death = true_pars[["phi_death"]],
    alpha0 = true_pars[["alpha0"]],
    alpha1 = true_pars[["alpha1"]],
    kappa = true_pars[["kappa"]])

n_params <- length(parameters)
par_initializer = function() {
  priors$from_estimation_scale(priors$to_estimation_scale(parameters) + rnorm(n_params, 0, 0.1))
}


# Build Kernel ------------------------------------------------------------
# specify the kernel
mcmc_kern <-
  mcmc_kernel(
    parameter_blocks =
      list(parblock(
        pars_nat = c("beta", "gamma", "nu_early", "mu_rec", "mu_death", "rho_death", "phi_death", "alpha0", "alpha1", "kappa"),
        # par_est should have different names from pars_nat
        pars_est = c("R0_est", "dur_latent_est", "dur_early_est", "dur_progress_est", "ifr_est", "rho_death_est", "phi_death_est", "alpha0_est", "alpha1_est", "kappa_est"),
        priors = priors,
        # alg = "mvnss",
        alg = "mvnmh",
        sigma = diag(0.01, n_params),
        initializer = par_initializer,
        control =
          # mvnss_control(stop_adaptation = 1e2))),
          mvnmh_control(stop_adaptation = 150000,
                        scale_cooling = 0.85,
                        scale_constant = 1,
                        step_size = 0.25))),
    lna_ess_control = lna_control(bracket_update_iter = 5e3,
                                  joint_initdist_update = FALSE))


# Fit Model ---------------------------------------------------------------

n_chains <- 4
thinning_interval <- 100
iterations <- 350000

stem_fit_list <- foreach(chain = 1:n_chains,
                         .packages = "stemr",
                         .export = ls()) %dorng% {
                           fit_stem(stem_object = stem_object,
                                    method = "ode", # or "lna"
                                    mcmc_kern = mcmc_kern,
                                    iterations = iterations,
                                    thinning_interval = thinning_interval,
                                    print_progress = 1e3)
                         }

multi_chain_stem_fit <- list()
multi_chain_stem_fit$data <- dat
multi_chain_stem_fit$stem_fit_list <- stem_fit_list
multi_chain_stem_fit$n_iterations <- 2000
multi_chain_stem_fit$n_chains <- n_chains
multi_chain_stem_fit$thinning_interval <- thinning_interval

city_folder <- path(county_results_folder, str_replace_all(str_to_lower(city_name), "\\s", "-"))
dir_create(city_folder)

write_rds(multi_chain_stem_fit, path(city_folder, "original", ext = "rds"))

forecast_from_folder(results_folder = city_folder)
