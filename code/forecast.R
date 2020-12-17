forecast_from_folder <- function(results_folder) {

  # Setup
  first_day <- ymd(str_sub(results_folder, start = -21, end = -12))
  last_day <- ymd(str_sub(results_folder, start = -10))
  last_forecast_day <- last_day + days(30)

  deaths_at_t0 <- read_csv("data/oc_data.csv") %>%
    select(date, deaths) %>%
    filter(date <= first_day) %>%
    pull(deaths) %>%
    sum()


  multi_chain_stem_fit <- read_rds(path(results_folder, "original", ext = "rds"))

  init_states <- multi_chain_stem_fit$stem_fit_list[[1]]$dynamics$initdist_params %>%
    enframe() %>%
    mutate(name = str_sub(name, end = -3)) %>%
    deframe()

  dat_original <- multi_chain_stem_fit$data

  date_delta <- dat_original$end_date[2] - dat_original$end_date[1]

  dat <- dat_original %>%
    add_row(tibble(end_date = seq(last_day + date_delta, last_forecast_day, by = date_delta))) %>%
    mutate(start_date = end_date - (date_delta - 1)) %>%
    mutate(time = as.numeric((end_date - min(start_date) + 1) / 7)) %>%
    replace_na(replace = list(cases = 0, tests = round(mean(dat_original$tests)), deaths = 0))

  # Build Model
  strata <- NULL
  compartments <- c("S", "E", "Ie", "Ip", "R", "D")
  obs_times <- dat$time

  rates <-
    list(
      rate(
        rate = "beta * (Ie + 0.8 * Ip)", # individual level rate (unlumped)
        from = "S", # source compartment
        to = "E", # destination compartment
        incidence = F
      ), # compute incidence of S2I transitions, required for simulating incidence data
      rate(
        rate = "gamma",
        from = "E",
        to = "Ie",
        incidence = F
      ),
      rate(
        rate = "nu_early",
        from = "Ie",
        to = "Ip",
        incidence = T
      ),
      rate(
        rate = "mu_rec",
        from = "Ip",
        to = "R",
        incidence = F
      ),
      rate(
        rate = "mu_death", # individual level rate
        from = "Ip", # source compartment
        to = "D", # destination compartment
        incidence = TRUE
      )
    ) # compute incidence of I2R transitions (not required for simulating data)

  state_initializer <-
    list(stem_initializer(
      init_states = init_states, # must match compartment names
      fixed = T
    )) # initial state fixed for simulation, we'll change this later


  parameters <- numeric(10)
  names(parameters) <- c("beta", "gamma", "nu_early", "mu_rec", "mu_death", "rho_death", "phi_death", "alpha0", "alpha1", "kappa")
  constants <- c(t0 = 0)
  tcovar <- data.frame(
    time = obs_times,
    tests = dat$tests
  )
  tmax <- max(tcovar$time)

  emissions <-
    list(
      emission(
        meas_var = "cases", # transition or compartment being measured (S->I transitions)
        distribution = "betabinomial", # emission distribution
        emission_params =
          c(
            "tests",
            "kappa * (exp(alpha0) * (Ie2Ip / popsize) ^ alpha1) / (exp(alpha0) * (Ie2Ip / popsize) ^ alpha1 + (1 - Ie2Ip/popsize) ^ alpha1)",
            "kappa * ((1 - Ie2Ip/popsize) ^ alpha1) / (exp(alpha0) * (Ie2Ip / popsize) ^ alpha1 + (1 - Ie2Ip/popsize) ^ alpha1)"
          ), # distribution pars, here overdispersion and mean
        incidence = TRUE, # is the data incidence
        obstimes = obs_times
      ), # vector of observation times
      emission(
        meas_var = "deaths",
        distribution = "negbinomial",
        emission_params = c(
          "phi_death",
          "rho_death * Ip2D"
        ),
        incidence = T,
        obstimes = obs_times
      )
    ) # vector of observation times)
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
      compile_ode = T, # compile ODE functions
      compile_rates = F, # compile MJP functions for Gillespie simulation
      compile_lna = T, # compile LNA functions
      messages = F # don't print messages
    )

  measurement_process <-
    stem_measure(
      emissions = emissions,
      dynamics = dynamics,
      data = dat %>%
        mutate(t = obs_times) %>%
        select(t, cases, deaths)
    )

  stem_object <- make_stem(dynamics = dynamics, measurement_process = measurement_process)

  # Forecast

  sim_results_tbl <- imap_dfr(multi_chain_stem_fit$stem_fit_list, function(stem_fit, chain) {
    simulation_parameters_list <- split_along_dim(stem_fit$results$posterior$parameter_samples_nat, 1)
    init_dist_list <- split_along_dim(stem_fit$results$posterior$initdist_samples, 1)

    map2_dfr(simulation_parameters_list, init_dist_list, function(simulation_parameters, init_dist) {
      R_plus_D <- init_dist[["R_0"]] + init_dist[["D_0"]]
      ifr <- expit(log(simulation_parameters[["mu_death"]]) - log(simulation_parameters[["mu_rec"]]))
      init_dist[["R_0"]] <- R_plus_D * (1 - ifr)
      init_dist[["D_0"]] <- R_plus_D * ifr

      stem_object$dynamics$initdist_params <- stem_object$dynamics$initdist_priors <- stem_object$dynamics$initializer[[1]]$init_states <- stem_object$dynamics$initializer[[1]]$prior <- init_dist

      stem_object$dynamics$parameters <- simulation_parameters
      tibble(sim_result = list(unlist(simulate_stem(stem_object = stem_object, method = "ode"), recursive = F)))
    }) %>%
      mutate(
        .chain = chain,
        .iteration = 1:multi_chain_stem_fit$n_iterations,
        .draw = tidybayes:::draw_from_chain_and_iteration_(chain = .chain, iteration = .iteration)
      )
  }) %>%
    unnest_wider(sim_result)

  # Return
  forecast_obj <- list(
    forecast_results = sim_results_tbl,
    data = dat %>%
      add_row(tibble(
        time = 0,
        start_date = min(dat$start_date) - date_delta,
        end_date = min(dat$end_date) - date_delta),
        .before = 1)
    )

  write_rds(forecast_obj, path = path(results_folder, "forecast", ext = "rds"))
}
