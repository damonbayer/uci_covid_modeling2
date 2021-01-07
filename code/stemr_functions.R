`%notin%` <- Negate(`%in%`)

# to_human_scale ----------------------------------------------------------
to_human_scale = function(params_est) {
  c(R0 = exp(params_est[["R0_est"]]),
    dur_latent = exp(-params_est[["dur_latent_est"]]),
    dur_early = exp(-params_est[["dur_early_est"]]),
    dur_progress = exp(-params_est[["dur_progress_est"]]),
    ifr = expit(params_est[["ifr_est"]]),
    rho_death = expit(params_est[["rho_death_est"]]),
    phi_death = exp(-2 * params_est[["phi_death_est"]]),
    alpha0 = exp(params_est[["alpha0_est"]]),
    alpha1 = expit(params_est[["alpha1_est"]]),
    kappa = exp(-2 * params_est[["kappa_est"]]))
}


# extract_stem_parameter_posterior ----------------------------------------
extract_stem_parameter_posterior <- function(multi_chain_stem_fit, transform = "estimation") {
  map(multi_chain_stem_fit$stem_fit_list,
      function(stem_fit) {
        if (is.character(transform) && transform == "estimation") {
          parameter_samples <- stem_fit$results$posterior$parameter_samples_est
        } else if (is.character(transform) && transform == "natural") {
          parameter_samples <- stem_fit$results$posterior$parameter_samples_nat
        } else {
          parameter_samples <- t(apply(stem_fit$results$posterior$parameter_samples_est, 1, transform))
        }
        # mcmc(parameter_samples, start = 100, end = 200000, thin = 100)}) %>%
        mcmc(parameter_samples)}) %>%
    as.mcmc.list()
}


# extract_epi_curves ------------------------------------------------------
extract_epi_curves <- function(multi_chain_stem_fit, curve_type = "prevalence", tidy = F) {
  if (curve_type == "prevalence") curve_type <- "p"
  if (curve_type == "incidence") curve_type <- "i"
  if (curve_type %notin% c("i", "p")) stop('curve_type must be one of "i", "incidence", "p", or "prevalence"')
  if (!is.logical(tidy)) stop('tidy must be TRUE or FALSE')

  epi_curves <- imap_dfr(multi_chain_stem_fit$stem_fit_list, function(stem_fit, chain) {
    latent_paths_list <- split_along_dim(stem_fit$result$posterior$latent_paths, 3)

    if (stem_fit$dynamics$fixed_inits) {
      init_dist_sample_list <- rep(list(stem_fit$dynamics$initdist_params), multi_chain_stem_fit$n_iterations)
    } else {
      init_dist_sample_list <- split_along_dim(stem_fit$results$posterior$initdist_samples, 1)
    }

    if (curve_type == "i") {
      path_df <- map_dfr(latent_paths_list, as_tibble)
    } else {
      path_df <- map2(.x = latent_paths_list, .y = init_dist_sample_list, .f =
                        ~incidence2prevalence(path = .x,
                                              flow_matrix = stem_fit$dynamics$flow_matrix_ode,
                                              init_state = .y)) %>%
        map_dfr(as_tibble)
    }

    path_df <- mutate(path_df, .iteration = rep(1:multi_chain_stem_fit$n_iterations, each = nrow(path_df) / multi_chain_stem_fit$n_iterations))
    path_df
  } %>%
    mutate(.chain = chain)
  ) %>%
    mutate(.draw = tidybayes:::draw_from_chain_and_iteration_(chain = .chain, iteration = .iteration))

  if (tidy == T) {
    epi_curves <- epi_curves %>%
      pivot_longer(-c(time, .iteration, .chain, .draw)) %>%
      mutate(name = fct_inorder(name))
  }
  epi_curves
}


# plot_epi_curves ---------------------------------------------------------
# possible add an option to separate by chain
plot_epi_curves <- function(multi_chain_stem_fit, curve_type = "p") {
  get_epi_curves(multi_chain_stem_fit = multi_chain_stem_fit, curve_type = curve_type, tidy = T) %>%
    ggplot(aes(time, value)) +
    stat_lineribbon() +
    facet_wrap(. ~ name, scales = "free_y") +
    scale_fill_brewer()
}
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


gather_array <- function(a, value, ..., .id=NULL){
  qs <- rlang::quos(...)
  if (missing(value)) {
    evalue <- rlang::sym("var")}
  else {
    evalue <- rlang::enquo(value)
  }
  len <- length(qs)
  d <- dim(a)

  # Default Values
  if (len > 0) {
    dimnames <- purrr::map(qs, rlang::quo_name) %>%
      as_vector()
  } else {
    dimnames <- paste0("dim_", 1:length(d))
  }

  l <- list()
  for (i in 1:length(d)){
    l[[i]] <- 1:d[i]
  }
  names(l) <- dimnames
  tidy <- expand.grid(l)
  tidy[[rlang::quo_name(evalue)]] <- a[as.matrix(tidy)]
  if (!is.null(.id)) tidy[[.id]] <- rlang::expr_name(a)
  return(tidy)
}

split_along_dim <- function(a, n){
  setNames(lapply(split(a, arrayInd(seq_along(a), dim(a))[, n]),
                  array, dim = dim(a)[-n], dimnames(a)[-n]),
           dimnames(a)[[n]])
}

dataset_to_list <- function(dat) {
  dat <- as.matrix(dat)
  time <- dat[,1]
  dat_na <- is.na(dat[,-1])
  unique_na_structure <- unique(dat_na, MARGIN = 2)

  assignments <- apply(dat_na, 2, function(x) which(apply(unique_na_structure, 2, function(y) isTRUE(all.equal(x, y)))))

  dat_list <- list()
  for (i in 1:length(unique(assignments))) {
    dat_list[[i]] <- dat[!unique_na_structure[,i], c(1, which(assignments == i) + 1)]
  }
  dat_list
}

named_group_split <- function(.tbl, ...) {
  grouped <- group_by(.tbl, ...)
  names <- rlang::eval_bare(rlang::expr(paste(!!!group_keys(grouped), sep = " / ")))

  grouped %>%
    group_split(.keep = F) %>%
    rlang::set_names(names)
}
