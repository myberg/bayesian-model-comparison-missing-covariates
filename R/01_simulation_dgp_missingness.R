# -----------------------------------------------------------------------------
# 01_simulation_dgp_missingness.R
# DGP plus MCAR/MAR/MNAR mechanisms
# -----------------------------------------------------------------------------


# generate probit data
# setup is not indiviudally adjustable
# simulating routine discribed in paper
simulate_re_probit_data <- function(config, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  N <- config$N
  J <- config$J
  beta_true <- config$beta_true
  sigma2_a_true <- config$sigma2_a_true
  
  id <- sample(seq_len(J), size = N, replace = TRUE)
  
  # continuous primary covariates: X2, X3, X4
  mu <- c(0, 0, 0)
  Sigma <- matrix(
    c(
      1.00,  0.45, -0.15,
      0.45,  1.00,  0.25,
      -0.15,  0.25,  1.00
    ),
    nrow = 3,
    byrow = TRUE
  )
  
  X_cont <- MASS::mvrnorm(n = N, mu = mu, Sigma = Sigma)
  colnames(X_cont) <- c("x2", "x3", "x4")
  
  # categorical X5 with levels A, B, C
  x5_fac <- factor(
    sample(c("A", "B", "C"), size = N, replace = TRUE),
    levels = c("A", "B", "C")
  )
  
  x5B <- as.integer(x5_fac == "B")
  x5C <- as.integer(x5_fac == "C")
  
  # redundant continuous noise variable X7
  x7 <- rnorm(N, mean = 0, sd = 1)
  
  # categorical X8 with levels A, B, C
  x8_fac <- factor(
    sample(c("A", "B", "C"), size = N, replace = TRUE),
    levels = c("A", "B", "C")
  )
  
  x8B <- as.integer(x8_fac == "B")
  x8C <- as.integer(x8_fac == "C")
  
  # redundant continuous noise variable X10
  x10 <- rnorm(N, mean = 0, sd = 1)
  
  dat <- data.frame(
    id = factor(id),
    x2 = X_cont[, "x2"],
    x3 = X_cont[, "x3"],
    x4 = X_cont[, "x4"],
    x5B = x5B,
    x5C = x5C,
    x7 = x7,
    x8B = x8B,
    x8C = x8C,
    x10 = x10
  )
  
  # random effects with variance sigma2_a_true
  a_j <- rnorm(J, mean = 0, sd = sqrt(sigma2_a_true))
  
  X_true <- model.matrix(
    ~ x2 + x3 + x4 + x5B + x5C,
    data = dat
  )
  
  z <- as.numeric(X_true %*% beta_true) +
    a_j[id] +
    rnorm(N, mean = 0, sd = 1)
  
  y <- as.integer(z > 0)
  
  dat$y <- y
  dat$z <- z
  dat$a_true <- a_j[id]
  
  dat
}

# Help-Function: logit
#  use your own function
logit <- function(x) 1 / (1 + exp(-x))

# Help-Function: scale to rate
scale_to_rate <- function(p_raw, target_rate, eps = 0.001) {
  p <- p_raw * target_rate / mean(p_raw)
  pmin(pmax(p, eps), 0.95)
}


# MCAR
apply_mcar_missingness <- function(dat, config, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  out <- dat
  rates <- config$missing_rates_mcar
  
  for (v in names(rates)) {
    n_mis <- floor(nrow(out) * rates[[v]])
    idx_mis <- sample(seq_len(nrow(out)), size = n_mis, replace = FALSE)
    out[[v]][idx_mis] <- NA
  }
  
  attr(out, "missingness") <- "MCAR"
  out
}


# MAR
apply_mar_missingness <- function(dat, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  out <- dat
  N <- nrow(out)
  
  eta <- rnorm(N)
  
  # x3 missing depends on observed x2, x4, x7
  score_x3 <- plogis(0.4 * scale(out$x2)[, 1] -
                       0.3 * scale(out$x4)[, 1] +
                       0.2 * scale(out$x7)[, 1] +
                       eta)
  out$x3[score_x3 > quantile(score_x3, 0.90, na.rm = TRUE)] <- NA
  
  # x5B missing depends on observed x2 and x7
  score_x5B <- plogis(0.5 * scale(out$x2)[, 1] +
                        0.3 * scale(out$x7)[, 1] +
                        rnorm(N))
  out$x5B[score_x5B > quantile(score_x5B, 0.85, na.rm = TRUE)] <- NA
  
  # x5C missing depends on observed x4 and x7
  score_x5C <- plogis(-0.4 * scale(out$x4)[, 1] +
                        0.4 * scale(out$x7)[, 1] +
                        rnorm(N))
  out$x5C[score_x5C > quantile(score_x5C, 0.75, na.rm = TRUE)] <- NA
  
  # x7 missing depends on observed x2, x3, x4
  # Important: use the complete x3 from dat, not out$x3 after missingness
  score_x7 <- plogis(0.3 * scale(dat$x2)[, 1] -
                       0.3 * scale(dat$x3)[, 1] +
                       0.2 * scale(dat$x4)[, 1] +
                       rnorm(N))
  out$x7[score_x7 > quantile(score_x7, 0.80, na.rm = TRUE)] <- NA
  
  attr(out, "missingness") <- "MAR"
  out
}


# MNAR
apply_mnar_missingness <- function(data, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  out <- data
  N <- nrow(out)
  
  # MNAR for relevant covariate x2: 20% missing
  score_x2 <- 0.80 * as.numeric(scale(data$x2)) + rnorm(N)
  cutoff_x2 <- quantile(score_x2, probs = 0.80, na.rm = TRUE)
  out$x2[score_x2 > cutoff_x2] <- NA
  
  # MNAR for redundant covariate x7: 20% missing
  score_x7 <- 0.80 * as.numeric(scale(data$x7)) + rnorm(N)
  cutoff_x7 <- quantile(score_x7, probs = 0.80, na.rm = TRUE)
  out$x7[score_x7 > cutoff_x7] <- NA
  
  attr(out, "missingness") <- "MNAR"
  attr(out, "mnar_variables") <- c("x2", "x7")
  attr(out, "mnar_rates") <- c(x2 = 0.20, x7 = 0.20)
  
  out
}


# apply missingness
apply_missingness <- function(dat,
                              mechanism = c("MCAR", "MAR", "MNAR"),
                              config,
                              seed = NULL) {
  mechanism <- match.arg(mechanism)
  
  if (mechanism == "MCAR") {
    return(apply_mcar_missingness(dat, config = config, seed = seed))
  }
  
  if (mechanism == "MAR") {
    return(apply_mar_missingness(dat, seed = seed))
  }
  
  if (mechanism == "MNAR") {
    return(apply_mnar_missingness(dat, seed = seed))
  }
}


# simulated datasets
make_simulation_data <- function(n_rep,
                                          mechanisms = c("MCAR", "MAR", "MNAR"),
                                          sim_control = sim_config,
                                          config = sim_config) {
  out <- list()
  counter <- 1L
  
  mechanisms <- match.arg(
    mechanisms,
    choices = c("MCAR", "MAR", "MNAR"),
    several.ok = TRUE
  )
  
  seed_base <- config$seed_base
  if (is.null(seed_base) || length(seed_base) != 1L || is.na(seed_base)) {
    seed_base <- 10000L
  }
  seed_base <- as.integer(seed_base)
  
  for (rep_id in seq_len(n_rep)) {
    
    data_full <- simulate_re_probit_data(
      config = config,
      seed = seed_base + rep_id
    )
    
    for (mechanism in mechanisms) {
      
      data_mis <- apply_missingness(
        dat = data_full,
        mechanism = mechanism,
        config = config,
        seed = seed_base + 1000L + 100L * rep_id + match(mechanism, c("MCAR", "MAR", "MNAR"))
      )
      
      out[[counter]] <- list(
        rep_id = rep_id,
        mechanism = mechanism,
        data_full = data_full,
        data_mis = data_mis,
        complete_case_rate = mean(stats::complete.cases(data_mis))
      )
      
      counter <- counter + 1L
    }
  }
  
  out
}


# summary of the missings
missing_summary <- function(dat) {
  data.frame(
    variable = names(dat),
    n_missing = colSums(is.na(dat)),
    missing_rate = round(colMeans(is.na(dat)), 3)
  )
}


# Export simulated datasets as .CSV
export_simulation_data <- function(n_rep,
                                            mechanisms = c("MCAR", "MAR", "MNAR"),
                                            sim_control_arg = NULL,
                                            model_specs_arg = NULL,
                                            out_dir = "data/simulated") {
  if (is.null(sim_control_arg)) {
    sim_control_arg <- sim_config
  }
  
  if (is.null(model_specs_arg)) {
    model_specs_arg <- model_specs
  }
  
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  sim_data <- make_simulation_data(
    n_rep = n_rep,
    mechanisms = mechanisms,
    sim_control = sim_control_arg
  )
  
  index <- data.frame()
  
  for (i in seq_along(sim_data)) {
    rep_id <- sim_data[[i]]$rep_id
    mech <- sim_data[[i]]$mechanism
    
    file_full <- file.path(out_dir, sprintf("rep_%03d_%s_full.csv", rep_id, mech))
    file_mis  <- file.path(out_dir, sprintf("rep_%03d_%s_mis.csv", rep_id, mech))
    
    write.csv(sim_data[[i]]$data_full, file_full, row.names = FALSE, na = "")
    write.csv(sim_data[[i]]$data_mis,  file_mis,  row.names = FALSE, na = "")
    
    index <- rbind(index, data.frame(
      rep_id = rep_id,
      mechanism = mech,
      file_full = normalizePath(file_full, winslash = "/", mustWork = FALSE),
      file_mis = normalizePath(file_mis, winslash = "/", mustWork = FALSE),
      complete_case_rate = sim_data[[i]]$complete_case_rate
    ))
  }
  
  index_file <- file.path(out_dir, "simulation_index.csv")
  write.csv(index, index_file, row.names = FALSE)
  
  model_tab <- do.call(rbind, lapply(names(model_specs_arg), function(nm) {
    f <- model_specs_arg[[nm]]$fixed
    terms <- attr(stats::terms(f), "term.labels")
    
    data.frame(
      model = nm,
      terms = paste(terms, collapse = ","),
      stringsAsFactors = FALSE
    )
  }))
  
  model_file <- file.path(out_dir, "model_specs.csv")
  write.csv(model_tab, model_file, row.names = FALSE)
  
  true_beta_tab <- data.frame(
    term = names(sim_config$beta_true),
    true_value = as.numeric(sim_config$beta_true)
  )
  
  true_beta_file <- file.path(out_dir, "true_beta.csv")
  write.csv(true_beta_tab, true_beta_file, row.names = FALSE)
  
  invisible(list(
    index = index,
    index_file = normalizePath(index_file, winslash = "/", mustWork = TRUE),
    model_file = normalizePath(model_file, winslash = "/", mustWork = TRUE),
    true_beta_file = normalizePath(true_beta_file, winslash = "/", mustWork = TRUE)
  ))
  
  print(index)
}
