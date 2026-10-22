source("R/00_config.R")
load_required_packages()
source("R/02_cart_bb_update.R")
source("R/03_bayes_sampler.R")

set.seed(7301)
n <- 120L
J <- 8L
id <- rep(seq_len(J), each = n / J)
x2 <- rnorm(n)
x3 <- rbinom(n, 1, 0.45)
x4 <- rnorm(n)
a <- rnorm(J, 0, 0.7)
eta <- 0.2 - 0.9 * x2 + 0.5 * x3 + a[id]
y <- as.integer(rnorm(n) < eta)
dat <- data.frame(y = y, id = id, x2 = x2, x3 = x3, x4 = x4)
dat$x2[sample.int(n, 24L)] <- NA_real_
dat$x3[sample.int(n, 18L)] <- NA_integer_
dat$x4[sample.int(n, 20L)] <- NA_real_

control <- list(
  minbucket_range = 5:8,
  maxdepth = 3L,
  cp = 0.001,
  max_donors = 20L
)

reference_fit <- fit_bayes_re_probit_cart(
  data = dat,
  fixed_formula = y ~ x2 + x3,
  n_iter = 160L,
  burnin = 40L,
  cart_control = control,
  seed = 91L,
  verbose_every = 0L
)
medoid <- select_medoid_completed_data(reference_fit, dat, medoid_thin = 4L)
stopifnot(!anyNA(medoid$data[c("x2", "x3", "x4")]))

fixed_fit <- fit_bayes_re_probit_cart(
  data = medoid$data,
  fixed_formula = y ~ x2 + x3,
  n_iter = 220L,
  burnin = 60L,
  cart_control = control,
  seed = 92L,
  verbose_every = 0L
)
theta_star <- get_chib_fixed_points(fixed_fit)$posterior_mean
evidence <- compute_chib_logml(
  fixed_fit,
  medoid$data,
  theta_star,
  reduced_iter = 120L,
  reduced_burnin = 30L,
  reduced_thin = 2L,
  details = TRUE
)

stopifnot(is.finite(evidence$log_marginal_likelihood))
print(evidence)
