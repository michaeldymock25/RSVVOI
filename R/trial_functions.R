
#' @title gen_hyper_parms
#' @description Generates hyperparameters for the statistical model parameter prior distributions using a summary power prior method.
#' @param outcome Optional outcome of interest out of "MA" for medical attendances or "HP" for hospitalisations. Only required if use_lit == TRUE. Default is NULL.
#' @param use_lit Logical. If set to TRUE, the data from the literature will be used to inform the hyperparameters. Otherwise, user supplied sample sizes (n) and events (y) are required. Default is TRUE.
#' @param n Optional named vector of participant sample sizes. Names must be NS, MV and mAbs. Default is NULL.
#' @param y Optional named vector of events. Names must be NS, MV and mAbs. Default is NULL.
#' @param seed Optional seed for replication. Default is NULL.
#' @param N_sim Number of samples to draw to estimate the relative risk distributions for the moment-matching step. Default is 10,000,000.
#' @param path Path to the location of the Parameters folder. Only required if saving the output. Default is "." and may be used if the function is run inside the package.
#' @param save_output Logical. If set to TRUE, then the generated distributions will be saved into the Parameters folder. Default is FALSE.
#' @return Named matrix containing hyperparameters.
#' @rdname gen_hyper_parms
#' @export
gen_hyper_parms <- function(outcome = NULL, use_lit = TRUE, n = NULL, y = NULL, seed = NULL, N_sim = 10000000, path = ".", save_output = FALSE){
  if(!is.null(seed)) set.seed(seed)
  if(use_lit & (is.null(outcome) || !(outcome %in% c("MA", "HP")))) stop('Outcome must be "MA" or "HP" if use_lit == TRUE.')
  if(!use_lit & (!all(c("NS", "MV", "mAbs") %in% names(n)) | !all(c("NS", "MV", "mAbs") %in% names(y))))
    stop("Vectors n and y must contain values named NS, MV and mAbs.")
  if(!use_lit & any(y > n)) stop("Events in y must not exceed sample sizes in n.")
  if(use_lit){
    if(outcome == "MA"){
      n <- c(NS = 3563 + 1003, MV = 3585, mAbs = 2009)
      y <- c(NS = 132 + 54   , MV = 67  , mAbs = 24)
    } else if(outcome == "HP"){
      n <- c(NS = 3563 + 4019, MV = 3585, mAbs = 4038)
      y <- c(NS = 47 + 68    , MV = 21  , mAbs = 12)
    }
  }
  p_sim <- sapply(c("NS", "MV", "mAbs"), function(d) rbeta(N_sim, 1 + y[d], 1 + n[d] - y[d]))
  l_r_d_sim <- sapply(c("MV", "mAbs"), function(d) log(p_sim[,d]/p_sim[,"NS"]))
  hyper_parms <- list(NS   = c(a  =          unname(y["NS"]), b     = unname(n["NS"] - y["NS"])),
                      MV   = c(mu =   mean(l_r_d_sim[,"MV"]), sigma =      sd(l_r_d_sim[,"MV"])),
                      mAbs = c(mu = mean(l_r_d_sim[,"mAbs"]), sigma =    sd(l_r_d_sim[,"mAbs"])))
  if(save_output) saveRDS(hyper_parms, paste0(path, "/Parameters/hyper_parms.rds"))
  return(hyper_parms)
}

#' @title gen_prior_parms
#' @description Generates the (prior) parameter distributions for the statistical model parameters.
#' @param N_draw Number of parameter draws to characterise distributions.
#' @param hyper_parms Hyperparameters. Must be a named list ("NS", "MV", "mAbs") with named elements ("a", "b") for "NS" and ("mu", "sigma") for "MV" and "mAbs".
#' @param alpha Power prior contribution. Must be a scalar between zero and 1 (inclusive).
#' @param seed Optional seed for replication. Default is NULL.
#' @return Statistical model prior parameter distributions.
#' @rdname gen_prior_parms
#' @export
gen_prior_parms <- function(N_draw, hyper_parms, alpha, seed = NULL){
  if(!is.null(seed)) set.seed(seed)
  if(alpha < 0 | alpha > 1) stop("Power prior contribution alpha must be between zero and one.")
  prior_NS   <- rbeta(N_draw, 1 + alpha*(hyper_parms$NS["a"] - 1), 1 + alpha*(hyper_parms$NS["b"] - 1))
  if(alpha == 0){
    sds <- c(MV = 1, mAbs = 1)
  } else {
    sds <- c(MV   = unname(sqrt(hyper_parms$MV["sigma"]^2/alpha)),
             mAbs = unname(sqrt(hyper_parms$mAbs["sigma"]^2/alpha)))
  }
  prior_MV   <- rnorm(N_draw, hyper_parms$MV["mu"], sds["MV"])
  prior_mAbs <- rnorm(N_draw, hyper_parms$mAbs["mu"], sds["mAbs"])
  prior_parms <- cbind(NS = prior_NS, MV = prior_MV, mAbs = prior_mAbs)
  return(prior_parms)
}

#' @title gen_randomisation
#' @description Generates a randomisation list.
#' @param N Total sample size.
#' @param p_alloc Allocation probabilities for each arm. Must be a vector of values between zero and one, and must sum to one.
#' @param exact Logical. If set to TRUE, then the total sample size will be split across the arms according to the proportions defined by p_alloc (equivalent to a blocked design). If set to FALSE then allocations are generated dynamically (i.e., without blocking) and risks imbalance.
#' @return Vector of arm allocations of length N.
#' @examples
#' N <- 1000
#' p_alloc <- c(0.5, 0.5)
#'
#' rand_exact <- gen_randomisation(N, p_alloc, exact = TRUE)
#' table(rand_exact)
#'
#' rand_inexact <- gen_randomisation(N, p_alloc, exact = FALSE)
#' table(rand_inexact)
#' @rdname gen_randomisation
#' @export
gen_randomisation <- function(N, p_alloc, exact){
  if(sum(p_alloc) != 1 | !all(p_alloc >= 0 & p_alloc <= 1)) stop("p_alloc must be a vector of values between zero and one and must sum to one.")
  if(exact){
    alloc <- floor(N*p_alloc)
    remainder <- N - sum(alloc)
    if(remainder > 0) alloc[seq_len(remainder)] <- alloc[seq_len(remainder)] + 1
    rand <- sample(rep(1:length(p_alloc), times = alloc))
  } else {
    rand <- sample(1:length(p_alloc), size = N, replace = TRUE, prob = p_alloc)
  }
  rand <- factor(rand, labels = c("NS", "MV", "mAbs"))
  return(rand)
}

#' @title gen_outcomes
#' @description Generates Bernoulli outcome data.
#' @param rand Randomisation list generated by gen_randomisation().
#' @param p_d Vector of outcome probabilities for each arm. Values must be between zero and one.
#' @return Vector of outcomes of equal length to rand.
#' @examples
#' N <- 1000
#' p_alloc <- c(0.5, 0.5)
#' rand <- gen_randomisation(N, p_alloc, exact = TRUE)
#' p_d <- c(0.2, 0.3)
#'
#' out <- gen_outcomes(rand, p_d)
#' prop.table(table(rand, out), 1)
#' @rdname gen_outcomes
#' @export
gen_outcomes <- function(rand, p_d){
  if(!all(p_d >= 0 & p_d <= 1)) stop("p_d must be a vector of values between zero and one.")
  rbinom(n = length(rand), size = 1, prob = p_d[rand])
}

#' @title log_marginal_post
#' @description Computes the log marginal joint posterior distribution.
#' @param log_r Proposal vector for the log relative risks (length two).
#' @param n Vector of participant sample sizes.
#' @param y Vector of events.
#' @param hyper_parms Hyperparameters. Must be a named list ("NS", "MV", "mAbs") with named elements ("a", "b") for "NS" and ("mu", "sigma") for "MV" and "mAbs".
#' @param alpha Power prior contribution. Must be a scalar between zero and 1 (inclusive).
#' @param N_grid Number of grid points to sample from to evaluate the integral.
#' @return The log marginal joint posterior evaluated at log_r.
#' @rdname log_marginal_post
#' @export
log_marginal_post <- function(log_r, n, y, hyper_parms, alpha, N_grid){
  p_max <- min(c(1, exp(-log_r)))
  p_grid <- seq(0, p_max, length.out = N_grid + 1)[-1]
  r <- c(1, exp(log_r))
  p <- outer(p_grid, r)
  a <- 1 + alpha*(hyper_parms$NS["a"] - 1)
  b <- 1 + alpha*(hyper_parms$NS["b"] - 1)
  log_vals <- apply(p, 1, function(x){
    if(any(x >= 1)) return(-Inf)
    sum(dbinom(y, n, x, log = TRUE)) + dbeta(x[1], a, b, log = TRUE)
  })
  integral <- sum(exp(log_vals - max(log_vals)))*(p_grid[2] - p_grid[1])
  max(log_vals) + log(integral) +
    dnorm(log_r[1], hyper_parms$MV["mu"], sqrt(hyper_parms$MV["sigma"]^2/alpha), log = TRUE) +
    dnorm(log_r[2], hyper_parms$mAbs["mu"], sqrt(hyper_parms$mAbs["sigma"]^2/alpha), log = TRUE)
}

#' @title upt_posterior
#' @import MASS
#' @description Updates the posterior distribution using a Laplace approximation.
#' @param n Vector of participant sample sizes.
#' @param y Vector of events.
#' @param hyper_parms Hyperparameters. Must be a named list ("NS", "MV", "mAbs") with named elements ("a", "b") for "NS" and ("mu", "sigma") for "MV" and "mAbs".
#' @param alpha Power prior contribution. Must be a scalar between zero and 1 (inclusive).
#' @param init_vals Vector containing initial values for optimisation. Default is c(0, 0).
#' @param N_draw Number of samples to draw. Default is 10,000.
#' @param N_grid Number of grid points to sample from to evaluate the integral. Default is 10,000.
#' @return Named matrix of posterior draws for the relative risk parameters.
#' @rdname upt_posterior
#' @export
upt_posterior <- function(n, y, hyper_parms, alpha, init_vals = c(0, 0), N_draw = 10000, N_grid = 10000){
  opt <- optim(par = init_vals,
               fn = log_marginal_post,
               n = n,
               y = y,
               hyper_parms = hyper_parms,
               alpha = alpha,
               N_grid = N_grid,
               control = list(fnscale = -1),
               hessian = TRUE)
  log_r_hat <- opt$par
  log_r_Sigma <- solve(-opt$hessian)
  r_draws <- exp(mvrnorm(n = N_draw, mu = log_r_hat, Sigma = log_r_Sigma))
  colnames(r_draws) <- c("MV", "mAbs")
  return(r_draws)
  # p_NS_prop <- rbeta(N_draw, hpars["a", "NS"] + y["NS"], hpars["b", "NS"] + n["NS"] - y["NS"])
  # r_d_prop <- sapply(c("MV", "mAbs"), function(d) rbeta(N_draw, hpars["a", d], hpars["b", d]))
  # p_d_prop <- apply(r_d_prop, 2, function(x) pmin(pmax(p_NS_prop*x, 1e-12), 1 - 1e-12))
  # logw <- rowSums(sapply(c("MV", "mAbs"), function(d) lchoose(n[d], y[d]) + y[d]*log(p_d_prop[,d]) + (n[d] - y[d])*log1p(-p_d_prop[,d])))
  # logw <- logw - max(logw)
  # w <- exp(logw)
  # w <- w/sum(w)
  # u <- runif(1, 0, 1/N_draw) + (0:(N_draw - 1))/N_draw
  # idx <- pmin(findInterval(u, cumsum(w)) + 1, length(w))
  # p_NS_post <- p_NS_prop[idx]
  # r_d_post <- r_d_prop[idx,]
  # p_d_post <- apply(r_d_post, 2, function(x) p_NS_post*x)
  # rho_post <- 1 - r_d_post
  # return(list(p_NS = p_NS_post, r_d = r_d_post, p_d = p_d_post, rho = rho_post))
}

#' @title sim_trial
#' @description Simulates a trial by generating a randomisation list, generating the outcomes and updating the Beta-Binomial posterior distribution.
#' @param N Total sample size.
#' @param p_alloc Allocation probabilities for each arm. Must be a vector of values between zero and one, and must sum to one.
#' @param p_k Vector of outcome probabilities for each arm. Values must be between zero and one.
#' @param priors Prior distribution hyperparameters. Must be a matrix with rows for the hyperparameters "a" and "b", and columns for the arms.
#' @param exact Logical. If set to TRUE, then the total sample size will be split across the arms according to the proportions defined by p_alloc (equivalent to a blocked design). If set to FALSE then allocations are generated dynamically (i.e., without blocking) and risks imbalance. Default is TRUE.
#' @return Matrix of updated hyperparameters with the same dimension names as the input priors.
#' @examples
#' N <- 1000
#' p_alloc <- c(0.5, 0.5)
#' p_k <- c(0.2, 0.3)
#' priors <- matrix(rep(1, 4), nrow = 2, ncol = 2, dimnames = list(c("a", "b"), c("NS", "MV")))
#'
#' sim_trial(N, p_alloc, p_k, priors)
#' @rdname sim_trial
#' @export
sim_trial <- function(N, p_alloc, p_k, priors, exact = TRUE){
  rand <- gen_randomisation(N = N, p_alloc = p_alloc, exact = exact)
  out <- gen_outcomes(rand = rand, p_k = p_k)
  pars <- upt_posterior(pars = priors, rand = rand, out = out)
  return(pars)
}

#' @title draw_beta
#' @description Draws samples from the Beta distribution given a set of hyperparameters.
#' @param pars Matrix of hyperparameters with rows for the hyperparameters "a" and "b", and columns for the arms.
#' @param N_draw Number of samples to draw.
#' @return Matrix containing N_draw samples (rows) for each arm (columns).
#' @examples
#' pars <- matrix(c(1, 2, 1, 3), nrow = 2, ncol = 2, dimnames = list(c("a", "b"), c("NS", "MV")))
#' N_draw <- 10000
#'
#' draws <- draw_beta(pars, N_draw)
#' colMeans(draws)
#' @rdname draw_beta
#' @export
draw_beta <- function(pars, N_draw) apply(pars, 2, function(arm_pars) rbeta(N_draw, arm_pars["a"], arm_pars["b"]))

#' @title gen_eff_distributions
#' @description Generates uncertainty distributions for the strategy effectiveness parameters.
#' @param pars Matrix of hyperparameters with rows for the hyperparameters "a" and "b", and columns for the arms.
#' @param N_draw Number of samples to draw. Default is 10,000.
#' @return Matrix containing N_draw samples drawn (rows) for each arm (columns).
#' @examples
#' pars <- matrix(c(2, 1, 1, 2, 1, 5), nrow = 2, ncol = 3, dimnames = list(c("a", "b"), c("NS", "MV", "IM")))
#'
#' draws <- gen_eff_distributions(pars)
#' colMeans(draws)
#' @rdname gen_eff_distributions
#' @export
gen_eff_distributions <- function(pars, N_draw = 10000){
  p_k <- draw_beta(pars = pars, N_draw = N_draw)
  cbind(MV = 1 - p_k[,"MV"]/p_k[,"NS"],
        IM = 1 - p_k[,"IM"]/p_k[,"NS"])
}

