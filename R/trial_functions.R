
#' @title gen_hyper_parms
#' @description Generates hyperparameters for the statistical model parameter prior distributions using a summary power prior method.
#' @param outcome Optional outcome of interest out of "MA" for medical attendances or "HP" for hospitalisations. Only required if use_lit == TRUE. Default is NULL.
#' @param use_lit Logical. If set to TRUE, the data from the literature will be used to inform the hyperparameters. Otherwise, user supplied sample sizes (n) and events (y) are required. Default is TRUE.
#' @param n Optional named vector of participant sample sizes. Names must be NS, MV and mAbs. Default is NULL.
#' @param y Optional named vector of events. Names must be NS, MV and mAbs. Default is NULL.
#' @param OR Logical. If set to TRUE, hyperparameters will be estimated for the odds ratio. Otherwise, hyperparameters will be estimated for the relatve risk. Default is TRUE.
#' @param seed Optional seed for replication. Default is NULL.
#' @param N_sim Number of samples to draw to estimate the relative risk distributions for the moment-matching step. Default is 10,000,000.
#' @param path Path to the location of the Parameters folder. Only required if saving the output. Default is "." and may be used if the function is run inside the package.
#' @param save_output Logical. If set to TRUE, then the generated distributions will be saved into the Parameters folder. Default is FALSE.
#' @return Named matrix containing hyperparameters.
#' @rdname gen_hyper_parms
#' @export
gen_hyper_parms <- function(outcome = NULL, use_lit = TRUE, n = NULL, y = NULL, OR = TRUE, seed = NULL, N_sim = 10000000, path = ".", save_output = FALSE){
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
  if(OR){
    l_par_sim <- sapply(c("MV", "mAbs"), function(d) log((p_sim[,d]*(1 - p_sim[,"NS"]))/((1 - p_sim[,d])*p_sim[,"NS"])))
  } else {
    l_par_sim <- sapply(c("MV", "mAbs"), function(d) log(p_sim[,d]/p_sim[,"NS"]))
  }
  hyper_parms <- list(NS   = c(a  =          unname(y["NS"]), b     = unname(n["NS"] - y["NS"])),
                      MV   = c(mu =   mean(l_par_sim[,"MV"]), sigma =      sd(l_par_sim[,"MV"])),
                      mAbs = c(mu = mean(l_par_sim[,"mAbs"]), sigma =    sd(l_par_sim[,"mAbs"])))
  if(save_output) saveRDS(hyper_parms, paste0(path, "/Parameters/hyper_parms.rds"))
  return(hyper_parms)
}

#' @title gen_prior_parms
#' @description Generates the (prior) parameter distributions for the statistical model parameters.
#' @param N_draw Number of parameter draws to characterise distributions.
#' @param hyper_parms Hyperparameters. Must be a named list ("NS", "MV", "mAbs") with named elements ("a", "b") for "NS" and ("mu", "sigma") for "MV" and "mAbs".
#' @param alpha Power prior contribution. Must be a scalar between zero and 1 (inclusive).
#' @param OR Logical. If set to TRUE, parameters will be transformed as if they were odds ratios. Otherwise, parameters will be transformed as if they were relative risks. Default is TRUE.
#' @param seed Optional seed for replication. Default is NULL.
#' @return Statistical model prior parameter distributions.
#' @rdname gen_prior_parms
#' @export
gen_prior_parms <- function(N_draw, hyper_parms, alpha, OR = TRUE, seed = NULL){
  if(!is.null(seed)) set.seed(seed)
  if(alpha < 0 | alpha > 1) stop("Power prior contribution alpha must be between zero and one.")
  p_NS <- rbeta(N_draw, 1 + alpha*(hyper_parms$NS["a"] - 1), 1 + alpha*(hyper_parms$NS["b"] - 1))
  if(alpha == 0){
    sds <- c(MV = 1, mAbs = 1)
  } else {
    sds <- c(MV   = unname(sqrt(hyper_parms$MV["sigma"]^2/sqrt(alpha))),
             mAbs = unname(sqrt(hyper_parms$mAbs["sigma"]^2/sqrt(alpha))))
  }
  pi_MV   <- rnorm(N_draw, hyper_parms$MV["mu"], sds["MV"])
  pi_mAbs <- rnorm(N_draw, hyper_parms$mAbs["mu"], sds["mAbs"])
  if(OR){
    p_MV <- plogis(qlogis(p_NS) + pi_MV)
    p_mAbs <- plogis(qlogis(p_NS) + pi_mAbs)
  } else {
    p_MV <- p_NS*exp(pi_MV)
    p_mAbs <- p_NS*exp(pi_mAbs)
  }
  r_MV <- p_MV/p_NS
  r_mAbs <- p_mAbs/p_NS
  prior_parms <- cbind(p_NS = p_NS, p_MV = p_MV, p_mAbs = p_mAbs, pi_MV = pi_MV, pi_mAbs = pi_mAbs, r_MV = r_MV, r_mAbs = r_mAbs)
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
  rand <- factor(factor(rand, labels = c("NS", "MV", "mAbs")[sort(unique(rand))]), levels = c("NS", "MV", "mAbs"))
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

#' @title gen_data
#' @description Generates data from a series of trials.
#' @param N_draw Number of parameter draws to characterise distributions.
#' @param N Total trial sample size.
#' @param p_alloc Allocation probabilities for each arm. Must be a vector of values between zero and one, and must sum to one.
#' @param p_d Named list ("MA", "HP") of matrices containing draws of p_d with columns for p_NS, p_MV and p_mAbs.
#' @param exact Logical. If set to TRUE, then the total sample size will be split across the arms according to the proportions defined by p_alloc (equivalent to a blocked design). If set to FALSE then allocations are generated dynamically (i.e., without blocking) and risks imbalance. Default is TRUE.
#' @return Arrays containing simulated data for each trial.
#' @rdname gen_data
#' @export
gen_data <- function(N_draw, N, p_alloc, p_d, exact = TRUE){
  rand <- lapply(1:N_draw, function(i) gen_randomisation(N = N, p_alloc = p_alloc, exact = exact))
  data_arr <- array(NA, dim = c(2, 3, N_draw, 2), dimnames = list(c("n", "y"), c("NS", "MV", "mAbs"), trial = 1:N_draw, c("MA", "HP")))
  for(q in c("MA", "HP")){
    data_tmp <- lapply(1:N_draw, function(i){
      outcomes <- gen_outcomes(rand = rand[[i]], p_d = p_d[[q]][i,])
      n <- tapply(outcomes, rand[[i]], length)
      y <- tapply(outcomes, rand[[i]], sum)
      n[is.na(n)] <- 0
      y[is.na(y)] <- 0
      rbind(n = n, y = y)
    })
    data_arr[,,,q] <- unlist(data_tmp)
  }
  return(data_arr)
}

#' @title upt_posterior
#' @import cmdstanr
#' @description Updates the posterior distribution using Markov Chain Monte-Carlo.
#' @param n Vector of participant sample sizes.
#' @param y Vector of events.
#' @param alpha Power prior contribution. Must be a scalar between zero and 1 (inclusive).
#' @param hyper_parms Hyperparameters. Must be a named list ("NS", "MV", "mAbs") with named elements ("a", "b") for "NS" and ("mu", "sigma") for "MV" and "mAbs".
#' @param OR Logical. If set to TRUE, the odds ratio model will be used. Otherwise, the relative risk model will be used.
#' @param N_draw Number of samples to draw.
#' @param MCMC Logical. If set to TRUE, posterior will be estimated using Markov Chain Monte-Carlo. Otherwise, a Laplace approximation will be used (not currently available). Default is TRUE.
#' @param nchains Number of chains to run in parallel.
#' @return Named matrix of posterior draws.
#' @rdname upt_posterior
#' @export
upt_posterior <- function(n, y, alpha, hyper_parms, OR, N_draw, MCMC = TRUE, nchains = 8){
  if(MCMC){
    mod <- cmdstan_model("Models/RSV_statistical_model.stan")
    stan_dat <- list(n = n,
                     y = y,
                     alpha = alpha,
                     a_NS = hyper_parms$NS["a"],
                     b_NS = hyper_parms$NS["b"],
                     mu = c(hyper_parms$MV["mu"], hyper_parms$mAbs["mu"]),
                     sigma = c(hyper_parms$MV["sigma"], hyper_parms$mAbs["sigma"]),
                     OR = OR)
    fit <- mod$sample(data = stan_dat, chains = nchains, parallel_chains = nchains, iter_sampling = ceiling(N_draw/nchains))
    par_draws <- apply(fit$draws(c("p_NS", "p_d", "pi", "r")), 3, as.vector)
    colnames(par_draws) <- c("p_NS", "p_MV", "p_mAbs", "pi_MV", "pi_mAbs", "r_MV", "r_mAbs")
  } else {
    stop("Only the MCMC implementation is currently available.")
  }
  return(par_draws)
}
