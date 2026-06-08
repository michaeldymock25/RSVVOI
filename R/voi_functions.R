
#' @title evpi
#' @description Computes the expected value of perfect information.
#' @param INMB Named list containing incremental net monetary benefit for each willingness-to-pay threshold.
#' @return Named vector containing the expected value of perfect information for each willingness-to-pay threshold.
#' @rdname evpi
#' @export
evpi <- function(INMB) sapply(INMB, function(x) mean(apply(x, 1, max)) - max(colMeans(x)))

#' @title evppi
#' @import mgcv
#' @import RhpcBLASctl
#' @import stats
#' @description Computes the expected value of partial perfect information using a non-parametric regression approximation method.
#' @param INMB Named list containing incremental net monetary benefit for each willingness-to-pay threshold.
#' @param psa Named matrix of parameter draws for parameters of interest.
#' @param model Generalised additive regression model specification (formula). Should be a function of the parameters of interest.
#' @return Named list containing named vector with the expected value of partial perfect information for each willingness-to-pay threshold in addition to a named list of the partial incremental net benefit.
#' @rdname evppi
#' @export
evppi <- function(INMB, psa, model){
  RhpcBLASctl::blas_set_num_threads(1)
  N_draw <- unique(sapply(INMB, nrow))
  D <- unique(sapply(INMB, ncol))
  INMB_partial <- lapply(1:length(INMB), function(i) matrix(data = NA, nrow = N_draw, ncol = D))
  EVPPI <- vector(length = length(INMB))
  for(i in 1:length(INMB)){
    INMB_partial[[i]][,1] <- 0
    for(d in 2:D) INMB_partial[[i]][,d] <- gam(update(formula(INMB[[i]][,d] ~ .), formula(paste(".~", model))), data = as.data.frame(psa))$fitted
    EVPPI[i] <- mean(apply(INMB_partial[[i]], 1, max)) - max(colMeans(INMB_partial[[i]]))
  }
  names(INMB_partial) <- names(INMB)
  names(EVPPI) <- names(INMB)
  return(list(EVPPI = EVPPI, INMB_partial = INMB_partial))
}

#' @title prep_dat_evsi
#' @description Prepares data for input into the evsi function.
#' @import EVSI
#' @param p_d Named list ("MA", "HP") of matrices containing draws of p_d with columns for p_NS, p_MV and p_mAbs.
#' @param N Total trial sample size.
#' @param p_alloc Allocation probabilities for each arm. Must be a vector of values between zero and one, and must sum to one.
#' @param MM Logical. If set to TRUE, then data will be simulated in preparation for the moment matching approximation method to be used in the EVSI estimation. Default is FALSE.
#' @param Q Number of model reruns to estimate the expected variance of the posterior net benefit. Only required if MM == TRUE. Default is NULL.
#' @param exact Logical. If set to TRUE, then the total sample size will be split across the arms according to the proportions defined by p_alloc (equivalent to a blocked design). If set to FALSE then allocations are generated dynamically (i.e., without blocking) and risks imbalance. Default is TRUE.
#' @return Expected value of sample information
#' @rdname prep_dat_evsi
#' @export
##
prep_dat_evsi <- function(p_d, N, p_alloc, MM = FALSE, Q = NULL, exact = TRUE){
  if(MM){
    p_d_tmp <- do.call(cbind, p_d)
    colnames(p_d_tmp) <- paste0(colnames(p_d_tmp), rep(c("_MA", "_HP"), each = 3))
    p_d_tmp <- gen.quantiles(parameter = colnames(p_d_tmp), param.mat = p_d_tmp, Q = Q)
    p_d <- list("MA" = as.matrix(p_d_tmp[,c("p_NS_MA", "p_MV_MA", "p_mAbs_MA")]), "HP" = as.matrix(p_d_tmp[,c("p_NS_HP", "p_MV_HP", "p_mAbs_HP")]))
  }
  dat <- gen_data(N_draw = unique(sapply(p_d, nrow)), N = N, p_alloc = p_alloc, p_d = p_d, exact = exact)
  y_tmp <- t(apply(dat["y",,,], 2, as.vector))
  n_tmp <- t(apply(dat["n",,,], 2, as.vector))
  p_tmp <- y_tmp/n_tmp
  dat <- cbind(y_tmp, n_tmp, p_tmp)
  colnames(dat) <- as.vector(sapply(c("y", "n", "p"), function(i) sapply(c("MA", "HP"), function(q) paste(i, c("NS", "MV", "mAbs"), q, sep = "_"))))
  return(dat)
}

#' @title evsi
#' @description Computes the expected value of sample information using a non-parametric regression or moment matching approximation method.
#' @import mgcv
#' @import RhpcBLASctl
#' @import stats
#' @param INMB Named list containing incremental net monetary benefit for each willingness-to-pay threshold.
#' @param dat Named matrix of simulated data summaries.
#' @param method Approximation method. Either NP for non-parametric regression or MM for moment matching. The moment matching method requires the evppi function to be run in advance using the non-parametric regression method to generate INMB_partial. Default is NP.
#' @param model Generalised additive regression model specification (formula). Only required for the non-parametric approximation method. Should be a function of variables in dat. Default is NULL.
#' @param alpha Power prior contribution. Must be a scalar between zero and 1 (inclusive). Only required for the moment matching method. Default is NULL.
#' @param hyper_parms Hyperparameters. Must be a named list ("NS", "MV", "mAbs") with named elements ("a", "b") for "NS" and ("mu", "sigma") for "MV" and "mAbs". Only required for the moment matching method. Default is NULL.
#' @param OR Logical. If set to TRUE, the odds ratio model will be used. Otherwise, the relative risk model will be used. Only required for the moment matching method. Default is TRUE.
#' @param MCMC Logical. If set to TRUE, posterior will be estimated using Markov Chain Monte-Carlo. Otherwise, a Laplace approximation will be used (not currently available).  Only required for the moment matching method. Default is TRUE.
#' @param nchains Number of chains to run in parallel. Only required for the moment matching method. Default is eight.
#' @param trans_parms Transmission model parameters. Only required for the moment matching method. Default is NULL.
#' @param y0_burn Array of initial values for the transmission models. Must be of dimension c(N_draw, 75 (age groups), 5 (SEIR + 1)). Only required for the moment matching method. Default is NULL.
#' @param years Number of years to run the transmission models. Only required for the moment matching method. Default is NULL.
#' @param cea_parms Health economic parameter distributions. Only required for the moment matching method. Default is NULL.
#' @param WTP Vector containing willingness-to-pay thresholds. Only required for the moment matching method. Default is NULL.
#' @param ref Reference group to use for the increments when computing the incremental net monetary benefit. Only required for the moment matching method. Default is "NS".
#' @param INMB_partial Samples of INMB for the parameters of interest generated from the evppi function using the non-parametric regression approximation method. Only required for the moment matching approximation method. Default is NULL.
#' @return Expected value of sample information
#' @rdname evsi
#' @export
evsi <- function(INMB, dat, method = "NP", model = NULL, alpha = NULL, hyper_parms = NULL, OR = TRUE, MCMC = TRUE, nchains = 8,
                 trans_parms = NULL, y0_burn = NULL, years = NULL, cea_parms = NULL, WTP = NULL, ref = "NS", INMB_partial = NULL){
  N_draw <- unique(sapply(INMB, nrow))
  D <- unique(sapply(INMB, ncol))
  if(method == "NP"){
    EVSI <- vector(length = length(INMB))
    RhpcBLASctl::blas_set_num_threads(1)
    for(i in 1:length(INMB)){
      g_hat <- matrix(data = NA, nrow = N_draw, ncol = D)
      g_hat[,1] <- 0
      for(d in 2:D) g_hat[,d] <- gam(update(formula(INMB[[i]][,d] ~ .), formula(paste(".~", model))), data = as.data.frame(dat))$fitted
      EVSI[i] <- mean(apply(g_hat, 1, max)) - max(colMeans(g_hat))
    }
  } else if(method == "MM"){
    var_est <- list()
    trans_parms_tmp <- trans_parms
    cea_parms_tmp <- cea_parms
    for(i in 1:nrow(dat)){
      postr_MA <- upt_posterior(n = dat[i, c("n_NS_MA", "n_MV_MA", "n_mAbs_MA")], y = dat[i, c("y_NS_MA", "y_MV_MA", "y_mAbs_MA")],
                                alpha = alpha, hyper_parms = hyper_parms[["MA"]], OR = OR, N_draw = N_draw, MCMC = MCMC, nchains = nchains)
      postr_HP <- upt_posterior(n = dat[i, c("n_NS_HP", "n_MV_HP", "n_mAbs_HP")], y = dat[i, c("y_NS_HP", "y_MV_HP", "y_mAbs_HP")],
                                alpha = alpha, hyper_parms = hyper_parms[["HP"]], OR = OR, N_draw = N_draw, MCMC = MCMC, nchains = nchains)
      trans_parms_tmp[c("rho_V", "rho_M")] <-  list("rho_V" = 1 - postr_MA[,"r_MV"], "rho_M" = 1 - postr_MA[,"r_mAbs"])
      inc_parms <- run_trans_models(N_draw = N_draw, y0 = y0_burn, trans_parms = trans_parms_tmp, years = years)
      cea_parms_tmp[c("p_HP_I_MV", "p_HP_I_mAbs")] <- list("p_HP_I_MV" = cea_parms$p_HP_I_NS*postr_HP[,"r_MV"], "p_HP_I_mAbs" = cea_parms$p_HP_I_NS*postr_HP[,"r_mAbs"])
      cea_out <- eval_cea(N_draw = N_draw, inc_parms = inc_parms, cea_parms = cea_parms_tmp, ages = trans_parms$age_years)
      INMB_tmp <- eval_INMB(cea_out, WTP = WTP, ref = ref)
      var_est[[i]] <- lapply(INMB_tmp, function(x) var(x[,-1]))
    }
    EVSI <- lapply(1:length(INMB), function(i){
               prior_var <- var(INMB[[i]][,-1])
               mu_mn <- colMeans(INMB[[i]][,-1])
               mu_var <- prior_var - Reduce("+", lapply(var_est, function(x) x[[i]]))/Q
               mu_var_eigen <- eigen(mu_var)
               mu_var_sqrt <- mu_var_eigen$vectors%*%diag(sqrt(mu_var_eigen$values))%*%t(mu_var_eigen$vectors)
               prior_var_eigen <- eigen(prior_var)
               prior_var_sqrt_inv <- chol2inv(chol(prior_var_eigen$vectors%*%diag(sqrt(prior_var_eigen$values))%*%t(prior_var_eigen$vectors)))
               INB_rescaled <- t(t(t(t(INMB_partial[[i]][,-1]) - mu_mn)%*%prior_var_sqrt_inv%*%mu_var_sqrt) + mu_mn)
               mean(apply(INB_rescaled, 1, function(x) max(0, x))) - max(colMeans(INMB[[i]]))
    })
  }
  names(EVSI) <- names(INMB)
  return(EVSI)
}
