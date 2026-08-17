
#' @title quants_to_pars
#' @import nleqslv
#' @description Computes parameters for either a Beta or Gamma distribution given an interval. Uses a normal approximation and method of moments to generate reasonable initial values then solves numerically.
#' @param dist Either "beta" or "gamma" to specify distribution.
#' @param lower Lower bound of interval.
#' @param upper Upper bound of interval.
#' @param pb Vector containing lower and upper probability bounds.
#' @param N_sim Number of simulations used to validate solution.
#' @return List containing parameter solution (shape and rate parameters if Gamma) with sampled mean and interval to validate solution.
#' @rdname quants_to_pars
#' @export
quants_to_pars <- function(dist, lower, upper, pb, N_sim){
  sigma <- (upper - lower)/(qnorm(pb[2]) - qnorm(pb[1]))
  mu <- lower - qnorm(pb[1])*sigma
  if(dist == "beta"){
    init_vals <- c(mu*((mu*(1 - mu))/sigma^2 - 1), (1 - mu)*((mu*(1 - mu))/sigma^2 - 1))
    f <- function(pars) c(pbeta(lower, pars[1], pars[2]) - pb[1], pbeta(upper, pars[1], pars[2]) - pb[2])
    sol <- nleqslv(init_vals, f)$x
    sam <- rbeta(N_sim, sol[1], sol[2])
  } else if(dist == "gamma"){
    init_vals <- c(mu^2/sigma^2, mu/sigma^2)
    f <- function(pars) c(pgamma(lower, pars[1], rate = pars[2]) - pb[1], pgamma(upper, pars[1], rate = pars[2]) - pb[2])
    sol <- nleqslv(init_vals, f)$x
    sam <- rgamma(N_sim, sol[1], rate = sol[2])
  }
  return(list(pars = sol, mean = mean(sam), int = quantile(sam, prob = pb)))
}

#' @title mean_and_quant_to_pars
#' @description Computes parameters for either a Beta or Gamma distribution given mean and single quantile. Uses a normal approximation and method of moments to generate reasonable initial values then solves numerically.
#' @param dist Either "beta" or "gamma" to specify distribution.
#' @param mn Distribution mean.
#' @param quant Distribution quantile.
#' @param pb Probability associated with quantile.
#' @param N_sim Number of simulations used to validate solution.
#' @return List containing parameter solution (shape and rate parameters if Gamma) with sampled mean and interval to validate solution.
#' @rdname mean_and_quant_to_pars
#' @export
mean_and_quant_to_pars <- function(dist, mn, quant, pb, N_sim){
  low_root <- 0.1
  upp_root <- 10
  if(dist == "beta"){
    f <- function(par) pbeta(quant, mn*par, (1 - mn)*par) - pb
    while(f(low_root)*f(upp_root) > 0){
      low_root <- low_root/10
      upp_root <- upp_root*10
      if(upp_root > 1e12) stop("Could not bracket root.")
    }
    phi <- uniroot(f, lower = low_root, upper = upp_root)$root
    sol <- c(mn*phi, (1 - mn)*phi)
    sam <- rbeta(N_sim, sol[1], sol[2])
  } else if(dist == "gamma"){
    f <- function(par) pgamma(quant, mn*par, rate = par) - pb
    while(f(low_root)*f(upp_root) > 0){
      low_root <- low_root/10
      upp_root <- upp_root*10
      if(upp_root > 1e12) stop("Could not bracket root.")
    }
    lambda <- uniroot(f, lower = low_root, upper = upp_root)$root
    sol <- c(mn*lambda, lambda)
    sam <- rgamma(N_sim, sol[1], rate = sol[2])
  } else if(dist == "gamma_shift"){
    sigma <- (quant[2] - quant[1])/(qnorm(pb[2]) - qnorm(pb[1]))
    alpha0 <- (mn/sigma)^2
    beta0 <- mn/sigma^2
    delta0 <- quant[1] - qgamma(pb[1], alpha0, rate = beta0)
    init_vals <- c(log(alpha0), log(beta0), delta0)
    f <- function(pars)
      pars[3] + c(exp(pars[1])/exp(pars[2]) - mn,
                  qgamma(pb[1], exp(pars[1]), rate = exp(pars[2])) - quant[1],
                  qgamma(pb[2], exp(pars[1]), rate = exp(pars[2])) - quant[2])
    sol <- nleqslv(init_vals, f)$x
    sol[1:2] <- exp(sol[1:2])
    sam <- sol[3] + rgamma(N_sim, sol[1], rate = sol[2])
  }
  return(list(pars = sol, mean = mean(sam), int = quantile(sam, prob = pb)))
}

#' @title comp_pars
#' @import nleqslv
#' @description Computes parameters for either a Beta, Gamma or Shifted-Gamma distribution using either a mean and single quantile or an interval. Uses a normal approximation and method of moments to generate reasonable initial values then solves numerically.
#' @param dist Either "beta", "gamma", or "gamma_shift" to specify distribution.
#' @param mn Distribution mean. If specified then only one of lower or upper may also be specified, unless fitting a Shifted-Gamma distribution. Default is NULL.
#' @param lower Lower bound of interval. If specified then only one of mean or upper may also be specified, unless fitting a Shifted-Gamma distribution. Default is NULL.
#' @param upper Upper bound of interval. If specified then only one of mean or lower may also be specified, unless fitting a Shifted-Gamma distribution. Default is NULL.
#' @param pb Probability bound/s. If using two quantiles then pb must be a vector of length two with the upper and lower probability bounds. If using the mean and a single quantile then must be a scalar with the corresponding quantile.
#' @param N_sim Number of simulations used to validate solution. Default is 100000.
#' @return List containing parameter solution (shape and rate parameters if Gamma or Shifted-Gamma) with sampled mean and interval to validate solution.
#' @rdname comp_pars
#' @export
comp_pars <- function(dist, mn = NULL, lower = NULL, upper = NULL, pb, N_sim = 100000){
  if(!(dist %in% c("beta", "gamma", "gamma_shift"))) stop("dist must be either beta, gamma, or gamma_shift.")
  if(dist %in% c("beta", "gamma") &
     (as.integer(!is.null(mn)) +
      as.integer(!is.null(lower)) +
      as.integer(!is.null(upper)) != 2)) stop("Only two out of mn, lower, and upper may be specified if fitting a Beta or Gamma distribution.")
  if(dist == "gamma_shift" &
     (as.integer(!is.null(mn)) +
      as.integer(!is.null(lower)) +
      as.integer(!is.null(upper)) != 3)) stop("Must specify mn, lower, and upper if fitting a Shifted-Gamma distribution.")
  if(dist %in% c("beta", "gamma") & !is.null(mn) & length(pb) > 1)
    stop("pb must be a scalar if fitting a Beta or Gamma distribution and mn is specified.")
  if(is.null(mn)){
    out <- quants_to_pars(dist = dist, lower = lower, upper = upper, pb = pb, N_sim = N_sim)
  } else {
    if(dist %in% c("beta", "gamma")) quant <- ifelse(!is.null(lower), lower, upper)
    if(dist == "gamma_shift") quant <- c(lower, upper)
    out <- mean_and_quant_to_pars(dist = dist, mn = mn, quant = quant, pb = pb, N_sim = N_sim)
  }
  return(out)
}

#' @title beta_mean
#' @description Computes the analytic mean of a Beta distribution.
#' @param x Vector of length two containing the shape parameters.
#' @return Scalar mean.
#' @rdname beta_mean
#' @export
beta_mean <- function(x) x[1]/(x[1] + x[2])

#' @title beta_prod_mean
#' @description Computes the analytic mean of the product of two independent Beta distributions.
#' @param x Vector of length two containing the shape parameters for the first Beta distribution.
#' @param y Vector of length two containing the shape parameters for the second Beta distribution.
#' @return Scalar mean.
#' @rdname beta_prod_mean
#' @export
beta_prod_mean <- function(x, y) beta_mean(x)*beta_mean(y)

#' @title beta_var
#' @description Computes the analytic variance of a Beta distribution.
#' @param x Vector of length two containing the shape parameters.
#' @return Scalar variance.
#' @rdname beta_var
#' @export
beta_var <- function(x) x[1]*x[2]/((x[1] + x[2])^2*(x[1] + x[2] + 1))

#' @title beta_prod_var
#' @description Computes the analytic variance of the product of two independent Beta distributions.
#' @param x Vector of length two containing the shape parameters for the first Beta distribution.
#' @param y Vector of length two containing the shape parameters for the second Beta distribution.
#' @return Scalar variance.
#' @rdname beta_prod_var
#' @export
beta_prod_var <- function(x, y) beta_var(x)*beta_var(y) + beta_var(x)*beta_mean(y)^2 + beta_var(y)*beta_mean(x)^2

#' @title derive_conditional
#' @description Derives the conditional distribution from the marginal distributions of a uni-directional probability tree by moment matching a Beta distribution.
#' @param pars Two x two matrix containing the Beta distribution parameters (columns) for the marginal distributions (rows). Specify either pars OR mu AND sigma2. Default is NULL.
#' @param mu Vector of length two containing the means of the marginal distributions. Specify either pars OR mu AND sigma2. Default is NULL.
#' @param sigma2 Vector of length two containing the variances of the marginal distributions. Specify either pars OR mu AND sigma2. Default is NULL.
#' @return Vector of length two containing the matched Beta distribution parameters for the conditional distribution.
#' @rdname derive_conditional
#' @export
derive_conditional <- function(pars = NULL, mu = NULL, sigma2 = NULL){
  if(( is.null(pars) &  (is.null(mu) | is.null(sigma2))) |
     (!is.null(pars) & !(is.null(mu) & is.null(sigma2)))) stop("Specify either pars OR mu AND sigma2")
  if(!is.null(pars)){
    mu <- apply(pars, 1, beta_mean)
    sigma2 <- apply(pars, 1, beta_var)
  }
  mu_cond <- mu[2]/mu[1]
  sigma2_cond <- (sigma2[2] + mu[2]^2)/(sigma2[1] + mu[1]^2) - mu_cond^2
  a_cond <- mu_cond*(mu_cond*(1 - mu_cond)/sigma2_cond - 1)
  b_cond <- (1 - mu_cond)*(mu_cond*(1 - mu_cond)/sigma2_cond - 1)
  return(c(a = a_cond, b = b_cond))
}
