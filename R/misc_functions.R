
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
  }
  return(list(pars = sol, mean = mean(sam), int = quantile(sam, prob = pb)))
}

#' @title comp_pars
#' @import nleqslv
#' @description Computes parameters for either a Beta or Gamma distribution using either a mean and single quantile or an interval. Uses a normal approximation and method of moments to generate reasonable initial values then solves numerically.
#' @param dist Either "beta" or "gamma" to specify distribution.
#' @param mn Distribution mean. If specified then only one of lower or upper may also be specified. Default is NULL.
#' @param lower Lower bound of interval. If specified then only one of mean or upper may also be specified. Default is NULL.
#' @param upper Upper bound of interval. If specified then only one of mean or lower may also be specified. Default is NULL.
#' @param pb Probability bound/s. If using two quantiles then pb must be a vector of length two with the upper and lower probability bounds. If using the mean and a single quantile then must be a scalar with the corresponding quantile.
#' @param N_sim Number of simulations used to validate solution. Default is 100000.
#' @return List containing parameter solution (shape and rate parameters if Gamma) with sampled mean and interval to validate solution.
#' @rdname comp_pars
#' @export
comp_pars <- function(dist, mn = NULL, lower = NULL, upper = NULL, pb, N_sim = 100000){
  if(!(dist %in% c("beta", "gamma"))) stop("dist must be either beta or gamma.")
  if(as.integer(!is.null(mn)) +
     as.integer(!is.null(lower)) +
     as.integer(!is.null(upper)) != 2) stop("Only two out of mn, lower and upper may be specified.")
  if(!is.null(mn) & length(pb) > 1) stop("pb must be a scalar if mn is specified.")
  if(is.null(mn)){
    out <- quants_to_pars(dist = dist, lower = lower, upper = upper, pb = pb, N_sim = N_sim)
  } else {
    quant <- ifelse(!is.null(lower), lower, upper)
    out <- mean_and_quant_to_pars(dist = dist, mn = mn, quant = quant, pb = pb, N_sim = N_sim)
  }
  return(out)
}

