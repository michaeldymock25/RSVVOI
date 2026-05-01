
#' @title gen_trans_parms
#' @description Generates the parameter distributions for the transmission model parameters. Excludes strategy effectiveness parameters. Requires access to the Data folder containing data for the Australian population (size by age group) and the estimated mixing matrix (computed via the RSVModels package).
#' @param N_draw Number of parameter draws to characterise distributions.
#' @param seed Optional seed for replication. Default is NULL.
#' @param path Path to the location of the Data and Parameters folders. Default is "." and may be used if the function is run inside the package.
#' @param save_output Logical. If set to TRUE, then the generated distributions will be saved into the Parameters folder. Default is FALSE.
#' @return List containing transmission model parameter distributions.
#' @rdname gen_trans_parms
#' @export
gen_trans_parms <- function(N_draw, seed = NULL, path = ".", save_output = FALSE){
  if(!is.null(seed)) set.seed(seed)
  pop_agg <- read.csv(paste0(path, "/Data/australian_population_aggregated.csv"))
  mixing <- read.csv(paste0(path, "/Data/mixing.csv"))
  trans_parms <- list(total_pop   = sum(pop_agg$population),
                      age_years   = c(seq(0, 5, 1/12), seq(10, 75, 5)),
                      age_months  = 12*c(seq(0, 5, 1/12), seq(10, 75, 5)),
                      size_months = c(rep(pop_agg$population[1]/60, 60), pop_agg$population[2:15], sum(pop_agg$population[16:17])),
                      mixing      = as.matrix(mixing),
                      b0          = 0.087,
                      b1          = -0.193,
                      phi         = 1.536,
                      omega       = c(rep(1.00, 59), rep(0.35, 16)),
                      delta       = 1/(rgamma(N_draw, 16, 4)/(365/12)),
                      gamma       = 1/(rgamma(N_draw, 81, 9)/(365/12)),
                      nu          = 1/(rgamma(N_draw, 529, 2.3)/(365/12)),
                      r_sigma     = runif(N_draw, 0.5, 1),
                      dur_V       = sample(3:9, size = N_draw, replace = TRUE),
                      kappa_V     = rbeta(N_draw, 10, 4.286),
                      dur_M       = sample(3:9, size = N_draw, replace = TRUE),
                      kappa_M     = rbeta(N_draw, 10, 1.765))
  if(save_output) saveRDS(trans_parms, paste0(path, "/Parameters/trans_parms.rds"))
  return(trans_parms)
}

#' @title burn_trans_model
#' @import RSVModels
#' @description Runs the base transmission model for a set burn-in period.
#' @param N_draw Number of simulated runs.
#' @param trans_parms Transmission model parameters.
#' @param seed Optional seed for replication. Default is NULL.
#' @param burn_time Number of months to burn-in model. Default is 6000 (equivalent to 500 years).
#' @param batch_size Size of each batch for parallel simulations. Default is 100.
#' @param ncores Number of cores to run in parallel. Default is one.
#' @param path Path to the location of the Data and Parameters folders. Default is "." and may be used if the function is run inside the package.
#' @param save_output Logical. If set to TRUE, then the generated distributions will be saved into the Parameters folder. Default is FALSE.
#' @return Array of estimated model states after burn-in to use as initial values.
#' @rdname burn_trans_model
#' @export
burn_trans_model <- function(N_draw, trans_parms, seed = NULL, burn_time = 6000, batch_size = 100, ncores = 1, path = ".", save_output = FALSE){
  if(!is.null(seed)) set.seed(seed)
  y0 <- initial_values(mod = "base", size_months = trans_parms$size_months, N_sim = N_draw)
  y0_burn <- mod_base(y0 = y0, max_time = burn_time, parms = trans_parms, N_sim = N_draw, batch_size = batch_size, ncores = ncores)
  y0_burn <- y0_burn[,burn_time,,]
  y0_burn[,,5] <- 0
  if(save_output) saveRDS(y0_burn, paste0(path, "/Parameters/y0_burn.rds"))
  return(y0_burn)
}

#' @title run_trans_model
#' @import RSVModels
#' @description Runs a transmission model and extracts the incidence outputs.
#' @param N_draw Number of simulated runs.
#' @param y0 Array of initial values. Must be of dimension c(N_draw, 75 (age groups), 4 SEIR)).
#' @param mod Transmission model to run. Either "base", "vax" or "mab".
#' @param trans_parms Transmission model parameters.
#' @param years Number of years to run the model.
#' @param batch_size Size of each batch for parallel simulations. Default is 100.
#' @param ncores Number of cores to run in parallel. Default is one.
#' @return Incidence aggregated by year.
#' @rdname run_trans_model
#' @export
run_trans_model <- function(N_draw, y0, mod, trans_parms, years, batch_size = 100, ncores = 1){
  inc_dim <- ifelse(mod == "base", 5, 6)
  if(mod %in% c("vax", "mab")){
    y0_old <- y0
    y0 <- array(data = 0, dim = c(N_draw, 75, 6))
    y0[,,c(1:4, 6)] <- y0_old
  }
  mod_f <- match.fun(paste0("mod_", mod))
  out <- mod_f(y0 = y0, max_time = years*12, parms = trans_parms, N_sim = N_draw, batch_size = batch_size, ncores = ncores)
  inc_over_ages <- apply(out[,,,inc_dim], c(1,2), sum)
  inc <- apply(inc_over_ages, 1, function(x) colSums(matrix(x, nrow = 12, ncol = years)))
  return(inc)
}
