
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
#' @param burn_time Number of months to burn-in model. Default is 6000 (equivalent to 500 years).
#' @param burn_batch_size Size of each burn_time batch to run in parallel. Default is 100. Reduce to manage memory workload.
#' @param ncores Number of cores to run in parallel. Default is one.
#' @param path Path to the location of the Data and Parameters folders. Default is "." and may be used if the function is run inside the package.
#' @param save_output Logical. If set to TRUE, then the generated distributions will be saved into the Parameters folder. Default is FALSE.
#' @return Array of estimated model states after burn-in to use as initial values.
#' @rdname burn_trans_model
#' @export
burn_trans_model <- function(N_draw, trans_parms, burn_time = 6000, burn_batch_size = 100, ncores = 1, path = ".", save_output = FALSE){
  N_burn_batch <- ceiling(burn_time/burn_batch_size)
  burn_batch_lens <- sapply(1:N_burn_batch, function(i) ifelse(i < N_burn_batch, burn_batch_size, burn_time - (N_burn_batch - 1)*burn_batch_size))
  y0_burn <- RSVModels::initial_values(mod = "base", size_months = trans_parms$size_months, N_sim = N_draw)
  for(i in 1:N_burn_batch){
    y0_burn <- RSVModels::mod_base(y0 = y0_burn, max_time = burn_batch_lens[i], parms = trans_parms, N_sim = N_draw, batch_size = N_draw/ncores, ncores = ncores)
    y0_burn <- y0_burn[,burn_batch_lens[i],,]
    gc()
  }
  y0_burn[,,5] <- 0
  if(save_output) saveRDS(y0_burn, paste0(path, "/Parameters/y0_burn.rds"))
  return(y0_burn)
}

#' @title run_trans_models
#' @import RSVModels
#' @description Runs the transmission models for each strategy and extracts the incidence outputs and administration outputs (where applicable).
#' @param N_draw Number of simulated runs.
#' @param y0 Array of initial values. Must be of dimension c(N_draw, 75 (age groups), 5 (SEIR + 1)).
#' @param trans_parms Transmission model parameters.
#' @param years Number of years to run the model.
#' @param ncores Number of cores to run in parallel. Default is one.
#' @return Incidence and administration aggregated by year.
#' @rdname run_trans_models
#' @export
run_trans_models <- function(N_draw, y0, trans_parms, years, ncores = 1){
  y0 <- y0[1:N_draw,,]
  mod_funs <- list(base = RSVModels::mod_base, vax = RSVModels::mod_vax, mab = RSVModels::mod_mab)
  out <- lapply(c("base", "vax", "mab"), function(mod){
    if(mod == "base"){
      y0_tmp <- y0
    } else {
      y0_old <- y0
      y0_tmp <- array(data = 0, dim = c(N_draw, 75, 6))
      y0_tmp[,,c(1:4, 6)] <- y0_old
    }
    mod_out <- mod_funs[[mod]](y0 = y0_tmp, max_time = years*12, parms = trans_parms, N_sim = N_draw, batch_size = ceiling(N_draw/ncores), ncores = ncores, minimal = TRUE)
    inc <- aperm(apply(mod_out[,,,"Incidence"], c(1,3), function(x) colSums(matrix(x, nrow = 12, ncol = years))), c(2, 1, 3))
    dimnames(inc) <- list("simulation" = dimnames(inc)[[1]], "year" = paste("year", 1:years), "age" = dimnames(inc)[[3]])
    if(mod == "base"){
      return(list(inc = inc))
    } else {
      admin <- aperm(apply(mod_out[,,,1], c(1,3), function(x) colSums(matrix(x, nrow = 12, ncol = years))), c(2, 1, 3))
      dimnames(admin) <- list("simulation" = dimnames(admin)[[1]], "year" = paste("year", 1:years), "age" = dimnames(admin)[[3]])
      return(list(inc = inc, admin = admin))
    }
  })
  names(out) <- c("NS", "MV", "mAbs")
  return(out)
}
