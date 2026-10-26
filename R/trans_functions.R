
#' @title sub_trans_parms
#' @description Subsets the parameter distributions for the transmission model parameters. Includes strategy effectiveness parameters.
#' @param trans_parms Full set of transmission parameter distributions.
#' @param start First draw to subset.
#' @param end Final draw to subset.
#' @return Subsetted list containing transmission model parameter distributions.
#' @rdname sub_trans_parms
#' @export
sub_trans_parms <- function(trans_parms, start, end){
  trans_parms_sub <- list(total_pop   = trans_parms$total_pop,
                          age_years   = trans_parms$age_years,
                          age_months  = trans_parms$age_months,
                          size_months = trans_parms$size_months,
                          mixing      = trans_parms$mixing,
                          b0          = trans_parms$b0,
                          b1          = trans_parms$b1,
                          phi         = trans_parms$phi,
                          omega       = trans_parms$omega,
                          delta       = trans_parms$delta[start:end],
                          gamma       = trans_parms$gamma[start:end],
                          nu          = trans_parms$nu[start:end],
                          r_sigma     = trans_parms$r_sigma[start:end],
                          dur_V       = trans_parms$dur_V[start:end],
                          kappa_V     = trans_parms$kappa_V[start:end],
                          dur_M       = trans_parms$dur_M[start:end],
                          kappa_M     = trans_parms$kappa_M[start:end],
                          rho_V       = trans_parms$rho_V[start:end],
                          rho_M       = trans_parms$rho_M[start:end])
  return(trans_parms_sub)
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
#' @return Array of estimated model states after burn-in to use as initial values.
#' @rdname burn_trans_model
#' @export
burn_trans_model <- function(N_draw, trans_parms, burn_time = 6000, burn_batch_size = 100, ncores = 1, path = "."){
  N_burn_batch <- ceiling(burn_time/burn_batch_size)
  burn_batch_lens <- sapply(1:N_burn_batch, function(i) ifelse(i < N_burn_batch, burn_batch_size, burn_time - (N_burn_batch - 1)*burn_batch_size))
  y0_burn_fnames <- list.files("Parameters", paste0("y0_burn_", N_draw/1000, "K_"), full.names = TRUE)
  if(length(y0_burn_fnames) == 0){
    st <- 1
  } else {
    st <- max(as.numeric(sapply(y0_burn_fnames, function(txt) strsplit(txt, "_")[[1]][4]))) + 1
  }
  for(i in st:N_burn_batch){
    cat("Starting batch", i, "of", N_burn_batch, "at", as.character(Sys.time()), "\n")
    if(i == 1){
      y0_burn <- RSVModels::initial_values(mod = "base", size_months = trans_parms$size_months, N_sim = N_draw)
    } else {
      y0_burn <- readRDS(paste0(path, "/Parameters/y0_burn_", N_draw/1000, "K_", i - 1, "_of_", N_burn_batch, ".rds"))
    }
    y0_burn <- RSVModels::mod_base(y0 = y0_burn, max_time = burn_batch_lens[i], parms = trans_parms, N_sim = N_draw, batch_size = N_draw/ncores, ncores = ncores)
    y0_burn <- y0_burn[,burn_batch_lens[i],,]
    y0_burn[,,5] <- 0
    saveRDS(y0_burn, paste0(path, "/Parameters/y0_burn_", N_draw/1000, "K_", i, "_of_", N_burn_batch, ".rds"))
    gc()
  }
  saveRDS(y0_burn, paste0(path, "/Parameters/y0_burn_", N_draw/1000, "K.rds"))
  return(y0_burn)
}

#' @title run_trans_models
#' @import RSVModels
#' @description Runs the transmission models for each strategy and extracts the incidence outputs and administration outputs (where applicable).
#' @param N_draw Number of simulated runs.
#' @param y0 Array of initial values. Must be of dimension c(N_draw, 75 (age groups), 5 (SEIR + 1)).
#' @param trans_parms Transmission model parameters.
#' @param years Number of years to run the model.
#' @param subsets Number of subsets to run the models. Default is one. Increase to manage memory workload. Must be a factor of N_draw.
#' @param batch_size Size of each batch to run in parallel. Default is NULL, which is then set to ceiling(N_draw/ncores). Reduce to manage memory workload.
#' @param ncores Number of cores to run in parallel. Default is one.
#' @return Incidence and administration aggregated by year.
#' @rdname run_trans_models
#' @export
run_trans_models <- function(N_draw, y0, trans_parms, years, subsets = 1, batch_size = NULL, ncores = 1){
  if(N_draw%%subsets != 0) stop("N_draw must be divisible by subsets.")
  arr_tmp <- array(NA, dim = c(N_draw, years, length(trans_parms$age_years)),
                   dimnames = list("simulation" = 1:N_draw, "year" = paste("year", 1:years), "age" = trans_parms$age_years))
  out <- list(NS   = list(inc = arr_tmp),
              MV   = list(inc = arr_tmp, inc_strat = arr_tmp, admin = arr_tmp[,,1]),
              mAbs = list(inc = arr_tmp, inc_strat = arr_tmp, admin = arr_tmp[,,1]))
  sub_mat <- matrix(matrix(1:N_draw, ncol = subsets)[c(1,N_draw/subsets),], ncol = subsets)
  if(is.null(batch_size)) batch_size <- ceiling(N_draw/subsets/ncores)
  mod_funs <- list(NS = RSVModels::mod_base, MV = RSVModels::mod_vax, mAbs = RSVModels::mod_mab)
  for(j in 1:subsets){
    cat("Generating incidence parameters for subset", j, "of", subsets, "at", as.character(Sys.time()), "\n")
    start <- sub_mat[1,j]
    end <- sub_mat[2,j]
    y0_j <- y0[start:end,,]
    trans_parms_j <- sub_trans_parms(trans_parms = trans_parms, start = start, end = end)
    for(mod in c("NS", "MV", "mAbs")){
      if(mod == "NS"){
        y0_tmp <- y0_j
      } else {
        y0_tmp <- array(data = 0, dim = c(N_draw/subsets, 75, 7))
        y0_tmp[,,c(1:4, 6)] <- y0_j
        mod_str <- ifelse(mod == "MV", "V", "M")
      }
      mod_out <- mod_funs[[mod]](y0 = y0_tmp, max_time = years*12, parms = trans_parms_j, N_sim = N_draw/subsets,
                                 batch_size = batch_size, ncores = ncores, minimal = TRUE)
      if(mod == "NS"){
        inc <- aperm(apply(mod_out[,,,"Incidence"], c(1,3), function(x) colSums(matrix(x, nrow = 12, ncol = years))), c(2, 1, 3))
        dimnames(inc) <- list("simulation" = dimnames(inc)[[1]], "year" = paste("year", 1:years), "age" = dimnames(inc)[[3]])
      } else {
        inc <- aperm(apply(mod_out[,,,"Incidence"], c(1,3), function(x) colSums(matrix(x, nrow = 12, ncol = years))), c(2, 1, 3))
        dimnames(inc) <- list("simulation" = dimnames(inc)[[1]], "year" = paste("year", 1:years), "age" = dimnames(inc)[[3]])
        inc_strat <- aperm(apply(mod_out[,,,paste0("Incidence_", mod_str)], c(1,3), function(x) colSums(matrix(x, nrow = 12, ncol = years))),
                           c(2, 1, 3))
        dimnames(inc_strat) <- list("simulation" = dimnames(inc_strat)[[1]], "year" = paste("year", 1:years), "age" = dimnames(inc_strat)[[3]])
        admin <- array(rep(trans_parms_j[[paste0("kappa_", mod_str)]]*trans_parms_j$size_months[1]*12, years), dim = c(N_draw/subsets, years),
                       dimnames = list("simulation" = dimnames(inc)[[1]], "year" = paste("year", 1:years)))
        out[[mod]]$inc_strat[start:end,,] <- inc_strat
        out[[mod]]$admin[start:end,] <- admin
      }
      out[[mod]]$inc[start:end,,] <- inc
    }
    gc()
  }
  return(out)
}
