
rm(list = ls())

## read in the libraries and source the functions

library(data.table)

source("R/psa_functions.R")
source("R/trans_functions.R")
source("R/cea_functions.R")
source("R/voi_functions.R")
source("R/trial_functions.R")

## set N_draw, alphas, the subset of years, the seed for replication, number of subsets for memory management,
## number of cores and willingness-to-pay thresholds
## note that simulations for N_draw = 10000 are provided. If a different N_draw is desired, then the model will need to be burned in again.

N_draw <- 10000
alphas <- c(0.005, 0.01, 0.1, 0.25, 0.5, 1)
years_sub <- c(1:5, 10, 15)
seed <- 35234
subsets <- 20
ncores <- 15
WTP <- c(20000, 50000)

## generate and/or read in the burned transmission model states and the hyperparameters

if(!file.exists("Data/mixing.csv")) source("Run/make_mixing_matrix.R")
if(!file.exists("Parameters/y0_burn_10K.rds")) source("Run/burn_trans_model.R")
if(!file.exists("Parameters/hyper_parms.rds")) source("Run/make_hyper_parms.R")
y0_burn <- readRDS("Parameters/y0_burn_10K.rds")
hyper_parms <- readRDS("Parameters/hyper_parms.rds")

## generate the strategy effectiveness parameters, the transmission model parameters, and the health economic model parameters

strat_eff_parms_l <- list()
trans_parms_l <- list()
cea_parms_l <- list()
trans_parms <- gen_trans_parms(N_draw = N_draw, seed = seed)
cea_parms <- gen_cea_parms(N_draw = N_draw, seed = seed)
for(i in 1:length(alphas)){
  strat_eff_parms_MA <- gen_strat_eff_parms(N_draw = N_draw, hyper_parms = hyper_parms[["MA"]], alpha = alphas[i], seed = seed)
  strat_eff_parms_HP <- gen_strat_eff_parms(N_draw = N_draw, hyper_parms = hyper_parms[["HP"]], alpha = alphas[i], seed = seed)
  r_mAbs_tmp <- r_MV_tmp <- matrix(1, nrow = N_draw, ncol = length(trans_parms$age_years))
  for(sim in 1:N_draw){
    r_MV_tmp[sim, 1:trans_parms$dur_V[sim]] <- strat_eff_parms_HP[sim,"r_MV"]
    r_mAbs_tmp[sim, 1:trans_parms$dur_M[sim]] <- strat_eff_parms_HP[sim,"r_mAbs"]
  }
  strat_eff_parms_l[[i]] <- list("MA" = strat_eff_parms_MA, "HP" = strat_eff_parms_HP)
  trans_parms_l[[i]] <- c(trans_parms, list("rho_V" = 1 - strat_eff_parms_MA[,"r_MV"], "rho_M" = 1 - strat_eff_parms_MA[,"r_mAbs"]))
  cea_parms_l[[i]] <- c(cea_parms, list("p_HP_ED_MV" = cea_parms$p_HP_ED_NS*r_MV_tmp, "p_HP_ED_mAbs" = cea_parms$p_HP_ED_NS*r_mAbs_tmp))
}
saveRDS(strat_eff_parms_l, "Parameters/strat_eff_parms_l_10K.rds")
saveRDS(trans_parms_l, "Parameters/trans_parms_l_10K.rds")
saveRDS(cea_parms_l, "Parameters/cea_parms_l_10K.rds")

## setup the design configurations

cfg_par <- CJ(alpha = alphas,
              sorted = FALSE)

cfg_voi <- CJ(alpha = alphas,
              design = c("A", "B", "C", "D"),
              N = c(1000, 5000, 10000, 20000, 50000),
              sorted = FALSE)

## run the simulations for each design configuration and save the output

cat("Starting parameter generation simulations at", as.character(Sys.time()), "\n")
for(i in 1:nrow(cfg_par)){
  cat("Starting configuration", i, "of", nrow(cfg_par), "at", as.character(Sys.time()), "\n")
  fname <- paste0("sims_par_cfg_", formatC(i, width = 2, flag = "0"), ".rds")
  if(fname %in% list.files("Simulations/10K")){
    cat("Skipping configuration as parameter generation simulations already exist.\n")
    next
  }
  inc_parms <- run_trans_models(N_draw = N_draw, y0 = y0_burn, trans_parms = trans_parms_l[[i]],
                                years = max(years_sub), subsets = subsets, ncores = ncores)
  saveRDS(list(cfg = cfg_par[i], inc_parms = inc_parms), paste0("Simulations/10K/", fname))
  rm(inc_parms)
  gc()
}

par_fnames <- paste0("Simulations/10K/sims_par_cfg_", formatC(1:nrow(cfg_par), width = 2, flag = "0"), ".rds")
sims_par <- vector("list", nrow(cfg_par))
for(i in 1:nrow(cfg_par)) if(file.exists(par_fnames[i])) sims_par[[i]] <- readRDS(par_fnames[i])

cat("Starting cost-effectiveness simulations at", as.character(Sys.time()), "\n")
for(i in 1:nrow(cfg_par)){
  cat("Starting configuration", i, "of", nrow(cfg_par), "at", as.character(Sys.time()), "\n")
  fname <- paste0("sims_cea_cfg_", formatC(i, width = 2, flag = "0"), ".rds")
  if(fname %in% list.files("Simulations/10K")){
    cat("Skipping configuration as cost-effectiveness simulations already exist.\n")
    next
  }
  cea <- eval_cea(N_draw = N_draw, inc_parms = sims_par[[i]]$inc_parms, cea_parms = cea_parms_l[[i]],
                  ages = trans_parms$age_years, minimal = TRUE)
  INMB <- lapply(years_sub, function(yy) eval_INMB(cea = cea, WTP = WTP, years = yy, max_age_ind = max(which(trans_parms$age_years < 5))))
  names(INMB) <- paste("years =", years_sub)
  saveRDS(list(cfg = cfg_par[i], cea = cea, INMB = INMB), paste0("Simulations/10K/", fname))
  rm(cea, INMB)
  gc()
}

cea_fnames <- paste0("Simulations/10K/sims_cea_cfg_", formatC(1:nrow(cfg_par), width = 2, flag = "0"), ".rds")
sims_cea <- vector("list", nrow(cfg_par))
for(i in 1:nrow(cfg_par)) if(file.exists(cea_fnames[i])) sims_cea[[i]] <- readRDS(cea_fnames[i])

cat("Starting EVPI/EVPPI value of information simulations at", as.character(Sys.time()), "\n")
for(i in 1:nrow(cfg_par)){
  cat("Starting configuration", i, "of", nrow(cfg_par), "at", as.character(Sys.time()), "\n")
  fname <- paste0("sims_voi_p_cfg_", formatC(i, width = 2, flag = "0"), ".rds")
  if(fname %in% list.files("Simulations/10K")){
    cat("Skipping configuration as EVPI/EVPPI value of information simulations already exist.\n")
    next
  }
  INMB <- sims_cea[[i]]$INMB
  cat("Estimating EVPI at", as.character(Sys.time()), "\n")
  EVPI <- lapply(INMB, function(inmb) evpi(INMB = inmb))
  names(EVPI) <- names(INMB)
  cat("Estimating EVPPI at", as.character(Sys.time()), "\n")
  psa <- cbind(rho_MV_MA     =                trans_parms_l[[i]]$rho_V, rho_mAbs_MA    =                  trans_parms_l[[i]]$rho_M,
               r_MV_HP       = strat_eff_parms_l[[i]][["HP"]][,"r_MV"], r_mAbs_HP      = strat_eff_parms_l[[i]][["HP"]][,"r_mAbs"],
               delta         =                       trans_parms$delta, gamma          =                         trans_parms$gamma,
               nu            =                          trans_parms$nu, r_sigma        =                       trans_parms$r_sigma,
               dur_V         =                       trans_parms$dur_V, dur_M          =                         trans_parms$dur_M,
               kappa_V       =                     trans_parms$kappa_V, kappa_M        =                       trans_parms$kappa_M,
               MV_dose_cost  =                  cea_parms$MV_dose_cost, mAbs_dose_cost =                  cea_parms$mAbs_dose_cost,
               p_waste       =                       cea_parms$p_waste, admin_cost     =                      cea_parms$admin_cost,
               GP_cost       =                       cea_parms$GP_cost, N_GP           =                            cea_parms$N_GP,
               Virology_cost =                 cea_parms$Virology_cost, NA_ED_cost     =                      cea_parms$NA_ED_cost,
               ED_cost       =                       cea_parms$ED_cost, HP_cost        =                         cea_parms$HP_cost,
               ICU_cost      =                      cea_parms$ICU_cost,
               GP_NA_ED_dis  =                  cea_parms$GP_NA_ED_dis, HP_dis         =                          cea_parms$HP_dis,
               ICU_dis       =                       cea_parms$ICU_dis, GP_NA_ED_dur   =                    cea_parms$GP_NA_ED_dur,
               HP_dur        =                        cea_parms$HP_dur, ICU_dur        =                         cea_parms$ICU_dur)
  models <- c(strategy_effectiveness = "s(rho_MV_MA) + s(rho_mAbs_MA) + s(r_MV_HP) + s(r_mAbs_HP)",
              transmission           = "s(delta) + s(gamma) + s(nu) + s(r_sigma)",
              strategy_coverage      = "s(dur_V) + s(dur_M) + s(kappa_V) + s(kappa_M)",
              strategy_costs         = "s(MV_dose_cost) + s(mAbs_dose_cost) + s(p_waste) + s(admin_cost)",
              health_costs           = "s(GP_cost) + s(N_GP) + s(Virology_cost) + s(NA_ED_cost) + s(ED_cost) + s(HP_cost) + s(ICU_cost)",
              health_utilities       = "s(GP_NA_ED_dis) + s(HP_dis) + s(ICU_dis) + s(GP_NA_ED_dur) + s(HP_dur) + s(ICU_dur)")
  EVPPI <- lapply(INMB, function(inmb){
              out <- lapply(models, function(mod)
                        tryCatch(evppi(INMB = inmb, psa = psa, model = mod)$EVPPI,
                                 error = function(e) evppi(INMB = inmb, psa = psa, model = gsub("\\)", ", k = 3)", mod))$EVPPI))
              names(out) <- names(models)
              out
  })
  names(EVPPI) <- names(INMB)
  saveRDS(list(cfg = cfg_par[i], EVPI = EVPI, EVPPI = EVPPI), paste0("Simulations/10K/", fname))
  rm(INMB, psa, EVPPI)
  gc()
}

cat("Starting EVSI value of information simulations at", as.character(Sys.time()), "\n")
for(i in 1:nrow(cfg_voi)){
  cat("Starting configuration", i, "of", nrow(cfg_voi), "at", as.character(Sys.time()), "\n")
  fname <- paste0("sims_voi_s_cfg_", formatC(i, width = 3, flag = "0"), ".rds")
  if(fname %in% list.files("Simulations/10K")){
    cat("Skipping configuration as EVSI value of information simulations already exist.\n")
    next
  }
  if(cfg_voi[i]$design == "A"){
    p_alloc <- c(0.5, 0.5, 0)
    evsi_mod <- "s(p_NS_MA) + s(p_MV_MA) + s(p_NS_HP) + s(p_MV_HP)"
  } else if(cfg_voi[i]$design == "B"){
    p_alloc <- c(0.5, 0, 0.5)
    evsi_mod <- "s(p_NS_MA) + s(p_mAbs_MA) + s(p_NS_HP) + s(p_mAbs_HP)"
  } else if(cfg_voi[i]$design == "C"){
    p_alloc <- c(0, 0.5, 0.5)
    evsi_mod <- "s(p_MV_MA) + s(p_mAbs_MA) + s(p_MV_HP) + s(p_mAbs_HP)"
  } else if(cfg_voi[i]$design == "D"){
    p_alloc <- c(1/3, 1/3, 1/3)
    evsi_mod <- "s(p_NS_MA) + s(p_MV_MA) + s(p_mAbs_MA) + s(p_NS_HP) + s(p_MV_HP) + s(p_mAbs_HP)"
  }
  alpha <- cfg_voi[i]$alpha
  N <- cfg_voi[i]$N
  strat_eff_parms_tmp <- strat_eff_parms_l[[which(cfg_par$alpha == alpha)]]
  INMB <- sims_cea[[which(cfg_par$alpha == alpha)]]$INMB
  cat("Estimating EVSI using the non-parametric regression method at", as.character(Sys.time()), "\n")
  p_d <- list("MA" = strat_eff_parms_tmp[["MA"]][,c("p_NS", "p_MV", "p_mAbs")],
              "HP" = strat_eff_parms_tmp[["HP"]][,c("p_NS", "p_MV", "p_mAbs")])
  dat_NP <- prep_dat_evsi(p_d = p_d, N = N, p_alloc = p_alloc, seed = seed)
  EVSI <- lapply(INMB, function(inmb)
             tryCatch(evsi(INMB = inmb, dat = dat_NP, method = "NP", model = evsi_mod),
                      error = function(e) evsi(INMB = inmb, dat = dat_NP, method = "NP", model = gsub("\\)", ", k = 5)", evsi_mod))))
  names(EVSI) <- names(INMB)
  saveRDS(list(cfg = cfg_voi[i], EVSI = EVSI), paste0("Simulations/10K/", fname))
  rm(strat_eff_parms_tmp, INMB, p_d, dat_NP)
  gc()
}
cat("Finished at", as.character(Sys.time()), "\n")
