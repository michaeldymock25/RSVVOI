
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

#' @title gen_strat_eff_parms
#' @description Generates the strategy effectiveness parameter distributions.
#' @param N_draw Number of parameter draws to characterise distributions.
#' @param hyper_parms Hyperparameters. Must be a named list ("NS", "MV", "mAbs") with named elements ("a", "b") for "NS" and ("mu", "sigma") for "MV" and "mAbs".
#' @param alpha Power prior contribution. Must be a scalar between zero and 1 (inclusive).
#' @param OR Logical. If set to TRUE, parameters will be transformed as if they were odds ratios. Otherwise, parameters will be transformed as if they were relative risks. Default is TRUE.
#' @param seed Optional seed for replication. Default is NULL.
#' @return Strategy effectiveness parameter distributions.
#' @rdname gen_strat_eff_parms
#' @export
gen_strat_eff_parms <- function(N_draw, hyper_parms, alpha, OR = TRUE, seed = NULL){
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
  strat_eff_parms <- cbind(p_NS = p_NS, p_MV = p_MV, p_mAbs = p_mAbs, pi_MV = pi_MV, pi_mAbs = pi_mAbs, r_MV = r_MV, r_mAbs = r_mAbs)
  return(strat_eff_parms)
}

#' @title gen_trans_parms
#' @description Generates the parameter distributions for the transmission model parameters. Excludes strategy effectiveness parameters. Requires access to the Data folder containing data for the Australian population (size by age group) and the estimated mixing matrix (computed via the RSVModels package).
#' @param N_draw Number of parameter draws to characterise distributions.
#' @param seed Optional seed for replication. Default is NULL.
#' @param Data_path Path to the location of the Data folder. Default is "." and may be used if the function is run inside the package.
#' @param save_output Logical. If set to TRUE, then the generated distributions will be saved into the Parameters folder. Default is FALSE.
#' @param save_path Path to the location of the Parameters folder (or other location where output is to be saved). Default is "Parameters/trans_parms.rds" and may be used if the function is run inside the package. Only required if saving the output.
#' @return List containing transmission model parameter distributions.
#' @rdname gen_trans_parms
#' @export
gen_trans_parms <- function(N_draw, seed = NULL, Data_path = ".", save_output = FALSE, save_path = "Parameters/trans_parms.rds"){
  if(!is.null(seed)) set.seed(seed)
  pop_agg <- read.csv(paste0(Data_path, "/Data/australian_population_aggregated.csv"))
  mixing <- read.csv(paste0(Data_path, "/Data/mixing.csv"))
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
  if(save_output) saveRDS(trans_parms, save_path)
  return(trans_parms)
}

#' @title cea_parms_age_helper
#' @description Helper function to generate the health economic parameter distributions. Allows for age-specific distributions.
#' @param N_draw Number of parameter draws to characterise distributions.
#' @param par_nm Name of the parameter to simulate.
#' @return Named matrix containing age-specific parameter distributions.
#' @rdname cea_parms_age_helper
#' @export
cea_parms_age_helper <- function(N_draw, par_nm){
  ages <- c(seq(0, 5, 1/12), seq(10, 75, 5))
  par_mat <- matrix(NA, nrow = N_draw, ncol = length(ages), dimnames = list(NULL, "age" = ages))
  if(par_nm == "p_SI_I"){
    pars <- list(rbeta(N_draw, 20, 10) ,
                 rbeta(N_draw, 30, 30) ,
                 rbeta(N_draw, 50, 150),
                 rbeta(N_draw, 30, 25) ,
                 rbeta(N_draw, 7 , 1))
    inds <- rep(1:5, times = c(60, 2, 5, 5, 3))
  } else if(par_nm == "p_MA_SI"){
    pars <- list(rbeta(N_draw, 20, 6)  ,
                 rbeta(N_draw, 40, 60) ,
                 rbeta(N_draw, 40, 70) ,
                 rbeta(N_draw, 40, 100),
                 rbeta(N_draw, 50, 200))
    inds <- rep(1:5, times = c(6, 18, 36, 2, 13))
  } else if(par_nm == "p_ED_MA"){
    pars <- list(rbeta(N_draw, 50, 50) ,
                 rbeta(N_draw, 40, 60) ,
                 rbeta(N_draw, 30, 60) ,
                 rbeta(N_draw, 20, 60) ,
                 rbeta(N_draw, 15, 90) ,
                 rbeta(N_draw, 25, 170),
                 rbeta(N_draw, 60, 550))
    inds <- rep(1:7, times = c(2, 1, 3, 6, 6, 6, 51))
  } else if(par_nm == "p_HP_ED_NS"){
    pars <- list(rbeta(N_draw, 17, 22),
                 rbeta(N_draw, 15, 35),
                 rbeta(N_draw, 15, 60),
                 rbeta(N_draw, 6 , 16),
                 rbeta(N_draw, 30, 95),
                 rbeta(N_draw, 10, 40),
                 rbeta(N_draw, 3 , 25),
                 rbeta(N_draw, 3 , 20),
                 rbeta(N_draw, 3 , 12))
    inds <- rep(1:9, times = c(2, 1, 3, 6, 6, 6, 36, 2, 13))
  } else if(par_nm == "p_D_ICU"){
    pars <- list(rbeta(N_draw, 60, 100000),
                 rbeta(N_draw, 60, 2500),
                 rbeta(N_draw, 60, 1000))
    inds <- rep(1:3, times = c(60, 12, 3))
  }
  for(j in 1:length(ages)) par_mat[,j] <- pars[[inds[j]]]
  return(par_mat)
}

#' @title gen_cea_parms
#' @description Generates the parameter distributions for the health economic model parameters. Requires access to the Data folder containing data for the Australian life tables (life expectancy by age group) and the quality adjusted life year weights (QALYs by age group).
#' @param N_draw Number of parameter draws to characterise distributions.
#' @param seed Optional seed for replication. Default is NULL.
#' @param Data_path Path to the location of the Data folder. Default is "." and may be used if the function is run inside the package.
#' @param save_output Logical. If set to TRUE, then the generated distributions will be saved into the Parameters folder. Default is FALSE.
#' @param save_path Path to the location of the Parameters folder (or other location where output is to be saved). Default is "Parameters/cea_parms.rds" and may be used if the function is run inside the package. Only required if saving the output.
#' @return List containing cost-effectiveness model parameter distributions.
#' @rdname gen_cea_parms
#' @export
gen_cea_parms <- function(N_draw, seed = NULL, Data_path = ".", save_output = FALSE, save_path = "Parameters/cea_parms.rds"){
  if(!is.null(seed)) set.seed(seed)
  lifetable <- read.csv(paste0(Data_path, "/Data/australian_life_tables_111125.csv"))
  lifetable$all_lx <- rowMeans(lifetable[, c("males_lx", "females_lx")])
  lifetable$all_Lx <- rowMeans(lifetable[, c("males_Lx", "females_Lx")])
  age_qalys <- read.csv(paste0(Data_path, "/Data/age_qalys.csv"))
  lifetable <- merge(lifetable, age_qalys, by = "age")
  cea_parms <- list(lifetable = lifetable,
                    p_SI_I = cea_parms_age_helper(N_draw, par_nm = "p_SI_I"),
                    p_MA_SI = cea_parms_age_helper(N_draw, par_nm = "p_MA_SI"),
                    p_ED_MA = cea_parms_age_helper(N_draw, par_nm = "p_ED_MA"),
                    p_HP_ED_NS = cea_parms_age_helper(N_draw, par_nm = "p_HP_ED_NS"),
                    p_ICU_HP = rbeta(N_draw, 15, 400),
                    p_D_ICU = cea_parms_age_helper(N_draw, par_nm = "p_D_ICU"),
                    MV_dose_cost = 65 + rgamma(N_draw, 2, 0.02),
                    mAbs_dose_cost = 232 + rgamma(N_draw, 0.4, 0.007),
                    p_waste = rbeta(N_draw, 5, 95),
                    admin_cost = rgamma(N_draw, 150, 4),
                    GP_cost = rgamma(N_draw, 16, 0.4),
                    N_GP = sample(1:3, N_draw, replace = TRUE),
                    Virology_cost = rgamma(N_draw, 15, 0.5),
                    NA_ED_cost = rgamma(N_draw, 15, 0.02),
                    ED_cost = rgamma(N_draw, 15, 0.01),
                    HP_cost = rgamma(N_draw, 15, 0.002),
                    ICU_cost = rgamma(N_draw, 15, 0.0006),
                    GP_NA_ED_dis = rbeta(N_draw, 15, 80),
                    HP_dis = rbeta(N_draw, 9, 13),
                    ICU_dis = rbeta(N_draw, 6, 4),
                    GP_NA_ED_dur = rgamma(N_draw, 15, 3)/365.25,
                    HP_dur = rgamma(N_draw, 27, 6)/365.25,
                    ICU_dur = rgamma(N_draw, 40, 4)/365.25,
                    discount = 0.05)
  if(save_output) saveRDS(cea_parms, save_path)
  return(cea_parms)
}
