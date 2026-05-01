

#' @title gen_cea_parms
#' @description Generates the parameter distributions for the cost-effectiveness model parameters. Requires access to the Data folder containing data for the Australian life tables (life expectancy by age group) and the quality adjusted life year weights (QALYs by age group).
#' @param N_draw Number of parameter draws to characterise distributions.
#' @param seed Optional seed for replication. Default is NULL.
#' @param path Path to the location of the Data and Parameters folders. Default is "." and may be used if the function is run inside the package.
#' @param save_output Logical. If set to TRUE, then the generated distributions will be saved into the Parameters folder. Default is FALSE.
#' @return List containing cost-effectiveness model parameter distributions.
#' @rdname gen_cea_parms
#' @export
gen_cea_parms <- function(N_draw, seed = NULL, path = ".", save_output = FALSE){
  if(!is.null(seed)) set.seed(seed)
  lifetable <- read.csv(paste0(path, "/Data/australian_life_tables_111125.csv"))
  lifetable$all_lx <- rowMeans(lifetable[, c("males_lx", "females_lx")])
  lifetable$all_Lx <- rowMeans(lifetable[, c("males_Lx", "females_Lx")])
  age_qalys <- read.csv(paste0(path, "/Data/age_qalys.csv"))
  lifetable <- merge(lifetable, age_qalys, by = "age")
  cea_parms <- list(lifetable = lifetable,
                    p_HP_I_NS = rbeta(N_draw, 30, 2000),
                    p_ED_I = rbeta(N_draw, 0.8, 10),
                    p_GP_I = rbeta(N_draw, 1, 5),
                    p_ICU_HP = rbeta(N_draw, 15, 500),
                    p_D_ICU = rbeta(N_draw, 4, 1000),
                    MV_dose_cost = runif(N_draw, 75, 370),
                    mAbs_dose_cost = runif(N_draw, 290, 660),
                    p_waste = rbeta(N_draw, 5, 95),
                    admin_cost = rgamma(N_draw, 150, 4),
                    GP_cost = rgamma(N_draw, 16, 0.4),
                    N_GP = sample(1:3, N_draw, replace = TRUE),
                    NA_ED_cost = rgamma(N_draw, 15, 0.02),
                    Virology_cost = rgamma(N_draw, 15, 0.5),
                    ED_cost = rgamma(N_draw, 15, 0.01),
                    H_cost = rgamma(N_draw, 15, 0.002),
                    ICU_cost = rgamma(N_draw, 15, 0.0006),
                    NA_dis = rbeta(N_draw, 15, 80),
                    H_dis = rbeta(N_draw, 9, 13),
                    ICU_dis = rbeta(N_draw, 6, 4),
                    NA_dur = rgamma(N_draw, 15, 3)/365.25,
                    H_dur = rgamma(N_draw, 27, 6)/365.25,
                    ICU_dur = rgamma(N_draw, 40, 4)/365.25,
                    discount = 0.05)
  if(save_output) saveRDS(cea_parms, paste0(path, "/Parameters/cea_parms.rds"))
  return(cea_parms)
}

costs <- function(d_trans, cea_parms){
  n_MV <- d_trans$n_MV
  n_mAbs <- d_trans$n_mAbs
  n_GP <- d_trans$n_GP
  n_NA_ED <- d_trans$n_NA_ED
  n_H <- d_trans$n_H
  n_ICU <- d_trans$n_ICU
  year <- d_trans$year

  MV_dose_cost <- cea_parms$MV_dose_cost
  mAbs_dose_cost <- cea_parms$mAbs_dose_cost
  p_waste <- cea_parms$p_waste
  admin_cost <- cea_parms$admin_cost
  discount <- cea_parms$discount
  GP_cost <- cea_parms$GP_cost
  N_GP <- cea_parms$N_GP
  NA_ED_cost <- cea_parms$NA_ED_cost
  Virology_cost <- cea_parms$Virology_cost
  ED_cost <- cea_parms$ED_cost
  H_cost <- cea_parms$H_cost
  ICU_cost <- cea_parms$ICU_cost

  strategy_cost <- (n_MV*MV_dose_cost + n_mAbs*mAbs_dose_cost)*(1 + p_waste) + (n_MV + n_mAbs)*admin_cost
  outpatient_cost <- n_GP*(N_GP*GP_cost + Virology_cost) + n_NA_ED*NA_ED_cost
  inpatient_cost <- n_H*(H_cost + ED_cost) + n_ICU*(ICU_cost + ED_cost)
  total_cost <- (strategy_cost + outpatient_cost + inpatient_cost)/(1 + discount)^(year - 1)
  return(total_cost)
}

qalys <- function(age, d_trans, cea_parms, lifetable){
  n_GP <- d_trans$n_GP
  n_NA_ED <- d_trans$n_NA_ED
  n_H <- d_trans$n_H
  n_ICU <- d_trans$n_ICU
  n_D <- d_trans$n_D
  year <- d_trans$year

  NA_dis <- cea_parms$NA_dis
  H_dis <- cea_parms$H_dis
  ICU_dis <- cea_parms$ICU_dis
  NA_dur <- cea_parms$NA_dur
  H_dur <- cea_parms$H_dur
  ICU_dur <- cea_parms$ICU_dur
  discount <- cea_parms$discount

  dead_qaly <- sapply(1:length(age), function(i){
    fut_Ql <- lifetable[lifetable$age >= age[i], "qaly"]
    fut_Lx <- lifetable[lifetable$age >= age[i], "all_Lx"]
    fut_age <- lifetable[lifetable$age >= age[i], "age"]
    curr_lx <- lifetable[lifetable$age == age[i], "all_lx"]
    -n_D[i]*sum(fut_Ql*fut_Lx/((1 + discount)^(fut_age - age[i])))/curr_lx
  })

  outpatient_qaly <- -((n_GP + n_NA_ED)*NA_dis*NA_dur)
  inpatient_qaly <- -(n_H*H_dis*H_dur + n_ICU*ICU_dis*ICU_dur)
  alive_qaly <- (outpatient_qaly + inpatient_qaly)/(1 + discount)^(year - 1)
  total_qaly <- alive_qaly + dead_qaly
  return(total_qaly)
}

## cea function to be run per scenario (d_trans and cea_parms are both uncertain with the same lengths)

cea <- function(ages, d_trans, cea_parms, lifetable){
  cea_l <- lapply(ages, function(age){
    cea_costs <- costs(d_trans = d_trans[d_trans$age == age,], cea_parms = cea_parms)
    cea_qalys <- qalys(age = age, d_trans = d_trans[d_trans$age == age,], cea_parms = cea_parms, lifetable = lifetable)
    data.frame(sim = 1:nrow(d_trans[d_trans$age == age,]), age = age, costs = cea_costs, qalys = cea_qalys)
  })
  cea_df <- do.call(rbind, cea_l)
  cea_df <- aggregate(cbind(costs, qalys) ~ sim, data = cea_df, FUN = sum)
  return(cea_df)
}

  # mutate(
  #   cost_0 = costAll[strategy == reference_strategy],
  #   qaly_0 = qalyAll[strategy == reference_strategy],
  #   yll_0  = yll[strategy == reference_strategy],
  #   iCost = costAll - cost_0,
  #   iQaly = qalyAll - qaly_0,
  #   icer  = iCost / iQaly)



