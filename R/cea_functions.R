
#' @title gen_cea_parms
#' @description Generates the parameter distributions for the health economic model parameters. Requires access to the Data folder containing data for the Australian life tables (life expectancy by age group) and the quality adjusted life year weights (QALYs by age group).
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
                    HP_cost = rgamma(N_draw, 15, 0.002),
                    ICU_cost = rgamma(N_draw, 15, 0.0006),
                    NA_dis = rbeta(N_draw, 15, 80),
                    HP_dis = rbeta(N_draw, 9, 13),
                    ICU_dis = rbeta(N_draw, 6, 4),
                    NA_dur = rgamma(N_draw, 15, 3)/365.25,
                    HP_dur = rgamma(N_draw, 27, 6)/365.25,
                    ICU_dur = rgamma(N_draw, 40, 4)/365.25,
                    discount = 0.05)
  if(save_output) saveRDS(cea_parms, paste0(path, "/Parameters/cea_parms.rds"))
  return(cea_parms)
}

#' @title eval_cea
#' @description Evaluates the health economic model.
#' @param N_draw Number of parameter draws to use. Should match the number of simulations used to estimate the incidence parameters.
#' @param inc_parms RSV incidence parameter distributions.
#' @param cea_parms Health economic parameter distributions.
#' @param ages Vector of ages in years used for transmission modelling.
#' @return List containing arrays for states (number of individuals in each health economic state), costs and QALYs.
#' @rdname eval_cea
#' @export
eval_cea <- function(N_draw, inc_parms, cea_parms, ages){
  states <- lapply(c("NS", "MV", "mAbs"), function(d){
    inc <- inc_parms[[d]]$inc[1:N_draw,,]
    HP <- inc*cea_parms[[paste0("p_HP_I_", d)]][1:N_draw]
    ED <- inc*cea_parms$p_ED_I[1:N_draw]
    GP <- inc*cea_parms$p_GP_I[1:N_draw]
    No_MA <- inc - HP - ED - GP
    ICU <- HP*cea_parms$p_ICU_HP[1:N_draw]
    No_ICU <- HP - ICU
    D <- ICU*cea_parms$p_D_ICU[1:N_draw]
    No_D <- ICU - D
    return(list(No_MA = No_MA, GP = GP, ED = ED, No_ICU = No_ICU, No_D = No_D, D = D))
  })
  names(states) <- c("NS", "MV", "mAbs")

  costs <- lapply(c("NS", "MV", "mAbs"), function(d){
    states_tmp <- states[[d]]
    if(d == "NS"){
      strategy_cost <- 0*inc_parms[[d]]$inc[1:N_draw,,]
    } else {
      strategy_cost <- inc_parms[[d]]$admin[1:N_draw,,]*(cea_parms[[paste0(d, "_dose_cost")]][1:N_draw]*(1 + cea_parms$p_waste[1:N_draw]) + cea_parms$admin_cost[1:N_draw])
    }
    GP_cost <- (states_tmp$GP + states_tmp$ED + states_tmp$No_ICU + states_tmp$No_D + states_tmp$D)*(cea_parms$N_GP[1:N_draw]*cea_parms$GP_cost[1:N_draw] + cea_parms$Virology_cost[1:N_draw])
    ED_cost <- states_tmp$ED*cea_parms$NA_ED_cost[1:N_draw]
    No_ICU_cost <- (states_tmp$No_ICU + states_tmp$No_D + states_tmp$D)*(cea_parms$ED_cost[1:N_draw] + cea_parms$HP_cost[1:N_draw])
    No_D_cost <- (states_tmp$No_D + states_tmp$D)*cea_parms$ICU_cost[1:N_draw]
    total_cost <- strategy_cost + GP_cost + ED_cost + No_ICU_cost + No_D_cost
    total_cost_disc <- aperm(apply(total_cost, c(1,3), function(x) x/(1 + cea_parms$discount)^(1:dim(total_cost)[2] - 1)), c(2, 1, 3))
    return(total_cost_disc)
  })
  names(costs) <- c("NS", "MV", "mAbs")

  qalys <- lapply(c("NS", "MV", "mAbs"), function(d){
    states_tmp <- states[[d]]
    outpatient_qaly <- -((states_tmp$GP + states_tmp$ED + states_tmp$No_ICU + states_tmp$No_D + states_tmp$D)*cea_parms$NA_dis[1:N_draw]*cea_parms$NA_dur[1:N_draw])
    No_ICU_qaly <- -((states_tmp$No_ICU + states_tmp$No_D + states_tmp$D)*cea_parms$HP_dis[1:N_draw]*cea_parms$HP_dur[1:N_draw])
    No_D_qaly <- -((states_tmp$No_D + states_tmp$D)*cea_parms$ICU_dis[1:N_draw]*cea_parms$ICU_dur[1:N_draw])
    alive_qaly <- outpatient_qaly + No_ICU_qaly + No_D_qaly
    alive_qaly_disc <- aperm(apply(alive_qaly, c(1,3), function(x) x/(1 + cea_parms$discount)^(1:dim(alive_qaly)[2] - 1)), c(2, 1, 3))
    D_age_disc <- sapply(1:length(ages), function(i){
      fut_Ql <- cea_parms$lifetable[cea_parms$lifetable$age >= ages[i], "qaly"]
      fut_Lx <- cea_parms$lifetable[cea_parms$lifetable$age >= ages[i], "all_Lx"]
      fut_age <- cea_parms$lifetable[cea_parms$lifetable$age >= ages[i], "age"]
      curr_lx <- cea_parms$lifetable[cea_parms$lifetable$age == round(ages[i]), "all_lx"]
      sum(fut_Ql*fut_Lx/(1 + cea_parms$discount)^(fut_age - ages[i]))/curr_lx
    })
    dead_qaly_disc <- aperm(aperm(states_tmp$D, c(3, 1, 2))*-D_age_disc, c(2, 3, 1))
    total_qaly_disc <- alive_qaly_disc + dead_qaly_disc
    return(total_qaly_disc)
  })
  names(qalys) <- c("NS", "MV", "mAbs")

  return(list(states = states, costs = costs, qalys = qalys))
}

#' @title eval_icer
#' @description Evaluates the incremental cost-effectiveness ratio (ICER) for the health economic model.
#' @param costs List containing arrays of costs (output from eval_cea() function).
#' @param qalys List containing arrays of QALYs (output from eval_cea() function).
#' @param curr Name of strategy to use as the comparator (e.g., the currently implemented strategy).
#' @param new Name of the strategy to use as the new strategy (e.g., to compare against the currently implemented strategy).
#' @param years Number of years to compute cumulative costs and QALYs. If NULL, then the maximum number of years will be used. Default is NULL.
#' @return List containing incremental costs, incremental QALYs and ICERs.
#' @rdname eval_icer
#' @export
eval_icer <- function(costs, qalys, curr, new, years = NULL){
  if(is.null(years)) years <- dim(costs[[1]])[2]
  cost_increm <- apply(costs[[new]][,1:years,] - costs[[curr]][,1:years,], 1, sum)
  qaly_increm <- apply(qalys[[new]][,1:years,] - qalys[[curr]][,1:years,], 1, sum)
  icer <- cost_increm/qaly_increm
  return(list(cost_increm = cost_increm, qaly_increm = qaly_increm, icer = icer))
}

#' @title plot_cea
#' @import ggplot2
#' @description Plots the distribution of the incremental costs against the incremental QALYs for the health economic model.
#' @param cea_increm List containing incremental costs and incremental QALYs for each comparison (output from eval_icer() function).
#' @param WTP Vector containing willingness-to-pay thresholds to plot. If NULL, then no thresholds will be included. Default is NULL.
#' @param sub Vector containing subset of comparisons to visualise. If NULL, then all comparisons will be included. Default is NULL.
#' @return ggplot object visualising the incremental costs against the incremental QALYs for each comparison.
#' @rdname plot_cea
#' @export
plot_cea <- function(cea_increm, WTP = NULL, sub = NULL){
  cost_increm <- sapply(cea_increm, function(x) x$cost_increm/1000000)
  qaly_increm <- sapply(cea_increm, function(x) x$qaly_increm)
  dat_vis <- data.frame(comp = factor(rep(colnames(cost_increm), each = nrow(cost_increm)), levels = colnames(cost_increm)),
                        cost_increm = as.vector(cost_increm),
                        qaly_increm = as.vector(qaly_increm))
  if(!is.null(sub)) dat_vis <- droplevels(dat_vis[dat_vis$comp %in% sub,])
  p <- ggplot(dat_vis, aes(x = qaly_increm, y = cost_increm)) +
    facet_wrap(~comp, scales = "free") +
    geom_point(alpha = 0.6) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    xlab("Incremental QALY") +
    ylab("Incremental cost (in million $AUD)")
  if(!is.null(WTP)){
    dat_slope <- data.frame(slope = WTP/1000000,
                            WTP = factor(paste0("$", WTP/1000, "K"), levels = paste0("$", WTP/1000, "K")))
    p <- p + geom_abline(data = dat_slope, mapping = aes(intercept = 0, slope = slope, colour = WTP), linetype = "dashed")
  }
  return(p)
}

#' @title eval_INMB
#' @description Evaluates the incremental net monetary benefit for the health economic model.
#' @param cea List containing costs and QALYs for each strategy (output from eval_cea() function).
#' @param WTP Vector containing willingness-to-pay thresholds.
#' @param ref Reference group to use for the increments. Default is "NS".
#' @param years Number of years to compute cumulative costs and QALYs. If NULL, then the maximum number of years will be used. Default is NULL.
#' @return Named list containing incremental net monetary benefit for each willingness-to-pay threshold.
#' @rdname eval_INMB
#' @export
eval_INMB <- function(cea, WTP, ref = "NS", years = NULL){
  if(is.null(years)) years <- dim(cea$costs[[1]])[2]
  costs <- sapply(cea$costs, function(x) apply(x[,1:years,], 1, sum))
  qalys <- sapply(cea$qalys, function(x) apply(x[,1:years,], 1, sum))
  INMB <- lapply(WTP, function(wtp){
    NB <- wtp*qalys - costs
    NB - NB[,ref]
  })
  names(INMB) <- paste("WTP =", WTP)
  return(INMB)
}
