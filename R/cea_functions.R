
#' @title eval_cea
#' @description Evaluates the health economic model.
#' @param N_draw Number of parameter draws to use. Should match the number of simulations used to estimate the incidence parameters.
#' @param inc_parms RSV incidence parameter distributions.
#' @param cea_parms Health economic parameter distributions.
#' @param ages Vector of ages in years used for transmission modelling.
#' @param minimal Logical indicator. If set to TRUE, then only minimal output is stored (i.e., no states). Default is FALSE.
#' @return List containing arrays for states (number of individuals in each health economic state), costs and QALYs.
#' @rdname eval_cea
#' @export
eval_cea <- function(N_draw, inc_parms, cea_parms, ages, minimal = FALSE){
  states <- lapply(c("NS", "MV", "mAbs"), function(d){
    inc <- inc_parms[[d]]$inc[1:N_draw,,]
    if(d == "NS"){
      D <- ICU <- HP <- ED <- MA <- SI <- array(NA, dim = dim(inc), dimnames = dimnames(inc))
      for(j in 1:length(ages)){
        SI[,,j]  <- inc[,,j]*cea_parms$p_SI_I[1:N_draw,j]
        MA[,,j]  <- SI[,,j]*cea_parms$p_MA_SI[1:N_draw,j]
        ED[,,j]  <- MA[,,j]*cea_parms$p_ED_MA[1:N_draw,j]
        HP[,,j]  <- ED[,,j]*cea_parms$p_HP_ED_NS[1:N_draw,j]
        ICU[,,j] <- HP[,,j]*cea_parms$p_ICU_HP[1:N_draw]
        D[,,j]   <- ICU[,,j]*cea_parms$p_D_ICU[1:N_draw,j]
      }
    } else {
      inc_strat <- inc_parms[[d]]$inc_strat[1:N_draw,,]
      inc_nstrat <- inc - inc_strat
      D_strat <- ICU_strat <- HP_strat <- ED_strat <- MA_strat <- SI_strat <- array(NA, dim = dim(inc), dimnames = dimnames(inc))
      D_nstrat <- ICU_nstrat <- HP_nstrat <- ED_nstrat <- MA_nstrat <- SI_nstrat <- array(NA, dim = dim(inc), dimnames = dimnames(inc))
      for(j in 1:length(ages)){
        SI_strat[,,j]   <- inc_strat[,,j]*cea_parms$p_SI_I[1:N_draw,j]
        SI_nstrat[,,j]  <- inc_nstrat[,,j]*cea_parms$p_SI_I[1:N_draw,j]
        MA_strat[,,j]   <- SI_strat[,,j]*cea_parms$p_MA_SI[1:N_draw,j]
        MA_nstrat[,,j]  <- SI_nstrat[,,j]*cea_parms$p_MA_SI[1:N_draw,j]
        ED_strat[,,j]   <- MA_strat[,,j]*cea_parms$p_ED_MA[1:N_draw,j]
        ED_nstrat[,,j]  <- MA_nstrat[,,j]*cea_parms$p_ED_MA[1:N_draw,j]
        HP_strat[,,j]   <- ED_strat[,,j]*cea_parms[[paste0("p_HP_ED_", d)]][1:N_draw,j]
        HP_nstrat[,,j]  <- ED_nstrat[,,j]*cea_parms$p_HP_ED_NS[1:N_draw,j]
        ICU_strat[,,j]  <- HP_strat[,,j]*cea_parms$p_ICU_HP[1:N_draw]
        ICU_nstrat[,,j] <- HP_nstrat[,,j]*cea_parms$p_ICU_HP[1:N_draw]
        D_strat[,,j]    <- ICU_strat[,,j]*cea_parms$p_D_ICU[1:N_draw,j]
        D_nstrat[,,j]   <- ICU_nstrat[,,j]*cea_parms$p_D_ICU[1:N_draw,j]
      }
      SI  <- SI_strat + SI_nstrat
      MA  <- MA_strat + MA_nstrat
      ED  <- ED_strat + ED_nstrat
      HP  <- HP_strat + HP_nstrat
      ICU <- ICU_strat + ICU_nstrat
      D   <- D_strat + D_nstrat
    }
    No_SI  <- inc - SI
    No_MA  <- SI - MA
    No_ED  <- MA - ED
    No_HP  <- ED - HP
    No_ICU <- HP - ICU
    No_D   <- ICU - D
    rm(inc, SI, MA, ED, HP, ICU)
    if(d != "NS") rm(inc_strat, inc_nstrat, SI_strat, SI_nstrat, MA_strat, MA_nstrat, ED_strat, ED_nstrat,
                     HP_strat, HP_nstrat, ICU_strat, ICU_nstrat, D_strat, D_nstrat)
    gc()
    return(list(No_SI = No_SI, No_MA = No_MA, No_ED = No_ED, No_HP = No_HP, No_ICU = No_ICU, No_D = No_D, D = D))
  })
  names(states) <- c("NS", "MV", "mAbs")

  costs <- lapply(c("NS", "MV", "mAbs"), function(d){
    states_tmp <- states[[d]]
    strategy_cost <- 0*inc_parms[[d]]$inc[1:N_draw,,]
    if(d != "NS"){
      strategy_cost[,,1] <- inc_parms[[d]]$admin[1:N_draw,]*(cea_parms[[paste0(d, "_dose_cost")]][1:N_draw]*(1 + cea_parms$p_waste[1:N_draw]) + cea_parms$admin_cost[1:N_draw])
    }
    GP_cost <- (states_tmp$No_ED + states_tmp$No_HP + states_tmp$No_ICU + states_tmp$No_D + states_tmp$D)*(cea_parms$N_GP[1:N_draw]*cea_parms$GP_cost[1:N_draw] + cea_parms$Virology_cost[1:N_draw])
    ED_cost <- states_tmp$No_HP*cea_parms$NA_ED_cost[1:N_draw] + (states_tmp$No_ICU + states_tmp$No_D + states_tmp$D)*cea_parms$ED_cost[1:N_draw]
    HP_cost <- (states_tmp$No_ICU + states_tmp$No_D + states_tmp$D)*cea_parms$HP_cost[1:N_draw]
    ICU_cost <- (states_tmp$No_D + states_tmp$D)*cea_parms$ICU_cost[1:N_draw]
    total_cost <- strategy_cost + GP_cost + ED_cost + HP_cost + ICU_cost
    total_cost_disc <- aperm(apply(total_cost, c(1,3), function(x) x/(1 + cea_parms$discount)^(1:dim(total_cost)[2] - 1)), c(2, 1, 3))
    return(total_cost_disc)
  })
  names(costs) <- c("NS", "MV", "mAbs")

  qalys <- lapply(c("NS", "MV", "mAbs"), function(d){
    states_tmp <- states[[d]]
    GP_NA_ED_qaly <- -((states_tmp$No_ED + states_tmp$No_HP + states_tmp$No_ICU + states_tmp$No_D + states_tmp$D)*cea_parms$GP_NA_ED_dis[1:N_draw]*cea_parms$GP_NA_ED_dur[1:N_draw])
    HP_qaly <- -((states_tmp$No_ICU + states_tmp$No_D + states_tmp$D)*cea_parms$HP_dis[1:N_draw]*cea_parms$HP_dur[1:N_draw])
    ICU_qaly <- -((states_tmp$No_D + states_tmp$D)*cea_parms$ICU_dis[1:N_draw]*cea_parms$ICU_dur[1:N_draw])
    alive_qaly <- GP_NA_ED_qaly + HP_qaly + ICU_qaly
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

  if(minimal){
    return(list(costs = costs, qalys = qalys))
  } else {
    return(list(states = states, costs = costs, qalys = qalys))
  }
}

#' @title eval_icer
#' @description Evaluates the incremental cost-effectiveness ratio (ICER) for the health economic model.
#' @param costs List containing arrays of costs (output from eval_cea() function).
#' @param qalys List containing arrays of QALYs (output from eval_cea() function).
#' @param curr Name of strategy to use as the comparator (e.g., the currently implemented strategy).
#' @param new Name of the strategy to use as the new strategy (e.g., to compare against the currently implemented strategy).
#' @param years Number of years to compute cumulative costs and QALYs. If NULL, then the maximum number of years will be used. Default is NULL.
#' @param max_age_ind Maximum age to compute cumulative costs and QALYs. If NULL, then the maximum number of ages will be used. Default is NULL.
#' @return List containing incremental costs, incremental QALYs and ICERs.
#' @rdname eval_icer
#' @export
eval_icer <- function(costs, qalys, curr, new, years = NULL, max_age_ind = NULL){
  if(is.null(years)) years <- dim(costs[[1]])[2]
  if(is.null(max_age_ind)) max_age_ind <- dim(costs[[1]])[3]
  cost_increm <- apply(costs[[new]][,1:years,1:max_age_ind] - costs[[curr]][,1:years,1:max_age_ind], 1, sum)
  qaly_increm <- apply(qalys[[new]][,1:years,1:max_age_ind] - qalys[[curr]][,1:years,1:max_age_ind], 1, sum)
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
#' @param max_age_ind Maximum age to compute cumulative costs and QALYs. If NULL, then the maximum number of ages will be used. Default is NULL.
#' @return Named list containing incremental net monetary benefit for each willingness-to-pay threshold.
#' @rdname eval_INMB
#' @export
eval_INMB <- function(cea, WTP, ref = "NS", years = NULL, max_age_ind = NULL){
  if(is.null(years)) years <- dim(cea$costs[[1]])[2]
  if(is.null(max_age_ind)) max_age_ind <- dim(cea$costs[[1]])[3]
  costs <- sapply(cea$costs, function(x) apply(x[,1:years,1:max_age_ind], 1, sum))
  qalys <- sapply(cea$qalys, function(x) apply(x[,1:years,1:max_age_ind], 1, sum))
  INMB <- lapply(WTP, function(wtp){
    NB <- wtp*qalys - costs
    NB - NB[,ref]
  })
  names(INMB) <- paste("WTP =", WTP)
  return(INMB)
}
