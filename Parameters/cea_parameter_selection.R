
## the below code was run to assist the selection of the parameters for the probability distributions of the cost-effectiveness analysis parameters

source("R/misc_functions.R")

## p_HP_I_NS
comp_pars(dist = "beta", lower = 0.01, upper = 0.02, pb = c(0.025, 0.975))

## p_ED_I
comp_pars(dist = "beta", lower = 0.001, upper = 0.25, pb = c(0.025, 0.975))

## p_GP_I
comp_pars(dist = "beta", lower = 0.01, upper = 0.50, pb = c(0.025, 0.975))

## p_ICU_HP
comp_pars(dist = "beta", lower = 0.02, upper = 0.05, pb = c(0.025, 0.975))

## p_D_ICU
comp_pars(dist = "beta", lower = 0.0012, upper = 0.0084, pb = c(0.025, 0.975))

## p_waste
comp_pars(dist = "beta", mn = 0.05, upper = 0.10, pb = 0.975)

## admin_cost
comp_pars(dist = "gamma", mn = 36.95, upper = 42.85, pb = 0.975)

## GP_cost
comp_pars(dist = "gamma", mn = 41.4, upper = 62.1, pb = 0.975)

## Virology_cost
comp_pars(dist = "gamma", mn = 28.65, upper = 42.975, pb = 0.975)

## NA_ED_cost
comp_pars(dist = "gamma", mn = 724, upper = 1086, pb = 0.975)

## ED_cost
comp_pars(dist = "gamma", mn = 1335, upper = 2000, pb = 0.975)

## HP_cost
comp_pars(dist = "gamma", mn = 7980, upper = 11970, pb = 0.975)

## ICU_cost
comp_pars(dist = "gamma", mn = 25490, upper = 38235, pb = 0.975)

## NA_dis
comp_pars(dist = "beta", mn = 0.16, upper = 0.24, pb = 0.975)

## HP_dis
comp_pars(dist = "beta", mn = 0.41, upper = 0.615, pb = 0.975)

## ICU_dis
comp_pars(dist = "beta", mn = 0.60, upper = 0.9, pb = 0.975)

## NA_dur
comp_pars(dist = "gamma", mn = 5, upper = 8, pb = 0.975)

## HP_dur
comp_pars(dist = "gamma", mn = 4.4, upper = 6.2, pb = 0.975)

## ICU_dur
comp_pars(dist = "gamma", mn = 9.8, upper = 13.6, pb = 0.975)
