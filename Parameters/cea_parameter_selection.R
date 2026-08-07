
## the below code was run to assist the selection of the parameters for the probability distributions of the cost-effectiveness analysis parameters

library(nleqslv)
source("R/misc_functions.R")

## p_SI_I
## taken directly from Nazareno (2025)

## p_MA_SI
## taken directly from Nazareno (2025)

## p_ED_MA
## formula from Carlin (2026) cbind(c(seq(0, 5, 1/12), seq(10, 75, 5)), 0.0972 + exp(-1.53*(c(seq(0, 5, 1/12), seq(10, 75, 5))+0.581)))
comp_pars(dist = "beta", mn = 0.5  , upper = 0.6  , pb = 0.975) ## 0-1 months
comp_pars(dist = "beta", mn = 0.4  , upper = 0.5  , pb = 0.975) ## 2 months
comp_pars(dist = "beta", mn = 0.35 , upper = 0.45 , pb = 0.975) ## 3-5 months
comp_pars(dist = "beta", mn = 0.25 , upper = 0.35 , pb = 0.975) ## 6-11 months
comp_pars(dist = "beta", mn = 0.15 , upper = 0.225, pb = 0.975) ## 12-17 months
comp_pars(dist = "beta", mn = 0.125, upper = 0.175, pb = 0.975) ## 18-23 months
comp_pars(dist = "beta", mn = 0.1  , upper = 0.125, pb = 0.975) ## +24 months

## p_HP_SI_NS
comp_pars(dist = "beta", lower = 0.10 , upper = 0.25 , pb = c(0.025, 0.975)) ## 0-1 months
comp_pars(dist = "beta", lower = 0.05 , upper = 0.15 , pb = c(0.025, 0.975)) ## 2 months
comp_pars(dist = "beta", lower = 0.03 , upper = 0.09 , pb = c(0.025, 0.975)) ## 3-5 months
comp_pars(dist = "beta", lower = 0.01 , upper = 0.05 , pb = c(0.025, 0.975)) ## 6-11 months
comp_pars(dist = "beta", lower = 0.01 , upper = 0.03 , pb = c(0.025, 0.975)) ## 12-17 months
comp_pars(dist = "beta", lower = 0.005, upper = 0.02 , pb = c(0.025, 0.975)) ## 18-23 months
comp_pars(dist = "beta", lower = 0.001, upper = 0.01 , pb = c(0.025, 0.975)) ## +24 months

## p_HP_ED_NS
## derived using the stick breaking and moment matching method
derive_conditional(mu     = c(beta_prod_mean(x = c(50, 50) , y = c(20, 6))  , beta_mean(x = c(15, 75))),  ## 0-1 months
                   sigma2 = c(beta_prod_var(x  = c(50, 50) , y = c(20, 6))  , beta_var(x  = c(15, 75))))
derive_conditional(mu     = c(beta_prod_mean(x = c(40, 60) , y = c(20, 6))  , beta_mean(x = c(12, 120))), ## 2 months
                   sigma2 = c(beta_prod_var(x  = c(40, 60) , y = c(20, 6))  , beta_var(x  = c(12, 120))))
derive_conditional(mu     = c(beta_prod_mean(x = c(30, 60) , y = c(20, 6))  , beta_mean(x = c(10, 200))), ## 3-5 months
                   sigma2 = c(beta_prod_var(x  = c(30, 60) , y = c(20, 6))  , beta_var(x  = c(10, 200))))
derive_conditional(mu     = c(beta_prod_mean(x = c(20, 60) , y = c(40, 60)) , beta_mean(x = c(6, 200))),  ## 6-11 months
                   sigma2 = c(beta_prod_var(x  = c(20, 60) , y = c(40, 60)) , beta_var(x  = c(6, 200))))
derive_conditional(mu     = c(beta_prod_mean(x = c(15, 90) , y = c(40, 60)) , beta_mean(x = c(10, 700))), ## 12-17 months
                   sigma2 = c(beta_prod_var(x  = c(15, 90) , y = c(40, 60)) , beta_var(x  = c(10, 700))))
derive_conditional(mu     = c(beta_prod_mean(x = c(25, 170), y = c(40, 60)) , beta_mean(x = c(8, 750))),  ## 18-23 months
                   sigma2 = c(beta_prod_var(x  = c(25, 170), y = c(40, 60)) , beta_var(x  = c(8, 750))))
derive_conditional(mu     = c(beta_prod_mean(x = c(60, 550), y = c(40, 70)) , beta_mean(x = c(3, 800))),  ## 24 months - 4 years
                   sigma2 = c(beta_prod_var(x  = c(60, 550), y = c(40, 70)) , beta_var(x  = c(3, 800))))
derive_conditional(mu     = c(beta_prod_mean(x = c(60, 550), y = c(40, 100)), beta_mean(x = c(3, 800))),  ## 5-14 years
                   sigma2 = c(beta_prod_var(x  = c(60, 550), y = c(40, 100)), beta_var(x  = c(3, 800))))
derive_conditional(mu     = c(beta_prod_mean(x = c(60, 550), y = c(50, 200)), beta_mean(x = c(3, 800))),  ## +15 years
                   sigma2 = c(beta_prod_var(x  = c(60, 550), y = c(50, 200)), beta_var(x  = c(3, 800))))

## p_ICU_HP
## taken directly from Nazareno (2025)

## p_D_HP
## taken directly from Nazareno (2025)

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
