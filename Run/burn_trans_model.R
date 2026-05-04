
## this script burns in the base transmission model for 500 years to find the stationary solution and saves the output

library(RSVVOI)
library(RSVModels)

N_draw <- 10000
ncores <- 16
trans_parms <- gen_trans_parms(N_draw = N_draw, seed = seed)
y0_burn <- burn_trans_model(N_draw = N_draw, trans_parms = trans_parms, ncores = ncores, save_output = TRUE)
rm(N_draw, ncores, trans_parms, y0_burn)
