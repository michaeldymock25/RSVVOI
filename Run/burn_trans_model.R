
## this script burns in the base transmission model for 500 years to find the stationary solution and saves the output

library(RSVModels)

source("R/trans_functions.R")

N_draw <- 100000
seed <- 35234
ncores <- 8
trans_parms <- gen_trans_parms(N_draw = N_draw, seed = seed)
y0_burn <- burn_trans_model(N_draw = N_draw, trans_parms = trans_parms, burn_batch_size = 20, ncores = ncores, save_output = TRUE)

rm(N_draw, seed, ncores, trans_parms, y0_burn)
gc()
