
## script to generate the statistical model hyperparameters

seed <- 35234

source("R/trial_functions.R")

hyper_parms <- list("MA" = gen_hyper_parms(outcome = "MA", seed = seed),
                    "HP" = gen_hyper_parms(outcome = "HP", seed = seed))
saveRDS(hyper_parms, "Parameters/hyper_parms.rds")

rm(list = ls())
