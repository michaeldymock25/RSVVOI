
## this script visualises some of the outputs for the manuscript

rm(list = ls())

## read in libraries and load functions and simulation outputs

library(ggplot2)
library(data.table)

source("R/psa_functions.R")
source("R/cea_functions.R")

N_draw <- 10000
years_sub <- c(1:5, 10, 15)
WTP <- c(20000, 50000)

sims_par <- lapply(list.files(path = "Simulations/10K", pattern = "sims_par", full.names = TRUE), readRDS)
sims_cea <- lapply(list.files(path = "Simulations/10K", pattern = "sims_cea", full.names = TRUE), readRDS)
sims_voi_p <- lapply(list.files(path = "Simulations/10K", pattern = "sims_voi_p", full.names = TRUE), readRDS)
sims_voi_s <- lapply(list.files(path = "Simulations/10K", pattern = "sims_voi_s", full.names = TRUE), readRDS)

## Main text

### Figure 1

nodes <- data.frame(x     = c(                             4,
                                                       3,
                                                   2,      4,
                                               1,      3,
                                           0,      2,
                                      -1,      1,
                                  -2,      0,
                              -3,     -1,
                                  -2),
                    y     = c(                             4,
                                                       3,
                                                   2,      2,
                                               1,      1,
                                           0,      0,
                                      -1,     -1,
                                  -2,     -2,
                              -3,     -3,
                                  -4),
                    shape = c(                            22,
                                                      21,
                                                  21,     22,
                                              21,     22,
                                          21,     22,
                                      21,     22,
                                  21,     22,
                              23,     22,
                                  22))
nodes$type <- factor(ifelse(nodes$shape == 22, "Terminal", ifelse(nodes$shape == 21, "Chance", "Decision")),
                     levels = c("Decision", "Chance", "Terminal"))
edges <- data.frame(x    = c( 3,  3,
                              2,  2,
                              1,  1,
                              0,  0,
                             -1, -1,
                             -2, -2,
                             -3, -3),
                    y    = c( 3,  3,
                              2,  2,
                              1,  1,
                              0,  0,
                             -1, -1,
                             -2, -2,
                             -3, -3),
                    xend = c( 4,  4,
                              3,  3,
                              2,  2,
                              1,  1,
                              0,  0,
                             -1, -1,
                             -2, -2),
                    yend = c( 4,  2,
                              3,  1,
                              2,  0,
                              1, -1,
                              0, -2,
                             -1, -3,
                             -2, -4),
                    type = c( 1,  1,
                              1,  1,
                              1,  1,
                              1,  1,
                              1,  1,
                              1,  1,
                              2,  2))
edges$type <- ifelse(edges$type == 1, "solid", "dashed")
node_text <- data.frame(text = c(                                                                     "D",
                                                                                            "ICU",
                                                                                   "HP",           "No D",
                                                                          "ED",          "No ICU",
                                                                 "MA",          "No HP",
                                                        "SI",          "No ED",
                                                "I",          "No MA",
                                 "Decision",         "No SI",
                                             "No I"),
                        x    = c(                                                                       4,
                                                                                                3,
                                                                                      2,                4,
                                                                             1,                 3,
                                                                    0,                2,
                                                          -1,                1,
                                                 -2,                0,
                                       -3.2,              -1,
                                                 -2),
                        y    = c(                                                                     4.4,
                                                                                              3.4,
                                                                                    2.4,              2.4,
                                                                           1.4,               1.4,
                                                                  0.4,              0.4,
                                                        -0.6,             -0.6,
                                               -1.6,             -1.6,
                                       -2.6,            -2.6,

                                               -3.6))
edge_text <- data.frame(
  text  = c(
                                                                                                                 "$p^{D|ICU}$",
                                                                                             "$p^{ICU|HP}$", "$1 - p^{D|ICU}$",
                                                                        "$p_d^{HP|ED}$", "$1 - p^{ICU|HP}$",
                                                     "$p^{ED|MA}$", "$1 - p_d^{HP|ED}$",
                                  "$p^{MA|SI}$", "$1 - p^{ED|MA}$",
                "$p^{SI|I}$", "$1 - p^{MA|SI}$",
            "$1 - p^{SI|I}$"),
  x     = c(                                                                                                               3.5,
                                                                                                        2.5,               3.4,
                                                                                    1.5,                2.4,
                                                               0.5,                 1.4,
                                           -0.5,               0.4,
                        -1.5,              -0.6,
                        -1.6),
  y     = c(                                                                                                               3.8,
                                                                                                        2.8,               2.3,
                                                                                    1.8,                1.3,
                                                               0.8,                 0.3,
                                           -0.2,              -0.8,
                        -1.2,              -1.8,
                        -2.8),
 angle = c(                                                                                                                 32,
                                                                                                         32,               -38,
                                                                                     32,                -38,
                                                                32,                 -38,
                                              32,              -38,
                         32,                 -38,
                        -38))

fig_1 <- ggplot() +
  geom_segment(data = edges[edges$type == "solid",], aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_segment(data = edges[edges$type == "dashed",], aes(x = x, y = y, xend = xend, yend = yend), linetype = "dashed") +
  geom_point(data = nodes, aes(x = x, y = y, shape = shape, colour = type, fill = type), size = 4) +
  geom_text(data = edge_text, aes(x = x, y = y, label = latex2exp::TeX(text), angle = angle), size = 2.6, parse = TRUE) +
  geom_text(data = node_text, aes(x = x, y = y, label = text), size = 2.6) +
  scale_shape_identity() +
  scale_fill_manual("", values = c("#88CCEE", "#117733", "#661100")) +
  scale_colour_manual("", values = c("#88CCEE", "#117733", "#661100")) +
  theme_void() +
  theme(legend.position = "bottom")

tiff("Figure 1.tiff", width = 6, height = 5, units = "in", res = 1000)
fig_1
dev.off()

png("Figure 1.png", width = 6, height = 5, units = "in", res = 200)
fig_1
dev.off()

### Figure 2

dat_par_vis <- rbindlist(lapply(sims_par, function(par_sim)
  data.table(alpha     = paste("alpha ==", par_sim$cfg$alpha),
             incidence = unlist(lapply(par_sim$inc_parms, function(x) x$inc)),
             strategy  = rep(c("NS", "MV", "mAbs"), each = N_draw*max(years_sub)*75),
             sim       = rep(1:N_draw, max(years_sub)*75*3),
             years_n   = rep(rep(1:max(years_sub), each = N_draw), 75*3),
             age       = as.numeric(rep(rep(dimnames(par_sim$inc_parms$NS$inc)[[3]], each = N_draw*max(years_sub)), 3)))))
dat_par_vis <- dcast.data.table(dat_par_vis, alpha + sim + years_n + age ~ strategy, value.var = "incidence")
#dat_par_vis_abs <- copy(dat_par_vis)
dat_par_vis_rel <- copy(dat_par_vis)
#dat_par_vis_abs[, `:=`(MV = NS - MV, mAbs = NS - mAbs)]
dat_par_vis_rel[, `:=`(MV = (NS - MV)/NS, mAbs = (NS - mAbs)/NS)]
#dat_par_vis_abs <- melt.data.table(dat_par_vis_abs, id.vars = c("alpha", "sim", "years_n", "age"), measure.vars = c("MV", "mAbs"))
dat_par_vis_rel <- melt.data.table(dat_par_vis_rel, id.vars = c("alpha", "sim", "years_n", "age"), measure.vars = c("MV", "mAbs"))
#dat_par_vis_abs <- dat_par_vis_abs[, .(value = mean(value)), by = .(alpha, sim, age, variable)]
dat_par_vis_rel <- dat_par_vis_rel[, .(value = mean(value)), by = .(alpha, sim, age, variable)]
#dat_par_vis_abs <- dat_par_vis_abs[, .(med = median(value), low = quantile(value, 0.025), upp = quantile(value, 0.975)),
#                                   by = .(alpha, age, variable)]
dat_par_vis_rel <- dat_par_vis_rel[, .(med = median(value), low = quantile(value, 0.025), upp = quantile(value, 0.975)),
                                   by = .(alpha, age, variable)]

# fig_2a <- ggplot(dat_par_vis_abs[age <= 1.5], aes(x = age, y = med, ymin = low, ymax = upp, colour = variable)) +
#   facet_wrap(~ alpha, labeller = labeller(alpha = label_parsed), ncol = 2) +
#   geom_pointrange(position = position_dodge(width = 0.07)) +
#   geom_hline(yintercept = 0) +
#   geom_vline(xintercept = 2/12, linetype = "dashed") +
#   geom_vline(xintercept = 9/12, linetype = "dotted") +
#   scale_x_continuous("Age in months", breaks = seq(0, 1.5, by = 1/12), labels = 0:18) +
#   scale_y_continuous("Average annual infections averted compared to status quo") +
#   scale_colour_manual("Strategy", values = c("#88CCEE", "#CC6677")) +
#   theme(legend.position = "bottom")

# tiff("Figure 2A.tiff", width = 10, height = 5.5, units = "in", res = 1000)
# fig_2a
# dev.off()

# png("Figure 2A.png", width = 10, height = 5.5, units = "in", res = 200)
# fig_2a
# dev.off()

fig_2 <- ggplot(dat_par_vis_rel[age <= 1.5 & alpha == "alpha == 1"], aes(x = age, y = med, ymin = low, ymax = upp, colour = variable)) +
  geom_pointrange(position = position_dodge(width = 0.07)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 2/12, linetype = "dashed") +
  geom_vline(xintercept = 9/12, linetype = "dotted") +
  scale_x_continuous("Age in months", breaks = seq(0, 1.5, by = 1/12), labels = 0:18) +
  scale_y_continuous("Average annual reduction in infections compared to status quo", labels = scales::label_percent(), n.breaks = 8) +
  scale_colour_manual("Strategy", values = c("#88CCEE", "#CC6677")) +
  theme(legend.position = "bottom")

tiff("Figure 2.tiff", width = 10, height = 5.5, units = "in", res = 1000)
fig_2
dev.off()

png("Figure 2.png", width = 10, height = 5.5, units = "in", res = 200)
fig_2
dev.off()

### Figure 3

dat_icer_vis <- rbindlist(lapply(sims_cea, function(cea_sim){
  cea_increm_full <- list("MV vs SQ"   = eval_icer(costs = cea_sim$cea$costs, qalys = cea_sim$cea$qalys, curr = "NS", new = "MV"),
                          "mAbs vs SQ" = eval_icer(costs = cea_sim$cea$costs, qalys = cea_sim$cea$qalys, curr = "NS", new = "mAbs"),
                          "mAbs vs MV" = eval_icer(costs = cea_sim$cea$costs, qalys = cea_sim$cea$qalys, curr = "MV", new = "mAbs"))
  cea_increm_rest <- list("MV vs SQ"   = eval_icer(costs = cea_sim$cea$costs, qalys = cea_sim$cea$qalys, curr = "NS", new = "MV",
                                                   max_age_ind = 60),
                          "mAbs vs SQ" = eval_icer(costs = cea_sim$cea$costs, qalys = cea_sim$cea$qalys, curr = "NS", new = "mAbs",
                                                   max_age_ind = 60),
                          "mAbs vs MV" = eval_icer(costs = cea_sim$cea$costs, qalys = cea_sim$cea$qalys, curr = "MV", new = "mAbs",
                                                   max_age_ind = 60))
  data.table(alpha = paste("alpha ==", cea_sim$cfg$alpha),
             rest  = factor(rep(c("All ages", "Age < 5 years"), each = N_draw*3), levels = c("All ages", "Age < 5 years")),
             comp  = factor(rep(rep(names(cea_increm_full), each = N_draw), 2), levels = names(cea_increm_full)),
             costs = c(as.vector(sapply(cea_increm_full, function(x) x$cost_increm/1000000)),
                       as.vector(sapply(cea_increm_rest, function(x) x$cost_increm/1000000))),
             qalys = c(as.vector(sapply(cea_increm_full, function(x) x$qaly_increm)),
                       as.vector(sapply(cea_increm_rest, function(x) x$qaly_increm))))
}))
dat_slope <- data.frame(slope = WTP/1000000,
                        WTP = factor(paste0("$", WTP/1000, "K"), levels = paste0("$", WTP/1000, "K")))

fig_3 <- ggplot(dat_icer_vis[rest == "Age < 5 years"], aes(x = qalys/max(years_sub), y = costs/max(years_sub), colour = alpha)) +
    facet_grid(~comp) +
    geom_point(alpha = 0.05) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    geom_abline(data = dat_slope, mapping = aes(intercept = 0, slope = slope, linetype = WTP)) +
    scale_x_continuous("Incremental quality-adjusted life years", limits = c(-120, 300)) +
    scale_y_continuous("Incremental cost (in million $AUD)", limits = c(-180, 150)) +
    scale_colour_manual("",
                        labels = function(x) parse(text = x),
                        values = c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499")) +
    scale_linetype_manual("Willingness-to-pay", values = c("dashed", "dotted")) +
    guides(colour = guide_legend(order = 1, override.aes = list(alpha = 1))) +
    theme(legend.position = "bottom")

tiff("Figure 3.tiff", width = 10, height = 5, units = "in", res = 1000)
fig_3
dev.off()

png("Figure 3.png", width = 10, height = 5, units = "in", res = 200)
fig_3
dev.off()

#### results to report in text

dat_icer_vis[rest == "Age < 5 years",
             .(qaly_summ = paste0(round(median(qalys/max(years_sub))),
                                  " (", round(quantile(qalys/max(years_sub), 0.025)), ",", round(quantile(qalys/max(years_sub), 0.975)), ")"),
               cost_summ = paste0(round(median(costs/max(years_sub))),
                                  " (", round(quantile(costs/max(years_sub), 0.025)), ",", round(quantile(costs/max(years_sub), 0.975)), ")")),
             by = comp]

dat_inmb_vis <- rbindlist(lapply(sims_cea, function(cea_sim)
  data.table(alpha    = paste("alpha ==", cea_sim$cfg$alpha),
             INMB     = unlist(cea_sim$INMB),
             sim      = rep(1:N_draw, 3*length(WTP)*length(years_sub)),
             strategy = rep(rep(c("NS", "MV", "mAbs"), each = N_draw), length(WTP)*length(years_sub)),
             WTP      = factor(rep(rep(paste0("WTP = $", WTP/1000, "K"), each = N_draw*3), length(years_sub))),
             years_f  = factor(rep(paste("years =", years_sub), each = N_draw*3*length(WTP)), levels = paste("years =", years_sub)),
             years_n  = rep(years_sub, each = N_draw*3*length(WTP)))))
dat_inmb_vis[, alpha := factor(alpha)]
dat_inmb_vis <- dcast.data.table(dat_inmb_vis, alpha + WTP + years_f + sim ~ strategy, value.var = "INMB")
dat_inmb_vis[, max_strategy := ifelse(NS > MV & NS > mAbs, "NS", ifelse(MV > NS & MV > mAbs, "MV", ifelse(mAbs > NS & mAbs > MV, "mAbs", NA)))]
dat_inmb_vis <- dat_inmb_vis[, .(NS = mean(max_strategy == "NS"), MV = mean(max_strategy == "MV"), mAbs = mean(max_strategy == "mAbs")),
                             by = .(alpha, WTP, years_f)]
dat_inmb_vis <- melt.data.table(dat_inmb_vis, id.vars = c("alpha", "WTP", "years_f"), measure.vars = c("NS", "MV", "mAbs"))

dat_inmb_vis[, .(paste0(format(round(median(value), 2), nsmall = 2), " (",
                        format(round(quantile(value, 0.025), 2), nsmall = 2), ",",
                        format(round(quantile(value, 0.975), 2), nsmall = 2), ")")),
             by = variable]

### Figure 4

dat_voi_p_vis <- rbindlist(lapply(sims_voi_p, function(voi_p_sim)
  data.table(alpha    = paste("alpha ==", voi_p_sim$cfg$alpha),
             WTP      = factor(rep(paste0("WTP = $", WTP/1000, "K"), (1 + 6)*length(years_sub))),
             model    = factor(c(rep("-", length(WTP)*length(years_sub)),
                                 rep(rep(c("strategy_effectiveness", "transmission", "strategy_coverage",
                                           "strategy_costs", "health_costs", "health_utilities"), each = length(WTP)), length(years_sub))),
                               levels = c("-", "strategy_effectiveness", "strategy_coverage", "strategy_costs",
                                          "transmission", "health_costs", "health_utilities"),
                               labels = c("-", "Strategy Effectiveness", "Strategy Coverage", "Strategy Cost",
                                          "Transmission Dynamics", "Healthcare Cost", "Healthcare Utility")),
             years_f  = factor(c(rep(paste("years =", years_sub), each = length(WTP)),
                                 rep(paste("years =", years_sub), each = length(WTP)*6)),
                               levels = paste("years =", years_sub)),
             years_n  = c(rep(years_sub, each = length(WTP)),
                          rep(years_sub, each = length(WTP)*6)),
             variable = rep(c("EVPI", "EVPPI"), times = length(WTP)*length(years_sub)*c(1, 6)),
             value    = c(unlist(voi_p_sim$EVPI), unlist(voi_p_sim$EVPPI))/1000000)))
dat_voi_p_vis[, alpha := factor(alpha)]

fig_4 <- ggplot(dat_voi_p_vis[variable == "EVPI"], aes(x = years_n, y = value, colour = alpha, shape = WTP, group = WTP)) +
  geom_point(position = position_dodge(width = 0.6)) +
  scale_x_continuous("Decision horizon (years)", breaks = years_sub) +
  ylab("Value (in million $AUD)") +
  scale_colour_manual("", labels = function(x) parse(text = x), values = c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499")) +
  scale_shape_discrete("") +
  guides(colour = guide_legend(order = 1), shape = guide_legend(order = 2)) +
  theme(legend.position = "bottom")

tiff("Figure 4.tiff", width = 6, height = 5, units = "in", res = 1000)
fig_4
dev.off()

png("Figure 4.png", width = 6, height = 5, units = "in", res = 200)
fig_4
dev.off()

#### results to report in text

dat_voi_p_vis[variable == "EVPI" & years_n == 15, .(min = min(value), max = max(value))]
dat_voi_p_vis[variable == "EVPPI" & years_n == 15, .(min = min(value), max = max(value)), by = model]

### Figure 5

fig_5 <- ggplot(dat_voi_p_vis[variable == "EVPPI" & model %in% c("Strategy Effectiveness", "Strategy Coverage", "Strategy Cost",
                                                                 "Transmission Dynamics", "Healthcare Cost", "Healthcare Utility")],
                aes(x = years_n, y = value, colour = alpha, shape = WTP, group = WTP)) +
  facet_wrap(~model) +
  geom_point(position = position_dodge(width = 0.6)) +
  scale_x_continuous("Decision horizon (years)", breaks = years_sub) +
  ylab("Value (in million $AUD)") +
  scale_colour_manual("", labels = function(x) parse(text = x), values = c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499")) +
  scale_shape_discrete("") +
  guides(colour = guide_legend(order = 1), shape = guide_legend(order = 2)) +
  theme(legend.position = "bottom")

tiff("Figure 5.tiff", width = 7, height = 5, units = "in", res = 1000)
fig_5
dev.off()

png("Figure 5.png", width = 7, height = 5, units = "in", res = 200)
fig_5
dev.off()

### Figure 6

dat_voi_s_vis <- rbindlist(lapply(sims_voi_s, function(voi_sim)
  data.table(alpha   = paste("alpha ==", voi_sim$cfg$alpha),
             design  = voi_sim$cfg$design,
             N_f     = factor(paste("N =", voi_sim$cfg$N), levels = c("N = 1000", "N = 5000", "N = 10000", "N = 20000", "N = 50000")),
             N_n     = voi_sim$cfg$N,
             WTP     = factor(rep(paste0("WTP = $", WTP/1000, "K"), length(years_sub))),
             years_f = factor(rep(paste("years =", years_sub), each = length(WTP)), levels = paste("years =", years_sub)),
             years_n = rep(years_sub, each = length(WTP)),
             EVSI = unlist(voi_sim$EVSI)/1000000)))
dat_voi_s_vis[, alpha := factor(alpha)]

fig_6 <- ggplot(dat_voi_s_vis[design == "C"], aes(x = years_n, y = EVSI, colour = alpha, shape = WTP, group = WTP)) +
  facet_wrap(~ N_f) +
  geom_point(position = position_dodge(width = 0.6)) +
  scale_x_continuous("Decision horizon (years)", breaks = years_sub) +
  ylab("Value (in million $AUD)") +
  scale_colour_manual("", labels = function(x) parse(text = x), values = c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499")) +
  scale_shape_discrete("") +
  guides(colour = guide_legend(order = 1), shape = guide_legend(order = 2)) +
  theme(legend.position = "bottom")

tiff("Figure 6.tiff", width = 7, height = 5, units = "in", res = 1000)
fig_6
dev.off()

png("Figure 6.png", width = 7, height = 5, units = "in", res = 200)
fig_6
dev.off()

#### results to report in text

dat_voi_s_vis[years_n == 15, .(min = min(EVSI), max = max(EVSI)), by = .(design, N_n)]

## Supplementary Information

### Figure S1

nodes <- data.frame(text = c("S[i]", "E[i]", "I[i]", "R[i]"),
                    x    = c(     0,      1,      2,      3),
                    y    = c(     0,      0,      0,      0))
edges_straight <- data.frame(x    = c(  0,   1,   2),
                             y    = c(  0,   0,   0),
                             xend = c(0.9, 1.9, 2.9),
                             yend = c(  0,   0,   0))
edges_curved <- data.frame(x    = c(   3),
                           y    = c(   0),
                           xend = c(0.08),
                           yend = c(   0))

fig_s1 <- ggplot() +
  geom_segment(data = edges_straight, aes(x = x, y = y, xend = xend, yend = yend), arrow = arrow(length = unit(0.2, "cm"))) +
  geom_curve(data = edges_curved, aes(x = x, y = y, xend = xend, yend = yend), arrow = arrow(length = unit(0.2, "cm")), curvature = 0.15) +
  geom_point(data = nodes, aes(x = x, y = y), shape = 21, fill = "white", size = 10) +
  geom_text(data = nodes, aes(x = x, y = y, label = text), size = 4, parse = TRUE) +
  annotate("text", x = 0.6, y = 0, label = expression(lambda[i]*sigma[i]), size = 4, vjust = 1.5) +
  annotate("text", x = 1.5, y = 0, label = expression(delta), size = 4, vjust = 1.5) +
  annotate("text", x = 2.5, y = 0, label = expression(gamma), size = 4, vjust = 1.5) +
  annotate("text", x = 1.5, y = 0, label = expression(nu), size = 4, vjust = -3) +
  theme_void()

tiff("Figure S1.tiff", width = 6, height = 3, units = "in", res = 1000)
fig_s1
dev.off()

png("Figure S1.png", width = 6, height = 3, units = "in", res = 200)
fig_s1
dev.off()

### Figure S2

nodes <- data.frame(text = c(     "V[1]",      "S[1]",      "E[1]",      "I[1]",      "R[1]",
                                  "V[2]",      "S[2]",      "E[2]",      "I[2]",      "R[2]",
                             "V[dur[V]]", "S[dur[V]]", "E[dur[V]]", "I[dur[V]]", "R[dur[V]]",
                                  "S[i]",      "E[i]",      "I[i]",      "R[i]"),
                    x    = c(          0,           1,           2,           3,           4,
                                       0,           1,           2,           3,           4,
                                       0,           1,           2,           3,           4,
                                       1,           2,           3,           4),
                    y    = c(          2,           2,           2,           2,           2,
                                       1,           1,           1,           1,           1,
                                      -1,          -1,          -1,          -1,          -1,
                                      -2,          -2,          -2,          -2))
edges_straight <- data.frame(x      = c(  1,   2,   3,
                                          1,   2,   3,
                                          1,   2,   3,
                                          1,   2,   3),
                             y      = c(  2,   2,   2,
                                          1,   1,   1,
                                         -1,  -1,  -1,
                                         -2,  -2,  -2),
                             xend   = c(1.9, 2.9, 3.9,
                                        1.9, 2.9, 3.9,
                                        1.9, 2.9, 3.9,
                                        1.9, 2.9, 3.9),
                             yend   = c(  2,   2,   2,
                                          1,   1,   1,
                                         -1,  -1,  -1,
                                         -2,  -2,  -2))
edges_curved_one <- data.frame(x    = c(  4,   4,    4,    4),
                               y    = c(  2,   1,   -1,   -2),
                               xend = c(1.1, 1.1,  1.1,  1.1),
                               yend = c(2.1, 1.1, -0.9, -1.9))
edges_curved_two <- data.frame(x    = c(  0,   0,    0),
                               y    = c(  2,   1,   -1),
                               xend = c(1.9, 1.9,  1.9),
                               yend = c(2.1, 1.1, -0.9))
edges_dashed <- data.frame(x    = c(   0,    1,
                                       0,    1,    2,    3,    4,
                                       0,    1,    2,    3,    4,
                                       0,    1,    2,    3,    4,
                                       0,    1,    2,    3,    4,
                                       1,    2,    3,    4),
                           y    = c( 2.7,  2.7,
                                     1.8,  1.8,  1.8,  1.8,  1.8,
                                     0.8,  0.8,  0.8,  0.8,  0.8,
                                    -0.2, -0.2, -0.2, -0.2, -0.2,
                                    -1.2, -1.2, -1.2, -1.2, -1.2,
                                    -2.2, -2.2, -2.2, -2.2),
                           xend = c(   0,    1,
                                       0,    1,    2,    3,    4,
                                       0,    1,    2,    3,    4,
                                       0,    1,    2,    3,    4,
                                     0.9,    1,    2,    3,    4,
                                       1,    2,    3,    4),
                           yend = c( 2.2,  2.2,
                                     1.2,  1.2,  1.2,  1.2,  1.2,
                                     0.2,  0.2,  0.2,  0.2,  0.2,
                                    -0.8, -0.8, -0.8, -0.8, -0.8,
                                      -2, -1.8, -1.8, -1.8, -1.8,
                                    -2.7, -2.7, -2.7, -2.7))
edges_angle <- data.frame(x    = c( -0.1,   0.9,   1.9,   2.9,   3.9,
                                    -0.1,   0.9,   1.9,   2.9,   3.9),
                          y    = c(-0.15, -0.15, -0.15, -0.15, -0.15,
                                   -0.05, -0.05, -0.05, -0.05, -0.05),
                          xend = c(  0.1,   1.1,   2.1,   3.1,   4.1,
                                     0.1,   1.1,   2.1,   3.1,   4.1),
                          yend = c( 0.05,  0.05,  0.05,  0.05,  0.05,
                                    0.15,  0.15,  0.15,  0.15,  0.15))
edge_text <- data.frame(text = c(                                 "kappa[V]",                 "1 - kappa[V]",
                                           "lambda[1]*sigma[1]*(1 - rho[V])",           "lambda[1]*sigma[1]", "delta", "gamma",  "nu",
                                           "lambda[2]*sigma[2]*(1 - rho[V])",           "lambda[2]*sigma[2]", "delta", "gamma",  "nu",
                                 "lambda[dur[V]]*sigma[dur[V]]*(1 - rho[V])", "lambda[dur[V]]*sigma[dur[V]]", "delta", "gamma",  "nu",
                                                                                        "lambda[i]*sigma[i]", "delta", "gamma",  "nu"),
                        x    = c(                                       0.18,                           1.18,
                                                                        0.55,                            1.6,     2.5,     3.4,   2.5,
                                                                        0.55,                            1.6,     2.5,     3.4,   2.5,
                                                                        0.55,                            1.6,     2.5,     3.4,   2.5,
                                                                                                         1.6,     2.5,     3.4,   2.5),
                        y    = c(                                        2.5,                            2.6,
                                                                         1.9,                           1.82,    1.82,    1.82,  2.35,
                                                                         0.9,                           0.82,    0.82,    0.82,  1.35,
                                                                        -1.1,                          -1.17,   -1.17,   -1.17, -0.65,
                                                                                                       -2.17,   -2.17,   -2.17, -1.65))

fig_s2 <- ggplot() +
  geom_segment(data = edges_straight, aes(x = x, y = y, xend = xend, yend = yend), arrow = arrow(length = unit(0.2, "cm"))) +
  geom_curve(data = edges_curved_one, aes(x = x, y = y, xend = xend, yend = yend), arrow = arrow(length = unit(0.2, "cm")), curvature = 0.15) +
  geom_curve(data = edges_curved_two, aes(x = x, y = y, xend = xend, yend = yend), arrow = arrow(length = unit(0.2, "cm")), curvature = -0.15) +
  geom_segment(data = edges_dashed, aes(x = x, y = y, xend = xend, yend = yend), arrow = arrow(length = unit(0.2, "cm")), linetype = "dashed") +    geom_segment(data = edges_angle, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_point(data = nodes, aes(x = x, y = y), shape = 21, fill = "white", size = 10) +
  geom_text(data = nodes, aes(x = x, y = y, label = text), size = 3, parse = TRUE) +
  geom_text(data = edge_text, aes(x = x, y = y, label = text), size = 3, parse = TRUE) +
  coord_cartesian(ylim = c(-3, 3)) +
  theme_void()

tiff("Figure S2.tiff", width = 6, height = 7, units = "in", res = 1000)
fig_s2
dev.off()

png("Figure S2.png", width = 6, height = 7, units = "in", res = 200)
fig_s2
dev.off()

### Figure S3

nodes <- data.frame(text = c(     "M[1]",      "S[1]",      "E[1]",      "I[1]",      "R[1]",
                                  "M[2]",      "S[2]",      "E[2]",      "I[2]",      "R[2]",
                             "M[dur[M]]", "S[dur[M]]", "E[dur[M]]", "I[dur[M]]", "R[dur[M]]",
                                  "S[i]",      "E[i]",      "I[i]",      "R[i]"),
                    x    = c(          0,           1,           2,           3,           4,
                                       0,           1,           2,           3,           4,
                                       0,           1,           2,           3,           4,
                                       1,           2,           3,           4),
                    y    = c(          2,           2,           2,           2,           2,
                                       1,           1,           1,           1,           1,
                                      -1,          -1,          -1,          -1,          -1,
                                      -2,          -2,          -2,          -2))
edges_straight <- data.frame(x      = c(  1,   2,   3,
                                          1,   2,   3,
                                          1,   2,   3,
                                          1,   2,   3),
                             y      = c(  2,   2,   2,
                                          1,   1,   1,
                                         -1,  -1,  -1,
                                         -2,  -2,  -2),
                             xend   = c(1.9, 2.9, 3.9,
                                        1.9, 2.9, 3.9,
                                        1.9, 2.9, 3.9,
                                        1.9, 2.9, 3.9),
                             yend   = c(  2,   2,   2,
                                          1,   1,   1,
                                         -1,  -1,  -1,
                                         -2,  -2,  -2))
edges_curved_one <- data.frame(x    = c(  4,   4,    4,    4),
                               y    = c(  2,   1,   -1,   -2),
                               xend = c(1.1, 1.1,  1.1,  1.1),
                               yend = c(2.1, 1.1, -0.9, -1.9))
edges_curved_two <- data.frame(x    = c(  0,   0,    0),
                               y    = c(  2,   1,   -1),
                               xend = c(1.9, 1.9,  1.9),
                               yend = c(2.1, 1.1, -0.9))
edges_dashed <- data.frame(x    = c(   0,    1,
                                       0,    1,    2,    3,    4,
                                       0,    1,    2,    3,    4,
                                       0,    1,    2,    3,    4,
                                       0,    1,    2,    3,    4,
                                       1,    2,    3,    4),
                           y    = c( 2.7,  2.7,
                                     1.8,  1.8,  1.8,  1.8,  1.8,
                                     0.8,  0.8,  0.8,  0.8,  0.8,
                                    -0.2, -0.2, -0.2, -0.2, -0.2,
                                    -1.2, -1.2, -1.2, -1.2, -1.2,
                                    -2.2, -2.2, -2.2, -2.2),
                           xend = c(   0,    1,
                                       0,    1,    2,    3,    4,
                                       0,    1,    2,    3,    4,
                                       0,    1,    2,    3,    4,
                                     0.9,    1,    2,    3,    4,
                                       1,    2,    3,    4),
                           yend = c( 2.2,  2.2,
                                     1.2,  1.2,  1.2,  1.2,  1.2,
                                     0.2,  0.2,  0.2,  0.2,  0.2,
                                    -0.8, -0.8, -0.8, -0.8, -0.8,
                                      -2, -1.8, -1.8, -1.8, -1.8,
                                    -2.7, -2.7, -2.7, -2.7))
edges_angle <- data.frame(x    = c( -0.1,   0.9,   1.9,   2.9,   3.9,
                                    -0.1,   0.9,   1.9,   2.9,   3.9),
                          y    = c(-0.15, -0.15, -0.15, -0.15, -0.15,
                                   -0.05, -0.05, -0.05, -0.05, -0.05),
                          xend = c(  0.1,   1.1,   2.1,   3.1,   4.1,
                                     0.1,   1.1,   2.1,   3.1,   4.1),
                          yend = c( 0.05,  0.05,  0.05,  0.05,  0.05,
                                    0.15,  0.15,  0.15,  0.15,  0.15))
edge_text <- data.frame(text = c(                                 "kappa[M]",                 "1 - kappa[M]",
                                           "lambda[1]*sigma[1]*(1 - rho[M])",           "lambda[1]*sigma[1]", "delta", "gamma",  "nu",
                                           "lambda[2]*sigma[2]*(1 - rho[M])",           "lambda[2]*sigma[2]", "delta", "gamma",  "nu",
                                 "lambda[dur[M]]*sigma[dur[M]]*(1 - rho[M])", "lambda[dur[M]]*sigma[dur[M]]", "delta", "gamma",  "nu",
                                                                                        "lambda[i]*sigma[i]", "delta", "gamma",  "nu"),
                        x    = c(                                       0.18,                           1.18,
                                                                        0.55,                            1.6,     2.5,     3.4,   2.5,
                                                                        0.55,                            1.6,     2.5,     3.4,   2.5,
                                                                        0.55,                            1.6,     2.5,     3.4,   2.5,
                                                                                                         1.6,     2.5,     3.4,   2.5),
                        y    = c(                                        2.5,                            2.6,
                                                                         1.9,                           1.82,    1.82,    1.82,  2.35,
                                                                         0.9,                           0.82,    0.82,    0.82,  1.35,
                                                                        -1.1,                          -1.17,   -1.17,   -1.17, -0.65,
                                                                                                       -2.17,   -2.17,   -2.17, -1.65))

fig_s3 <- ggplot() +
  geom_segment(data = edges_straight, aes(x = x, y = y, xend = xend, yend = yend), arrow = arrow(length = unit(0.2, "cm"))) +
  geom_curve(data = edges_curved_one, aes(x = x, y = y, xend = xend, yend = yend), arrow = arrow(length = unit(0.2, "cm")), curvature = 0.15) +
  geom_curve(data = edges_curved_two, aes(x = x, y = y, xend = xend, yend = yend), arrow = arrow(length = unit(0.2, "cm")), curvature = -0.15) +
  geom_segment(data = edges_dashed, aes(x = x, y = y, xend = xend, yend = yend), arrow = arrow(length = unit(0.2, "cm")), linetype = "dashed") +    geom_segment(data = edges_angle, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_point(data = nodes, aes(x = x, y = y), shape = 21, fill = "white", size = 10) +
  geom_text(data = nodes, aes(x = x, y = y, label = text), size = 3, parse = TRUE) +
  geom_text(data = edge_text, aes(x = x, y = y, label = text), size = 3, parse = TRUE) +
  coord_cartesian(ylim = c(-3, 3)) +
  theme_void()

tiff("Figure S3.tiff", width = 6, height = 7, units = "in", res = 1000)
fig_s3
dev.off()

png("Figure S3.png", width = 6, height = 7, units = "in", res = 200)
fig_s3
dev.off()

### Figure S4

hyper_parms <- readRDS("Parameters/hyper_parms.rds")
alpha <- c(0.005, 0.01, 0.1, 0.25, 0.5, 1)
N_draw <- 10000
seed <- 35234
priors <- lapply(alpha, function(alph)
  list(gen_strat_eff_parms(N_draw = N_draw, hyper_parms = hyper_parms[["MA"]], alpha = alph, seed = seed),
       gen_strat_eff_parms(N_draw = N_draw, hyper_parms = hyper_parms[["HP"]], alpha = alph, seed = seed)))
dat_vis <- data.frame(x = unlist(priors),
                      alpha = rep(paste("alpha ==", alpha), each = N_draw*7*2),
                      Outcome = factor(rep(rep(c("Medical Attendance", "Hospitalisation"), each = N_draw*7), length(priors)),
                                       levels = c("Medical Attendance", "Hospitalisation")),
                      Strategy = factor(rep(rep(c("SQ", "MV", "mAbs", "MV", "mAbs", "MV", "mAbs"), each = N_draw), 2*length(priors)),
                                        levels = c("SQ", "MV", "mAbs")),
                      Parameter = rep(rep(c("prob", "prob", "prob", "lOR", "lOR", "RR", "RR"), each = N_draw), 2*length(priors)))

fig_s4 <- ggplot(dat_vis[dat_vis$Parameter == "prob",], aes(x = x, colour = Strategy)) +
  facet_grid(alpha ~ Outcome, labeller = labeller(alpha = label_parsed, Outcome = label_value)) +
  geom_density() +
  scale_colour_manual(values = c("black", "#88CCEE", "#CC6677")) +
  scale_x_continuous(expression(p[d]^q), limits = c(0, 0.08), n.breaks = 5) +
  ylab("Density") +
  theme(legend.position = "bottom")

tiff("Figure S4.tiff", width = 6, height = 6, units = "in", res = 1000)
fig_s4
dev.off()

png("Figure S4.png", width = 6, height = 6, units = "in", res = 200)
fig_s4
dev.off()

### Figure S5

fig_s5 <- ggplot(dat_vis[dat_vis$Parameter == "RR",], aes(x = x, colour = Strategy)) +
  facet_grid(alpha ~ Outcome, labeller = labeller(alpha = label_parsed, Outcome = label_value)) +
  geom_density() +
  scale_colour_manual(values = c("#88CCEE", "#CC6677")) +
  scale_x_continuous(expression(r[d]^q), n.breaks = 5, limits = c(0, 1)) +
  ylab("Density") +
  theme(legend.position = "bottom")

tiff("Figure S5.tiff", width = 6, height = 6, units = "in", res = 1000)
fig_s5
dev.off()

png("Figure S5.png", width = 6, height = 6, units = "in", res = 200)
fig_s5
dev.off()

## Figure S6

fig_s6 <- ggplot(dat_vis[dat_vis$Parameter == "lOR",], aes(x = x, colour = Strategy)) +
  facet_grid(alpha ~ Outcome, labeller = labeller(alpha = label_parsed, Outcome = label_value)) +
  geom_density() +
  scale_colour_manual(values = c("#88CCEE", "#CC6677")) +
  scale_x_continuous(expression(pi[d]^q), n.breaks = 5, limits = c(-5, 2)) +
  ylab("Density") +
  theme(legend.position = "bottom")

tiff("Figure S6.tiff", width = 6, height = 6, units = "in", res = 1000)
fig_s6
dev.off()

png("Figure S6.png", width = 6, height = 6, units = "in", res = 200)
fig_s6
dev.off()

## Figure S7

fig_s7 <- ggplot(dat_par_vis_rel[age <= 1.5], aes(x = age, y = med, ymin = low, ymax = upp, colour = variable)) +
  facet_wrap(~ alpha, labeller = labeller(alpha = label_parsed), ncol = 2) +
  geom_pointrange(position = position_dodge(width = 0.07)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 2/12, linetype = "dashed") +
  geom_vline(xintercept = 9/12, linetype = "dotted") +
  scale_x_continuous("Age in months", breaks = seq(0, 1.5, by = 1/12), labels = 0:18) +
  scale_y_continuous("Average annual reduction in infections compared to status quo", labels = scales::label_percent(), n.breaks = 8) +
  scale_colour_manual("Strategy", values = c("#88CCEE", "#CC6677")) +
  theme(legend.position = "bottom")

tiff("Figure S7.tiff", width = 10, height = 5.5, units = "in", res = 1000)
fig_s7
dev.off()

png("Figure S7.png", width = 10, height = 5.5, units = "in", res = 200)
fig_s7
dev.off()

## Figure S8

fig_s8 <- ggplot(dat_par_vis_rel[age > 1.5 & age < 5], aes(x = age, y = med, ymin = low, ymax = upp, colour = variable)) +
  facet_wrap(~ alpha, labeller = labeller(alpha = label_parsed), ncol = 2) +
  geom_pointrange(position = position_dodge(width = 0.07)) +
  geom_hline(yintercept = 0) +
  scale_x_continuous("Age in years", breaks = seq(1, 5, by = 0.5)) +
  scale_y_continuous("Average annual reduction in infections compared to status quo", labels = scales::label_percent(), n.breaks = 8) +
  scale_colour_manual("Strategy", values = c("#88CCEE", "#CC6677")) +
  theme(legend.position = "bottom")

tiff("Figure S8.tiff", width = 10, height = 5.5, units = "in", res = 1000)
fig_s8
dev.off()

png("Figure S8.png", width = 10, height = 5.5, units = "in", res = 200)
fig_s8
dev.off()

## Figure S9

fig_s9 <- ggplot(dat_par_vis_rel[age >= 5], aes(x = age, y = med, ymin = low, ymax = upp, colour = variable)) +
  facet_wrap(~ alpha, labeller = labeller(alpha = label_parsed), ncol = 2) +
  geom_pointrange(position = position_dodge(width = 2)) +
  geom_hline(yintercept = 0) +
  scale_x_continuous("Age in years", breaks = seq(5, 75, by = 5)) +
  scale_y_continuous("Average annual reduction in infections compared to status quo", labels = scales::label_percent(), n.breaks = 8) +
  scale_colour_manual("Strategy", values = c("#88CCEE", "#CC6677")) +
  theme(legend.position = "bottom")

tiff("Figure S9.tiff", width = 10, height = 5.5, units = "in", res = 1000)
fig_s9
dev.off()

png("Figure S9.png", width = 10, height = 5.5, units = "in", res = 200)
fig_s9
dev.off()

## Figure S10

fig_s10 <- ggplot(dat_voi_s_vis[design == "A"], aes(x = years_n, y = EVSI, colour = alpha, shape = WTP, group = WTP)) +
  facet_wrap(~ N_f) +
  geom_point(position = position_dodge(width = 0.6)) +
  scale_x_continuous("Decision horizon (years)", breaks = years_sub) +
  ylab("Value (in million $AUD)") +
  scale_colour_manual("", labels = function(x) parse(text = x), values = c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499")) +
  scale_shape_discrete("") +
  guides(colour = guide_legend(order = 1), shape = guide_legend(order = 2)) +
  theme(legend.position = "bottom")

tiff("Figure S10.tiff", width = 7, height = 5, units = "in", res = 1000)
fig_s10
dev.off()

png("Figure S10.png", width = 7, height = 5, units = "in", res = 200)
fig_s10
dev.off()

## Figure S11

fig_s11 <- ggplot(dat_voi_s_vis[design == "B"], aes(x = years_n, y = EVSI, colour = alpha, shape = WTP, group = WTP)) +
  facet_wrap(~ N_f) +
  geom_point(position = position_dodge(width = 0.6)) +
  scale_x_continuous("Decision horizon (years)", breaks = years_sub) +
  ylab("Value (in million $AUD)") +
  scale_colour_manual("", labels = function(x) parse(text = x), values = c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499")) +
  scale_shape_discrete("") +
  guides(colour = guide_legend(order = 1), shape = guide_legend(order = 2)) +
  theme(legend.position = "bottom")

tiff("Figure S11.tiff", width = 7, height = 5, units = "in", res = 1000)
fig_s11
dev.off()

png("Figure S11.png", width = 7, height = 5, units = "in", res = 200)
fig_s11
dev.off()

## Figure S12

fig_s12 <- ggplot(dat_voi_s_vis[design == "D"], aes(x = years_n, y = EVSI, colour = alpha, shape = WTP, group = WTP)) +
  facet_wrap(~ N_f) +
  geom_point(position = position_dodge(width = 0.6)) +
  scale_x_continuous("Decision horizon (years)", breaks = years_sub) +
  ylab("Value (in million $AUD)") +
  scale_colour_manual("", labels = function(x) parse(text = x), values = c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499")) +
  scale_shape_discrete("") +
  guides(colour = guide_legend(order = 1), shape = guide_legend(order = 2)) +
  theme(legend.position = "bottom")

tiff("Figure S12.tiff", width = 7, height = 5, units = "in", res = 1000)
fig_s12
dev.off()

png("Figure S12.png", width = 7, height = 5, units = "in", res = 200)
fig_s12
dev.off()
