############################################
## This generates plots for the supplementary materials of the simulation paper.
## This is for simulation 10.10 using the bryophyte template (large community)
##
## This script is meant to be run after sourcing analyze_performance.R
##
## inputs:  * bryophyte_template_evals_13Sep2019.csv 
##          * bryophyte_template_butterfly_bias_evals_12Sep2019.csv
##
## 
## author: Willson Gaul wgaul@hotmail.com
## created: 10 Sep 2019 based on code from plots.R
## last modified: 1 Nov 2019
############################################

library(wgutil)
#library(skimr)
library(Hmisc)
library(rgdal)
library(raster)
library(parallel)
library(pROC)
library(akima)
library(tidyverse)
library(virtualspecies)
library(simulator) # this file was created under simulator version 0.2.0

setwd("~/Documents/Data_Analysis/UCD/simulation/sims_10.10_bryophyte/")

t_size <- 19
line_size <- 1.5

sim_name <- "Bryophyte Simulation"

bias_names <- c("no_bias", "extreme_mothNBDC", "moderate_bryBBS",
                "butterflyNBDC", "least_birdNBDC")

## optionally, read in a .csv with auc results from simulation
## large community simulation results
ev_df <- read_csv("./bryophyte_template_evals_24Sep2019.csv")
df_butt <- read_csv("../sims_10.10_bryophyte_butterflyBias/bryophyte_template_butterfly_bias_evals_12Sep2019.csv")
ev_df <- bind_rows(ev_df, df_butt)

# subset to only block CV results
ev_df_cv <- ev_df[which(ev_df$Method %in% c("glm_poly_wg_block_cv", 
                                            "random_forest_wg_block_cv", 
                                            "brt_wg_block_cv", 
                                            "idw_interp_wg_block_cv")), ]

# subset to non-block CV results
ev_df_noCV <- ev_df[which(ev_df$Method %nin% c("glm_poly_wg_block_cv", 
                                               "random_forest_wg_block_cv", 
                                               "brt_wg_block_cv", 
                                               "idw_interp_wg_block_cv", 
                                               "occ_det_wg_block_cv")), ]

## small community simulation results
od_ev_df <- read_csv("../sims_10.10_odonata/odonata_template_evals.csv")
# od_butBias_df <- read_csv("../sims_10.10_odonata_butterflyBias/odonata_template_butterfly_bias_evals_25Aug2019.csv")

# remove methods that did not fully run (see simulation_tracking_status.ods)
od_ev_df <- od_ev_df[od_ev_df$Method %in% c("glm_poly", "glm_poly_wg_block_cv", 
                                            "idw_interp", 
                                            "idw_interp_wg_block_cv"), ]
# get only CV results
od_ev_df_cv <- od_ev_df[od_ev_df$Method %in% c("glm_poly_wg_block_cv",
                                               "idw_interp_wg_block_cv"), ]


## load predictor variables
load("~/Documents/Data_Analysis/UCD/predictor_variables/ETOPO1/elevation_hec_ETOPO1.RData")
load("~/Documents/Data_Analysis/UCD/predictor_variables/eobs/annual_precip_hectad.RData")
load("~/Documents/Data_Analysis/UCD/predictor_variables/eobs/summer_tx_hectad.RData")
load("~/Documents/Data_Analysis/UCD/predictor_variables/eobs/winter_tn_hectad.RData")
load("~/Documents/Data_Analysis/UCD/predictor_variables/eobs/mean_pp_hectad.RData")
load("~/Documents/Data_Analysis/UCD/predictor_variables/CORINE/corine_label_1_hectad.RData")
# Load Ireland coastline
ir <- readOGR(dsn='../../mapping/data/', layer='ireland_coastline')
ir_TM75 <- spTransform(ir, CRS("+init=epsg:29903"))
rm(ir)

pred_rast_brick <- brick(list(
  "minimum temperature" = resample(krg_mean_tn_rast, krg_mean_rr_rast), 
  "maximum temperature" = resample(krg_mean_tx_rast, krg_mean_rr_rast), 
  "annual precipitation" = krg_mean_rr_rast, 
  "atmospheric pressure" = resample(krg_mean_pp_rast, krg_mean_rr_rast), 
  "agricultural areas" = resample(agricultural_l1_rast, krg_mean_rr_rast), 
  "artificial surfaces" = resample(artificial_surfaces_l1_rast, krg_mean_rr_rast), 
  "forest/semi-natural" = resample(forest_seminatural_l1_rast, krg_mean_rr_rast),
  "wetlands" = resample(wetlands_l1_rast, krg_mean_rr_rast), 
  "water bodies" = resample(water_l1_rast, krg_mean_rr_rast), 
  "elevation" = resample(elev_hec, krg_mean_rr_rast)))

rm(krg_mean_pp, krg_mean_pp_rast, krg_mean_rr_predict, krg_mean_rr_rast, 
   krg_mean_tn_predict, krg_mean_tn_rast, krg_mean_tx_predict, krg_mean_tx_rast, 
   agricultural_l1_rast, artificial_surfaces_l1_rast, forest_seminatural_l1_rast, 
   wetlands_l1_rast, water_l1_rast, elev_hec)

## map predictor variables --------------------------------------------------
plot(mask(pred_rast_brick, ir_TM75))


##############################################################################
### block cross-validation GAM plots
## GAM smooth of AUC by bias
# use predictions from the four fitted GAM models (one model for each SDM method)

# make data for rug plot (with factor levels matching those in preds_auc_brt)
rug_df <- data.frame(method = ev_df_cv$Method, 
                     n_obs_per_sp = ev_df_cv$n_obs_per_sp, 
                     bias.name = ev_df_cv$bias.name, 
                     pred_auc = as.numeric(NA), 
                     pred_rmse = as.numeric(NA))
rug_df$bias.name <- gsub(".*_|NBDC|BBS", "", rug_df$bias.name)
rug_df$bias.name <- gsub("bias", "no bias", rug_df$bias.name)
rug_df$method <- gsub("_.*", "", rug_df$method)
rug_df$method <- gsub("idw", "idwi", rug_df$method)
rug_df$method <- gsub("random", "rf", rug_df$method)

auc_smooth <- ggplot(data = preds_auc_gam, 
                   aes(x = n_obs_per_sp, y = pred_auc, 
                       group = factor(bias.name,
                                      levels = c("no bias", "bird", "bry", 
                                                 "butterfly", "moth"),
                                      labels = c(
                                        "\nnone\n", "\nlow\n(Birds)", 
                                        "\nmoderate\n(Bryophytes)\n",
                                        "\nmedian\n(Butterflies)\n", 
                                        "\nsevere\n(Moths)\n"),
                                      ordered = TRUE))) + 
  geom_line(aes(color = factor(bias.name,
                               levels = c("no bias", "bird", "bry", 
                                          "butterfly", "moth"),
                               labels = c(
                                 "\nnone\n", "\nlow\n(Birds)", 
                                 "\nmoderate\n(Bryophytes)\n",
                                 "\nmedian\n(Butterflies)\n", 
                                 "\nsevere\n(Moths)\n"),
                               ordered = TRUE), 
                linetype = factor(bias.name,
                                  levels = c("no bias", "bird", "bry", 
                                             "butterfly", "moth"),
                                  labels = c(
                                    "\nnone\n", "\nlow\n(Birds)", 
                                    "\nmoderate\n(Bryophytes)\n",
                                    "\nmedian\n(Butterflies)\n", 
                                    "\nsevere\n(Moths)\n"),
                                  ordered = TRUE)), size = line_size) + 
  # ggtitle("Large Community Simulation") + #sim_name
  geom_rug(data = rug_df, aes(x = n_obs_per_sp), show.legend = F, sides = "b") +
  facet_wrap( ~ factor(method, 
                       levels = c("glm", "rf", "brt", "idwi"), 
                       labels = c("GLM", "Random\nForest", 
                                  "Boosted\nRegression\nTree", 
                                  "Inverse\ndistance-weighted\ninterpolation"))) + 
  xlab("Average number of observations \nper species in data set") + 
  ylab("AUC\n(block cross-validated)") +
  scale_linetype_discrete(name = "Spatial\nSampling\nBias") + 
  scale_color_viridis_d(name = "Spatial\nSampling\nBias", option = "viridis") +
  theme_bw() + 
  theme(text = element_text(size = t_size), 
        legend.key.width = unit(2*t_size, "points"))
auc_smooth

## GAM smooth of RMSE by bias
# use predictions from the four fitted GAM models (one model for each SDM method)
rmse_smooth <- ggplot(data = preds_rmse_gam, 
                    aes(x = n_obs_per_sp, y = pred_rmse, 
                        group = factor(bias.name,
                                       levels = c("no bias", "bird", "bry", 
                                                  "butterfly", "moth"),
                                       labels = c(
                                         "\nnone\n", "\nlow\n(Birds)", 
                                         "\nmoderate\n(Bryophytes)\n",
                                         "\nmedian\n(Butterflies)\n", 
                                         "\nsevere\n(Moths)\n"),
                                       ordered = TRUE))) + 
  geom_line(aes(color = factor(bias.name,
                               levels = c("no bias", "bird", "bry", 
                                          "butterfly", "moth"),
                               labels = c(
                                 "\nnone\n", "\nlow\n(Birds)", 
                                 "\nmoderate\n(Bryophytes)\n",
                                 "\nmedian\n(Butterflies)\n", 
                                 "\nsevere\n(Moths)\n"),
                               ordered = TRUE), 
                linetype = factor(bias.name,
                                  levels = c("no bias", "bird", "bry", 
                                             "butterfly", "moth"),
                                  labels = c(
                                    "\nnone\n", "\nlow\n(Birds)", 
                                    "\nmoderate\n(Bryophytes)\n",
                                    "\nmedian\n(Butterflies)\n", 
                                    "\nsevere\n(Moths)\n"),
                                  ordered = TRUE)), size = line_size) + 
  geom_rug(data = rug_df, aes(x = n_obs_per_sp), show.legend = F, sides = "b") +
  # ggtitle("Large Community Simulation") + #sim_name
  facet_wrap( ~ factor(method, 
                       levels = c("glm", "rf", "brt", "idwi"), 
                       labels = c("GLM", "Random\nForest", 
                                  "Boosted\nRegression\nTree", 
                                  "Inverse\ndistance-weighted\ninterpolation"))) + 
  xlab("Average number of observations \nper species in data set") + 
  ylab("RMSE\n(block cross-validated)") +
  scale_linetype_discrete(name = "Spatial\nSampling\nBias") + 
  scale_color_viridis_d(name = "Spatial\nSampling\nBias", option = "viridis") +
  theme_bw() + 
  theme(text = element_text(size = t_size), 
        legend.key.width = unit(2*t_size, "points"))
rmse_smooth


############################################################################
### no cross-validation plots
### plots for large community simulation --------------------------------------
## AUC plots 
## plot AUC by avg. num of observations per species 
## use fitted GAM models

# make data for rug plot (with factor levels matching those in preds_auc)
rug_df <- data.frame(method = ev_df_noCV$Method, 
                     n_obs_per_sp = ev_df_noCV$n_obs_per_sp, 
                     bias.name = ev_df_noCV$bias.name, 
                     pred_auc = as.numeric(NA), 
                     pred_rmse = as.numeric(NA))
rug_df$bias.name <- gsub(".*_|NBDC|BBS", "", rug_df$bias.name)
rug_df$bias.name <- gsub("bias", "no bias", rug_df$bias.name)
rug_df$method <- gsub("_.*", "", rug_df$method)
rug_df$method <- gsub("idw", "idwi", rug_df$method)
rug_df$method <- gsub("random", "rf", rug_df$method)
rug_df <- rug_df[rug_df$bias.name != "butterfly" & 
                   rug_df$method %in% c("glm", "idwi", "rf", "brt"), ] # drop butterfly bias


auc_smooth_noCV <- ggplot(data = preds_auc_gam_noCV, 
                         aes(x = n_obs_per_sp, y = pred_auc, 
                             group = factor(bias.name,
                                            levels = c("no bias", "bird", "bry", 
                                                       "moth"),
                                            labels = c(
                                              "\nnone\n",
                                              "\nlow\n(Birds)", 
                                              "\nmoderate\n(Bryophytes)\n", 
                                              "\nsevere\n(Moths)\n"),
                                            ordered = TRUE))) +
  geom_line(aes(color = factor(bias.name,
                               levels = c("no bias", "bird", "bry", "moth"),
                               labels = c(
                                 "\nnone\n", "\nlow\n(Birds)", 
                                 "\nmoderate\n(Bryophytes)\n",
                                 "\nsevere\n(Moths)\n"),
                               ordered = TRUE), 
                linetype = factor(bias.name,
                                  levels = c("no bias", "bird", "bry", "moth"),
                                  labels = c(
                                    "\nnone\n", "\nlow\n(Birds)", 
                                    "\nmoderate\n(Bryophytes)\n",
                                    "\nsevere\n(Moths)\n"),
                                  ordered = TRUE)), size = line_size) + 
  geom_rug(data = rug_df, aes(x = n_obs_per_sp), show.legend = F, sides = "b") +
  facet_wrap( ~ factor(method, 
                       levels = c("glm", "rf", "idwi", "brt"), 
                       labels = c("GLM", "Random\nForest", 
                                  "Boosted\nRegression\nTree", 
                                  "Inverse\ndistance-weighted\ninterpolation"))) + 
  xlab("Average number of observations \nper species in data set") + 
  ylab("AUC\n(no cross-validation)") +
  ylim(c(0.5, 1)) + 
  scale_linetype_discrete(name = "Spatial\nSampling\nBias") + 
  scale_color_viridis_d(name = "Spatial\nSampling\nBias", option = "viridis") +
  theme_bw() + 
  theme(text = element_text(size = t_size), 
        legend.key.width = unit(2*t_size, "points"))
auc_smooth_noCV

## plot AUC boxplots sample size and bias 
auc_boxplot_noCV <- ggplot(data = ev_df_noCV[ev_df_noCV$Method %in% 
                                               c("glm_poly", "random_forest", 
                                                 "brt", "idw_interp"), ], 
                           aes(x = factor(bias.name,
                                          levels = c("no_bias",
                                                     "least_birdNBDC",
                                                     "moderate_bryBBS",
                                                     "butterflyNBDC", 
                                                     "extreme_mothNBDC"),
                                          labels = c("none", "low", "moderate", 
                                                     "median", "severe"),
                                          ordered = TRUE), y = auc, 
                               fill = factor(round(n_obs_per_sp, 0), 
                                             levels = c("0", "1", "2", "4", "8", 
                                                        "20", "79", "158"), 
                                             labels = c("< 1", "1", "2", "4", "8", 
                                                        "20", "79", "158")))) +
  geom_boxplot(varwidth = TRUE) +
  facet_wrap(~ factor(Method, 
                      levels = c("glm_poly", 
                                 "random_forest", 
                                 "brt", 
                                 "idw_interp"), 
                      labels = c("GLM", "Random Forest", 
                                 "Boosted\nRegression\nTrees", 
                                 "Inverse\ndistance-weighted\ninterpolation"))) + 
  xlab("Sampling Bias") + 
  ylab("AUC\n(no cross-validation)") + 
  scale_fill_viridis_d(name = "Average\nnumber of\nrecords\nper species", 
                       option = "magma") + 
  theme_bw() + 
  theme(text = element_text(size = t_size), 
        axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1)) 
auc_boxplot_noCV

# RMSE boxplots block CV results
rmse_cv_boxplot <- ggplot(data = ev_df_cv, 
                          aes(x = factor(bias.name,
                                         levels = c("no_bias",
                                                    "least_birdNBDC",
                                                    "moderate_bryBBS",
                                                    "butterflyNBDC", 
                                                    "extreme_mothNBDC"),
                                         labels = c("none",
                                                    "low",
                                                    "moderate", 
                                                    "median", 
                                                    "severe"),
                                         ordered = TRUE), 
                              y = rmse, 
                              fill = factor(round(n_obs_per_sp, 0), 
                                            levels = c("0", "1", "2", "4", "8", 
                                                       "20", "79", "158"), 
                                            labels = c("< 1", "1", "2", "4", "8", 
                                                       "20", "79", "158")))) +
  geom_boxplot(varwidth = TRUE) +
  facet_wrap(~ factor(Method, 
                      levels = c("glm_poly_wg_block_cv", 
                                 "random_forest_wg_block_cv", 
                                 "brt_wg_block_cv", 
                                 "idw_interp_wg_block_cv"), 
                      labels = c("GLM", "Random\nForest", 
                                 "Boosted\nRegression\nTree", 
                                 "Inverse\ndistance-weighted\ninterpolation"))) + 
  # ggtitle("Large Community Simulation") + #sim_name
  xlab("Sampling Bias") + 
  ylab("RMSE\n(spatial block cross-validated)") + 
  scale_fill_viridis_d(name = "Average\nnumber of\nrecords\nper species", 
                       option = "magma") + 
  theme_bw() + 
  theme(text = element_text(size = t_size), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
rmse_cv_boxplot

## plot prediction performance by prevalence
# gather prediction performance values into better format for plotting
prev_df <- select(ev_df_cv, bias.name, Method, auc, rmse, prevalence, 
                  n_obs_per_sp) %>%
  gather(key = "metric", value = "value", auc:rmse)

prevalence_auc_rmse_plot <- ggplot(data = prev_df, 
                                   aes(x = prevalence, y = value)) + 
  geom_point() + 
  facet_wrap(~factor(metric, levels = c("auc", "rmse"), 
                     labels = c("AUC", "RMSE"))) + 
  xlab("Prevalence of simulated species\n(proportion of grid cells occupied)") + 
  ylab("") +
  theme_bw() + 
  theme(text = element_text(size = t_size))
prevalence_auc_rmse_plot

## plot histogram of prevalence
prev_hist <- ggplot(data = prev_df[prev_df$Method == "glm_poly_wg_block_cv" & 
                                     prev_df$bias.name == "no_bias" & 
                                     prev_df$n_obs_per_sp > 157 & 
                                     prev_df$metric == "auc", ], 
                    aes(x = prevalence)) + 
  geom_histogram() + 
  xlab("Prevalence of simulated species\n(proportion of grid cells occupied)") + 
  ylab("Number of species") + 
  theme_bw()
prev_hist

rm(prev_df)
### end plots for large community simulation ----------------------------------






### plots for small community simulation --------------------------------------
auc_smooth_all_od <- ggplot(data = preds_auc_gam_od, 
                            aes(x = n_obs_per_sp, y = pred_auc, 
                                group = factor(bias.name,
                                               levels = c("no bias", "bird", 
                                                          "bry", "moth"),
                                               labels = c(
                                                 "\nnone\n",
                                                 "\nlow\n(Birds)", 
                                                 "\nmoderate\n(Bryophytes)\n", 
                                                 "\nsevere\n(Moths)\n"),
                                               ordered = TRUE))) +
  geom_line(aes(color = factor(bias.name,
                               levels = c("no bias", "bird", "bry", "moth"),
                               labels = c(
                                 "\nnone\n", "\nlow\n(Birds)", 
                                 "\nmoderate\n(Bryophytes)\n",
                                 "\nsevere\n(Moths)\n"),
                               ordered = TRUE), 
                linetype = factor(bias.name,
                                  levels = c("no bias", "bird", "bry", "moth"),
                                  labels = c(
                                    "\nnone\n", "\nlow\n(Birds)", 
                                    "\nmoderate\n(Bryophytes)\n",
                                    "\nsevere\n(Moths)\n"),
                                  ordered = TRUE)), size = line_size) + 
  facet_wrap( ~ factor(method, 
                       levels = c("glm", "idwi"), 
                       labels = c("GLM", "Inverse\ndistance-weighted\ninterpolation"))) + 
  xlab("Average number of observations \nper species in data set") + 
  ylab("AUC\n(block cross-validation)") +
  scale_linetype_discrete(name = "Spatial\nSampling\nBias") + 
  scale_color_viridis_d(name = "Spatial\nSampling\nBias", option = "viridis") +
  theme_bw() + 
  theme(text = element_text(size = t_size), 
        axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1), 
        legend.key.width = unit(2*t_size, "points"))
auc_smooth_all_od

### end small community sim plots --------------------------------------------


### Tables --------------------------------------------------------------------
## variable importance for RMSE
rmse_summary




#### print pdf -------------------------------------------------------------
pdf("plots_supplementary.pdf")
print(auc_smooth)
print(rmse_smooth)
print(auc_smooth_noCV)
print(auc_boxplot_noCV)
print(rmse_cv_boxplot)
print(rmse_var_importance)
print(auc_smooth_all_od)
print(prevalence_auc_rmse_plot)
print(## evenness metrics correlation 
  pairs(evenness[, colnames(evenness) %in% 
                   c("camargo", "shannon_evenness", 
                     "simpson_evenness")]))
print(prev_hist)
dev.off()

ggsave("S1.svg", auc_smooth, width = 25, height = 25, units = "cm", device = "svg")
ggsave("S2.svg", rmse_smooth, width = 25, height = 25, units = "cm", device = "svg")
ggsave("S3.svg", auc_smooth_noCV, width = 25, height = 25, units = "cm", device = "svg")
ggsave("S4.svg", auc_boxplot_noCV, width = 25, height = 25, units = "cm", device = "svg")
ggsave("S5.svg", rmse_cv_boxplot, width = 25, height = 25, units = "cm", device = "svg")
ggsave("S6.svg", rmse_var_importance, width = 25, height = 25, units = "cm", device = "svg")
ggsave("S7.svg", auc_smooth_all_od, width = 25, height = 25, units = "cm", device = "svg")
ggsave("S8.svg", prevalence_auc_rmse_plot, width = 25, height = 25, units = "cm", device = "svg")
ggsave("S9.svg", prev_hist, width = 25, height = 25, units = "cm")
svg("S10.svg")
print(## evenness metrics correlation 
  pairs(evenness[, colnames(evenness) %in% 
                   c("camargo", "shannon_evenness", 
                     "simpson_evenness")]))
dev.off()
svg("S11.svg")
print(## predictor variables
  plot(mask(pred_rast_brick, ir_TM75)))
dev.off()

ggsave("S1.jpg", auc_smooth, width = 25, height = 25, units = "cm", device = "jpg")
ggsave("S2.jpg", rmse_smooth, width = 25, height = 25, units = "cm", device = "jpg")
ggsave("S3.jpg", auc_smooth_noCV, width = 25, height = 25, units = "cm", device = "jpg")
ggsave("S4.jpg", auc_boxplot_noCV, width = 25, height = 25, units = "cm", device = "jpg")
ggsave("S5.jpg", rmse_cv_boxplot, width = 25, height = 25, units = "cm", device = "jpg")
ggsave("S6.jpg", rmse_var_importance, width = 25, height = 25, units = "cm", device = "jpg")
ggsave("S7.jpg", auc_smooth_all_od, width = 25, height = 25, units = "cm", device = "jpg")
ggsave("S8.jpg", prevalence_auc_rmse_plot, width = 25, height = 25, units = "cm", device = "jpg")
ggsave("S9.jpg", prev_hist, width = 25, height = 25, units = "cm", device = "jpg")
