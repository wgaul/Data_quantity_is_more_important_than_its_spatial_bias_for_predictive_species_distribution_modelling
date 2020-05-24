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
## last modified: 22 May 2020
############################################

library(wgutil)
#library(skimr)
library(Hmisc)
library(rgdal)
library(raster)
library(parallel)
library(pROC)
library(tidyverse)
library(virtualspecies)
library(simulator) # this file was created under simulator version 0.2.0

t_size <- 19
line_size <- 1.5

sim_name <- "Bryophyte Simulation"

bias_names <- c("no_bias", "extreme_mothNBDC", 
                "butterflyNBDC", "least_birdNBDC")

## optionally, read in a .csv with auc results from simulation
## large community simulation results
ev_df <- read_csv("./bryophyte_template_evals_10April2020.csv")
# drop sample sizes that are 0.08 and 0.8 records/sp because of difficulty 
# communicating how there can be less than 1 record per species.  This does not
# change conclusions in any way.
ev_df <- ev_df[ev_df$n_obs_per_sp > 1, ]

# subset to only block CV results
ev_df_cv <- ev_df[which(ev_df$Method %in% c("glm_poly_wg_block_cv", 
                                            "brt_wg_block_cv", 
                                            "idw_interp_wg_block_cv")), ]

# subset to non-block CV results
ev_df_noCV <- ev_df[which(ev_df$Method %nin% c("glm_poly_wg_block_cv",  
                                               "brt_wg_block_cv", 
                                               "idw_interp_wg_block_cv")), ]

## small community simulation results
od_ev_df <- read_csv("./butterfly_template_evals_17April.csv")

# remove methods that did not fully run (see simulation_tracking_status.ods)
od_ev_df <- od_ev_df[od_ev_df$Method %in% c("glm_poly", "glm_poly_wg_block_cv", 
                                            "idw_interp", 
                                            "idw_interp_wg_block_cv"), ]
# get only CV results
od_ev_df_cv <- od_ev_df[od_ev_df$Method %in% c("glm_poly_wg_block_cv",
                                               "idw_interp_wg_block_cv"), ]

## load number of predictors per model results 
nvar_glm <- read_csv("nvar_df_glm_poly_wg_block_cv.csv")
nvar_idwi <- read_csv("nvar_df_idw_interp_wg_block_cv.csv")
nvar_brt <- read_csv("nvar_df_brt_wg_block_cv.csv")
# nvar_df <- bind_rows(nvar_glm, nvar_rf, nvar_idwi, nvar_brt)

# add n_obs_per_sp and AUC to nvar data frames
nvar_glm$n_obs_per_sp <- nvar_glm$nobs/community_list_lengths$bry_BBS$community_size
nvar_glm <- data.frame(left_join(nvar_glm, 
                      ev_df_cv[ev_df_cv$Method == "glm_poly_wg_block_cv", 
                               c("n.obs", "bias.name", "Draw", "auc", "rmse", 
                                 "prevalence")], 
                      by = c("nobs" = "n.obs", "bias" = "bias.name", 
                             "sp_name" = "Draw")))
nvar_brt$n_obs_per_sp <- nvar_brt$nobs/community_list_lengths$bry_BBS$community_size
nvar_brt <- data.frame(left_join(
  nvar_brt, 
  ev_df_cv[ev_df_cv$Method == "brt_wg_block_cv", 
           c("n.obs", "bias.name", "Draw", "auc", "rmse", "prevalence")], 
  by = c("nobs" = "n.obs", "bias" = "bias.name", 
         "sp_name" = "Draw")))
nvar_idwi$n_obs_per_sp <- nvar_idwi$nobs/community_list_lengths$bry_BBS$community_size
nvar_idwi <- data.frame(left_join(
  nvar_idwi, 
  ev_df_cv[ev_df_cv$Method == "idw_interp_wg_block_cv", 
           c("n.obs", "bias.name", "Draw", "auc", "rmse", "prevalence")], 
  by = c("nobs" = "n.obs", "bias" = "bias.name", 
         "sp_name" = "Draw")))

### load environmental predictor data ------------------------------------------
elev_hec <- readRDS("elevation_hec_ETOPO1.rds")
krg_mean_rr_rast <- readRDS("annual_precip_hectad.rds")
krg_mean_tx_rast <- readRDS("summer_tx_hectad.rds")
krg_mean_tn_rast <- readRDS("winter_tn_hectad.rds")
mean_pp_rast <- readRDS("mean_pp_hectad.rds")
load("corine_label_1_hectad.RData")
# Load Ireland coastline
ir <- readOGR(dsn='./data/', layer='ireland_coastline')
ir_TM75 <- spTransform(ir, CRS("+init=epsg:29903"))
rm(ir)

pred_rast_brick <- brick(list(
  "minimum temperature" = resample(krg_mean_tn_rast, krg_mean_rr_rast),
  "maximum temperature" = resample(krg_mean_tx_rast, krg_mean_rr_rast),
  "annual precipitation" = krg_mean_rr_rast,
  "atmospheric pressure" = resample(mean_pp_rast, krg_mean_rr_rast),
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


## plot methods by bias - AUC --------------------------------------------------
# use only block CV results
methods_bias_line <- ggplot(data = preds_auc_brt, aes(
  x = as.numeric(as.character(n_obs_per_sp)), 
  y = pred_auc, 
  group = factor(method, 
                 levels = c("glm", "brt", "idwi"), 
                 labels = c("GLM", 
                            "Boosted\nRegression\nTree", 
                            "Inverse distance-weighted\ninterpolation")))) +
  geom_line(aes(
    linetype = factor(method, 
                      levels = c("glm", "brt", "idwi"), 
                      labels = c("GLM",  
                                 "\nBoosted\nRegression\nTree", 
                                 "\nInverse\ndistance-weighted\ninterpolation")), 
    color = factor(method, levels = c("glm", "brt", "idwi"), 
                   labels = c("GLM",  
                              "\nBoosted\nRegression\nTree", 
                              "\nInverse\ndistance-weighted\ninterpolation"))), 
    size = line_size) +
  geom_rug(data = rug_df, aes(x = n_obs_per_sp), show.legend = F, sides = "b") +
  facet_wrap(~ factor(bias.name,
                      levels = c("no bias", "bird", 
                                 "butterfly", "moth"),
                      labels = c("no bias", "low bias\n(Birds)",
                                 "median bias\n(Butterflies)", 
                                 "severe bias\n(Moths)"),
                      ordered = TRUE)) + 
  # ggtitle("Large Community Simulation") + #simulation_name
  ylab("AUC\n(block cross-validated)") + 
  xlab("Number of observations\nper species") + 
  scale_linetype_discrete(name = "SDM\nModelling\nMethod") +
  scale_color_viridis_d(name = "SDM\nModelling\nMethod", option = "C") + 
  theme_bw() + 
  theme(text = element_text(size = t_size), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        legend.key.width = unit(2*t_size, "points"))
methods_bias_line
## end plot methods by bias - AUC -------------------------------------------


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
                       levels = c("glm", "brt", "idwi"), 
                       labels = c("GLM", 
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
                                       levels = c("no bias", "bird", 
                                                  "butterfly", "moth"),
                                       labels = c(
                                         "\nnone\n", "\nlow\n(Birds)", 
                                         "\nmedian\n(Butterflies)\n", 
                                         "\nsevere\n(Moths)\n"),
                                       ordered = TRUE))) + 
  geom_line(aes(color = factor(bias.name,
                               levels = c("no bias", "bird",
                                          "butterfly", "moth"),
                               labels = c(
                                 "\nnone\n", "\nlow\n(Birds)", 
                                 "\nmedian\n(Butterflies)\n", 
                                 "\nsevere\n(Moths)\n"),
                               ordered = TRUE), 
                linetype = factor(bias.name,
                                  levels = c("no bias", "bird", 
                                             "butterfly", "moth"),
                                  labels = c(
                                    "\nnone\n", "\nlow\n(Birds)", 
                                    "\nmedian\n(Butterflies)\n", 
                                    "\nsevere\n(Moths)\n"),
                                  ordered = TRUE)), size = line_size) + 
  geom_rug(data = rug_df, aes(x = n_obs_per_sp), show.legend = F, sides = "b") +
  # ggtitle("Large Community Simulation") + #sim_name
  facet_wrap( ~ factor(method, 
                       levels = c("glm", "brt", "idwi"), 
                       labels = c("GLM", 
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
rug_df <- rug_df[rug_df$bias.name != "butterfly" & 
                   rug_df$method %in% c("glm", "idwi", "brt"), ] # drop butterfly bias

preds_auc_gam$eval_method <- "Block CV"
preds_auc_gam_noCV$eval_method <- "no CV"
all_preds_gam <- bind_rows(preds_auc_gam, preds_auc_gam_noCV)

auc_smooth_cvNoCv <- ggplot(data = all_preds_gam, 
                            aes(x = n_obs_per_sp, y = pred_auc)) + 
  geom_line(aes(color = factor(bias.name,
                               levels = c("no bias", "bird", "butterfly", 
                                          "moth"),
                               labels = c(
                                 "\nnone\n", "\nlow\n(Birds)", 
                                 "\nmedian\n(Butterfly)\n",
                                 "\nsevere\n(Moths)\n"),
                               ordered = TRUE), 
                linetype = factor(eval_method,
                                  levels = c("Block CV", "no CV"),
                                  labels = c("Spatial\nBlock\nCV", 
                                             "\nNo CV"),
                                  ordered = TRUE)), size = line_size) + 
  geom_rug(data = rug_df, aes(x = n_obs_per_sp), show.legend = F, sides = "b") +
  facet_wrap( ~ factor(method, 
                       levels = c("glm", "brt", "idwi"), 
                       labels = c("A", 
                                  "B", 
                                  "C"))) + 
  xlab("Average number of observations \nper species in data set") + 
  ylab("AUC") +
  ylim(c(0.5, 1)) + 
  scale_linetype_discrete(name = "Model\nevaluation\nmethod") + 
  scale_color_viridis_d(name = "Spatial\nsampling\nbias", option = "viridis") +
  theme_bw() + 
  theme(text = element_text(size = t_size), 
        legend.key.width = unit(2*t_size, "points"))
auc_smooth_cvNoCv

# ev_df$eval_method <- "no CV"
# ev_df$eval_method[grepl(".*wg_blo.*", ev_df$Method)] <- "Block CV"
# ev_df$Method <- gsub("_wg_block_cv", "", ev_df$Method)
# 
# rmse_cvNoCv <- ggplot(data = ev_df[ev_df$Method %in% c("glm_poly", "brt", 
#                                                        "idw_interp") & 
#                                      ev_df$bias.name != "butterflyNBDC", ], 
#                             aes(x = n_obs_per_sp, y = rmse)) + 
#   geom_smooth(aes(color = factor(bias.name,
#                                levels = c("no_bias", "least_birdNBDC", 
#                                           "moderate_bryBBS", "extreme_mothNBDC"),
#                                labels = c(
#                                  "\nnone\n", "\nlow\n(Birds)", 
#                                  "\nmoderate\n(Bryophytes)\n",
#                                  "\nsevere\n(Moths)\n"),
#                                ordered = TRUE), 
#                 linetype = factor(eval_method,
#                                   levels = c("Block CV", "no CV"),
#                                   labels = c("Spatial\nBlock\nCV", 
#                                              "\nNo CV"),
#                                   ordered = TRUE)), 
#               method = "loess", se = F, size = line_size) + 
#   facet_wrap( ~ factor(Method, 
#                        levels = c("glm_poly", "brt", "idw_interp"), 
#                        labels = c("A", "B", "C"))) + 
#   # geom_rug(data = rug_df, aes(x = n_obs_per_sp), show.legend = F, sides = "b") +
#   xlab("Average number of observations \nper species in data set") + 
#   ylab("RMSE") +
#   scale_linetype_discrete(name = "Model\nevaluation\nmethod") + 
#   scale_color_viridis_d(name = "Spatial\nsampling\nbias", option = "viridis") +
#   theme_bw() + 
#   theme(text = element_text(size = t_size), 
#         legend.key.width = unit(2*t_size, "points"))
# rmse_cvNoCv


# ## plot AUC boxplots sample size and bias 
# auc_boxplot_noCV <- ggplot(data = ev_df_noCV[ev_df_noCV$Method %in% 
#                                                c("glm_poly",  
#                                                  "brt", "idw_interp"), ], 
#                            aes(x = factor(bias.name,
#                                           levels = c("no_bias",
#                                                      "least_birdNBDC",
#                                                      "moderate_bryBBS",
#                                                      "butterflyNBDC", 
#                                                      "extreme_mothNBDC"),
#                                           labels = c("none", "low", "moderate", 
#                                                      "median", "severe"),
#                                           ordered = TRUE), y = auc, 
#                                fill = factor(round(n_obs_per_sp, 0), 
#                                              levels = c("0", "1", "2", "4", "8", 
#                                                         "20", "79", "158"), 
#                                              labels = c("< 1", "1", "2", "4", "8", 
#                                                         "20", "79", "158")))) +
#   geom_boxplot(varwidth = TRUE) +
#   facet_wrap(~ factor(Method, 
#                       levels = c("glm_poly", 
#                                  "brt", 
#                                  "idw_interp"), 
#                       labels = c("GLM", 
#                                  "Boosted\nRegression\nTrees", 
#                                  "Inverse\ndistance-weighted\ninterpolation"))) + 
#   xlab("Sampling Bias") + 
#   ylab("AUC\n(no cross-validation)") + 
#   scale_fill_viridis_d(name = "Average\nnumber of\nrecords\nper species", 
#                        option = "magma") + 
#   theme_bw() + 
#   theme(text = element_text(size = t_size), 
#         axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1)) 
# auc_boxplot_noCV

###################################
## RMSE plots 
## BRT prediction of RMSE by bias
rmse_line <- ggplot(data = preds_rmse_evaluation_brt, 
                    aes(x = n_obs_per_sp, y = pred_rmse, 
                        group = factor(bias.name,
                                       levels = c("no bias", "bird",  
                                                  "butterfly", "moth"),
                                       labels = c(
                                         "\nnone\n", "\nlow\n", 
                                         "\nmedian\n", "\nsevere\n"),
                                       ordered = TRUE))) + 
  geom_line(aes(color = factor(bias.name,
                               levels = c("no bias", "bird", 
                                          "butterfly", "moth"),
                               labels = c(
                                 "\nnone\n", "\nlow\n", 
                                 "\nmedian\n", "\nsevere\n"),
                               ordered = TRUE)), size = line_size) + 
  geom_rug(data = rug_df, aes(x = n_obs_per_sp), 
           show.legend = F, sides = "b") +
  # ggtitle("Large Community Simulation") + #simulation_name
  facet_wrap( ~ factor(method, 
                       levels = c("glm", "brt", "idwi"), 
                       labels = c("A", "B", "C"))) + 
  xlab("Average number of records \nper species in data set") + 
  ylab("RMSE\n(block cross-validated)") +
  # scale_linetype_discrete(name = "Spatial\nSampling\nBias") + 
  scale_color_viridis_d(name = "Spatial\nSampling\nBias", option = "magma", 
                        begin = 0, end = 0.65) +
  theme_bw() + 
  theme(text = element_text(size = t_size*1.2), 
        legend.key.width = unit(1.8*t_size, "points"), 
        axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1))
rmse_line

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
                                     prev_df$n_obs_per_sp == 200 & 
                                     prev_df$metric == "auc", ], 
                    aes(x = prevalence)) + 
  geom_histogram() + 
  xlab("Prevalence of simulated species\n(proportion of grid cells occupied)") + 
  ylab("Number of species") + 
  theme_bw()  + 
  theme(text = element_text(size = t_size))
prev_hist

rm(prev_df)

### plot number of predictors per model ---------------------------------------
## GLM
nvar_glm_plot <- ggplot(data = nvar_glm[nvar_glm$n_obs_per_sp > 1, ], 
                    aes(x = n_obs_per_sp, y = n_predictors, 
                        color = factor(bias, 
                                       levels = c("no_bias", 
                                                  "least_birdNBDC", 
                                                  "median_butterflyNBDC", 
                                                  "extreme_mothNBDC"), 
                                       labels = c("none", "low", "median", 
                                                  "severe"), 
                                       ordered = TRUE))) + 
  geom_boxplot(aes(group = factor(as.character(nobs)))) + 
  geom_smooth(method = "loess", alpha = 0.3) + 
  xlab("Average number of records per species") + 
  ylab("Number of terms in model") + 
  scale_color_viridis_d(name = "Spatial\nSampling\nBias", option = "magma", 
                        begin = 0, end = 0.65) + 
  # scale_x_log10() +
  theme_bw() + 
  theme(text = element_text(size = t_size*1.2), 
        legend.key.width = unit(1.8*t_size, "points"))
nvar_glm_plot

nvar_auc_glm <- ggplot(data = nvar_glm[nvar_glm$n_obs_per_sp > 1 & 
                                         !is.na(nvar_glm$auc), ], 
                       aes(x = n_predictors, y = auc, 
                           colour = factor(bias, 
                                           levels = c("no_bias", 
                                                      "least_birdNBDC", 
                                                      "median_butterflyNBDC", 
                                                      "extreme_mothNBDC"), 
                                           labels = c("none", "low", "median", 
                                                      "severe"), 
                                           ordered = TRUE))) + 
  geom_boxplot(aes(group = factor(n_predictors)), varwidth = T, 
               color = "black") + 
   geom_smooth(method = "loess", se = F) +
  scale_color_viridis_d(name = "Spatial\nSampling\nBias", option = "magma", 
                        begin = 0, end = 0.65, alpha = 0.6) +
  facet_wrap(~factor(as.character(n_obs_per_sp), 
                     levels = c("2", "5", "10", "50", "100", "200"), 
                     labels = c("A", "B", "C", "D", "E", "F"))) + 
  xlab("Number of terms in model") + 
  ylab("AUC (spatial block CV)") + 
  # scale_x_log10() + 
  theme_bw() + 
  theme(text = element_text(size = t_size*1.1), 
        legend.key.width = unit(1.8*t_size, "points"), 
        axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1), 
        panel.spacing.x = unit(4, "mm"))
nvar_auc_glm

# auc_prevalence_glm <- ggplot(data = nvar_glm[nvar_glm$n_obs_per_sp > 1, ], 
#                              aes(x = prevalence, y = auc)) + 
#   geom_point() + 
#   geom_smooth(aes(colour = factor(round(n_obs_per_sp),
#                                   levels = c("2", "4", "8", "20", "79",
#                                              "158"))),
#               method = "loess") +
#   facet_wrap(~factor(bias, 
#                      levels = c("no_bias", "least_birdNBDC",
#                                 "moderate_bryBBS", "extreme_mothNBDC"), 
#                      labels = c("no bias", "low bias", 
#                                 "moderate bias", "severe bias"))) + 
#   xlab("Species prevalence") + 
#   ylab("AUC (spatial block CV)") + 
#   scale_x_log10() + 
#   scale_color_viridis_d(name = "Average\nnumber of\nrecords\nper species")  
# auc_prevalence_glm

## bin prevalence
# GLM
nvar_glm$prevalence_bin <- NA
nvar_glm$prevalence_bin[nvar_glm$prevalence < 0.15] <- "0 - 0.149" 
nvar_glm$prevalence_bin[nvar_glm$prevalence >= 0.15 & 
                          nvar_glm$prevalence < 0.3] <- "0.15 - 0.29" 
nvar_glm$prevalence_bin[nvar_glm$prevalence >= 0.3 & 
                          nvar_glm$prevalence < 0.45] <- "0.3 - 0.449" 
nvar_glm$prevalence_bin[nvar_glm$prevalence >= 0.45 & 
                          nvar_glm$prevalence < 0.7] <- "0.45 - 0.7" 
# BRT
nvar_brt$prevalence_bin <- NA
nvar_brt$prevalence_bin[nvar_brt$prevalence < 0.15] <- "0 - 0.149" 
nvar_brt$prevalence_bin[nvar_brt$prevalence >= 0.15 & 
                          nvar_brt$prevalence < 0.3] <- "0.15 - 0.29" 
nvar_brt$prevalence_bin[nvar_brt$prevalence >= 0.3 & 
                          nvar_brt$prevalence < 0.45] <- "0.3 - 0.449" 
nvar_brt$prevalence_bin[nvar_brt$prevalence >= 0.45 & 
                          nvar_brt$prevalence < 0.7] <- "0.45 - 0.7"
# idwi
nvar_idwi$prevalence_bin <- NA
nvar_idwi$prevalence_bin[nvar_idwi$prevalence < 0.15] <- "0 - 0.149" 
nvar_idwi$prevalence_bin[nvar_idwi$prevalence >= 0.15 & 
                          nvar_idwi$prevalence < 0.3] <- "0.15 - 0.29" 
nvar_idwi$prevalence_bin[nvar_idwi$prevalence >= 0.3 & 
                           nvar_idwi$prevalence < 0.45] <- "0.3 - 0.449" 
nvar_idwi$prevalence_bin[nvar_idwi$prevalence >= 0.45 & 
                          nvar_idwi$prevalence < 0.7] <- "0.45 - 0.7"

nvar_all <- bind_rows(nvar_glm[, colnames(nvar_brt) %nin% 
                                 c("elapsed_minutes", "n_brt_trees")], 
                      nvar_brt[, colnames(nvar_brt) %nin% 
                                 c("elapsed_minutes", "n_brt_trees")], 
                      nvar_idwi[, colnames(nvar_brt) %nin% 
                                  c("elapsed_minutes", "n_brt_trees")])

auc_prevalence_bin <- ggplot(
  data = nvar_all[nvar_all$n_obs_per_sp > 1 & !is.na(nvar_all$prevalence), ], 
  aes(x = factor(prevalence_bin, 
                 levels = c("0 - 0.149", 
                            "0.15 - 0.29", 
                            "0.3 - 0.449", 
                            "0.45 - 0.7")), 
      y = auc)) + 
  geom_boxplot() + 
  facet_wrap(~factor(method, 
                     levels = c("glm_poly_wg_block_cv", 
                                "brt_wg_block_cv", 
                                "idw_interp_wg_block_cv"), 
                     labels = c("A", "B", "C"))) + 
  xlab("Species prevalence") + 
  ylab("AUC (spatial block CV)") +   
  theme_bw() + 
  theme(text = element_text(size = t_size*1.2), 
        legend.key.width = unit(1.8*t_size, "points"),
        axis.text.x = element_text(angle = 35, vjust = 1, hjust = 1))
auc_prevalence_bin


# auc_prevalence_bin_glm <- ggplot(data = nvar_glm[nvar_glm$n_obs_per_sp > 1, ], 
#                              aes(x = factor(prevalence_bin, 
#                                             levels = c("0 - 0.149", 
#                                                        "0.15 - 0.29", 
#                                                        "0.3 - 0.449", 
#                                                        "0.45 - 0.7")), 
#                                  y = auc)) + 
#   geom_boxplot() + 
#   facet_wrap(~factor(bias, 
#                      levels = c("no_bias", "least_birdNBDC",
#                                 "moderate_bryBBS", "extreme_mothNBDC"), 
#                      labels = c("no bias", "low bias", 
#                                 "moderate bias", "severe bias"))) + 
#   xlab("Species prevalence") + 
#   ylab("AUC (spatial block CV)") +   
#   theme_bw() + 
#   theme(text = element_text(size = t_size*1.2), 
#         legend.key.width = unit(1.8*t_size, "points"),
#         axis.text.x = element_text(angle = 35, vjust = 1, hjust = 1))
# auc_prevalence_bin_glm

# auc_prevalence_bin_brt <- ggplot(data = nvar_brt[nvar_brt$n_obs_per_sp > 1, ], 
#                                 aes(x = factor(prevalence_bin, 
#                                                levels = c("0 - 0.149", 
#                                                           "0.15 - 0.29", 
#                                                           "0.3 - 0.449", 
#                                                           "0.45 - 0.7")), 
#                                     y = auc)) + 
#   geom_boxplot() + 
#   facet_wrap(~factor(bias, 
#                      levels = c("no_bias", "least_birdNBDC",
#                                 "moderate_bryBBS", "extreme_mothNBDC"), 
#                      labels = c("no bias", "low bias", 
#                                 "moderate bias", "severe bias"))) + 
#   xlab("Species prevalence") + 
#   ylab("AUC (spatial block CV)") +   
#   theme_bw() + 
#   theme(text = element_text(size = t_size*1.2), 
#         legend.key.width = unit(1.8*t_size, "points"),
#         axis.text.x = element_text(angle = 35, vjust = 1, hjust = 1))
# auc_prevalence_bin_brt
# 
# auc_prevalence_bin_idwi <- ggplot(
#   data = nvar_idwi[nvar_idwi$n_obs_per_sp > 1, ], 
#   aes(x = factor(prevalence_bin, 
#                  levels = c("0 - 0.149", 
#                             "0.15 - 0.29", 
#                             "0.3 - 0.449", 
#                             "0.45 - 0.7")), 
#       y = auc)) + 
#   geom_boxplot() + 
#   facet_wrap(~factor(bias, 
#                      levels = c("no_bias", "least_birdNBDC",
#                                 "moderate_bryBBS", "extreme_mothNBDC"), 
#                      labels = c("no bias", "low bias", 
#                                 "moderate bias", "severe bias"))) + 
#   xlab("Species prevalence") + 
#   ylab("AUC (spatial block CV)") +   
#   theme_bw() + 
#   theme(text = element_text(size = t_size*1.2), 
#         legend.key.width = unit(1.8*t_size, "points"),
#         axis.text.x = element_text(angle = 35, vjust = 1, hjust = 1))
# auc_prevalence_bin_idwi

# ## plot number of terminal nodes for RF
# n_term_node_rf_plot <- ggplot(
#   data = nvar_rf, 
#   aes(x = n_train_detections / (n_train_nonDets + n_train_detections), 
#       y = mean_n_term_nodes, colour = factor(bias))) + 
#   geom_point() + 
#   facet_wrap(~factor(round(n_obs_per_sp))) + 
#   xlab("Species prevalence in training data")
# n_term_node_rf_plot

## plot number of trees for BRT 
# this tells us whether we were being limited by computation time
n_brt_tree_hist <- ggplot(data = nvar_brt[which(nvar_brt$n_obs_per_sp >1), ], 
                          aes(x = n_brt_trees)) + 
  geom_histogram() + 
  facet_wrap(~factor(bias, 
                     levels = c("no_bias", "least_birdNBDC", 
                                "median_butterflyNBDC", "extreme_mothNBDC"), 
                     labels = c("no bias", "low bias", 
                                "median bias", "severe bias"))) + 
  xlab("Number of trees") +
  theme_bw()
n_brt_tree_hist

n_brt_trees <- ggplot(data = nvar_brt[which(nvar_brt$n_obs_per_sp >1), ], 
                      aes(x = n_obs_per_sp, y = n_brt_trees)) + 
  # geom_point() + 
  geom_boxplot(aes(group = factor(nobs)), varwidth = T) + 
  facet_wrap(~factor(bias, 
                     levels = c("no_bias", "least_birdNBDC", 
                                "median_butterflyNBDC", "extreme_mothNBDC"), 
                     labels = c("no bias", "low bias", 
                                "median bias", "severe bias"))) + 
  # geom_smooth(method = "loess") + 
  xlab("Number of records per species") + 
  ylab("Number of trees") + 
  scale_x_log10() +
  theme_bw() + 
  theme(text = element_text(size = t_size*1.2), 
        legend.key.width = unit(1.8*t_size, "points"),
        axis.text.x = element_text(angle = 35, vjust = 1, hjust = 1))
n_brt_trees

auc_ntree <- ggplot(data = nvar_brt, 
                    aes(x = n_brt_trees, y = auc)) + 
  geom_boxplot(aes(group = factor(n_brt_trees))) + 
  ylab("AUC") + xlab("Number of trees used") + 
  theme_bw() + 
  theme(text = element_text(size = t_size*1.2), 
        legend.key.width = unit(1.8*t_size, "points"),
        axis.text.x = element_text(angle = 35, vjust = 1, hjust = 1))
auc_ntree

### Plot performance by number of detections in training set
nvar_glm$ntrain_bins <- NA
nvar_glm$ntrain_bins[nvar_glm$n_train_detections <= 119] <- "0 - 119"
nvar_glm$ntrain_bins[nvar_glm$n_train_detections > 119 & 
                       nvar_glm$n_train_detections <= 239] <- "120 - 239"
nvar_glm$ntrain_bins[nvar_glm$n_train_detections > 239 & 
                       nvar_glm$n_train_detections <= 359] <- "239 - 359"
nvar_glm$ntrain_bins[nvar_glm$n_train_detections > 359 & 
                       nvar_glm$n_train_detections <= 479] <- "360 - 479"
nvar_glm$ntrain_bins[nvar_glm$n_train_detections > 479 & 
                       nvar_glm$n_train_detections <= 600] <- "480 - 600"

auc_ntrain_glm <- ggplot(data = nvar_glm[!is.na(nvar_glm$ntrain_bins), ], 
                         aes(x = ntrain_bins, y = auc)) + 
  geom_boxplot(aes(group = factor(ntrain_bins))) + 
  facet_wrap(~factor(bias)) + 
  xlab("Number of detections in training set") + ylab("AUC") + 
  theme_bw() + 
  theme(axis.text = element_text(angle = 35, hjust = 1, vjust = 1))
auc_ntrain_glm

nvar_brt$ntrain_bins <- NA
nvar_brt$ntrain_bins[nvar_brt$n_train_detections <= 119] <- "0 - 119"
nvar_brt$ntrain_bins[nvar_brt$n_train_detections > 119 & 
                       nvar_brt$n_train_detections <= 239] <- "120 - 239"
nvar_brt$ntrain_bins[nvar_brt$n_train_detections > 239 & 
                       nvar_brt$n_train_detections <= 359] <- "239 - 359"
nvar_brt$ntrain_bins[nvar_brt$n_train_detections > 359 & 
                       nvar_brt$n_train_detections <= 479] <- "360 - 479"
nvar_brt$ntrain_bins[nvar_brt$n_train_detections > 479 & 
                       nvar_brt$n_train_detections <= 600] <- "480 - 600"

auc_ntrain_brt <- ggplot(data = nvar_brt[!is.na(nvar_brt$ntrain_bins), ], 
                         aes(x = ntrain_bins, y = auc)) + 
  geom_boxplot(aes(group = factor(ntrain_bins))) + 
  facet_wrap(~factor(bias)) + 
  xlab("Number of detections in training set") + ylab("AUC") + 
  theme_bw() + 
  theme(axis.text = element_text(angle = 35, hjust = 1, vjust = 1))
auc_ntrain_brt

### Plot performance by number of detections in test set
nvar_all$ntest_bins <- NA
nvar_all$ntest_bins[nvar_all$n_test_detections <= 59] <- "0 - 59"
nvar_all$ntest_bins[nvar_all$n_test_detections > 59 & 
                      nvar_all$n_test_detections <= 119] <- "60 - 119"
nvar_all$ntest_bins[nvar_all$n_test_detections > 119 & 
                      nvar_all$n_test_detections <= 179] <- "120 - 179"
nvar_all$ntest_bins[nvar_all$n_test_detections > 179 & 
                      nvar_all$n_test_detections <= 239] <- "180 - 239"
nvar_all$ntest_bins[nvar_all$n_test_detections > 239 & 
                      nvar_all$n_test_detections <= 300] <- "240 - 300"

auc_ntest_all <- ggplot(data = nvar_all[!is.na(nvar_all$ntest_bins), ], 
                         aes(x = factor(ntest_bins, 
                                        levels = c("0 - 59", "60 - 119", 
                                                   "120 - 179", "180 - 239", 
                                                   "240 - 300")), 
                             y = auc)) + 
  geom_boxplot(aes(group = factor(ntest_bins))) + 
  facet_wrap(~factor(method, 
                     levels = c("glm_poly_wg_block_cv", 
                                "brt_wg_block_cv", 
                                "idw_interp_wg_block_cv"), 
                     labels = c("A", "B", "C"))) + 
  xlab("Number of detections in test set") + ylab("AUC (spatial block CV)") + 
  theme_bw() + 
  theme(text = element_text(size = t_size*1.2), 
        axis.text = element_text(angle = 35, hjust = 1, vjust = 1))
auc_ntest_all
### end plots for large community simulation ----------------------------------





### plots for small community simulation --------------------------------------
auc_smooth_all_od <- ggplot(data = preds_auc_gam_od, 
                            aes(x = n_obs_per_sp, y = pred_auc, 
                                group = factor(bias.name,
                                               levels = c("no bias", "bird", 
                                                          "butterfly", "moth"),
                                               labels = c(
                                                 "\nnone\n",
                                                 "\nlow\n(Birds)", 
                                                 "\nmedian\n(Butterfly)\n", 
                                                 "\nsevere\n(Moths)\n"),
                                               ordered = TRUE))) +
  geom_line(aes(color = factor(bias.name,
                               levels = c("no bias", "bird", "butterfly", 
                                          "moth"),
                               labels = c(
                                 "\nnone\n", "\nlow\n(Birds)", 
                                 "\nmedian\n(Butterfly)\n",
                                 "\nsevere\n(Moths)\n"),
                               ordered = TRUE), 
                linetype = factor(bias.name,
                                  levels = c("no bias", "bird", "butterfly",
                                             "moth"),
                                  labels = c(
                                    "\nnone\n", "\nlow\n(Birds)", 
                                    "\nmedian\n(Butterfly)\n",
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
## Adjusted R2 of GAMs for AUC
var_imp_r2 # from analyze_performance_GAM.R script

## variable importance for RMSE
rmse_summary



#### print plots to files ----------------------------------------------------
ggsave("S3.jpg", prevalence_auc_rmse_plot, width = 25, height = 25, 
       units = "cm", device = "jpg")
ggsave("S4.jpg", auc_ntest_all, width = 25, height = 25, units = "cm", 
       device = "jpg")
ggsave("S5.jpg", prev_hist, width = 25, height = 25, units = "cm", device = "jpg")
ggsave("S6.jpg", rmse_line, width = 25, height = 25, 
       units = "cm", device = "jpg")
ggsave("S7.jpg", auc_smooth_cvNoCv, width = 25, height = 25, units = "cm", 
       device = "jpg")
ggsave("S8.jpg", nvar_glm_plot, width = 25, height = 25, 
       units = "cm", device = "jpg")
ggsave("S9.jpg", nvar_auc_glm, width = 25, height = 25, units = "cm", device = "jpg")
ggsave("S10.jpg", n_brt_trees, width = 25, height = 25, units = "cm", 
       device = "jpg")
ggsave("S11.jpg", auc_ntree, width = 25, height = 25, units = "cm", 
       device = "jpg")
ggsave("S12.jpg", auc_smooth_all_od, width = 25, height = 25, units = "cm", 
       device = "jpg")



## each figure as a pdf ---------------------------------------------------
pdf("S1.pdf")
print(## predictor variables
  plot(mask(pred_rast_brick, ir_TM75), axes = F))
dev.off()
pdf("S3.pdf")
rmse_line
dev.off()
pdf("S4.pdf")
prevalence_auc_rmse_plot
dev.off()
pdf("S5.pdf")
auc_prevalence_bin
dev.off()
pdf("S6.pdf")
auc_ntest_all
dev.off()
pdf("S7.pdf")
prev_hist
dev.off()
pdf("S8.pdf")
auc_smooth_cvNoCv
dev.off()
pdf("S9.pdf")
nvar_glm_plot
dev.off()
pdf("S10.pdf")
nvar_auc_glm
dev.off()
pdf("S11.pdf")
n_brt_trees
dev.off()
pdf("S12.pdf")
auc_ntree
dev.off()
pdf("S13.pdf")
auc_smooth_all_od
dev.off()

