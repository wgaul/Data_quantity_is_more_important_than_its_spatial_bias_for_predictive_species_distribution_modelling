############################################
## This analyzes prediction performance of models from the simulation
## This is for non-block-CV methods and others for supplementary materials
## This is for simulation 10.10 
## This script is meant to be sourced from within main.R
##
## inputs:  workspace from main.R
## 
## author: Willson Gaul wgaul@hotmail.com
## created: 24 Sep 2019
## last modified: 22 May 2020
############################################

# no need to set working directory or load packages as that is done in main.R
# library(mgcv)

k <- 5

### load data -----------------------------------------------------------------
# read in a .csv with auc results from simulation
## large community simulation results
ev_df <- read_csv("./bryophyte_template_evals_10April2020.csv")

# keep only relevant columns
ev_df <- ev_df[, colnames(ev_df) %in% c("n.obs", "bias.name", "Method", "Draw", 
                                        "auc", "rmse", "time", 
                                        "prevalence", "n_obs_per_sp")]
# relabel bias factor levels
ev_df$bias.name <- factor(ev_df$bias.name, 
                          levels = c("no_bias", "least_birdNBDC", 
                                     "median_butterflyNBDC", 
                                     "extreme_mothNBDC"), 
                          labels = c("no bias", "bird",
                                     "butterfly", "moth"))
ev_df <- ev_df[ev_df$Method %nin% c("random_forest_wg_block_cv", 
                                    "random_forest", "occ_det", 
                                    "occ_det_wg_block_cv"), ]

## small community simulation results
od_ev_df <- read_csv("./butterfly_template_evals_17April.csv")

# keep only relevant columns
od_ev_df <- od_ev_df[, colnames(od_ev_df) %in% c("n.obs", "bias.name", "Method", 
                                                 "Draw", "auc", "rmse", "time", 
                                                 "prevalence", "n_obs_per_sp")]

# relabel bias factor levels
od_ev_df$bias.name <- factor(od_ev_df$bias.name, 
                             levels = c("no_bias", "least_birdNBDC", 
                                        "median_butterflyNBDC", 
                                        "extreme_mothNBDC"), 
                             labels = c("no bias", "bird", "butterfly", "moth"))

# remove methods that did not fully run (if any)
od_ev_df <- od_ev_df[od_ev_df$Method %in% c("glm_poly", "glm_poly_wg_block_cv", 
                                            "idw_interp", 
                                            "idw_interp_wg_block_cv"), ]
# get only CV results
od_ev_df_cv <- od_ev_df[od_ev_df$Method %in% c("glm_poly_wg_block_cv",
                                               "idw_interp_wg_block_cv"), ]
### end load data -------------------------------------------------------------


### large community sim no CV models ------------------------------------------
### fit each method separately so I do not have to model method in 
### interaction, and I can then model nobs*bias 
## fit only glm -----------------------------------------------------
auc_glm_noCV <- gam(auc ~ bias.name + 
                      s(n_obs_per_sp, 
                        by = bias.name, 
                        k = k, bs = "cr"), 
                    data = ev_df[!is.na(ev_df$auc) & 
                                   ev_df$Method == "glm_poly", ], 
                    family = betar(link = "logit"), 
                    select = TRUE)
auc_glm_noCV
gam.check(auc_glm_noCV) # some fanning in fitted vs. residuals plot.  edf and k-index ok.
plot(auc_glm_noCV, pages = 1, all.terms = T)

## fit only idwi -----------------------------------------------------
auc_idwi_noCV <- gam(auc ~ bias.name + 
                  s(n_obs_per_sp, 
                    by = bias.name, 
                    k = k, bs = "cr"), 
                data = ev_df[!is.na(ev_df$auc) & 
                               ev_df$Method == "idw_interp", ], 
                family = betar(link = "logit"))
auc_idwi_noCV
gam.check(auc_idwi_noCV) # a few big residuals. edf and k-index ok
plot(auc_idwi_noCV, pages = 1, all.terms = T)


## fit only brt -----------------------------------------------------
auc_brt_noCV <- gam(auc ~ bias.name + 
                 s(n_obs_per_sp, 
                   by = bias.name, 
                   k = k, bs = "cr"), 
               data = ev_df[!is.na(ev_df$auc) & 
                              ev_df$Method == "brt", ], 
               family = betar(link = "logit"))
auc_brt_noCV
gam.check(auc_brt_noCV) # minor fanning in resids v fitteds, and a few big residuals.  
# edf and k-index ok.  
plot(auc_brt_noCV, pages = 1, all.terms = T)


### generate AUC predictions to use in plotting
nobs_range <- 2:200 # set up range of n_obs_per_sp values
newdata <- data.frame(n_obs_per_sp = rep(nobs_range, 4), 
                      bias.name = c(rep("no bias", length(nobs_range)), 
                                    rep("bird", length(nobs_range)), 
                                    rep("butterfly", length(nobs_range)), 
                                    rep("moth", length(nobs_range))))
newdata$bias.name <- factor(newdata$bias.name, 
                            levels = c("no bias", "bird", "butterfly", "moth"), 
                            labels = c("no bias", "bird", "butterfly", "moth"))

# predictions when using a GLM SDM
preds_glm_noCV <- newdata
preds_glm_noCV$pred_auc <- predict(auc_glm_noCV, newdata = newdata, 
                                   type = "response")
preds_glm_noCV$method <- "glm"

# predictions when using an inverse distance-weighted interpolation SDM
preds_idwi_noCV <- newdata
preds_idwi_noCV$pred_auc <- predict(auc_idwi_noCV, newdata = newdata, 
                                    type = "response")
preds_idwi_noCV$method <- "idwi"

# predictions when using a boosted classification tree SDM
preds_brt_noCV <- newdata
preds_brt_noCV$pred_auc <- predict(auc_brt_noCV, newdata = newdata, 
                                   type = "response")
preds_brt_noCV$method <- "brt"

# bind all preditions into a single data frame
preds_auc_gam_noCV <- bind_rows(preds_glm_noCV, preds_idwi_noCV, 
                            preds_brt_noCV)

# add simpson evenness values to predictions
# temp_evenness should be in workspace, it was created in analyze_performance.R
preds_auc_gam_noCV <- left_join(preds_auc_gam_noCV, temp_evenness, 
                            by = c("bias.name" = "taxon"))
### end large community sim no CV models ---------------------------------------

### small community sim models (for contour plot) ------------------------------
## fit only glm -----------------------------------------------------
auc_glm_od <- gam(auc ~ bias.name + 
                    s(n_obs_per_sp, 
                      by = bias.name, 
                      k = k, bs = "cr"), 
                  data = data.frame(od_ev_df_cv[!is.na(od_ev_df_cv$auc) & 
                                                  od_ev_df_cv$Method == "glm_poly_wg_block_cv", ]), 
                  family = betar(link = "logit"), 
                  select = TRUE)
auc_glm_od
gam.check(auc_glm_od) # minor fanning in fitted vs. residuals plot.  
# edf high but raising to k = 8 does not change shape of smooths much, though it 
# increases width of CI, k-index ok, p-values not good
plot(auc_glm_od, pages = 1, all.terms = T, rug = T)


## fit only idwi -----------------------------------------------------
auc_idwi_od <- gam(auc ~ bias.name + 
                     s(n_obs_per_sp, 
                       by = bias.name, 
                       k = k, bs = "cr"), 
                   data = od_ev_df_cv[!is.na(od_ev_df_cv$auc) & 
                                        od_ev_df_cv$Method == "idw_interp_wg_block_cv", ], 
                   family = betar(link = "logit"))
auc_idwi_od
gam.check(auc_idwi_od) # minor fanning in fitted vs. resids plot.  
# edf high but raising to k = 8 does not change shape of smooths much, though it 
# increases width of CI, k-index ok
plot(auc_idwi_od, pages = 1, all.terms = T, rug = T)


### generate AUC predictions to use in plotting
nobs_range <- 2:max(od_ev_df_cv$n_obs_per_sp) # set up range of n_obs_per_sp values
newdata <- data.frame(n_obs_per_sp = rep(nobs_range, 4), 
                      bias.name = c(rep("no bias", length(nobs_range)), 
                                    rep("bird", length(nobs_range)), 
                                    rep("butterfly", length(nobs_range)), 
                                    rep("moth", length(nobs_range))))
newdata$bias.name <- factor(newdata$bias.name, 
                            levels = c("no bias", "bird", 
                                       "butterfly", "moth"), 
                            labels = c("no bias", "bird",  
                                       "butterfly", "moth"))

# predictions when using a GLM SDM
preds_glm_od <- newdata
preds_glm_od$pred_auc <- predict(auc_glm_od, newdata = newdata, 
                                 type = "response")
preds_glm_od$method <- "glm"

# predictions when using an inverse distance-weighted interpolation SDM
preds_idwi_od <- newdata
preds_idwi_od$pred_auc <- predict(auc_idwi_od, newdata = newdata, 
                                  type = "response")
preds_idwi_od$method <- "idwi"

# bind all preditions into a single data frame
preds_auc_gam_od <- bind_rows(preds_glm_od, preds_idwi_od)

# add simpson evenness values to predictions
# temp_evenness should be in workspace, it was created in analyze_performance.R
preds_auc_gam_od <- left_join(preds_auc_gam_od, temp_evenness, 
                              by = c("bias.name" = "taxon"))

##### Numbers and stats for paper ----------------------------------------------

