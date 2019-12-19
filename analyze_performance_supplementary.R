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
## last modified: 25 Sep 2019
############################################

# no need to set working directory or load packages as that is done in main.R
# library(mgcv)

k <- 7

### load data -----------------------------------------------------------------
# read in a .csv with auc results from simulation
## large community simulation results
ev_df <- read_csv("./bryophyte_template_evals_24Sep2019.csv")
df_butt <- read_csv("../sims_10.10_bryophyte_butterflyBias/bryophyte_template_butterfly_bias_evals_12Sep2019.csv")
ev_df <- bind_rows(ev_df, df_butt)

# keep only relevant columns
ev_df <- ev_df[, colnames(ev_df) %in% c("n.obs", "bias.name", "Method", "Draw", 
                                        "auc", "rmse", "time", 
                                        "prevalence", "n_obs_per_sp")]
# relabel bias factor levels
ev_df$bias.name <- factor(ev_df$bias.name, 
                          levels = c("no_bias", "least_birdNBDC", 
                                     "moderate_bryBBS", "butterflyNBDC", 
                                     "extreme_mothNBDC"), 
                          labels = c("no bias", "bird", "bry", 
                                     "butterfly", "moth"))

## small community simulation results
od_ev_df <- read_csv("../sims_10.10_odonata/odonata_template_evals_25Sep.csv")

# for now do not look at butterfly bias results b/c idwi has not been run completely
# od_butBias_df <- read_csv("../sims_10.10_odonata_butterflyBias/odonata_template_butterfly_bias_evals_25Aug2019.csv")
# od_ev_df <- bind_rows(od_ev_df, od_butBias_df)

# keep only relevant columns
od_ev_df <- od_ev_df[, colnames(od_ev_df) %in% c("n.obs", "bias.name", "Method", 
                                                 "Draw", "auc", "rmse", "time", 
                                                 "prevalence", "n_obs_per_sp")]

# relabel bias factor levels
od_ev_df$bias.name <- factor(od_ev_df$bias.name, 
                             levels = c("no_bias", "least_birdNBDC", 
                                        "moderate_bryBBS", 
                                        "extreme_mothNBDC"), 
                             labels = c("no bias", "bird", "bry", "moth"))

# remove methods that did not fully run (see simulation_tracking_status.ods)
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

## fit only random forest ----------------------------------------------------
auc_rf_noCV <- gam(auc ~ bias.name + 
                s(n_obs_per_sp, 
                  by = bias.name, 
                  k = k, bs = "cr"), 
              data = ev_df[!is.na(ev_df$auc) & 
                             ev_df$Method == "random_forest", ], 
              family = betar(link = "logit"))
auc_rf_noCV
gam.check(auc_rf_noCV) # fanning in fitted v. resids plot.  efs and k-index ok.  
plot(auc_rf_noCV, pages = 1, all.terms = T)


### generate AUC predictions to use in plotting
nobs_range <- 1:158 # set up range of n_obs_per_sp values
newdata <- data.frame(n_obs_per_sp = rep(nobs_range, 4), 
                      bias.name = c(rep("no bias", length(nobs_range)), 
                                    rep("bird", length(nobs_range)), 
                                    rep("bry", length(nobs_range)), 
                                    rep("moth", length(nobs_range))))
newdata$bias.name <- factor(newdata$bias.name, 
                            levels = c("no bias", "bird", "bry", "moth"), 
                            labels = c("no bias", "bird", "bry", "moth"))

# predictions when using a GLM SDM
preds_glm_noCV <- newdata
preds_glm_noCV$pred_auc <- predict(auc_glm_noCV, newdata = newdata, 
                                   type = "response")
preds_glm_noCV$method <- "glm"

# predictions when using a Random Forest SDM
preds_rf_noCV <- newdata
preds_rf_noCV$pred_auc <- predict(auc_rf_noCV, newdata = newdata, 
                                  type = "response")
preds_rf_noCV$method <- "rf"

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
preds_auc_gam_noCV <- bind_rows(preds_glm_noCV, preds_rf_noCV, preds_idwi_noCV, 
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
# increases width of CI, k-index ok, p-values not good
plot(auc_idwi_od, pages = 1, all.terms = T, rug = T)


### generate AUC predictions to use in plotting
nobs_range <- 3:max(od_ev_df_cv$n_obs_per_sp) # set up range of n_obs_per_sp values
newdata <- data.frame(n_obs_per_sp = rep(nobs_range, 4), 
                      bias.name = c(rep("no bias", length(nobs_range)), 
                                    rep("bird", length(nobs_range)), 
                                    rep("bry", length(nobs_range)), 
                                    rep("moth", length(nobs_range))))
newdata$bias.name <- factor(newdata$bias.name, 
                            levels = c("no bias", "bird", "bry", 
                                       "butterfly", "moth"), 
                            labels = c("no bias", "bird", "bry", 
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
# change in no-CV AUC when sample size goes from 20 to 80  
preds_auc_gam_noCV$pred_auc[preds_auc_gam_noCV$method == "rf" & 
                              preds_auc_gam_noCV$n_obs_per_sp == 80 & 
                              preds_auc_gam_noCV$bias.name == "no bias"] - 
  preds_auc_gam_noCV$pred_auc[preds_auc_gam_noCV$method == "rf" & 
                                preds_auc_gam_noCV$n_obs_per_sp == 20 & 
                                preds_auc_gam_noCV$bias.name == "no bias"]
