############################################
## This analyzes prediction performance of models from the simulation
## This is for simulation 10.10 
## This script is meant to be sourced from within main.R
##
## inputs:  workspace from main.R
## 
## author: Willson Gaul wgaul@hotmail.com
## created: 5 Sep 2019
## last modified: 10 April 2020
############################################

# no need to set working directory or load packages as that is done in main.R
# library(mgcv)
# library(DHARMa)

k <- 5

### load data -----------------------------------------------------------------
# optionally, read in a .csv with auc results from simulation
ev_df <- read_csv("./bryophyte_template_evals_10April2020.csv")

# keep only relevant columns
ev_df <- ev_df[, colnames(ev_df) %in% 
                 c("n.obs", "bias.name", "Method", "Draw", 
                   "auc", "rmse", "time", 
                   "prevalence", "n_obs_per_sp")]

# subset to only block CV results
ev_df_cv <- ev_df[which(ev_df$Method %in% c("glm_poly_wg_block_cv", 
                                            "brt_wg_block_cv", 
                                            "idw_interp_wg_block_cv")), ]

ev_df_cv$bias.name <- factor(ev_df_cv$bias.name, 
                             levels = c("no_bias", "least_birdNBDC", 
                                        "median_butterflyNBDC", 
                                        "extreme_mothNBDC"), 
                             labels = c("no bias", "bird",  
                                        "butterfly", "moth"))
### end load data -------------------------------------------------------------

##### model AUC -------------------------------------------------------------
# fit a GAM to model AUC as a function of overall sample size (n_obs_per_sp),
# spatial sampling bias (bias.name), and modeling method (Method)

# from looking at boxplots, I expect AUC to vary by Method, nobs, Method*nobs, 
# bias, bias*Method (e.g. bias matters more for brt than idwi), and perhaps
# bias*nobs*method (e.g. RF high nobs matters at low bias but not high bias)

# because I can't fit a 3-way interaction with mgcv, I will fit each method 
# separately so I do not have to model method in an interaction, and I can then 
# model nobs*bias 
## fit only glm -----------------------------------------------------
auc_glm <- gam(auc ~ bias.name + 
                 s(n_obs_per_sp, 
                   by = bias.name, 
                   k = k, bs = "cr"), 
               data = ev_df_cv[!is.na(ev_df_cv$auc) & 
                                 ev_df_cv$Method == "glm_poly_wg_block_cv", ], 
               family = betar(link = "logit"), 
               select = TRUE)
auc_glm
gam.check(auc_glm) # some fanning in fitted vs. residuals plot.  edf and k-index ok.
plot(auc_glm, pages = 1, all.terms = T)

## fit only idwi -----------------------------------------------------
auc_idwi <- gam(auc ~ bias.name + 
                 s(n_obs_per_sp, 
                   by = bias.name, 
                   k = k, bs = "cr"), 
               data = ev_df_cv[!is.na(ev_df_cv$auc) & 
                                 ev_df_cv$Method == "idw_interp_wg_block_cv", ], 
               family = betar(link = "logit"))
auc_idwi
gam.check(auc_idwi) # minor fanning in fitted vs. resids plot.  
# edf and k-index ok, despite significant p-values.  Nothing much changes if I 
# increase k to 6.  
plot(auc_idwi, pages = 1, all.terms = T)


## fit only brt -----------------------------------------------------
auc_gam_brt <- gam(auc ~ bias.name + 
                  s(n_obs_per_sp, 
                    by = bias.name, 
                    k = k, bs = "cr"), 
                data = ev_df_cv[!is.na(ev_df_cv$auc) & 
                                  ev_df_cv$Method == "brt_wg_block_cv", ], 
                family = betar(link = "logit"))
auc_gam_brt
gam.check(auc_gam_brt) # edf and k-index ok.  
plot(auc_gam_brt, pages = 1, all.terms = T)

### generate AUC predictions to use in plotting
nobs_range <- 2:200 # set up range of n_obs_per_sp values
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
preds_glm <- newdata
preds_glm$pred_auc <- predict(auc_glm, newdata = newdata, 
                              type = "response")
preds_glm$method <- "glm"

# predictions when using an inverse distance-weighted interpolation SDM
preds_idwi <- newdata
preds_idwi$pred_auc <- predict(auc_idwi, newdata = newdata, 
                             type = "response")
preds_idwi$method <- "idwi"

# predictions when using a boosted classification tree SDM
preds_brt <- newdata
preds_brt$pred_auc <- predict(auc_gam_brt, newdata = newdata, 
                              type = "response")
preds_brt$method <- "brt"

# bind all preditions into a single data frame
preds_auc_gam <- bind_rows(preds_glm, preds_idwi, preds_brt)

# add simpson evenness values to predictions
temp_evenness <- evenness[evenness$taxon %in% c("insect - moth", 
                                                "insect - butterfly", 
                                                "bird"), 
                          colnames(evenness) %in% c("taxon", 
                                                    "simpson_evenness")]
temp_evenness[nrow(temp_evenness) + 1, ] <- c("no bias", 1)
temp_evenness$taxon[temp_evenness$taxon == "insect - moth"] <- "moth"
temp_evenness$taxon[temp_evenness$taxon == "insect - butterfly"] <- "butterfly"
preds_auc_gam <- left_join(preds_auc_gam, temp_evenness, 
                           by = c("bias.name" = "taxon"))


## construct adjusted R2 variable importance table
var_imp_r2 <- data.frame(variable = c("full_model", "-sample_size", 
                                      "-spatial_bias"), 
                         R2_glm = NA, R2_idwi = NA, R2_brt = NA)
# fill in values for full models (all methods)
var_imp_r2[1, 2:ncol(var_imp_r2)] <- round(c(summary(auc_glm)$r.sq, 
                                       summary(auc_idwi)$r.sq, 
                                       summary(auc_gam_brt)$r.sq), digits = 2)

# fit models without sample size
auc_glm_noSampSize <- gam(auc ~ bias.name, 
                          data = ev_df_cv[!is.na(ev_df_cv$auc) & 
                                            ev_df_cv$Method == "glm_poly_wg_block_cv", ], 
                          family = betar(link = "logit"), 
                          select = TRUE)
auc_idwi_noSampSize <- gam(auc ~ bias.name, 
                          data = ev_df_cv[!is.na(ev_df_cv$auc) & 
                                            ev_df_cv$Method == "idw_interp_wg_block_cv", ], 
                          family = betar(link = "logit"), 
                          select = TRUE)
auc_brt_noSampSize <- gam(auc ~ bias.name, 
                          data = ev_df_cv[!is.na(ev_df_cv$auc) & 
                                            ev_df_cv$Method == "brt_wg_block_cv", ], 
                          family = betar(link = "logit"), 
                          select = TRUE)
# fit models without bias
auc_glm_noBias <- gam(auc ~ s(n_obs_per_sp, k = k, bs = "cr"), 
                      data = ev_df_cv[!is.na(ev_df_cv$auc) & 
                                        ev_df_cv$Method == "glm_poly_wg_block_cv", ], 
                      family = betar(link = "logit"), 
                      select = TRUE)
auc_idwi_noBias <- gam(auc ~ s(n_obs_per_sp, k = k, bs = "cr"), 
                      data = ev_df_cv[!is.na(ev_df_cv$auc) & 
                                        ev_df_cv$Method == "idw_interp_wg_block_cv", ], 
                      family = betar(link = "logit"), 
                      select = TRUE)
auc_brt_noBias <- gam(auc ~ s(n_obs_per_sp, k = k, bs = "cr"), 
                      data = ev_df_cv[!is.na(ev_df_cv$auc) & 
                                        ev_df_cv$Method == "brt_wg_block_cv", ], 
                      family = betar(link = "logit"), 
                      select = TRUE)

## fit null models
auc_glm_null <- gam(auc ~ 1, 
                 data = ev_df_cv[!is.na(ev_df_cv$auc) & 
                                   ev_df_cv$Method == "glm_poly_wg_block_cv", ], 
                 family = betar(link = "logit"), 
                 select = TRUE)
auc_idwi_null <- gam(auc ~ 1, 
                     data = ev_df_cv[!is.na(ev_df_cv$auc) & 
                                       ev_df_cv$Method == "idw_interp_wg_block_cv", ], 
                     family = betar(link = "logit"), 
                     select = TRUE)
auc_brt_null <- gam(auc ~ 1, 
                    data = ev_df_cv[!is.na(ev_df_cv$auc) & 
                                      ev_df_cv$Method == "brt_wg_block_cv", ], 
                    family = betar(link = "logit"), 
                    select = TRUE)

# calculate proportion of deviance explained by smooth of sample size
# See Dina's email of 26 Nov 2019.  (Null_deviance - Reduced_deviance) *100/Null_deviance

prop_dev_nobs_glm <- (deviance(auc_glm_null) - deviance(auc_glm_noBias)) * 
  100/deviance(auc_glm_null)
# prop_dev_nobs_rf
# prop_dev_nobs_idwi
# prop_dev_nobs_brt
# # calculate proportion of deviance explained by spatial bias variable
# prop_dev_bias_glm
# prop_dev_bias_rf
# prop_dev_bias_idwi
# prop_dev_bias_brt



# fill in adjusted R^2 values for models without terms
var_imp_r2[var_imp_r2$variable == "-sample_size", 2:ncol(var_imp_r2)] <- 
  round(c(summary(auc_glm_noSampSize)$r.sq, 
          summary(auc_idwi_noSampSize)$r.sq, 
          summary(auc_brt_noSampSize)$r.sq), digits = 2)
var_imp_r2[var_imp_r2$variable == "-spatial_bias", 2:ncol(var_imp_r2)] <- 
  round(c(summary(auc_glm_noBias)$r.sq, 
          summary(auc_idwi_noBias)$r.sq, 
          summary(auc_brt_noBias)$r.sq), digits = 2)
var_imp_r2 # removing the sample size term results in much less variance explained.
rm(auc_glm_noSampSize, auc_idwi_noSampSize, auc_brt_noSampSize, 
   auc_glm_noBias, auc_idwi_noBias, auc_brt_noBias)

###########################################################
### model RMSE
## fit only glm -----------------------------------------------------
# rmse_glm <- gam(rmse ~ factor(bias.name, 
#                               levels = c("no_bias", "least_birdNBDC", 
#                                          "moderate_bryBBS", "butterflyNBDC", 
#                                          "extreme_mothNBDC"), 
#                               labels = c("no bias", "bird", "bry", 
#                                          "butterfly", "moth")) + 
#                   s(n_obs_per_sp, 
#                     by = factor(bias.name, 
#                                 levels = c("no_bias", "least_birdNBDC", 
#                                            "moderate_bryBBS", "butterflyNBDC", 
#                                            "extreme_mothNBDC"), 
#                                 labels = c("no bias", "bird", "bry", 
#                                            "butterfly", "moth")), 
#                     k = k, bs = "cr") + 
#                   s(prevalence, k = k, bs = "cr"), 
#                 data = ev_df_cv[!is.na(ev_df_cv$auc) & 
#                                   ev_df_cv$Method == "glm_poly_wg_block_cv", ], 
#                 family = gaussian(link = "identity"))
# rmse_glm
# gam.check(rmse_glm)
# #qq.gam(rmse_glm)
# plot(rmse_glm, pages = 1, all.terms = T)
# plot(rmse_glm, pages = 1, residuals = T)

## fit only glm -----------------------------------------------------
# thin plate spline does not cause big changes in coefficient estimates, edf, or significance of terms
rmse_glm <- gam(rmse ~  s(prevalence, k = k, bs = "cr") + bias.name + 
                 s(n_obs_per_sp, 
                   by = bias.name, 
                   k = k, bs = "cr"), 
               data = ev_df_cv[!is.na(ev_df_cv$rmse) & 
                                 ev_df_cv$Method == "glm_poly_wg_block_cv", ], 
               family = betar(link = "logit"))
rmse_glm
gam.check(rmse_glm) # edf and k-index look bad-ish, but increasing k gives smooths that wiggle too much at medium-high prevalences.
plot(rmse_glm, pages = 1, all.terms = T)

## fit only inverse distance-weighted interpolation ---------------------------
# tp spline causes no major changes in model
rmse_idwi <- gam(rmse ~  s(prevalence, k = k, bs = "cr") + bias.name + 
                 s(n_obs_per_sp, by = bias.name, k = k, bs = "cr"), 
               data = ev_df_cv[!is.na(ev_df_cv$rmse) & 
                                 ev_df_cv$Method == "idw_interp_wg_block_cv", ], 
               family = betar(link = "logit"))
rmse_idwi
gam.check(rmse_idwi) # a few big residuals.  K issue is the same, so ok.  Removing residuals with abs() > 5 does not change model in any meaningful way.  
plot(rmse_idwi, pages = 1, all.terms = T)

## fit only boosted classification trees ---------------------------------------
rmse_brt <- gam(rmse ~  s(prevalence, k = k, bs = "cr") + bias.name + 
                   s(n_obs_per_sp, by = bias.name, k = k, bs = "cr"), 
                 data = ev_df_cv[!is.na(ev_df_cv$rmse) & 
                                   ev_df_cv$Method == "brt_wg_block_cv", ], 
                 family = betar(link = "logit"))
rmse_brt
gam.check(rmse_brt) # some big (negative) residuals.  K issue maybe better, ok.  Removing residuals > 5 does not change model meaningfully except making smooth of n_obs by butterfly and moth biases not significant.  
plot(rmse_brt, pages = 1, all.terms = T)


### generate RMSE predictions to use in plotting -----------------------------
# get median prevalence 
med_prev <- median(ev_df_cv$prevalence[ev_df_cv$Method == "glm_poly_wg_block_cv" & 
                                         ev_df_cv$bias.name == "no bias" & 
                                         ev_df_cv$n_obs_per_sp == 200], 
                   na.rm = T)
newdata <- data.frame(n_obs_per_sp = rep(nobs_range, 4), 
                      bias.name = c(rep("no bias", length(nobs_range)), 
                                    rep("bird", length(nobs_range)), 
                                    rep("butterfly", length(nobs_range)), 
                                    rep("moth", length(nobs_range))), 
                      prevalence = med_prev)
# predictions when using a GLM SDM
preds_rmse_gam_glm <- newdata
preds_rmse_gam_glm$pred_rmse <- predict(rmse_glm, newdata = newdata, 
                              type = "response")
preds_rmse_gam_glm$method <- "glm"

# predictions when using an inverse distance-weighted interpolation SDM
preds_rmse_gam_idwi <- newdata
preds_rmse_gam_idwi$pred_rmse <- predict(rmse_idwi, newdata = newdata, 
                               type = "response")
preds_rmse_gam_idwi$method <- "idwi"

# predictions when using a boosted classification tree SDM
preds_rmse_gam_brt <- newdata
preds_rmse_gam_brt$pred_rmse <- predict(rmse_brt, newdata = newdata, 
                              type = "response")
preds_rmse_gam_brt$method <- "brt"

# bind all preditions into a single data frame
preds_rmse_gam <- bind_rows(preds_rmse_gam_glm, 
                            preds_rmse_gam_idwi, preds_rmse_gam_brt)










################################################################################
### Stats and numbers for text of paper
################################################################################
# AUC GLM no bias v. butterfly, n.obs/sp = 158
preds_auc_gam$pred_auc[preds_auc_gam$method == "glm" & 
                     preds_auc_gam$n_obs_per_sp == 158 & 
                     preds_auc_gam$bias.name == "no bias"] - 
  preds_auc_gam$pred_auc[preds_auc_gam$method == "glm" & 
                       preds_auc_gam$n_obs_per_sp == 158 & 
                       preds_auc_gam$bias.name == "butterfly"] # 0.049
# auc w/ no bias, glm, max nobs
preds_auc_gam$pred_auc[preds_auc_gam$method == "glm" & 
                     preds_auc_gam$n_obs_per_sp == 158 & 
                     preds_auc_gam$bias.name == "no bias"] # 0.736
# auc w/ butterfly bias, glm, max nobs
preds_auc_gam$pred_auc[preds_auc_gam$method == "glm" & 
                     preds_auc_gam$n_obs_per_sp == 158 & 
                     preds_auc_gam$bias.name == "butterfly"] # 0.6868

# AUC BRT no bias v. butterfly, n.obs/sp = 158
preds_auc_gam$pred_auc[preds_auc_gam$method == "brt" & 
                     preds_auc_gam$n_obs_per_sp == 158 & 
                     preds_auc_gam$bias.name == "no bias"] - 
  preds_auc_gam$pred_auc[preds_auc_gam$method == "brt" & 
                       preds_auc_gam$n_obs_per_sp == 158 & 
                       preds_auc_gam$bias.name == "butterfly"] # 0.022

# AUC idwi no bias v. butterfly, n.obs/sp = 158
preds_auc_gam$pred_auc[preds_auc_gam$method == "idwi" & 
                     preds_auc_gam$n_obs_per_sp == 158 & 
                     preds_auc_gam$bias.name == "no bias"] - 
  preds_auc_gam$pred_auc[preds_auc_gam$method == "idwi" & 
                       preds_auc_gam$n_obs_per_sp == 158 & 
                       preds_auc_gam$bias.name == "butterfly"] # 0.026


## AUC median vs no bias, sample size vs. spatial bias
preds_auc_gam$pred_auc[preds_auc_gam$method == "glm" & 
                     preds_auc_gam$n_obs_per_sp == 158 & 
                     preds_auc_gam$bias.name == "no bias"] - 
  preds_auc_gam$pred_auc[preds_auc_gam$method == "glm" & 
                       preds_auc_gam$n_obs_per_sp == 158 & 
                       preds_auc_gam$bias.name == "butterfly"] # 0.049

# samp size reduction to get AUC decrease equivalent to switching to 
# butterfly bias (decrease of 0.04902 (above))
target <- preds_auc_gam$pred_auc[preds_auc_gam$method == "glm" & 
                               preds_auc_gam$n_obs_per_sp == 158 & 
                               preds_auc_gam$bias.name == "butterfly"]
preds_auc_gam[preds_auc_gam$method == "glm" & preds_auc_gam$bias.name == "no bias" & 
            abs(preds_auc_gam$pred_auc[preds_auc_gam$method == "glm" & 
                                 preds_auc_gam$bias.name == "no bias"] - target) == 
            min(abs(preds_auc_gam$pred_auc[preds_auc_gam$method == "glm" & 
                  preds_auc_gam$bias.name == "no bias"] - target)), ] # nobs/sp = 33

## delta AUC w/nobs_per_sp = 10 and spat bias changing from median to none
preds_auc_gam$pred_auc[preds_auc_gam$method == "glm" & 
                     preds_auc_gam$n_obs_per_sp == 10 & 
                     preds_auc_gam$bias.name == "no bias"] - 
  preds_auc_gam$pred_auc[preds_auc_gam$method == "glm" & 
                       preds_auc_gam$n_obs_per_sp == 10 & 
                       preds_auc_gam$bias.name == "butterfly"] # 0.017
# increase in sample size to get AUC improvement equivalent to switching from
# butterfly to no bias (above)
target <- preds_auc_gam$pred_auc[preds_auc_gam$method == "glm" & 
                               preds_auc_gam$n_obs_per_sp == 10 & 
                               preds_auc_gam$bias.name == "no bias"]
preds_auc_gam[preds_auc_gam$method == "glm" & preds_auc_gam$bias.name == "butterfly" & 
            abs(preds_auc_gam$pred_auc[preds_auc_gam$method == "glm" & 
                                     preds_auc_gam$bias.name == "butterfly"] - target) == 
            min(abs(preds_auc_gam$pred_auc[preds_auc_gam$method == "glm" & 
                                         preds_auc_gam$bias.name == "butterfly"] - target)), ] # n_obs_per_sp = 15


## RMSE for GLM severe bias - RMSE for GLM no bias
preds_rmse_gam$pred_rmse[preds_rmse_gam$method == "glm" & 
                       preds_rmse_gam$n_obs_per_sp == 158 & 
                       preds_rmse_gam$bias.name == "moth"] - 
  preds_rmse_gam$pred_rmse[preds_rmse_gam$method == "glm" & 
                         preds_rmse_gam$n_obs_per_sp == 158 & 
                         preds_rmse_gam$bias.name == "no bias"] # 0.008
# GLMs trained with severely spatially biased training data did have higher RMSE than GLMs trained with less spatially biased data, but the size of this effect was small (difference in predicted RMSE between GLMs trained with severely biased data and GLMs trained with unbiased data was 0.0083 when using the largest sample size and fixing species prevalence at the median value). 

## evennes metrics
evenness$simpson_evenness[evenness$taxon == "insect - moth"] # 0.021
evenness$simpson_evenness[evenness$taxon == "insect - butterfly"] # 0.126
evenness$simpson_evenness[evenness$taxon == "bird"] # 0.762
simpson_even(bias_rasters$no_bias$layer[complete.cases(
  bias_rasters$no_bias$layer)]) # 1 (no bias)

# quantiles of evenness values
even_ecdf <- ecdf(evenness$simpson_evenness)
even_ecdf(evenness$simpson_evenness[evenness$taxon == "insect - moth"]) # 0.06
even_ecdf(evenness$simpson_evenness[evenness$taxon == "insect - butterfly"]) # 0.33
even_ecdf(evenness$simpson_evenness[evenness$taxon == "bird"]) # 1

