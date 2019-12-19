############################################
## This analyzes prediction performance of models from the simulation
## using boosted regression trees
## This is for simulation 10.10 
## This script is meant to be sourced from within main.R
##
## inputs:  workspace from main.R
## 
## author: Willson Gaul wgaul@hotmail.com
## created: 4 Oct 2019
## last modified: 6 Nov 2019
############################################

# no need to set working directory or load packages as that is done in main.R
# library(dismo)

set.seed(04102019) # date I started this script

### load data -----------------------------------------------------------------
# optionally, read in a .csv with auc results from simulation
ev_df <- read_csv("./bryophyte_template_evals_24Sep2019.csv")
df_butt <- read_csv("../sims_10.10_bryophyte_butterflyBias/bryophyte_template_butterfly_bias_evals_12Sep2019.csv")
ev_df <- bind_rows(ev_df, df_butt)

# keep only relevant columns
ev_df <- ev_df[, colnames(ev_df) %in% c("n.obs", "bias.name", "Method", "Draw", 
                                        "auc", "rmse", "time", 
                                        "prevalence", "n_obs_per_sp")]

# subset to only block CV results
ev_df_cv <- ev_df[which(ev_df$Method %in% c("glm_poly_wg_block_cv", 
                                            "random_forest_wg_block_cv", 
                                            "brt_wg_block_cv", 
                                            "idw_interp_wg_block_cv")), ]

ev_df_cv$bias.name <- factor(ev_df_cv$bias.name, 
                             levels = c("no_bias", "least_birdNBDC", 
                                        "moderate_bryBBS", "butterflyNBDC", 
                                        "extreme_mothNBDC"), 
                             labels = c("no bias", "bird", "bry", 
                                        "butterfly", "moth"))
ev_df_cv$Method <- factor(ev_df_cv$Method, 
                          levels = c("glm_poly_wg_block_cv", "brt_wg_block_cv",  
                          "idw_interp_wg_block_cv", "random_forest_wg_block_cv"), 
                          labels = c("glm", "brt", "idwi", "rf"))
### end load data -------------------------------------------------------------

##### model AUC -------------------------------------------------------------
# fit a BRT to model AUC as a function of overall sample size (n_obs_per_sp),
# spatial sampling bias (bias.name), and modeling method (Method)
# 
# from looking at boxplots, I expect AUC to vary by Method, nobs, Method*nobs, 
# bias, bias*Method (e.g. bias matters more for brt than idwi), and perhaps
# bias*nobs*method (e.g. RF high nobs matters at low bias but not high bias)

auc_brt <- gbm.step(data = data.frame(ev_df_cv[complete.cases(ev_df_cv), ]), 
                    gbm.x = which(colnames(ev_df_cv) %in% 
                                    c("bias.name", "Method", "n_obs_per_sp")), 
                    gbm.y = which(colnames(ev_df_cv) == "auc"), 
                    tree.complexity = 5, learning.rate = 0.01, 
                    family = "gaussian", max.trees = 10000, plot.main = T)

ints <- gbm.interactions(auc_brt)

auc_brt
plot(auc_brt)
auc_summary <- summary(auc_brt)
gbm.plot(auc_brt)
ints$interactions # look at pairwise interactions

### generate AUC predictions to use in plotting
nobs_range <- 1:158 # set up range of n_obs_per_sp values

newdata <- data.frame(expand.grid(nobs_range, 
                                  unique(as.character(ev_df_cv$bias.name)), 
                                  unique(as.character(ev_df_cv$Method))))
colnames(newdata) <- c("n_obs_per_sp", "bias.name", "Method")
newdata$bias.name <- factor(newdata$bias.name, 
                            levels = c("no bias", "bird", "bry", 
                                       "butterfly", "moth"), 
                            labels = c("no bias", "bird", "bry", 
                                       "butterfly", "moth"))
newdata$Method <- factor(newdata$Method, 
                         levels = c("glm", "brt", "idwi", "rf"), 
                         labels = c("glm", "brt", "idwi", "rf"))

preds_auc_brt <- newdata
preds_auc_brt$pred_auc <- predict(auc_brt, newdata = newdata, 
                              n.trees = auc_brt$gbm.call$best.trees, 
                              type = "response")
colnames(preds_auc_brt)[colnames(preds_auc_brt) == "Method"] <- "method" # to match plotting code

# add simpson evenness values to predictions
temp_evenness <- evenness[evenness$taxon %in% c("insect - moth", 
                                                "insect - butterfly", 
                                                "bry_BBS", "bird"), 
                          colnames(evenness) %in% c("taxon", "simpson_evenness")]
temp_evenness[nrow(temp_evenness) + 1, ] <- c("no bias", 1)
temp_evenness$taxon[temp_evenness$taxon == "insect - moth"] <- "moth"
temp_evenness$taxon[temp_evenness$taxon == "insect - butterfly"] <- "butterfly"
temp_evenness$taxon[temp_evenness$taxon == "bry_BBS"] <- "bry"
preds_auc_brt <- left_join(preds_auc_brt, temp_evenness, by = c("bias.name" = "taxon"))
### end model AUC -------------------------------------------------------------

##### model RMSE -------------------------------------------------------------
rmse_evaluation_brt <- gbm.step(data = data.frame(ev_df_cv[complete.cases(ev_df_cv), ]), 
                    gbm.x = which(colnames(ev_df_cv) %in% 
                                    c("bias.name", "Method", "n_obs_per_sp", 
                                      "prevalence")), 
                    gbm.y = which(colnames(ev_df_cv) == "rmse"), 
                    tree.complexity = 5, learning.rate = 0.01, 
                    family = "gaussian", max.trees = 10000, plot.main = T)

rmse_ints <- gbm.interactions(rmse_evaluation_brt)

rmse_evaluation_brt
plot(rmse_evaluation_brt)
rmse_summary <- summary(rmse_evaluation_brt)
gbm.plot(rmse_evaluation_brt)
rmse_ints$interactions # look at pairwise interactions

### generate AUC predictions to use in plotting
nobs_range <- 1:158 # set up range of n_obs_per_sp values
# get median prevalence 
med_prev <- median(ev_df_cv$prevalence[ev_df_cv$Method == "glm" & 
                                         ev_df_cv$bias.name == "no bias" & 
                                         ev_df_cv$n.obs == 200000], na.rm = T)

newdata <- data.frame(expand.grid(nobs_range, 
                                  unique(as.character(ev_df_cv$bias.name)), 
                                  unique(as.character(ev_df_cv$Method))))
colnames(newdata) <- c("n_obs_per_sp", "bias.name", "Method")
newdata$bias.name <- factor(newdata$bias.name, 
                            levels = c("no bias", "bird", "bry", 
                                       "butterfly", "moth"), 
                            labels = c("no bias", "bird", "bry", 
                                       "butterfly", "moth"))
newdata$Method <- factor(newdata$Method, 
                         levels = c("glm", "brt", "idwi", "rf"), 
                         labels = c("glm", "brt", "idwi", "rf"))
newdata$prevalence <- med_prev

preds_rmse_evaluation_brt <- newdata
preds_rmse_evaluation_brt$pred_rmse <- predict(rmse_evaluation_brt, newdata = newdata, 
                              n.trees = rmse_evaluation_brt$gbm.call$best.trees, 
                              type = "response")
colnames(preds_rmse_evaluation_brt)[colnames(preds_rmse_evaluation_brt) == "Method"] <- "method" # to match plotting code
### end model RMSE -------------------------------------------------------------



################################################################################
### Stats and numbers for text of paper
################################################################################
# expected AUC for GLMs with max n.obs and no bias
preds_auc_brt$pred_auc[preds_auc_brt$n_obs_per_sp == max(preds_auc_brt$n_obs_per_sp) & 
                preds_auc_brt$bias.name == "no bias" & 
                preds_auc_brt$method == "glm"]
# expected AUC for GLMs with max n.obs and median bias
preds_auc_brt$pred_auc[preds_auc_brt$n_obs_per_sp == max(preds_auc_brt$n_obs_per_sp) & 
                         preds_auc_brt$bias.name == "butterfly" & 
                         preds_auc_brt$method == "glm"]
# difference in AUC for GLMs between no and median bias (with max n.obs)
preds_auc_brt$pred_auc[preds_auc_brt$n_obs_per_sp == max(preds_auc_brt$n_obs_per_sp) & 
                         preds_auc_brt$bias.name == "no bias" & 
                         preds_auc_brt$method == "glm"] - 
  preds_auc_brt$pred_auc[preds_auc_brt$n_obs_per_sp == max(preds_auc_brt$n_obs_per_sp) & 
                           preds_auc_brt$bias.name == "butterfly" & 
                           preds_auc_brt$method == "glm"]
# difference in AUC for BRTs between no and median bias (with max n.obs)
preds_auc_brt$pred_auc[preds_auc_brt$n_obs_per_sp == max(preds_auc_brt$n_obs_per_sp) & 
                         preds_auc_brt$bias.name == "no bias" & 
                         preds_auc_brt$method == "brt"] - 
  preds_auc_brt$pred_auc[preds_auc_brt$n_obs_per_sp == max(preds_auc_brt$n_obs_per_sp) & 
                           preds_auc_brt$bias.name == "butterfly" & 
                           preds_auc_brt$method == "brt"]
# difference in AUC for inverse distance-weighted interpolation between no and median bias (with max n.obs)
preds_auc_brt$pred_auc[preds_auc_brt$n_obs_per_sp == max(preds_auc_brt$n_obs_per_sp) & 
                         preds_auc_brt$bias.name == "no bias" & 
                         preds_auc_brt$method == "idwi"] - 
  preds_auc_brt$pred_auc[preds_auc_brt$n_obs_per_sp == max(preds_auc_brt$n_obs_per_sp) & 
                           preds_auc_brt$bias.name == "butterfly" & 
                           preds_auc_brt$method == "idwi"]


# variable importance values
auc_summary

