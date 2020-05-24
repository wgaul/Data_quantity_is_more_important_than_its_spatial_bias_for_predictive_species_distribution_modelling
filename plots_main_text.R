############################################
## This generates plots for the main text of the simulation paper.
## This is for simulation 10.11 using the bryophyte template (large community)
## 
## This script is meant to be run after sourcing analyze_performance.R
##
## inputs:  * bryophyte_template_evals_9Sep2019.csv 
##          * bryophyte_template_butterfly_bias_evals_9Sep2019.csv
##
## 
## author: Willson Gaul wgaul@hotmail.com
## created: 10 Sep 2019 based on code from plots.R
## last modified: 27 April 2020
############################################

library(wgutil)
library(Hmisc)
library(rgdal)
library(raster)
library(parallel)
library(pROC)
library(tidyverse)
library(plotly)

simulation_name <- "Bryophyte Simulation"

bias_names <- c("no_bias", "extreme_mothNBDC", 
                "median_butterflyNBDC", "least_birdNBDC")

## optionally, read in a .csv with auc results from simulation
## Large community simulation results
ev_df <- read_csv("./bryophyte_template_evals_10April2020.csv")

# subset to only block CV results
ev_df_cv <- ev_df[which(ev_df$Method %in% c("glm_poly_wg_block_cv", 
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

## load truth maps
truth_maps <- readRDS("truth_maps.rds")
## load example observations
example_obs <- readRDS("example_obs_n63400.rds")
## load example block CV fold assignments
example_folds <- readRDS("example_folds.rds")

t_size <- 19
line_size <- 1.5
c_f = 1 # coord fixed for plots


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

### plot sampling biases ------------------------------------------------------
# load sample bias rasters
spat_nrec_hec <- readRDS("spat_nrec_hec.rds")

bias_df <- bias_rasters
# set value for no bias map to 1 in each cell. 
bias_df$no_bias$layer[!is.na(bias_df$no_bias$layer)] <- 1 # 0.05
bias_df <- mapply(bias_df, names(bias_df), SIMPLIFY = F, USE.NAMES = T,
                  FUN = function(d, nm) {
                    d$bias <- rep(nm, nrow(d))
                    d})
bias_df <- bind_rows(bias_df)
# transform for better visualization?
# bias_df$layer[!is.na(bias_df$layer)] <- sqrt(bias_df$layer[!is.na(bias_df$layer)])
# bias_df$layer[!is.na(bias_df$layer)] <- log(bias_df$layer[!is.na(bias_df$layer)] + 0.00001)

bias_df$bias <- factor(bias_df$bias, levels = c("no_bias", 
                                                "least_birdNBDC", 
                                                "median_butterflyNBDC", 
                                                "extreme_mothNBDC"), 
                       labels = c("A", "B", "C", "D"))

annot <- data.frame(x1 = c(265000, 310000, 60000, 60000), 
                    x2 = c(365000, 310000, 60000, 60000), 
                    y1 = c(60000, 40000, 400000, 380000), 
                    y2 = c(60000, 40000, 455000, 380000),
                    label = c(NA, "100 km", NA, "N"), 
                    bias = "D")

sampling_bias_maps <- ggplot() + 
  geom_raster(data = bias_df[complete.cases(bias_df), ], 
              aes(x = x, y = y, fill = layer)) + 
  coord_fixed(c_f) + 
  # scale_fill_gradient(name = "Relative\nRecording\nEffort\n", 
  #                     low = "light grey", high = "black") + 
  scale_fill_gradient(name = "Relative\nRecording\nEffort\n", 
                      low = "light grey", high = "black") + 
  facet_wrap(~bias) +
  geom_segment(data = annot[1, ], aes(x = x1, xend = x2, y = y1, yend = y2)) + 
  geom_text(data = annot[c(2, 4), ], aes(x = x1, y = y1, label = label)) + 
  geom_segment(data = annot[3, ], aes(x = x1, xend = x2, y = y1, yend = y2), 
               arrow = arrow(length = unit(0.1, "npc"))) + 
  theme_bw() + 
  theme(text = element_text(size = t_size), 
        axis.title = element_blank(),
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        strip.text = element_text(hjust = -0.01))
sampling_bias_maps
### end plot sampling biases --------------------------------------------------


### plot example species distributions -----------------------------------------
# plot(subset(truth_maps, c(98, 14, 32, 47, 63, 100, 102)))
# sp 1
sp_rast_1 <- data.frame(rasterToPoints(subset(truth_maps, 98)))
colnames(sp_rast_1)[3] <- "present"

p_dist_1 <- ggplot(data = sp_rast_1, aes(x = x, y = y), cex = 2) + 
  geom_raster(aes(fill = factor(as.character(present), 
                                levels = c(0, 1), 
                                labels = c("Absent", "Present")))) + 
  coord_fixed(c_f) + 
  scale_fill_manual(name = element_blank(),
                    values = c("light grey", "dark green")) +
  ggtitle("A") + 
  theme_bw() + 
  theme(text = element_text(size = t_size), 
        axis.title = element_blank(),
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        legend.position = "none")

# sp 2
sp_rast_2 <- data.frame(rasterToPoints(subset(truth_maps, 32)))
colnames(sp_rast_2)[3] <- "present"

p_dist_2 <- ggplot(data = sp_rast_2, aes(x = x, y = y), cex = 2) + 
  geom_raster(aes(fill = factor(as.character(present), 
                                levels = c(0, 1), 
                                labels = c("Absent", "Present")))) + 
  coord_fixed(c_f) + 
  scale_fill_manual(name = element_blank(),
                    values = c("light grey", "dark green")) +
  ggtitle("B") + 
  theme_bw() + 
  theme(text = element_text(size = t_size), 
        axis.title = element_blank(),
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        legend.position = "none")

# sp 3
sp_rast_3 <- data.frame(rasterToPoints(subset(truth_maps, 47)))
colnames(sp_rast_3)[3] <- "present"

p_dist_3 <- ggplot(data = sp_rast_3, aes(x = x, y = y), cex = 2) + 
  geom_raster(aes(fill = factor(as.character(present), 
                                levels = c(0, 1), 
                                labels = c("Absent", "Present")))) + 
  coord_fixed(c_f) + 
  scale_fill_manual(name = element_blank(),
                    values = c("light grey", "dark green")) +
  ggtitle("C") + 
  theme_bw() + 
  theme(text = element_text(size = t_size), 
        axis.title = element_blank(),
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        legend.position = "none")

# sp 4
sp_rast_4 <- data.frame(rasterToPoints(subset(truth_maps, 63)))
colnames(sp_rast_4)[3] <- "present"

p_dist_4 <- ggplot(data = sp_rast_4, aes(x = x, y = y), cex = 2) + 
  geom_raster(aes(fill = factor(as.character(present), 
                                levels = c(0, 1), 
                                labels = c("Absent", "Present")))) + 
  coord_fixed(c_f) + 
  scale_fill_manual(name = element_blank(),
                    values = c("light grey", "dark green")) +
  ggtitle("D") + 
  theme_bw() + 
  theme(text = element_text(size = t_size), 
        axis.title = element_blank(),
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        legend.position = "none")
### end plot example distributions --------------------------------------------


### plot number of models fit -------------------------------------------------
# count number of models producing RMSE measures
n_rmse_fits <- ev_df_cv[!is.na(ev_df_cv$rmse), ]
n_rmse_fits <- data.frame(table(n_rmse_fits$bias.name, n_rmse_fits$n.obs, 
                                n_rmse_fits$Method))
colnames(n_rmse_fits) <- c("bias", "n_obs", "method", "n_fits")
print(n_rmse_fits)

n_rmse_point <- ggplot(data = n_rmse_fits, 
                       aes(x = as.numeric(as.character(n_obs))/1267, 
                           y = n_fits, 
                           shape = factor(
                             bias, 
                             levels = c("no_bias", 
                                        "least_birdNBDC", 
                                        "median_butterflyNBDC", 
                                        "extreme_mothNBDC"), 
                             labels = c("\nno bias\n", 
                                        "\nlow\n",
                                        "\nmedian\n",  
                                        "\nsevere\n")),
                           color = factor(
                             bias, 
                             levels = c("no_bias", 
                                        "least_birdNBDC", 
                                        "median_butterflyNBDC", 
                                        "extreme_mothNBDC"), 
                             labels = c("\nno bias\n", 
                                        "\nlow\n",
                                        "\nmedian\n",  
                                        "\nsevere\n")))) + 
  geom_line(size = line_size) + 
  facet_wrap(~ factor(method, 
                      levels = c("glm_poly_wg_block_cv", 
                                 "brt_wg_block_cv", 
                                 "idw_interp_wg_block_cv"), 
                      labels = c("A", "B", "C"))) + 
  # ggtitle("Number of species for which models were sucessfully fitted") + 
  ylab("Number of species fit") + 
  xlab("Average number of records\nper species in dataset") + 
  scale_shape_discrete(name = "Spatial Bias", solid = F) +  
  scale_colour_viridis_d(option = "magma", name = "Spatial Bias", 
                         begin = 0, end = 0.65) +
  scale_y_continuous(breaks = seq(10, 110, by = 20)) + 
  theme_bw() + 
  theme(text = element_text(size = t_size*1.2), 
        legend.key.width = unit(2*t_size, "points"), 
        axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1))
n_rmse_point
### end plot number of models fit ---------------------------------------------


###################################
## AUC plots 
## Boosted regression tree model of AUC by bias

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

auc_line <- ggplot(data = preds_auc_brt, 
                     aes(x = n_obs_per_sp, y = pred_auc, 
                         group = factor(bias.name,
                                        levels = c("no bias", "bird",  
                                                   "butterfly", "moth"),
                                        labels = c(
                                          "\nnone\n", "\nlow\n", 
                                          "\nmedian\n", 
                                          "\nsevere\n"),
                                        ordered = TRUE))) + 
  geom_line(aes(color = factor(bias.name,
                               levels = c("no bias", "bird",
                                          "butterfly", "moth"),
                               labels = c(
                                 "\nnone\n", "\nlow\n", 
                                 "\nmedian\n", "\nsevere\n"),
                               ordered = TRUE)), 
            size = line_size) + 
  # ggtitle("Large Community Simulation") + #simulation_name
  geom_rug(data = rug_df, aes(x = n_obs_per_sp), 
           show.legend = F, sides = "b") +
  facet_wrap( ~ factor(method, 
                       levels = c("glm", "brt", "idwi"), 
                       labels = c("A", "B", "C"))) + 
  xlab("Average number of records \nper species in data set") + 
  ylab("AUC\n(block cross-validated)") +
  # scale_linetype_discrete(name = "Spatial\nSampling\nBias") + 
  scale_color_viridis_d(name = "Spatial\nSampling\nBias", option = "magma", 
                        begin = 0, end = 0.65) +
  theme_bw() + 
  theme(text = element_text(size = t_size*1.2), 
        legend.key.width = unit(1.8*t_size, "points"), 
        axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1))
auc_line

auc_line_presentation <- ggplot(data = preds_auc_brt, 
                                aes(x = n_obs_per_sp, y = pred_auc, 
                                    group = factor(bias.name,
                                                   levels = c("no bias", 
                                                              "bird", 
                                                              "butterfly", 
                                                              "moth"),
                                                   labels = c(
                                                     "\nnone\n", "\nlow\n", 
                                                     "\nmedian\n", 
                                                     "\nsevere\n"),
                                                   ordered = TRUE))) + 
  geom_line(aes(color = factor(bias.name,
                               levels = c("no bias", "bird",
                                          "butterfly", "moth"),
                               labels = c(
                                 "\nnone\n", "\nlow\n", 
                                 "\nmedian\n", "\nsevere\n"),
                               ordered = TRUE)), 
  size = 1.2*line_size) + 
  # ggtitle("Large Community Simulation") + #simulation_name
  geom_rug(data = rug_df, aes(x = n_obs_per_sp), show.legend = F, sides = "b") +
  facet_wrap( ~ factor(method, 
                       levels = c("glm", "brt", "idwi"), 
                       labels = c("GLM",  
                                  "Boosted\nregression\ntree", 
                                  "Inverse\ndistance-\nweighted\ninterpolation"))) + 
  xlab("Average number of\nrecords per species") + 
  ylab("AUC\n(block cross-validated)") +
  # scale_linetype_discrete(name = "Spatial\nSampling\nBias") + 
  scale_color_viridis_d(name = "Spatial\nSampling\nBias", option = "magma", 
                        begin = 0, end = 0.65, direction = -1) +
  theme_bw() + 
  theme(text = element_text(size = t_size*1.5), 
        legend.key.width = unit(1.8*t_size, "points"), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
auc_line_presentation

# AUC boxplots
auc_boxplot <- ggplot(data = ev_df_cv, 
                      aes(x = factor(bias.name,
                                     levels = c("no_bias", "least_birdNBDC", 
                                                "median_butterflyNBDC", 
                                                "extreme_mothNBDC"),
                                     labels = c("none",
                                                "low",
                                                "median", 
                                                "severe"),
                                     ordered = TRUE), 
                          y = auc, 
                          fill = factor(n_obs_per_sp, 
                                        levels = c("2", "5", "10", 
                                                   "50", "100", "200"), 
                                        labels = c("2", "5", "10", 
                                                   "50", "100", "200")))) +
  geom_boxplot(varwidth = TRUE) +
  facet_wrap(~ factor(Method, 
                      levels = c("glm_poly_wg_block_cv", 
                                 "brt_wg_block_cv", 
                                 "idw_interp_wg_block_cv"), 
                      labels = c("A", "B", "C"))) + 
  # ggtitle("Large Community Simulation") + #simulation_name
  xlab("Spatial Sampling Bias") + 
  ylab("AUC\n(block cross-validated)") + 
  scale_fill_viridis_d(name = "Average\nnumber\nof\nrecords\nper\nspecies", 
                       option = "viridis", begin = 0, end = 0.85) + 
  theme_bw() + 
  theme(text = element_text(size = t_size*1.1), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
auc_boxplot
## end AUC boxplots -----------------------------------------------------------

## RMSE boxplots block CV results ---------------------------------------------
rmse_cv_boxplot <- ggplot(data = ev_df_cv,
                          aes(x = factor(bias.name,
                                         levels = c("no_bias",
                                                    "least_birdNBDC",
                                                    "median_butterflyNBDC",
                                                    "extreme_mothNBDC"),
                                         labels = c("none",
                                                    "low",
                                                    "median",
                                                    "severe"),
                                         ordered = TRUE),
                              y = rmse,
                              fill = factor(as.character(n_obs_per_sp),
                                            levels = c("2", "5", "10", "50", 
                                                       "100", "200")))) +
  geom_boxplot(varwidth = TRUE) +
  facet_wrap(~ factor(Method,
                      levels = c("glm_poly_wg_block_cv",
                                 "brt_wg_block_cv",
                                 "idw_interp_wg_block_cv"),
                      labels = c("A", "B", "C"))) +
  # ggtitle("Large Community Simulation") + #sim_name
  xlab("Sampling Bias") +
  ylab("RMSE\n(spatial block cross-validated)") +
  scale_fill_viridis_d(name = "Average\nnumber of\nrecords\nper species",
                       option = "viridis", begin = 0, end = 0.85) +
  theme_bw() +
  theme(text = element_text(size = t_size),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
rmse_cv_boxplot
## end RMSE boxplot ----------------------------------------------------------
##############################################################################



##################
### example training and test data
# using species r20.3, which I plotted a truth map of also.
# add jitter to coordinates
example_obs$no_bias$x_jit <- example_obs$no_bias$x + 
  runif(min = 500, max = 8000, 
        n = nrow(example_obs$no_bias))
example_obs$no_bias$y_jit <- example_obs$no_bias$y + 
  runif(min = 500, max = 8000, 
        n = nrow(example_obs$no_bias))
example_obs$median_bias$x_jit <- example_obs$median_bias$x + 
  runif(min = 500, max = 8000, 
        n = nrow(example_obs$median_bias))
example_obs$median_bias$y_jit <- example_obs$median_bias$y + 
  runif(min = 500, max = 8000, 
        n = nrow(example_obs$median_bias))

# add en column to data frames
sp_rast_1$en <- paste0(sp_rast_1$x, "_", sp_rast_1$y)
sp_rast_1$fold <- NA
sp_rast_1$fold[sp_rast_1$en %in% example_folds[[1]]$train_en] <- "training"
sp_rast_1$fold[sp_rast_1$en %in% example_folds[[1]]$test_en] <- "test"
example_obs <- lapply(example_obs, FUN = function(x) {
  x$en <- paste0(as.character(round(x$x)), "_", as.character(round(x$y)))
  x})

example_sp <- ggplot(data = sp_rast_1, aes(x = x, y = y), cex = 2) + 
  geom_raster(aes(fill = factor(as.character(present), 
                                levels = c(0, 1), 
                                labels = c("Absent", "Present")))) + 
  coord_fixed(c_f) + 
  scale_fill_manual(name = element_blank(), # "True\nvirtual\nspecies\nstatus"
                    values = c("light grey", "dark green")) +
  ggtitle("A") + 
  # geom_segment(data = annot[1, ], aes(x = x1, xend = x2, y = y1, yend = y2)) + 
  # geom_text(data = annot[c(2), ], aes(x = x1, y = y1, label = label)) +
  theme_bw() + 
  theme(text = element_text(size = t_size*3), 
        legend.position = "none", # legend.key.width = unit(1.8*t_size, "points")
        axis.ticks = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank()) 
example_sp

noBias_obs_map <- ggplot() + 
  geom_raster(data = sp_rast_1, 
              aes(x = x, y = y, fill = fold)) + 
  geom_point(
    data = example_obs$no_bias[example_obs$no_bias$list_length > 40 & 
                                       example_obs$no_bias$en %in%
                                       example_folds[[1]]$train_en, ], 
             aes(x = x_jit, y = y_jit, 
                 color = factor(as.character(sp3), 
                                levels = c("0", "1"), 
                                labels = c("No", "Yes"))), 
    size = t_size/4) + 
  coord_fixed(c_f) + 
  scale_color_discrete(name = "Species\nrecorded", h = c(110, 30), 
                       l = c(10, 65)) +
  scale_fill_manual(name = "Cross\nvalidation\nfold", 
                    values = c("dark grey", "light grey")) +
  # geom_segment(data = annot[1, ], aes(x = x1, xend = x2, y = y1, yend = y2)) + 
  # geom_text(data = annot[c(2), ], aes(x = x1, y = y1, label = label)) +
  # xlab("Eastings") + ylab("Northings") + 
  ggtitle("B") + 
  theme_bw() + 
  theme(text = element_text(size = t_size*3), 
        legend.position = "none",
        # legend.key.width = unit(1.8*t_size, "points"), 
        axis.ticks = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank())
noBias_obs_map

median_obs_map <- ggplot() + 
  geom_raster(data = sp_rast_1,
              aes(x = x, y = y, fill = fold)) +
  geom_point(
    data = example_obs$median_bias[
      example_obs$median_bias$list_length > 40 & 
        example_obs$median_bias$en %in%
        example_folds[[1]]$train_en, ], 
    aes(x = x_jit, y = y_jit, 
        color = factor(as.character(sp3), 
                       levels = c("0", "1"), 
                       labels = c("No", "Yes"))), 
    size = t_size/4) + 
  coord_fixed(c_f) + 
  # geom_segment(data = annot[1, ], aes(x = x1, xend = x2, y = y1, yend = y2)) +
  # geom_text(data = annot[c(2), ], aes(x = x1, y = y1, label = label)) +
  scale_color_discrete(name = "Species\nrecorded", h = c(110, 30), 
                       l = c(10, 65)) +
  scale_fill_manual(name = "Cross\nvalidation\nfold", 
                    values = c("dark grey", "light grey")) +
  xlab("Eastings") + ylab("Northings") + 
  ggtitle("C") + 
  theme_bw() + 
  theme(text = element_text(size = t_size*3), 
        legend.position = "none", 
        # legend.key.width = unit(1.8*t_size, "points"), 
        axis.ticks = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank())
median_obs_map

test_points_map <- ggplot(data = sp_rast_1, cex = 2) + 
  geom_raster(aes(x = x, y = y, fill = fold)) +
  geom_point(data = sp_rast_1[sp_rast_1$en %in% example_folds[[1]]$test_en, ], 
             aes(x = x, y = y, 
                 color = factor(as.character(present), 
                                levels = c("0", "1"), 
                                labels = c("Absent", "Present"))), 
             size = t_size/4) + 
  coord_fixed(c_f) + 
  # geom_segment(data = annot[1, ], aes(x = x1, xend = x2, y = y1, yend = y2)) +
  # geom_text(data = annot[c(2), ], aes(x = x1, y = y1, label = label)) +
  scale_color_discrete(name = "Species\npresent", h = c(110, 30), 
                       l = c(10, 65)) +
  scale_fill_manual(name = "Cross\nvalidation\nfold", 
                    values = c("dark grey", "light grey")) +
  # xlab("Eastings") + ylab("Northings") + 
  ggtitle("D") + 
  theme_bw() + 
  theme(text = element_text(size = t_size*3), 
        legend.position = "none", 
        # legend.key.width = unit(1.8*t_size, "points"), 
        axis.ticks = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank())
test_points_map



#############################
### contour map of AUC by nobs and bias  
# code from Dina email 10 Sep 2019
# large community simulation

f1 <- list(family = "Arial, sans-serif",
           size = t_size*1.5) 
m <- list(l = 90, b = 100, t = 50)
x_axis <- list(title = "Mean number of records per species",
               showticklabels = TRUE, showgrid = F, 
               titlefont = f1, tickfont = f1, exponentformat = "E")
y_axis <- list(title = "Simpson evenness",
               showticklabels = TRUE, showgrid = F,
               titlefont = f1, tickfont = f1, exponentformat = "E")
y_axis_small_community <- list(showticklabels = TRUE, showgrid = F,
                               titlefont = f1, tickfont = f1, 
                               exponentformat = "E")
label_font <- list(size = t_size*1.3)

# select NBDC datasets to plot for comparison
plot_taxa <- c("insect - moth", "insect - butterfly", "bird")
# get values for plotting NBDC dataset points
nbdc <- evenness[evenness$taxon %in% plot_taxa, ] 

# make a matrix to use in plotly - rows are y axis
auc_mat <- spread(preds_auc_gam[preds_auc_gam$method == "glm", 
                            colnames(preds_auc_gam) %in% c("n_obs_per_sp", 
                                                       "pred_auc", 
                                                       "simpson_evenness")], 
                  key = "simpson_evenness", value = "pred_auc")
nobs_labs <- auc_mat$n_obs_per_sp
simps <- round(as.numeric(as.character(colnames(auc_mat)[2:ncol(auc_mat)])),
               digits = 4) # get simpson values to use on axis label

# get matrix of values tested in sim
simulated_values <- expand.grid(unique(
  ev_df_cv$n_obs_per_sp[ev_df_cv$Method == "glm_poly_wg_block_cv"]), simps)
colnames(simulated_values) <- c("n_obs_per_sp", "simpson_evenness")


# taxa labels to go up from dot in large community plot
large_com_taxa_up <- c("insect - moth")
# large_com_taxa_down_r <- c()
# separate into taxa for which the arrow goes up or down
large_com_taxa_up <- nbdc[nbdc$taxon %in% large_com_taxa_up, ]
# large_com_taxa_down_r <- nbdc[nbdc$taxon %in% large_com_taxa_down_r, ]

# shorten taxon group names
large_com_taxa_up$taxon <- gsub(" \\(.*\\)", "", large_com_taxa_up$taxon)
large_com_taxa_up$taxon <- gsub("insect - ", "", large_com_taxa_up$taxon)
# large_com_taxa_down_r$taxon <- gsub("insect - moth", "", large_com_taxa_down_r$taxon)

contour_large_community <- plot_ly(width = 700, height = 800) %>% # width = 1500, height = 1500
  add_trace(x = simulated_values$n_obs_per_sp,
            y = simulated_values$simpson_evenness, type = "scatter",
            mode = "markers", showlegend = F, color = I("black"), 
            symbol = I("circle-open"), marker = list(size = t_size*0.5)) %>%
  add_trace(x = nbdc$n_obs_per_sp[nbdc$n_obs_per_sp <= 200],
            y = nbdc$simpson_evenness[nbdc$n_obs_per_sp <= 200], 
            type = "scatter", mode = "markers", 
            showlegend = F, color = I("black"), 
            symbol = I("circle"), marker = list(size = t_size*0.5)) %>%
  add_annotations(x = large_com_taxa_up$n_obs_per_sp, 
                  y = large_com_taxa_up$simpson_evenness,
                  text = large_com_taxa_up$taxon, 
                  xref = "x", yref = "y", showarrow = TRUE, arrowhead = 1,
                  arrowsize = .5, xanchor = "right", ax = -30, ay = -15, 
                  font = label_font) %>%
  # add_annotations(x = large_com_taxa_down_r$n_obs_per_sp,
  #                 y = large_com_taxa_down_r$simpson_evenness,
  #                 text = large_com_taxa_down_r$taxon,
  #                 xref = "x", yref = "y", showarrow = TRUE, arrowhead = 1,
  #                 arrowsize = .5, ax = -5, ay = 40,
  #                 font = label_font) %>%
  add_trace(y = simps, z = matrix(auc_mat[, 2:ncol(auc_mat)]),
            type = "contour",
            contours = list(start = 0.52, end = 0.74, showlabels = TRUE, 
                            labelfont = list(size = t_size)), 
            autocolorscale = FALSE, 
            colorscale = list(c(0, 'rgb(245,245,245)'), 
                              c(1, 'rgb(0,0,200')), #c(0, 'rgb(245,245,245)'), c(1, 'rgb(30,30,30)')
            reversescale = TRUE, 
            showlegend = FALSE, showscale = FALSE) %>%
  # layout(title = "(a)", titlefont = label_font) %>% # Large Community Simulation\nGLM performance (modelled by GAM)
  layout(title = "A",  titlefont = list(size = t_size*1.8),
         xaxis = x_axis, yaxis = y_axis,
         showlegend = FALSE, margin = m)
contour_large_community


### contour map of AUC by nobs and bias - small community simulation
# make a matrix to use in plotly - rows are y axis
auc_mat_od <- spread(
  preds_auc_gam_od[preds_auc_gam_od$method == "glm", 
                   colnames(preds_auc_gam_od) %in% c("n_obs_per_sp", 
                                                     "pred_auc", 
                                                     "simpson_evenness")], 
  key = "simpson_evenness", value = "pred_auc")
nobs_labs_od <- auc_mat_od$n_obs_per_sp
simps_od <- round(as.numeric(as.character(colnames(auc_mat_od)[2:ncol(auc_mat_od)])),
                  digits = 4) # get simpson values to use on axis label

# get matrix of values tested in sim
simulated_values_od <- expand.grid(unique(
  od_ev_df_cv$n_obs_per_sp[od_ev_df_cv$Method == "glm_poly_wg_block_cv"]), 
  simps_od)
colnames(simulated_values_od) <- c("n_obs_per_sp", "simpson_evenness")



# labels that go up from dot for small community plot:
text_taxa_up <- c("insect - moth") # "insect - hymenopteran","insect - stonefly (Plecoptera)" "insect - mayfly (Ephemeroptera)", "flowering plant", 
# labels that go down from dot: 
text_taxa_down <- c("insect - butterfly", "bird") # "millipede", "insect - dragonfly (Odonata)"

# separate into taxa for which the arrow goes up or down
nbdc_text_up <- nbdc[nbdc$taxon %in% text_taxa_up, ]
nbdc_text_down <- nbdc[nbdc$taxon %in% text_taxa_down, ]

# shorten taxon group names
nbdc_text_up$taxon <- gsub(" \\(.*\\)", "", nbdc_text_up$taxon)
nbdc_text_up$taxon <- gsub("insect - ", "", nbdc_text_up$taxon)
nbdc_text_down$taxon <- gsub(" \\(.*\\)", "", nbdc_text_down$taxon)
nbdc_text_down$taxon <- gsub("insect - ", "", nbdc_text_down$taxon)

contour_small_community <- plot_ly(width = 800, height = 800) %>% # width = 1500, height = 1500
  add_trace(x = simulated_values_od$n_obs_per_sp,
            y = simulated_values_od$simpson_evenness, type = "scatter",
            mode = "markers", showlegend = F, color = I("black"), 
            symbol = I("circle-open"), marker = list(size = t_size*0.5)) %>%
  add_trace(x = nbdc$n_obs_per_sp,
            y = nbdc$simpson_evenness, type = "scatter",
            mode = "markers", showlegend = F, color = I("black"), 
            symbol = I("circle"), marker = list(size = t_size*0.5)) %>%
  # add_text(x = nbdc_text$n_obs_per_sp, y = nbdc_text$simpson_evenness,
  #          text = nbdc_text$taxon, textposition = "top right") %>%
  add_annotations(x = nbdc_text_down$n_obs_per_sp, 
                  y = nbdc_text_down$simpson_evenness,
                  text = nbdc_text_down$taxon, 
                  xref = "x", yref = "y", showarrow = TRUE, arrowhead = 1, 
                  arrowsize = .5, xanchor = "right", ax = -20, ay = -20, 
                  font = label_font) %>%
  add_annotations(x = nbdc_text_up$n_obs_per_sp, 
                  y = nbdc_text_up$simpson_evenness,
                  text = nbdc_text_up$taxon, 
                  xref = "x", yref = "y", showarrow = TRUE, arrowhead = 1,
                  arrowsize = .5, xanchor = "left", ax = 15, ay = -18, 
                  font = label_font) %>%
  add_trace(y = simps_od, z = matrix(auc_mat_od[, 2:ncol(auc_mat_od)]),
            type = "contour",
            contours = list(start = 0.52, end = 0.74, showlabels = TRUE, 
                            labelfont = list(size = t_size)), 
            autocolorscale = FALSE, 
            colorscale = list(c(0, 'rgb(245,245,245)'), c(1, 'rgb(0,0,200)')),#c(1, 'rgb(30,30,30)')
            reversescale = TRUE, 
            showlegend = TRUE, 
            colorbar = list(title = "AUC", tickfont = label_font, 
                            titlefont = label_font)) %>%
  layout(title = "B", titlefont = list(size = t_size*1.8), 
         xaxis = x_axis, yaxis = y_axis_small_community,
         showlegend = FALSE, margin = m)
contour_small_community


### Tables --------------------------------------------------------------------
## Table of predictor variables
pred_table <- data.frame(
  variable = c("annual minimum temperature (deg C)", 
               "annual maximum temperature (deg C)",
               "annual precipitation (mm)", 
               "average daily sea level atmospheric pressure (hecto Pascals)", 
               "agricultural areas", 
               "artificial surfaces", 
               "forest and semi-natural areas", 
               "water bodies", 
               "wetlands",  
               "elevation"), 
  description = c("2% quantile of annual temperatures in each grid cell averaged over the years 1995-2016", 
                  "98% quantile of annual temperatures in each grid cell averaged over the years 1995-2016",
                  "Average total annual precipitation in each grid cell over the years 1995-2016 (excluding 2010-2012)", 
                  "Average daily sea level atmospheric pressure over the years 1995-2016",
                  "Proportion of each grid cell classified as agricultural areas", 
                  "Proportion of each grid cell classified as artificial surfaces", 
                  "Proportion of each grid cell classified as forest and semi-natural areas", 
                  "Proportion of each grid cell classified as water bodies",
                  "Proportion of each grid cell classified as wetlands", 
                  "Average elevation in each grid cell"), 
  data_source = c("E-OBS", 
                  "E-OBS", 
                  "E-OBS", 
                  "E-OBS", 
                  "CORINE Land Cover Database", 
                  "CORINE Land Cover Database", 
                  "CORINE Land Cover Database", 
                  "CORINE Land Cover Database", 
                  "CORINE Land Cover Database", 
                  "ETOPO1"), 
  Moran_I = NA)

pred_table$Moran_I[grepl(".*minimum.temp.*", 
                         pred_table$variable)] <- round(Moran(
                           mask(pred_rast_brick$minimum.temperature, ir_TM75)), 
                           digits = 2)
pred_table$Moran_I[grepl(".*maximum.temp.*", 
                         pred_table$variable)] <- round(Moran(
                           mask(pred_rast_brick$maximum.temperature, ir_TM75)), 
                           digits = 2)
pred_table$Moran_I[grepl(".*annual.prec.*", 
                         pred_table$variable)] <- round(Moran(
                           mask(pred_rast_brick$annual.precipitation, ir_TM75)), 
                           digits = 2)
pred_table$Moran_I[grepl(".*atmosp.*", 
                         pred_table$variable)] <- round(Moran(
                           mask(pred_rast_brick$atmospheric.pressure, ir_TM75)), 
                           digits = 2)
pred_table$Moran_I[grepl(".*agricul.*", 
                         pred_table$variable)] <- round(Moran(
                           mask(pred_rast_brick$agricultural.areas, ir_TM75)), 
                           digits = 2)
pred_table$Moran_I[grepl(".*artific.*", 
                         pred_table$variable)] <- round(Moran(
                           mask(pred_rast_brick$artificial.surfaces, ir_TM75)), 
                           digits = 2)
pred_table$Moran_I[grepl(".*forest.*", 
                         pred_table$variable)] <- round(Moran(
                           mask(pred_rast_brick$forest.semi.natural, ir_TM75)), 
                           digits = 2)
pred_table$Moran_I[grepl(".*wetlands.*", 
                         pred_table$variable)] <- round(Moran(
                           mask(pred_rast_brick$wetlands, ir_TM75)), 
                           digits = 2)
pred_table$Moran_I[grepl(".*water.*", 
                         pred_table$variable)] <- round(Moran(
                           mask(pred_rast_brick$water.bodies, ir_TM75)), 
                           digits = 2)
pred_table$Moran_I[grepl(".*elevation.*", 
                         pred_table$variable)] <- round(Moran(
                           mask(pred_rast_brick$elevation, ir_TM75)), 
                           digits = 2)
pred_table


# # find range of spatial autocorrelation
# gs_minTemp <- gstat(id = "min_temp", 
#                  formula = minimum.temperature ~ 1,
#                  data = data.frame(as(mask(pred_rast_brick, ir_TM75), 
#                                       "SpatialGridDataFrame")), 
#                  locations = ~ s1 + s2)
# plot(variogram(gs_minTemp))
# 
# gs_elev <- gstat(id = "elev", 
#                  formula = elevation ~ 1,
#                  data = data.frame(as(mask(pred_rast_brick, ir_TM75), 
#                                       "SpatialGridDataFrame")), 
#                  locations = ~ s1 + s2)
# var_elev <- variogram(gs_elev)
# f_vg_elev <- fit.variogram(var_elev, 
#                            vgm(c("Nug", "Exp", "Sph", "Gau", "Exc", "Mat", "Ste")))
# plot(variogramLine(f_vg_elev, max(var_elev$dist)), type = 'l')
# points(var_elev[, 2:3], pch = 20, col = 'red')


## Table of variable importance
auc_summary

## Simpson's evenness
nbdc # this shows evenness values for NBDC template datasets
# this shows evenness for the "no bias" template
simpson_even(bias_df$layer[bias_df$bias == "A" & !is.na(bias_df$layer)])

#### print pdf -------------------------------------------------------------
contour_large_community # save this manually or use orca()
contour_small_community # save this manually or use orca()

ggsave("Fig1.svg", sampling_bias_maps, width = 25, height = 25, units = "cm", 
       device = "svg")
ggsave("Fig3.svg", multiplot(p_dist_1, p_dist_3, p_dist_2, p_dist_4, 
                             cols = 2), 
       width = 25, height = 25, units = "cm", 
       device = "svg")

ggsave("Fig4.svg", n_rmse_point, width = 25, height = 25, units = "cm", 
       device = "svg")
ggsave("Fig5.svg", auc_line, width = 25, height = 25, units = "cm", 
       device = "svg")
ggsave("Fig6.svg", auc_boxplot, width = 25, height = 25, units = "cm", 
       device = "svg")
ggsave("Fig7.svg", rmse_cv_boxplot, width = 25, height = 25, units = "cm", 
       device = "svg")
ggsave("example_sp.svg", example_sp, width = 25, height = 25, units = "cm", 
       device = "svg")
ggsave("noBias_obs_map.svg", noBias_obs_map, width = 25, height = 25, units = "cm", 
       device = "svg")
ggsave("median_obs_map.svg", median_obs_map, width = 25, height = 25, units = "cm", 
       device = "svg")
ggsave("test_points_map.svg", test_points_map, width = 25, height = 25, units = "cm", 
       device = "svg")


ggsave("Fig1.jpg", sampling_bias_maps, width = 25, height = 25, units = "cm", 
       device = "jpg")
ggsave("Fig3.jpg", multiplot(p_dist_1, p_dist_3, p_dist_2, p_dist_4, 
                             cols = 2), 
       width = 25, height = 25, units = "cm", 
       device = "jpg")
ggsave("Fig4.jpg", n_rmse_point, width = 25, height = 25, units = "cm", 
       device = "jpg")
ggsave("Fig5.jpg", auc_line, width = 25, height = 25, units = "cm", 
       device = "jpg")
ggsave("Fig5_slideshow.jpg", auc_line_presentation, width = 25, height = 25, 
       units = "cm", device = "jpg")
ggsave("Fig6.jpg", auc_boxplot, width = 25, height = 25, units = "cm", 
       device = "jpg")
ggsave("Fig7.jpg", rmse_cv_boxplot, width = 25, height = 25, units = "cm", 
       device = "jpg")
ggsave("example_sp.jpg", example_sp, width = 25, height = 25, units = "cm", 
       device = "jpg")
ggsave("noBias_obs_map.jpg", noBias_obs_map, width = 25, height = 25, units = "cm", 
       device = "jpg")
ggsave("median_obs_map.jpg", median_obs_map, width = 25, height = 25, units = "cm", 
       device = "jpg")
ggsave("test_points_map.jpg", test_points_map, width = 25, height = 25, units = "cm", 
       device = "jpg")

write_csv(pred_table, "predictor_variable_table.csv")


### save plots as pdf
pdf("Fig1.pdf")
sampling_bias_maps
dev.off()
pdf("Fig3.pdf")
multiplot(p_dist_1, p_dist_3, p_dist_2, p_dist_4, 
          cols = 2)
dev.off()
pdf("Fig4.pdf")
n_rmse_point
dev.off()
pdf("Fig5.pdf")
auc_line
dev.off()
pdf("Fig6.pdf")
auc_boxplot
dev.off()
pdf("Fig7.pdf")
rmse_cv_boxplot
dev.off()
