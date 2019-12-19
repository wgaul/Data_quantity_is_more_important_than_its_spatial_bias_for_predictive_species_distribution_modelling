############################################
## This generates plots for the main text of the simulation paper.
## This is for simulation 10.10 using the bryophyte template (large community)
## 
## This script is meant to be run after sourcing analyze_performance.R
##
## inputs:  * bryophyte_template_evals_9Sep2019.csv 
##          * bryophyte_template_butterfly_bias_evals_9Sep2019.csv
##
## 
## author: Willson Gaul wgaul@hotmail.com
## created: 10 Sep 2019 based on code from plots.R
## last modified: 22 Nov 2019
############################################

library(wgutil)
#library(skimr)
library(Hmisc)
library(rgdal)
library(raster)
library(parallel)
library(pROC)
# library(akima)
library(tidyverse)
library(plotly)

setwd("~/Documents/Data_Analysis/UCD/simulation/sims_10.10_bryophyte/")

sim_name <- "Bryophyte Simulation"

bias_names <- c("no_bias", "extreme_mothNBDC", "moderate_bryBBS",
                "butterflyNBDC", "least_birdNBDC")

## optionally, read in a .csv with auc results from simulation
## Large community simulation results
ev_df <- read_csv("./bryophyte_template_evals_13Sep2019.csv")
df_butt <- read_csv("../sims_10.10_bryophyte_butterflyBias/bryophyte_template_butterfly_bias_evals_12Sep2019.csv")
ev_df <- bind_rows(ev_df, df_butt)

# subset to only block CV results
ev_df_cv <- ev_df[which(ev_df$Method %in% c("glm_poly_wg_block_cv", 
                                            "random_forest_wg_block_cv", 
                                            "brt_wg_block_cv", 
                                            "idw_interp_wg_block_cv")), ]

## small community simulation results
od_ev_df <- read_csv("../sims_10.10_odonata/odonata_template_evals.csv")
od_butBias_df <- read_csv("../sims_10.10_odonata_butterflyBias/odonata_template_butterfly_bias_evals_25Aug2019.csv")
od_ev_df <- bind_rows(od_ev_df, od_butBias_df)

# remove methods that did not fully run (see simulation_tracking_status.ods)
od_ev_df <- od_ev_df[od_ev_df$Method %in% c("glm_poly", "glm_poly_wg_block_cv", 
                                            "idw_interp", 
                                            "idw_interp_wg_block_cv"), ]
# get only CV results
od_ev_df_cv <- od_ev_df[od_ev_df$Method %in% c("glm_poly_wg_block_cv",
                                               "idw_interp_wg_block_cv"), ]

t_size <- 19
line_size <- 1.5
c_f = 1 # coord fixed for plots

### plot sampling biases ------------------------------------------------------
# load sample bias rasters
load("~/Documents/Data_Analysis/UCD/sampling_distribution/saved_objects/spat_nrec_hec.RData")
bias_rasters$butterfly <- as.data.frame(spat_nrec_hec$`insect - butterfly`$rast, 
                                        xy = TRUE)
bias_df <- bias_rasters
# set value for no bias map to 0.05 in each cell.  This is arbitrary but needed 
# to turn cells purple so it is clear there is sampling happening and it's even.
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
                                                "moderate_bryBBS", 
                                                "butterfly", 
                                                "extreme_mothNBDC"), 
                       labels = c("(a)", "(b)", "(c)", "(d)", "(e)"))

annot <- data.frame(x1 = c(265000, 310000, 60000, 60000), 
                    x2 = c(365000, 310000, 60000, 60000), 
                    y1 = c(60000, 40000, 400000, 380000), 
                    y2 = c(60000, 40000, 455000, 380000),
                    label = c(NA, "100 km", NA, "N"), 
                    bias = "(e)")

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
  # annotate("segment", x = 265000, xend = 365000, y = 60000, yend = 60000) + 
  # annotate("text", x = 310000, y = 40000, label = "100 km") + 
  theme_bw() + 
  theme(text = element_text(size = t_size), 
        axis.title = element_blank(),
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        strip.text = element_text(hjust = -0.01), 
        legend.position = c(0.84, 0.27))
sampling_bias_maps
### end plot sampling biases --------------------------------------------------


### plot number of models fit --------------------------------------------------
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
                                        "moderate_bryBBS", 
                                        "butterflyNBDC", 
                                        "extreme_mothNBDC"), 
                             labels = c("\nno bias\n", 
                                        "\nlow\n",
                                        "\nmoderate\n", 
                                        "\nmedian\n",  
                                        "\nsevere\n")),
                           color = factor(
                             bias, 
                             levels = c("no_bias", 
                                        "least_birdNBDC", 
                                        "moderate_bryBBS", 
                                        "butterflyNBDC", 
                                        "extreme_mothNBDC"), 
                             labels = c("\nno bias\n", 
                                        "\nlow\n",
                                        "\nmoderate\n", 
                                        "\nmedian\n",  
                                        "\nsevere\n")))) + 
  geom_point() + 
  geom_line(size = line_size) + 
  facet_wrap(~ factor(method, 
                      levels = c("glm_poly_wg_block_cv", 
                                 "random_forest_wg_block_cv", 
                                 "brt_wg_block_cv", 
                                 "idw_interp_wg_block_cv"), 
                      labels = c("(a)", "(b)", "(c)", "(d)"))) + 
  # ggtitle("Number of species for which models were sucessfully fitted") + 
  ylab("Number of species fit") + 
  xlab("Average number of records\nper species in dataset") + 
  scale_shape_discrete(name = "Spatial Bias", solid = F) +  
  scale_colour_viridis_d(option = "magma", name = "Spatial Bias", 
                         begin = 0, end = 0.65) + 
  theme_bw() + 
  theme(text = element_text(size = t_size*1.2), 
        legend.key.width = unit(2*t_size, "points"))
n_rmse_point
### end plot number of models fit ----------------------------------------------


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
rug_df$method <- gsub("random", "rf", rug_df$method)

auc_line <- ggplot(data = preds_auc_brt, 
                     aes(x = n_obs_per_sp, y = pred_auc, 
                         group = factor(bias.name,
                                        levels = c("no bias", "bird", "bry", 
                                                   "butterfly", "moth"),
                                        labels = c(
                                          "\nnone\n", "\nlow\n", 
                                          "\nmoderate\n",
                                          "\nmedian\n", 
                                          "\nsevere\n"),
                                        ordered = TRUE))) + 
  geom_line(aes(color = factor(bias.name,
                               levels = c("no bias", "bird", "bry", 
                                          "butterfly", "moth"),
                               labels = c(
                                 "\nnone\n", "\nlow\n", 
                                 "\nmoderate\n",
                                 "\nmedian\n", 
                                 "\nsevere\n"),
                               ordered = TRUE) #, 
                # linetype = factor(bias.name,
                #                   levels = c("no bias", "bird", "bry", 
                #                              "butterfly", "moth"),
                #                   labels = c(
                #                     "\nnone\n", "\nlow\n(Birds)", 
                #                     "\nmoderate\n(Bryophytes)\n",
                #                     "\nmedian\n(Butterflies)\n", 
                #                     "\nsevere\n(Moths)\n"),
                #                   ordered = TRUE)
                ), 
            size = line_size) + 
  # ggtitle("Large Community Simulation") + #sim_name
  geom_rug(data = rug_df, aes(x = n_obs_per_sp), show.legend = F, sides = "b") +
  facet_wrap( ~ factor(method, 
                       levels = c("glm", "rf", "brt", "idwi"), 
                       labels = c("(a)", "(b)", 
                                  "(c)", 
                                  "(d)"))) + 
  xlab("Mean number of records \nper species in data set") + 
  ylab("AUC\n(block cross-validated)") +
  # scale_linetype_discrete(name = "Spatial\nSampling\nBias") + 
  scale_color_viridis_d(name = "Spatial\nSampling\nBias", option = "magma", 
                        begin = 0, end = 0.65) +
  theme_bw() + 
  theme(text = element_text(size = t_size*1.2), 
        legend.key.width = unit(1.8*t_size, "points"))
auc_line

auc_line_presentation <- ggplot(data = preds_auc_brt, 
                                aes(x = n_obs_per_sp, y = pred_auc, 
                                    group = factor(bias.name,
                                                   levels = c("no bias", "bird", 
                                                              "bry", 
                                                              "butterfly", 
                                                              "moth"),
                                                   labels = c(
                                                     "\nnone\n", "\nlow\n", 
                                                     "\nmoderate\n",
                                                     "\nmedian\n", 
                                                     "\nsevere\n"),
                                                   ordered = TRUE))) + 
  geom_line(aes(color = factor(bias.name,
                               levels = c("no bias", "bird", "bry", 
                                          "butterfly", "moth"),
                               labels = c(
                                 "\nnone\n", "\nlow\n", 
                                 "\nmoderate\n",
                                 "\nmedian\n", 
                                 "\nsevere\n"),
                               ordered = TRUE) #, 
                # linetype = factor(bias.name,
                #                   levels = c("no bias", "bird", "bry", 
                #                              "butterfly", "moth"),
                #                   labels = c(
                #                     "\nnone\n", "\nlow\n(Birds)", 
                #                     "\nmoderate\n(Bryophytes)\n",
                #                     "\nmedian\n(Butterflies)\n", 
                #                     "\nsevere\n(Moths)\n"),
                #                   ordered = TRUE)
                ), 
  size = 1.2*line_size) + 
  # ggtitle("Large Community Simulation") + #sim_name
  geom_rug(data = rug_df, aes(x = n_obs_per_sp), show.legend = F, sides = "b") +
  facet_wrap( ~ factor(method, 
                       levels = c("glm", "rf", "brt", "idwi"), 
                       labels = c("GLM", "Random\nforest", 
                                  "Boosted\nregression\ntree", 
                                  "Inverse\ndistance-\nweighted\ninterpolation"))) + 
  xlab("Mean number of\nrecords per species") + 
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
                                                "moderate_bryBBS", 
                                                "butterflyNBDC", 
                                                "extreme_mothNBDC"),
                                     labels = c("none",
                                                "low",
                                                "moderate", 
                                                "median", 
                                                "severe"),
                                     ordered = TRUE), 
                          y = auc, 
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
                      labels = c("(a)", "(b)", "(c)","(d)"))) + 
  # ggtitle("Large Community Simulation") + #sim_name
  xlab("Spatial Sampling Bias") + 
  ylab("AUC\n(block cross-validated)") + 
  scale_fill_viridis_d(name = "Mean\nnumber of\nrecords\nper species", 
                       option = "viridis", begin = 0, end = 0.85) + 
  theme_bw() + 
  theme(text = element_text(size = t_size*1.2), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
auc_boxplot
## end AUC boxplots -----------------------------------------------------------

## plot methods by bias - AUC --------------------------------------------------
# use only block CV results
methods_bias_line <- ggplot(data = preds_auc_brt, aes(
  x = as.numeric(as.character(n_obs_per_sp)), 
  y = pred_auc, 
  group = factor(method, 
                 levels = c("glm", "rf", "brt", "idwi"), 
                 labels = c("GLM", "Random\nForest", 
                            "Boosted\nRegression\nTree", 
                            "Inverse distance-weighted\ninterpolation")))) +
  geom_line(aes(
    linetype = factor(method, 
                      levels = c("glm", "rf", "brt", "idwi"), 
                      labels = c("GLM", "\nRandom\nForest", 
                                 "\nBoosted\nRegression\nTree", 
                                 "\nInverse\ndistance-weighted\ninterpolation")), 
    color = factor(method, levels = c("glm", "rf", "brt", "idwi"), 
                   labels = c("GLM", "\nRandom\nForest", 
                              "\nBoosted\nRegression\nTree", 
                              "\nInverse\ndistance-weighted\ninterpolation"))), 
    size = line_size) +
  geom_rug(data = rug_df, aes(x = n_obs_per_sp), show.legend = F, sides = "b") +
  facet_wrap(~ factor(bias.name,
                      levels = c("no bias", "bird", "bry", 
                                 "butterfly", "moth"),
                      labels = c("no bias",
                                 "low bias\n(Birds)",
                                 "moderate bias\n(Bryophytes)", 
                                 "median bias\n(Butterflies)", 
                                 "severe bias\n(Moths)"),
                      ordered = TRUE)) + 
  # ggtitle("Large Community Simulation") + #sim_name
  ylab("AUC\n(block cross-validated)") + 
  xlab("Number of observations\nper species") + 
  scale_linetype_discrete(name = "SDM\nModelling\nMethod") +
  scale_color_viridis_d(name = "SDM\nModelling\nMethod", option = "C") + 
  theme_bw() + 
  theme(text = element_text(size = t_size), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        legend.key.width = unit(2*t_size, "points"))
methods_bias_line
## end plot methods by bias - AUC ----------------------------------------------

## end AUC plots
##############################################################################



###################################
## RMSE plots 
## GAM smooth of RMSE by bias
# use predictions from the four fitted GAM models (one model for each SDM method)
rmse_line <- ggplot(data = preds_rmse_evaluation_brt, 
                     aes(x = n_obs_per_sp, y = pred_rmse, 
                         group = factor(bias.name,
                                        levels = c("no bias", "bird", "bry", 
                                                   "butterfly", "moth"),
                                        labels = c(
                                          "\nnone\n", "\nlow\n", 
                                          "\nmoderate\n", "\nmedian\n", 
                                          "\nsevere\n"),
                                        ordered = TRUE))) + 
  geom_line(aes(color = factor(bias.name,
                               levels = c("no bias", "bird", "bry", 
                                          "butterfly", "moth"),
                               labels = c(
                                 "\nnone\n", "\nlow\n", "\nmoderate\n", 
                                 "\nmedian\n", "\nsevere\n"),
                               ordered = TRUE)), size = line_size) + 
  geom_rug(data = rug_df, aes(x = n_obs_per_sp), show.legend = F, sides = "b") +
  # ggtitle("Large Community Simulation") + #sim_name
  facet_wrap( ~ factor(method, 
                       levels = c("glm", "rf", "brt", "idwi"), 
                       labels = c("(a)", "(b)", "(c)", "(d)"))) + 
  xlab("Mean number of records \nper species in data set") + 
  ylab("RMSE\n(block cross-validated)") +
  # scale_linetype_discrete(name = "Spatial\nSampling\nBias") + 
  scale_color_viridis_d(name = "Spatial\nSampling\nBias", option = "magma", 
                        begin = 0, end = 0.65) +
  theme_bw() + 
  theme(text = element_text(size = t_size*1.2), 
        legend.key.width = unit(1.8*t_size, "points"))
rmse_line


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
plot_taxa <- c("insect - moth", "insect - butterfly", 
               "insect - dragonfly (Odonata)", "bry_BBS", "bird")
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
large_com_taxa_up <- c("bry_BBS")
# large_com_taxa_down_r <- c()
# separate into taxa for which the arrow goes up or down
large_com_taxa_up <- nbdc[nbdc$taxon %in% large_com_taxa_up, ]
# large_com_taxa_down_r <- nbdc[nbdc$taxon %in% large_com_taxa_down_r, ]

# shorten taxon group names
large_com_taxa_up$taxon[large_com_taxa_up$taxon == "bry_BBS"] <- "bryophyte"
large_com_taxa_up$taxon <- gsub(" \\(.*\\)", "", large_com_taxa_up$taxon)
large_com_taxa_up$taxon <- gsub("insect - ", "", large_com_taxa_up$taxon)
# large_com_taxa_down_r$taxon <- gsub("insect - moth", "", large_com_taxa_down_r$taxon)

contour_large_community <- plot_ly(width = 700, height = 800) %>% # width = 1500, height = 1500
  add_trace(x = simulated_values$n_obs_per_sp,
            y = simulated_values$simpson_evenness, type = "scatter",
            mode = "markers", showlegend = F, color = I("black"), 
            symbol = I("circle-open"), marker = list(size = t_size*0.5)) %>%
  add_trace(x = nbdc$n_obs_per_sp[nbdc$n_obs_per_sp <= 160],
            y = nbdc$simpson_evenness[nbdc$n_obs_per_sp <= 160], 
            type = "scatter", mode = "markers", 
            showlegend = F, color = I("black"), 
            symbol = I("circle"), marker = list(size = t_size*0.5)) %>%
  add_annotations(x = large_com_taxa_up$n_obs_per_sp, 
                  y = large_com_taxa_up$simpson_evenness,
                  text = large_com_taxa_up$taxon, 
                  xref = "x", yref = "y", showarrow = TRUE, arrowhead = 1,
                  arrowsize = .5, xanchor = "right", ax = -4, ay = -40, 
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
  layout(title = "(a)",  titlefont = list(size = t_size*1.8),
         xaxis = x_axis, yaxis = y_axis,
         showlegend = FALSE, margin = m)
contour_large_community


### contour map of AUC by nobs and bias - small community simulation
# make a matrix to use in plotly - rows are y axis
auc_mat_od <- spread(preds_auc_gam_od[preds_auc_gam_od$method == "glm", 
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
text_taxa_up <- c("insect - moth",  
               "bry_BBS") # "insect - hymenopteran","insect - stonefly (Plecoptera)" "insect - mayfly (Ephemeroptera)", "flowering plant", 
# labels that go down from dot: 
text_taxa_down <- c("insect - butterfly", "bird") # "millipede", "insect - dragonfly (Odonata)"

# separate into taxa for which the arrow goes up or down
nbdc_text_up <- nbdc[nbdc$taxon %in% text_taxa_up, ]
nbdc_text_down <- nbdc[nbdc$taxon %in% text_taxa_down, ]

# shorten taxon group names
nbdc_text_up$taxon[nbdc_text_up$taxon == "bry_BBS"] <- "bryophyte"
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
                  arrowsize = .5, xanchor = "left", ax = 20, ay = 20, 
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
  layout(title = "(b)", titlefont = list(size = t_size*1.8), 
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



#### print pdf -------------------------------------------------------------
pdf("plots_main_text.pdf")
print(sampling_bias_maps)
print(n_rmse_point)
print(auc_line)
print(auc_boxplot)
print(rmse_line)
dev.off()
contour_large_community # save this manually or use orca()
contour_small_community # save this manually or use orca()

ggsave("Fig1.svg", sampling_bias_maps, width = 25, height = 25, units = "cm", 
       device = "svg")
ggsave("Fig2.svg", n_rmse_point, width = 25, height = 25, units = "cm", 
       device = "svg")
ggsave("Fig3.svg", auc_line, width = 25, height = 25, units = "cm", 
       device = "svg")
ggsave("Fig4.svg", auc_boxplot, width = 25, height = 25, units = "cm", 
       device = "svg")
ggsave("Fig5.svg", rmse_line, width = 25, height = 25, units = "cm", 
       device = "svg")


ggsave("Fig1.jpg", sampling_bias_maps, width = 25, height = 25, units = "cm", 
       device = "jpg")
ggsave("Fig2.jpg", n_rmse_point, width = 25, height = 25, units = "cm", 
       device = "jpg")
ggsave("Fig3.jpg", auc_line, width = 25, height = 25, units = "cm", 
       device = "jpg")
ggsave("Fig3_slideshow.jpg", auc_line_presentation, width = 25, height = 25, 
       units = "cm", device = "jpg")
ggsave("Fig4.jpg", auc_boxplot, width = 25, height = 25, units = "cm", 
       device = "jpg")
ggsave("Fig5.jpg", rmse_line, width = 25, height = 25, units = "cm", 
       device = "jpg")

write_csv(pred_table, "predictor_variable_table.csv")
