############################################
## This is the main simulator file for the species simulation
## This is for simulation 10.10 
## Using bryophytes as template
##
## inputs:  * ireland_coastline shapefile
##          * elevation_hec_ETOPO1.RData - elevation predictor variable 
##              interpolated to hectad scale
##          * annual_precip_hectad.RData - predictor variable interpolated from 
##              EOBS
##          * winter_tn_hectad.RData - minimum winter temperature (2nd quantile)
##              interpolated from EOBS
##          * summer_tx_hectad.RData - maximum summer temperature (98th 
##              quantile) interpolated from EOBS
##          * mean_pp_hectad.RData - mean daily average sea level pressure 
##              (units = hecto Pascals)
##          * corine_label_1_hectad.RData - rasters with proportions of CORINE
##              LABEL 1 land cover classes in hectads
##          * spat_nrec_hec.RData - rasters giving sampling bias weights based
##              on number of records per hectad in NBDC datasets
##
## TODO: figure out where to set seed
## 
## author: Willson Gaul wgaul@hotmail.com
## created: 24 Sep 2018
## last modified: 19 Dec 2019
############################################

rm(list = ls())
library(wgutil)
library(Hmisc)
library(rgdal)
library(raster)
library(parallel)
library(randomForest) 
library(dismo)
library(pROC)
library(gstat)
library(unmarked)
library(tidyverse)
library(virtualspecies)
library(simulator) # this file was created under simulator version 0.2.0

new_draws <- T  # run the simulation again because parameters have been changed
new_methods <- T
new_evals <- T
add_to_sim <- F # add draws or models to an existing simulation
analyze_performance <- T
print_plots <- T
dbg <- F # debug

n_cores <- 1
seed <- 10061983 + 12202018 + 1300  # wg bday + today's date + current time

options("scipen" = 100, "digits" = 4) # so that file names don't contain scientific notation

setwd("~/Documents/Data_Analysis/UCD/simulation/sims_10.10_bryophyte/")
sim_dest_dir <- "./"

# load data (this will need to be modified by each user)
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
# load sample bias rasters
load("./spat_nrec_hec.RData")
# load data for defining list lengths and community size
load("./community_size_and_list_lengths.RData")


source("./simulation_functions.R")
source("./utility_functions.R")
source("./model_functions.R")
source("./method_functions.R")
source("./eval_functions.R")


## set simulation parameters --------------------------------------------------
sim_name <- "bryophyte_template"
lab <- "Bryophyte template"
n_obs <- list(100, 1000, 2500, 5000, 10000, 25000, 100000, 200000) 
nsim <- 240 # number of species to simulate.  Ideally divisible by n_cores.
# set community size
community_size <- community_list_lengths$bry_BBS$community_size
list_length <- community_list_lengths$bry_BBS$list_lengths_std

## end set simulation params --------------------------------------------------


## Function definitions -----------------------------------------------------
calc_coef_from_vertex <- function(vert_x, vert_y, a) {
  # calculate coefficients based on location of vertex of a parabola
  # ax^2 + bx + c = a(x - H)^2 + K  where H is v_x and K is v_y
  # so, b = -a2H  and  c = a(H^2) + K
  # from https://www.geogebra.org/m/FhP89QWZ#material/RVF4GYcX
  # ARGS: vert_x - x coordinate of vertex
  #       vert_x - y coordinate of vertex
  #       a - coefficient of the squared term
  # RETURN: a vector giving the a, b and c coefficients
  b <- -a*2*vert_x
  c <- a*(vert_x^2) + vert_y
  results <- c(a = a, b = b, c = c)
  results
}
## end function definitions -----------------------------------------------

## specify bias rasters -----------------------------------------------------
## TODO: make additional bias rasters where sampling is correlated with PC1
# make a no-bias raster using one of the existing rasters as a template
no_bias_rast <- spat_nrec_hec$bird$rast
no_bias_rast <- calc(no_bias_rast, fun = function(x) {
  n_cell <- length(x[!is.na(x)])
  x[!is.na(x)] <- 1/n_cell
  x
})

bias_rasters <- list(no_bias = no_bias_rast, 
                     extreme_mothNBDC = spat_nrec_hec$`insect - moth`$rast, 
                     moderate_bryBBS = spat_nrec_hec$bry_BBS$rast, 
                     least_birdNBDC = spat_nrec_hec$bird$rast)
## TODO: 28 Sep.  Moved this next bit into simulation function for now as it 
## needs to be done after I re-create the rasters in the function.  Do this more
## elegantly if I decide to keep it that way.  Alternatively, un-comment this
## to have it done here again.
# bias_rasters <- lapply(bias_rasters, FUN = function(x, ir_TM75) {
#   raster::mask(x, mask = ir_TM75)}, ir_TM75 = ir_TM75)
## end specify bias rasters -------------------------------------------------

## combine all predictors into a single df ------------------------------------
# initiate dataframe
## TODO: masking by ir_TM75 throws out some hectads that do have records in the
## real NBDC data I think (reducing total hectads from around 1003 to 840).  Do
## I want to somehow keep all 1003 grid cells?
raw_pred_df <- full_join(data.frame(as(mask(krg_mean_rr_rast, ir_TM75), 
                                       "SpatialGridDataFrame")), 
                         data.frame(as(mask(krg_mean_tn_rast, ir_TM75), 
                                       "SpatialGridDataFrame")), 
                         by = c("s1", "s2"))
raw_pred_df$s1 <- round(raw_pred_df$s1) # fix number representation issues
raw_pred_df$s2 <- round(raw_pred_df$s2)

# put remaining predictors in a list
pred_l <- c(krg_mean_tx_rast, krg_mean_pp_rast, elev_hec, agricultural_l1_rast, 
            artificial_surfaces_l1_rast, forest_seminatural_l1_rast, 
            water_l1_rast, wetlands_l1_rast)

# add each remaining predictor to df
for (i in 1:length(pred_l)) {
  df <- data.frame(as(mask(pred_l[[i]], ir_TM75), 
                      "SpatialGridDataFrame"))
  df$s1 <- round(df$s1) # fix number representation issues
  df$s2 <- round(df$s2)
  raw_pred_df <- full_join(raw_pred_df, df, 
                           by = c("s1", "s2"))
}
rm(df)

# reorder columns to put "s1" and "s2" (grid coordinate) columns first
new_order <- names(raw_pred_df)
new_order <- new_order[c(2:3, 1, 4:length(new_order))]
raw_pred_df <- raw_pred_df[, new_order]

## make predictor df to use for model fitting
# center and square predictors for using in model fitting
sq_pred_df <- raw_pred_df
sq_pred_df$mean_annual_rr_c <- sq_pred_df$mean_annual_rr - 
  mean(sq_pred_df$mean_annual_rr)

center <- function(x) {x - mean(x)} # function to center variables

sq_pred_df <- sq_pred_df %>%
  mutate(mean_annual_rr_c = center(mean_annual_rr), 
         mean_winter_tn_c = center(mean_winter_tn), 
         mean_tx_c = center(mean_tx), 
         mean_pp_c = center(mean_pp), 
         elevation.pred_c = center(elevation.pred), 
         Agricultural.areas_c = center(Agricultural.areas), 
         Artificial.surfaces_c = center(Artificial.surfaces), 
         Forest.and.semi.natural.areas_c = 
           center(Forest.and.semi.natural.areas), 
         Water.bodies_c = center(Water.bodies), 
         Wetlands_c = center(Wetlands)) %>%
  mutate(mean_annual_rr_c_sq = mean_annual_rr_c^2, 
         mean_winter_tn_c_sq = mean_winter_tn_c^2, 
         mean_tx_c_sq = mean_tx_c^2, 
         mean_pp_c_sq = mean_pp_c^2, 
         elevation.pred_c_sq = elevation.pred_c^2, 
         Agricultural.areas_c_sq = Agricultural.areas_c^2, 
         Artificial.surfaces_c_sq = Artificial.surfaces_c^2, 
         Forest.and.semi.natural.areas_c_sq = 
           Forest.and.semi.natural.areas_c^2, 
         Water.bodies_c_sq = Water.bodies_c^2, 
         Wetlands_c_sq = Wetlands_c^2)


# remove original predictors to clean workspace
rm(year_tot_rr, year_winter_tn, year_summer_tx, elev_hec, krg_mean_rr_predict, 
   krg_mean_rr_rast, krg_mean_tn_predict, krg_mean_tn_rast, krg_mean_tx_predict,
   krg_mean_tx_rast, krig_mean_rr_map, krig_mean_tn_map, krig_mean_tx_map, 
   agricultural_l1_rast, artificial_surfaces_l1_rast, clc_l1_props_hecs, 
   forest_seminatural_l1_rast, water_l1_rast, wetlands_l1_rast, pred_l, 
   new_order, krg_mean_pp_rast, krg_mean_pp, krig_mean_pp_map)
## end combine predictors in to single df -----------------------------------

## calculate metrics for datasets -------------------------------------------
evenness <- mapply(FUN = function(x, nm) {
  data.frame(taxon = as.character(nm), camargo = camargo(x$spat_df$nrec), 
             shannon_evenness = shannon_even(x$spat_df$nrec), 
             simpson_evenness = simpson_even(x$spat_df$nrec))
}, spat_nrec_hec, names(spat_nrec_hec), SIMPLIFY = F)

evenness <- bind_rows(evenness)
evenness <- evenness[order(evenness$camargo), ]

# add original sample size to evenness dataframe
evenness <- evenness[evenness$taxon %in% names(community_list_lengths), ]
evenness$n_obs_per_sp <- NA
for(i in 1:nrow(evenness)) {
  if(evenness$taxon[i] == "bry_BBS") {
    evenness$n_obs_per_sp[i] <- sum(spat_nrec_hec$bry_BBS$spat_df$nrec) / 
      community_list_lengths$bry_BBS$community_size
  } else {
    tname <- evenness$taxon[i]
    try(
      evenness$n_obs_per_sp[i] <- tryCatch({
        sum(spat_nrec_hec[[tname]]$spat_df$nrec) / 
          community_list_lengths[[tname]]$community_size}, 
        error = function(x) NA))
  }
}
## end calculate metrics for datasets ----------------------------------------
rm(spat_nrec_hec)
rm(no_bias_rast)

## @knitr init

## @knitr main ----------------------------------------------------------------
## TODO: 28 Sep. Converting bias rasters to data frames so sha1() can deal with 
## them. This is just a quick fix to let the simulator process run.  If I end 
## up keeping this change, will probably need to pass all the information for
## converting it back to a raster as arguments so the dfs can be turned back in
## to rasters inside the simulation functions. 
bias_rasters <- lapply(bias_rasters, as.data.frame, xy = TRUE)

## run/load simulation ------------------------------------------------------
# set parameters for species creation (these still get passed to 
# make_bio_recs_model() in order to record their values in the simulation file
# name for reproducability and to trigger a re-run of the simulation when these
# values change).
shape <- "hump"
polynomial.resp <- c(1, 2) 
polynomial.coef.min <- 0 
polynomial.coef.max <- 1.3 # 0.1 and 1.3 work. Setting lower limit to 0 lets there be essentially straight-line sp also, right?
sl.intercept.min <- NA 
sl.intercept.max <- NA # 0.1 and 0.7 work
sl.coef.min <- NA
sl.coef.max <- NA # TODO: refine slope. 0.7 and 1 might work?  values 
pca <- TRUE 
error.prob <- 0
env.predictors <- raw_pred_df


# create simulation if it doesn't yet exist.  Otherwise, load it.
if(new_draws) {
  ## simulate data
  set.seed(seed, "L'Ecuyer")
  # set number of simulations to draw within each index chunk
  n_sim_per_index <- rep(ceiling(nsim/n_cores), n_cores)
  # reduce index chunk sizes until the number of simulations equals nsim
  while(sum(n_sim_per_index) > nsim) { 
    n_sim_per_index[[which(n_sim_per_index == max(n_sim_per_index))[1]]] <- 
      n_sim_per_index[[which(n_sim_per_index == max(n_sim_per_index))[1]]] - 1
    if(any(n_sim_per_index == 0)) {
      n_sim_per_index <- n_sim_per_index[-which(n_sim_per_index == 0)]
      n_cores <- n_cores - 1
    }
  }
  
  sp_sim <- new_simulation(
    name = sim_name, 
    label = lab, 
    dir = sim_dest_dir) %>%
    generate_model(make_bio_recs_model, 
                   seed = seed,
                   community.size = community_size, 
                   n.obs = n_obs,
                   n.obs.reference = "community", 
                   shape = shape, 
                   polynomial.resp = polynomial.resp, 
                   polynomial.coef.min = polynomial.coef.min, 
                   polynomial.coef.max = polynomial.coef.max, 
                   sl.intercept.min = sl.intercept.min, 
                   sl.intercept.max = sl.intercept.max, 
                   sl.coef.min = sl.coef.min, sl.coef.max = sl.coef.max, 
                   pca = pca, 
                   bias.rasters = bias_rasters, 
                   bias.name = as.list(names(bias_rasters)), 
                   list.lengths = list_length, 
                   env.predictors = env.predictors, 
                   squared.predictors = sq_pred_df, 
                   error.prob = 0, 
                   randomPoints.replace = TRUE, 
                   on.sonic = on_sonic, 
                   vary_along = c("n.obs", "bias.name")) %>%
    # nsim will generate that many draws from each of n indices, and each chunk
    # specified in the index argument will get a separate RNG stream.  So the 
    # total number of unique simulated species will by nsim*n indices.  
    simulate_from_model(nsim = n_sim_per_index, 
                        index = 1:n_cores,
                        parallel = list(socket_names = n_cores,
                                        libraries = c("wgutil", "Hmisc", "rgdal",
                                                      "raster", "tidyverse",
                                                      "virtualspecies",
                                                      "simulator"))
    ) 
  save(sp_sim, file = paste0(sim_dest_dir, "files/", sim_name, 
                                "_draws.Rdata"))
} else sp_sim <- load_simulation(name = sim_name, 
                                    dir = sim_dest_dir)

if(new_methods) {
  ## loop through n_obs and run methods in a newly-named simulation just for 
  ## that value of n.obs.  Print the value of n.obs that has been run to the
  ## standard output.  This will run methods incrementally so that if wall time
  ## is exceeded I will know which draws the methods have been run on.  I can
  ## then load a simulation with references to all files in the 'files' 
  ## directory, which will give me the methods for all draws that have been run.
  #target_n <- which(n_obs == 100000) # specify n.obs to run
  for(i in 1:length(n_obs)) { # target_n
    print(paste0("\nStarting methods with n.obs = ", n_obs[[i]]))
    # run methods on draws for a single n.obs value
    sp_sim_sub <- subset_simulation(sp_sim, n.obs == n_obs[[i]]) %>%
      simulator::rename(paste0("sp_sim_methods_nobs", n_obs[[i]])) %>%
      run_method(list(glm_poly, glm_poly + wg_block_cv, 
                      rf, rf + wg_block_cv, 
                      brt, brt + wg_block_cv, 
                      idw_interp, idw_interp + wg_block_cv),
                 parallel = list(socket_names = n_cores,
                                 libraries = c("wgutil", "Hmisc", "rgdal",
                                               "raster", "tidyverse",
                                               "virtualspecies", 
                                               "randomForest",
                                               "dismo", 
                                               "pROC", 
                                               "gstat", 
                                               "unmarked", 
                                               "simulator"))) 
    print(paste0("\nFinished methods with n.obs = ", n_obs[[i]]))
  }
  # make the sp_sim object contain all methods that have been run
  sp_sim <- get_simulation_with_all_files(dir = sim_dest_dir) %>%
    simulator::rename(name = sim_name)
  save_simulation(sp_sim)
  save(sp_sim, file = paste0(sim_dest_dir, "files/", sim_name, 
                                "_outputs.Rdata"))
} else sp_sim <- load_simulation(name = sim_name, 
                                 dir = sim_dest_dir)
  
if(new_evals) {
  sp_sim <- evaluate(sp_sim, list(auc, rmse)) # evaluate(sp_sim, list(rmse))
} else sp_sim <- load_simulation(name = sim_name, 
                                 dir = sim_dest_dir)

## end run/load simulation --------------------------------------------------


### analyze performance --------------------------------------------------------
if(analyze_performance) {
  library(mgcv)
  library(dismo)
  source("./analyze_performance_BRT.R")
  source("./analyze_performance_GAM.R")
  source("./analyze_performance_supplementary.R")
}

if(print_plots) {
  source("plots_main_text.R")
  source("plots_supplementary.R")
}
### end analyze performance ----------------------------------------------------

options("scipen" = 0, "digits" = 7)
print(sessionInfo())