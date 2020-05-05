############################################
## This is the main simulator file for the species simulation
## This is for simulation 10.11
## Using bryophytes as template
## In this version (10.11) I have fixed the variable selection for GLMs
## (see brown notebook 24 March 2020)
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
## 
## author: Willson Gaul wgaul@hotmail.com
## created: 24 Sep 2018
## last modified: 5 May 2020
############################################

rm(list = ls())
library(wgutil)
library(Hmisc)
library(rgdal)
library(raster)
#library(blockCV)
library(parallel)
# library(MuMIn) # (needed for step in occu) conflicts with randomForest
# library(randomForest) # conflicts with MuMIn
library(dismo)
library(pROC)
library(gstat)
library(unmarked)
library(tidyverse)
library(virtualspecies)
library(simulator) # this file was created under simulator version 0.2.0

new_draws <- F  # run the simulation again because parameters have been changed
new_methods <- F # run the SDM modelling methods 
new_evals <- F # calculate AUC and RMSE for SDMs
add_to_sim <- F # add draws or models to an existing simulation
extract_truth_maps <- F # get truth maps out of simulated object and save as a raster brick
analyze_performance <- T # analyze the effect of sampling bias, sample size, etc.
print_plots <- T # print plots to .jpg and .pdf files
run_tests <- F # print some diagnostic tests to the screen
dbg <- F # run with some debugging options
on_laptop <- T # is this running on WG's laptop?
on_sonic <- F # is this running on sonic?

n_cores <- 1
seed <- 10061983 + 12202018 + 1300  # wg bday + today's date + current time

options("scipen" = 100, "digits" = 4) # so that file names don't contain scientific notation

if(on_laptop) {
  sim_dest_dir <- "./"
  
  ### prepare and load environmental data ---------------------------------------
  if(!all(file.exists("elevation_hec_ETOPO1.rds") & 
          file.exists("annual_precip_hectad.rds") & 
          file.exists("summer_tx_hectad.rds") & 
          file.exists("winter_tn_hectad.rds") &
          file.exists("mean_pp_hectad.rds") &
          file.exists("corine_label_1_hectad.RData"))) {
    source("interp_elev_hec_etopo1.R")
    source("eobs_annual_precipitation_hectad.R")
    source("eobs_max_summer_temp_hec.R")
    source("eobs_min_winter_temp_hec.R")
    source("eobs_pp_hectad.R")
    source("prep_corine.R")
  }
  
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
  
  if(!all(file.exists("spat_nrec_hec.rds") & 
          file.exists("community_size_and_list_lengths.rds"))) {
    # item to use from here is tab_nrec_hec
    # This uses NBDC data that is not stored in this simulation directory
    # After this is run once to create the .rds file, the simulation code should
    # run without needing access to the NBDC data (which is restricted)
    source("format_for_sampling_distribution_plotting.R")
    source("raster_nrec_sampling_weights.R")
    source("define_community_size_and_list_lengths.R")
  }
  # load sample bias rasters
  spat_nrec_hec <- readRDS("spat_nrec_hec.rds")
  # load data for defining list lengths and community size
  community_list_lengths <- readRDS("community_size_and_list_lengths.rds")
}
if(on_sonic) {
  setwd("~/scratch/sims_10.11/")
  sim_dest_dir <- "./"
  
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
  # load sample bias rasters
  spat_nrec_hec <- readRDS("spat_nrec_hec.rds")
  # load data for defining list lengths and community size
  community_list_lengths <- readRDS("community_size_and_list_lengths.rds")
  
  if(dbg) print(paste0("number of cores: ", n_cores))
}

source("./simulation_functions.R")
source("./utility_functions.R")
source("./model_functions.R")
source("./method_functions.R")
source("./eval_functions.R")


## set simulation parameters --------------------------------------------------
sim_name <- "bryophyte_template"
lab <- "Bryophyte template"
community_size <- community_list_lengths$bry_BBS$community_size
list_length <- community_list_lengths$bry_BBS$list_lengths_std
n_obs <- c(2, 5, 10, 50, 100, 200, 1000, 5000) * community_size #
nsim <- 110 # number of species to simulate.  Ideally divisible by n_cores.
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
# make a no-bias raster using one of the existing rasters as a template
no_bias_rast <- spat_nrec_hec$bird$rast
no_bias_rast <- calc(no_bias_rast, fun = function(x) {
  n_cell <- length(x[!is.na(x)])
  x[!is.na(x)] <- 1/n_cell
  x
})

bias_rasters <- list(
  no_bias = no_bias_rast, 
  extreme_mothNBDC = spat_nrec_hec$`insect - moth`$rast, 
  median_butterflyNBDC = spat_nrec_hec$`insect - butterfly`$rast, 
  least_birdNBDC = spat_nrec_hec$bird$rast)
## TODO: 28 Sep.  Moved this next bit into simulation function for now as it 
## needs to be done after I re-create the rasters in the function.  Do this more
## elegantly if I decide to keep it that way.  Alternatively, un-comment this
## to have it done here again.
# bias_rasters <- lapply(bias_rasters, FUN = function(x, ir_TM75) {
#   raster::mask(x, mask = ir_TM75)}, ir_TM75 = ir_TM75)
## end specify bias rasters -------------------------------------------------

## combine all predictors into a single df ------------------------------------
## TODO: masking by ir_TM75 throws out some hectads that do have records in the
## real NBDC data I think (reducing total hectads from around 1003 to 840).  Do
## I want to somehow keep all 1003 grid cells?
pred_brick <- brick(list(
  elev = resample(elev_hec, krg_mean_rr_rast), 
  mean_rr = krg_mean_rr_rast, 
  mean_tx = resample(krg_mean_tx_rast, krg_mean_rr_rast), 
  mean_tn = resample(krg_mean_tn_rast, krg_mean_rr_rast), 
  mean_pp = resample(mean_pp_rast, krg_mean_rr_rast), 
  agricultural = resample(agricultural_l1_rast, krg_mean_rr_rast), 
  artificial_surfaces = resample(artificial_surfaces_l1_rast, krg_mean_rr_rast), 
  forest_seminatural = resample(forest_seminatural_l1_rast, krg_mean_rr_rast),
  water = resample(water_l1_rast, krg_mean_rr_rast), 
  wetlands = resample(wetlands_l1_rast, krg_mean_rr_rast)
  ))
# mask pred brick by ir_TM75 to get only Irish land cells
pred_brick <- mask(pred_brick, ir_TM75)
# scale and centre environmental predictors over study extent
pred_brick <- scale(pred_brick, center = TRUE, scale = TRUE)

pred_brick_sq <- pred_brick^2 # square predictors
names(pred_brick_sq) <- paste0(names(pred_brick), "_sq")

raw_pred_df <- left_join(data.frame(as(pred_brick, "SpatialGridDataFrame")),
                         data.frame(as(pred_brick_sq, "SpatialGridDataFrame")), 
                         by = c("s1", "s2"))

raw_pred_df$s1 <- round(raw_pred_df$s1) # fix number representation issues
raw_pred_df$s2 <- round(raw_pred_df$s2)

# reorder columns to put "s1" and "s2" (grid coordinate) columns first
new_order <- names(raw_pred_df)
new_order <- new_order[c(which(new_order %in% c("s1", "s2")), 
                         which(new_order %nin% c("s1", "s2")))]
raw_pred_df <- raw_pred_df[, new_order]

# remove original predictors to clean workspace
rm(elev_hec, krg_mean_rr_rast, krg_mean_tn_rast, krg_mean_tx_rast, 
   agricultural_l1_rast, artificial_surfaces_l1_rast, clc_l1_props_hecs, 
   forest_seminatural_l1_rast, water_l1_rast, wetlands_l1_rast, 
   new_order, mean_pp_rast, pred_brick, pred_brick_sq)
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

## @knitr main ---------------------------------------------------------------
## TODO: 28 Sep. Converting bias rasters to data frames so sha1() can deal 
## with them. This is just a quick fix to let the simulator process run.  
## If I end up keeping this change, will probably need to pass all the 
## information for converting it back to a raster as arguments so the dfs 
## can be turned back in to rasters inside the simulation functions. 
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
                   n.obs = as.list(n_obs),
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
                   env.predictors = raw_pred_df[, which(
                     !grepl(".*_sq", colnames(raw_pred_df))), ], 
                   predictors.for.models = raw_pred_df, 
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
                                        libraries = c("wgutil", "Hmisc", 
                                                      "rgdal",
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
      run_method(list(brt, brt + wg_block_cv),
                 parallel = list(socket_names = n_cores,
                                 libraries = c("wgutil", "Hmisc", "rgdal",
                                               "raster", "tidyverse",
                                               "virtualspecies", 
                                               #"randomForest", 
                                               "dismo", 
                                               "pROC", 
                                               "gstat", 
                                               "unmarked", 
                                               #"MuMIn", # masks from randomForest
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
  # save(sp_sim, file = paste0(sim_dest_dir, "files/", sim_name, 
  #                               "_evals.Rdata"))
  # # write a .csv of performance results to use for making plots
  # evals_df <- as.data.frame(evals(sp_sim))
  # print(paste0("Size of model performance data frame: ", 
  #       pryr::object_size(evals_df)))
  # write_csv(evals_df, paste0("model_performance_results_", 
  #                            sim_name, ".csv"))
} else sp_sim <- load_simulation(name = sim_name, 
                                 dir = sim_dest_dir)

## end run/load simulation --------------------------------------------------

if(extract_truth_maps) {
  ## Get truth maps
  # get draws
  if(on_laptop) {
    sp_sim <- load_simulation(name = sim_name, 
                              dir = "/media/willson/Backup Plus/simulation_data/sims_10.11")
    sp_draws <- draws(sp_sim, n.obs == 6340)
  }
  if(on_sonic) {
    sp_draws <- draws(sp_sim, n.obs == 6340)
  }
  
  truth_maps <- lapply(sp_draws[[1]]@draws, 
                       FUN = function(x) x$truth_maps[[1]])
  truth_maps <- brick(truth_maps)
  try(saveRDS(truth_maps, "truth_maps.rds"))
  
  ## Get an example sample
  sp_draws <- draws(sp_sim, n.obs == 126800)
  # get butterfly bias draws
  bb_sim <- load_simulation(name = "butterfly_template",
                            dir = "/media/willson/Backup Plus/simulation_data/sims_10.11_butterfly/")
  
  example_obs_n63400 <- list(no_bias = sp_draws[[1]]@draws$r20.3$observations, 
                           median_bias = sp_draws[[3]]@draws$r20.3$observations,
                           extreme_bias = sp_draws[[2]]@draws$r20.3$observations)
  try(saveRDS(example_obs_n63400, "example_obs_n63400.rds"))
  try(rm(truth_maps, example_obs_n63400))
  
  ## Get block cv fold assignments
  ex_refs <- unlist(output(bb_sim, methods = "glm_poly_wg_block_cv", 
                           reference = T))
  ex_outs <- output(bb_sim, subset = c(ex_refs[[1]]@model_name), 
                     methods = "glm_poly_wg_block_cv", index = 1)
  try(saveRDS(ex_outs@out$r1.1$folds_en, "example_folds.rds"))
}



### testing code -------------------------------------------------------------
if(run_tests) {
  # make sure all n.obs finished running - each draws ref should have the
  # same length (which should be the number of cores)
  for(i in 1:length(sim@draws_refs)) {
    print(paste0(i, " ", length(sim@draws_refs[[i]])))
  }
  
  # using subset is same as just using evals function, but perhaps faster
  identical(evals(subset_simulation(sp_sim, n.obs == 1000)), 
            evals(sp_sim, n.obs == 1000))
  
  # expect nrow to be nsim*length(bias_rasters)*number of n.obs in this subset 
  # (may be slightly higher if nsim wasn't evenly divisible by n_cores)
  dim(test_eval_df) 
  names(test_evals[[1]]@evals$glm_straight_line)
  test_eval_df[c(1, 13), ] 
  identical(test_eval_df$Model[1], test_eval_df$Model[13])
  
  draw_test <- draws(sp_sim, n.obs == 6340)
  length(draw_test)
  # are the objects with same draw name the same virtualspecies?  They should be
  # if the 4 draws are the 4 biases (can tell by printing draw_test and looking
  # at the name)
  draw_test[[1]]@draws$r1.1 
  draw_test[[2]]@draws$r1.1
  
  # are the objects with the same number within each index iteration the same
  # virtualspecies?  They should NOT be.
  draw_test[[1]]@draws$r1.1
  draw_test[[1]]@draws$r2.1
  
  draw_test_all <- draws(sp_sim_test)
  length(draw_test_all)
  draw_test_all[[1]]@draws$r1.1
  draw_test_all[[7]]@draws$r1.1
  draw_test_all[[5]]@draws$r1.1
  draw_test_all[[5]]@draws$r2.1
  
}
### end testing code ---------------------------------------------------------

### analyze performance -------------------------------------------------------
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
