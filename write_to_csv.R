#####################################
## This script is for use on sonic to write results from the simulation to .csv 
## 
## author: Willson Gaul wgaul@hotmail.com
## created: 14 Nov 2018
## last modified: 27 March 2020
#####################################

rm(list = ls())
library(wgutil)
#library(skimr)
library(Hmisc)
library(rgdal)
library(raster)
library(parallel)
library(tidyverse)
library(virtualspecies)
library(simulator) # this file was created under simulator version 0.2.0

add_to_sim <- F # add draws or models to an existing simulation
dbg <- T
on_sonic <- T
print_plots <- F
n_cores <- 1

options("scipen" = 100, "digits" = 4) # so that file names don't contain scientific notation

if(!on_sonic) {
  sim_dest_dir <- "."
}
if(on_sonic) {
  sim_dest_dir <- "./"
}

source("./simulation_functions.R")
source("./utility_functions.R")
source("./model_functions.R")
source("./method_functions.R")
source("./eval_functions.R")

sim_name <- "bryophyte_template"

sim <- load_simulation(name = sim_name, dir = sim_dest_dir)
#sim <- subset_simulation(sim, methods = "glm_poly_cv")


bias_names <- c("no_bias", "extreme_mothNBDC", "median_butterflyNBDC",
                "least_birdNBDC")

# extract community size for this simulation (for calculating nobs per species)
community.size <- sim@model_refs[[1]]@name
community.size <- gsub(".*community.size_", "", community.size)
community.size <- gsub("/env.predi.*", "", community.size)
community.size <- as.numeric(community.size)

n_obs <- c(2, 5, 10, 50, 100, 200) * community.size

### get data frame of AUC values ----------------------------------------
ev_df <- sim %>% evals %>% as.data.frame
model_df <- sim %>% model %>% as.data.frame
ev_df <- dplyr::right_join(model_df, ev_df, by = c("name" = "Model"))

# get prevalence of each species
ds <- draws(sim, n.obs == n_obs[1])
prev_df <- data.frame(matrix(nrow = length(ds[[1]]@draws), ncol = 2))
colnames(prev_df) <- c("Draw", "prevalence")
for(i in 1:length(ds[[1]]@draws)) {
  prev_df$Draw[i] <- names(ds[[1]]@draws)[i]
  tb <- table(ds[[1]]@draws[[i]]$truth_maps[[1]]@data@values)
  prev_df$prevalence[i] <- tb["1"] / sum(tb)
}

ev_df <- left_join(ev_df, prev_df, by = "Draw")

# get number of observations per species.  This should be an average obtained by
# dividing n.obs for the taxonomic group by the number of species.  NOT a sum of
# giving the number of positive observations for each species (number of 
# presences).  This is because the glm method and others assume absences and so
# effectively have data points for each species for each checklist.  So the
# number of observations per species is the number of checklists, each of which 
# can give either a 1 or a 0 for each species.  

# n_bias_df$n_obs_per_sp <- n_bias_df$n_obs/community.size

ev_df$n_obs_per_sp <- ev_df$n.obs/community.size

## end make df of AUC values -------------------------------------------------
print(pryr::object_size(ev_df))
write_csv(ev_df, paste0(sim_name, "_evals.csv"))
