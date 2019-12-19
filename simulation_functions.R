##############################
## Define functions for simulating species and drawing samples from the 
## simulated distributions.  
## Simulation 10.10
## 
##
## notes: - because I need list length, but don't want to be able to simulate
##            large numbers of species without making the total possible number
##            of species unrealistically large (for purposes of sampling by 
##            list length), I need to create multiple smaller communities to get
##            the total number of simulated species that I desire. 
##
## author: Willson Gaul wgaul@hotmail.com
## created: 24 Sep 2018
## last modified: 1 Nov 2018
##############################

simulate_fun <- function(community.size = NULL, n.obs = NULL, 
                         n.obs.reference = NULL, 
                         shape = NULL, polynomial.resp = NULL, 
                         polynomial.coef.min = NULL, polynomial.coef.max = NULL, 
                         sl.intercept.min = NULL, sl.intercept.max = NULL, 
                         sl.coef.min = NULL, sl.coef.max = NULL, 
                         pca = TRUE, bias.rasters = NULL, bias.name = NULL, 
                         list.lengths = NULL, 
                         env.predictors = NULL, squared.predictors = NULL, 
                         error.prob = 0, randomPoints.replace = NULL, 
                         on.sonic = NULL, nsim = NULL) {
  # this function must return a list of length nsim where each element is a df 
  # of all checklists for the community, with 
  # a row for each checklist and columns for grid cell, checklist ID, 
  # list length, and the observed/not observed status of the focal species (so
  # a different species for each element of the list).  This will still preserve
  # the list length and assumed non-detection info but without duplicating the
  # observations for all species inside the result returned for each individual
  # species.  This should save a lot of memory. 

  # tif path for use as irish_hec_raster in simulate_fun
  if(on_sonic) tif.path <- "../data/eobs/IE_10km_hecs_raster.tif" 
  if(!on_sonic) tif.path <- "~/Documents/Data_Analysis/UCD/predictor_variables/eobs/IE_10km_hecs_raster.tif"

  ## restore bias.raster to raster
  bias.raster <- bias.rasters[[which(names(bias.rasters) == bias.name)]]
  
  if(any(!is.raster(bias.raster))) {
    warning("bias.rasters are not actually rasters because of hash issue.  Quick fix implemented here, but if I keep this I should pass all the relevant parts of this as arguments into simulate_fun().")
    # template raster
    irish_hec_raster <- raster::raster(tif.path)
    irish_hec_raster <- projectRaster(from = irish_hec_raster, 
                                      crs = CRS("+init=epsg:29903"))
    bias.raster <- make_spatial(bias.raster, temp_rast = irish_hec_raster, 
                                field = "layer")
    bias.raster <- raster::mask(bias.raster, mask = ir_TM75)
    ## TODO: if I keep it this way, need to pass in ir_TM75
  }
  
  ## create species list - this should be different each time -----------------
  n_sp <- nsim
  while(n_sp %% community.size != 0) {
    n_sp <- n_sp + 1
  }
  
  # make all species
  # this creates as many species as needed to make full communities containing
  # at least nsim species total
  sp_list <- make_all_species(n.sp = n_sp, shape = shape, 
                              polynomial.resp = polynomial.resp,
                              polynomial.coef.min = polynomial.coef.min, 
                              polynomial.coef.max = polynomial.coef.max,
                              sl.intercept.min = sl.intercept.min, 
                              sl.intercept.max = sl.intercept.max, 
                              sl.coef.min = sl.coef.min, 
                              sl.coef.max = sl.coef.max, pca = pca, 
                              env.predictors = env.predictors,
                              error.prob = error.prob, 
                              community.size = community.size)
  
  ## group species into communities
  sp_list <- group_communities(community.size = community.size, 
                               sp.list = sp_list) 
  sp.list <- sp_list
  ## end create species list -------------------------------------------------
  
  
  
  ## prepare bias raster
  bias.raster <- raster::crop(bias.raster, sp.list[[1]]$truth_maps[[1]])
  # calculate probability of sampling for each cell based on the current values 
  # which are number of records per hectad normalized from 0 to 1 for the dataset.
  calc_proportion <- function(x) {
    x[!is.na(x)] <- x[!is.na(x)]/sum(x[!is.na(x)])
    x
  }
  bias.raster <- calc(bias.raster, fun = calc_proportion)
  
  ## sample each community to make a df of obs for that community
  community_dfs <- sample_all_communities(sp.list = sp.list, 
                                          community.size = community.size, 
                                          n.obs = n.obs, 
                                          n.obs.reference = n.obs.reference, 
                                          list.lengths = list.lengths, 
                                          bias.raster = bias.raster, 
                                          randomPoints.replace = 
                                            randomPoints.replace)

  # select only the number of species desired for simulation
  sp.list <- sp.list[1:nsim] 
  # extract from the community observation dfs the unique checklist IDs, the
  # list lengths, and the observed/not observed status for each species
  for(i in 1:length(sp.list)) {
    d <- community_dfs[[sp.list[[i]]$community]]
    d_ll <- count(d, checklist_ID) # calculate list length
    colnames(d_ll)[2] <- "list_length"
    d <- tryCatch({
      left_join(d, d_ll, by = "checklist_ID") %>% # add list length to sp df
        mutate(present = 1) %>% 
        spread(sp_name, value = present, fill = 0) %>%
        select(x, y, checklist_ID, list_length, eval(names(sp.list)[i]))}, 
      error = function(x) NA)
    sp.list[[i]]$observations <- d
  }
  sp.list # return the list of length nsim with a df for each sp.
}

make_all_species <- function(shape = shape, n.sp = NULL, ...) {
  # This function calls make_sp_response() as many times as necessary to create
  # complete communities that hold at least nsim species.
  #
  # ARGS  n.sp - integer giving the number of species to create
  #       shape - character vector of length n.sp giving the response shapes for
  #               all species that are to be created
  #       ... - (optional) values to be passed to make_sp_response
  # 
  # VALUE a list of virtualspecies objects
  list2env(list(...), envir = environment())
  
  if(is.null(shape) || !is.character(shape)) stop("A character vector must be provided to 'shape' in make_all_species().")
  if(is.null(n.sp) || !is.numeric(n.sp)) stop("An integer value must be provided to 'n.sp' in make_all_species().")
  
  ## make all species
  sp_names <- paste0(rep("sp", n.sp), 1:n.sp) # make species names
  sp <- lapply(sp_names, FUN = make_sp_response, 
               env.predictors = env.predictors, 
               pca = pca, shape = shape, polynomial.resp = polynomial.resp, 
               polynomial.coef.min = polynomial.coef.min, 
               polynomial.coef.max = polynomial.coef.max, 
               sl.intercept.min = sl.intercept.min, 
               sl.intercept.max = sl.intercept.max, 
               sl.coef.min = sl.coef.min, 
               sl.coef.max = sl.coef.max, 
               community.size = community.size)
  names(sp) <- sp_names
  
  ## make true distributions
  # only use mclapply if simulator is not using parallel
  # sp <- mclapply(sp, FUN = call_make_truth_map, n = 1,
  #                mc.cores = detectCores())
  sp <- lapply(sp, FUN = call_make_truth_map, n = 1)
  sp # return list of virtual species
}

group_communities <- function(community.size, sp.list) {
  # This function groups a list of virtual species into communities.
  # ARGS  community.size: The total number of species in the taxonomic group
  #       sp.list: a list of virtualspecies objects defining each species
  #
  # VALUE a list of virtualspecies objects (one for each species), with an 
  #       additional element named "community" that holds a character string
  #       identifying which community the species has been assigned to.

  ## group species into communities
  communities <- c() # make community labels
  for(j in 1:(length(sp.list)/community.size)) {
    communities <- c(communities, paste0(rep("community_", community.size), j))
  }
  # assign community labels to species
  for(j in 1:length(sp.list)) {
    sp.list[[j]]$community <- communities[j]
  }
  sp.list # return species list with community groups
}

make_sp_response <- function(sp.name = NULL, env.predictors = NULL, pca = TRUE, 
                             shape = NULL, polynomial.resp = NULL, 
                             polynomial.coef.min = NULL, 
                             polynomial.coef.max = NULL, 
                             sl.intercept.min = NULL, 
                             sl.intercept.max = NULL, 
                             sl.coef.min = NULL, 
                             sl.coef.max = NULL, ...) {
  # Function to make an individual species
  # Response functions are given here in terms of the linear predictors, so the
  # coefficients will be the change in the logit(probability.of.occurence) for a 
  # 1-unit increase in predicictor x.  The output raster will therefore also be
  # in the units of the linear predictors, in this case logit(prob), and will
  # need to be transformed with logistic(raster) to get a raster giving 
  # probabilities of occurrence in each cell.
  #
  # So, the equation is: ??? TODO ??? 
  # ? logit(prob.of.occurrence) = intercept + beta1*env_var1[i] + beta2*env_var2[i] ??
  #
  # ARGS: sp.name - a single species name
  #       env.predictors - data frame holding raw predictor variable values for 
  #         each grid cell.  This will be used for PCA if pca = TRUE.
  #       pca - logical indicating whether to calculate 2 principal components
  #         from env.predictors to use for the 2 predictor variables.  If pca = F, 
  #         the first two predictor variables in env.predictors will be used as the
  #         predictors.
  #       shape - character string. "hump" or "straight line" indicating the 
  #         response shape for the species
  #       polynomial.resp - numeric vector giving the positions of the predictors
  #         for which a polynomial (squared) response should be included. e.g. 
  #         polynomial.resp = c(1, 2) will produce a polynomial response for
  #         both predictor variables.  This argument is only valid
  #         if the sp is a hump species.  
  #       polynomial.coef.min - (optional) numeric, the minimum value for the 
  #               polynomial coefficients
  #       polynomial.coef.max - (optional) the maximum value for the polynomial 
  #               coefficients
  #       sl.intercept.min - (optional) the minimum value for the straight-line 
  #               intercept
  #       sl.intercept.max - (optional) the maximum value for the straight-line 
  #               intercept
  #       sl.coef.min - (optional) the minimum value for the straight-line 
  #               coefficients
  #       sl.coef.max - (optional) the maximum value for the straight-line 
  #               coefficients
  #       
  # VALUE: object of class "virtualspecies" containing rasters with values for
  #         both logit(prob) and probability of occurrence.
  list2env(list(...), envir = environment())
  
  ## error checking --------------------------------------------------------
  if(is.null(env.predictors)) stop("Must provide a data frame with at least one predictor variable.")
  if(!is.null(polynomial.resp) && is.null(polynomial.coef.min) || 
     is.null(polynomial.coef.max)) {
    warning("Argument provided to polynomial.resp but no polynomial coefficients provided.")
  }
  if(!is.null(polynomial.resp)) {
    if(!is.numeric(polynomial.resp)) {
      stop("polynomial.resp must be either numeric.")
    }
    if(any(polynomial.resp > 2)) {
      stop("polynomial.resp must be 1, 2, or c(1, 2).  Species with more than two predictors are not implemented right now.")
    }
  }
  # end error checking -------------------------------------------------
  
  pred_list <- list(pred_1 = list(pred = NULL, coef = NULL, intercept = NULL), 
                    pred_2 = list(pred = NULL, coef = NULL, intercept = NULL), 
                    pred1_sq = list(pred = NULL, coef = NULL), 
                    pred2_sq = list(pred = NULL, coef = NULL))
  
  if(pca == T) {
    ## run PCA
    # randomly drop 1/3 of preditor variables each time
    drop_cols <- c(1, 2, sample(3:ncol(env.predictors), 
                                size = trunc((ncol(env.predictors)-2) / 3)))
    pred_pca <- ade4::dudi.pca(env.predictors[, -drop_cols], 
                               center = T, scale = T, 
                               scannf = F, nf = 2)
    
    # ade4::s.arrow(pred_pca$c1, lab = names(pred_pca$tab))
    # ade4::s.corcircle(pred_pca$co, lab = names(pred_pca$tab), full = T, box = T)
    # ade4::scatter(pred_pca)
    
    # make spatialDataFrame with principal component values
    pca_df <- data.frame(pred_pca$li)
    coordinates(pca_df) <- as.matrix(env.predictors[, c("s1", "s2")])
    proj4string(pca_df) <- CRS("+init=epsg:29903")
    gridded(pca_df) <- TRUE # promote to SpatialPixelsDataFrame
    pca_df <- as(pca_df, "SpatialGridDataFrame")
    
    # Make rasters of each PC
    pca1_rast <- raster::raster(pca_df["Axis1"])
    pca2_rast <- raster::raster(pca_df["Axis2"])
  } else {
    stop("Have not yet implemented option to use the first 2 columns of env.predictors as the 2 predictors.")
    # here I would assign the 2 first columns of env.predictors to pca1_rast and
    # pca2_rast
  }
  
  pred_list$pred_1$pred <- pca1_rast # store pc 1 raster
  pred_list$pred_2$pred <- pca2_rast # store pc 2 raster
  
  theoretical_prevalence <- 0
  
  while(theoretical_prevalence < 0.01) {
    if(shape == "hump") {
      ## TODO: 1 Oct. Do I need to center pca variables to use for generating 
      ## responses?
      # make squared predictor raster(s)
      if(1 %in% polynomial.resp) {
        pred_list$pred1_sq$pred <- pca1_rast^2
      } 
      ## TODO: 25 Sep. bug - if only 1 variabel is to have a polynomical term, 
      ## then this doesn't create the
      ## second raster and then the raster stack has the wrong number of layers.
      if(2 %in% polynomial.resp) { 
        pred_list$pred2_sq$pred <- pca2_rast^2
      } else warning("if only 1 variabel is to have a polynomical term, then this doesn't create the second raster and then the raster stack has the wrong number of layers.")
    }
    
    ## set coefficients for each predictor term  --------------------------
    if(shape == "hump") {
      ## set squared term coefficients
      # For squared term, set randomly.  For intercept and 1st degree term, set
      # the location of the vertex randomly within a defined range, then calculate
      # intercept and 1st degree coefficient.
      # Bigger (>1) coefficients give narrower hump, smaller (<1) coefficients 
      # give wide hump.
      # Use only negative coefficients for hump species, as a positive coefficient
      # would turn the hump into a U
      # TODO: refin min and max values to get reasonable coefficients
      
      pred_list$pred1_sq$coef <- -1*runif(n = 1, 
                                          min = polynomial.coef.min, 
                                          max = polynomial.coef.max)
      pred_list$pred2_sq$coef <- -1*runif(n = 1, 
                                          min = polynomial.coef.min, 
                                          max = polynomial.coef.max)
      
      # set vertex coordinates
      v_xmin <- quantile(pred_list$pred_1$pred@data@values, probs = 0.1, 
                         na.rm = T) #TODO: get this right - I don't want 
      # quantile of the values, I want quantile of the range of values
      v_xmax <- quantile(pred_list$pred_1$pred@data@values, probs = 0.9, 
                         na.rm = T)
      v_ymin <- logit(0.3) 
      v_ymax <- logit(0.8)
      v_y <- runif(n= 1, min = v_ymin, max = v_ymax) # set y location of vertex
      v_x <- runif(n = 1, min = v_xmin, max = v_xmax) # set x location of vertex
      
      ## set coefficients for 1st order terms (based on coefs for squared terms 
      # and vertex location that were set above)
      # calculate coefficients based on vertex location
      coefs <- as.list(calc_coef_from_vertex(vert_x = v_x, vert_y = v_y, 
                                             a = pred_list$pred1_sq$coef))
      pred_list$pred_1$coef <- coefs$b # set coefficient for pred1 term
      pred_list$pred_1$intercept <- coefs$c # set intercept for pred1 response
      
      if(!is.null(pred_list$pred_2$pred)) {
        ## coefficients for response to predictor 2
        # determine vertex x-coordinate
        v_xmin <- min(pred_list$pred_2$pred@data@values, na.rm = T)
        v_xmax <- max(pred_list$pred_2$pred@data@values, na.rm = T)
        v_ymin <- logit(0.2) 
        v_ymax <- logit(0.7)
        v_y <- runif(n= 1, min = v_ymin, max = v_ymax) # set y location of vertex
        v_x <- runif(n = 1, min = v_xmin, max = v_xmax) # set x location of vertex
        
        # calculate coefficients based on vertex location
        coefs <- as.list(calc_coef_from_vertex(vert_x = v_x, vert_y = v_y, 
                                               a = pred_list$pred2_sq$coef))
        pred_list$pred_2$coef <- coefs$b # set coefficient for pred2 term
        pred_list$pred_2$intercept <- coefs$c # set intercept for pred2 response
      }
    }
    
    if(shape == "straight line") {
      ## set intercept (baseline probability of occurrence)
      # TODO: relate these prob_occ intercept values to the range for these values
      #       in hump sp.  Perhaps median of that range?
      # TODO: confirm that lm recovers the correct parameter values
      prob_occ <- runif(min = sl.intercept.min, 
                        max = sl.intercept.max, n = 1)
      
      pred_list$pred_1$intercept <- logit(prob_occ)
      pred_list$pred_2$intercept <- logit(prob_occ)
      
      ## set slopes
      pred_list$pred_1$coef <- sample(c(1, -1), size = 1)*
        runif(n = 1, min = sl.coef.min, max = sl.coef.max)
      pred_list$pred_2$coef <- sample(c(1, -1), size = 1)*
        runif(n = 1, min = sl.coef.min, max = sl.coef.max)
    }
    
    # make raster stack
    to_stack <- list() # list to hold rasters
    for(i in 1:length(pred_list)) {
      if(!is.null(pred_list[[i]]$pred)) {
        to_stack[[i]] <- pred_list[[i]]$pred
        names(to_stack)[i] <- names(pred_list)[i]
      }
    }
    keep <- !sapply(to_stack, FUN = is.null)
    to_stack <- to_stack[keep]
    var_stack <- stack(to_stack)
    names(var_stack) <- names(to_stack)
    
    if(shape == "hump") {
      # set intercept to 0 for all but the first variable, because all the linear
      # functions are added to create the final response.
      s <- generateSpFromFun(var_stack, 
                             parameters = list(
                               pred_1 = list(fun = "linearFun", 
                                             args = list(
                                               a = pred_list$pred_1$coef, 
                                               b = pred_list$pred_1$intercept)), 
                               pred_2 = list(fun = "linearFun",
                                             args = list(
                                               a = pred_list$pred_2$coef,
                                               # TODO: here 5 Feb!!!!!!  Need to
                                               # carefully evaluate whether adding a different 
                                               # intercept here changes 
                                               # parameter values for coefficients or only changes 
                                               # the overall intercept for 
                                               # the combined response.  So far it looks like glm 
                                               # recovers the parameter 
                                               # values exactly as I set them, even with this 
                                               # combined intercept.
                                               b = pred_list$pred_2$intercept)), 
                               pred1_sq = list(fun = "linearFun", 
                                               args = list(
                                                 a = pred_list$pred1_sq$coef, 
                                                 b = 0)), 
                               pred2_sq = list(fun = "linearFun",
                                               args = list(
                                                 a = pred_list$pred2_sq$coef,
                                                 b = 0))),
                             rescale = F, species.type = "additive", 
                             rescale.each.response = F, plot = F)
    }
    if(shape == "straight line") {
      s <- generateSpFromFun(var_stack, 
                             parameters = list(
                               pred_1 = list(fun = "linearFun", 
                                             args = list(
                                               a = pred_list$pred_1$coef, 
                                               b = pred_list$pred_1$intercept)), 
                               pred_2 = list(fun = "linearFun", 
                                             args = list(
                                               a = pred_list$pred_2$coef,
                                               b = 0))), 
                             rescale = F, species.type = "additive", 
                             rescale.each.response = F, plot = F)
    }
    
    s$prob.raster <- calc(s$suitab.raster, fun = logistic)
    # calculate likely prevalence by summing probabilities and dividing by
    # number of cells
    s$theoretical_prevalence <- sum(s$prob.raster@data@values, na.rm = T) / 
      length(which(!is.na(s$prob.raster@data@values) == T))
    theoretical_prevalence <- s$theoretical_prevalence # for while loop
    
    s$wg_params <- pred_list
    s$pca_dudi <- pred_pca
    s$pca_df <- pca_df
    s$community_size <- community.size
    s$sp_name <- sp.name
    # TODO: can only set prob_occ if I use the same intercept for both predictors
    # s$wg_params$prob_occ <- logistic(pred_list$pred_1$intercept)
  } # end while loop
  
  s # return the final object
}
### end make_sp_response -------------------------------------------------------


# list2env(sim_params, envir = environment()) # TODO what is this doing and where does it belong?

## Functions to Create True Distribution Maps -------------------------------
make_truth_map <- function(prob_occ_rast) {
  ## make a raster of realized p/a values based on probability of occurrence
  # ARGS: prob_occ_rast - a raster with each cell giving a probability of 
  #         occurrence for a species in that cell
  # RETURN: a raster giving a realized species p/a (1/0) for each cell.
  tm <- calc(prob_occ_rast, fun = function(x) {
    x[!is.na(x)] <- rbinom(n = length(which(!is.na(x))), size = 1, 
                           prob = x[!is.na(x)])
    x})
  tm
}

call_make_truth_map <- function(sp_list_element, n = 1) {
  # ARGS: sp_list_element - a virtualspecies object from the list called "sp" in
  #                         my simulation
  #       n - integer number of truth maps to create
  # RETURN: sp_list_element with an additional elemtent, which is a list of 
  #         n different realizations of truth maps.
  sp_list_element$truth_maps <- list()
  for (i in 1:n) {
    sp_list_element$truth_maps[[i]] <- make_truth_map(sp_list_element$prob.raster)
  }
  sp_list_element
}
## end functions to create truth maps ----------------------------------------

#### sample communities -------------------------------------------------------
## function to be called inside sample_all_communities() to draw observations
draw_samps <- function(n.obs = NULL, list.lengths = NULL, sp = NULL, 
                       sp_prev = NULL, counter = NULL, ...) {
  ## Function to make rows of a data frame with samples from a community.
  ## This is to be called in sample_all_communities() repeatedly until the 
  ## community samples data frame has the same number of rows as n.obs. 
  ##
  ## ARGS:  n.obs - integer.  Number of observations to draw
  ##        list.lengths - a list of list lengths
  ##        sp - a list of virtualspecies objects for species in the community
  ##        sp_prev - prevalence (0 to 1) of each species in sp
  ##        counter - a integer incremented in sample_all_communities used here
  ##                  for generating unique checklist IDs.
  list2env(list(...), envir = environment())
  # randomly resample list lengths
  list.lengths <- sample(list.lengths, size = n.obs, 
                        replace = TRUE)
  # determine how many checklists are needed to get n.obs
  ll <- c()
  j <- 1
  while(sum(as.numeric(ll)) < n.obs) {
    ll <- c(ll, list.lengths[j])
    j <- j + 1
  }
  
  ## determine sampling locations for this community
  # select locations
  locs <- data.frame(wg_randomPoints(bias.raster, n = length(ll), 
                                     prob = TRUE, tryf = 1, 
                                     replaceCells = randomPoints.replace))
  # draw samples and bind into a dataframe
  samps <- bind_rows(mapply(FUN = function(ll, location, list.id, 
                                                locs, sp, sp_prev) {
    coords <- locs[location, c(1, 2)]
    # subset to species present at this site
    present_sp <- sapply(sp, FUN = function(x, coords) {
      raster::extract(x$truth_maps[[1]], y = coords)}, coords = coords)
    present_sp <- present_sp[which(present_sp > 0)]
    l <- min(ll, length(present_sp)) # if not enough sp present, l = n sp.
   # if(l < 1) warning("l is less than 0 in sample_all_communitites")
   # if(l < ll) warning("Number of present species is less than list length")
    
    # get prevalences for present species
    present_sp <- sp_prev[which(names(sp_prev) %in% names(present_sp))]
    
    # Sample 
    # 
    # Probabilities of sampling each species are defined by the
    # twentieth root of prevalence.  This is to make the probability of 
    # recording rare species higher than using prevalence alone would do 
    # (in comparison to common species).
    obs <- tryCatch(sample(present_sp, size = l, replace = F, 
             prob = twentiethrt(present_sp)), 
             error = function(x) NA) 
    # put this sample into a data frame that will be joined by bind_rows
    samp <- data.frame(matrix(nrow = l, ncol = 0))
    samp$sp_name <- names(obs)
    samp$x <- rep(coords$x, l)
    samp$y <- rep(coords$y, l)
    samp$checklist_ID <- rep(paste0("l", counter, ".", list.id), l)
    data.frame(samp)
  }, ll = ll, location = 1:length(ll), list.id = 1:length(ll), 
  MoreArgs = list(locs = locs, sp = sp, sp_prev = sp_prev), SIMPLIFY = F))
  
  samps # return data frame of samples
}

## function to sample communities --------------------------------------------
sample_all_communities <- function(sp.list = NULL, n.obs = NULL, 
                                   n.obs.reference = NULL, 
                                   list.lengths = NULL, ...) {
  # Function to sample a dataset from the species distributions for each 
  # community.  This function produces a dataset that looks like a biological
  # records dataset.
  #
  # ARGS: sp.list - a list of virtualspecies objects that each include an 
  #         element named "community", which is a character string giving
  #         a community identifier (e.g. "community_1")
  #       n.obs - the number of observations to draw
  #       n.obs.reference - "community" or "target species", indicating whether
  #         n.obs refers to the number of records to be drawn for the entire
  #         taxonomic group / community, or the number of records required for
  #         a single target species
  #       list.lengths - a numeric vector of list lengths, given as the 
  #         proportion of all possible species that appear on each list
  #       target.species - (optional) if n.obs.reference = "target species", a 
  #         character vector giving the name of the species that n.obs refers to
  
  # list.lengths gives the proportion of total possible species that are on 
  # each checklist.
  list2env(list(...), envir = environment())
  # set list length for all checklists
  list.lengths <- ceiling(list.lengths * community.size)
  
  # list unique communities
  comms <- as.list(unique(sapply(sp.list, FUN = function(x) {x$community})))
  names(comms) <- comms
  
  if(n.obs.reference == "community") {
    for(i in 1:length(comms)) {
      # select only species in the i^th community
      sp <- sp.list[sapply(sp.list, FUN = function(x) {x$community == comms[i]})]
      # set relative probability of observing each species based on prevalence
      sp_prev <- sapply(sp, FUN = function(x) {x$theoretical_prevalence})
      
      # initiate df for storing samples
      comm_samps <- data.frame(matrix(nrow = 0, ncol = 4))
      colnames(comm_samps) <- c("sp_name", "x", "y", "checklist_ID")
      
      counter <- 0
      while(nrow(comm_samps) < n.obs) {
        n_to_draw <- n.obs - nrow(comm_samps)
        new_samps <- draw_samps(list.lengths = list.lengths,
                                n.obs = n_to_draw, 
                                sp = sp, sp_prev = sp_prev, 
                                counter = counter,
                                bias.raster = bias.raster, 
                                randomPoints.replace = 
                                  randomPoints.replace) 
        comm_samps <- bind_rows(comm_samps, new_samps)
        counter <- counter + 1
      }
      
      # set comm_samps to be exactly the size specified by n.obs
      if(nrow(comm_samps) > n.obs) comm_samps <- comm_samps[1:n.obs, ] 
      comms[[i]] <- comm_samps
    }
  }
  if(n.obs.reference == "target species") {
    stop("n.obs.reference = 'target species' is not yet implemented in sample_all_communities().")
  }
  comms # return list with a biological records df for each community
}
## end function to sample communities ----------------------------------------


## define function to make rasters ------------------------------------------
# turn tab_nrec_hec into spatialGridDataFrame and then raster
make_spatial <- function(x, temp_rast = NULL, field = NULL) {
  coord_cols <- which(colnames(x) %in% c("x", "y"))
  if(any(is.na(x$layer))) {
    # remove blocks with no coordinate values (probably ocean-only blocks)
    x <- x[-which(is.na(x$layer)), ] 
  }
  coordinates(x) <- as.matrix(x[, coord_cols])
  proj4string(x) <- CRS("+init=epsg:29903")
  rast <- rasterize(x, temp_rast, field = field)
  # result <- list(spat_df = x, rast = rast)
  # result
  rast
}
## end function to make rasters -----------------------------------------------





#### values for testing -----------------------------------------------
# 
# seed = 10061983
# community.size = 50
# n.obs = 700
# n.obs.reference = "community"
# shape = "hump"
# polynomial.resp = c(1, 2) # there is a bug if I only want 1 hump, but I think that's what I'll want
# polynomial.coef.min = 0.1
# polynomial.coef.max = 1.3 # 0.1 and 1.3 work
# sl.intercept.min = NULL
# sl.intercept.max = NULL # 0.1 and 0.7 work
# sl.coef.min = NULL
# sl.coef.max = NULL
# pca = TRUE
# bias.raster = bias_rasters[[3]]
# list.lengths = list_length
# env.predictors = raw_pred_df
# error.prob = 0
# randomPoints.replace = TRUE
# vary_along = c("n.obs", "bias.raster")
# 
# nsim = 60
# list.length = list_length
# 
# sp.list = sp_list
# 
# #sp.names <- paste0(rep("sp", n_sp), 1:n_sp)
# 
# ## graphs for testing
# plot(s)
