############################################
## This defines the species distribution modeling methods to be tested.
## Simulation 10.10 (the first version using the simulator package)
## 
## author: Willson Gaul wgaul@hotmail.com
## created: 24 Sep 2018
## last modified: 2 July 2019
############################################

## @knitr methods

make_random_folds <- function(n, nfolds) {
  ## Make folds for random cross validation
  # Divides the indices 1:n into nfolds random folds of about the same size.
  # ARGS:   n - sample size
  #         nfolds - number of folds
  nn <- round(n / nfolds) # calculate number of cells needed per fold
  sizes <- rep(nn, nfolds) # set size of each fold
  # add additional cells to the last fold to get the correct total number of
  # cells accross the sum of all folds (this is needed because founding above 
  # if the total number of cells is not evenly divisible by nfolds).
  sizes[nfolds] <- sizes[nfolds] + n - nn * nfolds 
  b <- c(0, cumsum(sizes)) # vector of fold boundaries to use in for loop below
  ii <- sample(n) # get vector of indices in random order
  folds <- list()
  # assign data indices to folds
  for (i in seq(nfolds))
    folds[[i]] <- ii[seq(b[i] + 1, b[i + 1])]
  folds # return the list of data indices for each fold
}

make_block_folds <- function(sp_data, nfolds, block_size, xOffset = 0, 
                             yOffset = 0) {
  ## Make folds for block cross validation
  # ARGS:   sp_data - Spatial Points dataframe (or something that can be 
  #           coerced to SpatialPointsDataFrame) with species observations
  #         nfolds - number of folds
  #         block_size - integer size of blocks in meters.  This should be based
  #           on autocorrelation range of data (either predictor data or model
  #           residuals of a fitted model.  See vignette for blockCV package).
  #         xOffset; yOffset - integer of offset for blocks in meters
  
  # make spatial blocking by specified range with random assignment
  sb <- spatialBlock(speciesData = sp_data, 
                     species = "Species", 
                     theRange = block_size, # size of blocks in meters
                     k = nfolds, 
                     selection = "random", 
                     iteration = 100, # find evenly dispersed folds
                     biomod2Format = TRUE, 
                     xOffset = xOffset, yOffset = yOffset, 
                     showBlocks = FALSE, progress = FALSE)
  sb
}

make_wg_block_folds <- function(sp_data, all_cells, nfolds, block_size) {
  ## Make folds for block cross validation without using any spatial packages.
  # This is so that it can run on sonic, which has an issue with gdal that 
  # prevents using the blockCV package (as of January 2019).
  # ARGS:   sp_data - a data frame with species observations and 
  #                   eastings/northings sampled cells
  #         all_cells - data frame with e/n values for all cells in study extent
  #         nfolds - integer number of folds to produce
  #         block_size - integer size of blocks in meters
  # VALUE:  all_cells - data frame with easting/northing of all cells and 
  #             columns indicating to which block and fold the cell belongs

  # get original resolution of observation grid cells (to use for offsetting 
  # block boundaries)
  res <- abs(unique(all_cells$x)[order(unique(all_cells$x))][1] - 
    unique(all_cells$x)[order(unique(all_cells$x))][2])
  # get unique northings and eastings values
  norths <- unique(all_cells$y)
  easts <- unique(all_cells$x)
  
  # Ensure that there are at least two folds with data
  # This is a minimal level of optimization of block placement.  Ideally I will
  # do this better somehow to get an even distribution of samples accross folds,
  # but I need to be careful that it doesn't just always put fold boundaries 
  # right through high data-density areas, because this would mean that in most
  # cases there would still be large amounts of test and training data from 
  # directly adjacent areas.  
  folds_w_data <- 0
  iter <- 0
  while(folds_w_data < 2 && iter < 1000) {
  # set the ggrid cell that will be the top left corner of the first block
  top_l_start <- sample_n(all_cells, size = 1)
  # initiate vectors to hold e/n values of grid cells in left, right, top, and 
  # bottom rows of each block

  # make list of which norths values are above the cells in top rows of blocks
  top_rows <- norths[which((norths - top_l_start$y) %% 
                             block_size == 0)]
  # make maximum and minimum row bounds
  top_rows <- c(top_rows, max(top_rows) + block_size, 
                min(top_rows) - block_size)
  # make list of which easts values are above the cells in left column of blocks
  left_cols <- easts[which((easts - top_l_start$x) %% 
                             (block_size) == 0)]
  # make maximum and minimum column bounds
  left_cols <- c(left_cols, max(left_cols) + block_size, 
                 min(left_cols) - block_size)
  
  top_rows <- top_rows - 0.5*res # offset block bounds to above cells
  top_rows <- top_rows[order(top_rows, decreasing = T)]
  left_cols <- left_cols - 0.5*res # offset block bounts to left of cells
  left_cols <- left_cols[order(left_cols, decreasing = F)]
  
  # check map (diagnostic)
  # plot(all_cells$y ~ all_cells$x)
  # for(i in 1:length(left_cols)) {abline(v = left_cols[i])}
  # for(i in 1:length(top_rows)) {abline(h = top_rows[i])}
  
  assign_block <- function(x, top_rows, left_cols) {
    for(j in 1:(length(top_rows)-1)) {
      if(abs(as.numeric(as.character(x[["y"]])) - top_rows[j]) <
             block_size && 
         abs(as.numeric(as.character(x[["y"]])) - top_rows[j + 1]) < 
         block_size) {
        for(k in 1:(length(left_cols)-1)) {
          if(abs(as.numeric(as.character(x[["x"]])) - left_cols[k]) < 
              block_size && 
             abs(as.numeric(as.character(x[["x"]])) - left_cols[k + 1]) < 
              block_size) {
            block <- paste0("bl_", j, ".", k)
          }
        }
      }
    }
    block
  }
  
  blocks <- apply(all_cells, MARGIN = 1, FUN = assign_block, 
                  top_rows = top_rows, left_cols = left_cols)
  all_cells$block <- blocks
  
  # visually check folds (debugging)
  # ggplot(data = all_cells, aes(x = x, y = y, color = block)) + 
  #   geom_point()

  # assign blocks to folds
  blocks <- sample(unique(all_cells$block), replace = F)
  folds <- rep(1:5, length.out = length(blocks))
  folds <- data.frame(block = blocks, fold = folds)
  
  block_df <- left_join(all_cells, folds, by = "block")

  # count how many folds contain samples
  samps_folds <- left_join(sp_data, 
                           block_df[, which(colnames(block_df) %in% 
                                              c("en", "block", "fold"))], 
                           by = c("en")) %>%
    group_by(fold) %>%
    summarise(n = n())
  
  folds_w_data <- length(which(samps_folds$n > 0)) 
  iter <- iter + 1
}
  
  # visually check folds (debugging)
  # ggplot(data = block_df, aes(x = x, y = y, color = as.factor(fold))) + 
  #   geom_point()
  
  block_df
} # end make_wg_block_folds()


### random cross-validation extension ---------------------------------------
cv <- new_method_extension(
  "cv", "random cross validation",
  method_extension = function(model, draw, out,
                              base_method) {
    warning("TODO: re-write random cv method to include folds_en in the output, with train_en and test_en (see block_cv method).")
    nfolds <- 5
    if(is.data.frame(draw$observations)) {
      draw$observations$en <- paste0(round(draw$observations$x), "_", 
                                     round(draw$observations$y))
    }
    #  err <- matrix(NA, ncol(out$beta), nfolds)
    # define possible grid cells
    all_cells <- data.frame(rasterToPoints(draw$prob.raster))
    all_cells$en <- paste0(all_cells$x, "_", all_cells$y)
    ii <- tryCatch({make_random_folds(nrow(all_cells), 
                               nfolds)}, 
                   error = function(x) NA)
    fits <- list() # initialize list for results
    preds <- list() # initialize list for results
    folds_en <- list() # initialize list for cv fold en values
    
    for (i in seq_along(ii)) {
      train <- model
      # define training data cells
      train_cells <- all_cells[-ii[[i]], ]
      train_draw <- draw
      train_draw$observations <- tryCatch(
        {train_draw$observations[which(train_draw$observations$en %in% 
                                         train_cells$en), ]}, 
        error = function(x) NA)
      
      test <- model
      test_cells <- all_cells[ii[[i]], ]
      test_draw <- draw
      test_draw$observations <- tryCatch({
        test_draw$observations[which(test_draw$observations$en %in% 
                                       test_cells$en), ]},
        error = function(x) NA)
      
      m <- tryCatch({base_method@method(
        model = train, draw = train_draw)}, 
        error = function(c) paste0(
          c, ". CV model fitting failed with species ", 
          draw$sp_name))
      fits[i] <- tryCatch({m$fit}, 
                            error = function(x) NA)
      preds[i] <- tryCatch({m$preds}, 
                             error = function(x) NA)
      folds_en[[i]] <- tryCatch({list(train_en = train_cells$en, 
                                      test_en = test_cells$en)}, 
                                error = function(x) NA)
      
      # yhat <- test$x %*% fit$beta
      # ll <- seq(ncol(yhat))
      # err[ll, i] <- colMeans((yhat - test_draw)^2)
    }
    # m <- rowMeans(err)
    # se <- apply(err, 1, sd) / sqrt(nfolds)
    # imin <- which.min(m)
    # ioneserule <- max(which(m <= m[imin] + se[imin]))
    # list(err = err, m = m, se = se, imin = imin,
    #      ioneserule = ioneserule,
    #      beta = out$beta[, imin],
    #      yhat = model$x %*% out$beta[, imin])
    list(fit = fits, preds = preds, folds_en = folds_en) 
  }
)
### end random cross-validation extension ------------------------------------

### block cv extension -----------------------------------------------------
## Using blockCV package (requires gdal-config)
block_cv <- new_method_extension(
  "block_cv", "block cross validation",
  method_extension = function(model, draw, out,
                              base_method) {
    nfolds <- 5
    block_size <- 100000
    yOffset <- 0
    xOffset <- 0
    train_model <- model
    
    if(is.data.frame(draw$observations)) {
      draw$observations$en <- paste0(round(draw$observations$x), "_", 
                                     round(draw$observations$y))
    }
    # convert obsevations to spatial points df
    sp_data <- tryCatch({
      SpatialPointsDataFrame(
        coords = draw$observations[, c("x", "y")],
        data = draw$observations[ , c("checklist_ID", "list_length", "sp1")],
        proj4string = CRS("+init=epsg:29903"))
    }, 
    error = function(x) NA)
    if(class(sp_data) == "SpatialPointsDataFrame") {
      sp_data$en <- paste0(round(sp_data$x, 0), "_", round(sp_data$y))
    }
    
    # define possible grid cells
    truth_cells <- data.frame(rasterToPoints(draw$truth_maps[[1]]))
    truth_cells$en <- paste0(truth_cells$x, "_", truth_cells$y)
    truth_cells <- SpatialPointsDataFrame(
      coords = truth_cells[, c("x", "y")], 
      data = truth_cells[, c("en", "layer")], 
      proj4string = crs(draw$prob.raster)
    )
    
    sb <- tryCatch({  
      # make spatial blocking by specified range with random assignment
      spatialBlock(speciesData = sp_data, 
                   species = "sp1", 
                   theRange = block_size, # size of blocks in meters
                   k = nfolds, 
                   selection = "random", 
                   rasterLayer = draw$truth_maps[[1]], 
                   maskBySpecies = FALSE, 
                   iteration = 100, # find evenly dispersed folds
                   biomod2Format = TRUE, 
                   xOffset = xOffset, yOffset = yOffset, 
                   showBlocks = FALSE, progress = FALSE)}, 
      error = function(x) NA)
    folds <- tryCatch({sb$folds}, error = function(x) NA)
    
    # optionally plot blocks with checklist points:
    # sb$plots + 
    #   geom_point(data = as.data.frame(coordinates(sp_data)), 
    #              aes(x = x, y = y))
    # foldExplorer(sb, rasterLayer = draw$truth_maps[[1]], speciesData = sp_data)
    
    fits <- list() # initialize list for results
    preds <- list() # initialize list for results
    folds_en <- list() # initialize list for cv fold en values

    for(k in 1:length(folds)) {
      # training and test data must stay within the format of draws so that
      # they can be used by the base model fitting methods.  So, I need
      # to subset the draw$observations to get training data, and subset
      # draw$truth_maps[[1]] to get true sp p/a values for testing
      
      #TODO 19 Dec.  (On 4 Jan I am wondering if I have already done this and
      # forgot to remove this note?  I guess I need to re-check at some point...
      # This sets train_en using en values from sp_data indices, but
      # there are no-data cells that SHOULD be in training data but aren't in 
      # sp_data.  This really matters because test data is set by en values not
      # in train_en.  So I need to find a way to get ALL cells that would be
      # within the test & training blocks.  But feeding sp_data into the 
      # blockCV function is nice because it tries to get even-ish data in all
      # blocks.  Need to figure out if I can get the box limits of blocks
      # from the blockCV function.
      # USE 'blocks' from spatialBlocks output to mask all_cells to train_en
      train_ind <- unlist(folds[[k]][1]) # extract the training set indices
      test_ind <- unlist(folds[[k]][2]) # extract test set indices
      # get en values for grid cells in training and test sets
      test_en <- tryCatch({
        truth_cells[sb$blocks[which(sb$blocks$folds == k), ], ]$en
        }, 
        error = function(x) NA)
      if(!is.na(test_en)) {
        train_en <- unique(truth_cells$en[which(truth_cells$en %nin% test_en)])
      } else train_en <- NA
      # The above selects all grid cells that do not
      # contain any training data.  TODO: 18 Dec check that this works in 
      # extreme bias cases - does spatialBlock always make sure there is at 
      # least one observation in every block?  If so, then this way of selecting
      # test data should work.  If there are sometimes blocks with no 
      # observations, then this method will results in some test data containing
      # blocks that are technically within the "training" set (even there is no
      # data from those blocks being used to train model because there are no
      # observations from those blocks).  Not sure if that's a problem. Maybe
      # not since no data from there is being used to train the model.
      
      train_draw <- draw # copy the draw object
      # subset observations to only those in the traning set
      train_draw$observations <- tryCatch( 
        {train_draw$observations[which(train_draw$observations$en %in% 
                                         train_en), ]}, 
        error = function(x) NA)
      
      # optionally plot test and training grid cells, and observation points
      # plot(ir_TM75)
      # points(truth_cells[which(truth_cells$en %in% train_en), ])
      # points(truth_cells[which(truth_cells$en %in% test_en), ], col = "red")
      # points(x = train_draw$observations$x, y = train_draw$observations$y, 
      #        col = "green")
      # points(sp_data[sb$folds[[k]][[1]], ], col = "blue", pch = 3)
      
      # fit base model with only training data
      m <- tryCatch({base_method@method(
        model = train_model, draw = train_draw)}, 
        error = function(c) paste0(
          c, ". Block CV model fitting failed with species ", 
          draw$sp_name))
      
      # store results for this fold in lists that will be returned
      fits[k] <- tryCatch({m$fit}, 
                          error = function(x) NA)
      preds[k] <- tryCatch({m$preds}, 
                           error = function(x) NA)
      folds_en[[k]] <- tryCatch({list(train_en = train_en, 
                                      test_en = test_en)}, 
                                error = function(x) NA)
    }
    
    list(fit = fits, preds = preds, folds_en = folds_en) 
  }
)
### end block cv extension ----------------------------------------------------

### wg block cv extension -----------------------------------------------------
# Do block cross-validation, but make blocks without using spatial packages 
#(so it can run on Sonic)
wg_block_cv <- new_method_extension(
  "wg_block_cv", "wg block cross validation",
  method_extension = function(model, draw, out,
                              base_method) {
    # set parameters
    nfolds <- 5
    block_size <- 100000 # this should be based on autocorrelation range of resids
    train_model <- model
    
    if(is.data.frame(draw$observations)) {
      draw$observations$en <- paste0(round(draw$observations$x), "_", 
                                     round(draw$observations$y))
    }
    # define possible grid cells
    all_cells <- data.frame(rasterToPoints(draw$truth_maps[[1]]))
    all_cells$en <- paste0(all_cells$x, "_", all_cells$y)
    
    # make spatial blocking and assign blocks to folds
    # This could potentially be repeated to find optimal blocking in terms of
    # even number of observations in each fold
    sb <- tryCatch({  
      # make spatial blocking by specified range with random assignment
      make_wg_block_folds(sp_data = draw$observations, all_cells = all_cells,
                          nfolds = nfolds, 
                          block_size = block_size)
    }, 
    error = function(x) NA)
    folds <- tryCatch({unique(sb$fold)}, error = function(x) NA)

    # optionally plot blocks with checklist points:
    # ggplot() + 
    #   geom_point(data = sb, aes(x = x, y = y, color = factor(fold))) + 
    #   geom_point(data = draw$observations, aes(x = x, y = y), shape = "x", size = 3)
    
    fits <- list() # initialize list for results
    preds <- list() # initialize list for results
    folds_en <- list() # initialize list for cv fold en values
    
    for(k in 1:length(folds)) {
      # training and test data must stay within the format of draws so that
      # they can be used by the base model fitting methods.  So, I need
      # to subset the draw$observations to get training data, and subset
      # draw$truth_maps[[1]] to get true sp p/a values for testing
      
      # get en values for grid cells in training and test sets
      test_en <- tryCatch({
        sb$en[which(sb$fold == k)]
      }, 
      error = function(x) NA)
      
      if(!is.na(test_en)) {
        train_en <- sb$en[which(sb$fold != k)]
      } else train_en <- NA
      
      train_draw <- draw # copy the draw object
      # subset observations to only those in the traning set
      train_draw$observations <- tryCatch( 
        {train_draw$observations[which(train_draw$observations$en %in% 
                                         train_en), ]}, 
        error = function(x) NA)
      
      # optionally plot test and training grid cells, and observation points
      # plot(ir_TM75)
      # points(sb[which(sb$en %in% train_en), ])
      # points(sb[which(sb$en %in% test_en), ], col = "red")
      # points(x = train_draw$observations$x, y = train_draw$observations$y, 
      #        col = "green")
      # points(draw$observations[which(draw$observations$en %in% test_en), ], 
      #        col = "blue", pch = 3)
      
      # fit base model with only training data
      m <- tryCatch({base_method@method(
        model = train_model, draw = train_draw)}, 
        error = function(c) paste0(
          c, ". Block CV model fitting failed with species ", 
          draw$sp_name))
      
      # store results for this fold in lists that will be returned
      fits[k] <- tryCatch({m$fit}, 
                          error = function(x) NA)
      preds[k] <- tryCatch({m$preds}, 
                           error = function(x) NA)
      folds_en[[k]] <- tryCatch({list(train_en = train_en, 
                                      test_en = test_en)}, 
                                error = function(x) NA)
      }
    
    list(fit = fits, preds = preds, folds_en = folds_en) 
  }
)

### end wg block cv extension

##############   
#### end method extensions -------------------------------------------------
##############

## GLM with polynomial terms -------------------------------------------------
glm_poly <- new_method(
  "glm_poly", "Logistic Regression with polynomial terms", 
  method = function(model, draw) {
    ## outputs: - a list with 
    ##              * a glm object for each species trained using only 
    ##                training data
    ##              * a (raster or data frame) with the predicted probability
    ##                of species presence in each cell for each species
    ##
    ## TODO:  
    ##        - or automatic variable selection with AUC VIM (Neuman et al 2016)

    sp_name <- draw$sp_name
    chkl_df <- draw$observations 
    predictors <- model@params$env.predictors
    predictors$en <- paste(round(predictors$s1), round(predictors$s2), 
                           sep = "_")
    sp_df <- tryCatch({
      gather(chkl_df, key = "sp_name", value = "present", eval(sp_name))
    }, error = function(x) NA)
    
    sp_df$en <- tryCatch({
      paste(round(sp_df$x), round(sp_df$y), sep = "_")
    }, error = function(x) NA)
    
    sp_df <- tryCatch({
      # join predictor variables to species observation df
      sp_df <- left_join(sp_df, predictors, by = "en")
    }, error = function(x) NA)
    
    ## fit glm
    # make formula using a random subset of predictors each time  
    var_names <- names(predictors)[which(names(predictors) %nin% c("s1", "s2", 
                                                                   "en"))]
    var_names <- sample(var_names, size = 0.5*length(var_names), replace = F)
    form <- as.formula(paste0("present ~ ", paste(var_names, collapse = " + "), 
                              " + ", 
                              paste("poly(", var_names, ", degree = 2)",
                                    collapse = " + ")))
    
    # fit full glm model
    m <- tryCatch({
      do.call("glm", list(form, 
              data = as.name("sp_df"), family = "binomial", model = FALSE,
              x = FALSE))
    }, 
    error = function(c) {
      c$message <- paste0(c$message, 
                          " glm do.call failing while fitting glm to species ", 
                          sp_name, ".") 
    })
    
    # Optimize glm model
    m <- tryCatch({
      do.call("step", list(m, direction = "both", trace = 0))
    }, 
    error = function(c) {
      c$message <- paste0(c$message, 
                          " GLM step selection failing with species ", 
                          sp_name, ".")
    })
    
    # predict probabilities of presence in each grid cell from fitted model
    predictions_df <- predictors
    predictions_df$list_length <- rep(model@params$community.size, 
                                      nrow(predictions_df))
    
    predictions_df$predicted_prob_occ <- tryCatch({
      do.call("predict", list(m, newdata = as.name("predictions_df"), 
                              type = "response"))
    }, 
    error = function(c) {
      c$message <- paste0(c$message, " predict do.call failing while predicting occurrences for species ", 
                          sp_name, ".")
    })
    
    predictions_df <- predictions_df[, which(names(predictions_df) %in% 
                                               c("s1", "s2", "predicted_prob_occ"))]
    predictions_df$en <- paste0(round(predictions_df$s1), "_", 
                                round(predictions_df$s2))
    
    # return the model object and the data frame with predictions for each cell
    # using this fitted model.
    list(fit = list(m), preds = list(predictions_df))
  }
)
## end polynomial GLM --------------------------------------------------------


## Random Forest ------------------------------------------------------------
rf <- new_method(
  "random_forest", "Random Forest", 
  method = function(model, draw) {
    ## outputs: - a list with 
    ##              * a randomForest object for each species trained using only 
    ##                training data
    ##              * a (raster or data frame) with the predicted probability
    ##                of species presence in each cell for each species
    sp_name <- draw$sp_name
    chkl_df <- draw$observations 
    predictors <- model@params$env.predictors
    predictors$en <- paste(round(predictors$s1), round(predictors$s2), 
                           sep = "_")
    
    sp_df <- tryCatch({ # gather the column of observations for this sp.
      gather(chkl_df, key = "sp_name", value = "present", eval(sp_name))
    }, error = function(x) NA)
    
    sp_df$en <- tryCatch({ # make east/north character string column
      paste(round(sp_df$x), round(sp_df$y), sep = "_")
    }, error = function(x) NA)
    
    sp_df <- tryCatch({
      # join predictor variables to species observation df
      sp_df <- left_join(sp_df, predictors, by = "en")
    }, error = function(x) NA)
    
    # use random subset of variables for each species
    var_names <- names(predictors)[which(names(predictors) %nin% c("s1", "s2", 
                                                                   "en"))]
    var_names <- sample(var_names, size = 0.5*length(var_names), replace = F)
    
    ## fit random forest model
    # Note that the "predicted", "votes", and error metrics in the object 
    # produced by the randomForest call are predictions and error metrics for 
    # the training data, which in this case is CHECKLISTS, not grid cells.  So
    # for biased sampling, there is not test in the model of how the model 
    # predicts agains non-sampled grid cells.  Therefore these metrics are not
    # of interest to me in evaluating how well the random forest actually works.
    # For that I need to predict to the grid cells (as done below).  The 
    # prediction and preformance results returned in the randomForest object are
    # however of interest because this is what would be used to tune the model
    # in a real-world setting (e.g. variable selection, etc).  

    # TODO: optimize mtry here using rfcv to avoid over fitting?  
    # Probably only needed if I have many predictor variables, which I 
    # currently do not. 
    
    # optimize number of variables to use at each node
    mtry_err <- c() # make list to hold errors for the mtry tests
    nvar <- length(var_names) # number of variables
    # mtry values to test - default, half that, and twice that
    mtry_tests <- c(ceiling(sqrt(nvar)/2), ceiling(sqrt(nvar)), 
                    ceiling(sqrt(nvar)*2)) 
    
    for(k in 1:length(mtry_tests) ) {
      m_k <- tryCatch({
        randomForest(
          x = sp_df[, colnames(sp_df) %in% var_names],
          y = factor(sp_df[, colnames(sp_df) == "present"]),
          ntree = 1000, 
          mtry = mtry_tests[k], 
          nodesize = 1, 
          replace = TRUE, classwt = NULL, 
          importance = FALSE, 
          keep.forest = FALSE)},
        error = function(c) {
          c$message <- paste0(c$message, 
                              " RF failing with species ", sp_name, ".")
        })
      
      # if error from this model is lowest so far, keep this model
      if(class(m_k) == "randomForest" && 
         !is.na(m_k$err.rate[nrow(m_k$err.rate), "OOB"])) {
        mtry_err[k] <- m_k$err.rate[nrow(m_k$err.rate), "OOB"] # error for this mtry
      }
    }
    
    # get best mtry value.  If multiply mtry values tied for best error rate,
    # use the one that is closest to the square root of the number of variables
    if(length(which(mtry_err == min(mtry_err))) > 1) {
      # calculate the distance of each tied "best" mtry value from the default
      dist_from_default <- tryCatch({
        abs(mtry_tests[mtry_err == min(mtry_err)] - sqrt(nvar))}, 
        error = function(x) NA)
      # keep the mtry value that is closest to the default
      mtry_best <- tryCatch({
        mtry_tests[mtry_err == min(mtry_err)][dist_from_default == 
                                                min(dist_from_default)]}, 
        error = function(x) NA)
    } else mtry_best <- tryCatch(mtry_tests[mtry_err == min(mtry_err)], 
                                 error = function(x) NA)
    
    # fit model with optimum mtry 
    # use 1000 trees which is hopefully high enough to get stable variable
    # importance if I want it (see manual linked in help documentation)
    m <- tryCatch({
      randomForest(
        x = sp_df[, which(colnames(sp_df) %in% var_names)],
        y = factor(sp_df[, which(colnames(sp_df) == "present")]),
        ntree = 1000, 
        mtry = mtry_best, 
        nodesize = 1, 
        replace = TRUE, classwt = NULL, 
        importance = FALSE, 
        keep.forest = TRUE)}, error = function(x) NA)
    
    m$data <- sp_df # store data for auc calculation in evals()
    ## predict probabilities of presence in each cell from rf model
    predictions_df <- predictors
    predictions_df$predicted_prob_occ <- tryCatch({
      predict(m, newdata = predictions_df, type = "prob")[, "1"]}, 
    error = function(c) {
      c$message <- paste0(c$message, 
                          " predict do.call failing while predicting occurrences for species ", 
                          sp_name, ".")})
    
    predictions_df <- predictions_df[, which(names(predictions_df) %in% 
                                               c("s1", "s2", 
                                                 "predicted_prob_occ"))]
    predictions_df$en <- paste0(round(predictions_df$s1), "_", 
                                round(predictions_df$s2))
    
    ## return the model object and the data frame with predictions for each cell
    # using this fitted model.
    list(fit = list(m), preds = list(predictions_df))
  }
)
## end random forest --------------------------------------------------------

## Boosted regression tree -------------------------------------------------
brt <- new_method(
  "brt", "Boosted Regression Tree", 
  method = function(model, draw) {
    ## outputs: - a list with 
    ##              * a gbm object for each species trained using only 
    ##                training data
    ##              * a (raster or data frame) with the predicted probability
    ##                of species presence in each cell for each species
    ##            
    ## notes: - response and predictor variables are indicated by column 
    ##          position, not name.      
    ## starting with:
    #   - tree complexity of 5
    #   - learning rate 0.005
    #   - maximum n.tree of 5000
    #   - stepping function to determine number of trees given lr
    sp_name <- draw$sp_name
    chkl_df <- draw$observations 
    predictors <- model@params$env.predictors
    predictors$en <- paste(round(predictors$s1), round(predictors$s2), 
                           sep = "_")
    
    # predictive performance is generally better with smaller learning rate, 
    # but computation and memory are increased (because of more iterations). 
    # See vignette('gbm')
    # I start with learning rate small (which will result in large n.trees, 
    # long computation times, and large memory).  I set the maximum n.trees I 
    # want to allow based on computation time.  I then try to fit a model with
    # the small learning rate.  If the model fits with fewer than the 
    # max n.trees, I keep it.  If the model exceeds max n.trees, I discard it,
    # increase the learning rate, and try again.  I try this for learning rates
    # between 0.001 (best) and 0.01 (worst).  If no model fits with a learning
    # rate of <= 0.01, I give up fitting a model for this case.
    
    lr <- 0.001 # try this learning rate first 
    max_ntree <- 10000 # maximum number of trees to allow
    tc <- 5 # tree complexity

    chkl_df <- tryCatch({ # gather the column of observations for this sp.
      gather(chkl_df, key = "sp_name", value = "present", eval(sp_name))
    }, error = function(x) NA)
    
    chkl_df$en <- tryCatch({ # make east/north character string column
      paste(round(chkl_df$x), round(chkl_df$y), sep = "_")
    }, error = function(x) NA)
    
    chkl_df <- tryCatch({
      # join predictor variables to species observation df
      chkl_df <- left_join(chkl_df, predictors, by = "en")
    }, error = function(x) NA)
    
    # use random subset of variables for each species
    var_names <- names(predictors)[which(names(predictors) %nin% c("s1", "s2", 
                                                                   "en"))]
    var_names <- sample(var_names, size = 0.5*length(var_names), replace = F)

    fit_model <- TRUE # indicator for while loop
    # while an adequate model has not been fit, try to fit a model
    while (fit_model == TRUE) { 
      m1 <- tryCatch({gbm.step(chkl_df, 
                               gbm.x = var_names,
                               gbm.y = "present", 
                               tree.complexity = 5, 
                               learning.rate = lr, 
                               bag.fraction = 0.5, 
                               max.trees = max_ntree,
                               family = "bernoulli", 
                               plot.main = FALSE)
      }, error = function(x) NA)
      # If model fit successfully with fewer than max.trees, keep this model and
      # do not try any other learning rates, as lower learning rates are 
      # generally prefered as I understand it (see vignette('gbm')). 
      if(!is.null(m1) && !is.na(m1) && m1$n.trees < max_ntree) {
        m <- m1
        fit_model <- FALSE
      }

      # if model fitting stopped because max.trees was reached, then discard
      # this model and fit a new one, unless the learning rate is higher than 
      # 0.01, in which case give up. 
      if(fit_model == TRUE) {
        lr <- lr + 0.001
        if(lr > 0.02) {
          fit_model <- FALSE
          m <- "Not brt model fit with fewern than the maximum number of trees and a learning rate less than 0.01 (see vignette('gbm'))."
        }
      }
    } # end while loop
    
    ## predict probabilities of presence in each cell from rf model
    predictions_df <- predictors
    predictions_df$predicted_prob_occ <- tryCatch({
      predict(m, newdata = predictions_df, n.trees = m$n.trees, 
              type = "response")
    }, 
    error = function(c) {
      c$message <- paste0(c$message, ". Prediction do.call failing while predicting occurrences for species ", 
                          sp_name, ".")
    })

    if(class(m) == "gbm") m$Terms <- NA # remove terms to save memory
    ## return the model object and the data frame with predictions for each cell
    # using this fitted model.
    list(fit = list(m), preds = list(predictions_df))
  }
)
## end boosted regression tree -------------------------------------------------

## biomod2 ensemble -----------------------------------------------------------
biomod_ens <- new_method(
  "biomod_ens", "biomod2 Ensemble", 
  method = function(model, draw) {
    ## outputs: - a list with 
    ##              * a biomod2 object for each species trained using only 
    ##                training data
    ##              * a (raster or data frame) with the predicted probability
    ##                of species presence in each cell for each species
    sp_name <- draw$sp_name
    chkl_df <- draw$observations 
    predictors <- model@params$env.predictors
    predictors$en <- paste(round(predictors$s1), round(predictors$s2), 
                           sep = "_")
    chkl_df <- tryCatch({ # gather the column of observations for this sp.
      gather(chkl_df, key = "sp_name", value = "present", eval(sp_name))
    }, error = function(x) NA)
    
    chkl_df$en <- tryCatch({ # make east/north character string column
      paste(round(chkl_df$x), round(chkl_df$y), sep = "_")
    }, error = function(x) NA)
    
    chkl_df <- tryCatch({
      # join predictor variables to species observation df
      chkl_df <- left_join(chkl_df, predictors, by = "en")
    }, error = function(x) NA)
    
    # use random subset of variables for each species
    var_names <- names(predictors)[which(names(predictors) %nin% c("s1", "s2", 
                                                                   "en"))]
    var_names <- sample(var_names, size = 0.5*length(var_names), replace = F)
    
    # make a directory to hold biomod data and temporarily set it to be working
    # directory.  This is because biomod2 writes objects to the working 
    # directory and I can't find a way to specify a differnt location to write
    # them to.  But I want to keep my main working directory for the simulation
    # clean (not have a new directory created for every species).
    if(!dir.exists("./biomod2_outs")) {dir.create("biomod2_outs")}
    orig_dir <- getwd()
    setwd(paste0(getwd(), "/biomod2_outs"))
    bm_data <- BIOMOD_FormatingData(resp.var = chkl_df$present, 
                                    resp.xy = chkl_df[, c(
                                      which(colnames(chkl_df) == "x"), 
                                      which(colnames(chkl_df) == "y"))], 
                                    resp.name = sp_name, 
                                    expl.var = chkl_df[, which(
                                      colnames(chkl_df) %in% var_names)])
    # specify modelling options
    # TODO 12 Feb decide what options I want to specify for each method
    warning("Biomod2 currently using default options.  WG should specify these.")
    browser()
    bm_opts <- BIOMOD_ModelingOptions()
    
    mod <- BIOMOD_Modeling(bm_data, 
                           models = c("GLM", "RF"), 
                           models.options = bm_opts, 
                           NbRunEval = 1, 
                           DataSplit = 100, 
                           Yweights = NULL, 
                           Prevalence = NULL, 
                           VarImport = 0, 
                           models.eval.meth = c("TSS", "ROC"), 
                           SaveObj = TRUE, 
                           rescal.all.models = FALSE, 
                           do.full.models = TRUE)
    
    # TODO here 12 Feb.  Need to make ensemble next.
    # TODO 13 Feb look at EnsembleModelingAssembly vignette
    ens_mod <- BIOMOD_EnsembleModeling(
      modeling.output = mod, 
      chosen.models = "all", 
      em.by = "all", 
      eval.metric = c("ROC"), 
      eval.metric.quality.threshold = 0.6, 
      models.eval.meth = c("ROC", "TSS"), 
      prob.mean = T,
      prob.cv = T, 
      prob.ci = T, 
      prob.ci.alpha = 0.05, 
      prob.median = T, 
      committee.averaging = T, 
      prob.mean.weight = T, 
      prob.mean.weight.decay = "proportional")
    
    ens_pred <- BIOMOD_Projection(
      modeling.output = ens_mod, 
      new.env = #TODO: 13 Feb put data to predict to here.
    )
    
    setwd(orig_dir) # reset working directory to what it should be
  }
)
## end biomod2 ensemble -------------------------------------------------------

## krigging -------------------------------------------------------------------
krig <- new_method(
  "krig", "Krigging proportion of lists", 
  method = function(model, draw) {
    ## outputs: - a list with 
    ##              * a TODO: ??krig model?? for each species trained using only 
    ##                training data
    ##              * a (raster or data frame) with the predicted probability
    ##                of species presence in each cell for each species
    sp_name <- draw$sp_name
    chkl_df <- draw$observations 
    chkl_df <- tryCatch({ # gather the column of observations for this sp.
      gather(chkl_df, key = "sp_name", value = "present", eval(sp_name))
    }, error = function(x) NA)
    
    chkl_df$en <- tryCatch({ # make east/north character string column
      paste(round(chkl_df$x), round(chkl_df$y), sep = "_")
    }, error = function(x) NA)
    
    # function to calculate prop. of checklists on which the species appears
    calc_prop <- function(x, list_id) {sum(x) / length(unique(list_id))}
    
    chkl_df <- tryCatch({
      # summarise observations to the grid cell level by calculating the 
      # proportion of all checklists in the grid cell that have the species.
      # This should produce a df with a single row for each grid cell that has
      # a checklist. No value will be produced for grid cells with no checklists
      group_by(chkl_df, en) %>%
        summarise(prob_present = calc_prop(x = present, list_id = checklist_ID), 
                  x = unique(x), y = unique(y), sp_name = unique(sp_name))
    }, error = function(x) NA)

    # create empirical variogram for prob_present
    # TODO 3 March - I could use universal krigging to add list length or nrec
    # or n lists to formula on right-hand side
    gs_prop <- gstat(id = "prob_occ", 
                     formula = prob_present ~ 1,
                     data = chkl_df, 
                     locations = ~ x + y)
    vg_prop <- variogram(gs_prop)
    #vgm(as.character(vgm()[, 1])) # show vgm options
    # TODO: 3 March This does krigging, but I think just leaves values in cells 
    # as-is if there is any data in the cell.  Make sure it's at least filling 
    # in some blank cells.  If not, what's the point?
    warning("krigging method is likely doing dumn stuff right now.")
    f_vg_prop <- tryCatch({
      fit.variogram(vg_prop, 
                    vgm(c("Nug", "Exp", "Sph", "Gau", "Exc", "Mat", "Ste")))
    }, error = function(x) NA)
     plot(variogramLine(f_vg_prop, max(vg_prop$dist)), type = 'l')
     points(vg_prop[, 2:3], pch = 20, col = 'red')
    
    # use variogram in krigging interpolation
     krg_prop <- gstat(formula = prob_present ~ 1, 
                       data = chkl_df, 
                       model = f_vg_prop, 
                       locations = ~ x + y)
     # get locations to predict at from env.predictors (which has all grid cells)
     pred_locs <- model@params$env.predictors[, which(
       colnames(model@params$env.predictors) %in% c("s1", "s2"))]
     colnames(pred_locs) <- c("x", "y") # give column names to match "locations"
     pred_locs$en <- paste0(pred_locs$x, "_", pred_locs$y)
     pred_locs <- pred_locs[which(pred_locs$en %nin% chkl_df$en), ]
     
     # predict probability of presence in each cell by ordinary krigging
     predictions_df <- tryCatch({
       predict(krg_prop, newdata = pred_locs)
       colnames(predictions_df)[which(colnames(predictions_df) == "var1.pred")] <-
         "prob_present"
       predictions_df$en <- paste0(predictions_df$x, "_", predictions_df$y)
       predictions_df <- tryCatch({full_join(chkl_df, predictions_df)}, 
                                  error = function(x) NA)
       }, 
       error = function(c) {
         c$message <- paste0(c$message, 
                             ". Krigging prediction failing with species ", 
                             sp_name, ". ")
       })

     # print(ggplot(data = predictions_df, 
     #              aes(x = x, y = y, color = prob_present)) + 
     #   geom_point())
     # print(ggplot(data = predictions_df, aes(x = x, y = y, color = var1.var)) + 
     #   geom_point())
     
     ## return the gstat object and the data frame with predictions for each cell
     # using this fitted model.
     list(fit = list(krg_prop), preds = list(predictions_df))
  }
)
## end krigging ---------------------------------------------------------------

## inverse distance weighted interpolation ------------------------------------
idw_interp <- new_method(
  "idw_interp", "Inverse Distance Interpolation - Proportion of Lists", 
  method = function(model, draw) {
    ## outputs: - a list with 
    ##              * a gstat object with the IDW interpolation parameters
    ##              * a (raster or data frame) with the predicted probability
    ##                of species presence in each cell for each species
    sp_name <- draw$sp_name
    chkl_df <- draw$observations 
    chkl_df <- tryCatch({ # gather the column of observations for this sp.
      gather(chkl_df, key = "sp_name", value = "present", eval(sp_name))
    }, error = function(x) NA)
    chkl_df$en <- tryCatch({ # make east/north character string column
      paste(round(chkl_df$x), round(chkl_df$y), sep = "_")
    }, error = function(x) NA)
    
    # function to calculate prop. of checklists on which the species appears
    calc_prop <- function(x, list_id) {sum(x) / length(unique(list_id))}
    
    chkl_df <- tryCatch({
      # summarise observations to the grid cell level by calculating the 
      # proportion of all checklists in the grid cell that have the species.
      # This should produce a df with a single row for each grid cell that has
      # a checklist. No value will be produced for grid cells with no checklists
      group_by(chkl_df, en) %>%
        summarise(prob_present = calc_prop(x = present, list_id = checklist_ID), 
                  x = unique(x), y = unique(y), sp_name = unique(sp_name))
    }, error = function(x) NA)
    
    # get locations for prediction from env.predictors (which has all grid cells)
    pred_locs <- model@params$env.predictors[, which(
      colnames(model@params$env.predictors) %in% c("s1", "s2"))]
    colnames(pred_locs) <- c("x", "y") # give column names to match "locations"
    pred_locs$en <- paste0(pred_locs$x, "_", pred_locs$y)
    pred_locs <- tryCatch({pred_locs[which(pred_locs$en %nin% chkl_df$en), ]}, 
                          error = function(x) NA)
    
    ## choose best weighting power and number of points for inverse distance 
    ## weigthed interpolation
    ## using 3-fold cross-validation testing using the observed data (not the
    ## truth).  I think this might be how a practitioner would chose the 
    ## distance power and number of points
    # make data frame to hold powers, number of points, and AUC for each combo
    # test powers from 0 to 10 in increments of 0.5
    # test number of points from 1 to all points in increments of 2

    performance <- tryCatch({
      expand.grid(power = round(seq(from = 0, to = 10, by = 0.5), digits = 2), 
                  npoints = seq(from = 1, to = nrow(chkl_df), by = 2), 
                  RMSE = NA)}, 
      error = function(x) {
        data.frame(power = NA, npoints = NA, auc = NA, RMSE = NA)})
    
    # make data folds 
    folds <- tryCatch({
      sample(rep(1:3, ceiling(nrow(chkl_df)/length(1:3))), 
             size = nrow(chkl_df), replace = F)}, 
      error = function(x) NA)
    
    # choose optimal power and number of points using RMSE 
    # store AUC also in case I want to do comparison between RMSE and AUC later
    for(i in 1:nrow(performance)) {
      p <- performance$power[i]
      ndat <- performance$npoints[i]
      cv_aucs <- c() # initialize object to hold AUC values for each fold
      rmspe <- c() # initialize object to hold RMSPE values for each fold
      for(j in 1:length(unique(folds))) {
        temp_prop <-tryCatch({
          gstat(formula = prob_present ~ 1, 
                data = chkl_df[folds != j, ], # train with data not in the jth fold
                nmin = ndat,
                nmax = ndat,
                force = T, # use at least nmin points
                locations = ~ x + y, 
                set = list(idp = p))}, 
          error = function(x) NA)

        # predict probability of presence in each holdout cell by interpolation
        temp_locs <- tryCatch({chkl_df[folds == j, ]}, 
                              error = function(x) NA)
        temp_pred_df <- tryCatch({
          predict(temp_prop, newdata = temp_locs)
        }, 
        error = function(x) NA)
        
        try({
          # change column name to match that returned by other methods
          colnames(temp_pred_df)[which(
            colnames(temp_pred_df) == "var1.pred")] <- "predicted_prob_occ"
          temp_pred_df$en <- paste0(round(temp_pred_df$x), "_", 
                                    round(temp_pred_df$y))
          temp_pred_df$x <- round(temp_pred_df$x)
          temp_pred_df$y <- round(temp_pred_df$y)}
        )
        
        try(if(!identical(temp_pred_df$en, temp_locs$en)) warning("Test data and prediction locations do not match in idwi power optimization loop."))
        
        # test performance with random cross-validation
        rmspe[j] <- tryCatch(
          # square root of mean (predicted minus real)^2
          # for out-of-sample locations
          sqrt(mean((temp_pred_df$predicted_prob_occ -
                       temp_locs$prob_present)^2)),
          error = function(x) NA)
      }
      performance$RMSE[i] <- tryCatch(mean(rmspe, na.rm = T), 
                                      error = function(x) NA)
    }

    # select power that gave best RMS in test with observed data
    # if multiple powers performed equally well, use the largest power
    pwr <- tryCatch({
      max(performance$power[which(performance$RMSE == min(performance$RMSE, 
                                                               na.rm = T))], 
          na.rm = T)}, 
      error = function(x) NA)
    # select number of points that gave best RMSE in test with observed data
    # if multiple n points were equally good, use the smallest for fast computing
    npts <- tryCatch({
      min(performance$npoints[which(
        performance$RMSE == min(performance$RMSE, na.rm = T))], na.rm = T)}, 
      error = function(x) NA)
    
    ## make gstat object using inverse distance interpolation using power 
    ## selected above
    int_prop <-tryCatch({
      gstat(formula = prob_present ~ 1, 
            data = chkl_df, 
            nmin = min(npts, nrow(chkl_df)),
            nmax = npts,
            force = T, # use at least nmin points
            locations = ~ x + y, 
            set = list(idp = pwr))}, 
      error = function(x) NA)
    
    # predict probability of presence in each cell by interpolation
    predictions_df <- tryCatch({
      predict(int_prop, newdata = pred_locs)
    }, 
    error = function(c) {
      c$message <- paste0(c$message, 
                          ". Inverse distance interpolation failing with species ", 
                          sp_name, ". ")
    })
    
    try({
      colnames(predictions_df)[which(colnames(predictions_df) == "var1.pred")] <-
        "prob_present"
      predictions_df$en <- paste0(round(predictions_df$x), "_", 
                                  round(predictions_df$y))
      predictions_df <- tryCatch({full_join(chkl_df, predictions_df)}, 
                                 error = function(x) NA)
      predictions_df$x <- round(predictions_df$x)
      predictions_df$y <- round(predictions_df$y)
      # change column name to match that returned by other methods
      colnames(predictions_df)[which(colnames(predictions_df) == 
                                       "prob_present")] <- "predicted_prob_occ"
      predictions_df <- predictions_df[, -which(
        colnames(predictions_df) == "var1.var")]
    })
    
    ## plot for debugging:
    # print(ggplot(data = predictions_df, 
    #              aes(x = x, y = y, color = prob_present)) + 
    #         geom_point())
    
    ## return the gstat object and the data frame with predictions for each cell
    # using this fitted model.
    list(fit = list(int_prop), preds = list(predictions_df))
  }
)
## end inverse distance interpolation -----------------------------------------


## occupancy-detection model --------------------------------------------------
occ_det <- new_method(
  "occ_det", "Hierarchical occupancy-detection model", 
  method = function(model, draw) {
    ## outputs: - a list with 
    ##              * an unmarkedFit model object
    ##              * a (raster or data frame) with the predicted probability
    ##                of species presence in each cell for each species
    sp_name <- draw$sp_name
    chkl_df <- draw$observations 
    predictors <- model@params$env.predictors
    predictors$en <- paste(round(predictors$s1), round(predictors$s2), 
                           sep = "_")

    sp_df <- tryCatch({ # gather the column of observations for this sp.
      gather(chkl_df, key = "sp_name", value = "present", eval(sp_name))
    }, error = function(x) NA)
    sp_df$en <- tryCatch({ # make east/north character string column
      paste(round(sp_df$x), round(sp_df$y), sep = "_")
    }, error = function(x) NA)
    sp_df <- tryCatch({
      # join predictor variables to species observation df
      sp_df <- left_join(sp_df, predictors, by = "en")
    }, error = function(x) NA)
    
    # use random subset of variables for each species
    var_names <- names(predictors)[which(names(predictors) %nin% c("s1", "s2", 
                                                                   "en"))]
    var_names <- sample(var_names, size = 0.5*length(var_names), replace = F)
    
    ## fit full occupancy/detection model
    # unmarked needs to be fed a df with sites as rows and each column a 
    # replicate visit to that site, with the p/a value for a single species.
    # (This goes into the 'y' argument) 
    # Wrangel sp_df so that there is a single row for each site and a single
    # column for each visit

    # Assign each visit to a site a visit number to make spreading easier, 
    # then spread to make a df with a row for each site and a column for each
    # visit up to the maximum number of visits for any site.
    # I will keep a site names column here and remove it in the call to occu()
    sp_obs <- tryCatch(select(sp_df, present, en) %>%
                         group_by(en) %>%
                         mutate(visit = row_number(en)) %>%
                         mutate(visit = paste("v", visit, sep = ".")) %>%
                         spread(key = visit, value = present),
                       error = function(x) NA)
    # get rows in order to match s_covs
    sp_obs <- tryCatch(sp_obs[order(sp_obs$en), ], error = function(x) NA)
    
    # make site covariates df to feed into occu()
    s_covs <- tryCatch(select(sp_df, c(en, var_names)) %>%
                         distinct(), error = function(x) NA)
    # get rows in order to match sp_obs
    s_covs <- tryCatch(s_covs[order(s_covs$en), ], error = function(x) NA)
    
    # this basic occ_det method will not include any observation (visit-level) 
    # covariates, so I do not need to make a dataframe of observation covariates
    det_vars <- NA

    # make double right-hand side formula with covariates for detection 
    # and occupancy
    form <- tryCatch(as.formula(paste0("~ ", paste(" 1", names(det_vars), 
                                                   collapse = " + "), " ~ ", 
                                       paste(var_names, collapse = " + "), 
                                       " + ", 
                                       paste("poly(", var_names, ", degree = 2)",
                                             collapse = " + "))), 
                     error = function(x) NA)
    if(tryCatch(!identical(sp_obs$en, s_covs$en), error = function(x) FALSE)) {
      warning("sp_obs and s_covs do not have matching site rows in occ_det().")}
    
    # make unmarked frame
    ufo <- tryCatch(unmarkedFrameOccu(y = sp_obs[, 2:ncol(sp_obs)], 
                                      siteCovs = s_covs[, 2:ncol(s_covs)]),
                    error = function(c) {
                      c$message <- paste0(c$message, " in unmarkedFrameOccu().")
                    })

    # fit model
    # m <- tryCatch(occu(form, data = ufo, engine = "C" 
    #                      #, starts = c(0.5, 0.5, 0, 0, 0, 0)
    #                      # in sim_09b more models are fit without my starting values
    # ),
    # error = function(c) {c$message <- paste0(c$message, " in occu() call")
    # })

    all_vars <- c(var_names, paste0("poly(", var_names, ", degree = 2)"))
    m_base <- tryCatch(occu(~ 1 ~ 1, data = ufo, engine = "C"),
                            error = function(c) {c$message <- paste0(c$message, " in occu() call")
                            })
    
    m <- tryCatch(f.AICc.occu.sig(start.model = m_base, detocc = 1, 
                                  blocks = all_vars, max.iter = 30, 
                                  ufo = ufo), 
                  error = function(x) NA)


    # predict probabilities of presence in each grid cell from fitted model
    predictions_df <- predictors
    predictions_df$list_length <- rep(model@params$community.size, 
                                      nrow(predictions_df))
    preds <- tryCatch({
      do.call("predict", list(m, newdata = as.name("predictions_df"), 
                              type = "state"))
    }, 
    error = function(c) {
      c$message <- paste0(c$message, " predict do.call failing while predicting occurrences for species ", 
                          sp_name, ".")
    })
    try(colnames(preds) <- c("predicted_prob_occ", "SE", "pred_lower", 
                             "pred_upper"))
    if(is.data.frame(preds)){
      predictions_df <- tryCatch(bind_cols(predictions_df, preds), 
                                 error = function(x) NA)}
    
    predictions_df <- tryCatch(
      predictions_df[, names(predictions_df) %in% 
                       c("s1", "s2", "en", "predicted_prob_occ", "SE", 
                         "pred_lower", "pred_upper")], 
      error = function(x) NA)
    
    ## return the unmarkedFit model object and the data frame with predictions 
    # for each cell using this fitted model.
    list(fit = list(m), preds = list(predictions_df))
  }
)
## end occupancy-detection model ----------------------------------------------