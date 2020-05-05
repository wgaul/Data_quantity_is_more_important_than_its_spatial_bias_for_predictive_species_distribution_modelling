############################################
## This defines the species distribution modeling methods to be tested.
## Simulation 10.10 (the first version using the simulator package)
## 
## author: Willson Gaul wgaul@hotmail.com
## created: 24 Sep 2018
## last modified: 30 March 2020
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
    sp_name <- draw$sp_name
    chkl_df <- draw$observations 
    predictors <- model@params$predictors.for.models
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
    # Use a random subset of half of the predictors each time  
    var_names <- names(predictors)[which(names(predictors) %nin% c("s1", "s2", 
                                                                   "en"))]
    # remove squared terms so I can select half the predictors
    var_names <- gsub("_sq", "", var_names) 
    var_names <- unique(var_names)
    var_names <- sample(var_names, size = 0.5*length(var_names), replace = F)
    
    # add squared terms back in
    var_names <- c(var_names, paste0(var_names, "_sq"))
    
    # make all possible combinations of terms
    n0 <- length(var_names)
    ind1 <- unlist(lapply(1:n0, function(x) combn(1:n0, x, simplify = F)), 
                   recursive = F)
    all_mods <- lapply(ind1, function(x) {
      paste(var_names[x], collapse = " + ")})
    
    # remove models that have a 2nd degree term but no 1st degree term for any
    # variable.  These are not valide models.
    invalid_mods <- rep_len(FALSE, length.out = length(all_mods))
    for(trm in var_names[grepl(".*_sq", var_names)]) {
      t1 <- gsub("_sq", " ", trm) # get the 1st degree term for this variable
      for(j in 1:length(all_mods)) {
        if(grepl(trm, all_mods[[j]]) & !grepl(t1, all_mods[[j]])){
          invalid_mods[j] <- TRUE}
      }
    }
    all_mods <- all_mods[!invalid_mods] # this is now a list of all valid models
  
    ## remove models that have more terms than we have data to fit.  
    # I will use the rule of thumb that each term needs at least 10 events or 
    # non-events (whichever is smaller).  
   
    # find smaller of number of events or non-events
    nevents <-tryCatch(min(as.numeric(table(sp_df$present))), 
                       error = function(x) NA)
    # find which models have number of terms < number of events / 10
    keep_mods <- sapply(all_mods, function(x, nevents = nevents) {
      nterms <- tryCatch(length(strsplit(x, "+", fixed = TRUE)[[1]]), 
                         error = function (x) NA)
      tryCatch({if(nterms <= nevents/10) {TRUE} else FALSE}, 
               error = function(x) FALSE)
    }, nevents = nevents)
    
    # Only fit models if more than 5 events/nonevents in the training data
    if(!is.na(nevents) & nevents > 5) {
      # keep only models with few enough predictors to fit given n events
      all_mods <- tryCatch(all_mods[keep_mods], error = function(x) NA)
      if(length(all_mods) > 0) {
        all_mods <- lapply(all_mods, function(x) {
          tryCatch(paste0("present ~ 1 + ", x), error = function(x) NA)})
      }
      
      # add the intercept-only model
      try(all_mods[[length(all_mods) + 1]] <- "present ~ 1") 
      
      # make vector to hold aic values for all models
      aic_vec <- tryCatch(rep_len(NA, length.out = length(all_mods)), 
                          error = function(x) NA)
      
      # Fit all models and record AIC, keeping the best (lowest AIC) model
      for(mi in 1:length(all_mods)) {
        form <- tryCatch(as.formula(all_mods[[mi]]), 
                         error = function (x) NA)
        # fit glm model
        this_mod <- tryCatch({
          do.call("glm", list(form, data = as.name("sp_df"), family = "binomial", 
                              model = FALSE, x = FALSE))}, 
          error = function(c) {
            c$message <- paste0(c$message, 
                                " glm do.call failing while fitting glm to species ", 
                                sp_name, ".")})
        # put AIC for this model in vector
        aic_vec[mi] <- tryCatch(AIC(this_mod), error = function(x) NA)
        if(!is.na(aic_vec[mi]) & aic_vec[mi] == min(aic_vec, na.rm = T)) {
          m <- this_mod}
      }
      
      if(!exists("m")) m <- "No model fit."
      
      # predict probabilities of presence in each grid cell from best fitted model
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
    } else {
      m <- "No model fit because there were fewer than 6 detections."
      predictions_df <- "No model fit because there were fewer than 6 detections."
      all_mods <- "No model fit because there were fewer than 6 detections."
      aic_vec <- "No model fit because there were fewer than 6 detections."
    }
    
    # return the model object and the data frame with predictions for each cell
    # using this fitted model.
    list(fit = list(m), preds = list(predictions_df), 
         models_tested = list(models = all_mods, aics = aic_vec))
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
    
    # find smaller of number of events or non-events
    nevents <-tryCatch(min(as.numeric(table(sp_df$present))), 
                       error = function(x) NA)
    
    # Only fit models if more than 5 events/nonevents in the training data
    if(!is.na(nevents) & nevents > 5) {
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
    } else {
      m <- "No model fit because there were fewer than 6 detections."
      predictions_df <- "No model fit because there were fewer than 6 detections."
    }
    
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
    max_ntree <- 30000 # maximum number of trees to allow
    tc1 <- 2 # tree complexity
    tc2 <- 5 # alternate tree complexity to try

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
    
    # find smaller of number of events or non-events
    nevents <-tryCatch(min(as.numeric(table(chkl_df$present))), 
                       error = function(x) NA)
    
    # Only fit models if more than 5 events/nonevents in the training data
    if(!is.na(nevents) & nevents > 5) {
      fit_model <- TRUE # indicator for while loop
      # while an adequate model has not been fit, try to fit a model
      while (fit_model == TRUE) { 
        m1 <- tryCatch({gbm.step(chkl_df, 
                                 gbm.x = var_names,
                                 gbm.y = "present", 
                                 tree.complexity = tc1, 
                                 learning.rate = lr, 
                                 bag.fraction = 0.5, 
                                 n.trees = 50,
                                 step.size = 50, 
                                 max.trees = max_ntree,
                                 family = "bernoulli", 
                                 n.folds = 10, 
                                 prev.stratify = TRUE, 
                                 plot.main = FALSE)
        }, error = function(x) NA)
        # If model fit successfully with more than 1000 trees and fewer than 
        # max.trees, keep this model and
        # do not try any other learning rates, as lower learning rates are 
        # generally prefered as I understand it (see vignette('gbm')). 
        if(!is.null(m1) && !is.na(m1) && m1$n.trees >= 1000 && 
           m1$n.trees <= max_ntree) {
          fit_model <- FALSE
        }
        
        # if model fitting stopped because max.trees was reached, then discard
        # this model and fit a new one, unless the learning rate is higher than
        # 0.1, in which case give up. 
        if(fit_model == TRUE) {
          # if too few trees (< 1000) were used, make learning rate smaller
          # (until learning rates as small as 0.00001 have been tested, after
          # which point if the model still has not fit then abandon it.)
          if(!is.null(m1) && !is.na(m1) && m1$n.trees < 1000 && lr > 0.00001) {
            lr <- lr - lr*0.5
          } else { 
            if(!is.null(m1) && !is.na(m1)) {
              # presumably more than max trees were used if I'm in this if 
              # clause, in which case learning rate needs to be made bigger.
              lr <- lr + 0.002
              if(lr > 0.1) {
                fit_model <- FALSE
                m <- "Brt model not fit with fewer than the maximum number of trees and a learning rate less than 0.1 (see vignette('gbm'))."
              }
            } else {
              # if model is NULL, try reducing learning rate to see if it will
              # fit, unless learning rate is already below 0.00001, in which
              # case abandon this attempt.
              if(is.null(m1) && lr > 0.00001) {
                lr <- lr - lr*0.5
              } else {
                fit_model <- FALSE
                m <- "Brt model not fit, perhaps because learning rate and step size were never made small enough."
              }
            }
          }
        }
      } # end while loop
      
      ## Try tree complexity of 5
      lr <- 0.001 # reset learning rate to smallest to try 1st 
      fit_model <- TRUE # indicator for while loop
      # while an adequate model has not been fit, try to fit a model
      while (fit_model == TRUE) { 
        m2 <- tryCatch({gbm.step(chkl_df, 
                                 gbm.x = var_names,
                                 gbm.y = "present", 
                                 tree.complexity = tc2, 
                                 learning.rate = lr, 
                                 bag.fraction = 0.5, 
                                 n.trees = 50,
                                 step.size = 50, 
                                 max.trees = max_ntree,
                                 family = "bernoulli", 
                                 n.folds = 10, 
                                 prev.stratify = TRUE, 
                                 plot.main = FALSE)
        }, error = function(x) NA)
        # If model fit successfully with more than 1000 trees and fewer than 
        # max.trees, keep this model and
        # do not try any other learning rates, as lower learning rates are 
        # generally prefered as I understand it (see vignette('gbm')). 
        if(!is.null(m2) && !is.na(m2) && m2$n.trees >= 1000 && 
           m2$n.trees <= max_ntree) {
          fit_model <- FALSE
        }
        
        # if model fitting stopped because max.trees was reached, then discard
        # this model and fit a new one, unless the learning rate is higher than 
        # 0.1, in which case give up. 
        if(fit_model == TRUE) {
          # if too few trees (< 1000) were used, make learning rate smaller
          # (until learning rates as small as 0.00001 have been tested, after
          # which point if the model still has not fit then abandon it.)
          if(!is.null(m2) && !is.na(m2) && m2$n.trees < 1000 && lr > 0.00001) {
            lr <- lr - lr*0.5
          } else { 
            if(!is.null(m2) && !is.na(m2)) {
              # presumably more than max trees were used if I'm in this if 
              # clause, in which case learning rate needs to be made bigger.
              lr <- lr + 0.002
              if(lr > 0.1) {
                fit_model <- FALSE
                m <- "Brt model not fit with fewer than the maximum number of trees and a learning rate less than 0.1 (see vignette('gbm'))."
              }
            } else {
              # if model is NULL, try reducing learning rate to see if it will
              # fit, unless learning rate is already below 0.00001, in which
              # case abandon this attempt.
              if(is.null(m2) && lr > 0.00001) {
                lr <- lr - lr*0.5
              } else {
                fit_model <- FALSE
                m <- "Brt model not fit, perhaps because learning rate and step size were never made small enough."
              }
            }
          }
        }
      } # end second while loop (for tree complexity = 5)
      
      ## Use the model with the tree complexity that was best
      if(!is.null(m1) && !is.na(m1) && !is.null(m2) && !is.na(m2)) {
        # if both models fit, if m1 has lower CV deviance than m2, keep m1.
        if(m1$cv.statistics$deviance.mean < m2$cv.statistics$deviance.mean) {
          m <- m1
        } else m <- m2 # otherwise, keep m2 (higher tree complexity)
      } else {
        # if only one of the two models fit, use the one that fit
        if(!is.null(m1) && !is.na(m1)) {m <- m1}
        if(!is.null(m2) && !is.na(m2)) {m <- m2}
      }
      
      # If no model has been selected as the model m for some reason, put
      # an error message in m.  
      if(!exists("m")) m <- "No brt model fit for some reason."
      
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
    } else {
      m <- "No model fit because there were fewer than 6 detections."
      predictions_df <- "No model fit because there were fewer than 6 detections."
    }
    
    ## return the model object and the data frame with predictions for each cell
    # using this fitted model.
    list(fit = list(m), preds = list(predictions_df))
  }
)
## end boosted regression tree -------------------------------------------------


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
    
    # find smaller of number of events or non-events
    nevents <-tryCatch(min(as.numeric(table(chkl_df$present))), 
                       error = function(x) NA)
    
    # Only fit models if more than 5 events/nonevents in the training data
    if(!is.na(nevents) & nevents > 5) {
      # function to calculate prop. of checklists on which the species appears
      calc_prop <- function(x, list_id) {sum(x) / length(unique(list_id))}
      
      chkl_df <- tryCatch({
        # summarise observations to the grid cell level by calculating the 
        # proportion of all checklists in the grid cell that have the species.
        # This should produce a df with a single row for each grid cell that has
        # a checklist. No value will be produced for grid cells with no checklists
        group_by(chkl_df, en) %>%
          summarise(prob_present = calc_prop(x = present, 
                                             list_id = checklist_ID), 
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
    } else {
      int_prop <- "No model fit because there were fewer than 6 detections."
      predictions_df <- "No model fit because there were fewer than 6 detections."
    }
    
    ## return the gstat object and the data frame with predictions for each cell
    # using this fitted model.
    list(fit = list(int_prop), preds = list(predictions_df))
  }
)
## end inverse distance interpolation -----------------------------------------
