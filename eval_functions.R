############################################
## This defines the evaluation methods to be used for testing model performance.
## Simulation 10.10 (the first version using the simulator package)
## 
## author: Willson Gaul wgaul@hotmail.com
## created: 24 Sep 2018
## last modified: 12 Sep 2019
############################################

## @knitr metrics

## AUC (can take random or block cross-validation, or a single data set with
## no test/training split)
auc <- new_metric(
  "auc", "AUC", 
  metric = function(model, draw, out) {
    cv_aucs <- c() # initialize object to hold AUC values for each fold
    test_data <- as.data.frame(draw$truth_maps[[1]], xy = TRUE)
    test_data <- test_data[complete.cases(test_data), ]
    test_data$en <- paste0(round(test_data$x), "_", 
                           round(test_data$y))

    # if multiple folds exist, evaluate models against test data set
    if(length(out$fit) > 1) {
      for(i in 1:length(out$fit)) {
        # subset test data to points in test set for this fold
        tdat <- tryCatch({test_data[which(test_data$en %in% 
                                            out$folds_en[[i]]$test_en), ]}, 
                         error = function(x) NA)
        # subset predictions to points in test set for this fold
        out$preds[[i]] <- tryCatch({
          out$preds[[i]][which(out$preds[[i]]$en %in% 
                                 out$folds_en[[i]]$test_en), ]}, 
          error = function(x) NA)
        if(is.data.frame(out$preds[[i]])) {
          names(out$preds[[i]])[which(names(out$preds[[i]]) == "s1")] <- "x"
          names(out$preds[[i]])[which(names(out$preds[[i]]) == "s2")] <- "y"
          tdat <- tdat[order(tdat$en), ]
          out$preds[[i]] <- out$preds[[i]][order(out$preds[[i]]$en), ]
          try(if(!identical(tdat$x, out$preds[[i]]$x) | 
                 !identical(tdat$y, out$preds[[i]]$y)) {
            warning("test_data and preds dataframes don't have matching coordinates in auc calculation.")
          })
        }
        cv_aucs[[i]] <- tryCatch(
          roc(response = factor(tdat[, "layer"], 
                                levels = c("0", "1")),
              predictor = out$preds[[i]]$predicted_prob_occ,
              auc = T)$auc[1], 
          error = function(x) NA)
      }
    }
    # if only one fold exists (no test-training split), then evaluate models
    # against all data
    if(length(out$fit) == 1) {
      for(i in 1:length(out$fit)) {
        # use all data as test data, and all predictions for comparison 
        tdat <- tryCatch({test_data}, error = function(x) NA)
        if(is.data.frame(out$preds[[i]])) {
          names(out$preds[[i]])[which(names(out$preds[[i]]) == "s1")] <- "x"
          names(out$preds[[i]])[which(names(out$preds[[i]]) == "s2")] <- "y"
          tdat <- tdat[order(tdat$en), ]
          out$preds[[i]] <- out$preds[[i]][order(out$preds[[i]]$en), ]
          try(if(!identical(as.character(tdat$x), 
                            as.character(out$preds[[i]]$x)) | 
                 !identical(as.character(tdat$y), 
                            as.character(out$preds[[i]]$y))) {
            warning("test_data and preds dataframes don't have matching coordinates in auc calculation.")
          })
        }
        cv_aucs[[i]] <- tryCatch(
          roc(response = factor(tdat[, "layer"], 
                                levels = c("0", "1")),
              predictor = out$preds[[i]]$predicted_prob_occ,
              auc = T)$auc[1], 
          error = function(x) NA)
      }
    }
    # return mean auc for all folds (even if only one fold)
    mean_auc <- tryCatch({mean(cv_aucs, na.rm = TRUE)}, error = function(x) NA)
    mean_auc
  }
)


## RMSE (can take random or block cross-validation, or a single data set with
## no test/training split)
rmse <- new_metric(
  "rmse", "RMSE", 
  metric = function(model, draw, out) {
    cv_rmses <- c() # initialize object to hold RMSE values for each fold
    test_data <- as.data.frame(draw$prob.raster, xy = TRUE)
    test_data <- test_data[complete.cases(test_data), ]
    test_data$en <- paste0(round(test_data$x), "_", 
                           round(test_data$y))
    
    # if multiple folds exist, evaluate models against test data set
    if(length(out$fit) > 1) {
      for(i in 1:length(out$fit)) {
        # subset test data to points in test set for this fold
        tdat <- tryCatch({test_data[which(test_data$en %in% 
                                            out$folds_en[[i]]$test_en), ]}, 
                         error = function(x) NA)
        # subset predictions to points in test set for this fold
        out$preds[[i]] <- tryCatch({
          out$preds[[i]][which(out$preds[[i]]$en %in% 
                                 out$folds_en[[i]]$test_en), ]}, 
          error = function(x) NA)
        if(is.data.frame(out$preds[[i]])) {
          names(out$preds[[i]])[which(names(out$preds[[i]]) == "s1")] <- "x"
          names(out$preds[[i]])[which(names(out$preds[[i]]) == "s2")] <- "y"
          tdat <- tdat[order(tdat$en), ]
          out$preds[[i]] <- out$preds[[i]][order(out$preds[[i]]$en), ]
          try(if(!identical(tdat$x, out$preds[[i]]$x) | 
                 !identical(tdat$y, out$preds[[i]]$y)) {
            warning("test_data and preds dataframes don't have matching coordinates in rmse calculation.")
          })
        }
        errs <- tryCatch(out$preds[[i]]$predicted_prob_occ - 
                           as.numeric(as.character(tdat[, "layer"])), 
                   error = function(x) NA)
        errs_sq <- tryCatch(errs^2, error = function(x) NA)
        
        cv_rmses[[i]] <- tryCatch(
          sqrt(mean(errs_sq, na.rm = T)), error = function(x) NA)
      }
    }
    # if only one fold exists (no test-training split), then evaluate models
    # against all data
    if(length(out$fit) == 1) {
      for(i in 1:length(out$fit)) {
        # use all data as test data, and all predictions for comparison 
        tdat <- tryCatch({test_data}, error = function(x) NA)
        if(is.data.frame(out$preds[[i]])) {
          names(out$preds[[i]])[which(names(out$preds[[i]]) == "s1")] <- "x"
          names(out$preds[[i]])[which(names(out$preds[[i]]) == "s2")] <- "y"
          tdat <- tdat[order(tdat$en), ]
          out$preds[[i]] <- out$preds[[i]][order(out$preds[[i]]$en), ]
          try(if(!identical(as.character(tdat$x), 
                            as.character(out$preds[[i]]$x)) | 
                 !identical(as.character(tdat$y), 
                            as.character(out$preds[[i]]$y))) {
            warning("test_data and preds dataframes don't have matching coordinates in rmse calculation.")
          })
        }
        errs <- tryCatch(out$preds[[i]]$predicted_prob_occ - 
                           as.numeric(as.character(tdat[, "layer"])), 
                         error = function(x) NA)
        errs_sq <- tryCatch(errs^2, error = function(x) NA)
        
        cv_rmses[[i]] <- tryCatch(
          sqrt(mean(errs_sq, na.rm = T)), error = function(x) NA)
      }
    }
    # return mean auc for all folds (even if only one fold)
    mean_rmse <- tryCatch({mean(cv_rmses, na.rm = TRUE)}, error = function(x) NA)
    mean_rmse
  }
)
