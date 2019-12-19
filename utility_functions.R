##############################
## Define utility functions for simulating species and drawing samples from the 
## simulated distributions.  
## Simulation 10.10
##
## author: Willson Gaul wgaul@hotmail.com
## created: 26 Sep 2018
## last modified: 17 Sep 2019
##############################

## this is the function from virtualspecies, but it isn't exported from that 
## package so I need to create it here.
wg_randomPoints <- function(mask, n, p, ext=NULL, extf=1.1, excludep=TRUE, 
                          prob=FALSE, cellnumbers=FALSE, tryf=3, warn=2, 
                          lonlatCorrection=TRUE,
                          replaceCells = FALSE) {
  
  if (nlayers(mask) > 1) { 
    mask <- raster(mask, 1)	
  }
  
  tryf <- max(round(tryf[1]), 1)
  
  if (missing(p)) { 
    excludep <- FALSE
  } else {
    if (inherits(p, 'SpatialPoints')) {
      p <- sp::coordinates(p)
    }
  }
  
  if (class(ext)=='character') {
    if (! ext %in% c('points')) { 
      stop("if ext is a character variable it should be 'points'") 
    } else if (missing(p)) { 
      warning("if p is missing, 'ext=points' is meaningless") 
      ext <- extent(mask)  
    } else {
      ext <- extent(min(p[,1]), max(p[,1]), min(p[,2]), max(p[,2]))
    }
  } 
  
  if (! is.null(ext)) {
    ext <- extent(ext)
    ext <- ext * extf
    ext <- intersect(ext, extent(mask))
    mask2 <- crop(raster(mask), ext)
  }  else {
    mask2 <- raster(mask)
  }
  
  nn <- n * tryf
  nn <- max(nn, 10)
  
  if (prob) {
    stopifnot(hasValues(mask))
    cells <- crop(mask, mask2)
    cells <- try( stats::na.omit(cbind(1:ncell(cells), getValues(cells))))
    if (class(cells) == 'try-error') {
      stop("the raster is too large to be used with 'prob=TRUE'")
    }
    prob <- cells[,2]
    cells <- cells[,1]
    if (couldBeLonLat(mask)) {
      rows <- rowFromCell(mask2, cells)
      y <- yFromRow(mask2, rows)
      dx <- pointDistance(cbind(0, y), cbind(xres(mask2), y), longlat=TRUE)  
      prob <- prob * dx
    }
    
    cells <- sample(cells, nn, prob=prob, replace = replaceCells)
    xy <- xyFromCell(mask2, cells)
    cells <- cellFromXY(mask, xy)
    
  } else 	if (canProcessInMemory(mask2)) {
    
    cells <- crop(mask, mask2)
    if (hasValues(cells)) {
      cells <- which(! is.na(getValues(cells)) )
    } else {
      cells <- 1:ncell(cells)
    }
    nn <- min(length(cells), nn)
    
    if (lonlatCorrection & couldBeLonLat(mask)) {
      # which rows are that?
      rows <- rowFromCell(mask2, cells)
      # what is the latitude?
      y <- yFromRow(mask2, rows)
      # what is the 'width' of a cell?
      dx <- pointDistance(cbind(0, y), cbind(xres(mask2), y), longlat=TRUE)  
      cells <- sample(cells, nn, prob=dx, replace = replaceCells)
      
    } else {
      cells <- sample(cells, nn, replace = replaceCells)
    }
    xy <- xyFromCell(mask2, cells)
    cells <- cellFromXY(mask, xy)
    
  } else {
    
    nn <- min(ncell(mask2), nn)
    if (couldBeLonLat(mask2)) {
      
      cells <- .randomCellsLonLat(mask2, nn)
      
    } else {
      if (nn >= ncell(mask2)) {
        cells <- 1:ncell(mask2)
      } else {
        cells <- sampleInt(ncell(mask2), nn, replace = replaceCells)
      }
      
    }
    xy <- xyFromCell(mask2, cells)
    cells <- cellFromXY(mask, xy)
    
    if (hasValues(mask)) {
      
      vals <- cbind(cells, extract(mask, cells))
      cells <- stats::na.omit(vals)[,1]
    }
  }
  
  if (excludep) {	
    pcells <- cellFromXY(mask, p)
    cells <- cells[ ! (cells %in% pcells) ] 	
  }
  
  if (length(cells) > n) { 
    
    cells <- sample(cells, n) 
    
  } else if (length(cells) < n) {
    
    frac <- length(cells) / n
    if (frac < 0.1) {
      stop("generated random points = ", frac," times requested number; Use a higher value for tryf" )
    }
    if (frac < 0.5  & warn==1) {
      warning("generated random points = ", frac," times requested number; Use a higher value for tryf" )
    } else if (warn > 1) {
      warning("generated random points = ", frac," times requested number")
    }
  }
  
  if (cellnumbers) {
    return(cells)
  } else {
    return(xyFromCell(mask, cells))
  }
}

### twentieth root function ---------------------------------------------------
twentiethrt <- function(x) {
  # calculate a tenth root
  sign(x) * abs(x)^(1/20)
}
### end twentieth root function ----------------------------------------------


### stepwise model selection for unmakred models ------------------------------
f.AICc.occu.sig <- function(start.model, blocks, max.iter=NULL, detocc = 1, AICcut = 1, print.log=TRUE, ...){
  # A function for "stepwise" regression using occupancy models of package unmarked 
  # 
  # ARGS:   * start.model: initial model e.g. occu(~1~1, UMF)
  #         * detocc: if set to 1 (default) runs the function on the occupancy side; 
  #                   2 does the detectability
  #         * blocks: variables to add incrementally to the model
  #         * max.iter: don't know exactly what this is for
  #         * AICcut: don't know exactly what this is for
  #         * print.log: logical?  Presumably this has to do with printing 
  #         
  # VALUE:  * The best model (lowest AIC) as selected by this forward stepwise
  #           process.
  #         
  # NOTES: 
  # This code by R. Schuster downloaded on 24 May 2019 by Willson Gaul
  # For details see https://groups.google.com/forum/#!topic/unmarked/XuI8kdxkqQo
  # 
  # See also Schuster & Arcese (2012) in Ecography
  # Parts based on Forward.lmer by Rense Nieuwenhuis 
  # (http://www.rensenieuwenhuis.nl/r-sessions-32/)
  # with some additions by for Nick Isaac 
  #
  # Author: Richard Schuster (mail@richard-schuster.com 
  # 08 October 2012
  # last modified: 27 May 2019 W. Gaul wgaul@hotmail.com
  # 
  list2env(list(...), envir = environment())
 
  # start.model <- eval(start.model)
  modlst <- c(start.model)
  x <- 2
  if (detocc == 1) coeff <- length(start.model@estimates@estimates$state@estimates)
  else coeff <- length(start.model@estimates@estimates$det@estimates)                                                                               
  best <- FALSE
  model.basis <- start.model
  keep <- list(start.model)
  AICmin <- AIClst <- AICc(start.model)
  #cutoff for when to exclude values
  cutoff <- 20	
  
  # Maximum number of iterations cannot exceed number of blocks, but this is also the default
  if (is.null(max.iter) | max.iter > length(blocks)) max.iter <- length(blocks)
  
  # Setting up the outer loop
  for(ii in 1:max.iter){
    models <- list()
    coeff <- coeff + 1
    cnt <- 1
    for (jj in 1:length(keep)) {  
      # Iteratively updating the model with addition of one block of variable(s)
      for(kk in 1:length(blocks)) {
        
        if (detocc == 1) form <- as.formula(paste("~. ~. + ", blocks[kk]))
        else form <- as.formula(paste("~. + ", blocks[kk], "~. "))
        
        if (class(dummy <- try(update(keep[[jj]], form))) == "unmarkedFitOccu") {
          flag <- 0
          
          #check if there is any NAN's in the SE's
          if (detocc == 1) {
            for(dd in 1:length(sqrt(diag (vcov(dummy@estimates@estimates$state))))) {
              if(diag (vcov(dummy@estimates@estimates$state))[[dd]] < 0
                 || sqrt(diag (vcov(dummy@estimates@estimates$state))[[dd]]) > cutoff)
                flag <- 1
              break
            }
          }
          else {
            for(dd in 1:length(sqrt(diag (vcov(dummy@estimates@estimates$det))))) {
              if(diag (vcov(dummy@estimates@estimates$det))[[dd]] < 0
                 || sqrt(diag (vcov(dummy@estimates@estimates$det))[[dd]]) > cutoff)
                flag <- 1
              break
            }
          } 
          if (flag == 0) { 
            for (bb in 1:length(AIClst)) {
              if (round(AICc(dummy),digits=6) == round(AIClst[bb], digits=6)) {
                flag <- 1
                break
              }
            }
          }          
          
        }
        else {flag <- 1}
        
        if (flag == 0) {
          models[[cnt]] <- dummy
          modlst[[x]] <- models[[cnt]]
          AIClst <- c(AIClst, AICc(models[[cnt]]))
          x <- x + 1
          cnt <- cnt + 1
        } 
      }
    }
    
    if(length(LL <- unlist(lapply(models, function(x) {AICc(x)}))) == 0) break
    keep <- list()
    k = 1
    cont <- 0
    #check for improvement in AIC, if none stop loop
    for (mm in order(LL, decreasing=FALSE)) {
      if (LL[mm]< AICmin + AICcut) {
        if (detocc == 1) {
          # 27 May 2019.  WG changed the following if statement to be true
          # if there are more estimates than coeff as well as if the number of
          # estimates is the same as coeff.  This is because coeff in the loops
          # above gets only "1" added to it even when a "poly()" term is the 
          # term added from blocks, but the "poly()" term adds two estimates to
          # the fitted model.  
          # I do not know why this if statement was here in the original
          # function.  I have not read through it closely enough to understand
          # what conditions would lead to the number of estimates being 
          # different from coef.  So maybe I have messed things up here, I don't
          # really know. 
          if(length(models[[mm]]@estimates@estimates$state@estimates) >= coeff) {
            keep[[k]]<- models[[mm]]
            k <- k + 1
            if (LL[mm]<AICmin) {
              AICmin <- LL[mm]
              icont <- 1
            }
          }
        }
        else {
          stop("As of 27 May 2019. This is not fixed to accomodate poly() because I am not optimizing the detetion portion of the model.  See above for how to deal with poly.  wg.")
          if (length(models[[mm]]@estimates@estimates$det@estimates) == coeff) {
            keep[[k]]<- models[[mm]]
            k <- k + 1
            if (LL[mm]<AICmin) {
              AICmin <- LL[mm]
              cont <- 1
            }
          }
        }
      }
      else break
    }
    
    if (length(keep) == 0) break                                      
  }  
  ## Get best model
  aics <- sapply(modlst, FUN = function(x) x@AIC) # get AIC values from models
  best_mod <- modlst[aics == min(aics)][[1]] # select model with smallest AIC 
  ## Return the best model
  return(best_mod)
}
### end stepwise model selection for unmakred models ---------------------------

############ eveness ---------------------------------------------------------
###
### Camargo evenness (Magurran & McGill 2011 Box 5.1 p. 57)
## "The highest possible evenness is when p_i = p_j = 1/S, so Camargo et al. 
## (1993) suggested a direct measurement of deviation from this ideal.  
## E = 1 - sum over(|p_i - p_j|/S) where the sum is taken over i = 1...S and 
## j = i+1...S
## For my situation, instead of abundance and species, I will have number of 
## records and hectads.  So, 
## E = 1 -sum over(|p_i - p_i+1|/nhec) where p_i is the proportion of records
## in hectad i (nrec_i/total_nrec)

camargo <- function(x) {
  # Calculate Camargo evenness (Camargo et al 1993) 
  # E = 1 - (sum over(|p_i - p_i+1|/nhec))
  # 
  # ARGS: x - a vector of abundances (or number of records in my case)
  N = length(x) 
  x = x/sum(x) # convert abundances to proportions
  # calculate absolute value of differences between hectad i and hec i+1
  difs = matrix(nrow = N, ncol = N)
  for(i in 1:N) {
    for(j in 1:N) {
      difs[i, j] <- abs(x[i] - x[j])
    }
  }
  E = 1 - (sum(difs)/(2*N))
  E
}

### Shannon evenness (Magurran & McGill 2011 Box 5.1 p. 57)
### "The highest value of Shannon Diversity when all species are equally 
### abundant [is] ln(S) so dividing by ln(s) will give an index from 0 to 1
### E = Shannon Diversity / ln(S)
### Shannon Diversity (entropy H) = -sum over p_i*ln(p_i)  
### OR  the finite population size version "Brillouin's index".  
### Brillouin's index = (1/N)log[(N!)/product of all (N_i!)]
### where N = sum of all abundances N_i

shannon_even <- function(x) {
  # Calculate Shannon Evenness 
  # E = Shannon Diversity / ln(S)  
  # where S is number of species (or hectads in my case) and 
  # Shannon Diversity (H) = -sum over p_i*ln(p_i)
  # OR substitute Brillouin's index for finite population size in place of 
  # Shannon Diversity (H)
  # E = (1/N)log[(N!)/product of all (N_i!)] / ln(S)
  # where N = sum of all abundances N_i
  # 
  # 
  # ARGS: x - a vector of abundances (or number of records in my case)
  S = length(x)
  N = sum(x)
  # with high abundances (number of records) N! gets too big for R I think.
  if(is.finite(suppressWarnings(factorial(N)))) {
    prod_Ni_fact = prod(sapply(x, factorial)) # product of all N_i!
    brillouin_d = (1/N) * log(factorial(N) / (prod_Ni_fact))
    E = brillouin_d / log(S)
  } else {
    x_p = x/sum(x)
    # When N_i is 0, instead of calculating the log(0) for shannon H,
    # set H to 0 since that is the limit of H is - as p goes to 0
    shannon_H = -sum(sapply(x_p, 
                            FUN = function(p) {
                              p[p == 0] = rep(0, length(which(p == 0))) 
                              p[p != 0] = p[p != 0] *log(p[p != 0])
                              p}))
    E = shannon_H / log(S)
  }
  E
} 

simpson_even <- function(x) {
  ### Simpson diversity (Magurran & McGill 2011 Box 5.1 p. 57)
  ### "D = sum over p_i^2 gave the probability that two individuals drawn at 
  ### random from an infinite community would belong to the same species.  ...
  ### As such, D is the inverse of diversity.... The most common way of converting
  ### homogeneity into diversity is D = 1/D."  Simpson evenness is scaled 0 to 1
  ### by dividing Simpson diversity by species richness (or number of hectads in 
  ### my case).
  # Simpson_e = (1 / sumOver(p_i)^2) / S
  # where S is number of observed species (or number of hectads) and p_i is the
  # proportion of abundance for species i (p_i = N_i/N)
  S = length(x)
  x_p = x/sum(x)
  simpsons_d = 1 / sum(sapply(x_p, function(x) {x^2}))
  E = simpsons_d / S
  E
}