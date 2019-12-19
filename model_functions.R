############################################
## This creates the model function that defines how to simulate data.
## This is for simulation 10.10 (the first version using the simulator package)
##
## author: Willson Gaul wgaul@hotmail.com
## created: 24 Sep 2018
## last modified: 5 Oct 2018
############################################

## @knitr models

# make_my_model <- function(n, prob) {
#   new_model(name = "contaminated-normal",
#             label = sprintf("Contaminated normal (n = %s, prob = %s)", n, prob),
#             params = list(n = n, mu = 2, prob = prob),
#             simulate = function(n, mu, prob, nsim) {
#               # this function must return a list of length nsim
#               contam <- runif(n * nsim) < prob
#               x <- matrix(rep(NA, n * nsim), n, nsim)
#               x[contam] <- rexp(sum(contam))
#               x[!contam] <- rnorm(sum(!contam))
#               x <- mu + x # true mean is mu
#               return(split(x, col(x))) # make each col its own list element
#             })
# }

make_bio_recs_model <- function(sp.list, community.size, n.obs, n.obs.reference, 
                                shape, polynomial.resp, 
                                polynomial.coef.min, polynomial.coef.max, 
                                sl.intercept.min, sl.intercept.max, 
                                sl.coef.min, sl.coef.max, pca, bias.rasters, 
                                bias.name, 
                                list.lengths, env.predictors, 
                                error.prob = 0, ...) {

  list2env(list(...), envir = environment())
  
  ## generate model
  new_model(name = "species_model", 
            label = sprintf("Species model (comm.size = %s, n.obs = %s, bias = %s, n.obs.ref = %s)", 
                            community.size, n.obs, bias.name, n.obs.reference), 
            params = list(community.size = community.size, 
                          n.obs = n.obs, 
                          n.obs.reference = n.obs.reference, 
                          shape = shape, 
                          polynomial.resp = polynomial.resp, 
                          polynomial.coef.min = polynomial.coef.min, 
                          polynomial.coef.max = polynomial.coef.max,
                          sl.intercept.min = sl.intercept.min, 
                          sl.intercept.max = sl.intercept.max, 
                          sl.coef.min = sl.coef.min, 
                          sl.coef.max = sl.coef.max,
                          pca = pca, 
                          bias.rasters = bias.rasters, 
                          bias.name = bias.name, 
                          list.lengths = list.lengths, 
                          env.predictors = env.predictors,
                          squared.predictors = squared.predictors, 
                          randomPoints.replace = randomPoints.replace, 
                          error.prob = error.prob, 
                          on.sonic = on.sonic), 
            simulate = simulate_fun
  )
}
