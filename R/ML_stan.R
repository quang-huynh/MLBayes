#' Bayesian mean length estimator
#'
#' @param MLZ_data An object with mean length data. Slots \code{Year}, \code{MeanLength}, \code{ss}, and \code{Lc} are used.
#' @param MLZ_prior An object with describing model structure (including priors and the number of change points).
#' @param prior_only Logical. If TRUE, only the priors will be sampled in the MCMC.
#' @param chains The number of MCMC chains.
#' @param iter The number of MCMC iterations.
#' @param warmup The number of warmup iterations.
#' @param thin The number for which each \code{thin}-th iteration in the chain is saved.
#' @param seed Integer for random number generation.
#' @param cores The number of cores for running the chains in parallel.
#' @param ... Other arguments to pass to \link[rstan]{sampling}.
#' 
#' @note By default, uses the No U-Turn Sampling algorithm for MCMC sampling. The default settings for the number of 
#' total iterations (\code{iter}), warmup iterations (\code{warmup}), and thinning (\code{thin}) are much more conservative
#' than the rstan defaults. More iterations and thinning may still be needed. Always check for convergence.
#' @return An object of class \linkS4class{stanfit}. 
#' @seealso \link[rstan]{sampling}
#' @importMethodsFrom MLZ summary 
#' @export
#' @examples
#' 
#' \dontrun{ 
#' # Create MLZ_data object with mean length data 
#' # Use utility functions in MLZ package
#' libray(MLZ)
#' data(Goosefish)
#' 
#' # Create an object with priors for a model with 2 change points in mortality (ncp = 2)
#' uninformative_priors <- new("MLZ_prior", ncp = 2)
#' 
#' # Run the MCMC (calls rstan::sampling)
#' res <- ML_stan(Goosefish, uninformative_priors)
#' 
#' # All generics for stanfit object from rstan package are available.
#' stan_trace(res)
#' stan_dens(res, separate_chains = TRUE)
#' 
#' # Test priors only
#' res2 <- ML_stan(Goosefish, uninformative, prior_only = TRUE)
#' stan_dens(res2, separate_chains = TRUE)
#' }
#' 
ML_stan <- function(MLZ_data, MLZ_prior, prior_only = FALSE, chains = 2, iter = 6e4, warmup = 1e4, thin = 20, seed = 20, cores = chains, ...) {
  ncp <- as.integer(MLZ_prior@ncp)
  data_summary <- summary(MLZ_data)
  if(MLZ_prior@Z_dist == "uniform") Z_dist <- 0L
  if(MLZ_prior@Z_dist == "lognormal") Z_dist <- 1L
  if(MLZ_prior@sigma_dist == "uniform") sigma_dist <- 0L
  if(MLZ_prior@sigma_dist == "lognormal") sigma_dist <- 1L
  
  ss <- data_summary$ss
  if(prior_only) ss <- rep(0, length(ss))
  if (ncp == 0L) {
    stan_obj <- stanmodels$MLeq
    stan_data <- list(count = nrow(data_summary), Lobs = data_summary$MeanLength, ss = ss,
                      Lc = MLZ_data@Lc, Linf = MLZ_data@vbLinf, K = MLZ_data@vbK, Z_dist = Z_dist, Z_par = MLZ_prior@Z_par,
                      sigma_dist = sigma_dist, sigma_par = MLZ_prior@sigma_par, prior_only = prior_only)
  } else {
    stan_obj <- stanmodels$ML
    stan_data <- list(count = nrow(data_summary), nbreaks = ncp, Lobs = data_summary$MeanLength, ss = ss,
                      Lc = MLZ_data@Lc, Linf = MLZ_data@vbLinf, K = MLZ_data@vbK, Z_dist = Z_dist, Z_par = MLZ_prior@Z_par,
                      alpha_dirichlet = MLZ_prior@alpha, sigma_dist = sigma_dist, sigma_par = MLZ_prior@sigma_par, prior_only = prior_only)
  }
  
  res <- sampling(stan_obj, data = stan_data, chains = chains, iter = iter, warmup = warmup, thin = thin, 
                  seed = seed, cores = cores, ...)
  
  # Convert yearZ to calendar year
  Year1 <- data_summary$Year[1]
  new_samples <- lapply(res@sim$samples, add_yearZ, Year1)
  res@sim$samples <- new_samples
  
  return(res)
}

add_yearZ <- function(x, Year1) {
  ind <- grep("yearZ", names(x))
  x[ind] <- lapply(x[ind], function(y) y + Year1 - 1)
  return(x)
}

