#' Bayesian mean length estimator
#'
#' @param MLZ_data An object with mean length data. Only slots \code{Year}, \code{ML}, and \code{ss} are used.
#' @param MLZ_prior An object with describing model structure (including priors and the number of change points).
#' @param chains The number of MCMC chains.
#' @param iter The number of MCMC iterations.
#' @param warmup The number of warmup iterations.
#' @param thin The number for which each \code{thin}-th iteration in the chain is saved.
#' @param cores The number of cores for running the chains in parallel.
#' @param ... Arguments to pass to \link[rstan]{sampling}.
#' 
#' @note By default, uses the No U-Turn Sampling algorithm for MCMC sampling. The default settings for the number of 
#' total iterations (\code{iter}), warmup iterations (\code{warmup}), and thinning (\code{thin}) are much more conservative
#' than the rstan defaults. More iterations and thinning may still be needed. Always check for convergence.
#' @seealso \link[rstan]{sampling}
#' @export
ML_stan <- function(MLZ_data, MLZ_prior, chains = 2, iter = 6e4, warmup = 1e4, thin = 20, cores = chains, ...) {
  
  ncp <- as.integer(MLZ_prior@ncp)
  if (ncp == 0L) {
    stan_obj <- stanmodels$MLeq
  } else {
    stan_obj <- stanmodels$ML
  }
  
  res <- sampling(stan_obj, data = NULL, chains = chains, iter = iter, warmup = warmup, thin = thin, cores = cores, ...)
  return(res)
}
