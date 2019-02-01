
#' @export
setMethod("plot", signature(x = "stanfit", y = "missing"), 
          function(x, subplot = c("par", "ts"), interval = 0.95, ...) {
            subplot <- match.arg(subplot, several.ok = TRUE)
            
            old_par <- par(no.readonly = TRUE)
            on.exit(par(old_par))
            
            if(any(subplot == "par")) make_parameter_plot(x)
            if(any(subplot == "ts")) make_ts_plot(x, interval)
            invisible()
            
          })

#' Plot function for stanfit objects in MLZBayes
#' 
#' This function produces up two plots. The first plot draws density plots for priors and posteriors
#' of model parameters (dotted lines indicate median). The second plot is a time series plot of observed (black) and predicted 
#' (red, with median in bolded lines and confidence intervals in dotted lines) mean lengths and estimated Z.
#' 
#' @param x A stanfit object from \link{ML_stan}.
#' @param y An optional stanfit object for the same model as \code{x} but only sampled the priors. See example.
#' @param subplot A character vector of plots to draw. (\code{"par"} for parameter priors/posteriors, 
#' \code{"ts"} for time series of mean length and Z.
#' @param interval The probability coverage interval for the time series plots.
#' @seealso \link{ML_stan}
#' @importClassesFrom rstan stanfit
#' @examples 
#' \dontrun{
#' library(MLZ)
#' data(Goosefish)
#' 
#' # Create an object with priors for a model with 2 change points in mortality (ncp = 2)
#' uninformative_priors <- new("MLZ_prior", ncp = 2)
#' 
#' # Run the MCMC (calls rstan::sampling)
#' res <- ML_stan(Goosefish, uninformative_priors)
#' plot(res)
#' 
#' # Test priors only
#' res2 <- ML_stan(Goosefish, uninformative, prior_only = TRUE)
#' plot(res2)
#' 
#' # Overlay posteriors on top of priors
#' plot(res, res2) 
#' }
#' @export
setMethod("plot", signature(x = "stanfit", y = "stanfit"), 
          function(x, y, subplot = c("par", "ts"), interval = 0.95, ...) {
            subplot <- match.arg(subplot, several.ok = TRUE)
            
            old_par <- par(no.readonly = TRUE)
            on.exit(par(old_par))
            
            if(any(subplot == "par")) make_parameter_plot(x, y)
            if(any(subplot == "ts")) make_ts_plot(x, interval)
            invisible()
            
          })

generate_plot_data <- function(x, var = "Lpred", interval = NULL) {
  probs <- c(0.5 - 0.5 * interval, 0.5, 0.5 + 0.5 * interval)
  data_subset <- lapply(x, function(z) z[grep(var, names(z))])
  n_data <- length(data_subset[[1]])
  
  reorg_data <- list()
  for(i in 1:n_data) reorg_data[[i]] <- lapply(data_subset, getElement, i)
  reorg_data <- lapply(reorg_data, function(z) do.call(c, z))
  
  if(is.null(interval)) return(reorg_data) else {
    reorg_data <- lapply(reorg_data, quantile, probs = probs)
  }
  
  return(do.call(cbind, reorg_data))
}

make_parameter_plot <- function(x, y = NULL) {
  MLZ_data <- get("MLZ_data", envir = x@.MISC)
  MLZ_prior <- get("MLZ_prior", envir = x@.MISC)
  ncp <- MLZ_prior@ncp
  
  Z_post <- generate_plot_data(x@sim$samples, "Z")[1:(ncp+1)]
  p_post <- generate_plot_data(x@sim$samples, "p")[1:(ncp+1)]
  D_post <- generate_plot_data(x@sim$samples, "D")
  sigma_post <- generate_plot_data(x@sim$samples, "sigma")
  post_list <- c(Z_post, p_post, D_post, sigma_post)
  
  warmup_ind <- 1:(x@sim$warmup %/% x@sim$thin)
  post_dens <- lapply(post_list, function(z) density(z[-warmup_ind], from = min(z), to = max(z)))
  
  pr_dens <- new("list", rep(1, length(post_dens)))
  if(!is.null(y)) {
    Z_pr <- generate_plot_data(y@sim$samples, "Z")[1:(ncp+1)]
    p_pr <- generate_plot_data(y@sim$samples, "p")[1:(ncp+1)]
    D_pr <- generate_plot_data(y@sim$samples, "D")
    sigma_pr <- generate_plot_data(y@sim$samples, "sigma")
    pr_list <- c(Z_pr, p_pr, D_pr, sigma_pr)
    pr_dens <- lapply(pr_list, function(z) density(z[-warmup_ind], from = min(z), to = max(z)))
  }
  
  npar <- length(post_dens)
  nrow_par <- min(3, ceiling(npar/3))
  xlab <- c(paste0("Z[", 1:(ncp+1), "]"), paste0("p[", 1:(ncp+1), "]"),
            paste0("D[", 1:ncp, "]"), "sigma")
  
  par(mfrow = c(nrow_par, 3), mar = c(4, 3, 1, 1),
      oma = c(ifelse(is.null(y), 0, 2), 2, 0, 0))
  for(i in 1:npar) {
    plot(post_dens[[i]]$x, post_dens[[i]]$y, typ = "n", ylab = "", xlab = xlab[i])
    abline(h = 0, col = "grey")
    lines(post_dens[[i]]$x, post_dens[[i]]$y, lwd = 3, col = "red")
    abline(v = median(post_list[[i]]), col = "red", lty = 2)
    
    if(!is.null(y)) {
      lines(pr_dens[[i]]$x, pr_dens[[i]]$y, lwd = 3)
      abline(v = median(pr_list[[i]]), lty = 2)
    }
  }
  mtext("Density", side = 2, outer = TRUE)
  
  if(!is.null(y)) {
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    
    plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
    legend("bottom", legend = c("Prior", "Posterior"), col = c("black", "red"), lwd = 3, bty = "n", horiz = TRUE)
  }
  
  invisible()
}

make_ts_plot <- function(x, interval) {
  MLZ_data <- get("MLZ_data", envir = x@.MISC)
  
  Lpred <- generate_plot_data(x@sim$samples, "Lpred", interval)
  ylab <- ifelse(nchar(MLZ_data@length.units) > 0, paste0("Mean Length (", MLZ_data@length.units, ")"),
                 "Mean Length")
  ylim <- c(0.9, 1.1) * range(c(Lpred, MLZ_data@MeanLength), na.rm = TRUE)
  
  par(mfrow = c(2, 1), mar = c(3, 4, 1, 1), oma = c(2, 0, 0, 0))
  plot(MLZ_data@Year, MLZ_data@MeanLength, typ = "o", pch = 16,
       ylab = ylab, ylim = ylim)
  lines(MLZ_data@Year, Lpred[2, ], col = "red", lwd = 3)
  lines(MLZ_data@Year, Lpred[1, ], col = "red", lty = 3)
  lines(MLZ_data@Year, Lpred[3, ], col = "red", lty = 3)
  
  Z_yr <- generate_plot_data(x@sim$samples, "Z_yr", interval)
  ylim <- c(0, 1.1) * max(Z_yr)
  plot(MLZ_data@Year, Z_yr[2, ], typ = "l", col = "red", lwd = 3,
       ylab = "Total Mortality Z", ylim = ylim)
  lines(MLZ_data@Year, Z_yr[1, ], col = "red", lty = 3)
  lines(MLZ_data@Year, Z_yr[3, ], col = "red", lty = 3)
  
  mtext("Year", side = 1, outer = TRUE)
  
  invisible()
  
}