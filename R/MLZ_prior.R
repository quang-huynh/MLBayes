


setClassUnion("numeric_integer", c("integer", "numeric"))

#' Class-\code{MLZ_prior}
#'
#' An S4 class that contains priors and model structure input.
#'
#' @name MLZ_prior-class
#' @docType class
#'
#' @slot ncp The number of change points in mortality.
#' @slot Z_dist The probability distribution for the mortality rate. Currently supports lognormal.
#' @slot Z_par Parameters describing the priors for the mortality rate.
#' @slot cp_par The alpha vector for the change points (uses Dirichlet formulation, see vignette).
#' @slot sigma_dist The probability distribution for the standard deviation. Currently supports lognormal.
#' @slot sigma_par Parameters describing the priors for the standard deviation.
#' @author Q. Huynh
#' @export
#' @exportClass MLZ_prior
#' @importClassesFrom MLZ MLZ_data
MLZ_prior <- setClass("MLZ_prior",
                      slots = c(ncp = "numeric_integer", Z_dist = "character", Z_par = "matrix",
                                cp_par = "numeric", sigma_dist = "character", sigma_par = "numeric"))

setMethod("initialize", "MLZ_prior", 
          function(.Object, ...) {
            .Object@Z_dist <- "uniform"
            .Object@sigma_dist <- "uniform"
            dots <- list(...)
            if (length(dots) > 0) {
              for (i in 1:length(dots)) slot(.Object, names(dots)[[i]]) <- dots[[i]]
            }
            
            # Default prior
            if (!is.null(.Object@ncp)) { 
              if (length(.Object@Z_par == 0)) {
                .Object@Z_par <- matrix(NA, nrow = .Object@ncp + 1, ncol = 2)
                if (.Object@Z_dist == "lognormal") {
                  .Object@Z_par[, 1] <- 0.5
                  .Object@Z_par[, 2] <- 1
                }
                if (.Object@Z_dist == "uniform") {
                  .Object@Z_par[, 1] <- 0
                  .Object@Z_par[, 2] <- 3
                }
              }
              if (.Object@ncp > 0 && length(.Object@cp_par == 0)) .Object@cp_par <- rep(1e-3, .Object@ncp + 1)
            }
            if (length(.Object@sigma_par) == 0) {
              if (.Object@sigma_dist == "lognormal") .Object@sigma_par <- c(1, 1)
              if (.Object@sigma_dist == "uniform") .Object@sigma_par <- c(0, 1e5)
            }
            
            # Error check
            #if (length(.Object@ncp) > 1) warning("slot ncp should be an integer.")
            #if (!is.null(.Object@ncp)) {
            #  if(length)
            #}
            
            .Object
          })

