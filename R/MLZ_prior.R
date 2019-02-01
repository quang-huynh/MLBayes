#' @rdname MLZBayes-package
#' @name MLZBayes-package
#' @aliases MLZBayes-package MLZBayes
#' @title \packageTitle{MLZBayes}
#' @description \packageDescription{MLZBayes}
#' @details An overview of how the package works
#' @import stats graphics


setClassUnion("numeric_integer", c("integer", "numeric"))


#' Class-\code{MLZ_prior}
#'
#' An S4 class that contains priors and model structure input. See example for template code.
#'
#' @name MLZ_prior-class
#' @docType class
#'
#' @slot ncp The number of change points in mortality.
#' @slot Z_dist Character. The probability distribution for the mortality rate. Currently supports \code{lognormal} or \code{uniform}.
#' Defaults to \code{uniform}.
#' @slot Z_par A matrix with \code{ncp+1} rows and 2 columns. Parameters describing the priors for the mortality rate. 
#' If lognormal, the first column is the mean and the second column is the standard deviation of the prior (in normal space).
#' If uniform, the first column is the minimum and the second column is the maximum of the distribution. Defaults to \code{uniform} with bounds 0 and 3.
#' @slot alpha The alpha vector of the Dirichlet distribution of length \code{ncp+1} for parameterization the priors for the change points 
#' (when \code{ncp > 0}). Ignored if \code{ncp = 0}. Defaults to \code{rep(1, ncp + 1)} which is a relatively uninformative prior.
#' @slot sigma_dist The probability distribution for the standard deviation. Currently supports \code{lognormal} or \code{uniform}.
#' Defaults to \code{uniform}.
#' @slot sigma_par A vector of length 2. Parameters describing the priors for the standard deviation. 
#' If lognormal, the mean and standard deviation of the prior (in normal space).
#' If uniform, the first column is the minimum and the second column is the maximum of the distribution. Defaults to \code{uniform} with bounds 0 and 1e5.
#' @author Q. Huynh
#' @export
#' @exportClass MLZ_prior
#' @importClassesFrom MLZ MLZ_data
MLZ_prior <- setClass("MLZ_prior",
                      slots = c(ncp = "numeric_integer", Z_dist = "character", Z_par = "matrix",
                                alpha = "numeric", sigma_dist = "character", sigma_par = "numeric"))

setMethod("initialize", "MLZ_prior", 
          function(.Object, ...) {
            .Object@Z_dist <- "uniform"
            .Object@sigma_dist <- "uniform"
            dots <- list(...)
            if (length(dots) > 0) {
              for (i in 1:length(dots)) {
                slot(.Object, names(dots)[[i]]) <- dots[[i]]
              }
            }
            
            # Default prior
            if (!is.null(.Object@ncp)) { 
              if (length(.Object@Z_par) == 0) {
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
              if (.Object@ncp > 0 && length(.Object@alpha) == 0) {
                .Object@alpha <- rep(1, .Object@ncp + 1)
              }
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

