
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
#' @slot cp_par The alpha for the change points (uses Dirichlet formulation, see vignette).
#' @slot sigma_dist The probability distribution for the standard deviation. Currently supports lognormal.
#' @slot sigma_par Parameters describing the priors for the standard deviation.
#' @author Q. Huynh
#' @export
#' @exportClass MLZ_prior
#' @importClassesFrom MLZ MLZ_data
MLZ_prior <- setClass("MLZ_prior",
                      slots = c(ncp = "integer", Z_dist = "character", Z_par = "matrix",
                                cp_par = "numeric", sigma_dist = "character", sigma_par = "numeric"))
