#' Bayesian mean length estimator
#'
#' @param ... Arguments to pass through
#'
#' @export
ML_stan <- function(...) {
  res <- sampling(stanmodels$ML, ...)
  return(res)
}
