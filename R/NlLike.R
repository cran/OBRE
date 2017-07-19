#' Negative Log-Likelihood
#'
#' The function compute the Negative Log-Likelihood value that has to be used for optimization in MLE function.
#' @keywords internal
#'
#' @param nvTheta Parameters of the distribution.
#' @param nvData The vector of the data.
#' @param lDensityExpr The distribution density espressions.
#'
#' @return Negative log likelihood value.
#'

NlLike = function (nvTheta, nvData, lDensityExpr) {

# assigning parameters ####
  nTheta1 = nvTheta[1]
  nTheta2 = nvTheta[2]
  nvDensity = eval(expr = lDensityExpr$eDensityFun)
  nvLogDensity = suppressWarnings(log(nvDensity))
  nvNegLogLH = try(-sum(nvLogDensity), silent = TRUE)
  if (is.na(nvNegLogLH) || is.nan(nvNegLogLH)) {
    nvNegLogLH = .Machine$double.xmax
  }
  return(nvNegLogLH)
}
