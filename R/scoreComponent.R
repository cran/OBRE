#' First component of the score function.
#'
#' The function evaluates the formula used to compute the first component of the score function. The
#' missing elements are imputed with 0.
#' @keywords internal
#'
#' @param nvData The vector of data.
#' @param nTheta1 The first parameter.
#' @param nTheta2 The second parameter.
#' @param lDensityExpr The list of symbolic expressions of density, cumulative and derivatives.
#' @param nParIndex Which component parameter needs to be calculated.
#'
#' @return The first component of the score function.
#'

scoreComponent = function(nvData, nTheta1, nTheta2, lDensityExpr, nParIndex) {

# Computation ####
  nvScore1 = ((eval(expr = lDensityExpr$eDerivDensityFunTheta[nParIndex])) / eval(expr = lDensityExpr$eDensityFun))
  # Subtitute the NA values with zero
  nvScore1[is.na(nvScore1)] = 0
  return(nvScore1)
}
