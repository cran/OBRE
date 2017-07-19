#' Element [2, 2] of matrix Q.
#'
#' Function computing argument element [2, 2] of matrix Q of asymptotic covariance matrix V.
#' @keywords internal
#'
#' @param nvData The vector of data.
#' @param nTheta1 The first parameter.
#' @param nTheta2 The second parameter.
#' @param lDensityExpr List of symbolic expressions of density, cumulative and derivatives.
#' @param nCParOBRE OBRE c parameter.
#' @param matA Matrix A.
#' @param nvA Vector a.
#'

OBREMatVMatQEl22 = function(nvData, nTheta1, nTheta2, lDensityExpr, nCParOBRE, matA, nvA) {
  # compute score functions components
  nvScore1 = scoreComponent(nvData = nvData, nTheta1 = nTheta1, nTheta2 = nTheta2,
                            lDensityExpr = lDensityExpr, nParIndex = 1)
  nvScore2 = scoreComponent(nvData = nvData, nTheta1 = nTheta1, nTheta2 = nTheta2,
                            lDensityExpr = lDensityExpr, nParIndex = 2)
  # compute weights function
  nvWeights = OBREWeightsFun(nvData = nvData, nTheta1 = nTheta1, nTheta2 = nTheta2, lDensityExpr = lDensityExpr,
                             nCParOBRE = nCParOBRE, matA = matA, nvA = nvA)
  # compute the part of the element in common to both shifted and truncated distributions
  nvEl22matQ = ((nvScore2 - nvA[2]) * nvWeights) ^ 2
  nvEl22matQ = nvEl22matQ * eval(expr = lDensityExpr$eDensityFun)
  # substitute NA values with zero values
  nvEl22matQ[is.na(nvEl22matQ)] = 0
  return(nvEl22matQ)
}
