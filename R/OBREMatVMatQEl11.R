#' Element [1, 1] of matrix Q.
#'
#' Function computing argument element [1, 1] of matrix Q of asymptotic covariance matrix V.
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

OBREMatVMatQEl11 = function(nvData, nTheta1, nTheta2, lDensityExpr, nCParOBRE, matA, nvA) {
  # compute score functions components
  nvScore1 = scoreComponent(nvData = nvData, nTheta1 = nTheta1, nTheta2 = nTheta2,
                            lDensityExpr = lDensityExpr, nParIndex = 1)
  nvScore2 = scoreComponent(nvData = nvData, nTheta1 = nTheta1, nTheta2 = nTheta2,
                            lDensityExpr = lDensityExpr, nParIndex = 2)
  # compute weights function
  nvWeights = OBREWeightsFun(nvData = nvData, nTheta1 = nTheta1, nTheta2 = nTheta2, lDensityExpr = lDensityExpr,
                             nCParOBRE = nCParOBRE, matA = matA, nvA = nvA)
  # compute the part of the element in common to both shifted and truncated distributions
  nvEl11matQ = ((nvScore1 - nvA[1]) * nvWeights ) ^ 2
  nvEl11matQ = nvEl11matQ * eval(expr = lDensityExpr$eDensityFun)
  # substitute NA values with zero values
  nvEl11matQ[is.na(nvEl11matQ)] = 0
  return(nvEl11matQ)
}
