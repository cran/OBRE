#' Element [1, 2] of matrix M.
#'
#' Function computing element [1, 2] of matrix M, for the computation of asymptotic covariance matrix V.
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

OBREMatVMatMEl12 = function(nvData, nTheta1, nTheta2, lDensityExpr, nCParOBRE, matA, nvA) {

  # compute score functions components
  nvScore1 = scoreComponent(nvData = nvData, nTheta1 = nTheta1, nTheta2 = nTheta2,
                            lDensityExpr = lDensityExpr, nParIndex = 1)
  nvScore2 = scoreComponent(nvData = nvData, nTheta1 = nTheta1, nTheta2 = nTheta2,
                            lDensityExpr = lDensityExpr, nParIndex = 2)
  # compute weights function
  nvWeights = OBREWeightsFun(nvData = nvData, nTheta1 = nTheta1, nTheta2 = nTheta2, lDensityExpr = lDensityExpr,
                             nCParOBRE = nCParOBRE, matA = matA, nvA = nvA)
  # compute the part of the element in common to both shifted and truncated distributions
  nvEl12matM = (nvScore1 - nvA[1]) * nvWeights * (nvScore2)
  nvEl12matM = nvEl12matM * eval(expr = lDensityExpr$eDensityFun)
  # substitute NA values with zero values
  nvEl12matM[is.na(nvEl12matM)] = 0
  return(nvEl12matM)
}
