#' Denominator for nvA
#'
#' Function computing denominator for OBRE numeric vector nvA evaluation.
#'
#' @param nvData The vector of data.
#' @param nTheta1 The first parameter.
#' @param nTheta2 The second parameter.
#' @param lDensityExpr List of symbolic expressions of density, cumulative and derivatives.
#' @param nCParOBRE OBRE c parameter.
#' @param matA Matrix A.
#' @param nvA Vector a.
#'

OBREnvADen = function(nvData, nTheta1, nTheta2, lDensityExpr, nCParOBRE, matA, nvA) {

  nvDen = OBREWeightsFun(nvData = nvData, nTheta1 = nTheta1, nTheta2 = nTheta2, lDensityExpr = lDensityExpr,
                         nCParOBRE = nCParOBRE, matA = matA, nvA = nvA) *
    eval(expr = lDensityExpr$eDensityFun)
  nvDen[is.na(nvDen)] = 0
  return(nvDen)
}
