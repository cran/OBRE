#' Argument B for OBRE matrix M integrals.
#'
#' Function computing argument B for OBRE matrix M integrals.
#'
#' @param nvData The vector of data.
#' @param nTheta1 The first parameter.
#' @param nTheta2 The second parameter.
#' @param lDensityExpr List of symbolic expressions of density, cumulative and derivatives.
#' @param nCParOBRE OBRE c parameter.
#' @param matA Matrix A.
#' @param nvA Vector a.
#' @param nK Exponent which differentiate M_1 from M_2.
#'

OBREmatMArgumentB = function(nvData, nTheta1, nTheta2, lDensityExpr, nCParOBRE, matA, nvA, nK) {

  nvArgB = ((eval(expr = lDensityExpr$eDerivDensityFunTheta[2]) / eval(expr = lDensityExpr$eDensityFun) - nvA[2])^2) *
    ((OBREWeightsFun(nvData = nvData, nTheta1 = nTheta1, nTheta2 = nTheta2, lDensityExpr = lDensityExpr,
                     nCParOBRE = nCParOBRE, matA = matA, nvA = nvA)) ^ nK) *
    eval(expr = lDensityExpr$eDensityFun)
  nvArgB[is.na(nvArgB)] = 0
  return(nvArgB)
}
