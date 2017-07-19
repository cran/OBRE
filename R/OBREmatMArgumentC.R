#' Argument C for OBRE matrix M integrals.
#'
#' Function computing argument C for OBRE matrix M integrals.
#' @keywords internal
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

OBREmatMArgumentC = function(nvData, nTheta1, nTheta2, lDensityExpr, nCParOBRE, matA, nvA, nK) {

  nvArgC = (eval(expr = lDensityExpr$eDerivDensityFunTheta[1]) / eval(expr = lDensityExpr$eDensityFun) - nvA[1]) *
    (eval(expr = lDensityExpr$eDerivDensityFunTheta[2]) / eval(expr = lDensityExpr$eDensityFun) - nvA[2]) *
    ((OBREWeightsFun(nvData = nvData, nTheta1 = nTheta1, nTheta2 = nTheta2, lDensityExpr = lDensityExpr,
                     nCParOBRE = nCParOBRE, matA = matA, nvA = nvA)) ^ nK) *
    eval(expr = lDensityExpr$eDensityFun)
  nvArgC[is.na(nvArgC)] = 0
  return(nvArgC)
}
