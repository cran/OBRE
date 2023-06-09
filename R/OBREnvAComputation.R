#' OBRE vector a.
#'
#' The function evaluates integrals used to compute the components of OBRE a vector.
#'
#' @param nvData The vector of data.
#' @param nTheta1 The first parameter.
#' @param nTheta2 The second parameter.
#' @param lDensityExpr The list of symbolic expressions of density, cumulative and derivatives.
#' @param nCParOBRE OBRE c parameter.
#' @param matA OBRE matrix A.
#' @param nvA OBRE vector a.
#'
#' @return The OBRE a vector.
#'

OBREnvAComputation = function(nvData, nTheta1, nTheta2, lDensityExpr, nCParOBRE, matA, nvA) {

  nUpBound = Inf
  nLowBound = 0

# compute numerator for the component 1 ####
  nIntNum1 = quadinf(f = OBREnvANum1, xa = nLowBound, xb = nUpBound, nTheta1 = nTheta1, nTheta2 = nTheta2,
                     lDensityExpr = lDensityExpr, nCParOBRE = nCParOBRE, matA = matA,
                     nvA = nvA, tol = .Machine$double.eps)$Q

# compute numerator for the component 2 ####
  nIntNum2 = quadinf(f = OBREnvANum2, xa = nLowBound, xb = nUpBound, nTheta1 = nTheta1, nTheta2 = nTheta2,
                     lDensityExpr = lDensityExpr, nCParOBRE = nCParOBRE, matA = matA,
                     nvA = nvA, tol = .Machine$double.eps)$Q

# compute denominator ####
  nIntDen = quadinf(f = OBREnvADen, xa = nLowBound, xb = nUpBound, nTheta1 = nTheta1, nTheta2 = nTheta2,
                    lDensityExpr = lDensityExpr, nCParOBRE = nCParOBRE, matA = matA,
                    nvA = nvA, tol = .Machine$double.eps)$Q

# compure components ####
  nComponent1 = nIntNum1 / nIntDen
  nComponent2 = nIntNum2 / nIntDen

# generate output ####
  nvA = c(nComponent1, nComponent2)
  return(nvA)
}
