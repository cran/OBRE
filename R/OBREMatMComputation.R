#' Function computing the OBRE matrix M.
#'
#' The function evaluates integrals used to compute the M_1 and M_2 OBRE matrices. Element (1,1) uses
#' argument (A,B,F); element (1,2) uses argument (B,D,E,F); elements (2,2) uses arguments (C,D,F).
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
#' @return OBRE M matrix (M_1 if nK = 1; M_2 if nK = 2).
#'

OBREMatMComputation = function(nvData, nTheta1, nTheta2, lDensityExpr, nCParOBRE, matA, nvA, nK){

  nUpBound = Inf
  nLowBound = 0

# components for mat M computation ####
  nIntA = quadinf(f = OBREmatMArgumentA, xa = nLowBound, xb = nUpBound, tol = .Machine$double.eps,
                  nTheta1 = nTheta1, nTheta2 = nTheta2, lDensityExpr = lDensityExpr,
                  nCParOBRE = nCParOBRE, matA = matA, nvA = nvA, nK = nK)$Q

  nIntB = quadinf(f = OBREmatMArgumentB, xa = nLowBound, xb = nUpBound, tol = .Machine$double.eps,
                  nTheta1 = nTheta1, nTheta2 = nTheta2, lDensityExpr = lDensityExpr,
                  nCParOBRE = nCParOBRE, matA = matA, nvA = nvA, nK = nK)$Q

  nIntC = quadinf(f = OBREmatMArgumentC, xa = nLowBound, xb = nUpBound, tol = .Machine$double.eps,
                  nTheta1 = nTheta1, nTheta2 = nTheta2, lDensityExpr = lDensityExpr,
                  nCParOBRE = nCParOBRE, matA = matA, nvA = nvA, nK = nK)$Q
  # computation of element [1, 1]
  nEl11 = nIntA
  # computation of element [2, 2]
  nEl22 = nIntB
  # computation of element [1, 2]
  nEl12 = nIntC

# build the output matrix ####
  matM = matrix(c(nEl11, nEl12, nEl12, nEl22), nrow = 2, ncol = 2)

  return(matM)
}
