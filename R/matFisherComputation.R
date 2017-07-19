#' Fisher information matrix
#'
#' Function calculating the Fisher information matrix.
#' @keywords internal
#'
#' @param nTheta1 First parameter.
#' @param nTheta2 Second parameter.
#' @param lDensityExpr List of symbolic expressions of density, cumulative and derivatives.
#'
#' @return The Fisher information matrix.
#'

matFisherComputation = function(nTheta1, nTheta2, lDensityExpr) {

# inizialization of bounds ####
  nUpBound = Inf
  nLowBound = 0

# computation of integrals ####
  # Integral for element (1,1)
  nFisherInt11Part1 = quadinf(f = fisherEl11Part1, xa = nLowBound, xb = nUpBound,
                							nTheta1 = nTheta1, nTheta2 = nTheta2,
                							lDensityExpr = lDensityExpr, tol = .Machine$double.eps)$Q

  nFisherInt11Part2 = quadinf(f = fisherEl11Part2, xa = nLowBound, xb = nUpBound,
                							nTheta1 = nTheta1, nTheta2 = nTheta2,
                							lDensityExpr = lDensityExpr, tol = .Machine$double.eps)$Q
  # Integral for element (2,2)
  nFisherInt22Part1 = quadinf(f = fisherEl22Part1, xa = nLowBound, xb = nUpBound,
                							nTheta1 = nTheta1, nTheta2 = nTheta2,
                							lDensityExpr = lDensityExpr, tol = .Machine$double.eps)$Q

  nFisherInt22Part2 = quadinf(f = fisherEl22Part2, xa = nLowBound, xb = nUpBound,
                							nTheta1 = nTheta1, nTheta2 = nTheta2,
                							lDensityExpr = lDensityExpr, tol = .Machine$double.eps)$Q
  # Integral for element (1,2)
  nFisherInt12Part1 = quadinf(f = fisherEl12Part1, xa = nLowBound, xb = nUpBound,
                              nTheta1 = nTheta1, nTheta2 = nTheta2,
                              lDensityExpr = lDensityExpr, tol = .Machine$double.eps)$Q

  nFisherInt12Part2 = quadinf(f = fisherEl12Part2, xa = nLowBound, xb = nUpBound,
                              nTheta1 = nTheta1, nTheta2 = nTheta2,
                              lDensityExpr = lDensityExpr, tol = .Machine$double.eps)$Q

# computation of matrix elements ####
  nFisherEl11 = -(nFisherInt11Part1 + nFisherInt11Part2)
  nFisherEl22 = -(nFisherInt22Part1 + nFisherInt22Part2)
  nFisherEl12 = -(nFisherInt12Part1 + nFisherInt12Part2)

# build the output matrix ####
  matFisher = matrix(data = c(nFisherEl11, nFisherEl12, nFisherEl12, nFisherEl22),
                     nrow = 2, ncol = 2)
  return(matFisher)
}
