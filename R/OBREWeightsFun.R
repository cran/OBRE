#' OBRE weights.
#'
#' Function for computing OBRE weights. The function computes the score function for both parameters and
#' build the score matrix. The score matrix is then modified using OBRE parameters A matrix and a vector and
#' an euclidean norm is derived. The weights are finally found as the minimum between the normalized nCParOBRE
#' and 1.
#'
#' @param nvData The vector of data.
#' @param nTheta1 The first parameter.
#' @param nTheta2 The second parameter.
#' @param lDensityExpr The list of symbolic expressions of density, cumulative and derivatives.
#' @param nCParOBRE OBRE c parameter.
#' @param matA OBRE matrix A.
#' @param nvA OBRE vector a.
#'
#' @return A numeric vector containing OBRE weights.
#'

OBREWeightsFun = function(nvData, nTheta1, nTheta2, lDensityExpr, nCParOBRE, matA, nvA) {

# Compute score components ####
  nvScoreComponent1 = scoreComponent(nvData = nvData, nTheta1 = nTheta1, nTheta2 = nTheta2,
                                     lDensityExpr = lDensityExpr, nParIndex = 1)
  nvScoreComponent2 = scoreComponent(nvData = nvData, nTheta1 = nTheta1, nTheta2 = nTheta2,
                                     lDensityExpr = lDensityExpr, nParIndex = 2)
  # Gather the two components in one matrix
  matScore = matrix( data = c(nvScoreComponent1, nvScoreComponent2), nrow = 2,
                     ncol = length(nvScoreComponent1), byrow = TRUE)

# Compute the argument of the norm and then the euclidean norm ####
  nvArgumentNorm = matA %*% (matScore - nvA)
  nvNorm = sqrt(nvArgumentNorm[1, ] ^ 2 + nvArgumentNorm[2, ] ^ 2)

# Find the weight function as the minimum between the normalized OBRE c parameter and 1 ####
  nvWeightFun = pmin(1, nCParOBRE / nvNorm)
  return(nvWeightFun)
}
