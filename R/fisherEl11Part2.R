#' Part 2 of element [1, 1] for Fisher Information matrix
#'
#' Function computing part 2 of element [1, 1] for Fisher Information matrix computation. The Fisher Information matrix is
#' splitted in the four elements ([1, 1], [1, 2], [2, 1], [2, 2]). Each element is split in part 1 and part 2
#'
#' @param nvData The vector of data.
#' @param nTheta1 The first parameter.
#' @param nTheta2 The second parameter.
#' @param lDensityExpr List of symbolic expressions of density, cumulative and derivatives.
#'

fisherEl11Part2 = function(nvData, nTheta1, nTheta2, lDensityExpr) {
  # computation of the argument
  nvArg11Part2 = eval(expr = lDensityExpr$eDeriv2DensityFunTheta1Theta1)
  # substitute NA values with zero values
  nvArg11Part2[is.na(nvArg11Part2)] = 0
  return(nvArg11Part2)
}
