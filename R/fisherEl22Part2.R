#' Part 2 of element [2, 2] for Fisher Information matrix
#'
#' Function computing part 2 of element [2, 2] for Fisher Information matrix computation. The Fisher Information matrix is
#' splitted in the four elements ([1, 1], [1, 2], [2, 1], [2, 2]). Each element is split in part 1 and part 2
#'
#' @param nvData The vector of data.
#' @param nTheta1 The first parameter.
#' @param nTheta2 The second parameter.
#' @param lDensityExpr List of symbolic expressions of density, cumulative and derivatives.
#'

fisherEl22Part2 = function(nvData, nTheta1, nTheta2, lDensityExpr) {
  # computation of the argument
  nvArg22Part2 = eval(expr = lDensityExpr$eDeriv2DensityFunTheta2Theta2)
  # substitute NA values with zero values
  nvArg22Part2[is.na(nvArg22Part2)] = 0
  return(nvArg22Part2)
}
