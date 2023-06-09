#' Check if OBRE matrix A and vector a are final.
#'
#' The function compute the relative distance from the past to the current iteration of matrix A, with
#' respect to the relative tolerance if at the current iteration matrix A is not null. Otherwise the absolute
#' error is checked. Then the vector a is checked in the same way.
#'
#' @param matANew Matrix A at the current iteration.
#' @param matAOld Matrix A at the past iteration.
#' @param nvANew Vector a at the current iteration.
#' @param nvAOld Vector a at the past iteration.
#' @param nRelTol Relative tolerance.
#' @param nAbsTol Absolute tolerance.
#'
#' @return A flag indicating if condition on matrix A and vector a are both satisfied.
#'

OBRECheckTolParameters = function (matANew, matAOld, nvANew, nvAOld, nRelTol, nAbsTol) {

# initialize the flag ####
  flagDone = TRUE

# control on vector a ####
  # using the relative error
  if(nvANew[1] != 0 && nvANew[2] != 0) {
    if(abs(nvANew[1] - nvAOld[1]) / nvANew[1] > nRelTol || abs(nvANew[2] - nvAOld[2]) / nvANew[2] > nRelTol) {
      flagDone = FALSE
    }
  # using the absolute error
  } else {
    if(abs(nvANew[1] - nvAOld[1]) > nAbsTol || abs(nvANew[2] - nvAOld[2]) > nAbsTol ) {
      flagDone = FALSE
    }
  }

# control on matrix A ####
  # using the relative error
  if(matANew[1, 1] != 0 && matANew[1, 2] != 0 && matANew[2, 2] != 0) {
    if(abs(matANew[1, 1] - matAOld[1, 1]) / matANew[1, 1] > nRelTol ||
         abs(matANew[1, 2] - matAOld[1, 2]) / matANew[1, 2] > nRelTol ||
         abs(matANew[2, 2] - matAOld[2, 2]) / matANew[2, 2] > nRelTol) {
      flagDone = FALSE
    }
  # using the absolute error
  } else {
    if(abs(matANew[1, 1] - matAOld[1, 1]) > nAbsTol || abs(matANew[1, 2] - matAOld[1, 2]) > nAbsTol ||
         abs(matANew[2, 2] - matAOld[2, 2]) > nAbsTol) {
      flagDone = FALSE
    }
  }
  return(flagDone)
}
