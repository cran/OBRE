#' Generic summary method
#'
#' @param object \dots
#' @export
summary <- function(object) UseMethod("summary")

#' Function that summarize the results contained in an OBREresult object.
#'
#' The function shows the estimated parameters, the OBRE tuning parameter,
#' the proportion of data weighted and the relative efficiency with respect to MLE of an OBREresult object.
#'
#' @param object The OBREresult object (output of OBRE function) that has to be plotted.
#'
#' @return The summary an OBREresult obect with the estimated parameters, the OBRE tuning parameter,
#' the proportion of data weighted and the relative efficiency with respect to MLE.
#' @examples
#' \donttest{try({# Generates the Normal distribution input for OBRE
#' distrForOBRE <- densityExpressions(strDistribution = "normal")
#' # Generates input data
#' simData = c(rnorm(100, 12, 1), rnorm(10, 10, 10))
#' # Estimates OBREresult object
#' estOBRE <- OBRE(nvData = simData, strDistribution = distrForOBRE, nCParOBRE = 3)
#' # Summary of the results
#' summary(estOBRE)})}
#'
#' @export

summary.OBREresult = function(object) {
  if(!is(object, "OBREresult")) {
    return(cat("Input is not an OBREresult object"))
  }
 cat(paste(object$strMess, "\n"))
 cat("Estimated parameters", "\n",object$nvTheta, "\n")
 cat("OBRE tuning parameter equal to", object$nCParOBRE, "\n")
 cat("Proportion of data weighted:", length(object$nvWeights[object$nvWeights!=1])/length(object$nvWeights), "\n")
 lOBRECov = OBRECovarianceMatrix(object)
 cat("Relative efficiency with respect to MLE:", lOBRECov$nReiOBRE)
}
