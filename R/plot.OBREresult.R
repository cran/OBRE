
#'
#'
#' Function that plot an OBREresult object.
#'
#' The function computes the plot of the OBRE computation
#'
#' @importFrom methods is
#' @importFrom graphics par
#' @param x The OBREresult object (output of OBRE function) that has to be plotted.
#' @param ... Added argument for consistency with the plot generic function.
#'
#' @return A graphical representation of an OBREresult obect. The plot is composed by four plots: the value of input data in logaritmic scale, the values of score function evaluated in the input data, the OBRE weights, the values of OBRE components.
#' @examples
#' \donttest{try({# Generates the Normal distribution input for OBRE
#' distrForOBRE <- densityExpressions(strDistribution = "normal")
#' # Generates input data
#' simData = c(rnorm(100, 12, 1), rnorm(10, 10, 10))
#' # Estimates OBREresult object
#' estOBRE = OBRE(nvData = simData, strDistribution = "normal", nCParOBRE = 3)
#' plot(estOBRE)})}
#'
#' @export
plot.OBREresult = function(x, ...) {
  if (!is(x, "OBREresult")) {
    cat("The input is not a OBREresult object.")
    return("Input error.")
  }
  par(mfrow = c(2, 2))
  plot(x$nvData, main = "Input data", ylab = "x")
  nvScore1 = scoreComponent(nvData = x$nvData, nTheta1 = x$nvTheta[1], nTheta2 = x$nvTheta[2],
                         lDensityExpr = x$lDensityExpr, nParIndex = 1)
  nvScore2 = scoreComponent(nvData = x$nvData, nTheta1 = x$nvTheta[1], nTheta2 = x$nvTheta[2],
                         lDensityExpr = x$lDensityExpr, nParIndex = 2)
  plot(x = nvScore1, y = nvScore2, xlab = "First parameter",
       ylab = "Second parameter", main = "Values of the score function")
  plot(x$nvWeights, main = "OBRE Weights", ylab = "Weights")
  plot(x = (nvScore1 - x$nvA[1]) * x$matA[1,1] * x$nvWeights,
       y = (nvScore2 - x$nvA[2]) * x$matA[2,2] * x$nvWeights,
       xlab = "First parameter", ylab = "Second parameter",
       main = "Values of OBRE components", xlim = c(-3, 3),
       ylim = c(-2,4))
  par(mfrow = c(1, 1))
}
