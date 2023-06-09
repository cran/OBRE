#' Numerical Maximum Likelihood Estimator
#'
#' The parameters Maximum Likelihood Estimation is obtained by numerical optimization.
#'
#' @importFrom stats optim var
#' @importFrom methods is
#' @param nvData The vector of the data.
#' @param strDistribution The distribution name.
#' @param lDensityExpr The distribution expression,
#'
#' @return A list with distribution name, distribution parameters, value of the objective function
#' corresponding to the parameters, additional information returned by the optimizer, convergence of
#' the algorithm.
#'


MLE = function (nvData, strDistribution, lDensityExpr) {

  nTheta1 = mean(nvData)
  nTheta2 = sqrt(var(nvData))

# Likelihood parameter optimization ####
  nvTheta = c(nTheta1, nTheta2)
  # optimize using Nelder-Mead
  nlminbOut = try(optim(par = nvTheta, fn = NlLike, method = "Nelder-Mead",
                        control = list(maxit = 1e+06, reltol = 10^-10), nvData = nvData,
                        lDensityExpr = lDensityExpr), silent = TRUE)
  # did not work: optimize with BFGS
  if (is(nlminbOut, "try-error")) {
    nlminbOut1 = try(optim(par = nvTheta, fn = NlLike, method = "BFGS",
                           control = list(maxit = 1e+06, reltol = 10^-12, ndeps = rep(10^-10, 3)),
                           nvData = nvData, lDensityExpr = lDensityExpr), silent = TRUE)
    if (is(nlminbOut1, "try-error")) {
      objective = NlLike(nvTheta = nvTheta, nvData = nvData, lDensityExpr = lDensityExpr)
      message = "NM and BFGS went in error. Initial values returned"
      nlminbOut = list()
      nlminbOut$par = nvTheta
      nlminbOut$value = objective
      nlminbOut$counts = c(-1, -1)
      nlminbOut$convergence = -1
      nlminbOut$message = message
    } else if (!is(nlminbOut1, "try-error")) {
      nlminbOut = nlminbOut1
    }
  } else if (!is(nlminbOut, "try-error")) {
    # NM worked, refine with BFGS if possible
    nlminbOut1 = try(optim(par = nlminbOut$par, fn = NlLike, method = "BFGS",
                           control = list(maxit = 1e+06, reltol = 10^-12, ndeps = rep(10^-10, 3)),
                           nvData = nvData, lDensityExpr = lDensityExpr), silent = TRUE)
    if (!is(nlminbOut1, "try-error")) {
      nlminbOut = nlminbOut1
    }
  }
  nvTheta = c(nlminbOut$par[1 : 2])
  objective = nlminbOut$value
  message = nlminbOut$message

# Output ####
  out = list(strDistribution = strDistribution, nvTheta = nvTheta, objective = objective,
             message = message)
  invisible(out)
}
