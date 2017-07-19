#' Optimal B-Robust Estimator
#'
#' Function for obtaining the Optimal B-Robust Estimates starting by a vector of data and a
#' two parameters distribution.
#'
#' @param nvData The vector of data.
#' @param strDistribution The distribution name between "normal" (Normal distribution), "logNormal" (logNormal distribution),
#' "weibull" (Weibull distribution), "logLogistic" (logLogistic distribution), "gpd2" (Generalized Pareto
#' Distribution with two parameters) or "custom" if the distribution is written by the user as an input of "eDensityFun" parameter.
#' Alternatively, the input of "strDistribution" can be an object of class "OBREdist", obtained using function densityExpressions.
#' @param nCParOBRE OBRE robustness parameter.
#' @param dfParOBRE A data frame containing oprimization parameters, i.e. nEta, the precision
#' between two parameters optimization, nMaxIterLoopWc and nMaxIterLoopA, the number of iterations in the
#' optimization proceture, nRelTol and nAbsTol, the relative and absolute tolerances.
#' @param nTheta1Init First parameter for the beginning of the computation.
#' @param nTheta2Init Second parameter for the beginning of the computation.
#' @param eDensityFun The density of a two parameters distribution. To be inserted if in strDistribution the "custom" option is chosen. This should be an expression object, the
#' two parameters should be called "nTheta1" and "nTheta2", the data "nvData" and its formula should be derivable
#'
#' @return A list with the vector containing the final parameters, the exit OBRE message, the values of vector a and matrix A,
#' the OBRE tuning parameter c, the initial values of the parameters (if unspecified by the user, the values of MLE are reported),
#' the vector of data, the density expression.
#'
#' @export
#' @examples
#' # Using the densityExpressions function for initialize the distribution
#' distrForOBRE <- densityExpressions(strDistribution = "normal")
#' simData = c(rnorm(1000, 12, 2),200,150)
#' \donttest{estOBRE <- OBRE(nvData = simData, strDistribution = distrForOBRE, nCParOBRE = 3)
#' # Launching the generation of the density expression directly from OBRE
#' simData = c(rnorm(1000, 12, 2),200,150)
#' estOBRE <- OBRE(nvData = simData, strDistribution = "normal", nCParOBRE = 3)
#' # Using the "custom" option and using the normal distribution
#' simData = c(rnorm(1000, 12, 2),200,150)
#' estOBRE <- OBRE(nvData = simData, strDistribution = "custom", nCParOBRE = 3,
#' eDensityFun = expression((exp( -((nvData - nTheta1)^2) / (2 * nTheta2^2)) /
#' (sqrt(2 * pi) * nTheta2))))}
#'
#' @references Bellio, R. (2007). Algorithms for bounded-influence estimation. Comput. Stat. Data Anal. 51, 2531-2541.
#' @references Hampel F (1968). Contributions to the theory of robust estimation. University of California.
#' @references Hampel, F., Ronchetti, E., Rousseeuw, P. & Stahel, W. (1985). Robust Statistics. The approach based on influence function. John Wiley and Sons Ltd., Chichester, UK.
#' @references Victoria-Feser, M.P. & Ronchetti, E. (1994). Robust methods for personal-income distribution models. Canadian Journal of Statistics 22, 247-258.
#'

OBRE = function(nvData, strDistribution, nCParOBRE, dfParOBRE = data.frame(nEta = 1e-6,
                                                                           nMaxIterLoopWc = 10,
                                                                           nMaxIterLoopA = 10,
                                                                           nRelTol = 0.001,
                                                                           nAbsTol = 0.5,
                                                                           stringsAsFactors = FALSE),
                nTheta1Init = NA, nTheta2Init = NA, eDensityFun = NA) {

# STEP 1: initialize starting values ####
  # Calculate the A matrix from initial matrix matFisher
  # starting matrix A
  if(class(strDistribution)=="OBREdist"){
    lDensityExpr =  strDistribution
  }else{
    lDensityExpr = densityExpressions(strDistribution = strDistribution, eDensityFun = eDensityFun)
  }
  if (is.na(nTheta1Init) || is.na(nTheta2Init)) {
    lMLE = MLE(nvData = nvData, strDistribution = strDistribution, lDensityExpr = lDensityExpr)
    nTheta1Init = lMLE$nvTheta[1]
    nTheta2Init = lMLE$nvTheta[2]
  }
  matFisherMLE = try(matFisherComputation(nTheta1 = nTheta1Init, nTheta2 = nTheta2Init,
                                          lDensityExpr = lDensityExpr), silent = TRUE)
  invFisher = try(solve(matFisherMLE), silent = TRUE)
  eInvFisher = try(eigen(invFisher), silent = TRUE)
  # check on matrix singularity
  if(class(invFisher) == "try-error") {
    strMess = "Fisher Matrix is singular"
    lOutOBRE = list(nvTheta = c(NA, NA), strMess = strMess, nvA = c(NA, NA),
                    matA = matrix(NA, nrow = 2, ncol = 2), nCParOBRE = nCParOBRE,
                    nvThetaInit = c(nTheta1Init, nTheta2Init),
                    nvWeights = rep(NA, length(nvData)),
                    nvData = nvData,
                    lDensityExpr = lDensityExpr)
    return(lOutOBRE)
  }
  # check on the positiveness of the matrix
  if(min(eInvFisher$values) < 0) {
    strMess = "Fisher Matrix is not positive definite"
    lOutOBRE = list(nvTheta = c(NA, NA), strMess = strMess, nvA = c(NA, NA),
                    matA = matrix(NA, nrow = 2, ncol = 2), nCParOBRE = nCParOBRE,
                    nvThetaInit = c(nTheta1Init, nTheta2Init),
                    nvWeights = rep(NA, length(nvData)),
                    nvData = nvData,
                    lDensityExpr = lDensityExpr)
    return(lOutOBRE)
  }
  # vector a starting value
  nvA = c(0, 0)
  matA = eInvFisher$vectors %*% diag(sqrt(eInvFisher$values)) %*% t(eInvFisher$vectors)

# STEP 2: solve a and A ####
  # Initial value of deltaTheta
  nTheta1 = nTheta1Init
  nTheta2 = nTheta2Init
  nvDeltaTheta = c(nTheta1Init, nTheta2Init)
  nIterLoopWc = 0
  # LOOP B on delta theta
  while(max(abs(nvDeltaTheta / c(nTheta1, nTheta2))) > dfParOBRE$nEta && nIterLoopWc < dfParOBRE$nMaxIterLoopWc) {
    cat("#")
    # Initialise while loop stopping parameters
    nIterLoopA = 0
    flagDone = FALSE
    # loop B: finds vector a and matrix A
    while(!flagDone && nIterLoopA < dfParOBRE$nMaxIterLoopA ) {
      cat(".")
      # Update iteration counter for loop A
      nIterLoopA = nIterLoopA + 1
      # Assign matA and nvA of the previous iteration
      matAOld = matA
      nvAOld = nvA
      nvA = OBREnvAComputation(nvData = nvData, nTheta1 = nTheta1, nTheta2 = nTheta2,
                               lDensityExpr = lDensityExpr, nCParOBRE = nCParOBRE,
                               matA = matA, nvA = nvA)
      # computation of matrix matM2
      matM2 = OBREMatMComputation(nvData = nvData, nTheta1 = nTheta1, nTheta2 = nTheta2,
                                  lDensityExpr = lDensityExpr, nCParOBRE = nCParOBRE,
                                  matA = matA, nvA = nvA, nK = 2)
      invM2 = try(solve(matM2), silent = TRUE)
      # check on the matrix
      if(class(invM2) == "try-error") {
        strMess = "OBRE matrix M2 is singular"
        lOutOBRE = list(nvTheta = c(NA, NA), strMess = strMess, nvA = nvA,
                        matA = matA, nCParOBRE = nCParOBRE,
                        nvThetaInit = c(nTheta1Init, nTheta2Init),
                        nvWeights = rep(NA, length(nvData)),
                        nvData = nvData,
                        lDensityExpr = lDensityExpr)
        return(lOutOBRE)
      }
      eInvM2 = eigen(invM2)
      # check on the positiveness of the matrix
      if(min(eInvM2$values) < 0) {
        strMess = "OBRE matrix M2 is not positive definite"
        lOutOBRE = list(nvTheta = c(NA, NA), strMess = strMess, nvA = nvA,
                        matA = matA, nCParOBRE = nCParOBRE,
                        nvThetaInit = c(nTheta1Init, nTheta2Init),
                        nvWeights = rep(NA, length(nvData)),
                        nvData = nvData,
                        lDensityExpr = lDensityExpr)
        return(lOutOBRE)
      }
      # Use Cholesky decomposition to find matrix A
      matA = chol(invM2)
      # check that values of vector a and matrix A don't exceed given tolerance
      flagDone = OBRECheckTolParameters(matANew = matA, matAOld = matAOld, nvANew = nvA, nvAOld = nvAOld,
                                        nRelTol = dfParOBRE$nRelTol, nAbsTol = dfParOBRE$nAbsTol)
    }

# STEP 3: compute delta theta and update theta ####
    # calculate matrix M1
    matM1 = OBREMatMComputation(nvData = nvData, nTheta1 = nTheta1, nTheta2 = nTheta2,
                                lDensityExpr = lDensityExpr, nCParOBRE = nCParOBRE,
                                matA = matA, nvA = nvA, nK = 1)
    # calculate deltaTheta
    nvWeightFunLoop = OBREWeightsFun(nvData = nvData, nTheta1 = nTheta1, nTheta2 = nTheta2,
                                     lDensityExpr = lDensityExpr, nCParOBRE = nCParOBRE,
                                     matA = matA, nvA = nvA)
    nMeanComponent1 = mean((scoreComponent(nvData = nvData, nTheta1 = nTheta1, nTheta2 = nTheta2,
                                           lDensityExpr = lDensityExpr, nParIndex = 1) - nvA[1]) * nvWeightFunLoop)
    nMeanComponent2 = mean((scoreComponent(nvData = nvData, nTheta1 = nTheta1, nTheta2 = nTheta2,
                                           lDensityExpr = lDensityExpr, nParIndex = 2) - nvA[2]) * nvWeightFunLoop)
    nvDeltaTheta = solve(matM1) %*% c(nMeanComponent1, nMeanComponent2)

    # Update theta1 and theta2
    nTheta1 = nTheta1 + nvDeltaTheta[1]
    nTheta2 = nTheta2 + nvDeltaTheta[2]
    nIterLoopWc = nIterLoopWc + 1
  }
  strMess = "OBRE converged in less than maxiterLoopWc iterations"
  if (nIterLoopWc >= dfParOBRE$nMaxIterLoopWc) {
    strMess = "OBRE terminates with maximum number of iteration"
  }

# Output ####
  lOutOBRE = list(nvTheta = c(nTheta1, nTheta2), strMess = strMess, nvA = nvA,
                  matA = matA, nCParOBRE = nCParOBRE,
                  nvThetaInit = c(nTheta1Init, nTheta2Init),
                  nvWeights = nvWeightFunLoop,
                  nvData = nvData,
                  lDensityExpr = lDensityExpr)
  class(lOutOBRE) = "OBREresult"
  return(lOutOBRE)
}
