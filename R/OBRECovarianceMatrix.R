#' Function that computes the OBRE covariance matrix.
#'
#' The function computes matrices M (Jacobian) and Q (Variability) and
#' uses them to evaluate the covariance matrix V.
#'
#' @import pracma
#' @param lOBRE List of all the variables resulting from the OBRE computation.
#'
#' @export
#' @return A list containing Jacobian of the estimate function, variability and asymptotic covariance
#' matrices, as well as the relative efficiency with respect to Maximum Likelihood Estimator
#'
#' @examples
#' \donttest{distrForOBRE <- densityExpressions(strDistribution = "normal")
#' simData = c(rnorm(1000, 12, 2),200,150)
#' estOBRE <- OBRE(nvData = simData, strDistribution = distrForOBRE, nCParOBRE = 3)
#' lOBRECov = OBRECovarianceMatrix(estOBRE)}
#'
#' @references Hampel, F., Ronchetti, E., Rousseeuw, P. & Stahel, W. (1985). Robust Statistics. The approach based on influence function. John Wiley and Sons Ltd., Chichester, UK.
#' @references Heritier S, Cantoni E, Copt S, Victoria-Feser M (2011). Robust Methods in Biostatistics. John Wiley and Sons Ltd., Chichester, UK.
#'


OBRECovarianceMatrix = function(lOBRE) {

#Inizialisation of the variables ####
  nTheta1 = lOBRE$nvTheta[1]
  nTheta2 = lOBRE$nvTheta[2]
  lDensityExpr = lOBRE$lDensityExpr
  nCParOBRE = lOBRE$nCParOBRE
  matA = lOBRE$matA
  nvA = lOBRE$nvA

# Computation of Jacobian of the estimate function matrix (M) ####
  nMatMEl11 = quadinf(f = OBREMatVMatMEl11, xa = 0, xb = Inf, tol = .Machine$double.eps,
                      nTheta1 = nTheta1, nTheta2 = nTheta2, lDensityExpr = lDensityExpr,
                      nCParOBRE = nCParOBRE, matA = matA, nvA = nvA)$Q
  nMatMEl12 = quadinf(f = OBREMatVMatMEl12, xa = 0, xb = Inf, tol = .Machine$double.eps,
                      nTheta1 = nTheta1, nTheta2 = nTheta2, lDensityExpr = lDensityExpr,
                      nCParOBRE = nCParOBRE, matA = matA, nvA = nvA)$Q
  nMatMEl21 = quadinf(f = OBREMatVMatMEl21, xa = 0, xb = Inf, tol = .Machine$double.eps,
                      nTheta1 = nTheta1, nTheta2 = nTheta2, lDensityExpr = lDensityExpr,
                      nCParOBRE = nCParOBRE, matA = matA, nvA = nvA)$Q
  nMatMEl22 = quadinf(f = OBREMatVMatMEl22, xa = 0, xb = Inf, tol = .Machine$double.eps,
                      nTheta1 = nTheta1, nTheta2 = nTheta2, lDensityExpr = lDensityExpr,
                      nCParOBRE = nCParOBRE, matA = matA, nvA = nvA)$Q
  matJacobianEstFun = matrix(data = c(nMatMEl11, nMatMEl12, nMatMEl21, nMatMEl22), nrow = 2, ncol = 2)

# Computation of variability matrix (Q) ####
  nMatQEl11 = quadinf(f = OBREMatVMatQEl11, xa = 0, xb = Inf, tol = .Machine$double.eps,
                      nTheta1 = nTheta1, nTheta2 = nTheta2, lDensityExpr = lDensityExpr,
                      nCParOBRE = nCParOBRE, matA = matA, nvA = nvA)$Q
  nMatQEl12 = quadinf(f = OBREMatVMatQEl12, xa = 0, xb = Inf, tol = .Machine$double.eps,
                      nTheta1 = nTheta1, nTheta2 = nTheta2, lDensityExpr = lDensityExpr,
                      nCParOBRE = nCParOBRE, matA = matA, nvA = nvA)$Q
  nMatQEl22 = quadinf(f = OBREMatVMatQEl22, xa = 0, xb = Inf, tol = .Machine$double.eps,
                      nTheta1 = nTheta1, nTheta2 = nTheta2, lDensityExpr = lDensityExpr,
                      nCParOBRE = nCParOBRE, matA = matA, nvA = nvA)$Q
  matVariability = matrix(data = c(nMatQEl11, nMatQEl12, nMatQEl12, nMatQEl22), nrow = 2, ncol = 2)

# Computation of the asymptotic covariance matrix (V) ####
  matAsymptCovariance = solve(matJacobianEstFun) %*% matVariability %*% solve(t(matJacobianEstFun))

# Computation of the OBRE efficiency indicator
  matFisherOBRE = try(matFisherComputation(nTheta1 = nTheta1, nTheta2 = nTheta2,
                                           lDensityExpr = lDensityExpr), silent = TRUE)
  invFisherOBRE = try(solve(matFisherOBRE), silent = TRUE)
  if(class(invFisherOBRE) == "try-error") {
    nReiOBRE = NA
  } else {
    nArgumentREIOBRE = det(invFisherOBRE) / det(matAsymptCovariance)
    if (nArgumentREIOBRE < 0) {
      nReiOBRE = 0
    } else {
      nReiOBRE = sqrt(nArgumentREIOBRE)
    }
  }

# Generating output ####
  lOut = list(matJacobianEstFun = matJacobianEstFun, matVariability = matVariability,
              matAsymptCovariance = matAsymptCovariance, nReiOBRE = nReiOBRE)
  return(lOut)
}
