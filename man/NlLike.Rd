% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NlLike.R
\name{NlLike}
\alias{NlLike}
\title{Negative Log-Likelihood}
\usage{
NlLike(nvTheta, nvData, lDensityExpr)
}
\arguments{
\item{nvTheta}{Parameters of the distribution.}

\item{nvData}{The vector of the data.}

\item{lDensityExpr}{The distribution density espressions.}
}
\value{
Negative log likelihood value.
}
\description{
The function compute the Negative Log-Likelihood value that has to be used for optimization in MLE function.
}
