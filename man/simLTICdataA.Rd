\name{simLTICdataA}
\alias{simLTICdataA}
\docType{data}
\title{
A Generated Set with LTIC Data in Simulation Scenario (a) of Chen et al. (2013)
}
\description{
A generated set with left-truncated and interval-censored data in simulation scenario (a) of Chen et al. (2013) can be used in the package MixtureRegLTIC for fitting mixture regression models.}
\usage{data(simLTICdataA)}
\format{
There are 500 simulated observations with the following variables:
\describe{
\item{\code{time1}}{left time-endpoint of interval censored data}
\item{\code{time2}}{right time-endpoint of interval censored data}
\item{\code{status}}{status indicator with 0=right censored, 1=exact event, 2=left censored, and 3=interval censored.}
\item{\code{entry}}{time at study entry}
\item{\code{X1}}{a binary covariate}
}
}
\seealso{\code{\link{NPMLEsurv}}, \code{\link{MixtureLogitAFT}}}
\references{
Chen CH, Tsay YC, Wu YC and Horng CF. Logistic-AFT location-scale mixture regression models with nonsusceptibility for left-truncated and general interval-censored data. Statistics in Medicine, 2013; 32:4285–4305.
}
\keyword{datasets}
