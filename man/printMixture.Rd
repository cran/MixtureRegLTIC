\name{printMixture}
\alias{printMixture}
\title{Print the Summary Table of Regression Results}
\description{Print the summary table of regression results from the output of mixture regression models}
\usage{
printMixture(fit, digits = 3, file = NULL)
}
\arguments{
\item{fit}{the output object from the fitted mixture regression model.}

\item{digits}{an integer indicating the number of decimal places after rounding.}

\item{file}{a character string of output filename for saving the table of regression results.}
}
\seealso{\code{\link{MixtureLogitAFT}}}
\examples{
data(simLTICdataA)

##### fit the logistic-AFT location-scale model for LTIC data
fit=MixtureLogitAFT(formula=Surv(time1,time2,status)~1,
    eventprobreg=~X1,locationreg=~X1,scalereg=~X1,
    var.entry="entry",data=simLTICdataA)

##### print regression results of the fitted model
printMixture(fit)

}
