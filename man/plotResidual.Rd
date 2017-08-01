\name{plotResidual}
\alias{plotResidual}
\title{Estimate and Plot Empirical Distributions of Residuals in the Fitted Regression Model}
\description{
A function to estimate and plot estimated empirical distribution functions of the residuals in the fitted regression model.
}
\usage{
plotResidual(fit, xlab = NULL, ylab = NULL, main = NULL,
col = NULL, lty = NULL, lwd = 1, axes = T)
}
\arguments{

\item{fit}{the output object from the fitted mixture regression model.}

\item{xlab}{the title for x axis.}

\item{ylab}{the title for y axis.}

\item{main}{the main title.}

\item{col}{a vector of colors.}

\item{lty}{a vector of line types.}

\item{lwd}{a numeric value specifies the line width.}

\item{axes}{a logical value specifies whether axes should be drawn. If axes=FALSE, both x and y axes are not shown .}
}
\seealso{\code{\link{MixtureLogitAFT}}}
\references{
Chen CH, Tsay YC, Wu YC and Horng CF. Logistic-AFT location-scale mixture regression models with nonsusceptibility for left-truncated and general interval-censored data. Statistics in Medicine, 2013; 32:4285–4305.
}
\examples{
data(simLTICdataA)

##### fit the logistic-AFT location-scale model for LTIC data
fit=MixtureLogitAFT(formula=Surv(time1,time2,status)~1,
    eventprobreg=~X1,locationreg=~X1,scalereg=~X1,
    var.entry="entry",data=simLTICdataA)

##### plot the empirical distribution of residuals
plot.res=plotResidual(fit)
legend(-9.5,1,legend=plot.res$legend,col=plot.res$col,lty=plot.res$lty,
       title=" Strata (Case / Total)")

}
