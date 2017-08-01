\name{NPMLEsurv}
\alias{NPMLEsurv}
\title{Estimate Overall and Conditional Survival/Event Curve(s) Based on the NPMLE for LTIC Data}
\description{
A function to estimate the nonparametric maximum likelihood estimator (NPMLE) of the distribution function for LTIC data based on the Turnbull (1976) method corrected by Frydman (1994).
}
\usage{
NPMLEsurv(formula, var.entry, var.weight = NULL, data, time.origin = 0)
}
\arguments{

\item{formula}{A formula object specifies the regression response on the left of a ~ operator using a survival object of the same form as from the Surv function, and fitted covariates on the right. The status indicator is 0=right censored, 1=event at time, 2=left censored, and 3=interval censored.}

\item{var.entry}{A variable specifies each study subject's left-truncated time at the entry in the follow-up study. The var.entry = NULL corresponds to no left-truncation.}

\item{var.weight}{A numeric variable specifies the weight for each observation contributed differently to the log-likelihood. The default var.weight = NULL corresponds to equal contributions.}

\item{data}{A data.frame contains the variables named in the formula, var.entry and var.weight.}

\item{time.origin}{A numeric value specifies the time origin. The default time.origin = 0.}
}
\references{
Frydman H. A note on nonparametric estimation of the distribution function from interval-censored and truncated observations. Journal of the Royal Statistical Society, Series B, 1994; 56:71–74.

Turnbull BW. The empirical distribution function with arbitrarily grouped, censored and truncated data. Journal of the Royal Statistical Society, Series B, 1976; 38:290–295.
}
\seealso{\code{\link{MixtureLogitAFT}}, \code{\link{plotNPMLEsurv}}}
\examples{
data(simLTICdataA)

##### estimate the NPMLE
est=NPMLEsurv(formula=Surv(time1,time2,status)~X1,var.entry="entry",data=simLTICdataA)

### plot estimated survival curves with NPMLE
#win.graph(width=18,height=10)
#par(mfrow=c(1,2))
plot.NPMLE=plotNPMLEsurv(est,lty=c(2,2),col=c("red","blue"))
legend(20,0.4,legend=plot.NPMLE$legend,col=plot.NPMLE$col,lty=plot.NPMLE$lty,
       title=" Strata (Case / Total)")

plotD.NPMLE=plotNPMLEsurv(est,dist="cond",lty=c(2,2),col=c("red","blue"))
legend(3,0.4,legend=plotD.NPMLE$legend,col=plotD.NPMLE$col,lty=plotD.NPMLE$lty,
       title=" Strata (Case / Total)")

}

