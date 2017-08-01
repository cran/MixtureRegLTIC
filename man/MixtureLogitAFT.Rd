\name{MixtureLogitAFT}
\alias{MixtureLogitAFT}
\title{Fit Logistic-AFT Location-Scale Mixture Regression Models}
\description{
A function to fit parametric logistic-AFT location-scale mixture regression models with nonsusceptibility for LTIC data.
}
\usage{
MixtureLogitAFT(formula, eventprobreg = ~1, locationreg = ~1, scalereg = ~1,
var.entry, var.mixturetype = NULL, var.weight = NULL,
data, time.origin = 0, shape = NULL)
}
\arguments{
\item{formula}{A formula object specifies the regression response using a survival object of the same form as from the Surv function. The status indicator is 0=right censored, 1=event at time, 2=left censored, and 3=interval censored.}

\item{eventprobreg}{The formula for covariates of the logistic regression part for the event (susceptibility/non-cure) probability in the mixture model. The eventprobreg = NULL corresponds to the one-component model.}

\item{locationreg}{The formula for covariates of the location regression part in the AFT submodel.}

\item{scalereg}{The formula for covariates of the scale regression part in the AFT submodel.}

\item{var.entry}{A variable specifies each study subject's left-truncated time at the entry in the follow-up study. The var.entry = NULL corresponds to no left-truncation.}

\item{var.mixturetype}{A variable specifies the number of fitted component(s) in the mixture regression model for each study subject. The defined value of 1 is for "1-component", and the value of 2 is for “2-component with cure”. All subjects in the same stratum must have the same defined values. If all strata have either one-component or two-component models, then var.mixturetype needs not to be included in the data.frame by setting var.mixturetype = NULL.}

\item{var.weight}{A numeric variable specifies the weight for each observation contributed differently to the log-likelihood. The default var.weight = NULL corresponds to equal contributions.}

\item{data}{A data.frame contains the variables named in the formula, eventprobreg, locationreg, scalereg, var.entry, var.mixturetype and var.weight.}

\item{time.origin}{A numeric value specifies the time origin in the AFT submodel. The default time.origin = 0.}

\item{shape}{
The default shape = NULL indicates the shape parameter of the generalized log-gamma distribution to be estimated from the data. Otherwise,  a numeric value assigns the fixed shape parameter of a specified generalized log-gamma distribution of the error variate to be fitted in the AFT submodel, or a character string "logistic" stands for the logistic distribution of the error variate to be fitted in the AFT submodel.}
}
\details{
This function fits the logistic-AFT location-scale mixture regression models with nonsusceptibility/cure for LTIC data. The event time is assumed following either a generalized gamma distribution or a log-logistic distribution, i.e., the error variate in the AFT submodel has a generalized log-gamma distribution or a logistic distribution. The family of generalized gamma distributions for the event time includes many important distributions as its special cases: Weibull (shape = 1), lognormal (shape = 0), and reciprocal Weibull (shape = -1) distributions. Ordinarily, the shape parameter of the generalized gamma distribution is estimated jointly with the logistic, location and scale regression parameters. If a special error distribution is preferred, then the argument shape should be specified.

To obtain the maximum likelihood estimators of the parameters in the logistic-AFT location-scale mixture regression model, R function optim() is applied using the negative log-likelihood function and its first and second derivatives.

As indicated in the scenario (e) of Figure 1 in Chen et al. (2013), some strata may present with nonsusceptibility/cure and the other strata without that. This function can simultaneously handle situations in which both one-component and two-component location-scale regression models jointly emerge as in the above-mentioned scenario. The usage of var.mixturetype facilitates this kind of analysis. However, a consideration of the joint one-component and two-component mixture regression models needs special caution in selecting the covariates and interpreting the corresponding parameter estimates in the logistic regression part.
}
\note{To obtain the standard errors and confidence intervals of the estimated regression parameters, formulas of the associated first and second derivatives were given in Appendix A of Chen et al. (2013) by the chain rule and the Leibnitz rule. An unexported R function HGNlogL() calculates the negative log-likelihood and the corresponding first and second derivatives with several unexported R functions. In addition, a FORTRAN program with algorithm and codes in Moore (1982) were used to compute the incomplete gamma integral.}

\references{
Chen CH, Tsay YC, Wu YC and Horng CF. Logistic-AFT location-scale mixture regression models with nonsusceptibility for left-truncated and general interval-censored data. Statistics in Medicine, 2013; 32:4285–4305.

Moore RJ. Algorithm AS. 187: derivatives of the incomplete gamma integral. Applied Statistics - Journal of the Royal Statistical Society, Series C, 1982; 31:330–333.
}
\seealso{\code{\link{NPMLEsurv}}, \code{\link{plotMixture}}, \code{\link{plotNPMLEsurv}}, \code{\link{plotResidual}}, \code{\link{printMixture}}}
\examples{
data(simLTICdataA)

##### fit the logistic-AFT location-scale model for LTIC data
fit=MixtureLogitAFT(formula=Surv(time1,time2,status)~1,
    eventprobreg=~X1,locationreg=~X1,scalereg=~X1,
    var.entry="entry",data=simLTICdataA)

##### print regression results of the fitted regression model
printMixture(fit)

##### plot estimated survival curves
#win.graph(width=18,height=10)
#par(mfrow=c(1,2))
plot.fit=plotMixture(fit)
legend(20,0.4,legend=plot.fit$legend,col=plot.fit$col,lty=plot.fit$lty,
       title=" Strata (Case / Total)")

plotD.fit=plotMixture(fit,dist="cond")
legend(3,0.4,legend=plotD.fit$legend,col=plotD.fit$col,lty=plotD.fit$lty,
       title=" Strata (Case / Total)")

##### estimate the NPMLE
est=NPMLEsurv(formula=Surv(time1,time2,status)~X1,var.entry="entry",data=simLTICdataA)

##### plot estimated event curves with both the regression model and NPMLE
#win.graph(width=18,height=10)
#par(mfrow=c(1,2))
plot.fit=plotMixture(fit,curve="event",col=c("red","blue"))
legend(20,1,legend=plot.fit$legend,col=plot.fit$col,lty=plot.fit$lty,
       title=" Strata (Case / Total)")
par(new=TRUE)
plot.NPMLE=plotNPMLEsurv(est,curve="event",lty=c(2,2),col=c("red","blue"))

plotD.fit=plotMixture(fit,curve="event",dist="cond",col=c("red","blue"))
legend(3,1,legend=plotD.fit$legend,col=plotD.fit$col,lty=plotD.fit$lty,
       title=" Strata (Case / Total)")
par(new=TRUE)
plotD.NPMLE=plotNPMLEsurv(est,dist="cond",curve="event",lty=c(2,2),col=c("red","blue"))

}
