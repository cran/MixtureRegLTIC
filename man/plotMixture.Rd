\name{plotMmixture}
\alias{plotMixture}
\title{Plot Survival/Event Curve(s) from the Fitted Regression Model}
\description{
A function to plot estimated overall and conditional survival/event curve(s) from the fitted regression models.
}
\usage{
plotMixture(fit, dist = "overall", curve = "survival",
xlab = NULL, ylab = NULL, main = NULL, col = NULL, lty = NULL, lwd = 1, axes = T)
}
\arguments{
\item{fit}{the output object from the fitted mixture regression model.}

\item{dist}{a character string specified with either "overall" or "cond". The default dist="overall" plots the overall distribution of the event time, and dist="cond" plots the conditional distribution for the logarithm of the event time for susceptible/non-cured subjects.}

\item{curve}{a character string specifies the type of desired curves to be plotted. The default curve="survival" plots the survival curves, and curve="event" plots event curves.}

\item{xlab}{the title for x axis.}

\item{ylab}{the title for y axis.}

\item{main}{the main title of the plot.}

\item{col}{a vector of colors.}

\item{lty}{a vector of line types.}

\item{lwd}{a numeric value specifies the line width.}

\item{axes}{a logical value specifies whether axes should be drawn. If axes=FALSE, both x and y axes are not shown .}
}
\seealso{\code{\link{MixtureLogitAFT}}, \code{\link{plotNPMLEsurv}}}
\examples{
data(simLTICdataE)

##### fit the logistic-AFT location-scale model for LTIC data
fit=MixtureLogitAFT(formula=Surv(time1,time2,status)~1,
    eventprobreg=~X1,locationreg=~X1+X2,scalereg=~X1+X2,
    var.entry="entry",var.mixturetype="mtype",data=simLTICdataE)

##### print regression results of the fitted model
printMixture(fit)

##### plot estimated event curves of the fitted model
#win.graph(width=18,height=10)
#par(mfrow=c(1,2))
plot.fit=plotMixture(fit,curve="event",col=c("red","blue","deeppink"))
legend(55,0.95,legend=plot.fit$legend,col=plot.fit$col,lty=plot.fit$lty,
       title=" Strata (Case / Total)",cex=0.85)

plotD.fit=plotMixture(fit,dist="cond",curve="event",col=c("red","blue","deeppink"))
legend(3,0.95,legend=plotD.fit$legend,col=plotD.fit$col,lty=plotD.fit$lty,
       title=" Strata (Case / Total)",cex=0.85)

}
