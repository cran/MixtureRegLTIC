\name{plotNPMLEsurv}
\alias{plotNPMLEsurv}
\title{Plot Survival/Event Curve(s) from the NPMLE}
\description{
A function to plot overall and conditional survival/event curve(s) from nonparametric maximum likelihood estimators (NPMLEs) of distribution functions.
}
\usage{
plotNPMLEsurv(est, dist = "overall", curve = "survival", type = "s",
xlab = NULL, ylab = NULL, main = NULL, col = NULL, lty = NULL, lwd = 1, axes = T)
}
\arguments{
\item{est}{the output object from the estimated NPMLE}

\item{dist}{a character string specified with either "overall" or "cond". The default dist="overall" plots the overall distribution of the event time, and dist="cond" plots the conditional distribution for the logarithm of the event time for susceptible/non-cured subjects.}

\item{curve}{a character string specifies the type of desired curves to be plotted. The default curve="survival" plots the survival curves, and. curve="event" plots event curves.}

\item{type}{a character string specifies the type of plot, as in plot.}

\item{xlab}{the title for x axis.}

\item{ylab}{the title for y axis.}

\item{main}{the main title of the plot.}

\item{col}{a vector of colors.}

\item{lty}{a vector of line types.}

\item{lwd}{a numeric value specifies the line width.}

\item{axes}{a logical value specifies whether axes should be drawn. If axes=FALSE, both x and y axes are not shown .}
}
\seealso{\code{\link{NPMLEsurv}}, \code{\link{plotMixture}}}
\examples{
data(simLTICdataA)

##### estimate the NPMLE #####
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
