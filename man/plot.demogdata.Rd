\name{plot.demogdata}
\alias{plot.demogdata}
\alias{lines.demogdata}
\alias{points.demogdata}
\title{Plot age-specific demographic functions}
\usage{\method{plot}{demogdata}(x, series=ifelse(!is.null(x$rate),names(x$rate)[1],names(x$pop)[1]), 
    datatype=ifelse(!is.null(x$rate),"rate","pop"),
    years=x$year, ages=x$age, max.age=max(x$age), 
    transform=(x$type=="mortality"), 
    plot.type=c("functions","time","depth","density"), 
    type="l", main=NULL, xlab, ylab,...)
\method{lines}{demogdata}(x, series=ifelse(!is.null(x$rate),names(x$rate)[1],names(x$pop)[1]), 
    datatype=ifelse(!is.null(x$rate),"rate",""),
    years=x$year, ages=x$age, max.age=max(x$age), 
    transform=(x$type=="mortality"), 
    plot.type=c("functions","time","depth","density"), ...)
\method{points}{demogdata}(..., pch = 1)
}


\arguments{
\item{x}{Demogdata object such as created using \code{\link{read.demogdata}} or \code{\link{smooth.demogdata}}.}
\item{series}{Name of series to plot. Default: the first matrix within \code{datatype}.}
\item{datatype}{Name of demogdata object which contains series. Default \dQuote{rate}. Alternative: \dQuote{pop}.}
\item{years}{Vector indicating which years to plot. Default: all available years.}
\item{ages}{Vector indicating which ages to plot. Default: all available ages.}
\item{max.age}{Maximum age to plot. Default: all available ages.}
\item{transform}{Should a transformation of the data be plotted? Default is TRUE if the object contains mortality data and datatype="rate", and FALSE otherwise.}
\item{plot.type}{Type of plot: either \dQuote{functions} or \dQuote{time}.}
\item{type}{What type of plot should be drawn. See \code{\link[graphics]{plot}} for possible types.}
\item{main}{Main title for the plot.}
\item{xlab}{Label for x-axis.}
\item{ylab}{Label for y-axis.}
\item{pch}{Plotting character.}
\item{...}{Other plotting parameters. In \code{points.demogdata}, all arguments are passed to \code{lines.demogdata.}}
}

\description{If \code{plot.type="functions"}, then years are plotted using a rainbow palette so the
earliest years are red, followed by orange, yellow, green, blue
and indigo with the most recent years plotted in violet.

If \code{plot.type="time"}, then each age is shown as a separate time series in a time plot.}

\value{None. Function produces a plot}

\author{Rob J Hyndman}
\examples{
plot(fr.mort)
par(mfrow=c(1,2))
plot(aus.fert,plot.type="time")
plot(aus.fert,plot.type="functions")
}
\keyword{hplot}
