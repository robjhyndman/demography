\name{mean.demogdata}
\alias{mean.demogdata}
\alias{median.demogdata}
\title{Mean and median functions for data of class demogdata}
\usage{\method{mean}{demogdata}(x, series = names(x$rate)[1], transform = TRUE, na.rm = TRUE, ...)
\method{median}{demogdata}(x, na.rm=FALSE, series = names(x$rate)[1], transform = TRUE,
    method = c("hossjercroux", "coordinate"), ...)
}

\arguments{
\item{x}{Demogdata object such as created using \code{\link{read.demogdata}} or \code{\link{smooth.demogdata}}.}
\item{series}{Name of demogdata series to plot..}
\item{transform}{Should transform of data be taken first?}
\item{na.rm}{a logical value indicating whether NA values should be stripped before the computation proceeds.}
\item{method}{Method for computing the median. Either "coordinate" for a coordinate-wise median, or "hossjercroux" for the
L1-median using the Hossjer-Croux algorithm.}
\item{...}{Other arguments.}
}

\description{Computes mean or median of demographic rates for each age level.}

\value{A list containing \code{x}=ages and \code{y}=mean or median rates.}

\author{Rob J Hyndman}

\references{Hossjer, O., and Croux, C. (1995) Generalized univariate signed rank statistics for testing
and estimating a multivariate location parameter. \emph{Nonparametric Statistics}, \bold{4}, 293-308.
}

\examples{
plot(fr.mort)
lines(mean(fr.mort),lwd=2)
lines(median(fr.mort),lwd=2,col=2)
}
\keyword{models}
