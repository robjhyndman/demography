\name{coherentfdm}
\alias{coherentfdm}
\title{Coherent functional demographic model for grouped data}
\usage{coherentfdm(data, order1=6, order2=6,  ...)
}

\arguments{
\item{data}{demogdata object containing at least two groups.}
\item{order1}{Number of basis functions to fit to the model for the geometric mean.}
\item{order2}{Number of basis functions to fit to the models for each ratio.}
\item{...}{Extra arguments passed to \code{\link{fdm}}.}
}

\description{Fits a coherent functional model to demographic data as described in Hyndman, Booth & Yasmeen (2012). If two of the series
in \code{data} are named \code{male} and \code{female}, then it will use these two groups. Otherwise it will use all available groups.
}

\value{A list (of class \code{fdmpr}) consisting of two objects: \code{product} (an \code{\link{fdm}} object containing a 
model for the geometric mean of the data) and \code{ratio} (a list of \code{\link{fdm}} objects, being the models for the 
ratio of each series with the geometric mean).}


\author{Rob J Hyndman.}

\references{Hyndman, R.J., Booth, H., and Yasmeen, F. (2012) Coherent mortality forecasting: the product-ratio method with functional time series models.
\emph{Demography}, to appear. \url{http://robjhyndman.com/papers/coherentfdm}
}

\seealso{\code{\link{fdm}}, \code{\link{forecast.fdmpr}}}

\examples{
fr.short <- extract.years(fr.sm,1950:2006)
fr.fit <- coherentfdm(fr.short)
summary(fr.fit)
plot(fr.fit$product, components=3)
}
\keyword{models}
