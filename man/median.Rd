\name{median}
\alias{median}
\alias{median.default}
\title{Median}
\usage{median(x, na.rm=FALSE, ...)
}

\arguments{
\item{x}{an object for which a method has been defined, or a numeric vector containing the values whose median is to be computed.}
\item{na.rm}{a logical value indicating whether NA values should be stripped before the computation proceeds.}
\item{...}{Arguments passed to specific methods.}
}

\description{Generic function for the median.}

\details{This is necessary because the \code{\link[stats]{median}} function in the \code{stats} package does not allow additional parameters.}

\value{Refer to specific methods. For numeric vectors, see the \code{\link[stats]{median}} function in the \code{stats} package.}

\author{Rob J Hyndman}

\seealso{\code{\link{median.demogdata}}}

\keyword{models}
