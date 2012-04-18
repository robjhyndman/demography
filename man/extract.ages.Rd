\name{extract.ages}
\alias{extract.ages}
\title{Extract some ages from a demogdata object}
\usage{extract.ages(data, ages, combine.upper=TRUE)
}

\arguments{
\item{data}{Demogdata object such as created using \code{\link{read.demogdata}} or \code{\link{smooth.demogdata}}.}
\item{ages}{Vector of ages to extract from data}
\item{combine.upper}{If TRUE, ages beyond the maximum of \code{ages} are combined into the upper age group.}
}

\description{Creates subset of demogdata object.
}

\value{Demogdata object with same components as \code{data} but with a subset of ages.
}

\author{Rob J Hyndman}
\examples{
france.teens <- extract.ages(fr.mort,13:19,FALSE)
plot(france.teens)
}
\keyword{manip}
