\name{extract.years}
\alias{extract.years}
\title{Extract some years from a demogdata object}
\usage{extract.years(data, years)
}

\arguments{
\item{data}{Demogdata object such as created using \code{\link{read.demogdata}} or \code{\link{smooth.demogdata}}.}
\item{years}{Vector of years to extract from data}
}

\description{Creates subset of demogdata object.
}

\value{Demogdata object with same components as \code{data} but with a subset of years.
}

\author{Rob J Hyndman}
\examples{
france.1918 <- extract.years(fr.mort,1918)
}
\keyword{manip}
