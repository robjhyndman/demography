\name{set.upperage}
\alias{set.upperage}
\title{Combine the upperages of a demogdata object.}
\usage{set.upperage(data, max.age)
}

\arguments{
\item{data}{Demogdata object such as created using \code{\link{read.demogdata}} or \code{\link{smooth.demogdata}}.}
\item{max.age}{Upper age group. Ages beyond this are combined into the upper age group.}
}

\description{Computes demographic rates by combining age groups.
}

\value{Demogdata object with same components as \code{data} but with a subset of ages.
}

\author{Rob J Hyndman}
\examples{
france.short <- set.upperage(fr.mort, 85)
}
\keyword{manip}
