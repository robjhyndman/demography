\name{sex.ratio}
\alias{sex.ratio}
\title{Compute sex ratios from mortality rates}
\usage{
sex.ratio(data)
}

\arguments{
\item{data}{Demogdata object of type \dQuote{mortality} such as obtained from \code{\link{read.demogdata}}, 
or an object of class \code{fmforecast} such as the output from  \code{\link{forecast.fdm}} or \code{\link{forecast.lca}}.}
}

\description{Calculates the Male/Female ratios from historical or forecasted mortality rates.}

\value{Functional time series of sex ratios.}

\author{Rob J Hyndman}
\examples{
plot(sex.ratio(fr.mort),ylab="Sex ratios (M/F)")
}
\keyword{models}
