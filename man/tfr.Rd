\name{tfr}
\alias{tfr}
\title{Compute total fertility rate from fertility rates}
\usage{
tfr(data, PI=FALSE, nsim=500, ...)
}

\arguments{
\item{data}{Demogdata object of type \code{"fertility"} such as obtained from \code{\link{read.demogdata}}, \code{\link{forecast.fdm}}}.
\item{PI}{If TRUE, produce a prediction interval.}
\item{nsim}{Number of simulations to use when computing a prediction interval.}
\item{...}{Other arguments passed to \code{simulate} when producing prediction intervals.}
}

\value{If data are of class \code{demogdata}, the function returns a time series of fertility rates.
If data are from \code{\link{forecast.fdm}}, the function returns an object of class \code{forecast} containing
point forecasts and (optionally) prediction intervals.
}

\description{Compute total fertility rates from age-specific fertility rates contained in a \code{demogdata} object.}


\seealso{\code{\link{fdm}}}

\author{Rob J Hyndman}


\examples{
plot(tfr(aus.fert))
ausfert.fcast <- forecast(fdm(aus.fert))
plot(tfr(ausfert.fcast,PI=TRUE,nsim=400))
}
\keyword{models}
