\name{models}
\alias{models}
\alias{models.fmforecast}
\alias{models.fmforecast2}
\title{Show model information for the forecast coefficients in FDM models.}
\usage{
\method{models}{fmforecast}(object, select=0, ...)
\method{models}{fmforecast2}(object, ...)
}

\arguments{
\item{object}{Output from \code{\link{forecast.fdm}} or \code{\link{forecast.fdmpr}}.}
\item{select}{Indexes of coefficients to display. If select=0, all coefficients are displayed.}
\item{...}{Other arguments.}
}

\description{The models for the time series coefficients used in forecasting fdm models are shown.
}

\value{Nothing is returned.}


\author{Rob J Hyndman.}

\seealso{\code{\link{forecast.fdm}}, \code{\link{forecast.fdmpr}}.}

\examples{
\dontrun{
fr.short <- extract.years(fr.sm,1950:2006)
fr.fit <- fdm(fr.short,series="male")
fr.fcast <- forecast(fr.fit)
models(fr.fcast)

fr.fit <- coherentfdm(fr.short)
fr.fcast <- forecast(fr.fit)
models(fr.fcast,select=1:3)
}}
\keyword{models}
