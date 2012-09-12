\name{forecast.fdmpr}
\alias{forecast.fdmpr}
\title{Forecast coherent functional demographic model.}
\usage{
\method{forecast}{fdmpr}(object, h=50, level=80, K=100, drange=c(0.0,0.5), ...)
}

\arguments{
\item{object}{Output from \code{\link{coherentfdm}}.}
\item{h}{Forecast horizon.}
\item{level}{Confidence level for prediction intervals.}
\item{K}{Maximum number of years to use in forecasting coefficients for ratio components.}
\item{drange}{Range of fractional differencing parameter for the ratio coefficients.}
\item{...}{Other arguments as for \code{\link{forecast.fdm}}.}
}

\description{The product and ratio models from \code{\link{coherentfdm}} are forecast, and the results combined to give forecasts for each group in the original data.}

\value{Object of class \code{fmforecast2} containing a list of objects each of class \code{fmforecast}. The forecasts for each group in the original data are given first. Then the forecasts from the product model, and finally a list of forecasts from each of the ratio models.}


\author{Rob J Hyndman.}

\seealso{\code{\link{coherentfdm}}, \code{\link{forecast.fdm}}.}

\examples{
fr.short <- extract.years(fr.sm,1950:2006)
fr.fit <- coherentfdm(fr.short)
fr.fcast <- forecast(fr.fit)
plot(fr.fcast$male)
plot(fr.fcast$ratio$male, plot.type='component', components=3)
models(fr.fcast)
}
\keyword{models}
