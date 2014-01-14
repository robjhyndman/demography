\name{update.fmforecast}
\alias{update.fmforecast}
\alias{update.fmforecast2}
\title{Updating functional demographic models and coherent functional demographic models.}
\usage{\method{update}{fmforecast}(object, ...)
\method{update}{fmforecast2}(object, ...)
}

\arguments{
\item{object}{Output from either \code{\link{fdm}} or \code{\link{coherentfdm}}.}
\item{...}{Extra arguments currently ignored.}
}

\description{\code{update.fmforecast()} updates \code{fdm} forecasts. The argument \code{object} is the output from \code{\link{forecast.fdm}} which has been subsequently modified with new coefficient forecasts. These new forecasts are used when re-calculating the forecast of the mortality or fertility rates, or net migration numbers.

\code{update.fmforecast2()} updates \code{fdmpr} forecasts. The argument \code{object} is the output from \code{\link{forecast.fdmpr}} which has been subsequently modified with new coefficient forecasts. 
}

\value{A list of the same class as \code{object}.
}

\author{Rob J Hyndman.}

\seealso{\code{\link{forecast.fdm}}, \code{\link{forecast.fdmpr}}}

\examples{
\dontrun{
france.fit <- fdm(fr.mort,order=2)
france.fcast <- forecast(france.fit,50)
# Replace first coefficient model with ARIMA(0,1,2)+drift
france.fcast$coeff[[2]] <- forecast(Arima(france.fit$coeff[,2], 
  order=c(0,1,2), include.drift=TRUE), h=50, level=80)
france.fcast <- update(france.fcast)

fr.short <- extract.years(fr.sm,1950:2006)
fr.fit <- coherentfdm(fr.short)
fr.fcast <- forecast(fr.fit)
par(mfrow=c(1,2))
plot(fr.fcast$male)
# Replace first coefficient model in product component with a damped ETS model:
fr.fcast$product$coeff[[2]] <- forecast(ets(fr.fit$product$coeff[,2], damped=TRUE), 
  h=50, level=80)
fr.fcast <- update(fr.fcast)
plot(fr.fcast$male)
}
}
\keyword{models}
