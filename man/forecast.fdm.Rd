\name{forecast.fdm}
\alias{forecast.fdm}
\title{Forecast functional demographic model.}
\usage{
\method{forecast}{fdm}(object, h = 50, level=80, jumpchoice = c("fit", "actual"), method =
				 "arima", warnings=FALSE, ...)
}

\arguments{
\item{object}{Output from \code{\link{fdm}}.}
\item{h}{Forecast horizon.}
\item{level}{Confidence level for prediction intervals.}
\item{jumpchoice}{If "actual", the forecasts are bias-adjusted by the difference between the fit and the last year of observed data.
	Otherwise, no adjustment is used.}
\item{method}{Forecasting method to be used.}
\item{warnings}{If TRUE, warnings arising from the forecast models for coefficients will be shown. Most of these can be ignored, so the default is \code{warnings=FALSE}.}
\item{...}{Other arguments as for \code{\link[ftsa]{forecast.ftsm}}.}
}

\description{The coefficients from the fitted object are forecast using a univariate time series model.
The forecast coefficients are then multiplied by the basis functions to
obtain a forecast demographic rate curve.}

\value{Object of class \code{fmforecast} with the following components:
\item{label}{Name of region from which the data are taken.}
\item{age}{Ages from \code{lcaout} object.}
\item{year}{Years from \code{lcaout} object.}
\item{rate}{List of matrices containing forecasts, lower bound and upper bound of prediction intervals.
	Point forecast matrix takes the same name as the series that has been forecast.}
\item{error}{Matrix of one-step errors for historical data}
\item{fitted}{Matrix of one-step forecasts for historical data}
\item{coeff}{List of objects of type \code{forecast} containing the coefficients and their forecasts.}
\item{coeff.error}{One-step errors for each of the coefficients.}
\item{var}{List containing the various components of variance: model, error, mean, total and coeff.}
\item{model}{Fitted model in \code{obj}.}
\item{type}{Type of data: \dQuote{mortality}, \dQuote{fertility} or \dQuote{migration}.}
}

\author{Rob J Hyndman.}

\seealso{\code{\link{fdm}}, \code{\link{forecast.lca}}, \code{\link[ftsa]{forecast.ftsm}}.}
\examples{
france.fit <- fdm(fr.mort,order=2)
france.fcast <- forecast(france.fit,50)
plot(france.fcast)
models(france.fcast)
}
\keyword{models}
