\name{isfe}
\alias{isfe}
\alias{isfe.demogdata}
\title{Integrated Squared Forecast Error for models of various orders}
\usage{
\method{isfe}{demogdata}(data, series=names(data$rate)[1], max.order=N-3, N=10, h=5:10,
    ages=data$age, max.age=max(ages), method=c("classical", "M", "rapca"),
    fmethod = c("arima","ar","arfima","ets","ets.na","struct","rwdrift","rw"), 
    lambda=3, ...)
}

\arguments{
\item{data}{demogdata object.}
\item{series}{name of series within data holding rates (1x1)}
\item{ages}{Ages to include in fit.}
\item{max.age}{Maximum age to fit.}
\item{max.order}{Maximum number of basis functions to fit.}
\item{N}{Minimum number of functional observations to be used in fitting a model.}
\item{h}{Forecast horizons over which to average.}
\item{method}{Method to use for principal components decomposition. Possibilities are \dQuote{M}, \dQuote{rapca} and \dQuote{classical}.}
\item{fmethod}{Method used for forecasting. Current possibilities are \dQuote{ets}, \dQuote{arima}, \dQuote{ets.na},
	\dQuote{struct}, \dQuote{rwdrift} and \dQuote{rw}.}
\item{lambda}{Tuning parameter for robustness when \code{method="M"}.}
\item{...}{Additional arguments control the fitting procedure.}
}
\description{Computes ISFE values for functional time series models of various orders.
}

\value{Numeric matrix with \code{(max.order+1)} rows and \code{length(h)} columns
containing ISFE values for models of orders 0:max.order.}

\author{Rob J Hyndman.}

\references{Hyndman, R.J., and Ullah, S. (2007) Robust forecasting of mortality and
fertility rates: a functional data approach. \emph{Computational Statistics & Data Analysis}, 
\bold{51}, 4942-4956. \url{http://robjhyndman.com/papers/funcfor}
}

\seealso{\code{\link{fdm}}, \code{\link{forecast.fdm}}.}

\keyword{models}
