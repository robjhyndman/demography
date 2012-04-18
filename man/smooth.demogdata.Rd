\name{smooth.demogdata}
\alias{smooth.demogdata}
\title{Create smooth demogdata functions}
\usage{smooth.demogdata(data, method = switch(data$type, mortality = "mspline",
    fertility = "cspline", migration = "loess"), age.grid,
    power = switch(data$type, mortality = 0.4, fertility = 1, migration = 1),
    b = 65, k = 30, span = 0.2, lambda = 1e-10, interpolate = FALSE,
    weight = data$type != "migration", obs.var = "empirical")
}

\arguments{
\item{data}{Demogdata object such as created using \code{\link{read.demogdata}}.}
\item{method}{Method of smoothing. Possibilities: \code{"mspline"} (monotonic regression splines),
    \code{"cspline"} (concave regression splines),
    \code{"spline"} (unconstrained regression splines),
    \code{"loess"} (local quadratic using \code{\link{loess}}).}
\item{age.grid}{Ages to use for smoothed curves. Default is single years over a slightly greater range than the unsmoothed data.}
\item{power}{Power transformation for age variable before smoothing. Default is 0.4 for mortality data and 1 (no transformation) for fertility or migration data.}
\item{b}{Lower age for monotonicity if \code{method=="mspline"}. Above this, the smooth curve
is assumed to be monotonically increasing.}
\item{k}{Number of knots to use for penalized regression spline estimate. Ignored if \code{method=="loess"}.}
\item{span}{Span for loess smooth if \code{method=="loess"}.}
\item{lambda}{Penalty for constrained regression spline if \code{method=="cspline"}.}
\item{interpolate}{If \code{interpolate==TRUE}, a linear interpolation is used instead of smoothing.}
\item{weight}{If TRUE, uses weighted smoothing.}
\item{obs.var}{Method for computing observational variance. Possible values: \dQuote{empirical} or \dQuote{theoretical}.}
}

\description{Smooth demogdata data using one of four methods depending on the value of \code{method}}

\details{The value of \code{method} determines the type of smoothing used.
\describe{
\item{method="mspline"}{Weighted penalized regression splines with a monotonicity constraint. The curves are monotonically
    increasing for age greater than b. Smoothness controlled by \code{k}. Methodology based on Wood (1994). Code calls \code{\link[mgcv]{gam}} for the basic
    computations.}
\item{method="cspline"}{Weighted regression B-splines with a concavity constraint. Smoothness controlled by \code{lambda}.
Methodology based on He and Ng (1999). Code calls \code{\link[cobs]{cobs}} for the basic computations.}
\item{method="spline"}{Unconstrained weighted penalized regression splines. Equivalent to "mspline" but with \code{b=Inf}.}
\item{method="loess"}{Weighted locally quadratic regression. Smoothness controlled by span. Code calls
\code{\link{loess}} for the basic computations.}
}
}

\value{Demogdata object identical to \code{data} except all
rate matrices are replaced with smooth versions and pop matrices are replaced with disaggregated population estimates obtained
using monotonic spline interpolation applied to the cumulative population data.
Weight
matrices are also added to the object showing the inverse
variances of the estimated smooth curves. }

\author{Rob J Hyndman}
\examples{
france.sm <- smooth.demogdata(extract.years(fr.mort,1980:1997))
plot(france.sm)
plot(fr.mort,years=1980,type="p",pch=1)
lines(france.sm,years=1980,col=2)
}
\keyword{smooth}
