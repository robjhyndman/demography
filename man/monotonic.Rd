\name{cm.spline}
\alias{cm.spline}
\alias{cm.splinefun}
\title{Monotonic interpolating splines}

\description{Perform cubic spline monotonic interpolation of given data points, returning either a list of points obtained by the interpolation or a function performing the interpolation. The splines are constrained to be monotonically increasing (i.e., the slope is never negative).}

\usage{cm.splinefun(x, y = NULL, ...)
cm.spline(x, y = NULL, n = 3*length(x), 
       xmin = min(x), xmax = max(x), ...)
}

\arguments{
\item{x,y}{vectors giving the coordinates of the points to be interpolated. Alternatively a single plotting structure can be specified: see \code{\link[grDevices]{xy.coords}}.}
\item{n}{interpolation takes place at n equally spaced points spanning the interval [\code{xmin}, \code{xmax}].}
\item{xmin}{left-hand endpoint of the interpolation interval.}
\item{xmax}{right-hand endpoint of the interpolation interval.}
\item{...}{Other arguments are ignored.}
}

\details{From version 1.14 in R 2.15.2 or later, these are simply wrappers to the \code{\link[stats]{spline}} and \code{\link[stats]{splinefun}} functions from the 
stats package.}

\value{
\item{cm.spline}{returns a list containing components \code{x} and \code{y} which give the ordinates where interpolation took place and the interpolated values.}
\item{cm.splinefun}{returns a function which will perform cubic spline interpolation of the given data points. This is often more useful than \code{spline}.}
}


\references{
Forsythe, G. E., Malcolm, M. A. and Moler, C. B. (1977) \emph{Computer Methods for Mathematical Computations}.

Hyman (1983) \emph{SIAM J. Sci. Stat. Comput.} \bold{4}(4):645-654.

Dougherty, Edelman and Hyman 1989 \emph{Mathematics of Computation}, \bold{52}: 471-494.
}

\author{Rob J Hyndman}
\examples{
x <- seq(0,4,l=20)
y <- sort(rnorm(20))
plot(x,y)
lines(spline(x, y, n = 201), col = 2) # Not necessarily monotonic
lines(cm.spline(x, y, n = 201), col = 3) # Monotonic
}
\keyword{smooth}
