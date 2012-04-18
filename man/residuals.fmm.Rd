\name{residuals.fdm}
\alias{residuals.fdm}
\alias{residuals.lca}
\alias{fitted.fdm}
\alias{fitted.lca}
\title{Compute residuals and fitted values from functional demographic model or Lee-Carter model}
\usage{\method{residuals}{fdm}(object,...)
\method{residuals}{lca}(object,...)
\method{fitted}{fdm}(object,...)
\method{fitted}{lca}(object,...)}

\description{After fitting a Lee-Carter model or functional
demographic model, it is useful to inspect the residuals or plot
the fitted values. These functions extract the relevant
information from the fit object.}

\arguments{
\item{object}{Output from \code{\link{fdm}} or \code{\link{lca}}.}
\item{...}{Other arguments.}
}

\value{\code{residuals.fdm} and \code{residuals.lca} produce an object of class \dQuote{fmres}
containing the residuals from the model.

\code{fitted.fdm} and \code{fitted.lca} produce an object of class \dQuote{fts} containing the fitted values
from the model.}

\author{Rob J Hyndman.}

\seealso{\code{\link{fdm}}, \code{\link{lca}}, \code{\link{bms}}}

\examples{
fit1 <- lca(fr.mort)
plot(residuals(fit1))
plot(fitted(fit1))
}
\keyword{models}
