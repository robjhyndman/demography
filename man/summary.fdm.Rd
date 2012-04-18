\name{summary.fdm}
\alias{summary.fdm}
\alias{summary.lca}
\title{Summary for functional demographic model or Lee-Carter model}
\usage{\method{summary}{fdm}(object,...)
\method{summary}{lca}(object,...)
}

\arguments{
\item{object}{Output from \code{\link{fdm}} or \code{\link{lca}}.}
\item{...}{Other arguments.}
}

\description{Summarizes a basis function model fitted to age-specific demographic rate data.
	It returns various measures of goodness-of-fit.}

\value{None.}

\author{Rob J Hyndman.}

\seealso{\code{\link{fdm}}, \code{\link{lca}}, \code{\link{bms}}, \code{\link{compare.demogdata}}}

\examples{
fit1 <- lca(fr.mort)
fit2 <- bms(fr.mort,breakmethod="bai")
fit3 <- fdm(fr.mort)
summary(fit1)
summary(fit2)
summary(fit3)
}
\keyword{models}
