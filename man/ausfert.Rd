\name{aus.fert}
\alias{aus.fert}

\docType{data}

\title{Australian fertility data}

\description{Age-specific fertility rates and female child-bearing population for Australia.}

\details{Australian fertility rates and populations (1921-2002) for age groups (<20, 20-24, 25-29, 30-34, 35-39, 40-44, 45+).
Data taken from v3.2b of the Australian Demographic Data Bank released 10 February 2005.
}

\format{Object of class \code{demogdata} containing the following components:
\describe{\item{year}{Vector of years}
\item{age}{Vector of ages}
\item{rate}{List containing one matrix with one age group per row and one column per year.}
\item{pop}{Population data in same form as \code{rate}.}
\item{type}{Type of object. In this case, \dQuote{fertility}.}
\item{label}{Character string giving area from which data are taken. In this case, \dQuote{Australia}.}
}}

\source{The Australian Demographic Data Bank (courtesy of Len Smith).}

\author{Rob J Hyndman}

\examples{
plot(aus.fert)
}
\keyword{data}
