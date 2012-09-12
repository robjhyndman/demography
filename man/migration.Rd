\name{netmigration}
\alias{netmigration}
\title{Calculate net migration from mortality and fertility data}
\usage{netmigration(mort, fert, mfratio = 1.05)
}

\arguments{
\item{mort}{Demogdata object of type \code{"mortality"}}
\item{fert}{Demogdata object of type \code{"fertility"}}
\item{mfratio}{Male-female ratio to be used in simulating births.}
}

\description{Function to compute the net number of migrants in each year and for each age, 
based on the total population numbers, deaths and births in each year.}

\value{Object of class \dQuote{demogdata} with the following components:
\item{year}{Vector of years}
\item{age}{Vector of ages}
\item{rate}{List containing matrices of net migration numbers (not "rates") with with one age group per 
    row and one column per year. Names of matrices are the same as for \code{mort$rate}.}
\item{pop}{List containing matrices of populations in same form as \code{rate} and containing population numbers.}
\item{type}{Type of object. In this case, \dQuote{migration}.}
\item{label}{label from \code{mort$label}}
}

\seealso{\code{\link{demogdata}}}

\author{Rob J Hyndman}

\examples{
\dontrun{
require(addb)
aus.mig <- netmigration(australia,aus.fertility)
plot(aus.mig)}
}
 
\keyword{manip}
