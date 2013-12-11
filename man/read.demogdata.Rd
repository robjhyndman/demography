\name{read.demogdata}
\alias{read.demogdata}
\title{Read demographic data and construct demogdata object}
\usage{read.demogdata(file, popfile, type, label, max.mx = 10, skip = 2,
                 popskip = skip, lambda, scale=1)
}

\arguments{
\item{file}{Filename containing demographic rates.}
\item{popfile}{Filename containing population numbers.}
\item{type}{Character string showing type of demographic series:
    either \dQuote{mortality}, \dQuote{fertility} or \dQuote{migration}.}
\item{label}{Name of area from which the data are taken.}
\item{max.mx}{Maximum allowable value for demographic rate. All values greater than max.mx will be set to max.mx.}
\item{skip}{Number of lines to skip at the start of \code{file}.}
\item{popskip}{Number of lines to skip at the start of \code{popfile}.}
\item{lambda}{Box-Cox transformation parameter to be used in modelling and plotting. If missing, default values are 0 (for mortality), 0.4 (for fertility) and 1 (for migration).}
\item{scale}{Number of people in the rate definition. \code{scale=1} indicates the rates are per person; \code{scale=1000} indicates the rates are per 1000 people.}
}

\description{Read data from text files and construct a demogdata object
suitable for plotting using \code{\link{plot.demogdata}} and
fitting an LC or BMS model using \code{\link{lca}} or an FDA model using \code{\link{fdm}}.}

\details{All data are assumed to be tab-delimited text files with the first column
containing the year of observation and the second column containing
the age level. All remaining columns are assumed to be demographic
rates for sections of the population. The first row of the text file
is assumed to contain the names of each column.
Population data are assumed to have the same format but with population numbers in place
of  rates.  The columns names in the two files should be identical. Note that this format is what is used by the Human
Mortality Database \url{http://www.mortality.org}. If \code{popfile} contains the
Exposures and \code{file} contains the Mx rates from the HMD, then everything will work seamlessly.}

\value{Object of class \dQuote{demogdata} with the following components:
\item{year}{Vector of years}
\item{age}{Vector of ages}
\item{rate}{A list containing one or more rate matrices with one age group per row and one column per year.}
\item{pop}{A list of the same form as \code{rate} but containing population numbers instead of demographic rates.}
\item{type}{Type of object: \dQuote{mortality}, \dQuote{fertility} or \dQuote{migration}.}
\item{label}{label}
}

\seealso{\code{\link{demogdata}}}

\examples{
\dontrun{
norway <- read.demogdata("Mx_1x1.txt", "Exposures_1x1.txt",
    type="mortality", label="Norway")}
}

\author{Rob J Hyndman}
\keyword{manip}
