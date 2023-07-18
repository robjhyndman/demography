#' Read data from HMD and construct a mortality demogdata object
#'
#' \code{hmd.mx} reads "Mx" (1x1) data from the Human Mortality Database (HMD
#' \url{https://www.mortality.org}) and constructs a demogdata object suitable
#' for plotting using \code{\link{plot.demogdata}} and fitting an LC or BMS
#' model using \code{\link{lca}} or an FDA model using \code{\link{fdm}}.
#' \code{hmd.pop} reads "Population" (1x1) data from the HMD and constructs a
#' demogdata object suitable for plotting using \code{\link{plot.demogdata}}.
#' \code{hmd.e0} reads life expectancy at birth from the HMD and returns the
#' result as a \code{ts} object.
#'
#' In order to read the data, users are required to create their account via the HMD website (\url{https://www.mortality.org}),
#' and obtain a valid username and password. 
#'
#' @param country Directory abbreviation from the HMD. For instance, Australia =
#'   "AUS".
#' @param username HMD username (case-sensitive)
#' @param password HMD password (case-sensitive)
#' @param label Character string giving name of country from which the data are
#'   taken.
#'
#' @return \code{hmd.mx} returns an object of class \code{demogdata} with the following components:
#' \item{year}{Vector of years}
#' \item{age}{Vector of ages}
#' \item{rate}{A list containing one or more rate matrices with one age group per row and one column per year.}
#' \item{pop}{A list of the same form as \code{rate} but containing population numbers instead of demographic rates.}
#' \item{type}{Type of object: \dQuote{mortality}, \dQuote{fertility} or \dQuote{migration}.}
#' \item{label}{label}
#' \code{hmd.pop} returns a similar object but without the \code{rate} component.
#' \code{hmd.e0} returns an object of class \code{ts} with columns \code{male}, \code{female} and \code{total}.
#'
#' @seealso \code{\link{demogdata}},\code{\link{read.demogdata}},\code{\link{plot.demogdata}}, \code{\link{life.expectancy}}
#' @author Rob J Hyndman
#' @examples
#' \dontrun{
#' norway <- hmd.mx("NOR", username, password, "Norway")
#' summary(norway)
#' }
#' @keywords manip
#' @name hmd
#' @export
hmd.mx <- function(country, username, password, label = country) {
  # Read raw MX and Exposure data
  mx <- HMDHFDplus::readHMDweb(country, item = "Mx_1x1",
    username = username, password = password, fixup = TRUE)
  pop <- HMDHFDplus::readHMDweb(country, item = "Exposures_1x1",
    username = username, password = password, fixup = TRUE)

  # Construct output
  obj <- list(type = "mortality", label = label, lambda = 0)
  obj$year <- sort(unique(mx[, "Year"]))
  n <- length(obj$year)
  m <- length(unique(mx[, "Age"]))
  obj$age <- mx[seq(m), "Age"]
  mnames <- names(mx)[-c(1:2, NCOL(mx))]
  n.mort <- length(mnames)
  obj$rate <- obj$pop <- list()
  for (i in seq(n.mort)) {
    obj$rate[[i]] <- matrix(mx[, i + 2], nrow = m, ncol = n)
    obj$rate[[i]][obj$rate[[i]] < 0] <- NA
    obj$pop[[i]] <- matrix(pop[, i + 2], nrow = m, ncol = n)
    obj$pop[[i]][obj$pop[[i]] < 0] <- NA
    dimnames(obj$rate[[i]]) <- dimnames(obj$pop[[i]]) <- list(obj$age, obj$year)
  }
  names(obj$pop) <- names(obj$rate) <- tolower(mnames)

  return(structure(obj, class = "demogdata"))
}

#' @rdname hmd
#' @export
hmd.e0 <- function(country, username, password) {
  # Read raw e0 data
  lt <- HMDHFDplus::readHMDweb(country, item = "E0per",
    username = username, password = password, fixup = TRUE)
  # Convert to a ts object
  ts(lt[, -1], start = lt[1, 1], frequency = 1)
}


#' @rdname hmd
#' @export
hmd.pop <- function(country, username, password, label = country) {
  # Read raw data
  pop <- HMDHFDplus::readHMDweb(country, item = "Population",
    username = username, password = password, fixup = FALSE)

  # Only keep 1 January populations
  pop <- pop[, !grepl("2$", colnames(pop))]
  mnames <- names(pop)
  mnames <- sub("1$", "", mnames)
  names(pop) <- mnames

  # Construct output
  obj <- list(type = "population", label = label, lambda = 0)
  obj$year <- sort(unique(pop[, "Year"]))
  n <- length(obj$year)
  m <- length(unique(pop[, "Age"]))
  obj$age <- pop[seq(m), "Age"]
  pop <- pop[, !(mnames %in% c("Year","Age","OpenInterval"))]
  n.pop <- NCOL(pop)
  obj$pop <- list()
  for (i in seq(n.pop)) {
    obj$pop[[i]] <- matrix(pop[, i], nrow = m, ncol = n)
    obj$pop[[i]][obj$pop[[i]] < 0] <- NA
    dimnames(obj$pop[[i]]) <- list(obj$age, obj$year)
  }
  names(obj$pop) <- tolower(colnames(pop))

  return(structure(obj, class = "demogdata"))
}
