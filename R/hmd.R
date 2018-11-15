# Function to construct a mortality demogdata object from HMD

#' Summary for functional demographic model or Lee-Carter model
#'
#' \code{hmd.mx} reads "Mx" (1x1) data from the Human Mortality Database (HMD
#' \url{http://www.mortality.org}) and constructs a demogdata object suitable
#' for plotting using \code{\link{plot.demogdata}} and fitting an LC or BMS
#' model using \code{\link{lca}} or an FDA model using \code{\link{fdm}}.
#' \code{hmd.pop} reads "Population" (1x1) data from the HMD and constructs a
#' demogdata object suitable for plotting using \code{\link{plot.demogdata}}.
#' \code{hmd.e0} reads life expectancy at birth from the HMD and returns the
#' result as a \code{ts} object.
#'
#' In order to read the data, users are required to create their account via the HMD website (\url{http://www.mortality.org}),
#' and obtain a valid username and password.
#' 
#' The country codes (as at 23 December 2016) are as follows.
#' \tabular{ll}{
#'   Australia \tab AUS\cr
#'   Austria \tab AUT\cr
#'   Belarus \tab BLR\cr
#'   Belgium \tab BEL\cr
#'   Bulgaria \tab BGR\cr
#'   Canada \tab CAN\cr
#'   Chile \tab CHL\cr
#'   Czech Republic \tab CZE\cr
#'   Denmark \tab DNK\cr
#'   Estonia \tab EST\cr
#'   Finland \tab FIN\cr
#'   France\cr
#'   -- France total population \tab FRATNP\cr
#'   -- France civilian population \tab FRACNP\cr
#'   Germany\cr
#'   -- Germany total population \tab DEUTNP\cr
#'   -- West Germany \tab DEUTFRG\cr
#'   -- East Germany \tab DEUTGDR\cr
#'   Greece \tab GRC\cr
#'   Hungary \tab HUN\cr
#'   Iceland \tab ISL\cr
#'   Ireland \tab IRL\cr
#'   Israel \tab ISR\cr
#'   Italy \tab ITA\cr
#'   Japan \tab JPN\cr
#'   Latvia \tab LVA\cr
#'   Lithuania \tab LTU\cr
#'   Luxembourg \tab LUX\cr
#'   Netherlands \tab NLD\cr
#'   New Zealand	\cr
#'   -- NZ total population \tab NZL_NP\cr
#'   -- NZ Maori \tab NZL_MA\cr
#'   -- NZ non-Maori \tab NZL_NM\cr
#'   Norway \tab NOR\cr
#'   Poland \tab POL\cr
#'   Portugal \tab PRT\cr
#'   Russia \tab RUS\cr
#'   Slovakia \tab SVK\cr
#'   Slovenia \tab SVN\cr
#'   Spain \tab ESP\cr
#'   Sweden \tab SWE\cr
#'   Switzerland \tab CHE\cr
#'   Taiwan \tab TWN\cr
#'   United Kingdom\cr
#'   -- UK Total Population \tab GBR_NP\cr
#'   -- England & Wales Total Population \tab GBRTENW\cr
#'   -- England & Wales Civilian Population \tab GBRCENW\cr
#'   -- Scotland \tab GBR_SCO\cr
#'   -- Northern Ireland \tab GBR_NIR\cr
#'   U.S.A. \tab USA\cr
#'   Ukraine \tab UKR\cr
#' }
#' 
#' 
#' 
#' @param country Directory abbreviation from the HMD. For instance, Australia =
#'   "AUS". See below for other countries.
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
#' summary(norway)}
#' @keywords manip
#' @name hmd
#' @export
hmd.mx <- function(country, username, password, label=country)
{
    path <- paste("https://www.mortality.org/hmd/", country, "/STATS/", "Mx_1x1.txt", sep = "")
    userpwd <- paste(username, ":", password, sep = "")
    txt <- RCurl::getURL(path, userpwd = userpwd)
    con <- textConnection(txt)
    mx <- try(utils::read.table(con, skip = 2, header = TRUE, na.strings = "."),TRUE)
    close(con)
    if(class(mx)=="try-error")
        stop("Connection error at www.mortality.org. Please check username, password and country label.")

    path <- paste("https://www.mortality.org/hmd/", country, "/STATS/", "Exposures_1x1.txt", sep = "")
    userpwd <- paste(username, ":", password, sep = "")
    txt <- RCurl::getURL(path, userpwd = userpwd)
    con <- textConnection(txt)
    pop <- try(utils::read.table(con, skip = 2, header = TRUE, na.strings = "."),TRUE)
    close(con)
    if(class(pop)=="try-error")
        stop("Exposures file not found at www.mortality.org")

    obj <- list(type="mortality",label=label,lambda=0)

    obj$year <- sort(unique(mx[, 1]))
    #obj$year <- ts(obj$year, start=min(obj$year))
    n <- length(obj$year)
    m <- length(unique(mx[, 2]))
    obj$age <- mx[1:m, 2]
    mnames <- names(mx)[-c(1, 2)]
    n.mort <- length(mnames)
    obj$rate <- obj$pop <- list()
    for (i in 1:n.mort)
    {
        obj$rate[[i]] <- matrix(mx[, i + 2], nrow = m, ncol = n)
        obj$rate[[i]][obj$rate[[i]] < 0] <- NA
        obj$pop[[i]] <- matrix(pop[, i + 2], nrow = m, ncol = n)
        obj$pop[[i]][obj$pop[[i]] < 0] <- NA
        dimnames(obj$rate[[i]]) <- dimnames(obj$pop[[i]]) <- list(obj$age, obj$year)
    }
    names(obj$pop) = names(obj$rate) <- tolower(mnames)

    obj$age <- as.numeric(as.character(obj$age))
    if (is.na(obj$age[m]))
        obj$age[m] <- 2 * obj$age[m - 1] - obj$age[m - 2]
    return(structure(obj, class = "demogdata"))
}

#' @rdname hmd
#' @export
hmd.e0 <- function(country, username, password)
{
    path <- paste("https://www.mortality.org/hmd/", country, "/STATS/", "E0per.txt", sep = "")
    userpwd <- paste(username, ":", password, sep = "")
    txt <- RCurl::getURL(path, userpwd = userpwd)
    con <- textConnection(txt)
    lt <- try(utils::read.table(con, skip = 2, header = TRUE, na.strings = "."), TRUE)
    close(con)
    if (class(lt) == "try-error")
        stop("Life expectancy file not found at www.mortality.org")
	lt <- ts(lt[,-1],start=lt[1,1],frequency=1)
    return(lt)
}


#' @rdname hmd
#' @export
hmd.pop <- function(country, username, password, label=country)
{
    path <- paste("https://www.mortality.org/hmd/", country, "/STATS/", "Population.txt", sep = "")
    userpwd <- paste(username, ":", password, sep = "")
    txt <- RCurl::getURL(path, userpwd = userpwd)
    con <- textConnection(txt)
    pop <- try(utils::read.table(con, skip = 2, header = TRUE, na.strings = "."),TRUE)
    close(con)
    if(class(pop)=="try-error")
        stop("Population file not found at www.mortality.org")

    obj <- list(type="population",label=label,lambda=0)

    obj$year = sort(unique(pop[, 1]))
    #obj$year <- ts(obj$year, start=min(obj$year))
    n <- length(obj$year)
    m <- length(unique(pop[, 2]))
    obj$age <- pop[1:m, 2]
    mnames <- names(pop)[-c(1, 2)]
    n.pop <- length(mnames)
    obj$pop <- list()
    for (i in 1:n.pop)
    {
        obj$pop[[i]] <- matrix(pop[, i + 2], nrow = m, ncol = n)
        obj$pop[[i]][obj$pop[[i]] < 0] <- NA
        dimnames(obj$pop[[i]]) <- list(obj$age, obj$year)
    }
    names(obj$pop) <- tolower(mnames)

    obj$age <- as.numeric(as.character(obj$age))
    if (is.na(obj$age[m]))
        obj$age[m] <- 2 * obj$age[m - 1] - obj$age[m - 2]
    return(structure(obj, class = "demogdata"))
}

