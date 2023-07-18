#####################################################
### FUNCTIONS FOR HANDLING DEMOGRAPHIC RATE DATA ###
####################################################

## BASIC DEMOGDATA STRUCTURE:
# year - time series
# age - numeric vector of mid points of age groups
# rate - main rate (mortality, fertility or migration, etc)
#          matrices: female, male and total.
# pop - population matrices matching rates
# label (e.g., country, state, etc.) - character string
# type - mortality, fertility, migration, etc. character string.

#' Create demogdata object from raw data matrices
#'
#' Create demogdata object suitable for plotting using \code{\link{plot.demogdata}} and
#' fitting an LC or BMS model using \code{\link{lca}} or an FDA model using \code{\link{fdm}}.
#'
#' @param data Matrix of data: either mortality rates or fertility rates
#' @param pop Matrix of population values of same dimension as data.
#' These are population numbers as at 30 June of each year (i.e., the "exposures").
#' So, for example, the number of deaths is data*pop if data contains mortality rates.
#' @param ages Vector of ages corresponding to rows of \code{data}.
#' @param years Vector of years corresponding to columns of \code{data}.
#' @param type Character string showing type of demographic series:
#' either \dQuote{mortality}, \dQuote{fertility} or \dQuote{migration}.
#' @param label Character string of the name of area from which the data are taken.
#' @param name Name of series: usually male, female or total.
#' @param lambda Box-Cox transformation parameter.
#'
#' @return Object of class \dQuote{demogdata} with the following components:
#' \item{year}{Vector of years}
#' \item{age}{Vector of ages}
#' \item{rate}{A list containing one or more rate matrices with one age group per row and one column per year.}
#' \item{pop}{A list of the same form as \code{rate} but containing population numbers instead of demographic rates.}
#' \item{type}{Type of object: \dQuote{mortality}, \dQuote{fertility} or \dQuote{migration}.}
#' \item{label}{label}
#' \item{lambda}{lambda}
#'
#' @seealso \code{\link{read.demogdata}}
#' @author Rob J Hyndman
#' @keywords manip
#' @export
demogdata <- function(data, pop, ages, years, type, label, name, lambda)
{
  p <- nrow(data)
  n <- ncol(data)
  if(nrow(pop) != p | ncol(pop) != n)
    stop("data and pop are of different size")
  if(length(ages) != p)
    stop("Number of ages doesn't match data")
  if(length(years) != n)
    stop("Number of years doesn't match data")

  types <- c("mortality","fertility","migration","population")
  idx <- pmatch(type,types)
  if(is.na(idx))
    warning("Unknown type")
  else
    type <- types[idx]

  if(missing(lambda))
  {
    if(type=="mortality")
      lambda <- 0
    else if(type=="fertility")
      lambda <- 0.4
    else
      lambda <- 1
  }

  if(type=="population")
  {
    obj <- list(year=years, age=ages, pop=list(as.matrix(pop)), type=type,
        label=label, lambda=lambda)
    dimnames(obj$pop[[1]]) <- list(ages,years)
  }
  else
  {
    obj <- list(year=years, age=ages, rate=list(as.matrix(data)), pop=list(as.matrix(pop)), type=type,
        label=label, lambda=lambda)
    dimnames(obj$rate[[1]]) <- dimnames(obj$pop[[1]]) <- list(ages,years)
    names(obj$rate) <- names(obj$pop) <- name
  }
  return(structure(obj,class="demogdata"))
}

# Function to read demographic rates from file
# Assumed input format:
#  col 1: year
#  col 2: age
#  col 3+: rates for males, females and total (any order but labels required)
# Assume age cycles within year and all in order.
# Population file of same format.
# Output format:  object of class demogdata
# scale indicates the rates. scale=1 means per person. scale=1000 means per 1000 people.

#' Read demographic data and construct demogdata object
#'
#' Read data from text files and construct a demogdata object suitable for
#' plotting using \code{\link{plot.demogdata}} and fitting an LC or BMS model
#' using \code{\link{lca}} or an FDA model using \code{\link{fdm}}.
#'
#' All data are assumed to be tab-delimited text files with the first column
#' containing the year of observation and the second column containing the age
#' level. All remaining columns are assumed to be demographic rates for sections
#' of the population. The first row of the text file is assumed to contain the
#' names of each column. Population data are assumed to have the same format but
#' with population numbers in place of  rates.  The columns names in the two
#' files should be identical. Note that this format is what is used by the Human
#' Mortality Database \url{http://www.mortality.org}. If \code{popfile} contains
#' the Exposures and \code{file} contains the Mx rates from the HMD, then
#' everything will work seamlessly.
#'
#' @param file Filename containing demographic rates.
#' @param popfile Filename containing population numbers.
#' @param type Character string showing type of demographic series:
#'     either \dQuote{mortality}, \dQuote{fertility} or \dQuote{migration}.
#' @param label Name of area from which the data are taken.
#' @param max.mx Maximum allowable value for demographic rate. All values greater than max.mx will be set to max.mx.
#' @param skip Number of lines to skip at the start of \code{file}.
#' @param popskip Number of lines to skip at the start of \code{popfile}.
#' @param lambda Box-Cox transformation parameter to be used in modelling and plotting. If missing, default values are 0 (for mortality), 0.4 (for fertility) and 1 (for migration).
#' @param scale Number of people in the rate definition. \code{scale=1} indicates the rates are per person; \code{scale=1000} indicates the rates are per 1000 people.
#'
#' @return Object of class \dQuote{demogdata} with the following components:
#'   \item{year}{Vector of years} \item{age}{Vector of ages} \item{rate}{A list
#'   containing one or more rate matrices with one age group per row and one
#'   column per year.} \item{pop}{A list of the same form as \code{rate} but
#'   containing population numbers instead of demographic rates.}
#'   \item{type}{Type of object: \dQuote{mortality}, \dQuote{fertility} or
#'   \dQuote{migration}.} \item{label}{label}
#'
#' @seealso \code{\link{demogdata}}
#'
#' @examples
#'
#'   \dontrun{ norway <- read.demogdata("Mx_1x1.txt",
#'   "Exposures_1x1.txt", type="mortality", label="Norway")}
#' @author Rob J Hyndman
#' @keywords manip
#'
#' @export
read.demogdata <- function(file,popfile,type,label,max.mx=10,skip=2,popskip=skip,lambda, scale=1)
{
  if(missing(lambda))
  {
    if(type=="mortality")
      lambda <- 0
    else if(type=="fertility")
      lambda <- 0.4
    else
      lambda <- 1
  }

  mfile <- !missing(file)
  mpopfile <- !missing(popfile)

  obj <- list(type=type,label=label,lambda=lambda)

  if(mfile)
  {
    tmp1 <- utils::read.table(file,header=TRUE,na.strings=".",skip=skip)
    obj$year=sort(unique(tmp1[,1]))
    n <- length(obj$year)
    m <- length(unique(tmp1[,2]))
    obj$age <- tmp1[1:m,2]
    mnames <- names(tmp1)[-c(1,2)]
    n.mort <- length(mnames)
    obj$rate <- list()
    for(i in 1:n.mort)
    {
      obj$rate[[i]] = matrix(tmp1[,i+2], nrow=m, ncol=n)
      # Check bounds
      obj$rate[[i]][obj$rate[[i]] < 0] <- NA
      obj$rate[[i]][obj$rate[[i]] > max.mx] <- max.mx
      if(scale > 1)
        obj$rate[[i]] <- obj$rate[[i]] / scale
      dimnames(obj$rate[[i]]) <- list(obj$age,obj$year)
    }
    names(obj$rate) = tolower(mnames)
  }

  if(mpopfile)
  {
    tmp2 <- utils::read.table(popfile,header=TRUE,na.strings=".",skip=popskip)
    obj$year=sort(unique(tmp2[,1]))
    n <- length(obj$year)
    m <- length(unique(tmp2[,2]))
    obj$age <- tmp2[1:m,2]
    pnames <- names(tmp2)[-c(1,2)]
    if(mfile)
    {
      if(sum(pnames==mnames) != length(pnames))
      {
        warning("Population names different from rates names")
        if(length(pnames) <- length(mnames))
          pnames <- mnames
      }
      if(n!=ncol(obj$rate[[1]]) | m != nrow(obj$rate[[1]]))
        warning("Population matrices different size from rates matrices")
    }
    p.mort <- length(pnames)
    obj$pop <- list()
    for(i in 1:p.mort)
    {
      obj$pop[[i]] = matrix(tmp2[,i+2], nrow=m, ncol=n)
      # Check bounds
      obj$pop[[i]][obj$pop[[i]] < 0] <- NA
      dimnames(obj$pop[[i]]) <- list(obj$age,obj$year)
    }
    names(obj$pop) = tolower(pnames)
  }

  junk <- options(warn=-1)
  obj$age <- as.numeric(as.character(obj$age))
  options(warn=junk$warn)
  if(is.na(obj$age[m]))
    obj$age[m] <- 2*obj$age[m-1] - obj$age[m-2]

  return(structure(obj,class="demogdata"))
}


#' Plot age-specific demographic functions
#'
#'
#' If \code{plot.type="functions"}, then years are plotted using a rainbow palette so the
#' earliest years are red, followed by orange, yellow, green, blue
#' and indigo with the most recent years plotted in violet.
#' If \code{plot.type="time"}, then each age is shown as a separate time series in a time plot.
#'
#' @param x Demogdata object such as created using \code{\link{read.demogdata}} or \code{\link{smooth.demogdata}}.
#' @param series Name of series to plot. Default: the first matrix within \code{datatype}.
#' @param datatype Name of demogdata object which contains series. Default \dQuote{rate}. Alternative: \dQuote{pop}.
#' @param years Vector indicating which years to plot. Default: all available years.
#' @param ages Vector indicating which ages to plot. Default: all available ages.
#' @param max.age Maximum age to plot. Default: all available ages.
#' @param transform Should a transformation of the data be plotted? Default is TRUE if the object contains mortality data and datatype="rate", and FALSE otherwise.
#' @param plot.type Type of plot: either \dQuote{functions} or \dQuote{time}.
#' @param type What type of plot should be drawn. See \code{\link[base]{plot}} for possible types.
#' @param main Main title for the plot.
#' @param xlab Label for x-axis.
#' @param ylab Label for y-axis.
#' @param pch Plotting character.
#' @param ... Other plotting parameters. In \code{points.demogdata}, all arguments are passed to \code{lines.demogdata}.
#'
#' @return None. Function produces a plot
#' @author Rob J Hyndman
#' @examples
#' plot(fr.mort)
#' par(mfrow=c(1,2))
#' plot(aus.fert,plot.type="time")
#' plot(aus.fert,plot.type="functions")
#' @keywords hplot
#' @export
plot.demogdata <- function(x, series=ifelse(!is.null(x$rate),names(x$rate)[1],names(x$pop)[1]),
    datatype=ifelse(!is.null(x$rate),"rate","pop"),
    years=x$year, ages=x$age, max.age=max(x$age), transform=(x$type=="mortality"),
    plot.type= c("functions", "time", "depth", "density"), type="l", main=NULL, xlab, ylab,...)
{
  plot.type <- match.arg(plot.type)
  series <- tolower(series)
  ages <- ages[ages <= max.age]
  data <- extract.ages(extract.years(x,years),ages,FALSE)
  if(x$type == "population")
    datatype <- "pop"

  # Extract data matrix
  if(!is.element(datatype,names(data)))
    stop(paste("Data type",datatype,"not found"))
  tmp <- data[[pmatch(datatype,names(data))]]
  y <- get.series(tmp,series)

  # Transformation
  if(data$lambda > 1-1e-6)
	 transform <- FALSE
  if(transform)
  {
    if(datatype=="rate")
      y <- BoxCox(y,data$lambda)
    else # Population
      y <- log(y)
    y[abs(y)==Inf] <- NA
  }

  # Not sure why we needed this, so removed for now
  # remove NA/NaN/Inf values
  #if(any(is.na(y)|is.nan(y)|is.infinite(y)))
  #{
  #  y <- apply(y, 2, na.interp)
  #}

  # Choose appropriate y axis label
  if(missing(ylab))
  {
    if(datatype=="pop")
      ylab <- "Population"
    else if(data$type=="mortality")
      ylab <- "Death rate"
    else if(data$type=="fertility")
      ylab <- "Fertility rate"
    else if(data$type=="migration")
      ylab <- "Net migration"
    else if(data$type=="population")
      ylab <- "Population"
    else
      stop("This shouldn't happen!")
    if(transform)
    {
      if(data$lambda==0 | datatype=="pop")
        ylab <- paste("Log",tolower(ylab))
      else
        ylab <- paste("Transformed",tolower(ylab))
    }
  }

  # Choose appropriate axis title
  if(is.null(main))
  {
    if(data$type=="fertility" & series=="female")
      main <- data$label
    else
      main <- paste(data$label,": ",series,sep="")
    if(datatype=="pop")
      main <- paste(main,"population")
    else if(data$type=="mortality")
      main <- paste(main,"death rates")
    else if(data$type=="fertility")
      main <- paste(main,"fertility rates")
    else if(data$type=="migration")
      main <- paste(main,"net migration")
    if(length(years)==1)
      main <- paste(main,"  (",years,")",sep="")
    else
      main <- paste(main,"  (",min(years),"-",max(years),")",sep="")
  }

  # Produce plot
  if(length(data$age)==1)
  {
    main <- paste(main,"  Age:",data$age)
    if(missing(xlab))
      xlab <- "Year"
  }
  else if(missing(xlab))
    xlab <- "Age"

  plot(fts(data$age,y,start=years[1],frequency=1,yname="",xname=""),plot.type=plot.type,xlab=xlab, ylab=ylab,main=main,type=type,...)
}

#' @rdname plot.demogdata
#' @export
lines.demogdata <- function(x, series=ifelse(!is.null(x$rate),names(x$rate)[1],names(x$pop)[1]),
    datatype=ifelse(!is.null(x$rate),"rate",""),
    years=x$year, ages=x$age, max.age=max(x$age), transform=(x$type=="mortality"),
    plot.type= c("functions", "time", "depth", "density"), ...)
{
  plot.type <- match.arg(plot.type)

  series <- tolower(series)
  ages <- ages[ages <= max.age]
  data <- extract.ages(extract.years(x,years),ages,FALSE)

  if(x$type == "population")
    datatype <- "pop"

  # Extract data matrix
  if(!is.element(datatype,names(data)))
    stop(paste("Data type",datatype,"not found"))
  tmp <- data[[pmatch(datatype,names(data))]]
  y <- get.series(tmp,series)

  # Transformation
  if(transform)
  {
    if(datatype=="rate")
      y <- BoxCox(y,data$lambda)
    else # Population
      y <- log(y)
}

  # Set other arguments to appropriate values
  lines(fts(data$age,y,start=years[1],frequency=1),plot.type=plot.type, ...)
}

#' @rdname plot.demogdata
#' @export
points.demogdata <- function(...,pch=1)
{
  lines.demogdata(...,type="p",pch=pch)
}
#' @export
print.demogdata <- function(x,...)
{
  Type <- x$type
  substr(Type,1,1) <- toupper(substr(Type,1,1))
  cat(paste(Type,"data for",x$label))
  cat("\n    Series: ")
  if(!is.null(x$rate))
    cat(names(x$rate))
  else
    cat(names(x$pop))
  cat(paste("\n    Years:",min(x$year),"-",max(x$year)))
  minx <- ifelse(min(x$age) < 0,"B",min(x$age))
  cat(paste("\n    Ages: ",minx,"-",max(x$age),"\n"))
}
#' @export
summary.demogdata <- function(object, ...)
{
	print(object)
}


#' Extract some years from a demogdata object
#'
#' Creates subset of demogdata object.
#'
#' @param data Demogdata object such as created using \code{\link{read.demogdata}} or \code{\link{smooth.demogdata}}.
#' @param years Vector of years to extract from data.
#'
#' @return Demogdata object with same components as \code{data} but with a subset of years.
#'
#' @author Rob J Hyndman
#' @examples
#' france.1918 <- extract.years(fr.mort,1918)
#' @keywords manip
#' @export
extract.years <- function(data,years)
{
  idx <- match(years,data$year)
  idx <- idx[!is.na(idx)]
	if(length(idx)==0)
		stop("No data available for those years")
  no.pop <- is.null(data$pop)
  no.rate <- is.null(data$rate)
  if(!no.rate)
  {
    nn <- length(data$rate)
    dname <- dimnames(data$rate[[1]])
    dname <- list(dname[[1]],dname[[2]][idx])
  }
  else if(!no.pop)
  {
    nn <- length(data$pop)
    dname <- dimnames(data$pop[[1]])
    dname <- list(dname[[1]],dname[[2]][idx])
  }
  else
    stop("No data!")

  for(j in 1:nn)
  {
    if(!no.rate)
    {
      data$rate[[j]] <- matrix(data$rate[[j]][,idx],ncol=length(idx))
      dimnames(data$rate[[j]]) <- dname
    }
    if(!no.pop)
    {
      data$pop[[j]] <- matrix(data$pop[[j]][,idx],ncol=length(idx))
      dimnames(data$pop[[j]]) <- dname
    }
    if(!is.null(data$obs.var))
    {
      if(length(dim(data$obs.var[[j]]))==2)
        data$obs.var[[j]] <- data$obs.var[[j]][,idx]
      else
        data$obs.var[[j]] <- data$obs.var[[j]][idx]
    }
  }
  data$year <- data$year[idx]
  return(data)
}



#' Extract some ages from a demogdata object
#'
#' Creates subset of demogdata object.
#'
#' @param data Demogdata object such as created using \code{\link{read.demogdata}} or \code{\link{smooth.demogdata}}.
#' @param ages Vector of ages to extract from data.
#' @param combine.upper If TRUE, ages beyond the maximum of \code{ages} are combined into the upper age group.
#' @return Demogdata object with same components as \code{data} but with a subset of ages.
#' @author Rob J Hyndman
#' @examples
#' france.teens <- extract.ages(fr.mort,13:19,FALSE)
#' plot(france.teens)
#' @keywords manip
#' @export
extract.ages <- function(data,ages,combine.upper=TRUE)
{
	no.pop <- is.null(data$pop)
  no.rate <- is.null(data$rate)
  if(combine.upper)
	{
		if(no.pop)
    {
      if(max(data$age) > max(ages))
        warning("No population data available for combining upper ages")
    }
		else
			data <- set.upperage(data,max(ages))
	}
  idx <- match(ages,data$age)
  idx <- idx[!is.na(idx)]
  no.pop <- is.null(data$pop)
  no.rate <- is.null(data$rate)
  if(!no.rate)
  {
    nn <- length(data$rate)
    dname <- dimnames(data$rate[[1]])
    dname <- list(dname[[1]][idx],dname[[2]])
  }
  else if(!no.pop)
  {
    nn <- length(data$pop)
    dname <- dimnames(data$pop[[1]])
    dname <- list(dname[[1]][idx],dname[[2]])
  }
  else
    stop("No data!")

  for(j in 1:nn)
  {
    if(!no.rate)
    {
      data$rate[[j]] <- matrix(data$rate[[j]][idx,],nrow=length(idx))
      dimnames(data$rate[[j]]) <- dname
    }
    if(!no.pop)
    {
      data$pop[[j]] <- matrix(data$pop[[j]][idx,],nrow=length(idx))
      dimnames(data$pop[[j]]) <- dname
    }
    if(!is.null(data$obs.var))
    {
      if(length(dim(data$obs.var[[j]]))==2)
        data$obs.var[[j]] <- data$obs.var[[j]][idx,]
      else
        data$obs.var[[j]] <- data$obs.var[[j]][idx]
    }
  }
  data$age <- data$age[idx]

  return(data)
}

#' Combine the upperages of a demogdata object.
#'
#' Computes demographic rates by combining age groups.
#'
#' @param data Demogdata object such as created using \code{\link{read.demogdata}} or \code{\link{smooth.demogdata}}.
#' @param max.age Upper age group. Ages beyond this are combined into the upper age group.
#'
#' @return Demogdata object with same components as \code{data} but with a subset of ages.
#' @author Rob J Hyndman
#' @examples
#' france.short <- set.upperage(fr.mort, 85)
#' @keywords manip
#' @export
set.upperage <- function(data, max.age)
{
  if(max(data$age) < max.age)
    stop("max.age too large")
  else if(max(data$age) == max.age)
    return(data)

  if(is.null(data$pop))
    stop("This procedure needs the population data")

  no.rate <- is.null(data$rate)
  multiple.series <- is.list(data$pop)

  if(multiple.series)
  {
    nn <- length(data$pop)
    fred <- data
    for(j in 1:nn)
    {
      if(!no.rate)
        fred$rate <- data$rate[[j]]
      fred$pop <- data$pop[[j]]
      tmp <- set.upperage(fred,max.age)
      if(!no.rate)
        data$rate[[j]] <- tmp$rate
      data$pop[[j]] <- tmp$pop
    }
    data$age <- tmp$age
  }
  else
  {
    idx <- data$age >= max.age
    if(sum(!idx) > 0)
    {
		  age <- data$age[!idx]
		  rnames <- rownames(data$pop)[!idx]
		  if(max(age) < max.age)
		  {
        age <- c(age,max.age)
        rnames <- c(rnames,NA)
		  }
    }
    else
    {
		  age <- max.age
		  rnames <- ""
    }
    rnames[length(rnames)] <- paste(max.age,"+",sep="")
    upper.pop <- data$pop[idx,]
    pop <- apply(matrix(upper.pop,nrow=sum(idx)),2,sum,na.rm=TRUE)
    data$pop <- rbind(matrix(data$pop[!idx,],ncol=ncol(data$pop)),pop)
    rownames(data$pop) <- rnames
    colnames(data$pop) <- data$year
    if(!no.rate)
    {
      upper.rate <- data$rate[idx,]
      actuals = apply(matrix(upper.rate*upper.pop,nrow=sum(idx)),2,sum,na.rm=TRUE)
      data$rate <- rbind(matrix(data$rate[!idx,],ncol=ncol(data$rate)),actuals/pop)
      rownames(data$rate) <- rnames
      colnames(data$rate) <- data$year
    }
    data$age <- age
    if(is.null(data$obs.var))
    {
      if(length(dim(data$obs.var))==2)
        data$obs.var <- data$obs.var[data$age <= max.age,]
      else
        data$obs.var <- data$obs.var[data$age <= max.age]
    }
  }
  return(data)
}


# Case insensitive version of is.element
is.el <- function(el,set)
{
  is.element(toupper(el),toupper(set))
}

get.series <- function(data,series)
{
  if(!is.el(series,names(data)))
    stop(paste("Series",series,"not found"))
  i <- match(toupper(series),toupper(names(data)))
  return(as.matrix(data[[i]]))
}



#' Combine two demogdata objects into one demogdata object
#'
#' Function to combine demogdata objects containing
#' different years but the same age structure into one demogdata
#' object. The standard use for this function will be combining
#' historical data with forecasts. The objects must be of the same type.
#'
#' @param obj1 First demogdata object (e.g., historical data).
#' @param obj2 Second demogdata object (e.g., forecasts).
#'
#' @return Object of class \dQuote{demogdata} with the following components:
#' \item{year}{Vector of years}
#' \item{age}{Vector of ages}
#' \item{rate}{Matrix of rates with with one age group per row and one column per year.}
#' \item{pop}{Matrix of populations in same form as \code{rate} and containing population numbers. This is only
#'   produced when both objects contain a \code{pop} component.}
#' \item{type}{Type of object: \dQuote{mortality}, \dQuote{fertility} or \dQuote{migration}.}
#' \item{label}{Name of area from which the data are taken.}
#'
#' @seealso \code{\link{demogdata}}
#' @author Rob J Hyndman
#' @examples
#' fit <- fdm(fr.mort)
#' fcast <- forecast(fit, h=50)
#' france2 <- combine.demogdata(fr.mort,fcast)
#' plot(france2)
#' plot(life.expectancy(france2))
#' lines(rep(max(fr.mort$year)+0.5,2),c(0,100),lty=3)
#' @keywords manip
#' @export
combine.demogdata <- function(obj1, obj2)
{
  if(!is.element("demogdata",class(obj1))  | !is.element("demogdata",class(obj2)))
    stop("Not demogdata objects")
  if(obj1$type != obj2$type)
    stop("Objects not of the same type")
  if(min(obj1$year) > min(obj2$year))
  {
    tmp <- obj2
    obj2 <- obj1
    obj1 <- tmp
  }
  if(max(obj1$year) > min(obj2$year))
    stop("Years overlap")

  pop <- (!is.null(obj1$pop) & !is.null(obj2$pop))
  idx <- match(obj1$age,obj2$age)
  idx <- idx[!is.na(idx)]
  m.obj <- list(label=obj1$label,age=obj1$age[idx],year=c(obj1$year,obj2$year),rate=list(),lambda = obj1$lambda)
  if(pop)
    m.obj$pop = list()
  if(!is.list(obj1$rate))
  {
    obj1$rate <- list(obj1$rate)
    if(pop)
      obj1$pop <- list(obj1$pop)
  }
  nn <- length(obj1$rate)
  k <- 0
  for (j in 1:nn)
  {
    i <- match(names(obj1$rate)[j], names(obj2$rate))
    if(is.na(i))
      i <- match(names(obj1$rate)[j], obj2$name)
    if(!is.na(i))
    {
      k <- k + 1
      m.obj$rate[[k]] <- cbind(obj1$rate[[j]][idx,],obj2$rate[[i]][idx,])
      names(m.obj$rate)[k] <- names(obj1$rate)[j]
      colnames(m.obj$rate[[k]]) <- m.obj$year
      rownames(m.obj$rate[[k]]) <- m.obj$age
      if(pop)
      {
        m.obj$pop[[k]] <- cbind(obj1$pop[[j]][idx,],obj2$pop[[i]][idx,])
        names(m.obj$pop)[k] <- names(obj1$pop)[j]
        colnames(m.obj$pop[[k]]) <- m.obj$year
        rownames(m.obj$pop[[k]]) <- m.obj$age
      }
    }
  }
  m.obj$type <- obj1$type
  return(structure(m.obj,class="demogdata"))
}



#' Mean and median functions for data of class demogdata
#'
#' Computes mean or median of demographic rates for each age level.
#'
#' @param x Demogdata object such as created using \code{\link{read.demogdata}} or \code{\link{smooth.demogdata}}.
#' @param series Name of demogdata series to plot..
#' @param transform Should transform of data be taken first?
#' @param na.rm a logical value indicating whether NA values should be stripped before the computation proceeds.
#' @param method Method for computing the median. Either "coordinate" for a coordinate-wise median, or "hossjercroux" for the
#'   L1-median using the Hossjer-Croux algorithm.
#' @param ... Other arguments.
#'
#' @return A list containing \code{x}=ages and \code{y}=mean or median rates.
#' @author Rob J Hyndman
#' @references
#' Hossjer, O., and Croux, C. (1995) Generalized univariate signed rank statistics for testing
#' and estimating a multivariate location parameter. \emph{Nonparametric Statistics}, \bold{4}, 293-308.
#'
#' @examples
#' plot(fr.mort)
#' lines(mean(fr.mort),lwd=2)
#' lines(median(fr.mort),lwd=2,col=2)
#' @keywords models
#' @export
mean.demogdata <- function(x,series=names(x$rate)[1],transform=TRUE,na.rm=TRUE,...)
{
  mx <- get.series(x$rate,series)
  # Transformation
  if(transform)
    mx <- BoxCox(mx,x$lambda)

  mx[mx < -1e9] <- NA
  loc <- rowMeans(mx,na.rm=na.rm)
  return(list(x=x$age,y=loc))
}

#' @rdname mean.demogdata
#' @export
median.demogdata <- function(x,  na.rm=FALSE, series=names(x$rate)[1],
    transform=TRUE,method=c("hossjercroux","coordinate"),...)
{
  method = match.arg(method)
  mx <- get.series(x$rate,series)
  if(transform)
    mx <- BoxCox(mx,x$lambda)
  mx[mx < -1e9] <- NA
  loc <- L1median(t(mx),method=method)
  return(list(x=x$age,y=loc))
}


# Sex ratios
#' Compute sex ratios from mortality rates
#'
#' Calculates the Male/Female ratios from historical or forecasted mortality rates.
#'
#' @param data Demogdata object of type \dQuote{mortality} such as obtained from \code{\link{read.demogdata}},
#' or an object of class \code{fmforecast} such as the output from  \code{\link{forecast.fdm}} or \code{\link{forecast.lca}}.
#'
#' @return Functional time series of sex ratios.
#'
#' @author Rob J Hyndman
#' @examples plot(sex.ratio(fr.mort),ylab="Sex ratios (M/F)")
#' @keywords models
#' @export
sex.ratio <- function(data)
{
  if (inherits(data, "demogdata"))
    rate.sr <- fts(x=data$age,y=data$rate$male/data$rate$female,start=min(data$year),frequency=1,xname="Age",yname="Sex ratio (M/F)")
  else if(inherits(data, "fmforecast2"))
    rate.sr <- fts(x=data$male$age,data$male$rate$male/data$female$rate$female,start=min(data$male$year),frequency=1,xname="Age",yname="Sex ratio (M/F)")
  else
    stop("Unknown class of data")
  return(rate.sr)
}
