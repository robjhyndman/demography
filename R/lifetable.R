# Lifetable functions
# Produce lifetable from mortality rates

#' Construct lifetables from mortality rates
#' 
#' Computes period and cohort lifetables from mortality rates for multiple years.
#' 
#' For period lifetables, all years and all ages specified are included in the tables. For cohort lifetables, 
#' if \code{ages} takes a scalar value, then the cohorts are taken to be of that age in each year contained in \code{years}. 
#' But if \code{ages} is a vector of values, then the cohorts are taken to be of those ages in the first year contained in \code{years}.
#' 
#' For example, if \code{ages=0} then lifetables of the birth cohorts for all years in \code{years} are computed. On the other hand, 
#' if \code{ages=0:100} and \code{years=1950:2010}, then lifetables of each age cohort in 1950 are computed.
#' 
#' In all cases, \eqn{q_x = m_x/(1+[(1-a_x)m_x])}{qx = mx/(1 + ((1-ax) * mx))} as per Chiang (1984).
#' 
#' Warning: the code has only been tested for data based on single-year age groups.
#' 
#' @param data Demogdata object such as obtained from \code{\link{read.demogdata}}, 
#' \code{\link{forecast.fdm}} or \code{\link{forecast.lca}}.
#' @param series Name of series to use.  Default is the first series in data\$rate.
#' @param years Vector indicating which years to include in the tables.
#' @param ages Vector indicating which ages to include in table.
#' @param max.age Age for last row. Ages beyond this are combined.
#' @param type Type of lifetable: \code{period} or \code{cohort}.
#' 
#' @return Object of class \dQuote{lifetable} containing the following components:
#' \item{label}{Name of region from which data are taken.}
#' \item{series}{Name of series}
#' \item{age}{Ages for lifetable}
#' \item{year}{Period years or cohort years}
#' \item{mx}{Death rate at age x.}
#' \item{qx}{The probability that an individual of exact age x will die before exact age x+1.}
#' \item{lx}{Number of survivors to exact age x.  The radix is 1.}
#' \item{dx}{The number of deaths between exact ages x and x+1.}
#' \item{Lx}{Number of years lived between exact age x and exact age x+1.}
#' \item{Tx}{Number of years lived after exact age x.}
#' \item{ex}{Remaining life expectancy at exact age x.}
#' Note that the lifetables themselves are not returned, only their components. However, there is a print method that constructs (and returns)
#' the lifetables from the above components.
#' 
#' @seealso \code{\link{life.expectancy}}
#' @author Heather Booth, Leonie Tickle, Rob J Hyndman, John Maindonald and Timothy Miller
#' @references 
#' Chiang CL. (1984) \emph{The life table and its applications}. Robert E Krieger Publishing Company: Malabar.
#' 
#' Keyfitz, N, and Caswell, H. (2005) \emph{Applied mathematical demography}, Springer-Verlag: New York.
#' 
#' Preston, S.H., Heuveline, P., and Guillot, M. (2001) \emph{Demography: measuring and modeling population processes}. Blackwell
#' 
#' @examples 
#' france.lt <- lifetable(fr.mort)
#' plot(france.lt)
#' lt1990 <- print(lifetable(fr.mort,year=1990))
#' 
#' france.LC <- lca(fr.mort)
#' france.fcast <- forecast(france.LC)
#' france.lt.f <- lifetable(france.fcast)
#' plot(france.lt.f)
#' 
#' # Birth cohort lifetables, 1900-1910
#' france.clt <- lifetable(fr.mort,type="cohort",age=0, years=1900:1910)
#' 
#' # Partial cohort lifetables for 1950
#' lifetable(fr.mort,type="cohort",years=1950)
#' @keywords models
#' @export
lifetable <- function(data, series=names(data$rate)[1], years=data$year, ages=data$age,
  max.age=min(100,max(data$age)), type=c("period","cohort"))
{
  if(!is.element("demogdata",class(data)))
    stop("data must be a demogdata object")
  if(data$type != "mortality")
    stop("data must be a mortality object")
  type <- match.arg(type)
  if(!is.el(series,names(data$rate)))
    stop(paste("Series",series,"not found"))
  if(is.na(sum(match(years,data$year))))
    stop("Years not present in data")
  sex <- series
  if(sex!="female" & sex!="male" & sex!="total")
  {
    if(is.element("model",names(data)))
      sex <- names(data$model)[4]
  }
  na <- length(data$age)
  if(na > 4)
    agegroup <- data$age[na-1]-data$age[na-2]
  else
    stop("Insufficient age groups")

  if(type=="period")
  {
    max.age <- min(max.age,max(ages))
    data <- extract.ages(data,ages,combine.upper=FALSE)
    if(max.age < max(ages))
      data <- extract.ages(data,min(ages):max.age,combine.upper=TRUE)
    data <- extract.years(data,years=years)
    mx <- get.series(data$rate,series)
    n <- length(years)
    p <- nrow(mx)
    rownames(mx) <- ages <- data$age
    colnames(mx) <- years
    qx <- lx <- dx <- Lx <- Tx <- ex <- rx <- mx*NA
    rownames(rx) <- ages - 1
    rownames(rx)[ages==0] <- "B"
    for(i in 1:n)
    {
      ltable <- lt(mx[,i],min(ages),agegroup,sex)
      nn <- length(ltable$qx)
      qx[1:nn,i] <- ltable$qx
      lx[1:nn,i] <- ltable$lx
      dx[1:nn,i] <- ltable$dx
      Lx[1:nn,i] <- ltable$Lx
      Tx[1:nn,i] <- ltable$Tx
      ex[1:nn,i] <- ltable$ex
      rx[1:nn,i] <- ltable$rx
    }
  }
  else if(type=="cohort" & length(ages)>1) # multiple ages, single year.
  {
    data <- extract.ages(data,min(ages):max.age,combine.upper=TRUE)
    data <- extract.years(data,years=seq(min(years),max(data$year),by=1))
    n <- length(data$year)
    p <- length(data$age)
    cmx <- matrix(NA,p,p)
    rownames(cmx) <- data$age
    colnames(cmx) <- paste(min(years)," age ",data$age,sep="")
    qx <- dx <- Tx <- lx <- Lx <- ex <- cmx
    cohort <- match(ages,data$age)
    cohort <- cohort[!is.na(cohort)]
    if(length(cohort)==0)
      stop("No data available")
    if(min(data$age[cohort]+n-1) < max.age)
      warning("Insufficient data for other lifetables. Try reducing max.age")
    for(coh in cohort)
    {
      if(data$age[coh]+n-1 > max.age)
      {
        subdata <- extract.ages(data,data$age[coh]:max.age,combine.upper=TRUE)
        mx <- get.series(subdata$rate,series)
        p <- nrow(mx)
        for (j in 1:p)
          cmx[coh+j-1,coh] <- mx[j,j]
        ltable <- lt(cmx[coh+(1:p)-1,coh], data$age[coh], agegroup, sex=sex)
        p <- length(ltable$lx)
        lx[coh+(1:p)-1,coh] <- ltable$lx
        Lx[coh+(1:p)-1,coh] <- ltable$Lx
        ex[coh+(1:p)-1,coh] <- ltable$ex
        qx[coh+(1:p)-1,coh] <- ltable$qx
        dx[coh+(1:p)-1,coh] <- ltable$dx
        Tx[coh+(1:p)-1,coh] <- ltable$Tx
      }
    }
    mx <- cmx
    # Retain columns in required cohort
    mx <- mx[,cohort,drop=FALSE]
    lx <- lx[,cohort,drop=FALSE]
    Lx <- Lx[,cohort,drop=FALSE]
    ex <- ex[,cohort,drop=FALSE]
    qx <- qx[,cohort,drop=FALSE]
    dx <- dx[,cohort,drop=FALSE]
    Tx <- Tx[,cohort,drop=FALSE]
    rx <- NULL
  }
  else #single age, multiple years.
  {
    data <- extract.years(data,years=seq(min(years),max(data$year),by=1))
    data <- extract.ages(data,ages:max.age,combine.upper=TRUE)
    n <- length(data$year)
    p <- length(data$age)
    ny <- length(years)
    cmx <- matrix(NA,p,ny)
    rownames(cmx) <- data$age
    colnames(cmx) <- paste(years," age ",ages,sep="")
    qx <- dx <- Tx <- lx <- Lx <- ex <- cmx
    minage <- max.age
    for(i in 1:ny)
    {
      subdata <- extract.years(data,years=seq(years[i],max(data$year),by=1))
      upperage <- min(ages+length(subdata$year)-1)
      minage <- min(minage,upperage)
      if(upperage >= max.age)
      {
        mx <- get.series(subdata$rate,series)
        p <- nrow(mx)
        for (j in 1:p)
          cmx[j,i] <- mx[j,j]
        ltable <- lt(cmx[,i],ages, agegroup, sex=sex)
        p <- length(ltable$lx)
        lx[1:p,i] <- ltable$lx
        Lx[1:p,i] <- ltable$Lx
        ex[1:p,i] <- ltable$ex
        qx[1:p,i] <- ltable$qx
        dx[1:p,i] <- ltable$dx
        Tx[1:p,i] <- ltable$Tx
      }
    }
    mx <- cmx
    rx <- NULL
  }

  return(structure(list(age=ages,year=years, mx=mx,qx=qx,lx=lx,dx=dx,Lx=Lx,Tx=Tx,ex=ex,rx=rx,
        series=series, type=type, label=data$label),class="lifetable"))
}



lt <- function (mx, startage = 0, agegroup = 5, sex)
{
  # Omit missing ages
  if (is.na(mx[1]))
    mx[1] <- 0
  firstmiss <- (1:length(mx))[is.na(mx)][1]
  if (!is.na(firstmiss))
    mx <- mx[1:(firstmiss - 1)]
  nn <- length(mx)
  if (nn < 1)
    stop("Not enough data to proceed")

  # Compute width of each age group
  if (agegroup == 1)
    nx <- c(rep(1, nn - 1), Inf)
  else if (agegroup == 5) # First age group 0, then 1-4, then 5-year groups.
    nx <- c(1, 4, rep(5, nn - 3), Inf)
  else stop("agegroup must be either 1 or 5")

  if (agegroup == 5 & startage > 0 & startage < 5)
    stop("0 < startage < 5 not supported for 5-year age groups")

  if (startage == 0) # for single year data and the first age (0) in 5-year data
  {
    if (sex == "female")
    {
      if (mx[1] < 0.107)
        a0 <- 0.053 + 2.8 * mx[1]
      else a0 <- 0.35
    }
    else if (sex == "male")
    {
      if (mx[1] < 0.107)
        a0 <- 0.045 + 2.684 * mx[1]
      else a0 <- 0.33
    }
    else # if(sex == "total")
    {
      if (mx[1] < 0.107)
        a0 <- 0.049 + 2.742 * mx[1]
      else a0 <- 0.34
    }
  }
  else if (startage > 0)
    a0 <- 0.5
  else
    stop("startage must be non-negative")
  if (agegroup == 1)
  {
    if (nn > 1)
      ax <- c(a0, rep(0.5, nn - 2), Inf)
    else
      ax <- Inf
  }
  else if (agegroup == 5 & startage == 0)
  {
    if (sex == "female")
    {
      if (mx[1] < 0.107)
        a1 <- 1.522 - 1.518 * mx[1]
      else
        a1 <- 1.361
    }
    else if (sex == "male")
    {
      if (mx[1] < 0.107)
        a1 <- 1.651 - 2.816 * mx[1]
      else
        a1 <- 1.352
    }
    else # sex == "total"
    {
      if (mx[1] < 0.107)
        a1 <- 1.5865 - 2.167 * mx[1]
      else a1 <- 1.3565
    }
    ax <- c(a0, a1, rep(2.6, nn - 3), Inf)
    ### ax=2.5 known to be too low esp at low levels of mortality
  }
  else # agegroup==5 and startage > 0
  {
    ax <- c(rep(2.6, nn - 1), Inf)
    nx <- c(rep(5, nn))
  }
  qx <- nx * mx/(1 + (nx - ax) * mx)
  # age <- startage + cumsum(nx) - 1
  # if (max(age) >= 75) {
  #    idx <- (age >= 75)
  #   ax[idx] <- (1/mx + nx - nx/(1 - exp(-nx * mx)))[idx]
  #  qx[idx] <- 1 - exp(-nx * mx)[idx]
  #    }
  #qx[qx > 1] <- 1  ################  NOT NEEDED IN THEORY

  #plot(qx)  #### TO CHECK RESULT RE QX>1

  qx[nn] <- 1
  if (nn > 1)
  {
    lx <- c(1, cumprod(1 - qx[1:(nn - 1)]))
    dx <- -diff(c(lx, 0))
  }
  else
    lx <- dx <- 1
  Lx <- nx * lx - dx * (nx - ax)
  Lx[nn] <- lx[nn]/mx[nn]
  Tx <- rev(cumsum(rev(Lx)))
  ex <- Tx/lx
  if (nn > 2)
    rx <- c(Lx[1]/lx[1], Lx[2:(nn - 1)]/Lx[1:(nn - 2)], Tx[nn]/Tx[nn-1])
  else if (nn == 2)
    rx <- c(Lx[1]/lx[1], Tx[nn]/Tx[nn - 1])
  else
    rx <- c(Lx[1]/lx[1])
  if (agegroup == 5)
    rx <- c(0, (Lx[1] + Lx[2])/5 * lx[1], Lx[3]/(Lx[1]+Lx[2]),
            Lx[4:(nn - 1)]/Lx[3:(nn - 2)], Tx[nn]/Tx[nn-1])
  result <- data.frame(ax = ax, mx = mx, qx = qx, lx = lx,
      dx = dx, Lx = Lx, Tx = Tx, ex = ex, rx = rx, nx = nx)
  return(result)
}


#' Plot life expectancy from lifetable
#' 
#' plots life expectancy for each age and each year as functional time series.
#' 
#' @param x Output from \code{\link{lifetable}}.
#' @param years Years to plot. Default: all available years.
#' @param main Main title.
#' @param xlab Label for x-axis.
#' @param ylab Label for y-axis.
#' @param ... Additional arguments passed to \code{\link[rainbow]{plot.fds}}.
#' 
#' @seealso \code{\link{life.expectancy}}, \code{\link{lifetable}}.
#' @author Rob J Hyndman
#' @examples 
#' france.lt <- lifetable(fr.mort)
#' plot(france.lt)
#' 
#' france.LC <- lca(fr.mort)
#' france.fcast <- forecast(france.LC)
#' france.lt.f <- lifetable(france.fcast)
#' plot(france.lt.f,years=2010)
#' @export
#' @keywords models
plot.lifetable <- function(x,years=x$year,main,xlab="Age",ylab="Expected number of years left",...)
{
  if(x$type != "period")
    stop("Currently only period lifetables can be plotted.")
  # Extract years
  idx <- match(years,x$year)
  idx <- idx[!is.na(idx)]
  idx <- idx[idx <= ncol(x$ex)]
  if(length(idx)==0)
    stop("Year not available")
  years <- x$year[idx]
  ny <- length(years)

  if(missing(main))
  {
    main <- paste("Life expectancy:",x$label,x$series)
    if(ny>1)
      main <- paste(main," (",min(years),"-",max(years),")",sep="")
    else
      main <- paste(main," (",years,")",sep="")
  }

  plot(fts(x$age,x$ex[,idx],start=years[1],frequency=1),main=main,ylab=ylab,xlab=xlab,...)
}

#' @rdname plot.lifetable
#' @export
lines.lifetable <- function(x,years=x$year,...)
{
  if(x$type != "period")
    stop("Currently only period lifetables can be plotted.")
  # Extract years
  idx <- match(years,x$year)
  idx <- idx[!is.na(idx)]
  idx <- idx[idx <= ncol(x$ex)]
  if(length(idx)==0)
    stop("Year not available")
  years <- x$year[idx]
  ny <- length(years)

  lines(fts(x$age,x$ex[,idx],start=x$year[1],frequency=1),...)
}


#' @export
print.lifetable <- function(x,digits=4,...)
{
  ny <- ncol(x$ex)
  outlist <- vector(length=ny,mode="list")
  for(i in 1:ny)
  {
    idx2 <- !is.na(x$mx[,i])
    if(sum(idx2)>0)
    {
      outlist[[i]] <- data.frame(x$mx[,i],x$qx[,i],x$lx[,i],x$dx[,i],x$Lx[,i],x$Tx[,i],x$ex[,i])[idx2,]
      rownames(outlist[[i]]) <- rownames(x$ex)[idx2]
      colnames(outlist[[i]]) <- c("mx","qx","lx","dx","Lx","Tx","ex")
    }
  }
  if(x$type=="period")
  {
    names(outlist) = x$year
    cat("Period ")
  }
  else
  {
    names(outlist) <- colnames(x$ex)
    cat("Cohort ")
  }
  cat(paste("lifetable for",x$label,":",x$series,"\n\n"))
  for(i in 1:ny)
  {
    if(!is.null(outlist[[i]]))
    {
      if(x$type=="period")
        cat(paste("Year:",names(outlist)[i],"\n"))
      else
        cat(paste("Cohort:",names(outlist)[i],"\n"))
      print(round(outlist[[i]],digits=digits))
      cat("\n")
    }
  }
  invisible(outlist)
}


# Compute expected age from single year mortality rates
get.e0 <- function(x,agegroup,sex,startage=0)
{
  lt(x, startage, agegroup, sex)$ex[1]
}

# Compute expected ages for multiple years
#' Estimate life expectancy from mortality rates
#' 
#' All three functions estimate life expectancy from \code{lifetable}. 
#' The function \code{flife.expectancy} is primarily designed for forecast life expectancies and will optionally 
#' produce prediction intervals. Where appropriate, it will package the results as a forecast object
#' which makes it much easier to product nice plots of forecast life expectancies.
#' The \code{e0} function is a shorthand wrapper for \code{flife.expectancy} with \code{age=0}.
#' 
#' @return Time series of life expectancies (one per year), or a forecast object of life expectancies (one per year).
#'   
#' @param data Demogdata object of type \dQuote{mortality} such as obtained from \code{\link{read.demogdata}}, 
#' or an object of class \code{fmforecast} such as the output from  \code{\link{forecast.fdm}} or \code{\link{forecast.lca}}, 
#' or an object of class \code{fmforecast2} such as the output from \code{\link{forecast.fdmpr}}.
#' @param series Name of mortality series to use. Default is the first demogdata series in data.
#' @param years Vector indicating which years to use.
#' @param type Either \code{period} or \code{cohort}.
#' @param age Age at which life expectancy is to be calculated.
#' @param max.age Maximum age for life table calculation.
#' @param PI If TRUE, produce a prediction interval.
#' @param nsim Number of simulations to use when computing a prediction interval.
#' @param ... Other arguments passed to \code{simulate} when producing prediction intervals.
#' 
#' @seealso \code{\link{lifetable}}
#' @author Rob J Hyndman
#' @examples 
#' plot(life.expectancy(fr.mort),ylab="Life expectancy")
#' 
#' france.LC <- lca(fr.mort,adjust="e0",years=1950:1997)
#' france.fcast <- forecast(france.LC,jumpchoice="actual")
#' france.e0.f <- life.expectancy(france.fcast)
#' 
#' france.fdm <- fdm(extract.years(fr.mort,years=1950:2006))
#' france.fcast <- forecast(france.fdm)
#' \dontrun{
#'   e0.fcast <- e0(france.fcast,PI=TRUE,nsim=200)
#'   plot(e0.fcast)}
#' 
#' life.expectancy(fr.mort,type='cohort',age=50)
#' 
#' @keywords models
#' @export
life.expectancy <- function(data,series=names(data$rate)[1],years=data$year,
    type=c("period","cohort"), age=min(data$age), max.age=min(100,max(data$age)))
{
  type <- match.arg(type)
  if(!is.el(series,names(data$rate)))
    stop(paste("Series",series,"not found"))
  if(is.null(max.age))
    max.age <- min(100,max(data$age))
  if(age > max.age | age > max(data$age))
    stop("age is greater than maximum age")
  else if(age < min(data$age))
    stop("age is less than minimum age")
  if(type=="period")
    data.lt <- lifetable(data,series,years,type=type,max.age=max.age)$ex
  else
   data.lt <- lifetable(data,series,years,type=type,ages=age,max.age=max.age)$ex
  idx <- match(age,rownames(data.lt))
  #if(sum(is.na(data.lt[idx,]))>0
  # max(data.lt[idx,]) > 1e9)
  #    warning("Some missing or infinite values in the life table calculation.\n  These can probably be avoided by setting max.age to a lower value.")

  return(ts(data.lt[idx,],start=years[1],frequency=1))
}

#' @rdname life.expectancy
#' @export
flife.expectancy <- function(data, series=NULL, years=data$year,
    type=c("period","cohort"), age, max.age=NULL, PI=FALSE, nsim=500, ...)
{
  type <- match.arg(type)
  if(is.element("fmforecast",class(data)))
  {
    if(data$type != "mortality")
       stop("data not a mortality object")
    hdata <- list(year=data$model$year, age=data$model$age,
         type=data$type, label=data$model$label, lambda=data$lambda)
    hdata$rate <- list(data$model[[4]])
    if(min(hdata$rate[[1]],na.rm=TRUE) < 0 | !is.null(data$model$ratio)) # Transformed
      hdata$rate <- list(InvBoxCox(hdata$rate[[1]],data$lambda))
    if(type=="cohort")
    {
      hdata$year <- c(hdata$year,data$year)
      hdata$rate <- list(cbind(hdata$rate[[1]],data$rate[[1]]))
    }
    names(hdata$rate) <- names(data$model)[4]
    if(!is.null(data$model$pop))
    {
      hdata$pop = list(data$model$pop)
      names(hdata$pop) <- names(hdata$rate)
      if(type=="cohort") # Add bogus population for future years
      {
        n <- ncol(hdata$pop[[1]])
        h <- length(hdata$year)-n
        hdata$pop[[1]] <- cbind(hdata$pop[[1]],matrix(rep(hdata$pop[[1]][,n],h),nrow=nrow(hdata$pop[[1]]),ncol=h))
      }
    }
    class(hdata) <- "demogdata"
    # Fix missing values. Why are they there?
    hdata$rate[[1]][is.na(hdata$rate[[1]])] <- 1-1e-5
    if(is.null(max.age))
      max.age <- max(data$age)
    if(missing(age))
      age <- min(hdata$age)

    x <- stats::window(life.expectancy(hdata,type=type,age=age,max.age=max.age),end=max(data$model$year))
    xf <- life.expectancy(data,years=years,type=type,age=age,max.age=max.age)
    if(type=="cohort")
    {
      xf <- ts(c(stats::window(x,start=max(data$model$year)-max.age+age+1, extend=TRUE),xf),end=max(time(xf)))
      if(sum(!is.na(xf)) > 0)
        xf <- stats::na.omit(xf)
      else
        xf <- stop("Not enough data to continue")
      if(min(time(x)) > max(data$model$year)-max.age+age)
        x <- NULL#ts(NA,end=min(time(xf))-1)
      else
        x <- stats::window(x,end=max(data$model$year)-max.age+age)
    }

    out <- structure(list(x=x,mean=xf,method="FDM model"),class="forecast")
    if(is.element("lca",class(data$model)))
      out$method = "LC model"
    else if(!is.null(data$product))
      out$method = "Coherent FDM model"
    if(PI) # Compute prediction intervals
    {
      e0calc <- (!is.element("product",names(data$rate)) & !is.element("ratio",names(data$rate)) & is.null(data$model$ratio))
      if(is.null(data$product) & is.null(data$var) & is.null(data$kt.f))
        warning("Incomplete information. Possibly this is from a coherent\n  model and you need to pass the entire object.")
      else
      {
        sim <- simulate(data,nsim,adjust.model.var=FALSE,...)
        if(type=="cohort") # Add actual rates for first few years
        {
          usex <- length(x)*any(!is.na(x))
          ny <- length(data$model$year) - usex
          sim2 <- array(NA,c(dim(sim)[1],dim(sim)[2]+ny,dim(sim)[3]))
          sim2[,(ny+1):dim(sim2)[2],] <- sim
          hrates <- hdata$rate[[1]][,usex + (1:ny)]
          sim2[,1:ny,] <- array(rep(hrates,dim(sim)[2]),c(dim(sim)[1],ny,dim(sim)[3]))
          sim <- sim2
          rm(sim2)
        }
        if(e0calc)
        {
          e0sim <- matrix(NA,dim(sim)[2],dim(sim)[3])
          simdata <- data
          if(type=="cohort")
            simdata$year <- min(time(out$mean))-1 + 1:dim(sim)[2]
          for(i in 1:dim(sim)[3])
          {
            simdata$rate <- list(as.matrix(sim[,,i]))
            names(simdata$rate) <- names(data$rate)[1]
            e0sim[,i] <- life.expectancy(simdata,type=type,age=age,max.age=max.age)
          }
          e0sim <- e0sim[1:length(xf), , drop = FALSE]
          if(is.element("lca",class(data$model)))
            out$level <- data$kt.f$level
          else
            out$level <- data$coeff[[1]]$level
          out$lower <- ts(apply(e0sim,1,quantile,prob=0.5 - out$level/200,na.rm=TRUE))
          out$upper <- ts(apply(e0sim,1,quantile,prob=0.5 + out$level/200,na.rm=TRUE))
          stats::tsp(out$lower) <- stats::tsp(out$upper) <- stats::tsp(out$mean)
        }
        out$sim <- sim
      }
    }
    return(out)
  }
  else if(is.element("fmforecast2",class(data)))
  {
    if(data[[1]]$type != "mortality")
      stop("data not a mortality object")
    if(is.null(series))
      series <- names(data)[1]
    if(is.element("product",names(data))) # Assume coherent model
    {
      if(missing(age))
        age <- min(data[["product"]]$age)
      if(is.null(max.age))
        max.age <- max(data[["product"]]$age)
      if(series=="total")
      {
        # Compute forecast total mortality rates
        # Assume sex ratios of populations stay the same as last k years
        ny <- length(data[["product"]]$model$year)
        h <- length(data[["product"]]$year)
        k <- 5
        sex.ratio <- (data$female$model$pop / data$male$model$pop)[,ny-(1:k)+1]
        sex.ratio <- apply(sex.ratio, 1, median)
        frate <- data$female$rate$female
        mrate <- data$male$rate$male
        totalrate <- frate/(1+1/sex.ratio) + mrate/(1+sex.ratio)
        # Construct demogdata object with estimated total rates
        # Use last year of population as approximation
        total <- demogdata(totalrate,
                   pop=matrix(data$female$model$pop[,ny] +
                        data$male$model$pop[,ny],
                              nrow=max.age+1,ncol=h),
                   ages=data[["product"]]$age,
                   years=data[["product"]]$year,
                   type="mortality",
                   label=data$product$model$label,
                   name="total",
                   lambda=0)
        # Compute total life expectancy
        out <- flife.expectancy(total, PI=FALSE, age=age, max.age=max.age, type=type)
        if(PI)
        {
          # Create forecast object
          if(is.element("ts",class(out)))
          {
            out <- structure(list(mean=out, method="Coherent FDM model"), class='forecast')
            frate <- exp(data$female$model$female)
            mrate <- exp(data$male$model$male)
            htotalrate <- frate/(1+1/sex.ratio) + mrate/(1+sex.ratio)
            historictotal <- demogdata(htotalrate,
                   pop=data$female$model$pop + data$male$model$pop,
                   ages=data$female$model$age,
                   years=data$female$model$year,
                   type="mortality",
                   label=data$product$model$label,
                   name="total",
                   lambda=0)
            out$x <- life.expectancy(historictotal)
          }
          # Generate simulated products and ratios
          prodsim <- flife.expectancy(data$product,nsim=nsim,PI=TRUE,age=age,max.age=max.age,type=type)
          data$ratio[["female"]]$model$ratio <- TRUE # To avoid PI calculation on flife.expectancy
          ratiosim <- flife.expectancy(data$ratio[["female"]],nsim=nsim,PI=TRUE,age=age,max.age=max.age,type=type)
          # Simulated future rates
          fsim <- prodsim$sim * ratiosim$sim
          msim <- prodsim$sim / ratiosim$sim
          sim <- fsim/(1+1/sex.ratio) + msim/(1+sex.ratio)
          dimsim <- dim(sim)
          simdata <- total
          nnn <- dim(total$rate[[1]])
          useyears <- (1:nnn[2]) + (dimsim[2]-nnn[2])
          e0sim <- matrix(NA,nnn[2],dimsim[3])
          for(i in 1:dimsim[3])
          {
            simdata$rate <- list(as.matrix(sim[,useyears,i]))
            names(simdata$rate) <- names(total$rate)[1]
            fl <- flife.expectancy(simdata,type=type,age=age,max.age=max.age)
            if(class(fl)=="ts")
              e0sim[,i] <- fl
            else if(class(fl)=="forecast")
              e0sim[,i] <- fl$mean
            else
              stop("No idea what's going on here.")
          }
          out$level <- data$product$coeff[[1]]$level
          out$lower <- ts(apply(e0sim,1,quantile,prob=0.5 - out$level/200))
          out$upper <- ts(apply(e0sim,1,quantile,prob=0.5 + out$level/200))
          stats::tsp(out$lower) <- stats::tsp(out$upper) <- stats::tsp(out$mean)
        }
      }
      else
      {
        if(missing(age))
          age <- min(data[[series]]$age)
        if(is.null(max.age))
          max.age <- max(data[[series]]$age)
        out <- flife.expectancy(data[[series]],PI=FALSE,age=age,max.age=max.age,type=type)
        out$method <- "Coherent FDM model"
        if(PI)
        {
          # Generate simulated products and ratios
          prodsim <- flife.expectancy(data$product,nsim=nsim,PI=TRUE,age=age,max.age=max.age,type=type)
          data$ratio[[series]]$model$ratio <- TRUE # To avoid PI calculation on flife.expectancy
          ratiosim <- flife.expectancy(data$ratio[[series]],nsim=nsim,PI=TRUE,age=age,max.age=max.age,type=type)
          # Simulated future rates
          sim <- prodsim$sim * ratiosim$sim
          dimsim <- dim(sim)
          simdata <- data[[series]]
          nnn <- dim(data[[series]]$rate[[1]])
          useyears <- (1:nnn[2]) + (dimsim[2]-nnn[2])
          e0sim <- matrix(NA,nnn[2],dimsim[3])
          for(i in 1:dim(sim)[3])
          {
            simdata$rate <- list(as.matrix(sim[,useyears,i]))
            names(simdata$rate) <- names(data[[series]]$rate)[1]
            fl <- flife.expectancy(simdata,type=type,age=age,max.age=max.age)
            if(class(fl)=="ts")
              e0sim[,i] <- fl
            else if(class(fl)=="forecast")
              e0sim[,i] <- fl$mean
            else
              stop("No idea what's going on here.")
          }
          #browser()
          out$level <- data$product$coeff[[1]]$level
          out$lower <- ts(apply(e0sim,1,quantile,prob=0.5 - out$level/200))
          out$upper <- ts(apply(e0sim,1,quantile,prob=0.5 + out$level/200))
          stats::tsp(out$lower) <- stats::tsp(out$upper) <- stats::tsp(out$mean)
        }
      }
    }
    else
      out <- flife.expectancy(data[[series]],PI=PI,nsim=nsim,max.age=max.age,type=type,age=age)
    return(out)
  }
  else
  {
    if(!is.element("demogdata",class(data)))
      stop("data must be a demogdata object")
    if(data$type != "mortality")
      stop("data must be a mortality object")
    if(is.null(series))
      series <- names(data$rate)[1]
    if(missing(age))
      age <- min(data$age)
    return(life.expectancy(data,series=series,years=years,type=type,age=age,max.age=max.age))
  }
}

#' @rdname life.expectancy
#' @export
e0 <- function(data, series=NULL, years=data$year,
    type=c("period","cohort"), max.age=NULL,PI=FALSE, nsim=500, ...)
{
  flife.expectancy(data, series=series, years=years,age=0,
    type=type, max.age=max.age,PI=PI,nsim=nsim,...)
}
