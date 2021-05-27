## Last updated 29 May 2003 by RJH
## Added plot.lca and print.lca
## Changed way of limiting ages in lca
## Fixed bug which arose occasionally when using method "dt"
## Made "dt" the default as in original LC paper.

#' Model mortality or fertility data using Lee-Carter approach
#'
#' Lee-Carter model of mortality or fertility rates. \code{lca} produces a
#' standard Lee-Carter model by default, although many other options are
#' available. \code{bms} is a wrapper for \code{lca} and returns a model based
#' on the Booth-Maindonald-Smith methodology.
#'
#'
#' All mortality or fertility data are assumed to be in matrices of
#' mortality or fertility rates within \code{data$rate}. Each row is one age group
#' (assumed to be single years). Each column is one year. The
#' function produces a model for the \code{series} mortality or fertility rate matrix
#' within \code{data$rate}. Forecasts from this model can be obtained using \code{\link{forecast.lca}}.
#'
#' @param data demogdata object of type \dQuote{mortality} or
#'   \dQuote{fertility}. Output from read.demogdata.
#' @param series name of series within data containing mortality or fertility
#'   values (1x1)
#' @param years years to include in fit. Default: all available years.
#' @param ages ages to include in fit. Default: all available ages up to
#'   \code{max.age}.
#' @param max.age upper age to include in fit. Ages beyond this are collapsed
#'   into the upper age group.
#' @param adjust method to use for adjustment of coefficients \eqn{k_t kt}.
#'   Possibilities are \dQuote{dxt} (BMS method), \dQuote{dt} (Lee-Carter
#'   method), \dQuote{e0} (method based on life expectancy) and \dQuote{none}.
#'   Defaults are \dQuote{dxt} for \code{bms()} and \dQuote{dt} for
#'   \code{lca()}.
#' @param chooseperiod If TRUE, it will choose the best fitting period.
#' @param minperiod Minimum number of years to include in fitting period if
#'   chooseperiod=TRUE.
#' @param breakmethod method to use for identifying breakpoints if
#'   chooseperiod=TRUE. Possibilities are \dQuote{bai} (Bai's method computed
#'   using \code{\link[strucchange]{breakpoints}} in the strucchange package)
#'   and \dQuote{bms} (method based on mean deviance ratios described in BMS).
#' @param scale If TRUE, it will rescale bx and kt so that kt has drift
#'   parameter = 1.
#' @param restype method to use for calculating residuals. Possibilities are
#'   \dQuote{logrates}, \dQuote{rates} and \dQuote{deaths}.
#' @param interpolate If TRUE, it will estimate any zero mortality or fertility
#'   rates using the same age group from nearby years.
#'
#' @return Object of class \dQuote{lca} with the following components:
#' \item{label}{Name of region}
#' \item{age}{Ages from \code{data} object.}
#' \item{year}{Years from \code{data} object.}
#' \item{<series>}{Matrix of mortality or fertility data as contained in \code{data}. It takes the name given by the series argument.}
#' \item{ax}{Average deathrates across fitting period}
#' \item{bx}{First principal component in Lee-Carter model}
#' \item{kt}{Coefficient of first principal component}
#' \item{residuals}{Functional time series of residuals.}
#' \item{fitted}{Functional time series containing estimated mortality or fertility rates from model}
#' \item{varprop}{Proportion of variance explained by model.}
#' \item{y}{The data stored as a functional time series object.}
#' \item{mdev}{Mean deviance of total and base lack of fit, as described in Booth, Maindonald and Smith.}
#'
#' @references Booth, H., Maindonald, J., and Smith, L. (2002) Applying Lee-Carter
#' under conditions of variable mortality decline. \emph{Population Studies}, \bold{56}, 325-336.
#'
#' Lee, R.D., and Carter, L.R. (1992) Modeling and forecasting US mortality. \emph{Journal of
#'   the American Statistical Association}, \bold{87}, 659-671.
#'
#' @author Heather Booth, Leonie Tickle, John Maindonald and Rob J Hyndman.
#'
#' @seealso \code{\link{forecast.lca}}, \code{\link{fdm}}
#' @examples
#' \dontrun{
#' france.LC1 <- lca(fr.mort, adjust="e0")
#' plot(france.LC1)
#' par(mfrow=c(1,2))
#' plot(fr.mort,years=1953:2002,ylim=c(-11,1))
#' plot(forecast(france.LC1,jumpchoice="actual"),ylim=c(-11,1))
#'
#' france.bms <- bms(fr.mort, breakmethod="bai")
#' fcast.bms <- forecast(france.bms)
#' par(mfrow=c(1,1))
#' plot(fcast.bms$kt)
#' }
#' @keywords models
#' @export
lca <-  function(data,series=names(data$rate)[1],years=data$year, ages=data$age,
    max.age=100, adjust=c("dt","dxt","e0","none"),
    chooseperiod=FALSE, minperiod=20, breakmethod=c("bai","bms"),
    scale=FALSE, restype=c("logrates","rates","deaths"), interpolate=FALSE)
{
  if (class(data) != "demogdata") {
    stop("Not demography data")
  }
  if (!any(data$type == c("mortality", "fertility"))) {
    stop("Neither mortality nor fertility data")
  }

    adjust <- match.arg(adjust)
    restype <- match.arg(restype)
    breakmethod <- match.arg(breakmethod)

    data <- extract.ages(data,ages,combine.upper=FALSE)
    if(max.age < max(ages))
        data <- extract.ages(data,min(ages):max.age,combine.upper=TRUE)
    startage <- min(data$age)

    # Extract mortality rates and population numbers
    mx <- get.series(data$rate,series)
    pop <- get.series(data$pop,series)

    # Truncate years
    startyear <- min(years)
    stopyear <- max(years)
    if(startyear > max(data$year) | stopyear < min(data$year))
        stop("Year not found")
    startyear <- max(startyear,min(data$year))
    if(!is.null(stopyear))
        stopyear <- min(stopyear,max(data$year))
    else
        stopyear <- max(data$year)
    id2 <- stats::na.omit(match(startyear:stopyear,data$year))

    mx <- mx[,id2]
    pop <- pop[,id2]
    year <- data$year[id2]
    deltat <- year[2]-year[1]
    ages <- data$age
    n <- length(ages)
    m <- length(year)
    mx <- matrix(mx,nrow=n,ncol=m)
    colnames(mx) <- year
    rownames(mx) <- ages

    # Interpolate where rates are zero
    if(interpolate)
    {
        # Remove missing values
        mx[is.na(mx)] <- 0
        if(sum(abs(mx)<1e-9,na.rm=TRUE) > 0)
        {
            warning("Replacing zero values with estimates")
            for(i in 1:n)
                mx[i,] <- fill.zero(mx[i,])
        }
    }

    # Transpose data and get deaths and logrates
    mx <- t(mx)
    mx[mx==0] <- NA
    logrates <- log(mx)

    pop <- t(pop)
    deaths <- pop*mx

    # Do SVD
    ax <- apply(logrates,2,mean, na.rm=TRUE) # ax is mean of logrates by column
    if(sum(ax < -1e9) > 0)
        stop(sprintf("Some %s rates are zero.\n Try reducing the maximum age or setting interpolate=TRUE.", data$type))
    clogrates <- sweep(logrates,2,ax) # central log rates (with ax subtracted) (dimensions m*n)
    svd.mx <- svd(clogrates)

    # Extract first principal component
    sumv <- sum(svd.mx$v[,1])
    bx <- svd.mx$v[,1]/sumv
    kt <- svd.mx$d[1] * svd.mx$u[,1] * sumv

    # Adjust kt
    ktadj <- rep(0,m)
    logdeathsadj <- matrix(NA,n,m)
    z <- log(t(pop))+ax

    # Use regression to guess suitable range for root finding method
    x <- 1:m
    ktse <- stats::predict(stats::lm(kt ~ x),se.fit=TRUE)$se.fit
    ktse[is.na(ktse)] <- 1
    agegroup = ages[4]-ages[3]

    if (adjust=="dxt")
    {
        options(warn=-1) # Prevent warnings on non-integer population values
        for(i in 1:m)
        {
            y <- as.numeric(deaths[i,])
            zi <- as.numeric(z[,i])
            weight <- as.numeric(zi > -1e-8)  # Avoid -infinity due to zero population
            yearglm <- stats::glm(y ~ offset(zi)-1+bx, family=stats::poisson, weights=weight)
            ktadj[i] <- yearglm$coef[1]
            logdeathsadj[,i] <- z[,i]+bx*ktadj[i]
        }
        options(warn=0)
    }
    else if(adjust=="dt")
    {
        FUN <- function(p,Dt,bx,ax,popi) {Dt - sum(exp(p*bx+ax)*popi)}
        for(i in 1:m)
        {
            if(i==1)
                guess <- kt[1]
            else
                guess <- mean(c(ktadj[i-1],kt[i]))
            ktadj[i] <- findroot(FUN, guess=guess, margin = 3*ktse[i],ax=ax,bx=bx,popi=pop[i,],Dt=sum(as.numeric(deaths[i,])))
            logdeathsadj[,i] <- z[,i]+bx*ktadj[i]
        }
    }
    else if(adjust=="e0")
    {
        e0 <- apply(mx,1,get.e0,agegroup=agegroup,sex=series,startage=startage)
        FUN2 <- function(p,e0i,ax,bx,agegroup,series,startage){e0i - estimate.e0(p,ax,bx,agegroup,series,startage)}
        for (i in 1:m)
        {
            if(i==1)
                guess <- kt[1]
            else
                guess <- mean(c(ktadj[i-1],kt[i]))
            ktadj[i] <- findroot(FUN2, guess=guess, margin = 3*ktse[i],e0i=e0[i],ax=ax,bx=bx,agegroup=agegroup,series=series,startage=startage)
        }
    }
    else if(adjust=="none")
        ktadj <- kt
    else
        stop("Unknown adjustment method")

    kt <- ktadj

    # Find linear section of kt and refit
    if(chooseperiod)
    {
        if(breakmethod=="bai")
        {
            x <- 1:m
            # Find breakpoints
            bp <- strucchange::breakpoints(kt ~ x)$breakpoints
            # Omit breakpoints less than minperiod from end
            bp <- bp[bp <= (m-minperiod)]
            bestbreak <- max(bp)
            return(lca(data,series,year[(bestbreak+1):m],ages=ages,max.age=max.age,
                adjust=adjust,interpolate=interpolate,chooseperiod=FALSE,scale=scale))
        }
        else
        {
            RS <- devlin <- devadd <- numeric(m-2)
            for(i in 1:(m-2))
            {
                tmp <- lca(data,series,year[i:m],ages=ages,max.age=max.age,adjust=adjust,chooseperiod=FALSE,interpolate=interpolate,scale=scale)
                devlin[i] <- tmp$mdev[2]
                devadd[i] <- tmp$mdev[1]
                RS[i] <- (tmp$mdev[2]/tmp$mdev[1])
            }
            bestbreak <- order(RS[1:(m-minperiod)])[1]-1
            out <- lca(data,series,year[(bestbreak+1):m],ages=ages,max.age=max.age,
                adjust=adjust,chooseperiod=FALSE,interpolate=interpolate,scale=scale)
            out$mdevs <- ts(cbind(devlin,devadd,RS),start=startyear,deltat=deltat)
            dimnames(out$mdevs)[[2]] <- c("Mean deviance total","Mean deviance base","Mean deviance ratio")
            return(out)
        }
    }

    # Estimate rates from fitted values and get residuals
    logfit <- fitmx(kt,ax,bx,transform=TRUE)
    colnames(logfit) <- ages
    rownames(logfit) <- year
    if(restype=="logrates")
    {
        fit <- logfit
        res <- logrates - fit
    }
    else if(restype=="rates")
    {
        fit <- exp(logfit)
        res <- exp(logrates) - fit
    }
    else if(restype=="deaths")
    {
        fit <- exp(logfit)*pop
        res <- deaths - fit
    }
    residuals <- fts(ages,t(res),frequency=1/deltat,start=years[1],xname="Age",
                     yname=paste("Residuals", data$type, "rate"))
    fitted <- fts(ages,t(fit),frequency=1/deltat,start=years[1],xname="Age",
                  yname=paste("Fitted", data$type, "rate"))

    names(ax) <- names(bx) <- ages

    # Rescaling bx, kt
    if(scale)
    {
        avdiffk <- -mean(diff(kt))
        bx <- bx*avdiffk
        kt <- kt/avdiffk
    }

    # Compute deviances
    deathsadjfit <- exp(logfit)*pop
    drift <- mean(diff(kt))
    ktlinfit <- mean(kt) + drift* (1:m - (m+1)/2)
    deathslinfit <- fitmx(ktlinfit,ax,bx,transform=FALSE)*pop
    dflogadd <- (m-2)*(n-1)
    mdevlogadd <- 2/dflogadd*sum(deaths*log(deaths/deathsadjfit)-(deaths-deathsadjfit))
    dfloglin <- (m-2)*n
    mdevloglin <- 2/dfloglin*sum(deaths*log(deaths/deathslinfit)-(deaths-deathslinfit))
    mdev <- c(mdevlogadd,mdevloglin)

    #Return
    output <- list(label=data$label,age=ages,year=year, mx=t(mx),
        ax=ax, bx=bx, kt=ts(kt,start=startyear,deltat=deltat), residuals=residuals, fitted=fitted,
        varprop=svd.mx$d[1]^2/sum(svd.mx$d^2),
        y=fts(ages,t(mx),start=years[1],frequency=1/deltat,xname="Age",
              yname=ifelse(data$type == "mortality", "Mortality", "Fertility")),
        mdev=mdev)
    names(output)[4] <- series
    output$call <- match.call()
    names(output$mdev) <- c("Mean deviance base","Mean deviance total")
    output$adjust <- adjust
    output$type <- data$type
    return(structure(output,class="lca"))
}

#' @rdname lca
#' @export
bms <-  function(data,series=names(data$rate)[1],years=data$year, ages=data$age,
    max.age=100, minperiod=20, breakmethod=c("bms","bai"), scale=FALSE, restype=c("logrates","rates","deaths"),
    interpolate=FALSE)
{
    restype <- match.arg(restype)
    breakmethod <- match.arg(breakmethod)
    out <- lca(data,series=series, years=years, ages=ages, max.age=max.age, adjust="dxt",
        chooseperiod=TRUE, minperiod=minperiod, scale=scale, restype=restype, breakmethod=breakmethod,
        interpolate=interpolate)
    out$call <- match.call()
    return(out)
}


estimate.e0 <- function(kt,ax,bx,agegroup,series,startage=0)
{
    if(length(kt)>1)
        stop("Length of kt greater than 1")
    mx <- c(fitmx(kt,ax,bx))
    return(get.e0(mx,agegroup,series,startage=startage))
}

fitmx <- function (kt,ax,bx,transform=FALSE)
{
    # Derives mortality rates from kt mortality index,
    # per Lee-Carter method
    clogratesfit <- outer(kt, bx)
    logratesfit <- sweep(clogratesfit,2,ax,"+")
    if(transform)
        return(logratesfit)
    else
        return(exp(logratesfit))
}

#' @rdname plot.fmforecast
#' @export
plot.lca <- function(x,...)
{
    x$basis <- cbind(x$ax,x$bx)
    x$coeff <- cbind(rep(1,length(x$kt)),x$kt)
    colnames(x$basis) <- c("mean","bx")
    if(x$adjust != "none")
        xlab <- "kt (adjusted)"
    else
        xlab <- "kt"
    ftsa::plot.ftsm(x,1,xlab1="Age",ylab1="bx",xlab2="Year",ylab2=xlab,mean.lab="ax",...)
}

#' @export
print.lca <- function(x,...)
{
    cat("Lee-Carter analysis\n")
    cat(paste("\nCall:",deparse(x$call),"\n"))
    cat(paste("\nAdjustment method:",x$adjust))
    cat(paste("\nRegion:"),x$label)
    cat(paste("\nYears in fit:",min(x$year),"-",max(x$year)))
    cat(paste("\nAges in fit:",min(x$age),"-",max(x$age),"\n"))
    cat(paste("\nPercentage variation explained: ",round(x$varprop*100,1),"%\n",sep=""))
}

#' @rdname summary.fdm
#' @export
summary.lca <- function(object,...)
{
    print(object)

    cat(sprintf("\nERROR MEASURES BASED ON %s RATES\n", toupper(object$type)))
    printout(fdmMISE(object[[4]],exp(object$fitted$y),age=object$y$x,years=object$year))

    cat(sprintf("\nERROR MEASURES BASED ON LOG %s RATES\n", toupper(object$type)))
    printout(fdmMISE(log(object[[4]]),object$fitted$y,age=object$y$x,years=object$year))
}


printout <- function(output)
{
    junk1 <- cbind(output$ME,output$MSE,output$MPE,output$MAPE)
    rownames(junk1) <- output$age
    colnames(junk1) <- c("ME","MSE","MPE","MAPE")
    junk2 <- cbind(output$MIE,output$MISE,output$MIPE,output$MIAPE)
    rownames(junk2) = output$year
    colnames(junk2) = c("IE","ISE","IPE","IAPE")
    cat(paste("\nAverages across ages:\n"))
    print(round(apply(junk1,2,mean),5))
    cat(paste("\nAverages across years:\n"))
    print(round(apply(junk2,2,mean),5))
    cat("\n")
}

# Function performs predictions of k and life expectancy based on leecarter results (in lcaout)

#' Forecast demogdata data using Lee-Carter method.
#'
#' The kt coefficients are forecast using a random walk with drift.
#' The forecast coefficients are then multiplied by bx to
#' obtain a forecast demographic rate curve.
#'
#' @param object Output from \code{\link{lca}}.
#' @param h Number of years ahead to forecast.
#' @param se Method used for computation of standard error. Possibilities: \dQuote{innovdrift} (innovations and drift) and \dQuote{innovonly} (innovations only).
#' @param jumpchoice Method used for computation of jumpchoice. Possibilities: \dQuote{actual} (use actual rates from final year) and \dQuote{fit} (use fitted rates).
#' @param level Confidence level for prediction intervals.
#' @param ... Other arguments.
#'
#' @return Object of class \code{fmforecast} with the following components:
#' \item{label}{Region from which the data are taken.}
#' \item{age}{Ages from \code{object}.}
#' \item{year}{Years from \code{object}.}
#' \item{rate}{List of matrices containing forecasts, lower bound and upper bound of prediction intervals.
#'   Point forecast matrix takes the same name as the series that has been forecast.}
#' \item{fitted}{Matrix of one-step forecasts for historical data}
#' Other components included are
#' \item{e0}{Forecasts of life expectancies (including lower and upper bounds)}
#' \item{kt.f}{Forecasts of coefficients from the model.}
#' \item{type}{Data type.}
#' \item{model}{Details about the fitted model}
#'
#' @author Rob J Hyndman
#' @examples
#' france.lca <- lca(fr.mort, adjust="e0")
#' france.fcast <- forecast(france.lca, 50)
#' plot(france.fcast)
#' plot(france.fcast,'c')
#' @keywords models
#' @export
forecast.lca <- function(object, h=50, se=c("innovdrift","innovonly"), jumpchoice=c("fit","actual"), level=80, ...)
{
    se <- match.arg(se)
    jumpchoice <- match.arg(jumpchoice)

    # Section 1 Read in data from object
    jumpyear <- max(object$year)
    nyears <- length(object$year)
    nages <- length(object$age)

    # Find jumprates
    if(jumpchoice=="actual")
        jumprates <- object[[4]][,nyears]
    else if(jumpchoice=="fit")
        jumprates <- exp(object$ax + object$bx*object$kt[nyears])
    else
        stop(paste("Unknown jump choice:",jumpchoice))
    object$kt <- object$kt - object$kt[nyears]

    # Time series estimation of kt as Random walk with drift
    fit <- forecast::rwf(object$kt, drift=TRUE)
    kt.drift <- fit$model$par$drift
    sec <- fit$model$par$drift.se
    see <- sqrt(fit$model$sigma2)

    # Project kt
    x <- 1:h
    zval <- stats::qnorm(0.5 + 0.005*level)
    kt.forecast <- object$kt[nyears] + (x * kt.drift)

    # Calculate standard errors of forecast kt
    if (se=="innovdrift")
        kt.stderr <- sqrt(x*(see^2) + (x*sec)^2)
    else if(se=="innovonly")
        kt.stderr <- sqrt(x*(see^2))
    kt.lo.forecast <- kt.forecast - (zval*kt.stderr)
    kt.hi.forecast <- kt.forecast + (zval*kt.stderr)
    kt.f <- data.frame(kt.forecast,kt.lo.forecast,kt.hi.forecast)
    names(kt.f) <- c("kt forecast","kt lower","kt upper")
    deltat <- object$year[2] - object$year[1]
    kt.f <- ts(kt.f,start=object$year[nyears]+deltat,deltat=deltat)

    # Calculate expected life and mx forecasts
    e0.forecast <- rep(0,h)
    mx.forecast <- matrix(0,nrow=nages,ncol=h)
    colnames(mx.forecast) <- seq(h)
    rownames(mx.forecast) <- object$age
    mx.lo.forecast <- mx.hi.forecast <- mx.forecast
    logjumprates <- log(jumprates)
    series <- names(object)[4]
    agegroup <- object$age[4]-object$age[3]
    for (cnt in 1:h)
    {
        mx.forecast[,cnt] <- fitmx(kt.f[cnt,1], logjumprates, object$bx)
        mx.lo.forecast[,cnt] <- fitmx(kt.f[cnt,2], logjumprates, object$bx)
        mx.hi.forecast[,cnt] <- fitmx(kt.f[cnt,3], logjumprates, object$bx)
        e0.forecast[cnt] <- get.e0(mx.forecast[,cnt],agegroup,series,startage=min(object$age))
    }
    kt.f <- data.frame(kt.forecast,kt.lo.forecast,kt.hi.forecast)
    names(kt.f) <- c("kt forecast","kt lower","kt upper")
    kt.f <- ts(kt.f,start=object$year[nyears]+deltat,deltat=deltat)

    output = list(label=object$label,age=object$age,year=object$year[nyears] + x*deltat,
            rate=list(forecast=mx.forecast,lower=mx.lo.forecast,upper=mx.hi.forecast),
            fitted=object$fitted,
            e0=ts(e0.forecast,start=object$year[nyears]+deltat,deltat=deltat),
            kt.f=structure(list(mean=kt.f[,1],lower=kt.f[,2],upper=kt.f[,3],level=level,x=object$kt,
                                method="Random walk with drift"),class="forecast"),
                type = object$type,lambda=0)
    names(output$rate)[1] = names(object)[4]
    output$model <- object
    output$model$jumpchoice <- jumpchoice
		output$model$jumprates <- jumprates
    output$call <- match.call()
    output$name <- names(object)[4]
    return(structure(output,class=c("fmforecast","demogdata")))
}

#' @rdname residuals.fdm
#' @export
fitted.lca <- function(object,...)
{
    object$fitted
}

#' @rdname residuals.fdm
#' @export
residuals.lca <- function(object,...)
{
    return(structure(list(x=object$year,y=object$age,z=t(object$residuals$y)),class="fmres"))
}


findroot <- function(FUN,guess,margin,try=1,...)
{
    # First try in successively larger intervals around best guess
    for(i in 1:5)
    {
        rooti <- try(stats::uniroot(FUN,interval=guess+i*margin/3*c(-1,1),...),silent=TRUE)
        if(class(rooti) != "try-error")
            return(rooti$root)
    }
    # No luck. Try really big intervals
    rooti <- try(stats::uniroot(FUN,interval=guess+10*margin*c(-1,1),...),silent=TRUE)
    if(class(rooti) != "try-error")
        return(rooti$root)

    # Still no luck. Try guessing root using quadratic approximation
    if(try<3)
    {
        root <- try(quadroot(FUN,guess,10*margin,...),silent=TRUE)
        if(class(root)!="try-error")
            return(findroot(FUN,root,margin,try+1,...))
        root <- try(quadroot(FUN,guess,20*margin,...),silent=TRUE)
        if(class(root)!="try-error")
            return(findroot(FUN,root,margin,try+1,...))
    }

    # Finally try optimization
    root <- try(newroot(FUN,guess,...),silent=TRUE)
    if(class(root)!="try-error")
        return(root)
    else
        root <- try(newroot(FUN,guess-margin,...),silent=TRUE)
    if(class(root)!="try-error")
        return(root)
    else
        root <- try(newroot(FUN,guess+margin,...),silent=TRUE)
    if(class(root)!="try-error")
        return(root)
    else
        stop("Unable to find root")
}

quadroot <- function(FUN,guess,margin,...)
{
    x1 <- guess-margin
    x2 <- guess+margin
    y1 <- FUN(x1,...)
    y2 <- FUN(x2,...)
    y0 <- FUN(guess,...)
    if(is.na(y1) | is.na(y2) | is.na(y0))
        stop("Function not defined on interval")
    b <- 0.5*(y2-y1)/margin
    a <- (0.5*(y1+y2)-y0)/(margin^2)
    tmp <- b^2 - 4*a*y0
    if(tmp < 0)
        stop("No real root")
    tmp <- sqrt(tmp)
    r1 <- 0.5*(tmp-b)/a
    r2 <- 0.5*(-tmp-b)/a
    if(abs(r1) < abs(r2))
        root <- guess+r1
    else
        root <- guess+r2
    return(root)
}

# Try finding root using minimization
newroot <- function(FUN,guess,...)
{
    fred <- function(x,...){(FUN(x,...)^2)}
    junk <- stats::nlm(fred,guess,...)
    if(abs(junk$minimum)/fred(guess,...) > 1e-6)
        warning("No root exists. Returning closest")
    return(junk$estimate)
}


# Replace zeros with interpolated values
fill.zero <- function(x,method="constant")
{
    tt <- 1:length(x)
    zeros <- abs(x) < 1e-9
    xx <- x[!zeros]
    tt <- tt[!zeros]
    x <- stats::approx(tt,xx,1:length(x),method=method,f=0.5,rule=2)
    return(x$y)
}

