####################################################
### FUNCTIONS FOR FUNCTIONAL DEMOGRAPHIC MODELS  ###
####################################################

fdm <- function(data, series=names(data$rate)[1], order=6, ages=data$age, max.age=max(ages),
    method=c("classical","M","rapca"), lambda=3, mean=TRUE, level=FALSE, transform=TRUE,...)
{
    series <- tolower(series)
    method <- match.arg(method)
    data <- extract.ages(data,ages,combine.upper=FALSE)
    if(max.age < max(ages))
        data <- extract.ages(data,min(ages):max.age,combine.upper=TRUE)
    ages <- data$age

    if(data$type=="mortality")
        yname <- "Mortality rate"
    else if(data$type=="fertility")
        yname <- "Fertility rate"
    else if(data$type=="migration")
        yname <- "Net migration"
    else
        stop("Not sure what to do with this data type")
    if(transform)
    {
        mx <- BoxCox(get.series(data$rate,series),data$lambda)
        mx[mx < -1e9] <- NA
    }
    else
        mx <- get.series(data$rate,series)

    data.fts <- fts(ages,mx,start=data$year[1],xname="Age",yname=yname)
    fit <- ftsm(data.fts, order=order, method=method, mean=mean, level=level, lambda=lambda, ...)

	# Adjust signs of output so that basis functions are primarily positive. (Easier to interpret.)
	nx <- length(data$age)
	for(i in 1+(1:order))
	{
		if(sum(stats::na.omit(fit$basis[,i]) > 0) < nx/2)
		{
			fit$basis[,i] <- -fit$basis[,i]
			fit$coeff[,i] <- -fit$coeff[,i]
		}
	}

    if(is.element("obs.var",names(data)))
        ov <- data$obs.var[[match(series,tolower(names(data$rate)))]]
    else
        ov <- mx*0

    output <- list(label=data$label,age=fit$y$x,year=data$year,mx=mx,
        pop=get.series(data$pop,series),
        fitted=fit$fitted,
        residuals=fit$residuals,
        basis=fit$basis,coeff=fit$coeff,
        mean.se=fit$mean.se,varprop=fit$varprop, weights=fit$wt,wt=fit$wt,
        obs.var=ov,
        v=fit$v,type=data$type,
        y=data.fts,
        basis2 = fit$basis2,
        coeff2 = fit$coeff2,
        lambda = data$lambda,
        call=match.call())
    names(output)[4] <- series
    class(output)=c("fdm","ftsm")
    return(output)
}


print.fdm <- function(x,...)
{
    cat("Functional demographic model\n")
    cat(paste("\nCall:",deparse(x$call,200),"\n"))
    cat(paste("\nRegion:"),x$label)
    cat(paste("\nData type:"),x$type)
    cat(paste("\nYears in fit:",min(x$year),"-",max(x$year)))
    cat(paste("\nAges in fit:",min(x$age),"-",max(x$age),"\n"))
    ord <- ncol(x$basis)-1
    cat(paste("\nOrder:",ord))
    if(ord>1)
    {
        cat("\nPercentage variation due to basis functions: ")
        cat(paste(formatC(x$varprop*100,1,format="f"),"%",sep=""))
    }
    cat("\n")
}

summary.fdm <- function(object,...)
{
    print(object)
    junk <- fdmMISE(object[[4]],object$fitted$y,age=object$age,years=object$year)
    junk1 <- cbind(junk$ME,junk$MSE,junk$MPE,junk$MAPE)
    rownames(junk1) <- object$age
    colnames(junk1) <- c("ME","MSE","MPE","MAPE")
    junk2 <- cbind(junk$MIE,junk$MISE,junk$MIPE,junk$MIAPE)
    rownames(junk2) = object$year
    colnames(junk2) = c("IE","ISE","IPE","IAPE")
    cat(paste("\nAverages across ages:\n"))
    print(round(colMeans(junk1),5))
    cat(paste("\nAverages across years:\n"))
    print(round(colMeans(junk2),5))
    cat("\n")
}

summary.fdmpr <- function(object, ...)
{
	cat("*** PRODUCT MODEL ***\n")
	summary(object$product)
	cat("\n*** RATIO MODELS ***\n")
	for(i in 1:length(object$ratio))
	{
		cat(paste("\n",toupper(names(object$ratio))[i],"\n",sep=""))
		summary(object$ratio[[i]])
	}
}
residuals.fdm <- function(object,...)
{
    return(structure(list(x=object$year,y=object$age,z=t(object$residuals$y)),class="fmres"))
}

forecast.fdm <- function(object, h=50, level=80, jumpchoice=c("fit","actual"),
    method="arima", warnings=FALSE, ...)
{
    jumpchoice <- match.arg(jumpchoice)

    if(sum(object$weights < 0.1)/length(object$weights) > 0.2) # Probably exponential weights for fitting. Can be ignored for forecasting
        object$weights[object$weights > 0] <- 1

    if(warnings)
      fcast <- forecast.ftsm(object,h=h,level=level,jumpchoice=jumpchoice,method=method,...)
    else
      suppressWarnings(fcast <- forecast.ftsm(object,h=h,level=level,jumpchoice=jumpchoice,method=method,...))

    # Compute observational variance
    # Update to deal with general transformations
#    fred <- InvBoxCox(object[[4]],object$lambda)
#    if(object$type != "migration")
#        fred <- pmax(fred,0.000000001)
#    if(object$type == "mortality") # Use binomial distribution
#        s <- sqrt((1-fred)*fred^(2*object$lambda-1)/pmax(object$pop,1))
#    else if(object$type=="fertility") # Use Poisson distribution
#        s <- sqrt((fred/1000)^(2*object$lambda-1)/pmax(object$pop,1))
#    else
#      {
        s <- sqrt(abs(object$obs.var))
#      }
#    browser()
    ysd <- s*NA
    for(i in 1:ncol(ysd))
       ysd[,i] <- stats::spline(object$age, s[,i], n=nrow(fcast$var$total))$y
    ysd <- rowMeans(ysd)

#    browser()

    # Add observational variance to total variance
    fcast$var$observ <- ysd^2
    fcast$var$total <- sweep(fcast$var$total,1,fcast$var$observ,"+")

    # Correct prediction intervals and reverse transform
    fse <- stats::qnorm(.5 + fcast$coeff[[1]]$level/200) * sqrt(fcast$var$total)
    fcast$upper$y <- InvBoxCox(fcast$mean$y + fse,object$lambda)
    fcast$lower$y <- InvBoxCox(fcast$mean$y - fse,object$lambda)
    fcast$mean$y <- InvBoxCox(fcast$mean$y,object$lambda)
    if(object$type != "migration")
    {
        fcast$mean$y <- pmax(fcast$mean$y, 0.000000001)
        fcast$lower$y <- pmax(fcast$lower$y, 0.000000001)
        fcast$lower$y[is.na(fcast$lower$y)] <- 0
        fcast$upper$y <- pmax(fcast$upper$y, 0.000000001)
    }
#    if(object$type != "migration")
#    {
#        fcast$mean$y[is.na(fcast$mean$y)] <- 0
#        fcast$lower$y[is.na(fcast$upper$y)] <- 0
#        fcast$upper$y[is.na(fcast$lower$y)] <- 0
#    }
    output <- list(
        label=object$label,
        age=object$age,
        year=max(object$year)+(1:h)/stats::frequency(object$year),
        rate=list(forecast=fcast$mean$y,
                  lower=fcast$lower$y,
                  upper=fcast$upper$y),
        error=fcast$error,
        fitted=fcast$fitted,
        coeff=fcast$coeff,
        coeff.error=fcast$coeff.error,
        var=fcast$var,
        model=fcast$model,
        type=object$type,
        lambda=object$lambda)
    names(output$rate)[1] = names(object)[4]
    output$call <- match.call()
    return(structure(output,class=c("fmforecast","demogdata")))
}


print.fmforecast <- function(x,...)
{
    cat(paste("Forecasts for",x$label))
    cat(paste("\nData type:"),x$type,"\n\n")
    cat(paste("  Call:"),deparse(x$call))
    cat("\n  Based on model:",deparse(x$model$call))
    if(is.element("order",names(x$model)))
        cat(paste("\n  Model order:",x$model$order))
    if(is.element("adjust",names(x$model)))
        cat(paste("\n  Adjustment method:",x$model$adjust))
    if(is.element("jumpchoice",names(x$model)))
        cat(paste("\n  Jump-off method:",x$model$jumpchoice))
    cat(paste("\n\n    Years:",min(x$year),"-",max(x$year)))
    cat(paste("\n    Ages: ",min(x$age),"-",max(x$age),"\n"))
}

plot.fmforecast <- function(x,plot.type=c("function","component","variance"),vcol=1:4,mean.lab="Mean",
    xlab2="Year",h=1,...)
{
    plot.type=match.arg(plot.type)
    if(plot.type=="function")
        plot.demogdata(x,...)
    else if(plot.type=="variance")
    {
        ylim = range(x$var$model[,h],x$var$mean,x$var$error,x$var$observ,na.rm=TRUE)
        plot(x$age,x$var$model[,h],type="l",xlab="Age",ylab="Variance",col=vcol[1],ylim=ylim,...)
        abline(0,0,lty=2)
        lines(x$age,x$var$mean,col=vcol[2])
        lines(x$age,x$var$error,col=vcol[3])
        lines(x$age,x$var$observ,col=vcol[4])
    }
    else
    {
        if(is.element("ax",names(x$model))) #output from lca()
        {
            x$model$basis <- cbind(x$model$ax,x$model$bx)
            x$modelcoeff <- cbind(rep(1,length(x$model$kt)),x$model$kt)
            x$coeff <- list(NULL,x$kt.f)
            colnames(x$model$basis) <- c("mean","bx")
            if(x$model$adjust != "none")
                xlab <- "kt (adjusted)"
            else
                xlab <- "kt"
            plot.ftsf(x, "component", xlab2=xlab2, ylab2=xlab, mean.lab="ax", ...)
        }
        else
            plot.ftsf(x, "component", xlab2=xlab2, mean.lab=mean.lab, ...)
    }
}

models <- function(object, ...)
UseMethod("models")


models.fmforecast <- function(object, select=0, ...)
{
	if(!select[1])
		select <- 1:(length(object$coeff)-1)
	for(i in select)
	{
		cat("\n-- Coefficient",i,"--\n")
        mod <- object$coeff[[i+1]]$model$model 
        if(is.null(mod))
            mod <- object$coeff[[i+1]]$model
		print(mod)
	}
}

models.fmforecast2 <- function(object, ...)
{
	cat("\n************* PRODUCT MODEL *************\n")
	models(object$product,...)
	cat("\n\n\n************* RATIO MODELS *************")
	for(i in 1:length(object$ratio))
	{
		cat("\n\n\n***********", toupper(names(object$ratio)[i]),"***********\n\n")
		models(object$ratio[[i]],...)
	}
}

### fdmMISE
## Inputs:
##   actual = observed data
##   estimate = fitted or forecast data
## These must be matrices of the same order
## Assumed that each column is one year

fdmMISE <- function(actual,estimate,age=NULL,years=NULL,neval=1000)
{
    p <- nrow(actual)
    n <- ncol(actual)
    if(is.null(age))
        age <- 0:(p-1)
    if(is.null(years))
        years <- 1:n
    if(p != nrow(estimate) | n != ncol(estimate) | n != length(years))
        stop("Dimensions of inputs don't match")
    p <- length(age)
    actual = actual[1:p,]
    estimate = estimate[1:p,]
    out <- ftsaMISE(fts(age,actual,start=years[1],frequency=1),fts(age,estimate,start=years[1],frequency=1),neval=neval)
    out$age <- age
    out$years <- years
    return(out)
}

# Following function adapted from MISE borrowed from the ftsa package
# Originally written by RJH a long time ago before Han took over ftsa package
# Not exported by ftsa, and so included here to avoid using ::: in a call
ftsaMISE <- function (actual, estimate, neval = 1000) 
{
  m <- stats::frequency(actual$time)
  s <- stats::start(actual$time)
  x <- actual$x
  p <- length(x)
  n <- ncol(actual$y)
  if (length(estimate$x)!=p | nrow(actual$y)!=p | p!=nrow(estimate$y) | n!=ncol(estimate$y)) 
    stop("Dimensions of inputs don't match")
  if (max(abs(actual$x - estimate$x)) > 1e-05) 
    stop("Different x variables")
  e <- estimate$y - actual$y
  e.big <- pe.big <- matrix(NA, nrow = neval, ncol = n)
  pe <- e/actual$y
  pe.big <- matrix(NA, nrow = neval, ncol = n)
  for (i in 1:n) 
  {
    if (sum(is.na(e[, i])) == 0) 
    {
      idx <- (abs(e[,i]) == Inf) | is.na(e[,i])
      e.big[,i] <- stats::spline(x[!idx], e[!idx, i], n = neval, method = "natural")$y
      idx <- (abs(pe[,i]) == Inf) | is.na(pe[,i])
      pe.big[,i] <- stats::spline(x[!idx], pe[!idx,i], n = neval, method = "natural")$y
    }
  }
  delta <- (max(x) - min(x))/neval
  out <- list(x = x, error = e)
  out$MIE <- ts(colSums(e.big, na.rm = TRUE) * delta, start = s, frequency = m)
  out$MIAE <- ts(colSums(abs(e.big), na.rm = TRUE) * delta, start = s, frequency = m)
  out$MISE <- ts(colSums(e.big^2, na.rm = TRUE) * delta, start = s, frequency = m)
  out$ME <- rowMeans(e, na.rm = TRUE)
  out$MAE <- rowMeans(abs(e), na.rm = TRUE)
  out$MSE <- rowMeans(e^2, na.rm = TRUE)
  out$MIPE <- ts(colSums(pe.big, na.rm = TRUE) * delta,  start = s, frequency = m)
  out$MIAPE <- ts(colSums(abs(pe.big), na.rm = TRUE) * delta, start = s, frequency = m)
  out$MPE <- rowMeans(pe, na.rm = TRUE)
  out$MAPE <- rowMeans(abs(pe), na.rm = TRUE)
  return(out)
}

compare.demogdata <- function(data, forecast, series=names(forecast$rate)[1],
    ages = data$age, max.age=min(max(data$age),max(forecast$age)), years=data$year,
    interpolate=FALSE)
{
    years <- years[sort(match(forecast$year,years))]
    ages <- ages[ages <= max.age]
    if(length(years)==0)
        stop("No common years between data and forecasts")
    subdata <- extract.ages(extract.years(data,years=years),ages)
    forecast <- extract.ages(extract.years(forecast,years=years),ages)
    ages <- subdata$age
    mx <- get.series(subdata$rate,series)
    fmx <- get.series(forecast$rate,series)
    n <- nrow(mx)
    log.mx <- BoxCox(mx,data$lambda)
    if (interpolate)
    {
        if (sum(abs(mx) < 1e-09) > 0)
            warning("Replacing zero values with estimates")
        for (i in 1:n)
            log.mx[i, ] <- BoxCox(fill.zero(mx[i, ]),data$lambda)
    }
    junka <- fdmMISE(log.mx,BoxCox(fmx,forecast$lambda),ages,years)
    junkb <- fdmMISE(mx,fmx,ages,years)
    junk1 <- cbind(junka$ME,junka$MAE,junka$MSE,junkb$MPE,junkb$MAPE)
    rownames(junk1) <- ages
    colnames(junk1) <- c("ME","MAE","MSE","MPE","MAPE")
    junk2 <- cbind(junka$MIE,junka$MIAE,junka$MISE,junkb$MIPE,junkb$MIAPE)
    rownames(junk2) = years
    colnames(junk2) = c("IE","IAE","ISE","IPE","IAPE")
    fred <- list(label=data$label,age=ages,year=years,error=junkb$error,terror=junka$error,
        mean.error=junk1,int.error=ts(junk2,start=years[1],frequency=1))
    names(fred)[4] <- paste(series,"error")
    names(fred)[5] <- paste(series,"transformed error")

    # Add life expectancies
    if(data$type=="mortality")
    {
        actual.le <- life.expectancy(subdata,series=series)
        forecast.le <- life.expectancy(forecast,series=series)
        fred$life.expectancy <- cbind(actual.le,forecast.le,forecast.le-actual.le)
        dimnames(fred$life.expectancy)[[2]] <- c("Actual","Forecast","Error")
    }

    return(structure(fred,class="errorfdm"))
}

fitted.fdm <- function(object,...)
{
    object$fitted
}

print.errorfdm <- function(x,...)
{
    cat(paste("Demographic comparison for",x$label,"\n"))
    cat(paste("  Years: ",min(x$year),"-",max(x$year),"\n"))
    cat(paste("  Ages: ",min(x$age),"-",max(x$age),"\n"))
    fred1 <- apply(x$mean.error,2,mean)
    fred2 <- apply(x$int.error,2,mean)
    cat("\nTransformed data: Errors averaged across time and ages\n")
    print(fred1[1:3])
    cat("\nTransformed data: Errors averaged across time and integrated across ages\n")
    print(fred2[1:3])
    cat("\nRaw data: Percentage errors averaged across time and ages\n")
    print(fred1[4:5])
    cat("\nRaw data: Percentage errors averaged across time and integrated across ages\n")
    print(fred2[4:5])
    if(is.element("life.expectancy",names(x)))
    {
        cat("\nLife expectancy\n")
        x$life.expectancy <- rbind(x$life.expectancy,colMeans(x$life.expectancy))
        dimnames(x$life.expectancy)[[1]] <- c(paste(x$year),"Mean")
        print(x$life.expectancy)
    }
}

plot.errorfdm <- function(x,transform=TRUE,...)
{
    i <- ifelse(transform,5,4)
    plot(fts(x=x$age,y=x[[i]],start=x$year[1],frequency=1,xname="Age",yname=names(x)[i]),...)
}

isfe <- function(...) UseMethod("isfe")

isfe.demogdata <- function(data,series=names(data$rate)[1],max.order=N-3,N=10,h=5:10,
        ages=data$age, max.age=max(ages),
        method=c("classical","M","rapca"), fmethod=c("arima", "ar", "arfima", "ets","ets.na","struct","rwdrift","rw"),
        lambda=3, ...)
{
    series <- tolower(series)
    method <- match.arg(method)
    fmethod <- match.arg(fmethod)

    data <- extract.ages(data,ages,combine.upper=FALSE)
    if(max.age < max(ages))
        data <- extract.ages(data,min(ages):max.age,combine.upper=TRUE)
    ages <- data$age
    mx <- BoxCox(get.series(data$rate,series),data$lambda)
    mx[mx < -1e9] <- NA
    data.fts <- fts(ages,mx,start=data$year[1],xname="Age",yname="")
    return(isfe(data.fts,max.order=max.order,N=N,h=h,method=method,fmethod=fmethod,lambda=lambda,...))
}


summary.fmforecast <- function(object,...)
{
    print(object)

    cat(sprintf("\nERROR MEASURES BASED ON %s RATES\n", toupper(object$type)))
    printout(fdmMISE(object$model$y$y,exp(object$fitted$y),age=object$age,years=object$model$year))

    cat(sprintf("\nERROR MEASURES BASED ON LOG %s RATES\n", toupper(object$type)))
    printout(fdmMISE(log(object$model$y$y),object$fitted$y,age=object$age,years=object$model$year))
}

summary.fmforecast2 <- function(object,...)
{
	if(is.element("product",names(object))) # Assume coherent model
	{
		summary(object$product)
		for(i in 1:length(object$ratio))
			summary(object$ratio[[i]])
	}
	else
	{
		for(i in 1:length(object))
			summary(object[[i]])
	}
}
