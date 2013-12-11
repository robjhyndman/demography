## SMOOTH FTS
smooth.fts <- function(data, k=-1, xgrid=data$x, se.fit=FALSE, w=rep(1,nrow(data$y)))
{
    result <- seresult <- data
    x <- data$x
    result$x <- seresult$x <- xgrid
    ny <- ncol(data$y)
    result$y <- seresult$y <- matrix(NA,ncol=ny,nrow=length(xgrid))
    kvec <- numeric(ny)
    meany <- mean(data)$y
    data$y <- sweep(data$y,1,meany)
    if(is.null(dim(w)))
        w <- matrix(rep(w,ny),ncol=ny)

    # Find optimal k
    if(k<0)
    {
        for(j in 1:ny)
        {
            fit <- mgcv::gam(data$y[,j] ~ s(x,k=k), weights=w[,j])
            kvec[j] <- sum(fit$edf)
        }
        k <- round(median(kvec)+.5)
    }
    # Smooth using chosen k
    for(j in 1:ny)
    {
        fit <- mgcv::gam(data$y[,j] ~ s(x,k=k), weights=w[,j])
        smooth.fit <- mgcv::predict.gam(fit,newdata=data.frame(x=xgrid),se.fit=se.fit)
        if(se.fit)
        {
            result$y[,j] <- smooth.fit$fit
            seresults$y[,j] <- smooth.fit$se.fit
        }
        else
            result$y[,j] <- smooth.fit
    }
    interp <- spline(x,meany,n=500)
    interp <- approx(interp$x,interp$y,xout=xgrid)$y
    result$y <- sweep(result$y,1,interp,"+")
    if(se.fit)
        return(list(result,seresult))
    else
        return(result)
}

## Function to smooth mortality curves
## Divides age into three sections: 0-a, a-b and b+
## Will interpolate first period (0-a)
## Will smooth second period (a-b)
## For third period (b+) it uses montonically increasing smooths if monotonic TRUE
## Number of knots for smoothing set by k

smooth.demogdata <- function(data,method=switch(data$type,mortality="mspline",fertility="cspline",migration="loess"), age.grid,
    power=switch(data$type,mortality=0.4,fertility=1,migration=1),
    b=65, k=30, span=0.2, lambda=1e-10, interpolate=FALSE, weight=data$type != "migration",
    obs.var="empirical")
{
    method <- c("mspline","msplinecobs","cspline","spline","loess")[pmatch(method,c("mspline","msplinecobs","cspline","spline","loess"))]
    if(is.na(method))
        stop("Unknown smoothing method")
    obs.var <- c("empirical","theoretical")[pmatch(obs.var,c("empirical","theoretical"))]
    if(is.na(obs.var))
        stop("Unknown method for observational variance")

    # Smooth logged fertility and mortality data
    dlambda <- data$lambda
    if(data$type != "migration")
        data$lambda <- 0 # Ensure smooth curves are positive.

    minx <- min(data$age)
    maxx <- max(data$age)
    nx <- length(data$age)
    delta1 <- data$age[2]-data$age[1]
    delta2 <- data$age[nx] - data$age[nx-1]

    ## Construct upper ages for each group
    xmin <- minx + 0.5 - 0.5*delta1
    xmax <- maxx + 0.5 + 0.5*delta2
    upperage <- c(xmin, 0.5*(data$age[1:(nx-1)] + data$age[2:nx] + 1), xmax)-1

    # Construct age.grid if missing
    if(missing(age.grid))
    {
        xmin <- minx - 0.5*delta1 + 0.5
        xmax <- maxx + 0.5*delta2 - 0.5
    }
    else
    {
        xmin <- min(age.grid)
        xmax <- max(age.grid)
    }
    xmin <- floor(xmin)
    xmax <- ceiling(xmax)
    age.grid <- seq(xmin,xmax,by=1)

    if(method=="spline")
    {
        method <- "mspline"
        b <- 1000
        a <- -1000
    }

    n <- length(data$rate)

    if(data$type == "migration" | interpolate)
        weight <- FALSE
    if(weight)
        data$wt <- use.weight(data)
    if(obs.var=="theoretical")
        datawt <- use.weight(data,FALSE)
    fred <- matrix(1,nrow=length(age.grid),ncol=ncol(data$rate[[1]]))
    data$serate <- data$rate
    x <- data$age^power
    data$obs.var <- list()
    for(i in 1:n)
    {
        y <- BoxCox(as.matrix(data$rate[[i]]),data$lambda)
        if(weight)
            y[y< -1e20] <- -10
        data$obs.var[[i]] <- y*NA

        ny <- ncol(y)
        p <- nrow(y)
        newpop <- newy <- se.y <- matrix(NA,ncol=ny,nrow=length(age.grid))
        err <- y*NA
        for(j in 1:ny)
        {
            if(weight)
                w <- as.matrix(data$wt[[i]])[,j]
            else 
                w <- ww <- rep(1,p)
            w[y[,j] < -1e20] <- 0
            xx <- x[w>0]
            yy <- y[w>0,j]
            ww <- w[w>0]
            if(sum(is.na(y[,j]) | y[,j]==0) == p)
                smooth.fit <- list(fit=rep(NA,p),se=rep(NA,p))
			else if(interpolate)
				smooth.fit <- list(fit=approx(xx,yy,xout=age.grid^power,rule=1)$y,se=rep(0,length(age.grid)))
            else if(method=="loess")
            {
                fit <- loess(yy ~ xx,span=span,degree=2,weights=ww,surface="direct")
                smooth.fit <- predict(fit,newdata=data.frame(xx=age.grid^power),se=TRUE)
            }
            else if(method=="mspline")
                smooth.fit <- smooth.monotonic(xx, yy, b^power, max(min(round(length(xx)*.8),k),4), ww, age.grid^power)
            else if(method=="msplinecobs")
                smooth.fit <- smooth.monotonic.cobs(xx, yy, b^power, lambda=lambda, w=ww, newx=age.grid^power, nknots=k)
            else if(method=="cspline")
                smooth.fit <- fert.curve(xx,yy,ww,age.grid^power,lambda=lambda,interpolate=interpolate,tlambda=data$lambda)
            newy[,j] <- smooth.fit$fit
            se.y[,j] <- smooth.fit$se
            newpop[,j] <- diff(cm.spline(upperage,cumsum(c(0,data$pop[[i]][,j])),xmin=xmin-1,xmax=xmax,n=xmax-xmin+2)$y)
            if(sum(!is.na(newy[,j]))>2)
                err[,j] <- y[,j] - approx(age.grid,newy[,j],xout=data$age)$y
            if(obs.var=="theoretical")
                fred[,j] <- approx(data$age,datawt[[i]][,j],xout=age.grid,rule=2)$y
        }
        dimnames(newy) <- dimnames(newpop) <- dimnames(se.y) <- list(age.grid,data$year)
        data$rate[[i]] <- InvBoxCox(newy,data$lambda)
		if(interpolate)
			data$serate[[i]] <- data$obs.var[[i]] <- matrix(0,ncol=ny,nrow=length(age.grid))
        else if(obs.var=="theoretical")
        {
            data$serate[[i]] <- se.y/data$rate[[i]] # Needs fixing if lambda != 0
            data$obs.var[[i]] <- 1/fred
            data$obs.var[[i]][abs(data$obs.var[[i]])>1e9] <- max(data$obs.var[[i]][data$obs.var[[i]]<1e9])
        }
        else
        {
            data$serate[[i]] <- se.y
            y <- InvBoxCox(y,data$lambda)
            yy <- y
            for(j in 1:ny)
                yy[,j] <- approx(age.grid,newy[,j],xout=data$age)$y
            yy <- InvBoxCox(yy,data$lambda)
            xx <- data$age
            ov <- exp(predict(loess(log(rowMeans((y-yy)^2)) ~ xx,span=2/sqrt(length(data$age)),degree=2,surface="direct"),newdata=data.frame(xx=age.grid)))
            data$obs.var[[i]] <- matrix(ov,length(age.grid),ny)
        }
        dimnames(data$obs.var[[i]]) <- list(age.grid,data$year)
        data$pop[[i]] <- newpop
        data$err[[i]] <- err
    }
    names(data$obs.var) <- names(data$rate)
    data$age <- age.grid
    data$lambda <- dlambda
    return(data)
}


# Smooth interpolation of data

fert.curve <- function(x,y,w,age.grid,lambda=1,interpolate=TRUE,tlambda,...)
{
#    if(min(age.grid) < min(x))
#    {
#        x <- c(13,x)
#        y <- c(BoxCox(0.001,tlambda),y)
#        w <- c(max(w),w)
#    }
#    if(max(age.grid) > max(x))
#    {
#        x <- c(x,52)
#        y <- c(y,BoxCox(0.001,tlambda))
#        w <- c(w,max(w))
#    }
    w <- w/sum(w)
    # Transformation
    xx <- x#sign(x-30)*abs(x-30)^1.8

    oldwarn <- options(warn=-1)
    # Unweighted smoothing as there seems to be a problem with the cobs function when weights specified
    if(interpolate)
        fred <- predict(cobs(xx,y,constraint="concave",pointwise=cbind(rep(0,length(xx)),xx,y),
            lambda=lambda,print.warn=FALSE,print.mesg=FALSE,maxiter=1e4),interval="conf",nz=200)
    else
        fred <- predict(cobs(xx,y,constraint="concave",
            lambda=lambda,print.warn=FALSE,print.mesg=FALSE,maxiter=1e4),interval="conf",nz=200)
    options(warn=oldwarn$warn)

#    fred[,1] <- abs(fred[,1])^(5/9)*sign(fred[,1])+30

    fit <- approx(fred[,1],fred[,2],xout=age.grid,rule=1)$y
    se <- approx(fred[,1],(fred[,4]-fred[,3])/2/1.96, xout=age.grid,rule=1)$y
    return(list(fit=fit,se=se))
}

## Function to calculate the weight

standardize <- function(x,sumx=1)
{
    return(x/sum(x,na.rm=TRUE)*sumx)
}

use.weight <- function(data,standardize=TRUE)
{
    w <- list()
    n <- length(data$rate)
    for (i in 1:n)
    {
        # Extract rate and population matrices
        rate.dim <- dim(data$rate[[i]])
        rate <- data$rate[[i]]
        pop <- data$pop[[i]]
        if(standardize)
            pop <- pop/max(pop,na.rm=TRUE)
        if(data$type=="fertility")
        {
            rate <- rate/1000
            w[[i]] <- pop*rate^(1-2*data$lambda)
        }
        else if(data$type=="mortality")
            w[[i]] <- pop*rate^(1-2*data$lambda)
        else
            stop("I shouldn't be here!")
		#if(mean(w[[i]],na.rm=TRUE) > 1)
		#	stop("There's a problem. It looks like your rates are all too large.")
		#else 
		if(mean(w[[i]],na.rm=TRUE) < 0)
			stop("There's a problem. Do you have negative rates?")
		w[[i]][w[[i]] > 1e9] <- 0
        w[[i]][w[[i]] < 0] <- 0
        #w[[i]][rate < 1e-9] <- 0
        #w[[i]][BoxCox(rate,data$lambda) < -1e9] <- 0
        if(data$type=="mortality")
            w[[i]][log(rate) > -1e-9] <- 0
        if(standardize)
            w[[i]] <- apply(w[[i]],2,standardize,sumx=rate.dim[1])
        w[[i]][is.na(w[[i]])] <- 0
        colnames(w[[i]]) <- colnames(pop)
        rownames(w[[i]]) <- rownames(pop)
    }
    names(w) <- names(data$rate)
    return(w)
}
