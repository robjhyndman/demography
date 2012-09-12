# Parts of this file contain modifications of code taken from 
# the file src/library/stats/R/spline.R and src/library/stats/R/splinefun.R
# which are part of the R package, http://www.R-project.org
# The original code was part of the spline() and splinefun() functions
# and was modified by Simon Wood to allow for monotonic spline interpolation.
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

spl.coef.conv <- function(z)
# takes an object z containing equally lengthed arrays
# z$x, z$y, z$b, z$c, z$d defining a cubic spline interpolating 
# z$x, z$y and forces z$c and z$d to be consistent with z$y and
# z$b (gradient of spline). This is intended for use in conjunction
# with Hyman's monotonicity filter.
# Note that R's spline routine has s''(x)/2 as c and s'''(x)/6 as d.
# (c) Simon N. Wood 
{ 
    n <- length(z$x)
    h <- z$x[2:n]-z$x[1:(n-1)]
    y0 <- z$y[1:(n-1)];y1 <- z$y[2:n]
    b0 <- z$b[1:(n-1)];b1 <- z$b[2:n]
    cc <- -(3*(y0-y1)+(2*b0+b1)*h)/h^2 
    c1 <- (3*(y0[n-1]-y1[n-1])+(b0[n-1]+2*b1[n-1])*h[n-1])/h[n-1]^2
    dd <- (2*(y0-y1)/h+b0+b1)/h^2
    d1 <- dd[n-1]
    z$c <- c(cc,c1);z$d <- c(dd,d1)
    z 
}


hyman.filter <- function(z) 
# Filters cubic spline function to yield co-monotonicity in accordance
# with Hyman (1983) SIAM J. Sci. Stat. Comput. 4(4):645-654, z$x is knot
# position z$y is value at knot z$b is gradient at knot. See also
# Dougherty, Edelman and Hyman 1989 Mathematics of Computation
# 52: 471-494. (c) Simon N. Wood
{
    n <- length(z$x)
    ss <- (z$y[2:n]-z$y[1:(n-1)])/(z$x[2:n]-z$x[1:(n-1)])
    S0 <- c(ss[1],ss)
    S1 <- c(ss,ss[n-1])
    sig <- z$b
    ind <- (S0*S1>0)
    sig[ind] <- S1[ind]
    ind <- (sig>=0)
    if(sum(ind)) 
        z$b[ind] <- pmin(pmax(0,z$b[ind]),3*pmin(abs(S0[ind]),abs(S1[ind])))
    ind <- !ind
    if(sum(ind)) 
        z$b[ind] <- pmax(pmin(0,z$b[ind]),-3*pmin(abs(S0[ind]),abs(S1[ind])))
    z
}

cm.splinefun <- function(x, y = NULL, method = "fmm", gulim=0) 
# modification of stats package splinefun to produce co-monotonic 
#interpolant by Hyman Filtering. if gulim!=0 then it is taken as the upper
#limit on the gradient, and intrpolant is gradient limited rather than 
# monotonic. Modifications from R stats package splinefun() 
# (c) Simon N. Wood 2002 
{ 
    x <- xy.coords(x, y)
    y <- x$y
    x <- x$x
    n <- length(x)
    method <- match(method, c("periodic", "natural", "fmm"))
    if (is.na(method)) 
        stop("splinefun: invalid interpolation method")
    if (any(diff(x) < 0)) 
    {
        z <- order(x)
        x <- x[z]
        y <- y[z]
    }
    if(any(diff(y)<0))
        stop("Data are not monotonic")
    if (method == 1 && y[1] != y[n]) 
    {
        warning("first and last y values differ in spline - using y[1] for both")
        y[n] <- y[1]
    }
    z <- .C("spline_coef", method = as.integer(method), n = n, 
            x = as.double(x), y = as.double(y), b = double(n), c = double(n),
            d = double(n), e = double(if (method == 1) n else 0), PACKAGE = "stats")
    z$y <- z$y-z$x*gulim  # trick to impose upper
    z$b <- z$b-gulim      # limit on interpolator gradient
    z <- hyman.filter(z)  # filter gradients for co-monotonicity
    z$y <- z$y+z$x*gulim  # undo trick 
    z$b <- z$b+gulim      # transformation
    z <- spl.coef.conv(z) # force other coefficients to consistency
    rm(x, y, n, method)
    function(x) 
    {
        .C("spline_eval", z$method, length(x), x = as.double(x), y = double(length(x)), 
                z$n, z$x, z$y, z$b, z$c, z$d, PACKAGE = "stats")$y
    }
}

cm.spline <- function (x, y = NULL, n = 3 * length(x), method = "fmm", xmin = min(x), 
    xmax = max(x),gulim=0) 
# modification of stats package spline to produce co-monotonic 
#interpolant by Hyman Filtering. if gulim!=0 then it is taken as the upper
#limit on the gradient, and intrpolant is gradient limited rather than 
# monotonic. Modifications from R stats package spline() 
{
    x <- xy.coords(x, y)
    y <- x$y
    x <- x$x
    nx <- length(x)
    method <- match(method, c("periodic", "natural", "fmm"))
    if (is.na(method)) 
        stop("spline: invalid interpolation method")
    dx <- diff(x)
    if (any(dx < 0)) 
    {
        o <- order(x)
        x <- x[o]
        y <- y[o]
    }
    if(any(diff(y)<0))
        stop("Data are not monotonic")
    if (method == 1 && y[1] != y[nx]) 
    {
        warning("spline: first and last y values differ - using y[1] for both")
        y[nx] <- y[1]
    }
    z <- .C("spline_coef", method = as.integer(method), n = nx, 
        x = x, y = y, b = double(nx), c = double(nx), d = double(nx), 
        e = double(if (method == 1) nx else 0), PACKAGE = "stats")
    z$y <- z$y-z$x*gulim  # trick to impose upper
    z$b <- z$b-gulim      # limit on interpolator gradient
    z <- hyman.filter(z)  # filter gradients for co-monotonicity
    z$y <- z$y+z$x*gulim  # undo trick 
    z$b <- z$b+gulim      # transformation
    z <- spl.coef.conv(z) # force other coefficients to consistency
    u <- seq(xmin, xmax, length.out = n)
    .C("spline_eval", z$method, nu = length(u), x = u, y = double(n), 
        z$n, z$x, z$y, z$b, z$c, z$d, PACKAGE = "stats")[c("x","y")]
}

# Function to do cubic smoothing spline fit to y ~ x
# with constraint of monotonic increasing for x>b.
# Based on code provided by Simon Wood
# Last updated: 2 June 2004 to work with mgcv 1.0
#
smooth.monotonic <- function(x,y,b,k=-1,w=NULL,newx=x)
{
    require(mgcv)
    weight <- !is.null(w)
    if(k<3 & k!= -1)
        stop("Inappropriate value of k")
    # Unconstrained smooth.
    miss <- is.na(y)
    if(weight)
        miss <- miss | w < 1e-9
    yy <- y[!miss]
    xx <- x[!miss]
    if (weight)
    {
        w <- w[!miss]
        w <- w/sum(w)*length(w)
        f.ug <- gam(yy~s(xx,k=k),weights=w)
        assign("w",w,pos=1)
    }
    else
        f.ug <- gam(yy~s(xx,k=k))

    if(max(xx) <= b)
        return(predict.gam(f.ug,newdata=data.frame(xx=newx),se.fit=TRUE))

    # Create Design matrix, constraints etc. for monotonic spline....
    gam(yy~s(xx,k=k),data=data.frame(xx=xx,yy=yy),fit=FALSE) -> G
    if(weight)
        G$w <- w
    nc <- 200                  # number of constraints
    xc <- seq(b,max(xx),l=nc+1)# points at which to impose constraints
    A0 <- predict.gam(f.ug,data.frame(xx=xc),type="lpmatrix")
                               # A0%*%p will evaluate spline at the xc points
    A1 <- predict.gam(f.ug,data.frame(xx=xc+1e-6),type="lpmatrix")
    A <- (A1-A0)/1e-6          # approximate constraint matrix
                               #(A%%p is  -ve gradient of spline at points xc)

    G$Ain <- A                 # constraint matrix
    G$bin <- rep(0,nc+1)       # constraint vector
    G$sp <- f.ug$sp            # use smoothing parameters from un-constrained fit
    k <- G$smooth[[1]]$df+1
    G$p <- rep(0,k)
    G$p[k] <- 0.1              # get monotonic starting parameters, by 
                               # setting coefficiants of polynomial part of term
    G$p[k-1] <- -mean(0.1*xx)  # must ensure that gam side conditions are
                               # met so that sum of smooth over x's is zero
#    G$p <- rep(0,k+1)
#    G$p[k+1] <- 0.1
#    G$p[k] <- -mean(0.1*xx)
    G$y <- yy
    G$off <- G$off -1          # indexing inconsistency between pcls and internal gam
    p <- pcls(G)               # fit spline (using s.p. from unconstrained fit)

    # now modify the gam object from unconstrained fit a little, to use it
    # for predicting and plotting constrained fit.
    f.ug$coefficients <- p
    return(predict.gam(f.ug,newdata=data.frame(xx=newx),se.fit=TRUE))
}

smooth.monotonic.cobs <- function(x,y,b,lambda=0,w=NULL,newx=x,nknots=50)
{
    require(cobs)
    oldwarn <- options(warn=-1)

    weight <- !is.null(w)
    
    miss <- is.na(y)
    if(weight)
        miss <- miss | w < 1e-9
    yy <- y[!miss]
    xx <- x[!miss]
    if (weight)
    {
        w <- w[!miss]
        w <- w/sum(w)*length(w)
        f.ug <- cobs(xx,yy,w=w,print.warn=FALSE,print.mesg=FALSE,lambda=lambda,nknots=nknots)
    }
    else
        f.ug <- cobs(xx,yy,print.warn=FALSE,print.mesg=FALSE,lambda=lambda,nknots=nknots)

    fred <- predict(f.ug,interval="conf",nz=200)
    fit <- approx(fred[,1],fred[,2],xout=newx)$y
    se <- approx(fred[,1],(fred[,4]-fred[,3])/2/1.96, xout=newx)$y

    if(max(xx) > b)
    {
        delta <- (max(xx)-min(xx))/10
        xxx <- xx[xx>(b-delta)]
        yyy <- yy[xx>(b-delta)]
        if(weight)
            f.mono <- cobs(xxx,yyy,constraint="increase",w=w[xx>(b-delta)],print.warn=FALSE,print.mesg=FALSE,lambda=lambda,nknots=nknots)
        else
            f.mono <- cobs(xxx,yyy,constraint="increase",print.warn=FALSE,print.mesg=FALSE,lambda=lambda,nknots=nknots)
        fred <- predict(f.mono,interval="conf",nz=200)
        newfit <- approx(fred[,1],fred[,2],xout=newx[newx>(b-delta)])$y
        newse <- approx(fred[,1],(fred[,4]-fred[,3])/2/1.96,xout=newx[newx>(b-delta)])$y
        preb <- sum(newx <= (b-delta))
        newfit <- c(rep(0,preb),newfit)
        newse <- c(rep(0,preb),newse)
        postb <- sum(newx > b)
        n <- length(newx)
        cc <- c(rep(0,preb),seq(0,1,length=n-preb-postb),rep(1,postb))
        fit <- (1-cc)*fit + cc*newfit
        se <- (1-cc)*se + cc*newse
    }
    options(warn=oldwarn$warn)
    return(list(fit=fit,se=se))
}
    
    
        
