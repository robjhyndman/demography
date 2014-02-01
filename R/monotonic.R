cm.splinefun <- function(x, y = NULL, ...) 
# wrapper for splinefun()
# Function retained for backwards compatibility
{ 
    splinefun(x, y, method="hyman")
}

cm.spline <- function (x, y = NULL, n = 3 * length(x), xmin = min(x), xmax = max(x), ...) 
# wrapper for spline()
# Function retained for backwards compatibility
{
    spline(x, y, n=n, xmin=xmin, xmax=xmax, method="hyman")
}

# Function to do cubic smoothing spline fit to y ~ x
# with constraint of monotonic increasing for x>b.
# Based on code provided by Simon Wood
# Last updated: 1 February 2014

smooth.monotonic <- function(x,y,b,k=-1,w=NULL,newx=x)
{
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
        f.ug <- mgcv::gam(yy~s(xx,k=k),weights=w)
#        assign("w",w,pos=1)
    }
    else
        f.ug <- mgcv::gam(yy~s(xx,k=k))

    if(max(xx) <= b)
        return(mgcv::predict.gam(f.ug,newdata=data.frame(xx=newx),se.fit=TRUE))

    # Create Design matrix, constraints etc. for monotonic spline....
    mgcv::gam(yy~s(xx,k=k),data=data.frame(xx=xx,yy=yy),fit=FALSE) -> G
    if(weight)
        G$w <- w
    nc <- 200                  # number of constraints
    xc <- seq(b,max(xx),l=nc+1)# points at which to impose constraints
    A0 <- mgcv::predict.gam(f.ug,data.frame(xx=xc),type="lpmatrix")
                               # A0%*%p will evaluate spline at the xc points
    A1 <- mgcv::predict.gam(f.ug,data.frame(xx=xc+1e-6),type="lpmatrix")
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
    G$C <- matrix(0,0,0)       # fixed constraint matrix (there are none)
    p <- mgcv::pcls(G)         # fit spline (using s.p. from unconstrained fit)

    # now modify the gam object from unconstrained fit a little, to use it
    # for predicting and plotting constrained fit.
    f.ug$coefficients <- p
    return(mgcv::predict.gam(f.ug,newdata=data.frame(xx=newx),se.fit=TRUE))
}

smooth.monotonic.cobs <- function(x,y,b,lambda=0,w=NULL,newx=x,nknots=50)
{
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
        f.ug <- cobs::cobs(xx,yy,w=w,print.warn=FALSE,print.mesg=FALSE,lambda=lambda,nknots=nknots)
    }
    else
        f.ug <- cobs::cobs(xx,yy,print.warn=FALSE,print.mesg=FALSE,lambda=lambda,nknots=nknots)

    fred <- predict(f.ug,interval="conf",nz=200)
    fit <- approx(fred[,1],fred[,2],xout=newx)$y
    se <- approx(fred[,1],(fred[,4]-fred[,3])/2/1.96, xout=newx)$y

    if(max(xx) > b)
    {
        delta <- (max(xx)-min(xx))/10
        xxx <- xx[xx>(b-delta)]
        yyy <- yy[xx>(b-delta)]
        if(weight)
            f.mono <- cobs::cobs(xxx,yyy,constraint="increase",w=w[xx>(b-delta)],print.warn=FALSE,print.mesg=FALSE,lambda=lambda,nknots=nknots)
        else
            f.mono <- cobs::cobs(xxx,yyy,constraint="increase",print.warn=FALSE,print.mesg=FALSE,lambda=lambda,nknots=nknots)
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

