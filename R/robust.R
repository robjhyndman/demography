#median <- function(...) UseMethod("median")

# median.default <- function (x, na.rm = FALSE)
# {
#     if (is.factor(x) || mode(x) != "numeric")
#         stop("need numeric data")
#     if (na.rm)
#         x <- x[!is.na(x)]
#     else if (any(is.na(x)))
#         return(NA)
#     n <- length(x)
#     if (n == 0)
#         return(NA)
#     half <- (n + 1)/2
#     if (n%%2 == 1) {
#         sort(x, partial = half)[half]
#     }
#     else {
#         sum(sort(x, partial = c(half, half + 1))[c(half, half +
#             1)])/2
#     }
# }


# L1MEDIAN calculates the multivariate L1 median
# X is the data matrix
# tol is the convergence criterium; the iterative proces stops when ||m_k - m_{k+1}|| < tol.
#
# Ref: Hossjer and Croux (1995) "Generalizing Univariate Signed Rank Statistics for Testing
# and Estimating a Multivariate Location Parameter", Non-parametric Statistics, 4, 293-308.
# Assume columnwise location wanted.

L1median <- function(X,tol=1e-6,maxstep=200,na.rm=TRUE,method=c("hossjercroux","coordinate"))
{
    method <- match.arg(method)
    # Coordinatewise median
    if(method=="coordinate")
        return(apply(X,2,stats::median.default,na.rm=na.rm))
    # Gower's algorithm
#    else if(method=="gower")
#        return(gower(X,tol=tol,maxstep=maxstep,na.rm=na.rm))
    # Otherwise use Hossjer and Croux.
    else
        return(hossjercroux(X,tol=tol,maxstep=maxstep,na.rm=na.rm))
}

hossjercroux <- function(X, tol=1e-6, maxstep=100, na.rm=TRUE)
{
    n <- nrow(X)
    p <- ncol(X)
    m=apply(X,2,stats::median.default,na.rm=na.rm)
    hctol <- max(1,min(abs(m),na.rm=na.rm)) * tol
    for(k in 1:maxstep)
    {
        mold <- m
        XX <- sweep(X,2,m)
        dx <- norme(XX)
        if(min(abs(dx))>tol)
            w <- 1/dx
        else
        {
            w <- rep(0,n)
            w[dx>tol] <- 1/dx[dx>tol]
        }
        delta <- colSums(XX*repmat(w/sum(w),1,p),na.rm=na.rm)
        nd <- sqrt(sum(delta^2))
        maxhalf <- ifelse(nd < hctol, 0, log2(nd/hctol))
        m <- mold + delta
        nstep <- 0
        oldmobj <- mrobj(X,mold)
        while((mrobj(X,m) > oldmobj) & (nstep<=maxhalf))
        {
            nstep <- nstep+1
            m <- mold+delta/(2^nstep)
        }
        if (nstep>maxhalf)
            return(mold)
    }
#    warning("Iteration failed")
    return(mold)
}

#NORME calculates the euclidian norm of matrix X
# the output is a columnvector containing the norm of each row

norme <- function(X)
{
    return(sqrt(rowSums(X^2,na.rm=TRUE)))
}

#MROBJ computes objective function in m based on X

mrobj <- function(X,m)
{
    return(sum(norme(sweep(X,2,m))))
}

#repmat replicates the matrix A in an mxn block matrix.
repmat <- function(A,m,n=m)
{
    A <- as.matrix(A)
    tmp <- matrix(rep(t(A),m),nrow=m*nrow(A),byrow=TRUE)
    return(matrix(rep(tmp,n),ncol=n*ncol(tmp)))
}


Qn <- function(x)
{
    n <- length(x)
    diffs <- outer(x, x, "-")
    diffs <- diffs[!lower.tri(diffs,diag=TRUE)]
    qn <- 2.2219*quantile(abs(diffs),0.25)
    if(n==2)
        dn <- 0.399
    else if(n==3)
        dn <- 0.994
    else if(n==4)
        dn <- 0.512
    else if(n==5)
        dn <- 0.844
    else if(n==6)
        dn <- 0.611
    else if(n==7)
        dn <- 0.857
    else if(n==8)
        dn <- 0.669
    else if(n==9)
        dn <- 0.872
    else if(n %% 2 == 1)
        dn <- n/(n+1.4)
    else
        dn <- n/(n+3.8)
    return(dn*qn)
}


MAD <- function(x)
{
    med.x <- median(x)
    return(1.486*median(abs(x-med.x)))
}



# ROBUST PCA USING REFLECTION
# BASED ON HUBERT, ROUSSEEUW AND VERBOVEN

rstep <- function(x,FUN=Qn,order=4,r=matrix.rank(x),mean=TRUE)
{
    if(order<1)
        stop("Order must be positive")

    # transpose to be consistent with Hubert et al.
    X <- t(x)
    p <- ncol(X)
    n <- nrow(X)

    p1 <- min(order,r,floor(n/2))
    S <- numeric(p1)
    Bnorm <- numeric(n)
    V <- eig <- matrix(0,p,p1)
    Transfo <- diag(p)

    if(mean)
    {
        med <- L1median(X,method="hoss")
        xxx <- xx <- sweep(X,2,med)
    }
    else
        xxx <- xx <- X
    for(l in 1:p1)
    {
        B <- xxx
        for(i in 1:n)
            Bnorm[i] <- norm(B[i,],2)
        # Eliminate constant rows
        Bnormr <- Bnorm[Bnorm > 1e-12]
        B <- B[Bnorm > 1e-12,]
        A <- diag(1/Bnormr) %*% B
        Y <- xxx %*% t(A) #projected points in columns
        s <- colQn(Y)
        j <- order(s,decreasing=TRUE)[1]
        S[l] <- s[j]
        V[l:p,l] <- A[j,]

        # EigenVectors = columns of V
        # Constructing Transformation
        Base <- diag(p-l+1)
        ndiff <- norm(Base[,1]-V[l:p,l],Inf) #max norm of the normal vector
        if(ndiff > 1e-12)
        {
            if (sum(V[l:p,l] * Base[,1]) < 0)
                V[l:p,l] <- -V[l:p,l]
            u <- matrix(Base[,1]-V[l:p,l],ncol=1) / c(norm(Base[,1]-V[l:p,l]))
            U <- Base - 2 * repmat(t(u)%*%Base,p-l+1,1) * repmat(u,1,p-l+1)
        }
        else
            U <- Base

        # Transforming eigenvectors to the original pxp dimensional space
        eig[,l] <- Transfo %*% V[,l]
        if(l<p1)
        {
            Edge <- diag(p)
            Edge[l:p,l:p] <- U
            Transfo <- Transfo %*% Edge
            xxx <- xxx %*% U       #Reflection of data
            xxx <- as.matrix(xxx[,-1])
        }
    }
    coef <- xx %*% eig
    if(mean)
    {
        basis <- cbind(med,eig)
        coef <- cbind(rep(1,n),coef)
    }
    else
        basis <- eig

    return(list(basis=basis,coeff=coef,X=xx))
}

# RAPCA ALGORITHM for ROBUST PCA
# BASED ON HUBERT, ROUSSEEUW AND VERBOVEN

rapca <- function(x,FUN=Qn,order=4,mean=TRUE)
{
    if(order<1)
        stop("Order must be positive")

    X <- t(x)
    n <- nrow(X)
    p <- ncol(X)

    # First Step: classical SVD on the data
    # This step reduces the data space to the affine subspace
    # spanned by r=min(n-1,p) observations.
    if(mean)
    {
        med <- colMeans(X)
        xx <- sweep(X,2,med)
    }
    else
        xx <- X
    tmp <- La.svd(xx)
    r = sum(tmp$d > (max(n,p)*max(tmp$d)*1e-12)) # Approx rank
    P <- t(tmp$vt)[,1:r]

    # Second Step: Rstep on X
    # computes the robust eigenvectors and eigenvalues
    tmp2 <- rstep(t(xx %*% P), order=order, r=tmp$r, mean=mean)
    tmp <- P %*% tmp2$basis

    # Retransforming the robust location to the original space
    if(mean)
    {
        med <- c(med + tmp[,1])
        xx <- sweep(X,2,med)
        basis <- cbind(med, tmp[,(1:order)+1])
        coef <- cbind(rep(1,n),xx %*% basis[,-1])
    }
    else
    {
        basis <- tmp
        coef <- xx %*% basis
    }

    return(list(basis=basis,coeff=coef, X=xx))
}

matrix.rank <- function(X)
{
    X.sv <- abs(La.svd(X)$d)
    return(sum((X.sv/max(X.sv)) > 1e-9))
}

norm <- function(A,p=2)
{
    A <- as.matrix(A)
    if(min(dim(A))==1)
        A <- t(A)
    if(p==1)
        return(as.matrix(max(colSums(abs(A)))))
    else if(p==2)
    {
        A.sv <- La.svd(A)$d
        return(as.matrix(max(A.sv)))
    }
    else if(p>1e9)
        return(as.matrix(max(rowSums(abs(A)))))
    else
        stop("Unknown norm")
}

# Qn along columns of Z

colQn <- function(Z)
{
    return(apply(Z,2,Qn))
}
