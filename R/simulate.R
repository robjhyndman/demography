# Function to simulate future sample paths of functional data
# given output from forecast.fdm

simulate.fmforecast <- function(object,nsim=100,seed=NULL,bootstrap=FALSE, adjust.modelvar=TRUE,...)
{
  n <- length(object$model$year)
  p <- length(object$age)
  h <- length(object$year)

	# Fix lca objects to be able to use this function.
	if(is.element("lca",class(object$model)))
	{
		object$model$basis <- cbind(log(object$model$jumprates), object$model$bx)
		object$model$coeff <- ts(cbind(rep(1,n),object$model$kt))
		stats::tsp(object$model$coeff) <- stats::tsp(object$model$kt)
		zval <- stats::qnorm(0.5 + 0.005*object$kt.f$level)
		# refit model using Arima so it can be simulated more easily.
		ktmod <- Arima(object$model$kt,order=c(0,1,0),include.drift=TRUE)
		# Find stdev of kt from kt.f (to include coefficient error)
		ktmod$sigma2 <- ((object$kt.f$upper[1] - object$kt.f$lower[1])/zval/2)^2
		object$coeff <- list(rwf(rep(1,n),level=object$kt.f$level),forecast(ktmod,level=object$kt.f$level))
		colnames(object$model$basis) <- c("mean","bx")
		rownames(object$model$basis) <- names(object$model$bx)
		adjust.modelvar <- FALSE
	}
	
  nb <- length(object$coeff)

  ridx <- (1:n)#[!is.na(colSums(object$model$residuals$y))]
  # Set residuals to zero for the simulations
  resids <- object$model$residuals$y
  resids[is.na(resids)] <- 0

  fmean <- BoxCox(object$rate[[1]],object$lambda)

  fcoeff <- matrix(1,nrow=h,ncol=nb)
  output <- array(0,c(p,h,nsim))
  for(i in 1:nsim)
  {
    output[,,i] <- object$model$basis[,1]
    if(nb > 1)
    {
      for(j in 2:nb)
      {
        mod <- object$coeff[[j]]$model
        if(class(mod)=="forecast")
          mod <- mod$model
        output[,,i] <- output[,,i] + object$model$basis[,j] %*% matrix(simulate(mod,nsim=h,bootstrap=bootstrap,future=TRUE),nrow=1)
      }
    }
    output[,,i] <- output[,,i] + resids[,sample(ridx,h,replace=TRUE)]
    if(adjust.modelvar)
      output[,,i] <- fmean + sweep(output[,,i]-fmean,1,sqrt(object$var$adj.factor),"*")
  }
  dimnames(output) <- list(object$age,object$year,1:nsim)
  output <- InvBoxCox(output,object$lambda)
  output[is.na(output)] <- 0
  if(object$type != "migration")
    output[output<0] <- 0
  return(output)
}


simulate.fmforecast2 <- function (object, ...) 
{
	is.mortality <- (object$ratio[[1]]$type == "mortality")
	output <- unclass(object) # Just to retain same list structure
	if(is.element("product",names(object))) # Assume coherent model
	{
		output$product <- simulate(object$product,...)
		for(i in 1:length(object$ratio))
		{
			output$ratio[[i]] <- simulate(object$ratio[[i]],...)
			if(is.mortality)
				output[[i]] <- output$product * output$ratio[[i]]
			else
				output[[i]] <- output$product + output$ratio[[i]]
		}
	}
	else
	{
		for(i in 1:length(object))
			output[[i]] <- simulate(object[[i]],...)
	}
    return(output)
}
