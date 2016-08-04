
# COHERENT FORECASTING FOR MULTIPLE GROUPS

coherentfdm <- function(data, order1=6, order2=6, ...) 
{
	# Check if missing data
	
  # Use male and female if available
  gps <- names(data$rate)
  if(is.element("male",gps) & is.element("female",gps))
  {
    notneeded <- (1:length(gps))[-match(c("male","female"),gps)]
    for(i in notneeded)
      data$rate[[i]] <- data$pop[[i]] <- NULL
  }

  J <- length(data$rate)
  rate.ratio <- fdm.ratio <- list()

  # Construct ratio and product objects
	is.mortality <- (data$type=="mortality")
  y <- as.numeric(is.mortality)
  pop.total <- 0
  for (j in 1:J) 
  {
  	if(is.mortality)
  		y <- y * data$rate[[j]]^(1/J)
  	else
  		y <- y + data$rate[[j]]/J
    pop.total <- pop.total + data$pop[[j]]
  }
  rate.product <- demogdata(y, pop=pop.total, data$age, data$year, type=data$type, lambda=data$lambda,
      label=data$label, name="product")
  for (j in 1:J)
	{
		if(is.mortality)
		{
			rate.ratio[[j]] <- demogdata(data$rate[[j]]/rate.product$rate[[1]], pop=pop.total, data$age, 
				data$year, type=data$type, label=data$label, name=names(data$rate)[[j]],lambda=data$lambda)
		}
		else
		{
			rate.ratio[[j]] <- demogdata(data$rate[[j]]-rate.product$rate[[1]], pop=pop.total, data$age, 
				data$year, type=data$type, label=data$label, name=names(data$rate)[[j]],lambda=data$lambda)
		}
		# Set infinite rates to missing
		infrates <- rate.ratio[[j]]$rate[[1]] > 1e9
		rate.ratio[[j]]$rate[[1]][infrates] <- NA
	}
    
  # GM model
  fdm.mean <- fdm(rate.product, series="product", order=order1, ...)
  
  # Ratio model
  for (j in 1:J) 
    fdm.ratio[[j]] <- fdm(rate.ratio[[j]], series=names(data$rate)[j], order=order2, ...)
  names(fdm.ratio) <- names(data$rate)
  
  return(structure(list(product=fdm.mean, ratio=fdm.ratio), class="fdmpr"))
}


forecast.fdmpr <- function(object, h=50, level=80, K=100, drange=c(0.0,0.5), ...) 
{
  fcast.ratio <- fc <- totalvar.r <- list()
  J <- length(object$ratio)
  ny <- length(object$ratio[[1]]$year)
  K <- min(K,ny)
	
  # GM model
	fcast.mean <- forecast(object$product, method="arima", h=h, level=level, ...)
  # Make sure first coefficient is not I(1) with drift.
  #mod <- auto.arima(object$product$coeff[,2],d=2)
  #fcast.mean$coeff[[2]] <- forecast(mod, h=h, level=level, ...)
  #fcast.mean <- update(fcast.mean)
  
  # Obtain forecasts for each group
	is.mortality <- (object$product$type=="mortality")
  y <- as.numeric(is.mortality) #=1 for mortality and 0 for migration
  for (j in 1:J) 
  {
    # Use all available data other than last K years
    # As ARFIMA can't handle missing values
    object$ratio[[j]]$weights <- 0*object$ratio[[j]]$weights +1
    if(K < ny)
      object$ratio[[j]]$weights[1:(ny-K)] <- 0
    fcast.ratio[[j]] <- forecast(object$ratio[[j]], h=h, level=level, method="arfima", estim="mle", drange=drange,...)
    fc[[j]] <- fcast.mean
    if(is.mortality)
      fc[[j]]$rate[[1]] <- fcast.mean$rate$product * fcast.ratio[[j]]$rate[[1]]
    else
      fc[[j]]$rate[[1]] <- fcast.mean$rate$product + fcast.ratio[[j]]$rate[[1]]
    names(fc[[j]]$rate)[1] <- names(object$ratio)[j]
    fc[[j]]$coeff <- fc[[j]]$coeff.error <- fc[[j]]$call <- fc[[j]]$var <- NULL
    if(is.mortality)
      y <- y * fc[[j]]$rate[[1]]
    else
      y <- y + fc[[j]]$rate[[1]]
  }

  # Adjust forecasts so they multiply appropriately.
	if(is.mortality)
	{
    y <- y^(1/J)/fcast.mean$rate$product
		for(j in 1:J)
			fc[[j]]$rate[[1]] <- fc[[j]]$rate[[1]]/y
	}
	else
	{
    y <- y/J - fcast.mean$rate$product
		for(j in 1:J)
			fc[[j]]$rate[[1]] <- fc[[j]]$rate[[1]]-y
	}
  # Variance of forecasts
  qconf <- 2 * stats::qnorm(0.5 + fcast.mean$coeff[[1]]$level/200)
  for (j in 1:J) 
  {
    vartotal <- fcast.mean$var$total + fcast.ratio[[j]]$var$total
    tmp <- qconf * sqrt(vartotal)
    fc[[j]]$rate$lower <- InvBoxCox(BoxCox(fc[[j]]$rate[[1]],object$product$lambda) - tmp, object$product$lambda)
    fc[[j]]$rate$upper <- InvBoxCox(BoxCox(fc[[j]]$rate[[1]],object$product$lambda) + tmp, object$product$lambda)
    if(is.mortality)
    {
      fc[[j]]$model[[4]] <- BoxCox(InvBoxCox(object$product[[4]],object$product$lambda) * 
                                   InvBoxCox(object$ratio[[j]][[4]], object$product$lambda),object$product$lambda)
    }
    else
    {
      fc[[j]]$model[[4]] <- BoxCox(InvBoxCox(object$product[[4]],object$product$lambda) + 
                                     InvBoxCox(object$ratio[[j]][[4]], object$product$lambda),object$product$lambda)
    }
    names(fc[[j]]$model)[4] <- names(object$ratio)[j]
    fc[[j]]$coeff <- list(list(level= fcast.mean$coeff[[1]]$level))
  }
    
  names(fc) <- names(fcast.ratio) <- names(object$ratio)
  fc$product <- fcast.mean
  fc$ratio <- fcast.ratio
  
  return(structure(fc, class="fmforecast2"))
}
