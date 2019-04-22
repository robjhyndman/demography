
#' Updating functional demographic models and coherent functional demographic models.
#'
#' \code{update.fmforecast()} updates \code{fdm} forecasts. The argument \code{object} is the output from \code{\link{forecast.fdm}} which has been subsequently modified with new coefficient forecasts. These new forecasts are used when re-calculating the forecast of the mortality or fertility rates, or net migration numbers.
#' \code{update.fmforecast2()} updates \code{fdmpr} forecasts. The argument \code{object} is the output from \code{\link{forecast.fdmpr}} which has been subsequently modified with new coefficient forecasts.
#'
#' @param object Output from either \code{\link{fdm}} or \code{\link{coherentfdm}}.
#' @param ... Extra arguments currently ignored.
#'
#' @return A list of the same class as \code{object}.
#' @author Rob J Hyndman.
#' @seealso \code{\link{forecast.fdm}}, \code{\link{forecast.fdmpr}}
#' @examples
#' \dontrun{
#' france.fit <- fdm(fr.mort,order=2)
#' france.fcast <- forecast(france.fit,50)
#' # Replace first coefficient model with ARIMA(0,1,2)+drift
#' france.fcast$coeff[[2]] <- forecast(Arima(france.fit$coeff[,2],
#'                                           order=c(0,1,2), include.drift=TRUE), h=50, level=80)
#' france.fcast <- update(france.fcast)
#'
#' fr.short <- extract.years(fr.sm,1950:2006)
#' fr.fit <- coherentfdm(fr.short)
#' fr.fcast <- forecast(fr.fit)
#' par(mfrow=c(1,2))
#' plot(fr.fcast$male)
#' # Replace first coefficient model in product component with a damped ETS model:
#' fr.fcast$product$coeff[[2]] <- forecast(ets(fr.fit$product$coeff[,2], damped=TRUE),
#'                                         h=50, level=80)
#' fr.fcast <- update(fr.fcast)
#' plot(fr.fcast$male)
#' }
#' @keywords models
#' @name update
NULL

# Update fmforecast object
# Original object from forecast.fdm
# Assumed that the coefficient forecasts have been subsequently changed
# Object needs to be updated to reflect those changes.
#' @rdname update
#' @export
update.fmforecast <- function(object, ...)
{
  if(!is.element("fmforecast",class(object)))
    stop("object must be of class fmforecast")
  h <- length(object$year)
  nb <- ncol(object$model$basis)
  adjust <- length(object$var$adj.factor) > 1
  objnames <- dimnames(object$rate[[1]])

  # Update in-sample fitted values and errors.
  fitted <- matrix(NA,length(object$model$year),nb)
  meanfcast <- varfcast <- matrix(NA, nrow = h, ncol = nb)
  qconf <- stats::qnorm(0.5 + object$coeff[[2]]$level[1]/200)
  fitted[,1] <- 1
  meanfcast[,1] <- 1
  varfcast[,1] <- 0
  for(i in 2:nb)
  {
    fitted[,i] <- fitted(object$coeff[[i]]$model)
    meanfcast[,i] <- object$coeff[[i]]$mean
    varfcast[,i] <- ((object$coeff[[i]]$upper - object$coeff[[i]]$lower)/(2*qconf))^2
  }
  object$fitted$y <- object$model$basis %*% t(fitted)
  object$error$y <- object$model$y$y - object$fitted$y
  object$coeff.error <- object$model$coeff - fitted

  # Update point forecasts
  object$rate[[1]] <- object$model$basis %*% t(meanfcast)
  dimnames(object$rate[[1]]) <- objnames

  # Update forecast variances
  # Only model variance should have changed
  modelvar <- object$model$basis^2 %*% t(varfcast)
  totalvar <- sweep(modelvar, 1, object$var$error + object$var$mean, "+")
  if (adjust & nb > 1)
  {
    object$var$adj.factor <- rowMeans(object$error$y^2, na.rm = TRUE)/totalvar[,1]
    totalvar <- sweep(totalvar, 1, object$var$adj.factor, "*")
  }
  # Add observational variance to total variance
  object$var$total <- sweep(totalvar,1,object$var$observ,"+")

  # Update forecast intervals
  # Only parametric intervals computed here.
  tmp <- qconf * sqrt(object$var$total)
  object$rate$lower <- InvBoxCox(object$rate[[1]] - tmp, object$lambda)
  object$rate$upper <- InvBoxCox(object$rate[[1]] + tmp, object$lambda)
  object$rate[[1]] <- InvBoxCox(object$rate[[1]], object$lambda)
  if(object$type != "migration")
  {
    object$rate[[1]] <- pmax(object$rate[[1]],0.000000001)
    object$rate$lower <- pmax(object$rate$lower,0.000000001)
    object$rate$lower[is.na(object$rate$lower)] <- 0
    object$rate$upper <- pmax(object$rate$upper,0.000000001)
  }

  # Return updated object
  return(object)
}


# Function to combine product and ratio forecasts
# object is output from forecast.fdmpr, but with modified forecasts
# This function simply recombines them again.
#' @rdname update
#' @export
update.fmforecast2 <- function(object, ...)
{
  if(!is.element("fmforecast2",class(object)))
    stop("object must be of class fmforecast2")

  J <- length(object$ratio)
  ny <- length(object$ratio[[1]]$year)

  # GM model
  object$product <- update(object$product)

  # Obtain forecasts for each group
	is.mortality <- (object$product$type=="mortality")
  y <- as.numeric(is.mortality) #=1 for mortality and 0 for migration
  for (j in 1:J)
  {
    object$ratio[[j]] <- update(object$ratio[[j]])
    if(is.mortality)
      object[[j]]$rate[[1]] <- object$product$rate$product * object$ratio[[j]]$rate[[1]]
    else
      object[[j]]$rate[[1]] <- object$product$rate$product + object$ratio[[j]]$rate[[1]]
    if(is.mortality)
      y <- y * object[[j]]$rate[[1]]
    else
      y <- y + object[[j]]$rate[[1]]
  }

  # Adjust forecasts so they multiply appropriately.
	if(is.mortality)
	{
    y <- y^(1/J)/object$product$rate$product
		for(j in 1:J)
			object[[j]]$rate[[1]] <- object[[j]]$rate[[1]]/y
	}
	else
	{
    y <- y/J - object$product$rate$product
		for(j in 1:J)
			object[[j]]$rate[[1]] <- object[[j]]$rate[[1]]-y
	}
  # Variance of forecasts
  qconf <- 2 * stats::qnorm(0.5 + object$product$coeff[[1]]$level/200)
  for (j in 1:J)
  {
    vartotal <- object$product$var$total + object$ratio[[j]]$var$total
    tmp <- qconf * sqrt(vartotal)
    object[[j]]$rate$lower <- InvBoxCox(BoxCox(object[[j]]$rate[[1]],object$product$lambda) - tmp, object$product$lambda)
    object[[j]]$rate$upper <- InvBoxCox(BoxCox(object[[j]]$rate[[1]],object$product$lambda) + tmp, object$product$lambda)
  }

  return(object)
}
