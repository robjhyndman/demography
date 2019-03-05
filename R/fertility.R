
#' Compute total fertility rate from fertility rates
#'
#' Compute total fertility rates from age-specific fertility rates contained in
#' a \code{demogdata} object.
#'
#' @param data Demogdata object of type \code{"fertility"} such as obtained from
#'   \code{\link{read.demogdata}}, \code{\link{forecast.fdm}}.
#' @param PI If TRUE, produce a prediction interval.
#' @param nsim Number of simulations to use when computing a prediction
#'   interval.
#' @param ... Other arguments passed to \code{simulate} when producing
#'   prediction intervals.
#'
#' @return If data are of class \code{demogdata}, the function returns a time
#'   series of fertility rates. If data are from \code{\link{forecast.fdm}}, the
#'   function returns an object of class \code{forecast} containing point
#'   forecasts and (optionally) prediction intervals.
#'
#' @author Rob J Hyndman
#' @seealso \code{\link{fdm}}
#' @examples 
#' plot(tfr(aus.fert))
#' ausfert.fcast <- forecast(fdm(aus.fert))
#' plot(tfr(ausfert.fcast,PI=TRUE,nsim=400))
#' 
#' @keywords models
#' @export
tfr <- function(data, PI=FALSE, nsim=500, ...)
{
    if(!is.element("demogdata",class(data)))
        stop("data must be a demogdata object")
    if(data$type != "fertility")
        stop("data must be a fertility object")

    agegroup <- data$age[3]-data$age[2]
    n <- length(data$rate)
    tfr.mat <- matrix(NA,ncol=n,nrow=length(data$year))
    for(j in 1:n)
        tfr.mat[,j] <- colSums(data$rate[[j]],na.rm=TRUE)*agegroup
    out <- ts(tfr.mat,start=data$year[1],frequency=1)
    if(is.element("fmforecast",class(data)))
    {
			hdata <- data$model
			hdata$rate <- list(InvBoxCox(hdata[[4]],data$lambda))
			names(hdata$rate) <- names(hdata)[4]
			if(!is.null(data$model$pop))
			{
				hdata$pop = list(data$model$pop)
				names(hdata$pop) <- names(hdata$rate)
			}
			class(hdata) <- "demogdata"
			out <- structure(list(x=tfr(hdata),mean=out[,1], method="FDM model"),class="forecast")
			if(PI) # Compute prediction intervals
			{
				sim <- simulate(data,nsim,...)
				tfrsim <- matrix(NA,dim(sim)[2],dim(sim)[3])
				simdata <- data
				for(i in 1:dim(sim)[3])
				{
					simdata$rate[[1]] <- as.matrix(sim[,,i])
					tfrsim[,i] <- tfr(simdata,PI=FALSE)$mean
				}
				#browser()
				out$level <- data$coeff[[1]]$level
				out$lower <- apply(tfrsim,1,quantile,prob=0.5 - out$level/200)
				out$upper <- apply(tfrsim,1,quantile,prob=0.5 + out$level/200)
				out$sim <- sim
			}
    }
    return(out)
}
