#' Forecasting mortality and fertility data
#' 
#' Functions for demographic analysis including lifetable calculations,
#' Lee-Carter modelling and functional data analysis of mortality rates.
#' 
#' 
#' @author Rob J Hyndman with contributions from Heather Booth, Leonie Tickle, John Maindonald, Simon Wood and the R Core Team.  
#' @author Maintainer: <Rob.Hyndman@monash.edu>
#' 
#' @import forecast 
#' @import ftsa 
#' @import rainbow 
#' @import mgcv 
#' @import cobs 
#' @importFrom graphics abline lines plot 
#' @importFrom stats simulate update approx frequency glm lm loess na.omit nlm
#' @importFrom stats poisson predict qnorm rbinom rpois spline splinefun start time 
#' @importFrom stats ts tsp tsp<- uniroot window median 
#' @importFrom utils read.table 
#' 
#' @docType package
#' 
#' @name demography-package
#' @aliases demography
#' @keywords package
NULL
