# Forecast coherent functional demographic model.

The product and ratio models from [`coherentfdm`](coherentfdm.md) are
forecast, and the results combined to give forecasts for each group in
the original data.

## Usage

``` r
# S3 method for class 'fdmpr'
forecast(object, h = 50, level = 80, K = 100, drange = c(0, 0.5), ...)
```

## Arguments

- object:

  Output from [`coherentfdm`](coherentfdm.md).

- h:

  Forecast horizon.

- level:

  Confidence level for prediction intervals.

- K:

  Maximum number of years to use in forecasting coefficients for ratio
  components.

- drange:

  Range of fractional differencing parameter for the ratio coefficients.

- ...:

  Other arguments as for [`forecast.fdm`](forecast.fdm.md).

## Value

Object of class `fmforecast2` containing a list of objects each of class
`fmforecast`. The forecasts for each group in the original data are
given first. Then the forecasts from the product model, and finally a
list of forecasts from each of the ratio models.

## See also

[`coherentfdm`](coherentfdm.md), [`forecast.fdm`](forecast.fdm.md).

## Author

Rob J Hyndman

## Examples

``` r
fr.short <- extract.years(fr.sm, 1950:2006)
fr.fit <- coherentfdm(fr.short)
fr.fcast <- forecast(fr.fit)
plot(fr.fcast$male)

plot(fr.fcast$ratio$male, plot.type = "component", components = 3)

models(fr.fcast)
#> 
#> ************* PRODUCT MODEL *************
#> 
#> -- Coefficient 1 --
#> Series: xx[, i] 
#> ARIMA(1,1,0) with drift 
#> 
#> Coefficients:
#>           ar1    drift
#>       -0.4250  -0.2098
#> s.e.   0.1206   0.0228
#> 
#> sigma^2 = 0.06039:  log likelihood = 0.05
#> AIC=5.89   AICc=6.36   BIC=11.97
#> 
#> -- Coefficient 2 --
#> Series: xx[, i] 
#> ARIMA(0,1,0) 
#> 
#> sigma^2 = 0.01882:  log likelihood = 31.78
#> AIC=-61.56   AICc=-61.48   BIC=-59.53
#> 
#> -- Coefficient 3 --
#> Series: xx[, i] 
#> ARIMA(1,0,0) with zero mean 
#> 
#> Coefficients:
#>          ar1
#>       0.9518
#> s.e.  0.0449
#> 
#> sigma^2 = 0.02409:  log likelihood = 24.63
#> AIC=-45.27   AICc=-45.04   BIC=-41.18
#> 
#> -- Coefficient 4 --
#> Series: xx[, i] 
#> ARIMA(0,0,0) with zero mean 
#> 
#> sigma^2 = 0.03996:  log likelihood = 10.89
#> AIC=-19.78   AICc=-19.71   BIC=-17.74
#> 
#> -- Coefficient 5 --
#> Series: xx[, i] 
#> ARIMA(1,0,2) with zero mean 
#> 
#> Coefficients:
#>          ar1      ma1     ma2
#>       0.8425  -0.6319  0.3379
#> s.e.  0.1173   0.1584  0.1432
#> 
#> sigma^2 = 0.01901:  log likelihood = 33.08
#> AIC=-58.15   AICc=-57.38   BIC=-49.98
#> 
#> -- Coefficient 6 --
#> Series: xx[, i] 
#> ARIMA(2,0,2) with zero mean 
#> 
#> Coefficients:
#>          ar1      ar2      ma1     ma2
#>       1.8493  -0.9078  -1.3057  0.4781
#> s.e.  0.0705   0.0663   0.1613  0.1487
#> 
#> sigma^2 = 0.005993:  log likelihood = 66.08
#> AIC=-122.16   AICc=-120.99   BIC=-111.95
#> 
#> 
#> 
#> ************* RATIO MODELS *************
#> 
#> 
#> *********** FEMALE ***********
#> 
#> 
#> -- Coefficient 1 --
#> 
#> Call:
#>   arfima(y = xx[, i], drange = ..2, estim = "mle") 
#> 
#> *** Warning during (fdcov) fit: unable to compute correlation matrix; maybe change 'h'
#> 
#> Coefficients:
#>          d     ma.ma1     ma.ma2 
#>  0.4090936 -1.1513350 -0.5306586 
#> sigma[eps] = 0.1303427 
#> a list with components:
#>  [1] "log.likelihood"  "n"               "msg"             "d"              
#>  [5] "ar"              "ma"              "covariance.dpq"  "fnormMin"       
#>  [9] "sigma"           "stderror.dpq"    "correlation.dpq" "h"              
#> [13] "d.tol"           "M"               "hessian.dpq"     "length.w"       
#> [17] "residuals"       "fitted"          "call"            "x"              
#> [21] "series"         
#> 
#> -- Coefficient 2 --
#> 
#> Call:
#>   arfima(y = xx[, i], drange = ..2, estim = "mle") 
#> 
#> *** Warning during (fdcov) fit: unable to compute correlation matrix; maybe change 'h'
#> 
#> Coefficients:
#>             d        ma.ma1        ma.ma2 
#>  4.583013e-05 -1.622329e-01  5.425798e-01 
#> sigma[eps] = 0.1459614 
#> a list with components:
#>  [1] "log.likelihood"  "n"               "msg"             "d"              
#>  [5] "ar"              "ma"              "covariance.dpq"  "fnormMin"       
#>  [9] "sigma"           "stderror.dpq"    "correlation.dpq" "h"              
#> [13] "d.tol"           "M"               "hessian.dpq"     "length.w"       
#> [17] "residuals"       "fitted"          "call"            "x"              
#> [21] "series"         
#> 
#> -- Coefficient 3 --
#> 
#> Call:
#>   arfima(y = xx[, i], drange = ..2, estim = "mle") 
#> 
#> *** Warning during (fdcov) fit: unable to compute correlation matrix; maybe change 'h'
#> 
#> Coefficients:
#>            d       ar.ar1       ma.ma1 
#> 4.583013e-05 9.798979e-01 4.156993e-01 
#> sigma[eps] = 0.06410902 
#> a list with components:
#>  [1] "log.likelihood"  "n"               "msg"             "d"              
#>  [5] "ar"              "ma"              "covariance.dpq"  "fnormMin"       
#>  [9] "sigma"           "stderror.dpq"    "correlation.dpq" "h"              
#> [13] "d.tol"           "M"               "hessian.dpq"     "length.w"       
#> [17] "residuals"       "fitted"          "call"            "x"              
#> [21] "series"         
#> 
#> -- Coefficient 4 --
#> 
#> Call:
#>   arfima(y = xx[, i], drange = ..2, estim = "mle") 
#> 
#> *** Warning during (fdcov) fit: unable to compute correlation matrix; maybe change 'h'
#> 
#> Coefficients:
#>          d     ar.ar1     ar.ar2     ma.ma1     ma.ma2 
#>  0.3913578  1.6048157 -0.7526423  1.6155005 -0.9998107 
#> sigma[eps] = 0.06159422 
#> a list with components:
#>  [1] "log.likelihood"  "n"               "msg"             "d"              
#>  [5] "ar"              "ma"              "covariance.dpq"  "fnormMin"       
#>  [9] "sigma"           "stderror.dpq"    "correlation.dpq" "h"              
#> [13] "d.tol"           "M"               "hessian.dpq"     "length.w"       
#> [17] "residuals"       "fitted"          "call"            "x"              
#> [21] "series"         
#> 
#> -- Coefficient 5 --
#> 
#> Call:
#>   arfima(y = xx[, i], drange = ..2, estim = "mle") 
#> 
#> Coefficients:
#>            d       ar.ar1 
#> 4.583013e-05 7.006761e-01 
#> sigma[eps] = 0.06151039 
#> a list with components:
#>  [1] "log.likelihood"  "n"               "msg"             "d"              
#>  [5] "ar"              "ma"              "covariance.dpq"  "fnormMin"       
#>  [9] "sigma"           "stderror.dpq"    "correlation.dpq" "h"              
#> [13] "d.tol"           "M"               "hessian.dpq"     "length.w"       
#> [17] "residuals"       "fitted"          "call"            "x"              
#> [21] "series"         
#> 
#> -- Coefficient 6 --
#> 
#> Call:
#>   arfima(y = xx[, i], drange = ..2, estim = "mle") 
#> 
#> *** Warning during (fdcov) fit: unable to compute correlation matrix; maybe change 'h'
#> 
#> Coefficients:
#>            d       ar.ar1       ma.ma1 
#> 4.583013e-05 8.123545e-01 3.577327e-01 
#> sigma[eps] = 0.05375933 
#> a list with components:
#>  [1] "log.likelihood"  "n"               "msg"             "d"              
#>  [5] "ar"              "ma"              "covariance.dpq"  "fnormMin"       
#>  [9] "sigma"           "stderror.dpq"    "correlation.dpq" "h"              
#> [13] "d.tol"           "M"               "hessian.dpq"     "length.w"       
#> [17] "residuals"       "fitted"          "call"            "x"              
#> [21] "series"         
#> 
#> 
#> 
#> *********** MALE ***********
#> 
#> 
#> -- Coefficient 1 --
#> 
#> Call:
#>   arfima(y = xx[, i], drange = ..2, estim = "mle") 
#> 
#> *** Warning during (fdcov) fit: unable to compute correlation matrix; maybe change 'h'
#> 
#> Coefficients:
#>          d     ma.ma1     ma.ma2 
#>  0.4236553 -1.1435894 -0.5225683 
#> sigma[eps] = 0.1209133 
#> a list with components:
#>  [1] "log.likelihood"  "n"               "msg"             "d"              
#>  [5] "ar"              "ma"              "covariance.dpq"  "fnormMin"       
#>  [9] "sigma"           "stderror.dpq"    "correlation.dpq" "h"              
#> [13] "d.tol"           "M"               "hessian.dpq"     "length.w"       
#> [17] "residuals"       "fitted"          "call"            "x"              
#> [21] "series"         
#> 
#> -- Coefficient 2 --
#> 
#> Call:
#>   arfima(y = xx[, i], drange = ..2, estim = "mle") 
#> 
#> *** Warning during (fdcov) fit: unable to compute correlation matrix; maybe change 'h'
#> 
#> Coefficients:
#>             d        ma.ma1        ma.ma2 
#>  4.583013e-05 -1.622329e-01  5.425798e-01 
#> sigma[eps] = 0.1425088 
#> a list with components:
#>  [1] "log.likelihood"  "n"               "msg"             "d"              
#>  [5] "ar"              "ma"              "covariance.dpq"  "fnormMin"       
#>  [9] "sigma"           "stderror.dpq"    "correlation.dpq" "h"              
#> [13] "d.tol"           "M"               "hessian.dpq"     "length.w"       
#> [17] "residuals"       "fitted"          "call"            "x"              
#> [21] "series"         
#> 
#> -- Coefficient 3 --
#> 
#> Call:
#>   arfima(y = xx[, i], drange = ..2, estim = "mle") 
#> 
#> *** Warning during (fdcov) fit: unable to compute correlation matrix; maybe change 'h'
#> 
#> Coefficients:
#>            d       ar.ar1       ma.ma1 
#> 4.583013e-05 9.798979e-01 4.156993e-01 
#> sigma[eps] = 0.06414137 
#> a list with components:
#>  [1] "log.likelihood"  "n"               "msg"             "d"              
#>  [5] "ar"              "ma"              "covariance.dpq"  "fnormMin"       
#>  [9] "sigma"           "stderror.dpq"    "correlation.dpq" "h"              
#> [13] "d.tol"           "M"               "hessian.dpq"     "length.w"       
#> [17] "residuals"       "fitted"          "call"            "x"              
#> [21] "series"         
#> 
#> -- Coefficient 4 --
#> 
#> Call:
#>   arfima(y = xx[, i], drange = ..2, estim = "mle") 
#> 
#> Coefficients:
#>          d     ar.ar1     ar.ar2     ma.ma1     ma.ma2 
#>  0.3797002  1.6131684 -0.7587267  1.6171315 -0.9999822 
#> sigma[eps] = 0.06278872 
#> a list with components:
#>  [1] "log.likelihood"  "n"               "msg"             "d"              
#>  [5] "ar"              "ma"              "covariance.dpq"  "fnormMin"       
#>  [9] "sigma"           "stderror.dpq"    "correlation.dpq" "h"              
#> [13] "d.tol"           "M"               "hessian.dpq"     "length.w"       
#> [17] "residuals"       "fitted"          "call"            "x"              
#> [21] "series"         
#> 
#> -- Coefficient 5 --
#> 
#> Call:
#>   arfima(y = xx[, i], drange = ..2, estim = "mle") 
#> 
#> Coefficients:
#>            d       ar.ar1 
#> 4.583013e-05 7.006761e-01 
#> sigma[eps] = 0.06202872 
#> a list with components:
#>  [1] "log.likelihood"  "n"               "msg"             "d"              
#>  [5] "ar"              "ma"              "covariance.dpq"  "fnormMin"       
#>  [9] "sigma"           "stderror.dpq"    "correlation.dpq" "h"              
#> [13] "d.tol"           "M"               "hessian.dpq"     "length.w"       
#> [17] "residuals"       "fitted"          "call"            "x"              
#> [21] "series"         
#> 
#> -- Coefficient 6 --
#> 
#> Call:
#>   arfima(y = xx[, i], drange = ..2, estim = "mle") 
#> 
#> *** Warning during (fdcov) fit: unable to compute correlation matrix; maybe change 'h'
#> 
#> Coefficients:
#>            d       ar.ar1       ma.ma1 
#> 4.583013e-05 8.123545e-01 3.577327e-01 
#> sigma[eps] = 0.05376532 
#> a list with components:
#>  [1] "log.likelihood"  "n"               "msg"             "d"              
#>  [5] "ar"              "ma"              "covariance.dpq"  "fnormMin"       
#>  [9] "sigma"           "stderror.dpq"    "correlation.dpq" "h"              
#> [13] "d.tol"           "M"               "hessian.dpq"     "length.w"       
#> [17] "residuals"       "fitted"          "call"            "x"              
#> [21] "series"         
```
