# Updating functional demographic models and coherent functional demographic models.

`update.fmforecast()` updates `fdm` forecasts. The argument `object` is
the output from [`forecast.fdm`](forecast.fdm.md) which has been
subsequently modified with new coefficient forecasts. These new
forecasts are used when re-calculating the forecast of the mortality or
fertility rates, or net migration numbers. `update.fmforecast2()`
updates `fdmpr` forecasts. The argument `object` is the output from
[`forecast.fdmpr`](forecast.fdmpr.md) which has been subsequently
modified with new coefficient forecasts.

## Usage

``` r
# S3 method for class 'fmforecast'
update(object, ...)

# S3 method for class 'fmforecast2'
update(object, ...)
```

## Arguments

- object:

  Output from either [`fdm`](fdm.md) or [`coherentfdm`](coherentfdm.md).

- ...:

  Extra arguments currently ignored.

## Value

A list of the same class as `object`.

## See also

[`forecast.fdm`](forecast.fdm.md), [`forecast.fdmpr`](forecast.fdmpr.md)

## Author

Rob J Hyndman.

## Examples

``` r
if (FALSE) { # \dontrun{
france.fit <- fdm(fr.mort, order = 2)
france.fcast <- forecast(france.fit, 50)
# Replace first coefficient model with ARIMA(0,1,2)+drift
france.fcast$coeff[[2]] <- forecast(Arima(france.fit$coeff[, 2],
  order = c(0, 1, 2), include.drift = TRUE
), h = 50, level = 80)
france.fcast <- update(france.fcast)

fr.short <- extract.years(fr.sm, 1950:2006)
fr.fit <- coherentfdm(fr.short)
fr.fcast <- forecast(fr.fit)
par(mfrow = c(1, 2))
plot(fr.fcast$male)
# Replace first coefficient model in product component with a damped ETS model:
fr.fcast$product$coeff[[2]] <- forecast(ets(fr.fit$product$coeff[, 2], damped = TRUE),
  h = 50, level = 80
)
fr.fcast <- update(fr.fcast)
plot(fr.fcast$male)
} # }
```
