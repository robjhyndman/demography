# Compute total fertility rate from fertility rates

Compute total fertility rates from age-specific fertility rates
contained in a `demogdata` object.

## Usage

``` r
tfr(data, PI = FALSE, nsim = 500, ...)
```

## Arguments

- data:

  Demogdata object of type `"fertility"` such as obtained from
  [`read.demogdata`](read.demogdata.md),
  [`forecast.fdm`](forecast.fdm.md).

- PI:

  If TRUE, produce a prediction interval.

- nsim:

  Number of simulations to use when computing a prediction interval.

- ...:

  Other arguments passed to `simulate` when producing prediction
  intervals.

## Value

If data are of class `demogdata`, the function returns a time series of
fertility rates. If data are from [`forecast.fdm`](forecast.fdm.md), the
function returns an object of class `forecast` containing point
forecasts and (optionally) prediction intervals.

## See also

[`fdm`](fdm.md)

## Author

Rob J Hyndman

## Examples

``` r
plot(tfr(aus.fert))

ausfert.fcast <- forecast(fdm(aus.fert))
plot(tfr(ausfert.fcast, PI = TRUE, nsim = 400))
```
