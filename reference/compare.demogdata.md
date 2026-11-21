# Evaluation of demographic forecast accuracy

Computes mean forecast errors and mean square forecast errors for each
age level. Computes integrated squared forecast errors and integrated
absolute percentage forecast errors for each year.

## Usage

``` r
compare.demogdata(
  data,
  forecast,
  series = names(forecast$rate)[1],
  ages = data$age,
  max.age = min(max(data$age), max(forecast$age)),
  years = data$year,
  interpolate = FALSE
)
```

## Arguments

- data:

  Demogdata object such as created using
  [`read.demogdata`](read.demogdata.md) containing actual demographic
  rates.

- forecast:

  Demogdata object such as created using
  [`forecast.fdm`](forecast.fdm.md) or
  [`forecast.lca`](forecast.lca.md).

- series:

  Name of series to use. Default: the first matrix within
  `forecast$rate`.

- ages:

  Ages to use for comparison. Default: all available ages.

- max.age:

  Upper age to use for comparison.

- years:

  Years to use in comparison. Default is to use all available years that
  are common between data and forecast.

- interpolate:

  If TRUE, all zeros in data are replaced by interpolated estimates when
  computing the error measures on the log scale. Error measures on the
  original (rate) scale are unchanged.

## Value

Object of class "errorfdm" with the following components:

- label:

  Name of region from which data taken.

- age:

  Ages from `data` object.

- year:

  Years from `data` object.

- \<error\>:

  Matrix of forecast errors on rates.

- \<logerror\>:

  Matrix of forecast errors on log rates.

- mean.error:

  Various measures of forecast accuracy averaged across years.
  Specifically ME=mean error, MSE=mean squared error, MPE=mean
  percentage error and MAPE=mean absolute percentage error.

- int.error:

  Various measures of forecast accuracy integrated across ages.
  Specifically IE=integrated error, ISE=integrated squared error,
  IPE=integrated percentage error and IAPE=integrated absolute
  percentage error.

- life.expectancy:

  If `data$type="mortality"`, function returns this component which is a
  matrix containing actual, forecast and actual-forecast for life
  expectancies.

Note that the error matrices have different names indicating if the
series forecast was male, female or total.

## See also

[forecast.fdm](forecast.fdm.md),[plot.errorfdm](plot.errorfdm.md)

## Author

Rob J Hyndman

## Examples

``` r
fr.test <- extract.years(fr.sm, years = 1921:1980)
fr.fit <- fdm(fr.test, order = 2)
fr.error <- compare.demogdata(fr.mort, forecast(fr.fit, 20))
#> Warning: No population data available for combining upper ages
plot(fr.error)

par(mfrow = c(2, 1))
plot(fr.error$age, fr.error$mean.error[, "ME"],
  type = "l", xlab = "Age", ylab = "Mean Forecast Error"
)
plot(fr.error$int.error[, "ISE"],
  xlab = "Year", ylab = "Integrated Square Error"
)
```
