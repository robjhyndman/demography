# Forecast demogdata data using Lee-Carter method.

The kt coefficients are forecast using a random walk with drift. The
forecast coefficients are then multiplied by bx to obtain a forecast
demographic rate curve.

## Usage

``` r
# S3 method for class 'lca'
forecast(
  object,
  h = 50,
  se = c("innovdrift", "innovonly"),
  jumpchoice = c("fit", "actual"),
  level = 80,
  ...
)
```

## Arguments

- object:

  Output from [`lca`](lca.md).

- h:

  Number of years ahead to forecast.

- se:

  Method used for computation of standard error. Possibilities:
  “innovdrift” (innovations and drift) and “innovonly” (innovations
  only).

- jumpchoice:

  Method used for computation of jumpchoice. Possibilities: “actual”
  (use actual rates from final year) and “fit” (use fitted rates).

- level:

  Confidence level for prediction intervals.

- ...:

  Other arguments.

## Value

Object of class `fmforecast` with the following components:

- label:

  Region from which the data are taken.

- age:

  Ages from `object`.

- year:

  Years from `object`.

- rate:

  List of matrices containing forecasts, lower bound and upper bound of
  prediction intervals. Point forecast matrix takes the same name as the
  series that has been forecast.

- fitted:

  Matrix of one-step forecasts for historical data

Other components included are

- e0:

  Forecasts of life expectancies (including lower and upper bounds)

- kt.f:

  Forecasts of coefficients from the model.

- type:

  Data type.

- model:

  Details about the fitted model

## See also

[`lca`](lca.md), [`plot.fmforecast`](plot.fmforecast.md)

## Author

Rob J Hyndman

## Examples

``` r
france.lca <- lca(fr.mort, adjust = "e0")
france.fcast <- forecast(france.lca, 50)
plot(france.fcast)

plot(france.fcast, "c")
```
