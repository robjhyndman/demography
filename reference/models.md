# Show model information for the forecast coefficients in FDM models.

The models for the time series coefficients used in forecasting fdm
models are shown.

## Usage

``` r
models(object, ...)

# S3 method for class 'fmforecast'
models(object, select = 0, ...)

# S3 method for class 'fmforecast2'
models(object, ...)
```

## Arguments

- object:

  Output from [`forecast.fdm`](forecast.fdm.md) or
  [`forecast.fdmpr`](forecast.fdmpr.md).

- ...:

  Other arguments.

- select:

  Indexes of coefficients to display. If select=0, all coefficients are
  displayed.

## See also

[`forecast.fdm`](forecast.fdm.md),
[`forecast.fdmpr`](forecast.fdmpr.md).

## Author

Rob J Hyndman

## Examples

``` r
if (FALSE) { # \dontrun{
fr.short <- extract.years(fr.sm, 1950:2006)
fr.fit <- fdm(fr.short, series = "male")
fr.fcast <- forecast(fr.fit)
models(fr.fcast)

fr.fit <- coherentfdm(fr.short)
fr.fcast <- forecast(fr.fit)
models(fr.fcast, select = 1:3)
} # }
```
