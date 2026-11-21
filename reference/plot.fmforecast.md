# Plot forecasts from a functional demographic modell

Type of plot depends on value of `plot.type`:

- `plot.type="function"`:

  produces a plot of the forecast functions;

- `plot.type="components"`:

  produces a plot of the basis functions and coefficients with forecasts
  and prediction intervals for each coefficient;

- `plot.type="variance"`:

  produces a plot of the variance components.

## Usage

``` r
# S3 method for class 'fmforecast'
plot(
  x,
  plot.type = c("function", "component", "variance"),
  vcol = 1:4,
  mean.lab = "Mean",
  xlab2 = "Year",
  h = 1,
  ...
)

# S3 method for class 'lca'
plot(x, ...)
```

## Arguments

- x:

  Output from
  [`forecast.ftsm`](https://rdrr.io/pkg/ftsa/man/forecast.ftsm.html),
  [`forecast.fdm`](forecast.fdm.md) or [`lca`](lca.md).

- plot.type:

  Type of plot. See details.

- vcol:

  Colors to use if `plot.type="variance"`.

- mean.lab:

  Label for mean component.

- xlab2:

  x-axis label for coefficient time series.

- h:

  If `plot.type="variance"`, h gives the forecast horizon for which the
  variance is plotted.

- ...:

  Other arguments are passed to [`plot.demogdata`](plot.demogdata.md)
  (if `plot.type=="function"`),
  [`plot`](https://rdrr.io/r/base/plot.html) (if
  `plot.type=="variance"`) or
  [`plot.ftsf`](https://rdrr.io/pkg/ftsa/man/plot.ftsf.html) (if
  `plot.type=="component"`).

## Value

None. Function produces a plot

## See also

[fdm](fdm.md), [lca](lca.md), [forecast.fdm](forecast.fdm.md)

## Author

Rob J Hyndman

## Examples

``` r
france.fcast <- forecast(fdm(fr.mort))
plot(france.fcast)

plot(france.fcast, "c")

plot(france.fcast, "v")
```
