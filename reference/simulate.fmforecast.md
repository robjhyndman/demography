# Simulate future sample paths from functional demographic model forecasts.

This function will simulate future sample paths given forecasting models
from a functional demographic model such as those obtained using
[`forecast.fdm`](forecast.fdm.md) or
[`forecast.fdmpr`](forecast.fdmpr.md).

## Usage

``` r
# S3 method for class 'fmforecast'
simulate(
  object,
  nsim = 100,
  seed = NULL,
  bootstrap = FALSE,
  adjust.modelvar = TRUE,
  ...
)

# S3 method for class 'fmforecast2'
simulate(object, ...)
```

## Arguments

- object:

  Object of class `fmforecast`. Typically, this is output from
  [`forecast.fdm`](forecast.fdm.md).

- nsim:

  Number of sample paths to simulate.

- seed:

  Either NULL or an integer that will be used in a call to set.seed
  before simulating the time seriers. The default, NULL will not change
  the random generator state.

- bootstrap:

  If TRUE, simulation uses resampled errors rather than normally
  distributed errors.

- adjust.modelvar:

  If TRUE, will adjust the model variance by the ratio of the empirical
  and theoretical variances for one-step forecasts.

- ...:

  Other arguments passed to `simulate.fmforecast`.

## Value

An array containing the future simulated values (in the case of a
`fmforecast` object), or a list of arrays containing the future
simulated values (in the case of a `fmforecast2` object).

## See also

[`forecast.fdm`](forecast.fdm.md), [`forecast.lca`](forecast.lca.md),
[`forecast.ftsm`](https://rdrr.io/pkg/ftsa/man/forecast.ftsm.html).

## Author

Rob J Hyndman

## Examples

``` r
if (FALSE) { # \dontrun{
france.fit <- fdm(fr.mort, order = 2)
france.fcast <- forecast(france.fit, 50, method = "ets")
france.sim <- simulate(france.fcast, nsim = 100)

france.fit2 <- coherentfdm(fr.sm)
france.fcast2 <- forecast(france.fit2, 50)
france.sim2 <- simulate(france.fcast2, nsim = 100)
} # }
```
