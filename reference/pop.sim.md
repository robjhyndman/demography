# Population simulation

Simulate future sample paths of a population using functional models for
mortality, fertility and migration.

## Usage

``` r
pop.sim(
  mort,
  fert = NULL,
  mig = NULL,
  firstyearpop,
  N = 100,
  mfratio = 1.05,
  bootstrap = FALSE
)
```

## Arguments

- mort:

  Forecasts of class `fmforecast2` for mortality.

- fert:

  Forecasts of class `fmforecast` for female fertility.

- mig:

  Forecasts of class `fmforecast2` for net migration.

- firstyearpop:

  Population for first year of simulation.

- N:

  Number of sample paths to simulate.

- mfratio:

  Male-female ratio used in distributing births.

- bootstrap:

  If TRUE, simulation uses resampled errors rather than normally
  distributed errors.

## Value

A list of two arrays containing male and female future simulated
population values. The arrays are of dimension (p,h,N) where p is the
number of age groups, h is the forecast horizon and N is the number of
simulated sample paths.

## See also

[`simulate.fmforecast`](simulate.fmforecast.md),
[`simulate.fmforecast2`](simulate.fmforecast.md).

## Author

Rob J Hyndman

## Examples

``` r
if (FALSE) { # \dontrun{
require(addb)
# Construct data objects
mort.sm <- smooth.demogdata(set.upperage(extract.years(australia, 1950:2002), 100))
fert.sm <- smooth.demogdata(extract.years(aus.fertility, 1950:2002))
aus.mig <- netmigration(set.upperage(australia, 100), aus.fertility, mfratio = 1.0545)
# Fit models
mort.fit <- coherentfdm(mort.sm)
fert.fit <- fdm(fert.sm)
mig.fit <- coherentfdm(aus.mig)
# Produce forecasts
mort.fcast <- forecast(mort.fit)
fert.fcast <- forecast(fert.fit)
mig.fcast <- forecast(mig.fit)
# Simulate
aus.sim <- pop.sim(mort.fcast, fert.fcast, mig.fcast, australia)
} # }
```
