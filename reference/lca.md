# Model mortality or fertility data using Lee-Carter approach

Lee-Carter model of mortality or fertility rates. `lca` produces a
standard Lee-Carter model by default, although many other options are
available. `bms` is a wrapper for `lca` and returns a model based on the
Booth-Maindonald-Smith methodology.

## Usage

``` r
lca(
  data,
  series = names(data$rate)[1],
  years = data$year,
  ages = data$age,
  max.age = 100,
  adjust = c("dt", "dxt", "e0", "none"),
  chooseperiod = FALSE,
  minperiod = 20,
  breakmethod = c("bai", "bms"),
  scale = FALSE,
  restype = c("logrates", "rates", "deaths"),
  interpolate = FALSE
)

bms(
  data,
  series = names(data$rate)[1],
  years = data$year,
  ages = data$age,
  max.age = 100,
  minperiod = 20,
  breakmethod = c("bms", "bai"),
  scale = FALSE,
  restype = c("logrates", "rates", "deaths"),
  interpolate = FALSE
)
```

## Arguments

- data:

  demogdata object of type “mortality” or “fertility”. Output from
  read.demogdata.

- series:

  name of series within data containing mortality or fertility values
  (1x1)

- years:

  years to include in fit. Default: all available years.

- ages:

  ages to include in fit. Default: all available ages up to `max.age`.

- max.age:

  upper age to include in fit. Ages beyond this are collapsed into the
  upper age group.

- adjust:

  method to use for adjustment of coefficients \\k_t kt\\. Possibilities
  are “dxt” (BMS method), “dt” (Lee-Carter method), “e0” (method based
  on life expectancy) and “none”. Defaults are “dxt” for `bms()` and
  “dt” for `lca()`.

- chooseperiod:

  If TRUE, it will choose the best fitting period.

- minperiod:

  Minimum number of years to include in fitting period if
  chooseperiod=TRUE.

- breakmethod:

  method to use for identifying breakpoints if chooseperiod=TRUE.
  Possibilities are “bai” (Bai's method computed using
  [`breakpoints`](https://rdrr.io/pkg/strucchange/man/breakpoints.html)
  in the strucchange package) and “bms” (method based on mean deviance
  ratios described in BMS).

- scale:

  If TRUE, it will rescale bx and kt so that kt has drift parameter = 1.

- restype:

  method to use for calculating residuals. Possibilities are “logrates”,
  “rates” and “deaths”.

- interpolate:

  If TRUE, it will estimate any zero mortality or fertility rates using
  the same age group from nearby years.

## Value

Object of class “lca” with the following components:

- label:

  Name of region

- age:

  Ages from `data` object.

- year:

  Years from `data` object.

- \<series\>:

  Matrix of mortality or fertility data as contained in `data`. It takes
  the name given by the series argument.

- ax:

  Average deathrates across fitting period

- bx:

  First principal component in Lee-Carter model

- kt:

  Coefficient of first principal component

- residuals:

  Functional time series of residuals.

- fitted:

  Functional time series containing estimated mortality or fertility
  rates from model

- varprop:

  Proportion of variance explained by model.

- y:

  The data stored as a functional time series object.

- mdev:

  Mean deviance of total and base lack of fit, as described in Booth,
  Maindonald and Smith.

## Details

All mortality or fertility data are assumed to be in matrices of
mortality or fertility rates within `data$rate`. Each row is one age
group (assumed to be single years). Each column is one year. The
function produces a model for the `series` mortality or fertility rate
matrix within `data$rate`. Forecasts from this model can be obtained
using [`forecast.lca`](forecast.lca.md).

## References

Booth, H., Maindonald, J., and Smith, L. (2002) Applying Lee-Carter
under conditions of variable mortality decline. *Population Studies*,
**56**, 325-336.

Lee, R.D., and Carter, L.R. (1992) Modeling and forecasting US
mortality. *Journal of the American Statistical Association*, **87**,
659-671.

## See also

[`forecast.lca`](forecast.lca.md), [`plot.lca`](plot.fmforecast.md),
[`summary.lca`](summary.fdm.md), [`fdm`](fdm.md)

## Author

Heather Booth, Leonie Tickle, John Maindonald and Rob J Hyndman.

## Examples

``` r
if (FALSE) { # \dontrun{
france.LC1 <- lca(fr.mort, adjust = "e0")
plot(france.LC1)
par(mfrow = c(1, 2))
plot(fr.mort, years = 1953:2002, ylim = c(-11, 1))
plot(forecast(france.LC1, jumpchoice = "actual"), ylim = c(-11, 1))

france.bms <- bms(fr.mort, breakmethod = "bai")
fcast.bms <- forecast(france.bms)
par(mfrow = c(1, 1))
plot(fcast.bms$kt)
} # }
```
