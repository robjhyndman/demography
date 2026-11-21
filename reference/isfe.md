# Integrated Squared Forecast Error for models of various orders

Computes ISFE values for functional time series models of various
orders.

## Usage

``` r
isfe(...)

# S3 method for class 'demogdata'
isfe(
  data,
  series = names(data$rate)[1],
  max.order = N - 3,
  N = 10,
  h = 5:10,
  ages = data$age,
  max.age = max(ages),
  method = c("classical", "M", "rapca"),
  fmethod = c("arima", "ar", "arfima", "ets", "ets.na", "struct", "rwdrift", "rw"),
  lambda = 3,
  ...
)
```

## Arguments

- ...:

  Additional arguments control the fitting procedure.

- data:

  demogdata object.

- series:

  name of series within data holding rates (1x1)

- max.order:

  Maximum number of basis functions to fit.

- N:

  Minimum number of functional observations to be used in fitting a
  model.

- h:

  Forecast horizons over which to average.

- ages:

  Ages to include in fit.

- max.age:

  Maximum age to fit.

- method:

  Method to use for principal components decomposition. Possibilities
  are “M”, “rapca” and “classical”.

- fmethod:

  Method used for forecasting. Current possibilities are “ets”, “arima”,
  “ets.na”, “struct”, “rwdrift” and “rw”.

- lambda:

  Tuning parameter for robustness when `method="M"`.

## Value

Numeric matrix with `(max.order+1)` rows and `length(h)` columns
containing ISFE values for models of orders 0:max.order.

## References

Hyndman, R.J., and Ullah, S. (2007) Robust forecasting of mortality and
fertility rates: a functional data approach. *Computational Statistics &
Data Analysis*, **51**, 4942-4956.
<https://robjhyndman.com/publications/funcfor/>

## See also

[`fdm`](fdm.md), [`forecast.fdm`](forecast.fdm.md).

## Author

Rob J Hyndman
