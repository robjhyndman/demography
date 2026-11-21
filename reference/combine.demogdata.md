# Combine two demogdata objects into one demogdata object

Function to combine demogdata objects containing different years but the
same age structure into one demogdata object. The standard use for this
function will be combining historical data with forecasts. The objects
must be of the same type.

## Usage

``` r
combine.demogdata(obj1, obj2)
```

## Arguments

- obj1:

  First demogdata object (e.g., historical data).

- obj2:

  Second demogdata object (e.g., forecasts).

## Value

Object of class “demogdata” with the following components:

- year:

  Vector of years

- age:

  Vector of ages

- rate:

  Matrix of rates with with one age group per row and one column per
  year.

- pop:

  Matrix of populations in same form as `rate` and containing population
  numbers. This is only produced when both objects contain a `pop`
  component.

- type:

  Type of object: “mortality”, “fertility” or “migration”.

- label:

  Name of area from which the data are taken.

## See also

[`demogdata`](demogdata.md)

## Author

Rob J Hyndman

## Examples

``` r
fit <- fdm(fr.mort)
fcast <- forecast(fit, h = 50)
france2 <- combine.demogdata(fr.mort, fcast)
plot(france2)

plot(life.expectancy(france2))
#> Warning: No population data available for combining upper ages
lines(rep(max(fr.mort$year) + 0.5, 2), c(0, 100), lty = 3)
```
