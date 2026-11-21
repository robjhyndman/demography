# Compute sex ratios from mortality rates

Calculates the Male/Female ratios from historical or forecasted
mortality rates.

## Usage

``` r
sex.ratio(data)
```

## Arguments

- data:

  Demogdata object of type “mortality” such as obtained from
  [`read.demogdata`](read.demogdata.md), or an object of class
  `fmforecast` such as the output from [`forecast.fdm`](forecast.fdm.md)
  or [`forecast.lca`](forecast.lca.md).

## Value

Functional time series of sex ratios.

## Author

Rob J Hyndman

## Examples

``` r
plot(sex.ratio(fr.mort), ylab = "Sex ratios (M/F)")
```
