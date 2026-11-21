# Extract some years from a demogdata object

Creates subset of demogdata object.

## Usage

``` r
extract.years(data, years)
```

## Arguments

- data:

  Demogdata object such as created using
  [`read.demogdata`](read.demogdata.md) or
  [`smooth.demogdata`](smooth.demogdata.md).

- years:

  Vector of years to extract from data.

## Value

Demogdata object with same components as `data` but with a subset of
years.

## Author

Rob J Hyndman

## Examples

``` r
france.1918 <- extract.years(fr.mort, 1918)
```
