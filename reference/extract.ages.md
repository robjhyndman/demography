# Extract some ages from a demogdata object

Creates subset of demogdata object.

## Usage

``` r
extract.ages(data, ages, combine.upper = TRUE)
```

## Arguments

- data:

  Demogdata object such as created using
  [`read.demogdata`](read.demogdata.md) or
  [`smooth.demogdata`](smooth.demogdata.md).

- ages:

  Vector of ages to extract from data.

- combine.upper:

  If TRUE, ages beyond the maximum of `ages` are combined into the upper
  age group.

## Value

Demogdata object with same components as `data` but with a subset of
ages.

## Author

Rob J Hyndman

## Examples

``` r
france.teens <- extract.ages(fr.mort, 13:19, FALSE)
plot(france.teens)
```
