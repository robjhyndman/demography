# Combine the upperages of a demogdata object.

Computes demographic rates by combining age groups.

## Usage

``` r
set.upperage(data, max.age)
```

## Arguments

- data:

  Demogdata object such as created using
  [`read.demogdata`](read.demogdata.md) or
  [`smooth.demogdata`](smooth.demogdata.md).

- max.age:

  Upper age group. Ages beyond this are combined into the upper age
  group.

## Value

Demogdata object with same components as `data` but with a subset of
ages.

## Author

Rob J Hyndman

## Examples

``` r
france.short <- set.upperage(fr.mort, 85)
```
