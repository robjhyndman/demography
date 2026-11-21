# Mean and median functions for data of class demogdata

Computes mean or median of demographic rates for each age level.

## Usage

``` r
# S3 method for class 'demogdata'
mean(x, series = names(x$rate)[1], transform = TRUE, na.rm = TRUE, ...)

# S3 method for class 'demogdata'
median(
  x,
  na.rm = FALSE,
  series = names(x$rate)[1],
  transform = TRUE,
  method = c("hossjercroux", "coordinate"),
  ...
)
```

## Arguments

- x:

  Demogdata object such as created using
  [`read.demogdata`](read.demogdata.md) or
  [`smooth.demogdata`](smooth.demogdata.md).

- series:

  Name of demogdata series to plot..

- transform:

  Should transform of data be taken first?

- na.rm:

  a logical value indicating whether NA values should be stripped before
  the computation proceeds.

- ...:

  Other arguments.

- method:

  Method for computing the median. Either "coordinate" for a
  coordinate-wise median, or "hossjercroux" for the L1-median using the
  Hossjer-Croux algorithm.

## Value

A list containing `x`=ages and `y`=mean or median rates.

## References

Hossjer, O., and Croux, C. (1995) Generalized univariate signed rank
statistics for testing and estimating a multivariate location parameter.
*Nonparametric Statistics*, **4**, 293-308.

## Author

Rob J Hyndman

## Examples

``` r
plot(fr.mort)
lines(mean(fr.mort), lwd = 2)
lines(median(fr.mort), lwd = 2, col = 2)
```
