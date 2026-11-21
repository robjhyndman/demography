# Monotonic interpolating splines

Perform cubic spline monotonic interpolation of given data points,
returning either a list of points obtained by the interpolation or a
function performing the interpolation. The splines are constrained to be
monotonically increasing (i.e., the slope is never negative).

## Usage

``` r
cm.spline(x, y = NULL, n = 3 * length(x), xmin = min(x), xmax = max(x), ...)

cm.splinefun(x, y = NULL, ...)
```

## Arguments

- x, y:

  vectors giving the coordinates of the points to be interpolated.
  Alternatively a single plotting structure can be specified: see
  [`xy.coords`](https://rdrr.io/r/grDevices/xy.coords.html).

- n:

  interpolation takes place at n equally spaced points spanning the
  interval \[`xmin`, `xmax`\].

- xmin:

  left-hand endpoint of the interpolation interval.

- xmax:

  right-hand endpoint of the interpolation interval.

- ...:

  Other arguments are ignored.

## Value

- cm.spline:

  returns a list containing components `x` and `y` which give the
  ordinates where interpolation took place and the interpolated values.

- cm.splinefun:

  returns a function which will perform cubic spline interpolation of
  the given data points. This is often more useful than `spline`.

## Details

These are simply wrappers to the
[`splinefun`](https://rdrr.io/r/stats/splinefun.html) function family
from the stats package.

## References

Forsythe, G. E., Malcolm, M. A. and Moler, C. B. (1977) *Computer
Methods for Mathematical Computations*. Hyman (1983) *SIAM J. Sci. Stat.
Comput.* **4**(4):645-654. Dougherty, Edelman and Hyman 1989
*Mathematics of Computation*, **52**: 471-494.

## Author

Rob J Hyndman

## Examples

``` r
x <- seq(0, 4, l = 20)
y <- sort(rnorm(20))
plot(x, y)
lines(spline(x, y, n = 201), col = 2) # Not necessarily monotonic
lines(cm.spline(x, y, n = 201), col = 3) # Monotonic
```
