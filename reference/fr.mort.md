# French mortality data

Age-specific mortality rates and population for France.

## Details

`fr.mort` contains French mortality rates and populations (1899-2005)
for ages 0-110. Data taken from the Human Mortality Database on 20
February 2008. `fr.sm` contains a smoothed version of `fr.mort` obtained
using the [`smooth.demogdata`](smooth.demogdata.md) function.

## Format

Object of class `demogdata` containing the following components:

- year:

  Vector of years

- age:

  Vector of ages

- rate:

  List of matrices containing rates with with one age group per row and
  one column per year. Matrices: `total`, `female`, `male`.

- pop:

  Population data in same form as `rate`.

- type:

  Type of object. In this case, “mortality”.

- label:

  Character string giving area from which data are taken. In this case,
  “France”.

## Source

The Human Mortality Database (<http://www.mortality.org>).

## Author

Rob J Hyndman

## Examples

``` r
plot(fr.mort,years=1950:1997)


plot(fr.mort,years=1990,type='p',pch=1)
lines(fr.sm,years=1990)
```
