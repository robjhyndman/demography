# Create demogdata object from raw data matrices

Create demogdata object suitable for plotting using
[`plot.demogdata`](plot.demogdata.md) and fitting an LC or BMS model
using [`lca`](lca.md) or an FDA model using [`fdm`](fdm.md).

## Usage

``` r
demogdata(data, pop, ages, years, type, label, name, lambda)
```

## Arguments

- data:

  Matrix of data: either mortality rates or fertility rates

- pop:

  Matrix of population values of same dimension as data. These are
  population numbers as at 30 June of each year (i.e., the "exposures").
  So, for example, the number of deaths is data\*pop if data contains
  mortality rates.

- ages:

  Vector of ages corresponding to rows of `data`.

- years:

  Vector of years corresponding to columns of `data`.

- type:

  Character string showing type of demographic series: either
  “mortality”, “fertility” or “migration”.

- label:

  Character string of the name of area from which the data are taken.

- name:

  Name of series: usually male, female or total.

- lambda:

  Box-Cox transformation parameter.

## Value

Object of class “demogdata” with the following components:

- year:

  Vector of years

- age:

  Vector of ages

- rate:

  A list containing one or more rate matrices with one age group per row
  and one column per year.

- pop:

  A list of the same form as `rate` but containing population numbers
  instead of demographic rates.

- type:

  Type of object: “mortality”, “fertility” or “migration”.

- label:

  label

- lambda:

  lambda

## See also

[`read.demogdata`](read.demogdata.md)

## Author

Rob J Hyndman
