# Read data from HMD and construct a mortality demogdata object

`hmd.mx` reads "Mx" (1x1) data from the Human Mortality Database (HMD
<https://www.mortality.org>) and constructs a demogdata object suitable
for plotting using [`plot.demogdata`](plot.demogdata.md) and fitting an
LC or BMS model using [`lca`](lca.md) or an FDA model using
[`fdm`](fdm.md). `hmd.pop` reads "Population" (1x1) data from the HMD
and constructs a demogdata object suitable for plotting using
[`plot.demogdata`](plot.demogdata.md). `hmd.e0` reads life expectancy at
birth from the HMD and returns the result as a `ts` object.

## Usage

``` r
hmd.mx(country, username, password, label = country)

hmd.e0(country, username, password)

hmd.pop(country, username, password, label = country)
```

## Arguments

- country:

  Directory abbreviation from the HMD. For instance, Australia = "AUS".

- username:

  HMD username (case-sensitive)

- password:

  HMD password (case-sensitive)

- label:

  Character string giving name of country from which the data are taken.

## Value

`hmd.mx` returns an object of class `demogdata` with the following
components:

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

`hmd.pop` returns a similar object but without the `rate` component.
`hmd.e0` returns an object of class `ts` with columns `male`, `female`
and `total`.

## Details

In order to read the data, users are required to create their account
via the HMD website (<https://www.mortality.org>), and obtain a valid
username and password.

## See also

[`demogdata`](demogdata.md),[`read.demogdata`](read.demogdata.md),[`plot.demogdata`](plot.demogdata.md),
[`life.expectancy`](life.expectancy.md)

## Author

Rob J Hyndman

## Examples

``` r
if (FALSE) { # \dontrun{
norway <- hmd.mx("NOR", username, password, "Norway")
summary(norway)
} # }
```
