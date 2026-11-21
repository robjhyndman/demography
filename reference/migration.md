# Calculate net migration from mortality and fertility data

Function to compute the net number of migrants in each year and for each
age, based on the total population numbers, deaths and births in each
year.

## Usage

``` r
netmigration(mort, fert, startyearpop=mort, mfratio = 1.05)
```

## Arguments

- mort:

  Demogdata object of type `"mortality"`

- fert:

  Demogdata object of type `"fertility"`

- startyearpop:

  Demogdata object containing population data for first year of
  calculation.

- mfratio:

  Male-female ratio to be used in simulating births.

## Value

Object of class “demogdata” with the following components:

- year:

  Vector of years

- age:

  Vector of ages

- rate:

  List containing matrices of net migration numbers (not "rates") with
  with one age group per row and one column per year. Names of matrices
  are the same as for `mort$rate`.

- pop:

  List containing matrices of populations in same form as `rate` and
  containing population numbers.

- type:

  Type of object. In this case, “migration”.

- label:

  label from `mort$label`

## See also

[`demogdata`](demogdata.md)

## Author

Rob J Hyndman

## Examples

``` r
if (FALSE) { # \dontrun{
require(addb)
aus.mig <- netmigration(australia,aus.fertility)
plot(aus.mig)} # }
```
