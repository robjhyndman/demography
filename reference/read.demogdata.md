# Read demographic data and construct demogdata object

Read data from text files and construct a demogdata object suitable for
plotting using [`plot.demogdata`](plot.demogdata.md) and fitting an LC
or BMS model using [`lca`](lca.md) or an FDA model using
[`fdm`](fdm.md).

## Usage

``` r
read.demogdata(
  file,
  popfile,
  type,
  label,
  max.mx = 10,
  skip = 2,
  popskip = skip,
  lambda,
  scale = 1
)
```

## Arguments

- file:

  Filename containing demographic rates.

- popfile:

  Filename containing population numbers.

- type:

  Character string showing type of demographic series: either
  “mortality”, “fertility” or “migration”.

- label:

  Name of area from which the data are taken.

- max.mx:

  Maximum allowable value for demographic rate. All values greater than
  max.mx will be set to max.mx.

- skip:

  Number of lines to skip at the start of `file`.

- popskip:

  Number of lines to skip at the start of `popfile`.

- lambda:

  Box-Cox transformation parameter to be used in modelling and plotting.
  If missing, default values are 0 (for mortality), 0.4 (for fertility)
  and 1 (for migration).

- scale:

  Number of people in the rate definition. `scale=1` indicates the rates
  are per person; `scale=1000` indicates the rates are per 1000 people.

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

## Details

All data are assumed to be tab-delimited text files with the first
column containing the year of observation and the second column
containing the age level. All remaining columns are assumed to be
demographic rates for sections of the population. The first row of the
text file is assumed to contain the names of each column. Population
data are assumed to have the same format but with population numbers in
place of rates. The columns names in the two files should be identical.
Note that this format is what is used by the Human Mortality Database
<http://www.mortality.org>. If `popfile` contains the Exposures and
`file` contains the Mx rates from the HMD, then everything will work
seamlessly.

## See also

[`demogdata`](demogdata.md)

## Author

Rob J Hyndman

## Examples

``` r
  if (FALSE) norway <- read.demogdata("Mx_1x1.txt",
  "Exposures_1x1.txt", type="mortality", label="Norway") # \dontrun{}
```
