# Coherent functional demographic model for grouped data

Fits a coherent functional model to demographic data as described in
Hyndman, Booth & Yasmeen (2012). If two of the series in `data` are
named `male` and `female`, then it will use these two groups. Otherwise
it will use all available groups.

## Usage

``` r
coherentfdm(data, order1 = 6, order2 = 6, ...)
```

## Arguments

- data:

  demogdata object containing at least two groups.

- order1:

  Number of basis functions to fit to the model for the geometric mean.

- order2:

  Number of basis functions to fit to the models for each ratio.

- ...:

  Extra arguments passed to [`fdm`](fdm.md).

## Value

A list (of class `fdmpr`) consisting of two objects: `product` (an
[`fdm`](fdm.md) object containing a del for the geometric mean of the
data) and `ratio` (a list of [`fdm`](fdm.md) objects, being the models
for the ratio of each series with the geometric mean).

## References

Hyndman, R.J., Booth, H., and Yasmeen, F. (2012) Coherent mortality
forecasting: the product-ratio method with functional time series
models. *Demography*, to appear.
<https://robjhyndman.com/publications/coherentfdm/>

## See also

[`fdm`](fdm.md), [`forecast.fdmpr`](forecast.fdmpr.md)

## Author

Rob J Hyndman

## Examples

``` r
fr.short <- extract.years(fr.sm, 1950:2006)
fr.fit <- coherentfdm(fr.short)
summary(fr.fit)
#> *** PRODUCT MODEL ***
#> Functional demographic model
#> 
#> Call: fdm(data = rate.product, series = "product", order = order1) 
#> 
#> Region: FRATNP
#> Data type: mortality
#> Years in fit: 1950 - 2006
#> Ages in fit: 0 - 110 
#> 
#> Order: 6
#> Percentage variation due to basis functions: 95.1% 2.4% 1.3% 0.4% 0.3% 0.2%
#> 
#> Averages across ages:
#>       ME      MSE      MPE     MAPE 
#>  0.00012  0.00043 -0.00341  0.02145 
#> 
#> Averages across years:
#>       IE      ISE      IPE     IAPE 
#>  0.01507  0.04514 -0.16071  2.07038 
#> 
#> 
#> *** RATIO MODELS ***
#> 
#> FEMALE
#> Functional demographic model
#> 
#> Call: fdm(data = rate.ratio[[j]], series = names(data$rate)[j], order = order2) 
#> 
#> Region: FRATNP
#> Data type: mortality
#> Years in fit: 1950 - 2006
#> Ages in fit: 0 - 110 
#> 
#> Order: 6
#> Percentage variation due to basis functions: 79.2% 6.9% 4.9% 2.6% 1.7% 1.1%
#> 
#> Averages across ages:
#>      ME     MSE     MPE    MAPE 
#> 0.00000 0.00015 0.00182 0.07358 
#> 
#> Averages across years:
#>       IE      ISE      IPE     IAPE 
#> -0.00041  0.01629  0.25759  8.19290 
#> 
#> 
#> MALE
#> Functional demographic model
#> 
#> Call: fdm(data = rate.ratio[[j]], series = names(data$rate)[j], order = order2) 
#> 
#> Region: FRATNP
#> Data type: mortality
#> Years in fit: 1950 - 2006
#> Ages in fit: 0 - 110 
#> 
#> Order: 6
#> Percentage variation due to basis functions: 79.2% 6.9% 4.9% 2.6% 1.7% 1.1%
#> 
#> Averages across ages:
#>      ME     MSE     MPE    MAPE 
#> 0.00000 0.00015 0.00182 0.07358 
#> 
#> Averages across years:
#>      IE     ISE     IPE    IAPE 
#> 0.00041 0.01629 0.25759 8.19290 
#> 
plot(fr.fit$product, components = 3)
```
