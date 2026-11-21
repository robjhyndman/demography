# Summary for functional demographic model or Lee-Carter model

Summarizes a basis function model fitted to age-specific demographic
rate data. It returns various measures of goodness-of-fit.

## Usage

``` r
# S3 method for class 'fdm'
summary(object, ...)

# S3 method for class 'lca'
summary(object, ...)
```

## Arguments

- object:

  Output from [`fdm`](fdm.md) or [`lca`](lca.md).

- ...:

  Other arguments.

## See also

[`fdm`](fdm.md), [`lca`](lca.md), [`bms`](lca.md),
[`compare.demogdata`](compare.demogdata.md)

## Author

Rob J Hyndman

## Examples

``` r
fit1 <- lca(fr.mort)
fit2 <- bms(fr.mort, breakmethod = "bai")
fit3 <- fdm(fr.mort)
summary(fit1)
#> Lee-Carter analysis
#> 
#> Call: lca(data = fr.mort) 
#> 
#> Adjustment method: dt
#> Region: FRATNP
#> Years in fit: 1816 - 2006
#> Ages in fit: 0 - 100 
#> 
#> Percentage variation explained: 96.1%
#> 
#> ERROR MEASURES BASED ON MORTALITY RATES
#> 
#> Averages across ages:
#>       ME      MSE      MPE     MAPE 
#> -0.00051  0.00067  0.03237  0.14321 
#> 
#> Averages across years:
#>       IE      ISE      IPE     IAPE 
#> -0.04567  0.05467  3.24188 14.25819 
#> 
#> 
#> ERROR MEASURES BASED ON LOG MORTALITY RATES
#> 
#> Averages across ages:
#>      ME     MSE     MPE    MAPE 
#> 0.00556 0.05439 0.00639 0.05591 
#> 
#> Averages across years:
#>      IE     ISE     IPE    IAPE 
#> 0.55675 5.42655 0.58512 5.05281 
#> 
summary(fit2)
#> Lee-Carter analysis
#> 
#> Call: bms(data = fr.mort, breakmethod = "bai") 
#> 
#> Adjustment method: dxt
#> Region: FRATNP
#> Years in fit: 1979 - 2006
#> Ages in fit: 0 - 100 
#> 
#> Percentage variation explained: 91.2%
#> 
#> ERROR MEASURES BASED ON MORTALITY RATES
#> 
#> Averages across ages:
#>      ME     MSE     MPE    MAPE 
#> 0.00005 0.00001 0.00517 0.04500 
#> 
#> Averages across years:
#>      IE     ISE     IPE    IAPE 
#> 0.00514 0.00103 0.51670 4.37285 
#> 
#> 
#> ERROR MEASURES BASED ON LOG MORTALITY RATES
#> 
#> Averages across ages:
#>       ME      MSE      MPE     MAPE 
#>  0.00301  0.00433 -0.00046  0.00948 
#> 
#> Averages across years:
#>       IE      ISE      IPE     IAPE 
#>  0.30106  0.40835 -0.04694  0.90695 
#> 
summary(fit3)
#> Functional demographic model
#> 
#> Call: fdm(data = fr.mort) 
#> 
#> Region: FRATNP
#> Data type: mortality
#> Years in fit: 1816 - 2006
#> Ages in fit: 0 - 110 
#> 
#> Order: 6
#> Percentage variation due to basis functions: 93.8% 2.3% 1.6% 0.6% 0.4% 0.3%
#> 
#> Averages across ages:
#>       ME      MSE      MPE     MAPE 
#> -0.09913  0.05686     -Inf      Inf 
#> 
#> Averages across years:
#>       IE      ISE      IPE     IAPE 
#>  0.45441  0.73512 -1.66272  8.08544 
#> 
```
