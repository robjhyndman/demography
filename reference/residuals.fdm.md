# Compute residuals and fitted values from functional demographic model or Lee-Carter model

After fitting a Lee-Carter model or functional demographic model, it is
useful to inspect the residuals or plot the fitted values. These
functions extract the relevant information from the fit object.

## Usage

``` r
# S3 method for class 'fdm'
residuals(object, ...)

# S3 method for class 'fdm'
fitted(object, ...)

# S3 method for class 'lca'
fitted(object, ...)

# S3 method for class 'lca'
residuals(object, ...)
```

## Arguments

- object:

  Output from [`fdm`](fdm.md) or [`lca`](lca.md).

- ...:

  Other arguments.

## Value

`residuals.fdm` and `residuals.lca` produce an object of class “fmres”
containing the residuals from the model. `fitted.fdm` and `fitted.lca`
produce an object of class “fts” containing the fitted values from the
model.

## See also

[`fdm`](fdm.md), [`lca`](lca.md), [`bms`](lca.md)

## Author

Rob J Hyndman.

## Examples

``` r
fit1 <- lca(fr.mort)
plot(residuals(fit1))

plot(fitted(fit1))
```
