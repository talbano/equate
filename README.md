# equate

The R package *equate* provides methods for linking and equating measurements scales.

Install from CRAN with:
```R
install.packages("equate")
```

Or install the latest development version from Github using the *devtools* package:
```R
install.packages("devtools")
devtools::install_github("talbano/equate")
```

The general process of equating or linking is divided into three main steps:

1. Data prep, using the `freqtab` and `presmoothing` functions.
2. Linking and equating score scales, using the `equate` and `composite` functions.
3. Evaluating results, using `summary` methods, `plot` methods, and `bootstrap` standard errors. 

Some basic instructions on using the package are included below. For a detailed intro, see the main package vignette with `vignette("equatevignette", package = "equate")`. 

## 1. Data Prep

Data prep starts with the `freqtab` function, which converts vectors and matrices of raw total scores into frequency tables of class `freqtab`. These are simply undimensional or multidimensional arrays with `dimnames` set to the full score scale of each variable.

Data can be structured according to the single-group, equivalent-groups, and nonequivalent-groups with anchor test(s) designs.

The `presmoothing` function fits polynomial loglinear models to frequency distributions, with support for automatic creation of polynomial terms, model fit comparisons, and plotting.

## 2. Linking and Equating

The `equate` function does all the work. Give it your data, the `type` of equating function to apply, and the `method` for handling nonequivalent groups if necessary. You'll get descriptive statistics, conversion tables, synthetic score distributions, standard errors, and other goodies.

Equating types include identity, mean, linear, general linear, equipercentile, circle-arc, and composites of these. Methods for nonequivalent groups include synthetic, nominal weights, Tucker, Levine observed score, Levine true score, Braun/Holland, frequency estimation, and chained equating.

## 3. Evaluating Results
A `boostrap` function facilitates estimation of empirical and parametric bootstrap standard errors.

`summary` and `plot` methods are available for equating and bootstrap output.
