Simple simulated data set and example
================
Leo Frach
2023-12-21

## Citation

The following script provides an example on how to run Gsens, a
genetically informed sensitivity analysis, in R. Please cite the
following papers where the concept was first proposed and then
implemented:

1.  Pingault, J.-B., O’Reilly, P. F., Schoeler, T., Ploubidis, G. B.,
    Rijsdijk, F., & Dudbridge, F. (2018). Using genetic data to
    strengthen causal inference in observational research. Nature
    Reviews Genetics, 19(9), 566–580.
    <https://doi.org/10.1038/s41576-018-0020-3>

2.  Pingault, J. B., Rijsdijk, F., Schoeler, T., Choi, S. W., Selzam,
    S., Krapohl, E., … & Dudbridge, F. (2021). Genetic sensitivity
    analysis: Adjusting for genetic confounding in epidemiological
    associations. PLoS genetics, 17(6), e1009590.
    <https://doi.org/10.1371/journal.pgen.1009590>

3.  Frach, L., Rijsdijk, F., Dudbridge, F., & Pingault, J. B. (2022,
    November). Adjusting for Genetic Confounding Using Polygenic Scores
    Within Structural Equation Models. In BEHAVIOR GENETICS (Vol. 52,
    No. 6, pp. 359-359).

## Usage

We strongly recommend using the updated `gsensY()` function, since it
has been optimised, e.g., by including raw data as input.<br> As for all
lavaan models, we recommend *against* using a correlation
matrix as input but rather a covariance matrix if only summary data is
available.<br> Lastly, we recommend standardising the polygenic score
before regressing out potential batch effects (e.g., genotyping array,
genetic PCs).<br> All phenotypic variables should **not** be
standardised. However, if you want to get standardised estimates, the
additional `lavaan` argument `std.all = TRUE` can be passed
on to `gsensY()`, or output can be customized using the `parameterEstimates()` function with the argument `standardized = TRUE` on the output of the `gsensY()` function.

## Help

Run `?gsensY()` in your command line to read the function and argument
description.

## Simulate population model with two exposures, one outcome and one genetic factor

### Load packages for simulation

``` r
library(lavaan)
library(stringr)
library(simstandard)
library(dplyr)
```

### Simulation parameters for the population model

``` r
n = 1e4 # sample size
b1 = 0.10 # effect of X1 on outcome
b2 = 0.05 # effect of X2 on outcome
a1 = 0.10 # effect of PGS_outcome on X1
a2 = 0.08 # effect of PGS_outcome on X2
c = 0.6 # effect of PGS_outcome on outcome
pme = 0.80 # measurement error of G   

h <- a1*b1 + a2*b2 + c # h^2 = heritability 
```

### Create the population model

``` r
# Simulate an underlying model, where the polygenic score (G) is a noisy measure of the true genetic factor (GF), which comes with measurement error (pme)

population.modelGF = str_glue('  # use library stringr to pass the estimates to the model
             #paths and loading
             Y ~ {b1}*X1 + {b2}*X2 + {c}*GF  
             X1 ~ {a1}*GF 
             X2 ~ {a2}*GF
             GF =~ sqrt(1-{pme})*G  
          ')

# Simulate data accordingly
set.seed(123)
myData <- sim_standardized(population.modelGF, n = n,
                           latent = TRUE,
                           errors = TRUE)
```

### Checks

``` r
dim(myData)
```

    ## [1] 10000     9

``` r
head(myData)
```

    ## # A tibble: 6 × 9
    ##        Y     X1     X2      G     GF    e_Y   e_X1   e_X2     e_G
    ##    <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>   <dbl>
    ## 1 -0.304 -0.216  1.56   0.121  0.129 -0.438 -0.229  1.55   0.0631
    ## 2  1.05   0.414 -1.30  -0.814 -0.446  1.34   0.459 -1.26  -0.614 
    ## 3  0.671  0.302  0.355 -0.150 -0.556  0.957  0.358  0.399  0.0990
    ## 4  1.06   0.448 -2.00   0.416 -0.473  1.40   0.495 -1.96   0.627 
    ## 5 -1.29  -0.279 -1.07  -0.931 -0.625 -0.834 -0.217 -1.02  -0.652 
    ## 6 -0.457  0.959  0.253 -0.457  1.25  -1.32   0.834  0.153 -1.02

``` r
# Create covariance matrix
(cov1 = cov(myData)) #check all variances and covariances
```

    ##                Y           X1           X2             G            GF
    ## Y     1.00280986  0.143847452  0.097657243  0.2610850728  0.6081840000
    ## X1    0.14384745  0.990252901  0.004381535  0.0349070363  0.0844595533
    ## X2    0.09765724  0.004381535  0.994164814  0.0436810085  0.0717134770
    ## G     0.26108507  0.034907036  0.043681008  0.9832349333  0.4249479805
    ## GF    0.60818400  0.084459553  0.071713477  0.4249479805  0.9951226518
    ## e_Y   0.61863185 -0.006072647  0.004482762  0.0004415305 -0.0009212203
    ## e_X1  0.08302905  0.981806946 -0.002789812 -0.0075877618 -0.0150527119
    ## e_X2  0.04900252 -0.002375229  0.988427736  0.0096851700 -0.0078963351
    ## e_G  -0.01090308 -0.002864424  0.011609767  0.7931924191 -0.0200843986
    ##                e_Y          e_X1         e_X2           e_G
    ## Y     0.6186318545  0.0830290519  0.049002523 -0.0109030805
    ## X1   -0.0060726469  0.9818069457 -0.002375229 -0.0028644242
    ## X2    0.0044827621 -0.0027898125  0.988427736  0.0116097666
    ## G     0.0004415305 -0.0075877618  0.009685170  0.7931924191
    ## GF   -0.0009212203 -0.0150527119 -0.007896335 -0.0200843986
    ## e_Y   0.6195677132 -0.0059805249  0.004556460  0.0008535127
    ## e_X1 -0.0059805249  0.9833122169 -0.001585596 -0.0008559844
    ## e_X2  0.0045564598 -0.0015855955  0.989059443  0.0132165184
    ## e_G   0.0008535127 -0.0008559844  0.013216518  0.8021744352

### Create/load gsensY() function from [gsens.R](R/gsens.R)

``` r
# Install and load package
devtools::install_github("LeonardFrach/Gsens")
library(Gsens)
```

### Run Gsens

``` r
# Using raw data
gsensY(myData, h2 = h^2, exposures = c("X1", "X2"), pgs = "G", outcome = "Y") # this should correspond to the population model
```

    ## Using raw data as input.

    ##                            est    se     z   pvalue ci.lower ci.upper
    ## Adjusted Bx1y            0.095 0.015 6.383 1.74e-10    0.066    0.124
    ## Adjusted Bx2y            0.036 0.015 2.349 1.88e-02    0.006    0.065
    ## Mediation m1             0.008 0.001 5.256 1.47e-07    0.005    0.011
    ## Mediation m2             0.004 0.001 3.432 5.99e-04    0.002    0.006
    ## Total mediation          0.011 0.002 6.197 5.75e-10    0.008    0.015
    ## Genetic confounding Bx1y 0.050 0.014 3.602 3.16e-04    0.023    0.077
    ## Genetic confounding Bx2y 0.063 0.014 4.417 1.00e-05    0.035    0.091
    ## Genetic overlap x1y      0.050 0.014 3.575 3.50e-04    0.023    0.078
    ## Genetic overlap x2y      0.063 0.014 4.428 9.51e-06    0.035    0.091

``` r
# Using covariance matrix + sample size
gsensY(cov1, sample.nobs = n, h2 = h^2, exposures = c("X1", "X2"), pgs = "G", outcome = "Y")
```

    ## Using covariance matrix as input.

    ##                            est    se     z   pvalue ci.lower ci.upper
    ## Adjusted Bx1y            0.095 0.015 6.383 1.74e-10    0.066    0.124
    ## Adjusted Bx2y            0.036 0.015 2.349 1.88e-02    0.006    0.065
    ## Mediation m1             0.008 0.001 5.256 1.47e-07    0.005    0.011
    ## Mediation m2             0.004 0.001 3.432 5.99e-04    0.002    0.006
    ## Total mediation          0.011 0.002 6.197 5.75e-10    0.008    0.015
    ## Genetic confounding Bx1y 0.050 0.014 3.602 3.16e-04    0.023    0.077
    ## Genetic confounding Bx2y 0.063 0.014 4.417 1.00e-05    0.035    0.091
    ## Genetic overlap x1y      0.050 0.014 3.575 3.50e-04    0.023    0.078
    ## Genetic overlap x2y      0.063 0.014 4.428 9.51e-06    0.035    0.091
