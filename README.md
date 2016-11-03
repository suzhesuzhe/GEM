pirate
------

Aimed at generating optimal treatment effect modifiers for RCT, pre-treatment scalar patient characteristics are linearly combined as a Generated Effect Modifier (GEM). This package gives functions for

-   constructing the GEM

-   simulating three types of data set

-   calculating the population average benefit of a GEM model for a data set

-   calculating the effect size of a single moderator

-   calculating the permutation p-value of a GEM

See more detail in: E Petkova, T Tarpey, Z Su, and RT Ogden. Generated effect modifiers (GEMs) in randomized clinical trials. Biostatistics, (First published online: July 27, 2016). doi: 10.1093/biostatistics/kxw035.

Installation
------------

You can install:

-   the latest released version from CRAN with

    ``` r
    install.packages("pirate")
    ```

-   the latest development version from github with

    ``` r
    if (packageVersion("devtools") < 1.6) {
      install.packages("devtools")
    }
    devtools::install_github("suzhesuzhe/GEM")
    ```

Usage
-----

`library(pirate)` will load the core packages.

-   For fitting the GEM model with your own data set, please use `gem_fit` function and make sure the data frame is organized with first column as the treatment index, second column as the outcome, and the remaining columns as the covariates. One of the three methods could be choose and the default method is F-statistics.

    ``` r
    model <- gem_fit(dat = dat, method = "nu")
    ```

    You could get the permutation p-value of the GEM model by:

    ``` r
    permute_pvalue(dat,method = "nu")
    ```

-   For simulating data set, please refer to the 'data\_generators' help page to get detailed information about each data generator.

    ``` r
    #constructing the covariance matrix
    co <- matrix(0.2, 30, 30)
    diag(co) <- 1
    #simulate GEM type data set
    dataEx <- data_generator1(d = 0.3, R2 = 0.5, v2 = 1, n = 3000, 
                            co = co, beta1 = rep(1,30), inter = c(0,0))

    #simulate unconstrained data set with obeservation under both treatment and without error
    bigData <- data_generator3(n = 10000, co = co, bet = dataEx[[2]], inter = c(0,0))
    ```

-   For calculating the population average benefit under a GEM fit, the `gemObject` should be the second element of the ouput from the `gem_fit` function or a list with the same structure. please refer to the 'gem\_test' help page to get detailed information about two test functions.

    ``` r
    #dat has outcome under only one treatment
    gem_test_sample(dat, model[[2]])
    #dat with outcome under both treatment 
    gem_test_simsample(bigData[[1]], bigData[[2]], bigData[[3]], model[[2]])
    ```

Contact
-------

For questions and other discussion, please email <zhe.su@nyumc.org>.
