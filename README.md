# Nonparametric hazard regression
This is a barebones R package for nonparametric **haz**ard **r**egression (`hazr`). It is written entirely in C++, and given an `R` interface using the [`Rcpp`](https://github.com/RcppCore/Rcpp) package.

The main function in this package, `hazr::coha`, provides a nonparametric conditional hazard estimate at a specified time and univariate marker value. I created this package because another package I am developing requires estimates of the conditional hazard surface at many points, and I found existing solutions to be too slow.

The marginal hazard can also be calculated by just setting the marker bandwidth to be very wide. For a package that does this with a much more robust feature set (e.g. automated bandwidth selection, local/global bandwidth choice, kernel choice), see [`muhaz`](http://cran.r-project.org/web/packages/muhaz/).

## Example
Using the Mayo PBC data from the `survival` package, estimate the hazard rate at 1000 days for 50 year olds:

    library(dplyr)
    library(survival)
    library(hazr)
    
    # Slightly edit pbc data to prepare for feeding into coha(), save as `exdat`
    exdat <- pbc %>%
      select(time, status, age) %>%
      mutate(status = (status==2)*1) %>%
      as.matrix
      
    # t: evaluation time; m: marker value to evaluate at
    # bm: marker bandwidth; bt: time bandwidth
    # dat: data matrix with columns in following order: time, censoring indicator, marker
    coha(t = 1000, m = 50, bm = 5, bt = 365, dat = exdat)


## References
- [McKeague and Utikal (1990, Annals of Statistics)](http://projecteuclid.org/euclid.aos/1176347745)
- [Wang (2005, Encyclopedia of Biostatistics)](http://anson.ucdavis.edu/~wang/paper/hazardeob4.pdf)
