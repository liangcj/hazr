# Nonparametric hazard regression
**To install**, make sure you have the `devtools` library, then run: `devtools::install_github("liangcj/hazr")`

This is a barebones R package for nonparametric **haz**ard **r**egression (`hazr`). It is written entirely in C++, and given an `R` interface using the [`Rcpp`](https://github.com/RcppCore/Rcpp) package.

The main function in this package, `hazr::hazr`, provides a nonparametric conditional hazard estimate at a single time and univariate marker value. I created this package because another package I am developing requires estimates of the conditional hazard surface at many points, and I found existing solutions inadequate (the `muhaz` package should be more than enough for most other applications though - more on that later).

The details for the conditional hazard estimator are in [McKeague and Utikal (1990)](http://projecteuclid.org/euclid.aos/1176347745), but it basically involves calculating a conditional Nelson-Aalen estimator and kernel smoothing the result. Here is a quick graphical depiction of how the estimator works:

![coha](/pictures/coha.png)

The top is a plot of univariate survival data (taken from Mayo PBC dataset), with time on the x-axis and marker on the y-axis. To estimate the hazard function at marker=7, select the strip of everyone with marker values close to 7, then calculate the Nelson-Aalen cumuative hazard estimator (black line). To get the hazard function (blue line), simply do a kernel smooth of the cumulative hazard.

The marginal hazard can also be calculated by just setting the marker bandwidth to be very wide. For a really nice package that does this with a much more robust feature set (e.g. automated bandwidth selection, local/global bandwidth choice, boundary correction, kernel choice), see [`muhaz`](http://cran.r-project.org/web/packages/muhaz/).

`hazr` mostly replicates some of the core functionality in `muhaz`. The main differences are that `hazr` 1) is written in C++; 2) is slightly more flexible in that it can estimate the hazard at just a single time rather than an entire vector of times, saving computation (this can also be achieved by modifying the `muhaz` FORTRAN code, but I am much more comfortable with C++ so this was more efficient for me); and 3) has a logo based on a dated cat meme.

## Example
Using the Mayo PBC data from the `survival` package, estimate the hazard rate at 1000 days for 50 year olds:

    library(dplyr)
    library(survival)
    library(hazr)
    
    # Slightly edit pbc data to prepare for feeding into hazr::hazr(), save as `exdat`
    exdat <- pbc %>%
      slice(1:312) %>%
      select(time, status, age) %>%
      mutate(status = (status==2)*1) %>%
      arrange(time) %>%
      as.matrix
    
    # t: evaluation time; m: marker value to evaluate at
    # bm: marker bandwidth; bt: time bandwidth
    # dat: data matrix with columns in following order: time, censoring indicator, marker
    hazr(t = 1000, m = 50, bm = 5, bt = 730, dat = exdat)

What if we want to estimate the entire hazard curve for 50 year olds?

    library(dplyr)
    library(survival)
    library(hazr)
    library(foreach)
    
    # Slightly edit pbc data to prepare for feeding into hazr::hazr(), save as `exdat`
    exdat <- pbc %>%
      slice(1:312) %>%
      select(time, status, age) %>%
      mutate(status = (status==2)*1) %>%
      arrange(time) %>%
      as.matrix
    
    hazards <- foreach(i = 1:4500) %do% hazr(t = i, m = 50, bm = 5, bt = 730, dat = exdat)[1] %>%
      unlist
    plot(1:4500, hazards, type = "n")
    lines(1:4500, hazards)

## More about this package
This package was built for SPEED, rather than robustness. So please be careful if using it for your own applications. Any suggestions are always welcome.

For an **interactive JavaScript version**, go to my [interactive hazards repository](https://github.com/liangcj/interactivehazard) where there is also a link to a live demo.

## References
- [McKeague and Utikal (1990, Annals of Statistics)](http://projecteuclid.org/euclid.aos/1176347745)
- [Wang (2005, Encyclopedia of Biostatistics)](http://anson.ucdavis.edu/~wang/paper/hazardeob4.pdf)
