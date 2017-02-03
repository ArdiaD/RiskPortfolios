---
title: 'RiskPortfolios: Computation of risk-based portfolios in R'
date: "3 February 2017"
tags:
- risk
- portfolio
- optimization
- mean-variance
- minimum variance
- inverse-volatility
- equal-risk-contribution
- maximum diversification
- risk-efficient
authors:
- affiliation: Institute of Financial Analysis - University of Neuchâtel
  name: David Ardia
  orcid: 0000-0003-2823-782X
- affiliation: Solvay Business School - Vrije Universiteit Brussel
  name: Kris Boudt
  orcid: null
- affiliation: PSP Investments
  name: Jean-Philippe Gagnon-Fleury
  orcid: null
---

[![status](http://joss.theoj.org/papers/b52ded01411ff8f9f007b84a27e4d6d9/status.svg)](http://joss.theoj.org/papers/b52ded01411ff8f9f007b84a27e4d6d9)
[![CRAN](http://www.r-pkg.org/badges/version/RiskPortfolios)](https://cran.r-project.org/package=RiskPortfolios) 
[![Downloads](http://cranlogs.r-pkg.org/badges/RiskPortfolios?color=brightgreen)](http://www.r-pkg.org/pkg/RiskPortfolios)
[![Downloads](http://cranlogs.r-pkg.org/badges/grand-total/RiskPortfolios?color=brightgreen)](http://www.r-pkg.org/pkg/RiskPortfolios)

# RiskPortfolios
`RiskPortfolios` ([Ardia et al. (2017)](http://dx.doi.org/10.21105/joss.00171)) is an R package for constructing risk-based portfolios dedicated to portfolio managers 
and quantitative analysts. It provides a set of functionalities to build mean-variance, minimum variance, inverse-volatility weighted, 
equal-risk-contribution, maximum diversification, and risk-efficient portfolios. As risk-based portfolios are
mainly based on covariances, the package also provides a large set of covariance matrix estimators. See [Ardia et al. (2017)](http://dx.doi.org/10.21105/joss.00171) for details. A Monte Carlo study relying on `RiskPortfolios` is presented in [Ardia et al. (2016)](http://dx.doi.org/10.2139/ssrn.2650644).

The latest stable version of `RiskPortfolios` is available at 'https://cran.r-project.org/package=RiskPortfolios'.

The latest development version of `RiskPortfolios` is available at 'https://github.com/ArdiaD/RiskPortfolios'.

# Installation
To install the latest stable version of `RiskPortfolios`, run the following commands in R:

    R> install.packages("RiskPortfolios")

To install the development version of `RiskPortfolios`, run the following commands in R:

    R> install.packages("devtools")

    R> devtools::install_github("ArdiaD/RiskPortfolios")

Then check the help of the various files and run the examples:

    R> library("RiskPortfolios")

    R> ?RiskPortfolios
    
    
Please cite `RiskPortfolios` in publications:

    R> citation("RiskPortfolios")
    
    
# References
Ardia, D., Bolliger, G., Boudt, K., Gagnon-Fleury, J.-P. (2016).  
_The impact of covariance misspecification in risk-based portfolios_.  
Working paper. 
http://dx.doi.org/10.2139/ssrn.2650644

Ardia, D., Boudt, K., Gagnon-Fleury, J.-P. (2017).  
RiskPortfolios: Computation of risk-based portfolios in R.  
_Journal of Open Source Software_ **10**(2).
http://dx.doi.org/10.21105/joss.00171