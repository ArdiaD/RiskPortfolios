---
title: 'RiskPortfolios: Computation of Risk-Based Portfolios in R'
date: "10 January 2017"
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

# Introduction
`RiskPortfolios` is an R package for constructing risk-based portfolios dedicated to portfolio managers 
and quantitative analysts. It provides a set of functionalities to build mean-variance, minimum variance, inverse-volatility weighted, 
equal-risk-contribution, maximum diversification, and risk-efficient portfolios. As risk-based portfolios are
mainly based on covariances, the package also provides a large set of covariance matrix estimators.


# Installation
To install the package, run the following commands in R:

R> install.packages("devtools")

R> devtools::install_github("ArdiaD/RiskPortfolios", dependencies = TRUE)

Then check the help of the various files, and run the examples:

R> library("RiskPortfolios")

R> ?RiskPortfolios