---
title: 'RiskPortfolios: Computation of Risk-Based Portfolios in R'
bibliography: paper.bib
date: "22 December 2016"
tags:
- risk
- portfolio
- optimization
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

# Summary

`RiskPortfolios` is an R package (@R) for constructing risk-based portfolios. It provides a set of 
functionalities to build mean-variance, minimum variance, inverse-volatility weighted (@LeoteEtAl2012),
equal-risk-contribution (@MaillardEtAl2010), maximum diversification (@Choueifaty2008), and
risk-efficient (@AmencEtAl2011) portfolios. Optimization is achieved with the R packages `quadprog` (@quadprog) and `nloptr` (@nloptr).
Long or gross constraints can be added to the optimizer. 
As risk-based portfolios are mainly based on covariances, the package also provides a large set of covariance matrix estimators. A simulation study relying on the package is described in @ArdiaEtAl2016. The latest version of the package is 
available at 'https://github.com/ArdiaD/RiskPortfolios}'.


# References