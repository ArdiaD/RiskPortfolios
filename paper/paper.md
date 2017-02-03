---
title: 'RiskPortfolios: Computation of Risk-Based Portfolios in R'
bibliography: paper.bib
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
- affiliation: 1
  name: David Ardia
  orcid: 0000-0003-2823-782X
- affiliation: 2
  name: Kris Boudt
- affiliation: 3
  name: Jean-Philippe Gagnon-Fleury
affiliations:
- name: Institute of Financial Analysis - University of Neuchâtel
  index: 1
- name: Solvay Business School - Vrije Universiteit Brussel
  index: 2
- name: PSP Investments
  index: 3

---

# Summary

`RiskPortfolios` is an R package (@R) for constructing risk-based portfolios. It provides a set of
functionalities to build mean-variance, minimum variance, inverse-volatility weighted (@LeoteEtAl2012),
equal-risk-contribution (@MaillardEtAl2010), maximum diversification (@Choueifaty2008), and
risk-efficient (@AmencEtAl2011) portfolios. Optimization is achieved with the R packages `quadprog` (@quadprog) and `nloptr` (@nloptr).
Long or gross constraints can be added to the optimizer.
As risk-based portfolios are mainly based on covariances, the package also provides a large set of covariance matrix estimators. A simulation study relying on the package is described in @ArdiaEtAl2016. The latest version of the package is
available at 'https://github.com/ArdiaD/RiskPortfolios'.


# References
