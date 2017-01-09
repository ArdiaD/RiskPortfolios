---
title: 'RiskPortfolios: Computation of Risk-Based Portfolios in R'
date: "09 January 2017"
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

# Installation
To install the package, run the following commands in R:

R> install.packages("devtools")

R> devtools::install_github("ArdiaD/RiskPortfolios", dependencies = TRUE)

Then check the help of the various files, and run the examples:

R> library("RiskPortfolios")

R> ?optimalPortfolio

R> ?covEstimation

R> ?meanEstimation

R> ?semidevEstimation