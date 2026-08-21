## Submission

This is a bug-fix release for RiskPortfolios 2.1.7, currently on CRAN.

It corrects errors that made several estimators and portfolio optimizers return
incorrect numbers without any indication. The principal ones are:

* the equal-risk-contribution portfolio could return weights whose risk
  contributions were far from equal;
* the Bayes-Stein mean used the number of assets in place of the number of
  observations in its shrinkage intensity;
* the maximum-decorrelation portfolio was computed from the covariance matrix
  rather than the correlation matrix under the default constraint;
* the risk aversion parameter had no effect on the unconstrained mean-variance
  portfolio;
* the exponentially weighted covariance matrix was not normalized;
* bounds supplied by the caller were, in some cases, imposed after the
  optimization rather than within it.

Because these corrections change the values that are returned, results obtained
with 2.1.7 will differ from those obtained with 2.1.8. Every change is itemised
in NEWS, and the package now ships numerical regression tests that compare the
returned values against independently computed references rather than only
checking their shape.

There are no reverse dependencies on CRAN, so no other package is affected;
this was verified with `tools::package_dependencies(reverse = TRUE)` over
`Depends`, `Imports`, `LinkingTo` and `Suggests`.

## Test environments

* local: macOS 26.6.1, aarch64-apple-darwin20, R 4.5.2
* win-builder: R-release (R 4.6.1, x86_64-w64-mingw32, Windows Server 2022)
* win-builder: R-devel (submitted; please see the note below)

## R CMD check results

win-builder R-release: Status: OK. There were no ERRORs, WARNINGs or NOTEs.

win-builder R-devel: submitted at the same time as R-release; the result had not
yet been returned when these comments were written.

Locally there were no ERRORs or WARNINGs, and one NOTE:

```
* checking for future file timestamps ... NOTE
  unable to verify current time
```

This NOTE reflects the local machine being unable to reach the time service used
by that check; win-builder reports the same check as OK.

## Downstream dependencies

There are none.
