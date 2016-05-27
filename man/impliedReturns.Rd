\name{impliedReturns}
\alias{impliedReturns}
\title{
Implied returns estimation
}
\description{
This function computes the implied returns.
}
\usage{
impliedReturns(rets = NULL, mu = NULL, Sigma = NULL, semiDev = NULL, w = NULL, control = list())
}
\arguments{
  \item{rets}{a \eqn{(T \times N)}{(T x N)} matrix of returns.}
  \item{mu}{a \eqn{(N \times 1)}{(N x 1)} vector of mean (expected returns).}
  \item{Sigma}{a \eqn{(N \times N)}{(N x N)} matrix of covariances.}
  \item{semiDev}{a \eqn{(N \times 1)}{(N x 1)} vector of semideviation.}
  \item{w}{a \eqn{(N \times 1)}{(N x 1)} vector of portfolio weights.}
  \item{control}{control parameters (see *Details*).}
}
\details{
  The argument \code{control} is a list that can supply any of the following components:

  \describe{%
    \item{\code{type}}{method used to compute the implied returns, 
    among \code{'regression'}, \code{robust}, \code{constraint} and \code{'rl'} where:
    
   \code{'regression'} is use to compute the implied expected returns. The implied expected returns 
   are solved by fitting the following regression:
   
   \deqn{\tilde{\mu} = a + b[\hat{\Sigma}\hat{w}]}{\tilde{\mu} = a + b[\hat{\Sigma}\hat{w}]}
   
  The option \code{'reg'} can be used with \code{'regression'} if you want a robust estimation or a 
  constrained estimation. See \code{'reg'}.

   \code{'bl'} is used to compute the Black-Litterman implied expected returns. The Black-Litterman implied 
   expected return is estimated by fitting the following regression:
   
   \deqn{\mu = l \iota + \gamma \Sigma w.}
   }

   \item{\code{reg}}{is used when type = \code{'regression'}. \code{reg} can supply either \code{'robust'} or \code{'constraint'}.
    
   \code{'robust'} is used to compute the robust estimation of the implied expected returns.

   \code{'constraint'} is used to compute the implied expected returns with the restriction that the coefficient is larger than zero. 
   
   The constraint implied expected return estimate \eqn{\tilde{\mu}} is given by the following equation:
   
   \deqn{\tilde{\mu} = a + b(\hat{\Sigma}\hat{w})}{\tilde{\mu} = a + b(\hat{\Sigma}\hat{w})}
   %%
   } 
   \item{\code{gamma}}{is used for the computation of the Black-Litterman implied expected returns. Default: \code{gamma = 0.89}.}
  }
%%
First, if \eqn{w} is a mean-variance efficient portfolio with respect to a universe of assets with
known expected returns vector \eqn{\mu} and covariance matrix \eqn{\Sigma}, then there exists a linear relation between \eqn{\mu} 
and \eqn{\Sigma w}. Second, covariances can be much more accurately estimated than expected returns. 
Hence,if a mean-variance efficient weight vector is known and a covariance matrix has been 
accurately estimated, the above-mentioned linear relation can be exploited to create
implied expected returns which have a higher accuracy than sample means. The idea of reverse engineering 
for expected returns was first advocated by Best and Grauer (1985), and almost always applied with 
the market capitalization portfolio, which is well known to be mean-variance efficient under the 
Capital Asset Pricing Model (CAPM) assumptions. Third, in realistic markets, stock prices may deviate 
from their fundamental values and market capitalization portfolios tend to invest disproportionately 
in stocks that are overpriced and under-invest in stocks that are under-priced (Arnott et al., 2005; Treynor, 2005). 
Hence, we advocate that the market capitalization weighted portfolio is only one possible proxy for a 
mean-variance efficient portfolio, and that other proxies may lead to more accurate expected returns.
%%
Our main hypothesis is that there exist other (expected) return-insensitive proxies for mean-variance 
efficient portfolio weights that lead to more accurate expected return predictions. In particular, 
we will test this hypothesis for the smoothed market capitalization portfolio (Chen et al., 2007),
the equally-weighted portfolio (DeMiguel et al., 2009), the equal-risk-contribution 
portfolio (Maillard et al., 2010), the maximum diversification portfolio (Choueifaty and Coignard, 2008) and 
the risk-efficient portfolio (Amenc et al., 2011). See e.g. Lee (2011) for a recent overview.
%%
We consider a market with \eqn{N} risky securities and denote a generic portfolio in this 
market by the \eqn{(N \times 1)}{(N x 1)} vector \eqn{w}. 
%%
The expected arithmetic returns (in excess of the risk-free rate) at the desired holding horizon are denoted by 
the \eqn{(N \times 1)}{(N x 1)} vector \eqn{\mu} and the corresponding \eqn{(N \times N)}{(N x N)} covariance 
matrix of arithmetic returns is denoted by \eqn{\Sigma}. We denote by \eqn{\iota} the \eqn{(N \times 1)}{(N x 1)} 
vector of ones and by \eqn{0} the \eqn{(N \times 1)}{(N x 1)} vector of zeros.
%%
Our analysis builds on the assumption of mean--variance preferences. Let \eqn{0 < \gamma < \infty} be the 
risk aversion parameter. The mean-variance optimization problem is:
%%
\deqn{w^* = {argmax}\left\{\mu' w - \frac{1}{2} \gamma w' \Sigma w \right \}}{w* = argmax\{\mu' w - 1/2\gamma w'\Sigma w\}}
%%
where the optimization is done under the full-investment constraint. The corresponding Lagrangian is:
%%
\deqn{ L(w, l) = w'\mu -  \frac{\gamma}{2} \, w' \Sigma w - l (w' \iota - 1)\,}{L(w,l) = w' \mu - \gamma / 2 w'\Sigma w - l(w'\iota -1 )}
%%
The corresponding first order conditions are:
%%
\deqn{\mu - \gamma \Sigma w - l\iota = 0}
\deqn{w'\iota = 1}
%%
We note the linear relationship between \eqn{\mu} and \eqn{\Sigma w}:
%%
\deqn{\mu = l \iota + \gamma \Sigma w}
%%
}
\value{
A \eqn{(N \times 1)}{(N x 1)} vector of implied expected returns.
}
\references{
  Ardia, D., Boudt, K. (2013). 
  Implied expected returns and the choice of a mean-variance efficient portfolio proxy.
  \url{http://papers.ssrn.com/sol3/papers.cfm?abstract_id=2215042}.

  Best, M. J., Grauer, R. R., (1985). 
  Capital asset pricing compatible with observed market value weights. 
  \emph{Journal of Finance} \bold{40}(5), pp.85--103.
  
  Black, F., Litterman, R., (September-October 1992). 
  Global portfolio optimization. 
  \emph{Financial Analyst Journal} \bold{48}(5), pp.28--43.
}
\author{
David Ardia <david.ardia@unifr.ch> and Jean-Philippe Gagnon Fleury.
}
\examples{
# For the examples, we simply generate a 100 x 25 random matrix.
set.seed(3214)
T = 100
N = 25
rets = matrix(rnorm(T*N), nrow = T, ncol = N)

mu    = meanEstimation(rets)
Sigma = covEstimation(rets)
w     = rep(1, N) / N

#Computes the implied expected returns by Black-Litterman with gamma = 0.89.
impliedReturns(mu = mu, Sigma = Sigma, w = w, control = list(type = "bl"))

#Computes the implied expected returns by Black-Litterman with gamma = 1.
impliedReturns(mu = mu, Sigma = Sigma, w = w, control = list(type = "bl", gamma = 1))

#Computes the impled expected returns.
impliedReturns(mu = mu, Sigma = Sigma, w = w, control = list(type = "regression"))

#Compute the robust implied expected returns.
#impliedReturns(mu = mu, Sigma = Sigma, w = w, control = list(type = "regression", reg = "robust"))

#Compute the constraint implied expected returns.
#impliedReturns(mu = mu, Sigma = Sigma, w = w, control = list(type = "regression", reg = "constraint"))
}
\keyword{htest}