## Internal helpers shared by the estimation and optimization functions.
## Not exported, not documented.

.checkLambda <- function(lambda, allowOne = FALSE) {
  ## Decay parameter of the exponential weighting schemes. lambda = 0 gives a
  ## degenerate (all-zero) weight vector everywhere. lambda = 1 is the no-decay
  ## limit: it is a legitimate request for the mean and the semideviation, where
  ## it yields equal weights, but it divides by zero in the finite-sample
  ## normalisation (1 - lambda^T) of the ewma covariance.
  ok <- is.numeric(lambda) && length(lambda) == 1L && is.finite(lambda) &&
    lambda > 0 && (lambda < 1 || (allowOne && lambda == 1))
  if (!ok) {
    stop("'lambda' must be a single number in (0, 1", if (allowOne) "]" else ")")
  }
  invisible(lambda)
}

.checkReturns <- function(rets) {
  ## Shared validation of the (T x N) return matrix. Non-finite entries used to
  ## propagate silently into an all-NA covariance matrix and on into the
  ## optimizers.
  if (!is.matrix(rets)) {
    stop("rets must be a (T x N) matrix")
  }
  if (!is.numeric(rets)) {
    stop("rets must be a numeric matrix")
  }
  if (nrow(rets) < 2L || ncol(rets) < 1L) {
    stop("rets must have at least two rows and one column")
  }
  if (any(!is.finite(rets))) {
    stop("rets must not contain missing or infinite values")
  }
  invisible(rets)
}

.checkPSD <- function(Sigma) {
  ## A covariance matrix must be positive semidefinite. Without this, an
  ## indefinite matrix is accepted and, for instance, "minimum variance" is
  ## reported for a problem whose objective is unbounded below.
  ## Returns TRUE when Sigma is positive definite, FALSE when it is singular
  ## but admissible; the Cholesky attempt makes the usual case cheap.
  if (!inherits(try(chol(Sigma), silent = TRUE), "try-error")) {
    return(TRUE)
  }
  ev <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values
  scale <- max(abs(ev))
  tol <- 8 * ncol(Sigma) * max(scale, 1) * .Machine$double.eps
  if (min(ev) < -tol) {
    stop("Sigma is not positive semidefinite (smallest eigenvalue ",
         format(min(ev), digits = 3), "); it cannot be a covariance matrix")
  }
  return(FALSE)
}

.requirePD <- function(Sigma, what) {
  ## Objectives of the form (w'a) / sqrt(w'Sigma w) are unbounded over an
  ## unbounded feasible set as soon as Sigma is singular: any direction z with
  ## Sigma z = 0 drives the ratio to infinity. Refuse rather than return the
  ## enormous weights the optimizer wanders off to.
  if (inherits(try(chol(Sigma), silent = TRUE), "try-error")) {
    stop(what, " is unbounded when Sigma is singular; supply 'LB'/'UB' or the ",
         "'gross' constraint to bound the problem, or use a non-singular ",
         "covariance estimate")
  }
  invisible(TRUE)
}

.solveOrStop <- function(Sigma) {
  ## solve() reports a raw LAPACK message on a singular sample covariance
  ## matrix, which does not tell the caller what to do about it
  out <- try(solve(Sigma), silent = TRUE)
  if (inherits(out, "try-error")) {
    stop("the sample covariance matrix is singular or ill-conditioned and ",
         "cannot be inverted; the assets are collinear or the sample is too short")
  }
  return(out)
}

.shrinkage <- function(num, gamma, t) {
  ## Shrinkage intensity of the Ledoit-Wolf estimators, truncated to [0, 1].
  ## gamma is the squared Frobenius distance between the sample estimate and
  ## the prior; it is zero when the two coincide (a single asset), in which
  ## case the intensity is irrelevant and is set to zero rather than left NaN.
  if (!is.finite(gamma) || gamma <= 0) {
    return(0)
  }
  return(max(0, min(1, (num/gamma)/t)))
}

.projectBox <- function(w, LB = NULL, UB = NULL, tol = 1e-12, maxit = 200L) {
  ## Euclidean projection of w onto {LB <= w <= UB, w'1 = 1}. The solution of
  ##   min ||v - w||^2  s.t.  LB <= v <= UB, v'1 = 1
  ## has the form v(theta) = min(max(w - theta, LB), UB) with theta the
  ## multiplier of the equality constraint, and sum(v(theta)) is continuous and
  ## non-increasing in theta, so theta is found by bracketing and bisection.
  ## (Clipping into the box first and then spreading the residual equally over
  ## the free coordinates does restore feasibility, but is NOT this projection.)
  n <- length(w)
  lb <- if (is.null(LB)) rep(-Inf, n) else LB
  ub <- if (is.null(UB)) rep(Inf, n) else UB
  clip <- function(theta) pmin(pmax(w - theta, lb), ub)

  ## an infeasible box has no projection; return the closest point in the box
  ## and leave it to the caller to report the problem
  sumLB <- sum(lb)
  sumUB <- sum(ub)
  if (is.finite(sumLB) && sumLB > 1 + tol) {
    return(lb)
  }
  if (is.finite(sumUB) && sumUB < 1 - tol) {
    return(ub)
  }
  ## degenerate boxes for which theta is not finite
  if (is.finite(sumLB) && abs(sumLB - 1) <= tol) {
    return(lb)
  }
  if (is.finite(sumUB) && abs(sumUB - 1) <= tol) {
    return(ub)
  }

  lo <- -1
  while (sum(clip(lo)) < 1 && is.finite(lo)) {
    lo <- 2 * lo
  }
  hi <- 1
  while (sum(clip(hi)) > 1 && is.finite(hi)) {
    hi <- 2 * hi
  }
  for (it in seq_len(maxit)) {
    mid <- 0.5 * (lo + hi)
    s <- sum(clip(mid))
    if (abs(s - 1) < tol || hi - lo < tol) {
      break
    }
    if (s > 1) {
      lo <- mid
    } else {
      hi <- mid
    }
  }
  return(clip(0.5 * (lo + hi)))
}

.startingValues <- function(ctr, Sigma) {
  ## Starting values handed to slsqp. The objectives of the risk-based
  ## portfolios are not convex, so a single starting value is not enough: the
  ## optimizer is run from each of them and the best solution is kept.
  n <- dim(Sigma)[1]
  sig <- sqrt(diag(Sigma))
  invvol <- if (all(sig > 0)) sig^(-1)/sum(sig^(-1)) else rep(1, n)/n
  pts <- list(ctr$w0, rep(1, n)/n, invvol)
  if (!is.null(ctr$LB) && !is.null(ctr$UB)) {
    pts <- c(pts, list(0.5 * (ctr$LB + ctr$UB)))
  }
  pts <- lapply(pts, .projectBox, LB = ctr$LB, UB = ctr$UB)
  return(unique(pts))
}

.slsqpBest <- function(starts, fn, gr = NULL, hin = NULL, ctr) {
  ## Run slsqp from every starting value and keep the best feasible solution.
  ## The convergence status is inspected, which slsqp callers must do: the
  ## optimizer returns $par whether or not it succeeded.
  best <- NULL
  for (x0 in starts) {
    out <- try(nloptr::slsqp(x0 = x0, fn = fn, gr = gr, hin = hin,
                             heq = .eqConstraint, lower = ctr$LB, upper = ctr$UB,
                             nl.info = FALSE, control = ctr$ctr.slsqp,
                             deprecatedBehavior = FALSE), silent = TRUE)
    if (inherits(out, "try-error") || !is.finite(out$value)) {
      next
    }
    if (is.null(best) || out$value < best$value) {
      best <- out
    }
  }
  if (is.null(best)) {
    stop("the optimizer failed for every starting value")
  }
  ## nloptr codes: 1 success, 2 stopval, 3 ftol and 4 xtol reached are genuine
  ## convergence; anything negative is a failure; 5 maxeval and 6 maxtime mean
  ## the budget ran out before a stopping criterion was met. Only 1-4 may be
  ## accepted silently. This is informative only because the default control
  ## carries an attainable ftol_rel: with xtol_rel = 1e-18 alone every
  ## gross-constrained solve ended on code 5 whether or not it had converged.
  if (!(best$convergence %in% 1:4)) {
    warning("slsqp did not converge (code ", best$convergence, "): ", best$message,
            if (best$convergence %in% c(5L, 6L))
              " -- increase control$ctr.slsqp$maxeval")
  }
  return(best$par)
}

.finalizeWeights <- function(w, ctr) {
  ## Enforce the bounds and the summability constraint on the solution
  ## returned by the optimizer.
  w <- as.numeric(w)
  if (any(!is.finite(w))) {
    stop("the optimizer returned non-finite weights")
  }

  if (is.null(ctr$LB) && is.null(ctr$UB)) {
    s <- sum(w)
    if (abs(s) < sqrt(.Machine$double.eps)) {
      stop("the optimal weights do not sum to a non-zero value")
    }
    if (abs(s - 1) > 1e-06) {
      warning("the optimal weights summed to ", format(s),
              "; they have been rescaled to sum to one")
    }
    return(w/s)
  }

  ## With bounds, clipping and then rescaling can push a weight back outside
  ## the box it was just clipped into (LB = c(0.6, 0), w = c(0.6, 0.6) would
  ## be returned as c(0.5, 0.5)). Project instead: it satisfies both the box
  ## and the summability constraint by construction.
  wp <- .projectBox(w, ctr$LB, ctr$UB)
  if (max(abs(wp - w)) > 1e-06) {
    warning("the optimal weights did not satisfy the bounds and the summability ",
            "constraint; they have been projected onto them")
  }
  return(wp)
}

.checkGross <- function(w, ctr) {
  ## The gross exposure is a nonlinear constraint handed to the optimizer, and
  ## nothing downstream re-imposes it. A solution that leaves the optimizer
  ## violating it -- for instance because the run was truncated -- must not be
  ## returned as if it were admissible.
  if (ctr$constraint[1] != "gross") {
    return(invisible(w))
  }
  gross <- sum(abs(w))
  if (gross > ctr$gross.c * (1 + 1e-06) + 1e-08) {
    stop("the optimizer returned a portfolio with gross exposure ", format(gross, digits = 6),
         ", above the requested gross.c = ", format(ctr$gross.c),
         "; the optimization did not converge to an admissible portfolio")
  }
  invisible(w)
}
