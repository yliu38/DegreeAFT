#' Soft-thresholding operator
#'
#' Internal helper used by the proximal gradient updates.
#'
#' @keywords internal
soft <- function(u, t) sign(u) * pmax(abs(u) - t, 0)

#' Gradient of the (negative) log-likelihood in a Weibull AFT model
#'
#' @keywords internal
gradient.Weibull <- function(X, Y, eta, delta, sigma, n) {
  res <- (Y - eta) / sigma
  -crossprod(X, exp(res) - delta) / (n * sigma)
}

#' Inverse-Degree penalised Weibull AFT (Main method)
#'
#' Fits a sparse Weibull AFT model with an \eqn{\ell_1} penalty
#' plus an **inverse-degree** network penalty:
#'
#' @param time   numeric vector of survival times.
#' @param delta  numeric 0/1 censoring indicator (1 = event).
#' @param X      `n × p` design matrix (samples × genes).
#' @param DegInv sparse **p × p** inverse-degree matrix (class `dgCMatrix`).
#' @param l1     non-negative scalar – L1 penalty weight.
#' @param l2     non-negative scalar – network penalty weight.
#' @param itmax  maximum iterations for proximal gradient (default 2000).
#' @param tol    convergence tolerance on the max‐coef change (default 1e-6).
#'
#' @return numeric vector of length *p* with the fitted coefficients.
#' @export
prox_deg <- function(time, delta, X, DegInv, l1, l2,
                     itmax = 2000, tol = 1e-6) {
  n <- length(time); y <- log(time); p <- ncol(X)
  beta <- numeric(p); sigma <- 1
  for (it in seq_len(itmax)) {
    eta <- X %*% beta
    L   <- (max(exp((y - eta) / sigma)) / n) * norm(X, "2")^2 +
      l2 * max(diag(DegInv))
    g   <- gradient.Weibull(X, y, eta, delta, sigma, n) + l2 * (DegInv %*% beta)
    b   <- soft(beta - g / L, l1 / L)
    if (max(abs(b - beta)) < tol) break
    beta <- b
  }
  beta
}

#' Laplacian penalised Weibull AFT (Reference method)
#'
#' Identical interface to [prox_deg()] but uses a graph **Laplacian**
#'
#' @inheritParams prox_deg
#' @param Lap sparse **p × p** unnormalised Laplacian matrix.
#' @return numeric vector of length *p* with the fitted coefficients.
#' @export
prox_lap <- function(time, delta, X, Lap, l1, l2,
                     itmax = 2000, tol = 1e-6) {
  n <- length(time); y <- log(time); p <- ncol(X)
  beta <- numeric(p); sigma <- 1
  for (it in seq_len(itmax)) {
    eta <- X %*% beta
    L   <- (max(exp((y - eta) / sigma)) / n) * norm(X, "2")^2 +
      l2 * max(rowSums(abs(Lap)))
    g   <- gradient.Weibull(X, y, eta, delta, sigma, n) + l2 * (Lap %*% beta)
    b   <- soft(beta - g / L, l1 / L)
    if (max(abs(b - beta)) < tol) break
    beta <- b
  }
  beta
}
