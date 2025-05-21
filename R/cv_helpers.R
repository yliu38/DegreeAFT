#' Uno-AUC + C-index for one validation fold
#'
#' @keywords internal
metric_fold <- function(time_tr, status_tr, lp_tr,
                        time_va, status_va, lp_va) {

  Surv_tr <- survival::Surv(time_tr, status_tr)
  Surv_va <- survival::Surv(time_va, status_va)

  ## 1. Harrell's C-index
  cidx <- survival::concordance(Surv_va ~ lp_va)$concordance

  ## 2. Integrated Uno-AUC
  tau  <- seq(0, max(time_tr), length.out = 100)
  iauc <- survAUC::AUC.uno(
            Surv.rsp     = Surv_tr,
            Surv.rsp.new = Surv_va,
            lpnew        = lp_va,
            times        = tau)$iauc

  c(cidx = cidx, iauc = iauc)
}

#' K-fold cross-validation for any penalised AFT fit function
#'
#' @param time    survival times (numeric).
#' @param status  0/1 event indicator.
#' @param X       centred & scaled design matrix.
#' @param penMat  penalty matrix passed straight to `fitfun`.
#' @param fitfun  fitting function (*e.g.* [prox_deg()] or [prox_lap()]).
#' @param grid    `data.frame` with columns `lam1`, `lam2` (penalty grid).
#' @param K       number of folds (default 5).
#' @param seed    random seed for reproducibility.
#'
#' @return A list with:  
#' * **best_lambda** – single-row `data.frame` with the chosen \eqn{\lambda_1,\lambda_2};  
#' * **cv_metrics** – named numeric vector (`cidx`, `iauc`) averaged over folds.
#' @export
cv_select <- function(time, status, X, penMat, fitfun,
                      grid, K = 5, seed = 42) {
  set.seed(seed)
  n     <- length(time)
  folds <- sample(rep(1:K, length.out = n))

  score_one <- function(f, idx) {                # fold f, λ‑row idx
    tr <- folds != f;  va <- folds == f

    ## --------- NEW: compute train statistics and scale ----------
    mu <- colMeans(X[tr, , drop = FALSE])
    sd <- apply(X[tr, , drop = FALSE], 2, sd)
    sd[sd == 0] <- 1                    # avoid divide‑by‑zero

    X_tr <- sweep(sweep(X[tr, ], 2, mu,  `-`), 2, sd, `/`)
    X_va <- sweep(sweep(X[va, ], 2, mu,  `-`), 2, sd, `/`)
    ## ------------------------------------------------------------

    beta <- fitfun(time[tr], status[tr], X_tr, penMat,
                   grid$lam1[idx], grid$lam2[idx])

    lp_tr <- as.vector(X_tr %*% beta)
    lp_va <- as.vector(X_va %*% beta)

    metric_fold(time[tr], status[tr], lp_tr,
                time[va], status[va], lp_va)
  }


  res_grid <- future_sapply(1:nrow(grid), function(idx) {
    ## matrix  (2 rows = metrics, K columns = folds)
    fold_stats <- future_sapply(1:K, score_one, idx)

    ## average over folds – preserves row names ("cidx", "iauc")
    rowMeans(fold_stats, na.rm = TRUE)
  })

  ##  make sure row names exist
  if (is.null(rownames(res_grid)))
    rownames(res_grid) <- c("cidx", "iauc")

  best <- order(res_grid["cidx", ], res_grid["iauc", ], decreasing = TRUE)[1]
  list(best_lambda = grid[best, ],
       cv_metrics  = res_grid[, best])
}
