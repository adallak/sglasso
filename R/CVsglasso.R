
#' Cross-Validation for Simplified Graphical Lasso.
#'
#' The main funcitonality is borrowed from `CVglasso` package.
#'
#' @param X n x p data matrix
#' @param lambda  list of tuning parameters, default is `10^seq(-2, 2, 0.2)`
#' @param niter   maximum number of outer itarations; default is 30
#' @param pen.diag whether to penalize diagonal elements; `TRUE` by default
#' @param tolerance tolerance value for convergence; default is "ave" which the relative difference in the Fronenius norm of the precision matrices across two successive itarations
#' @param thresh  threshold value
#' @param K   numer of fold; default is 5
#' @param crit.cv CV criterian; default is `loglik`
#' @param start algorithm starting criterion; warm by default
#' @param ... other `sglasso()` parameters
#' @return `theta` p x p estimated precision matrix
#' @return `niter` number of iterations until convergence
#' @export
CVsglasso = function(X, lambda = 10^seq(-2, 2, 0.2), niter = 30,
                     tolerance = 1e-04,  thresh = 1e-4,
              pen.diag = FALSE, path = FALSE,
              adjmaxit = NULL, K = 5, crit.cv = c("loglik", "AIC", "BIC"),
              start = c("warm", "cold"), cores = 1, trace = c("progress",
                                                              "print", "none"), ...) {

      # match values
      crit.cv = match.arg(crit.cv)
      start = match.arg(start)
      trace = match.arg(trace)
      lambda = sort(lambda)

      # initialize
      Path = NULL
      initmaxit = niter
      S = cov(X)
      S.train = S.valid = S
      CV_errors = array(0, c(length(lambda), K))

      # set progress bar
      if (trace == "progress") {
            progress = txtProgressBar(max = K * length(lambda), style = 3)
         }

      # no need to create folds if K = 1
      if (K == 1) {

            # set sample size
            n = nrow(S)

            # initialize Path, if necessary
            if (path) {
                  Path = array(0, c(ncol(S), ncol(S), length(lambda)))
            }

      } else {

            # designate folds and shuffle -- ensures randomized folds
            n = nrow(X)
            ind = sample(n)

      }

      # parse data into folds and perform CV
      for (k in 1:K) {
            if (K > 1) {

                  # training set
                  leave.out = ind[(1 + floor((k - 1) * n/K)):floor(k *
                                                                         n/K)]
                  X.train = X[-leave.out, , drop = FALSE]
                  X_bar = apply(X.train, 2, mean)
                  X.train = scale(X.train, center = X_bar, scale = FALSE)

                  # validation set
                  X.valid = X[leave.out, , drop = FALSE]
                  X.valid = scale(X.valid, center = X_bar, scale = FALSE)
                  n = nrow(X.valid)

                  # sample covariances
                  S.train = crossprod(X.train)/(dim(X.train)[1])
                  S.valid = crossprod(X.valid)/(dim(X.valid)[1])

            }

            # re-initialize values for each fold
            maxit = initmaxit
            init = S.train
            initOmega = diag(ncol(S.train))

            # initial sigma
            if (!pen.diag) {

                  # provide estimate that is pd and dual feasible
                  Sminus = S.train
                  diag(Sminus) = 0
                  alpha = min(c(lambda[1]/max(abs(Sminus)), 1))
                  init = (1 - alpha) * S.train
                  diag(init) = diag(S.train)

            }


            # loop over all tuning parameters
            for (i in 1:length(lambda)) {

                  # set temporary tuning parameter
                  lam_ = lambda[i]

                  # update diagonal elements of init, if necessary
                  if (pen.diag) {
                        diag(init) = diag(S.train) + lam_
                  }

                  # compute the penalized likelihood precision matrix
                  # estimator
                  sglasso = sglasso(S = S, lambda = lam_, thr = tolerance,
                                  niter = niter, pen.diag = pen.diag,
                                 init_theta = initOmega,
                                  ...)

                  if (start == "warm") {

                        # option to save initial values for warm starts
                        initOmega = sglasso$theta
                       # maxit = adjmaxit

                  }

                  # compute the observed negative validation loglikelihood
                  # (close enoug)
                  CV_errors[i, k] = (nrow(X)/2) * (sum(sglasso$theta *
                                                             S.valid) - determinant(sglasso$theta, logarithm = TRUE)$modulus[1])

                  # update for crit.cv, if necessary
                  if (crit.cv == "AIC") {
                        CV_errors[i, k] = CV_errors[i, k] + sum(sglasso$theta !=
                                                                      0)
                  }
                  if (crit.cv == "BIC") {
                        CV_errors[i, k] = CV_errors[i, k] + sum(sglasso$theta !=
                                                                      0) * log(nrow(X))/2
                  }

                  # save estimate if path = TRUE
                  if (path) {
                        Path[, , i] = sglasso$theta
                  }

                  # update progress bar
                  if (trace == "progress") {
                        setTxtProgressBar(progress, i + (k - 1) *
                                                length(lambda))

                        # if not quiet, then print progress lambda
                  } else if (trace == "print") {
                        cat("\nFinished lam = ", paste(lam_, sep = ""))
                  }
            }

            # if not quiet, then print progress kfold
            if (trace == "print") {
                  cat("\nFinished fold ", paste(k, sep = ""))
            }
      }

      # determine optimal tuning parameters
      AVG = apply(CV_errors, 1, mean)
      best_lam = lambda[which.min(AVG)]
      error = min(AVG)

      cat("\n")
      # return best lam and alpha values
      return(list("lambda" = best_lam, path = Path, min.error = error,
                  avg.error = AVG, cv.error = CV_errors))

}
