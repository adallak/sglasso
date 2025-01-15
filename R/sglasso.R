#' Simplified Graphical Lasso
#'
#' This function implements GLasso using modifications proposed in
#' Dallakyan and Pourahmadi (2024)
#'
#' @param X n x p data matrix
#' @param init_theta p x p initial positive definite precision matrix
#' @param lambda  tuning parameter, default is 0.1
#' @param pen.diag whether to penalize diagonal elements; `TRUE` by default
#' @param niter   maximum number of outer itarations; default is 30
#' @param inner_niter maximum number of inner iterations; default is 1000
#' @param tolerance tolerance value for convergence;
#' @param inner.tolerance tolerance value for convergence of the inner loop;
#' @param tol.type default is "max" which the relative difference of the maximum of precision matrices across two successive iterations
#' @param thresh  threshold value
#' @param trace.objfnct whether to compute objective function. Increases computation time
#' @return `theta` p x p estimated precision matrix
#' @return `niter` number of iterations until convergence
#' @export
sglasso <-function(S, init_theta = NULL, lambda = 0.1,
                   pen.diag = TRUE, niter = 30, inner_niter = 1000,
                   tolerance=1e-5, inner.tolerance = 1e-7,
                   thresh = 1e-4, tol.type = c("max", "ave", "F"),
                   trace.objfnct = FALSE){
   tol.type = match.arg(tol.type)
   p = ncol(S)
   if(p != nrow(S)){
      stop("matrix S should be symmetric")
   }
   if (is.null(init_theta)){
      theta = diag(p)
      diag(theta) = 1/(diag(S) + lambda)
   } else{
      theta = init_theta
   }
   converge = FALSE
   iter = 0
   result = sglasso_R(S = S, lambda = lambda, init_theta = init_theta,
                         pen.diag = pen.diag, niter = niter, inner_niter = inner_niter,
                         tolerance=tolerance, tol.type = tol.type,
                         inner.tolerance = inner.tolerance,
                         thresh = thresh, trace.objfnct = trace.objfnct)
   if (niter == result$niter){
      warning("Algorithm doesn't converge")
   }
   ind_thresh = which(abs(result$theta) < thresh)
   result$theta[ind_thresh] = 0
   if(isTRUE(trace.objfnct)){
      return(list("theta" = theta, "error" = t(result$error),
                  "niter" = result$niter, "trace.objfnct" = result$trace.objfnct))
   } else {
      return(list("theta" = result$theta, "error" = t(result$error),
                  "niter" = result$niter))
   }
}




sglasso_R <- function(X = NULL, S = NULL,  lambda = 0.1,
                              init_theta = NULL, niter = 30, inner_niter = 1000,
                              tolerance=1e-5, tol.type = c("max","ave", "F"),
                              inner.tolerance = 1e-7,
                              pen.diag = TRUE, standardize = FALSE,
                              thresh = 1e-5, trace.objfnct = FALSE)
{
      if (is.null(X) && is.null(S)){
            stop("one of X or S should be non-null")
      }
      if (!is.null(X) && !is.null(S)){
            stop("only one of X or S should be provided")
      }
      if(!is.null(X)){
            if (isTRUE(standardize)){
                  X = scale(X, center = TRUE, scale = TRUE)
            }
            S = cov(X)
      }
      tol.type = match.arg(tol.type)
      p = ncol(S)
      if(p !=  nrow(S)) {
         stop("matrix S should be symmetric")
      }
      if(lambda < 0) {
         stop("lambda should be positive")
      }
      if(niter < 0 || inner_niter < 0) {
         stop("number of iterations should be positive")
      }
      if(tolerance < 0 || inner.tolerance < 0){
         stop("tolerance should be positive")
      }
      if (is.null(init_theta)){
            theta = diag(p)
            diag(theta) = 1/(diag(S) + lambda)
      } else{
         if (nrow(init_theta) != ncol(init_theta)){
            stop("initial matrix theta should be symmetric")
         }
            theta = init_theta

      }
      converge = FALSE
      iter = 0
      error_store = numeric(niter)
      loglik.trace = c()
      U.mat = diag(rep(lambda,p))
      while((converge == FALSE) && (iter <= niter)){
            theta_old = theta
            for(i in 1:p){
                  s_12 = S[-i, i]
                  s_22 = S[i, i]

                  converge_beta = FALSE
                  beta = matrix(1, nrow = p-1, ncol = 1)
                  if (isTRUE(pen.diag)){
                        denom = s_22 + lambda
                  } else {
                        denom = s_22
                  }
                  gamma = 1 / denom
                  res = box_qp_f(theta[-i, -i], U.mat[-i,i], s_12, lambda,
                                       inner_niter,
                                       inner.tolerance)
                  G = res$grad_vec
                  u = res$u
                  beta = -G/ denom
                  ## Update theta
                  theta[-i, i] = beta
                  theta[i, -i] = beta
                  theta[i, i] = gamma - sum((s_12 + u) * beta) / denom

            }
            err = compute_stopping(theta, theta_old, tolerance, S, tol.type)
            error_store[iter] = err$val
            check = err$check
            if(isTRUE(trace.objfnct)){
                  loglik.trace[iter] = primal_obj(theta, S, lambda)
            }
            if(check){
                  converge = TRUE
            } else{
                  iter = iter + 1
            }
      }
      if (iter == niter){
            warning("Algorithm doesn't converge")
      }
      ind_thresh = which(abs(theta) < thresh)
      theta[ind_thresh] = 0
      ind = which(error_store != 0)
      error_store = error_store[ind]
      niter = length(error_store)
      if(isTRUE(trace.objfnct)){
            return(list("theta" = theta, "error" = error_store,
                        "niter" = niter, "trace.objfnct" = loglik.trace))
      }else{
            return(list("theta" = theta, "error" = error_store, "niter" = niter))
      }
}

compute_theta_22 <- function(gamma, beta, theta_11_inv){
      return(gamma + t(beta) %*% theta_11_inv %*% beta)
}

soft <- function(a, b){
      return(sign(a) * pmax((abs(a) - b), 0 ))
}

compute_stopping <- function(new, old, tol, S, type = c("max","F", "ave")){
      crit = match.arg(type)
      if(crit == "F"){
         val = norm(new - old)/ norm(old)
         tr = (val < tol)
      } else if (crit == "max") {
         val = max(abs(new - old)) / max(abs(old))
         tr = (val < tol)
      } else {
         diag(new) = 0
         diag(old) = 0
         val = mean(abs(new - old))
         S_off = S
         diag(S_off) = 0
         tr = (val < (tol * mean(abs(S_off))))
      }
  #    err = primal_obj(new, S, lambda) - primal_obj(old, S, lambda)
      return(list("check" = tr, "val" = val))
}


dual_update<- function(Theta, s, lambda, tol , inner_niter) {
      p = length(s)
      converge = FALSE
      iter = 0
      u = matrix(1, nrow = p, ncol = 1)
      b = Theta %*% s
      min = -lambda
      max = lambda
      while((converge == FALSE) && (iter <= inner_niter)){
         u_old = u
         for(i in 1:p){
            z_i = (-b[i] - sum(Theta[-i,i] * u[-i])) / Theta[i,i]
            u[i] = proj_oper(min, max, z_i)
         }
         if(sqrt(sum(abs(u - u_old)^2)) < tol){
            converge = TRUE
         } else{
            iter = iter + 1
         }
      }
      if (iter == inner_niter){
         warning("Inner update doesn't converge")
      }
      return(u)
}

dual_obj <- function(theta, s, u){
      return(1/2 * t(s + u) %*% theta %*% (s + u))
}

primal_obj <- function(Theta, S, lambda){
      return(-log(det(Theta)) + sum(diag(S %*% Theta)) + lambda * sum(abs(Theta)))
}

proj_oper <- function(min, max, z){
   if(z < min){
      return(min)
   }
   if (z > max){
      return(max)
   }
   return(z)
}



