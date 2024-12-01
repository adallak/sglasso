#' Quadratic programming with box constrains
#'
#' This function is borrowed from dpglasso package
#'
#' @param Q symmetric p x p  matrix
#' @param u optimization variable, a vector of length p
#' @param b   is a vector of length p, this is a problem parameter.
#' @param rho tuning parameter
#' @param Maxiter   maximum number of itarations; default is 10^3
#' @param tol    convergence torlerance; default is 10^-4
#' @useDynLib sglasso, .registration = TRUE
#' @export
box_qp_f <- function (Q, u,b,rho,Maxiter=10^3,tol=10^-4) {
  qq=nrow(Q);
  obj_vals<-rep(0,Maxiter);

  mode(Q)="double";mode(u)="double";mode(b)="double";
  mode(Maxiter)="integer";mode(tol)="double";
  mode(rho)="double";
  mode(obj_vals)="double";
  mode(qq)="integer";

  junk<-.Fortran("box_qp_f",Q,uu=u,b,rho,Maxiter,tol,qq,grad_vec=double(qq))

  return(list(grad_vec=junk$grad_vec,u=junk$uu))

}
