#' Utility function to compare true and estimated matrices
#'
#'
#' @param est p x p estimated square matrix
#' @param true p x p true square matrix
#' @param report  report result as a table; `TRUE` by default
#' @export
compareMat <- function(est, true, report = TRUE){
  est = makeAdj(est)
  true = makeAdj(true)
  if((nrow(est) != ncol(est)) || nrow(true) != ncol(true)){
    stop("input matrices shold be square")
  }
  if(!isSymmetric(est) || !isSymmetric(true)){
    stop("input matrices should be symmetric")
  }
  if(nrow(est) != ncol(true)){
    stop("input matrices should have the same dimensions")
  }
  diffm = est - true
  numTrueGaps = (sum(true == 0)) / 2
  if (numTrueGaps == 0){
    fpr = 1
  } else {
    fpr = (sum(diffm >0) / 2) / numTrueGaps
  }
  diffm2 = true - est
  nmbTrueEdges = sum(true == 1) / 2
  if (nmbTrueEdges == 0)
  {
    tpr = 0
  }else{
    tpr = 1 - (sum(diffm2 > 0) / 2) / nmbTrueEdges
  }
  trueEstEdges = nmbTrueEdges - sum(diffm2 > 0) / 2
  if (nmbTrueEdges == 0)
  {
    if(trueEstEdges == 0)
    {
      tdr = 1
    }else{
      tdr = 0
    }
  }else{
    tdr = trueEstEdges / (sum(est == 1) / 2)
  }
  table = matrix(0, nrow = 1, ncol = 5)
  table[1,1] = tpr
  table[1,2] = fpr
  table[1,3] = tdr
  table[1,4] = sum(true)
  table[1,5] = sum(est)
  colnames(table) = c("TPR","FPR", "TDR", "True non-zero", "Est non-zero")
  if(isTRUE(report)){
    print(table)
  }
  return(list("tpr" = tpr, "fpr" = fpr, "tdr" = tdr,
              "true_nzero" = sum(true), "est_nzero"=sum(est)))
}

makeAdj <- function(mat){
  mat = abs(mat)
  mat[mat > 0] = 1
  diag(mat) = 0
  return(mat)
}
