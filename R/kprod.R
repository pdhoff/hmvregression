#' Kronecker Product
#' 
#' The function \code{kprod} calculates the Kronecker product of two 
#' or more matrices. 
#' 
#' @param \ldots Either a sequence of matrices or a single list of matrices. 
#' 
#' @return A matrix with dimensions equal to the product of the input dimensions. 
#' 
#' @examples
#' A<-rsan(c(2,3)) ; B<-rsan(c(3,4)) ; C<-rsan(c(4,5))
#' ABC<-kprod(A,B,C) 
#' dim(ABC)
#' dim(A)*dim(B)*dim(C) 
#'  
#' M<-list(A,B,C)
#' kprod(M) 
#' 
#' @author Peter Hoff
#' 
#' @export 
kprod<-function(...)
{
  M<-list(...)
  if(is.list(M[[1]])) { M<-M[[1]] }
  KM<-1 ; for(k in 1:length(M)){ KM<-kronecker(KM,M[[k]]) }
  KM
}

