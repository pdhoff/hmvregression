#' OLS Multivariate Linear Regression Estimates
#' 
#' This functions returns the OLS estimate of the regression 
#' matrix in the multivariate regression model. 
#' 
#' @param Y An \code{m} by \code{n} matrix. 
#' @param X A \code{p} by \code{n} matrix.
#' 
#' @return A \code{m} by {p} matrix of OLS regression estimates 
#' for the model Y=BX + E. 
#' 
#' @examples
#' m<-5 ; p<-3 ; n<-10
#' B<-rsan(c(m,p))
#' X<-rsan(c(p,n)) 
#' E<-rsan(c(m,n)) 
#' Y<-B%*%X+E 
#' 
#' t( lm( t(Y)~ -1 + t(X))$coef )  
#' 
#' mvreg_ols(Y,X) 
#' 
#' @author Peter Hoff
#'
#' @export
mvreg_ols<-function(Y,X)
{
  tcrossprod(Y,X) %*% solve(tcrossprod(X))
}




