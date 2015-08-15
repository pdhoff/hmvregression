#' Bayes Multivariate Linear Regression Posterior Simulation
#'
#' This functions returns a simulation from the posterior distribution
#' of the matrix of multivarate linear regression coefficients.
#'
#' @param Y An \code{m} by \code{n} matrix.
#' @param X A \code{p} by \code{n} matrix.
#' @param B0 An \code{m} by \code{p} matrix.
#' @param V0 A \code{mp} by \code{mp} matrix.
#' @param VE An \code{m} by \code{m} matrix.
#' @param YXt The value of \code{tcrossprod(Y,X)}.
#' @param XXt The value of \code{tcrossprod(X)}.
#' @param iV0b0 The value of \code{ solve(V0)\%*\%c(B0) }.
#' @param iV0 The value of \code{ solve(V0) }.
#' @param iVE The value of \code{ solve(VE) }.
#'
#' @return An \code{m} by \code{p} matrix simulated 
#' from the posterior distribution.
#'
#' @author Peter Hoff
#'
#' @export
#'
#' @examples
#' m<-5 ; p<-3 ; n<-10
#' B<-rsan(c(m,p))
#' X<-rsan(c(p,n))
#' E<-rsan(c(m,n))
#' Y<-B%*%X+E
#'
#' rB_fc(Y,X,matrix(0,m,p),diag(nrow=m*p),diag(m))
#'
#' mvreg_ols(Y,X)
#'
rB_fc<-function(Y,X,B0,V0,VE,
  YXt=tcrossprod(Y,X),
  XXt=tcrossprod(X),
  iV0b0=solve(V0)%*%c(B0),
  iV0=solve(V0),
  iVE=solve(VE))
{

  Vb<-solve( kprod(XXt,iVE)+iV0 )
  Eb<-Vb%*%( c( iVE%*%YXt ) + iV0b0 )
  Sb<-Eb + t(chol(Vb))%*%rnorm(nrow(Y)*nrow(X))
  # could possibly speed up with chol, chol2inv

  matrix(Sb,nrow(Y),nrow(X))
}





#' Bayes Hierarchical Multivariate Linear Regression Posterior Simulation
#'
#' This functions returns a simulation from the posterior distribution
#' of the array of matrix of multivarate linear regression coefficients.
#'
#' @param YA An \code{r} by \code{m} by \code{n} array.
#' @param XA An \code{r} by \code{p} by \code{n} array.
#' @param B0 An \code{m} by \code{p} matrix.
#' @param V0 An \code{mp} by \code{mp} matrix.
#' @param VE An \code{m} by \code{m} matrix.
#'
#' @return An \code{r} by \code{m} by \code{p} array, simulated
#' from the posterior distribution.
#'
#' @author Peter Hoff
#' 
#' @export
#'
rBA_fc<-function(YA,XA,B0,V0,VE)
{
  iVE<-solve(VE) ; iV0<-solve(V0) ; iV0b0<-iV0%*%c(B0)

  BA<-array(dim=c(dim(YA)[1],dim(YA)[2],dim(XA)[2]))
  for(i in 1:dim(YA)[1])
  {
    BA[i,,]<-rB_fc(array(YA[i,,],dim=dim(YA)[-1]),
                   array(XA[i,,],dim=dim(XA)[-1]),
                   B0,iV0b0=iV0b0,iV0=iV0,iVE=iVE)
  }

  BA
}


#' Bayes Hierarchical Multivariate Linear Regression Posterior Simulation
#'
#' This functions returns a simulation from the posterior distribution
#' of the shared error covariance matrix.
#'
#' @param YA An \code{r} by \code{m} by \code{n} array.
#' @param XA An \code{r} by \code{p} by \code{n} array.
#' @param BA An \code{r} by \code{m} by \code{p} array.
#'
#' @return An \code{m} by \code{m} matrix simulated
#' from the posterior distribution.
#'
#' @author Peter Hoff
#'
#' @export
#'
rVE_fc<-function(YA,XA,BA)
{
  S0<-diag( nrow=dim(YA)[2] )

  nu0<-dim(YA)[2]+1

  Sn<-matrix(0,dim(YA)[2],dim(YA)[2])  
  for(i in 1:dim(YA)[1]){ Sn<-Sn+tcrossprod(YA[i,,] - BA[i,,]%*%XA[i,,]) } 

  nun<-prod(dim(YA)[-2])

  solve( rwish( solve(S0+Sn), nu0+nun ) ) 
}


#' Bayes Hierarchical Multivariate Linear Regression Posterior Simulation
#'
#' This functions returns a simulation from the posterior distribution
#' of the across-group regression coefficient covariance.
#'
#' @param BA An \code{r} by \code{m} by \code{p} array.
#' @param B0 An \code{m} by \code{p} matrix. 
#'
#' @return An \code{mp} by \code{mp} matrix simulated
#' from the posterior distribution.
#'
#' @author Peter Hoff
#'
#' @export
#'
rV0_fc<-function(BA,B0)
{ 
  S0<-diag( nrow=prod(dim(B0)) )

  nu0<-prod(dim(B0))+1 

  Sn<-tcrossprod( apply( sweep(BA,c(2,3),B0,"-"),1,c) ) 

  nun<-dim(BA)[1] 

  solve( rwish( solve(S0+Sn), nu0+nun ) )
}

#' Bayes Hierarchical Multivariate Linear Regression Posterior Simulation
#'
#' This functions returns a simulation from the posterior distribution
#' of the across-group mean regression coefficient matrix. 
#'
#' @param BA An \code{r} by \code{m} by \code{p} array.
#' @param V0 An \code{mp} by \code{mp} matrix. 
#' @param iV0 The inverse of \code{V0}. 
#'
#' @return An \code{m} by \code{p} matrix simulated
#' from the posterior distribution.
#'
#' @author Peter Hoff
#'
#' @export
#'
rB0_fc<-function(BA,V0,iV0=solve(V0))
{ 

  Vb0<-solve( iV0*dim(BA)[1] + diag(nrow=nrow(iV0))/dim(BA)[1] ) 
  Eb0<-Vb0%*%iV0%*%apply(apply(BA,1,c),1,sum)

  matrix(Eb0+t(chol(Vb0))%*%rnorm(dim(BA)[2]*dim(BA)[3]),dim(BA)[2],dim(BA)[3]) 
}


