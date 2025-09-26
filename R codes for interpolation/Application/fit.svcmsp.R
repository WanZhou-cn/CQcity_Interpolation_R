#' Fitting spatially varying coefficient models with sign-preservation (SVCM-SP)
#'
#' This function is used to fit the SMVC-SP.
#'
#' @import BPST 
#' @import MGLM 
#' @import osqp
#'
#' @param B The spline basis function of dimension \code{n} by \code{nT*{(d+1)(d+2)/2}}, where n is the number of observationed points, \code{nT} is the number of triangles in the given triangulation, and \code{d} is the degree of the spline. If some points do not fall in the triangulation, the generation of the spline basis will not take those points into consideration.
#' \cr
#' @param Q2 The Q2 matrix after QR decomposition of the smoothness matrix H.
#' \cr
#' @param K The thin-plate energy function.
#' 
#' @param lambda The vector of the candidates of penalty parameter.
#' \cr
#' @param X The design matrix of dimension \code{n} by \code{p}, with an intercept. Each row is an observation vector.
#' \cr
#' @param Y The response of dimension \code{n} by one, where \code{n} is the number of observations.
#' \cr
#' @param sign.coeff Sign constraint on bivaraite coefficients. 1, -1 represent positive and negative preservation for each bivariate function, respectively. 0 is used for function without constraint -- default is NULL. For example, assume that there are two bivariate functions. For first and second functions with positive and negative sign-preservation, then \code{c(1,-1)} can be assigned for sign constraints.
#' \cr
#' @param NSVCM The indicator of whether the SVCM-SP is used or not -- \code{TRUE} for the SVCM-SP and \code{FALSE} for the SVCM without sign-constraints for bivariate functions.
#' \cr
#' @param H.accuracy The accuracy of tuning parameters -- \code{TRUE} is used for more accurate stopping criterior for the ADMM method. 
#' 
#' @return The function returns a list with the following items:
#' \item{beta}{The estimated coefficient functions.}
#' \item{gamma}{The estimated spline coefficient functions.}
#' \item{sse}{Sum of squared errors.}
#' \item{gcv}{Generalized cross-validation (GCV).}
#' \item{df}{Effective degree of freedom.}
#' \item{lambdac}{Selected tuning parameter for bivariate penalized spline based on GCV.}

fit.svcmsp = function(B, Q2, K, lambda, X, Y, sign.coeff = NULL, 
                      NSVCM = TRUE, H.accuracy = TRUE)
{
  start.time=Sys.time()
  n=length(Y)
  np=ncol(X)
  J=ncol(Q2)
  
  BQ2=B%*%Q2
  W=kr(X,BQ2,byrow=TRUE)
  WW=t(W)%*%W
  P=t(Q2)%*%K%*%Q2
  
  lambda=as.matrix(lambda)
  nl=nrow(lambda)
  if(ncol(lambda)==1){
    lambda=matrix(rep(lambda,times=np),nl,np)
  }
  
  gamma_all=c()
  beta_all=c()
  sse_all=c()
  df_all=c()
  gcv_all=c()
  bic_all=c()
  
  for(il in 1:nl){
    Lam=diag(lambda[il,])
    Dlam=kronecker(Lam,P)
    
    QQ = as.matrix(WW) + as.matrix(Dlam)
    qq = as.vector(-t(W)%*%Y)
    
    if (is.null(sign.coeff)){
      QQ2 = as.matrix(kronecker(diag(1, np), Q2))
    } else if (!is.null(sign.coeff)){
      QQ2 = as.matrix(kronecker(diag(sign.coeff, np), Q2))
    }
    
    QQ = Matrix(QQ, sparse = TRUE)
    QQ2 = Matrix(QQ2, sparse = TRUE)
    
    if (NSVCM == TRUE){
      if (H.accuracy == TRUE){
        settings = osqpSettings(sigma = 1e-06, alpha = 1.6, adaptive_rho = 1L, 
                                eps_abs= 1e-06, eps_rel = 1e-06, verbose = FALSE)
      } else {
        settings = osqpSettings(sigma = 1e-06, alpha = 1.6, adaptive_rho = 1L, verbose = FALSE)
      }
      GG.C = solve_osqp(P = QQ, q = qq, A = QQ2, l = rep(0,dim(QQ2)[1]), u = NULL,
                       pars = settings)

    theta = GG.C$x
    } else {
    GG.C = solve_osqp(P = QQ, q = qq, A = NULL, l = NULL, u = NULL,
                     pars = osqpSettings())
    theta = GG.C$x
    }
    
    lhs=WW+Dlam
    tmp=solve(lhs)%*%t(W)
    Slam=W%*%tmp
    theta.mtx=matrix(theta,J,np)
    gamma=Q2%*%theta.mtx
    gamma_all=cbind(gamma_all,gamma)
    beta=B%*%gamma
    beta_all=cbind(beta_all,beta)
    Yhat=rowSums(X*beta)
    sse=sum((Y-Yhat)^2)
    sse_all=c(sse_all,sse)
    df=sum(diag(Slam))
    df_all=c(df_all,df)
    gcv=n*sse/(n-df)^2
    gcv_all=c(gcv_all,gcv)
    bic=log(sse/n)+df*log(n)/n
    bic_all=c(bic_all,bic)
  }
  j=which.min(gcv_all)
  lambdac=lambda[j,]

  gamma=gamma_all[,(np*(j-1)+1):(np*j)]
  beta=beta_all[,(np*(j-1)+1):(np*j)]
  sse=sse_all[j];
  gcv=gcv_all[j];
  bic=bic_all[j];	
  df=df_all[j]
  
  end.time=Sys.time()
  time=end.time-start.time
  
  list(beta=beta, gamma=gamma, sse=sse, gcv=gcv, df=df,
       lambdac = matrix(lambdac, ncol=np))
}

