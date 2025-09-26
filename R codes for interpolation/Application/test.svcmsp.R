#' Test global or individual function stationarity.
#'
#' This function is used to test the goodness-of-fit or hypothesises for coefficient functions of the SVCM-SP.
#'
#' @import BPST 
#' @import MGLM 
#' @import osqp
#'

#' @param S The cooridinates of dimension \code{n} by two. Each row is the coordinates of an observation.
#' \cr
#' @param B The spline basis function of dimension \code{n} by \code{nT*{(d+1)(d+2)/2}}, where n is the number of observationed points, \code{nT} is the number of triangles in the given triangulation, and \code{d} is the degree of the spline. If some points do not fall in the triangulation, the generation of the spline basis will not take those points into consideration.
#' \cr
#' @param Q2 The Q2 matrix after QR decomposition of the smoothness matrix H.
#' \cr
#' @param K The thin-plate energy function.
#' \cr
#' @param lambda The vector of the candidates of penalty parameter -- a chosen lambda from the fitted SVCM is used.
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
#' \cr
#' @param nB The number of bootstrap samples -- default is 100.
#' \cr
#' @param test A type of test -- \code{global} for the goodness-of-fit test and \code{individual} for the individual test for stationarity of coefficient functions.
#' \cr
#' @param ind.test An index for individual test for coefficient functions.
#' \cr
#' @return The function returns a list with the following items:
#' \item{b.test}{The observed value of the bootstrap test statistics.}
#' \item{obs.test}{The observed value of the sample test statistics.}
#' \item{pvalue}{The estimated \code{p}-value for the selected test.}

test.svcmsp =
  function(mfit, S,
           B, Q2, K, lambda, X, Y, 
           sign.coeff = NULL, NSVCM = TRUE, H.accuracy = TRUE, 
           nB = 100, test = "global", ind.test = NULL) 
{
  if(!(test == "global" | test == "individual")){
    stop("There are two types of tests. Choose either global test or individulal test.")
  } else {
    
    
    y_bpst =  rowSums(X*mfit$beta)

    yi = Y
    Xi = X
    Si = S
    n = length(yi)
    np = ncol(Xi)

    if(test == "global"){
      # Global Test of Stationarity
      T_boot = NULL
      T_obs = NULL
      pvalue = NULL

      # lse with nonnegative
      QQL = as.matrix(t(Xi) %*% Xi) 
      qqL = as.vector(-t(Xi)%*%Y)
      
      if (NSVCM == TRUE){
        if (H.accuracy == TRUE){
          settings = osqpSettings(sigma = 1e-06, alpha = 1.6, adaptive_rho = 1L, eps_abs= 1e-06, eps_rel = 1e-06, verbose = FALSE)
        } else {
          settings = osqpSettings(sigma = 1e-06, alpha = 1.6, adaptive_rho = 1L, verbose = FALSE)
        }
        QQL = Matrix(QQL, sparse = TRUE)
        QQ2 = Matrix(diag(x = sign.coeff), sparse =TRUE)
        
        GGl.C = solve_osqp(P = QQL, q = qqL, A = QQ2 , l = rep(0,dim(QQ2)[2]), u = NULL,
                          pars = settings)  # check A
        beta_ls = GGl.C$x 
      } else {
        GG.C = solve_osqp(P = QQL, q = qqL, A = NULL, l = NULL, u = NULL,
                          pars = osqpSettings(verbose = FALSE))
        beta_ls = GGl.C$x
      }
      
      y_ls = Xi %*% beta_ls

      sse1 = sum((yi - y_bpst)^2)
      sse2 = sum((yi - y_ls)^2)
      T_obs = (sse2 - sse1) / sse1

      # Step 1
      Ei = yi - y_bpst
      Ei = Ei - mean(Ei)

      for(ib in 1:nB){
        # Step 2
        indp = sample(n, n, replace = TRUE)
        Eb = Ei[indp]
        ystar = y_ls + Eb
        ystar = as.vector(ystar)

        # Step 3
        fit.boot = fit.svcmsp(B, Q2, K, lambda, Xi, ystar,
                              sign.coeff = sign.coeff, NSVCM = NSVCM, H.accuracy = H.accuracy)
        Wb = as.matrix(kr(Xi, B, byrow = TRUE))
        y_bpst2 = Wb %*% as.vector(fit.boot$gamma)

        QQL = as.matrix(t(Xi) %*% Xi) 
        qqL = as.vector(-t(Xi)%*%ystar)
        
        if (NSVCM == TRUE){
          if (H.accuracy == TRUE){
            settings = osqpSettings(sigma = 1e-06, alpha = 1.6, adaptive_rho = 1L, eps_abs= 1e-06, eps_rel = 1e-06, verbose = FALSE)
          } else {
            settings = osqpSettings(sigma = 1e-06, alpha = 1.6, adaptive_rho = 1L, verbose = FALSE)
          }
          QQL = Matrix(QQL, sparse = TRUE)
          QQ2 = Matrix(diag(x = sign.coeff), sparse =TRUE)

          GGl.C2 = solve_osqp(P = QQL, q = qqL, A = QQ2 , l = rep(0,dim(QQ2)[2]), u = NULL,
                             pars = settings)
          betab2 = GGl.C2$x
        } else {
          GG.C2 = solve_osqp(P = QQL, q = qqL, A = NULL, l = NULL, u = NULL,
                            pars = osqpSettings(verbose = FALSE))
          betab2 = GGl.C2$x
        }
        
        y_ls2 = Xi %*% betab2


        sseb1 = sum((ystar - y_bpst2)^2)
        sseb2 = sum((ystar - y_ls2)^2)
        T_boot = c(T_boot, (sseb2 - sseb1) / sseb1)
      }

      # Step 4
      pvalue = mean(T_boot > T_obs)
      list(b.test = T_boot, obs.test = T_obs, pvalue = pvalue)

    } else if (test == "individual") {
      if(is.null(ind.test)){
        warnings("Put ind.test. It is set to 1.")
        ind.test = 1
      }
      T_boot = NULL
      T_obs = NULL
      pvalue = NULL
      
      QQL = as.matrix(t(Xi) %*% Xi) 
      qqL = as.vector(-t(Xi)%*%Y)
      
      
      X.l = X[,ind.test]
      X.nl = X[,-ind.test]
      red.sign.coeff = c(sign.coeff[ind.test],sign.coeff[-ind.test])
      red.lambda =  matrix(lambda[,-ind.test], ncol=np-1) ### hj / CHECK
      
      #reduced coeff:
      Red.result = reduce.fit.svcmsp(B, Q2, K, lambda = red.lambda, X.l, X.nl, ### hj
                                Y, 
                                sign.coeff = red.sign.coeff, 
                                NSVCM = NSVCM, H.accuracy = H.accuracy)
    
      y_ls = Red.result$Yhat
      
      
      
      sse1 = sum((yi - y_bpst)^2)
      sse2 = sum((yi - y_ls)^2)
      T_obs = (sse2 - sse1) / sse1
      
      # Step 1
      Ei = yi - y_bpst
      Ei = Ei - mean(Ei)
      
      for(ib in 1:nB){
        # Step 2
        indp = sample(n, n, replace = TRUE)
        Eb = Ei[indp]
        ystar = y_ls + Eb
        ystar = as.vector(ystar)
        
        # Step 3
        fit.boot = fit.svcmsp(B, Q2, K, lambda, Xi, ystar, 
                              sign.coeff = sign.coeff, NSVCM = NSVCM, H.accuracy = H.accuracy)
        Wb = as.matrix(kr(Xi, B, byrow = TRUE))
        y_bpst2 = Wb %*% as.vector(fit.boot$gamma)
    

        #reduced coeff:
        Red.result2 = reduce.fit.svcmsp(B, Q2, K, lambda = red.lambda, X.l, X.nl,  #hj
                                        ystar,
                                        sign.coeff = red.sign.coeff, 
                                        NSVCM = NSVCM, H.accuracy = H.accuracy)

        
        y_ls2 = Red.result2$Yhat
      
        sseb1 = sum((ystar - y_bpst2)^2)
        sseb2 = sum((ystar - y_ls2)^2)
        T_boot = c(T_boot, (sseb2 - sseb1) / sseb1)
      }
      
      # Step 4
      pvalue = mean(T_boot > T_obs)
      list(b.test = T_boot, obs.test = T_obs, pvalue = pvalue, ind.test = ind.test)
    }
  }
}



reduce.fit.svcmsp = function(B, Q2, K, lambda, X.l, X.nl, 
                             Y,  sign.coeff = NULL, NSVCM = TRUE, H.accuracy = TRUE
){
  
  if(!is.matrix(X.l)){
    warning("The explanatory variable, X.l, should be a matrix.")
    X.l = as.matrix(X.l)
  }
  if(!is.matrix(X.nl)){
    warning("The explanatory variable, X.nl, should be a matrix.")
    X.nl = as.matrix(X.nl)
  }
  
  start.time=Sys.time()
  n=length(Y)
  # np=ncol(X)
  n.l=ncol(X.l)
  n.nl=ncol(X.nl)
  
  np = n.l+n.nl    
  
  J=ncol(Q2)
  P=t(Q2)%*%K%*%Q2
  BQ2=B%*%Q2
  
  if(n.l==0 | is.null(X.l)){
    n.l=0
    X.l=NULL
  }
  
  if(n.nl==0 | is.null(X.nl)){
    n.nl=0
    X.nl=NULL
    XB=NULL
  }else{
    XB=kr(X.nl,BQ2,byrow=TRUE)
  }
  
  lambda=as.matrix(lambda)
  nl=nrow(lambda)
  
  if(ncol(lambda)==1){
    lambda=matrix(rep(lambda,times=np-1),nl,np-1)
  } #hj 
  
  
  Z = as.matrix(cbind(X.l, XB))
  D = matrix(0,(n.l+n.nl*J),(n.l+n.nl*J))
  
  W = Z
  WW = crossprod(W)
  
  
  gamma_all=c()
  beta_all.l=c()
  beta_all.nl=c()
  
  sse_all=c()
  df_all=c()
  gcv_all=c()
  bic_all=c()
  
  
  # il=1
  for(il in 1:nl){
    D[(n.l+1):(n.l+n.nl*J),(n.l+1):(n.l+n.nl*J)]=as.matrix(kronecker(diag(as.vector(lambda[il,]),n.nl),P)) #hj
    Dlam = D # hj
    
    QQ = as.matrix(WW) + as.matrix(Dlam)
    qq = as.vector(-t(W)%*%Y)
  
    if (is.null(sign.coeff)){
      QQ2 = bdiag(1, as.matrix(kronecker(diag(1, np-1), Q2)))
    } else if (!is.null(sign.coeff)){
      QQ2 = bdiag(sign.coeff[1], as.matrix(kronecker(diag(sign.coeff[-1], np-1), Q2)))
    }
    
    QQ = Matrix(QQ, sparse = TRUE)
    QQ2 = Matrix(QQ2, sparse = TRUE)
    
    if (NSVCM == TRUE){
      if (H.accuracy == TRUE){
        settings = osqpSettings(sigma = 1e-06, alpha = 1.6, adaptive_rho = 1L, eps_abs= 1e-06, eps_rel = 1e-06, verbose = FALSE)
      } else {
        settings = osqpSettings(sigma = 1e-06, alpha = 1.6, adaptive_rho = 1L, verbose = FALSE)
      }
      
      
      GG.C = solve_osqp(P = QQ, q = qq, A = QQ2, l = rep(0,dim(QQ2)[1]), u = NULL,
                        pars = settings)
      theta = GG.C$x
    } else {
      GG.C = solve_osqp(P = QQ, q = qq, A = NULL, l = NULL, u = NULL,
                        pars = osqpSettings(verbose = FALSE))
      theta = GG.C$x
    }
    
    lhs=WW+Dlam
    tmp=solve(lhs)%*%t(W)
    Slam=W%*%tmp
    # theta=tmp%*%Y
    
    beta.l= theta[1]
    theta.mtx.nl=matrix(theta[-1],J,np-1)
    gamma=Q2%*%theta.mtx.nl
    
    
    gamma_all=cbind(gamma_all,gamma)
    beta.nl=B%*%gamma
    
    beta_all.l=cbind(beta_all.l,beta.l)
    beta_all.nl=cbind(beta_all.nl,beta.nl)
    
    
    Yhat = X.l*beta.l+ rowSums(X.nl*beta.nl)
    # hist(Y-Yhat)
    # hist(Y)
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
  dim(beta_all.nl)  
  
  gamma.nl=gamma_all[,((np-1)*(j-1)+1):((np-1)*j)]
  beta.nl=beta_all.nl[,((np-1)*(j-1)+1):((np-1)*j)] 
  beta.l = beta_all.l[j] 
  Yhat = X.l*beta.l+ rowSums(X.nl*beta.nl)
  
  sse=sse_all[j];
  gcv=gcv_all[j];
  bic=bic_all[j];	
  df=df_all[j]
  
  end.time=Sys.time()
  time=end.time-start.time
  
  list(beta.l=beta.l,beta.nl, beta.nl, gamma.nl=gamma.nl,
       sse=sse,gcv=gcv,bic=bic,df=df,lambdac=lambdac,time=time, Yhat = Yhat)
}



