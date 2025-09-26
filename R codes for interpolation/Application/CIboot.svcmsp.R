#' Bootstrap confidence intervals for spatially varying coefficient models with sign-preservation (SVCM-SP)
#'
#' This function is used to obtain boostrap confidence intervals for the SMVC-SP.
#'
#' @import BPST 
#' @import MGLM 
#' @import osqp
#'
#' @param mfit The fitted result from \code{fit.svcmsp}.
#' \cr
#' @param y The response of dimension \code{n} by one, where \code{n} is the number of observations.
#' \cr
#' @param X The design matrix of dimension \code{n} by \code{p}, with an intercept. Each row is an observation vector.
#' \cr
#' @param B The spline basis function of dimension \code{n} by \code{nT*{(d+1)(d+2)/2}}, where n is the number of observationed points, \code{nT} is the number of triangles in the given triangulation, and \code{d} is the degree of the spline. If some points do not fall in the triangulation, the generation of the spline basis will not take those points into consideration.
#' \cr
#' @param Q2 The Q2 matrix after QR decomposition of the smoothness matrix H.
#' \cr
#' @param B0.pop The spline basis function for population locations. 
#' \cr 
#' @param K The thin-plate energy function.
#' \cr
#' @param lambda The vector of the candidates of penalty parameter -- a chosen lambda from the fitted SVCM is used.
#' \cr
#' @param sign.coeff Sign constraint on bivaraite coefficients. 1, -1 represent positive and negative preservation for each bivariate function, respectively. 0 is used for function without constraint -- default is NULL. For example, assume that there are two bivariate functions. For first and second functions with positive and negative sign-preservation, then \code{c(1,-1)} can be assigned for sign constraints.
#' \cr
#' @param NSVCM The indicator of whether the SVCM-SP is used or not -- \code{TRUE} for the SVCM-SP and \code{FALSE} for the SVCM without sign-constraints for bivariate functions.
#' \cr
#' @param H.accuracy The accuracy of tuning parameters -- \code{TRUE} is used for more accurate stopping criterior for the ADMM method. 
#' \cr
#' @param adj.lam The adjustment for smoothing parameter -- 0.1 is default. 
#' 
#' @return The function returns a list with the following items:
#' \item{Lower.CB}{The lower bound of confidance intervals based on population location.}
#' \item{Upper.CB}{The upper bound of confidance intervals based on population location.}
#' \item{SLower.CB}{The lower bound of confidance intervals based on sample location.}
#' \item{SUpper.CB}{The upper bound of confidance intervals based on sample location.}


CIboot.svcmsp = function(mfit, y, X, B, Q2, B0.pop, K, lambda, nB = 100, sign.coeff = NULL, NSVCM = TRUE, H.accuracy = TRUE, adj.lam = 0.1){
  
  lambda = as.matrix(adj.lam*lambda) 
  nl=nrow(lambda)

  if(ncol(lambda)==1){
    lambda=matrix(rep(lambda,times=np),nl,np)
  }

  BQ2 = as.matrix(B%*%Q2)
  np = ncol(X)
  J = dim(BQ2)[2]

  gamma = as.vector(mfit$gamma)
  
  W = as.matrix(kr(X,B,byrow=TRUE))
  
  star_beta_all_all=list(B)
  star_pop_beta_all_all=list(B)
  
  for(biter in 1:nB){
  
  # Step 1:
    yhat = W%*%gamma

    beta_all = c()
    pop_beta_all = c()
    for(p in 1:np){
      thetamat = matrix(gamma, ncol=np) 
      
      beta = B%*%thetamat[,p]
      beta_pop = B0.pop%*%thetamat[,p]
      
      beta_all = cbind(beta_all, beta)
      pop_beta_all = cbind(pop_beta_all, beta_pop)
    }
    tilte_beta = beta_all
    pop_tilte_beta = pop_beta_all
    hat_eps = y-yhat

    # Step 2:
    set.seed(biter)
    star_eps = sample(hat_eps, length(y), replace=TRUE)
    star_y = yhat + star_eps
    
    # Step 3:
    boot.result = fit.svcmsp(B,Q2,K,lambda, X, star_y,
                             sign.coeff = sign.coeff, NSVCM = NSVCM, H.accuracy = H.accuracy)
    

    star_beta_all = c()
    star_pop_beta_all = c()
    
    for(p in 1:np){
      star_thetamat = boot.result$gamma 
      star_beta = B%*%star_thetamat[,p]
      star_pop_beta = B0.pop%*%star_thetamat[,p]
      
      star_beta_all = cbind(star_beta_all, star_beta)
      star_pop_beta_all = cbind(star_pop_beta_all, star_pop_beta)
    }
    star_tilte_beta = star_beta_all
    star_pop_tilte_beta = star_pop_beta_all
    
    star_beta_all_all[[biter]] = star_tilte_beta 
    star_pop_beta_all_all[[biter]] = star_pop_tilte_beta 
  }

  result_beta_all =c() 
  result_pop_beta_all =c() 
  
  for(biter in 1:nB){
    result_beta_all = cbind(result_beta_all, as.vector(star_beta_all_all[[biter]]))
    result_pop_beta_all = cbind(result_pop_beta_all, as.vector(star_pop_beta_all_all[[biter]]))
  }
  

  dist.beta = sqrt(length(y))*(result_beta_all - as.vector(tilte_beta))
  dist.pop.beta = sqrt(length(y))*(result_pop_beta_all - as.vector(pop_tilte_beta))

  
  lbound = apply(dist.beta, 1, quantile, probs = c(0.025))
  ubound = apply(dist.beta, 1, quantile, probs = c(0.975))
  
  pop.lbound = apply(dist.pop.beta, 1, quantile, probs = c(0.025))
  pop.ubound = apply(dist.pop.beta, 1, quantile, probs = c(0.975))
  
  
  BLB = as.vector(tilte_beta)-1/sqrt(length(y))*ubound
  BUB = as.vector(tilte_beta)-1/sqrt(length(y))*lbound

  pop.BLB = as.vector(pop_tilte_beta)-1/sqrt(length(y))*pop.ubound
  pop.BUB = as.vector(pop_tilte_beta)-1/sqrt(length(y))*pop.lbound
  
  
  
  SLower.CB = matrix(BLB, ncol=np)
  SUpper.CB = matrix(BUB, ncol=np)
  Lower.CB = matrix(pop.BLB, ncol=np)
  Upper.CB = matrix(pop.BUB, ncol=np)

  

  for (indlower in 1:np){
      if((sign.coeff %in% 1)[indlower] == TRUE){
        SLower.CB[SLower.CB[,indlower] < 0,indlower] = 0        
        Lower.CB[Lower.CB[,indlower] < 0,indlower] = 0
        SUpper.CB[SUpper.CB[,indlower] < 0,indlower] = 0     
        Upper.CB[Upper.CB[,indlower] < 0,indlower] = 0   
      }
    }

  for (indupper in 1:np){
      if((sign.coeff %in% (-1))[indupper]  == TRUE){
        SLower.CB[SLower.CB[,indupper] > 0,indupper] = 0         
        Lower.CB[Lower.CB[,indupper] > 0,indupper] = 0   
        SUpper.CB[SUpper.CB[,indupper] > 0,indupper] = 0        
        Upper.CB[Upper.CB[,indupper] > 0,indupper] = 0
      }
    }


  list(Lower.CB = Lower.CB, Upper.CB = Upper.CB, SLower.CB = SLower.CB, SUpper.CB = SUpper.CB)
}

