#' Choosing the roughness parameter for spatially varying coefficient models with sign-preservation (SVCM-SP)
#'
#' This function is used to select roughness parameters of coefficient functions of the SMVC-SP.
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
#' @param nlambda The number of candidates for common roughness parameter
#' \cr
#' @param lambda_start The starting value of the sequence for common roughness parameter.
#' \cr
#' @param lambda_end The end value of the sequence for common roughness parameter.
#' \cr
#' @param nlambda2 The number of candidates for each roughness parameter of coefficient functions.
#' \cr
#' @param scale.lambda2 The value to adjust the starting and end values of the sequence for each roughness parameter of coefficient functions. 
#' \cr    
#' @param verbose \code{TRUE} is used to write out progress.
#' \cr
#' @return The chosen lambda based on the GCV.
#' 
svcmsp.sel=function(B, Q2, K, X, Y, sign.coeff = NULL, 
                    NSVCM = TRUE, H.accuracy = TRUE, nlambda = 5, 
                    lambda_start = log10(10^(-4)), lambda_end = log10(10^(1)),
                    nlambda2 = 5, scale.lambda2 = 10, verbose = FALSE)
{
# nlambda = 5  
# lambda_start = log10(10^(-4))
# lambda_end = log10(10^(1))
lambda = 10^(seq(lambda_start,lambda_end,length.out=nlambda))

lamb.Result = fit.svcmsp(B, Q2, K, lambda, X, Y, sign.coeff, NSVCM, H.accuracy) 
new.lambda = lamb.Result$lambdac

for (iter.X in 1:dim(X)[2]){
  # nlambda2 = 5
  lambda_start2 = log10(new.lambda[,iter.X]/scale.lambda2)
  lambda_end2 = log10(new.lambda[,iter.X]*scale.lambda2)
  lambda2 = 10^(seq(lambda_start2,lambda_end2,length.out=nlambda2))
  new.lambda = matrix(rep(new.lambda,nlambda2), nlambda2, byrow = T)
  new.lambda[,iter.X] = lambda2
  lamb.Result2 = fit.svcmsp(B, Q2, K, lambda = new.lambda, 
                            X, Y, sign.coeff, NSVCM, H.accuracy) 
  new.lambda = lamb.Result2$lambdac
  if (verbose == TRUE){
    cat("Lambda_Iter_",iter.X,"_lambda:",new.lambda,"\n")
  }
}


return(new.lambda)

}