#' Choosing the roughness parameter for spatially varying coefficient models (SVCM)
#'
#' This function is used to select roughness parameters of coefficient functions of the SMVC.
#'
#' @import BPST 
#' @import MGLM 
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
svcm.sel=function(B,Q2,K,X,Y, nlambda = 5, 
                    lambda_start = log10(10^(-4)), lambda_end = log10(10^(1)),
                    nlambda2 = 5, scale.lambda2 = 10, verbose = FALSE)
{

  
lambda = 10^(seq(lambda_start,lambda_end,length.out=nlambda))

lamb.Result = plsfitGCV_full(B,Q2,K,lambda,X,Y) 
new.lambda = lamb.Result$lambdac

for (iter.X in 1:dim(X)[2]){
  lambda_start2 = log10(new.lambda[,iter.X]/scale.lambda2)
  lambda_end2 = log10(new.lambda[,iter.X]*scale.lambda2)
  lambda2 = 10^(seq(lambda_start2,lambda_end2,length.out=nlambda2))
  new.lambda = matrix(rep(new.lambda,nlambda2), nlambda2, byrow=T)
  new.lambda[,iter.X] = lambda2
  lamb.Result2 =  plsfitGCV_full(B,Q2,K,new.lambda,X,Y) 
  new.lambda = lamb.Result2$lambdac
  if (verbose == TRUE){
    cat("Lambda_Iter_",iter.X,"_lambda:",new.lambda,"\n")
  }
}

return(new.lambda)

}