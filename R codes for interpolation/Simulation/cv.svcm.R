#' k-fold cross-validation MSPE for spatially varying coefficient regression (SVCM)
#'
#' This function implements k-fold cross-validation MSPE for spartially varying coefficient regression, and returns the mean squared prediction error.
#'
#' @importFrom BPST basis
#' @param y The response of dimension \code{n} by one, where \code{n} is the number of observations.
#' \cr
#' @param X The design matrix of dimension \code{n} by \code{p}, with an intercept. Each row is an observation vector.
#' \cr
#' @param S The sample location \code{n} by two.
#' \cr 
#' @param B The spline basis function of dimension \code{n} by \code{nT*{(d+1)(d+2)/2}}, where n is the number of observationed points, \code{nT} is the number of triangles in the given triangulation, and \code{d} is the degree of the spline. If some points do not fall in the triangulation, the generation of the spline basis will not take those points into consideration.
#' \cr
#' @param Q2 The Q2 matrix after QR decomposition of the smoothness matrix H.
#' \cr
#' @param lambda The vector of the candidates of penalty parameter -- default is grid points of 10 to the power of a sequence from -6 to 6 by 0.5.
#' \cr
#' @param nfolds The number of folds -- default is 10. Although \code{nfold} can be as large as the sample size (leave-one-out CV), it is not recommended for large datasets. Smallest value allowable for \code{nfolds} is 3.
#' \cr
#' @param initial The seed used for cross-validation -- default is 123.
#' \cr
#' @return The mean squared prediction error (MSPE) based on k-fold cross-validation
#' @details
#'
#'
#' @export
#'
cv.svcm =
function(y, X, S,  B, Q2, K,
         lambda = 10^seq(-6, 6, by = 0.5), 
         nfold = 10, initial = 123)
{
  if(nfold < 3){
    warning("The number of folds in CV is too small. Instead, the default 10-fold CV is used.")
    nfold = 10
  }
  if(!is.vector(y)){
    warning("The response variable, y, should be a vector.")
    y = as.vector(y)
  }
  if(!is.matrix(X)){
    warning("The explanatory variable, X, should be a matrix.")
    X = as.matrix(X)
  }
  if(!is.matrix(S)){
    warning("The coordinates, S, should be a matrix.")
    Z = as.matrix(S)
  }

  n = length(y)
  sfold = round(n / nfold)
  set.seed(initial)
  Test = sample(1:n)
  cv.error = c()

  for(ii in 1:nfold){
    if(ii < nfold){
      Test.set = sort(Test[((ii - 1) * sfold + 1):(ii * sfold)])
    }
    if(ii == nfold){
      Test.set = sort(Test[((ii - 1) * sfold + 1):n])
    }
    Train.set = setdiff(1:n, Test.set)

    if(is.vector(X) == 1){
      X.test = as.matrix(X[Test.set])
      X.train = as.matrix(X[Train.set])
    } else {
      X.test = X[Test.set, ]
      X.train = X[Train.set, ]
    }
    B.test = B[Test.set, ]
    B.train = B[Train.set, ]

    y.test = y[Test.set]
    y.train = y[Train.set]


    mfit.ii = plsfitGCV_full(B.train, Q2, K, lambda, X.train, y.train)
    
    W.test = as.matrix(kr(X.test, B.test, byrow = TRUE)) 
    ypred.ii = W.test %*% as.vector(mfit.ii$gamma)
    pred.error = mean((y.test - ypred.ii)^2)  
    cv.error = c(cv.error, pred.error)
  }
  return(cv.error)
}
