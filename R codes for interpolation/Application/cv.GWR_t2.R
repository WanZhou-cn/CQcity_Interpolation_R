#' k-fold cross-validation MSPE for geographically weighted regression (GWR)
#'
#' This function implements k-fold cross-validation MSPE forgeographically weighted regression, and returns the mean squared prediction error.
#'
#' @importFrom BPST basis
#' @param y The response of dimension \code{n} by one, where \code{n} is the number of observations.
#' \cr
#' @param X The design matrix of dimension \code{n} by \code{p}, with an intercept. Each row is an observation vector.
#' \cr
#' @param S The sample location \code{n} by two.
#' \cr
#' @param nfolds The number of folds -- default is 10. Although \code{nfold} can be as large as the sample size (leave-one-out CV), it is not recommended for large datasets. Smallest value allowable for \code{nfolds} is 3.
#' \cr
#' @param initial The seed used for cross-validation -- default is 123.
#' \cr
#' @return The mean squared prediction error (MSPE) based on k-fold cross-validation
cv.GWR =
function(y, X, S, iTest,
         nfold = 10, initial = 123)
{
 
  if(nfold < 3){
    warning("The number of folds in CV is too small. Instead, the default 10-fold CV is used.")
    nfold = 10
  }

  if(!is.matrix(X)){
    warning("The explanatory variable, X, should be a matrix.")
    X = as.matrix(X)
  }
  if(!is.matrix(S)){
    warning("The coordinates, S, should be a matrix.")
    S = as.matrix(S)
  }
  
  n = length(y)

  sfold = round(n / nfold)
  # set.seed(initial)
  # Test = sample(1:n)
  Test = iTest
  cv.error = c()
  list_out = list()
  mat_val0 = matrix(NA, nrow=n, ncol=4)
  
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
      S.test = as.matrix(S[Test.set])
      S.train = as.matrix(S[Train.set])
    } else {
      X.test = X[Test.set, ]
      X.train = X[Train.set, ]
      S.test = S[Test.set, ]
      S.train = S[Train.set, ]
    }

    y.test = y[Test.set]
    y.train = y[Train.set]
    
    # General formula:
    cov.names = paste0("X",1:dim(X)[2])
    tmp = paste0("y ~ -1")
    formula.GWR = as.formula(paste0(tmp, 
                                    paste0('+', cov.names, collapse = '')))
    cat("formula = ",as.character(formula.GWR),"\n")
    
    
    dat=data.frame(y.train,S.train,X.train)
    colnames(dat)=c('y','u','v',cov.names)
    coordinates(dat)=c("u","v")
    
    bw=gwr.sel(formula.GWR,data=dat)
    # gwr model fitting
    model.fit=gwr(formula.GWR,data=dat,bandwidth=bw, hatmatrix=T)
    new=data.frame(X.test, S.test)
    names(new)=c(cov.names,'u','v') 
    coordinates(new)=c("u","v")
    
    model.pred=gwr(formula.GWR,data=dat,bandwidth=bw,fit.points=new,
                   predictions=TRUE,se.fit=T,fittedGWRobject=model.fit)
    
    y.pred=model.pred$SDF$pred

    pred.error = mean((y.test - y.pred)^2) 
    cv.error = c(cv.error, pred.error)
    
    ##### record of Y values
    mat_val0[Test.set, 1] = y.test
    mat_val0[Test.set, 2] = y.pred
    mat_val0[Test.set, 3] = y.test - y.pred
    mat_val0[Test.set, 4] = ii
  }
  list_out[[1]] = cv.error
  list_out[[2]] = mat_val0
  return(list_out) 
}

