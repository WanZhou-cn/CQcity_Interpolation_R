#' Fitting spatially varying coefficient models (SVCM)
#'
#' This function is used to fit the SMVC.
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
#' @param lambda The vector of the candidates of penalty parameter.
#' \cr
#' @param X The design matrix of dimension \code{n} by \code{p}, with an intercept. Each row is an observation vector.
#' \cr
#' @param Y The response of dimension \code{n} by one, where \code{n} is the number of observations.
#' \cr
#' @return The function returns a list with the following items:
#' \item{beta}{The estimated coefficient functions.}
#' \item{gamma}{The estimated spline coefficient functions.}
#' \item{sse}{Sum of squared errors.}
#' \item{gcv}{Generalized cross-validation (GCV).}
#' \item{df}{Effective degree of freedom.}
#' \item{lambdac}{Selected tuning parameter for bivariate penalized spline based on GCV.}

plsfitGCV_full=function(B,Q2,K,lambda,X,Y){
  start.time=Sys.time()
	n=length(Y)
	np=ncol(X)
	J=ncol(Q2)

	BQ2=B%*%Q2
	W=kr(X,BQ2,byrow=TRUE)
	WW=t(W)%*%W
	#rhs=t(W)%*%Y
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
		lhs=WW+Dlam
		tmp=solve(lhs)%*%t(W)
		Slam=W%*%tmp
		theta=tmp%*%Y
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
	
	list(beta=beta,gamma=gamma,sse=sse,gcv=gcv,df=df,lambdac=matrix(lambdac, ncol=np))
}

