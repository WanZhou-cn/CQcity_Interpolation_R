#' Data generating function. 
#'
#' This function generates data in the simulation.

library(mgcv)

data.svcmsp <- function(iter=123){
  fsb=list(fs.boundary())
  
  uu=seq(-1,3.5,0.02)
  vv=seq(-1,1,0.02)
  n1=length(uu)
  n2=length(vv)
  u=rep(uu,n2)
  v=rep(vv,rep(n1,n2))
  m1=fs.test(u,v,b=1)
  N=length(m1)
  m1=m1- min(m1[is.na(m1)!=1])
  m2=m1
  m2[!is.na(m1)] = 2*0.05*pi*(u[!is.na(m1)]^2+v[!is.na(m1)]^2)
  
  set.seed(iter)
  eps = rnorm(N,0,1)
  x1=rep(1,N)
  x2=runif(N,0,1)

  m1.mtx=matrix(m1,n1,n2)
  m2.mtx=matrix(m2,n1,n2)
  y=m1*x1+m2*x2+eps

  pop=cbind(y,m1,m2,x1,x2,u,v)
  return(pop)
}



