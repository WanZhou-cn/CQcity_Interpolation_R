rm(list=ls())

# Library:
  # install.packages("devtools",repos="http://cran.r-project.org")
  # library(devtools)
  # install_github("funstatpackages/BPST")
  # install_github("funstatpackages/Triangulation")


library('MGLM') 
library('BPST')
library('Triangulation')
library('spgwr') 


# Source functions
source("plsfitGCV_full.R")
source("cv.svcm.R")
source("cv.GWR.R")
source("svcm.sel.R")

S.pop = as.matrix(read.csv('data/pop_location_d025.csv',header=T))
bb = as.matrix(read.csv('data/usa_bb.csv', header=FALSE))

# Triangulation:
Tr=as.matrix(read.csv('data/Tr_usa.csv',header=FALSE))
V=as.matrix(read.csv('data/V_usa.csv',header=FALSE))


# Data:
data = as.data.frame(read.csv('data/AirTempData.csv',header=TRUE))
Y = data$at
X = cbind(data$plst, data$dem)
X = scale(X, scale=TRUE) # standardized covariates
X = cbind(1, X)          # add intercept
S = cbind(data$x,data$y)


u= unique(S.pop[,1])
v= unique(S.pop[,2])
n1=length(u)
n2=length(v)

# Basis
d=3; r=1;

B0.pop=basis(V, Tr, d, r, S.pop)
Q2=B0.pop$Q2
B.pop=B0.pop$B
BQ2.pop=as.matrix(B.pop%*%Q2)
K=B0.pop$K
P=t(Q2)%*%K%*%Q2

# Sample basis:
B0 = basis(V, Tr, d, r, S)
B = B0$B


# SVCM:
lambda.svcm = svcm.sel(B, Q2, K, X, Y)
Result.svcm = plsfitGCV_full(B, Q2, K, lambda.svcm, X, Y) 
yhat.svcm = rowSums(X*Result.svcm$beta)
beta.hat.pop.svcm = B.pop%*%Result.svcm$gamma
mhat_all.svcm = beta.hat.pop.svcm


# 10CV-MSPE (SVCM):
CV.Result.svcm = cv.svcm(y = Y, X = X, S = S, B = B, Q2 = Q2, K = K, lambda = Result.svcm$lambdac)
CV.Result.svcm
mean(CV.Result.svcm)



# GWR:
dat = cbind(Y,X,S)
dat = data.frame(dat)
names(dat) = c("Y","X1","X2","X3","u","v")
cov.names = paste0("X",1:dim(X)[2])
tmp = paste0("Y ~ -1")
formula.GWR = as.formula(paste0(tmp, 
                                paste0('+', cov.names, collapse = '')))
cat("formula = ",as.character(formula.GWR),"\n")
coordinates(dat)=c("u","v")
bw=gwr.sel(formula.GWR, data=dat)
model.fit = gwr(formula.GWR, data=dat,bandwidth=bw, hatmatrix=T)
yhat.GWR = model.fit$SDF$pred

# 10CV-MSPE (GWR):
CV.Result.GWR = cv.GWR(y = Y, X = X, S = S)
mean(CV.Result.GWR)

# save.image("Final_app_result_SVCM_GWR.RData")

