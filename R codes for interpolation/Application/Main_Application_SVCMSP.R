rm(list=ls())

# Library:
  # install.packages("devtools",repos="http://cran.r-project.org")
  # library(devtools)
  # install_github("funstatpackages/BPST")
  # install_github("funstatpackages/Triangulation")


library('MGLM') 
library('BPST')
library('Triangulation')
library('osqp')

# Source functions
source('fit.svcmsp.R')
source("CIboot.svcmsp.R")
source("cv.svcmsp.R")
source("svcmsp.sel.R")

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

# Population basis:
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

# Fitting the model:
lambda = svcmsp.sel(B, Q2, K, X, Y, sign.coeff = c(0, 1, -1))
Result = fit.svcmsp(B, Q2, K, lambda, X, Y, sign.coeff = c(0, 1, -1)) 
beta.hat.pop = B.pop%*%Result$gamma
mhat_all = beta.hat.pop


# # # 10-fold MSPE of SVCM-SP:
# CV.Result = cv.svcmsp(y = Y, X = X, B = B, Q2 = Q2, K = K, lambda = Result$lambdac,
#                       sign.coeff = c(0, 1, -1))
# CV.Result
# mean(CV.Result)
# 
# 
# 
# # # Confidence Intervals of SVCM-SP:
# result_CB = CIboot.svcmsp(mfit = Result, y = Y, X = X, B = B, Q2 = Q2, B0.pop = B0.pop$B,
#                           K = K, lambda = Result$lambdac, sign.coeff = c(0, 1, -1))
# 


# source("test.svcmsp.R")
# # Global test:
# Test_global_result = test.svcmsp(mfit=Result, S=S,
#                                  B=B, Q2=Q2, K=K, lambda=Result$lambdac, X, Y,
#                                  sign.coeff = c(0,1,-1), nB = 100, NSVCM = TRUE, H.accuracy = FALSE,
#                                  test = "global")
# Test_global_result

# # Individual tests: ind.test = 1,2,3 for three individual test for bivariate functions
# Test_individual_result = test.svcmsp(mfit=Result, S=S, B=B, Q2=Q2, K=K, lambda=Result$lambdac,
#                                      X, Y, sign.coeff = c(0,1,-1), nB = 100, NSVCM = TRUE, H.accuracy = FALSE,
#                                      test = "individual", ind.test = 1)
# Test_individual_result
