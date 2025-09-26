rm(list=ls())

# Library:
  # install.packages("devtools",repos="http://cran.r-project.org")
  # library(devtools)
  # install_github("funstatpackages/BPST")
  # install_github("funstatpackages/Triangulation")

# library('BPST')
library('Triangulation')
library('mgcv')
library('MGLM') 
library('osqp')


# Source Functions:
source('fit.svcmsp.R')
source("CIboot.svcmsp.R")
source("svcmsp.sel.R")
source("cv.svcmsp.R")
source('data.svcmsp.R') 

source("cv.svcm.R")
source("svcm.sel.R")
source("plsfitGCV_full.R")

source("B0_Generator.R")
source("basis.R")
source("BPST.est.ho.R")
source("BPST.est.pc.R")
source("cv.BPST.R")
source("data.BPST.R")
source("findbdt.R")
source("fit.BPST.R")
source("inVT.R")
source("plagCV.R")
source("plot.BPST.R")
source("predict.BPST.R")
source("qrH.R")
source("RcppExports.R")
source("seval.R")
source("smoothness.R")
source("TArea.R")
source("tdata.R")
source("Tr1.R")
source("Tr2.R")
source("V1.R")
source("V2.R")


# Simulation Setup:
d = 2; r = 1;  # degree d = 2 or 3 // r = 1
n = 1000       # sample size n = 1000 or 2000
Tr_dist = 2    # Tr_dist = 2 / 3.65 / 5 
               # (They correspond to Tr1, Tr2, and Tr3 in the manuscript.) 

run.cv  = 0    # 1: run 10-cv MSPE code 
               # 0: skip 10-cv MSPE results 

run.bootstrap = 0 # 1: run bootstrap procedure for confidence interval 
                  # 0: skip bootstrap results

# Number of Simulation:
nsim = 500

# Triangulation:
data("horseshoe")
VT = TriMesh(horseshoe, n = Tr_dist)
V = VT$V
Tr = VT$Tr

# Population:
pop.all = data.svcmsp()
N.all=nrow(pop.all)
ind1=inVT(V,Tr,pop.all[,6],pop.all[,7])
ind1=ind1$ind.inside
ind2=(1:nrow(pop.all))[!is.na(pop.all[,1])]
ind=sort(intersect(ind1,ind2))
pop.r=pop.all[ind,]
Y.pop=pop.r[,1]
beta.pop=pop.r[,2:3]
X.pop=pop.r[,4:5]
S.pop=round(pop.r[,6:7],2)
Npop=nrow(pop.r)
np=ncol(X.pop)

u=unique(round(pop.all[,6],2))
v=unique(round(pop.all[,7],2))
n1=length(u)
n2=length(v)

# Population Basis:
B0.pop = basis(V,Tr,d,r,S.pop)
Q2 = B0.pop$Q2
B.pop = B0.pop$B
BQ2.pop = as.matrix(B.pop%*%Q2)
K = B0.pop$K
P = t(Q2)%*%K%*%Q2



mise.beta.all = c(); neg_number.all =c()
mise.beta.all2 = c(); neg_number.all2 =c()


target_Beta1.all = c(); target_Beta2.all =c()

CV.Result.all=c()
CV.Result.all2=c()


#main <- function(iter){
for(iter in 1:nsim){
  cat("\nIteration:",iter,"\n")
  set.seed(iter) 
  
  # Sample:
  ind.s = sort(sample(Npop,n,replace=FALSE))
  dat = as.matrix(pop.r[ind.s,])
  Y = dat[,1]
  beta0 = dat[,2:3]
  X = as.matrix(dat[,4:5])
  S = dat[,6:7]
  B = B.pop[ind.s,]
  
  # Fit the SVCM-SP: 
  lambda = svcmsp.sel(B, Q2, K, X, Y, sign.coeff = c(1,1))
  Result = fit.svcmsp(B, Q2, K, lambda, X, Y, sign.coeff = c(1,1)) 

  beta.hat.pop = B.pop%*%Result$gamma
  mise.beta = colMeans((beta.pop - beta.hat.pop)^2)
  mise.beta.all = rbind(mise.beta.all,mise.beta)
  neg_number.all = rbind(neg_number.all, round(colMeans(cbind(beta.hat.pop[,1] < 0, beta.hat.pop[,2] < 0)), 2))
  
  # # Fit the SVCM: 
  lambda2 = svcm.sel(B, Q2, K, X, Y) 
  Result2 = plsfitGCV_full(B, Q2, K, lambda2, X, Y) 
  
  beta.hat.pop2 = B.pop%*%Result2$gamma # SVCM
  mise.beta2 = colMeans((beta.pop - beta.hat.pop2)^2) 
  mise.beta.all2 = rbind(mise.beta.all2, mise.beta2) 
  neg_number.all2 = rbind(neg_number.all2, round(colMeans(cbind(beta.hat.pop2[,1] < 0, beta.hat.pop2[,2] < 0)), 2)) 
  

  
  # Confidence Intervals of SVCM-SP:
  if (run.bootstrap == 1){
    result_CB = CIboot.svcmsp(mfit = Result, y=Y, X, B, Q2, B0.pop = B0.pop$B, K = K,
                              lambda = Result$lambda, nB = 100, sign.coeff = c(1,1))
    LocTest = matrix(c(0.5, 0.5, 1.5, 1.5, 2.5, 2.5,
                       0.5,-0.5,0.5,-0.5,0.5,-0.5), ncol=2, byrow=F)
    Ind_Test = match(data.frame(t(S.pop)),data.frame(t(LocTest)))
    Ind_Test = which(is.na(Ind_Test)==0)
    Ind_Test = Ind_Test[c(4,1,5,2,6,3)] # locations S1 to S6.

    target_Beta1.all = cbind(target_Beta1.all, (result_CB$Lower.CB[Ind_Test,1]<beta.pop[Ind_Test,1])&(result_CB$Upper.CB[Ind_Test,1]>beta.pop[Ind_Test,1]))
    target_Beta2.all = cbind(target_Beta2.all, (result_CB$Lower.CB[Ind_Test,2]<beta.pop[Ind_Test,2])&(result_CB$Upper.CB[Ind_Test,2]>beta.pop[Ind_Test,2]))
  }
  
  # 10-fold MSPE:
  if (run.cv == 1){
    CV.Result = cv.svcmsp(y = Y, X = X, B = B, Q2 = Q2, K = K, lambda = Result$lambda, sign.coeff = c(1, 1))
    CV.Result.all = c(CV.Result.all, mean(CV.Result))
    
    CV.Result2 = cv.svcm(y = Y, X = X, S = S, B = B, Q2 = Q2, K = K, lambda = Result2$lambda) 
    CV.Result.all2 = c(CV.Result.all2, mean(CV.Result2)) 
  }
  
}

# AMISE and Proposion of Negative Values (SVCM-SP):
Result.all = round(c(colMeans(mise.beta.all),colMeans(neg_number.all)), 3)
names(Result.all) = c("amise.beta0","amise.beta1","Prop.Neg.beta0","Prop.Neg.beta1")
Result.all

# AMISE and Proposion of Negative Values (SVCM):
Result.all2 = round(c(colMeans(mise.beta.all2),colMeans(neg_number.all2)), 3) # SVCM
names(Result.all2) = c("amise.beta0","amise.beta1","Prop.Neg.beta0","Prop.Neg.beta1") 
Result.all2 # SVCM

# Coverage of 95% Confidence Intervals:
if (run.bootstrap == 1){
  CI.result.name = c("beta0", "beta1")
  CI.result =
    data.frame(rbind(round(rowMeans(target_Beta1.all), 3),    
                     round(rowMeans(target_Beta2.all), 3)))   
  CI.result = cbind(CI.result.name, CI.result)
  names(CI.result) = c("beta", "S1","S2","S3","S4","S5","S6")
  CI.result
}


# 10-CV MSPE:
if (run.cv == 1){
  cat("10-CV SVCM-SP: ", round(mean(CV.Result.all), 3),  "  ", 
      "10-CV SVCM: ", round(mean(CV.Result.all2), 3))
}


##################################################################
####################### Only for the Plot ########################
##################################################################
library(plot3D)
########################## True ##########################
beta.rr1 =  beta.pop[,1]
beta.rr2 =  beta.pop[,2]

beta.r1.all=matrix(NA,n1*n2,1)
beta.r1.all[ind,]= as.matrix(beta.rr1)
beta.r1.mtx=matrix(beta.r1.all,n1,n2)

beta.r2.all=matrix(NA,n1*n2,1)
beta.r2.all[ind,]= as.matrix(beta.rr2)
beta.r2.mtx=matrix(beta.r2.all,n1,n2)

image2D(x=u,y=v, z=beta.r1.mtx)
image2D(x=u,y=v, z=beta.r2.mtx)

########################## SVCM-SP ##########################
beta.hat.pop=B.pop%*%Result$gamma
beta.r1 = beta.hat.pop[,1]
beta.r2 = beta.hat.pop[,2]
mhat_all=cbind(beta.r1,beta.r2)

beta.rr1 = beta.r1
beta.rr2 = beta.r2

beta.r1.all=matrix(NA,n1*n2,1)
beta.r1.all[ind,]= as.matrix(beta.rr1)
beta.r1.mtx=matrix(beta.r1.all,n1,n2)

beta.r2.all=matrix(NA,n1*n2,1)
beta.r2.all[ind,]= as.matrix(beta.rr2)
beta.r2.mtx=matrix(beta.r2.all,n1,n2)

image2D(x=u,y=v, z=beta.r1.mtx)
image2D(x=u,y=v, z=beta.r2.mtx)
