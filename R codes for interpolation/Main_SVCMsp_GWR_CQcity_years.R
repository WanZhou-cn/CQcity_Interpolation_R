##########################################################################
##### 10-fold cross validation in Chongqing city for multiple years
##### 更改点：采用月均Landsat LST（30m）分布，插补每日气温
##########################################################################
rm(list=ls())
gc()
# Library:
# install.packages("devtools",repos="http://cran.r-project.org")
# library(devtools) 
# install_github("funstatpackages/BPST")
# install_github("funstatpackages/Triangulation")
### run in R version 4.0.2

library('MGLM') 
library('BPST')
library('Triangulation')
library('spgwr') 
library('graphics')
library('osqp')
library('tictoc')
library("doParallel")
library("parallel")
# library("foreign")

# 设置基础目录
basic_dir <- "D:/CQcity_Interpolation_R"

# 加载函数文件Source functions
source(file.path(basic_dir,"./3_code/Application/cv.svcmsp_t2.R"))
source(file.path(basic_dir,"./3_code/Application/cv.GWR_t2.R"))

source(file.path(basic_dir, "3_code/Application/fit.svcmsp.R"))
source(file.path(basic_dir, "3_code/Application/svcmsp.sel.R"))

# 设置输入输出路径（基于 basic_dir）
dir_raw    <- file.path(basic_dir, "0_raw")
dir_input  <- file.path(basic_dir, "1_input_Landsat")
dir_output <- file.path(basic_dir, "2_output_Landsat")


v_years = c(2024)         
v_days = c(1:366)
d = 3; r = 1;  # degree d = 2 or 3 // r = 1
Tr_dist = 8    # Tr_dist = 2 / 3.65 / 5 
initial = 123
# v_times = c("TMAX","TMIN")
v_times = c("TMIN")

# 判断是否为闰年函数
is_leap_year <- function(year) {
  (year %% 4 == 0 & year %% 100 != 0) | (year %% 400 == 0)
}

# 每月的天数向量（平年）
days_in_month_normal <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
days_in_month_leap <- c(31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)


##### locations of 0.05 degree grid points
f_S.pop = paste0(dir_input, "/pop_location_d005_CQcity.csv")
S.pop = as.matrix(read.csv(f_S.pop, header=T))
u= unique(S.pop[,1])
v= unique(S.pop[,2])
n1=length(u)
n2=length(v)

# Triangulation:
f_Tr0 = paste0(dir_input, "/Tr_CQcity", ".csv")
f_V0 = paste0(dir_input, "/V_CQcity", ".csv")
if(!file.exists(f_Tr0)|!file.exists(f_V0)) {
  f_bound = paste0(dir_raw, "/boundary_CQcity.csv")
  bound0 = read.csv(f_bound)
  VT = TriMesh(bound0, n = Tr_dist)
  V = VT$V
  Tr = VT$Tr
  write.csv(Tr, f_Tr0, row.names = FALSE, col.names = FALSE)
  write.csv(V, f_V0, row.names = FALSE, col.names = FALSE)
} else {
  Tr=as.matrix(read.csv(f_Tr0))
  V=as.matrix(read.csv(f_V0))
}

##### Basis of 0.05 degree grids
B0.pop=basis(V, Tr, d, r, S.pop)
Q2=B0.pop$Q2
B.pop=B0.pop$B
BQ2.pop=as.matrix(B.pop%*%Q2)
K=B0.pop$K
P=t(Q2)%*%K%*%Q2

for(year0 in v_years) {
  for(time0 in v_times) {
    
    dir_outcv = paste0(dir_output, "/", year0, time0,"_Month")
    if(!dir.exists(dir_outcv)) {dir.create(dir_outcv, recursive = TRUE)}
    
    f_mspe = paste0(dir_outcv, "/3_rmse_", year0, time0, ".csv")
    # if(file.exists(f_mspe)) {next}
    
    ##### read Data for specific year:
    f_st1 = paste0(dir_input, "/station_CQcity_", year0, ".csv")
    f_Ts1 = paste0(dir_input, "/Ts.station.month.", year0,  ".csv")
    f_Ta1 = paste0(dir_input, "/Ta.station.", year0, time0, ".csv")
    
    if(!file.exists(f_st1)|!file.exists(f_Ts1)|!file.exists(f_Ta1)) {
      ## 如果任意一个文件不存在，则执行这里的代码
      ##### data pre-processing
      f_st0 = paste0(dir_raw, "/station_CQcity.csv")
      mat_st0 = as.matrix(read.csv(f_st0))
      colnames(mat_st0)[1] <- "name"
      
      f_Ts0 = paste0(dir_raw, "/Ts.station.month.", year0, ".csv")
      mat_Ts0 = as.matrix(read.csv(f_Ts0))
      colnames(mat_Ts0)[1] <- "name"
      
      f_Ta0 = paste0(dir_raw, "/Ta.station.cn.", year0, time0, ".csv")
      mat_Ta0 = as.matrix(read.csv(f_Ta0))
      colnames(mat_Ta0)[1] <- "name" 
      
      v_st0 = mat_st0[, 1]
      v_sts = mat_Ts0[, 1]
      v_sta = mat_Ta0[, 1]
      
      ##### keep stations with data
      v_stv = intersect(v_st0, v_sts)
      v_stv = intersect(v_stv, v_sta)
      idx_st0 = match(v_stv, v_st0)
      idx_sts = match(v_stv, v_sts)
      idx_sta = match(v_stv, v_sta)
      
      mat_st0 = mat_st0[idx_st0, ]
      mat_Ts0 = mat_Ts0[idx_sts, ]
      mat_Ta0 = mat_Ta0[idx_sta, ]
      
      ##### remove stations with na records
      idx_nona = which(!is.na(mat_Ts0[, 2]))
      mat_st1 = mat_st0[idx_nona, ]
      mat_Ts1 = mat_Ts0[idx_nona, ]
      mat_Ta1 = mat_Ta0[idx_nona, ]
      
      ##### remove outliers
      idx_nona1 = c(1:length(mat_st1[, 1]))
      for(i0 in c(1:length(mat_st1[, 1]))) {
        v_value = c(mat_Ta1[i0, 2:366], mat_Ts1[i0, 2:366])
        v_00 = abs(v_value)
        if(max(v_00, na.rm=TRUE) > 100) {
          idx_nona1 = setdiff(idx_nona1, i0)
        }
      }
      mat_st1 = mat_st1[idx_nona1, ]
      mat_Ts1 = mat_Ts1[idx_nona1, ]
      mat_Ta1 = mat_Ta1[idx_nona1, ]
      
      write.csv(mat_st1, f_st1, row.names = FALSE)
      write.csv(mat_Ts1, f_Ts1, row.names = FALSE)
      write.csv(mat_Ta1, f_Ta1, row.names = FALSE)
    } else {
      mat_st1 = read.csv(f_st1)
      mat_Ts1 = read.csv(f_Ts1)
      mat_Ta1 = read.csv(f_Ta1)
      colnames(mat_st1)[1] <- "name"
      colnames(mat_Ts1)[1] <- "name"
      colnames(mat_Ta1)[1] <- "name"
    }
    
    ##### calculate basic parameters
    data1 = as.data.frame(mat_st1)
    S = cbind(data1$x, data1$y)
    ##### Sample basis:
    B0 = basis(V, Tr, d, r, S)
    B = B0$B
    
    # 定义天数序列
    leap = is_leap_year(year0)                # 判断闰年/平年
    days_in_month = if (leap) days_in_month_leap else days_in_month_normal
    total_days = sum(days_in_month)
    v_days = 1:total_days
    
    ##### estimation for each day
    # registerDoParallel(16)
    # foreach(i0 = c(1:length(v_days))) %dopar% {
    for(i0 in c(1:length(v_days))) {
      
      day0 = v_days[i0]
      
      # 计算当前 day0 属于哪一个月
      cum_days = cumsum(days_in_month)
      month_index = which(cum_days >= day0)[1]   # 当前天属于的月份索引（1~12）
      month_str <- sprintf("%02d", month_index)  # 格式化为两位数
      
      # 打印处理信息
      print(paste0("year", year0, " ", "-> month", month_index,"-> Type: ", time0, "-> Day: ", day0))
      
      f_gwr.cv = paste0(dir_outcv, "/2_gwr_cv_", year0, "_DOY", day0, time0, ".csv")
      if(!file.exists(f_gwr.cv)) {
        
        # 构造数据
        data = data.frame(
          name = data1$name,
          x=data1$x, 
          y=data1$y,
          lat=data1$y,
          plst = mat_Ts1[[month_index + 1]],  # 月LST列（跳过name列）
          dem=data1$dem, 
          at = mat_Ta1[[day0 + 1]]  # 日气温数据（跳过name列）
        )
        
        Y = data$at
        X = cbind(data$plst, data$dem)
        X = scale(X, scale=TRUE) # standardized covariates
        X = cbind(1, X)          # add intercept
        
        ###############################################################################################
        # SVCMsp: intercept, LST, DEM
        lambda.svcm = svcmsp.sel(B, Q2, K, X, Y, sign.coeff = c(0, 1, -1))
        Result.svcm = fit.svcmsp(B, Q2, K, lambda.svcm, X, Y, sign.coeff = c(0, 1, -1))
        yhat.svcm = rowSums(X*Result.svcm$beta)
        beta.hat.pop.svcm = B.pop%*%Result.svcm$gamma
        mhat_all.svcm = beta.hat.pop.svcm
        
        ##### output beta
        f_svcm = paste0(dir_outcv, "/1_svcmsp_beta_", year0, "_DOY", day0, time0, ".csv")
        mat_beta = as.matrix(Result.svcm$beta)
        write.csv(mat_beta, f_svcm, row.names = FALSE)
        
        # 10CV-MSPE (SVCMsp):
        set.seed(initial)
        Test0 = sample(1:length(Y))
        CV.Result.svcm = cv.svcmsp(y = Y, X = X, B = B, Q2 = Q2, K = K, iTest=Test0,
                                   lambda = Result.svcm$lambdac, sign.coeff = c(0, 1, -1), nfold = 10)
        ##### The mean squared prediction error (MSPE) based on k-fold cross-validation
        mspe0 = mean(sqrt(CV.Result.svcm[[1]]))
        
        ##### output cv results
        f_svcm.cv = paste0(dir_outcv, "/1_svcmsp_cv_", year0, "_DOY", day0, time0, ".csv")
        mat_cv0 = CV.Result.svcm[[2]]
        colnames(mat_cv0) = c("y_true", "y_pred", "bias","cvID")
        write.csv(mat_cv0, f_svcm.cv, row.names = FALSE)
        
        ################################################################################################
        ##### GWR:
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
        
        ##### output SDF
        f_gwr_sdf = paste0(dir_outcv, "/2_gwr_sdf_", year0, "_DOY", day0, time0, ".csv")
        mat_sdf = as.data.frame(model.fit$SDF)
        write.csv(mat_sdf, f_gwr_sdf, row.names = FALSE)
        
        # 10CV-MSPE (GWR):
        CV.Result.GWR = cv.GWR(y = Y, X = X, S = S, iTest=Test0, nfold = 10)
        mspe1 = mean(sqrt(CV.Result.GWR[[1]]))
        
        ##### output cv results
        mat_cv1 = CV.Result.GWR[[2]]
        colnames(mat_cv1) = c("y_true", "y_pred", "bias", "cvID")
        write.csv(mat_cv1, f_gwr.cv, row.names = FALSE)
      }
      
    }
    
    
  }
  
}

