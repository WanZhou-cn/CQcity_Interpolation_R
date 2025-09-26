##### calculate parameters in Chongqing city
rm(list=ls())
gc()
# Library:
# install.packages("devtools",repos="http://cran.r-project.org")
# library(devtools)
# install_github("funstatpackages/BPST")
# install_github("funstatpackages/Triangulation")

library('MGLM')
library('BPST')
library('Triangulation')
library('spgwr')
library('robustHD')
library('raster')
library('rgdal')
library('tictoc')
library('osqp')
# library("doParallel")
# library("parallel")
# library('envirem')
# library(SpaDES)

changeCoord <- function(v_lon, v_lat, proj_in, proj_out) {
  dat0 = data.frame(lon=v_lon, lat=v_lat)
  coordinates(dat0) = c("lon", "lat")
  proj4string(dat0) = proj_in
  dat1 = spTransform(dat0, proj_out)
  return(dat1)
}

# 设置基础目录
basic_dir <- "D:/CQcity_Interpolation_R"

# 加载函数文件Source functions
source(file.path(basic_dir, "3_code/Application/fit.svcmsp.R"))
source(file.path(basic_dir, "3_code/Application/svcmsp.sel.R"))

# 设置输入输出路径（基于 basic_dir）
dir_raw    <- file.path(basic_dir, "0_raw")
dir_input  <- file.path(basic_dir, "1_input")
dir_output <- file.path(basic_dir, "2_output")


v_years = c(2019)
# v_times = c("TMAX","TMIN")
v_times = c("TMAX")
# v_times = c("TMIN")
v_days = c(1:365)
d = 3; r = 1;  # degree d = 2 or 3 // r = 1
Tr_dist = 8    # Tr_dist = 2 / 3.65 / 5 
initial = 123
Ds=c(paste0('00',1:9),paste0('0',10:99),100:365)

##### read shapefile of stations
f_stations = paste0(dir_raw, "/pst_stations.shp")
shp_stations = shapefile(f_stations)
proj1 = shp_stations@proj4string

##### locations of grid points
f_S.pop = paste0(dir_input, "/pop_location_d005_CQcity.csv")
S.pop = as.matrix(read.csv(f_S.pop, header=T))
u= unique(S.pop[,1])
v= unique(S.pop[,2])
n1=length(u)
n2=length(v)

# Triangulation:
f_Tr0 = paste0(dir_input, "/Tr_CQcity.csv")
f_V0 = paste0(dir_input, "/V_CQcity.csv")
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

# 创建DEM 参数文件输出目录
dir_outdem = paste0(dir_input, "/dem_param")
if(!dir.exists(dir_outdem)){dir.create(dir_outdem, recursive = T)}

##### read dem data
f_dem1 = paste0(dir_input, "/dem_CQcity1.tif")
if(!file.exists(f_dem1)) {
  f_dem0 = paste0(dir_input, "/dem_CQcity0.tif")
  ras_dem0 = raster(f_dem0)
  ras_dem1 = raster::trim(ras_dem0, padding=0, values=NA)
  writeRaster(ras_dem1, f_dem1, overwrite=T, dataType="INT2S")
} else {
  ras_dem1 = raster(f_dem1)
}

# 分块处理 DEM 数据
len_1km = length(ras_dem1[])
v_start = seq(1, len_1km, by=5000)
len_tile = length(v_start)
v_start1 = v_start[2:len_tile] - 1
v_end = c(v_start1, len_1km)

# 获取 DEM 的投影信息
proj2 = ras_dem1@crs

#### get coordinates
# 获取 DEM 栅格的经纬度坐标
log=setValues(ras_dem1,xFromCell(ras_dem1,1:ncell(ras_dem1)))
lat=setValues(ras_dem1,yFromCell(ras_dem1,1:ncell(ras_dem1)))

#### calculate dem related parameters
# 分块处理 DEM 数据
# registerDoParallel(16)
# foreach(i = c(1:len_tile)) %dopar% {
# # for(i in c(1:len_tile)) {
#   f_bpred0 = paste0(dir_outdem, "/bpred0_", i, ".rds")
#   f_predgrid0 = paste0(dir_outdem, "/predgrid0_", i, ".rds")
#   if(!file.exists(f_bpred0)) {
#     ista = v_start[i]
#     iend = v_end[i]
#     ##### change projection of the coordinates
#     dat_coord = changeCoord(log@data@values[ista:iend], lat@data@values[ista:iend], proj2, proj1)
#     predgrid0=data.frame(id=c(ista:iend),xx=dat_coord@coords[, 1],yy=dat_coord@coords[, 2],dem=ras_dem1[ista:iend])
#     S.pred0 = cbind(predgrid0$xx, predgrid0$yy)
#     B0.pred0=basis(V, Tr, d, r, S.pred0)
#     B.pred0=B0.pred0$B
#     
#     saveRDS(predgrid0, f_predgrid0)
#     saveRDS(B.pred0, f_bpred0)
#   }
#   
# }


# 检查并创建目录，用于存储基函数参数文件（basis_param 文件夹）
dir_basis = paste0(dir_input, "/basis_param")
if(!dir.exists(dir_basis)){dir.create(dir_basis, recursive = T)}

for(time0 in v_times) {
  
  # 设置变量名称标识（日间/夜间）
  if(time0=="TMAX") {
    s_t0 = "Day"
  } else if(time0=="TMIN"){
    s_t0 = "Nit"
  }
  
  # 按年份处理数据:针对每一年份的数据进行循环处理，并为每个年份创建输出目录
  for(year0 in v_years) {
    
    # 输出结果系数图像的目录
    dir_outCoef = paste0(dir_output, "/Coef", year0)
    if(!dir.exists(dir_outCoef)){dir.create(dir_outCoef, recursive = T)}
    
    ##### read Data for specific year:
    f_st1 = paste0(dir_input, "/station_CQcity_", year0, ".csv")
    f_Ts1 = paste0(dir_input, "/Ts.station.cn.", year0, time0, ".csv")
    f_Ta1 = paste0(dir_input, "/Ta.station.cn.", year0, time0, ".csv")
    
    if(!file.exists(f_st1)|!file.exists(f_Ts1)|!file.exists(f_Ta1)) {
      ##### data pre-processing
      # 读取原始站点信息与气象数据
      f_st0 = paste0(dir_raw, "/stations_CQcity.csv")
      mat_st0 = as.matrix(read.csv(f_st0))
      colnames(mat_st0)[1] <- "name"
      
      f_Ts0 = paste0(dir_raw, "/Ts.station.cn.", year0, time0, ".csv")
      mat_Ts0 = as.matrix(read.csv(f_Ts0))
      colnames(mat_Ts0)[1] <- "name"
      
      f_Ta0 = paste0(dir_raw, "/Ta.station.cn.", year0, time0, ".csv")
      mat_Ta0 = as.matrix(read.csv(f_Ta0))
      colnames(mat_Ta0)[1] <- "name"
      
      # v_为站点名
      v_st0 = mat_st0[, 1]
      v_sts = mat_Ts0[, 1]
      v_sta = mat_Ta0[, 1]
      
      ##### keep stations with data
      v_stv = intersect(v_st0, v_sts)
      v_stv = intersect(v_stv, v_sta)
      # v_stv为交叉站点（三类数据中都存在的有效站点）
      idx_st0 = match(v_stv, v_st0)
      idx_sts = match(v_stv, v_sts)
      idx_sta = match(v_stv, v_sta)
      
      mat_st0 = mat_st0[idx_st0, ]
      mat_Ts0 = mat_Ts0[idx_sts, ]
      mat_Ta0 = mat_Ta0[idx_sta, ]
      
      ##### remove stations with na records
      # 删除缺失值站点
      idx_nona = which(!is.na(mat_Ts0[, 2]))
      mat_st1 = mat_st0[idx_nona, ]
      mat_Ts1 = mat_Ts0[idx_nona, ]
      mat_Ta1 = mat_Ta0[idx_nona, ]
      
      ##### remove outliers
      # 去除异常值（大于100℃）
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
    
    print("Cal. basis for samples")
    ##### calculate basic parameters
    data1 = as.data.frame(mat_st1)
    S = cbind(data1$x, data1$y)
    ##### Basis
    f_B0pop = paste0(dir_basis,"/B0pop_", year0, time0,".rds")
    if(!file.exists(f_B0pop)) {
      B0.pop=basis(V, Tr, d, r, S.pop)
      saveRDS(B0.pop, f_B0pop)
    } else {
      B0.pop = readRDS(f_B0pop)
    }
    Q2=B0.pop$Q2
    B.pop=B0.pop$B
    BQ2.pop=as.matrix(B.pop%*%Q2)
    K=B0.pop$K
    P=t(Q2)%*%K%*%Q2
    
    ##### Sample basis:
    f_B0 = paste0(dir_basis, "/B0_", year0, time0, ".rds")
    if(!file.exists(f_B0)) {
      B0 = basis(V, Tr, d, r, S)
      saveRDS(B0, f_B0)
    } else {
      B0 = readRDS(f_B0)
    }
    B = B0$B
    
    ##### estimation for each day
    # registerDoParallel(16)
    # foreach(i0 = c(161:240)) %dopar% {
    # 内层处理 具体的某一天（第200天、第194天）
    for(i0 in c(5,6)) {
      
      day0 = v_days[i0]
      f_coefDEM_svcm = paste0(dir_outCoef, "/CQcity_coefDEM_SVCMsp_",year0,time0,"_",Ds[day0],".tif")
      if(file.exists(f_coefDEM_svcm)){next}
      f_coefLST_svcm = paste0(dir_outCoef, "/CQcity_coefLST_SVCMsp_",year0,time0,"_",Ds[day0],".tif")
      if(file.exists(f_coefLST_svcm)){next}
      
      # 开始计时
      tic(i0)
      data = data.frame(name=data1$name, plst=mat_Ts1[, day0+1], dem=data1$dem, at=mat_Ta1[, day0+1])
      Y = data$at
      X = cbind(data$plst, data$dem)
      mean1 = mean(X[, 1])
      mean2 = mean(X[, 2])
      sd1 = sd(X[, 1])
      sd2 = sd(X[, 2])
      X = scale(X, scale=TRUE) # standardized covariates
      X = cbind(1, X)          # add intercept
      
      f_result.svcm = paste0(dir_basis, "/svcm_param_", year0, time0, Ds[day0], ".rds")
      if(!file.exists(f_result.svcm)) {
        lambda.svcm = svcmsp.sel(B, Q2, K, X, Y, sign.coeff = c(0, 1, -1))
        Result.svcm = fit.svcmsp(B, Q2, K, lambda.svcm, X, Y, sign.coeff = c(0, 1, -1)) 
        saveRDS(Result.svcm, f_result.svcm)
      } else {
        Result.svcm = readRDS(f_result.svcm)
      }
      
      ##### read LST
      # dir_lst = paste0("../GlobalGPFILL_statis/2_output/Final/GF_LST/",year0,"mosaicbf")
      # dir_lst = paste0("E:/", year0, s_t0,"_ChongqingCity")
      dir_lst = paste0("E:/",year0,"mosaicbf")
      dir_lst = paste0("D:/MODIS LST/", year0, s_t0,"_ChongqingCity")
      f_lst = paste0(dir_lst, "/CQcity_",s_t0,year0,"_",Ds[day0],".tif")
      ras_lst = raster(f_lst)
      
      ##### crop LST
      ras_lst = raster::crop(ras_lst, ras_dem1) / 100
      
      ras_coefDEM = ras_coefLST = ras_lst
      ras_coefDEM[] = ras_coefLST[] = NA
      
      ##### prediction Ta tile by tile
      for(i in c(1:len_tile)) {
        print(i)
        ista = v_start[i]
        iend = v_end[i]
        f_bpred0 = paste0(dir_outdem, "/bpred0_", i, ".rds")
        f_predgrid0 = paste0(dir_outdem, "/predgrid0_", i, ".rds")
        f_b0pred0 = paste0(dir_outdem, "/b0pred0_", i, ".rds")
        if(!file.exists(f_b0pred0)) {
          ##### change projection of the coordinates
          dat_coord = changeCoord(log@data@values[ista:iend], lat@data@values[ista:iend], proj2, proj1)
          predgrid0=data.frame(id=c(ista:iend),xx=dat_coord@coords[, 1],yy=dat_coord@coords[, 2],dem=ras_dem1[ista:iend])
          S.pred0 = cbind(predgrid0$xx, predgrid0$yy)
          B0.pred0=basis(V, Tr, d, r, S.pred0)
          B.pred0=B0.pred0$B
          
          saveRDS(predgrid0, f_predgrid0)
          saveRDS(B0.pred0, f_b0pred0)
          # saveRDS(B.pred0, f_bpred0)
        } else {
          predgrid0 = readRDS(f_predgrid0)
          # B.pred0 = readRDS(f_bpred0)
          B0.pred0 = readRDS(f_b0pred0)
          B.pred0 = B0.pred0$B
        }
        idx_in0 = B0.pred0$Ind.inside
        predgrid0$plst = ras_lst[ista:iend]
        predgrid0 = predgrid0[idx_in0, ]
        idx_comp = complete.cases(predgrid0)
        idx_sample = which(idx_comp==TRUE)
        predgrid0 = predgrid0[idx_sample, ]
        
        ##### create X and S
        X.pred = cbind(predgrid0$plst, predgrid0$dem)
        ##### standardlize and add intercept
        X.pred = scale(X.pred, center=c(mean1, mean2), scale=c(sd1, sd2))
        X.pred = cbind(1, X.pred)
        
        beta.hat.pred = B.pred0[idx_sample, ]%*%Result.svcm$gamma
        
        # yhat.pred = rowSums(X.pred*beta.hat.pred)
        mat_beta = as.matrix(beta.hat.pred)
        ras_coefLST[predgrid0$id] = mat_beta[, 2]
        ras_coefDEM[predgrid0$id] = mat_beta[, 3]
      }
      
      writeRaster(ras_coefDEM, f_coefDEM_svcm, overwrite=TRUE)
      writeRaster(ras_coefLST, f_coefLST_svcm, overwrite=TRUE)
      toc()
    }
  }
}



