# -*- coding: utf-8 -*-
#### 案例：重庆市中心城区气象要素插值-日最高温
## 更改点：采用月均Landsat LST（30m）分布，插补每日气温

rm(list=ls())
gc()
library('MGLM')
library('BPST')
library('Triangulation')
library('spgwr')
library('robustHD')
library('raster')
library('rgdal')
library('tictoc')
library('osqp')
library("doParallel")
library("parallel")


changeCoord <- function(v_lon, v_lat, proj_in, proj_out) {
  dat0 = data.frame(lon=v_lon, lat=v_lat)
  coordinates(dat0) = c("lon", "lat")
  proj4string(dat0) = proj_in
  dat1 = spTransform(dat0, proj_out)
  return(dat1)
}


# 设置基础目录
basic_dir <- "/data/wanzhou/CQcity_Interpolation_R/"

# 加载函数文件Source functions
source(file.path(basic_dir, "3_code/Application/fit.svcmsp.R"))
source(file.path(basic_dir, "3_code/Application/svcmsp.sel.R"))

# 设置输入输出路径（基于 basic_dir）
dir_raw    <- file.path(basic_dir, "0_raw")
dir_input  <- file.path(basic_dir, "1_input")
dir_output <- file.path(basic_dir, "2_output")


v_years = c(2019)            # 定义年份
# v_times = c("TMAX","TMIN")
# v_times = c("RHU")           # 定义时间变量湿度RHU
v_times = c("TMAX")
# v_times = c("TMIN")
v_days = c(1:365)           # 定义年内天数范围

# 判断是否为闰年函数
is_leap_year <- function(year) {
  (year %% 4 == 0 & year %% 100 != 0) | (year %% 400 == 0)
}

# 每月的天数向量（平年）
days_in_month_normal <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
days_in_month_leap <- c(31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

# 样条拟合的参数
d = 3; r = 1;     # 样条参数 degree d = 2 or 3 // r = 1
Tr_dist = 8       # 三角剖分距离 Tr_dist = 2 / 3.65 / 5 
initial = 123     # 随机种子，确保结果的可重复性 
Ds=c(paste0('00',1:9),paste0('0',10:99),100:365)  # 定义日期向量


##### read shapefile of stations
f_stations = paste0(dir_raw, "/pst_stations.shp")
# f_stations = paste0(dir_raw, "/Station_CQcity.shp")
shp_stations = shapefile(f_stations)
proj1 = shp_stations@proj4string    # 获取 shapefile 的投影信息（即坐标参考系，CRS）


##### locations of grid points(仅用于可视化模拟结果？)
f_S.pop = paste0(dir_input, "/pop_location_d005_CQcity.csv")
S.pop = as.matrix(read.csv(f_S.pop, header=T))
u = unique(S.pop[,1])
v = unique(S.pop[,2])
n1 =length(u)
n2 =length(v)


# Triangulation:(生成或加载三角剖分（Triangulation）结果)
f_Tr0 = paste0(dir_input, "/Tr_CQcity.csv")
f_V0 = paste0(dir_input, "/V_CQcity.csv")
if(!file.exists(f_Tr0)|!file.exists(f_V0)) {
  # 如果文件不存在，根据重庆市边界生成三角剖分
  f_bound = paste0(dir_raw, "/boundary_CQcity.csv")
  bound0 = read.csv(f_bound)
  VT = TriMesh(bound0, n = Tr_dist)
  V = VT$V
  Tr = VT$Tr
  write.csv(Tr, f_Tr0, row.names = FALSE, col.names = FALSE)
  write.csv(V, f_V0, row.names = FALSE, col.names = FALSE)
} else {
  # 将三角剖分结果中的顶点信息存储到 V，三角形信息存储到 Tr
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
  writeRaster(ras_dem1, f_dem1, overwrite=T, datatype="INT2S")
} else {
  ras_dem1 = raster(f_dem1)
}

# 分块处理 DEM 数据
len_30m = length(ras_dem1[])                   # 提取栅格数量
v_start = seq(1, len_30m, by=20000)
len_tile = length(v_start)
v_start1 = v_start[2:len_tile] - 1
v_end = c(v_start1, len_30m)

# 获取 DEM 的投影信息（单位：米）
proj2 = ras_dem1@crs

#### get coordinates
# 获取 DEM 栅格的经纬度坐标
log = setValues(ras_dem1,xFromCell(ras_dem1,1:ncell(ras_dem1)))
lat = setValues(ras_dem1,yFromCell(ras_dem1,1:ncell(ras_dem1)))

#### calculate dem related parameters
# 分块处理 DEM 数据(更新参数后需删除相关文件)
# B0.pred0（预测点上的基函数矩阵）只在 for 循环外部计算一次
for(i in c(1:len_tile)) {
  f_b0pred0 = paste0(dir_outdem, "/b0pred0_", i, ".rds")
  f_predgrid0 = paste0(dir_outdem, "/predgrid0_", i, ".rds")
  
  if(!file.exists(f_b0pred0)) {
    # 获取当前分块的索引范围
    ista = v_start[i]
    iend = v_end[i]
    
    ##### change projection of the coordinates
    # 坐标投影转换：MODIS → WGS84
    # 代码报错：经纬度坐标超出有效范围（修改函数）
    dat_coord = changeCoord(log@data@values[ista:iend], lat@data@values[ista:iend], proj2, proj1)
    
    # 构建预测网格
    predgrid0 = data.frame(
      id  = c(ista:iend),
      xx  = dat_coord@coords[, 1],
      yy  = dat_coord@coords[, 2],
      dem = ras_dem1[ista:iend]
    )
    
    ### ???? 添加调试代码：查看 predgrid0 的结构和 dem 值是否正常
    cat("---- Tile", i, "----\n")
    print(head(predgrid0))                # 查看前几行
    cat("DEM summary:\n")
    print(summary(predgrid0$dem))        # 查看 dem 列的统计信息
    cat("Number of NA in DEM:", sum(is.na(predgrid0$dem)), "\n")
    cat("Range of X (xx):", range(predgrid0$xx), "\n")
    cat("Range of Y (yy):", range(predgrid0$yy), "\n")
    cat("\n")
    
    # 基于三角剖分计算基函数参数，用于 DEM 数据的插值或拟合
    S.pred0 = cbind(predgrid0$xx, predgrid0$yy)
    B0.pred0=basis(V, Tr, d, r, S.pred0)
    B.pred0=B0.pred0$B
    
    # 保存结果：将预测网格和基函数结构保存为 .rds 文件，便于后续加载使用，无需重复计算
    saveRDS(predgrid0, f_predgrid0)
    saveRDS(B0.pred0, f_b0pred0)
  }
}


# 检查并创建目录，用于存储基函数参数文件（basis_param 文件夹）
dir_basis = paste0(dir_input, "/basis_param")
if(!dir.exists(dir_basis)){
  dir.create(dir_basis, recursive = T)
  message("目录 'basis_param' 已成功创建！")
} else {
  message("目录 'basis_param' 已存在，无需创建。")
}


# times(遍历处理变量)
for(time0 in v_times) {
  
  # 设置变量名称标识（日间/夜间）
  if(time0=="TMAX") {
    s_t0 = "Day"
  } else if(time0=="TMIN"){
    s_t0 = "Nit"
  }
  
  # 按年份处理数据:针对每一年份的数据进行循环处理，并为每个年份创建输出目录
  for(year0 in v_years) {
    
    # 生成输出文件地址
    dir_outTa = paste0(dir_output, "/Ta_month_CQcity", year0,time0)
    if(!dir.exists(dir_outTa)){dir.create(dir_outTa, recursive = T)}
    
    ##### read Data for specific year:
    # 加载当前年份的气象站数据、高程数据、地表温度（LST）观测数据，以及地面气温观测数据。
    f_st1 = paste0(dir_input, "/station_CQcity_", year0, ".csv")              # 气象站位置信息文件（包含站点的坐标、高程等）
    f_Ts1 = paste0(dir_input, "/Ts.station.month",  ".csv")  # 气象站对应LST数据
    f_Ta1 = paste0(dir_input, "/Ta.station.", year0, time0,".csv")                # 气象站对应气温数据
    
    # 读取数据
    mat_st1 <- read.csv(file(f_st1, encoding = "UTF-8-BOM"))
    mat_Ts1 <- read.csv(file(f_Ts1, encoding = "UTF-8-BOM"))
    mat_Ta1 <- read.csv(file(f_Ta1, encoding = "UTF-8-BOM"))
    mat_st1 = read.csv(f_st1)
    mat_Ts1 = read.csv(f_Ts1)
    mat_Ta1 = read.csv(f_Ta1)
    colnames(mat_st1)[1] <- "name"
    colnames(mat_Ts1)[1] <- "name"
    colnames(mat_Ta1)[1] <- "name"
    
    unique(mat_Ta1[[1]])
    
    print("Cal. basis for samples")
    ##### calculate basic parameters
    #  构建预测基函数（Basis）
    data1 = as.data.frame(mat_st1)
    # 强制重命名列名
    colnames(data1) <- c("name", "x", "y", "dem")
    S = cbind(data1$x, data1$y)
    
    ##### Basis
    # f_B0pop 为基函数文件路径
    f_B0pop = file.path(dir_basis, paste0("B0pop_", year0, time0, ".rds"))
    if (!file.exists(f_B0pop)) {
      B0.pop <- tryCatch({
        B0.pop = basis(V, Tr, d, r, S.pop)
      }, error = function(e) {
        stop("基函数计算失败：", e$message)
      })
      saveRDS(B0.pop, f_B0pop)
    } else {
      B0.pop <- readRDS(f_B0pop)
    }
    
    # 加载基函数结果
    Q2=B0.pop$Q2             # 表示列表 B0.pop 中的 "Q2" 元素
    B.pop=B0.pop$B           #表示列表 B0.pop 中的 "B" 元素
    BQ2.pop=as.matrix(B.pop%*%Q2)
    K=B0.pop$K
    P=t(Q2)%*%K%*%Q2
    
    
    ##### Sample basis:构建样本点的基函数
    # 计算或加载样本点的样条基函数（Sample Basis Functions）
    f_B0 = file.path(dir_basis, paste0("B0_", year0, time0, ".rds"))
    B0 <- tryCatch({
      if (!file.exists(f_B0)) {
        B0 = basis(V, Tr, d, r, S)
        saveRDS(B0, f_B0)
      } else {
        B0 <- readRDS(f_B0)
      }
    }, error = function(e) {
      stop("基函数计算或加载失败：", e$message)
    })
    B = B0$B        # 提取 B 矩阵
    
    # 定义天数序列
    leap = is_leap_year(year0)                # 判断闰年/平年
    days_in_month = if (leap) days_in_month_leap else days_in_month_normal
    total_days = sum(days_in_month)
    v_days = 1:total_days
    
    ##### estimation for each day
    # 逐日处理
    start_day <- 1
    for (i0 in seq(start_day, length(v_days))){
      day0 = v_days[i0]
      
      # 计算当前 day0 属于哪一个月
      cum_days = cumsum(days_in_month)
      month_index = which(cum_days >= day0)[1]   # 当前天属于的月份索引（1~12）
      month_str <- sprintf("%02d", month_index)  # 格式化为两位数
      
      # 打印处理信息
      print(paste("Processing day:", day0, "-> month:", month_index))
      
      # 构造数据
      data = data.frame(
        name = mat_st1$name,
        plst = mat_Ts1[[month_index + 1]],  # 月LST列（跳过name列）
        dem = mat_st1$dem,
        at = mat_Ta1[[day0 + 1]]  # 日气温数据（跳过name列）
      )
      
      # 因变量
      Y = data$at                     # 气象要素（如温度、湿度等）作为因变量
      
      # 标准化自变量
      X = cbind(data$plst, data$dem)
      mean1 = mean(X[, 1])
      mean2 = mean(X[, 2])
      sd1 = sd(X[, 1])
      sd2 = sd(X[, 2])
      X = scale(X, scale=TRUE)   # standardized covariates
      X = cbind(1, X)            # add intercept
      
      
      # SVCM 模型估算：拟合 SVCM 空间变系数模型
      f_result.svcm = paste0(dir_basis, "/svcm_param_month", year0, time0, Ds[day0], ".rds")
      if(!file.exists(f_result.svcm)) {
        lambda.svcm = svcmsp.sel(B, Q2, K, X, Y, sign.coeff = c(0, 1, -1))
        Result.svcm = fit.svcmsp(B, Q2, K, lambda.svcm, X, Y, sign.coeff = c(0, 1, -1)) 
        saveRDS(Result.svcm, f_result.svcm)
      } else {
        Result.svcm = readRDS(f_result.svcm)
      }
      
      ##### read Landsat LST 读取地表温度（Landsat LST）
      # land surface temperature (LST)
      
      # 构造月 LST 栅格文件路径
      dir_lst_month <- paste0("/data/wanzhou/Landsat LST/", "Month2014_2024_RF/")
      f_lst_month <- paste0(dir_lst_month, "CQcity_LST_Month_", month_str, ".tif")
      
      ##### crop LST 读取并裁剪
      ras_lst <- raster(f_lst_month)
      ras_lst <- raster::crop(ras_lst, ras_dem1)   
      print(paste("Processed Landsat LST for day", day0, "-> month", month_str))
      
      
      # 创建空的预测结果栅格模板
      ras_out_svcm = ras_lst         # 继承 ras_lst 的所有空间信息（分辨率、范围、投影等）
      ras_out_svcm[] = NA            # 清空所有值，以便填入预测结果
      
      
      ##### prediction Ta tile by tile
      # 逐块预测地面气象要素
      for(i in c(1:len_tile)) {
        ista = v_start[i]
        iend = v_end[i]
        f_b0pred0 = paste0(dir_outdem, "/b0pred0_", i, ".rds")
        f_predgrid0 = paste0(dir_outdem, "/predgrid0_", i, ".rds")
        if(!file.exists(f_b0pred0)) {
          
          ##### change projection of the coordinates
          # 将 raster 栅格的经纬度（）转为 shapefile 的投影坐标（WGS84），用于构建预测位置 xx、yy
          dat_coord = changeCoord(log@data@values[ista:iend], lat@data@values[ista:iend], proj2, proj1)
          
          # 获取当前 tile 的栅格经纬度坐标（log 和 lat 是 raster 对象）
          lon_vals <- log@data@values[ista:iend]
          lat_vals <- lat@data@values[ista:iend]
          # 创建 SpatialPoints 对象，用于提取 DEM 值（仍然用原始经纬度）
          coords_ll = data.frame(x = lon_vals, y = lat_vals)
          pts_ll = SpatialPoints(coords_ll, proj4string = proj2)
          # 使用 extract() 从 DEM 栅格中提取高程值
          dem_vals = raster::extract(ras_dem1, pts_ll)
          
          # 加载或计算 tile 的基函数与坐标
          # predgrid0=data.frame(id=c(ista:iend),xx=dat_coord@coords[, 1],yy=dat_coord@coords[, 2],dem=ras_dem1[ista:iend])
          predgrid0 = data.frame(
            id  = c(ista:iend),
            xx  = dat_coord@coords[, 1],
            yy  = dat_coord@coords[, 2],
            dem = dem_vals
          )
          S.pred0 = cbind(predgrid0$xx, predgrid0$yy)
          B0.pred0=basis(V, Tr, d, r, S.pred0)
          B.pred0=B0.pred0$B
          
          saveRDS(predgrid0, f_predgrid0)
          saveRDS(B0.pred0, f_b0pred0)
        } else {
          predgrid0 = readRDS(f_predgrid0)
          B0.pred0 = readRDS(f_b0pred0)
          B.pred0=B0.pred0$B
        }
        
        
        # 筛选有效像元
        idx_in0 = B0.pred0$Ind.inside
        predgrid0$plst = ras_lst[ista:iend]
        predgrid0 = predgrid0[idx_in0, ]
        idx_comp = complete.cases(predgrid0)
        idx_sample = which(idx_comp==TRUE)
        predgrid0 = predgrid0[idx_sample, ]
        
        ##### create X and S预测当日气象要素：气温
        # # 预测点的自变量矩阵X.pred，其维度为n_samples x 2
        X.pred = cbind(predgrid0$plst, predgrid0$dem)
        
        ##### standardize and add intercept
        # 对预测数据的协变量 X.pred 进行标准化处理
        X.pred = scale(X.pred, center=c(mean1, mean2), scale=c(sd1, sd2)) 
        
        # 添加一列常数 1，表示模型中的截距项：n_samples x 3的矩阵，第一列为常数1
        X.pred = cbind(1, X.pred)
        
        beta.hat.pred = B.pred0[idx_sample, ]%*%Result.svcm$gamma
        yhat.pred = rowSums(X.pred*beta.hat.pred)
        ras_out_svcm[predgrid0$id] = yhat.pred
      }
      
      # 保存预测结果：输出预测气温栅格
      f_Ta_svcm = paste0(dir_outTa, "/CQ_Ta_SVCMsp_",year0,time0,"_",Ds[day0],".tif")
      writeRaster(ras_out_svcm, f_Ta_svcm, overwrite=TRUE)
      
    }
  }
}
