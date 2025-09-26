import os
import glob
import numpy as np
import rasterio
import math
import geopandas as gpd
import matplotlib.pyplot as plt
from rasterio.plot import show
from rasterio.mask import mask

# 设置字体为 Times New Roman
plt.rcParams['font.family'] = 'Times New Roman'

# ========== 配置路径 ==========
lst_folder = r"D:\0 DataBase\20_Landsat LST\2015-2024"  # 存放 LST GeoTIFF 的文件夹


output_dir = r"D:\0 DataBase\20_Landsat LST\2015-2024\Figures"
os.makedirs(output_dir, exist_ok=True)  # 确保输出文件夹存在

# ========== 加载重庆市边界 ==========
CQ_shp_path = r"D:\0 DataBase\0 Chongqin Database\1 Boundary\Chongqing.shp"
chongqing = gpd.read_file(CQ_shp_path)
chongqing = chongqing.to_crs(epsg=4326)  # 确保坐标一致

CQcity_shp_path = r"D:\0 DataBase\0 Chongqin Database\1 Boundary\CQcity_Boundary.shp"
CQcity = gpd.read_file(CQcity_shp_path)
CQcity = CQcity.to_crs(epsg=4326)  # 确保坐标一致


# ========== 绘图函数 ==========
def plot_lst_month(tif_path, month, boundary_gdf1, boundary_gdf2, save_dir):
    with rasterio.open(tif_path) as src:
        lst_data = src.read(1)
        lst_data = np.where(lst_data == src.nodata, np.nan, lst_data)

        # 掩膜 NaN 区域
        masked_lst = np.ma.masked_invalid(lst_data)

        # 获取图像地理范围
        extent = [src.bounds.left, src.bounds.right, src.bounds.bottom, src.bounds.top]

        # 获取真实 min/max，向下/上取整
        real_min = np.nanmin(lst_data)
        real_max = np.nanmax(lst_data)
        vmin = math.floor(real_min)
        vmax = math.ceil(real_max)

    # 创建图像
    fig, ax = plt.subplots(figsize=(10, 8))
    im = ax.imshow(
        masked_lst,
        cmap='plasma',
        extent=extent,
        vmin=vmin,
        vmax=vmax,
        origin='upper'
    )

    # 叠加边界线（保持完整底图）
    boundary_gdf1.boundary.plot(ax=ax, edgecolor='black', linewidth=1)
    boundary_gdf2.boundary.plot(ax=ax, edgecolor='black', linewidth=1.2)

    # 添加标题与标签
    ax.set_title(f"LST in Month {month} (°C)", fontsize=17)
    ax.set_xlabel("Longitude(°)", fontsize=15)
    ax.set_ylabel("Latitude(°)", fontsize=15)
    ax.grid(True, linestyle='--', linewidth=0.5)
    ax.tick_params(labelsize=12)

    # ✅ 颜色条设置真实范围
    cbar = plt.colorbar(im, ax=ax, shrink=0.9, pad=0.05)
    cbar.set_label("LST (°C)", fontsize=14)
    cbar.ax.tick_params(labelsize=12)
    cbar.set_ticks(np.linspace(vmin, vmax, num=6))

    # 保存图像
    save_path = os.path.join(save_dir, f"LST_Chongqing_Month{month:02d}.png")
    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    plt.close()
    print(f"✅ 已保存：{save_path} | 范围: {vmin}°C ~ {vmax}°C")


# ========== 循环 1~12月绘图 ==========
for month in range(1, 4):
    tif_file = os.path.join(lst_folder, f"Mean_LST_Month_{month}.tif")
    if os.path.exists(tif_file):
        plot_lst_month(tif_file, month, chongqing, CQcity, output_dir)
    else:
        print(f"⚠️ 文件不存在：{tif_file}")