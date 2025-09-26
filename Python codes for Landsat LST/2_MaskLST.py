import os
import math
import numpy as np
import rasterio
from rasterio.mask import mask
import geopandas as gpd
import matplotlib.pyplot as plt
from shapely.geometry import mapping
import cartopy.crs as ccrs
from matplotlib import ticker

# 函数说明：根据CQcity的边界shp数据提取得到各月份栅格数据
# --------------------------------------------------------

# 设置字体为 Times New Roman
plt.rcParams['font.family'] = 'Times New Roman'

# 月份映射（用于英文标题）
month_names = {
    "01": "January", "02": "February", "03": "March", "04": "April",
    "05": "May", "06": "June", "07": "July", "08": "August",
    "09": "September", "10": "October", "11": "November", "12": "December"
}

# ========== 配置路径 ==========
lst_folder = r"D:\0 DataBase\20_Landsat LST\2015-2024"  # 原始 LST 数据文件夹
output_tif_dir = os.path.join(lst_folder, "CQcity")  # 掩膜后 tif 保存路径
output_fig_dir = os.path.join(lst_folder, "Figures")  # 图像保存路径
os.makedirs(output_tif_dir, exist_ok=True)
os.makedirs(output_fig_dir, exist_ok=True)

# ========== 加载重庆市边界 ==========
CQcity_shp_path = r"D:\0 DataBase\0 Chongqin Database\1 Boundary\CQcity_Boundary.shp"
CQcity = gpd.read_file(CQcity_shp_path)
CQcity = CQcity.to_crs(epsg=4326)  # 与 tif 文件 CRS 保持一致

# ==== 读取边界 ====
boundary = gpd.read_file(CQcity_shp_path).to_crs("EPSG:4326")
minx, miny, maxx, maxy = boundary.total_bounds
delta = 0.02
extent = [minx-delta , maxx+delta , miny-delta , maxy+delta ]

# ========== 掩膜并保存函数 ==========
def mask_and_save_tif(tif_path, boundary_gdf, save_folder):
    with rasterio.open(tif_path) as src:
        boundary_geom = [mapping(geom) for geom in boundary_gdf.geometry]
        out_image, out_transform = mask(src, boundary_geom, crop=True)
        out_meta = src.meta.copy()
        out_meta.update({
            "driver": "GTiff",
            "height": out_image.shape[1],
            "width": out_image.shape[2],
            "transform": out_transform
        })

        # 构建保存路径
        filename = os.path.basename(tif_path)
        save_path = os.path.join(save_folder, filename)

        # 保存掩膜后的 tif
        with rasterio.open(save_path, "w", **out_meta) as dest:
            dest.write(out_image)
    return save_path


# ========== 绘图函数 ==========
def plot_lst_month(tif_path, month, boundary_gdf, save_dir):
    with rasterio.open(tif_path) as src:
        # 将边界转为 GeoJSON 格式
        boundary_geom = [mapping(geom) for geom in boundary_gdf.geometry]

        # ✅ 使用 filled=False 获取掩膜数组（边界外为 masked）
        out_image, out_transform = mask(src, boundary_geom, crop=True, filled=False)
        lst_data = out_image[0].astype(np.float32)

        # ✅ 转换为 masked 数组
        masked_lst = np.ma.masked_array(lst_data, mask=np.isnan(lst_data) | (lst_data < -1e10))

        # ✅ 检查是否全为 NaN（masked）
        if masked_lst.mask.all():
            print(f"⚠️ 文件 {tif_path} 所有像元为 NaN，跳过绘图")
            return

        # ✅ 真实 min/max
        lst_min = masked_lst.min()
        lst_max = masked_lst.max()
        vmin = math.floor(lst_min)
        vmax = math.ceil(lst_max)

        # ✅ 获取绘图范围
        extent = [
            out_transform.c,
            out_transform.c + out_transform.a * masked_lst.shape[1],
            out_transform.f + out_transform.e * masked_lst.shape[0],
            out_transform.f
        ]


    # ✅ 开始绘图
    fig, ax = plt.subplots(figsize=(10, 8), subplot_kw={'projection': ccrs.PlateCarree()})
    ax.set_extent(extent, crs=ccrs.PlateCarree())

    im = ax.imshow(
        masked_lst,
        cmap='coolwarm',
        extent=extent,
        vmin=vmin,
        vmax=vmax,
        origin='upper'
    )

    # ✅ 绘制边界（确保有线）
    boundary_gdf.boundary.plot(ax=ax, edgecolor='black', linewidth=1)

    # ✅ 标题和坐标
    month_str = month_names.get(month, month)
    title = f"LST in {month_str}"
    ax.set_title(title, fontsize=17)
    ax.grid(True, linestyle='--', linewidth=0.5)
    ax.tick_params(labelsize=12)

    # ✅ 添加颜色条
    cbar = plt.colorbar(im, ax=ax, shrink=0.95, pad=0.05)
    cbar.set_label("LST (°C)", fontsize=14)
    cbar.locator = ticker.MultipleLocator(5)  # ✅ 设置间隔为 1°C
    cbar.update_ticks()
    cbar.ax.tick_params(labelsize=12)

    # ✅ 左上角文字
    text_str = f"LST max = {lst_max:.2f}°C\nLST min = {lst_min:.2f}°C"
    ax.text(0.02, 0.98, text_str,
            transform=ax.transAxes,
            fontsize=16,
            verticalalignment='top',
            horizontalalignment='left',
            # bbox=dict(facecolor='lightgray', alpha=0.7)
            )

    # ✅ 设置经纬度刻度为整数，间隔为 0.2°
    xticks = np.round(np.arange(extent[0]+0.1, extent[1] , 0.2), 1)
    yticks = np.round(np.arange(extent[2]+0.1, extent[3] , 0.2), 1)
    ax.set_xticks(xticks, crs=ccrs.PlateCarree())
    ax.set_yticks(yticks, crs=ccrs.PlateCarree())
    ax.set_xticklabels([f"{x:.1f}°" for x in xticks], fontsize=13)
    ax.set_yticklabels([f"{y:.1f}°" for y in yticks], fontsize=13)

    ax.set_xlabel("Longitude (°)", fontsize=14)
    ax.set_ylabel("Latitude (°)", fontsize=14)

    # ✅ 保存
    save_path = os.path.join(save_dir, f"LST_CQcity_Month{month:02d}.png")
    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    plt.close()
    print(f"✅ 已保存图像：{save_path} | 范围: {vmin}°C ~ {vmax}°C")


# ========== 循环处理每个月 ==========
for month in range(1, 12):  # 修改为全年
    tif_file = os.path.join(lst_folder, f"Mean_LST_Month_{month}.tif")
    if os.path.exists(tif_file):
        print(f"📂 正在处理：{tif_file}")

        # 掩膜提取并保存
        masked_tif_path = mask_and_save_tif(tif_file, CQcity, output_tif_dir)

        # 绘图
        plot_lst_month(masked_tif_path, month, CQcity, output_fig_dir)
    else:
        print(f"⚠️ 文件不存在：{tif_file}")