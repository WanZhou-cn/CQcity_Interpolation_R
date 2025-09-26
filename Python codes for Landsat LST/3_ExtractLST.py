import pandas as pd
import rasterio
from rasterio.transform import from_origin
from tqdm import tqdm
import matplotlib.pyplot as plt
from pyproj import Transformer
import os

# 设置字体为 Times New Roman
plt.rcParams['font.family'] = 'Times New Roman'

def extract_monthly_lst_to_csv(
    station_txt_path,
    lst_folder_path,
    output_csv_path,
    month_start=1,
    month_end=12
):
    # 读取站点信息
    stations = pd.read_csv(station_txt_path, sep='\t')
    station_names = stations['Name'].astype(str).tolist()
    lons = stations['Longitude'].tolist()
    lats = stations['Latitude'].tolist()

    # 自动识别年份（可选）
    folder_name = os.path.basename(lst_folder_path.strip("/\\"))
    year = ''.join(filter(str.isdigit, folder_name))

    # 初始化结果字典
    data_dict = {'name': station_names}

    # 遍历每个月
    for month in tqdm(range(month_start, month_end + 1), desc=f"Extracting LST"):
        month_str = f"{month:01d}"
        tif_name = f"Mean_LST_Month_{month_str}.tif"
        tif_path = os.path.join(lst_folder_path, tif_name)

        if not os.path.exists(tif_path):
            print(f"[WARNING] 未找到文件：{tif_path}")
            data_dict[f'm{month_str}'] = [None] * len(station_names)
            continue

        # 打开 tif 并提取 LST 值
        with rasterio.open(tif_path) as src:
            values = []
            for lon, lat in zip(lons, lats):
                try:
                    row, col = src.index(lon, lat)
                    value = src.read(1)[row, col]
                    if src.nodata is not None and value == src.nodata:
                        values.append(None)
                    else:
                        values.append(value)  # 直接使用摄氏度值
                except Exception as e:
                    values.append(None)

        data_dict[f'm{month_str}'] = values

    # 保存为 CSV
    df = pd.DataFrame(data_dict)
    os.makedirs(os.path.dirname(output_csv_path), exist_ok=True)
    df.to_csv(output_csv_path, index=False, float_format="%.2f")
    print(f"[INFO] 成功保存：{output_csv_path}")



# 主程序
# -----------------------------------------------------
month_start = 1
month_end = 12

# 提取站点对应LST数据
extract_monthly_lst_to_csv(
    station_txt_path=r"D:\0 DataBase\20_Landsat LST\CQcity\Station_CQcity.txt",
    lst_folder_path=r"D:\0 DataBase\20_Landsat LST\2014-2024\ChongqingBox",
    output_csv_path=r"D:\0 DataBase\20_Landsat LST\2014-2024\lst.station.month.csv",
    month_start=month_start,
    month_end=month_end
)