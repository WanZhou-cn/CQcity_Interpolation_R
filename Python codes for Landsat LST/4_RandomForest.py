import os
import numpy as np
import rasterio
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error

# ========== 路径设置 ==========
dem_path = r"D:\0 DataBase\0 Chongqin Database\1 Boundary\DEM_30m_CQcity.tif"
lst_folder = r"D:\Landsat LST\Month2014_2024"

output_folder = r"D:\Landsat LST\Month2014_2024_RF1"
os.makedirs(output_folder, exist_ok=True)


def calculate_nan_ratio(array, reference=None, name=""):
    """
    计算研究区域内的缺测率。

    参数：
    - array: 要评估的 LST 数据（2D数组）
    - reference: 掩膜数据（如 DEM），若为 None，自动使用 array 自身的非 NaN 区域
    - name: 标签
    """
    if reference is None:
        ref_mask = ~np.isnan(array)
    else:
        # 如果 reference 是整数型（如 DEM），0 作为 nodata
        if np.issubdtype(reference.dtype, np.integer):
            ref_mask = reference > 0
        else:
            ref_mask = ~np.isnan(reference)

    valid_pixels = array[ref_mask]
    total_valid = valid_pixels.size
    nan_count = np.isnan(valid_pixels).sum()

    ratio = nan_count / total_valid if total_valid > 0 else 0
    print(f"📊 [{name}] 区域内缺测: {nan_count} / {total_valid} => 缺测率: {ratio:.2%}")
    return ratio


def evaluate_model(model, X_test, y_test):

    """评估模型性能"""
    y_pred = model.predict(X_test)
    r2 = r2_score(y_test, y_pred)
    rmse = np.sqrt(mean_squared_error(y_test, y_pred))
    mae = mean_absolute_error(y_test, y_pred)

    print(f"OOB R² Score: {model.oob_score_:.4f}" if hasattr(model, "oob_score_") else "")
    print(f"Test R² Score: {r2:.4f}")
    print(f"RMSE: {rmse:.4f}")
    print(f"MAE: {mae:.4f}")
    return r2, rmse, mae


# ========== 读取 DEM ==========
with rasterio.open(dem_path) as dem_src:
    dem_data = dem_src.read(1)
    dem_meta = dem_src.meta

rows, cols = dem_data.shape

# ========== 读取全部 LST ==========
lst_data_dict = {}
for month in range(1, 13):
    file_path = os.path.join(lst_folder, f'CQcity_LST_Month_{month}.tif')
    with rasterio.open(file_path) as src:
        lst_data_dict[month] = src.read(1).astype(np.float32)


# ========== 插值主流程：插值并保存 ==========
# -----------------------------------------------------
for m in range(1, 13):
    output_path = os.path.join(output_folder, f'RF_LST_Month_{m:02d}.tif')

    # ✅ 已插值文件存在则跳过
    if os.path.exists(output_path):
        print(f"⏩ 文件已存在，跳过：{output_path}")
        continue

    print(f"\n📌 正在插值：{m}月")
    prev_m = 12 if m == 1 else m - 1
    next_m = 1 if m == 12 else m + 1

    lst_cur = lst_data_dict[m]
    lst_prev = lst_data_dict[prev_m]
    lst_next = lst_data_dict[next_m]

    # 使用 DEM 作为参考掩膜，只统计研究区域内有效数据
    calculate_nan_ratio(lst_cur, reference=dem_data, name=f"Month {m} (Before)")

    X_train = []
    y_train = []
    X_pred = []
    pred_indices = []

    for i in range(rows):
        for j in range(cols):
            dem = dem_data[i, j]
            cur = lst_cur[i, j]
            prev = lst_prev[i, j]
            nxt = lst_next[i, j]

            if np.isnan(dem) or dem <= 0:
                continue

            # 训练样本：当前月非NaN，前后月也非NaN
            if not np.isnan(cur) and not np.isnan(prev) and not np.isnan(nxt):
                X_train.append([dem, prev, nxt, i, j])
                y_train.append(cur)

            # 预测样本：当前月 NaN，前后月不为 NaN
            elif np.isnan(cur) and not np.isnan(prev) and not np.isnan(nxt):
                X_pred.append([dem, prev, nxt, i, j])
                pred_indices.append((i, j))

    lst_filled = lst_cur.copy()

    # ===== 阶段一：预测当前月 NaN，前后月有值 =====
    if len(X_train) > 100 and len(X_pred) > 0:
        # 拆分训练集用于评估
        X_train_part, X_test_part, y_train_part, y_test_part = train_test_split(
            X_train, y_train, test_size=0.2, random_state=42
        )
        rf = RandomForestRegressor(
            n_estimators=100, max_depth=15, random_state=42, oob_score=True, n_jobs=-1
        )
        rf.fit(X_train_part, y_train_part)

        print(f"📈 阶段一模型评估（月 {m}）")
        evaluate_model(rf, X_test_part, y_test_part)

        y_pred = rf.predict(X_pred)
        for idx, (i, j) in enumerate(pred_indices):
            lst_filled[i, j] = y_pred[idx]
    else:
        print(f"⚠️ 月 {m} 阶段一样本不足，跳过。")

    # ===== 阶段二：当前月仍为 NaN，DEM 存在 =====
    X2_train = []
    y2_train = []
    X2_pred = []
    pred2_indices = []

    for i in range(rows):
        for j in range(cols):
            dem = dem_data[i, j]
            cur = lst_cur[i, j]

            if np.isnan(dem) or dem <= 0:
                continue

            if not np.isnan(cur):
                X2_train.append([dem, i, j])
                y2_train.append(cur)
            elif np.isnan(lst_filled[i, j]):
                X2_pred.append([dem, i, j])
                pred2_indices.append((i, j))

    if len(X2_train) > 100 and len(X2_pred) > 0:
        rf2 = RandomForestRegressor(
            n_estimators=100, max_depth=10, random_state=42, oob_score=True, n_jobs=-1
        )
        rf2.fit(X2_train, y2_train)

        print(f"📈 阶段二模型评估（月 {m}）")
        X2_tr, X2_te, y2_tr, y2_te = train_test_split(X2_train, y2_train, test_size=0.2, random_state=42)
        rf_eval = RandomForestRegressor(n_estimators=100, random_state=42)
        rf_eval.fit(X2_tr, y2_tr)
        evaluate_model(rf_eval, X2_te, y2_te)

        y2_pred = rf2.predict(X2_pred)
        for idx, (i, j) in enumerate(pred2_indices):
            lst_filled[i, j] = y2_pred[idx]
    else:
        print(f"⚠️ 第二阶段样本不足，无法进一步填补。")

    # 使用 DEM 作为参考掩膜，只统计研究区域内有效数据
    calculate_nan_ratio(lst_filled, reference=dem_data, name=f"Month {m} (After)")

    # ✅ 屏蔽边界外像元
    lst_filled[dem_data <= 0] = np.nan

    # 保存插值结果
    output_path = os.path.join(output_folder, f'RF_LST_Month_{m:02d}.tif')
    out_meta = dem_meta.copy()
    out_meta.update({
        'dtype': 'float32',
        'nodata': np.nan,
        'count': 1
    })

    with rasterio.open(output_path, 'w', **out_meta) as dst:
        dst.write(lst_filled, 1)

    print(f"✅ 已保存插值结果：{output_path}")