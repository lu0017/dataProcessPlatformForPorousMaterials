# ----------------------------
# 使用origin软件进行dual site langmuir拟合，需要不同温度下的压强坐标一致
# 改代码可以导入不同温度数据，生成统一压强坐标的数据
# ----------------------------
import os
import sys
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from common import *
import constantsAndName as const

import T1fileSource.fileOperation as fl
import T1dataProcessSource.dataOperation as dop

def runInterpolateDsl(df, pressureUnit="kPa", plotFlag=True):
    """默认单位 kPa"""
    # 2. 提取三组数据
    #kPa
    scale_P = 1
    p_lists = {}
    q_lists = {}
    cleaned_p = {}
    cleaned_q = {}

    P_data, uptake_data, T_list = dop.build_P_and_uptake_data(df)

    for T in T_list:
        p_lists[T] = P_data[T]
        q_lists[T] = uptake_data[T]

    # #Pa
    if pressureUnit == "Pa":
        scale_P = 1000

    for T in T_list:
        cleaned_p[T], cleaned_q[T] = dop.clean_list(p_lists[T], q_lists[T], scale_x=scale_P)


    # 计算所有 T 的 p 最小值的最大值，和最大值的最小值
    p_min = max([min(cleaned_p[T]) for T in T_list])
    p_max = min([max(cleaned_p[T]) for T in T_list])

    # 生成公共 P 点
    p_common = np.linspace(p_min, p_max, 50)  # 50个点，可改

    # 4. 插值函数
    f_funcs = {}  # 字典存储每个 T 对应的插值函数

    # 5. 生成插值结果
    for T in T_list:
        f_funcs[T] = PchipInterpolator(cleaned_p[T], cleaned_q[T], extrapolate=False) #Hermite 插值，效果不错
        # f_funcs[T] = CubicSpline(cleaned_p[T], cleaned_q[T], extrapolate=False)  #样条曲线插值，效果最顺滑
        # f_funcs[T] = interp1d(cleaned_p[T], cleaned_q[T], kind="linear")

    # 初始化字典存储列
    data_dict = {"Pressure_common": p_common}

    # 遍历 T_list，动态生成列名和插值数据
    for T in T_list:
        data_dict[f"T{int(T + const.K2C)}_interp"] = f_funcs[T](p_common)

    # 创建 DataFrame
    df_interp = pd.DataFrame(data_dict)

    if plotFlag:
        plt.figure(figsize=(6, 4))

        # 绘制实验数据点
        for T in T_list:
            plt.scatter(cleaned_p[T], cleaned_q[T], label=f"T={T}°C exp")

        # 绘制插值曲线
        for T in T_list:
            plt.plot(p_common, f_funcs[T](p_common), label=f"T={T}°C Interp")

        plt.xlabel("P (kPa)")
        plt.ylabel("q (mmol/g)")
        plt.title("Interpolated and experiment")
        plt.legend()
        plt.grid(True)
        plt.show(block=True)

    return df_interp

def main(file_path=None):
    singleFile = True

    if singleFile:
        # # 选择文件
        sheet_name = "CC-Hy-550_60_5-650_15_5-1"
        file = fl.getFile()
        listData, out_path,_ = fl.readFileBySheet(file, sheet_name, expand = "Interp")
        df_interp = runInterpolateDsl(listData)
        # 6. 保存结果
        fl.export_to_excel_auto(df_interp, filename=out_path)

    else:
        file_path = fl.getFile()
        xls = pd.ExcelFile(file_path)
        fileName = os.path.splitext(os.path.basename(file_path))[0]
        out_path = fl.get_expanded_name(file_path, fileName, expand="Interp", expandPos=True, type="xlsx")
        for sheet in xls.sheet_names:
            if fl.should_skip(sheet): 
                continue

            listData, _, _ = fl.readFileBySheet(file_path, sheet)

            if isinstance(listData, dict):
                df0 = pd.DataFrame(list(listData.values())[0])  # 取第一个 sheet

            if df0.empty:
                print("DataFrame is empty")
                continue

            df_interp = runInterpolateDsl(listData, plotFlag=False)
            fl.export_to_excel_auto(df_interp, filename=out_path, sheet_name=sheet)

if __name__ == "__main__":
    f = sys.argv[1] if len(sys.argv) > 1 else None
    main(f)
