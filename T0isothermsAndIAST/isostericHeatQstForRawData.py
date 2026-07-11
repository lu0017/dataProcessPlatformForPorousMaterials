##每一个子平台固定开头，用于找到依赖
import os
import sys
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from common import *
import constantsAndName as const

import T1fileSource.fileOperation as fl
import T1dataProcessSource.dataOperation as dop

# ----------------------------
# 直接对实验数据通过Clausius–Clapeyron来计算Qst, 不使用拟合模型
# ----------------------------
# --------------------------
# 绘图
# --------------------------
def plotQstByRawData(Qst):
    # -----------------------------
    # 绘图
    # -----------------------------
    q = Qst["q_mmol_g"]
    qst = Qst["Qst_J_mol"]
    plt.figure()
    plt.plot(q, qst * const.J2kJ, marker="o", label="Qst")
    plt.xlabel("q (mmol/g)")
    plt.ylabel("Qst (kJ/mol)")
    plt.title("Isosteric Heat of Adsorption (Qst)")
    plt.grid(True)
    plt.legend()
    plt.show(block=True)

def runCalculateQstByRawData(df):
    P_data, uptake_data, T_list = dop.build_P_and_uptake_data(df)

    q_step = 0.1   # 吸附量步长 (mmol/g)

    q_min = max([np.nanmin(uptake_data[T]) for T in uptake_data])
    q_max = min([np.nanmax(uptake_data[T]) for T in uptake_data])
    print(f"公共吸附量范围: q_min = {q_min:.3f}, q_max = {q_max:.3f}")

    q_range = np.arange(q_min, q_max, q_step)
    # --------------------------
    # 计算 Qst
    # --------------------------
    Qst_list = []
    valid_q = []
    results_full = []

    for q_target in q_range:
        P_at_q = []
        P_dict = {}   #用于导出excel
        skip = False
        for T in T_list:
            uptake = uptake_data[T]
            P = P_data[T]
            P_clean, uptake_clean = dop.clean_list(P, uptake)

            # 检查长度是否足够插值
            if len(uptake_clean) < 2:
                skip = True
                break

            # 排序，确保单调递增
            P_sorted, uptake_sorted = dop.sortData(P_clean, uptake_clean)

            # 插值 此处x = uptake, y = Pressure
            # f_interp = PchipInterpolator(uptake_sorted, P_sorted, extrapolate=False) #Hermite 插值，效果不错
            #f_interp = CubicSpline(uptake_sorted, P_sorted, extrapolate=False)  #样条曲线插值，效果最顺滑
            f_interp = interp1d(uptake_sorted, P_sorted, bounds_error=True)  #线性插值，低吸附段较好，需要使用大的q_step
            P_q = f_interp(q_target)
            P_at_q.append(P_q)

            col_name = f"P_T{int(T + const.K2C)}" if float(T).is_integer() else f"P_T{int(T+ const.K2C)}"
            P_dict[col_name] = P_q
        if skip:
            continue
        
        P_at_q = np.array(P_at_q)
        lnP = np.log(P_at_q)  # ln(P)
        invT = 1 / T_list  # 1/T

        # Clausius–Clapeyron 拟合
        slope, intercept, r_value, p_value, std_err = linregress(invT, lnP) #线性拟合求斜率

        # Qst 计算
        Qst = -slope * const.R 
        Qst_list.append(Qst)
        valid_q.append(q_target)
        row = {
            "q_mmol_g": q_target,
            "Qst_J_mol": Qst,
            "Qst_kJ_mol": Qst/1000.0,
            "R2": r_value**2,
            "slope": slope,
            "intercept": intercept
        }
        row.update(P_dict)
        results_full.append(row)
    results = pd.DataFrame({"q_mmol_g": valid_q, "Qst_J_mol": Qst_list})
    # return results
    return pd.DataFrame(results_full)

def main(file_path=None):
    singleFile = True

    if singleFile:
        # 选择文件
        sheet_name = "uio-66 Co 10"
        file = fl.getFile()
        listData, out_path, _ = fl.readFileBySheet(file, sheet_name, expand = const.QST_EXP_SHEET)
        Qst_exp = runCalculateQstByRawData(listData)
        newPath = fl.check_or_create_folder(out_path, sub_folder_name=const.QST_FOLDER)
        newFilePath = fl.get_expanded_name(newPath, sheet_name, expand=const.QST_EXP_SHEET)
        fl.export_to_excel_auto(Qst_exp, filename=newFilePath, sheet_name=sheet_name)
        plotQstByRawData(Qst_exp)

    else:
        file_path = fl.getFile()
        xls = pd.ExcelFile(file_path)
        fileName = os.path.splitext(os.path.basename(file_path))[0]
        out_path = fl.get_expanded_name(file_path, fileName, expand=const.QST_EXP_SHEET, expandPos=True, type="xlsx")
        for sheet in xls.sheet_names:
            if fl.should_skip(sheet): 
                continue
            listData, _ , _ = fl.readFileBySheet(file_path, sheet)

            if isinstance(listData, dict):
                listData0 = pd.DataFrame(list(listData.values())[0])  # 取第一个 sheet

            if listData0.empty:
                print("DataFrame is empty")
                continue

            Qst_exp = runCalculateQstByRawData(listData)
            fl.export_to_excel_auto(Qst_exp, filename=out_path, sheet_name=sheet)
            plotQstByRawData(Qst_exp)

if __name__ == "__main__":
    f = sys.argv[1] if len(sys.argv) > 1 else None
    main(f)