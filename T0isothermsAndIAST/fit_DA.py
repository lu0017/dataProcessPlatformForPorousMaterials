# -*- coding: utf-8 -*-
"""
Dubinin–Astakhov (D-A) function fitting  + Qst/H/S/G 计算
"""
##平台固定开头，用于找到依赖
import os
import sys
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from common import *
import constantsAndName as const

import T1fileSource.fileOperation as fl
import T1dataProcessSource.dataOperation as dop

import modelAndFit as mf


# ==========================================
# 工具函数
# ==========================================

def poltDA(fits_data, rawData, fixedFlag = False):

    grouped_data = {T: gdf for T, gdf in rawData.groupby("T_K")}
    T_list = sorted(grouped_data.keys())
    plt.figure(figsize=(6,4))
    for T in T_list:
        gdf=grouped_data[T]
        plt.scatter(gdf["P_Pa"]/1000, gdf["q_mmol_g"], label=f"T={T-273.15:.1f}°C exp")
    for T in T_list:
        q_fit = fits_data[fits_data["T_K"]==T]["q_fit_mmol_g"].values
        P_fit = fits_data[fits_data["T_K"]==T]["P_Pa"].values
        plt.plot(P_fit/1000,q_fit,label=f"T={T-273.15}°C D-A")
    plt.xlabel("P (kPa)")
    plt.ylabel("q (mmol/g)")
    plt.title("D-A and experiment")
    plt.legend()
    plt.grid(True)
    if fixedFlag:
        plt.show(block=True)
    else:
        plt.show(block=False)


def export_DA_params_custom(DApara_or_dict, file_name="results.xlsx", sheet_name="DA_params", singleFlag=False):
    """
    导出 D-A 拟合参数到 Excel，自定义列名和单位：
        W0(cm3/g), E(J/mol), n(-)
    参数：
        DApara_or_dict: 
            - singleFlag=False 时，为 lmfit MinimizerResult（全局拟合）
            - singleFlage=True 时，为 dict[T] = MinimizerResult（单温度拟合）
        file_name: Excel 文件名
        sheet_name: 工作表名
        singleFlag: 是否为单温度拟合
    """
    rows = []

    if not singleFlag:
        # 全局拟合
        W0, E, n = getParaForDA(DApara_or_dict.params)
        rows.append({
            "Temperature_K": "Global",
            "W0_cm3_g": W0,
            "E_J_mol": E,
            "n": n
        })
    else:
        # 单温度拟合
        for T, daPara in DApara_or_dict.items():
            W0, E, n = getParaForDA(daPara.params)
            rows.append({
                "Temperature_K": T,
                "W0_cm3_g": W0,
                "E_J_mol": E,
                "n": n
            })

    df_export = pd.DataFrame(rows)

    # 使用已有自动导出函数
    fl.export_to_excel_auto(df_export, file_name, sheet_name=sheet_name)

def getCO2LiquidPhaseDensityForT(T_K):
    Vboiling = 1 / (const.rho_b * const.Pa2kPa)
    Vt = Vboiling * np.exp(const.alph * (T_K - const.Tb_K))
    rho1 = 1 / (Vt * const.M_CO2) * const.kPa2Pa # mmol/cm3
    rho = const.rho_b / (const.M_CO2 * np.exp(const.alph * (T_K - const.Tb_K)))
    return rho

def getAandP0dsorptionPotentialForT(T_K, p):
    P0 = const.get_P0(T_K) * const.kPa2Pa
    A = const.R * T_K * np.log(P0 / p)
    return A, P0

def build_DAparams(params):
    """
    构建 lmfit Parameters
    """
    # params = Parameters()
    # 添加参数
    params.add("W0", value=0.349, min=0, max=5)
    params.add("E", value=10209, min=1000, max=50000)
    params.add("n", value=1.92, min=0.2, max=6)
        
    # 保存单位
    params["W0"].user_data = {"unit": "cm3/g"}
    params["E"].user_data = {"unit": "J/mol"}
    params["n"].user_data = {"unit": "-"}

    return params

def getParaForDA(params):
    W0 = params["W0"].value
    E = params["E"].value
    n = params["n"].value

    return W0, E, n

def DA_model0( W0, A, E, n):
    W =  W0 * np.exp( - (A / E)**n )
    return W

def DA_model(rho, W0, A, E, n):
    W = DA_model0( W0, A, E, n)
    q = rho * W
    return q, W

def residual_DA_single(params, data, T):

    # 读取 P (Pa) 和实验 q (mmol/g)
    P = data["P_Pa"].values
    q_exp = data["q_mmol_g"].values

    W0 = params["W0"].value
    E = params["E"].value
    n = params["n"].value

    rho = getCO2LiquidPhaseDensityForT(T)
    A, _ = getAandP0dsorptionPotentialForT(T, P)

    q_fit, _ = DA_model(rho, W0, A, E, n)

    return q_exp - q_fit  # lmfit 要求返回一维数组

def residual_DA(params, grouped_data):
    """
    Residual function for D-A global fitting across multiple temperatures.
    grouped_data: dict[T_K] = df(P, q)
    """
    # 取 D-A 参数
    W0 = params["W0"].value
    E  = params["E"].value
    n  = params["n"].value

    resid_list = []
    # --- 循环每个温度 ---
    for T_K, gdf in grouped_data.items():

        # 读取 P (Pa) 和实验 q (mmol/g)
        P = gdf["P_Pa"].values
        q_obs = gdf["q_mmol_g"].values

        # Adsorbed-phase density (mmol/cm3)
        rho = getCO2LiquidPhaseDensityForT(T_K)

        # Adsorption potential A (J/mol)
        A, P0 = getAandP0dsorptionPotentialForT(T_K, P)

        # 计算 D-A 拟合值
        q_fit, _ = DA_model(rho, W0, A, E, n)

        # 记录残差
        resid_list.append(q_obs - q_fit)
        #加权版本
        # resid_list.append((q_obs - q_fit)/np.max(q_obs))

    # lmfit 要求一维残差向量
    return np.concatenate(resid_list)

def fit_DA01(T_list, grouped_data, sheet=None, file=None):
    """
    对每个温度独立拟合 D-A 模型（W0, E, n 都单独拟合）。
    输入：
        - T_list: 温度列表
        - grouped_data: dict[T_K] = DataFrame(P_Pa, q_mmol_g)
    输出：
        - dict[T_K] = MinimizerResult 对象
        - fits_data DataFrame，包含所有温度的拟合结果
    """
    params0 = Parameters()
    validSheet = False

    all_DA_results = {}
    fits_rows = []

    for T in T_list:
        gdf = grouped_data[T]

        # 构建初始拟合参数

        if sheet: ##origin fitting data
            params = build_DAparams(params0)
            #params, validSheet = get_DAparams_from_origin(params1, T_list, sheet, file) #未开发
        else: 
            params = build_DAparams(params0)

        # 读取 P (Pa) 和实验 q (mmol/g)
        P = gdf["P_Pa"].values
        q_exp = gdf["q_mmol_g"].values

        # Adsorbed-phase density (mmol/cm3) & Adsorption potential (J/mol)
        rho = getCO2LiquidPhaseDensityForT(T)
        A, P0 = getAandP0dsorptionPotentialForT(T, P)

        # 拟合
        daPara = minimize(residual_DA_single, params, args=(gdf, T), method="leastsq", max_nfev=20000)
        all_DA_results[T] = daPara
        print(f"Temperature {T} K fit:")
        print(report_fit(daPara))

        # 构建拟合曲线
        W0, E, n = getParaForDA(daPara.params)
        q_fit, W = DA_model(rho, W0, A, E, n)

        # 计算拟合指标
        residuals, relative_residuals, r2, rmsd = mf.compute_fit_metrics(
            fit_x=P, fit_y=q_fit, exp_x=P, exp_y=q_exp
        )

        # 保存结果
        fits_rows.append(pd.DataFrame({
            "T_K": T,
            "P_Pa": P,
            "q_exp_mmol_g": q_exp,
            "q_fit_mmol_g": q_fit,
            "rho_mmol_cm3": rho,
            "W0_cm3_g": W0,
            "W_cm3_g": W,
            "E_J_mol": E,
            "n": n,
            "A_J_mol": A,
            "P0_kPa": P0 * const.Pa2kPa,
            "R2": r2,
            "residuals": residuals,
            "relative_residuals": relative_residuals,
            "RMSD": rmsd
        }))

    
    fits_data = pd.concat(fits_rows, ignore_index=True)


    return all_DA_results, fits_data, validSheet

def fit_DA0(T_list,  grouped_data, sheet = None, file = None):
        # 构建初始拟合参数，全局拟合策略
    params0 = Parameters()
    validSheet = False
    if sheet: ##origin fitting data
        params = build_DAparams(params0)
        #params, validSheet = get_DAparams_from_origin(params1, T_list, sheet, file) #未开发
    else: 
        params = build_DAparams(params0)

    # 拟合 D-A
    daPara = minimize(residual_DA, params, args=(grouped_data,), method="leastsq", max_nfev=20000)
    print(report_fit(daPara))

    fits_rows = []
    all_P = []
    all_q_exp = []
    all_q_fit = []

    for T in T_list:
        gdf = grouped_data[T]

        # 获取拟合参数（保持原名，单位保存在属性里）
        P = gdf["P_Pa"].values
        q_exp = gdf["q_mmol_g"].values
        rho = getCO2LiquidPhaseDensityForT(T)
        A, P0 = getAandP0dsorptionPotentialForT(T, P)
        W0, E, n = getParaForDA(daPara.params)

        # 拟合曲线
        q_fit, W = DA_model(rho, W0, A, E, n)

         # 计算 R², RMSD
        residuals, relative_residuals, r2, rmsd = mf.compute_fit_metrics(
            fit_x=P, fit_y=q_fit, exp_x=P, exp_y=q_exp
            )

        all_P.append(P)
        all_q_exp.append(q_exp)
        all_q_fit.append(q_fit)

        # 构建 DataFrame，列名保持原始名称
        fits_rows.append(pd.DataFrame({
            "T_K": T,
            "P_Pa": P,
            "q_exp_mmol_g": gdf["q_mmol_g"].values,
            "q_fit_mmol_g": q_fit,
            "rho_mmol_cm3": rho,
            "W0_cm3_g": W0,
            "W_cm3_g": W,
            "E_J_mol": E,
            "n": n,
            "A_J_mol": A,
            "P0_kPa": P0 * const.Pa2kPa,
            "R2": r2,
            "residuals": residuals,
            "relative_residuals": relative_residuals,
            "RMSD": rmsd
        }))

    P_all = np.concatenate(all_P)
    q_exp_all = np.concatenate(all_q_exp)
    q_fit_all = np.concatenate(all_q_fit)

    global_residuals, global_relative_residuals, global_r2, global_rmsd = mf.compute_fit_metrics(
            fit_x=P_all,
            fit_y=q_fit_all,
            exp_x=P_all,
            exp_y=q_exp_all,
            use_interp=False   # 全局
        )

    fits_data = pd.concat(fits_rows, ignore_index=True)
    
    fits_data["global_R2"] = global_r2
    fits_data["global_RMSD"] = global_rmsd
    fits_data["global_residuals"] = global_residuals
    fits_data["global_relative_residuals"] = global_relative_residuals

    return daPara, fits_data, validSheet

def fit_DA(rawData, out_path = None, sheet_name=None):
    """
    拟合 D-A 模型（Dubinin–Astakhov）到原始数据
    输入 rawData：已转换为基础单位的 DataFrame，包含列 [P_Pa, q_mmol_g, T_K]
    输出：
        - D-A 拟合结果（MinimizerResult）
        - fits_data DataFrame，包含列 ["T_K","P_Pa","q_exp","q_fit","W0","E","n","A","rho","alph","Vb","Tb"],
        - q1","q2" unit: mmol_g
        - "b1","b2" unit: 1/Pa
          单位信息保存在 dsl.params[param].unit 中
    """
    # 按温度分组
    grouped_data = {T: gdf for T, gdf in rawData.groupby("T_K")}
    T_list = sorted(grouped_data.keys())

     #全局拟合
    DAPara, fits_data, validSheet = fit_DA0(T_list,  grouped_data, sheet = sheet_name, file = out_path)

    #单独拟合每个温度
    DAPara_single, fits_data_single, validSheet_single = fit_DA01(T_list,  grouped_data, sheet = sheet_name, file = out_path)

    return DAPara, fits_data, DAPara_single, fits_data_single

def runDAsimu(listData, out_path, plotFlag=True, sheet_name = None):
    # 宽表转长表
    rawData0= dop.reshape_wide_to_long(listData)
    # 单位转换
    rawData = dop.convert_units(rawData0)
    #拟合D-A
    DAPara, fits_data, DAPara_single, fits_data_single = fit_DA(rawData, out_path, sheet_name)

    newPath = fl.check_or_create_folder(out_path, sub_folder_name=const.DA_FOLDER)
    newFilePath = fl.get_expanded_name(newPath, sheet_name, expand=const.DA_FILE)
    export_DA_params_custom(DAPara, file_name=newFilePath, sheet_name=const.DA_PARA_SHEET, singleFlag=False)
    fl.export_to_excel_auto(fits_data, newFilePath, sheet_name=const.DA_FITS_SHEET)

    export_DA_params_custom(DAPara_single, file_name=newFilePath, sheet_name=const.DA_BYT_PARA_SHEET, singleFlag=True)
    fl.export_to_excel_auto(fits_data_single, newFilePath, sheet_name=const.DA_BYT_FITS_SHEET)

    if plotFlag:
        poltDA(fits_data, rawData)
        poltDA(fits_data_single, rawData, fixedFlag = True)

# ==========================================
# 主函数
# ==========================================
def main(file_path=None):

    singleFile = True

    if singleFile:
        # # 选择文件
        sheet_name = "CC-Hy-550_60_5-650_15_5-1"
        file = fl.getFile()
        listData, out_path, _ = fl.readFileBySheet(file, sheet_name)
        runDAsimu(listData, out_path, plotFlag=True, sheet_name = sheet_name)

    else:
        file = fl.getFile()
        xls = pd.ExcelFile(file)
        for sheet_name in xls.sheet_names:
            if fl.should_skip(sheet_name): 
                continue

            listData, out_path, _  = fl.readFileBySheet(file, sheet_name)

            if isinstance(listData, dict):
                listData0 = pd.DataFrame(list(listData.values())[0])  # 取第一个 sheet

            if listData0.empty:
                print("DataFrame is empty")
                continue

            runDAsimu(listData, out_path, plotFlag=False, sheet_name = sheet_name)

if __name__ == "__main__":
    f = sys.argv[1] if len(sys.argv) > 1 else None
    main(f)