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

import T0isothermsAndIAST.fit_DSL_simul as dsl
import T0isothermsAndIAST.strongDSL as sl


# -----------------------------
# 定义函数
# -----------------------------
def save_selectivity_to_excel(select, select1, out_path):
    """
    将两个 DataFrame 保存到同一个 Excel 文件的不同 sheet。
    文件名固定为 selectivity.xlsx，保存到 folder_path。
    """
    fl.export_to_excel_auto(select, out_path, sheet_name="P-S")
    fl.export_to_excel_auto(select1, out_path, sheet_name="y-S")

def selec(q, y):
    sele = (q[0] / y[0]) / (q[1] / y[1])
    return sele
    
def getMaxP(df_CO2, df_N2, y_co2):
    Pmax_CO2 = df_CO2["P_bar"].max()
    Pmax_N2 = df_N2["P_bar"].max()

    total_pressure_co2 = Pmax_CO2 / y_co2
    total_pressure_n2 = Pmax_N2 / (1 - y_co2)
    total_pressure_max = min(total_pressure_n2, total_pressure_co2)
    return total_pressure_max

def cal_henry_from_slope(isotherm, P_max):
    """
    :param isotherm: Description
    :param P_max: better less 0.05 bar
    """
    iso_low = isotherm[isotherm["P_bar"] < P_max]
    if len(iso_low) < 3:
        raise ValueError("Not enough low-pressure points for Henry fitting.")
    K_H_slope, intercept = np.polyfit(iso_low["P_bar"], iso_low["q_mmol_g"], 1)
    print(f"K_H = {K_H_slope:.6f}")
    print(f"Intercept = {intercept:.6e}")
    return K_H_slope, intercept

def calHenrySecletivity(df_CO2, df_N2, co2_isotherm, n2_isotherm):
    """
    Calculate CO2/N2 Henry selectivity from fitted isotherm parameters
    (P -> 0 limit).
    Validation logic:
    - Henry constant from model:
        CO2 (DSL): K_H = M1*K1 + M2*K2
        N2  (Langmuir): K_H = M*K
    - Henry constant from data:
        Linear slope of q vs P in low-pressure region (P < P_max)

    Relative deviation definition:
        deviation = |K_model - K_slope| / K_slope * 100%
    
    Evaluation criteria:
        | Relative deviation | Interpretation                     |
        | ------------------ | ---------------------------------- |
        | < 5%               | Excellent agreement                |
        | 5 - 10%            | Acceptable                         |
        | 10 - 15%           | Needs explanation                  |
        | > 15%              | Henry region or model questionable |
    """

    # numerical Henry (slope)
    P_max = 0.05
    co2_slope_low, co2_intercept = cal_henry_from_slope(df_CO2, P_max)
    n2_slope_low, n2_intercept = cal_henry_from_slope(df_N2, P_max)
    sele_henry_numerical = co2_slope_low / n2_slope_low
    
    # model Henry (P → 0)
    co2_qA = co2_isotherm.params["M1"]
    co2_qB = co2_isotherm.params["M2"]
    co2_bA = co2_isotherm.params["K1"]
    co2_bB = co2_isotherm.params["K2"]
    n2_q = n2_isotherm.params["M"]
    n2_b = n2_isotherm.params["K"]

    co2_henry = co2_qA * co2_bA + co2_qB * co2_bB
    n2_henry = n2_q * n2_b
    sele_henry_model = co2_henry / n2_henry

    s_rel_error = abs(sele_henry_model - sele_henry_numerical) / sele_henry_numerical
    co2_rel_error = abs(co2_henry - co2_slope_low) / co2_slope_low
    n2_rel_error = abs(n2_henry - n2_slope_low) / n2_slope_low

    select_henry = pd.DataFrame({
        "sele_henry_model": sele_henry_model,
        "sele_henry_numerical": sele_henry_numerical,
        "co2_henry_model": co2_henry,
        "n2_henry_model": n2_henry,
        "co2_henry_slope": co2_slope_low,
        "co2_henry_intercept": co2_intercept,
        "n2_henry_slope": n2_slope_low,
        "n2_henry_intercept": n2_intercept,
        "s_rel_error": s_rel_error,
        "co2_rel_error": co2_rel_error,
        "n2_rel_error": n2_rel_error
    }, index=[0])
    return select_henry


def calSecletivityAtFixedY(co2_isotherm, n2_isotherm, total_pressure_max, y_co2, p_num):

    p_arange = np.linspace(10e-4, total_pressure_max, p_num)
    y = np.array([y_co2, 1-y_co2])  # gas mole fractions

    q_fixY = []
    sel_fixY = []
    x_co2 = []

    for p_here in p_arange:
        q = pyiast.iast(p_here * y, [co2_isotherm, n2_isotherm], verboseflag=True)
        sel = selec(q, y)
        x = q[0] / (q[0] + q[1]) * 100
        q_fixY.append(q)
        sel_fixY.append(sel)
        x_co2.append(x)

    q_co2 = [q[0] for q in q_fixY]
    q_n2 = [q[1] for q in q_fixY]
    select = pd.DataFrame({
        "P (kPa)": p_arange * 1e2,
        "y(CO2)": y_co2,
        "Uptake": [c + n for c, n in zip(q_co2, q_n2)],
        "x(CO2)": x_co2,
        "x(N2)": [100.0 - xi for xi in x_co2],
        "CO2_uptake（mmol/g）": q_co2,
        "N2_uptake（mmol/g）": q_n2,
        "CO2/N2 selectivity": sel_fixY
    })
    return select

def calSecletivityAtFixedP(co2_isotherm, n2_isotherm, p, y_max, y_num):
    y_arange = np.linspace(0.001, y_max, y_num)
    
    q_fixP = []
    sel_fixP = []
    x_co2 = []

    for y_here in y_arange:
        y = np.array([y_here, 1-y_here])  # gas mole fractions
        q = pyiast.iast(p * y, [co2_isotherm, n2_isotherm], verboseflag=True)
        sel = selec(q, y)
        x = q[0] / (q[0] + q[1]) * 100
        q_fixP.append(q)
        sel_fixP.append(sel)
        x_co2.append(x)

    q_co2 = [q[0] for q in q_fixP]
    q_n2 = [q[1] for q in q_fixP]
    select = pd.DataFrame({
        "P (kPa)": p * 1e2,
        "y(CO2)": y_arange,
        "Uptake": [c + n for c, n in zip(q_co2, q_n2)],
        "x(CO2)": x_co2,
        "x(N2)": [100.0 - xi for xi in x_co2],
        "CO2_uptake（mmol/g）": q_co2,
        "N2_uptake（mmol/g）": q_n2,
        "CO2/N2 selectivity": sel_fixP
    })
    return select

def poltP2Sele(select):
    plt.figure(figsize=(6,4))
    plt.plot(select["P (kPa)"], select["CO2/N2 selectivity"], marker="o", label="Selectivity")

    plt.xlabel("Pressure (kPa)")
    plt.ylabel("Selectivity (-)")
    plt.title("P–Selectivity Curve")

    # 在图上写 y_co2
    y_value = select["y(CO2)"].iloc[0]   # 取第一个（如果是常数列）
    plt.text(
        x=0.7*select["P (kPa)"].max(),  # 横坐标：选在右上角
        y=0.9*select["CO2/N2 selectivity"].max(),  # 纵坐标：选在上方
        s=f"y_CO2 = {y_value:.2f}", 
        fontsize=10,
        bbox=dict(facecolor="white", alpha=0.6, edgecolor="gray")
    )
    plt.legend()
    plt.grid(True)  
    plt.show(block=False)

def poltY2Sele(select):
    plt.figure(figsize=(6,4))
    plt.plot(select["y(CO2)"], select["CO2/N2 selectivity"], marker="o", label="Selectivity")
    plt.xlabel("y_co2 (-)")
    plt.ylabel("Selectivity (-)")
    plt.title("y_co2 – Selectivity Curve")
    # 在图上写 pressure
    p_value = select["P (kPa)"].iloc[0]   # 取第一个（如果是常数列）
    plt.text(
        x=0.7*select["y(CO2)"].max(),  # 横坐标：选在右上角
        y=0.9*select["CO2/N2 selectivity"].max(),  # 纵坐标：选在上方
        s=f"P = {p_value:.2f} kPa",          # 带单位
        fontsize=10,
        bbox=dict(facecolor="white", alpha=0.6, edgecolor="gray")
    )
    plt.legend()
    plt.grid(True)
    plt.show(block=False)

def poltY2X(select):
    plt.figure(figsize=(6,4))
    plt.plot(select["y(CO2)"], select["x(CO2)"], marker="o", label="Purity")
    plt.xlabel("y_co2 (-)")
    plt.ylabel("x_CO2 (%)")
    plt.title("y_co2 – x_CO2 Curve")
    # 在图上写 pressure
    p_value = select["P (kPa)"].iloc[0]   # 取第一个（如果是常数列）
    plt.text(
        x=0.7*select["y(CO2)"].max(),  # 横坐标：选在右上角
        y=0.9*select["x(CO2)"].max(),  # 纵坐标：选在上方
        s=f"P = {p_value:.2f} kPa",          # 带单位
        fontsize=10,
        bbox=dict(facecolor="white", alpha=0.6, edgecolor="gray")
    )
    # 添加垂直虚线 x=0.04 和 x=0.15
    plt.axvline(x=0.04, color='red', linestyle='--', label='x=0.04')
    plt.axvline(x=0.15, color='green', linestyle='--', label='x=0.15')

    # 添加水平虚线 y=90
    plt.axhline(y=90, color='blue', linestyle='--', label='y=90')
    plt.legend()
    plt.grid(True)
    plt.show(block=True)

def outputTxt(T, list, path, fileName):

    CO2_PName = dop.getName(T, unitChange = 0, prefix=None, suffix = "Absolute Pressure (kPa)")
    N2_PName = dop.getName(T, unitChange = 0, prefix="N2_", suffix = "Absolute Pressure (kPa)")
    CO2_qName = dop.getName(T, unitChange = 0, prefix=None, suffix = "uptake (mmol/g)")
    N2_qName = dop.getName(T, unitChange = 0, prefix="N2_", suffix = "uptake (mmol/g)")

    co2_25 = dop.cols_to_clean_df(list, CO2_PName, CO2_qName)
    n2_25 = dop.cols_to_clean_df(list, N2_PName, N2_qName)

    newPath = fl.check_or_create_folder(path, sub_folder_name=const.SELE_FOLDER_IAST_SOFT)
    out_path_CO2 = fl.get_expanded_name(newPath, fileName, expand="CO2_25C", expandPos=False, type="txt")
    out_path_N2 = fl.get_expanded_name(newPath, fileName, expand="N2_25C", expandPos=False, type="txt")

    co2_25.to_csv(out_path_CO2, sep=" ", index=False, header=False)
    n2_25.to_csv(out_path_N2, sep=" ", index=False, header=False)

def safe_value(x):
    """从 lmfit Parameter 或数值中安全取值"""
    return x.value if hasattr(x, "value") else float(x)

def saveModelParamForIAST(T_SELE, co2_isotherm, n2_isotherm, out_path = None, sheet_name="Model_params"):
    rows = []

    if n2_isotherm.model == "Henry":
        rows = [
        {
            "gas": "CO2",
            "model": co2_isotherm.model,
            "T(°C)": T_SELE,
            "qA(mmol/g)": safe_value(co2_isotherm.params["M1"]),
            "qB(mmol/g)": safe_value(co2_isotherm.params["M2"]),
            "bA(1/kPa)": safe_value(co2_isotherm.params["K1"]) * const.kPa2bar,
            "bB(1/kPa)": safe_value(co2_isotherm.params["K2"]) * const.kPa2bar,
        },
        {
            "gas": "N2",
            "model": n2_isotherm.model,
            "T(°C)": T_SELE,
            "KH(mmol/g*kPa)": safe_value(n2_isotherm.params["KH"]) * const.kPa2bar,
        }
        ]
    elif n2_isotherm.model == "Langmuir":    
        rows = [
        {
            "gas": "CO2",
            "model": co2_isotherm.model,
            "T(°C)": T_SELE,
            "qA(mmol/g)": safe_value(co2_isotherm.params["M1"]),
            "qB(mmol/g)": safe_value(co2_isotherm.params["M2"]),
            "bA(1/kPa)": safe_value(co2_isotherm.params["K1"]) * const.kPa2bar,
            "bB(1/kPa)": safe_value(co2_isotherm.params["K2"]) * const.kPa2bar,
        },
        {
            "gas": "N2",
            "model": n2_isotherm.model,
            "T(°C)": T_SELE,
            "q(mmol/g)": safe_value(n2_isotherm.params["M"]),
            "b(1/kPa)": safe_value(n2_isotherm.params["K"]) * const.kPa2bar,
        }
        ]

    df_export = pd.DataFrame(rows)
    
    # 使用已有自动导出函数
    fl.export_to_excel_auto(df_export, out_path, sheet_name=sheet_name)

def getParaByUserDefine(listData, co2_isotherm, n2_isotherm, out_path = None, sheet_name=None):

    co2_dsl = copy.deepcopy(co2_isotherm)
    n2_sl = copy.deepcopy(n2_isotherm)

    rawData0= dop.reshape_wide_to_long(listData)
        # 单位转换
    rawData = dop.convert_units(rawData0)
        #拟合DSL
    co2_dsl0, _ = dsl.fit_dsl(rawData, out_path, sheet_name)

    grouped_data = {T: gdf for T, gdf in rawData.groupby("T_K")}
    T_max = max(grouped_data.keys())
    b1_name = dsl.pname("b1", T_max)
    b2_name = dsl.pname("b2", T_max)

    from lmfit import Parameter
    co2_dsl.params["M1"] = Parameter(name="M1", value=co2_dsl0.params["q1"].value, min=0)
    co2_dsl.params["M2"] = Parameter(name="M2", value=co2_dsl0.params["q2"].value, min=0)
    co2_dsl.params["K1"] = Parameter(name="K1", value=co2_dsl0.params[b1_name].value * const.bar2Pa, min=0)
    co2_dsl.params["K2"] = Parameter(name="K2", value=co2_dsl0.params[b2_name].value * const.bar2Pa, min=0)

    T_C = T_max + const.K2C
    N2_PName = dop.getName(T_C, unitChange = 0, prefix="N2_", suffix = "Absolute Pressure (kPa)")
    N2_qName = dop.getName(T_C, unitChange = 0, prefix="N2_", suffix = "uptake (mmol/g)")
    list_N2  = dop.extract_PairData(N2_PName, N2_qName, listData, x_Scale = const.kPa2Pa)
    df_N2 = pd.DataFrame(list_N2, columns=["P_Pa", "q_mmol_g"])

    if n2_isotherm.model == "Langmuir": 
        n2_sl0 = sl.fit_single_langmuir(df_N2)
        
        n2_sl.params["M"] = Parameter(name="M", value=n2_sl0.params["q_sat"].value, min=0)
        n2_sl.params["K"] = Parameter(name="K", value=n2_sl0.params["b"].value * const.bar2Pa, min=0)
    elif n2_isotherm.model == "Henry":
        n2_sl0 = sl.fit_henry_model(df_N2)
        n2_sl.params["KH"] = Parameter(name="KH", value=n2_sl0.params["KH"].value * const.bar2Pa, min=0)
    return co2_dsl, n2_sl

def getDataForIAST(listData, T = 25):
    """
    支持 DataFrame 或字典(list) 形式输入，返回 CO2 与 N2 数据的 numpy 数组。
    输出：
        df_CO2, df_N2 : 均为二维 numpy 数组 [Pressure_bar, Uptake_mmol_g]
    """
    CO2_PName = dop.getName(T, unitChange = 0, prefix=None, suffix = "Absolute Pressure (kPa)")
    N2_PName = dop.getName(T, unitChange = 0, prefix="N2_", suffix = "Absolute Pressure (kPa)")
    CO2_qName = dop.getName(T, unitChange = 0, prefix=None, suffix = "uptake (mmol/g)")
    N2_qName = dop.getName(T, unitChange = 0, prefix="N2_", suffix = "uptake (mmol/g)")

    list_CO2 = dop.extract_PairData(CO2_PName, CO2_qName, listData, x_Scale = const.kPa2bar)
    list_N2  = dop.extract_PairData(N2_PName, N2_qName, listData, x_Scale = const.kPa2bar)
    
    df_CO2 = pd.DataFrame(list_CO2, columns=["P_bar", "q_mmol_g"])
    df_N2 = pd.DataFrame(list_N2, columns=["P_bar", "q_mmol_g"])
    return df_CO2, df_N2

def build_co2_iast_shell():
    # 构造一个“必定可拟合”的假 DSL 数据
    fake_df = pd.DataFrame({
        "P_bar": [0.01, 0.1, 1.0],
        "q_mmol_g": [0.01, 0.02, 0.03]
    })

    iso = pyiast.ModelIsotherm(
        fake_df,
        pressure_key="P_bar",
        loading_key="q_mmol_g",
        model="Langmuir"   # 注意：先用 Langmuir，最稳
    )

    # 标记为 DSL（供你后续逻辑判断）
    iso.model = "DSLangmuir"

    # 清空参数（你后面会全部覆盖）
    iso.params = {}

    return iso

def build_n2_iast_shell(model="Langmuir"):
    fake_df = pd.DataFrame({
        "P_bar": [0.01, 0.1, 1.0],
        "q_mmol_g": [0.001, 0.002, 0.003]
    })

    iso = pyiast.ModelIsotherm(
        fake_df,
        pressure_key="P_bar",
        loading_key="q_mmol_g",
        model="Langmuir"
    )

    iso.model = model
    iso.params = {}

    return iso

def build_co2_isotherm_safe(df_CO2):
    try:
        co2_isotherm = pyiast.ModelIsotherm(
            df_CO2,
            loading_key="q_mmol_g",
            pressure_key="P_bar",
            model="DSLangmuir",
            optimization_method="trust-constr" #trust-constr // CG
        )
        return co2_isotherm

    except Exception as e:
        print("[WARNING] CO2 DSL fitting failed, fallback to IAST shell.")
        print("          Reason:", e)

        return build_co2_iast_shell()
    
def build_n2_isotherm_safe(df_N2):
    try:
        n2_isotherm = pyiast.ModelIsotherm(
            df_N2,
            loading_key="q_mmol_g",
            pressure_key="P_bar",
            model="Langmuir",
            optimization_method="CG"
        )
        return n2_isotherm

    except Exception as e:
        print("[WARNING] N2 Langmuir fitting failed, fallback to IAST shell.")
        print("          Reason:", e)

        return build_n2_iast_shell(model="Langmuir")

def runPyIast(listData, T_C = 25, plotFlag = False, out_path = None, sheet_name=None):
    if listData is not None:
        # 创建新的 listData，只包含 25°C 的压力和吸附量
        df_CO2, df_N2 = getDataForIAST(listData, T_C)
    
    co2_isotherm = build_co2_isotherm_safe(df_CO2)
    n2_isotherm  = build_n2_isotherm_safe(df_N2)

    co2_dsl, n2_sl = getParaByUserDefine(listData, co2_isotherm, n2_isotherm, out_path, sheet_name)

    sele_henry = calHenrySecletivity(df_CO2, df_N2, co2_dsl, n2_sl)
    
    y_co2 = 0.15 #co2 组分比例
    p_num = 100
    total_pressure_max = getMaxP(df_CO2, df_N2, y_co2)
    P_max = max(total_pressure_max, 1)  #bar = 1
    select = calSecletivityAtFixedY(co2_dsl, n2_sl, P_max, y_co2, p_num)

    p = 1  #bar = 1
    y_max = 0.2
    y_num = 100
    select1 = calSecletivityAtFixedP(co2_dsl, n2_sl, p, y_max, y_num)

    newPath = fl.check_or_create_folder(out_path, sub_folder_name=const.SELE_FOLDER_IAST_PY)
    newOutPath = fl.get_expanded_name(newPath, fileName = sheet_name, expand=const.SELE_FILE_IAST_PY)

    saveModelParamForIAST(T_C, co2_dsl, n2_sl, out_path = newOutPath)
    save_selectivity_to_excel(select, select1, newOutPath)
    fl.export_to_excel_auto(sele_henry, filename = newOutPath, sheet_name="S-HENRY")
    if plotFlag:
        poltP2Sele(select)
        poltY2Sele(select1)
        poltY2X(select1)
# -----------------------------
# main函数*************************************************************************************************************
# -----------------------------
def main():
    singleFile = True
    T_C = 25

    if singleFile:
        sheet = "CC-Hy-550_60_5-650_15_5-1"
        file = fl.getFile()
        listData, out_path, _ = fl.readFileBySheet(file, sheet, expand = const.SELE_FILE_IAST_SOFT)
        outputTxt(T_C, listData, out_path, fileName = sheet)
        runPyIast(listData, T_C, plotFlag=True, out_path = out_path, sheet_name=sheet)

    else:
        file = fl.getFile()
        xls = pd.ExcelFile(file)
        for sheet in xls.sheet_names:
            if fl.should_skip(sheet): 
                continue

            listData, out_path, _ = fl.readFileBySheet(file, sheet, expand = const.SELE_FILE_IAST_SOFT)

            if isinstance(listData, dict):
                listData0 = pd.DataFrame(list(listData.values())[0])  # 取第一个 sheet

            if listData0.empty:
                print("DataFrame is empty")
                continue

            outputTxt(T_C, listData, out_path, fileName = sheet)
            runPyIast(listData, T_C, plotFlag=False, out_path = out_path, sheet_name=sheet)
    
# -----------------------------
# Python入口
# -----------------------------
if __name__ == "__main__":
    main()