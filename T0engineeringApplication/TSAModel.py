import os
import sys
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from common import *
import constantsAndName as const

import T1fileSource.fileOperation as fl
import T1dataProcessSource.dataOperation as dop

import T0isothermsAndIAST.fit_DSL_simul as DSL

# ----------------------------
# 通用 q_ads 函数
# SL / DSL 通用
# 高温 ΔH 修正可开关 (apply_dH_correction=True/False)
# T ≤ 60℃ 时不修正 ΔH
# 支持 拟合温度列表 计算参考温度 T_low
# 向量化压力 / 温度，适合 TSA 循环模拟
# 单次 TSA：apply_cycle=False，只算一次 q_ads、q_des、q_working、Qsorption 等
# 循环 TSA：apply_cycle=True，n_cycles=3，输出每个循环结果
# ΔH 修正开关：apply_dH_correction=True
# 使用 T_fit_SL 计算参考温度 T_low 进行 ΔH 高温修正： Kirchhoff 方程更为完美，当前只使用线性近似，是 Kirchhoff 方程的简化版
# ----------------------------


def K_from_HS(T, DeltaH, DeltaS):    #"""    计算无量纲的平衡常数:    """
    return np.exp((DeltaS / const.R) - (DeltaH / (const.R * T)))

def DeltaH_correct(DeltaH, T, T_low=298.15, k=0.0):     #温度超 60C 修正
    T_C = T - 273.15
    return np.where(T_C > 60.0, DeltaH - k*(T - T_low), DeltaH)

def q_model(P, T, params, apply_dH_correction=False, kA=0.0, kB=0.0, T_fit=None):
    """
    计算平衡吸附量 q(P, T)。
    参数:
    - P: 压力 (kPa)
    - T: 温度 (K)
    - params:
        SL: {"q": q_sat (mol/g), "dH": ΔH (J/mol), "dS": ΔS (J/mol/K)}
        DSL: {"qA": , "dHA": , "dSA": ,
              "qB": , "dHB": , "dSB": }  (同上, qA/qB 单位 mol/g)
    - apply_dH_correction: 是否应用 ΔH 修正 (True/False)
    - kA, kB: ΔCp 参数 (J/mol/K)，Kirchhoff 方程修正用; 具体值还需要确认
    - T_fit: 拟合数据的参考温度 (K)，修正时需要, 具体值对应拟合DSL时吸附曲线的具体温度

    返回:
    - q (mol/g)
    """
    P0 = 100   #kPa:100   Pa:1e5  bar: 1
    P = np.atleast_1d(P).astype(float)
    T = np.atleast_1d(T).astype(float)
    T_low = np.mean(T_fit) if T_fit is not None else 298.15

    # 单位点 Langmuir
    if "q" in params:
        DeltaH_eff = DeltaH_correct(params["dH"], T, T_low, kA) if apply_dH_correction else params["dH"]
        K = K_from_HS(T, DeltaH_eff, params["dS"])
        alpha = K / P0
        return float(params["q"] * alpha * P / (1.0 + alpha * P))
    
    # 双位点 Langmuir
    else:
        DeltaH_A = DeltaH_correct(params["dHA"], T, T_low, kA) if apply_dH_correction else params["dHA"]
        DeltaH_B = DeltaH_correct(params["dHB"], T, T_low, kB) if apply_dH_correction else params["dHB"]
        KA = K_from_HS(T, DeltaH_A, params["dSA"])
        KB = K_from_HS(T, DeltaH_B, params["dSB"])
        bA = KA / P0
        bB = KB / P0
        qA = params["qA"] * bA * P / (1.0 + bA * P)
        qB = params["qB"] * bB * P / (1.0 + bB * P)
        return float(qA + qB)

# ----------------------------
# 单次 TSA 模拟（支持 Qst）
# ----------------------------
def TSA_single(P_ads, P_des, T_ads, T_des, params, m=1.0, cp=1.0,
               apply_dH_correction=False, kA=0.0, kB=0.0, T_fit=None,
               n_points=100, Qst_func=None):
    """
    参数:
    - P_ads: 吸附压力 (kPa)
    - P_des: 解吸压力 (kPa)
    - T_ads: 吸附温度 (K)
    - T_des: 解吸温度 (K)
    - params: 吸附模型参数字典 (SL/DSL)
    - m: 吸附剂质量 (g)
    - cp: 比热容 (J/g*K)
    - apply_dH_correction: 是否修正 ΔH
    - kA, kB: ΔCp 参数 (J/mol/K)，Kirchhoff 修正用
    - T_fit: 拟合数据对应温度 (K)
    - n_points: 积分点数
    - Qst_func: 可选，返回 (-Qst(q)) 的函数 (J/mol)，随 q 变化

    返回:
    - dict，包含 q_ads, q_des, q_working, Qsorption, Qtemp, Qregen
      以及 P_ads, P_des, T_ads, T_des 和单位
    """
    results = []

    q_ads_val = q_model(P_ads, T_ads, params, apply_dH_correction, kA, kB, T_fit)
    q_des_val = q_model(P_des, T_des, params, apply_dH_correction, kA, kB, T_fit)
    q_working = q_ads_val - q_des_val

    if "q" in params:
        dH_avg = params["dH"]
    else:
        dH_avg = (params["qA"]*params["dHA"] + params["qB"]*params["dHB"]) / (params["qA"]+params["qB"])

    # ----------------------------
    # Qsorption 计算
    # ----------------------------
    q_vals = np.linspace(q_des_val, q_ads_val, n_points).ravel()
    if Qst_func is not None:
        Qst_vals = np.array(Qst_func(q_vals)).ravel()  # 输入 Qst(q) 数组  # J/mol
        Qsorption = -np.trapezoid(Qst_vals, q_vals) / 1000  # kJ/g， 梯形积分
        # Qsorption00 = (-dH_avg / 1000) * q_working 
        # print("Qst_vals shape:", np.shape(Qst_vals))
        # print("Qsorption raw:", Qsorption)
        # print("Qsorption00:", Qsorption00)
        # print("Qst_vals:", Qst_vals)
        # print("q_vals:", q_vals)
    else:
        Qsorption = (-dH_avg / 1000) * q_working  # kJ/g 平均 ΔH 近似

    # 温度贡献
    Qtemp = m * cp * (T_des - T_ads) / 1000  # kJ
    if m != 1.0:
        Qsorption_total = Qsorption * m       # kJ
        Qregen = Qsorption_total + Qtemp
    else:
        Qregen = Qsorption + Qtemp            # kJ/g, 默认 m=1g

    T_ads_C = dop.kelvin_to_celsius(T_ads, unit="K")
    T_des_C = dop.kelvin_to_celsius(T_des, unit="K")
    # return pd.DataFrame([{
    #     "q_ads (mol/g)": q_ads_val,
    #     "q_des (mol/g)": q_des_val,
    #     "q_working (mol/g)": q_working,
    #     "q_working (mmol/g)": q_working * 1000.0,
    #     "Qsorption (kJ/g)": Qsorption,
    #     "Qtemp (kJ/g)": Qtemp,
    #     "Qregen (kJ/g)": Qregen,
    #     "P_ads (kPa)": P_ads,
    #     "P_des (kPa)": P_des,
    #     "T_ads (K)": T_ads,
    #     "T_ads_C (C)": T_ads_C,
    #     "T_des (K)": T_des,
    #     "T_des_C (C)": T_des_C,
    #     "q_working/Qregen (mmol/kJ)": (q_working * 1000.0) / Qregen
    # }])
    row = {
        "q_ads (mol/g)": q_ads_val,
        "q_des (mol/g)": q_des_val,
        "q_working (mol/g)": q_working,
        "q_working (mmol/g)": q_working * 1000.0,
        "Qsorption (kJ/g)": Qsorption,
        "Qtemp (kJ/g)": Qtemp,
        "Qregen (kJ/g)": Qregen,
        "P_ads (kPa)": P_ads,
        "P_des (kPa)": P_des,
        "T_ads (K)": T_ads,
        "T_ads_C (C)": T_ads_C,
        "T_des (K)": T_des,
        "T_des_C (C)": T_des_C,
        "q_working/Qregen (mmol/kJ)": (q_working * 1000.0) / Qregen
    }
    results.append(row)
    return pd.DataFrame(results)


# ----------------------------
# 基于变化 P/T 的循环 TSA
# ----------------------------
def TSA_cycle(P_ads, T_ads, params,
              P_des_range=None, T_des_range=None, num=1, m=1.0, cp=1.0,
              apply_dH_correction=False, kA=0.0, kB=0.0, T_exp=None, Qst_func=None):
    """
    灵活的 TSA 循环模拟函数，支持单变量扫描和笛卡尔积扫描。

    参数:
    - P_ads: 吸附压力 (kPa)，固定标量
    - P_des: 解吸压力 (kPa)，标量或列表/数组
    - T_ads: 吸附温度 (K)，固定标量
    - T_des: 解吸温度 (K)，标量或列表/数组
    - params: 吸附模型参数字典 (SL/DSL)，例如：
        SL: {"q": q_sat (mol/g), "dH": ΔH (J/mol), "dS": ΔS (J/mol/K)}
        DSL: {"qA": , "dHA": , "dSA": ,
              "qB": , "dHB": , "dSB": }
    - P_des_range: 解吸压强 [start, stop] (kPa)，若为 None，则使用 P_ads
    - T_des_range: 解吸温度 [start, stop] (K)，若为 None，则使用 T_ads
    - step: 解吸变量步长 (kPa 或 K)
    - num: 采样点数
    - m: 吸附剂质量 (g)
    - cp: 比热容 (J/g/K)
    - apply_dH_correction: 是否修正 ΔH (True/False)
    - kA, kB: ΔCp 参数 (J/mol/K)，用于 Kirchhoff 修正，默认为0
    - T_exp: 拟合吸附曲线时的温度数组 (K)，用于 ΔH 修正

    返回:
    - results: list，每个元素是 dict，对应每个循环结果，包含：
        "q_ads", "q_des", "q_working", "Qsorption", "Qtemp", "Qregen"
    """
    # 压力列表
    if P_des_range is None:
        P_des_list = [P_ads]
    elif isinstance(P_des_range, (int, float)):  # 单个指定值
        P_des_list = [P_des_range]
    else:  # [start, stop]
        start, stop = P_des_range
        if num == 1:
            P_des_list = [start]
        else:
            P_des_list = list(np.linspace(start, stop, num))

    # 温度列表
    if T_des_range is None:
        T_des_list = [T_ads]
    elif isinstance(T_des_range, (int, float)):  # 单个指定值
        T_des_list = [T_des_range]
    else:  # [start, stop]
        start, stop = T_des_range
        if num == 1:
            T_des_list = [start]
        else:
            T_des_list = list(np.linspace(start, stop, num))
    
    # ----------------------------
    # 笛卡尔积循环
    # ----------------------------
    results = []
    for P_des in P_des_list:
        for T_des in T_des_list:
            res = TSA_single(P_ads, P_des, T_ads, T_des, params, m, cp,
                             apply_dH_correction=apply_dH_correction,
                             kA=kA, kB=kB, T_fit=T_exp, Qst_func=Qst_func)
            results.append(res)
    df_results = pd.concat(results, ignore_index=True)
    return df_results

def runTsaModel(listData,out_path,sheet_name):
    # 宽表转长表
    rawData0=dop.reshape_wide_to_long(listData)
    # 单位转换
    rawData = dop.convert_units(rawData0)
    grouped_data={T: gdf for T,gdf in rawData.groupby("T_K")}
    T_list_K=sorted(grouped_data.keys())
    
    #拟合DSL
    dsl, _ = DSL.fit_dsl(rawData=rawData,out_path = out_path,sheet_name=sheet_name)
    param =dsl.params
    thermo_df, _ = DSL.compute_vant_hoff(param, T_list_K)
    Qst_exp, _ = DSL.getQst(dsl, thermo_df, rawData, 20)  
    Qst_func = interp1d(
        Qst_exp["q_mmol_g"]/1000,   # 转换为 mol/g
        -Qst_exp["Qst_J_mol"],      # 注意 TSA_single 里假设 Qst_func 返回的是 -Qst(q)
        kind="linear",
        bounds_error=False,
        fill_value="extrapolate"
        )
    params_DSL = {
    "qA": param['q1'].value /1000,     # 位点 A 饱和吸附量 (mol/g)
    "dHA": thermo_df.loc[thermo_df["site"] == "b1", "dH_Jmol"].values[0], # 位点 A 吸附焓 ΔH (J/mol)
    "dSA": thermo_df.loc[thermo_df["site"] == "b1", "dS_Jmol_K"].values[0],    # 位点 A 吸附熵 ΔS (J/mol/K)
    "qB": param['q2'].value / 1000,     # 位点 B 饱和吸附量 (mol/g)
    "dHB": thermo_df.loc[thermo_df["site"] == "b2", "dH_Jmol"].values[0], # 位点 B 吸附焓 ΔH (J/mol)
    "dSB": thermo_df.loc[thermo_df["site"] == "b2", "dS_Jmol_K"].values[0]     # 位点 B 吸附熵 ΔS (J/mol/K)
    }

    # T_exp_C = np.array([0, 15, 25]) #吸附曲线温度，单位 K
    # T_exp_K = to_kelvin(T_exp_C, unit="C")

    # 吸附条件固定
    P_ads = 15       # kPa
    T_ads_K = dop.celsius_to_kelvin(25, unit="C")  # 单个值

    # 解吸压力和温度扫描范围
    P_des_range = 100   # kPa,P_des_range = 100 或P_des_range = [100, 200]
    T_des_range_C = [100, 200]  # C -> 转 K
    T_des_range_K = dop.celsius_to_kelvin(T_des_range_C, unit="C")

    cycle_results = TSA_cycle(
        P_ads=P_ads,
        T_ads=T_ads_K,
        params=params_DSL,
        P_des_range=P_des_range,
        T_des_range=T_des_range_K,
        num=10,                 # 例如采样10个点
        apply_dH_correction=False,
        kA=10.0,
        T_exp=T_list_K,
        Qst_func=Qst_func
    )

    return cycle_results
    
# ==========================================
# 主函数
# ==========================================
def main(file_path=None):
    
    # 选择文件
    sheet_name = "CC-Hy-550_60_5-650_15_5-1"
    file = fl.getFile()
    listData, out_path, _ = fl.readFileBySheet(file, sheet_name, expand = const.TSA_FILE)
    Tsa = runTsaModel(listData, out_path, sheet_name)
    fl.export_to_excel_auto(Tsa, filename=out_path)

    # file_path = fl.getFile()
    # xls = pd.ExcelFile(file_path)
    # fileName = os.path.splitext(os.path.basename(file_path))[0]
    # out_path = fl.get_expanded_name(file_path, fileName, expand=const.TSA_FILE, expandPos=True, type="xlsx")
    # for sheet in xls.sheet_names:
    #     if fl.should_skip(sheet): 
    #         continue

    #     listData, _ , _= fl.readFileBySheet(file_path, sheet)

    #     if isinstance(listData, dict):
    #         listData0 = pd.DataFrame(list(listData.values())[0])  # 取第一个 sheet

    #     if listData0.empty:
    #         print("DataFrame is empty")
    #         continue

    #     Tsa = runTsaModel(listData, out_path, sheet_name)
    #     fl.export_to_excel_auto(Tsa, filename=out_path, sheet_name=sheet)
    
if __name__ == "__main__":
    f = sys.argv[1] if len(sys.argv) > 1 else None
    main(f)

# # ----------------------------
# # 示例：SL 循环
# # ----------------------------
# params_SL = {
#     "q": 5.0,      # 吸附剂饱和吸附量 (mol/g)
#     "dH": -25000,  # 吸附焓 ΔH (J/mol)
#     "dS": -50      # 吸附熵 ΔS (J/mol/K)
# }

# # ----------------------------
# # DSL 示例参数
# # ----------------------------
# params_DSL = {
#     "qA": 3.0,     # 位点 A 饱和吸附量 (mol/g)
#     "dHA": -20000, # 位点 A 吸附焓 ΔH (J/mol)
#     "dSA": -40,    # 位点 A 吸附熵 ΔS (J/mol/K)
#     "qB": 2.0,     # 位点 B 饱和吸附量 (mol/g)
#     "dHB": -25000, # 位点 B 吸附焓 ΔH (J/mol)
#     "dSB": -50     # 位点 B 吸附熵 ΔS (J/mol/K)
# }

# T_exp_C = np.array([0, 15, 25]) #吸附曲线温度，单位 K
# T_exp_K = to_kelvin(T_exp_C, unit="C")

# # 吸附条件固定
# P_ads = 50       # kPa
# T_ads_K = to_kelvin(40, unit="C")  # 单个值

# # 解吸压力和温度扫描范围
# P_des_range = 100   # kPa,P_des_range = 100 或P_des_range = [100, 200]
# T_des_range_C = [100, 200]  # C -> 转 K
# T_des_range_K = to_kelvin(T_des_range_C, unit="C")

# cycle_results = TSA_cycle(
#     P_ads=P_ads,
#     T_ads=T_ads_K,
#     params=params_DSL,
#     P_des_range=P_des_range,
#     T_des_range=T_des_range_K,
#     num=10,                 # 例如采样10个点
#     apply_dH_correction=True,
#     kA=10.0,
#     T_exp=T_exp_K
# )

# save_TSA_to_excel(cycle_results)
# # 输出结果
# for idx, r in enumerate(cycle_results):
#     print(f"Cycle {idx+1}:")
#     for k, v in r.items():
#         if isinstance(v, np.ndarray):
#             v_str = np.array2string(v, precision=4, separator=',')  # 数组打印
#         else:
#             v_str = f"{v:.4f}"  # 标量打印
#         print(f"  {k}: {v_str}")

