# -*- coding: utf-8 -*-
"""
Dual-Site Langmuir (DSL) 多温度拟合 + Qst/H/S/G 计算
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

import T0isothermsAndIAST.modelAndFit as mf
import T0isothermsAndIAST.strongDSL as sl
import T0isothermsAndIAST.isostericHeatQstForRawData as qstByRawdata

# ==========================================
# 工具函数
# ==========================================

def export_dsl_params_custom(dsl, file_name="results.xlsx", sheet_name="DSL_params"):
    """
    导出 DSL 拟合参数到 Excel，自定义列名和单位：
    T(K)	bA(1/kPa)	bB(1/kPa)	qA(mmol/g)	qB(mmol/g)
    b1,b2 单位从 1/Pa 转为 1/kPa
    """
    rows = []
    # 遍历所有 b1 参数
    for key in sorted([k for k in dsl.params.keys() if k.startswith("b1")]):
        # 解析温度
        T_str = key.split("_", 1)[1].replace("_", ".")
        T = float(T_str)

        b1 = dsl.params[key].value * 1e3  # 直接用 key，不再重新生成
        b2 = dsl.params[key.replace("b1", "b2", 1)].value * 1e3

        row = {
            "T(K)": T,
            "qA(mmol/g)": dsl.params["q1"].value,
            "qB(mmol/g)": dsl.params["q2"].value,
            "bA(1/kPa)": b1,
            "bB(1/kPa)": b2            
        }
        rows.append(row)
    
    df_export = pd.DataFrame(rows)
    
    # 使用已有自动导出函数
    fl.export_to_excel_auto(df_export, file_name, sheet_name=sheet_name)

def dsl_model(P, q1, b1, q2, b2):
    return q1 * (b1*P)/(1+b1*P) + q2*(b2*P)/(1+b2*P)

def pname(prefix, T):
    if float(T).is_integer():
        return f"{prefix}_{int(T)}"
    else:
        return f"{prefix}_{str(T).replace('.', '_')}"
    
def get_dslparams_from_origin(params, T_list, sheet, file):
    path = fl.getPathOnly(file)  # current folder
    file_name = "DslOriginFit"
    file = fl.get_expanded_name(path, file_name, type="xlsx")
    
    cols, _, validSheet = fl.readFileBySheet(file, sheet, header=[0,1])
    
    if validSheet:
        # 先构建 dict: 行名 -> 列数据索引
        # 这里假设 cols 的所有 list 长度相同
        n_rows = len(next(iter(cols.values())))

        # -----------------------------
        # 遍历温度 T_list，赋值 b1, b2
        # -----------------------------
        for i, T in enumerate(T_list):
            T_label = f"T{int(T - 273.15)}_interp"
            # 找到行索引
            if i >= n_rows:
                raise ValueError(f"Missing data for {T_label} in Excel")
            
            # bA 和 bB
            b1_name = pname("b1", T)
            b2_name = pname("b2", T)
            
            bA_col = ('bA', 'Value')
            bB_col = ('bB', 'Value')
            
            params[b1_name].value = float(cols[bA_col][i])
            params[b2_name].value = float(cols[bB_col][i])
        
        # -----------------------------
        # q1 和 q2 用平均值
        # -----------------------------
        q1_col = ('qA', 'Value')
        q2_col = ('qB', 'Value')
        
        params["q1"].value = float(np.mean(cols[q1_col]))
        params["q2"].value = float(np.mean(cols[q2_col]))
    
    return params, validSheet

def build_paramsBySL(T_list, grouped_data, flagForQ1less=True):
    params = Parameters()
    
    init_params, sl_summary = sl.get_initial_dsl_params_multitemp(T_list, grouped_data)

    if flagForQ1less:
        params.add("q1", value=init_params["q1_sat"], min=0)
        params.add("q2", value=init_params["q2_sat"], min=0)
    else:
        params.add("q1", value=init_params["q2_sat"], min=0)
        params.add("q2", value=init_params["q1_sat"], min=0)
    # 保存单位
    params["q1"].user_data = {"unit": "mmol/g"}
    params["q2"].user_data = {"unit": "mmol/g"}

    for T in T_list:
        b1_name = pname("b1", T)
        b2_name = pname("b2", T)
        if flagForQ1less:
            params.add(b1_name, value=init_params["b1"], min=0)
            params.add(b2_name, value=init_params["b2"], min=0)
        else:
            params.add(b1_name, value=init_params["b2"], min=0)
            params.add(b2_name, value=init_params["b1"], min=0)
        params[b1_name].user_data = {"unit": "1/Pa"}
        params[b2_name].user_data = {"unit": "1/Pa"}
    return params

def build_paramsByQmax(params, T_list, qmax, flagForQ1less=True):
    """
    构建 lmfit Parameters：
    - 全局：q1, q2
    - 每个温度：b1_T, b2_T
    并为每个参数添加单位属性
    """
    # params = Parameters()
    # 添加参数

    if flagForQ1less:
        params.add("q1", value=max(qmax*0.3, 0.1), min=0)
        params.add("q2", value=max(qmax*0.7, 0.1), min=0)
    else:
        params.add("q1", value=max(qmax*0.7, 0.1), min=0)
        params.add("q2", value=max(qmax*0.3, 0.1), min=0)
        
    # 保存单位
    params["q1"].user_data = {"unit": "mmol/g"}
    params["q2"].user_data = {"unit": "mmol/g"}

    for T in T_list:
        b1_name = pname("b1", T)
        b2_name = pname("b2", T)
        if flagForQ1less:
            params.add(b1_name, value=0.1, min=0)
            params.add(b2_name, value=0.05, min=0)
        else:
            params.add(b1_name, value=0.05, min=0)
            params.add(b2_name, value=0.1, min=0)
        params[b1_name].user_data = {"unit": "1/Pa"}
        params[b2_name].user_data = {"unit": "1/Pa"}

    return params

def build_params(T_list, qmax, flagForQ1less=True):
    params = Parameters()
    if flagForQ1less:
        params.add("q1", value=max(0.5, 0.1), min=0, max=qmax)
        params.add("q2", value=max(13.4, 0.1), min=0, max=qmax)
    else:
        params.add("q1", value=max(13.4, 0.1), min=0, max=qmax)
        params.add("q2", value=max(0.5, 0.1), min=0, max=qmax)

    for T in T_list:
        b1_name = pname("b1", T)
        b2_name = pname("b2", T)
        if flagForQ1less:
            params.add(b1_name, value=0.05,  min=1e-8, max=10)
            params.add(b2_name, value=0.1,  min=1e-8, max=10)
        else:
            params.add(b1_name, value=0.1,  min=1e-8, max=10)
            params.add(b2_name, value=0.05,  min=1e-8, max=10)
        params[b1_name].user_data = {"unit": "1/Pa"}
        params[b2_name].user_data = {"unit": "1/Pa"}
    return params

def build_initial_pool(T_list, qmax, flagForQ1less, sheet, file):
    pool = []

    # 随机 / 规则初值（qA > qB, qB > qA）
    pool.append(build_paramsByQmax(Parameters(), T_list, qmax, flagForQ1less))
    pool.append(build_paramsByQmax(Parameters(), T_list, qmax, not flagForQ1less))

    # Origin 初值（如果存在）
    if sheet:
        params0 = build_paramsByQmax(Parameters(), T_list, qmax, flagForQ1less)
        params_origin, valid = get_dslparams_from_origin(params0, T_list, sheet, file)
        if valid:
            pool.append(params_origin)

    return pool

def compute_sse(result, grouped_data):
    sse = 0.0
    for T, gdf in grouped_data.items():
        P = gdf["P_Pa"].values
        q_exp = gdf["q_mmol_g"].values

        p = getDSLpara(result, T)
        q_fit = dsl_model(P, p["q1"], p["b1"], p["q2"], p["b2"])

        sse += np.sum((q_exp - q_fit) ** 2)
    return sse

def dsl_diagnosis(result, T_list, site2_thresh=0.05):
    if not result.success:
        return False, "not_converged"

    if result.covar is None:
        return False, "singular_covar"

    for T in T_list:
        p = getDSLpara(result, T)

        q1, q2, b1, b2 = p["q1"], p["q2"], p["b1"], p["b2"]

        # 非物理
        if min(q1, q2, b1, b2) <= 0:
            return False, "nonphysical"

        # 第二 site 不可辨
        if q2 / (q1 + q2) < site2_thresh:
            return False, "site2_unidentifiable"

    return True, "ok"

def fit_with_constrained(initial_params, grouped_data):
    return minimize(
        residualForDSL,
        initial_params,
        args=(grouped_data,),
        method="trust-constr",
        max_nfev=20000
    )

def run_constrained_multistart(initial_pool, grouped_data, T_list):
    records = []

    for params in initial_pool:
        res = fit_with_constrained(params, grouped_data)
        ok, reason = dsl_diagnosis(res, T_list)

        records.append({
            "result": res,
            "ok": ok,
            "reason": reason,
            "sse": compute_sse(res, grouped_data) if ok else np.inf
        })

    valid = [r for r in records if r["ok"]]
    if not valid:
        return None, records

    best = min(valid, key=lambda x: x["sse"])
    return best["result"], records

def fit_with_lm(initial_params, grouped_data):
    return minimize(
        residualForDSL,
        initial_params,
        args=(grouped_data,),
        method="leastsq"   # lmfit 对应 LM
    )

def run_lm_multistart(initial_pool, grouped_data, T_list):
    records = []

    for params in initial_pool:
        res = fit_with_lm(params, grouped_data)
        ok, reason = dsl_diagnosis(res, T_list)

        records.append({
            "result": res,
            "ok": ok,
            "reason": reason,
            "sse": compute_sse(res, grouped_data) if ok else np.inf
        })

    valid = [r for r in records if r["ok"]]
    if not valid:
        return None, records

    best = min(valid, key=lambda x: x["sse"])
    return best["result"], records

def enforce_qA_less_qB(result, T_list):
    """
    规范化 DSL 的 site 标签，强制 q1 <= q2
    q1, q2 为全局参数
    b1, b2 为温度相关参数，需对所有 T 同步交换
    """

    q1 = result.params["q1"].value
    q2 = result.params["q2"].value

    # 若已经满足，不做任何事
    if q1 <= q2:
        return result

    # 1️⃣ 交换全局 q1 / q2
    result.params["q1"].value = q2
    result.params["q2"].value = q1

    # 2️⃣ 对所有 T 交换 b1(T) / b2(T)
    for T in T_list:
        name_b1 = pname("b1", T)
        name_b2 = pname("b2", T)

        b1 = result.params[name_b1].value
        b2 = result.params[name_b2].value

        result.params[name_b1].value = b2
        result.params[name_b2].value = b1

    return result

def compute_param_changes(dsl_orig, dsl_cropped, T_list):
    """
    计算 DSL 拟合参数在裁剪前后的变化百分比，直接使用 dsl 对象
    输入：
        dsl_orig: 原始拟合结果（MinimizerResult）
        dsl_cropped: 裁剪后拟合结果（MinimizerResult）
        T_list: 温度列表
    输出：
        param_changes: DataFrame，每行一个温度，每列参数变化百分比
    """
    
    records = []
    for T in T_list:
        orig = getDSLpara(dsl_orig, T)
        cropped = getDSLpara(dsl_cropped, T)
        record = {"T_K": T}
        for param in ["q1","q2","b1","b2"]:
            val_orig = orig[param]
            val_new = cropped[param]
            change = abs(val_new - val_orig)/val_orig if val_orig != 0 else float('nan')
            record[f"{param}_change"] = change
        records.append(record)
    
    param_changes = pd.DataFrame(records)
    return param_changes

def residualForDSL(params, grouped_data):
    q1 = params["q1"].value
    q2 = params["q2"].value
    res = []
    resid_all = []
    for T, gdf in grouped_data.items():
        b1 = params[pname("b1",T)].value
        b2 = params[pname("b2",T)].value
        P = gdf["P_Pa"].values
        q_obs = gdf["q_mmol_g"].values
        q_fit = dsl_model(P, q1, b1, q2, b2)
    #     res.append(q_obs - q_fit)
    # return np.concatenate(res)
        resid_all.append(q_obs - q_fit)
    return np.concatenate(resid_all)



def getDSLpara(dsl, T):
    """
    获取 DSL 拟合参数，返回字典（键名不带单位）
    单位信息保存在 dsl.params[param].unit 中
    """
    return {
        "q1": dsl.params["q1"].value,
        "q2": dsl.params["q2"].value,
        "b1": dsl.params[pname("b1", T)].value,
        "b2": dsl.params[pname("b2", T)].value
    }

def fit_dsl0(T_list, qmax, grouped_data, flagForQ1less=True, sheet = None, file = None):
    initial_pool = build_initial_pool(T_list, qmax, flagForQ1less, sheet, file)
    best_result = None

    best_result, lm_records = run_lm_multistart(initial_pool, grouped_data, T_list)

    # Step 3: fallback
    if best_result is None:
        best_result, tc_records = run_constrained_multistart(
            initial_pool, grouped_data, T_list
        )

    dsl = enforce_qA_less_qB(best_result, T_list)
    print(report_fit(dsl))


#     #############
#     # 构建初始拟合参数
#     params0 = Parameters()
#     validSheet = False
#     if sheet: ##origin fitting data
#         params1 = build_paramsByQmax(params0, T_list, qmax, flagForQ1less)
#         params, validSheet = get_dslparams_from_origin(params1, T_list, sheet, file)
#     else: 
#         params = build_paramsByQmax(params0, T_list, qmax, flagForQ1less)
#         # params = build_paramsBySL(T_list, grouped_data, flagForQ1less)
#         # params = build_params(T_list, qmax, flagForQ1less)          
# ######################
#     # 拟合 DSL  
#     dsl = minimize(residualForDSL, params, args=(grouped_data,), method="trust-constr", max_nfev=20000)
#     print(report_fit(dsl))

    fits_rows = []
    all_P = []
    all_q_exp = []
    all_q_fit = []
    for T in T_list:
        gdf = grouped_data[T]

        # 获取拟合参数（保持原名，单位保存在属性里）
        P = gdf["P_Pa"].values
        q_exp = gdf["q_mmol_g"].values
        params_dict = getDSLpara(dsl, T)
        q1 = params_dict["q1"]
        q2 = params_dict["q2"]
        b1 = params_dict["b1"]
        b2 = params_dict["b2"]

        # 拟合曲线
        q_fit = dsl_model(
            P,
            q1, b1,
            q2, b2
        )
        all_P.append(P)
        all_q_exp.append(q_exp)
        all_q_fit.append(q_fit)

         # 计算 R², RMSD
        residuals, relative_residuals, r2, rmsd = mf.compute_fit_metrics(
            fit_x=P, fit_y=q_fit, exp_x=P, exp_y=q_exp
            )

        # 构建 DataFrame，列名保持原始名称
        fits_rows.append(pd.DataFrame({
            "T_K": T,
            "P_Pa": P,
            "q_exp_mmol_g": q_exp,
            "q_fit_mmol_g": q_fit,
            "q1": q1,
            "q2": q2,
            "b1": b1,
            "b2": b2,
            "R2": r2,  # 可以重复写入每行，方便导出 Excel 查看,
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

    return dsl, fits_data
    # return dsl, fits_data, validSheet

def fit_dsl(rawData, out_path = None, sheet_name=None):
    """
    拟合 DSL 模型（双位点 Langmuir）到原始数据
    输入 rawData：已转换为基础单位的 DataFrame，包含列 [P_Pa, q_mmol_g, T_K]
    输出：
        - dsl 拟合结果（MinimizerResult）
        - fits_data DataFrame，包含列 ["T_K","P_Pa","q_obs","q_fit","q1","q2","b1","b2"],
        - q1","q2" unit: mmol_g
        - "b1","b2" unit: 1/Pa
          单位信息保存在 dsl.params[param].unit 中
    """

    # 按温度分组
    grouped_data = {T: gdf for T, gdf in rawData.groupby("T_K")}
    T_list = sorted(grouped_data.keys())

    # qmax 取整个数据的最大值
    qmax = rawData["q_mmol_g"].max()

    invalidOringinPar = True
    if sheet_name:##origin fitting para
        dsl, fits_data = fit_dsl0(T_list, qmax, grouped_data, sheet = sheet_name, file = out_path)
        # dsl, fits_data, validSheet = fit_dsl0(T_list, qmax, grouped_data, sheet = sheet_name, file = out_path)
    #     if validSheet == True:
    #         invalidOringinPar = False

    # if invalidOringinPar: 
    #     dsl1, fits_data1,_ = fit_dsl0(T_list, qmax, grouped_data, flagForQ1less=True)
    #     dsl2, fits_data2,_ = fit_dsl0(T_list, qmax, grouped_data, flagForQ1less=False)
        
    #     unique_T1 = fits_data1[["T_K", "R2"]].drop_duplicates().sort_values("T_K")
    #     r2_1 = unique_T1.iloc[0]["R2"]

    #     unique_T2 = fits_data2[["T_K", "R2"]].drop_duplicates().sort_values("T_K")
    #     r2_2 = unique_T2.iloc[0]["R2"]

    #     if r2_2 > r2_1:
    #         dsl = dsl2
    #         fits_data = fits_data2
    #     else:
    #         dsl = dsl1
    #         fits_data = fits_data1

    return dsl, fits_data



def find_P_for_q(q_target, T, params, P_min=1e-6, P_max=1e9):
    """数值解出 DSL 中给定 q 对应 P"""
    def f(P):
            b1_param = params.get(pname("b1", T))
            b2_param = params.get(pname("b2", T))
            if b1_param is None or b2_param is None:
                #print(f"Missing b1 or b2 at T={T}")
                return np.nan
            b1 = b1_param.value
            b2 = b2_param.value
            val = dsl_model(P, params["q1"].value, b1, params["q2"].value, b2) - q_target
            #print(f"P={P}, f(P)={val}")
            return val
    fmin, fmax = f(P_min), f(P_max)
    if fmin*fmax>0: raise ValueError("q_target outside achievable range")
    return brentq(f, P_min, P_max, maxiter=200)

def compute_HSG(params, T_list, R=8.314):
    """计算 DSL 每个位点的 H/S/G, 存在错误, 没计算斜率"""
    q1 = params["q1"].value
    q2 = params["q2"].value
    rows = []
    for T in T_list:
        b1 = params[pname("b1", T)].value
        b2 = params[pname("b2", T)].value
        # Van't Hoff: b = exp(-ΔH/RT + ΔS/R)
        # 解析 H/S/G
        # 这里采用 ln(b*R*T) = -H/(R*T) + S/R
        H1 = -R*T*np.log(b1*R*T)
        H2 = -R*T*np.log(b2*R*T)
        S1 = R*np.log(b1*R*T) + H1/T
        S2 = R*np.log(b2*R*T) + H2/T
        G1 = H1 - T*S1
        G2 = H2 - T*S2
        rows.append({"T":T,"q1":q1,"q2":q2,
                     "b1":b1,"b2":b2,
                     "H1_Jmol":H1,"H2_Jmol":H2,
                     "S1_JmolK":S1,"S2_JmolK":S2,
                     "G1_Jmol":G1,"G2_Jmol":G2})
    return pd.DataFrame(rows)

def compute_vant_hoff(params, T_list_K, P0=1e5, R=8.314462618):
    """
    根据拟合得到的 b_i(T)（params 中的 b1_*, b2_*），对每个位点做 Van't Hoff 回归。
    输入:
      - params: lmfit Parameters (含 b1_<T>, b2_<T>)
      - T_list_K: iterable of temperatures in K (exact same values used when creating params)
      - P0: 标准压强，Pa (默认 1e5 Pa = 1 bar)
      - R: 气体常数 J/mol/K
    输出:
      - thermo_df: DataFrame 每个位点一行，包含 dH (J/mol), dS (J/mol/K), R2, n_points, slope, intercept
      - dg_df: DataFrame 每行一个 (site, T, b, DeltaG_Jmol, DeltaG_kJmol)
    """
    sites = [1, 2]
    thermo_rows = []
        # 使用字典列表一次性存储列
    dg_dict = {
        "site": [],
        "T_K": [],
        "b": [],
        "DeltaG_Jmol": [],
        "DeltaG_kJmol": [],
        "DeltaG_Jmol_byHS": [],
        "DeltaG_kJmol_byHS": []
    }

    for i in sites:
        b_vals = []
        Ts_for_b = []
        # 提取每个温度的 b 值（跳过缺失或非正值）
        for T in T_list_K:
            key = pname(f"b{i}", T)
            p = params.get(key)
            if p is None:
                continue
            b = p.value
            if b is None or not np.isfinite(b) or b <= 0:
                continue
            b_vals.append(b)
            Ts_for_b.append(T)

        # 每个温度也计算 DeltaG（哪怕回归不可行），根据 Gibbs free energy–equilibrium constant relation: ΔG∘=−RTlnK
        for T in T_list_K:
            key = pname(f"b{i}", T)
            p = params.get(key)
            if p is None:
                continue
            b_here = p.value
            if (b_here is None) or (not np.isfinite(b_here)) or (b_here <= 0):
                dg = np.nan
            else:
                dg = -R * T * np.log(b_here * P0)   # J/mol
            # 先填入 dg_dict，DeltaG_byHS 先填 NaN，后续回归完成再更新
            dg_dict["site"].append(f"b{i}")
            dg_dict["T_K"].append(T)
            dg_dict["b"].append(b_here)
            dg_dict["DeltaG_Jmol"].append(dg)
            dg_dict["DeltaG_kJmol"].append(dg/1000.0 if np.isfinite(dg) else np.nan)
            dg_dict["DeltaG_Jmol_byHS"].append(np.nan)
            dg_dict["DeltaG_kJmol_byHS"].append(np.nan)

        if len(b_vals) < 2:
            # 不足以做 Van't Hoff 回归
            thermo_rows.append({
                "site": f"b{i}",
                "n": len(b_vals),
                "slope": np.nan,
                "intercept": np.nan,
                "R2": np.nan,
                "dH_Jmol": np.nan,
                "dH_kJmol": np.nan,
                "dS_Jmol_K": np.nan
            })
            continue

        x = 1.0 / np.array(Ts_for_b)    # 1/T (K^-1)
        y = np.log(np.array(b_vals) * P0)   # ln(b * P0), dimensionless

        lr = linregress(x, y)  # slope, intercept, r_value, p_value, stderr
        slope = lr.slope
        intercept = lr.intercept
        r2 = lr.rvalue**2

        dH = -slope * R            # J/mol
        dS = intercept * R         # J/mol*K

        # 更新 DeltaG_byHS 列
        for idx, T in enumerate(T_list_K):
            if T in Ts_for_b:
                dG_byHS = dH - T * dS  #delta_G  by: Gibbs–Helmholtz 方程的基本形式ΔG=ΔH−TΔS
                # key = pname(f"b{i}", T)
                # p = params.get(key)
                # b_here = p.value if (p is not None) else np.nan
                dg_dict["DeltaG_Jmol_byHS"][len(dg_dict["DeltaG_Jmol_byHS"]) - len(T_list_K) + idx] = dG_byHS
                dg_dict["DeltaG_kJmol_byHS"][len(dg_dict["DeltaG_kJmol_byHS"]) - len(T_list_K) + idx] = dG_byHS / 1000.0

        thermo_rows.append({
            "site": f"b{i}",
            "n": len(b_vals),
            "slope": slope,
            "intercept": intercept,
            "R2": r2,
            "dH_Jmol": dH,
            "dH_kJmol": dH / 1000.0,
            "dS_Jmol_K": dS,
            })

    thermo_df = pd.DataFrame(thermo_rows)
    dg_df = pd.DataFrame(dg_dict)

    return thermo_df, dg_df

def getQstByDslAndVanHoff(T_list_K, params, thermo_df):
    """DSL +Van't Hoff 计算Qst"""
    # Qst​(q,T) = (q1_sat​ * ΔH1​ * c1 * ​(1−c1​)+q2_sat​* ΔH2 * ​c2 *​(1−c2​)​) / (q1_sat​ * c1 * ​(1−c1​)+q2_sat​ * ​c2 *​(1−c2​)​)
    #ci​ = (bi​(T) * P​) / (1 + bi​(T) * P) 
    #总覆盖度 θtotal ​= (​q1​(T,P)+q2​(T,P)) / (q1_sat ​+ q2_sat) ​= (​q1_sat * ​θ1​+q2_sat * ​θ2​​) / (q1_sat ​+ q2_sat)
    #覆盖度分数 ​θi = (qi_sat * ​ci​​) / (q1_sat * ​c1​+q2_sat * ​c2​)
    
    if hasattr(params, "params"):  # 如果传入的是 MinimizerResult
        params = params.params
    q1_sat = params["q1"].value if hasattr(params["q1"], "value") else params["q1"]
    q2_sat = params["q2"].value if hasattr(params["q2"], "value") else params["q2"]
    dH1 = thermo_df.loc[thermo_df["site"]=="b1", "dH_Jmol"].values[0]
    dH2 = thermo_df.loc[thermo_df["site"]=="b2", "dH_Jmol"].values[0]

    P_max = 1e5 #unit: Pa
    P_grid=np.linspace(1e-6, P_max, 50) #使用点数量
    Qst_val = []

    for T in T_list_K:
        row = []
        
        b1 = params[pname("b1", T)].value
        b2 = params[pname("b2", T)].value
        # b1 = b0_1 * np.exp(DeltaS1/R) * np.exp(-DeltaH1/(R*T))
        # b2 = b0_2 * np.exp(DeltaS2/R) * np.exp(-DeltaH2/(R*T))

        for P in P_grid:
            c1 = b1*P/(1+b1*P)
            c2 = b2*P/(1+b2*P)

            Qst = -(q1_sat*dH1*c1*(1-c1) + q2_sat*dH2*c2*(1-c2)) / (q1_sat*c1*(1-c1) + q2_sat*c2*(1-c2))

            theta1 = q1_sat*c1/(q1_sat*c1 + q2_sat*c2)
            theta2 = q2_sat*c2/(q1_sat*c1 + q2_sat*c2)
            # theta = (q1_sat*theta1 + q2_sat*theta2) / (q1_sat + q2_sat)
            theta = (q1_sat * c1 + q2_sat * c2) / (q1_sat + q2_sat)
            q = q1_sat * c1 + q2_sat * c2
            row.append({
                "T_K": T,
                "P_Pa": P,
                "theta1": theta1,
                "theta2": theta2,
                "theta": theta,
                "q_mmol_g": q,
                "Qst_J_mol": Qst,
                "Qst_kJ_mol": Qst/1000
            })
        Qst_val.append(row)

    Qst_simul = [item for row in Qst_val for item in row]
    return pd.DataFrame(Qst_simul)

def compute_Qst(uptake_grid, T_list_K, params, R=8.314):
    """Clausius-Clapeyron Qst计算"""
    results = []
    T_list_K = list(T_list_K)  # 确保可索引
    for q_star in uptake_grid:
        P_list = []
        P_dict = {}   #用于导出excel
        skip = False
        for T in T_list_K:
            try:
                Pj = find_P_for_q(q_star, T, params)
                Pj = float(Pj)   # 强制转为标量 float
            except Exception as e:
                # 可以记录 e 或 q_star 以便调试
                skip = True
                break
            # 只接受有限且正的压力值
            if not np.isfinite(Pj) or Pj <= 0:
                skip = True
                break
            P_list.append(Pj)
            col_name = f"P_{int(T)}K" if float(T).is_integer() else f"P_{T}K"
            P_dict[col_name] = Pj

        if skip:
            continue

        # 确保 x, y 长度一致
        x = 1.0 / np.array(T_list_K)          # shape (n,)
        y = np.log(np.array(P_list))         # shape (n,)
        if x.shape[0] != y.shape[0]:
            continue
        if not np.all(np.isfinite(y)):
            continue

        slope, intercept, r_value, _, _ = linregress(x, y)
        Qst = -R * slope
        theta = q_star / (params["q1"].value + params["q2"].value)

        row = {
            "q_mmol_g": q_star,
            "theta": theta,
            "Qst_J_mol": Qst,
            "Qst_kJ_mol": Qst/1000.0,
            "R2": r_value**2,
            "slope": slope,
            "intercept": intercept
        }
        # 合并每个温度对应的压力列
        row.update(P_dict)
        results.append(row)
    return pd.DataFrame(results)


def getQst(dsl, thermo_df, rawData, uptake_num):  #目前只有DSL模型，rawdata拟合计算还没有开发

    grouped_data={T: gdf for T,gdf in rawData.groupby("T_K")}
    T_list_K=sorted(grouped_data.keys())

    uptake_min=max([gdf["q_mmol_g"].min() for gdf in grouped_data.values()])
    uptake_max=min([gdf["q_mmol_g"].max() for gdf in grouped_data.values()])
    uptake_grid=np.linspace(uptake_min,uptake_max,uptake_num) #使用点数量
    #uptake_grid = np.arange(uptake_min, uptake_min, uptake_step) #使用步长

    Qst_simuByCC=compute_Qst(uptake_grid, T_list_K, dsl.params)
    Qst_simul = getQstByDslAndVanHoff(T_list_K, dsl, thermo_df)

    return Qst_simuByCC, Qst_simul

def fitN2Isothem(T_max, listData, model="Henry"):
    T_C = T_max + const.K2C
    N2_PName = dop.getName(T_C, unitChange = 0, prefix="N2_", suffix = "Absolute Pressure (kPa)")
    N2_qName = dop.getName(T_C, unitChange = 0, prefix="N2_", suffix = "uptake (mmol/g)")
    list_N2  = dop.extract_PairData(N2_PName, N2_qName, listData, x_Scale = const.kPa2Pa)
    df_N2 = pd.DataFrame(list_N2, columns=["P_Pa", "q_mmol_g"])

    P = df_N2["P_Pa"].values
    q_exp = df_N2["q_mmol_g"].values

    params_dict = {}
    if model == "Langmuir": 
        n2_sl = sl.fit_single_langmuir(df_N2)

        # 拟合曲线
        q_sat=n2_sl.params["q_sat"].value
        b=n2_sl.params["b"].value
        q_fit = sl.langmuir_model(P, q_sat, b)
        params_dict = {"q_sat(mmol/g)": q_sat, "b(1/kPa)": b * const.kPa2Pa}

    elif model == "Henry":
        n2_sl = sl.fit_henry_model(df_N2)
        KH = n2_sl.params["KH"].value
        q_fit = sl.henry_model(P, KH)
        params_dict = {"KH(mmol/g*kPa)": KH * const.kPa2Pa}

    residuals, relative_residuals, r2, rmsd = mf.compute_fit_metrics(
            fit_x=P, fit_y=q_fit, exp_x=P, exp_y=q_exp
            )
            # 构建 DataFrame，列名保持原始名称

    result = {
    "T_C": T_C,
    "P_Pa": P,
    "q_exp_mmol_g": q_exp,
    "q_fit_mmol_g": q_fit,
    "model": model,
    "residuals": residuals,
    "relative_residuals": relative_residuals,
    "R2": r2,
    "RMSD": rmsd
    }
    result.update(params_dict)
    
    return pd.DataFrame(result)

def runDslSimu(listData, out_path, plotFlag=True, sheet_name=None):
    # 宽表转长表
        rawData0= dop.reshape_wide_to_long(listData)
        # 单位转换
        rawData = dop.convert_units(rawData0)
        #拟合DSL
        dsl,fits_data = fit_dsl(rawData, out_path, sheet_name)

        #单独拟合低压段数据
        # raw_low = dop.preprocess_low_pressure(rawData, P_max_Pa=5000)
        # dsl_low, fits_data_low = fit_dsl(raw_low, out_path, sheet_name)

        grouped_data={T: gdf for T,gdf in rawData.groupby("T_K")}
        T_list_K=sorted(grouped_data.keys())

    # 拟合N2吸附曲线
        T_max = max(grouped_data.keys())
        n2_fit = fitN2Isothem(T_max, listData, model="Henry")
    ##

        # #删除两个数据，用于比较b和q的变化是否敏感
        # cropped_df2 = dop.crop_low_points(rawData, n_remove=2)
        # dsl_cropped, _ = fit_dsl(cropped_df2, out_path, sheet_name)
        # param_changes = compute_param_changes(dsl, dsl_cropped, T_list_K)
        # print(param_changes)

        thermo_df, dg_df = compute_vant_hoff(dsl.params,T_list_K)
        
        # ther = compute_HSG(dsl.params,T_list_K) #该函数计算有问题，没计算斜率和截距，后续再修改 250930
        # Qst计算
        Qst_simuByCC, Qst_simul = getQst(dsl, thermo_df, rawData, 20)
        
        
        Qst_exp = qstByRawdata.runCalculateQstByRawData(listData)

        newPath = fl.check_or_create_folder(out_path, sub_folder_name=const.DSL_FOLDER)
        newFilePath = fl.get_expanded_name(newPath, sheet_name, expand=const.DSL_FILE)
        export_dsl_params_custom(dsl, file_name=newFilePath, sheet_name=const.DSL_PARA_SHEET)
        # export_dsl_params_custom(dsl_low, file_name=newFilePath, sheet_name=cn.DSL_LOW_PARA_SHEET)
        fl.export_to_excel_auto(fits_data, newFilePath, sheet_name=const.DSL_FITS_SHEET)
        # fl.export_to_excel_auto(fits_data_low, newFilePath, sheet_name=cn.DSL_LOW_FITS_SHEET)
        fl.export_to_excel_auto(thermo_df, newFilePath, sheet_name=const.THERMO_SHEET)
        fl.export_to_excel_auto(dg_df, newFilePath, sheet_name=const.GIBBS_SHEET)
        fl.export_to_excel_auto(Qst_simuByCC, newFilePath, sheet_name=const.QST_SIM_CC_SHEET)
        fl.export_to_excel_auto(Qst_simul, newFilePath, sheet_name=const.QST_SIM_SHEET)
        fl.export_to_excel_auto(Qst_exp, newFilePath, sheet_name=const.QST_EXP_SHEET)
        fl.export_to_excel_auto(n2_fit, newFilePath, sheet_name=const.N2_FITS_SHEET)

        if plotFlag:
            # 绘图excel
            plt.figure(figsize=(6,4))
            for T in T_list_K:
                gdf=grouped_data[T]
                plt.scatter(gdf["P_Pa"]/1000, gdf["q_mmol_g"], label=f"T={T-273.15:.1f}°C exp")
            for T in T_list_K:
                q_fit = fits_data[fits_data["T_K"]==T]["q_fit_mmol_g"].values
                P_fit = fits_data[fits_data["T_K"]==T]["P_Pa"].values
                plt.plot(P_fit/1000,q_fit,label=f"T={T-273.15}°C DSL")
            plt.xlabel("P (kPa)")
            plt.ylabel("q (mmol/g)")
            plt.title("DSL and experiment")
            plt.legend()
            plt.grid(True)
            plt.show(block=False)

            # Qst图
            plt.figure()
            plt.plot(Qst_simuByCC["q_mmol_g"],Qst_simuByCC["Qst_kJ_mol"],marker="*", label="Qst_simul")
            # 绘制实验点
            plt.scatter(Qst_exp["q_mmol_g"], Qst_exp["Qst_kJ_mol"], color="green", marker="o", s=50, label="Qst_exp")
            plt.xlabel("q (mmol/g)")
            plt.ylabel("Qst (kJ/mol)")
            plt.title("Qst (Clausius–Clapeyron)")
            plt.legend()
            plt.grid(True)
            plt.show(block=False)

            markers = ['o', 's', '^', 'v', 'D', '*']  # 不同标记
            colors = ['b', 'g', 'r', 'c', 'm', 'y'] # 不同颜色

            plt.figure(figsize=(6,4))

            for i, T in enumerate(T_list_K):
                df_T = Qst_simul[Qst_simul["T_K"] == T]
                plt.plot(df_T["theta"], df_T["Qst_kJ_mol"], marker=markers[i % len(markers)],
                    color=colors[i % len(colors)], label=f"T={T-273.15:.0f}°C")

            plt.xlabel("Coverage θ")
            plt.ylabel("Qst (kJ/mol)")
            plt.title("Coverage vs Qst at different temperatures")
            plt.legend()
            plt.grid(True)
            plt.show(block=True)

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
        runDslSimu(listData, out_path, plotFlag=True, sheet_name = sheet_name)

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

            runDslSimu(listData, out_path, plotFlag=False, sheet_name = sheet_name)

if __name__ == "__main__":
    f = sys.argv[1] if len(sys.argv) > 1 else None
    main(f)

