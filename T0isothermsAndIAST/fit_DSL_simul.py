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
import T1plotSource.plotOperation as myPlt
import T1dataProcessSource.residualBootstrap as redBoots
import T1dataProcessSource.modelAndFit as mf
import T0isothermsAndIAST.strongDSL as sl
import T0isothermsAndIAST.isostericHeatQstForRawData as qstByRawdata
# ==========================================
# 工具函数
# ==========================================
def plotDSLResults( grouped_data, T_list_K, fits_data, newFigPath=None):
    """Plot experimental isotherms and DSL fitting curves."""
    data, fits_by_model, _ = dop.prepareIsothermPlotData( grouped_data=grouped_data, T_list_K=T_list_K, 
                                                         fits_data={"DSL": fits_data}, )
    original_labels = {
        "0°C": "0°C Exp",
        "15°C": "15°C Exp",
        "25°C": "25°C Exp",
    }
    fit_labels = {
        "0°C": "0°C Fit",
        "15°C": "15°C Fit",
        "25°C": "25°C Fit",
    }
    ax = myPlt.plotCurve(
        data=data,
        fit=fits_by_model["DSL"],
        x="Pressure",
        y="Uptake",
        xlabel="Pressure (kPa)",
        ylabel=r"CO$_2$ uptake (mmol g$^{-1}$)",
        original_labels=original_labels,
        fit_label=fit_labels,
        marker=True,
        line=False,
        colors="PAPER1",
        markersize=45,
        fit_linewidth=2,
        fit_linestyle="-",
        figsize=(6, 4),
        textSize=12,
        legend=True,
        legendPosition="upper left",
        title=None,
        savepath=newFigPath,
    )
    plt.show(block=False)
    return ax
def plotQstComparison( Qst_simuByCC, Qst_exp, newFigPath=None):
    """Plot simulated and experimental Qst against uptake."""
    qst_data = {
        "Qst_simul": ( Qst_simuByCC["q_mmol_g"].to_numpy(),
            Qst_simuByCC["Qst_kJ_mol"].to_numpy(), ),
        "Qst_exp": ( Qst_exp["q_mmol_g"].to_numpy(),
            Qst_exp["Qst_kJ_mol"].to_numpy(), ), }
    ax = myPlt.plotCurve(
        data=qst_data,
        x="Uptake",
        y="Qst",
        xlabel=r"$q$ (mmol g$^{-1}$)",
        ylabel=r"$Q_{st}$ (kJ mol$^{-1}$)",
        original_labels={
            "Qst_simul": "Qst_simul",
            "Qst_exp": "Qst_exp", },
        marker=True,
        line=True,
        colors=["blue", "green"],
        markersize=50,
        legend=True,
        legendPosition="upper right",
        figsize=(6, 4),
        textSize=12,
        title="Qst (Clausius–Clapeyron)",
        savepath=newFigPath,
    )
    plt.show(block=False)
    return ax
def plotQstCoverage( Qst_simul, T_list_K, newFigPath=None ):
    """Plot Qst against coverage at different temperatures."""
    qst_coverage_data = {}
    temperature_labels = {}
    for T in T_list_K:
        temp_C = T - 273.15
        name = f"{temp_C:.0f}°C"
        df_T = Qst_simul[Qst_simul["T_K"] == T]
        qst_coverage_data[name] = (
            df_T["theta"].to_numpy(),
            df_T["Qst_kJ_mol"].to_numpy(), )
        temperature_labels[name] = f"T = {temp_C:.0f}°C"
    ax = myPlt.plotCurve(
        data=qst_coverage_data,
        x="Coverage",
        y="Qst",
        xlabel=r"Coverage $\theta$",
        ylabel=r"$Q_{st}$ (kJ mol$^{-1}$)",
        original_labels=temperature_labels,
        marker=True,
        line=True,
        colors="PAPER1",
        markersize=35,
        linewidth=2,
        legend=True,
        legendPosition="upper right",
        figsize=(6, 4),
        textSize=12,
        title="Coverage vs Qst at different temperatures",
        savepath=newFigPath,
    )
    plt.show(block=False)
    return ax
def makeDSLFitPredictFunc(
    T_list,
    *,
    qmax,
    fit_dsl0_func,
    getDSLpara_func,
    dsl_model_func,
    flagForQ1less=True,
    sheet=None,
    file=None,
):
    """
    Adapter for the existing multi-temperature DSL model.
    Shared parameters:
        q1, q2
    Temperature-dependent parameters:
        b1(T), b2(T)
    """
    temperatures = list(T_list)
    qmax = float(qmax)
    if not np.isfinite(qmax):
        raise ValueError("qmax must be finite.")
    def _fit_predict_func( grouped_data, *, initial_fit=None, verbose=False, ):
        # 1. Validate the temperature groups.
        missing = set(temperatures) - set(grouped_data)
        extra = set(grouped_data) - set(temperatures)
        if missing or extra:
            raise ValueError( f"Temperature-group mismatch: missing={missing}, extra={extra}." )
        # 2. Recover the original fit parameters for Bootstrap refitting.
        initial_params = None
        if initial_fit is not None:
            native_result = initial_fit.get("fit_result")
            if ( native_result is None or not hasattr(native_result, "params") ):
                raise ValueError( "The initial DSL fit does not expose optimizer parameters." )
            initial_params = native_result.params.copy()
        # Origin parameters are only relevant to the original fit.
        current_sheet = sheet if initial_fit is None else None
        current_file = file if initial_fit is None else None
        # 3. Fit all temperature groups jointly.
        dsl, fits_data = fit_dsl0_func( temperatures, qmax, grouped_data, flagForQ1less=flagForQ1less, 
                sheet=current_sheet, file=current_file, initial_params=initial_params, verbose=verbose, )
        # 4. Extract the shared parameters.
        q1 = float(dsl.params["q1"].value)
        q2 = float(dsl.params["q2"].value)
        parameters = { "q1_mmol_g": q1, 
                       "q2_mmol_g": q2, }
        fitted_values = {}
        # 5. Calculate predictions in the original row order.
        for T in temperatures:
            df = grouped_data[T]
            P = df["P_Pa"].to_numpy(dtype=float)
            para = getDSLpara_func(dsl, T)
            q1_T = float(para["q1"])
            q2_T = float(para["q2"])
            b1 = float(para["b1"])
            b2 = float(para["b2"])
            if not ( np.isclose(q1_T, q1) and np.isclose(q2_T, q2) ):
                raise ValueError( f"Shared DSL capacities are inconsistent at T={T!r}." )
            q_fit = np.asarray( dsl_model_func(P, q1, b1, q2, b2), dtype=float, ).reshape(-1)
            if len(q_fit) != len(df):
                raise ValueError( f"DSL prediction length mismatch at T={T!r}." )
            if not np.isfinite(q_fit).all():
                raise ValueError( f"DSL predictions contain non-finite values at T={T!r}." )
            fitted_values[T] = q_fit
            # Stable parameter names across all Bootstrap iterations.
            T_token = format(float(T), ".12g")
            parameters[f"b1_Pa_inv_T{T_token}K"] = b1
            parameters[f"b2_Pa_inv_T{T_token}K"] = b2
        # 6. Record boundary diagnostics.
        boundary_hits = []
        for name, param in dsl.params.items():
            value = float(param.value)
            tol = 1e-7 * max(1.0, abs(value))
            at_min = ( np.isfinite(param.min) and abs(value - float(param.min)) <= tol )
            at_max = ( np.isfinite(param.max) and abs(value - float(param.max)) <= tol )
            if at_min or at_max:
                boundary_hits.append(str(name))
        diagnostics = {
            "message": str(getattr(dsl, "message", "")),
            "nfev": getattr(dsl, "nfev", None),
            "chisqr": getattr(dsl, "chisqr", None),
            "redchi": getattr(dsl, "redchi", None),
            "boundary_hits": boundary_hits,
            "site_order_q1_le_q2": bool(q1 <= q2),
        }
        return {
            "fit_result": dsl,
            "fitted_values": fitted_values,
            "parameters": parameters,
            "success": bool(getattr(dsl, "success", True)),
            "diagnostics": diagnostics,
        }
    return _fit_predict_func
def buildDSLTemperatureContrastCITable(
    contrast,
    contrast_ci,
    ci_result,
    *,
    sample_name=None,
    parameter_key_a,
    parameter_key_b,
    metadata=None,
    unit_scale=1e3,
    output_unit="1/kPa",
    confidence_level=0.95,
    ci_method="percentile",
):
    """
    整理同一样品不同温度下 DSL 参数差值的 Bootstrap CI。
    Parameters
    ----------
    contrast : dict
        runBootstrapContrastAnalysis() 的返回结果。
    contrast_ci : dict
        差值的 CI 计算结果，至少包含 lower_bound 和 upper_bound。
    ci_result : dict
        runResidualBootstrapCIAnalysis() 的返回结果，
        用于提取两个温度下各参数的原始拟合估计值。
    sample_name : str or None
        样品名称。
    parameter_key_a, parameter_key_b : str
        两个需要比较的参数键，例如：
        "b1_Pa_inv_T273.15K"
        "b1_Pa_inv_T298.15K"
    metadata : dict or None
        可选的附加信息。
    unit_scale : float
        输出单位的换算系数。DSL 的 b 参数从 1/Pa 转为
        1/kPa 时使用 1000；若保持 1/Pa，则使用 1。
    output_unit : str
        输出参数和差值的单位标签。
    confidence_level : float
        备用置信水平。当 contrast_ci 和 ci_result 都没有
        confidence_level 时使用。
    ci_method : str
        备用 CI 方法。当 contrast_ci 中没有 method 时使用。
    Returns
    -------
    pd.DataFrame
        单行结果表，包含两个温度的参数估计值、参数差值、
        Bootstrap SE、95% CI 以及诊断信息。
    """
    # ==================================================
    # 1. 检查输入
    # ==================================================
    if not isinstance(contrast, dict):
        raise TypeError("contrast must be a dictionary.")
    if not isinstance(contrast_ci, dict):
        raise TypeError("contrast_ci must be a dictionary.")
    if not isinstance(ci_result, dict):
        raise TypeError("ci_result must be a dictionary.")
    statistics = ci_result.get("statistics")
    if not isinstance(statistics, dict):
        raise ValueError(
            "ci_result must contain a 'statistics' dictionary."
        )
    if (
        parameter_key_a not in statistics
        or parameter_key_b not in statistics
    ):
        raise KeyError(
            "Parameter keys not found in ci_result['statistics']: "
            f"{parameter_key_a}, {parameter_key_b}"
        )
    # ==================================================
    # 2. 解析参数名称和温度
    # ==================================================
    pattern = re.compile(
        r"^(b\d+)_.*_T([-+]?\d+(?:\.\d+)?)K$"
    )
    match_a = pattern.fullmatch(parameter_key_a)
    match_b = pattern.fullmatch(parameter_key_b)
    if match_a is None or match_b is None:
        raise ValueError(
            "Unable to parse DSL parameter keys. Expected names "
            "like 'b1_Pa_inv_T273.15K'."
        )
    parameter_a = match_a.group(1)
    parameter_b = match_b.group(1)
    # 当前函数处理同一位点在不同温度下的比较
    if parameter_a != parameter_b:
        raise ValueError(
            "Temperature contrast must compare the same DSL "
            f"parameter. Got {parameter_a} and {parameter_b}."
        )
    T_a = float(match_a.group(2))
    T_b = float(match_b.group(2))
    # ==================================================
    # 3. 提取原始参数估计值
    # ==================================================
    estimate_a = float(
        statistics[parameter_key_a]["estimate"]
    )
    estimate_b = float(
        statistics[parameter_key_b]["estimate"]
    )
    # 原始数据拟合得到的参数差值
    difference = float(contrast["estimate"])
    # 检查差值与两个参数估计值是否一致
    expected_difference = estimate_a - estimate_b
    if not np.isclose(
        difference,
        expected_difference,
        rtol=1e-6,
        atol=1e-12,
    ):
        raise ValueError(
            "The contrast estimate is inconsistent with the "
            "two original parameter estimates. Check the "
            "contrast calculation and parameter ordering."
        )
    # ==================================================
    # 4. 提取差值 CI 和 Bootstrap SE
    # ==================================================
    lower = float(contrast_ci["lower_bound"])
    upper = float(contrast_ci["upper_bound"])
    bootstrap_se = float(
        contrast["bootstrap_standard_error"]
    )
    if (
        not np.isfinite(lower)
        or not np.isfinite(upper)
        or lower > upper
    ):
        raise ValueError(
            "Invalid confidence interval bounds."
        )
    # ==================================================
    # 5. 兼容 CI 元数据：优先读取已有结果
    # ==================================================
    resolved_confidence_level = contrast_ci.get(
        "confidence_level",
        ci_result.get(
            "confidence_level",
            confidence_level,
        ),
    )
    resolved_ci_method = contrast_ci.get(
        "method",
        ci_method,
    )
    # ==================================================
    # 6. 单位转换
    # ==================================================
    unit_scale = float(unit_scale)
    if not np.isfinite(unit_scale) or unit_scale <= 0:
        raise ValueError(
            "unit_scale must be a positive finite number."
        )
    estimate_a_out = estimate_a * unit_scale
    estimate_b_out = estimate_b * unit_scale
    difference_out = difference * unit_scale
    bootstrap_se_out = bootstrap_se * unit_scale
    lower_out = lower * unit_scale
    upper_out = upper * unit_scale
    # ==================================================
    # 7. 整理单行输出
    # ==================================================
    row = dict(metadata or {})
    if sample_name is not None:
        row["Sample"] = sample_name
    row.setdefault("Model", "DSL")
    row.setdefault(
        "Analysis",
        "Within-sample temperature contrast",
    )
    row.update({
        "Parameter": parameter_a,
        "Parameter_Key_A": parameter_key_a,
        "Parameter_Key_B": parameter_key_b,
        "Temperature_A_K": T_a,
        "Temperature_A_C": T_a - 273.15,
        "Temperature_B_K": T_b,
        "Temperature_B_C": T_b - 273.15,
        f"Estimate_A_{output_unit}": estimate_a_out,
        f"Estimate_B_{output_unit}": estimate_b_out,
        f"Difference_A_minus_B_{output_unit}": difference_out,
        f"Bootstrap_SE_{output_unit}": bootstrap_se_out,
        f"CI_Lower_{output_unit}": lower_out,
        f"CI_Upper_{output_unit}": upper_out,
        "Confidence_Level": resolved_confidence_level,
        "CI_Method": resolved_ci_method,
        "Bootstrap_Method": ci_result.get("method"),
        "Contrast_Mode": contrast.get("mode"),
        "N_Contrast_Resamples": contrast.get("n_valid"),
        "N_Contrast_Input": contrast.get("n_input"),
        "Contrast_Random_State": contrast.get("random_state"),
        "N_Bootstrap_A": statistics[
            parameter_key_a
        ].get("n_valid"),
        "N_Bootstrap_B": statistics[
            parameter_key_b
        ].get("n_valid"),
        "CI_Excludes_Zero": (
            lower > 0 or upper < 0
        ),
    })
    return pd.DataFrame([row])
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
    # -----------------------------
    # 检查文件是否存在
    # -----------------------------
    if file is None:
        print(f"[Warning] File not found: {file}")
        return params, False
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
    return minimize( residualForDSL, initial_params, args=(grouped_data,), method="trust-constr", max_nfev=20000 )
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
    return minimize( residualForDSL, initial_params, args=(grouped_data,), method="leastsq" )  # lmfit 对应 LM )
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
def fit_dsl0( T_list, qmax, grouped_data, flagForQ1less=True, sheet=None, file=None, *, initial_params=None, verbose=True, ):
    if initial_params is None:
        # Original fit: retain the existing multistart strategy.
        initial_pool = build_initial_pool( T_list, qmax, flagForQ1less, sheet, file, )
        best_result, lm_records = run_lm_multistart( initial_pool, grouped_data, T_list, )
        if best_result is None:
            best_result, tc_records = run_constrained_multistart( initial_pool, grouped_data, T_list, ) 
    else:
        # Bootstrap refit: first try the original fit's optimum
        # as the starting point.
        best_result = None
        for fit_func in (fit_with_lm, fit_with_constrained):
            try:
                candidate = fit_func( initial_params.copy(), grouped_data, )
                ok, reason = dsl_diagnosis(candidate, T_list)
                if ok:
                    best_result = candidate
                    break
            except Exception:
                # Try the next fitting strategy.
                continue
        # Fallback: retain the existing multistart strategy.
        if best_result is None:
            fallback_pool = build_initial_pool( T_list, qmax, flagForQ1less, None, None, )
            best_result, lm_records = run_lm_multistart( fallback_pool, grouped_data, T_list, )
            if best_result is None:
                best_result, tc_records = run_constrained_multistart( fallback_pool, grouped_data, T_list, )
    if best_result is None:
        raise RuntimeError( "DSL fitting failed: no acceptable solution was found." )
    # Canonicalize site labels: q1 <= q2.
    dsl = enforce_qA_less_qB(best_result, T_list)
    if verbose:
        print(report_fit(dsl))
    # Keep your existing fits_rows / metrics / fits_data code below.
    # initial_pool = build_initial_pool(T_list, qmax, flagForQ1less, sheet, file)
    # best_result = None
    # best_result, lm_records = run_lm_multistart(initial_pool, grouped_data, T_list)
    # # Step 3: fallback
    # if best_result is None:
    #     best_result, tc_records = run_constrained_multistart( initial_pool, grouped_data, T_list )
    # dsl = enforce_qA_less_qB(best_result, T_list)
    # print(report_fit(dsl))
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
        q_fit = dsl_model( P, q1, b1, q2, b2 )
        all_P.append(P)
        all_q_exp.append(q_exp)
        all_q_fit.append(q_fit)
         # 计算 R², RMSD
        residuals, relative_residuals, r2, rmsd = mf.compute_fit_metrics( fit_x=P, fit_y=q_fit, exp_x=P, exp_y=q_exp )
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
    global_residuals, global_relative_residuals, global_r2, global_rmsd = mf.compute_fit_metrics( fit_x=P_all, 
            fit_y=q_fit_all, exp_x=P_all, exp_y=q_exp_all,
            use_interp=False   # 全局
        )
    fits_data = pd.concat(fits_rows, ignore_index=True)
    fits_data["global_R2"] = global_r2
    fits_data["global_RMSD"] = global_rmsd
    fits_data["global_residuals"] = global_residuals
    fits_data["global_relative_residuals"] = global_relative_residuals
    return dsl, fits_data
    # return dsl, fits_data, validSheet
def fit_dsl(rawData, out_path = None, sheet_name=None, run_bootstrap=False,n_bootstrap=1000,):
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
        # Fit DSL regardless of whether Origin data are available.
    # Origin parameters are optional initial guesses.
    dsl, fits_data = fit_dsl0( T_list, qmax, grouped_data, sheet=sheet_name, file=out_path, )
    #对全局拟合的参数进行residual bootstrap 操作计算置信区间confidence intervals
    dsl_bootstrap = None
    dsl_ci = None
    # modelCIcheck = True
    if run_bootstrap:
        # Create the DSL adapter.
        dsl_fit_predict = makeDSLFitPredictFunc( T_list, qmax=qmax, fit_dsl0_func=fit_dsl0, 
                getDSLpara_func=getDSLpara, dsl_model_func=dsl_model, flagForQ1less=True, 
                sheet=sheet_name, file=out_path, )
        # First, test the workflow with 100 iterations.
        dsl_bootstrap = redBoots.runResidualBootstrapAnalysis( grouped_data, dsl_fit_predict, 
                response_col="q_mmol_g", n_bootstrap=n_bootstrap, residual_centering="within_group", random_state=2026, )
        dsl_ci = redBoots.runResidualBootstrapCIAnalysis( dsl_bootstrap, 
                confidence_interval_func=dop.runConfidenceIntervalAnalysis, confidence_level=0.95, )
        # 4. Select the two parameters to compare.
        name_a = "b1_Pa_inv_T273.15K"
        name_b = "b1_Pa_inv_T298.15K"
        # Retrieve Bootstrap distributions.
        # Support either of the two result structures.
        dist_a = dsl_bootstrap["parameter_distributions"][name_a]
        dist_b = dsl_bootstrap["parameter_distributions"][name_b]
        # 5. Construct the Bootstrap distribution of the difference.
        # Paired mode preserves iteration-wise correspondence because
        # both parameters came from the same joint DSL Bootstrap run.
        b1_contrast = redBoots.runBootstrapContrastAnalysis( dist_a, dist_b,
            estimate_a=dsl_ci["statistics"][name_a]["estimate"],
            estimate_b=dsl_ci["statistics"][name_b]["estimate"],
            mode="paired", 
        )
        # 6. Calculate the 95% CI of the difference.
        b1_contrast_ci = dop.runConfidenceIntervalAnalysis( method="percentile", confidence_level=0.95,
                        estimate=b1_contrast["estimate"], bootstrap_distribution=b1_contrast["bootstrap_distribution"], )
        lower = b1_contrast_ci["lower_bound"]
        upper = b1_contrast_ci["upper_bound"]
        print(f"Contrast: {name_a} - {name_b}")
        print(f"Difference: {b1_contrast['estimate']:.6g}")
        print(f"95% CI: ({lower:.6g}, {upper:.6g})")
        print( f"Bootstrap SE: {b1_contrast['bootstrap_standard_error']:.6g}" )
        print( "Difference supported:", lower > 0 or upper < 0, )
        for name, result in dsl_ci["statistics"].items():
            lower, upper = result["confidence_interval"]
            print(
                f"{name}: "
                f"estimate={result['estimate']:.6g}, "
                f"95% CI=({lower:.6g}, {upper:.6g}), "
                f"SE={result['bootstrap_standard_error']:.6g}, "
                f"valid={result['n_valid']}, "
                f"invalid={result['n_invalid']}"
            )
        print("Successful fits:", dsl_bootstrap["n_successful_fits"])
        print("Failed fits:", dsl_bootstrap["n_failed_fits"])
        print("Failure reasons:", dsl_bootstrap["failure_reasons"])
    if run_bootstrap:
        return dsl, fits_data, dsl_bootstrap, dsl_ci, b1_contrast, b1_contrast_ci
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
def runDslSimu(listData, out_path, plotFlag=True, sheet_name=None, run_bootstrap=False,n_bootstrap=1000,):
    # 宽表转长表
    rawData0= dop.reshape_wide_to_long(listData)
    # 单位转换
    rawData = dop.convert_units(rawData0)
    #拟合DSL
    # 拟合 DSL
    if run_bootstrap:
        dsl, fits_data, dsl_bootstrap, dsl_ci, b1_contrast, b1_contrast_ci = fit_dsl( rawData, out_path, sheet_name=sheet_name, run_bootstrap=True, n_bootstrap=n_bootstrap,)
    else:
        dsl, fits_data = fit_dsl( rawData, out_path, sheet_name=sheet_name, run_bootstrap=False, )
        dsl_bootstrap = None
        dsl_ci = None
        b1_contrast = None
        b1_contrast_ci = None
    #单独拟合低压段数据
    # raw_low = dop.preprocess_low_pressure(rawData, P_max_Pa=5000)
    # dsl_low, fits_data_low = fit_dsl(raw_low, out_path, sheet_name)
    grouped_data={T: gdf for T,gdf in rawData.groupby("T_K")}
    T_list_K=sorted(grouped_data.keys())
# 拟合N2吸附曲线
    T_max = max(grouped_data.keys())
    n2_fit = fitN2Isothem(T_max, listData, model="Langmuir")
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
    if run_bootstrap:
        temperature_contrast_df = (
            buildDSLTemperatureContrastCITable( contrast=b1_contrast, contrast_ci=b1_contrast_ci, ci_result=dsl_ci, sample_name=sheet_name,
                parameter_key_a="b1_Pa_inv_T273.15K",
                parameter_key_b="b1_Pa_inv_T298.15K",
                metadata={ "Model": "DSL", "Analysis": "Within-sample temperature contrast", }, ) )
        extra_tables = { "DSL_TemperatureContrast_CI": temperature_contrast_df, }
        redBoots.exportBootstrapAnalysisResults( bootstrap_result=dsl_bootstrap, ci_result=dsl_ci, 
                file_name=newFilePath, metadata={ "Sample": sheet_name, "Model": "DSL", "Analysis": "Joint fit", }, 
                extra_tables=extra_tables, )
    if plotFlag:
        newFigPath = fl.get_expanded_name(newPath, sheet_name, expand=const.DSL_FILE, type="png",)
    # 1. DSL fitting
        plotDSLResults( grouped_data=grouped_data, T_list_K=T_list_K, fits_data=fits_data, newFigPath=newFigPath)
        # 2. Simulated vs experimental Qst
        plotQstComparison( Qst_simuByCC=Qst_simuByCC, Qst_exp=Qst_exp, )
        # 3. Qst vs coverage
        plotQstCoverage( Qst_simul=Qst_simul, T_list_K=T_list_K, )
    return {
    "sample_name": sheet_name,
    "dsl": dsl,
    "fits_data": fits_data,
    "dsl_bootstrap": dsl_bootstrap,
    "dsl_ci": dsl_ci,
}
# ==========================================
# 主函数
# ==========================================
def main(file_path=None):
    # =========================
    # 运行配置
    # =========================
    singleFile = True
    run_bootstrap = False
    n_bootstrap = 100
    file = ( file_path if file_path is not None else fl.getFile() )
    all_sample_results = {}
    # 用于保存样品间比较结果的输出目录
    comparison_out_path = None
    # =========================
    # 模式一：单样品拟合
    # =========================
    if singleFile:
        sheet_name = "CC-Hy-900-1-1"
        listData, out_path, _ = fl.readFileBySheet( file, sheet_name, )
        result = runDslSimu( listData, out_path, plotFlag=True, sheet_name=sheet_name, run_bootstrap=run_bootstrap, n_bootstrap=n_bootstrap, )
        # 单样品模式不执行样品间比较
        return result
    # =========================
    # 模式二：批量处理全部样品
    # =========================
    xls = pd.ExcelFile(file)
    for sheet_name in xls.sheet_names:
        if fl.should_skip(sheet_name):
            continue
        listData, out_path, _ = fl.readFileBySheet( file, sheet_name, )
        # -------------------------
        # 检查数据
        # -------------------------
        if isinstance(listData, dict):
            if not listData:
                print(f"Skip empty data: {sheet_name}")
                continue
            listData0 = pd.DataFrame( next(iter(listData.values())) )
        elif isinstance(listData, pd.DataFrame):
            listData0 = listData
        else:
            print(f"Skip unsupported data: {sheet_name}")
            continue
        if listData0.empty:
            print(f"DataFrame is empty: {sheet_name}")
            continue
        # 记录比较结果的保存目录
        if comparison_out_path is None:
            comparison_out_path = out_path
        # -------------------------
        # 拟合当前样品
        # -------------------------
        print(f"\nProcessing sample: {sheet_name}")
        result = runDslSimu( listData, out_path, plotFlag=False, sheet_name=sheet_name, run_bootstrap=run_bootstrap, n_bootstrap=n_bootstrap, )
        all_sample_results[sheet_name] = result
    # =========================
    # 模式三：不同样品之间的参数 CI
    # =========================
    contrast_df = pd.DataFrame()
    comparison_file = None
    if run_bootstrap:
        # ----------------------------------
        # 方式 A：仅比较指定样品对
        # 使用时保留这部分，注释掉方式 B
        # ----------------------------------
        sample_pairs = [
            ("CC-Hy-600-0.5-1", "CC-Hy-600-1-1"),
            ("CC-Hy-600-0.5-1", "CC-Hy-700-1-1"),
            ("CC-Hy-600-2-1", "CC-Hy-900-1-1"),
            # 按需添加其他样品对
        ]
        # ----------------------------------
        # 方式 B：全部有效样品两两比较
        # 使用时注释掉方式 A，取消下面这部分注释
        # ----------------------------------
        # sample_names = [
        #     name
        #     for name, result in all_sample_results.items()
        #     if (
        #         result.get("dsl_bootstrap") is not None
        #         and result.get("dsl_ci") is not None
        #     )
        # ]
        #
        # sample_pairs = list(
        #     combinations(sample_names, 2)
        # )
        # ----------------------------------
        # 执行 DSL b1 样品间 CI 比较
        # ----------------------------------
        contrast_df = redBoots.compareSampleBootstrapParameters( all_sample_results=all_sample_results, 
                parameter_name="b1_Pa_inv_T298.15K", sample_pairs=sample_pairs, bootstrap_key="dsl_bootstrap", 
                ci_key="dsl_ci", confidence_level=0.95, )
        if not contrast_df.empty:
            # 增加比较任务的元数据
            contrast_df = contrast_df.copy()
            contrast_df.insert( 0, "Model", "DSL", )
            contrast_df.insert( 1, "Analysis", "Between-sample pairwise contrast", )
            print("\nDSL b1 pairwise contrast results:")
            print( contrast_df.to_string( index=False, float_format=lambda x: f"{x:.10g}", ) )
            # ----------------------------------
            # 保存到独立的 Excel 文件
            # ----------------------------------
            if comparison_out_path is None:
                raise RuntimeError( "Cannot determine output directory for sample-pair comparison results." )
            comparison_dir = fl.check_or_create_folder( comparison_out_path, sub_folder_name=const.DSL_FOLDER, )
            comparison_file = os.path.join( comparison_dir, "DSL_SamplePairwise_CI.xlsx", )
            # 统一保存入口：仅传入需要保存的结果表
            redBoots.exportBootstrapAnalysisResults( file_name=comparison_file, 
                                                    extra_tables={ "DSL_b1_Pairwise_CI": contrast_df, }, )
            print( "\nSample-pair CI results saved to:\n" f"{comparison_file}" )
        else:
            print(
                "No valid sample-pair comparison results "
                "were generated. Nothing was saved."
            )
    return all_sample_results, contrast_df
if __name__ == "__main__":
    f = sys.argv[1] if len(sys.argv) > 1 else None
    main(f)
