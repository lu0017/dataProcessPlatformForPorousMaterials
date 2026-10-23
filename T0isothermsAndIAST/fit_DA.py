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

import T1dataProcessSource.modelAndFit as mf
import T1dataProcessSource.residualBootstrap as redBoots

# ==========================================
# 工具函数
# ==========================================
def printResidualBootstrapSummary(
    bootstrap_result: Mapping[str, Any],
    ci_result: Mapping[str, Any] | None = None,
    *,
    group_labels: Mapping[Hashable, str] | None = None,
    response_unit: str = "mmol/g",
    digits: int = 4,
) -> None:
    """Print a human-readable summary of a residual-bootstrap analysis.

    Parameters
    ----------
    bootstrap_result
        Result returned by ``runResidualBootstrapAnalysis``.
    ci_result
        Optional result returned by ``runResidualBootstrapCIAnalysis``.
    group_labels
        Optional display labels for group keys, e.g.
        ``{273.15: "0 °C (273.15 K)"}``.
    response_unit
        Unit displayed in the residual-diagnostic table.
    digits
        Number of decimal places for residual diagnostics.
    """
    def _fmt(value: Any, ndigits: int = digits) -> str:
        try:
            value = float(value)
        except (TypeError, ValueError):
            return str(value)
        if not np.isfinite(value):
            return "NA"
        return f"{value:.{ndigits}f}"

    def _group_label(key: Any) -> str:
        labels = group_labels or {}
        if key in labels:
            return str(labels[key])
        # Some diagnostic counters store group keys as strings.
        if str(key) in labels:
            return str(labels[str(key)])
        return str(key)

    n_bootstrap = int(bootstrap_result.get("n_bootstrap", 0))
    n_successful = int(bootstrap_result.get("n_successful_fits", 0))
    n_failed = int(bootstrap_result.get("n_failed_fits", max(0, n_bootstrap - n_successful)))
    success_rate = 100.0 * n_successful / n_bootstrap if n_bootstrap else float("nan")

    print("=" * 88)
    print("RESIDUAL BOOTSTRAP SUMMARY")
    print("=" * 88)
    print(f"Method                 : {bootstrap_result.get('method', 'NA')}")
    print(f"Bootstrap iterations   : {n_bootstrap:,}")
    print(f"Successful refits      : {n_successful:,} ({_fmt(success_rate, 2)}%)")
    print(f"Failed refits          : {n_failed:,}")
    print(f"Residual centering     : {bootstrap_result.get('residual_centering', 'NA')}")
    print(f"Random seed            : {bootstrap_result.get('random_state', 'NA')}")

    failure_reasons = bootstrap_result.get("failure_reasons", {}) or {}
    if failure_reasons:
        print("\nFailure reasons:")
        for reason, count in sorted(failure_reasons.items()):
            print(f"  - {reason}: {count}")
    elif n_failed == 0:
        print("Failure reasons        : none")

    # Residual diagnostics: display group-specific offsets before centering.
    residual_diagnostics = bootstrap_result.get("residual_diagnostics", {}) or {}
    if residual_diagnostics:
        print("\nResidual diagnostics by group")
        print(f"Residual unit: {response_unit}; SD uses ddof=1")
        headers = [
            "Group", "N", "Mean before", "Center applied", "Mean after",
            "SD before", "SD after", "Min centered", "Max centered",
        ]
        rows = []
        for key, d in residual_diagnostics.items():
            label = _group_label(key)
            rows.append([
                label,
                str(d.get("n_points", "NA")),
                _fmt(d.get("mean_residual_before_centering")),
                _fmt(d.get("centering_value")),
                _fmt(d.get("mean_residual_after_centering")),
                _fmt(d.get("residual_sd_before_centering")),
                _fmt(d.get("residual_sd_after_centering")),
                _fmt(d.get("residual_min_after_centering")),
                _fmt(d.get("residual_max_after_centering")),
            ])
        widths = [max(len(headers[i]), max(len(row[i]) for row in rows)) for i in range(len(headers))]
        print("  ".join(headers[i].ljust(widths[i]) for i in range(len(headers))))
        print("  ".join("-" * widths[i] for i in range(len(headers))))
        for row in rows:
            print("  ".join(row[i].ljust(widths[i]) for i in range(len(headers))))
        print("\nNote: 'Mean before' shows the original group-wise residual bias. 'Mean after' is near zero by construction when within-group centering is used.")

    # Negative simulated responses, if tracked by the engine.
    negative_by_group = bootstrap_result.get("n_negative_simulated_responses_by_group", {}) or {}
    if negative_by_group:
        total_negative = int(sum(negative_by_group.values()))
        total_points = int(bootstrap_result.get("n_simulated_points", 0))
        percent_negative = 100.0 * total_negative / total_points if total_points else 0.0
        print("\nNegative simulated responses")
        print(f"Total negative responses : {total_negative:,} / {total_points:,} ({_fmt(percent_negative, 3)}%)")
        for key, count in negative_by_group.items():
            label = _group_label(key)
            print(f"  - {label}: {int(count):,}")

    boundary_hits = bootstrap_result.get("boundary_hit_counts", {}) or {}
    if boundary_hits:
        print("\nParameter-boundary hits across successful refits")
        for name, count in sorted(boundary_hits.items()):
            print(f"  - {name}: {count:,}")
    else:
        print("\nParameter-boundary hits : none recorded")

    if ci_result is not None:
        statistics = ci_result.get("statistics", {}) or {}
        print("\nParameter uncertainty")
        print(f"CI method: {ci_result.get('method', 'NA')}; confidence level: {_fmt(100 * ci_result.get('confidence_level', float('nan')), 1)}%")
        headers = ["Parameter", "Estimate", "Bootstrap SE", "CI lower", "CI upper", "Valid", "Invalid"]
        rows = []
        for name, d in statistics.items():
            ci = d.get("confidence_interval", (float("nan"), float("nan")))
            rows.append([
                str(name),
                _fmt(d.get("estimate"), max(digits, 5)),
                _fmt(d.get("bootstrap_standard_error"), max(digits, 5)),
                _fmt(ci[0], max(digits, 5)),
                _fmt(ci[1], max(digits, 5)),
                str(d.get("n_valid", "NA")),
                str(d.get("n_invalid", "NA")),
            ])
        if rows:
            widths = [max(len(headers[i]), max(len(row[i]) for row in rows)) for i in range(len(headers))]
            print("  ".join(headers[i].ljust(widths[i]) for i in range(len(headers))))
            print("  ".join("-" * widths[i] for i in range(len(headers))))
            for row in rows:
                print("  ".join(row[i].ljust(widths[i]) for i in range(len(headers))))

    print("=" * 88)
FitPredictFunc = Callable[..., Mapping[str, Any]]
def makeDAGlobalFitPredictFunc(
    T_list,
    *,
    fit_DA0_func: Callable[..., Any],
    getParaForDA_func: Callable[..., Any],
) -> FitPredictFunc:
    """Build an adapter around the project's global ``fit_DA0`` function.

    ``fit_DA0_func`` must support the keyword arguments ``verbose`` and
    ``initial_params``. Its existing three-item return convention is retained:
    ``(da_result, fits_data, valid_sheet)``.

    The fit data must contain ``q_fit_mmol_g`` and ``P_Pa`` columns, and its
    rows must be concatenated in ``T_list`` order, matching the current
    implementation of ``fit_DA0``.
    """
    temperatures = list(T_list)

    def _fit_predict_func(
        grouped_data,
        *,
        initial_fit=None,
        verbose=False,
    ):
        initial_params = None
        if initial_fit is not None:
            native_result = initial_fit.get("fit_result")
            if native_result is None or not hasattr(native_result, "params"):
                raise ValueError(
                    "The initial D-A fit output does not expose optimizer parameters."
                )
            initial_params = native_result.params

        da_result, fits_data, valid_sheet = fit_DA0_func(
            temperatures,
            grouped_data,
            verbose=verbose,
            initial_params=initial_params,
        )

        fitted_values: dict[Hashable, np.ndarray] = {}
        cursor = 0
        for T in temperatures:
            if T not in grouped_data:
                raise KeyError(f"Temperature {T!r} is missing from grouped_data.")

            source_df = grouped_data[T]
            n_rows = len(source_df)
            fit_group = fits_data.iloc[cursor:cursor + n_rows]
            cursor += n_rows

            if len(fit_group) != n_rows:
                raise ValueError(
                    f"D-A fitting output for temperature {T!r} has an unexpected row count."
                )

            input_p = source_df["P_Pa"].to_numpy(dtype=float)
            fit_p = fit_group["P_Pa"].to_numpy(dtype=float)
            if not np.allclose(input_p, fit_p, rtol=1e-10, atol=0.0):
                raise ValueError(
                    f"D-A fitted pressures do not align with input rows for temperature {T!r}."
                )

            fitted_values[T] = fit_group["q_fit_mmol_g"].to_numpy(dtype=float)

        if cursor != len(fits_data):
            raise ValueError(
                "D-A fitting output contains extra rows beyond those in grouped_data."
            )

        W0, E, n = getParaForDA_func(da_result.params)
        boundary_hits = []
        native_params = getattr(da_result, "params", None)
        if native_params is not None:
            for name, param in native_params.items():
                value = float(param.value)
                tol = 1e-7 * max(1.0, abs(value))
                at_min = (
                    np.isfinite(param.min)
                    and abs(value - float(param.min)) <= tol
                )
                at_max = (
                    np.isfinite(param.max)
                    and abs(value - float(param.max)) <= tol
                )
                if at_min or at_max:
                    boundary_hits.append(str(name))

        diagnostics = {
            "message": str(getattr(da_result, "message", "")),
            "nfev": getattr(da_result, "nfev", None),
            "chisqr": getattr(da_result, "chisqr", None),
            "validSheet": valid_sheet,
            "boundary_hits": boundary_hits,
        }

        return {
            "fit_result": da_result,
            "fitted_values": fitted_values,
            "parameters": {
                "W0_cm3_g": float(W0),
                "E_J_mol": float(E),
                "n": float(n),
            },
            "success": bool(getattr(da_result, "success", False)),
            "diagnostics": diagnostics,
        }

    return _fit_predict_func

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

def fit_DA0(T_list,  grouped_data, sheet = None, file = None,verbose=True, initial_params=None,):
        # 构建初始拟合参数，全局拟合策略
    # params0 = Parameters()
    # validSheet = False
    # if sheet: ##origin fitting data
    #     params = build_DAparams(params0)
    #     #params, validSheet = get_DAparams_from_origin(params1, T_list, sheet, file) #未开发
    # else: 
    #     params = build_DAparams(params0)

    validSheet = False

    if initial_params is None:
        params0 = Parameters()
        params = build_DAparams(params0)
    else:
        # Reuse the original fit's parameter constraints and values
        # as the starting point for a bootstrap refit.
        params = initial_params.copy()

    # 拟合 D-A
    daPara = minimize(residual_DA, params, args=(grouped_data,), method="leastsq", max_nfev=20000)
    # print(report_fit(daPara))

    if verbose:
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
    modelCIcheck = True
    if modelCIcheck:
        # 1. Create the adapter for the existing global D-A fitter.
        da_fit_predict = makeDAGlobalFitPredictFunc( T_list, fit_DA0_func=fit_DA0, getParaForDA_func=getParaForDA, )
        # 2. Run a small test first.
        da_bootstrap_test = redBoots.runResidualBootstrapAnalysis( grouped_data, da_fit_predict, response_col="q_mmol_g", 
                n_bootstrap=100, residual_centering="within_group", random_state=2026, )
        printResidualBootstrapSummary( da_bootstrap_test, 
                    group_labels={ 273.15: "0 °C (273.15 K)", 288.15: "15 °C (288.15 K)", 298.15: "25 °C (298.15 K)", }, 
                    response_unit="mmol/g", )
        da_ci = redBoots.runResidualBootstrapCIAnalysis( da_bootstrap_test, confidence_interval_func=dop.runConfidenceIntervalAnalysis, 
                                confidence_level=0.95, )

        printResidualBootstrapSummary( da_bootstrap_test, ci_result=da_ci, 
                    group_labels={ 273.15: "0 °C (273.15 K)", 288.15: "15 °C (288.15 K)", 298.15: "25 °C (298.15 K)", }, 
                    response_unit="mmol/g", )

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
        sheet_name = "CC-Hy-900-1-1"
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