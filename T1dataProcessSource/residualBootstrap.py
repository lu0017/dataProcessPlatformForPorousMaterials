from common import *
import T1dataProcessSource.dataOperation as dop
import T1fileSource.fileOperation as fl
def getBootstrapDistributions(bootstrap_result):
    """
    统一读取不同分析模块的 Bootstrap 分布。
    优先使用通用字段 quantity_distributions，
    同时兼容现有模型参数分布字段。
    """
    if bootstrap_result is None:
        return {}
    if not isinstance(bootstrap_result, Mapping):
        raise TypeError( "bootstrap_result must be a mapping." )
    supported_keys = (
        "quantity_distributions",
        "statistic_distributions",
        "parameter_distributions",  # 兼容 D-A / DSL
        "bootstrap_distributions",
    )
    for key in supported_keys:
        distributions = bootstrap_result.get(key)
        if isinstance(distributions, Mapping):
            return distributions
    return {}
def resolveBootstrapCIMetadata( ci_result, confidence_level=0.95, ci_method="percentile", ):
    """
    优先从 ci_result 读取 CI 元数据；
    缺失时使用调用者传入的备用值。
    """
    if ci_result is None:
        ci_result = {}
    # 优先使用分析结果中已记录的值
    stored_confidence_level = ci_result.get( "confidence_level" )
    stored_method = ci_result.get("method")
    resolved_confidence_level = (
        stored_confidence_level
        if stored_confidence_level is not None
        else confidence_level
    )
    resolved_method = ( stored_method if stored_method is not None else ci_method )
    return {
        "confidence_level": resolved_confidence_level,
        "ci_method": resolved_method,
    }
def buildBootstrapDistributionTable( bootstrap_result, metadata=None, ):
    """
    将任意标准化 Bootstrap 参数/统计量分布整理成长表。
    Expected input
    --------------
    bootstrap_result["parameter_distributions"]:
        {
            "E_J_mol": [...],
            "b1_Pa_inv_T298.15K": [...],
            "K": [...],
            "g": [...],
            "corr_K_g": [...],
        }
    Parameters
    ----------
    bootstrap_result : dict
        Bootstrap 分析返回的结果。
    metadata : dict or None
        任意附加信息，例如 Sample、Model、Analysis。
    Returns
    -------
    pd.DataFrame
        每行代表一个统计量的一次 Bootstrap 结果。
    """
    distributions = getBootstrapDistributions( bootstrap_result )
    if not distributions:
        raise ValueError( "No supported Bootstrap distributions found." )
    if not isinstance(distributions, dict) or not distributions:
        raise ValueError( "bootstrap_result must contain a non-empty parameter_distributions' dictionary." )
    metadata = dict(metadata or {})
    rows = []
    for quantity_name, values in distributions.items():
        values = np.asarray(values, dtype=float).reshape(-1)
        for index, value in enumerate(values, start=1):
            row = {
                **metadata,
                "Quantity": quantity_name,
                "Distribution_Index": index,
                "Bootstrap_Value": value,
                "Is_Finite": bool(np.isfinite(value)),
            }
            rows.append(row)
    return pd.DataFrame(rows)
def buildBootstrapCISummaryTable( bootstrap_result, ci_result, metadata=None, confidence_level=0.95, ci_method="percentile", ):
    """
    整理任意 Bootstrap CI 的估计值、标准误、
    置信区间与诊断信息。
    不负责计算 CI，只保存已计算的 CI。
    """
    ci_metadata = resolveBootstrapCIMetadata( ci_result, confidence_level=confidence_level, ci_method=ci_method, )
    confidence_level = ci_metadata["confidence_level"]
    ci_method = ci_metadata["ci_method"]
    statistics = ci_result.get("statistics")
    if not isinstance(statistics, dict) or not statistics:
        raise ValueError( "ci_result must contain a non-empty 'statistics' dictionary." )
    distributions = bootstrap_result.get( "parameter_distributions", {} )
    metadata = dict(metadata or {})
    n_successful = bootstrap_result.get( "n_successful_fits" )
    n_failed = bootstrap_result.get( "n_failed_fits" )
    failure_reasons = bootstrap_result.get( "failure_reasons" )
    if failure_reasons is not None:
        failure_reasons = json.dumps( failure_reasons, ensure_ascii=False, default=str, )
    rows = []
    for quantity_name, stat in statistics.items():
        ci_bounds = stat.get("confidence_interval")
        if ci_bounds is None:
            lower, upper = np.nan, np.nan
        else:
            lower, upper = ci_bounds
        values = np.asarray( distributions.get(quantity_name, []), dtype=float, ).reshape(-1)
        n_finite = int(np.isfinite(values).sum())
        n_total = int(values.size)
        n_valid = stat.get("n_valid", n_finite)
        n_invalid = stat.get( "n_invalid", max(n_total - n_finite, 0), )
        row = {
            **metadata,
            "Quantity": quantity_name,
            "Estimate": stat.get("estimate", np.nan),
            "Bootstrap_SE": stat.get( "bootstrap_standard_error", np.nan ),
            "CI_Lower": lower,
            "CI_Upper": upper,
            "Confidence_Level": confidence_level,
            "CI_Method": ci_method,
            "N_Valid": n_valid,
            "N_Invalid": n_invalid,
            "N_Distribution_Entries": n_total,
            "N_Successful_Fits": n_successful,
            "N_Failed_Fits": n_failed,
            "Failure_Reasons": failure_reasons,
        }
        rows.append(row)
    return pd.DataFrame(rows)
def exportBootstrapAnalysisResults(
    bootstrap_result=None,
    ci_result=None,
    file_name=None,
    metadata=None,
    confidence_level=0.95,
    ci_method="percentile",
    distribution_sheet_name="Bootstrap_Distribution",
    summary_sheet_name="Bootstrap_CI_Summary",
    extra_tables=None,
):
    """
    通用 Bootstrap 结果保存入口。
    支持：
    1. Bootstrap 分布 + CI 汇总；
    2. 仅 Bootstrap 分布；
    3. 仅 CI 汇总；
    4. 仅保存 extra_tables，例如样品间比较结果；
    5. 上述结果的任意组合。
    """
    # ==================================================
    # 1. 检查输出路径
    # ==================================================
    if not file_name:
        raise ValueError("file_name must be provided.")
    metadata = dict(metadata or {})
    extra_tables = dict(extra_tables or {})
    tables_to_export = {}
    # ==================================================
    # 2. 只有传入 Bootstrap 结果时，才整理分布
    # ==================================================
    if bootstrap_result is not None:
        distribution_df = buildBootstrapDistributionTable(
            bootstrap_result,
            metadata=metadata,
        )
        tables_to_export[distribution_sheet_name] = (
            distribution_df
        )
    # ==================================================
    # 3. 只有传入 CI 结果时，才整理 CI 汇总
    # ==================================================
    if ci_result is not None:
        # 允许仅保存 CI，不要求同时传入 Bootstrap 分布
        bootstrap_context = (
            bootstrap_result
            if bootstrap_result is not None
            else {}
        )
        summary_df = buildBootstrapCISummaryTable(
            bootstrap_context,
            ci_result,
            metadata=metadata,
            confidence_level=confidence_level,
            ci_method=ci_method,
        )
        tables_to_export[summary_sheet_name] = summary_df
    # ==================================================
    # 4. 加入任意附加结果表
    # ==================================================
    for sheet_name, table_df in extra_tables.items():
        if not isinstance(table_df, pd.DataFrame):
            raise TypeError(
                f"Extra table '{sheet_name}' must be a DataFrame."
            )
        if sheet_name in tables_to_export:
            raise ValueError(
                f"Duplicate sheet name: '{sheet_name}'"
            )
        tables_to_export[sheet_name] = table_df
    # ==================================================
    # 5. 确认至少有一张表需要保存
    # ==================================================
    if not tables_to_export:
        raise ValueError(
            "No results provided. Supply bootstrap_result, "
            "ci_result, or extra_tables."
        )
    # ==================================================
    # 6. 统一导出 Excel
    # ==================================================
    for sheet_name, table_df in tables_to_export.items():
        fl.export_to_excel_auto(
            table_df,
            filename=file_name,
            sheet_name=sheet_name,
            index=False,
        )
    return tables_to_export
FitPredictFunc = Callable[..., Mapping[str, Any]]
def _validateResidualBootstrapData(
    grouped_data: Mapping[Hashable, pd.DataFrame],
    response_col: str,
) -> list[Hashable]:
    """Validate grouped input data and return group keys in stable order."""
    if not isinstance(grouped_data, Mapping) or not grouped_data:
        raise TypeError("grouped_data must be a non-empty mapping of group keys to DataFrames.")
    group_keys = list(grouped_data.keys())
    for key in group_keys:
        df = grouped_data[key]
        if not isinstance(df, pd.DataFrame):
            raise TypeError(f"Group {key!r} must be a pandas DataFrame.")
        if response_col not in df.columns:
            raise KeyError(f"Group {key!r} is missing response column {response_col!r}.")
        if df.empty:
            raise ValueError(f"Group {key!r} contains no observations.")
        try:
            response = df[response_col].to_numpy(dtype=float)
        except (TypeError, ValueError) as exc:
            raise TypeError(
                f"Response column {response_col!r} in group {key!r} must be numeric."
            ) from exc
        if not np.isfinite(response).all():
            raise ValueError(
                f"Response column {response_col!r} in group {key!r} contains non-finite values."
            )
    return group_keys
def _normalizeResidualBootstrapFitOutput(
    fit_output: Mapping[str, Any],
    group_keys: list[Hashable],
    group_sizes: Mapping[Hashable, int],
    *,
    context: str,
) -> dict[str, Any]:
    """Validate and normalize the model-adapter output."""
    if not isinstance(fit_output, Mapping):
        raise TypeError(f"{context}: fit_predict_func must return a mapping.")
    required = {"success", "fitted_values", "parameters"}
    missing = required - set(fit_output)
    if missing:
        raise KeyError(f"{context}: fit output is missing keys {sorted(missing)}.")
    if not bool(fit_output["success"]):
        details = fit_output.get("diagnostics", {})
        raise RuntimeError(f"{context}: model fit reported failure. Diagnostics: {details!r}")
    fitted_mapping = fit_output["fitted_values"]
    if not isinstance(fitted_mapping, Mapping):
        raise TypeError(f"{context}: 'fitted_values' must be a mapping by group key.")
    missing_groups = set(group_keys) - set(fitted_mapping)
    extra_groups = set(fitted_mapping) - set(group_keys)
    if missing_groups or extra_groups:
        raise ValueError(
            f"{context}: fitted-value group keys mismatch; missing={missing_groups!r}, "
            f"extra={extra_groups!r}."
        )
    fitted_values: dict[Hashable, np.ndarray] = {}
    for key in group_keys:
        try:
            values = np.asarray(fitted_mapping[key], dtype=float).reshape(-1)
        except (TypeError, ValueError) as exc:
            raise TypeError(f"{context}: fitted values for group {key!r} must be numeric.") from exc
        if len(values) != group_sizes[key]:
            raise ValueError(
                f"{context}: fitted values for group {key!r} have length {len(values)}, "
                f"expected {group_sizes[key]}. Ensure predictions follow the input row order."
            )
        if not np.isfinite(values).all():
            raise ValueError(f"{context}: fitted values for group {key!r} contain non-finite values.")
        fitted_values[key] = values
    parameters_input = fit_output["parameters"]
    if not isinstance(parameters_input, Mapping) or not parameters_input:
        raise TypeError(f"{context}: 'parameters' must be a non-empty mapping of scalar parameters.")
    parameters: dict[str, float] = {}
    for name, value in parameters_input.items():
        if isinstance(value, (bool, np.bool_)):
            raise TypeError(f"{context}: parameter {name!r} must be numeric, not boolean.")
        try:
            scalar = float(value)
        except (TypeError, ValueError) as exc:
            raise TypeError(f"{context}: parameter {name!r} must be a scalar number.") from exc
        if not np.isfinite(scalar):
            raise ValueError(f"{context}: parameter {name!r} is not finite.")
        parameters[str(name)] = scalar
    return {
        "fit_result": fit_output.get("fit_result"),
        "fitted_values": fitted_values,
        "parameters": parameters,
        "success": True,
        "diagnostics": dict(fit_output.get("diagnostics", {}) or {}),
        "adapter_output": fit_output,
    }
def _centerResidualPools(
    residuals_by_group: Mapping[Hashable, np.ndarray],
    *,
    method: str,
) -> tuple[dict[Hashable, np.ndarray], dict[Hashable, dict[str, float]]]:
    """Center residual pools using the selected centering rule."""
    allowed = {"within_group", "global", "none"}
    if method not in allowed:
        raise ValueError(f"residual_centering must be one of {sorted(allowed)}; got {method!r}.")
    if method == "global":
        pooled = np.concatenate([np.asarray(v, dtype=float) for v in residuals_by_group.values()])
        global_mean = float(np.mean(pooled))
    else:
        global_mean = 0.0
    centered: dict[Hashable, np.ndarray] = {}
    diagnostics: dict[Hashable, dict[str, float]] = {}
    for key, residuals in residuals_by_group.items():
        residuals = np.asarray(residuals, dtype=float)
        mean_before = float(np.mean(residuals))
        if method == "within_group":
            center_value = mean_before
        elif method == "global":
            center_value = global_mean
        else:
            center_value = 0.0
        pool = residuals - center_value
        centered[key] = pool
        diagnostics[key] = {
            "n_points": int(len(residuals)),
            "mean_residual_before_centering": mean_before,
            "centering_value": float(center_value),
            "mean_residual_after_centering": float(np.mean(pool)),
            "residual_sd_before_centering": float(np.std(residuals, ddof=1)) if len(residuals) > 1 else 0.0,
            "residual_sd_after_centering": float(np.std(pool, ddof=1)) if len(pool) > 1 else 0.0,
            "residual_min_after_centering": float(np.min(pool)),
            "residual_max_after_centering": float(np.max(pool)),
        }
    return centered, diagnostics
def compareSampleBootstrapParameters(
    all_sample_results,
    parameter_name,
    sample_pairs,
    bootstrap_key,
    ci_key,
    confidence_level=0.95,
    n_resamples = 20_000,
):
    """
    比较不同样品的同一拟合参数，计算参数差值的 Bootstrap CI。
    Parameters
    ----------
    all_sample_results : dict
        每个样品对应的分析结果。
    parameter_name : str
        要比较的参数名称。
        DSL 示例: "b1_Pa_inv_T298.15K"
        D-A 示例:  "E_J_mol"
    sample_pairs : iterable
        需要比较的样品对。
    bootstrap_key : str
        Bootstrap 结果对应的键。
        DSL: "dsl_bootstrap"
        D-A: "da_bootstrap"
    ci_key : str
        单样品参数 CI 对应的键。
        DSL: "dsl_ci"
        D-A: "da_ci"
    confidence_level : float
        差值 CI 的置信水平。
    Returns
    -------
    pd.DataFrame
        每个样品对的参数估计值、差值、SE 和置信区间。
    """
    rows = []
    for sample_a, sample_b in sample_pairs:
        # ----------------------------------
        # 1. 检查样品是否存在
        # ----------------------------------
        if sample_a not in all_sample_results:
            print(f"Skip: sample not found: {sample_a}")
            continue
        if sample_b not in all_sample_results:
            print(f"Skip: sample not found: {sample_b}")
            continue
        result_a = all_sample_results[sample_a]
        result_b = all_sample_results[sample_b]
        # ----------------------------------
        # 2. 提取两个样品的 Bootstrap 结果
        # ----------------------------------
        bootstrap_a = result_a.get(bootstrap_key)
        bootstrap_b = result_b.get(bootstrap_key)
        ci_a = result_a.get(ci_key)
        ci_b = result_b.get(ci_key)
        if any(
            x is None
            for x in (bootstrap_a, bootstrap_b, ci_a, ci_b)
        ):
            print(
                f"Skip: bootstrap/CI unavailable for "
                f"{sample_a} vs {sample_b}"
            )
            continue
        distributions_a = bootstrap_a["parameter_distributions"]
        distributions_b = bootstrap_b["parameter_distributions"]
        statistics_a = ci_a["statistics"]
        statistics_b = ci_b["statistics"]
        # ----------------------------------
        # 3. 检查参数名称
        # ----------------------------------
        if (
            parameter_name not in distributions_a
            or parameter_name not in distributions_b
            or parameter_name not in statistics_a
            or parameter_name not in statistics_b
        ):
            print(
                f"Skip: parameter '{parameter_name}' unavailable "
                f"for {sample_a} vs {sample_b}"
            )
            continue
        # ----------------------------------
        # 4. 获取 Bootstrap 参数分布
        # ----------------------------------
        dist_a = np.asarray(
            distributions_a[parameter_name],
            dtype=float,
        ).ravel()
        dist_b = np.asarray(
            distributions_b[parameter_name],
            dtype=float,
        ).ravel()
        # 去除无效值
        dist_a = dist_a[np.isfinite(dist_a)]
        dist_b = dist_b[np.isfinite(dist_b)]
        if len(dist_a) < 2 or len(dist_b) < 2:
            print(
                f"Skip: insufficient valid bootstrap values "
                f"for {sample_a} vs {sample_b}"
            )
            continue
        # 单样品的原始参数估计值
        estimate_a = statistics_a[parameter_name]["estimate"]
        estimate_b = statistics_b[parameter_name]["estimate"]
        # ----------------------------------
        # 5. 构建样品间参数差值分布
        # ----------------------------------
        contrast = runBootstrapContrastAnalysis(
            dist_a,
            dist_b,
            estimate_a=estimate_a,
            estimate_b=estimate_b,
            mode="independent",
            n_resamples = n_resamples,
        )
        # ----------------------------------
        # 6. 计算差值的置信区间
        # ----------------------------------
        contrast_ci = dop.runConfidenceIntervalAnalysis(
            method="percentile",
            confidence_level=confidence_level,
            estimate=contrast["estimate"],
            bootstrap_distribution=contrast[
                "bootstrap_distribution"
            ],
        )
        lower = contrast_ci["lower_bound"]
        upper = contrast_ci["upper_bound"]
        # ----------------------------------
        # 7. 汇总结果
        # ----------------------------------
        rows.append({
            "Sample_A": sample_a,
            "Sample_B": sample_b,
            "Parameter": parameter_name,
            "Estimate_A": estimate_a,
            "Estimate_B": estimate_b,
            "Difference_A_minus_B": contrast["estimate"],
            "Bootstrap_SE": contrast[
                "bootstrap_standard_error"
            ],
            "CI_lower": lower,
            "CI_upper": upper,
            "CI_excludes_zero": (
                lower > 0 or upper < 0
            ),
            "N_bootstrap_A": len(dist_a),
            "N_bootstrap_B": len(dist_b),
        })
    return pd.DataFrame(rows)
def runBootstrapContrastAnalysis(
    distribution_a,
    distribution_b,
    *,
    estimate_a: float,
    estimate_b: float,
    mode: str = "paired",
    n_resamples: Optional[int] = 10000,
    random_state: Optional[int] = 50,
) -> Dict[str, Any]:
    #mode="paired"对应同一次拟合下的两个变量，如DSL全局拟合，固定q，b随温度变化，那么计算b_0C 和b_25C间的关系
    #mode="independent" 是不同样品间同一个变量，如样品S1 和S5的D-A的E
    """
    Generate a Bootstrap distribution for the difference between
    two parameter estimates.
    Parameters
    ----------
    distribution_a, distribution_b :
        Bootstrap distributions of two parameters.
    estimate_a, estimate_b : float
        Original-data parameter estimates.
    mode : {"paired", "independent"}
        "paired": 
            Preserve iteration-wise correspondence. Use when both
            parameters come from the same joint Bootstrap run.
        "independent":
            Independently sample from the two marginal Bootstrap
            distributions. Use when the two estimates come from
            independent Bootstrap analyses.
    n_resamples : int or None
        Number of differences generated in independent mode.
        Defaults to the larger valid distribution size.
        Not used in paired mode.
    random_state : int or None
        Random seed for independent mode.
    Returns
    -------
    dict
        Original difference estimate, Bootstrap difference
        distribution, Bootstrap SE, and diagnostics.
    """
    # 1. Validate original estimates.
    estimate_a = float(estimate_a)
    estimate_b = float(estimate_b)
    if not np.isfinite(estimate_a) or not np.isfinite(estimate_b):
        raise ValueError("Original parameter estimates must be finite.")
    # 2. Validate Bootstrap distributions.
    dist_a = np.asarray(distribution_a, dtype=float).reshape(-1)
    dist_b = np.asarray(distribution_b, dtype=float).reshape(-1)
    if dist_a.size < 2 or dist_b.size < 2:
        raise ValueError( "Each Bootstrap distribution must contain at least two values." )
    # 3. Construct the difference distribution.
    if mode == "paired":
        if dist_a.size != dist_b.size:
            raise ValueError( "Paired Bootstrap distributions must have equal lengths." )
        # Keep the original iteration-wise pairing.
        valid = np.isfinite(dist_a) & np.isfinite(dist_b)
        valid_a = dist_a[valid]
        valid_b = dist_b[valid]
        if len(valid_a) < 2:
            raise ValueError( "Fewer than two valid paired Bootstrap estimates remain." )
        difference_distribution = valid_a - valid_b
        n_input = int(dist_a.size)
        n_valid = int(len(difference_distribution))
    elif mode == "independent":
        # The two estimates come from independent analyses.
        valid_a = dist_a[np.isfinite(dist_a)]
        valid_b = dist_b[np.isfinite(dist_b)]
        if len(valid_a) < 2 or len(valid_b) < 2:
            raise ValueError( "Each distribution must contain at least two finite values." )
        # Validate the user-specified number of resamples
        if n_resamples is None:
            n_resamples = 10_000
        if ( isinstance(n_resamples, (bool, np.bool_))
            or not isinstance(n_resamples, (int, np.integer))
            or n_resamples < 2 ):
            raise ValueError("n_resamples must be an integer >= 2.")
        # Automatically scale the number of contrast resamples
        # based on the larger valid bootstrap distribution.
        base_n = max(len(valid_a), len(valid_b))
        n_resamples = max( n_resamples, int(np.ceil(base_n * 1.5)), )
        rng = np.random.default_rng(random_state)
        # Draw independently from each empirical distribution.
        index_a = rng.integers( 0, len(valid_a), size=n_resamples )
        index_b = rng.integers( 0, len(valid_b), size=n_resamples )
        difference_distribution = ( valid_a[index_a] - valid_b[index_b] )
        n_input = int(dist_a.size + dist_b.size)
        n_valid = int(n_resamples)
    else:
        raise ValueError( "mode must be either 'paired' or 'independent'." )
    # 4. Summarize the difference distribution.
    difference_estimate = estimate_a - estimate_b
    return {
        "estimate": float(difference_estimate),
        "bootstrap_distribution": difference_distribution,
        "bootstrap_standard_error": float( np.std(difference_distribution, ddof=1) ),
        "mode": mode,
        "n_valid": n_valid,
        "n_input": n_input,
        "random_state": random_state,
    }
def runResidualBootstrapCIAnalysis(
    bootstrap_result: Mapping[str, Any],
    *,
    confidence_interval_func: Callable[..., Mapping[str, Any]],
    confidence_level: float = 0.95,
    method: str = "percentile",
) -> dict[str, Any]:
    """Calculate parameter CIs from a completed residual-bootstrap result.
    独立的置信区间层，通过你现有的 runConfidenceIntervalAnalysis() 计算 CI
    This is deliberately separate from ``runResidualBootstrapAnalysis``.
    ``confidence_interval_func`` is expected to be the project's existing
    ``runConfidenceIntervalAnalysis`` implementation.
    """
    if not callable(confidence_interval_func):
        raise TypeError("confidence_interval_func must be callable.")
    if not 0.0 < confidence_level < 1.0:
        raise ValueError("confidence_level must be between 0 and 1.")
    required = {"original_parameters", "parameter_distributions"}
    missing = required - set(bootstrap_result)
    if missing:
        raise KeyError(f"bootstrap_result is missing keys {sorted(missing)}.")
    results: dict[str, Any] = {}
    for name, estimate in bootstrap_result["original_parameters"].items():
        distribution = np.asarray(
            bootstrap_result["parameter_distributions"][name], dtype=float
        )
        valid_distribution = distribution[np.isfinite(distribution)]
        if len(valid_distribution) < 2:
            raise RuntimeError(
                f"Parameter {name!r} has fewer than two valid bootstrap estimates."
            )
        ci = confidence_interval_func(
            method=method,
            confidence_level=confidence_level,
            estimate=float(estimate),
            bootstrap_distribution=valid_distribution,
        )
        results[name] = {
            "estimate": float(estimate),
            "bootstrap_standard_error": float(np.std(valid_distribution, ddof=1)),
            "confidence_interval": (
                float(ci["lower_bound"]),
                float(ci["upper_bound"]),
            ),
            "confidence_level": confidence_level,
            "n_valid": int(len(valid_distribution)),
            "n_invalid": int(len(distribution) - len(valid_distribution)),
        }
    return {
        "statistics": results,
        "n_bootstrap": int(bootstrap_result["n_bootstrap"]),
        "n_successful_fits": int(bootstrap_result["n_successful_fits"]),
        "n_failed_fits": int(bootstrap_result["n_failed_fits"]),
        "confidence_level": confidence_level,
        "method": f"residual_bootstrap_{method}_ci",
    }
def runResidualBootstrapAnalysis(
    grouped_data: Mapping[Hashable, pd.DataFrame],
    fit_predict_func: FitPredictFunc,
    *,
    response_col: str,
    n_bootstrap: int = 10000,
    residual_centering: str = "within_group",
    random_state: int | None = None,
    verbose: bool = False,
) -> dict[str, Any]:
    """Run a model-agnostic, grouped residual bootstrap.
    通用核心：残差提取、中心化、组内有放回抽样、模拟数据生成、重复拟合、参数分布保存
    Parameters
    ----------
    grouped_data
        Mapping such as ``{temperature: DataFrame}``. Each group is resampled
        separately; for global D-A fitting, the groups should normally be the
        individual temperature datasets. For a single-temperature DSL fit, a
        one-group mapping is appropriate.
    fit_predict_func
        Callback with signature::
            fit_predict_func(data, *, initial_fit=None, verbose=False) -> dict
        It must return at least::
            {
                "success": bool,
                "fit_result": <native optimizer result, optional>,
                "fitted_values": {group_key: 1D predictions aligned to input rows},
                "parameters": {parameter_name: scalar estimate},
                "diagnostics": {},  # optional
            }
        The callback should not mutate its input. If ``initial_fit`` is given,
        it is the normalized output from the original fit and can be used to
        warm-start the optimizer.
    response_col
        Name of the measured response column (e.g. ``"q_mmol_g"``).
    n_bootstrap
        Number of bootstrap refitting iterations.
    residual_centering
        ``"within_group"`` (recommended when groups represent temperatures),
        ``"global"``, or ``"none"``. Residuals are always sampled within their
        original groups in this first implementation.
    random_state
        Seed for reproducibility.
    verbose
        Passed to the model adapter for the original fit; bootstrap fits are
        always called with ``verbose=False``.
    Returns
    -------
    dict
        Original estimates/residuals, per-parameter bootstrap distributions
        (failed replicates are NaN), and convergence diagnostics. Confidence
        intervals are deliberately calculated by a separate function.
    """
    if not callable(fit_predict_func):
        raise TypeError("fit_predict_func must be callable.")
    if ( isinstance(n_bootstrap, (bool, np.bool_)) 
        or not isinstance(n_bootstrap, (int, np.integer)) 
        or n_bootstrap < 2 ):
        raise ValueError("n_bootstrap must be an integer >= 2.")
    group_keys = _validateResidualBootstrapData(grouped_data, response_col)
    group_sizes = {key: len(grouped_data[key]) for key in group_keys}
    # Work on defensive copies to protect the caller's original data.
    original_data = {key: grouped_data[key].copy(deep=True) for key in group_keys}
    # 1. Fit the original data and validate the adapter contract.
    original_raw = fit_predict_func( {key: df.copy(deep=True) for key, df in original_data.items()}, 
                                    initial_fit=None, verbose=verbose, )
    original_fit = _normalizeResidualBootstrapFitOutput( original_raw, group_keys, group_sizes, context="Original fit", )
    # 2. Calculate residual pools: residual = experimental - fitted.
    residuals_by_group: dict[Hashable, np.ndarray] = {}
    for key in group_keys:
        observed = original_data[key][response_col].to_numpy(dtype=float)
        predicted = original_fit["fitted_values"][key]
        residuals_by_group[key] = observed - predicted
    residual_pools, residual_diagnostics = _centerResidualPools( residuals_by_group, method=residual_centering, )
    # 3. Prepare distribution arrays, using NaN to retain failed iterations.
    original_parameters = dict(original_fit["parameters"])
    parameter_names = list(original_parameters.keys())
    distributions = { name: np.full(n_bootstrap, np.nan, dtype=float) for name in parameter_names }
    failure_log: list[dict[str, Any]] = []
    failure_reasons: Counter[str] = Counter()
    boundary_hit_counts: Counter[str] = Counter()
    n_negative_simulated_by_group: Counter[str] = Counter()
    n_simulated_points = 0
    rng = np.random.default_rng(random_state)
    # 4. Draw residuals within each group, generate y*, and refit.
    for iteration in range(n_bootstrap):
        bootstrap_data: dict[Hashable, pd.DataFrame] = {}
        for key in group_keys:
            pool = residual_pools[key]
            n_group = group_sizes[key]
            sampled_residuals = rng.choice(pool, size=n_group, replace=True) #有放回抽样
            simulated_response = original_fit["fitted_values"][key] + sampled_residuals
            copied_df = original_data[key].copy(deep=True)
            copied_df[response_col] = simulated_response #只替换response_col指定的列
            bootstrap_data[key] = copied_df
            n_simulated_points += len(simulated_response)
            # Diagnostic only; the engine does not clip negative responses.
            n_negative_simulated_by_group[str(key)] += int(np.sum(simulated_response < 0))
        try:
            fit_raw = fit_predict_func( bootstrap_data, initial_fit=original_fit, verbose=False, )
            fit_normalized = _normalizeResidualBootstrapFitOutput( fit_raw, group_keys, group_sizes, context=f"Bootstrap iteration {iteration}", )
            current_parameters = fit_normalized["parameters"]
            if set(current_parameters) != set(parameter_names):
                raise ValueError(
                    "Parameter names returned by a bootstrap fit differ from the original fit: "
                    f"expected {sorted(parameter_names)}, got {sorted(current_parameters)}."
                )
            for name in parameter_names:
                distributions[name][iteration] = current_parameters[name]
            diagnostics = fit_normalized.get("diagnostics", {})
            boundary_hits = diagnostics.get("boundary_hits", ()) if isinstance(diagnostics, Mapping) else ()
            if isinstance(boundary_hits, Mapping):
                for name, is_hit in boundary_hits.items():
                    if bool(is_hit):
                        boundary_hit_counts[str(name)] += 1
            else:
                for name in boundary_hits:
                    boundary_hit_counts[str(name)] += 1
        except Exception as exc:
            failure_reasons[type(exc).__name__] += 1
            # Keep a failure record for diagnostics; successful iterations are
            # represented by their parameter distributions.
            failure_log.append({
                "iteration": int(iteration),
                "exception_type": type(exc).__name__,
                "message": str(exc),
            })
    n_valid_by_parameter = { name: int(np.isfinite(values).sum()) for name, values in distributions.items() }
    # A replicate is successful only if every reported parameter was stored.
    successful_mask = np.ones(n_bootstrap, dtype=bool)
    for values in distributions.values():
        successful_mask &= np.isfinite(values)
    n_successful = int(successful_mask.sum())
    return {
        "original_fit": original_fit,
        "original_parameters": original_parameters,
        "original_fitted_values": { key: values.copy() for key, values in original_fit["fitted_values"].items() },
        "original_residuals": { key: values.copy() for key, values in residuals_by_group.items() },
        "centered_residual_pools": { key: values.copy() for key, values in residual_pools.items() },
        "parameter_distributions": distributions,
        "residual_diagnostics": residual_diagnostics,
        "n_bootstrap": int(n_bootstrap),
        "n_successful_fits": n_successful,
        "n_failed_fits": int(n_bootstrap - n_successful),
        "n_valid_by_parameter": n_valid_by_parameter,
        "failure_reasons": dict(failure_reasons),
        "failure_log": failure_log,
        "boundary_hit_counts": dict(boundary_hit_counts),
        "n_negative_simulated_responses_by_group": dict(n_negative_simulated_by_group),
        "n_simulated_points": int(n_simulated_points),
        "response_col": response_col,
        "group_keys": group_keys,
        "residual_centering": residual_centering,
        "random_state": random_state,
        "method": "grouped_residual_bootstrap",
    }