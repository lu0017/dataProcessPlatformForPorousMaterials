from common import *
import T1dataProcessSource.dataOperation as dop

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
    if (
        isinstance(n_bootstrap, (bool, np.bool_))
        or not isinstance(n_bootstrap, (int, np.integer))
        or n_bootstrap < 2
    ):
        raise ValueError("n_bootstrap must be an integer >= 2.")

    group_keys = _validateResidualBootstrapData(grouped_data, response_col)
    group_sizes = {key: len(grouped_data[key]) for key in group_keys}

    # Work on defensive copies to protect the caller's original data.
    original_data = {key: grouped_data[key].copy(deep=True) for key in group_keys}

    # 1. Fit the original data and validate the adapter contract.
    original_raw = fit_predict_func(
        {key: df.copy(deep=True) for key, df in original_data.items()},
        initial_fit=None,
        verbose=verbose,
    )
    original_fit = _normalizeResidualBootstrapFitOutput(
        original_raw,
        group_keys,
        group_sizes,
        context="Original fit",
    )

    # 2. Calculate residual pools: residual = experimental - fitted.
    residuals_by_group: dict[Hashable, np.ndarray] = {}
    for key in group_keys:
        observed = original_data[key][response_col].to_numpy(dtype=float)
        predicted = original_fit["fitted_values"][key]
        residuals_by_group[key] = observed - predicted

    residual_pools, residual_diagnostics = _centerResidualPools(
        residuals_by_group,
        method=residual_centering,
    )

    # 3. Prepare distribution arrays, using NaN to retain failed iterations.
    original_parameters = dict(original_fit["parameters"])
    parameter_names = list(original_parameters.keys())
    distributions = {
        name: np.full(n_bootstrap, np.nan, dtype=float)
        for name in parameter_names
    }
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
            sampled_residuals = rng.choice(pool, size=n_group, replace=True)
            simulated_response = original_fit["fitted_values"][key] + sampled_residuals

            copied_df = original_data[key].copy(deep=True)
            copied_df[response_col] = simulated_response
            bootstrap_data[key] = copied_df

            n_simulated_points += len(simulated_response)
            # Diagnostic only; the engine does not clip negative responses.
            n_negative_simulated_by_group[str(key)] += int(np.sum(simulated_response < 0))

        try:
            fit_raw = fit_predict_func(
                bootstrap_data,
                initial_fit=original_fit,
                verbose=False,
            )
            fit_normalized = _normalizeResidualBootstrapFitOutput(
                fit_raw,
                group_keys,
                group_sizes,
                context=f"Bootstrap iteration {iteration}",
            )

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

    n_valid_by_parameter = {
        name: int(np.isfinite(values).sum())
        for name, values in distributions.items()
    }
    # A replicate is successful only if every reported parameter was stored.
    successful_mask = np.ones(n_bootstrap, dtype=bool)
    for values in distributions.values():
        successful_mask &= np.isfinite(values)
    n_successful = int(successful_mask.sum())

    return {
        "original_fit": original_fit,
        "original_parameters": original_parameters,
        "original_fitted_values": {
            key: values.copy() for key, values in original_fit["fitted_values"].items()
        },
        "original_residuals": {
            key: values.copy() for key, values in residuals_by_group.items()
        },
        "centered_residual_pools": {
            key: values.copy() for key, values in residual_pools.items()
        },
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