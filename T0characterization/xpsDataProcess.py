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
# -----------------------------------------------------------------------------
# Configuration: compatible with the user's existing C1s peak-config structure
# -----------------------------------------------------------------------------
XPS_PEAK_CONFIG = {
    "C1s": {
        "x_range": (280.0, 295.0),
        "background": "shirley",
        "peaks": [
            {
                "name": "C-C/C=C",
                "center": 284.8,
                "center_bounds": (284.2, 285.3),
                "fwhm": 1.2,
                "fwhm_bounds": (0.5, 2.0),
                "eta": 0.3,
                "eta_bounds": (0.0, 1.0),
            },
            {
                "name": "C-O",
                "center": 286.0,
                "center_bounds": (285.3, 286.8),
                "fwhm": 1.3,
                "fwhm_bounds": (0.5, 2.0),
                "eta": 0.3,
                "eta_bounds": (0.0, 1.0),
            },
            {
                "name": "C=O",
                "center": 287.5,
                "center_bounds": (286.8, 288.3),
                "fwhm": 1.3,
                "fwhm_bounds": (0.5, 2.0),
                "eta": 0.3,
                "eta_bounds": (0.0, 1.0),
            },
            {
                "name": "O-C=O",
                "center": 289.0,
                "center_bounds": (288.2, 290.0),
                "fwhm": 1.4,
                "fwhm_bounds": (0.5, 2.2),
                "eta": 0.3,
                "eta_bounds": (0.0, 1.0),
            },
        ],
        # The 3-peak candidate uses the first three entries by default.
        # Override candidate_peak_names if a different chemically justified
        # three-peak subset is required.
        "candidate_peak_names": {
            3: ["C-C/C=C", "C-O", "C=O"],
            4: ["C-C/C=C", "C-O", "C=O", "O-C=O"],
        },
        "fit_eta": False,  # False preserves the old fixed-eta behaviour.
        "maxfev": 50000,
    },
    "O1s": {
        "x_range": (528.0, 536.0),
        "background": "shirley",
        "peaks": [
            {
                "name": "C=O",
                "center": 531.0,
                "center_bounds": (530.3, 531.8),
                "fwhm": 1.5,
                "fwhm_bounds": (0.8, 2.5),
                "eta": 0.3,
                "eta_bounds": (0.0, 1.0),
            },
            {
                "name": "C-O",
                "center": 532.5,
                "center_bounds": (531.5, 533.5),
                "fwhm": 1.5,
                "fwhm_bounds": (0.8, 2.5),
                "eta": 0.3,
                "eta_bounds": (0.0, 1.0),
            },
            {
                "name": "O-C=O",
                "center": 533.5,
                "center_bounds": (532.5, 535.0),
                "fwhm": 1.5,
                "fwhm_bounds": (0.8, 2.5),
                "eta": 0.3,
                "eta_bounds": (0.0, 1.0),
            },
        ],
    },
}
# -----------------------------------------------------------------------------
# Data preparation
# -----------------------------------------------------------------------------
def prepareXPSData(df, config, x="Binding Energy", y="Intensity"):
    """Clean, range-filter, sort and de-duplicate one XPS spectrum."""
    if x not in df.columns or y not in df.columns:
        raise KeyError(f"DataFrame must contain columns {x!r} and {y!r}.")
    work = df[[x, y]].copy()
    work[x] = pd.to_numeric(work[x], errors="coerce")
    work[y] = pd.to_numeric(work[y], errors="coerce")
    work = work.replace([np.inf, -np.inf], np.nan).dropna().sort_values(x)
    xmin, xmax = config.get("x_range", (None, None))
    if xmin is not None:
        work = work[work[x] >= xmin]
    if xmax is not None:
        work = work[work[x] <= xmax]
    work = work.groupby(x, as_index=False)[y].mean()
    X = work[x].to_numpy(dtype=float)
    Y = work[y].to_numpy(dtype=float)
    if len(X) < 10:
        raise ValueError(f"Insufficient data points after filtering: {len(X)}")
    if not np.all(np.isfinite(X)) or not np.all(np.isfinite(Y)):
        raise ValueError("Spectrum contains non-finite values after cleaning.")
    if np.any(np.diff(X) <= 0):
        raise ValueError("Binding energy must be strictly increasing.")
    return X, Y
# -----------------------------------------------------------------------------
# Background estimation
# -----------------------------------------------------------------------------
def shirleyBackground(X, Y, max_iter=200, tol=1e-5):
    """Iterative Shirley background estimate for a sorted spectrum.
    This is a practical iterative implementation. Confirm endpoint selection and
    convergence visually for each spectrum, especially for sloping baselines.
    """
    X = np.asarray(X, dtype=float)
    Y = np.asarray(Y, dtype=float)
    if X.ndim != 1 or Y.ndim != 1 or len(X) != len(Y):
        raise ValueError("X and Y must be one-dimensional arrays of equal length.")
    if len(X) < 3:
        raise ValueError("At least three points are required for Shirley background.")
    # Endpoint levels initialize the background; the integral is evaluated over X.
    background = np.linspace(Y[0], Y[-1], len(Y))
    y_min = min(float(Y[0]), float(Y[-1]))
    signal = np.maximum(Y - background, 0.0)
    for _ in range(max_iter):
        old = background.copy()
        signal = np.maximum(Y - background, 0.0)
        total = trapezoid(signal, X)
        if not np.isfinite(total) or total <= np.finfo(float).eps:
            return np.linspace(Y[0], Y[-1], len(Y))
        # Cumulative area from the low-binding-energy end, normalized to [0, 1].
        cumulative = np.concatenate(([0.0], np.cumsum((signal[1:] + signal[:-1]) * np.diff(X) / 2.0)))
        fraction = np.clip(cumulative / total, 0.0, 1.0)
        background = Y[0] + (Y[-1] - Y[0]) * fraction
        if np.max(np.abs(background - old)) <= tol * max(1.0, np.max(np.abs(old))):
            break
    return background
def estimateXPSBackground(X, Y, method="shirley"):
    """Dispatch background estimation according to config['background']."""
    method = str(method).lower().strip()
    if method == "shirley":
        return shirleyBackground(X, Y)
    if method == "linear":
        return np.linspace(Y[0], Y[-1], len(Y))
    if method == "none":
        return np.zeros_like(Y, dtype=float)
    raise ValueError(f"Unsupported background method: {method!r}")
# -----------------------------------------------------------------------------
# Peak line shapes: all profiles are normalized by area over the full line shape
# -----------------------------------------------------------------------------
def gaussianPeak(X, area, center, fwhm):
    X = np.asarray(X, dtype=float)
    sigma = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    return area / (sigma * np.sqrt(2.0 * np.pi)) * np.exp(-0.5 * ((X - center) / sigma) ** 2)
def lorentzianPeak(X, area, center, fwhm):
    X = np.asarray(X, dtype=float)
    gamma = fwhm / 2.0
    return area / np.pi * gamma / ((X - center) ** 2 + gamma**2)
def pseudoVoigtPeak(X, area, center, fwhm, eta):
    """Area-normalized pseudo-Voigt; eta is the Lorentzian fraction."""
    eta = float(eta)
    if not 0.0 <= eta <= 1.0:
        raise ValueError("eta must be between 0 and 1.")
    return (1.0 - eta) * gaussianPeak(X, area, center, fwhm) + eta * lorentzianPeak(X, area, center, fwhm)
def evaluatePeak(X, area, center, fwhm, shape="pseudo_voigt", eta=0.3):
    """Evaluate one peak for the requested candidate line shape."""
    shape = str(shape).lower().replace("-", "_").strip()
    if shape in {"gaussian", "g"}:
        return gaussianPeak(X, area, center, fwhm)
    if shape in {"lorentzian", "lorentz", "l"}:
        return lorentzianPeak(X, area, center, fwhm)
    if shape in {"pseudo_voigt", "pseudovoigt", "pv"}:
        return pseudoVoigtPeak(X, area, center, fwhm, eta)
    raise ValueError(f"Unsupported peak shape: {shape!r}")
# -----------------------------------------------------------------------------
# Model parameter setup and evaluation
# -----------------------------------------------------------------------------
def validatePeakConfig(peak, fit_eta=False):
    required = ("name", "center", "center_bounds")
    missing = [key for key in required if key not in peak]
    if missing:
        raise KeyError(f"Peak config missing required keys {missing}: {peak}")
    center = float(peak["center"])
    center_lo, center_hi = map(float, peak["center_bounds"])
    fwhm = float(peak.get("fwhm", 1.2))
    fwhm_lo, fwhm_hi = map(float, peak.get("fwhm_bounds", (0.5, 2.0)))
    eta = float(peak.get("eta_fixed", peak.get("eta", 0.3)))
    eta_lo, eta_hi = map(float, peak.get("eta_bounds", (0.0, 1.0)))
    if not center_lo < center_hi or not center_lo <= center <= center_hi:
        raise ValueError(f"Invalid initial center/bounds for peak {peak['name']!r}.")
    if not 0.0 < fwhm_lo < fwhm_hi or not fwhm_lo <= fwhm <= fwhm_hi:
        raise ValueError(f"Invalid initial FWHM/bounds for peak {peak['name']!r}.")
    if not 0.0 <= eta <= 1.0:
        raise ValueError(f"Eta for peak {peak['name']!r} must be between 0 and 1.")
    if not 0.0 <= eta_lo < eta_hi <= 1.0:
        raise ValueError(f"eta_bounds for peak {peak['name']!r} must lie within [0, 1].")
    if fit_eta and not eta_lo <= eta <= eta_hi:
        raise ValueError(f"Initial eta for peak {peak['name']!r} is outside eta_bounds.")
    return {
        "name": str(peak["name"]),
        "center": center,
        "center_bounds": (center_lo, center_hi),
        "fwhm": fwhm,
        "fwhm_bounds": (fwhm_lo, fwhm_hi),
        "eta": eta,
        "eta_bounds": (eta_lo, eta_hi),
    }
def prepareFitParameters(X, corrected, peak_configs, shape, fit_eta=False):
    """Build initial values, bounds, and parameter labels for one model fit."""
    signal_scale = max(float(np.max(corrected)), float(np.std(corrected)), 1e-8)
    p0, lower, upper, param_names = [], [], [], []
    normalized_configs = []
    is_pv = str(shape).lower().replace("-", "_") in {"pseudo_voigt", "pseudovoigt", "pv"}
    for peak in peak_configs:
        p = validatePeakConfig(peak, fit_eta=fit_eta and is_pv)
        normalized_configs.append(p)
        area0 = max(signal_scale * p["fwhm"], 1e-8)
        p0.extend([area0, p["center"], p["fwhm"]])
        lower.extend([0.0, p["center_bounds"][0], p["fwhm_bounds"][0]])
        upper.extend([np.inf, p["center_bounds"][1], p["fwhm_bounds"][1]])
        param_names.extend([f"{p['name']}.Area", f"{p['name']}.Center", f"{p['name']}.FWHM"])
        if is_pv and fit_eta:
            p0.append(p["eta"])
            lower.append(p["eta_bounds"][0])
            upper.append(p["eta_bounds"][1])
            param_names.append(f"{p['name']}.Eta")
    return (np.asarray(p0, dtype=float), np.asarray(lower, dtype=float),
            np.asarray(upper, dtype=float), param_names, normalized_configs)
def unpackPeakParameters(params, peak_configs, shape, fit_eta=False):
    """Convert a flat optimizer vector into per-peak parameter dictionaries."""
    shape_key = str(shape).lower().replace("-", "_")
    is_pv = shape_key in {"pseudo_voigt", "pseudovoigt", "pv"}
    stride = 4 if is_pv and fit_eta else 3
    unpacked = []
    for i, peak in enumerate(peak_configs):
        start = stride * i
        area, center, fwhm = params[start:start + 3]
        eta = params[start + 3] if stride == 4 else peak["eta"]
        unpacked.append({"name": peak["name"], "area_parameter": float(area),
                         "center": float(center), "fwhm": float(fwhm), "eta": float(eta)})
    return unpacked
def evaluateMultiPeakModel(X, params, peak_configs, shape="pseudo_voigt", fit_eta=False):
    """Sum the configured component profiles."""
    peak_params = unpackPeakParameters(params, peak_configs, shape, fit_eta=fit_eta)
    total = np.zeros_like(np.asarray(X, dtype=float))
    for p in peak_params:
        total += evaluatePeak(X, p["area_parameter"], p["center"], p["fwhm"], shape, p["eta"])
    return total
def curveFitMultiPeakAdapter(X, *params, peak_configs, shape="pseudo_voigt", fit_eta=False):
    """Top-level adapter matching scipy.optimize.curve_fit's callable API."""
    return evaluateMultiPeakModel(X, params, peak_configs, shape=shape, fit_eta=fit_eta)
# -----------------------------------------------------------------------------
# Fit-quality statistics and diagnostics
# -----------------------------------------------------------------------------
def calculateXPSMetrics(observed, fitted, residual):
    observed = np.asarray(observed, dtype=float)
    fitted = np.asarray(fitted, dtype=float)
    residual = np.asarray(residual, dtype=float)
    ss_res = float(np.sum(residual**2))
    ss_tot = float(np.sum((observed - np.mean(observed))**2))
    r2 = float(1.0 - ss_res / ss_tot) if ss_tot > 0 else np.nan
    rmse = float(np.sqrt(np.mean(residual**2)))
    mae = float(np.mean(np.abs(residual)))
    n = len(observed)
    return {"R2": r2, "RMSE": rmse, "MAE": mae, "SSE": ss_res, "n_points": n}
def calculateInformationCriteria(metrics, n_parameters):
    """AIC/BIC based on Gaussian residual likelihood; lower is preferred."""
    n = int(metrics["n_points"])
    sse = max(float(metrics["SSE"]), np.finfo(float).tiny)
    k = int(n_parameters)
    return {
        "AIC": float(n * np.log(sse / n) + 2 * k),
        "BIC": float(n * np.log(sse / n) + k * np.log(n)),
    }
def checkParameterBounds( params, lower, upper, param_names, relative_tol=1e-4, ):
    """Return bound diagnostics for every fitted parameter."""
    rows = []
    for value, lo, hi, name in zip( params, lower, upper, param_names ):
        value = float(value)
        lo = float(lo)
        hi = float(hi)
        scale = max(1.0, abs(value))
        near_lower = ( np.isfinite(lo) and abs(value - lo) <= relative_tol * max(scale, abs(lo)) )
        near_upper = ( np.isfinite(hi) and abs(value - hi) <= relative_tol * max(scale, abs(hi)) )
        if near_lower and near_upper:
            status = "both"
        elif near_lower:
            status = "lower"
        elif near_upper:
            status = "upper"
        else:
            status = "none"
        rows.append({
            "Parameter": name,
            "Value": value,
            "Lower bound": lo,
            "Upper bound": hi,
            "Distance to lower": value - lo if np.isfinite(lo) else np.nan,
            "Distance to upper": hi - value if np.isfinite(hi) else np.nan,
            "Bound status": status,
            "At bound": status != "none",
        })
    return pd.DataFrame(rows)
def printXPSResidualDiagnostics(X, observed, fitted, residual):
    """Print concise residual summaries; this is called only when requested."""
    residual = np.asarray(residual, dtype=float)
    print("[Residual diagnostics] "
          f"mean={np.mean(residual):.4g}; std={np.std(residual):.4g}; "
          f"min={np.min(residual):.4g}; max={np.max(residual):.4g}")
# -----------------------------------------------------------------------------
# Single candidate fit: one shape + one explicit peak list only
# -----------------------------------------------------------------------------
def fitXPSRegion(df, config, x="Binding Energy", y="Intensity", mode=None, verbose=False, *,
                 shape="pseudo_voigt", peak_configs=None, fit_eta=None,
                 model_name=None):
    """Fit exactly one XPS model. Does not compare models or print by default.
    Parameters
    ----------
    config : dict
        Region config (e.g. XPS_PEAK_CONFIG['C1s']).
    shape : {'gaussian', 'lorentzian', 'pseudo_voigt'}
        Line shape used for every component in this single fit.
    peak_configs : list[dict], optional
        Explicit peaks to fit. Defaults to config['peaks'].
    fit_eta : bool, optional
        If True, optimize eta for pseudo-Voigt components using eta_bounds.
        Defaults to config.get('fit_eta', False). Otherwise eta is fixed, as in
        the original implementation. eta_bounds are validated but not optimized.
    verbose : bool
        Print metrics and residual diagnostics only when True.
    """
    # Backward compatibility: the old function accepted mode="components".
    # Envelope fitting is now intentionally handled by the separate fitXPSEnvelope().
    if mode is not None and str(mode).lower().strip() != "components":
        raise ValueError("fitXPSRegion now fits one component model only; use fitXPSEnvelope() for envelope diagnostics.")
    if peak_configs is None:
        peak_configs = config.get("peaks", [])
    if not peak_configs:
        raise ValueError("A component fit requires a non-empty peak_configs list.")
    if fit_eta is None:
        fit_eta = bool(config.get("fit_eta", False))
    if shape not in {"gaussian", "lorentzian", "pseudo_voigt", "pseudo-voigt", "pseudovoigt", "pv"}:
        raise ValueError(f"Unsupported shape: {shape!r}")
    X, Y = prepareXPSData(df, config, x=x, y=y)
    background = estimateXPSBackground(X, Y, config.get("background", "shirley"))
    corrected = Y - background
    p0, lower, upper, param_names, normalized_peaks = prepareFitParameters(
        X, corrected, peak_configs, shape, fit_eta=fit_eta
    )
    fit_eta_effective = bool(fit_eta and str(shape).lower().replace("-", "_") in {"pseudo_voigt", "pseudovoigt", "pv"})
    # partial binds the model settings while the top-level adapter accepts *params.
    fit_callable = partial(curveFitMultiPeakAdapter, peak_configs=normalized_peaks,
                           shape=shape, fit_eta=fit_eta_effective)
    popt, pcov = curve_fit(fit_callable, X, corrected, p0=p0,
                           bounds=(lower, upper), maxfev=int(config.get("maxfev", 50000)))
    fitted_corrected = evaluateMultiPeakModel(X, popt, normalized_peaks, shape, fit_eta=fit_eta_effective)
    fitted = background + fitted_corrected
    residual = Y - fitted
    metrics = calculateXPSMetrics(Y, fitted, residual)
    information = calculateInformationCriteria(metrics, len(popt))
    parameterBounds = checkParameterBounds(popt, lower, upper, param_names)
    boundHits = parameterBounds.loc[ parameterBounds["At bound"] ].reset_index(drop=True)
    unpacked = unpackPeakParameters(popt, normalized_peaks, shape, fit_eta=fit_eta_effective)
    components = {}
    rows = []
    for p in unpacked:
        component = evaluatePeak(X, p["area_parameter"], p["center"], p["fwhm"], shape, p["eta"])
        components[p["name"]] = component
        area_in_window = float(trapezoid(component, X))
        rows.append({"Peak": p["name"], "Center (eV)": p["center"],
                     "FWHM (eV)": p["fwhm"], "Eta (Lorentzian fraction)": p["eta"],
                     "Area parameter": p["area_parameter"], "Area in region": area_in_window})
    peak_table = pd.DataFrame(rows)
    total_area = float(peak_table["Area in region"].sum())
    peak_table["Area (%)"] = peak_table["Area in region"] / total_area * 100.0 if total_area > 0 else np.nan
    result = {
        "model_name": model_name or f"{shape}_{len(normalized_peaks)}peaks",
        "shape": shape,
        "x": X,
        "raw": Y,
        "background": background,
        "corrected": corrected,
        "components": components,
        "fitted_corrected": fitted_corrected,
        "fitted": fitted,
        "residual": residual,
        "peaks": peak_table,
        "covariance": pcov,
        "parameters": popt,
        "parameter_names": param_names,
        "metrics": metrics,
        "information_criteria": information,
        "success": True,
        "fit_mode": "components",
        "diagnostic_only": False,
        "fit_range": (float(X[0]), float(X[-1])),
        "parameter_bounds": parameterBounds,
        "bound_hits": boundHits,
        "n_parameters": len(popt),
    }
    if verbose:
        print(f"[{result['model_name']}] R2={metrics['R2']:.5f}; "
              f"RMSE={metrics['RMSE']:.5g}; AIC={information['AIC']:.5g}; BIC={information['BIC']:.5g}")
        printXPSResidualDiagnostics(X, Y, fitted, residual)
    return result
# -----------------------------------------------------------------------------
# Candidate selection and comparison
# -----------------------------------------------------------------------------
def selectPeakConfigs(config, n_peaks):
    """Select a declared peak subset for a candidate model."""
    peak_configs = config.get("peaks", [])
    candidate_names = config.get("candidate_peak_names", {}).get(n_peaks)
    if candidate_names is None:
        if len(peak_configs) < n_peaks:
            raise ValueError(f"Config has only {len(peak_configs)} peaks, cannot fit {n_peaks} peaks.")
        return peak_configs[:n_peaks]
    by_name = {p["name"]: p for p in peak_configs}
    missing = [name for name in candidate_names if name not in by_name]
    if missing:
        raise KeyError(f"candidate_peak_names[{n_peaks}] refers to undefined peaks: {missing}")
    if len(candidate_names) != n_peaks:
        raise ValueError(f"candidate_peak_names[{n_peaks}] must contain exactly {n_peaks} names.")
    return [by_name[name] for name in candidate_names]
def buildXPSCandidateSpecs(config, peak_counts=(3, 4), shapes=("gaussian", "lorentzian", "pseudo_voigt")):
    """Create the six default candidates: three shapes x three/four peaks."""
    specs = []
    for shape in shapes:
        for n_peaks in peak_counts:
            peak_configs = selectPeakConfigs(config, n_peaks)
            specs.append({"model_name": f"{shape}_{n_peaks}peaks", "shape": shape,
                          "n_peaks": n_peaks, "peak_configs": peak_configs})
    return specs
def collectXPSParameterBounds(modelResult):
    """
    Extract parameter-bound diagnostics from one model result.
    Preference:
        1. parameter_bounds: all fitted parameters
        2. bound_hits: parameters close to bounds
    Returns
    -------
    pandas.DataFrame
    """
    parameterBounds = modelResult.get("parameter_bounds")
    if isinstance(parameterBounds, pd.DataFrame):
        return parameterBounds.copy()
    boundHits = modelResult.get("bound_hits")
    if isinstance(boundHits, pd.DataFrame):
        return boundHits.copy()
    return pd.DataFrame()
def summarizeXPSModelComparison(results):
    """Build a ranking table; lower AIC/BIC are preferred, not R2 alone."""
    rows = []
    for name, result in results.items():
        rows.append({
            "Model": name,
            "Shape": result["shape"],
            "N peaks": len(result["peaks"]),
            "R2": result["metrics"]["R2"],
            "RMSE": result["metrics"]["RMSE"],
            "MAE": result["metrics"]["MAE"],
            "AIC": result["information_criteria"]["AIC"],
            "BIC": result["information_criteria"]["BIC"],
            "Parameters": result["n_parameters"],
            "Bound hits": len(result["bound_hits"]),
        })
    table = pd.DataFrame(rows)
    if not table.empty:
        table = table.sort_values(["BIC", "AIC"], ascending=True).reset_index(drop=True)
        table.insert(0, "BIC rank", np.arange(1, len(table) + 1))
    return table
def compareXPSModels(df, config, x="Binding Energy", y="Intensity", *,
                     candidate_specs=None, peak_counts=(3, 4),
                     shapes=("gaussian", "lorentzian", "pseudo_voigt"),
                     fit_eta=None, verbose=False, continue_on_error=True):
    """High-level runner: fit candidate models one by one and summarize results.
    Returns a dict with per-model results, a comparison table, and failures.
    Each candidate is fitted independently by fitXPSRegion().
    """
    if candidate_specs is None:
        candidate_specs = buildXPSCandidateSpecs(config, peak_counts=peak_counts, shapes=shapes)
    results, failures = {}, {}
    for spec in candidate_specs:
        name = spec.get("model_name", f"{spec['shape']}_{spec.get('n_peaks', len(spec['peak_configs']))}peaks")
        try:
            result = fitXPSRegion(
                df, config, x=x, y=y,
                shape=spec["shape"],
                peak_configs=spec["peak_configs"],
                fit_eta=fit_eta,
                model_name=name,
                verbose=verbose,
            )
            results[name] = result
        except Exception as exc:
            failures[name] = {"error_type": type(exc).__name__, "message": str(exc)}
            if not continue_on_error:
                raise
    summary = summarizeXPSModelComparison(results)
    return {"results": results, "summary": summary, "failures": failures}
# -----------------------------------------------------------------------------
# Optional envelope diagnostic, kept separate from component-model comparison
# -----------------------------------------------------------------------------
def asymmetricPseudoVoigtEnvelope(X, height, center, width_left, width_right, eta=0.3):
    """Asymmetric pseudo-Voigt envelope parameterized by peak height."""
    X = np.asarray(X, dtype=float)
    widths = np.where(X < center, width_left, width_right)
    widths = np.maximum(widths, np.finfo(float).eps)
    # Convert a height-like amplitude to a unit-height Gaussian/Lorentzian blend.
    z = (X - center) / widths
    gaussian = np.exp(-4.0 * np.log(2.0) * z**2)
    lorentzian = 1.0 / (1.0 + 4.0 * z**2)
    return height * ((1.0 - eta) * gaussian + eta * lorentzian)
def curveFitEnvelopeAdapter(X, height, center, width_left, width_right, eta=0.3):
    """Top-level scipy adapter for the asymmetric envelope fit."""
    return asymmetricPseudoVoigtEnvelope(X, height, center, width_left, width_right, eta=eta)
def fitXPSData( data, peakConfig, regions=("C1s",), x="Binding Energy", y="Intensity", modelSelection=True, shape="pseudo_voigt", fit_eta=None, verbose=False, ):
    """
    Batch-fit XPS spectra for multiple samples and regions.
    Parameters
    ----------
    data : dict
        {sample: {region: DataFrame}}
        Output from extractXPSRegions().
    peakConfig : dict
        Region-specific fitting configuration, e.g. XPS_PEAK_CONFIG.
    regions : tuple or list
        Regions to fit.
    modelSelection : bool
        True: compare the six candidate models.
        False: fit only one specified model.
    shape : str
        Used only when modelSelection=False.
    Returns
    -------
    dict
        {sample: {region: fitting_result}}
    """
    fitResults = {}
    for sample, sampleData in data.items():
        fitResults[sample] = {}
        for region in regions:
            if region not in peakConfig:
                fitResults[sample][region] = { "success": False, "error": f"No peak configuration found for {region}.", }
                continue
            df = sampleData.get(region)
            if not isinstance(df, pd.DataFrame) or df.empty:
                fitResults[sample][region] = { "success": False, "error": f"Missing or empty spectrum for {region}.", }
                continue
            config = peakConfig[region]
            try:
                if modelSelection:
                    # Compare Gaussian/Lorentzian/Pseudo-Voigt,
                    # each with 3-peak and 4-peak candidates.
                    result = compareXPSModels( df=df, config=config, x=x, y=y, fit_eta=fit_eta, verbose=verbose, continue_on_error=True, )
                    result["success"] = len(result["results"]) > 0
                else:
                    # Fit only one specified model.
                    result = fitXPSRegion( df=df, config=config, x=x, y=y, shape=shape, fit_eta=fit_eta, verbose=verbose, )
                fitResults[sample][region] = result
            except Exception as exc:
                fitResults[sample][region] = { "success": False, "error": f"{type(exc).__name__}: {exc}", }
    return fitResults
def collectXPSModelSummaries(fitResults, debug=True):
    """Collect XPS model summaries, peak parameters, bound diagnostics, and fit metrics."""
    rows = []
    for sample, sampleResults in fitResults.items():
        for region, regionResult in sampleResults.items():
            modelResults = regionResult["results"] if isinstance(regionResult, dict) and isinstance(regionResult.get("results"), dict) else ({regionResult.get("model_name", "model"): regionResult} if isinstance(regionResult, dict) else {})
            for modelName, modelResult in modelResults.items():
                if not isinstance(modelResult, dict):
                    continue
                metrics, info = modelResult.get("metrics", {}), modelResult.get("information_criteria", {})
                parameterBounds, boundHitsDf = modelResult.get("parameter_bounds"), modelResult.get("bound_hits")
                if isinstance(parameterBounds, pd.DataFrame):
                    selectedBounds, boundsSource = parameterBounds.copy(), "parameter_bounds"
                elif isinstance(boundHitsDf, pd.DataFrame):
                    selectedBounds, boundsSource = boundHitsDf.copy(), "bound_hits (fallback)"
                else:
                    selectedBounds, boundsSource = pd.DataFrame(), "unavailable"
                if not selectedBounds.empty:
                    boundRows = selectedBounds.loc[selectedBounds["At bound"].fillna(False).astype(bool)] if "At bound" in selectedBounds.columns else selectedBounds
                    boundHits = len(boundRows)
                    boundParameters = boundRows["Parameter"].astype(str).tolist() if "Parameter" in boundRows.columns else []
                    parameterBoundRecords = selectedBounds.to_dict(orient="records")
                else:
                    boundHits, boundParameters, parameterBoundRecords = 0, [], []
                peaks = modelResult.get("peaks", [])
                if isinstance(peaks, pd.DataFrame):
                    peaksDf = peaks.copy()
                elif isinstance(peaks, dict):
                    peaksDf = pd.DataFrame(peaks if all(isinstance(v, (dict, list, tuple, np.ndarray, pd.Series)) for v in peaks.values()) else [peaks])
                elif isinstance(peaks, (list, tuple)):
                    peaksDf = pd.DataFrame([p for p in peaks if isinstance(p, dict)])
                else:
                    peaksDf = pd.DataFrame()
                if not peaksDf.empty:
                    areaCol = next((c for c in peaksDf.columns if str(c).lower() in ("area", "peak area", "integrated area")), None)
                    if areaCol:
                        areas = pd.to_numeric(peaksDf[areaCol], errors="coerce")
                        totalArea = areas.sum()
                        peaksDf["Area (%)"] = areas / totalArea * 100 if totalArea > 0 else np.nan
                row = {
                    "Sample": sample, "Region": region, "Model": modelName,
                    "Success": modelResult.get("success", True),
                    "R2": metrics.get("R2", np.nan), "RMSE": metrics.get("RMSE", np.nan),
                    "MAE": metrics.get("MAE", np.nan), "AIC": info.get("AIC", np.nan),
                    "BIC": info.get("BIC", np.nan), "N peaks": len(peaksDf) if not peaksDf.empty else len(peaks) if isinstance(peaks, (list, tuple)) else 0,
                    "Bound hits": boundHits, "Bound parameters": ", ".join(boundParameters),
                    "Parameter bounds": parameterBoundRecords, "Error": modelResult.get("error"),
                }
                rows.append(row)
                if debug:
                    print(f"\n{'=' * 90}\n[XPS DEBUG] Sample={sample} | Region={region} | Model={modelName}\n{'-' * 90}")
                    print("Available result keys:", list(modelResult.keys()))
                    print(f"Fit metrics: R2={row['R2']:.5f} | RMSE={row['RMSE']:.4f} | MAE={row['MAE']:.4f} | AIC={row['AIC']:.3f} | BIC={row['BIC']:.3f}")
                    print(f"Parameter-bound source: {boundsSource}")
                    print(selectedBounds.to_string(index=False) if not selectedBounds.empty else "No parameter-bound information saved.")
                    print(f"Bound-hit count: {boundHits} | Bound-hit parameters: {boundParameters}")
                    print("Peak parameters:")
                    print(peaksDf.to_string(index=False, float_format=lambda v: f"{v:.4f}") if not peaksDf.empty else "No peak-parameter table found in modelResult['peaks'].")
                    print("=" * 90)
    summary = pd.DataFrame(rows)
    if debug and not summary.empty:
        displayColumns = ["Sample", "Region", "Model", "Success", "R2", "RMSE", "MAE", "AIC", "BIC", "N peaks", "Bound hits"]
        print("\n=== XPS Model Comparison Summary ===")
        print(summary[displayColumns].to_string(index=False, float_format=lambda value: f"{value:.5f}"))
    return summary
def fitXPSEnvelope(df, config, x="Binding Energy", y="Intensity", verbose=False):
    """Fit a single asymmetric envelope for diagnosis only, not chemical areas."""
    if "envelope_range" not in config:
        raise ValueError("Envelope mode requires config['envelope_range'].")
    X_all, Y_all = prepareXPSData(df, config, x=x, y=y)
    B_all = estimateXPSBackground(X_all, Y_all, config.get("background", "shirley"))
    C_all = Y_all - B_all
    env_min, env_max = config["envelope_range"]
    mask = (X_all >= env_min) & (X_all <= env_max)
    if np.count_nonzero(mask) < 10:
        raise ValueError("Insufficient data points in envelope_range.")
    X, Y, B, C = X_all[mask], Y_all[mask], B_all[mask], C_all[mask]
    idx = int(np.argmax(C))
    center0 = float(X[idx])
    center_lo, center_hi = config.get("envelope_center_bounds", (max(X[0], center0 - 0.8), min(X[-1], center0 + 0.8)))
    width_lo, width_hi = config.get("envelope_width_bounds", (0.5, 4.0))
    eta = float(config.get("envelope_eta_fixed", 0.3))
    if not 0.0 <= eta <= 1.0:
        raise ValueError("envelope_eta_fixed must be between 0 and 1.")
    p0 = [max(float(np.max(C)), 1e-8), float(np.clip(center0, center_lo, center_hi)),
          float(np.clip(1.5, width_lo + 1e-6, width_hi - 1e-6)),
          float(np.clip(1.5, width_lo + 1e-6, width_hi - 1e-6))]
    envelope_callable = partial(curveFitEnvelopeAdapter, eta=eta)
    lower = [0.0, center_lo, width_lo, width_lo]
    upper = [np.inf, center_hi, width_hi, width_hi]
    popt, pcov = curve_fit(envelope_callable, X, C, p0=p0, bounds=(lower, upper),
                           maxfev=int(config.get("maxfev", 50000)))
    fitted_corrected = envelope_callable(X, *popt)
    fitted = B + fitted_corrected
    residual = Y - fitted
    metrics = calculateXPSMetrics(Y, fitted, residual)
    height, center, width_left, width_right = popt
    peaks = pd.DataFrame([{"Feature": "Main envelope (diagnostic only)",
                           "Center (eV)": center, "Width left (eV)": width_left,
                           "Width right (eV)": width_right, "Eta fixed": eta,
                           "Peak height": height,
                           "Envelope area in fit window": float(trapezoid(fitted_corrected, X))}])
    if verbose:
        print(f"[Envelope diagnostic] R2={metrics['R2']:.5f}; RMSE={metrics['RMSE']:.5g}")
        printXPSResidualDiagnostics(X, Y, fitted, residual)
    return {"x": X, "raw": Y, "background": B, "corrected": C,
            "components": {"Main envelope (diagnostic only)": fitted_corrected},
            "fitted_corrected": fitted_corrected, "fitted": fitted,
            "residual": residual, "peaks": peaks, "covariance": pcov,
            "metrics": metrics, "success": True, "fit_mode": "envelope",
            "diagnostic_only": True, "fit_range": (float(X[0]), float(X[-1]))}
def saveXPS(xpsData, filename="XPS.xlsx"):
    """
    Save XPS data.
    One Region -> one sheet
    xpsData structure
    -----------------
    {
        Sample:{
            metadata:{...},
            Region: DataFrame,
            ...
        }
    }
    Excel format
    ------------
        Sample1                    Sample2
        Binding Energy Intensity   Binding Energy Intensity
        ...
    """
    # =========================
    # Get all regions
    # =========================
    regions = set()
    for sampleData in xpsData.values():
        regions.update( key for key, value in sampleData.items() if isinstance(value, pd.DataFrame) )
    regions = sorted(regions)
    # =========================
    # Export each region
    # =========================
    for region in regions:
        data_blocks = []
        header_row1 = []
        header_row2 = []
        # ---------------------
        # Collect every sample
        # ---------------------
        for sample, sampleData in xpsData.items():
            if region not in sampleData:
                continue
            df = sampleData[region]
            if not isinstance(df, pd.DataFrame):
                continue
            temp = df[ [ "Binding Energy", "Intensity", ] ].reset_index(drop=True)
            data_blocks.append(temp)
            header_row1.extend( [ sample, sample, ] )
            header_row2.extend( [ "Binding Energy", "Intensity", ] )
        if len(data_blocks) == 0:
            continue
        # ---------------------
        # Horizontal concatenate
        # ---------------------
        output = pd.concat( data_blocks, axis=1, )
        # ---------------------
        # Two header rows
        # ---------------------
        header = pd.DataFrame( [ header_row1, header_row2, ], columns=output.columns, )
        final = pd.concat( [ header, output, ], ignore_index=True, )
        # ---------------------
        # Export
        # ---------------------
        fl.export_to_excel_auto( final, filename=filename, sheet_name=region, )
        print(f"Export {region} -> {filename}")
def readXPS( folder, pattern="*.csv", keyName=None, parser_kwargs=None, ):
    """
    Batch read XPS spectra from CSV files.
    Returns
    -------
    dict
    {
        sample_name: {
            "metadata": metadata,
            "C1s": dataframe,
            "O1s": dataframe,
            ...
        }
    }
    """
    parser_kwargs = parser_kwargs or {}
    region_map = {
        "su": "Survey",
        "survey": "Survey",
        "c": "C1s",
        "c1s": "C1s",
        "o": "O1s",
        "o1s": "O1s",
        "n": "N1s",
        "n1s": "N1s",
        "s": "S2p",
        "s2p": "S2p",
        "p": "P2p",
        "p2p": "P2p",
    }
    xps = {}
    folder = Path(folder)
    for file in folder.glob(pattern):
        result = fl.parseGroupedTable(
            file=file,
            metadata_rows=(0, 2),
            metadata_names=("Sample", "Original"),
            data_start_row=4,
            group_size=2,
            data_names=("Binding Energy", "Intensity"),
            **parser_kwargs,
        )
        # Use CSV filename as the top-level key,
        # instead of the internal sample metadata such as "6".
        sample_name = file.stem
        metadata = result["metadata"]
        grouped_data = result["data"]
        if sample_name not in xps:
            xps[sample_name] = { "metadata": metadata, }
        for _, regions in grouped_data.items():
            for region_name, df in regions.items():
                region_key = str(region_name).strip().lower()
                standardized_region = region_map.get( region_key, str(region_name).strip(), )
                xps[sample_name][standardized_region] = df
    return xps
def plotXPSFit(fit_result, title=None, figsize=(7, 5), savepath=None, dpi=600):
    """Plot XPS spectrum, background, fitted components, and residual."""
    if not fit_result.get("success", False):
        raise ValueError("Cannot plot an unsuccessful fit.")
    X, Y = fit_result["x"], fit_result["raw"]
    background, fitted = fit_result["background"], fit_result["fitted"]
    fig, (ax, ax_res) = plt.subplots(2, 1, figsize=figsize, sharex=True,
                                     gridspec_kw={"height_ratios": [3, 1]})
    ax.scatter(X, Y, s=8, label="Data", zorder=3)
    ax.plot(X, background, "--", lw=1.2, label="Background")
    for name, component in fit_result["components"].items():
        ax.plot(X, background + component, lw=1.1, label=name)
    ax.plot(X, fitted, lw=2, label="Total fit")
    ax.set_ylabel("Intensity (a.u.)")
    ax.legend(fontsize=8, frameon=False)
    ax_res.axhline(0, color="gray", ls="--", lw=0.8)
    ax_res.plot(X, fit_result["residual"], lw=1)
    ax_res.set_xlabel("Binding Energy (eV)")
    ax_res.set_ylabel("Residual")
    ax_res.set_xlim(max(X), min(X))
    if title:
        ax.set_title(title)
    fig.tight_layout()
    if savepath:
        fig.savefig(savepath, dpi=dpi, bbox_inches="tight")
    plt.show(block=False)
    return fig, (ax, ax_res)
def plotSpectrumXPS(
    data,
    regions=None,
    x="Binding Energy",
    y="Intensity",
    offset=0,
    savepath=None,
    **kwargs,
):
    """
    Plot nested XPS spectra by sample and region.
    Parameters
    ----------
    data : dict
        {
            sampleName: {
                "metadata": metadata,
                "Survey": DataFrame,
                "C1s": DataFrame,
                "O1s": DataFrame,
                ...
            }
        }
    regions : list or str, optional
        Regions to plot, e.g. ["C1s", "O1s"].
        If None, plot all available DataFrame regions.
    x, y : str
        X and Y column names.
    savepath : str or Path, optional
        Output directory or filename prefix.
    **kwargs
        Other parameters passed to plotSpectrum().
    Returns
    -------
    dict
        {
            region: axes
        }
    """
    results = {}
    # Collect available regions
    if regions is None:
        regions = []
        for sample_data in data.values():
            for region, value in sample_data.items():
                if isinstance(value, pd.DataFrame) and region not in regions:
                    regions.append(region)
    elif isinstance(regions, str):
        regions = [regions]
    for region in regions:
        region_data = {
            sample: sample_data[region]
            for sample, sample_data in data.items()
            if region in sample_data
            and isinstance(sample_data[region], pd.DataFrame)
        }
        if not region_data:
            continue
        region_savepath = None
        if savepath is not None:
            savepath = Path(savepath)
            region_savepath = savepath.parent / f"{savepath.stem}_{region}{savepath.suffix}"
        results[region] = myPlt.plotSpectrum(
            region_data,
            x=x,
            y=y,
            xlabel="Binding Energy (eV)",
            ylabel="Intensity (a.u.)",
            reverse_x=True,
            offset=offset,
            savepath=region_savepath,
            **kwargs,
        )
    plt.show(block=False)
    return results
def extractXPSRegions(
    data,
    highResolutionFlag=False,
    regions=("C1s", "O1s", "N1s", "K2p"),
    energyRanges=None,
    surveyKey="Survey",
    x="Binding Energy",
):
    """
    Extract XPS regions from Survey spectra or use
    directly measured high-resolution spectra.
    Parameters
    ----------
    data : dict
        If highResolutionFlag=False:
            {
                sample: {
                    "metadata": ...,
                    "Survey": DataFrame
                }
            }
        If highResolutionFlag=True:
            {
                sample: {
                    "metadata": ...,
                    "C1s": DataFrame,
                    "O1s": DataFrame,
                    ...
                }
            }
    highResolutionFlag : bool
        False: extract regions from Survey.
        True: use directly measured high-resolution spectra.
    regions : tuple or list
        Target regions.
    energyRanges : dict, optional
        Binding-energy ranges for Survey extraction.
    surveyKey : str
        Key for the Survey spectrum.
    x : str
        Binding-energy column name.
    Returns
    -------
    dict
        {
            sample: {
                "C1s": DataFrame,
                "O1s": DataFrame,
                ...
            }
        }
    """
    if energyRanges is None:
        energyRanges = {
            "C1s": (280, 295),
            "O1s": (525, 540),
            "N1s": (395, 410),
            "K2p": (290, 305),
        }
    result = {}
    for sample, sampleData in data.items():
        result[sample] = {}
        for region in regions:
            # ------------------------------------------
            # Mode 1: Directly measured high-resolution
            # ------------------------------------------
            if highResolutionFlag:
                df = sampleData.get(region)
                if isinstance(df, pd.DataFrame):
                    result[sample][region] = df.copy()
            # ------------------------------------------
            # Mode 2: Extract from Survey
            # ------------------------------------------
            else:
                survey = sampleData.get(surveyKey)
                if not isinstance(survey, pd.DataFrame):
                    continue
                if region not in energyRanges:
                    continue
                xmin, xmax = energyRanges[region]
                mask = survey[x].between(xmin, xmax)
                result[sample][region] = ( survey.loc[mask].copy().reset_index(drop=True) )
    return result
def cropXPSData(data, xmin=None, xmax=None):
    """
    Crop XPS data for each sample and spectrum region.
    Parameters
    ----------
    data : dict
        {
            sampleName: {
                "metadata": metadata,
                "C1s": DataFrame,
                "O1s": DataFrame,
                ...
            }
        }
    xmin, xmax : float, optional
        Minimum and maximum X values.
    Returns
    -------
    dict
        Same nested structure, with DataFrames cropped.
    """
    result = {}
    for sample, sample_data in data.items():
        result[sample] = {}
        for region, value in sample_data.items():
            if isinstance(value, pd.DataFrame):
                cropped = dop.cropX(
                    {region: value},
                    xmin=xmin,
                    xmax=xmax,
                )
                result[sample][region] = cropped[region]
            else:
                # Preserve metadata and other non-DataFrame values
                result[sample][region] = value
    return result
def baselineCorrectionXPS( data, x="Binding Energy", y="Intensity", method="asls", ):
    """Baseline correction for nested XPS spectra."""
    result = {}
    for sample, sample_data in data.items():
        result[sample] = {}
        for region, value in sample_data.items():
            if isinstance(value, pd.DataFrame):
                corrected = dop.baselineCorrection(
                    {sample: value},
                    x=x,
                    y=y,
                    method=method,
                )
                result[sample][region] = corrected[sample]
            else:
                result[sample][region] = value
    return result
def normalizeSpectrumXPS( data, y="Intensity", method="vector", ):
    """Normalize nested XPS spectra."""
    result = {}
    for sample, sample_data in data.items():
        result[sample] = {}
        for region, value in sample_data.items():
            if isinstance(value, pd.DataFrame):
                normalized = dop.normalizeSpectrum( {sample: value}, y=y, method=method, )
                result[sample][region] = normalized[sample]
            else:
                result[sample][region] = value
    return result
def smoothXPSSpectrum(data, y="Intensity", regions=("Survey",), method="savgol", **kwargs):
    result = {}
    for sample, item in data.items():
        newitem = item.copy()
        for region in regions:
            if region not in item:
                continue
            df = item[region].copy()
            if y not in df.columns:
                raise KeyError(f"{sample}/{region}: missing '{y}', columns={df.columns.tolist()}")
            df[y] = dop.smooth(df[y].to_numpy(), method=method, **kwargs)
            newitem[region] = df
        result[sample] = newitem
    return result
def main(file_path=None):
    spectrum = "XPS"
    saveFile = True
    singleFile = False
    highResolutionFlag = False #False时将自动从Suvey中提取对应的区段谱
    charaData = {}
    if spectrum == "XPS":
        ########  FTIR 文件的设置  #########
        extension = ".csv"
        fileName = "XPS"
        keyFlag = False
        keyName = ["Binding Energy", "Intensity"]
        ncols = 2
        sampleMap = {
        "S1": ["600-05"],
        "S2": ["600-1"],
        "S3": ["700-1"],
        "S4": ["800-05"],
        "S5": ["800-1"],
        "S6": ["600-2"],
        "S7": ["900-1"],
        "S8": ["800-2"],
        "S9": ["600-4"],
        "S10": ["800-4"],
        }
    if singleFile:
        # 选择一个文件
        file = fl.getFile()
        sample = os.path.splitext(os.path.basename(file))[0]
        df = fl.readData( file=file, keyFlag=keyFlag, keyName=keyName, ncols=ncols, skiprows=4,)
        if not df.empty:
            charaData[sample] = df
    else:
        # 选择文件夹
        folder = fl.select_folder()
        out_path = fl.get_expanded_name(folder, fileName = fileName)
        charaData = readXPS( folder, pattern="*.csv", keyName = keyName, parser_kwargs=None, )
        charaData1 = dop.copySamples(charaData,sampleMap)
        charaData1 = dop.renameSamples(charaData,sampleMap)
        # charaData2= dop.swapDictLevels(charaData1)
    if saveFile:
        saveXPS( charaData1, filename=out_path )
    print(f"Read {len(charaData)} samples.")
    if spectrum == "XPS":
        # charaData1 = cropXPSData(charaData1, xmin=700, xmax=3500)
        # charaData2 = baselineCorrectionXPS( charaData1, x=keyName[0], y=keyName[1], method="airpls", )
        # charaData2 = smoothXPSSpectrum( charaData1, y=keyName[1], regions=("Survey",), method="savgol", window_length=11, polyorder=3, )
        # 2. Extract regions
        charaDataRegions = extractXPSRegions( charaData1, highResolutionFlag=False, 
                                regions=("C1s", "O1s"), x="Binding Energy", )  
        if saveFile:
                saveXPS( charaDataRegions, filename=out_path )
        # charaDataRegions = smoothXPSSpectrum( charaDataRegions, y=keyName[1], regions=("C1s",), method="savgol", window_length=11, polyorder=3, )
        # Batch fitting
        
        
        
        
        
        
        
        return
        fitResults = fitXPSData( data=charaDataRegions, peakConfig=XPS_PEAK_CONFIG, regions=("O1s",), modelSelection=True, verbose=False, )
        # Collect all sample/model results
        fitSummary = collectXPSModelSummaries(fitResults)
        # Plot selected models
        for sample, sampleResults in fitResults.items():
            regionResult = sampleResults.get("C1s", {})
            modelResults = regionResult.get("results", {})
            for modelName in ("lorentzian_3peaks", "lorentzian_4peaks"):
                if modelName in modelResults:
                    plotXPSFit(
                        modelResults[modelName],
                        title=f"{sample} C1s - {modelName}",
                        # savepath=f"{sample}_C1s_{modelName}.png",
                    )
 ##以上为调试
        charaData3 = normalizeSpectrumXPS( charaData2, y=keyName[1], method="vector", )
        figureName = fl.get_expanded_name(out_path, fileName = "FTIR", type="png")
        figureName = fl.get_expanded_name( out_path, fileName="XPS", type="png", )
        plotSpectrumXPS(
            charaDataRegions,
            # regions=["Su"],
            figsize=(8, 4),
            x=keyName[0],
            y=keyName[1],
            offset=100,
            show_yticks=False,
            colors="PAPER1",
            linewidth=1.0,
            # savepath=figureName,
        )
        plotSpectrumXPS(
                    charaData1,
                    # regions=["Su"],
                    figsize=(8, 4),
                    x=keyName[0],
                    y=keyName[1],
                    offset=3000,
                    show_yticks=False,
                    colors="PAPER1",
                    linewidth=1.0,
                    # savepath=figureName,
                )
    y =0
if __name__ == "__main__":
    f = sys.argv[1] if len(sys.argv) > 1 else None
    main(f)