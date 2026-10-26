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
    },
}
def shirleyBackground(x, y, max_iter=200, tol=1e-5):
    """
    Iterative Shirley background estimation.
    Parameters
    ----------
    x : array-like
        Binding energy, sorted in ascending order.
    y : array-like
        Measured intensity.
    max_iter : int, default=200
        Maximum number of iterations.
    tol : float, default=1e-5
        Absolute convergence tolerance.
    Returns
    -------
    background : ndarray
        Estimated Shirley background, not background-corrected intensity.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.ndim != 1 or y.ndim != 1 or len(x) != len(y):
        raise ValueError("x and y must be 1D arrays with the same length.")
    if len(x) < 3:
        raise ValueError("At least 3 data points are required.")
    if not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
        raise ValueError("x and y must contain only finite values.")
    dx = np.diff(x)
    if np.any(dx <= 0):
        raise ValueError(
            "x must be strictly increasing. "
            "Sort by binding energy and remove duplicate x values first."
        )
    # Initial linear background between the two endpoints
    background = np.linspace(y[0], y[-1], len(y))
    for _ in range(max_iter):
        # Background-corrected signal used for the Shirley integral
        signal = np.maximum(y - background, 0.0)
        # Trapezoidal integral over each adjacent interval
        interval_integral = (
            (signal[:-1] + signal[1:]) * dx / 2.0
        )
        # Integral from each point to the high-energy endpoint
        reverse_integral = np.zeros_like(signal)
        reverse_integral[:-1] = np.cumsum(
            interval_integral[::-1]
        )[::-1]
        total = reverse_integral[0]
        if not np.isfinite(total) or total <= 0:
            return np.linspace(y[0], y[-1], len(y))
        # Update Shirley background
        new_background = (
            y[-1]
            + (y[0] - y[-1]) * reverse_integral / total
        )
        # Check convergence
        if np.max(np.abs(new_background - background)) < tol:
            background = new_background
            break
        background = new_background
    return background
def pseudoVoigt(x, area, center, fwhm, eta):
    """
    Area-normalized pseudo-Voigt peak.
    Parameters
    ----------
    x : array-like
        Binding energy.
    area : float
        Total integrated peak area over the infinite domain.
    center : float
        Peak center.
    fwhm : float
        Full width at half maximum.
    eta : float
        Lorentzian fraction, from 0 to 1.
    """
    x = np.asarray(x, dtype=float)
    if fwhm <= 0:
        raise ValueError("FWHM must be positive.")
    if not 0.0 <= eta <= 1.0:
        raise ValueError("Eta must be between 0 and 1.")
    # Gaussian component
    sigma = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    gaussian = ( np.exp(-0.5 * ((x - center) / sigma) ** 2) / (sigma * np.sqrt(2.0 * np.pi)) )
    # Lorentzian component
    gamma = fwhm / 2.0
    lorentzian = ( gamma / ( np.pi * ((x - center) ** 2 + gamma ** 2) ) )
    return area * ( (1.0 - eta) * gaussian + eta * lorentzian )
def multiPeakModel(x, *params):
    """
    params per peak: area, center, fwhm, eta
    """
    y = np.zeros_like(np.asarray(x, dtype=float))
    for i in range(0, len(params), 4):
        y += pseudoVoigt(x, *params[i:i + 4])
    return y
def fitXPSRegion(
    df,
    config,
    x="Binding Energy",
    y="Intensity",
    mode=None,
    verbose=True,
):
    """
    Fit one XPS region using one of two modes.
    Parameters
    ----------
    df : pandas.DataFrame
        Experimental XPS spectrum.
    config : dict
        Region-specific configuration.
    x : str
        Binding-energy column name.
    y : str
        Intensity column name.
    mode : {"envelope", "components"}, optional
        "envelope":
            Diagnostic fit of the dominant peak envelope using an
            asymmetric pseudo-Voigt profile. Does not provide
            chemical-state area percentages.
        "components":
            Simultaneous fitting of configured peaks using fixed Eta.
    verbose : bool
        Print diagnostics and fitted parameters.
    Returns
    -------
    dict
        x, raw, background, corrected, components,
        fitted_corrected, fitted, residual, peaks,
        covariance, metrics, success, fit_mode
    """
    import numpy as np
    import pandas as pd
    from scipy.optimize import curve_fit
    from scipy.integrate import trapezoid
    # ======================================
    # 0. Select fitting mode
    # ======================================
    fit_mode = (
        mode
        if mode is not None
        else config.get("fit_mode", "components")
    )
    fit_mode = fit_mode.lower().strip()
    if fit_mode not in {"envelope", "components"}:
        raise ValueError(
            "mode must be 'envelope' or 'components'."
        )
    # ======================================
    # 1. Prepare experimental data
    # ======================================
    work = df[[x, y]].copy()
    work[x] = pd.to_numeric(work[x], errors="coerce")
    work[y] = pd.to_numeric(work[y], errors="coerce")
    work = work.replace(
        [np.inf, -np.inf],
        np.nan,
    )
    work = work.dropna().sort_values(x)
    xmin, xmax = config.get(
        "x_range",
        (None, None),
    )
    if xmin is not None:
        work = work[work[x] >= xmin]
    if xmax is not None:
        work = work[work[x] <= xmax]
    # Average duplicate binding-energy points
    work = work.groupby(
        x,
        as_index=False,
    )[y].mean()
    X = work[x].to_numpy(dtype=float)
    Y = work[y].to_numpy(dtype=float)
    if len(X) < 10:
        raise ValueError(
            f"Insufficient data points: {len(X)}"
        )
    if not np.all(np.isfinite(X)) or not np.all(np.isfinite(Y)):
        raise ValueError(
            "Spectrum contains non-finite values."
        )
    if np.any(np.diff(X) <= 0):
        raise ValueError(
            "Binding energy must be strictly increasing."
        )
    idx_max = int(np.argmax(Y))
    idx_min = int(np.argmin(Y))
    if verbose:
        print(
            f"[Raw spectrum] "
            f"Max: {X[idx_max]:.3f} eV, {Y[idx_max]:.2f}; "
            f"Min: {X[idx_min]:.3f} eV, {Y[idx_min]:.2f}"
        )
    # ======================================
    # 2. Estimate background
    # ======================================
    background_method = config.get(
        "background",
        "shirley",
    ).lower()
    if background_method == "shirley":
        background = shirleyBackground(X, Y)
    elif background_method == "linear":
        background = np.linspace(
            Y[0],
            Y[-1],
            len(Y),
        )
    elif background_method == "none":
        background = np.zeros_like(Y)
    else:
        raise ValueError(
            f"Unsupported background: {background_method}"
        )
    corrected = Y - background
    if verbose:
        idx_corr_max = int(np.argmax(corrected))
        print(
            f"[Peak diagnostic] "
            f"Raw max = {X[idx_max]:.3f} eV; "
            f"Corrected max = {X[idx_corr_max]:.3f} eV"
        )
    # ======================================
    # 3. Fit-quality helper
    # ======================================
    def calculateMetrics(raw, fitted, residual):
        ss_res = np.sum(residual ** 2)
        ss_tot = np.sum(
            (raw - np.mean(raw)) ** 2
        )
        r2 = (
            1.0 - ss_res / ss_tot
            if ss_tot > 0
            else np.nan
        )
        rmse = np.sqrt(
            np.mean(residual ** 2)
        )
        return {
            "R2": float(r2),
            "RMSE": float(rmse),
            "n_points": len(raw),
        }
    def printResidualDiagnostics(
        X_fit,
        Y_fit,
        fitted,
        residual,
    ):
        if not verbose:
            return
        idx_overfit = int(np.argmin(residual))
        idx_underfit = int(np.argmax(residual))
        print(
            f"[Largest overfit] "
            f"BE={X_fit[idx_overfit]:.3f} eV, "
            f"Raw={Y_fit[idx_overfit]:.2f}, "
            f"Fit={fitted[idx_overfit]:.2f}, "
            f"Residual={residual[idx_overfit]:.2f}"
        )
        print(
            f"[Largest underfit] "
            f"BE={X_fit[idx_underfit]:.3f} eV, "
            f"Raw={Y_fit[idx_underfit]:.2f}, "
            f"Fit={fitted[idx_underfit]:.2f}, "
            f"Residual={residual[idx_underfit]:.2f}"
        )
    def checkParameterBounds(
        params,
        lower,
        upper,
        param_names,
    ):
        bound_hits = []
        for i, name in enumerate(param_names):
            value = params[i]
            lo = lower[i]
            hi = upper[i]
            # Only positive/finite parameter values are expected here
            atol = 1e-3 * max(1.0, abs(value))
            if np.isclose(
                value,
                lo,
                rtol=0,
                atol=atol,
            ):
                bound_hits.append(
                    f"{name} at LOWER bound "
                    f"({value:.6g})"
                )
            if (
                np.isfinite(hi)
                and np.isclose(
                    value,
                    hi,
                    rtol=0,
                    atol=atol,
                )
            ):
                bound_hits.append(
                    f"{name} at UPPER bound "
                    f"({value:.6g})"
                )
        if verbose:
            if bound_hits:
                print(
                    "[Fit bounds] Parameters near bounds:"
                )
                for item in bound_hits:
                    print(f"  - {item}")
            else:
                print(
                    "[Fit bounds] "
                    "No parameters near bounds."
                )
        return bound_hits
    # ======================================
    # 4. Envelope diagnostic mode
    # ======================================
    if fit_mode == "envelope":
        # Explicit range is required to avoid accidentally
        # fitting the entire C1s region with one broad peak.
        envelope_range = config.get(
            "envelope_range",
            None,
        )
        if envelope_range is None:
            raise ValueError(
                "Envelope mode requires "
                "config['envelope_range'], e.g. (280, 286) "
                "for diagnosing the dominant C1s peak."
            )
        env_min, env_max = envelope_range
        env_mask = (
            (X >= env_min)
            & (X <= env_max)
        )
        if np.count_nonzero(env_mask) < 10:
            raise ValueError(
                "Insufficient data points in envelope_range."
            )
        X_fit = X[env_mask]
        Y_fit = Y[env_mask]
        B_fit = background[env_mask]
        C_fit = corrected[env_mask]
        # --------------------------------------
        # 4.1 Initial peak position and bounds
        # --------------------------------------
        idx_peak = int(np.argmax(C_fit))
        center0 = float(X_fit[idx_peak])
        # If no explicit center bounds are supplied,
        # use the observed corrected-signal maximum
        # with a local search interval around it.
        if "envelope_center_bounds" in config:
            center_lo, center_hi = config[
                "envelope_center_bounds"
            ]
        else:
            center_lo = max(
                float(X_fit[0]),
                center0 - 0.8,
            )
            center_hi = min(
                float(X_fit[-1]),
                center0 + 0.8,
            )
        if center_lo >= center_hi:
            raise ValueError(
                "Invalid envelope_center_bounds."
            )
        if not center_lo <= center0 <= center_hi:
            center0 = float(
                np.clip(
                    center0,
                    center_lo,
                    center_hi,
                )
            )
        # Width bounds apply separately to the
        # lower- and higher-binding-energy sides.
        width_lo, width_hi = config.get(
            "envelope_width_bounds",
            (0.5, 4.0),
        )
        if not (
            0 < width_lo < width_hi
        ):
            raise ValueError(
                "envelope_width_bounds must satisfy "
                "0 < lower < upper."
            )
        # Fixed line-shape mixing coefficient
        eta = float( config.get( "envelope_eta_fixed", 0.3, ) )
        if not 0.0 <= eta <= 1.0:
            raise ValueError( "envelope_eta_fixed must be between 0 and 1." )
        height0 = max( float(np.max(C_fit)), 1e-8, )
        width0 = float( np.clip( 1.5, width_lo + 1e-6, width_hi - 1e-6, ) )
        # --------------------------------------
        # 4.2 Asymmetric pseudo-Voigt envelope
        # --------------------------------------
        def asymmetricEnvelope( xx, height, center, width_left, width_right, ):
            """
            Diagnostic envelope model.
            The left and right sides have separate width
            parameters. This is a flexible line-shape
            diagnostic, not a chemical-state peak model.
            width_left:
                Width parameter on the lower-BE side.
            width_right:
                Width parameter on the higher-BE side.
            """
            xx = np.asarray( xx, dtype=float, )
            widths = np.where( xx < center, width_left, width_right, )
            z = (xx - center) / widths
            gaussian = np.exp( -4.0 * np.log(2.0) * z ** 2 )
            lorentzian = 1.0 / ( 1.0 + 4.0 * z ** 2 )
            return height * ( (1.0 - eta) * gaussian + eta * lorentzian )
        p0 = [ height0, center0, width0, width0, ]
        lower = [ 0.0, center_lo, width_lo, width_lo, ]
        upper = [ np.inf, center_hi, width_hi, width_hi, ]
        # --------------------------------------
        # 4.3 Fit envelope to corrected experimental data
        # --------------------------------------
        popt, pcov = curve_fit( asymmetricEnvelope, X_fit, C_fit, p0=p0, bounds=(lower, upper), maxfev=50000, )
        fitted_corrected = asymmetricEnvelope( X_fit, *popt, )
        fitted_total = B_fit + fitted_corrected
        residual = Y_fit - fitted_total
        height, center, width_left, width_right = popt
        # Diagnostic envelope integral only.
        # This is NOT a chemical-state area.
        envelope_area = trapezoid( fitted_corrected, X_fit, )
        peak_table = pd.DataFrame([{
            "Feature": "Main envelope (diagnostic only)",
            "Center (eV)": center,
            "Width left (eV)": width_left,
            "Width right (eV)": width_right,
            "Eta fixed": eta,
            "Peak height": height,
            "Envelope area in fit window": envelope_area,
        }])
        param_names = [ "Envelope height", "Envelope center", "Width left", "Width right", ]
        bound_hits = checkParameterBounds( popt, lower, upper, param_names, )
        metrics = calculateMetrics( Y_fit, fitted_total, residual, )
        printResidualDiagnostics( X_fit, Y_fit, fitted_total, residual, )
        if verbose:
            print(
                f"[Envelope fit] "
                f"R2={metrics['R2']:.4f}; "
                f"RMSE={metrics['RMSE']:.2f}"
            )
            print(
                f"[Envelope parameters] "
                f"Center={center:.3f} eV; "
                f"Width left={width_left:.3f} eV; "
                f"Width right={width_right:.3f} eV"
            )
        return {
            "x": X_fit,
            "raw": Y_fit,
            "background": B_fit,
            "corrected": C_fit,
            "components": {
                "Main envelope (diagnostic only)": (
                    fitted_corrected
                ),
            },
            "fitted_corrected": fitted_corrected,
            "fitted": fitted_total,
            "residual": residual,
            "peaks": peak_table,
            "covariance": pcov,
            "metrics": metrics,
            "success": True,
            "fit_mode": "envelope",
            "diagnostic_only": True,
            "fit_range": (
                float(X_fit[0]),
                float(X_fit[-1]),
            ),
            "bound_hits": bound_hits,
            "full_spectrum": {
                "x": X,
                "raw": Y,
                "background": background,
                "corrected": corrected,
            },
        }
    # ======================================
    # 5. Multi-component fitting mode
    # ======================================
    peak_configs = config.get( "peaks", [], )
    if not peak_configs:
        raise ValueError( "Component mode requires config['peaks']." )
    signal_scale = max( float(np.max(corrected)), float(np.std(corrected)), 1e-8, )
    p0 = []
    lower = []
    upper = []
    fixed_etas = []
    for peak in peak_configs:
        center = float(peak["center"])
        center_lo, center_hi = peak["center_bounds"]
        fwhm = float( peak.get("fwhm", 1.2) )
        fwhm_lo, fwhm_hi = peak.get( "fwhm_bounds", (0.5, 2.0), )
        eta = float( peak.get( "eta_fixed", peak.get("eta", 0.3), ) )
        if not 0.0 <= eta <= 1.0:
            raise ValueError( f"Eta for '{peak['name']}' must be between 0 and 1." )
        if not center_lo <= center <= center_hi:
            raise ValueError( f"Initial center for '{peak['name']}' is outside center_bounds." )
        if not fwhm_lo <= fwhm <= fwhm_hi:
            raise ValueError( f"Initial FWHM for '{peak['name']}' is outside fwhm_bounds." )
        fixed_etas.append(eta)
        area0 = max( signal_scale * fwhm, 1e-8, )
        p0.extend([ area0, center, fwhm, ])
        lower.extend([ 0.0, center_lo, fwhm_lo, ])
        upper.extend([ np.inf, center_hi, fwhm_hi, ])
    p0 = np.asarray(p0, dtype=float)
    lower = np.asarray(lower, dtype=float)
    upper = np.asarray(upper, dtype=float)
    # --------------------------------------
    # 5.1 Build fixed-Eta multi-peak model
    # --------------------------------------
    def model(xx, *params):
        xx = np.asarray( xx, dtype=float, )
        y_model = np.zeros_like(xx)
        for i, peak in enumerate(peak_configs):
            area, center, fwhm = params[ 3 * i: 3 * i + 3 ]
            eta = fixed_etas[i]
            y_model += pseudoVoigt( xx, area, center, fwhm, eta, )
        return y_model
    # --------------------------------------
    # 5.2 Simultaneous fitting
    # --------------------------------------
    popt, pcov = curve_fit( model, X, corrected, p0=p0, bounds=(lower, upper), maxfev=50000, )
    fitted_corrected = model( X, *popt, )
    fitted_total = background + fitted_corrected
    residual = Y - fitted_total
    # --------------------------------------
    # 5.3 Check parameter bounds
    # --------------------------------------
    param_names = []
    for peak in peak_configs:
        param_names.extend([ f"{peak['name']}.Area", f"{peak['name']}.Center", f"{peak['name']}.FWHM", ])
    bound_hits = checkParameterBounds( popt, lower, upper, param_names, )
    # --------------------------------------
    # 5.4 Peak areas
    # --------------------------------------
    components = {}
    peak_rows = []
    for i, peak in enumerate(peak_configs):
        area, center, fwhm = popt[ 3 * i: 3 * i + 3 ]
        eta = fixed_etas[i]
        component = pseudoVoigt( X, area, center, fwhm, eta, )
        components[peak["name"]] = component
        measured_area = trapezoid( component, X, )
        peak_rows.append({
            "Peak": peak["name"],
            "Center (eV)": center,
            "FWHM (eV)": fwhm,
            "Eta": eta,
            "Area parameter": area,
            "Area in region": measured_area,
        })
    peak_table = pd.DataFrame( peak_rows )
    total_area = peak_table[ "Area in region" ].sum()
    if total_area > 0:
        peak_table["Area (%)"] = ( peak_table["Area in region"] / total_area * 100.0 )
    else:
        peak_table["Area (%)"] = np.nan
    # --------------------------------------
    # 5.5 Fit quality and diagnostics
    # --------------------------------------
    metrics = calculateMetrics( Y, fitted_total, residual, )
    printResidualDiagnostics( X, Y, fitted_total, residual, )
    if verbose:
        print( f"[Component fit] R2={metrics['R2']:.4f}; RMSE={metrics['RMSE']:.2f}" )
        print("[Fitted peak centers]")
        for i, peak in enumerate(peak_configs):
            center = popt[3 * i + 1]
            print( f"  {peak['name']}: " f"{center:.3f} eV" )
    # ======================================
    # 6. Return results
    # ======================================
    return {
        "x": X,
        "raw": Y,
        "background": background,
        "corrected": corrected,
        "components": components,
        "fitted_corrected": fitted_corrected,
        "fitted": fitted_total,
        "residual": residual,
        "peaks": peak_table,
        "covariance": pcov,
        "metrics": metrics,
        "success": True,
        "fit_mode": "components",
        "diagnostic_only": False,
        "fit_range": (
            float(X[0]),
            float(X[-1]),
        ),
        "bound_hits": bound_hits,
    }
def fitXPSData( data, regions=("C1s",), peakConfig=None, x="Binding Energy", y="Intensity", ):
    """
    Fit selected XPS regions for all samples.
    Returns
    -------
    dict
        results[sample][region] = fit result
    """
    if peakConfig is None:
        peakConfig = XPS_PEAK_CONFIG
    results = {}
    for sample, sample_data in data.items():
        results[sample] = {}
        for region in regions:
            df = sample_data.get(region)
            if not isinstance(df, pd.DataFrame):
                print(f"[Skip] {sample} / {region}: no DataFrame")
                continue
            if region not in peakConfig:
                print(f"[Skip] {sample} / {region}: no peak config")
                continue
            try:
                results[sample][region] = fitXPSRegion( df=df, config=peakConfig[region], x=x, y=y, )
                print(
                    f"[Fit] {sample} / {region}: "
                    f"R²={results[sample][region]['metrics']['R2']:.4f}" )
            except (ValueError, RuntimeError, FloatingPointError) as exc:
                results[sample][region] = { "success": False, "error": str(exc), }
                print(f"[Failed] {sample} / {region}: {exc}")
    return results
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
        regions.update(
            key
            for key, value in sampleData.items()
            if isinstance(value, pd.DataFrame)
        )
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
            temp = df[
                [
                    "Binding Energy",
                    "Intensity",
                ]
            ].reset_index(drop=True)
            data_blocks.append(temp)
            header_row1.extend(
                [
                    sample,
                    sample,
                ]
            )
            header_row2.extend(
                [
                    "Binding Energy",
                    "Intensity",
                ]
            )
        if len(data_blocks) == 0:
            continue
        # ---------------------
        # Horizontal concatenate
        # ---------------------
        output = pd.concat(
            data_blocks,
            axis=1,
        )
        # ---------------------
        # Two header rows
        # ---------------------
        header = pd.DataFrame(
            [
                header_row1,
                header_row2,
            ],
            columns=output.columns,
        )
        final = pd.concat(
            [
                header,
                output,
            ],
            ignore_index=True,
        )
        # ---------------------
        # Export
        # ---------------------
        fl.export_to_excel_auto(
            final,
            filename=filename,
            sheet_name=region,
        )
        print(f"Export {region} -> {filename}")
def readXPS(
    folder,
    pattern="*.csv",
    keyName=None,
    parser_kwargs=None,
):
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
            xps[sample_name] = {
                "metadata": metadata,
            }
        for _, regions in grouped_data.items():
            for region_name, df in regions.items():
                region_key = str(region_name).strip().lower()
                standardized_region = region_map.get(
                    region_key,
                    str(region_name).strip(),
                )
                xps[sample_name][standardized_region] = df
    return xps
def plotXPSFit(
    fit_result,
    title=None,
    figsize=(7, 5),
    savepath=None,
    dpi=600,
):
    """
    Plot raw spectrum, background, individual peaks,
    total fitted spectrum, and residual.
    """
    if not fit_result.get("success", False):
        raise ValueError("Cannot plot an unsuccessful fit.")
    X = fit_result["x"]
    Y = fit_result["raw"]
    fig, (ax, ax_res) = plt.subplots(
        2, 1,
        figsize=figsize,
        sharex=True,
        gridspec_kw={"height_ratios": [3, 1]},
    )
    ax.scatter(
        X, Y, s=10, label="Data", zorder=3
    )
    ax.plot(
        X,
        fit_result["background"],
        linestyle="--",
        label="Background",
    )
    for name, component in fit_result["components"].items():
        ax.plot(
            X,
            fit_result["background"] + component,
            linewidth=1.2,
            label=name,
        )
    ax.plot(
        X,
        fit_result["fitted"],
        linewidth=2,
        label="Total fit",
    )
    ax.set_ylabel("Intensity (a.u.)")
    ax.legend(fontsize=8)
    ax_res.axhline(0, linestyle="--", linewidth=1)
    ax_res.plot(X, fit_result["residual"], linewidth=1)
    ax_res.set_xlabel("Binding Energy (eV)")
    ax_res.set_ylabel("Residual")
    # XPS convention: high binding energy on the left
    ax_res.set_xlim(X.max(), X.min())
    if title is not None:
        ax.set_title(title)
    fig.tight_layout()
    if savepath is not None:
        fig.savefig(
            savepath,
            dpi=dpi,
            bbox_inches="tight",
        )
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
                result[sample][region] = (
                    survey.loc[mask].copy().reset_index(drop=True)
                )
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
def baselineCorrectionXPS(
    data,
    x="Binding Energy",
    y="Intensity",
    method="asls",
):
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
def normalizeSpectrumXPS(
    data,
    y="Intensity",
    method="vector",
):
    """Normalize nested XPS spectra."""
    result = {}
    for sample, sample_data in data.items():
        result[sample] = {}
        for region, value in sample_data.items():
            if isinstance(value, pd.DataFrame):
                normalized = dop.normalizeSpectrum(
                    {sample: value},
                    y=y,
                    method=method,
                )
                result[sample][region] = normalized[sample]
            else:
                result[sample][region] = value
    return result
def main(file_path=None):
    spectrum = "XPS"
    saveFile = False
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
        charaData2= dop.swapDictLevels(charaData1)
    if saveFile:
        saveXPS( charaData1, filename=out_path )
    print(f"Read {len(charaData)} samples.")
    if spectrum == "XPS":
        # charaData1 = cropXPSData(charaData1, xmin=700, xmax=3500)
        # charaData2 = baselineCorrectionXPS( charaData1, x=keyName[0], y=keyName[1], method="airpls", )
        # 2. Extract regions
        charaDataRegions = extractXPSRegions( charaData1, highResolutionFlag=False, 
                                regions=("C1s", "O1s"), x="Binding Energy", )      
       #调试部分
        # 从原有配置复制 C1s 参数，不修改原始配置
        # =========================================================
        # 1. Configure envelope fitting
        # =========================================================
        envelopeConfig = copy.deepcopy(
            XPS_PEAK_CONFIG["C1s"]
        )
        # Main-peak diagnostic range (will be varied below)
        base_range = (280.0, 286.0)
        # Width constraints for the asymmetric envelope
        envelopeConfig["envelope_width_bounds"] = (
            0.5,
            4.0,
        )
        # Fixed pseudo-Voigt mixing parameter
        envelopeConfig["envelope_eta_fixed"] = 0.3
        # Test different upper limits of the envelope fitting range
        test_ranges = [
            (280.0, 286.0),
            (280.0, 286.5),
            (280.0, 287.0),
        ]
        # Samples to analyze
        samples = ["S1", "S3"]
        # =========================================================
        # 2. Fit all samples under all tested ranges
        # =========================================================
        # Store complete fitting results:
        # envelopeResults[sample][energy_range] = result
        envelopeResults = {}
        # Collect summary statistics for the comparison table
        summary_rows = []
        for sample in samples:
            print("\n" + "=" * 60)
            print(f"Sample: {sample} / C1s envelope")
            # Retrieve C1s data
            df = charaDataRegions[sample]["C1s"]
            envelopeResults[sample] = {}
            for energy_range in test_ranges:
                print(
                    f"\n--- Fitting range: "
                    f"{energy_range[0]:.1f}–"
                    f"{energy_range[1]:.1f} eV ---"
                )
                # Copy the base configuration so that each test
                # is independent of the previous fitting run
                testConfig = copy.deepcopy(envelopeConfig)
                testConfig["envelope_range"] = energy_range
                # Fit the envelope
                result = fitXPSRegion(
                    df=df,
                    config=testConfig,
                    mode="envelope",
                    verbose=False,
                )
                # Store the complete result
                envelopeResults[sample][energy_range] = result
                # Extract fitted envelope parameters
                peak = result["peaks"].iloc[0]
                metrics = result["metrics"]
                # Store diagnostics
                bound_hits = result.get(
                    "bound_hits",
                    [],
                )
                summary_rows.append({
                    "Sample": sample,
                    "Range (eV)": (
                        f"{energy_range[0]:.1f}–"
                        f"{energy_range[1]:.1f}"
                    ),
                    "R2": metrics["R2"],
                    "RMSE": metrics["RMSE"],
                    "Center (eV)": peak["Center (eV)"],
                    "Width left (eV)": peak["Width left (eV)"],
                    "Width right (eV)": peak["Width right (eV)"],
                    "Boundary hits": (
                        "; ".join(bound_hits)
                        if bound_hits
                        else "None"
                    ),
                })
                print(
                    f"R2 = {metrics['R2']:.4f}; "
                    f"RMSE = {metrics['RMSE']:.2f}"
                )
                print(
                    f"Center = {peak['Center (eV)']:.3f} eV; "
                    f"Width left = {peak['Width left (eV)']:.3f} eV; "
                    f"Width right = {peak['Width right (eV)']:.3f} eV"
                )
                print(
                    "Boundary hits:",
                    bound_hits if bound_hits else "None",
                )
        # =========================================================
        # 3. Display the range-sensitivity comparison table
        # =========================================================
        sensitivityTable = pd.DataFrame(summary_rows)
        print("\n" + "=" * 80)
        print("ENVELOPE FITTING RANGE SENSITIVITY")
        print("=" * 80)
        print(
            sensitivityTable.to_string(
                index=False,
                float_format=lambda value: f"{value:.4f}",
            )
        )
        # =========================================================
        # 4. Plot the baseline-range fit for S1 and S3
        # =========================================================
        for sample in samples:
            result = envelopeResults[sample][base_range]
            print("\n" + "=" * 60)
            print(
                f"Plotting {sample} / C1s envelope "
                f"({base_range[0]:.1f}–{base_range[1]:.1f} eV)"
            )
            plotXPSFit(
                result,
                title=(
                    f"{sample} - C1s envelope "
                    f"({base_range[0]:.1f}–{base_range[1]:.1f} eV)"
                ),
                # Uncomment the next line to save each figure
                # savepath=f"{sample}_C1s_envelope.png",
            )

        # Keep the original fitting ranges
        test_ranges = [
            (280.0, 286.0),
            (280.0, 286.5),
            (280.0, 287.0),
        ]

        # Test different upper limits for the right-side width
        right_width_limits = [4.0, 5.0, 6.0]

        samples = ["S1", "S3"]
        widthTestResults = []

        for sample in samples:

            df = charaDataRegions[sample]["C1s"]

            for energy_range in test_ranges:

                for width_upper in right_width_limits:

                    test_config = copy.deepcopy(envelopeConfig)

                    test_config["envelope_range"] = energy_range
                    test_config["envelope_width_bounds"] = (
                        0.5,
                        width_upper,
                    )

                    result = fitXPSRegion(
                        df=df,
                        config=test_config,
                        mode="envelope",
                        verbose=False,
                    )

                    peak = result["peaks"].iloc[0]
                    metrics = result["metrics"]

                    width_right = peak["Width right (eV)"]

                    width_at_upper = abs(
                        width_right - width_upper
                    ) < 1e-3

                    widthTestResults.append({
                        "Sample": sample,
                        "Range": (
                            f"{energy_range[0]:.1f}-"
                            f"{energy_range[1]:.1f}"
                        ),
                        "Width upper": width_upper,
                        "Width right": width_right,
                        "Center": peak["Center (eV)"],
                        "R2": metrics["R2"],
                        "RMSE": metrics["RMSE"],
                        "At upper bound": width_at_upper,
                    })

        widthTestTable = pd.DataFrame(widthTestResults)

        print(
            widthTestTable.to_string(
                index=False,
                float_format=lambda value: f"{value:.4f}",
            )
        )
# =========================================================
# 5. Optional: export the comparison table
# =========================================================
# Uncomment if you want to save the summary to Excel
# sensitivityTable.to_excel(
#     "C1s_envelope_range_sensitivity.xlsx",
#     index=False,
# )
        ##调试
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
        testPeakConfig = copy.deepcopy(XPS_PEAK_CONFIG)
        for peak in testPeakConfig["C1s"]["peaks"]:
            if peak["name"] == "C-C/C=C":
                peak["fwhm_bounds"] = (0.5, 3.5)
        fitResults_wide = fitXPSData(
            charaDataRegions,
            regions=("C1s",),
            peakConfig=testPeakConfig,
        )
        # fitResults = fitXPSData( charaDataRegions, regions=("C1s",), peakConfig=XPS_PEAK_CONFIG, )
        # 先关闭之前的图窗，避免把旧图误认为新结果
        for sample in ["S1", "S3"]:
            result = fitResults_wide[sample]["C1s"]
            mask = (
                (result["x"] >= 281.8)
                & (result["x"] <= 283.6)
            )
            diagnostic = pd.DataFrame({
                "BE (eV)": result["x"][mask],
                "Raw": result["raw"][mask],
                "Background": result["background"][mask],
                "Corrected": result["corrected"][mask],
                "Fitted components": result["fitted_corrected"][mask],
                "Total fit": result["fitted"][mask],
                "Residual": result["residual"][mask],
            })
            print(f"\n===== {sample}: 281.8–283.6 eV =====")
            print(diagnostic.round(2).to_string(index=False))
        plt.close("all")
        for sample in ["S1", "S3"]:
            result = fitResults_wide[sample]["C1s"]
            X = result["x"]
            Y = result["raw"]
            idx_max = np.argmax(Y)
            top_idx = np.argsort(Y)[-5:][::-1]
            print(f"\n{'=' * 40}")
            print(f"Sample: {sample}")
            print(f"Maximum: BE = {X[idx_max]:.3f} eV, "
                f"Intensity = {Y[idx_max]:.3f}")
            print(f"Intensity range: {Y.min():.3f} – {Y.max():.3f}")
            print("\nTop 5 intensity points:")
            print(pd.DataFrame({
                "Binding Energy (eV)": X[top_idx],
                "Intensity": Y[top_idx],
            }))
            plotXPSFit(
                result,
                title=f"{sample} - C1s (current fit)",
            )
        plt.show()
        # samples = ["S1", "S3"]
        # region = "C1s"
        # for sample in samples:
        #     result = fitResults[sample][region]
        #     if result["success"]:
        #         print(f"\n{'=' * 40}")
        #         print(f"Sample: {sample} / {region}")
        #         print(result["peaks"])
        #         print(result["metrics"])
        #         plotXPSFit(
        #             result,
        #             title=f"{sample} - {region}",
        #             # savepath=f"{sample}_{region}_fit.png",
        #         )
        #     else:
        #         print(f"[Warning] Fit failed: {sample} / {region}")
    # elif spectrum == "Raman": 
    #     charaData1 = dop.cropX(charaData1, xmin=500, xmax=2200)
    #     charaData2 = dop.smoothSpectrum( charaData1, y=keyName[1], method="gaussian", sigma=1.2, )
    #     charaData3 = dop.baselineCorrection( charaData2, x=keyName[0], y=keyName[1], method="rubberband", )
    #     figureName = fl.get_expanded_name(out_path, fileName = "Raman", type="png")
    #     myPlt.plotSpectrum( charaData1, figsize=(8, 4), x=keyName[0], y=keyName[1], 
    #                    reverse_x=True, offset=0, colors="PAPER1", linewidth=0.5, )
    #     myPlt.plotSpectrum( charaData2, figsize=(8, 4), x=keyName[0], y=keyName[1], 
    #                    offset=80, colors="PAPER1", linewidth=0.5, )
    #     myPlt.plotSpectrum( charaData3, figsize=(8, 4), x=keyName[0], y=keyName[1], 
    #                    offset=50, colors="PAPER1", linewidth=0.5, savepath=figureName)
    # plt.show(block=True)
    y =0
if __name__ == "__main__":
    f = sys.argv[1] if len(sys.argv) > 1 else None
    main(f)