from common import *
import T1dataProcessSource.dataOperation as dop
def calculateFitMetrics(y_true, y_pred):
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)
    if y_true.shape != y_pred.shape:
        raise ValueError( "y_true and y_pred must have the same shape." )
    residual = y_true - y_pred
    ss_res = np.sum(residual ** 2)
    ss_tot = np.sum( (y_true - np.mean(y_true)) ** 2 )
    if ss_tot > 0:
        r2 = 1 - ss_res / ss_tot
    else:
        r2 = np.nan
    rmse = np.sqrt( np.mean(residual ** 2) )
    y_mean = np.mean(y_true)
    rmse_rel = ( rmse / (np.abs(y_mean) + 1e-12) )
    y_range = ( np.max(y_true) - np.min(y_true) )
    rmse_nrm = ( rmse / (y_range + 1e-12) )
    mae = np.mean( np.abs(residual) )
    with np.errstate( divide="ignore", invalid="ignore" ):
        relative_residual = np.where( y_true != 0, residual / y_true, 0.0 )
    mape = np.mean( np.abs(relative_residual) ) * 100
    return {
        "Residual": residual,
        "RelativeResidual": relative_residual,
        "SSE": ss_res,
        "R2": r2,
        "RMSE": rmse,
        "RMSE_rel": rmse_rel,
        "RMSE_nrm": rmse_nrm,
        "MAE": mae,
        "MAPE": mape,
    }
def calculateModelSelectionMetrics( sse, n, k):
    """
    Calculate model selection metrics.
    Parameters
    ----------
    sse : float
        Sum of squared errors.
    n : int
        Number of observations.
    k : int
        Number of fitted parameters.
    Returns
    -------
    dict
        {
            "AIC": ...,
            "AICc": ...
        }
    """
    # ==========================================================
    # Validate
    # ==========================================================
    if not np.isfinite(sse):
        return { "AIC": np.nan, "AICc": np.nan }
    # ==========================================================
    # AIC
    # ==========================================================
    if sse > 0:
        aic = ( n * np.log(sse / n) + 2 * k )
    elif sse == 0:
        aic = -np.inf
    else:
        aic = np.nan
    # ==========================================================
    # AICc
    # ==========================================================
    if ( np.isfinite(aic) and n - k - 1 > 0 ):
        aicc = ( aic + 2 * k * (k + 1) / (n - k - 1) )
    elif aic == -np.inf:
        aicc = -np.inf
    else:
        aicc = np.inf
    return {
        "AIC": aic,
        "AICc": aicc
    }
def compute_fit_metrics(fit_x, fit_y, exp_x, exp_y, method='linear', use_interp=True,):
    # from sklearn.metrics import r2_score 旧版，新版是computeFitMetrics（）
    """
    计算拟合曲线与实验数据的拟合指标
    包括：
    - 残差 residuals
    - 相对残差 rems
    - R² r2
    自动处理拟合数据和实验数据维度不一致，通过插值对齐。
    参数:
    fit_x, fit_y : 拟合曲线数据
    exp_x, exp_y : 实验数据
    method : 插值方法，默认线性插值
    返回:
    residuals, relative_residuals, r2, rmsd
    """
    # 拟合曲线插值到实验点
    if use_interp:
        interp_func = interp1d(fit_x, fit_y, kind=method, fill_value="extrapolate")
        fit_y_interp = interp_func(exp_x)
    else:   ##计算全局数据时，不能使用插值，必须在调用函数前保证数据点数和顺序的准确性
        fit_y_interp = fit_y
    # 残差
    residuals = exp_y - fit_y_interp
    # 相对残差
    with np.errstate(divide='ignore', invalid='ignore'):
        relative_residuals = np.where(exp_y != 0, residuals / exp_y, 0.0)
    # R²
    r2 = r2_score(exp_y, fit_y_interp)
    rmsd = np.sqrt(np.mean((residuals)**2))
    return residuals, relative_residuals, r2, rmsd
def computeFitMetrics( fit_x, fit_y, exp_x, exp_y, method="linear", use_interp=True, allow_extrapolation=False):
    """
    Compute fit metrics between a fitted/predicted curve
    and experimental data.
    The fitted curve is aligned to the experimental pressure
    points by interpolation.
    Parameters
    ----------
    fit_x, fit_y : array-like
        Pressure and uptake of the fitted/predicted curve.
    exp_x, exp_y : array-like
        Pressure and uptake of the experimental data.
    method : str, default="linear"
        Interpolation method.
    use_interp : bool, default=True
        Whether to interpolate the fitted curve onto the
        experimental x coordinates.
    allow_extrapolation : bool, default=False
        Whether to allow extrapolation outside the overlapping
        x range.
    Returns
    -------
    dict
        Metrics returned by calculateFitMetrics().
    """
    fit_x = np.asarray(fit_x, dtype=float).ravel()
    fit_y = np.asarray(fit_y, dtype=float).ravel()
    exp_x = np.asarray(exp_x, dtype=float).ravel()
    exp_y = np.asarray(exp_y, dtype=float).ravel()
    # ======================================================
    # Basic validation
    # ======================================================
    if fit_x.size != fit_y.size:
        raise ValueError( "fit_x and fit_y must have the same length." )
    if exp_x.size != exp_y.size:
        raise ValueError( "exp_x and exp_y must have the same length." )
    if fit_x.size == 0 or exp_x.size == 0:
        raise ValueError( "Fit and experimental data must not be empty." )
    # ======================================================
    # Sort by x
    # ======================================================
    fit_idx = np.argsort(fit_x)
    fit_x = fit_x[fit_idx]
    fit_y = fit_y[fit_idx]
    exp_idx = np.argsort(exp_x)
    exp_x = exp_x[exp_idx]
    exp_y = exp_y[exp_idx]
    # ======================================================
    # Align x coordinates
    # ======================================================
    if use_interp:
        if allow_extrapolation:
            interp_func = interp1d( fit_x, fit_y, kind=method, fill_value="extrapolate", bounds_error=False )
            exp_x_eval = exp_x
            exp_y_eval = exp_y
            fit_y_aligned = interp_func( exp_x_eval )
        else:
            # --------------------------------------------------
            # Use only overlapping x range
            # --------------------------------------------------
            xmin = max( fit_x.min(), exp_x.min() )
            xmax = min( fit_x.max(), exp_x.max() )
            mask = ( (exp_x >= xmin) & (exp_x <= xmax) )
            exp_x_eval = exp_x[mask]
            exp_y_eval = exp_y[mask]
            if exp_x_eval.size < 2:
                raise ValueError( "Less than two overlapping points between fit and experimental data." )
            interp_func = interp1d( fit_x, fit_y, kind=method, bounds_error=False )
            fit_y_aligned = interp_func( exp_x_eval )
    else:
        if fit_y.shape != exp_y.shape:
            raise ValueError( "fit_y and exp_y must have the same shape when use_interp=False." )
        exp_x_eval = exp_x
        exp_y_eval = exp_y
        fit_y_aligned = fit_y
    # ======================================================
    # Calculate metrics
    # ======================================================
    return calculateFitMetrics( y_true=exp_y_eval, y_pred=fit_y_aligned )
def compute_global_fit_metrics(exp_y, fit_y):
    #旧版，新版是computeGlobalFitMetrics（）
    """
    计算全局拟合评价指标（global metrics）
    Parameters
    ----------
    exp_y : array-like
        所有实验数据（可由多温度、多等温线拼接而成）
    fit_y : array-like
        与 exp_y 一一对应的拟合数据
    Returns
    -------
    metrics : dict
        {
            "residuals": ndarray,
            "relative_residuals": ndarray,
            "R2": float,
            "RMSD": float
        }
    """
    import numpy as np
    from sklearn.metrics import r2_score
    exp_y = np.asarray(exp_y)
    fit_y = np.asarray(fit_y)
    if exp_y.shape != fit_y.shape:
        raise ValueError("exp_y and fit_y must have the same shape.")
    # residuals
    residuals = exp_y - fit_y
    # relative residuals
    with np.errstate(divide='ignore', invalid='ignore'):
        relative_residuals = np.where(exp_y != 0, residuals / exp_y, 0.0)
    # R²
    r2 = r2_score(exp_y, fit_y)
    # RMSD
    rmsd = np.sqrt(np.mean(residuals ** 2))
    return residuals, relative_residuals, r2, rmsd
def computeGlobalFitMetrics( exp_y, fit_y):
    """
        计算全局拟合评价指标（global metrics）
        Parameters
        ----------
        exp_y : array-like
            所有实验数据（可由多温度、多等温线拼接而成）
        fit_y : array-like
            与 exp_y 一一对应的拟合数据
    """
    exp_y = np.asarray(
        exp_y,
        dtype=float
    )
    fit_y = np.asarray(
        fit_y,
        dtype=float
    )
    if exp_y.shape != fit_y.shape:
        raise ValueError(
            "exp_y and fit_y must "
            "have the same shape."
        )
    return calculateFitMetrics(
        y_true=exp_y,
        y_pred=fit_y
    )
def buildInterpolationFunction( x, y, method="pchip", extrapolate=False):
    """
    Build interpolation function.
    Parameters
    ----------
    x : array-like
    y : array-like
    method : str
        "pchip", "linear", "cubic"
    extrapolate : bool
    """
    if method.lower() == "pchip":
        return PchipInterpolator(x, y, extrapolate=extrapolate)
    elif method.lower() == "linear":
        return interp1d(
            x,
            y,
            kind="linear",
            bounds_error=False,
            fill_value=np.nan
        )
    elif method.lower() == "cubic":
        return CubicSpline(
            x,
            y,
            extrapolate=extrapolate
        )
    else:
        raise ValueError(f"Unknown interpolation method: {method}")
def interpolateData(
        x,
        y,
        x_interp=None,
        n_points=50,
        method="pchip",
        extrapolate=False):
    """
    Interpolate one dataset.
    Parameters
    ----------
    x : array-like
        Original x values.
    y : array-like
        Original y values.
    x_interp : array-like, optional
        Interpolation x values. If None, uniformly generate
        n_points within the data range.
    n_points : int, default=50
        Number of interpolation points when x_interp is None.
    method : {"pchip", "linear", "cubic"}, default="pchip"
    Returns
    -------
    x_interp : ndarray
        Interpolation x values.
    y_interp : ndarray
        Interpolated y values.
    interp_func : callable
        Interpolation function.
    """
    if x_interp is None:
        x_interp = np.linspace(
            np.min(x),
            np.max(x),
            n_points
        )
    interp_func = buildInterpolationFunction(
        x,
        y,
        method=method,
        extrapolate=extrapolate
    )
    y_interp = interp_func(x_interp)
    return x_interp, y_interp, interp_func
def multiDataInterpolation(
        x_data,
        y_data,
        x_common=None,
        n_points=50,
        method="pchip"):
    """
    Interpolate multiple datasets onto a common x grid.
    Parameters
    ----------
    x_data : dict
        Dictionary of x arrays.
    y_data : dict
        Dictionary of y arrays.
    x_common : array-like, optional
        Common interpolation coordinates.
        If None, they are automatically generated.
    n_points : int, default=50
        Number of interpolation points when x_common is None.
    method : {"pchip", "linear", "cubic"}, default="pchip"
    Returns
    -------
    df_interp : pandas.DataFrame
        Interpolated datasets.
        The first column is 'x_common'.
    interp_funcs : dict
        Dictionary of interpolation functions.
    """
    if x_data.keys() != y_data.keys():
        raise ValueError("x_data and y_data must have identical keys.")
    if x_common is None:
        x_common = dop.generateCommonX(
            x_data,
            n_points=n_points
        )
    result = {
        "x_common": x_common
    }
    interp_funcs = {}
    for key in x_data:
        _, y_interp, interp_func = interpolateData(
            x=x_data[key],
            y=y_data[key],
            x_interp=x_common,
            method=method
        )
        result[key] = y_interp
        interp_funcs[key] = interp_func
    df_interp = pd.DataFrame(result)
    return df_interp, interp_funcs
def fitLinear(
        x,
        y,
        x_name="x",
        y_name="y",
        dropna=True):
    """
    Perform simple linear regression.
    Model
    -----
        y = slope * x + intercept
    """
    # ======================================================
    # 1. Convert to numpy arrays
    # ======================================================
    x = np.asarray(
        x,
        dtype=float
    )
    y = np.asarray(
        y,
        dtype=float
    )
    # ======================================================
    # 2. Remove invalid values
    # ======================================================
    if dropna:
        mask = (
            np.isfinite(x) &
            np.isfinite(y)
        )
        x = x[mask]
        y = y[mask]
    if len(x) < 2:
        raise ValueError(
            "Not enough valid data."
        )
    # ======================================================
    # 3. Linear regression
    # ======================================================
    reg = linregress(
        x,
        y
    )
    # ======================================================
    # 4. Parameters
    # ======================================================
    parameters = {
        "slope": reg.slope,
        "intercept": reg.intercept
    }
    parameter_stderr = {
        "slope": reg.stderr,
        "intercept": getattr(
            reg,
            "intercept_stderr",
            np.nan
        )
    }
    # ======================================================
    # 5. Fitted values
    # ======================================================
    y_fit = linearModel(
        x,
        slope=reg.slope,
        intercept=reg.intercept
    )
    # ======================================================
    # 6. Common fit metrics
    # ======================================================
    fit_metrics = calculateFitMetrics(
        y_true=y,
        y_pred=y_fit
    )
    # ======================================================
    # 7. Return
    # ======================================================
    return {
        "x_name": x_name,
        "y_name": y_name,
        "N": len(x),
        "X": x,
        "Y": y,
        "Y_fit": y_fit,
        # ==================================================
        # Model
        # ==================================================
        "FitFunc": linearModel,
        # Backward compatibility
        "FitFunction": linearModel,
        # ==================================================
        # Parameters
        # ==================================================
        "ParameterNames": [
            "slope",
            "intercept"
        ],
        "Parameters": parameters,
        "ParameterStdErr": parameter_stderr,
        "ParameterCovariance": None,
        # ==================================================
        # Fit metrics
        # ==================================================
        **fit_metrics,
    }
def fitNonlinear(
        x,
        y,
        func,
        p0=None,
        bounds=(-np.inf, np.inf),
        param_names=None,
        x_name="x",
        y_name="y",
        dropna=True,
        maxfev=10000):
    """
    Perform nonlinear curve fitting.
    Parameters
    ----------
    x, y : array-like
        Input data.
    func : callable
        Nonlinear fitting function.
    p0 : array-like or None
        Initial parameter estimates.
    bounds : 2-tuple
        Lower and upper bounds for parameters.
    param_names : list or None
        Names of fitted parameters.
    x_name, y_name : str
        Variable names.
    dropna : bool
        Remove NaN and Inf values.
    maxfev : int
        Maximum number of function evaluations.
    """
    # ======================================================
    # 1. Convert to numpy arrays
    # ======================================================
    x = np.asarray(x, dtype=float).ravel()
    y = np.asarray(y, dtype=float).ravel()
    # ======================================================
    # 2. Remove invalid values
    # ======================================================
    if dropna:
        mask = ( np.isfinite(x) & np.isfinite(y) )
        x = x[mask]
        y = y[mask]
    if len(x) < 2:
        raise ValueError( "Not enough valid data." )
    # ======================================================
    # 3. Debug information
    # ======================================================
    # print("\n" + "=" * 70)
    # print("Nonlinear fitting")
    # print("=" * 70)
    # print(f"N       : {len(x)}")
    # print(f"x range : {x.min()} ~ {x.max()}")
    # print(f"y range : {y.min()} ~ {y.max()}")
    # print(f"p0      : {p0}")
    # print(f"bounds  : {bounds}")
    # ======================================================
    # 4. Nonlinear fitting
    # ======================================================
    popt, pcov = curve_fit( func, x, y, p0=p0, bounds=bounds, maxfev=maxfev )
    # ======================================================
    # 5. Parameter names
    # ======================================================
    n_parameters = len(popt)
    if param_names is None:
        param_names = [ f"param_{i + 1}" for i in range(n_parameters) ]
    if len(param_names) != n_parameters:
        raise ValueError( "Length of param_names must match the number of fitted parameters." )
    # ======================================================
    # 6. Parameters
    # ======================================================
    parameters = dict( zip(param_names, popt) )
    # ======================================================
    # 7. Parameter uncertainty
    # ======================================================
    if pcov is not None:
        parameter_stderr = np.sqrt( np.diag(pcov) )
    else:
        parameter_stderr = np.full( n_parameters, np.nan )
    parameter_stderr_dict = dict( zip( param_names, parameter_stderr ) )
    # ======================================================
    # 8. Fitted values
    # ======================================================
    y_fit = func( x, *popt )
    # ======================================================
    # 9. Common fit metrics
    # ======================================================
    fit_metrics = calculateFitMetrics( y_true=y, y_pred=y_fit )
    # ======================================================
    # 10. Return
    # ======================================================
    return {
        "x_name": x_name,
        "y_name": y_name,
        "N": len(x),
        "X": x,
        "Y": y,
        "Y_fit": y_fit,
        # ==================================================
        # Model
        # ==================================================
        "FitFunc": func,
        # Backward compatibility
        "FitFunction": func,
        # ==================================================
        # Parameters
        # ==================================================
        "ParameterNames": param_names,
        "Parameters": parameters,
        "ParameterStdErr":
            parameter_stderr_dict,
        "ParameterCovariance":
            pcov,
        # ==================================================
        # Fit metrics
        # ==================================================
        **fit_metrics,
    }
def fittingResultToRecord00( result, parameter_name):
    """
    Convert a fitting result dictionary
    into one flat record.
    """
    record = {
        "Parameter": parameter_name,
        "N": result["N"],
        # ==================================================
        # Model
        # ==================================================
        "Model": getattr(
            result["FitFunc"],
            "model_name",
            result["FitFunc"].__name__
        ),
        "Model_Equation": getattr(
            result["FitFunc"],
            "equation",
            ""
        ),
        # ==================================================
        # Correlation
        # ==================================================
        "Pearson_r": result.get( "Pearson_r", np.nan ),
        "Pearson_p": result.get( "Pearson_p", np.nan ),
        "Spearman_r": result.get( "Spearman_r", np.nan ),
        "Spearman_p": result.get( "Spearman_p", np.nan ),
        # ==================================================
        # Goodness of fit
        # ==================================================
        "R2": result.get( "R2", np.nan ),
        "RMSE": result.get( "RMSE", np.nan ),
        "RMSE_rel": result.get( "RMSE_rel", np.nan ),
        "RMSE_nrm": result.get( "RMSE_nrm", np.nan ),
        "MAE": result.get( "MAE", np.nan ),
        "MAPE": result.get( "MAPE", np.nan ),
    }
    # ======================================================
    # Parameters
    # ======================================================
    for name, value in result.get( "Parameters", {} ).items():
        record[f"Fit_{name}"] = value
    # ======================================================
    # Parameter standard errors
    # ======================================================
    for name, value in result.get( "ParameterStdErr", {} ).items():
        record[ f"Fit_{name}_StdErr" ] = value
    return record
def fittingResultToRecord( result, parameter_name=None, extra=None, n_name="N"):
    """
    Convert an analysis/fitting result dictionary
    into one flat record.
    Compatible with:
        - correlationAnalysis()
        - fitLinear()
        - fitNonlinear()
    Parameters
    ----------
    result : dict
        Result dictionary.
    parameter_name : str, optional
        Optional parameter/analysis label.
    extra : dict, optional
        Additional fields to add to the record,
        e.g. pressure.
    n_name : str
        Name used for the sample-size field.
        Default: "N".
        Can be "N_pairs" for differential analysis.
    Returns
    -------
    dict
        Flat record suitable for DataFrame construction.
    """
    record = {}
    # ======================================================
    # Extra information
    # ======================================================
    if extra is not None:
        record.update(extra)
    # ======================================================
    # Parameter / analysis name
    # ======================================================
    if parameter_name is not None:
        record["Parameter"] = parameter_name
    # ======================================================
    # Sample size
    # ======================================================
    if "N" in result:
        record[n_name] = result["N"]
    # ======================================================
    # Correlation
    # ======================================================
    for name in [
        "Pearson_r",
        "Pearson_p",
        "Spearman_r",
        "Spearman_p",
    ]:
        if name in result:
            record[name] = result[name]
    # ======================================================
    # Goodness of fit
    # ======================================================
    for name in [
        "SSE",
        "R2",
        "RMSE",
        "RMSE_rel",
        "RMSE_nrm",
        "MAE",
        "MAPE",
    ]:
        if name in result:
            record[name] = result[name]
    # ======================================================
    # Model information
    # ======================================================
    fit_func = result.get( "FitFunc", result.get("FitFunction", None) )
    if fit_func is not None:
        record["Model"] = getattr(
            fit_func,
            "model_name",
            fit_func.__name__
        )
        record["Model_Equation"] = getattr(
            fit_func,
            "equation",
            ""
        )
    # ======================================================
    # Fitted parameters
    # ======================================================
    parameters = result.get( "Parameters", {} )
    for name, value in parameters.items():
        record[ f"Fit_{name}" ] = value
    # ======================================================
    # Parameter uncertainty
    # ======================================================
    parameter_stderr = result.get( "ParameterStdErr", {} )
    for name, value in parameter_stderr.items():
        record[ f"Fit_{name}_StdErr" ] = value
    return record
def calculateCorrelation(
        x,
        y,
        x_name="x",
        y_name="y",
        dropna=True):
    """
    Calculate Pearson and Spearman correlations.
    Parameters
    ----------
    x : array-like
        First variable.
    y : array-like
        Second variable.
    x_name : str
        Name of x.
    y_name : str
        Name of y.
    dropna : bool
        Remove NaN and infinite pairs.
    Returns
    -------
    dict
        Correlation results.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    # ======================================================
    # Remove invalid values
    # ======================================================
    if dropna:
        mask = np.isfinite(x) & np.isfinite(y)
        x = x[mask]
        y = y[mask]
    if len(x) < 2:
        raise ValueError("Not enough valid data.")
    # ======================================================
    # Pearson correlation
    # ======================================================
    pearson_r, pearson_p = pearsonr(x, y)
    # ======================================================
    # Spearman correlation
    # ======================================================
    spearman_r, spearman_p = spearmanr(x, y)
    return {
        "x_name": x_name,
        "y_name": y_name,
        "N": len(x),
        "X": x,
        "Y": y,
        "Pearson_r": pearson_r,
        "Pearson_p": pearson_p,
        "Spearman_r": spearman_r,
        "Spearman_p": spearman_p,
    }
def linearModel(x, slope, intercept):
    """
    Linear model.
    Mathematical form
    ------------------
    y = slope * x + intercept
    Parameters
    ----------
    x : array-like
        Independent variable.
    slope : float
        Slope of the linear relationship.
    intercept : float
        Intercept of the linear relationship.
    Returns
    -------
    y : float or np.ndarray
        Predicted value.
    """
    return slope * x + intercept
linearModel.model_name = "Linear"
linearModel.equation = ( "y = slope * x + intercept" )
def exponentialSaturation(x, a_inf, A, k):
    """
    Exponential saturation model.
    Mathematical form
    ------------------
    y = a_inf - A * exp(-k * x)
    Parameters
    ----------
    x : array-like
        Independent variable.
    a_inf : float
        Asymptotic value as x approaches infinity.
    A : float
        Initial deviation from the asymptotic value.
    k : float
        Saturation rate constant.
    Returns
    -------
    y : float or np.ndarray
        Predicted value.
    """
    return a_inf - A * np.exp(-k * x)
exponentialSaturation.model_name = "Exponential Saturation"
exponentialSaturation.equation = ( "y = a_inf - A * exp(-k * x)" )
def quadratic(x, a, b, c):
    """
    Quadratic model.
    Mathematical form
    ------------------
    y = a * x**2 + b * x + c
    Parameters
    ----------
    x : array-like
        Independent variable.
    a : float
        Quadratic coefficient.
    b : float
        Linear coefficient.
    c : float
        Intercept.
    Returns
    -------
    y : float or np.ndarray
        Predicted value.
    """
    return a * x**2 + b * x + c
quadratic.model_name = "Quadratic"
quadratic.equation = ( "y = a * x^2 + b * x + c" )
# ==========================================================
# 1. Model library
# ==========================================================
model_library = {
    "Linear": {
        "fit_method": "linear",
        "model": None,
        "p0": None,
        "bounds": (-np.inf, np.inf),
        "param_names": None
    },
    "Quadratic": {
        "fit_method": "nonlinear",
        "model": quadratic,
        "p0": [0.0, 0.0, 0.0],
        "bounds": (-np.inf, np.inf),
        "param_names": [
            "a",
            "b",
            "c"
        ]
    },
    "Exponential Saturation": {
        "fit_method": "nonlinear",
        "model": exponentialSaturation,
        "p0": [0.1, 0.1, 0.01],
        "bounds": (
            [-np.inf, -np.inf, 0.0],
            [np.inf, np.inf, np.inf]
        ),
        "param_names": [
            "a_inf",
            "A",
            "k"
        ]
    }
}