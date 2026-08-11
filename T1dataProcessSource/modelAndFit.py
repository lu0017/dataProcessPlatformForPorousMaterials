from common import *
def calculateFitMetrics(y_true, y_pred):
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)

    if y_true.shape != y_pred.shape:
        raise ValueError(
            "y_true and y_pred must have the same shape."
        )

    residual = y_true - y_pred

    ss_res = np.sum(residual ** 2)
    ss_tot = np.sum(
        (y_true - np.mean(y_true)) ** 2
    )

    if ss_tot > 0:
        r2 = 1 - ss_res / ss_tot
    else:
        r2 = np.nan

    rmse = np.sqrt(
        np.mean(residual ** 2)
    )

    y_mean = np.mean(y_true)

    rmse_rel = (
        rmse /
        (np.abs(y_mean) + 1e-12)
    )

    y_range = (
        np.max(y_true) -
        np.min(y_true)
    )

    rmse_nrm = (
        rmse /
        (y_range + 1e-12)
    )

    mae = np.mean(
        np.abs(residual)
    )

    with np.errstate(
        divide="ignore",
        invalid="ignore"
    ):
        relative_residual = np.where(
            y_true != 0,
            residual / y_true,
            0.0
        )

    mape = np.mean(
        np.abs(relative_residual)
    ) * 100

    return {
        "Residual": residual,
        "RelativeResidual": relative_residual,
        "R2": r2,
        "RMSE": rmse,
        "RMSE_rel": rmse_rel,
        "RMSE_nrm": rmse_nrm,
        "MAE": mae,
        "MAPE": mape,
    }
def compute_fit_metrics(fit_x, fit_y, exp_x, exp_y, method='linear', use_interp=True,):
    # from sklearn.metrics import r2_score
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
def computeFitMetrics(
        fit_x,
        fit_y,
        exp_x,
        exp_y,
        method="linear",
        use_interp=True):
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
    fit_x = np.asarray(fit_x, dtype=float)
    fit_y = np.asarray(fit_y, dtype=float)

    exp_x = np.asarray(exp_x, dtype=float)
    exp_y = np.asarray(exp_y, dtype=float)

    # ======================================================
    # Align x coordinates
    # ======================================================

    if use_interp:

        interp_func = interp1d(
            fit_x,
            fit_y,
            kind=method,
            fill_value="extrapolate"
        )

        fit_y_aligned = interp_func(exp_x)

    else:

        if fit_y.shape != exp_y.shape:
            raise ValueError(
                "fit_y and exp_y must have "
                "the same shape when "
                "use_interp=False."
            )

        fit_y_aligned = fit_y

    # ======================================================
    # Common metrics
    # ======================================================

    return calculateFitMetrics(
        y_true=exp_y,
        y_pred=fit_y_aligned
    )
def compute_global_fit_metrics(exp_y, fit_y):
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
        dropna=True):
    """
    Perform nonlinear curve fitting.
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
    # 3. Nonlinear fitting
    # ======================================================

    popt, pcov = curve_fit(
        func,
        x,
        y,
        p0=p0,
        bounds=bounds
    )

    # ======================================================
    # 4. Parameter names
    # ======================================================

    n_parameters = len(popt)

    if param_names is None:

        param_names = [
            f"param_{i + 1}"
            for i in range(
                n_parameters
            )
        ]

    if len(param_names) != n_parameters:

        raise ValueError(
            "Length of param_names must "
            "match the number of fitted "
            "parameters."
        )

    # ======================================================
    # 5. Parameters
    # ======================================================

    parameters = dict(
        zip(
            param_names,
            popt
        )
    )

    # ======================================================
    # 6. Parameter uncertainty
    # ======================================================

    if pcov is not None:

        parameter_stderr = np.sqrt(
            np.diag(pcov)
        )

    else:

        parameter_stderr = np.full(
            n_parameters,
            np.nan
        )

    parameter_stderr_dict = dict(
        zip(
            param_names,
            parameter_stderr
        )
    )

    # ======================================================
    # 7. Fitted values
    # ======================================================

    y_fit = func(
        x,
        *popt
    )

    # ======================================================
    # 8. Common fit metrics
    # ======================================================

    fit_metrics = calculateFitMetrics(
        y_true=y,
        y_pred=y_fit
    )

    # ======================================================
    # 9. Return
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

        "ParameterCovariance": pcov,

        # ==================================================
        # Fit metrics
        # ==================================================

        **fit_metrics,
    }

def fittingResultToRecord00(
        result,
        parameter_name):
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

        "Pearson_r": result.get(
            "Pearson_r",
            np.nan
        ),

        "Pearson_p": result.get(
            "Pearson_p",
            np.nan
        ),

        "Spearman_r": result.get(
            "Spearman_r",
            np.nan
        ),

        "Spearman_p": result.get(
            "Spearman_p",
            np.nan
        ),

        # ==================================================
        # Goodness of fit
        # ==================================================

        "R2": result.get(
            "R2",
            np.nan
        ),

        "RMSE": result.get(
            "RMSE",
            np.nan
        ),

        "RMSE_rel": result.get(
            "RMSE_rel",
            np.nan
        ),

        "RMSE_nrm": result.get(
            "RMSE_nrm",
            np.nan
        ),

        "MAE": result.get(
            "MAE",
            np.nan
        ),

        "MAPE": result.get(
            "MAPE",
            np.nan
        ),
    }

    # ======================================================
    # Parameters
    # ======================================================

    for name, value in result.get(
        "Parameters",
        {}
    ).items():

        record[f"Fit_{name}"] = value

    # ======================================================
    # Parameter standard errors
    # ======================================================

    for name, value in result.get(
        "ParameterStdErr",
        {}
    ).items():

        record[
            f"Fit_{name}_StdErr"
        ] = value

    return record
def fittingResultToRecord(
        result,
        parameter_name=None,
        extra=None,
        n_name="N"):
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

    fit_func = result.get(
        "FitFunc",
        result.get("FitFunction", None)
    )

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

    parameters = result.get(
        "Parameters",
        {}
    )

    for name, value in parameters.items():

        record[
            f"Fit_{name}"
        ] = value

    # ======================================================
    # Parameter uncertainty
    # ======================================================

    parameter_stderr = result.get(
        "ParameterStdErr",
        {}
    )

    for name, value in parameter_stderr.items():

        record[
            f"Fit_{name}_StdErr"
        ] = value

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
