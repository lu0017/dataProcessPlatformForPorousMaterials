from common import *
import constantsAndName as const
import T1dataProcessSource.modelAndFit as mf
def extract_temperature(col_name):
    """Extract numeric temperature from a column name like '15Absolute Pressure (kPa)'."""
    match = re.match(r"^(-?\d+\.?\d*)", col_name)
    if match:
        return float(match.group(1))
    return None
def getName(num, unitChange=0.0, prefix=None, suffix=None, decimals=None):
    """
    生成列名，支持前缀、后缀和小数点控制
    参数：
    - num: 数字
    - unitChange: 对 num 的偏移量
    - prefix: 列名前缀
    - suffix: 列名后缀
    - decimals: None 表示整数不保留小数位；整数 n 表示保留 n 位小数
    """
    value = num + unitChange
    # 根据 decimals 决定格式
    if decimals is None:
        # 整数时去掉小数
        if float(value).is_integer():
            value_str = f"{int(value)}"
        else:
            value_str = str(value)
    else:
        value_str = f"{value:.{decimals}f}"
    # 构造列名
    col_name = ""
    if prefix:
        col_name += f"{prefix}"
    col_name += value_str
    if suffix:
        col_name += f"{suffix}"
    return col_name
def crop_low_points(rawData, n_remove=1):
    """
    裁剪原始数据，删除每个温度下最低 n 个压力点
    输入：
        rawData: DataFrame，包含列 ["P_Pa","q_mmol_g","T_K"]
        n_remove: 每个温度下删除最低 pressure 点数量
    输出：
        cropped_data: DataFrame，裁剪后的数据
    """
    # 按温度分组，删除最低 n 个压力点
    cropped_list = []
    for T, df in rawData.groupby("T_K"):
        df_sorted = df.sort_values("P_Pa").reset_index(drop=True)
        if len(df_sorted) > n_remove:
            df_cropped = df_sorted.iloc[n_remove:].copy()
        else:
            df_cropped = df_sorted.copy()  # 数据点不足时不裁剪
        cropped_list.append(df_cropped)
    cropped_data = pd.concat(cropped_list, ignore_index=True)
    return cropped_data
def preprocess_low_pressure(rawData, P_max_Pa=5000, min_points=5):
    """
    从原始 rawData 中筛选低压段用于 DSL 拟合
    输入：
        rawData: DataFrame，含列 ["P_Pa","q_mmol_g","T_K"]
        P_max_Pa: 低压筛选的最大压力（Pa）
        min_points: 每个温度至少需要的点数，少于这个数会被剔除
    输出：
        lowP_df: 筛选后的 DataFrame
    """
    df = rawData.copy()
    # 1. 过滤低压数据
    df = df[df["P_Pa"] <= P_max_Pa]
    # 2. 删除异常值：负压力、负吸附量
    df = df[(df["P_Pa"] > 0) & (df["q_mmol_g"] >= 0)]
    # 3. 可能存在重复压力 -> 取平均
    df = (df.groupby(["T_K", "P_Pa"], as_index=False)
            .agg({"q_mmol_g":"mean"}))
    # 4. 每个温度至少需要 min_points 个点，否则 DSL 会拟合失败
    valid_T = df.groupby("T_K").filter(lambda x: len(x) >= min_points)
    if valid_T.empty:
        raise ValueError(f"低压范围( <= {P_max_Pa} Pa)内没有足够的数据用于拟合！")
    return valid_T.reset_index(drop=True)
def organize_isotherm_data(df_dict):
    """
    从你给的 df 字典（列名: 列数据）中自动按温度整理 Pressure 和 uptake。
    返回结构化数据:
    {
        T1: {"P": [...], "q": [...]},
        T2: {"P": [...], "q": [...]},
        ...
    }
    """
    result = {}
    for col_name, values in df_dict.items():
        T = extract_temperature(col_name)
        if T is None:
            continue
        if T not in result:
            result[T] = {"P": None, "q": None}
        if "Pressure" in col_name:
            result[T]["P"] = values
        elif "uptake" in col_name or "mmol" in col_name:
            result[T]["q"] = values
    T_C = sorted(list(result.keys()))
    return result, T_C
def build_P_and_uptake_data(df_dict):
    organized, T_C = organize_isotherm_data(df_dict)
    P_data = {}
    uptake_data = {}
    for T in T_C:
        # 将摄氏温度转换为开尔文
        T_K = T + 273.15  
        P_data[T_K] = organized[T]["P"]
        uptake_data[T_K] = organized[T]["q"]
    T_K = celsius_to_kelvin(T_C, unit="C")
    return P_data, uptake_data, T_K
def buildPeakModel( peaks, model="Lorentzian", ):
    """
    Build a combined peak model.
    Parameters
    ----------
    peaks : list of dict
        Peak definitions.
        Example
        -------
        peaks = [
            {"name": "D"},
            {"name": "G"},
        ]
    model : str
        "Lorentzian"
        "Gaussian"
        "Voigt"
        "PseudoVoigt"
    Returns
    -------
    lmfit.Model
        Combined model.
    """
    modelMap = {
        "Lorentzian": LorentzianModel,
        "Gaussian": GaussianModel,
        "Voigt": VoigtModel,
        "PseudoVoigt": PseudoVoigtModel,
    }
    if model not in modelMap:
        raise ValueError(
            f"Unknown peak model: {model}"
        )
    ModelClass = modelMap[model]
    peakModel = None
    for peak in peaks:
        prefix = peak["name"] + "_"
        m = ModelClass(
            prefix=prefix,
        )
        if peakModel is None:
            peakModel = m
        else:
            peakModel += m
    return peakModel
import numpy as np
def guessPeakParameters(
    peakModel,
    x,
    y,
    peaks,
):
    """
    Generate initial fitting parameters.
    Parameters
    ----------
    peakModel : lmfit.Model
        Peak model returned by buildPeakModel().
    x : ndarray
        X values.
    y : ndarray
        Y values.
    peaks : list of dict
        Peak definitions.
    Returns
    -------
    lmfit.Parameters
        Initial parameters.
    """
    x = np.asarray(x)
    y = np.asarray(y)
    # 利用模型创建 Parameters
    params = peakModel.make_params()
    ymax = np.max(y)
    for peak in peaks:
        prefix = peak["name"] + "_"
        # -----------------------------
        # center
        # -----------------------------
        center = peak["center"]
        params[prefix + "center"].set(
            value=center,
            min=peak.get("center_min", x.min()),
            max=peak.get("center_max", x.max()),
        )
        # -----------------------------
        # sigma
        # -----------------------------
        sigma = peak.get("sigma", 50)
        params[prefix + "sigma"].set(
            value=sigma,
            min=peak.get("sigma_min", 1),
            max=peak.get("sigma_max", 500),
        )
        # -----------------------------
        # amplitude (area)
        # -----------------------------
        amp = peak.get("amplitude")
        if amp is None:
            idx = np.argmin(np.abs(x - center))
            amp = abs(y[idx]) * sigma
            if amp <= 0:
                amp = ymax * sigma
        params[prefix + "amplitude"].set(
            value=amp,
            min=0,
        )
    return params
def calculatePeakProperties( result, peaks, ):
    """
    Extract peak properties from lmfit result.
    Parameters
    ----------
    result : lmfit.ModelResult
        Peak fitting result.
    peaks : list of dict
        Peak definitions.
    Returns
    -------
    dict
        Peak properties.
    """
    properties = {}
    params = result.params
    for peak in peaks:
        name = peak["name"]
        prefix = name + "_"
        p = {}
        # -----------------------
        # Center
        # -----------------------
        p["center"] = params[prefix + "center"].value
        # -----------------------
        # Sigma
        # -----------------------
        p["sigma"] = params[prefix + "sigma"].value
        # -----------------------
        # FWHM
        # -----------------------
        if prefix + "fwhm" in params:
            p["fwhm"] = params[prefix + "fwhm"].value
        else:
            p["fwhm"] = np.nan
        # -----------------------
        # Height
        # -----------------------
        if prefix + "height" in params:
            p["height"] = params[prefix + "height"].value
        else:
            p["height"] = np.nan
        # -----------------------
        # Area
        # -----------------------
        p["area"] = params[prefix + "amplitude"].value
        properties[name] = p
    return properties
def fitPeak( x, y, peaks, model="Lorentzian", method="leastsq", weights=None, fit_kws=None, ):
    """
    Fit multiple peaks.
    Parameters
    ----------
    x : ndarray
        X values.
    y : ndarray
        Y values.
    peaks : list of dict
        Peak definitions.
    model : str
        Peak model.
    method : str
        lmfit fitting method.
    weights : ndarray, optional
        Weight array.
    fit_kws : dict, optional
        Extra keyword arguments passed to lmfit.
    Returns
    -------
    ModelResult
        lmfit fitting result.
    """
    # -----------------------------
    # Build model
    # -----------------------------
    peakModel = buildPeakModel(
        peaks=peaks,
        model=model,
    )
    # -----------------------------
    # Initial parameters
    # -----------------------------
    params = guessPeakParameters(
    peakModel,
    x,
    y,
    peaks,
    )
    # -----------------------------
    # Fit
    # -----------------------------
    result = peakModel.fit(
        y,
        params=params,
        x=x,
        method=method,
        weights=weights,
        **({} if fit_kws is None else fit_kws),
    )
    return result
def calculateFitMetrics(xTrue, yTrue, xPred, yPred, allow_extrapolation=False):
    """
    Calculate R² and RMSE between experimental and predicted curves.
    If the x coordinates differ (including different numbers of points),
    the predicted curve is first aligned to the experimental x coordinates.
    By default, only the overlapping x range is used for evaluation.
    If ``allow_extrapolation=True``, linear extrapolation is performed
    outside the prediction range.
    Parameters
    ----------
    xTrue : array-like
        Experimental x values.
    yTrue : array-like
        Experimental y values.
    xPred : array-like
        Predicted x values.
    yPred : array-like
        Predicted y values.
    allow_extrapolation : bool, default=False
        Whether to allow linear extrapolation when the experimental
        x range extends beyond the prediction range.
        * False : evaluate only within the overlapping x range.
        * True  : linearly extrapolate the prediction to all
          experimental x values.
    Returns
    -------
    dict
        Dictionary containing:
            - R2 : coefficient of determination.
            - RMSE : root mean square error.
            - x : x values used for evaluation.
            - yTrue : experimental values used for evaluation.
            - yPred : predicted values (interpolated/extrapolated if needed).
            - n : number of points used for evaluation.
    """
    # -----------------------------
    # Convert to numpy arrays
    # -----------------------------
    xTrue = np.asarray(xTrue, dtype=float).ravel()
    yTrue = np.asarray(yTrue, dtype=float).ravel()
    xPred = np.asarray(xPred, dtype=float).ravel()
    yPred = np.asarray(yPred, dtype=float).ravel()
    # -----------------------------
    # Check lengths
    # -----------------------------
    if xTrue.size != yTrue.size:
        raise ValueError(
            f"xTrue ({xTrue.size}) and yTrue ({yTrue.size}) "
            "must have the same length."
        )
    if xPred.size != yPred.size:
        raise ValueError(
            f"xPred ({xPred.size}) and yPred ({yPred.size}) "
            "must have the same length."
        )
    if xTrue.size == 0:
        raise ValueError("Experimental data is empty.")
    if xPred.size == 0:
        raise ValueError("Predicted data is empty.")
    # -----------------------------
    # Check finite values
    # -----------------------------
    if not np.all(np.isfinite(xTrue)):
        raise ValueError("xTrue contains NaN or Inf.")
    if not np.all(np.isfinite(yTrue)):
        raise ValueError("yTrue contains NaN or Inf.")
    if not np.all(np.isfinite(xPred)):
        raise ValueError("xPred contains NaN or Inf.")
    if not np.all(np.isfinite(yPred)):
        raise ValueError("yPred contains NaN or Inf.")
    # -----------------------------
    # Sort by x
    # -----------------------------
    idx = np.argsort(xTrue)
    xTrue = xTrue[idx]
    yTrue = yTrue[idx]
    idx = np.argsort(xPred)
    xPred = xPred[idx]
    yPred = yPred[idx]
    # -----------------------------
    # Interpolate if necessary
    # -----------------------------
    if (
        xTrue.size != xPred.size
        or not np.allclose(xTrue, xPred)
    ):
        if allow_extrapolation:
            # Linear interpolation + extrapolation
            interp_func = interp1d(
                xPred,
                yPred,
                kind="linear",
                fill_value="extrapolate",
                bounds_error=False,
            )
            xEval = xTrue
            yTrueEval = yTrue
            yPredInterp = interp_func(xEval)
        else:
            # Use only the overlapping x range
            xmin = max(xTrue.min(), xPred.min())
            xmax = min(xTrue.max(), xPred.max())
            mask = (xTrue >= xmin) & (xTrue <= xmax)
            xEval = xTrue[mask]
            yTrueEval = yTrue[mask]
            if xEval.size < 2:
                raise ValueError(
                    "Less than two overlapping points between "
                    "experimental and predicted data."
                )
            yPredInterp = np.interp(xEval, xPred, yPred)
    else:
        xEval = xTrue
        yTrueEval = yTrue
        yPredInterp = yPred
    # -----------------------------
    # Calculate metrics
    # -----------------------------
    return {
        "R2": r2_score(yTrueEval, yPredInterp),
        "RMSE": np.sqrt(mean_squared_error(yTrueEval, yPredInterp)),
        "x": xEval,
        "yTrue": yTrueEval,
        "yPred": yPredInterp,
        "n": xEval.size,
    }
def copySamples(data, sampleMap):
    """
    根据 sampleMap 拷贝指定样品。
    Parameters
    ----------
    data : dict
        {
            sampleName: DataFrame
        }
    sampleMap : dict
        {
            newName: [oldName1, oldName2, ...]
        }
    Returns
    -------
    dict
    """
    result = {}
    for oldNames in sampleMap.values():
        for oldName in oldNames:
            if oldName in data:
                result[oldName] = data[oldName].copy()
    return result
def renameSamples(data, sampleMap):
    """
    根据 sampleMap 重命名样品。
    Parameters
    ----------
    data : dict
        {
            sampleName: DataFrame
        }
    sampleMap : dict
        {
            newName: [oldName1, oldName2, ...]
        }
    Returns
    -------
    dict
    """
    result = {}
    for newName, oldNames in sampleMap.items():
        for oldName in oldNames:
            if oldName in data:
                result[newName] = data[oldName].copy()
    return result
def swapDictLevels(data, keep_keys=None):
    """
    Swap the first and second levels of a nested dictionary.
    Parameters
    ----------
    data : dict
        Nested dictionary with two levels.
    keep_keys : list or tuple or set, optional
        Keys that should remain at the first level without swapping.
    Returns
    -------
    dict
        Dictionary with swapped first and second levels.
    Example
    -------
    Before
    {
        "S1": {
            "A": 1,
            "B": 2,
            "metadata": {...}
        },
        "S2": {
            "A": 3,
            "B": 4,
            "metadata": {...}
        }
    }
    After
    {
        "A": {
            "S1": 1,
            "S2": 3
        },
        "B": {
            "S1": 2,
            "S2": 4
        },
        "metadata": {
            "S1": {...},
            "S2": {...}
        }
    }
    """
    if keep_keys is None:
        keep_keys = []
    result = {}
    for outer_key, inner_dict in data.items():
        for inner_key, value in inner_dict.items():
            if inner_key in keep_keys:
                result.setdefault(inner_key, {})
                result[inner_key][outer_key] = deepcopy(value)
                continue
            result.setdefault(inner_key, {})
            result[inner_key][outer_key] = deepcopy(value)
    return result
def smooth( y, method="savgol", **kwargs, ):
    """
    Smooth one spectrum.
    Parameters
    ----------
    y : ndarray
        Spectrum intensity.
    method : str
        "savgol"
        "moving"
        "gaussian"
    **kwargs
        Method-specific parameters.
    Returns
    -------
    ndarray
        Smoothed spectrum.
    """
    y = np.asarray(y, dtype=float)
    # -----------------------------
    # Savitzky-Golay
    # -----------------------------
    if method == "savgol":
        window_length = kwargs.get("window_length", 9)
        polyorder = kwargs.get("polyorder", 3)
        # window 必须是奇数
        if window_length % 2 == 0:
            window_length += 1
        # window 不能超过数据长度
        if window_length >= len(y):
            window_length = len(y) - 1
            if window_length % 2 == 0:
                window_length -= 1
        return savgol_filter(
            y,
            window_length=window_length,
            polyorder=polyorder,
        )
    # -----------------------------
    # Moving average
    # -----------------------------
    elif method == "moving":
        window = kwargs.get("window", 5)
        kernel = np.ones(window) / window
        return np.convolve(
            y,
            kernel,
            mode="same",
        )
    # -----------------------------
    # Gaussian
    # -----------------------------
    elif method == "gaussian":
        sigma = kwargs.get("sigma", 1)
        return gaussian_filter1d(
            y,
            sigma=sigma,
        )
    else:
        raise ValueError(f"Unknown smoothing method: {method}")
def smoothSpectrum( data, y="Absorbance", method="savgol", **kwargs, ):
    """
    Smooth spectra.
    Parameters
    ----------
    data : dict
        {sample: DataFrame}
    y : str
        Y column name.
    method : str
        Smoothing method.
    Returns
    -------
    dict
        Smoothed spectra.
    """
    result = {}
    for sample, df in data.items():
        newdf = df.copy()
        newdf[y] = smooth(
            newdf[y].to_numpy(),
            method=method,
            **kwargs,
        )
        result[sample] = newdf
    return result
def normalize( y, method="vector", ):
    """
    Normalize one spectrum.
    Parameters
    ----------
    y : ndarray
    method : str
        "minmax"
        "max"
        "vector"
    Returns
    -------
    ndarray
    """
    y = np.asarray(y, dtype=float)
    if method == "minmax":
        ymin = np.min(y)
        ymax = np.max(y)
        if ymax > ymin:
            return (y - ymin) / (ymax - ymin)
        return np.zeros_like(y)
    elif method == "max":
        ymax = np.max(np.abs(y))
        if ymax != 0:
            return y / ymax
        return y
    elif method == "vector":
        norm = np.linalg.norm(y)
        if norm != 0:
            return y / norm
        return y
    else:
        raise ValueError(f"Unknown normalization method: {method}")
def normalizeSpectrum( data, y="Absorbance", method="vector", ):
    """
    Normalize spectra in a dictionary.
    """
    result = {}
    for sample, df in data.items():
        newdf = df.copy()
        newdf[y] = normalize(
            newdf[y].to_numpy(),
            method=method,
        )
        result[sample] = newdf
    return result
def baseline(x, y, method="rubberband", **kwargs):
    """
    Baseline correction.
    Parameters
    ----------
    x : ndarray
    y : ndarray
    method : str
        None
        "rubberband"
        "asls"
        "airpls"
    Returns
    -------
    ndarray
        Baseline corrected spectrum.
    """
    if method is None:
        return y
    fitter = Baseline(x)
    method = method.lower()
    if method == "rubberband":
        baseline, _ = fitter.rubberband(y, **kwargs)
    elif method == "asls":
        baseline, _ = fitter.asls(
            y,
            lam=1e6,
            p=0.001,
            **kwargs,
        )
    elif method == "airpls":
        baseline, _ = fitter.airpls(y, **kwargs)
    else:
        raise ValueError(
            f"Unknown baseline method: {method}"
        )
    return y - baseline
def baselineCorrection( data, x="Wavenumber", y="Absorbance", method="asls", ):
    """
    Baseline correction for spectra.
    """
    result = {}
    for sample, df in data.items():
        newdf = df.copy()
        newdf[y] = baseline(
            newdf[x].to_numpy(),
            newdf[y].to_numpy(),
            method=method,
        )
        result[sample] = newdf
    return result
def generateCommonX(x_data, n_points=50):
    """
    Generate common x coordinates for multiple datasets.
    The common interpolation range is determined by the overlap
    of all datasets.
    Parameters
    ----------
    x_data : dict
        Dictionary of x arrays.
        Example:
            {
                "25": x1,
                "35": x2,
                ...
            }
    n_points : int, default=50
        Number of interpolation points.
    Returns
    -------
    x_common : ndarray
        Common interpolation coordinates.
    """
    x_min = max(np.min(x) for x in x_data.values())
    x_max = min(np.max(x) for x in x_data.values())
    if x_min >= x_max:
        raise ValueError(
            "No overlapping x range exists among the datasets."
        )
    return np.linspace(x_min, x_max, n_points)
def cropX(data, xmin=None, xmax=None):
    """
    根据 X 范围截取数据。
    Parameters
    ----------
    data : dict
        {
            sampleName: DataFrame
        }
    xmin : float, optional
        最小 X 值。
    xmax : float, optional
        最大 X 值。
    Returns
    -------
    dict
    """
    result = {}
    for name, df in data.items():
        x = df.iloc[:, 0]
        mask = pd.Series(True, index=df.index)
        if xmin is not None:
            mask &= (x >= xmin)
        if xmax is not None:
            mask &= (x <= xmax)
        result[name] = df.loc[mask].reset_index(drop=True)
    return result
def clean_dataframe(df, x_col, y_col=None, scale_x=1.0, scale_y=1.0):
    """
    清理 DataFrame 数据，去掉 NaN 和 Inf。
    参数:
        df : pandas DataFrame
        x_col : str，x 列名（通常是压力）
        y_col : str，可选，y 列名（通常是吸附量）
        scale_x : float，x 列缩放系数，默认 1.0
        scale_y : float，y 列缩放系数，默认 1.0
    返回:
        如果 y_col 提供，返回 (x_clean, y_clean) numpy array
        如果 y_col 未提供，只返回 x_clean numpy array
    """
    if x_col not in df.columns:
        raise KeyError(f"列 {x_col} 不存在于 DataFrame 中")
    if y_col is not None and y_col not in df.columns:
        raise KeyError(f"列 {y_col} 不存在于 DataFrame 中")
    # 提取列并缩放
    x_array = df[x_col].to_numpy(dtype=float) * scale_x
    if y_col is not None:
        y_array = df[y_col].to_numpy(dtype=float) * scale_y
        # 按行过滤 NaN / Inf
        mask = np.isfinite(x_array) & np.isfinite(y_array)
        return x_array[mask], y_array[mask]
    else:
        mask = np.isfinite(x_array)
        return x_array[mask]
def extract_PairData(x_key, y_key, source, x_Scale = 1.0, y_Scale = 1.0):
    # 统一取出数据
    if isinstance(source, pd.DataFrame):
        x = source[x_key].values
        u = source[y_key].values
    else:  # 假定是字典形式
        x = np.array(source[x_key])
        u = np.array(source[y_key])
    # 转换单位 kPa → bar
    x_new = x * x_Scale
    y_new = u * y_Scale
    # 删除 NaN
    valid = ~np.isnan(x_new) & ~np.isnan(y_new)
    return np.column_stack((x_new[valid], y_new[valid]))
def clean_list(x_list, y_list=None, scale_x=1.0, scale_y=1.0):
    """
    清理 list 数据，去掉 NaN 和 Inf，并统一返回 (x_clean, y_clean)。
    参数:
        x_list : list 或可迭代，第一列数据（通常是压力）
        y_list : list 或可迭代，第二列数据（通常是吸附量），可选
        scale_x : float，x 列缩放系数，默认 1.0
        scale_y : float，y 列缩放系数，默认 1.0
    返回:
        x_clean : numpy array，过滤后的 x
        y_clean : numpy array，如果未提供 y_list，则返回 None
    """
    x_array = np.array(x_list, dtype=float) * scale_x
    if y_list is not None:
        y_array = np.array(y_list, dtype=float) * scale_y
        mask = np.isfinite(x_array) & np.isfinite(y_array)  # 按行过滤
        return x_array[mask], y_array[mask]
    else:
        mask = np.isfinite(x_array)
        return x_array[mask], None
def sortData(x, y=None):
    idx = np.argsort(x)
    x_sorted = x[idx]
    if y is not None:
        y_sorted = y[idx]
        return x_sorted, y_sorted
    else:
        return x_sorted, None
def naturalSortKey(value):
    """
    Generate a natural sort key.
    Examples
    --------
    S1, S2, S10
    A1, A2, A10
    Blank
    Reference
    """
    value = str(value)
    m = re.search(r"\d+", value)
    if m is None:
        # no number
        return (value.lower(), -1)
    prefix = value[:m.start()].lower()
    number = int(m.group())
    suffix = value[m.end():].lower()
    return (prefix, number, suffix)
def naturalSort(items):
    """
    Naturally sort iterable objects.
    Parameters
    ----------
    items : iterable
    Returns
    -------
    list
    """
    return sorted( list(items), key=naturalSortKey )
def naturalSortData(data, index=None, axis=0, level=0):
    """
    Sort supported data objects using natural order.
    Parameters
    ----------
    data : DataFrame, Series, Index, MultiIndex,
           list, tuple or ndarray
        Input object.
    index : iterable, optional
        External sorting labels. Only supported when `data`
        is a DataFrame. If None, the DataFrame index is used.
    axis : {0, 1}, default=0
        Axis to sort when `data` is a DataFrame.
        * axis=0 : sort row labels (index).
        * axis=1 : sort column labels.
    level : int, default=0
        MultiIndex level used when sorting labels.
    Returns
    -------
    Same type as input whenever possible.
    """
    axis = normalizeAxis(axis)
    # -----------------------------
    # DataFrame
    # -----------------------------
    if isinstance(data, pd.DataFrame):
        # Use external sorting labels
        if index is not None:
            order = natsort.index_natsorted(index)
            return data.iloc[order]
        labels = getSortedLabels(data, axis=axis, level=level)
        # sort index
        if axis == 0:
            return data.loc[labels]
        # sort columns
        if not isinstance(data.columns, pd.MultiIndex):
            return data.loc[:, labels]
        # MultiIndex columns
        columns = []
        for label in labels:
            columns.extend(
                [col for col in data.columns if col[level] == label]
            )
        return data.loc[:, columns]
    # index 参数仅支持 DataFrame
    if index is not None:
        raise TypeError(
            "'index' is only supported when data is a pandas.DataFrame."
        )
    # -----------------------------
    # Series
    # -----------------------------
    elif isinstance(data, pd.Series):
        labels = getSortedLabels(data)
        return data.loc[labels]
    # -----------------------------
    # Pandas Index
    # -----------------------------
    elif isinstance(data, pd.Index):
        return pd.Index(getSortedLabels(data))
    # -----------------------------
    # ndarray
    # -----------------------------
    elif isinstance(data, np.ndarray):
        return np.array(getSortedLabels(data))
    # -----------------------------
    # tuple
    # -----------------------------
    elif isinstance(data, tuple):
        return tuple(getSortedLabels(data))
    # -----------------------------
    # list / iterable
    # -----------------------------
    else:
        return getSortedLabels(data)
def naturalSortBy(data, by, reset_index=True):
    """
    Naturally sort a DataFrame by the values of a column.
    Parameters
    ----------
    data : pandas.DataFrame
        Input DataFrame.
    by : str
        Column name used as the natural sorting key.
    reset_index : bool, default=True
        Whether to reset the index after sorting.
    Returns
    -------
    pandas.DataFrame
        Naturally sorted DataFrame.
    """
    result = naturalSortData(data, index=data[by])
    if reset_index:
        result = result.reset_index(drop=True)
    return result
def normalizeAxis(axis):
    """
    Normalize axis specification.
    Parameters
    ----------
    axis : {0, 1, "index", "columns"}
    Returns
    -------
    int
        0 -> index
        1 -> columns
    Raises
    ------
    ValueError
        If axis is invalid.
    """
    if axis in (0, "index"):
        return 0
    if axis in (1, "columns"):
        return 1
    raise ValueError(
        "axis must be 0, 1, 'index' or 'columns'."
    )
def getLabels(data, axis=0, level=0, unique=True):
    """
    Extract labels from different data types.
    Parameters
    ----------
    data : DataFrame, Series, Index, MultiIndex,
           list, tuple or ndarray
        Input object.
    axis : {0, 1}, default=0
        Axis used when `data` is a DataFrame.
        * axis=0 : extract row labels (index).
        * axis=1 : extract column labels.
    level : int, default=0
        MultiIndex level used when extracting labels.
        Examples
        --------
        MultiIndex columns
              Sample          S1          S2
              Property     Width PSD   Width PSD
        level=0  ->  ["S1", "S2"]
        level=1  ->  ["Width", "PSD"]
    unique : bool, default=True
        Remove duplicated labels while preserving
        the original order.
    Returns
    -------
    list
        Extracted labels.
    """
    axis = normalizeAxis(axis)
    # -----------------------------
    # DataFrame
    # -----------------------------
    if isinstance(data, pd.DataFrame):
        if axis == 0:
            labels = list(data.index)
        elif axis == 1:
            if isinstance(data.columns, pd.MultiIndex):
                labels = list( data.columns.get_level_values(level) )
            else:
                labels = list(data.columns)
        else:
            raise ValueError("axis must be 0 or 1.")
    # -----------------------------
    # Series
    # -----------------------------
    elif isinstance(data, pd.Series):
        labels = list(data.index)
    # -----------------------------
    # Index
    # -----------------------------
    elif isinstance(data, pd.Index):
        labels = list(data)
    # -----------------------------
    # MultiIndex
    # -----------------------------
    elif isinstance(data, pd.MultiIndex):
        labels = list( data.get_level_values(level) )
    # -----------------------------
    # list / tuple / ndarray
    # -----------------------------
    else:
        labels = list(data)
    if unique:
        labels = list(dict.fromkeys(labels))
    return labels
def getSortedLabels(data, axis=0, level=0, unique=True):
    """
    Extract labels and sort them in natural order.
    Parameters
    ----------
    data : supported by getLabels()
    axis : {0, 1}, default=0
        Axis used for DataFrame.
    level : int, default=0
        MultiIndex level used when extracting labels.
    unique : bool, default=True
        Remove duplicated labels before sorting.
    Returns
    -------
    list
        Naturally sorted labels.
    Examples
    --------
    >>> getSortedLabels(metrics)
    ['S1', 'S2', 'S10']
    >>> getSortedLabels(metrics, axis=1)
    ['CompetitionIndex', 'HighLowRatio', 'Skewness']
    >>> getSortedLabels(psdData, axis=1)
    ['S1', 'S2', 'S10']
    >>> getSortedLabels(psdData, axis=1, level=1)
    ['PSD', 'Width']
    """
    axis = normalizeAxis(axis)
    labels = getLabels( data, axis=axis, level=level, unique=unique )
    return naturalSort(labels)
def buildRegions(splitPoints, names=None, digits=2):
    """
    Build region dictionary from split points.
    Parameters
    ----------
    splitPoints : sequence
        Boundary points.
        Example:
        [0.33, 0.50, 0.65, 0.80]
    names : list, optional
        Region names. If None, names will be generated automatically.
    digits : int
        Decimal places used in automatic names.
    Returns
    -------
    dict
    """
    if len(splitPoints) < 2:
        raise ValueError(
            "splitPoints must contain at least two values."
        )
    n = len(splitPoints) - 1
    if names is None:
        names = []
        for i in range(n):
            left = splitPoints[i]
            right = splitPoints[i + 1]
            if i == 0:
                names.append(f"<{right:.{digits}f}")
            else:
                names.append(
                    f"{left:.{digits}f}-{right:.{digits}f}"
                )
    elif len(names) != n:
        raise ValueError(
            "Number of names must equal len(splitPoints)-1."
        )
    return {
        name: (splitPoints[i], splitPoints[i + 1])
        for i, name in enumerate(names)
    }
def cols_to_clean_df(cols, x_col, y_col):
    """
    将 dict of list 转为 DataFrame，自动去掉 NaN/Inf
    """
    x, y = np.array(cols[x_col], dtype=float), np.array(cols[y_col], dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    df = pd.DataFrame({x_col: x[mask], y_col: y[mask]})
    return df
def reshape_wide_to_long(data):
    """
    将宽表格式（如 "5Absolute Pressure (kPa)", "5uptake (mmol/g)")
    转换成长表格式，支持：
      - DataFrame 格式（各列等长）
      - dict/list 格式（各列可不同长度）
    输出列：["T", "P", "q", "P_unit", "q_unit", "T_unit"]
    """
    long_data = []
    P_unit, q_unit, T_unit = None, None, "°C"  # 默认温度单位
    # 如果传入的是 DataFrame，则取其 columns，否则取 dict 的键
    cols = data.columns if hasattr(data, "columns") else list(data.keys())
    for col in cols:
        if "pressure" in col.lower():
            # 1️⃣ 提取温度
            T_match = re.match(r"(\d+\.?\d*)", col)
            if not T_match:
                continue
            T = float(T_match.group(1))
            # 2️⃣ 提取压力单位
            unit_match = re.search(r"\((.*?)\)", col)
            if unit_match:
                P_unit = unit_match.group(1).strip()
            # 3️⃣ 找对应 uptake 列
            uptake_col = None
            for c2 in cols:
                if c2.startswith(str(int(T))) and "uptake" in c2.lower():
                    uptake_col = c2
                    unit_match2 = re.search(r"\((.*?)\)", c2)
                    if unit_match2:
                        q_unit = unit_match2.group(1).strip()
                    break
            if uptake_col is None:
                continue
            # 4️⃣ 获取列数据（兼容 list / Series）
            P_vals0 = data[col] if not hasattr(data, "columns") else data[col].tolist()
            q_vals0 = data[uptake_col] if not hasattr(data, "columns") else data[uptake_col].tolist()
            # 转为 NumPy 并去除 NaN
            P_vals = np.array(P_vals0, dtype=float)
            q_vals = np.array(q_vals0, dtype=float)
            n = min(len(P_vals), len(q_vals))  # 对齐到最短长度
            mask = np.isfinite(P_vals[:n]) & np.isfinite(q_vals[:n])
            # 构建子表
            sub = pd.DataFrame({
                "T": [T] * np.sum(mask),
                "P": P_vals[:n][mask],
                "q": q_vals[:n][mask],
                "T_unit": [T_unit] * np.sum(mask),
                "P_unit": [P_unit] * np.sum(mask),
                "q_unit": [q_unit] * np.sum(mask) 
            })
            long_data.append(sub)
    if not long_data:
        raise ValueError("❌ 未识别到压力/吸附量列，请检查列名格式（例如 '25Absolute Pressure (kPa)')")
    df_long = pd.concat(long_data, ignore_index=True)
    return df_long
def kelvin_to_celsius(T, unit="K"):
    """
    将温度转为摄氏度 (°C)
    - T: 数值或数组
    - unit: "K" 表示输入为开尔文, "C" 表示输入已是摄氏度
    """
    T = np.array(T, dtype=float)
    if unit.upper() == "K":
        return T - 273.15
    elif unit.upper() == "C":
        return T
    else:
        raise ValueError("unit must be 'C' or 'K'")
def celsius_to_kelvin(T, unit="C"):
    """
    将温度转为开尔文 (K)
    - T: 数值或数组
    - unit: "C" 表示输入为摄氏度, "K" 表示输入已是开尔文
    """
    T = np.array(T, dtype=float)  # 保证是 numpy 数组，方便处理
    if unit.upper() == "C":
        return T + 273.15
    elif unit.upper() == "K":
        return T
    else:
        raise ValueError("unit must be 'C' or 'K'")
def convert_units(rawData):
    """
    输入: rawData (DataFrame), 包含列 [P, q, T, P_unit, q_unit, T_unit]
    输出: 转换为基础单位的 DataFrame
        - P_Pa [Pa]
        - T_K [K]
        - q_mmol_g [mmol/g]
        并保留原始数据列
    """
    df = rawData.copy()
    # ---- 压力单位转换 ----
    def convert_pressure(p, unit):
        unit = unit.lower()
        if unit == "pa":
            return p
        elif "kpa" in unit:
            return p * 1e3
        elif "bar" in unit:
            return p * 1e5
        elif "atm" in unit:
            return p * 101325
        elif "torr" in unit or "mmhg" in unit:
            return p * 133.322
        else:
            raise ValueError(f"未知的压力单位: {unit}")
    df["P_Pa"] = [convert_pressure(p, u) for p, u in zip(df["P"], df["P_unit"])]
    # ---- 温度单位转换 ----
    def convert_temperature(T, unit):
        unit = unit.lower()
        if unit == "k":
            return T
        elif "c" in unit:
            return T + 273.15
        else:
            raise ValueError(f"未知的温度单位: {unit}")
    df["T_K"] = [convert_temperature(T, u) for T, u in zip(df["T"], df["T_unit"])]
    # ---- 吸附量单位转换 ----
    def convert_uptake(q, unit):
        unit = unit.lower()
        if "mmol/g" in unit:
            return q
        elif "mol/kg" in unit:
            return q * 1000  # mol/kg = mmol/g
        else:
            raise ValueError(f"未知的吸附量单位: {unit}")
    df["q_mmol_g"] = [convert_uptake(q, u) for q, u in zip(df["q"], df["q_unit"])]
    # 只保留需要的三列
    df_new = df[["q_mmol_g", "T_K", "P_Pa", "q_unit", "T_unit", "P_unit"]].copy()
    # 重命名列名并加上单位
   # 覆盖单位信息，保证后续都是统一的
    df_new["P_unit"] = "Pa"
    df_new["T_unit"] = "K"
    df_new["q_unit"] = "mmol/g"
    return df_new
# =========================================================
# 1. Shared utility: validate and clean a 1D numeric array
# =========================================================
def _validateNumericArray( data, *, name: str = "data", nan_policy: str = "raise", ) -> tuple[np.ndarray, int]:
    """Validate a one-dimensional numeric array.
    Returns
    -------
    array : np.ndarray
        Valid finite values.
    n_excluded : int
        Number of non-finite values excluded.
    """
    if nan_policy not in {"raise", "omit"}:
        raise ValueError("nan_policy must be 'raise' or 'omit'.")
    if data is None:
        raise ValueError(f"{name} is required.")
    try:
        array = np.asarray(data, dtype=float)
    except (TypeError, ValueError) as exc:
        raise TypeError( f"{name} must contain numeric values." ) from exc
    if array.ndim != 1:
        raise ValueError( f"{name} must be a one-dimensional array." )
    if array.size == 0:
        raise ValueError(f"{name} must not be empty.")
    finite_mask = np.isfinite(array)
    n_excluded = int(np.sum(~finite_mask))
    if n_excluded > 0:
        if nan_policy == "raise":
            raise ValueError( f"{name} contains {n_excluded} non-finite value(s)." )
        array = array[finite_mask]
    if array.size == 0:
        raise ValueError( f"{name} contains no finite values after filtering." )
    return array, n_excluded
# =========================================================
# 2. Calculate Student's t confidence interval for a mean
# =========================================================
def calculateTConfidenceInterval( data, *, confidence_level: float = 0.95, nan_policy: str = "raise", ) -> Dict[str, Any]:
    """Calculate a one-sample Student's t CI for the population mean."""
    _validateConfidenceLevel(confidence_level)
    x, n_excluded = _validateNumericArray( data, name="data", nan_policy=nan_policy, )
    n = x.size
    if n < 2:
        raise ValueError( "At least two observations are required for a t interval." )
    estimate = float(np.mean(x))
    standard_error = float(stats.sem(x, ddof=1))
    degrees_of_freedom = n - 1
    alpha = 1.0 - confidence_level
    critical_value = float( stats.t.ppf( 1.0 - alpha / 2.0, df=degrees_of_freedom, ) )
    margin = critical_value * standard_error
    return {
        "method": "t",
        "confidence_level": float(confidence_level),
        "estimate": estimate,
        "lower_bound": estimate - margin,
        "upper_bound": estimate + margin,
        "interval_width": 2.0 * margin,
        "standard_error": standard_error,
        "sample_size": int(n),
        "degrees_of_freedom": int(degrees_of_freedom),
        "diagnostics": {
            "status": "success",
            "warnings": [],
            "n_excluded": n_excluded,
        },
    }
# =========================================================
# 3. Calculate normal-approximation confidence interval
# =========================================================
def calculateNormalConfidenceInterval( estimate: float, standard_error: float, *, confidence_level: float = 0.95, ) -> Dict[str, Any]:
    """Calculate a normal-approximation CI from estimate and SE."""
    _validateConfidenceLevel(confidence_level)
    estimate = _validateFiniteScalar(estimate, "estimate")
    standard_error = _validateFiniteScalar( standard_error, "standard_error", )
    if standard_error < 0:
        raise ValueError("standard_error must be non-negative.")
    alpha = 1.0 - confidence_level
    critical_value = float( stats.norm.ppf(1.0 - alpha / 2.0) )
    margin = critical_value * standard_error
    warnings = []
    if standard_error == 0:
        warnings.append( "The standard error is zero; the interval collapses to the point estimate." )
    return {
        "method": "normal",
        "confidence_level": float(confidence_level),
        "estimate": estimate,
        "lower_bound": estimate - margin,
        "upper_bound": estimate + margin,
        "interval_width": 2.0 * margin,
        "standard_error": standard_error,
        "sample_size": None,
        "degrees_of_freedom": None,
        "diagnostics": {
            "status": "success",
            "warnings": warnings,
            "n_excluded": 0,
        },
    }
# =========================================================
# 4. Calculate percentile confidence interval
# =========================================================
def calculatePercentileConfidenceInterval(
    sampling_distribution,
    *,
    confidence_level: float = 0.95,
    estimate: Optional[float] = None,
    nan_policy: str = "raise",
) -> Dict[str, Any]:
    """Calculate a percentile CI from an existing sampling distribution.
    This function does not generate bootstrap samples.
    """
    _validateConfidenceLevel(confidence_level)
    distribution, n_excluded = _validateNumericArray( sampling_distribution, name="sampling_distribution", nan_policy=nan_policy, )
    n = distribution.size
    if n < 2:
        raise ValueError( "At least two finite estimates are required for a percentile interval." )
    if estimate is None:
        point_estimate = float(np.mean(distribution))
    else:
        point_estimate = _validateFiniteScalar( estimate, "estimate", )
    alpha = 1.0 - confidence_level
    lower, upper = np.quantile( distribution, [alpha / 2.0, 1.0 - alpha / 2.0], )
    return {
        "method": "percentile",
        "confidence_level": float(confidence_level),
        "estimate": point_estimate,
        "lower_bound": float(lower),
        "upper_bound": float(upper),
        "interval_width": float(upper - lower),
        "standard_error": (
            float(np.std(distribution, ddof=1))
            if n >= 2
            else None
        ),
        "sample_size": int(n),
        "degrees_of_freedom": None,
        "diagnostics": {
            "status": "success",
            "warnings": [ "The percentile interval does not apply bias correction or acceleration correction." ],
            "n_excluded": n_excluded,
        },
    }
# =========================================================
# 5. Shared scalar and confidence-level validation
# =========================================================
def _validateConfidenceLevel(confidence_level: float) -> None:
    """Validate the confidence level."""
    try:
        confidence_level = float(confidence_level)
    except (TypeError, ValueError) as exc:
        raise TypeError( "confidence_level must be a numeric scalar." ) from exc
    if not np.isfinite(confidence_level) or not ( 0.0 < confidence_level < 1.0 ):
        raise ValueError( "confidence_level must be strictly between 0 and 1." )
def _validateFiniteScalar(value, name: str) -> float:
    """Validate a finite numeric scalar."""
    if not np.isscalar(value):
        raise TypeError(f"{name} must be a scalar.")
    try:
        value = float(value)
    except (TypeError, ValueError) as exc:
        raise TypeError( f"{name} must be numeric." ) from exc
    if not np.isfinite(value):
        raise ValueError(f"{name} must be finite.")
    return value
# =========================================================
# 6. Unified public entry point
# =========================================================
def runConfidenceIntervalAnalysis(
    data=None,
    *,
    method: str = "t",
    confidence_level: float = 0.95,
    estimate: Optional[float] = None,
    standard_error: Optional[float] = None,
    bootstrap_distribution=None,
    nan_policy: str = "raise",
) -> Dict[str, Any]:
    """Unified entry point for confidence interval calculations.
    Supported methods
    -----------------
    t:
        data -> one-sample t interval for the mean.
    normal:
        estimate + standard_error -> normal-approximation interval.
    percentile:
        bootstrap_distribution -> percentile interval.
    """
    _validateConfidenceLevel(confidence_level)
    if method not in {"t", "normal", "percentile"}:
        raise ValueError( "method must be 't', 'normal', or 'percentile'." )
    if nan_policy not in {"raise", "omit"}:
        raise ValueError("nan_policy must be 'raise' or 'omit'.")
    if method == "t":
        if ( estimate is not None
            or standard_error is not None
            or bootstrap_distribution is not None ):
            raise ValueError( "method='t' accepts data and confidence_level only." )
        return calculateTConfidenceInterval( data, confidence_level=confidence_level, nan_policy=nan_policy, )
    if method == "normal":
        if ( data is not None
            or bootstrap_distribution is not None ):
            raise ValueError( "method='normal' accepts estimate, standard_error, and confidence_level only." )
        if estimate is None or standard_error is None:
            raise ValueError( "Both estimate and standard_error are required for method='normal'." )
        return calculateNormalConfidenceInterval( estimate, standard_error, confidence_level=confidence_level, )
    # method == "percentile"
    if data is not None or standard_error is not None:
        raise ValueError( "method='percentile' accepts bootstrap_distribution, optional estimate, and confidence_level only." )
    if bootstrap_distribution is None:
        raise ValueError( "bootstrap_distribution is required for method='percentile'." )
    return calculatePercentileConfidenceInterval(
        bootstrap_distribution,
        confidence_level=confidence_level,
        estimate=estimate,
        nan_policy=nan_policy,
    )
def _validatePairedBootstrapData( data, nan_policy: str = "raise", ) -> np.ndarray:
    """Validate paired observations stored in rows."""
    if nan_policy not in {"raise", "omit"}:
        raise ValueError(
            "nan_policy must be 'raise' or 'omit'."
        )
    try:
        data = np.asarray(data, dtype=float)
    except (TypeError, ValueError) as exc:
        raise TypeError(
            "data must contain numeric values."
        ) from exc
    if data.ndim != 2:
        raise ValueError(
            "data must be a two-dimensional array."
        )
    if data.shape[0] < 2:
        raise ValueError(
            "At least two observations are required."
        )
    if data.shape[1] < 2:
        raise ValueError(
            "data must contain at least two variables."
        )
    valid_rows = np.isfinite(data).all(axis=1)
    if not valid_rows.all():
        if nan_policy == "raise":
            raise ValueError(
                "data contains NaN or infinite values."
            )
        data = data[valid_rows]
    if data.shape[0] < 2:
        raise ValueError(
            "At least two complete observations are required."
        )
    return data
def _calculatePairedBootstrapDistributions(
    data: np.ndarray,
    indices: np.ndarray,
    statistic_func: Callable,
) -> Dict[str, np.ndarray]:
    """
    Calculate multiple scalar statistics for paired bootstrap samples.
    statistic_func must return a dictionary of scalar statistics.
    Non-finite bootstrap estimates are recorded as NaN.
    """
    n_bootstrap = indices.shape[0]
    distributions = None
    expected_keys = None
    for i, sample_indices in enumerate(indices):
        bootstrap_sample = data[sample_indices]
        try:
            result = statistic_func(bootstrap_sample)
        except Exception as exc:
            raise ValueError(
                f"statistic_func failed during "
                f"bootstrap iteration {i + 1}."
            ) from exc
        if not isinstance(result, dict) or not result:
            raise TypeError(
                "statistic_func must return a non-empty dictionary."
            )
        if expected_keys is None:
            expected_keys = tuple(result.keys())
            distributions = {
                key: np.full(n_bootstrap, np.nan, dtype=float)
                for key in expected_keys
            }
        elif set(result.keys()) != set(expected_keys):
            raise ValueError(
                f"Inconsistent statistic keys during "
                f"bootstrap iteration {i + 1}."
            )
        for key in expected_keys:
            try:
                value = float(result[key])
            except (TypeError, ValueError) as exc:
                raise TypeError(
                    f"Statistic '{key}' must be scalar and numeric."
                ) from exc
            if np.isfinite(value):
                distributions[key][i] = value
    return distributions
def calculatePairedCorrelationStatistics(
    data,
    metric_names,
    target_names,
    statistic_names=("Pearson_r", "Spearman_r"),
):
    """
    Calculate selected correlation statistics between two groups
    of variables.
    Parameters
    ----------
    data : array-like, shape (n_samples, n_variables)
        Data matrix. Each row represents one sample.
    metric_names : dict
        Variable names and column indices for the first group.
        Example: {"Centroid": 0, "Skewness": 1}
    target_names : dict
        Variable names and column indices for the second group.
        Example: {"E": 3}
    statistic_names : sequence of str
        Statistics to retain from calculateCorrelation().
        Available:
        "Pearson_r", "Pearson_p", "Spearman_r", "Spearman_p".
    Returns
    -------
    dict
        Flat dictionary containing the selected statistics
        for each variable pair.
    """
    data = np.asarray(data, dtype=float)
    if data.ndim != 2:
        raise ValueError("data must be a two-dimensional array.")
    available_statistics = {
        "Pearson_r",
        "Pearson_p",
        "Spearman_r",
        "Spearman_p",
    }
    statistic_names = tuple(statistic_names)
    if not statistic_names:
        raise ValueError("statistic_names must not be empty.")
    invalid_statistics = (
        set(statistic_names) - available_statistics
    )
    if invalid_statistics:
        raise ValueError(
            f"Unsupported statistics: {sorted(invalid_statistics)}"
        )
    result = {}
    for group in (metric_names, target_names):
        for name, index in group.items():
            if (
                isinstance(index, (bool, np.bool_))
                or not isinstance(index, (int, np.integer))
                or not 0 <= index < data.shape[1]
            ):
                raise ValueError(
                    f"Invalid column index for '{name}': {index}"
                )
    for x_name, x_idx in metric_names.items():
        for y_name, y_idx in target_names.items():
            correlation = mf.calculateCorrelation(
                x=data[:, x_idx],
                y=data[:, y_idx],
                x_name=x_name,
                y_name=y_name,
                dropna=False,
            )
            pair_name = f"{x_name}_vs_{y_name}"
            for statistic_name in statistic_names:
                result[
                    f"{pair_name}_{statistic_name}"
                ] = correlation[statistic_name]
    return result
def runPairedBootstrapAnalysis(
    data,
    statistic_func: Callable,
    *,
    n_bootstrap: int = 10000,
    confidence_level: float = 0.95,
    random_state: Optional[int] = None,
    nan_policy: str = "raise",
) -> Dict[str, Any]:
    """
    Run paired bootstrap analysis for multiple scalar statistics.
    Each row represents one observation. All columns are resampled
    together to preserve the pairing between variables.
    statistic_func receives a two-dimensional sample and must return
    a dictionary of scalar statistics.
    one excample of Callable is def calculatePairedCorrelationStatistics()
    """
    if not callable(statistic_func):
        raise TypeError(
            "statistic_func must be callable."
        )
    if (
        isinstance(n_bootstrap, (bool, np.bool_))
        or not isinstance(n_bootstrap, (int, np.integer))
        or n_bootstrap < 2
    ):
        raise ValueError(
            "n_bootstrap must be an integer >= 2."
        )
    if (
        isinstance(confidence_level, (bool, np.bool_))
        or not isinstance(
            confidence_level,
            (int, float, np.integer, np.floating),
        )
        or not np.isfinite(confidence_level)
        or not 0.0 < confidence_level < 1.0
    ):
        raise ValueError(
            "confidence_level must be between 0 and 1."
        )
    # 1. Validate data
    data = _validatePairedBootstrapData(
        data=data,
        nan_policy=nan_policy,
    )
    # 2. Calculate statistics from the original data
    try:
        original_result = statistic_func(data)
    except Exception as exc:
        raise ValueError(
            "statistic_func failed on the original data."
        ) from exc
    if not isinstance(original_result, dict) or not original_result:
        raise TypeError(
            "statistic_func must return a non-empty dictionary."
        )
    estimates = {}
    for key, value in original_result.items():
        try:
            value = float(value)
        except (TypeError, ValueError) as exc:
            raise TypeError(
                f"Statistic '{key}' must be scalar and numeric."
            ) from exc
        if not np.isfinite(value):
            raise ValueError(
                f"Statistic '{key}' is non-finite on original data."
            )
        estimates[key] = value
    # 3. Generate one shared set of bootstrap indices
    indices = _generateBootstrapIndices(
        n_samples=data.shape[0],
        n_bootstrap=n_bootstrap,
        random_state=random_state,
    )
    # 4. Calculate all statistics using the same resamples
    distributions = _calculatePairedBootstrapDistributions(
        data=data,
        indices=indices,
        statistic_func=statistic_func,
    )
    if set(distributions) != set(estimates):
        raise ValueError(
            "Statistics returned for bootstrap samples do not match "
            "those returned for the original data."
        )
    # 5. Summarize each statistic and calculate its confidence interval
    statistics = {}
    for key, distribution in distributions.items():
        valid_distribution = distribution[
            np.isfinite(distribution)
        ]
        n_valid = len(valid_distribution)
        n_invalid = n_bootstrap - n_valid
        if n_valid < 2:
            raise ValueError(
                f"Statistic '{key}' has fewer than two valid "
                "bootstrap estimates; its uncertainty cannot "
                "be estimated reliably."
            )
        ci_result = runConfidenceIntervalAnalysis(
            method="percentile",
            confidence_level=confidence_level,
            estimate=estimates[key],
            bootstrap_distribution=valid_distribution,
        )
        # print("ci_result:", ci_result)
        statistics[key] = {
            "estimate": estimates[key],
            "bootstrap_distribution": distribution,
            "bootstrap_standard_error": float(
                np.std(valid_distribution, ddof=1)
            ),
            "confidence_interval": (
                ci_result["lower_bound"],
                ci_result["upper_bound"],
            ),
            "confidence_level": confidence_level,
            "n_valid": n_valid,
            "n_invalid": n_invalid,
            "invalid_fraction": n_invalid / n_bootstrap,
        }
    return {
        "statistics": statistics,
        "n_samples": int(data.shape[0]),
        "n_variables": int(data.shape[1]),
        "n_bootstrap": int(n_bootstrap),
        "confidence_level": confidence_level,
        "random_state": random_state,
        "method": "paired_bootstrap_percentile",
    }

def _validateBootstrapData(
    data,
    nan_policy: str = "raise",
) -> np.ndarray:
    """Validate and prepare one-dimensional bootstrap data."""
    if nan_policy not in {"raise", "omit"}:
        raise ValueError("nan_policy must be 'raise' or 'omit'.")
    try:
        x = np.asarray(data, dtype=float)
    except (TypeError, ValueError) as exc:
        raise TypeError(
            "data must contain numeric values."
        ) from exc
    if x.ndim != 1:
        raise ValueError(
            "data must be a one-dimensional array."
        )
    if x.size == 0:
        raise ValueError("data must not be empty.")
    finite_mask = np.isfinite(x)
    if not finite_mask.all():
        if nan_policy == "raise":
            raise ValueError(
                "data contains NaN or infinite values."
            )
        x = x[finite_mask]
    if x.size < 2:
        raise ValueError(
            "At least two finite observations are required."
        )
    return x
def _generateBootstrapIndices(
    n_samples: int,
    n_bootstrap: int,
    random_state: Optional[int] = None,
) -> np.ndarray:
    """
    Generate bootstrap indices by sampling with replacement.
    Returns
    -------
    np.ndarray
        Array with shape (n_bootstrap, n_samples).
    """
    rng = np.random.default_rng(random_state)
    return rng.integers(
        low=0,
        high=n_samples,
        size=(n_bootstrap, n_samples),
    )
def _calculateBootstrapDistribution(
    data: np.ndarray,
    indices: np.ndarray,
    statistic_func: Callable,
) -> np.ndarray:
    """
    Calculate a scalar statistic for each bootstrap resample.
    Parameters
    ----------
    data : np.ndarray
        Original data.
    indices : np.ndarray
        Bootstrap indices with shape (n_bootstrap, n_samples).
    statistic_func : callable
        Function that accepts one resampled dataset and returns
        a scalar statistic.
    Returns
    -------
    np.ndarray
        Bootstrap estimates with shape (n_bootstrap,).
    """
    bootstrap_estimates = np.empty(
        indices.shape[0],
        dtype=float,
    )
    for i, sample_indices in enumerate(indices):
        bootstrap_sample = data[sample_indices]
        try:
            estimate = float(statistic_func(bootstrap_sample))
        except Exception as exc:
            raise ValueError(
                f"statistic_func failed during bootstrap iteration {i + 1}."
            ) from exc
        if not np.isfinite(estimate):
            raise ValueError(
                f"statistic_func returned a non-finite estimate "
                f"during bootstrap iteration {i + 1}."
            )
        bootstrap_estimates[i] = estimate
    return bootstrap_estimates
def runBootstrapAnalysis(
    data,
    statistic_func: Callable,
    *,
    n_bootstrap: int = 10000,
    confidence_level: float = 0.95,
    random_state: Optional[int] = None,
    nan_policy: str = "raise",
) -> Dict[str, Any]:
    """
    Run bootstrap analysis for a scalar statistic.
    Parameters
    ----------
    data : array-like
        One-dimensional observations.
    statistic_func : callable
        Function that accepts a bootstrap sample and returns
        one scalar statistic.
    n_bootstrap : int, default=10000
        Number of bootstrap replicates.
    confidence_level : float, default=0.95
        Confidence level.
    random_state : int or None, default=None
        Random seed.
    nan_policy : {"raise", "omit"}, default="raise"
        Handling of NaN and infinite input values.
    Returns
    -------
    dict
        Original estimate, bootstrap distribution, standard
        error, confidence interval, and analysis metadata.
    """
    if not callable(statistic_func):
        raise TypeError( "statistic_func must be callable." )
    if (
        isinstance(n_bootstrap, (bool, np.bool_))
        or not isinstance(n_bootstrap, (int, np.integer))
        or n_bootstrap < 1
    ):
        raise ValueError( "n_bootstrap must be a positive integer." )
    if (
        isinstance(confidence_level, (bool, np.bool_))
        or not isinstance( confidence_level, (int, float, np.integer, np.floating), )
        or not np.isfinite(confidence_level)
        or not 0.0 < confidence_level < 1.0
    ):
        raise ValueError( "confidence_level must be between 0 and 1." )
    # 1. Validate and prepare the input data
    x = _validateBootstrapData( data=data, nan_policy=nan_policy, )
    # 2. Calculate the original-data estimate
    try:
        estimate = float(statistic_func(x))
    except Exception as exc:
        raise ValueError( "statistic_func failed on the original data." ) from exc
    if not np.isfinite(estimate):
        raise ValueError( "statistic_func returned a non-finite estimate for the original data." )
    # 3. Generate the bootstrap distribution
    indices = _generateBootstrapIndices( n_samples=x.size, n_bootstrap=n_bootstrap, random_state=random_state, )
    bootstrap_estimates = _calculateBootstrapDistribution( data=x, indices=indices, statistic_func=statistic_func, )
    # 4. Calculate the bootstrap standard error
    bootstrap_standard_error = float( np.std(bootstrap_estimates, ddof=1) )
    # 5. Calculate the percentile confidence interval
    ci_result = runConfidenceIntervalAnalysis( method="percentile", confidence_level=confidence_level, 
                                              estimate=estimate, bootstrap_distribution=bootstrap_estimates, )
    # 6. Assemble the results
    return {
        "estimate": estimate,
        "bootstrap_distribution": bootstrap_estimates,
        "bootstrap_standard_error": bootstrap_standard_error,
        "confidence_interval": ci_result["confidence_interval"],
        "confidence_level": confidence_level,
        "n_samples": int(x.size),
        "n_bootstrap": int(n_bootstrap),
        "method": "bootstrap_percentile",
        "random_state": random_state,
    }
def correlationAnalysis(
        x,
        y,
        x_name="x",
        y_name="y",
        dropna=True,
        fit_method="linear",
        fit_func=None,
        p0=None,
        bounds=(-np.inf, np.inf),
        param_names=None,
        weights=None):
    # ======================================================
    # Correlation
    # ======================================================
    correlation = mf.calculateCorrelation(
        x=x,
        y=y,
        x_name=x_name,
        y_name=y_name,
        dropna=dropna
    )
    # ======================================================
    # Fitting
    # ======================================================
    if fit_method == "linear":
        fitting = mf.fitLinear(
            x=x,
            y=y,
            x_name=x_name,
            y_name=y_name,
            dropna=dropna
        )
    elif fit_method == "nonlinear":
        if fit_func is None:
            raise ValueError(
                "fit_func must be provided "
                "when fit_method='nonlinear'."
            )
        fitting = mf.fitNonlinear(
            x=x,
            y=y,
            func=fit_func,
            p0=p0,
            bounds=bounds,
            param_names=param_names,
            x_name=x_name,
            y_name=y_name,
            dropna=dropna,
            weights=weights
        )
    elif fit_method is None:
        fitting = {}
    else:
        raise ValueError( f"Unknown fit_method: {fit_method}" )
    # ======================================================
    # Combine
    # ======================================================
    result = { **correlation, **fitting }
    return result
def correlationAnalysis0( x, y, x_name="x", y_name="y", dropna=True):
    #最初的版本，只满足线性拟合，主要用于计算PSD的三指标
    """
    Perform correlation and simple linear regression analysis.
    Parameters
    ----------
    x : array-like
        Independent variable.
    y : array-like
        Dependent variable.
    x_name : str
        Name of x.
    y_name : str
        Name of y.
    dropna : bool
        Remove NaN pairs.
    Returns
    -------
    dict
        Analysis results.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    # ----------------------------------------------------
    # Remove NaN
    # ----------------------------------------------------
    if dropna:
        mask = np.isfinite(x) & np.isfinite(y)
        x = x[mask]
        y = y[mask]
    if len(x) < 2:
        raise ValueError("Not enough valid data.")
    # ----------------------------------------------------
    # Linear regression
    # ----------------------------------------------------
    reg = linregress(x, y)
    y_fit = reg.slope * x + reg.intercept
    residual = y - y_fit
    rmse = np.sqrt(np.mean(residual ** 2))
    # ================================
    # NEW: scale-aware metrics
    # ================================
    y_mean = np.mean(y)
    y_range = np.max(y) - np.min(y)
    rmse_rel = rmse / (np.abs(y_mean) + 1e-12)
    rmse_nrm = rmse / (y_range + 1e-12)
    mae = np.mean(np.abs(y - y_fit))
    mape = np.mean(np.abs((y - y_fit) / y)) * 100
    # ----------------------------------------------------
    # Spearman
    # ----------------------------------------------------
    spearman_r, spearman_p = spearmanr(x, y)
    # ----------------------------------------------------
    # Output
    # ----------------------------------------------------
    result = {
        "x_name": x_name,
        "y_name": y_name,
        "N": len(x),
        "X": x,
        "Y": y,
        "Y_fit": y_fit,
        "Residual": residual,
        "Slope": reg.slope,
        "Intercept": reg.intercept,
        "Pearson_r": reg.rvalue,
        "Pearson_p": reg.pvalue,
        "R2": reg.rvalue ** 2,
        "StdErr": reg.stderr,
        # ------------------------
        # raw + normalized RMSE
        # ------------------------
        "RMSE": rmse,
        "RMSE_rel": rmse_rel,
        "RMSE_nrm": rmse_nrm,
        "MAE": mae,
        "MAPE": mape,
        "Spearman_r": spearman_r,
        "Spearman_p": spearman_p,
    }
    return result
def batchCorrelationAnalysis( xData, yData):
    """
    Batch correlation analysis.
    Parameters
    ----------
    xData : DataFrame
    yData : Series or DataFrame
    Returns
    -------
    summary : DataFrame
    results : dict
    """
    if isinstance(yData, pd.Series):
        yData = yData.to_frame()
    summary = []
    results = {}
    for x_name in xData.columns:
        results[x_name] = {}
        for y_name in yData.columns:
            result = correlationAnalysis( xData[x_name], yData[y_name], x_name, y_name )
            results[x_name][y_name] = result
            summary.append({
                "X": x_name,
                "Y": y_name,
                "Pearson_r": result["Pearson_r"],
                "Pearson_p": result["Pearson_p"],
                "Abs_Pearson_r": abs(result["Pearson_r"]),
                "Spearman_r": result["Spearman_r"],
                "Spearman_p": result["Spearman_p"],
                "R2": result["R2"],
                "RMSE": result["RMSE"],              # raw
                "RMSE_rel": result["RMSE_rel"],      # comparable
                "RMSE_nrm": result["RMSE_nrm"],      # robust scale-free
            })
    summary = pd.DataFrame(summary)
    return summary, results
def crossCorrelationAnalysis(dataX, dataY, xColumns=None, yColumns=None):
    """
    Prepare X/Y datasets and perform batch correlation analysis.
    Parameters
    ----------
    dataX : pandas.DataFrame
        DataFrame containing X variables.
    dataY : pandas.DataFrame
        DataFrame containing Y variables.
    xColumns : str or list[str], optional
        Columns selected from dataX.
        If None, all columns in dataX are used.
    yColumns : str or list[str], optional
        Columns selected from dataY.
        If None, all columns in dataY are used.
    Returns
    -------
    X : pandas.DataFrame
        Selected X variables.
    Y : pandas.DataFrame
        Selected Y variables.
    summary : pandas.DataFrame
        Summary of correlation analysis.
    results : dict
        Detailed correlation analysis results.
    """
    # -----------------------------
    # Convert column names to list
    # -----------------------------
    if xColumns is None:
        xColumns = list(dataX.columns)
    elif isinstance(xColumns, str):
        xColumns = [xColumns]
    if yColumns is None:
        yColumns = list(dataY.columns)
    elif isinstance(yColumns, str):
        yColumns = [yColumns]
    # -----------------------------
    # Sort samples
    # -----------------------------
    dataX = naturalSortData(dataX)
    dataY = naturalSortData(dataY)
    # -----------------------------
    # Select variables
    # -----------------------------
    X = dataX[xColumns]
    Y = dataY[yColumns]
    # -----------------------------
    # Batch correlation analysis
    # -----------------------------
    summary, results = batchCorrelationAnalysis(X, Y)
    return X, Y, summary, results
def matrixCorrelationAnalysis(dataX, dataY=None, columns=None):
    """
    Prepare data for correlation matrix analysis.
    Parameters
    ----------
    dataX : pandas.DataFrame
        First dataset.
    dataY : pandas.DataFrame
        Second dataset.
    columns : str or list[str], optional
        Variables to include.
        If None, all variables from both datasets are used.
    Returns
    -------
    X : pandas.DataFrame
        Variables used in correlation analysis.
    Y : pandas.DataFrame
        Same as X.
    summary : pandas.DataFrame
        Summary of correlation analysis.
    results : dict
        Detailed correlation analysis results.
    """
    # -----------------------------
    # Natural sort samples
    # -----------------------------
    dataX = naturalSortData(dataX)
    if dataY is None:
        data = dataX
    else:
        dataY = naturalSortData(dataY)
        data = pd.concat([dataX, dataY], axis=1)
    # -----------------------------
    # Select variables
    # -----------------------------
    if columns is None:
        columns = list(data.columns)
    elif isinstance(columns, str):
        columns = [columns]
    # -----------------------------
    # Cross correlation (X = Y)
    # -----------------------------
    X, Y, summary, results = crossCorrelationAnalysis( data, data, xColumns=columns, yColumns=columns )
    return X, Y, summary, results