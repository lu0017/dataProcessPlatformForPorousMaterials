from common import *
import constantsAndName as const
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
def smooth(
    y,
    method="savgol",
    **kwargs,
):
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
def naturalSortData(data, axis=0, level=0):
    """
    Sort supported data objects using natural order.
    Parameters
    ----------
    data : DataFrame, Series, Index, MultiIndex,
           list, tuple or ndarray
        Input object.
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
        labels = getSortedLabels( data, axis=axis, level=level )
        # sort index
        if axis == 0:
            return data.loc[labels]
        # sort columns
        if not isinstance(data.columns, pd.MultiIndex):
            return data.loc[:, labels]
        # MultiIndex columns
        columns = []
        for label in labels:
            columns.extend( [ col for col in data.columns if col[level] == label ] )
        return data.loc[:, columns]
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
def correlationAnalysis( x, y, x_name="x", y_name="y", dropna=True):
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
def matrixCorrelationAnalysis(dataX, dataY, columns=None):
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
    dataY = naturalSortData(dataY)
    # -----------------------------
    # Merge datasets
    # -----------------------------
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