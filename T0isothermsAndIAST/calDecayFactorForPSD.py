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
import T1dataProcessSource.modelAndFit as mf
def flattenNestedResult(obj, analysis_name):
    """
    将嵌套的 PCA / VIF 结果整理为一张长表。
    输出列：
        Analysis
        Result_Path
        Row_Key
        Column_Key
        Value
    """
    rows = []
    def add_value(path, value, row_key=None, column_key="Value"):
        # 转换 NumPy 标量，便于保存到 Excel
        if isinstance(value, np.generic):
            value = value.item()
        rows.append({
            "Analysis": analysis_name,
            "Result_Path": path or "Root",
            "Row_Key": (
                str(row_key) if row_key is not None else None
            ),
            "Column_Key": str(column_key),
            "Value": value,
        })
    def walk(value, path="Root"):
        if value is None:
            add_value(path, None)
            return
        # ----------------------------------
        # DataFrame：保留原始行列标签
        # ----------------------------------
        if isinstance(value, pd.DataFrame):
            for i, row_label in enumerate(value.index):
                for j, column_label in enumerate(value.columns):
                    cell_value = value.iat[i, j]
                    if isinstance(cell_value, (dict, list, tuple, np.ndarray,
                                               pd.Series, pd.DataFrame)):
                        cell_value = repr(cell_value)
                    add_value(
                        path,
                        cell_value,
                        row_key=row_label,
                        column_key=column_label,
                    )
            return
        # ----------------------------------
        # Series：保留原始索引
        # ----------------------------------
        if isinstance(value, pd.Series):
            for row_label, cell_value in value.items():
                if isinstance(cell_value, (dict, list, tuple, np.ndarray,
                                           pd.Series, pd.DataFrame)):
                    cell_value = repr(cell_value)
                add_value(
                    path,
                    cell_value,
                    row_key=row_label,
                    column_key=value.name or "Value",
                )
            return
        # ----------------------------------
        # dict：递归展开，但不创建新的 Sheet
        # ----------------------------------
        if isinstance(value, dict):
            for key, child in value.items():
                child_path = f"{path}.{key}"
                walk(child, child_path)
            return
        # ----------------------------------
        # ndarray：保留数组索引
        # ----------------------------------
        if isinstance(value, np.ndarray):
            if value.ndim == 0:
                walk(value.item(), path)
                return
            for index in np.ndindex(value.shape):
                index_label = ",".join(map(str, index))
                child_path = f"{path}[{index_label}]"
                walk(value[index], child_path)
            return
        # ----------------------------------
        # list / tuple：保留元素位置
        # ----------------------------------
        if isinstance(value, (list, tuple)):
            for i, child in enumerate(value):
                walk(child, f"{path}[{i}]")
            return
        # ----------------------------------
        # 标量
        # ----------------------------------
        if np.isscalar(value):
            add_value(path, value)
            return
        # 对不认识的对象保留表示形式，避免静默丢失
        add_value(path, repr(value))
    walk(obj)
    return pd.DataFrame(
        rows,
        columns=[
            "Analysis",
            "Result_Path",
            "Row_Key",
            "Column_Key",
            "Value",
        ],
    )
def savePSDAnalysisResults(
    statisticAnalysis,
    ci_sensitivity_result,
    filename,
    kineticdata=None,
    metrics=None,
    correlation_export=None,
):
    """
    Save PSD statistical analysis results to Excel.
    PCA and VIF nested results are each consolidated into one sheet.
    Other export behavior remains compatible with the original function.
    """
    used_sheets = set()
    # ==================================================
    # 1. Excel export helper
    # ==================================================
    def _exportDataFrame(df, sheet_name, index=False):
        if df is None:
            return
        if not isinstance(df, pd.DataFrame):
            df = pd.DataFrame(df)
        # Excel sheet name: maximum 31 characters
        base = re.sub(r"[\[\]:*?/\\]", "_", str(sheet_name))
        base = base.strip("'") or "Sheet"
        base = base[:31]
        safe_name = base
        counter = 1
        while safe_name.lower() in used_sheets:
            suffix = f"_{counter}"
            safe_name = f"{base[:31 - len(suffix)]}{suffix}"
            counter += 1
        used_sheets.add(safe_name.lower())
        fl.export_to_excel_auto(
            df,
            filename=filename,
            sheet_name=safe_name,
            index=index,
        )
    # ==================================================
    # 2. Input data
    # ==================================================
    if kineticdata is not None:
        _exportDataFrame(
            kineticdata,
            "kineticdata",
            index=True,
        )
    if metrics is not None:
        _exportDataFrame(
            metrics,
            "metrics",
            index=True,
        )
    combined_data = statisticAnalysis.get("combined_data")
    if combined_data is not None:
        _exportDataFrame(
            combined_data,
            "analysis_data",
            index=True,
        )
    # ==================================================
    # 3. Bootstrap summary
    # ==================================================
    bootstrap_result = statisticAnalysis.get("bootstrap_result")
    if bootstrap_result is not None:
        rows = []
        for name, stats in bootstrap_result.get(
            "statistics", {}
        ).items():
            ci = stats["confidence_interval"]
            rows.append({
                "Statistic": name,
                "Estimate": stats["estimate"],
                "CI_lower": ci[0],
                "CI_upper": ci[1],
                "Bootstrap_SE": stats["bootstrap_standard_error"],
                "Valid": stats["n_valid"],
                "Invalid": stats["n_invalid"],
            })
        if rows:
            _exportDataFrame(
                pd.DataFrame(rows),
                "bootstrap_summary",
            )
    # ==================================================
    # 4. Leave-one-out summary
    # ==================================================
    loo_result = statisticAnalysis.get("loo_result")
    if loo_result is not None:
        rows = []
        for name, stats in loo_result.get("summary", {}).items():
            absolute = stats.get("absolute_changes") or {}
            relative = stats.get("relative_changes") or {}
            rows.append({
                "Statistic": name,
                "Baseline": stats.get("baseline"),
                "Valid_LOO": stats.get("valid_count"),
                "Invalid_LOO": stats.get("invalid_count"),
                "LOO_min": stats.get("min"),
                "LOO_max": stats.get("max"),
                "Max_abs_change": absolute.get("max"),
                "Mean_abs_change": absolute.get("mean"),
                "Max_relative_change": relative.get("max"),
                "Mean_relative_change": relative.get("mean"),
                "Sign_consistency": stats.get("sign_consistency"),
            })
        if rows:
            _exportDataFrame(
                pd.DataFrame(rows),
                "loo_summary",
            )
    # ==================================================
    # 5. PCA: all nested results -> one sheet
    # ==================================================
    pca_result = statisticAnalysis.get("pca_result")
    if pca_result is not None:
        pca_df = flattenNestedResult(
            pca_result,
            analysis_name="PCA",
        )
        if not pca_df.empty:
            _exportDataFrame(
                pca_df,
                "PCA_Results",
            )
    # ==================================================
    # 6. VIF: all nested results -> one sheet
    # ==================================================
    vif_result = statisticAnalysis.get("vif_result")
    if vif_result is not None:
        vif_df = flattenNestedResult(
            vif_result,
            analysis_name="VIF",
        )
        if not vif_df.empty:
            _exportDataFrame(
                vif_df,
                "VIF_Results",
            )
    # ==================================================
    # 7. Partial correlation
    # ==================================================
    partial_result = statisticAnalysis.get(
        "partial_correlation_result"
    )
    if partial_result is not None:
        if isinstance(partial_result, pd.DataFrame):
            partial_df = partial_result
        elif isinstance(partial_result, pd.Series):
            partial_df = partial_result.rename(
                "Value"
            ).reset_index()
        elif isinstance(partial_result, dict):
            partial_df = pd.DataFrame(
                list(partial_result.items()),
                columns=["Statistic", "Value"],
            )
        else:
            partial_df = pd.DataFrame(partial_result)
        _exportDataFrame(
            partial_df,
            "partial_corr",
        )
    # ==================================================
    # 8. CompetitionIndex boundary sensitivity
    # ==================================================
    if ci_sensitivity_result is not None:
        ci_summary = ci_sensitivity_result.get("summary")
        ci_details = ci_sensitivity_result.get(
            "boundary_results"
        )
        if ci_summary is not None:
            _exportDataFrame(
                ci_summary,
                "CI_sensitivity",
            )
        if ci_details is not None:
            _exportDataFrame(
                ci_details,
                "CI_boundary_data",
            )
    # ==================================================
    # 9. Existing correlation export: unchanged
    # ==================================================
    if correlation_export is not None:
        fl.exportCorrelationExcel(
            filename=filename,
            **correlation_export,
        )
    print(f"\nAll analysis results exported to: {filename}")
def organizePSDCorrelation(metrics, kineticdata, xColumns, yColumns):
    """
    Organize X/Y variables and perform batch correlation analysis.
    Parameters
    ----------
    metrics : DataFrame
    kineticdata : DataFrame
    xColumns : list[str] or str
    yColumns : list[str] or str
    Returns
    -------
    X : DataFrame
    Y : DataFrame
    summary : DataFrame
    results : dict
    """
    # -----------------------------
    # Convert to list
    # -----------------------------
    if isinstance(xColumns, str):
        xColumns = [xColumns]
    if isinstance(yColumns, str):
        yColumns = [yColumns]
    metrics = dop.naturalSortData(metrics)
    kineticdata = dop.naturalSortData(kineticdata)
    # -----------------------------
    # Select variables
    # -----------------------------
    X = metrics[xColumns]
    Y = kineticdata[yColumns]
    # -----------------------------
    # Batch correlation
    # -----------------------------
    summary, results = dop.batchCorrelationAnalysis(X, Y)
    return X, Y, summary, results
def calculateRegionVolume(width, psd, lower, upper):
    """
    Integrate PSD within a pore-size range.
    Parameters
    ----------
    width : ndarray
        Pore size (nm)
    psd : ndarray
        dV/dD (cm3 g-1 nm-1)
    lower : float
        Lower pore size (nm)
    upper : float
        Upper pore size (nm)
    Returns
    -------
    float
        Integrated pore volume (cm3 g-1)
    """
    # interpolate PSD
    newWidth = np.linspace(
        width.min(),
        width.max(),
        5000
    )
    newPSD = np.interp(
        newWidth,
        width,
        psd
    )
    mask = (newWidth  >= lower) & (newWidth  <= upper)
    if np.sum(mask) < 2:
        return np.nan
    return np.trapezoid(newPSD[mask], newWidth[mask])
def calculateRegionVolume0(width, psd, regions, cumulative=False, cumulativeByUpper=False,):
    """
    Calculate integrated PSD volume for multiple pore-size regions.
    Parameters
    ----------
    width : ndarray
        Pore size (nm)
    psd : ndarray
        dV/dD (cm3 g-1 nm-1)
    regions : dict
        {
            region_name: (lower, upper)
        }
    Returns
    -------
    dict
    """
    valid = ~(np.isnan(width) | np.isnan(psd))
    width = width[valid]
    psd = psd[valid]
    result = {}
    cumulativeVolume = 0.0
    for name, (lower, upper) in regions.items():
        volume  = calculateRegionVolume( width, psd, lower, upper )
        result[name] = volume
        if cumulative:
            cumulativeVolume += volume
            if cumulativeByUpper:
                key = f"Cumulative_<{upper:.2f}"
            else:
                key = f"Cumulative_{name}"
            result[key] = cumulativeVolume
    result["Total"] = calculateRegionVolume( width, psd, width.min(), width.max() )
    return result
def calculatePSDVolumes( psdData, splitPoints, names=None, ):
    """
    Calculate segment volumes and cumulative volumes.
    Parameters
    ----------
    psdData : DataFrame
        PSD data with MultiIndex columns.
    splitPoints : sequence
        Region boundaries.
        Example:
        [0.33, 0.50, 0.65, 0.80]
    names : list, optional
        Region names.
    Returns
    -------
    DataFrame
    """
    regions = dop.buildRegions( splitPoints, names=names )
    samples = dop.naturalSort( psdData.columns.get_level_values(0).unique() )
    metrics = {}
    for sample in samples:
        sampleData = psdData[sample]
        width = sampleData.iloc[:,0].to_numpy()
        psd = sampleData.iloc[:,1].to_numpy()
        volumes = calculateRegionVolume0(
            width,
            psd,
            regions,
            cumulative=True,
            cumulativeByUpper=True
        )
        metrics[sample] = volumes
    metrics = pd.DataFrame.from_dict(
        metrics,
        orient="index"
    )
    metrics = dop.naturalSortData(metrics)
    return metrics
def calculateHighLowRatio(regionVolumes):
    """
    Calculate the ratio between high-energy
    and low-energy pore volumes.
    """
    VH = regionVolumes["High"]
    VL = regionVolumes["Low"]
    return {
        "HighLowRatio": VH / VL if VL != 0 else np.nan,
        "HighFraction": VH / regionVolumes["Total"],
        "LowFraction": VL / regionVolumes["Total"]
    }
def calculateCompetitionIndex(regionVolumes):
    """
    Competition between
    high-energy and low-energy pores.
    """
    VH = regionVolumes["High"]
    VL = regionVolumes["Low"]
    return {
        "CompetitionIndex":
        (VH - VL) / (VH + VL)
        if (VH + VL) != 0 else np.nan
    }
def calculateCentroid(width, psd):
    """
    Calculate PSD centroid.
    """
    total = np.trapezoid(psd, width)
    centroid = ( np.trapezoid(width * psd, width) / total )
    return {
        "Centroid": centroid
    }
def calculateDistributionMoment(width, psd):
    """
    Calculate statistical moments of the pore size distribution (PSD).
    Parameters
    ----------
    width : ndarray
        Pore size (nm)
    psd : ndarray
        dV/dD (cm3 g-1 nm-1)
    Returns
    -------
    dict
        {
            "Centroid": float,
            "Variance": float,
            "Std": float,
            "Skewness": float
        }
    """
    # Remove NaN
    valid = ~(np.isnan(width) | np.isnan(psd))
    width = width[valid]
    psd = psd[valid]
    # Total integrated pore volume
    total = np.trapezoid(psd, width)
    if total == 0:
        return {
            "Centroid": np.nan,
            "Variance": np.nan,
            "Std": np.nan,
            "Skewness": np.nan
        }
    # First moment (Centroid)
    centroid = np.trapezoid(width * psd, width) / total
    # Second central moment
    variance = ( np.trapezoid(((width - centroid) ** 2) * psd, width) / total )
    std = np.sqrt(variance)
    # Third standardized moment
    if std == 0:
        skewness = np.nan
    else:
        skewness = ( np.trapezoid(((width - centroid) ** 3) * psd, width) / total ) / (std ** 3)
    return {
        "Centroid": centroid,
        "Variance": variance,
        "Std": std,
        "Skewness": skewness
    }
def calculatePSDMetrics( psdData, regions, prefix=None ):
    """
    Calculate all PSD statistical descriptors.
    Parameters
    ----------
    psdData : DataFrame
        PSD data with MultiIndex columns.
    regions : dict
        Pore size regions.
    prefix : str, optional
        Prefix added to all descriptor names
        (e.g., "GCMC", "DFT").
    Returns
    -------
    metrics : DataFrame
    """
    samples = dop.naturalSort( psdData.columns.get_level_values(0).unique() )
    metrics = {}
    for sample in samples:
        sampleData = psdData[sample]
        width = sampleData.iloc[:, 0].to_numpy()
        psd = sampleData.iloc[:, 1].to_numpy()
        regionVolumes = calculateRegionVolume0( width, psd, regions, cumulative=True )
        highLow = calculateHighLowRatio(regionVolumes)
        competition = calculateCompetitionIndex(regionVolumes)
        moments = calculateDistributionMoment( width, psd )
        metrics[sample] = {}
        metrics[sample].update(regionVolumes)
        metrics[sample].update(highLow)
        metrics[sample].update(competition)
        metrics[sample].update(moments)
    metrics = pd.DataFrame.from_dict( metrics, orient="index" )
    metrics = dop.naturalSortData(metrics)
    if prefix is not None:
        metrics = metrics.add_prefix(f"{prefix}_")
    return metrics
def selectSamples(data, samples, sampleColumn="sample", axis=None):
    """
    Select specified samples from a DataFrame.
    Supports:
    1. Sample names stored in MultiIndex columns
    2. Sample names stored in a normal column
    3. Sample names stored in DataFrame index
    Parameters
    ----------
    data : pandas.DataFrame
        Input data.
    samples : list or None
        Samples to select.
        If None, all samples are retained.
    sampleColumn : str, optional
        Name of the sample column when samples are stored
        as a normal column.
    axis : int or None, optional
        None : automatically detect data structure
        0    : sample names are stored in rows/index
        1    : sample names are stored in MultiIndex columns
    Returns
    -------
    pandas.DataFrame
        Selected data.
    """
    # ==========================================
    # Select all samples
    # ==========================================
    if samples is None:
        return data.copy()
    # ==========================================
    # Case 1: MultiIndex columns
    # ==========================================
    if isinstance(data.columns, pd.MultiIndex):
        availableSamples = (
            data.columns
            .get_level_values(0)
            .unique()
        )
        missingSamples = [
            sample
            for sample in samples
            if sample not in availableSamples
        ]
        if missingSamples:
            print(
                f"Warning: samples not found: "
                f"{missingSamples}"
            )
        validSamples = [
            sample
            for sample in samples
            if sample in availableSamples
        ]
        return data.loc[
            :,
            validSamples
        ].copy()
    # ==========================================
    # Case 2: sample stored in a normal column
    # ==========================================
    if sampleColumn in data.columns:
        return data[
            data[sampleColumn].isin(samples)
        ].copy()
    # ==========================================
    # Case 3: sample stored in index
    # ==========================================
    if data.index.name == sampleColumn:
        return data.loc[
            data.index.isin(samples)
        ].copy()
    # ==========================================
    # Case 4: automatically detect sample index
    # ==========================================
    indexValues = data.index.astype(str)
    if all(
        str(sample) in indexValues
        for sample in samples
    ):
        return data.loc[
            data.index.astype(str).isin(
                [str(sample) for sample in samples]
            )
        ].copy()
    # ==========================================
    # Cannot identify sample structure
    # ==========================================
    raise ValueError(
        f"Cannot identify sample information.\n"
        f"Sample column '{sampleColumn}' not found.\n"
        f"Index name: {data.index.name}\n"
        f"Available columns: {list(data.columns)}"
    )
def _runLOOAnalysis( data, metric_names, target_names, statistic_names=("Pearson_r", "R2"), ):
    loo_result = dop.runPairedLeaveOneOutAnalysis(
        data=data,
        statistic_func=lambda sample: (
            dop.calculatePairedCorrelationStatistics(
                data=sample,
                metric_names=metric_names,
                target_names=target_names,
                statistic_names=statistic_names, ) ),
        nan_policy="raise", )
    # Print LOO results
    for name, stats in loo_result["summary"].items():
        print(f"\n{name}")
        print(f"  Baseline estimate: {stats['baseline']:.4f}")
        print(f"  Valid LOO runs: {stats['valid_count']}")
        print(f"  Invalid LOO runs: {stats['invalid_count']}")
        print( f"  LOO estimate range: [{stats['min']:.4f}, {stats['max']:.4f}]" )
        print( f"  Maximum absolute change: {stats['absolute_changes']['max']:.4f}" )
        print( f"  Mean absolute change: {stats['absolute_changes']['mean']:.4f}" )
        if stats["relative_changes"]["max"] is not None:
            print( f"  Maximum relative change: {stats['relative_changes']['max']:.2%}" )
        if stats["sign_consistency"] is not None:
            print( f"  Sign consistency: {stats['sign_consistency']:.1%}" )
    return loo_result
def _runBootstrapAnalysis(
    data,
    metric_names,
    target_names,
    statistic_names=("Pearson_r", "R2"),
    n_bootstrap=1000,
    confidence_level=0.95,
    random_state=42,
):
    """
    Run paired bootstrap analysis for Pearson correlation
    and selected goodness-of-fit statistics.
    Returns
    -------
    dict
        Original result returned by dop.runPairedBootstrapAnalysis().
    """
    # --------------------------------------------------
    # 1. Run Bootstrap
    # --------------------------------------------------
    bootstrap_result = dop.runPairedBootstrapAnalysis(
        data=data,
        statistic_func=lambda sample: (
            dop.calculatePairedCorrelationStatistics(
                data=sample,
                metric_names=metric_names,
                target_names=target_names,
                statistic_names=statistic_names,
            )
        ),
        n_bootstrap=n_bootstrap,
        confidence_level=confidence_level,
        random_state=random_state,
    )
    # --------------------------------------------------
    # 2. Print Bootstrap results
    # --------------------------------------------------
    for name, stats in bootstrap_result["statistics"].items():
        print(f"\n{name}")
        print(f"  Estimate: {stats['estimate']:.4f}")
        print(
            f"  {confidence_level:.0%} CI: "
            f"{stats['confidence_interval']}"
        )
        print(
            f"  Bootstrap SE: "
            f"{stats['bootstrap_standard_error']:.4f}"
        )
        print(f"  Valid: {stats['n_valid']}")
        print(f"  Invalid: {stats['n_invalid']}")
    # --------------------------------------------------
    # 3. Return the complete result
    # --------------------------------------------------
    return bootstrap_result
def _runPartialCorrelationAnalysis( data, metric_names, target_names, ):
    partial_result = dop.calculatePairedPartialCorrelationStatistics(
        data=data,
        metric_names=metric_names,
        target_names=target_names, )
    print("\nPartial Correlation Analysis")
    print("-" * 70)
    for name, value in partial_result.items():
        print(f"{name}: {value}")
    return partial_result
def staristicAnalysisForPSDdecay(
    psdData=None,
    metrics=None,
    kineticdata=None,
    metric_names=None,
    target_names=None,
    regions=None,
    statistic_names=("Pearson_r", "R2"),
    n_bootstrap=1000,
    confidence_level=0.95,
    random_state=42,
    target_data=None,
    check_rank=False,
):
    """
    Perform paired bootstrap correlation analysis between two groups
    of variables.
    Parameters
    ----------
    psdData : DataFrame or None
        PSD data used to calculate metrics when metrics is None.
    metrics : DataFrame or None
        DataFrame containing the first group of variables.
        If None, calculate PSD metrics from psdData.
    kineticdata : DataFrame or None
        DataFrame containing the target variables in the original
        two-DataFrame workflow. Used when target_data is None and
        metrics and kineticdata are separate DataFrames.
    metric_names : list, tuple, str, or None
        Column names of the first variable group.
        If None, use all columns in metrics.
    target_names : list, tuple, str, or None
        Column names of the second variable group.
        If None, use ["E (J/mol)"] when kineticdata is provided;
        otherwise, all columns in the selected target DataFrame
        must be specified explicitly or inferred as appropriate.
    regions : dict or None
        PSD regions used by calculatePSDMetrics when metrics is None.
    target_data : DataFrame or None
        Optional DataFrame containing the second variable group.
        If None, use kineticdata when provided; otherwise, use metrics.
        Set target_data=metrics to select both groups from the same
        DataFrame explicitly.
    check_rank : bool
        If True, print the ranks of selected variables for inspection.
    """
    # --------------------------------------------------
    # 1. Calculate PSD metrics if not already supplied
    # --------------------------------------------------
    if metrics is None:
        if psdData is None:
            raise ValueError( "Provide either metrics or psdData." )
        if regions is None:
            regions = {
                "High": (0.30, 0.50),
                "Middle": (0.50, 0.65),
                "Low": (0.65, 0.80),
            }
        metrics = calculatePSDMetrics( psdData, regions, )
    if not isinstance(metrics, pd.DataFrame):
        raise TypeError("metrics must be a pandas DataFrame.")
    # --------------------------------------------------
    # 2. Determine the source of target variables
    # --------------------------------------------------
    if target_data is not None:
        target_df = target_data
    elif kineticdata is not None:
        target_df = kineticdata
    else:
        target_df = metrics
    if not isinstance(target_df, pd.DataFrame):
        raise TypeError( "kineticdata or target_data must be a pandas DataFrame." )
    # --------------------------------------------------
    # 3. Normalize sample identifiers
    # --------------------------------------------------
    metrics = metrics.copy()
    target_df = target_df.copy()
    if "sample" in metrics.columns:
        metrics = metrics.set_index("sample")
    if "sample" in target_df.columns:
        target_df = target_df.set_index("sample")
    if not metrics.index.is_unique:
        raise ValueError( "metrics contains duplicate sample indices." )
    if not target_df.index.is_unique:
        raise ValueError( "The target DataFrame contains duplicate sample indices." )
    # Sort only when using the existing project's sorting function.
    metrics = dop.naturalSortData(metrics)
    target_df = dop.naturalSortData(target_df)
    # --------------------------------------------------
    # 4. Select variable names
    # --------------------------------------------------
    if metric_names is None:
        metric_names = list(metrics.columns)
    elif isinstance(metric_names, str):
        metric_names = [metric_names]
    else:
        metric_names = list(metric_names)
    if target_names is None:
        if "E (J/mol)" in target_df.columns:
            target_names = ["E (J/mol)"]
        elif target_df is metrics:
            raise ValueError( "Specify target_names when using the same DataFrame without an 'E (J/mol)' column." )
        else:
            raise ValueError( "Specify target_names to identify the target variables." )
    elif isinstance(target_names, str):
        target_names = [target_names]
    else:
        target_names = list(target_names)
    if not metric_names:
        raise ValueError("metric_names cannot be empty.")
    if not target_names:
        raise ValueError("target_names cannot be empty.")
    missing_metrics = [ name for name in metric_names if name not in metrics.columns ]
    missing_targets = [ name for name in target_names if name not in target_df.columns ]
    if missing_metrics:
        raise KeyError( f"Metric columns not found: {missing_metrics}" )
    if missing_targets:
        raise KeyError( f"Target columns not found: {missing_targets}" )
    # --------------------------------------------------
    # 5. Align observations by sample index
    # --------------------------------------------------
    common_samples = metrics.index[ metrics.index.isin(target_df.index) ]
    if len(common_samples) < 2:
        raise ValueError( "Fewer than two matching samples between the two variable groups." )
    selected_metrics = metrics.loc[ common_samples, metric_names ].copy()
    selected_targets = target_df.loc[ common_samples, target_names ].copy()
    # --------------------------------------------------
    # 6. Combine data into a DataFrame and numeric matrix
    # --------------------------------------------------
    duplicate_names = set(metric_names).intersection(target_names)
    if duplicate_names:
        raise ValueError(
            "The two variable groups contain overlapping column names: "
            f"{sorted(duplicate_names)}. Rename one group before analysis."
        )
    combined_df = pd.concat( [selected_metrics, selected_targets], axis=1, )
    try:
        combined_data = combined_df.to_numpy(dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError( "All selected variables must contain numeric data." ) from exc
    # --------------------------------------------------
    # 7. Map variable names to combined matrix column indices
    # --------------------------------------------------
    n_metrics = len(metric_names)
    combined_metric_names = { name: i for i, name in enumerate(metric_names) }
    combined_target_names = { name: n_metrics + i for i, name in enumerate(target_names) }
    # --------------------------------------------------
    # 8. Run Bootstrap analysis
    # --------------------------------------------------
    bootstrap_result = _runBootstrapAnalysis( data=combined_data, metric_names=combined_metric_names, 
            target_names=combined_target_names, statistic_names=statistic_names, n_bootstrap=n_bootstrap, 
            confidence_level=confidence_level, random_state=random_state, )
    # --------------------------------------------------
    # 9. Run Leave-One-Out analysis
    # --------------------------------------------------
    loo_result = _runLOOAnalysis( data=combined_data, metric_names=combined_metric_names, 
            target_names=combined_target_names, statistic_names=statistic_names, )
    # --------------------------------------------------
    # 10. Run PCA analysis
    # --------------------------------------------------
    pca_result = dop.runPCAAnalysis( data=selected_metrics, )
    # --------------------------------------------------
    # 11. Run VIF analysis
    # --------------------------------------------------
    vif_result = dop.runVIFAnalysis( data=selected_metrics, )
    # --------------------------------------------------
    # 12. Run Partial Correlation analysis
    # --------------------------------------------------
    partial_correlation_result = _runPartialCorrelationAnalysis( data=combined_data, 
            metric_names=combined_metric_names, target_names=combined_target_names, )
    return {
        "samples": list(common_samples),
        "metrics": selected_metrics,
        "kineticdata": selected_targets,
        "combined_data": combined_df,
        "bootstrap_result": bootstrap_result,
        "loo_result": loo_result,
        "pca_result": pca_result,
        "vif_result": vif_result,
        "partial_correlation_result": partial_correlation_result,
    }
def _prepareOrderedRegionLayout(regions):
    """
    Validate regions and return their pore-width-ordered layout.
    Returns
    -------
    region_names : list
        Region names ordered by increasing pore width.
    base_boundaries : numpy.ndarray
        Outer boundaries and all internal boundaries.
    """
    if not isinstance(regions, dict) or len(regions) < 2:
        raise ValueError(
            "regions must be a dictionary containing at least "
            "two contiguous pore-width regions."
        )
    region_items = []
    for name, bounds in regions.items():
        if len(bounds) != 2:
            raise ValueError(
                f"Region '{name}' must contain exactly two boundaries."
            )
        lower, upper = map(float, bounds)
        if not np.isfinite(lower) or not np.isfinite(upper):
            raise ValueError(
                f"Region '{name}' contains non-finite boundaries."
            )
        if lower >= upper:
            raise ValueError(
                f"Invalid boundaries for region '{name}': {bounds}"
            )
        region_items.append((name, lower, upper))
    # Sort by increasing pore width, independent of dict insertion order.
    region_items.sort(key=lambda item: item[1])
    for i in range(1, len(region_items)):
        previous_upper = region_items[i - 1][2]
        current_lower = region_items[i][1]
        if not np.isclose(
            previous_upper,
            current_lower,
            atol=1e-8,
            rtol=0,
        ):
            raise ValueError(
                "Regions must be contiguous and non-overlapping. "
                f"Gap or overlap detected between "
                f"'{region_items[i - 1][0]}' and "
                f"'{region_items[i][0]}'."
            )
    region_names = [item[0] for item in region_items]
    base_boundaries = np.array(
        [region_items[0][1]]
        + [item[2] for item in region_items],
        dtype=float,
    )
    return region_names, base_boundaries
def _buildRegionsFromBoundaries(region_names, boundaries):
    """Reconstruct contiguous regions while preserving their names."""
    return {
        name: (
            float(boundaries[i]),
            float(boundaries[i + 1]),
        )
        for i, name in enumerate(region_names)
    }
def _resolveCompetitionIndexColumn(metrics, ci_column=None):
    """
    Identify the CompetitionIndex column in calculatePSDMetrics output.
    """
    columns = list(metrics.columns)
    if ci_column is not None:
        if ci_column in columns:
            return ci_column
        # Permit a prefixed column such as GCMC_CompetitionIndex.
        requested = re.sub(
            r"[^a-z0-9]", "", str(ci_column).lower()
        )
        matches = [
            col for col in columns
            if re.sub(
                r"[^a-z0-9]", "", str(col).lower()
            ).endswith(requested)
        ]
        if len(matches) == 1:
            return matches[0]
        raise KeyError(
            f"Cannot uniquely identify CI column '{ci_column}'. "
            f"Available columns: {columns}"
        )
    normalized = {
        col: re.sub(r"[^a-z0-9]", "", str(col).lower())
        for col in columns
    }
    matches = [
        col for col, name in normalized.items()
        if "competitionindex" in name
    ]
    if len(matches) == 1:
        return matches[0]
    if len(matches) > 1:
        raise ValueError(
            "Multiple CompetitionIndex columns were found. "
            "Specify ci_column explicitly."
        )
    # Fallback for a function that names the output simply "CI".
    matches = [
        col for col, name in normalized.items()
        if name == "ci"
    ]
    if len(matches) == 1:
        return matches[0]
    raise KeyError(
        "Could not identify the CompetitionIndex column. "
        "Specify ci_column explicitly. "
        f"Available columns: {columns}"
    )
def _prepareSensitivityTargets(
    target_data,
    target_names,
    sample_index,
):
    """
    Align target data to the sample index of calculatePSDMetrics.
    """
    if isinstance(target_names, dict):
        names = list(target_names.keys())
    elif target_names is None:
        if isinstance(target_data, pd.DataFrame):
            names = list(target_data.columns)
        elif isinstance(target_data, pd.Series):
            names = [target_data.name or "Target"]
        else:
            raise ValueError(
                "target_names must be provided when target_data "
                "is not a DataFrame or Series."
            )
    elif isinstance(target_names, str):
        names = [target_names]
    else:
        names = list(target_names)
    if not names:
        raise ValueError("At least one target variable is required.")
    aligned = {}
    if isinstance(target_data, pd.DataFrame):
        for name in names:
            if name not in target_data.columns:
                raise KeyError(
                    f"Target '{name}' not found in target_data. "
                    f"Available columns: {list(target_data.columns)}"
                )
            aligned[name] = target_data[name].reindex(sample_index)
    elif isinstance(target_data, pd.Series):
        if len(names) != 1:
            raise ValueError(
                "A Series can be used with only one target name."
            )
        aligned[names[0]] = target_data.reindex(sample_index)
    elif isinstance(target_data, dict):
        for name in names:
            if name not in target_data:
                raise KeyError(
                    f"Target '{name}' not found in target_data."
                )
            values = target_data[name]
            if isinstance(values, pd.Series):
                aligned[name] = values.reindex(sample_index)
            else:
                values = np.asarray(values, dtype=float).reshape(-1)
                if len(values) != len(sample_index):
                    raise ValueError(
                        f"Target '{name}' has {len(values)} values, "
                        f"but {len(sample_index)} samples are expected."
                    )
                aligned[name] = pd.Series(
                    values, index=sample_index
                )
    else:
        values = np.asarray(target_data, dtype=float)
        if values.ndim == 1:
            values = values.reshape(-1, 1)
        if values.ndim != 2:
            raise ValueError(
                "target_data must be a 1D or 2D array, "
                "Series, DataFrame, or dictionary."
            )
        if values.shape != (len(sample_index), len(names)):
            raise ValueError(
                "Array target_data must have shape "
                "(n_samples, n_targets), with rows in the same "
                "order as the sample index of calculatePSDMetrics."
            )
        for j, name in enumerate(names):
            aligned[name] = pd.Series(
                values[:, j], index=sample_index
            )
    return aligned
def runCompetitionIndexSensitivityAnalysis(
    psdData,
    target_data,
    regions,
    target_names=None,
    boundary_offsets=(-0.05, -0.025, 0.0, 0.025, 0.05),
    prefix=None,
    ci_column=None,
):
    """
    Evaluate the sensitivity of CI-target correlations to changes
    in internal pore-size boundaries.
    Outer pore-width boundaries remain fixed. All valid combinations
    of the specified offsets are evaluated.
    Parameters
    ----------
    psdData : DataFrame
        PSD data with MultiIndex columns, as used by
        calculatePSDMetrics().
    target_data : DataFrame, Series, dict, or array-like
        Target data, e.g., D-A characteristic energy E.
    regions : dict
        Original contiguous pore-width regions.
    target_names : sequence of str, dict, or str, optional
        Target variable names.
    boundary_offsets : sequence of float
        Offsets in nm applied independently to each internal boundary.
        Must include zero to evaluate the original boundaries.
    prefix : str, optional
        Prefix passed to calculatePSDMetrics().
    ci_column : str, optional
        Explicit CompetitionIndex column name, if required.
    Returns
    -------
    dict
        boundary_results : DataFrame
            Pearson correlation and R2 for each boundary combination.
        summary : DataFrame
            Summary of boundary sensitivity.
        base_regions : dict
            Original region definitions.
    """
    # --------------------------------------------------
    # 1. Validate inputs
    # --------------------------------------------------
    if psdData is None:
        raise ValueError("psdData must not be None.")
    if target_data is None:
        raise ValueError("target_data must not be None.")
    region_names, base_boundaries = (
        _prepareOrderedRegionLayout(regions)
    )
    offsets = np.asarray(boundary_offsets, dtype=float)
    # Accept 1D arrays and single-row/single-column arrays.
    if offsets.ndim == 2 and 1 in offsets.shape:
        offsets = offsets.ravel()
    if offsets.ndim != 1 or offsets.size == 0:
        raise ValueError(
            "boundary_offsets must be a non-empty 1D sequence "
            "or a single-row/single-column array."
        )
    if not np.all(np.isfinite(offsets)):
        raise ValueError(
            "boundary_offsets must contain only finite values."
        )
    offsets = np.unique(offsets)
    if not np.any(np.isclose(offsets, 0.0)):
        raise ValueError(
            "boundary_offsets must include 0.0 so the original "
            "boundary configuration can be evaluated."
        )
    # --------------------------------------------------
    # 2. Store original and fixed outer boundaries
    # --------------------------------------------------
    outer_lower = base_boundaries[0]
    outer_upper = base_boundaries[-1]
    base_internal = base_boundaries[1:-1]
    base_regions = _buildRegionsFromBoundaries(
        region_names,
        base_boundaries,
    )
    # Calculate the original metrics to obtain the sample index.
    base_metrics = calculatePSDMetrics(
        psdData=psdData,
        regions=base_regions,
        prefix=prefix,
    )
    target_series = _prepareSensitivityTargets(
        target_data=target_data,
        target_names=target_names,
        sample_index=base_metrics.index,
    )
    # --------------------------------------------------
    # 3. Generate valid boundary combinations
    # --------------------------------------------------
    cutoff_combinations = {}
    n_internal = len(base_internal)
    for shifts in product(offsets, repeat=n_internal):
        candidate_internal = base_internal + np.asarray(
            shifts,
            dtype=float,
        )
        candidate_boundaries = np.concatenate([
            [outer_lower],
            candidate_internal,
            [outer_upper],
        ])
        # Skip invalid or overlapping regions.
        if np.any(np.diff(candidate_boundaries) <= 0):
            continue
        key = tuple(np.round(candidate_internal, 10))
        cutoff_combinations[key] = candidate_internal
    # Ensure that the original configuration is included.
    baseline_key = tuple(np.round(base_internal, 10))
    if baseline_key not in cutoff_combinations:
        cutoff_combinations[baseline_key] = (
            base_internal.copy()
        )
    # --------------------------------------------------
    # 4. Recalculate CI and evaluate each target
    # --------------------------------------------------
    records = []
    for candidate_internal in cutoff_combinations.values():
        candidate_boundaries = np.concatenate([
            [outer_lower],
            candidate_internal,
            [outer_upper],
        ])
        candidate_regions = _buildRegionsFromBoundaries(
            region_names,
            candidate_boundaries,
        )
        candidate_metrics = calculatePSDMetrics(
            psdData=psdData,
            regions=candidate_regions,
            prefix=prefix,
        )
        resolved_ci_column = _resolveCompetitionIndexColumn(
            candidate_metrics,
            ci_column=ci_column,
        )
        ci_values = candidate_metrics[
            resolved_ci_column
        ].reindex(base_metrics.index)
        is_baseline = np.allclose(
            candidate_internal,
            base_internal,
            atol=1e-9,
            rtol=0,
        )
        for target_name, target_series_i in target_series.items():
            y_values = target_series_i.reindex(
                candidate_metrics.index
            )
            # Use the existing correlation and linear-fitting workflow.
            analysis_result = dop.correlationAnalysis(
                x=ci_values.to_numpy(dtype=float),
                y=y_values.to_numpy(dtype=float),
                x_name=str(resolved_ci_column),
                y_name=str(target_name),
                dropna=True,
                fit_method="linear",
            )
            row = {
                "Target": target_name,
                "CI_column": resolved_ci_column,
                "IsBaseline": bool(is_baseline),
                "N": analysis_result["N"],
                "Pearson_r": analysis_result["Pearson_r"],
                "R2": analysis_result["R2"],
            }
            for i, boundary in enumerate(
                candidate_internal,
                start=1,
            ):
                row[f"Boundary_{i}_nm"] = float(boundary)
            records.append(row)
    boundary_results = pd.DataFrame(records)
    # --------------------------------------------------
    # 5. Summarize sensitivity relative to original boundaries
    # --------------------------------------------------
    summary_records = []
    for target_name, group in boundary_results.groupby(
        "Target",
        sort=False,
    ):
        baseline_rows = group[group["IsBaseline"]]
        if len(baseline_rows) != 1:
            raise RuntimeError(
                f"Expected one baseline result for '{target_name}', "
                f"but found {len(baseline_rows)}."
            )
        baseline_row = baseline_rows.iloc[0]
        r_values = group[
            "Pearson_r"
        ].dropna().to_numpy(dtype=float)
        r2_values = group[
            "R2"
        ].dropna().to_numpy(dtype=float)
        if len(r_values) == 0 or len(r2_values) == 0:
            raise ValueError(
                f"No valid Pearson_r or R2 values for target "
                f"'{target_name}'."
            )
        baseline_r = float(baseline_row["Pearson_r"])
        baseline_r2 = float(baseline_row["R2"])
        # Pearson correlation sensitivity.
        summary_row = {
            "Target": target_name,
            "N_boundary_combinations": len(group),
            "N_samples": int(baseline_row["N"]),
            "Baseline_Pearson_r": baseline_r,
            "Pearson_r_min": float(np.min(r_values)),
            "Pearson_r_max": float(np.max(r_values)),
            "Pearson_r_min_abs": float(
                np.min(np.abs(r_values))
            ),
            "Pearson_r_max_abs_change": float(
                np.max(np.abs(r_values - baseline_r))
            ),
            "Pearson_r_sign_consistency_pct": float(
                100.0 * np.mean(
                    np.sign(r_values) == np.sign(baseline_r)
                )
            ),
            # R2 sensitivity.
            "Baseline_R2": baseline_r2,
            "R2_min": float(np.min(r2_values)),
            "R2_max": float(np.max(r2_values)),
            "R2_max_abs_change": float(
                np.max(np.abs(r2_values - baseline_r2))
            ),
        }
        summary_records.append(summary_row)
    summary = pd.DataFrame(summary_records)
    # --------------------------------------------------
    # 6. Display results
    # --------------------------------------------------
    print("\nCompetitionIndex Boundary Sensitivity Analysis")
    print("-" * 70)
    print("\nBoundary sensitivity summary:")
    print(summary.to_string(index=False))
    print("\nResults for individual boundary combinations:")
    print(boundary_results.to_string(index=False))
    return {
        "boundary_results": boundary_results,
        "summary": summary,
        "base_regions": base_regions,
    }
##########################################################################
def main(file_path=None):
    # # 选择文件
    psd_sheet = "PSD-GCMC-CO2"
    kinetic_sheet = "thermalkinetic"
    regions_DFT = {
        "High": (0.3,0.65),
        "Middle": (0.65,0.9),
        "Low": (0.9,1)}
    regions_GCMC = {
        "High": (0.3,0.50),
        "Middle": (0.5,0.65),
        "Low": (0.65,0.8)}
    regions = regions_GCMC
    xColumns=[
        "HighLowRatio",
        "CompetitionIndex",
        "Skewness"
    ]
    yColumns=[
        "bA-T25",
        "SELE pyIAST",
        "E (J/mol)"
    ]
    file = fl.getFile()
    psdData, out_path, validFile, pdsMeta = fl.readFileBySheetWithMultiLevelHeader(file, psd_sheet, expand="correlation-GCMC_stati")
    kineticdata, kineticMeta = fl.readTableBySheet(file, kinetic_sheet)
    # # # ##################### 用于分析指定样品##############
    # # =========================
    # # Select samples
    # # =========================
    # samples = [
    #     "S4",
    #     "S5",
    #     "S8",
    #     "S10",
    # ]
    # psdData = selectSamples( psdData, samples, axis=1 )
    # kineticdata = selectSamples( kineticdata, samples, sampleColumn="Sample", axis=0 )
    # # # ##################### 用于分析指定样品##############
    psdData = dop.naturalSortData(psdData,axis=1)
    # plotData = {}
    # samples = psdData.columns.get_level_values(0).unique()
    # plotData = { sample: psdData[sample].copy() for sample in samples }
    # myPlt.plotCurve(
    #     data=plotData,
    #     x="Pore size(nm)",
    #     y="dV/dW (cm3/g·nm)",
    #     xlabel="Pore size(nm)",
    #     ylabel="dV/dW (cm3/g·nm)",
    #     marker= False,
    #     line=True,
    #     fit_label="PSD",
    #     # savepath="psd-DFT.png",
    # )
    kineticdata = dop.naturalSortData(kineticdata)
    metrics = calculatePSDMetrics(psdData, regions)
    # # ##################### 同时计算两个kernal的数据##############
    # psd_sheet_dft = "PSD-DFT-CO2"
    # psdData_dft, out_path_dft0, validFile, pdsMeta = fl.readFileBySheetWithMultiLevelHeader(file, psd_sheet_dft )
    # out_path = fl.get_expanded_name(out_path, fileName = "correlationMetrix", expand="all0809", expandPos=False)
    # psdData_dft = dop.naturalSortData(psdData_dft,axis=1)
    # metrics_GCMC = calculatePSDMetrics(psdData, regions_GCMC, prefix="GCMC")
    # metrics_dft = calculatePSDMetrics(psdData_dft, regions_DFT, prefix="DFT")
    # metrics = pd.concat( [metrics_GCMC, metrics_dft], axis=1 )
    # xColumns = [
    #     "GCMC_Centroid",
    #     "GCMC_CompetitionIndex",
    #     "GCMC_Skewness",
    #     "DFT_Centroid",
    #     "DFT_CompetitionIndex",
    #     "DFT_Skewness",
    #     ]
    # # ##################### 结束同时计算两个kernal的数据##############
    X, Y, summary, results = dop.crossCorrelationAnalysis( metrics, kineticdata,xColumns,yColumns=None )
    X_metrix, Y_metrix, summary_metrix, results_metrix = dop.matrixCorrelationAnalysis( metrics, kineticdata)
    # # ##################### 出传统porosity指标##############
    # porosity_sheet = "porosity"
    # include = ["BET", "Vultra (0.7)", "Vultra (1)", "Vmic(2)", "Vme", "Vt", "SELE pyIAST", "sel_henry"]
    # porosityData, porosityMeta = fl.readTableBySheet(file, porosity_sheet)
    # porosityData = dop.naturalSortData(porosityData,axis=1)
    # X_metrix_porspsity, Y_metrix_porspsity, summary_metrix_porspsity, results_metrix_porspsity = dop.matrixCorrelationAnalysis(porosityData, dataY=None, columns=None)
    # figureName = fl.get_expanded_name(out_path, fileName = "correlationMetrix-porosity", expand="", expandPos=True, type="png")
    # order = include
    # myPlt.plotCorrelogram(results = results_metrix_porspsity, order=order, include=include, cmap="RdBu_r",decimals=2, dpi=900, 
    #                       savePath = figureName
    #                       )
    # # ##################### 结束出传统porosity指标##############
    # ##################### 单独计算孔体积用##############
    # # X_metrix_po, Y_metrix_po, summary_metrix_po, results_metrix_po = dop.matrixCorrelationAnalysis( kineticdata)
    # # 不再使用3FLEX导出的累计体积
    # region_volume_GCMC = [0.3,0.45,0.5,0.55,0.6,0.65,0.7,0.8]
    # region_volume_DFT = [0.3,0.5,0.55,0.6,0.65,0.7,0.8,0.9]
    # volume_kernal=calculatePSDVolumes( psdData, region_volume_GCMC)
    # X_metrix_po, Y_metrix_po, summary_metrix_po, results_metrix_po = dop.matrixCorrelationAnalysis(volume_kernal, kineticdata)
    # out_path1= fl.get_expanded_name(out_path, fileName = "correlationMetrix-volume_gcmc", expand="", expandPos=True,)
    # # figureName = fl.get_expanded_name(out_path, fileName = "correlationMetrix-porosity2 r2", expand="", expandPos=True, type="png")
    # fl.export_to_excel_auto( volume_kernal, filename=out_path1, sheet_name="volume" )
    # fl.exportCorrelationExcel( out_path1,X_metrix_po, Y_metrix_po, X_matrix=None, Y_matrix=None, summary_matrix=summary_metrix_po,summary_matrix_SHHET="prosity2" )
    # myPlt.plotCorrelogram(results = results_metrix_po,textValue="R2", exclude=exclude, cmap="RdBu_r",decimals=2, dpi=900, 
    #                     #   savePath = figureName
    #                       )
    # ##################### 结束单独计算孔体积用##############
    # # ##################### 保存数据到excel##############
    # fl.export_to_excel_auto( kineticdata, filename=out_path, sheet_name="kineticdata" )
    # fl.export_to_excel_auto( metrics, filename=out_path, sheet_name="metrics" )
    # fl.exportCorrelationExcel( out_path, X, Y, summary, X_matrix=None, Y_matrix=None, summary_matrix=summary_metrix )
    # # ##################### 结束保存数据到excel##############
    ########################进行bootstrap leave-one-out PCA VIF Partial_correlation分析###################################
    # statisticAnalysis=staristicAnalysisForPSDdecay( psdData, metrics, kineticdata, metric_names=[ "Centroid", "CompetitionIndex", "Skewness" ], 
    #         target_names=["E (J/mol)"], regions = regions, statistic_names=("Pearson_r", "R2"), 
    #         n_bootstrap=10000, confidence_level=0.95, random_state=42, )
    # offsets = (-0.05, -0.025, 0.0, 0.025, 0.05)
    # ci_sensitivity_result = runCompetitionIndexSensitivityAnalysis( psdData=psdData, target_data=kineticdata, 
    #         regions=regions, target_names=["E (J/mol)"], boundary_offsets=offsets, prefix=None, )
    # savePSDAnalysisResults( statisticAnalysis=statisticAnalysis, ci_sensitivity_result=ci_sensitivity_result, kineticdata=kineticdata,
    #                        metrics=metrics, filename=out_path, correlation_export={ 
    #                            "X": X, "Y": Y, "summary": summary, "X_matrix": None, "Y_matrix": None, "summary_matrix": summary_metrix, }, )
    ########################结束进行bootstrap leave-one-out PCA VIF Partial_correlation分析###################################
    # # ##################### 输出图片部分##############
    # myPlt.plotBarByMetrics(metrics, columns=xColumns)
    # myPlt.plotBar(metrics.index, metrics["HighLowRatio"], xlabel="Sample",ylabel="HighLowRatio", figsize=(4.5, 3.5),gradientFlag=True,savepath="DFT-HighLowRatio")
    # myPlt.plotBar(metrics.index, metrics["CompetitionIndex"], xlabel="Sample",ylabel="CompetitionIndex", figsize=(4.5, 3.5),gradientFlag=True,savepath="DFT-CompetitionIndex")
    # myPlt.plotBar(metrics.index, metrics["Skewness"], xlabel="Sample",ylabel="Skewness", figsize=(4.5, 3.5),gradientFlag=True,savepath="DFT-Skewness")
    # myPlt.plotBar(metrics.index, metrics["Centroid"], xlabel="Sample",ylabel="Centroid (nm)", figsize=(4.5, 3.5),gradientFlag=True,showValue=True,valueRotation=90, valueInside=True, valuePosition=0.88, savepath="DFT-Centroid260907")
    # myPlt.plotBar(metrics.index, kineticdata["E (J/mol)"], xlabel="Sample",ylabel="E", figsize=(4, 3), gradientFlag=True,savepath="E")
    # myPlt.plotBar(metrics.index, kineticdata["bA-T25"], xlabel="Sample",ylabel="bA", figsize=(4, 3), gradientFlag=True,savepath="bA")
    # myPlt.plotBar(metrics.index, kineticdata["SELE pyIAST"], xlabel="Sample",ylabel="Selectivity (-)", gradientFlag=True,savepath="Sele")
    myPlt.plotSingleCorrelation(results_metrix["E (J/mol)"]["SELE pyIAST"], xlabel="E", ylabel="sele", text_position=(0.1, 0.95),figsize=(4.5, 3.5),
                                # savepath="E-sele_henry"
                                )
    myPlt.plotSingleCorrelation(results_metrix["bA-T25"]["b-N2"], xlabel="bA", ylabel="b-N2", text_position=(0.1, 0.95),figsize=(4.5, 3.5),
                                # savepath="bA-sele_henry"
                                )
    myPlt.plotSingleCorrelation(results_metrix["CompetitionIndex"]["E (J/mol)"], xlabel="CompetitionIndex", ylabel="E", text_position=(0.1, 0.95),figsize=(4.5, 3.5),
                                # savepath="DFT-E-CI"
                                )
    myPlt.plotSingleCorrelation(results_metrix["Centroid"]["E (J/mol)"], xlabel="Centroid", ylabel="E", text_position=(0.5, 0.95),figsize=(4.5, 3.5),
                                # savepath="DFT-E-Centroid"
                                )
    myPlt.plotSingleCorrelation(results_metrix["Skewness"]["E (J/mol)"], xlabel="Skewness", ylabel="E", text_position=(0.1, 0.95),figsize=(4.5, 3.5),
                                # savepath="DFT-E-Skewness"
                                )
    myPlt.plotSingleCorrelation(results_metrix["HighLowRatio"]["E (J/mol)"], xlabel="HighLowRatio", ylabel="E", text_position=(0.1, 0.95),figsize=(4.5, 3.5),savepath="DFT-E-HighLowRatio")
    myPlt.plotBatchCorrelation(results, topN=9)
    # # ##################### 结束输出图片部分##############
    exclude = [
        "Variance",
        "Std"
    ]
    include = [
        "O-EDS",
        "O(At)",
        "ID/IG",
        "co2_hen",
        "n2_hen",
        # "GCMC_High",
        # # "GCMC_HighLowRatio",
        # "GCMC_Centroid",
        # "GCMC_CompetitionIndex",
        # "GCMC_Skewness",
        # "DFT_High",
        # # "DFT_HighLowRatio",
        # "DFT_Centroid",
        # "DFT_CompetitionIndex",
        # "DFT_Skewness",
        "High",
        "HighLowRatio",
        "Centroid",
        "CompetitionIndex",
        "Skewness",
        "bA-T25",
        "b-N2",
        "E (J/mol)",
        "Eave (J/mol)",
        "SELE pyIAST",
        "sel_henry",
        #"sel_henry_numer",
    ]
    order = include
    figureName = fl.get_expanded_name(out_path, fileName = "correlationMetrix", expand="all", expandPos=True, type="png")
    myPlt.plotCorrelogram(results = results_metrix, textValue="Pearson_r",order= order, include=include,cmap="RdBu_r",decimals=2,
                        #   savePath = figureName
                          )
    # myPlt.plotCorrelogram(results = results_metrix, textValue="Pearson_r",cmap="RdBu_r",decimals=2,
    #                         #   savePath = figureName
    #                           )
    plt.show(block=True)
if __name__ == "__main__":
    f = sys.argv[1] if len(sys.argv) > 1 else None
    main(f)