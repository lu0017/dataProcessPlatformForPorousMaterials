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
def staristicAnalysisForPSDdecay(
    psdData=None,
    metrics=None,
    kineticdata=None,
    metric_names=None,
    target_names=None,
    regions=None,
    statistic_names=("Pearson_r", "Spearman_r"),
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
            raise ValueError(
                "Provide either metrics or psdData."
            )
        if regions is None:
            regions = {
                "High": (0.30, 0.50),
                "Middle": (0.50, 0.65),
                "Low": (0.65, 0.80),
            }
        metrics = calculatePSDMetrics(
            psdData,
            regions,
        )
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
    # 8. Run the existing paired bootstrap analysis
    # --------------------------------------------------
    bootstrap_result = dop.runPairedBootstrapAnalysis(
        data=combined_data,
        statistic_func=lambda sample: (
            dop.calculatePairedCorrelationStatistics(
                data=sample,
                metric_names=combined_metric_names,
                target_names=combined_target_names,
                statistic_names=statistic_names,
            )
        ),
        n_bootstrap=n_bootstrap,
        confidence_level=confidence_level,
        random_state=random_state,
    )
    # --------------------------------------------------
    # 9. Optional rank inspection
    # --------------------------------------------------
    if check_rank:
        print("\nVariable ranks:")
        print( combined_df.rank(method="average").to_string() )
    # --------------------------------------------------
    # 10. Print statistics
    # --------------------------------------------------
    for name, stats in bootstrap_result["statistics"].items():
        print(f"\n{name}")
        print(f"  Estimate: {stats['estimate']:.4f}")
        print(
            f"  {confidence_level:.0%} CI: "
            f"{stats['confidence_interval']}"
        )
        print(
            "  Bootstrap SE: "
            f"{stats['bootstrap_standard_error']:.4f}"
        )
        print(f"  Valid: {stats['n_valid']}")
        print(f"  Invalid: {stats['n_invalid']}")
    return {
        "samples": list(common_samples),
        "metrics": selected_metrics,
        "kineticdata": selected_targets,
        "combined_data": combined_df,
        "bootstrap_result": bootstrap_result,
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
    psdData, out_path, validFile, pdsMeta = fl.readFileBySheetWithMultiLevelHeader(file, psd_sheet, expand="correlation-GCMC")
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
    ########################进行bootstrap leave-one-out PCA VIF分析###################################
    staristicAnalysisForPSDdecay( psdData, metrics, kineticdata, metric_names=[ "Centroid", "CompetitionIndex", "Skewness" ], 
            target_names=["E (J/mol)"], regions = regions, statistic_names=("Pearson_r","Spearman_r"), 
            n_bootstrap=10000, confidence_level=0.95, random_state=42, )
    ########################结束进行bootstrap leave-one-out PCA VIF分析###################################
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
    # myPlt.plotSingleCorrelation(results_metrix["bA-T25"]["sel_henry"], xlabel="bA", ylabel="sele", text_position=(0.1, 0.95),figsize=(4.5, 3.5),savepath="bA-sele_henry")
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
        "ID/IG",
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