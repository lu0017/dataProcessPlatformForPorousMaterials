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
    mask = (width >= lower) & (width <= upper)
    if np.sum(mask) < 2:
        return np.nan
    return np.trapezoid(psd[mask], width[mask])

def calculateRegionVolume0(width, psd, regions):
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
    for name, (lower, upper) in regions.items():
        result[name] = calculateRegionVolume( width, psd, lower, upper )
    result["Total"] = calculateRegionVolume( width, psd, width.min(), width.max() )
    return result

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
def calculatePSDMetrics(psdData, regions):
    """
    Calculate all PSD statistical descriptors.
    Parameters
    ----------
    psdData : DataFrame
        PSD data with MultiIndex columns.
    regions : dict
    Returns
    -------
    metrics : DataFrame
    """
    samples = dop.naturalSort(
        psdData.columns.get_level_values(0).unique()
    )
    metrics = {}
    for sample in samples:
        sampleData = psdData[sample]
        width = sampleData.iloc[:, 0].to_numpy()
        psd = sampleData.iloc[:, 1].to_numpy()
        regionVolumes = calculateRegionVolume0(
            width,
            psd,
            regions
        )
        highLow = calculateHighLowRatio(regionVolumes)
        competition = calculateCompetitionIndex(regionVolumes)
        moments = calculateDistributionMoment(
            width,
            psd
        )
        metrics[sample] = {}
        metrics[sample].update(highLow)
        metrics[sample].update(competition)
        metrics[sample].update(moments)
    metrics = pd.DataFrame.from_dict(
        metrics,
        orient="index"
    )
    metrics = dop.naturalSortData(metrics)
    return metrics


def main(file_path=None):
    # # 选择文件
    psd_sheet = "PSD-GCMC-CO2"
    kinetic_sheet = "thermalkinetic"

    regions = {
    "High": (0.33,0.50),
    "Middle": (0.50,0.6),
    "Low": (0.6,0.80)}
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
    psdData, out_path, validFile, pdsMeta = fl.readFileBySheetWithMultiLevelHeader(file, psd_sheet, expand="correlation")
    kineticdata, kineticMeta = fl.readTableBySheet(file, kinetic_sheet)
    psdData = dop.naturalSortData(psdData)
    kineticdata = dop.naturalSortData(kineticdata)
    metrics = calculatePSDMetrics(psdData, regions)
    X, Y, summary, results = dop.crossCorrelationAnalysis( metrics, kineticdata,xColumns,yColumns )
    X_metrix, Y_metrix, summary_metrix, results_metrix = dop.matrixCorrelationAnalysis( metrics, kineticdata)
    fl.export_to_excel_auto( kineticdata, filename=out_path, sheet_name="kineticdata" )
    fl.export_to_excel_auto( metrics, filename=out_path, sheet_name="metrics" )
    fl.exportCorrelationExcel( out_path, X, Y, summary, X_matrix=None, Y_matrix=None, summary_matrix=summary_metrix )
    myPlt.plotBarByMetrics(metrics, columns=xColumns)
    myPlt.plotBatchCorrelation(results, topN=9)
    exclude = [
        "Variance",
        "Std"
    ]
    figureName = fl.get_expanded_name(out_path, fileName = "correlationMetrix", expand="", expandPos=True, type="png")
    myPlt.plotCorrelogram(results = results_metrix, exclude=exclude,cmap="RdBu_r",decimals=2,savePath = figureName)
    plt.show(block=True)

if __name__ == "__main__":
    f = sys.argv[1] if len(sys.argv) > 1 else None
    main(f)