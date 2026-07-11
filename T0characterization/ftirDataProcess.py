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
def exportRaman( resultDF, filename="Raman.xlsx", sheet_name="fitData", ):
    """
    Export Raman fitting results.
    Parameters
    ----------
    resultDF : DataFrame
        Output of fitRaman().
    filename : str
        Excel filename.
    sheet_name : str or None
        Excel sheet name.
    """
    columns = [
        "Sample",
        "D_center",
        "D_height",
        "G_center",
        "G_height",
        "ID/IG",
        "AD/AG",
    ]
    df = resultDF.loc[:, columns].copy()
    fl.export_to_excel_auto( df, filename=filename, sheet_name=sheet_name, )
def calculateRamanProperties( fitResult, peaks, ):
    """
    Extract Raman peak properties and calculate Raman ratios.
    Parameters
    ----------
    fitResult : lmfit.ModelResult
        Fitting result.
    peaks : list
        Peak definitions.
    Returns
    -------
    dict
        Raman properties.
    """
    params = fitResult.params
    row = {}
    # -----------------------------
    # Peak properties
    # -----------------------------
    for peak in peaks:
        prefix = peak["name"] + "_"
        row[prefix + "center"] = params[prefix + "center"].value
        row[prefix + "sigma"] = params[prefix + "sigma"].value
        # lmfit 的 amplitude 为峰面积
        row[prefix + "area"] = params[prefix + "amplitude"].value
        if prefix + "height" in params:
            row[prefix + "height"] = params[prefix + "height"].value
        else:
            row[prefix + "height"] = np.nan
        if prefix + "fwhm" in params:
            row[prefix + "fwhm"] = params[prefix + "fwhm"].value
        else:
            row[prefix + "fwhm"] = np.nan
    # -----------------------------
    # Raman ratios
    # -----------------------------
    if (
        "D_height" in row
        and "G_height" in row
        and row["G_height"] != 0
    ):
        row["ID/IG"] = row["D_height"] / row["G_height"]
    if (
        "D_area" in row
        and "G_area" in row
        and row["G_area"] != 0
    ):
        row["AD/AG"] = row["D_area"] / row["G_area"]
    # -----------------------------
    # Goodness of fit
    # -----------------------------
    row["redchi"] = fitResult.redchi
    row["chisqr"] = fitResult.chisqr
    row["AIC"] = fitResult.aic
    row["BIC"] = fitResult.bic
    return row
def fitRaman( data, x="Raman shift (cm-1)", y="Intensity (a.u.)", ):
    """
    Batch Raman peak fitting.
    Parameters
    ----------
    data : dict
        {sample: DataFrame}
    x : str
        X column name.
    y : str
        Y column name.
    Returns
    -------
    resultDF : DataFrame
        Raman fitting results.
    fitResults : dict
        {sample: ModelResult}
    """
    peaks = [
        {
            "name": "D",
            "center": 1350,
            "center_min": 1300,
            "center_max": 1400,
            "sigma": 60,
            "sigma_min": 20,
            "sigma_max": 150,
        },
        {
            "name": "G",
            "center": 1580,
            "center_min": 1540,
            "center_max": 1620,
            "sigma": 40,
            "sigma_min": 20,
            "sigma_max": 150,
        },
    ]
    model = "Lorentzian"
    fitResults = {}
    results = []
    data = dop.cropX(data, xmin=1020, xmax=1710)
    myPlt.plotSpectrum(
        data,
        figsize=(8, 4),
        x=x,
        y=y,
        reverse_x=False,
        offset=50,
        legend=True,
        colors="PAPER1",
        linewidth=0.5,
    )
    plt.show(block=False)
    for sample, df in data.items():
        X = df[x].to_numpy()
        Y = df[y].to_numpy()
        try:
            fitResult = dop.fitPeak( x=X, y=Y, peaks=peaks, model=model, )
            fitResults[sample] = fitResult
            row = calculateRamanProperties( fitResult, peaks, )
            row["Sample"] = sample
            results.append(row)
            print(fitResult.params.keys())
            print(sample)
            print(row)
        except Exception as e:
            print(f"{sample} fitting failed: {e}")
    resultDF = pd.DataFrame(results)
    return resultDF, fitResults
# ==========================================
# 主函数
# ==========================================
def main(file_path=None):
    spectrum = "Raman"
    saveFile = False
    singleFile = False
    charaData = {}
    if spectrum == "FTIR":
        ########  FTIR 文件的设置  #########
        extension = ".csv"
        fileName = "FTIR"
        keyFlag = False
        keyName = ["Wavenumber", "Absorbance"]
        ncols = 2
        sampleMap = {
        "S1": ["CC-Hy-600-.5-1(2)"],
        "S2": ["CC-Hy-600-1-1"],
        "S3": ["CC-Hy-700-1-1(2)"],
        "S4": ["CC-Hy-800-.5-1"],
        "S5": ["CC-Hy-800-1-1"],
        "S6": ["CC-Hy-600-2-1"],
        "S7": ["CC-Hy-900-1-1"],
        "S8": ["CC-Hy-800-2-1"],
        "S9": ["CC-Hy-600-4-1"],
        "S10": ["CC-Hy-800-4-1"],
        }
    elif spectrum == "Raman":
        ########  raman 文件的设置  #########
        extension = ".txt"
        fileName = "Raman"
        keyFlag = False
        keyName = ["Raman shift (cm-1)", "Intensity (a.u.)"]
        ncols = 2
        sampleMap = {
        "S1": ["CC-Hy-600-.5-1"],
        "S2": ["CC-Hy-600-1-1"],
        "S3": ["CC-Hy-700-1-1"],
        "S4": ["CC-Hy-800-.5-1"],
        "S5": ["CC-Hy-800-1-1"],
        "S6": ["CC-Hy-600-2-1"],
        "S7": ["CC-Hy-900-1-1"],
        "S8": ["CC-Hy-800-2-1"],
        "S9": ["CC-Hy-600-4-1"],
        "S10": ["CC-Hy-800-4-1"],
        }
    if singleFile:
        # 选择一个文件
        file = fl.getFile()
        sample = os.path.splitext(os.path.basename(file))[0]
        df = fl.readTable( file=file, keyFlag=keyFlag, keyName=keyName, ncols=ncols, )
        if not df.empty:
            charaData[sample] = df
    else:
        # 选择文件夹
        folder = fl.select_folder()
        out_path = fl.get_expanded_name(folder, fileName = fileName)
        charaData = fl.readFolderData( folder=folder, extension=extension, keyFlag=keyFlag, keyName=keyName, ncols=ncols, )
    if saveFile:
        for i, (sample, df) in enumerate(charaData.items()):
            fl.saveDataframe2Excel( df, sample, out_path, sheet_name=fileName)
    print(f"Read {len(charaData)} samples.")
    charaData1 = dop.copySamples(charaData,sampleMap)
    charaData1 = dop.renameSamples(charaData,sampleMap)
    if spectrum == "FTIR":
        charaData1 = dop.cropX(charaData1, xmin=200, xmax=3990)
        charaData2 = dop.baselineCorrection( charaData1, x=keyName[0], y=keyName[1], method="asls", )
        charaData3 = dop.normalizeSpectrum( charaData2, y=keyName[1], method="vector", )
    elif spectrum == "Raman": 
        charaData2 = dop.smoothSpectrum( charaData1, y=keyName[1], method="gaussian", sigma=1.2, )
        charaData3 = dop.baselineCorrection( charaData2, x=keyName[0], y=keyName[1], method="asls", )
        resultDF, fitResults = fitRaman(charaData3)
        exportRaman( resultDF=resultDF, filename=out_path)
    myPlt.plotSpectrum(
        charaData1,
        figsize=(8, 4),
        x=keyName[0],
        y=keyName[1],
        reverse_x=True,
        offset=0,
        legend=True,
        colors="PAPER1",
        linewidth=0.5,
    )
    myPlt.plotSpectrum(
        charaData2,
        figsize=(8, 4),
        x=keyName[0],
        y=keyName[1],
        reverse_x=True,
        offset=20,
        legend=True,
        colors="PAPER1",
        linewidth=0.5,
    )
    myPlt.plotSpectrum(
        charaData3,
        figsize=(8, 4),
        x=keyName[0],
        y=keyName[1],
        reverse_x=True,
        offset=50,
        legend=True,
        colors="PAPER1",
        linewidth=0.5,
    )
    plt.show(block=False)
    
    y =0
if __name__ == "__main__":
    f = sys.argv[1] if len(sys.argv) > 1 else None
    main(f)