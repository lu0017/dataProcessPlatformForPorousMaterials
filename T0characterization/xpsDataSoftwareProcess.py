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
import T0characterization.xpsDataProcess as xpsAuto
def readXPSFitDat(folder, pattern="*.dat"):
    regionMap = {
        "C 1S": "C1s", "C1S": "C1s",
        "O 1S": "O1s", "O1S": "O1s",
        "N 1S": "N1s", "N1S": "N1s",
        "K 2P": "K2p", "K2P": "K2p",
    }
    columns = ["B.E.(eV)", "Raw Intensity", "Peak Sum", "Background",
               "Peak 1", "Peak 2", "Peak 3", "difference"]
    xps = {}
    for file in Path(folder).glob(pattern):
        parts = file.stem.split()
        if len(parts) < 3:
            continue
        sample = parts[0]
        region = regionMap.get(" ".join(parts[1:3]).upper(), "".join(parts[1:3]))
        df = pd.read_csv(file, sep=r"\s+", skiprows=1, header=None, names=columns)
        df = df.apply(pd.to_numeric, errors="coerce")
        df = df.dropna(subset=["B.E.(eV)", "Raw Intensity", "Peak Sum"])
        xps.setdefault(sample, {"metadata": {"Sample": sample}})[region] = df.reset_index(drop=True)
    return xps
def readXPSPar(folder, pattern="*.par"):
    regionMap = {
        "C 1S": "C1s", "C1S": "C1s",
        "O 1S": "O1s", "O1S": "O1s",
        "N 1S": "N1s", "N1S": "N1s",
        "K 2P": "K2p", "K2P": "K2p",
    }
    columns = ["Peak", "Position (eV)", "Area", "FWHM (eV)", "%GL (%)"]
    xps = {}
    for file in Path(folder).glob(pattern):
        parts = file.stem.split()
        if len(parts) < 3:
            continue
        sample = parts[0]
        region = regionMap.get(" ".join(parts[1:3]).upper(), "".join(parts[1:3]))
        df = pd.read_csv(file, sep=r"\s+", skiprows=1, header=None, names=columns)
        df = df.apply(pd.to_numeric, errors="coerce").dropna(subset=["Peak"])
        xps.setdefault(sample, {"metadata": {"Sample": sample}})[region] = df.reset_index(drop=True)
    return xps
def exportXPSPeakSummary(xps, filename="XPS_Peak_Summary.xlsx"):
    import pandas as pd
    columns = ["sample", "region", "name", "Position (eV)", "Area",
               "FWHM (eV)", "Area (%)", "%GL (%)"]
    rows = []
    for sample, sampleData in xps.items():
        for region, df in sampleData.items():
            if not isinstance(df, pd.DataFrame) or "name" not in df.columns:
                continue
            for _, peak in df.iterrows():
                rows.append({
                    "sample": sample,
                    "region": region,
                    "name": peak["name"],
                    **{col: peak.get(col, None) for col in columns[3:]},
                })
    summary = pd.DataFrame(rows, columns=columns)
    fl.export_to_excel_auto(summary, filename=filename, sheet_name="Peak Summary", index=False)
    return summary
def saveXPSRegions(xpsData, path, regions=None,sub_folder=None):
    path = fl.getPathOnly(path)
    if sub_folder:
        path = fl.check_or_create_folder(path, sub_folder_name=sub_folder)
    os.makedirs(path, exist_ok=True)
    for sample, sampleData in xpsData.items():
        selected = regions or [k for k, v in sampleData.items() if isinstance(v, pd.DataFrame)]
        for region in selected:
            df = sampleData.get(region)
            if isinstance(df, pd.DataFrame):
                filename = os.path.join(path, f"{sample}_{region}.txt")
                df.to_csv(filename, sep="\t", index=False, header=False)
                print(f"Export {sample}_{region} -> {filename}")
def labelXPSPeaks(xps, peakConfig=xpsAuto.XPS_PEAK_CONFIG):
    for sampleData in xps.values():
        for region, df in sampleData.items():
            if not isinstance(df, pd.DataFrame) or region not in peakConfig:
                continue
            peaks = peakConfig[region]["peaks"]
            names, matches = [], []
            for pos in df["Position (eV)"]:
                candidates = [
                    p for p in peaks
                    if p["center_bounds"][0] <= pos <= p["center_bounds"][1]
                ]
                peak = min(candidates or peaks, key=lambda p: abs(pos - p["center"]))
                names.append(peak["name"])
                matches.append(bool(candidates))
            df["name"] = names
            df["matched"] = matches
    return xps
def processXPSPar(xps):
    for sampleData in xps.values():
        for region, df in sampleData.items():
            if not isinstance(df, pd.DataFrame) or "Area" not in df.columns:
                continue
            df.sort_values("Position (eV)", inplace=True)
            df["Area (%)"] = df["Area"] / df["Area"].sum() * 100
            sampleData[region] = df.reset_index(drop=True)
    return xps
def applyXPSCalibration(data, offsetMap, energyColumn="Binding Energy"):
    """Apply sample-specific binding-energy offsets to all XPS regions based on C 1S: C=C at 284.8eV ."""
    calibratedData = copy.deepcopy(data)
    for sample, regions in calibratedData.items():
        offset = offsetMap.get(sample, 0.0)
        for regionData in regions.values():
            if isinstance(regionData, pd.DataFrame) and energyColumn in regionData.columns:
                regionData[energyColumn] += offset
    return calibratedData
def main(file_path=None):
    STEP = 3 #1 原始.csv数据读取导出 .txt，， #2 XPS peak完成C 1S的offset确认后设置并重新导出txt，，， #3 读取.par 计算% 画图
    saveFile = True
    singleFile = False
    highResolutionFlag = False #False时将自动从Suvey中提取对应的区段谱
    charaData = {}
    extension = "*.csv"
    if STEP == 3:
        extension = "*.par"
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
    # C 1s calibration offsets (eV)
    offsetMap = {
        "S1": 0.8,
        "S2": 0.0,
        "S3": 0.0,
        "S4": 0.0,
        "S5": 0.0,
        "S6": 0.0,
        "S7": 0.0,
        "S8": 0.0,
        "S9": 0.0,
        "S10": 0.0,
    }
    if STEP == 3:
        folder = fl.select_folder()
        out_path = fl.get_expanded_name(folder, fileName = "XPS_Peak_Summary")
        peakResult = readXPSPar(folder, pattern="*.par")
        peakResult = dop.copySamples(peakResult,sampleMap)
        peakResult = dop.renameSamples(peakResult,sampleMap)
        peakpPecent = processXPSPar(peakResult)
        peakinfor = labelXPSPeaks(peakpPecent)
        # exportXPSPeakSummary(peakinfor, filename=out_path)

        xpsFit = readXPSFitDat(folder)
    else:
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
            charaData = xpsAuto.readXPS( folder, pattern=extension, keyName = keyName, parser_kwargs=None, )
            charaData1 = dop.copySamples(charaData,sampleMap)
            charaData1 = dop.renameSamples(charaData,sampleMap)
            charaDataRegions = xpsAuto.extractXPSRegions( charaData1, highResolutionFlag=False, 
                                        regions=("C1s", "O1s"), x="Binding Energy", ) 
        if STEP == 2:
            offsetCharaDataRegions = applyXPSCalibration(charaDataRegions, offsetMap)
        out_path2 = fl.get_expanded_name(folder, fileName = "OFFSET")
        if saveFile:
            # xpsAuto.saveXPS( charaData1, filename=out_path )
            # xpsAuto.saveXPS( charaDataRegions, filename=out_path )
            if STEP == 1:
                xpsAuto.saveXPS( charaDataRegions, filename=out_path2 )
                saveXPSRegions(charaDataRegions, out_path2, regions=None, sub_folder="original")
            elif STEP == 2:
                saveXPSRegions(offsetCharaDataRegions, out_path2, regions=None, sub_folder="offset")
            # saveXPSRegions(offsetCharaDataRegions, out_path2, regions=None)
        print(f"Read {len(charaData)} samples.")
    y = 0
if __name__ == "__main__":
    f = sys.argv[1] if len(sys.argv) > 1 else None
    main(f)
