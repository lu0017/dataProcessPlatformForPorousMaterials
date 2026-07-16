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
def readXPS( folder, pattern="*.csv", keyName=None, parser_kwargs=None, ):
    """
    Batch read XPS spectra from a folder.
    Returns
    -------
    dict
    {
        sample:{
            "metadata": metadata,
            "C1s": dataframe,
            "O1s": dataframe,
        }
    }
    """
    parser_kwargs = parser_kwargs or {}
    region_map = {
        "Su": "Survey",
        "Survey": "Survey",
        "C": "C1s",
        "C1s": "C1s",
        "O": "O1s",
        "O1s": "O1s",
        "N": "N1s",
        "N1s": "N1s",
        "S": "S2p",
        "S2p": "S2p",
        "P": "P2p",
        "P2p": "P2p",
    }
    xps = {}
    folder = Path(folder)
    for file in folder.glob(pattern):
        result = fl.parseGroupedTable(
            file=file,
            metadata_rows=(0,2),
            metadata_names=("Sample","Original"),
            data_start_row=4,
            group_size=2,
            data_names=("Binding Energy","Intensity"),
            **parser_kwargs,
        )
        metadata = result["metadata"]
        data = result["data"]
        samples = metadata["Sample"]
        regions = metadata["Original"]
        # 单谱情况
        if isinstance(samples, str):
            samples = [samples]
        if isinstance(regions, str):
            regions = [regions]
        for sample, regions in result["data"].items():
            if sample not in xps:
                xps[sample] = { "metadata": result["metadata"] }
            for region, df in regions.items():
                xps[sample][region] = df
    return xps
def readXPS0(
    folder,
    pattern="*.csv",
    keyName=None,
    parser_kwargs=None,
):
    """
    Batch read XPS spectra from a folder.
    Returns
    -------
    DataFrame
    Columns:
        Sample
        Region
        Binding Energy
        Intensity
    """
    parser_kwargs = parser_kwargs or {}
    region_map = {
        "Su": "Survey",
        "Survey": "Survey",
        "C": "C1s",
        "C1s": "C1s",
        "O": "O1s",
        "O1s": "O1s",
        "N": "N1s",
        "N1s": "N1s",
        "S": "S2p",
        "S2p": "S2p",
        "P": "P2p",
        "P2p": "P2p",
    }
    records = []
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
        metadata = result["metadata"]
        for sample, regions in result["data"].items():
            for region, df in regions.items():
                # 标准化region名称
                region = region_map.get(
                    region,
                    region
                )
                temp = df.copy()
                # 防止不同谱长度问题
                temp = temp.dropna(
                    axis=0,
                    how="any"
                )
                temp["Sample"] = sample
                temp["Region"] = region
                records.append(temp)
    if len(records) == 0:
        return pd.DataFrame()
    xps = pd.concat(
        records,
        ignore_index=True
    )
    # 调整列顺序
    cols = [
        "Sample",
        "Region",
        "Binding Energy",
        "Intensity",
    ]
    xps = xps[
        cols
    ]
    return xps
def main(file_path=None):
    spectrum = "XPS"
    saveFile = False
    singleFile = False
    charaData = {}
    if spectrum == "XPS":
        ########  FTIR 文件的设置  #########
        extension = ".csv"
        fileName = "XPS"
        keyFlag = False
        keyName = ["Binding Energy", "Intensity"]
        ncols = 2
        sampleMap = {
        "S1": ["CC-Hy-600-.5-1(2)"],
        "S2": ["6"],
        "S3": ["7"],
        "S4": ["CC-Hy-800-.5-1"],
        "S5": ["CC-Hy-800-1-1"],
        "S6": ["CC-Hy-600-2-1"],
        "S7": ["8"],
        "S8": ["CC-Hy-800-2-1"],
        "S9": ["9"],
        "S10": ["CC-Hy-800-4-1"],
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
    sample = next(iter(charaData))

    print(type(charaData[sample]["metadata"]))
    print(charaData[sample]["metadata"])
    if spectrum == "XPS":
        # charaData1 = dop.cropX(charaData1, xmin=700, xmax=3500)
        # charaData2 = dop.baselineCorrection( charaData1, x=keyName[0], y=keyName[1], method="airpls", )
        # charaData3 = dop.normalizeSpectrum( charaData2, y=keyName[1], method="vector", )
        # figureName = fl.get_expanded_name(out_path, fileName = "FTIR", type="png")
        myPlt.plotSpectrum( charaData2["N1s"], figsize=(8, 4), x=keyName[0], y=keyName[1], show_yticks = False,
                       reverse_x=True, offset=0, colors="PAPER1", linewidth=0.5, )
        # myPlt.plotSpectrum( charaData2, figsize=(8, 4), x=keyName[0], y=keyName[1], show_yticks = False,
        #                reverse_x=False, offset=-5, colors="PAPER1", linewidth=1.5, xbreak=(2100, 2700))
        # myPlt.plotSpectrum( charaData2, figsize=(8, 4), x=keyName[0], y=keyName[1], show_yticks = False,
        #                offset=-5, colors="PAPER1", linewidth=0.5, )
    elif spectrum == "Raman": 
        charaData1 = dop.cropX(charaData1, xmin=500, xmax=2200)
        charaData2 = dop.smoothSpectrum( charaData1, y=keyName[1], method="gaussian", sigma=1.2, )
        charaData3 = dop.baselineCorrection( charaData2, x=keyName[0], y=keyName[1], method="rubberband", )
        figureName = fl.get_expanded_name(out_path, fileName = "Raman", type="png")
        myPlt.plotSpectrum( charaData1, figsize=(8, 4), x=keyName[0], y=keyName[1], 
                       reverse_x=True, offset=0, colors="PAPER1", linewidth=0.5, )
        myPlt.plotSpectrum( charaData2, figsize=(8, 4), x=keyName[0], y=keyName[1], 
                       offset=80, colors="PAPER1", linewidth=0.5, )
        myPlt.plotSpectrum( charaData3, figsize=(8, 4), x=keyName[0], y=keyName[1], 
                       offset=50, colors="PAPER1", linewidth=0.5, savepath=figureName)
    plt.show(block=True)
    y =0
if __name__ == "__main__":
    f = sys.argv[1] if len(sys.argv) > 1 else None
    main(f)