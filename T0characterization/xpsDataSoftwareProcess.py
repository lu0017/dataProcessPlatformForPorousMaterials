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
XPS_DEFAULT_COLORS = ["#A000E8", "#4BA3F5", "#36A269", "#F39C45"]
XPS_PEAK_COLORS = {
    "C1s": {
        "C-C/C=C": "#4C78A8",
        "C-O":     "#F58518",
        "C=O":     "#54A24B",
        "O-C=O":   "#E45756",
    },
    "O1s": {
        "C=O":     "#4C78A8",
        "C-O":     "#F58518",
        "O-C=O":   "#54A24B",
    },
}
def getXPSPeakColumns(df):
    pattern = re.compile(r"^Peak\s*(\d+)$", re.I)
    peaks = []
    for c in df.columns:
        m = pattern.match(str(c).strip())
        if m:
            peaks.append((int(m.group(1)), c))
    return [c for _, c in sorted(peaks)]
def getXPSPeakColors(xpsFit, sample, region, defaultColors=None):
    peakNames = xpsFit.get(sample, {}).get(
        "metadata", {}
    ).get("peak_names", {}).get(region, {})
    regionColors = XPS_PEAK_COLORS.get(region, {})
    colors = {
        column: regionColors[name]
        for column, name in peakNames.items()
        if name in regionColors
    }
    peakColumns = getXPSPeakColumns(xpsFit[sample][region])
    fallbackColors = {
        column: XPS_DEFAULT_COLORS[i % len(XPS_DEFAULT_COLORS)]
        for i, column in enumerate(peakColumns)
    }
    return {**fallbackColors, **colors, **(defaultColors or {})}
def blendXPSColor(color, background="#FFFFFF", weight=1.0):
    c1, c2 = np.array(to_rgb(color)), np.array(to_rgb(background))
    return to_hex(np.clip(weight, 0, 1) * c1 + (1 - np.clip(weight, 0, 1)) * c2)
def plotXPSRawData(ax, df, color="#A8A8A8", markersize=3, linewidth=0.7, zorder=1):
    ax.plot(df["B.E.(eV)"], df["Raw Intensity"], "o", ms=markersize,
            mfc="none", mec=color, mew=linewidth, ls="none", zorder=zorder)
def plotXPSBackground(ax, df, color="black", linewidth=1.1, zorder=4):
    ax.plot(df["B.E.(eV)"], df["Background"], color=color,
            lw=linewidth, zorder=zorder)
def addGradientPeak(ax, x, y, baseline, color, n_layers=100, alpha=0.95, zorder=2):
    x, y = np.asarray(x, float), np.asarray(y, float)
    baseline = np.full_like(y, baseline, dtype=float) if np.isscalar(baseline) else np.asarray(baseline, float)
    valid = np.isfinite(x) & np.isfinite(y) & np.isfinite(baseline)
    x, y, baseline = x[valid], y[valid], baseline[valid]
    if len(x) < 2:
        return
    for i in range(n_layers):
        f0, f1 = i / n_layers, (i + 1) / n_layers
        y0, y1 = baseline + f0 * (y - baseline), baseline + f1 * (y - baseline)
        facecolor = blendXPSColor(color, weight=(i + 0.5) / n_layers)
        ax.fill_between(x, y0, y1, color=facecolor, alpha=alpha,
                        linewidth=0, zorder=zorder)
def plotXPSComponents(
    ax, df, colors=None, fill_base="background",
    n_layers=100, alpha=0.95, linewidth=0.8, zorder=2,
    xpsFit=None, sample=None, region=None
):
    peak_columns = getXPSPeakColumns(df)
    if xpsFit is not None and sample is not None and region is not None:
        color_map = getXPSPeakColors( xpsFit, sample, region, defaultColors=colors )
        peak_colors = [
            color_map.get(c, XPS_DEFAULT_COLORS[i % len(XPS_DEFAULT_COLORS)])
            for i, c in enumerate(peak_columns)
        ]
    else:
        palette = colors or XPS_DEFAULT_COLORS
        peak_colors = [ palette[i % len(palette)] for i in range(len(peak_columns)) ]
    x = df["B.E.(eV)"].to_numpy(float)
    if fill_base == "background":
        baseline = df["Background"].to_numpy(float)
    elif fill_base == "zero":
        baseline = np.zeros(len(df))
    else:
        raise ValueError("fill_base must be 'background' or 'zero'.")
    for i, (column, color) in enumerate(zip(peak_columns, peak_colors)):
        y = df[column].to_numpy(float)
        addGradientPeak( ax, x, y, baseline, color, n_layers, alpha, zorder + i * 0.01 )
        ax.plot( x, y, color=color, lw=linewidth, zorder=zorder + i * 0.01 + 0.001 )
    return peak_columns, peak_colors
def plotXPSTotalFit(ax, df, color="#C62828", linewidth=1.4, zorder=5):
    ax.plot(df["B.E.(eV)"], df["Peak Sum"], color=color, lw=linewidth, zorder=zorder)
def getXPSLegendLabels(xpsFit, sample, region):
    return xpsFit.get(sample, {}).get("metadata", {}).get( "peak_names", {} ).get(region, {})
def configureXPSAxes(ax, xlim=None, ylim=None, xticks=None, xtick_step=2,
                     reverse_x=True, show_x_ticks=True, show_y_ticks=False,
                     show_top_spine=True, show_right_spine=True,
                     spine_width=1, tick_length=5, tick_width=1,
                     tick_labelsize=11, axis_labelsize=13,
                     x_label="Binding energy (eV)", y_label="Intensity (a.u.)"):
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    if xticks is not None:
        ax.set_xticks(xticks)
    elif xtick_step is not None:
        lo, hi = sorted(xlim if xlim is not None else ax.get_xlim())
        ax.set_xticks(np.arange(np.ceil(lo / xtick_step) * xtick_step,
                                hi + xtick_step * 1e-8, xtick_step))
    if reverse_x:
        ax.set_xlim(sorted(ax.get_xlim(), reverse=True))
    ax.set_xlabel(x_label, fontsize=axis_labelsize)
    ax.set_ylabel(y_label, fontsize=axis_labelsize)
    ax.tick_params( axis="x", which="major", direction="in",
                    length=tick_length, width=tick_width,
                    labelsize=tick_labelsize,
                    bottom=show_x_ticks, top=False,
                    labelbottom=show_x_ticks, labeltop=False, )
    ax.tick_params(axis="y", which="both", left=show_y_ticks, labelleft=show_y_ticks)
    for side in ("left", "bottom", "top", "right"):
        ax.spines[side].set_linewidth(spine_width)
    ax.spines["top"].set_visible(show_top_spine)
    ax.spines["right"].set_visible(show_right_spine)
def buildXPSLegend(ax, peak_columns, peak_colors, labels=None, loc="upper right",
                   fontsize=10, frameon=True, show_raw=True,
                   show_background=True, show_total_fit=True, ncol=1):
    labels = labels or {}
    handles = []
    if show_raw:
        handles.append(Line2D([], [], color="#A8A8A8", marker="o", mfc="none",
                              mew=0.8, ms=4, ls="none",
                              label=labels.get("raw", "Extracted data")))
    if show_background:
        handles.append(Line2D([], [], color="black", lw=1.1,
                              label=labels.get("background", "Background")))
    for column, color in zip(peak_columns, peak_colors):
        handles.append(Patch(facecolor=blendXPSColor(color, weight=0.75),
                             edgecolor=color, lw=0.7,
                             label=labels.get(column, column)))
    if show_total_fit:
        handles.append(Line2D([], [], color="#C62828", lw=1.4,
                              label=labels.get("total", "Total fit")))
    return ax.legend(handles=handles, loc=loc, fontsize=fontsize,
                     frameon=frameon, ncol=ncol)
def plotXPSFit(xps, sample, region, colors=None, figsize=(7, 5),
               xlim=None, ylim=None, xticks=None, xtick_step=2,
               reverse_x=True, fill_base="background", show_legend=True,
               legend_loc="upper right", legend_fontsize=10, legend_labels=None,
               raw_markersize=3, peak_layers=100, peak_alpha=0.95,
               tick_labelsize=11, axis_labelsize=13, spine_width=1,
               title=None, save_path=None, dpi=600, show=True):
    df = xps[sample][region].copy()
    required = {"B.E.(eV)", "Raw Intensity", "Peak Sum", "Background"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")
    peak_columns = getXPSPeakColumns(df)
    if not peak_columns:
        raise ValueError(f"No peak columns found: {sample} / {region}")
    df = df.sort_values("B.E.(eV)").reset_index(drop=True)
    fig, ax = plt.subplots(figsize=figsize)
    plotXPSRawData(ax, df, markersize=raw_markersize)
    peak_columns, peak_colors = plotXPSComponents(
        ax, df, xpsFit=xps, sample=sample, region=region,
        colors=colors, fill_base=fill_base,
        n_layers=peak_layers, alpha=peak_alpha
    )
    plotXPSBackground(ax, df)
    plotXPSTotalFit(ax, df)
    configureXPSAxes(
        ax, xlim=xlim, ylim=ylim, xticks=xticks, xtick_step=xtick_step,
        reverse_x=reverse_x, tick_labelsize=tick_labelsize,
        axis_labelsize=axis_labelsize, spine_width=spine_width
    )
    if title:
        ax.set_title(title, fontsize=axis_labelsize)
    if show_legend:
        auto_labels = getXPSLegendLabels(xps, sample, region)
        legend_labels = {**auto_labels, **(legend_labels or {})}
        buildXPSLegend(
            ax, peak_columns, peak_colors, labels=legend_labels,
            loc=legend_loc, fontsize=legend_fontsize
        )
    fig.tight_layout()
    if save_path:
        path = Path(save_path)
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white")
    if show:
        plt.show(block=False)
    return fig, ax
def plotAllXPSFits(
    xpsFit, folder, figsize=(7, 1.5), fileType="png", **plotKwargs
):
    regionConfig = {
        "C1s": {"name": "C", "xlim": (280, 295), "xtick_step": 2},
        "O1s": {"name": "O", "xlim": (527.5, 537), "xtick_step": 2},
        "N1s": {"name": "N", "xlim": (395, 407), "xtick_step": 2},
        "K2p": {"name": "K", "xlim": (290, 300), "xtick_step": 2},
    }
    results = {}
    for sample, sampleData in xpsFit.items():
        results[sample] = {}
        for region, config in regionConfig.items():
            if region not in sampleData or not isinstance(sampleData[region], pd.DataFrame):
                continue
            kwargs = {
                "figsize": figsize,
                "xlim": config["xlim"],
                "xtick_step": config["xtick_step"],
                **plotKwargs,
            }
            figPath = fl.get_expanded_name(
                folder, fileName=f"{sample} {config['name']}", type=fileType
            )
            results[sample][region] = plotXPSFit(
                xpsFit, sample=sample, region=region,
                save_path=figPath, **kwargs
            )
    return results
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
def applyXPSPeakNames(xpsFit, peakinfor):
    for sample, sampleData in xpsFit.items():
        if sample not in peakinfor:
            continue
        for region, df in sampleData.items():
            peakInfo = peakinfor[sample].get(region)
            if not isinstance(df, pd.DataFrame) or not isinstance(peakInfo, pd.DataFrame):
                continue
            sampleData.setdefault("metadata", {}).setdefault(
                "peak_names", {}
            )[region] = {
                f"Peak {int(row['Peak']) + 1}": row["name"]
                for _, row in peakInfo.iterrows()
            }
    return xpsFit
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
def debugXPSPeakData(peakinfor, xpsFit):
    print("\n========== peakinfor ==========")
    print("Type:", type(peakinfor))
    if isinstance(peakinfor, dict):
        print("Keys:", list(peakinfor.keys())[:5])
        for key, value in list(peakinfor.items())[:2]:
            print(f"\nKey: {key!r}")
            print("Value type:", type(value))
            if isinstance(value, dict):
                print("Subkeys:", list(value.keys())[:10])
                for k, v in list(value.items())[:2]:
                    print(f"  {k!r}: type={type(v)}")
                    if isinstance(v, pd.DataFrame):
                        print(v.head(3).to_string())
                    else:
                        print(" ", repr(v)[:500])
            elif isinstance(value, pd.DataFrame):
                print("Columns:", value.columns.tolist())
                print(value.head(3).to_string())
            else:
                print(repr(value)[:500])
    print("\n========== xpsFit ==========")
    print("Type:", type(xpsFit))
    if isinstance(xpsFit, dict):
        print("Sample keys:", list(xpsFit.keys())[:5])
        for sample, regions in list(xpsFit.items())[:2]:
            print(f"\nSample: {sample!r}, type={type(regions)}")
            if isinstance(regions, dict):
                print("Subkeys:", list(regions.keys()))
                for region, data in list(regions.items())[:4]:
                    print(f"  {region!r}: type={type(data)}")
                    if isinstance(data, pd.DataFrame):
                        print("  Columns:", data.columns.tolist())
                        print(data.head(2).to_string(index=False))
                    else:
                        print(" ", repr(data)[:300])
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
        xpsFit0 = readXPSFitDat(folder)
        xpsFit0 = dop.copySamples(xpsFit0,sampleMap)
        xpsFit0 = dop.renameSamples(xpsFit0,sampleMap)
        debugXPSPeakData(peakinfor, xpsFit0)
        xpsFit=applyXPSPeakNames(xpsFit0, peakinfor)
        fig_path = fl.get_expanded_name(folder, fileName = "S1 O",type="png")
        plotXPSFit(
            xpsFit, sample="S1", region="O1s",figsize=(7, 1.5),
            # colors=["#A000E8", "#4BA3F5"],
            xlim=(527.5, 537), 
            xtick_step=2, show_legend=True,
            legend_loc="upper left",reverse_x=True,
            save_path=fig_path,
        )
        plotResults = plotAllXPSFits(
            xpsFit,
            folder=folder,
            figsize=(7, 1.5),
            show_legend=False,
            legend_loc="upper left",
            reverse_x=True,
        )
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
