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
def mapSampleNames(data, column="Sample", reverse=False):
    sample_map = {
        "S1": "CC-Hy-600-0.5-1",
        "S2": "CC-Hy-600-1-1",
        "S3": "CC-Hy-700-1-1",
        "S4": "CC-Hy-800-0.5-1",
        "S5": "CC-Hy-800-1-1",
        "S6": "CC-Hy-600-2-1",
        "S7": "CC-Hy-900-1-1",
        "S8": "CC-Hy-800-2-1",
        "S9": "CC-Hy-600-4-1",
        "S10": "CC-Hy-800-4-1",
    }
    mapping = {v: k for k, v in sample_map.items()} if not reverse else sample_map

    if isinstance(data, pd.DataFrame):
        result = data.copy()
        result[column] = result[column].map(mapping).fillna(result[column])
        return result

    if isinstance(data, str):
        return mapping.get(data, data)

    if isinstance(data, (list, tuple)):
        return type(data)(mapping.get(x, x) for x in data)

    return data
def plotParameterCI(
    df, parameter, samples=None, figsize=(6, 5),
    xlabel=None, colors="black", savepath=None,
    xlim=None, showSE=False, fontsize=12,
    labelsize=None, ticksize=None
):
    parameter_col = "Quantity" if "Quantity" in df.columns else "Parameter"
    data = df[df[parameter_col] == parameter].copy()

    if samples is None:
        samples = data["Sample"].drop_duplicates().tolist()

    sample_order = {s: i for i, s in enumerate(samples)}

    data["_order"] = data["Sample"].map(sample_order)
    data = data.sort_values("_order").dropna( subset=["Estimate", "CI_Lower", "CI_Upper"] )
    data = data.drop(columns="_order")

    labelsize = labelsize or fontsize
    ticksize = ticksize or fontsize

    fig, ax = plt.subplots(figsize=figsize)

    y = np.arange(len(data))
    estimate = data["Estimate"].to_numpy()
    lower = data["CI_Lower"].to_numpy()
    upper = data["CI_Upper"].to_numpy()
    xerr = np.vstack([estimate - lower, upper - estimate])

    ax.errorbar(
        estimate, y, xerr=xerr,
        fmt="o", markersize=6, capsize=3,
        elinewidth=1.3, linewidth=1.3,
        color=colors, markeredgecolor="black",
        markeredgewidth=0.5
    )

    ax.set_yticks(y)
    ax.set_yticklabels(data["Sample"].astype(str))
    ax.invert_yaxis()

    ax.set_xlabel(xlabel or parameter, fontsize=labelsize)
    ax.set_ylabel("Sample", fontsize=labelsize)

    if xlim is not None:
        ax.set_xlim(xlim)
    if parameter.startswith(("b1_", "b2_")):
        ax.ticklabel_format( axis="x", style="sci", scilimits=(0, 0), useMathText=True )
        ax.xaxis.get_offset_text().set_fontsize(ticksize)

    ax.grid(axis="x", linestyle="--", alpha=0.35)
    ax.tick_params(
        axis="both", which="major",
        labelsize=ticksize,
        top=False, right=False
    )

    ax.spines["top"].set_visible(True)
    ax.spines["right"].set_visible(True)

    if showSE and "Bootstrap_SE" in data.columns:
        for yi, (_, row) in zip(y, data.iterrows()):
            ax.annotate(
                f'SE={row["Bootstrap_SE"]:.3g}',
                (row["CI_Upper"], yi),
                xytext=(5, 0), textcoords="offset points",
                va="center", fontsize=fontsize - 2
            )

    fig.tight_layout()

    if savepath:
        fig.savefig(savepath, dpi=600, bbox_inches="tight")

    return ax
def plotBootstrapParameterCI(
    bootstrap, model, parameter=None, samples=None,
    xlabel=None, figsize=(6, 5), saveFlag=False, folder=None,
    fontsize=12, labelsize=None, ticksize=None,
    showSE=False, xlim=None, colors="black", **kwargs
):
    model = model.upper()
    data = collectBootstrapResults(bootstrap, model=model)
    summary = mapSampleNames(data["Summary"])

    if summary.empty:
        raise ValueError(f"No bootstrap summary found for model: {model}")

    parameters = [parameter] if parameter else summary["Quantity"].dropna().unique()
    results = {}

    if saveFlag:
        folder = Path(folder) if folder else Path.cwd()
        folder.mkdir(parents=True, exist_ok=True)

    for p in parameters:
        df = summary[summary["Quantity"] == p].copy()
        if df.empty:
            print(f"Skip {p}: no data found.")
            continue

        savepath = folder / f"{model.replace('-', '')}_{p}.png" if saveFlag else None
        parameter_xlabel = xlabel or PARAMETER_LABELS.get(p, p)

        results[p] = plotParameterCI(
            df=df,
            parameter=p,
            samples=samples,
            xlabel=parameter_xlabel,
            figsize=figsize,
            colors=colors,
            savepath=savepath,
            xlim=xlim,
            showSE=showSE,
            fontsize=fontsize,
            labelsize=labelsize,
            ticksize=ticksize,
            **kwargs
        )

    return results
def readBootstrapResults(folder, pattern="*.xlsx", readDetails=False):
    results = {}
    folder = Path(folder)

    for file in folder.glob(pattern):
        if file.name.startswith("~$"):
            continue

        try:
            results[file.stem] = {
                "Bootstrap_CI_Summary": pd.read_excel(
                    file, sheet_name="Bootstrap_CI_Summary"
                )
            }

            if readDetails:
                results[file.stem]["Bootstrap_Distribution"] = pd.read_excel(
                    file, sheet_name="Bootstrap_Distribution"
                )

        except Exception as e:
            print(f"Failed to read {file.name}: {e}")

    return results
def collectBootstrapResults(bootstrap, model=None, readDetails=False):
    summary, distribution = [], []

    for result in bootstrap.values():
        ci = result["Bootstrap_CI_Summary"]

        if model is not None:
            ci = ci[ci["Model"].astype(str).str.upper() == model.upper()]

        summary.append(ci)

        if readDetails:
            if "Bootstrap_Distribution" not in result:
                raise ValueError(
                    "Bootstrap details are not loaded. "
                    "Set readDetails=True in readBootstrapResults()."
                )

            dist = result["Bootstrap_Distribution"]

            if model is not None:
                dist = dist[dist["Model"].astype(str).str.upper() == model.upper()]

            distribution.append(dist)

    return {
        "Summary": pd.concat(summary, ignore_index=True) if summary else pd.DataFrame(),
        "Distribution": (
            pd.concat(distribution, ignore_index=True)
            if distribution else pd.DataFrame()
        )
    }
PARAMETER_LABELS = {
    # D-A
    "W0_cm3_g": r"Maximum adsorbed volume, $W_0$ (cm$^3$ g$^{-1}$)",
    "E_J_mol": r"Adsorption energy, $E$ (J mol$^{-1}$)",
    "n": r"Heterogeneity parameter, $n$ (-)",

    # DSL
    "q1_mmol_g": r"Adsorption capacity, $q_A$ (mmol g$^{-1}$)",
    "q2_mmol_g": r"Adsorption capacity, $q_B$ (mmol g$^{-1}$)",

    # "b1_Pa_inv_T273.15K": r"Affinity constant, $b_A$ (Pa$^{-1}$) at 273.15 K",
    # "b2_Pa_inv_T273.15K": r"Affinity constant, $b_B$ (Pa$^{-1}$) at 273.15 K",

    # "b1_Pa_inv_T288.15K": r"Affinity constant, $b_A$ (Pa$^{-1}$) at 288.15 K",
    # "b2_Pa_inv_T288.15K": r"Affinity constant, $b_B$ (Pa$^{-1}$) at 288.15 K",

    # "b1_Pa_inv_T298.15K": r"Affinity constant, $b_A$ (Pa$^{-1}$) at 298.15 K",
    # "b2_Pa_inv_T298.15K": r"Affinity constant, $b_B$ (Pa$^{-1}$) at 298.15 K",
    "b1_Pa_inv_T273.15K": r"Affinity constant, $b_A$ (Pa$^{-1}$) at 0 °C",
    "b2_Pa_inv_T273.15K": r"Affinity constant, $b_B$ (Pa$^{-1}$) at 0 °C",

    "b1_Pa_inv_T288.15K": r"Affinity constant, $b_A$ (Pa$^{-1}$) at 15 °C",
    "b2_Pa_inv_T288.15K": r"Affinity constant, $b_B$ (Pa$^{-1}$) at 15 °C",

    "b1_Pa_inv_T298.15K": r"Affinity constant, $b_A$ (Pa$^{-1}$) at 25 °C",
    "b2_Pa_inv_T298.15K": r"Affinity constant, $b_B$ (Pa$^{-1}$) at 25 °C",
}
def main(file_path=None):
    singleFile = False
    fileName="SamplePairwise_CI"
    model="DSL"
    readDetails=False
    if singleFile:
        # 选择一个文件
        file = fl.getFile()
        sample = os.path.splitext(os.path.basename(file))[0]
        # df = fl.readData( file=file, keyFlag=keyFlag, keyName=keyName, ncols=ncols, skiprows=4,)
        # if not df.empty:
        #     charaData[sample] = df
    else:
        # 选择文件夹
        folder = fl.select_folder()
        out_path = fl.get_expanded_name(folder, fileName = fileName)
    if model == "D-A":
        bootstrap = readBootstrapResults(folder,readDetails=readDetails)
        DA_bootstrap = collectBootstrapResults(bootstrap, model="D-A")
        DA_CI = plotBootstrapParameterCI( bootstrap=bootstrap, model="D-A", samples=[f"S{i}" for i in range(1, 11)],
                    folder=folder, figsize=(6, 5), saveFlag=True, fontsize=12, labelsize=14, ticksize=11, showSE=False )
    elif model == "DSL":
        bootstrap = readBootstrapResults(folder, readDetails=readDetails)
        DA_bootstrap = collectBootstrapResults(bootstrap, model="DSL")
        DA_CI = plotBootstrapParameterCI( bootstrap=bootstrap, model="DSL", samples=[f"S{i}" for i in range(1, 11)],
                     folder=folder, figsize=(6, 4), saveFlag=True, fontsize=16, labelsize=14, ticksize=11, showSE=False )
    
    plt.show()
    y=0
if __name__ == "__main__":
    f = sys.argv[1] if len(sys.argv) > 1 else None
    main(f)