##每一个子平台固定开头，用于找到依赖
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
import T0isothermsAndIAST.interpolateForDSL as intp
def select_folder():
    root = tk.Tk()
    root.withdraw()
    file_path = filedialog.askopenfilename(
        title="Select output file",
        filetypes=[("output files", "*.xlsx"), ("All files", "*.*")]
    )
    if not file_path:
        raise ValueError("No file selected!")
    print("Selected file:", file_path)
    return file_path
def readSimulationUptake(file):
    df = pd.read_excel(file, sheet_name="simulationUptake", header=None )
    # 找到 "Pore Size (nm)" 所在位置
    row_idx, col_idx = np.argwhere( df.values == "Pore Size (nm)" )[0]
    # 压力（标题行）
    pressures = ( df.iloc[row_idx, col_idx + 1:] .astype(float) .to_numpy() )
    # 孔径（第一列）
    poreSize = ( df.iloc[row_idx + 1:, col_idx] .astype(float) .to_numpy() )
    # 吸附数据 [cm^3 /g framework]
    uptake = ( df.iloc[row_idx + 1:, col_idx + 1:] .astype(float) .to_numpy() )
    return poreSize, pressures, uptake
def calculate_accessible_volume( helium_fraction, box_volume, framework_density):
    # Å³
    V_acc_A3 = helium_fraction * box_volume
    # cm³/g
    V_acc_cm3g = 1000.0 * helium_fraction / framework_density
    return V_acc_cm3g
def readSimulationDensity(file, HeliumFraction):
    # 不指定表头，全部读进来
    df = pd.read_excel( file, sheet_name="simulationDensity", header=None )
    # 找到 "Pore Size (nm)" 所在位置
    row_idx, col_idx = np.argwhere( df.values == "Pore Size (nm)" )[0]
    # 压力（标题行）
    pressures = ( df.iloc[row_idx, col_idx + 1:] .astype(float) .to_numpy() )
    # 孔径（第一列）
    poreSize = ( df.iloc[row_idx + 1:, col_idx] .astype(float) .to_numpy() )
    # 吸附数据 [cm^3 (STP)/cm^3 framework]
    density_cell = ( df.iloc[row_idx + 1:, col_idx + 1:] .astype(float) .to_numpy() )
    if HeliumFraction is not None:
        density = ( density_cell / HeliumFraction[:, None] ) # density_accessible [cm^3 (STP)/cm^3 framework]
    else:
        density = density_cell
    return poreSize, pressures, density
def readVolumeAndHeliumVoidFraction(file,sheet_name):
    # 不指定表头，全部读进来
    df = pd.read_excel( file, sheet_name, header=None )
    # 找到 "Pore Size (nm)" 所在位置
    row_idx, col_idx = np.argwhere( df.values == "Pore Size (nm)" )[0]
    # 孔径（第一列）
    poreSize = ( df.iloc[row_idx + 1:, col_idx] .astype(float) .to_numpy() )
    # HeliumVoidFraction
    HeliumFraction = ( df.iloc[row_idx + 1:, col_idx + 1] .astype(float) .to_numpy() )
    # BoxVolume (A^3)
    BoxVolume = ( df.iloc[row_idx + 1:, col_idx + 2] .astype(float) .to_numpy() )
    # Framework Density (kg/m^3)
    FrameworkDensity = ( df.iloc[row_idx + 1:, col_idx + 3] .astype(float) .to_numpy() )
    V_acc_cm3g = calculate_accessible_volume( HeliumFraction, BoxVolume, FrameworkDensity)
    return HeliumFraction, V_acc_cm3g
def readBoundary(file):
    df = pd.read_excel(file, sheet_name="boundary")
    poreSize = df.iloc[:, 0].to_numpy(dtype=float) / 10.0
    lower = df.iloc[:, 1].to_numpy(dtype=float) / 10.0
    upper = df.iloc[:, 2].to_numpy(dtype=float) / 10.0
    return poreSize, lower, upper
def findPSDHeader(df):
    for i in range(len(df)):
        row = df.iloc[i].astype(str).str.lower()
        if any("pore size" in item for item in row):
            return i
    raise ValueError("Cannot find PSD header row")
def mergePSD(sample):
    co2Pore = sample["co2Pore"]
    co2DV   = sample["co2DV"]
    n2Pore  = sample["n2Pore"]
    n2DV    = sample["n2DV"]
    maxPore = np.max(co2Pore)
    mask = n2Pore > maxPore
    pore = np.concatenate([co2Pore,n2Pore[mask]])
    dv = np.concatenate([co2DV,n2DV[mask]])
    order = np.argsort(pore)
    pore = pore[order]
    dv = dv[order]
    return pore,dv
def readExpUptake(file):
    df = pd.read_excel( file, sheet_name="exp", header=None )
    results = {}
    for col in range(df.shape[1]):
        sampleName = df.iloc[0, col]
        # 空单元格跳过
        if pd.isna(sampleName):
            continue
        sampleName = str(sampleName).strip()
        if sampleName == "":
            continue
        # 防止越界
        if col +1 >= df.shape[1]:
            break
        pressure = pd.to_numeric( df.iloc[2:, col], errors="coerce" ).dropna().to_numpy()
        expUptake = pd.to_numeric( df.iloc[2:, col+1], errors="coerce" ).dropna().to_numpy()
        results[sampleName] = {
            "pressure": pressure,
            "expUptake": expUptake
        }
    return results
def readPSD(file):
    df = pd.read_excel( file, sheet_name="PSD", header=None )
    results = {}
    for col in range(df.shape[1]):
        sampleName = df.iloc[0, col]
        # 空单元格跳过
        if pd.isna(sampleName):
            continue
        sampleName = str(sampleName).strip()
        if sampleName == "":
            continue
        # 防止越界
        if col + 3 >= df.shape[1]:
            break
        co2Pore = pd.to_numeric( df.iloc[3:, col], errors="coerce" ).dropna().to_numpy()
        co2DV = pd.to_numeric( df.iloc[3:, col+1], errors="coerce" ).dropna().to_numpy()
        n2Pore = pd.to_numeric( df.iloc[3:, col+2], errors="coerce" ).dropna().to_numpy()
        n2DV = pd.to_numeric( df.iloc[3:, col+3], errors="coerce" ).dropna().to_numpy()
        results[sampleName] = {
            "co2Pore": co2Pore,
            "co2DV": co2DV,
            "n2Pore": n2Pore,
            "n2DV": n2DV
        }
    return results
def calculatePSDVolumeAndWeight( psdPore, psdDV, boundary, ngrid=100): 
    interp = interp1d( psdPore, psdDV, bounds_error=False, fill_value=0.0 )
    volume = []
    for l, u in zip( boundary["lower"], boundary["upper"]):
        x = np.linspace(l, u, ngrid)
        area = np.trapezoid( interp(x), x )
        volume.append(area)
    volume = np.asarray(volume)
    total = np.sum(volume)
    if total > 0:
        weight = volume / total
    else:
        weight = np.zeros_like(volume)
    return weight, volume
def calculateIgnoredFraction(psdPoreOriginal,simPore):
    simMax = np.max(simPore)
    totalNum = len(psdPoreOriginal)
    exceedNum = np.sum(psdPoreOriginal > simMax)
    fraction = exceedNum / totalNum * 100
    print(f"Simulation max pore = {simMax:.3f} nm")
    print(f"PSD points outside range = {exceedNum}/{totalNum}")
    print(f"Outside fraction = {fraction:.2f}%")
    return fraction
def extendPSDToSimulationRange(psdPore,simPore):
    simMax = np.max(simPore)
    psdPore = psdPore.copy()
    num = np.sum(psdPore > simMax)
    if num > 0:
        print(f"{num} PSD points exceed simulation range")
        print(f"Assigning to {simMax:.3f} nm")
    psdPore[psdPore > simMax] = simMax
    return psdPore
def overlap(a1,a2,b1,b2):
    return max( 0.0, min(a2,b2)-max(a1,b1) )
def gcmcBoundary(pore):
    pore=np.sort(np.asarray(pore))
    edge=(pore[:-1]+pore[1:])/2
    lower=np.empty_like(pore)
    upper=np.empty_like(pore)
    lower[1:]=edge
    upper[:-1]=edge
    lower[0]=pore[0]-(pore[1]-pore[0])/2
    upper[-1]=pore[-1]+(pore[-1]-pore[-2])/2
    return {
        "lower": np.array(lower),
        "upper": np.array(upper)
    }
def buildBoundary(simPore,bdPore,bdLower,bdUpper):
    lower = []
    upper = []
    for pore in simPore:
        idx = np.where(np.isclose(pore,bdPore,atol=1e-6))[0]
        if len(idx) > 0:
            lower.append(bdLower[idx[0]])
            upper.append(bdUpper[idx[0]])
            continue
        pos = np.searchsorted(bdPore,pore)
        if pos == 0:
            l = pore - (bdPore[0] - pore) / 2
            u = (pore + bdPore[0]) / 2
        elif pos == len(bdPore):
            l = (bdPore[-1] + pore) / 2
            u = pore + (pore - bdPore[-1]) / 2
        else:
            left = bdPore[pos-1]
            right = bdPore[pos]
            l = (left + pore) / 2
            u = (pore + right) / 2
        lower.append(l)
        upper.append(u)
    return {
        "lower": np.array(lower),
        "upper": np.array(upper)
    }
# def calculateWeight( psdLower, psdUpper, psdVolume, simLower, simUpper):
#     weight=np.zeros(len(simLower))
#     for j in range(len(psdVolume)):
#         width=psdUpper[j]-psdLower[j]
#         for i in range(len(simLower)):
#             ov=overlap( psdLower[j], psdUpper[j], simLower[i], simUpper[i] )
#             if ov > 0:
#                 weight[i] += ( psdVolume[j] * ov / width )
#     return weight
def calculateUptakeByWeight(weight, adsorpUptake):
    contribution = weight[:, None] * adsorpUptake
    totalUptake = np.sum( contribution, axis=0 )   # mol/kg = mmol/g
    contributionPercent = np.divide( contribution, totalUptake[np.newaxis, :], 
                                    out=np.zeros_like(contribution), where=totalUptake[np.newaxis, :] != 0 )
    cumulative = np.cumsum( contributionPercent, axis=0 )
    return {
        "contribution": contribution,   # mol/kg
        "percent": contributionPercent,
        "cumulative": cumulative,
        "total": totalUptake           # mol/kg
    }
def calculateUptakeByDensity(volume, adsorpDensity):
    contribution = volume[:, None] * adsorpDensity
    totalUptakeSTP = np.sum(contribution, axis=0)  # cm3(STP)/g
    STP_MOLAR_VOLUME = 22.414  # cm3(STP)/mmol
    totalUptake = totalUptakeSTP / STP_MOLAR_VOLUME    # mmol/g
    contributionPercent = np.divide( contribution, totalUptakeSTP[np.newaxis, :], 
                                    out=np.zeros_like(contribution), where=totalUptakeSTP[np.newaxis, :] != 0 )
    cumulative = np.cumsum( contributionPercent, axis=0 )
    return {
        "contribution": contribution,      # cm3(STP)/g
        "percent": contributionPercent,
        "cumulative": cumulative,
        "totalUptakeSTP": totalUptakeSTP,       # cm3(STP)/g
        "total": totalUptake              # mmol/g
    }
def findThreshold(simPore, cumulative, contributionPercent, pressureIndex=-1):
    # 80/90/95%累计贡献孔径
    targets = [80, 90, 95]
    print("Cumulative contribution thresholds:")
    for target in targets:
        idx = np.where( cumulative[:, pressureIndex] * 100 >= target )[0]
        if len(idx) == 0:
            continue
        pore = simPore[idx[0]]
        print(f"{target}% contribution at pore size = {pore:.3f} nm")
    # 当前压力下贡献最大的三个孔径
    contrib = contributionPercent[:, pressureIndex]
    top_idx = np.argsort(contrib)[::-1][:3]
    print("\nTop 3 contributing pores:")
    for rank, i in enumerate(top_idx, start=1):
        print( f"{rank}. " f"{simPore[i]:.3f} nm " f"({contrib[i]*100:.2f}%)" )
def extractUptakeAtPressure(
        sampleResults,
        expData,
        pressure_kPa):

    records = []

    for sample in sampleResults:

        if sample not in expData:
            continue

        # ------------------------------------------------------
        # Experimental
        # ------------------------------------------------------

        _, _, exp_interp = dop.interpolateData(
            expData[sample]["pressure"],
            expData[sample]["expUptake"]
        )

        exp_uptake = float(
            exp_interp(pressure_kPa)
        )

        # ------------------------------------------------------
        # Simulation
        # ------------------------------------------------------

        _, _, sim_interp = dop.interpolateData(
            sampleResults[sample]["pressure"],
            sampleResults[sample]["uptake"]
        )

        sim_uptake = float(
            sim_interp(pressure_kPa)
        )

        if np.isnan(exp_uptake) or np.isnan(sim_uptake):
            continue

        records.append({
            "Sample": sample,
            "Pressure (kPa)": pressure_kPa,
            "Experimental": exp_uptake,
            "Simulation": sim_uptake,
        })

    return pd.DataFrame(records)
def analyzeNormalizedIsotherm(
        sampleResults,
        expData,
        pressure_kPa,
        plotflag=False):
    """
    Compare normalized uptake between simulation and experiment
    at a specified pressure.
    The normalized uptake is:
        q_norm = q / q_max
    where q_max is determined independently for each sample
    from the corresponding experimental/simulation isotherm.
    Parameters
    ----------
    sampleResults : dict
        {
            sample: {
                "pressure": ndarray,
                "uptake": ndarray,
            }
        }
    expData : dict
        Experimental isotherm data.
    pressure_kPa : float
        Target pressure.
    plotflag : bool
        Whether to plot the correlation.
    Returns
    -------
    df : pandas.DataFrame
        Normalized uptake comparison at the specified pressure.
    metrics : dict
        Correlation metrics across samples.
    """
    records = []
    for sample in sampleResults:
        if sample not in expData:
            continue
        # ==========================================================
        # 1. Experimental isotherm
        # ==========================================================
        exp_pressure = np.asarray(
            expData[sample]["pressure"]
        )
        exp_uptake = np.asarray(
            expData[sample]["expUptake"]
        )
        # ==========================================================
        # 2. Simulation isotherm
        # ==========================================================
        sim_pressure = np.asarray(
            sampleResults[sample]["pressure"]
        )
        sim_uptake = np.asarray(
            sampleResults[sample]["uptake"]
        )
        # Remove invalid values
        exp_mask = (
            np.isfinite(exp_pressure)
            & np.isfinite(exp_uptake)
        )
        sim_mask = (
            np.isfinite(sim_pressure)
            & np.isfinite(sim_uptake)
        )
        exp_pressure = exp_pressure[exp_mask]
        exp_uptake = exp_uptake[exp_mask]
        sim_pressure = sim_pressure[sim_mask]
        sim_uptake = sim_uptake[sim_mask]
        if len(exp_pressure) < 2 or len(sim_pressure) < 2:
            continue
        # ==========================================================
        # 3. Determine qmax independently
        # ==========================================================
        exp_qmax = np.max(exp_uptake)
        sim_qmax = np.max(sim_uptake)
        if exp_qmax <= 0 or sim_qmax <= 0:
            continue
        # ==========================================================
        # 4. Interpolate uptake at target pressure
        # ==========================================================
        if (
            pressure_kPa < np.min(exp_pressure)
            or pressure_kPa > np.max(exp_pressure)
            or pressure_kPa < np.min(sim_pressure)
            or pressure_kPa > np.max(sim_pressure)
        ):
            continue
        _, _, exp_interp = dop.interpolateData(
            exp_pressure,
            exp_uptake
        )
        _, _, sim_interp = dop.interpolateData(
            sim_pressure,
            sim_uptake
        )
        exp_uptake_p = float(
            exp_interp(pressure_kPa)
        )
        sim_uptake_p = float(
            sim_interp(pressure_kPa)
        )
        # ==========================================================
        # 5. Normalize
        # ==========================================================
        exp_norm = exp_uptake_p / exp_qmax
        sim_norm = sim_uptake_p / sim_qmax
        # ==========================================================
        # 6. Store
        # ==========================================================
        records.append({
            "Sample": sample,
            "Pressure (kPa)": pressure_kPa,
            "Experimental": exp_uptake_p,
            "Simulation": sim_uptake_p,
            "Experimental q/qmax": exp_norm,
            "Simulation q/qmax": sim_norm,
            "Difference": sim_norm - exp_norm,
            "ratio_exp/sim":
                exp_norm / sim_norm
                if sim_norm != 0 else np.nan,
        })
    # ==============================================================
    # Construct DataFrame
    # ==============================================================
    df = pd.DataFrame(records)
    min_samples = 3
    if len(df) < min_samples:
        print(
            f"Skip normalized isotherm at "
            f"{pressure_kPa:.2f} kPa: "
            f"only {len(df)} valid samples "
            f"(< {min_samples})."
        )
        return df, None
    # ==============================================================
    # Correlation
    # ==============================================================
    metrics = dop.correlationAnalysis(
        x=df["Experimental q/qmax"],
        y=df["Simulation q/qmax"],
        x_name="Experimental q/qmax",
        y_name="Simulation q/qmax",
    )
    if plotflag:
        myPlt.plotSingleCorrelation(
            metrics,
            xlabel="Experimental q/qmax",
            ylabel="Simulation q/qmax"
        )
    return df, metrics
def analyzeUptakeResponse(
        df,
        pressure_kPa,
        pair_mode="all", plotflag = False):
    """
    Analyze the differential uptake response between samples.
    Parameters
    ----------
    df : pandas.DataFrame
        Output DataFrame from analyzeUptakeAtPressure().
        Must contain:
            "Sample"
            "Experimental"
            "Simulation"
    pressure_kPa : float
        Pressure corresponding to the uptake data.
    pair_mode : str, optional
        "all"      : all possible sample pairs.
        "adjacent" : only adjacent samples in df.
    Returns
    -------
    df_delta : pandas.DataFrame
        Pairwise differential uptake data.
    metrics : dict or None
        Correlation analysis results for
        ΔSimulation vs ΔExperimental.
    """
    if df is None or len(df) < 2:
        return pd.DataFrame(), None
    records = []
    samples = df["Sample"].tolist()
    if pair_mode == "adjacent":
        pairs = zip(samples[:-1], samples[1:])
    elif pair_mode == "all":
        pairs = [
            (samples[i], samples[j])
            for i in range(len(samples))
            for j in range(i + 1, len(samples))
        ]
    else:
        raise ValueError(
            "pair_mode must be 'all' or 'adjacent'."
        )
    # Make lookup table
    df_lookup = df.set_index("Sample")
    for sample_i, sample_j in pairs:
        exp_i = df_lookup.loc[sample_i, "Experimental"]
        exp_j = df_lookup.loc[sample_j, "Experimental"]
        sim_i = df_lookup.loc[sample_i, "Simulation"]
        sim_j = df_lookup.loc[sample_j, "Simulation"]
        delta_exp = exp_j - exp_i
        delta_sim = sim_j - sim_i
        records.append({
            "Sample i": sample_i,
            "Sample j": sample_j,
            "Pressure (kPa)": pressure_kPa,
            "Experimental i": exp_i,
            "Experimental j": exp_j,
            "Simulation i": sim_i,
            "Simulation j": sim_j,
            "Delta Experimental": delta_exp,
            "Delta Simulation": delta_sim,
            "Delta Difference": delta_sim - delta_exp,
        })
    df_delta = pd.DataFrame(records)
    if len(df_delta) < 3:
        print(
            f"Skip Δresponse analysis at "
            f"{pressure_kPa:.2f} kPa: "
            f"only {len(df_delta)} valid pairs (< 3)."
        )
        return df_delta, None
    metrics = dop.correlationAnalysis(
        x=df_delta["Delta Experimental"],
        y=df_delta["Delta Simulation"],
        x_name="ΔExperimental uptake",
        y_name="ΔPSD-weighted uptake",
    )
    if plotflag:
        myPlt.plotSingleCorrelation( metrics, xlabel="ΔExperimental", ylabel="ΔSimulation", )
    return df_delta, metrics
def analyzeUptakeAtPressure(
        sampleResults,
        expData,
        pressure_kPa, plotflag = False):
    """
    Compare simulated and experimental uptake at a specified pressure.
    Parameters
    ----------
    sampleResults : dict
        {
            sample:{
                "pressure": ndarray,
                "uptake": ndarray,
            }
        }
    expData : dict
        Experimental isotherm data.
    pressure_kPa : float
        Target pressure.
    Returns
    -------
    df : pandas.DataFrame
    metrics : dict
        Correlation analysis results.
    """
    records = []
    for sample in sampleResults:
        # Experimental interpolation
        _, _, interp_func = dop.interpolateData(
            expData[sample]["pressure"],
            expData[sample]["expUptake"]
        )
        exp_uptake = float(interp_func(pressure_kPa))
        # Simulation interpolation (or direct value)
        _, _, sim_interp = dop.interpolateData(
            sampleResults[sample]["pressure"],
            sampleResults[sample]["uptake"]
        )
        sim_uptake = float(sim_interp(pressure_kPa))
        # Skip if either value is NaN
        if np.isnan(exp_uptake) or np.isnan(sim_uptake):
            continue
        # sim_uptake = sampleResults[sample]["uptake"][-1]
        records.append({
            "Sample": sample,
            "Pressure (kPa)": pressure_kPa,
            "Experimental": exp_uptake,
            "Simulation": sim_uptake,
            "Difference": sim_uptake - exp_uptake,
            "ratio_exp/sim": exp_uptake/sim_uptake,
        })
    df = pd.DataFrame(records)
    min_samples = 3
    if len(df) < min_samples:
        print(
            f"Skip {pressure_kPa:.2f} kPa: "
            f"only {len(df)} valid samples (< {min_samples})."
        )
        return df, None
    metrics = dop.correlationAnalysis(
        x=df["Experimental"],
        y=df["Simulation"],
        x_name="Experimental uptake",
        y_name="PSD-weighted uptake",
    )
    if plotflag:
        myPlt.plotSingleCorrelation(metrics,  xlabel="Experimental", ylabel="Simulation")
    return df, metrics
def analyzeSimulationValidation(
        sampleResults,
        expData,
        pressures,
        plotflag=False):
    """
    #Analyze simulation vs. experiment at all pressures.
    # ==========================================================
    # Validation relationships
    #
    # absolute:
    #     q_exp(P)  vs  q_sim(P)
    #
    # normalized:
    #     q_exp(P) / q_exp,max
    #     vs
    #     q_sim(P) / q_sim,max
    #
    # differential:
    #     Δq_exp(P)  vs  Δq_sim(P)
    # ==========================================================
    """
    compare_all = []
    metrics_all = []
    normalized_all = []
    normalized_metrics_all = []
    delta_all = []
    delta_metrics_all = []
    # ==========================================================
    # 1. Absolute uptake + differential response
    # ==========================================================
    for pressure in pressures / 1000:
        # ==========================================================
        # 1. Absolute uptake 每一个压强下：exp vs sim
        # ==========================================================
        df_compare, metrics = analyzeUptakeAtPressure(
            sampleResults=sampleResults,
            expData=expData,
            pressure_kPa=pressure,
            plotflag=plotflag
        )
        if not df_compare.empty:
            compare_all.append(df_compare)
        if metrics is not None:
            metrics_all.append({
                "Pressure (kPa)": pressure,
                "N": metrics["N"],
                "Pearson_r": metrics["Pearson_r"],
                "Pearson_p": metrics["Pearson_p"],
                "Spearman_r": metrics["Spearman_r"],
                "Spearman_p": metrics["Spearman_p"],
                "R2": metrics["R2"],
                "RMSE": metrics["RMSE"],
                "RMSE_rel": metrics["RMSE_rel"],
                "RMSE_nrm": metrics["RMSE_nrm"],
                "Slope": metrics["Slope"],
                "Intercept": metrics["Intercept"],
            })
        # ==========================================================
        # 2. Normalized uptake  每一个压强下：exp(P)/exp(MAX) vs sim(P)/sim(MAX)
        # ==========================================================
        df_normalized, normalized_metrics = \
            analyzeNormalizedIsotherm(
                sampleResults=sampleResults,
                expData=expData,
                pressure_kPa=pressure,
                plotflag=True
            )
        if not df_normalized.empty:
            normalized_all.append(df_normalized)
        if normalized_metrics is not None:
            normalized_metrics_all.append({
                "Pressure (kPa)": pressure,
                "N": normalized_metrics["N"],
                "Pearson_r": normalized_metrics["Pearson_r"],
                "Pearson_p": normalized_metrics["Pearson_p"],
                "Spearman_r": normalized_metrics["Spearman_r"],
                "Spearman_p": normalized_metrics["Spearman_p"],
                "R2": normalized_metrics["R2"],
                "RMSE": normalized_metrics["RMSE"],
                "RMSE_rel": normalized_metrics["RMSE_rel"],
                "RMSE_nrm": normalized_metrics["RMSE_nrm"],
                "Slope": normalized_metrics["Slope"],
                "Intercept": normalized_metrics["Intercept"],
            })
        # ==========================================================
        # 3. Differential uptake 每一个压强下，样品i,j之间 [exp(i)- exp(j)] vs [sim(i)- sim(j)]
        # ==========================================================
        if len(df_compare) >= 2:
            df_delta, delta_metrics = analyzeUptakeResponse(
                df=df_compare,
                pressure_kPa=pressure,
                pair_mode="all",
                plotflag=plotflag
            )
            if not df_delta.empty:
                delta_all.append(df_delta)
            if delta_metrics is not None:
                delta_metrics_all.append({
                    "Pressure (kPa)": pressure,
                    "N_pairs": delta_metrics["N"],
                    "Pearson_r": delta_metrics["Pearson_r"],
                    "Pearson_p": delta_metrics["Pearson_p"],
                    "Spearman_r": delta_metrics["Spearman_r"],
                    "Spearman_p": delta_metrics["Spearman_p"],
                    "R2": delta_metrics["R2"],
                    "RMSE": delta_metrics["RMSE"],
                    "RMSE_rel": delta_metrics["RMSE_rel"],
                    "RMSE_nrm": delta_metrics["RMSE_nrm"],
                    "Slope": delta_metrics["Slope"],
                    "Intercept": delta_metrics["Intercept"],
                })
    # ==========================================================
    # 2. Combine absolute results
    # ==========================================================
    if compare_all:
        df_absolute = pd.concat(
            compare_all,
            ignore_index=True
        )
    else:
        df_absolute = pd.DataFrame()
    df_absolute_metrics = pd.DataFrame(
        metrics_all
    )
    # ==========================================================
    # 3. Combine differential results
    # ==========================================================
    if delta_all:
        df_differential = pd.concat(
            delta_all,
            ignore_index=True
        )
    else:
        df_differential = pd.DataFrame()
    df_differential_metrics = pd.DataFrame(
        delta_metrics_all
    )
    # ==========================================================
    # 4. Combine normalized results
    # ==========================================================
    if normalized_all:
        df_normalized = pd.concat(
            normalized_all,
            ignore_index=True
        )
    else:
        df_normalized = pd.DataFrame()
    df_normalized_metrics = pd.DataFrame(
        normalized_metrics_all
    )
    # ==========================================================
    # Return validation results
    # ==========================================================
    validation = {
        "absolute": {
                "data": df_absolute,
                "metrics": df_absolute_metrics,
        },
        "normalized": {
                "data": df_normalized,
                "metrics": df_normalized_metrics,
        },
        "differential": {
                "data": df_differential,
                "metrics": df_differential_metrics,
        },
    }
    return validation
def plotContribution(simPore,contributionPercent,pressureIndex=-1):
    plt.figure(figsize=(6,4))
    plt.plot(simPore,contributionPercent[:,pressureIndex]*100,marker="o")
    plt.xlabel("Pore Size (nm)")
    plt.ylabel("Contribution (%)")
    plt.title("Pore Contribution")
    plt.tight_layout()
    plt.show()
def plotCumulative(simPore,cumulative,pressureIndex=-1):
    plt.figure(figsize=(6,4))
    plt.plot(simPore,cumulative[:,pressureIndex]*100,marker="o")
    plt.xlabel("Pore Size (nm)")
    plt.ylabel("Cumulative Contribution (%)")
    plt.title("Cumulative Contribution")
    plt.tight_layout()
    plt.show()
def plotPSDContribution(simPore,weight,contributionPercent,pressureIndex=-1):
    fig,ax1 = plt.subplots(figsize=(6,4))
    ax1.bar(simPore,weight,width=0.2)
    ax1.set_xlabel("Pore Size (nm)")
    ax1.set_ylabel("PSD Weight")
    ax2 = ax1.twinx()
    ax2.plot(simPore,contributionPercent[:,pressureIndex]*100,marker="o")
    ax2.set_ylabel("Contribution (%)")
    plt.title("PSD Weight vs Contribution")
    plt.tight_layout()
    plt.show()
def plotContributionHeatmap(simPore,pressures,contribution):
    plt.figure(figsize=(8,5))
    plt.imshow(contribution,aspect="auto",origin="lower")
    plt.colorbar(label="Contribution (mol/kg)")
    plt.xticks(np.arange(len(pressures)),[f"{int(p)}" for p in pressures],rotation=45)
    plt.yticks(np.arange(len(simPore)),[f"{p:.2f}" for p in simPore])
    plt.xlabel("Pressure (Pa)")
    plt.ylabel("Pore Size (nm)")
    plt.title("PSD-weighted Contribution")
    plt.tight_layout()
    plt.show()
def plotReconstructedIsotherm(pressures,totalUptake):
    plt.figure(figsize=(6,4))
    plt.plot(pressures,totalUptake,marker="o")
    plt.xscale("log")
    plt.xlabel("Pressure (Pa)")
    plt.ylabel("PSD-weighted Uptake (mol/kg)")
    plt.title("Reconstructed Isotherm")
    plt.tight_layout()
    plt.show()
def exportResult(simPore,boundary,pressures,weight,uptake,exp_interp_p, exp_interp_q,out_path):
    result = pd.DataFrame()
    result["Pore Size (nm)"] = simPore
    result["lower boundary (nm)"] = boundary["lower"]
    result["upper boundary (nm)"] = boundary["upper"]
    result["PSD Weight"] = weight
    contribution = uptake["contribution"]
    contributionPercent = uptake["percent"]
    cumulative = uptake["cumulative"]
    totalUptake = uptake["total"]
    result["Contribution"] = contribution[:,-1]
    result["Contribution (%)"] = contributionPercent[:,-1] * 100
    result["Cumulative (%)"] = cumulative[:,-1] * 100
    isotherm = pd.DataFrame({
    "Pressure (kPa)": pd.Series(pressures / 1000),
    "PSD-weighted Uptake (mol/kg)": pd.Series(totalUptake),
    "expPressure (kPa)": pd.Series(exp_interp_p),
    "expUptake (mol/kg)": pd.Series(exp_interp_q),
    })
    heatmap = pd.DataFrame(contribution)
    heatmap.index = simPore
    heatmap.columns = pressures
    with pd.ExcelWriter(out_path) as writer:
        result.to_excel(writer,sheet_name="Contribution",index=False)
        isotherm.to_excel(writer,sheet_name="Isotherm",index=False)
        heatmap.to_excel(writer,sheet_name="Heatmap")
def main(file_path=None):
    flagUsingDensity = False
    # ==== 输入参数 ====
    file_path = select_folder()   # 改成你的文件路径
    HeliumFraction, ModelVolumeAcc = readVolumeAndHeliumVoidFraction(file_path,sheet_name="volume")
    if flagUsingDensity:
        simPore,pressures,density_simu_acc = readSimulationDensity(file_path,HeliumFraction)
        simData = density_simu_acc
    else:
        simPore,pressures,uptake_simu_acc = readSimulationUptake(file_path)
        _,_,density_simu_acc = readSimulationDensity(file_path,HeliumFraction)
        simData = uptake_simu_acc
        simData_density = density_simu_acc
    boundary = gcmcBoundary(simPore)
    sampleAll = readPSD(file_path)
    expData =readExpUptake(file_path)
    sampleResults = {}
    for sample in sampleAll:
        # sample = "40-Bi-650-2-1"
        print(sample)
        psdPore,psdDV = mergePSD(sampleAll[sample])
        ignoredFraction = calculateIgnoredFraction(psdPore.copy(),simPore)
        psdPore = extendPSDToSimulationRange(psdPore,simPore)
        weight, volume = calculatePSDVolumeAndWeight(psdPore,psdDV,boundary)
        if flagUsingDensity:
            uptake = calculateUptakeByDensity(volume, simData)
        else:
            volumeWweight = volume/ModelVolumeAcc
            weight = volumeWweight
            uptake = calculateUptakeByWeight(weight, simData)
            uptake_density = calculateUptakeByDensity(volume, simData_density)
            # metrics = dop.calculateFitMetrics(xTrue=expData[sample]["pressure"], 
            #                                   yTrue=expData[sample]["expUptake"], xPred=pressures/1000, yPred=uptake["total"])
        sampleResults[sample] = {
            "pressure": pressures / 1000,
            "uptake": uptake["total"],
            }
        # print("ModelVolumeAcc= ", ModelVolumeAcc)
        # print("sum(volume)= ", np.sum(volume))
        # print("volume= ", volume)
        # print("volumeWweight= ", volume/ModelVolumeAcc)
        # print("weight= ", weight)
        # print("simPore= ",simPore)
        # print("HeliumFraction= ",HeliumFraction)
        # print("simData= ",np.max(simData,axis=1))
        findThreshold(simPore,uptake["cumulative"],uptake["percent"])
        # plotContribution(simPore,uptake["percent"])
        # plotCumulative(simPore,uptake["cumulative"])
        # plotPSDContribution(simPore,weight,uptake["percent"])
        # plotContributionHeatmap(simPore,pressures,uptake["contribution"])
        # plotReconstructedIsotherm(pressures,uptake["total"])
        # plt.ion()
        # myPlt.plotCurve( data={ sample: (expData[sample]["pressure"], expData[sample]["expUptake"]), }, fit={ sample: (pressures/1000, uptake["total"]), }, 
        #                             fit_label="PSD-weighted", marker=True, line=True, legendPosition="upper left")
        # if not flagUsingDensity:
        #     # plotCompareIsothermTwoMethod(sample, pressures, uptake, uptake_density)
        #     myPlt.plotCurve( data={ sample: (expData[sample]["pressure"], expData[sample]["expUptake"]), }, fit={ sample: (pressures/1000, uptake_density["total"]), }, 
        #                                         fit_label="PSD-weighted_density", marker=True, line=True, legendPosition="upper left")
        #     myPlt.plotCurve( data={ sample: (pressures/1000, uptake["total"]), }, fit={ sample: (pressures/1000, uptake_density["total"]), }, 
        #                                                     fit_label="PSD-weighted_density", marker=True, line=False, legendPosition="upper left")
            # plotCompareIsotherm(sample, pressures, uptake_density, expData, label_sim="PSD-weighted_density")
        out_path = fl.get_expanded_name(file_path, sample, expand="PSD_weighted", expandPos=True, type="xlsx")
        # exportResult(simPore,boundary,pressures,weight,uptake,exp_interp_p, exp_interp_q,out_path)
        # break
    validation = analyzeSimulationValidation( sampleResults=sampleResults, expData=expData, pressures=pressures, )
    out_path2 = fl.get_expanded_name(file_path, fileName="compare", expandPos=True, type="xlsx")
    # fl.export_to_excel_auto( df_compare_all, filename=out_path2, sheet_name="Comparison" )
    # fl.export_to_excel_auto( df_metrics, filename=out_path2, sheet_name="Correlation" )
    plt.show(block=True)
if __name__ == "__main__": 
    f = sys.argv[1] if len(sys.argv) > 1 else None
    main(f)