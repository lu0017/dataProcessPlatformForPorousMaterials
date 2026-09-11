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
import T1dataProcessSource.modelAndFit as mf
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
def debugPredictionResults(
        predictionResults,
        expDataCheck,
        printflag=True):
    """
    Debug prediction results by combining simulation,
    experimental, fitted parameters, and predicted uptake.
    Output columns
    --------------
    Sample
    Pressure (kPa)
    q_sim
    exp
    Slope
    Intercept
    q_pred
    Difference_pred_exp
    """
    debug_all = []
    for sample in predictionResults:
        result = predictionResults[sample]
        pressure = result["pressure"]
        q_sim = result["simUptake"]
        slope = result["slope"]
        intercept = result["intercept"]
        q_pred = result["predicted"]
        # --------------------------------------------------
        # Experimental data
        # --------------------------------------------------
        pressure_exp = expDataCheck[sample]["pressure"]
        uptake_exp = expDataCheck[sample]["expUptake"]
        # --------------------------------------------------
        # Build debug DataFrame
        # --------------------------------------------------
        df_sample = pd.DataFrame({
            "Sample": sample,
            "Slope": slope,
            "Intercept": intercept,
            "Pressure (kPa)": pressure,
            "q_sim": q_sim,
            "q_pred": q_pred,
        })
        # --------------------------------------------------
        # Match experimental uptake at pressure
        # --------------------------------------------------
        exp_interp = np.interp(
            pressure,
            pressure_exp,
            uptake_exp
        )
        df_sample["exp"] = exp_interp
        # Difference between predicted and experimental
        df_sample["Diff_q"] = (
            df_sample["q_pred"] - df_sample["exp"]
        )
        debug_all.append(df_sample)
    # ======================================================
    # Combine all samples
    # ======================================================
    df_debug = pd.concat(
        debug_all,
        ignore_index=True
    )
    # Reorder columns
    df_debug = df_debug[
        [
            "Sample",
            "Slope",
            "Intercept",
            "Pressure (kPa)",
            "q_sim",
            "exp",
            "q_pred",
            "Diff_q",
        ]
    ]
    # ======================================================
    # Print
    # ======================================================
    if printflag:
        for sample, df_sample in df_debug.groupby("Sample"):
            print("\n" + "=" * 100)
            print(f"Sample: {sample}")
            print("=" * 100)
            print(
                df_sample[
                    [
                        "Pressure (kPa)",
                        "Slope",
                        "Intercept",
                        "q_sim",
                        "exp",
                        "q_pred",
                        "Diff_q",
                    ]
                ].to_string(
                    index=False,
                    float_format=lambda x: f"{x:.6f}"
                )
            )
    return df_debug
def readSimulationUptake( file, pressure_min=None):
    # ======================================================
    # Read Excel
    # ======================================================
    df = pd.read_excel( file, sheet_name="simulationUptake", header=None )
    # ======================================================
    # Find "Pore Size (nm)"
    # ======================================================
    row_idx, col_idx = np.argwhere( df.values == "Pore Size (nm)" )[0]
    # ======================================================
    # Pressure
    # ======================================================
    pressures = ( df.iloc[ row_idx, col_idx + 1: ] .astype(float) .to_numpy() )
    # ======================================================
    # Pore size
    # ======================================================
    poreSize = ( df.iloc[ row_idx + 1:, col_idx ] .astype(float) .to_numpy() )
    # ======================================================
    # Uptake
    # ======================================================
    uptake = ( df.iloc[ row_idx + 1:, col_idx + 1: ] .astype(float) .to_numpy() )
    # ======================================================
    # Convert pressure
    # ======================================================
    pressures = pressures / 1000
    # ======================================================
    # Pressure filter
    # ======================================================
    if pressure_min is not None:
        mask = pressures >= pressure_min
        pressures = pressures[mask]
        # Uptake columns correspond to pressure
        uptake = uptake[:, mask]
    return ( poreSize, pressures, uptake )
def calculate_accessible_volume( helium_fraction, box_volume, framework_density):
    # Å³
    V_acc_A3 = helium_fraction * box_volume
    # cm³/g
    V_acc_cm3g = 1000.0 * helium_fraction / framework_density
    return V_acc_cm3g
def readSimulationDensity( file, HeliumFraction, pressure_min=None):
    # ======================================================
    # Read Excel
    # ======================================================
    df = pd.read_excel( file, sheet_name="simulationDensity", header=None )
    # ======================================================
    # Find "Pore Size (nm)"
    # ======================================================
    row_idx, col_idx = np.argwhere( df.values == "Pore Size (nm)" )[0]
    # ======================================================
    # Pressure
    # ======================================================
    pressures = ( df.iloc[ row_idx, col_idx + 1: ] .astype(float) .to_numpy() )
    # ======================================================
    # Pore size
    # ======================================================
    poreSize = ( df.iloc[ row_idx + 1:, col_idx ] .astype(float) .to_numpy() )
    # ======================================================
    # Density
    # ======================================================
    density_cell = ( df.iloc[ row_idx + 1:, col_idx + 1: ] .astype(float) .to_numpy() )
    # ======================================================
    # Helium correction
    # ======================================================
    if HeliumFraction is not None:
        density = ( density_cell / HeliumFraction[:, None] )
    else:
        density = density_cell
    # ======================================================
    # Pressure conversion
    # ======================================================
    pressures = pressures / 1000
    # ======================================================
    # Pressure filter
    # ======================================================
    if pressure_min is not None:
        mask = pressures >= pressure_min
        pressures = pressures[mask]
        # IMPORTANT:
        # density columns correspond to pressure
        density = density[:, mask]
    return ( poreSize, pressures, density )
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
def filterByMinimumPressure( data, pressure_min=None, pressure_col="Pressure (kPa)"):
    """
    Filter data by minimum pressure.
    Supports both:
        1. pandas.DataFrame containing a pressure column
        2. numpy.ndarray / list containing pressure values
    Parameters
    ----------
    data : pandas.DataFrame, numpy.ndarray, list
        Input data.
    pressure_min : float or None
        Minimum pressure. Values below this pressure are removed.
    pressure_col : str
        Pressure column name when data is a DataFrame.
    Returns
    -------
    filtered_data
        Filtered data with the same general data structure
        as the input.
    """
    if pressure_min is None:
        return data.copy()
    # ======================================================
    # DataFrame
    # ======================================================
    if isinstance(data, pd.DataFrame):
        if pressure_col not in data.columns:
            raise KeyError(
                f"Pressure column '{pressure_col}' not found."
            )
        return data[ data[pressure_col] >= pressure_min ].copy()
    # ======================================================
    # numpy array / list
    # ======================================================
    else:
        data = np.asarray(data)
        return data[data >= pressure_min]
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
def readExpUptake(file,sheet_name="exp"):
    df = pd.read_excel( file, sheet_name=sheet_name, header=None )
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
def readPSD(file,sheet_name="PSD"):
    df = pd.read_excel( file, sheet_name=sheet_name, header=None )
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
def preparePredictionAnalysisData( predictionResults, expData, modelfit=None, validation_type="absolute"):
    """
    Prepare prediction analysis data.
    The function:
        1. Keeps original experimental data.
        2. Keeps original prediction data.
        3. Builds a common pressure grid from both datasets.
        4. Interpolates experimental and predicted uptake onto
           the common pressure grid.
        5. Constructs point-wise prediction results.
        6. Collects fitting details if available.
    Returns
    -------
    predictionAnalysis : dict
        {
            "data_pred": DataFrame,
            "data_all": DataFrame,
            "fit_details": DataFrame
            "fit_metrics": DataFrame
        }
    """
    data_pred = []
    plot_data_all = []
    fit_details_all = []
    fit_metrics = []
    for sample in predictionResults:
        # ======================================================
        # 1. Experimental data
        # ======================================================
        pressure_exp = np.asarray( expData[sample]["pressure"], dtype=float )
        uptake_exp = np.asarray( expData[sample]["expUptake"], dtype=float )
        # Experimental interpolation
        _, _, exp_interp = mf.interpolateData( pressure_exp, uptake_exp, extrapolate=True )
        # ======================================================
        # 2. Prediction data
        # ======================================================
        pressure_pred = np.asarray( predictionResults[sample]["pressure"], dtype=float )
        simUptake_pred = np.asarray( predictionResults[sample]["simUptake"], dtype=float )
        uptake_pred = np.asarray( predictionResults[sample]["predicted"], dtype=float )
        # Prediction interpolation
        _, _, pred_interp = mf.interpolateData( pressure_pred, uptake_pred, extrapolate=True )
        metrics = mf.computeFitMetrics( fit_x=pressure_pred, fit_y=uptake_pred, exp_x=pressure_exp, exp_y=uptake_exp)
        fit_metrics.append({ "Sample": sample, **metrics, })
        # ======================================================
        # 3. Common pressure grid
        # ======================================================
        pressure_all = np.sort( np.unique( np.concatenate([ pressure_exp, pressure_pred ]) ) )
        # ======================================================
        # 4. Interpolated data
        # ======================================================
        uptake_exp_all = exp_interp( pressure_all )
        uptake_pred_all = pred_interp( pressure_all )
        # ======================================================
        # 5. Point-wise prediction results
        #
        #    Use prediction pressure as the original
        #    model calculation points.
        # ======================================================
        exp_at_pred = exp_interp( pressure_pred )
        difference = ( exp_at_pred - uptake_pred )
        ratio = np.divide( uptake_pred, exp_at_pred, 
                          out=np.full_like( uptake_pred, np.nan, dtype=float ), where=exp_at_pred != 0 )
        for p, sim_q, exp_q, pred_q, diff, rat in zip( pressure_pred, simUptake_pred, exp_at_pred, uptake_pred, difference, ratio):
            data_pred.append({
                "Sample": sample,
                "Validation type": validation_type,
                "Pressure (kPa)": p,
                "simUptake": sim_q,
                "Experimental": exp_q,
                "Predicted": pred_q,
                "Difference": diff,
                "Ratio": rat,
            })
        # ======================================================
        # 6. Plot data
        #    Dense common grid
        # ======================================================
        for p, exp_q, pred_q in zip(
                pressure_all,
                uptake_exp_all,
                uptake_pred_all):
            plot_data_all.append({
                "Sample": sample,
                "Pressure (kPa)": p,
                "Experimental": exp_q,
                "Predicted": pred_q,
            })
        # ======================================================
        # 7. Original prediction points
        #    Keep explicitly so they can be plotted as markers.
        # ======================================================
        for p, pred_q in zip( pressure_pred, uptake_pred):
            # Already contained in `data`,
            # but kept as metadata for clarity if needed.
            pass
        # ======================================================
        # 8. Fit details
        # ======================================================
        if modelfit is not None:
            fit_result = modelfit.get(sample, {})
            if isinstance(fit_result, dict):
                record = {
                    "Sample": sample,
                    "Validation type": validation_type,
                }
                record.update(fit_result)
                fit_details_all.append(record)
    # ==========================================================
    # 9. Construct DataFrames
    # ==========================================================
    df_data_pred = pd.DataFrame(data_pred)
    df_plot_data_all = pd.DataFrame(plot_data_all)
    df_fit_details = pd.DataFrame(fit_details_all)
    df_fit_metrics = pd.DataFrame(fit_metrics)
    return {
        "data_pred": df_data_pred,
        "data_all": df_plot_data_all,
        "fit_details": df_fit_details,
        "fit_metrics": df_fit_metrics,
    }
def calculatePSDVolume( psdPore, psdDV, boundary, ngrid=100): 
    interp = interp1d( psdPore, psdDV, bounds_error=False, fill_value=0.0 )
    volume = []
    for l, u in zip( boundary["lower"], boundary["upper"]):
        x = np.linspace(l, u, ngrid)
        area = np.trapezoid( interp(x), x )
        volume.append(area)
    volume = np.asarray(volume)
    return volume
def calculateWeightByVolume(volume, ModelVolumeAcc):
    return volume / ModelVolumeAcc
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
def extractUptakeAtPressure( sampleResults, expData, pressure_kPa, extrapolate=False):
    records = []
    for sample in sampleResults:
        if sample not in expData:
            continue
        # ------------------------------------------------------
        # Experimental
        # ------------------------------------------------------
        _, _, exp_interp = mf.interpolateData( expData[sample]["pressure"], 
                                              expData[sample]["expUptake"], extrapolate=extrapolate)
        exp_uptake = float( exp_interp(pressure_kPa) )
        # ------------------------------------------------------
        # Simulation
        # ------------------------------------------------------
        _, _, sim_interp = mf.interpolateData( sampleResults[sample]["pressure"], 
                                              sampleResults[sample]["simUptake"], extrapolate=extrapolate )
        sim_uptake = float( sim_interp(pressure_kPa) )
        if np.isnan(exp_uptake) or np.isnan(sim_uptake):
            continue
        records.append({
            "Sample": sample,
            "Pressure (kPa)": pressure_kPa,
            "q_exp": exp_uptake,
            "q_sim": sim_uptake,
        })
    return pd.DataFrame(records)
def extractUptakeAtAllPressures( sampleResults, expData, pressures):
    """
    Extract experimental and simulation uptake at all
    specified pressures.
    This function calls extractUptakeAtPressure() for each
    pressure and combines the results into one DataFrame.
    Parameters
    ----------
    sampleResults : dict
        Simulation isotherm data.
    expData : dict
        Experimental isotherm data.
    pressures : array-like
        Pressure points in Pa.
    Returns
    -------
    df : pandas.DataFrame
        Uptake data for all samples and pressures.
        Columns typically include:
            Sample
            Pressure (kPa)
            q_exp
            q_sim
            q_exp_max
            q_sim_max
            q_exp_norm
            q_sim_norm
    """
    all_data = []
    # ==========================================================
    # Extract uptake at each pressure
    # ==========================================================
    for pressure in pressures:      #  kPa
        df_pressure = extractUptakeAtPressure( sampleResults=sampleResults, expData=expData,
                                               pressure_kPa=pressure, extrapolate=True )
        if not df_pressure.empty:
            all_data.append(df_pressure)
    # ==========================================================
    # Combine all pressures
    # ==========================================================
    if all_data:
        df = pd.concat( all_data, ignore_index=True )
    else:
        df = pd.DataFrame()
    return df
def fitParameterBySingleModel(
        df_parameter,
        parameter_col,
        parameter_name,
        pressure_col="Pressure (kPa)",
        fit_method="nonlinear",
        model=None,
        p0=None,
        bounds=(-np.inf, np.inf),
        param_names=None,
        plotflag=False):
    """
    Fit the pressure dependence of a validation parameter
    using a single specified model.
    General form
    ------------
    Parameter(P) = f(P)
    Parameters
    ----------
    df_parameter : pandas.DataFrame
        DataFrame containing pressure and fitted parameter values.
    parameter_col : str
        Column containing the parameter to be fitted.
        For example:
            "Fit_slope"
            "Fit_intercept"
    parameter_name : str
        Name of the parameter for labeling and result storage.
    pressure_col : str
        Pressure column name.
    fit_method : str
        Fitting method used by correlationAnalysis().
        "linear" or "nonlinear".
    model : callable, optional
        Fitting model function.
        Required when fit_method="nonlinear".
    p0 : array-like, optional
        Initial parameters for nonlinear fitting.
    bounds : 2-tuple
        Parameter bounds.
    param_names : list, optional
        Names of model parameters.
    plotflag : bool
        Whether to plot the fitted pressure dependence.
    Returns
    -------
    fit_result : dict
        Complete fitting result returned by correlationAnalysis().
    """
    # ==========================================================
    # 1. Check required columns
    # ==========================================================
    required_columns = [ pressure_col, parameter_col ]
    missing_columns = [
        col for col in required_columns
        if col not in df_parameter.columns
    ]
    if missing_columns:
        raise KeyError(
            f"Missing columns for {parameter_name}: "
            f"{missing_columns}"
        )
    # ==========================================================
    # 2. Prepare fitting data
    # ==========================================================
    df_fit = df_parameter[ required_columns ].copy()
    df_fit = df_fit.replace( [np.inf, -np.inf], np.nan )
    df_fit = df_fit.dropna( subset=required_columns )
    if len(df_fit) < 2:
        return None
    # ==========================================================
    # 3. Perform fitting
    # ==========================================================
    fit_result = dop.correlationAnalysis(
        x=df_fit[pressure_col],
        y=df_fit[parameter_col],
        x_name=pressure_col,
        y_name=parameter_name,
        fit_method=fit_method,
        fit_func=model,
        p0=p0,
        bounds=bounds,
        param_names=param_names
    )
    # ==========================================================
    # 4. Plot
    # ==========================================================
    if plotflag:
        myPlt.plotSingleCorrelation( fit_result, xlabel=pressure_col, ylabel=parameter_name )
    return fit_result
def fitParameterByMultiModel(
        df_parameter,
        parameter_col,
        parameter_name,
        pressure_col="Pressure (kPa)",
        models=None,
        plotflag=False):
    """
    Fit the pressure dependence of a validation parameter
    using one or multiple candidate models.
    Model selection is based on AICc when multiple models
    are supplied.
    Parameters
    ----------
    df_parameter : pandas.DataFrame
        DataFrame containing pressure and fitted parameter values.
    parameter_col : str
        Column containing the parameter to be fitted.
    parameter_name : str
        Name of the parameter for labeling and result storage.
    pressure_col : str, default="Pressure (kPa)"
        Pressure column name.
    models : None, str, list, tuple, or dict, optional
        Model selection/configuration.
        None
            Use all default candidate models.
        str
            Use one fixed model, e.g.
                "Linear"
        list or tuple
            Compare selected models, e.g.
                ["Linear", "Quadratic"]
        dict
            Custom model configuration, e.g.
                {
                    "Linear": {
                        "fit_method": "linear",
                        "model": None,
                        "p0": None,
                        "bounds": (-np.inf, np.inf),
                        "param_names": None
                    }
                }
    plotflag : bool, default=False
        Whether to plot the best fitting model.
    Returns
    -------
    result : dict or None
        {
            "best_model": str,
            "best_fit": dict,
            "model_comparison": DataFrame,
            "fits": dict
        }
        Returns None if fitting is unsuccessful.
    """
    # ==========================================================
    # 1. Prepare fitting data
    # ==========================================================
    required_columns = [ pressure_col, parameter_col ]
    missing_columns = [ col for col in required_columns if col not in df_parameter.columns ]
    if missing_columns:
        raise KeyError(
            f"Missing columns for {parameter_name}: "
            f"{missing_columns}"
        )
    df_fit = df_parameter[ required_columns ].copy()
    df_fit = df_fit.replace( [np.inf, -np.inf], np.nan )
    df_fit = df_fit.dropna( subset=required_columns )
    # Sort by pressure
    df_fit = df_fit.sort_values( by=pressure_col ).reset_index(drop=True)
    if len(df_fit) < 3:
        print(
            f"Not enough data points to fit "
            f"{parameter_name}: N={len(df_fit)}"
        )
        return None
    # ==========================================================
    # 3. Resolve requested models
    # ==========================================================
    if models is None:
        # ------------------------------------------------------
        # Use all default candidate models
        # ------------------------------------------------------
        selected_models = mf.model_library.copy()
    elif isinstance(models, str):
        # ------------------------------------------------------
        # One fixed model
        # ------------------------------------------------------
        if models not in mf.model_library:
            raise ValueError(
                f"Unknown model: {models}. "
                f"Available models: "
                f"{list(mf.model_library.keys())}"
            )
        selected_models = { models: mf.model_library[models] }
    elif isinstance(models, (list, tuple)):
        # ------------------------------------------------------
        # Selected subset of candidate models
        # ------------------------------------------------------
        selected_models = {}
        for model_name in models:
            if model_name not in mf.model_library:
                raise ValueError(
                    f"Unknown model: {model_name}. "
                    f"Available models: "
                    f"{list(mf.model_library.keys())}"
                )
            selected_models[model_name] = ( mf.model_library[model_name] )
    elif isinstance(models, dict):
        # ------------------------------------------------------
        # Custom model configuration
        # ------------------------------------------------------
        selected_models = models
    else:
        raise TypeError(
            "`models` must be one of:\n"
            "  None\n"
            "  str\n"
            "  list / tuple\n"
            "  dict"
        )
    if not selected_models:
        raise ValueError(
            f"No models were specified for "
            f"{parameter_name}."
        )
    # ==========================================================
    # 4. Fit all selected models
    # ==========================================================
    fits = {}
    comparison_records = []
    n = len(df_fit)
    for model_name, config in selected_models.items():
        try:
            # --------------------------------------------------
            # Validate configuration
            # --------------------------------------------------
            required_config = [
                "fit_method",
                "model",
                "p0",
                "bounds",
                "param_names"
            ]
            missing_config = [ key for key in required_config if key not in config ]
            if missing_config:
                raise KeyError(
                    f"Missing model configuration for "
                    f"{model_name}: {missing_config}"
                )
            # --------------------------------------------------
            # Fit model
            # --------------------------------------------------
            fit_result = fitParameterBySingleModel(
                df_parameter=df_fit,
                parameter_col=parameter_col,
                parameter_name=parameter_name,
                pressure_col=pressure_col,
                fit_method=config["fit_method"],
                model=config.get("model"),
                p0=config.get("p0"),
                bounds=config.get( "bounds", (-np.inf, np.inf) ),
                param_names=config.get( "param_names" ),
                plotflag=False
            )
            if fit_result is None:
                continue
            # --------------------------------------------------
            # Convert fitting result to record
            # --------------------------------------------------
            record = mf.fittingResultToRecord( fit_result, parameter_name=parameter_name )
            # --------------------------------------------------
            # Number of fitted parameters
            # --------------------------------------------------
            if config["fit_method"] == "linear":
                k = 2
            else:
                param_names = config.get( "param_names" )
                if param_names is None:
                    raise ValueError(
                        f"`param_names` is required for "
                        f"nonlinear model: {model_name}"
                    )
                k = len(param_names)
            # --------------------------------------------------
            # Calculate SSE
            # --------------------------------------------------
            rmse = record.get( "RMSE", np.nan )
            sse = record.get( "SSE", np.nan )
            model_selection_metrics = mf.calculateModelSelectionMetrics( sse, n=n, k=k )
            fit_result.update({ "N": n, 
                               "N_parameters": k })
            fit_result.update(model_selection_metrics)
            # --------------------------------------------------
            # Store complete fitting result
            # --------------------------------------------------
            fits[model_name] = fit_result
            # ==========================================================
            # Store model comparison record
            # ==========================================================
            comparison_records.append({
                "Parameter": parameter_name,
                "Model": model_name,
                "N": fit_result.get("N", np.nan),
                "N_parameters": fit_result.get( "N_parameters", np.nan ),
                "R2": fit_result.get("R2", np.nan),
                "RMSE": fit_result.get("RMSE", np.nan),
                "SSE": fit_result.get("SSE", np.nan),
                "AIC": fit_result.get("AIC", np.nan),
                "AICc": fit_result.get("AICc", np.nan)
            })
        except Exception as e:
            # --------------------------------------------------
            # Keep failed model in comparison table
            # --------------------------------------------------
            comparison_records.append({
                "Parameter": parameter_name,
                "Model": model_name,
                "N": n,
                "N_parameters": np.nan,
                "R2": np.nan,
                "RMSE": np.nan,
                "AIC": np.nan,
                "AICc": np.nan,
                "Error": str(e)
            })
    # ==========================================================
    # 5. Build model comparison table
    # ==========================================================
    df_comparison = pd.DataFrame( comparison_records )
    if df_comparison.empty:
        print(
            f"No valid model was fitted for "
            f"{parameter_name}."
        )
        return None
    # ==========================================================
    # 6. Print model comparison
    # ==========================================================
    print( "\n================ MODEL COMPARISON ================" )
    print( df_comparison.to_string( index=False ) )
    print( "===================================================" )
    # ==========================================================
    # 7. Find valid models
    # ==========================================================
    valid = df_comparison[ df_comparison["AICc"].notna() ].copy()
    # Remove models whose fitting failed
    valid = valid[ valid["Model"].isin( fits.keys() ) ]
    if valid.empty:
        print( f"No valid fitted model for " f"{parameter_name}." )
        return None
    # ==========================================================
    # 8. Calculate Delta AICc
    # ==========================================================
    min_aicc = valid["AICc"].min()
    df_comparison["Delta_AICc"] = np.nan
    valid_indices = valid.index
    df_comparison.loc[ valid_indices, "Delta_AICc" ] = ( df_comparison.loc[ valid_indices, "AICc" ] - min_aicc )
    # ==========================================================
    # Add Delta AICc to each fitting result
    # ==========================================================
    for model_name, fit_result in fits.items():
        model_rows = df_comparison[ df_comparison["Model"] == model_name ]
        if not model_rows.empty:
            delta_aicc = model_rows.iloc[0]["Delta_AICc"]
            fit_result["Delta_AICc"] = delta_aicc
    # ==========================================================
    # 9. Select best model
    # ==========================================================
    best_idx = valid["AICc"].idxmin()
    best_model = df_comparison.loc[ best_idx, "Model" ]
    best_fit = fits[ best_model ]
    # ==========================================================
    # 10. Determine whether model was fixed
    # ==========================================================
    model_fixed = (
        isinstance(models, str)
        or ( isinstance(models, (list, tuple)) and len(models) == 1 )
    )
    # Add selection information
    df_comparison["Selected"] = False
    df_comparison.loc[ best_idx, "Selected" ] = True
    # ==========================================================
    # Add selection information to each fitting result
    # ==========================================================
    for model_name, fit_result in fits.items():
        fit_result["Selected"] = ( model_name == best_model )
    # ==========================================================
    # 11. Plot best model
    # ==========================================================
    if plotflag:
        myPlt.plotSingleCorrelation( best_fit, xlabel=pressure_col, ylabel=parameter_name, 
                                    # pos="bottom",xscale="log"
                                    )
    # ==========================================================
    # 12. Return
    # ==========================================================
    return {
        "best_model": best_model,
        "best_fit": best_fit,
        "model_comparison": df_comparison,
        "fits": fits,
        "model_fixed": model_fixed
    }
def fitValidationParameters(
        validation,
        validation_types=("absolute", "differential"),
        pressure_col="Pressure (kPa)",
        pressure_min=10,
        slope_models=None,
        intercept_models=None,
        plotflag=False):
    """
    Fit pressure dependence of validation Slope and Intercept.
    For each validation type:
        1. Fit Slope using multiple candidate models.
        2. Fit Intercept using multiple candidate models.
        3. Select the best model based on AICc.
    Parameters
    ----------
    validation : dict
        Output from validationSimAndExp().
    validation_types : tuple
        Validation types to analyze.
    pressure_col : str
        Pressure column name.
    pressure_min : float or None
        Minimum pressure used for parameter fitting.
    slope_models : dict, optional
        Candidate models for Slope(P).
    intercept_models : dict, optional
        Candidate models for Intercept(P).
    plotflag : bool
        Whether to plot the best fitting model.
    Returns
    -------
    parameter_fits : dict
    """
    # ==========================================================
    # 1. Prepare result container
    # ==========================================================
    parameter_fits = {}
    # ==========================================================
    # 2. Loop over validation types
    # ==========================================================
    for validation_type in validation_types:
        if validation_type not in validation:
            continue
        df_metrics = validation[ validation_type ]["metrics"]
        if df_metrics is None or df_metrics.empty:
            continue
        # ======================================================
        # 3. Required columns
        # ======================================================
        required_columns = [ pressure_col, "Fit_slope", "Fit_intercept" ]
        missing_columns = [
            col for col in required_columns
            if col not in df_metrics.columns
        ]
        if missing_columns:
            raise KeyError(
                f"Missing columns for {validation_type}: "
                f"{missing_columns}"
            )
        # ======================================================
        # 4. Prepare parameter data
        # ======================================================
        df_parameter = df_metrics[ required_columns ].copy()
        df_parameter = df_parameter.replace( [np.inf, -np.inf], np.nan )
        df_parameter = df_parameter.dropna( subset=required_columns )
        # ======================================================
        # 5. Pressure filter
        # ======================================================
        if pressure_min is not None:
            df_parameter = df_parameter[ df_parameter[pressure_col] >= pressure_min ].copy()
        if len(df_parameter) < 3:
            continue
        # ======================================================
        # 6. Fit Slope using multiple models
        # ======================================================
        slope_result = fitParameterByMultiModel(
            df_parameter=df_parameter,
            parameter_col="Fit_slope",
            parameter_name="Slope",
            pressure_col=pressure_col,
            models=slope_models,
            plotflag=plotflag
        )
        # ======================================================
        # 7. Fit Intercept using multiple models
        # ======================================================
        intercept_result = fitParameterByMultiModel(
            df_parameter=df_parameter,
            parameter_col="Fit_intercept",
            parameter_name="Intercept",
            pressure_col=pressure_col,
            models=intercept_models,
            plotflag=plotflag
        )
        if slope_result is None:
            continue
        if intercept_result is None:
            continue
        # Add best-fit Y values to parameter data
        df_parameter = df_parameter.copy()
        df_parameter["Slope_fit"] = slope_result["best_fit"]["Y_fit"]
        df_parameter["Intercept_fit"] = intercept_result["best_fit"]["Y_fit"]
        # ======================================================
        # 8. Build best-fit metrics
        # ======================================================
        slope_record = mf.fittingResultToRecord( slope_result["best_fit"], parameter_name="Slope" )
        intercept_record = mf.fittingResultToRecord( intercept_result["best_fit"], parameter_name="Intercept" )
        slope_record["Model"] = ( slope_result["best_model"] )
        intercept_record["Model"] = ( intercept_result["best_model"] )
        df_fit_metrics = pd.DataFrame([ slope_record, intercept_record ])
        # ======================================================
        # 9. Store results
        # ======================================================
        parameter_fits[validation_type] = {
            "data": df_parameter,
            "metrics": df_fit_metrics,
            "fits": {
                "Slope": slope_result["best_fit"],
                "Intercept": intercept_result["best_fit"]
            },
            "model_comparison": {
                "Slope": slope_result["model_comparison"],
                "Intercept": intercept_result["model_comparison"]
            },
            "all_fits": {
                "Slope": slope_result["fits"],
                "Intercept": intercept_result["fits"]
            }
        }
    return parameter_fits
def analyzeNormalizedIsotherm( df, df_uptake_all, plotflag=False):
    """
    Analyze normalized uptake at a specified pressure.
    Mathematical relationship
    --------------------------
        q_exp(P) / q_exp,max
        vs
        q_sim(P) / q_sim,max
    The normalization is performed independently for each sample.
    q_exp,max and q_sim,max are obtained from the complete isotherm
    of each sample.
    Parameters
    ----------
    df : pandas.DataFrame
        Uptake data at the current pressure.
        Required columns:
            "Sample"
            "q_exp"
            "q_sim"
    df_uptake_all : pandas.DataFrame
        Uptake data at all pressures.
        Required columns:
            "Sample"
            "q_exp"
            "q_sim"
    plotflag : bool
        Whether to plot the correlation.
    Returns
    -------
    df_valid : pandas.DataFrame
        Current-pressure data with normalized uptake columns added.
    metrics : dict or None
        Correlation analysis results.
    """
    df = df.copy()
    df_uptake_all = df_uptake_all.copy()
    # ==========================================================
    # 1. Calculate sample-specific qmax
    #
    # q_exp,max,i = max_P q_exp,i(P)
    # q_sim,max,i = max_P q_sim,i(P)
    # ==========================================================
    qmax = (
        df_uptake_all
        .groupby("Sample")
        .agg(
            q_exp_max=("q_exp", "max"),
            q_sim_max=("q_sim", "max")
        )
        .reset_index()
    )
    # ==========================================================
    # 2. Merge sample-specific qmax into current-pressure data
    # ==========================================================
    df = df.merge( qmax, on="Sample", how="left" )
    # ==========================================================
    # 3. Calculate normalized uptake
    # q_exp_norm = q_exp(P) / q_exp,max
    # q_sim_norm = q_sim(P) / q_sim,max
    # ==========================================================
    df["q_exp_norm"] = ( df["q_exp"] / df["q_exp_max"] )
    df["q_sim_norm"] = ( df["q_sim"] / df["q_sim_max"] )
    # ==========================================================
    # 4. Remove invalid values
    # ==========================================================
    valid = ( np.isfinite(df["q_exp_norm"]) & np.isfinite(df["q_sim_norm"]) )
    df_valid = df.loc[valid].copy()
    # ==========================================================
    # 5. Minimum sample requirement
    # ==========================================================
    min_samples = 3
    if len(df_valid) < min_samples:
        print(
            "Skip normalized analysis: "
            f"only {len(df_valid)} valid samples "
            f"(< {min_samples})."
        )
        return df_valid, None
    n_unique_exp = df_valid["q_exp_norm"].nunique()
    n_unique_sim = df_valid["q_sim_norm"].nunique()
    if n_unique_exp < 2 or n_unique_sim < 2:
        print( "Skip normalized analysis: normalized uptake has no variation at this pressure." )
        return df_valid, None
    # ==========================================================
    # 7. Correlation
    # ==========================================================
    metrics = dop.correlationAnalysis(
        x=df_valid["q_sim_norm"],
        y=df_valid["q_exp_norm"],
        x_name="Simulation q/qmax",
        y_name="Experimental q/qmax",
    )
    # ==========================================================
    # 8. Plot
    # ==========================================================
    if plotflag:
        myPlt.plotSingleCorrelation(
            metrics,
            xlabel=r"Simulation $q/q_{\max}$",
            ylabel=r"Experimental $q/q_{\max}$"
        )
    return df_valid, metrics
def analyzeUptakeResponse( df, pressure_kPa, pair_mode="all", plotflag=False):
    """
    Analyze pairwise differential uptake response between samples.
    Mathematical relationship
    --------------------------
        Δq_exp(P)
        =
        q_exp,i(P) - q_exp,j(P)
        vs
        Δq_sim(P)
        =
        q_sim,i(P) - q_sim,j(P)
    Parameters
    ----------
    df : pandas.DataFrame
        Output from extractUptakeAtPressure().
        Required columns:
            "Sample"
            "q_exp"
            "q_sim"
    pressure_kPa : float
        Current pressure.
    pair_mode : str
        "all" :
            Calculate all unique sample pairs.
    plotflag : bool
        Whether to plot the correlation.
    Returns
    -------
    df_delta : pandas.DataFrame
        Pairwise differential uptake.
    metrics : dict or None
        Correlation analysis results.
    """
    samples = df["Sample"].tolist()
    records = []
    # ==========================================================
    # Generate sample pairs
    # ==========================================================
    if pair_mode == "all":
        for i in range(len(samples)):
            for j in range(i + 1, len(samples)):
                sample_i = samples[i]
                sample_j = samples[j]
                row_i = df.loc[ df["Sample"] == sample_i ].iloc[0]
                row_j = df.loc[ df["Sample"] == sample_j ].iloc[0]
                # --------------------------------------------------
                # Differential experimental uptake
                # --------------------------------------------------
                delta_exp = ( row_i["q_exp"] - row_j["q_exp"] )
                # --------------------------------------------------
                # Differential simulation uptake
                # --------------------------------------------------
                delta_sim = ( row_i["q_sim"] - row_j["q_sim"] )
                records.append({
                    "Pressure (kPa)": pressure_kPa,
                    "Sample_i": sample_i,
                    "Sample_j": sample_j,
                    "q_exp_i": row_i["q_exp"],
                    "q_exp_j": row_j["q_exp"],
                    "q_sim_i": row_i["q_sim"],
                    "q_sim_j": row_j["q_sim"],
                    "Delta_q_exp": delta_exp,
                    "Delta_q_sim": delta_sim,
                })
    else:
        raise ValueError( f"Unsupported pair_mode: {pair_mode}" )
    # ==========================================================
    # Create DataFrame
    # ==========================================================
    df_delta = pd.DataFrame(records)
    if df_delta.empty:
        return df_delta, None
    # ==========================================================
    # Remove invalid values
    # ==========================================================
    valid = ( np.isfinite(df_delta["Delta_q_exp"]) & np.isfinite(df_delta["Delta_q_sim"]) )
    df_delta = df_delta.loc[valid].copy()
    # ==========================================================
    # Minimum pair requirement
    # ==========================================================
    min_pairs = 3
    if len(df_delta) < min_pairs:
        print(
            f"Skip differential analysis at "
            f"{pressure_kPa:.2f} kPa: "
            f"only {len(df_delta)} valid pairs "
            f"(< {min_pairs})."
        )
        return df_delta, None
    # ==========================================================
    # Correlation
    # ==========================================================
    metrics = dop.correlationAnalysis(
        x=df_delta["Delta_q_sim"],
        y=df_delta["Delta_q_exp"],
        x_name="Simulation Δq",
        y_name="Experimental Δq",
    )
    # ==========================================================
    # Plot
    # ==========================================================
    if plotflag:
        myPlt.plotSingleCorrelation( metrics, xlabel=r"Simulation $\Delta q$", ylabel=r"Experimental $\Delta q$" )
    return df_delta, metrics
def analyzeUptakeAtPressure( df, plotflag=False):
    if len(df) < 3:
        return df, None
    metrics = dop.correlationAnalysis( x=df["q_sim"], y=df["q_exp"], x_name="Simulation uptake", 
                                      y_name="Experimental uptake", )
    if plotflag:
        myPlt.plotSingleCorrelation( metrics, xlabel="Simulation", ylabel="Experimental" )
    return df, metrics
def validationSimAndExp( sampleResults, expData, pressures, plotflag=False):
    """
    Analyze simulation vs. experiment at all pressures.
    Validation relationships
    ------------------------
    absolute:
        q_exp(P) vs q_sim(P)
    normalized:
        q_exp(P) / q_exp,max
        vs
        q_sim(P) / q_sim,max
    differential:
        Δq_exp(P) vs Δq_sim(P)
    Returns
    -------
    validation : dict
        {   "absolute": {
                "data": DataFrame,
                "metrics": DataFrame,
            },
            "normalized": {
                "data": DataFrame,
                "metrics": DataFrame,
            },
            "differential": {
                "data": DataFrame,
                "metrics": DataFrame,
            },
        }
    """
    # ==========================================================
    # 1. Extract uptake at all pressures
    # ==========================================================
    df_uptake_all = extractUptakeAtAllPressures( sampleResults=sampleResults, expData=expData, pressures=pressures )
    # ==========================================================
    # 2. Prepare result containers
    # ==========================================================
    results = {
        "absolute": {
            "data": [],
            "metrics": []
        },
        "normalized": {
            "data": [],
            "metrics": []
        },
        "differential": {
            "data": [],
            "metrics": []
        }
    }
    # ==========================================================
    # 3. Analyze each pressure
    # ==========================================================
    for pressure, df_uptake in df_uptake_all.groupby( "Pressure (kPa)" ):
        # ======================================================
        # 3.1 Absolute uptake
        # q_exp(P) vs q_sim(P)
        # ======================================================
        ( df_absolute, absolute_metrics ) = analyzeUptakeAtPressure( df=df_uptake, plotflag=plotflag )
        if not df_absolute.empty:
            results["absolute"]["data"].append( df_absolute )
        if absolute_metrics is not None:
            results["absolute"]["metrics"].append(
                mf.fittingResultToRecord( absolute_metrics, extra={ "Pressure (kPa)": pressure } )
            )
        # ======================================================
        # 3.2 Normalized uptake
        # q_exp(P) / q_exp,max
        # vs
        # q_sim(P) / q_sim,max
        # ======================================================
        ( df_normalized, normalized_metrics ) = analyzeNormalizedIsotherm( df=df_uptake, 
                                                                          df_uptake_all=df_uptake_all, plotflag=plotflag )
        if not df_normalized.empty:
            results["normalized"]["data"].append( df_normalized )
        if normalized_metrics is not None:
            results["normalized"]["metrics"].append(
                mf.fittingResultToRecord( normalized_metrics, extra={ "Pressure (kPa)": pressure } )
            )
        # ======================================================
        # 3.3 Differential uptake
        # Δq_exp(P) vs Δq_sim(P)
        # ======================================================
        if len(df_uptake) >= 2:
            ( df_differential, differential_metrics ) = analyzeUptakeResponse( df=df_uptake, 
                                                                              pressure_kPa=pressure, pair_mode="all", 
                                                                              plotflag=plotflag )
            if not df_differential.empty:
                results["differential"]["data"].append( df_differential )
            if differential_metrics is not None:
                results["differential"]["metrics"].append(
                    mf.fittingResultToRecord(differential_metrics, extra={"Pressure (kPa)": pressure}, n_name="N_pairs")
                )
    # ==========================================================
    # 4. Combine results
    # ==========================================================
    validation = {}
    for validation_type, contents in results.items():
        # ------------------------------------------------------
        # Data
        # ------------------------------------------------------
        if contents["data"]:
            df_data = pd.concat( contents["data"], ignore_index=True )
        else:
            df_data = pd.DataFrame()
        # ------------------------------------------------------
        # Metrics
        # ------------------------------------------------------
        df_metrics = pd.DataFrame( contents["metrics"] )
        validation[validation_type] = { "data": df_data, "metrics": df_metrics }
    return validation
def modifyPoreVolume( volume, target_index, delta_volume, source_index=None):
    """
    Modify pore volume for perturbation analysis.
    Parameters
    ----------
    volume : array-like
        Original pore volume mapped to fixed GCMC pore windows.
    target_index : int
        Index of the pore window receiving additional volume.
    delta_volume : float
        Volume change.
        Addition:
            source_index = None
            V_target' = V_target + delta_volume
        Redistribution:
            source_index is not None
            V_source' = V_source - delta_volume
            V_target' = V_target + delta_volume
    source_index : int, optional
        Index of the pore window from which volume is removed.
        If None, this is a volume-addition perturbation.
    Returns
    -------
    modified_volume : np.ndarray
        Modified pore-volume distribution.
    """
    volume = np.asarray(volume, dtype=float)
    if delta_volume <= 0:
        raise ValueError( "delta_volume must be greater than zero." )
    if not (0 <= target_index < len(volume)):
        raise IndexError( f"target_index={target_index} is out of range." )
    modified_volume = volume.copy()
    # ======================================================
    # Volume addition
    # ======================================================
    if source_index is None:
        modified_volume[target_index] += delta_volume
    # ======================================================
    # Volume redistribution
    # ======================================================
    else:
        if not (0 <= source_index < len(volume)):
            raise IndexError( f"source_index={source_index} is out of range." )
        if source_index == target_index:
            raise ValueError( "source_index and target_index must be different." )
        if volume[source_index] < delta_volume:
            raise ValueError(
                f"Insufficient pore volume at source_index="
                f"{source_index}: "
                f"{volume[source_index]:.6g} < "
                f"{delta_volume:.6g}"
            )
        modified_volume[source_index] -= delta_volume
        modified_volume[target_index] += delta_volume
    # ======================================================
    # Numerical safety check
    # ======================================================
    if np.any(modified_volume < -1e-12):
        raise ValueError( "Modified pore volume contains negative values." )
    modified_volume[ np.abs(modified_volume) < 1e-12 ] = 0.0
    return modified_volume
def calculateSimByAddition( volume, target_index, delta_volume, simPore, simData, ModelVolumeAcc, flagUsingDensity=False, 
                           simData_density=None):
    """
    Calculate PSD-weighted GCMC uptake after adding pore volume
    to one fixed GCMC pore-size window.
    """
    # ------------------------------------------------------
    # 1. Modify pore volume
    # ------------------------------------------------------
    modified_volume = modifyPoreVolume( volume=volume, target_index=target_index, delta_volume=delta_volume )
    # ------------------------------------------------------
    # 2. Recalculate weighted simulation
    # ------------------------------------------------------
    result = calculateWeightedSimulationFromVolume( volume=modified_volume, simData=simData, 
                                                   ModelVolumeAcc=ModelVolumeAcc, flagUsingDensity=flagUsingDensity, 
                                                   simData_density=simData_density )
    # ------------------------------------------------------
    # 3. Store perturbation information
    # ------------------------------------------------------
    result["original_volume"] = np.asarray( volume, dtype=float ).copy()
    result["modified_volume"] = modified_volume
    result["target_index"] = target_index
    result["target_pore"] = simPore[target_index]
    result["delta_volume"] = delta_volume
    result["operation"] = "addition"
    return result
def plotContribution(simPore,contributionPercent,pressureIndex=-1):
    plt.figure(figsize=(6,4))
    plt.plot(simPore,contributionPercent[:,pressureIndex]*100,marker="o")
    plt.xlabel("Pore Size (nm)")
    plt.ylabel("Contribution (%)")
    plt.title("Pore Contribution")
    plt.tight_layout()
    plt.show(block=False)
def plotCumulative(simPore,cumulative,pressureIndex=-1):
    plt.figure(figsize=(6,4))
    plt.plot(simPore,cumulative[:,pressureIndex]*100,marker="o")
    plt.xlabel("Pore Size (nm)")
    plt.ylabel("Cumulative Contribution (%)")
    plt.title("Cumulative Contribution")
    plt.tight_layout()
    plt.show(block=False)
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
    plt.show(block=False)
def plotContributionHeatmap(simPore,pressures,contribution):
    plt.figure(figsize=(8,5))
    plt.imshow(contribution,aspect="auto",origin="lower")
    plt.colorbar(label="Contribution (mol/kg)")
    plt.xticks(np.arange(len(pressures)),[f"{int(p)}" for p in pressures],rotation=45)
    plt.yticks(np.arange(len(simPore)),[f"{p:.2f}" for p in simPore])
    plt.xlabel("Pressure (kPa)")
    plt.ylabel("Pore Size (nm)")
    plt.title("PSD-weighted Contribution")
    plt.tight_layout()
    plt.show(block=False)
def plotReconstructedIsotherm(pressures,totalUptake):
    plt.figure(figsize=(6,4))
    plt.plot(pressures,totalUptake,marker="o")
    plt.xscale("log")
    plt.xlabel("Pressure (kPa)")
    plt.ylabel("PSD-weighted Uptake (mol/kg)")
    plt.title("Reconstructed Isotherm")
    plt.tight_layout()
    plt.show(block=False)
def plotPredictionAnalysisByPressure( predictionAnalysis, sample=None, 
                           pressure_type="merge", fit_label="exp=f(sim)", legendPosition="upper left"):
    """
    Plot prediction analysis results.
    Parameters
    ----------
    predictionAnalysis : dict
        Output from preparePredictionAnalysisData().
    sample : str or None
        Sample to plot. If None, plot all samples.
    pressure_type : {"merge", "sim"}
        Pressure grid used for plotting.
        "merge":
            Use the common pressure grid from data_all.
        "sim":
            Use the original GCMC simulation pressure grid
            from data_pred.
    fit_label : str
        Label for prediction curve.
    legendPosition : str
        Position of the legend.
    """
        # ======================================================
    # 1. Select pressure data
    # ======================================================
    if pressure_type == "merge":
        df_plot = predictionAnalysis["data_all"]
        line = True
    elif pressure_type == "sim":
        df_plot = predictionAnalysis["data_pred"]
        line = False
    else:
        raise ValueError(
            "pressure_type must be either 'merge' or 'sim'."
        )
    if sample is None:
        samples = df_plot["Sample"].unique()
    else:
        samples = [sample]
    for sample_name in samples:
        # ======================================================
        # Dense interpolated curves
        # ======================================================
        df_sample = df_plot[ df_plot["Sample"] == sample_name ]
        myPlt.plotCurve(
            data={ sample_name: ( df_sample["Pressure (kPa)"].to_numpy(), df_sample["Experimental"].to_numpy() ) },
            fit={ sample_name: ( df_sample["Pressure (kPa)"].to_numpy(), df_sample["Predicted"].to_numpy() ) },
            fit_label=fit_label,
            marker=True,
            line=line,
            legendPosition=legendPosition
        )
    plt.show(block=False)
def plotPerturbationAddition(
        sample,
        perturbationResults,
        poreResult):
    """
    Plot original and modified results for one pore-volume perturbation.
    Three plots are generated:
        1. Pore-volume distribution (PSD)
        2. PSD-weighted simulation uptake (SIM)
        3. Predicted experimental uptake (PRE)
    Parameters
    ----------
    sample : str
        Sample name.
    perturbationResults : dict
        Complete perturbation result generated by
        calculatePerturbationByAddition().
    poreResult : dict
        Perturbation result for one pore-size window.
    Returns
    -------
    None
    """
    # ======================================================
    # 1. Extract original / baseline information
    # ======================================================
    simPore = np.asarray(
        perturbationResults["simPore"],
        dtype=float
    )
    pressures = np.asarray(
        perturbationResults["pressure"],
        dtype=float
    )
    original_volume = np.asarray(
        perturbationResults["original_volume"],
        dtype=float
    )
    original_sim = np.asarray(
        perturbationResults["original_simUptake"],
        dtype=float
    )
    original_pred = np.asarray(
        perturbationResults["original_predicted"],
        dtype=float
    )
    # ======================================================
    # 2. Extract perturbation information
    # ======================================================
    target_index = poreResult["target_index"]
    target_pore = poreResult["target_pore"]
    delta_volume = poreResult["delta_volume"]
    modified_volume = np.asarray(
        poreResult["modified_volume"],
        dtype=float
    )
    modified_sim = np.asarray(
        poreResult["modified_simUptake"],
        dtype=float
    )
    modified_pred = np.asarray(
        poreResult["modified_predicted"],
        dtype=float
    )
    # ======================================================
    # 3. PSD / pore-volume distribution
    # ======================================================
    myPlt.plotCurve(
        data={
            "Original": (
                simPore,
                original_volume
            ),
        },
        fit={
            "Original": (
                simPore,
                modified_volume
            ),
        },
        fit_label="Modified",
        x="Pore size",
        y="Volume",
        xlabel="Pore size (nm)",
        ylabel="Pore volume",
        title=(
            f"{sample} - PSD - "
            f"pore {target_pore:.2f} nm"
        ),
        marker=True,
        line=True,
        legendPosition="upper left",
    )
    # ======================================================
    # 4. PSD-weighted simulation uptake
    # ======================================================
    myPlt.plotCurve(
        data={
            "Original": (
                pressures,
                original_sim
            ),
        },
        fit={
            "Original": (
                pressures,
                modified_sim
            ),
        },
        fit_label="Modified",
        x="Pressure",
        y="Uptake",
        xlabel="Pressure (kPa)",
        ylabel="CO$_2$ uptake",
        title=(
            f"{sample} - SIM - "
            f"pore {target_pore:.2f} nm"
        ),
        marker=True,
        line=True,
        legendPosition="upper left",
    )
    # ======================================================
    # 5. Predicted experimental uptake
    # ======================================================
    myPlt.plotCurve(
        data={
            "Original": (
                pressures,
                original_pred
            ),
        },
        fit={
            "Original": (
                pressures,
                modified_pred
            ),
        },
        fit_label="Modified",
        x="Pressure",
        y="Uptake",
        xlabel="Pressure (kPa)",
        ylabel="Predicted CO$_2$ uptake",
        title=(
            f"{sample} - PRE - "
            f"pore {target_pore:.2f} nm"
        ),
        marker=True,
        line=True,
        legendPosition="upper left",
    )
    plt.show(block=False)
    # ======================================================
    # 6. Debug information
    # ======================================================
    print(
        f"[Perturbation Plot] {sample} | "
        f"index={target_index} | "
        f"pore={target_pore:.4f} nm | "
        f"deltaV={delta_volume}"
    )
def plotSensitivityHeatmap(
        perturbationResults,
        sensitivity_type="pred"):
    """
    Plot perturbation sensitivity heatmap.
    X-axis : Pressure (kPa)
    Y-axis : Pore size (nm)
    Color  : Sensitivity
    """
    if sensitivity_type == "sim":
        sensitivity_key = "sensitivity_sim"
        title = "PSD-weighted Simulation Sensitivity"
    elif sensitivity_type == "pred":
        sensitivity_key = "sensitivity_pred"
        title = "Predicted Experimental Sensitivity"
    else:
        raise ValueError(
            "sensitivity_type must be 'sim' or 'pred'."
        )
    pressures = np.asarray(
        perturbationResults["pressure"],
        dtype=float
    )
    poreResults = perturbationResults["pore_results"]
    simPore = np.array([
        result["target_pore"]
        for result in poreResults
    ])
    sensitivity = np.array([
        result[sensitivity_key]
        for result in poreResults
    ])
    plt.figure(
        figsize=(8, 5)
    )
    plt.imshow(
        sensitivity,
        aspect="auto",
        origin="lower"
    )
    plt.colorbar(
        label="Sensitivity"
    )
    plt.xticks(
        np.arange(len(pressures)),
        [f"{p:g}" for p in pressures],
        rotation=45
    )
    plt.yticks(
        np.arange(len(simPore)),
        [f"{p:.2f}" for p in simPore]
    )
    plt.xlabel(
        "Pressure (kPa)"
    )
    plt.ylabel(
        "Pore Size (nm)"
    )
    plt.title(
        f"{perturbationResults['sample']} - {title}"
    )
    plt.tight_layout()
    plt.show(block=False)
def exportParameterFittingToExcel( analysis, filename, prefix=None):
    """
    Export parameter fitting analysis and all candidate model
    fitting results to Excel.
    Parameters
    ----------
    analysis : dict
        Output of fitValidationParameters().
    filename : str
        Output Excel filename.
    prefix : str or None
        Optional prefix for sheet names.
    Notes
    -----
    DataFrame results are exported directly.
    The nested "all_fits" structure is automatically converted
    into a DataFrame with "Parameter" and "Model" columns.
    Other non-DataFrame results are not exported unless they
    are explicitly handled here.
    """
    exported_sheets = []
    # ==========================================================
    # 1. Loop over validation types
    # ==========================================================
    for analysis_type, contents in analysis.items():
        if not isinstance(contents, dict):
            continue
        # ======================================================
        # 2. Export normal DataFrame results
        # ======================================================
        for result_type, result in contents.items():
            # --------------------------------------------------
            # Skip all_fits here
            # --------------------------------------------------
            if result_type == "all_fits":
                continue
            if not isinstance(result, pd.DataFrame):
                continue
            # --------------------------------------------------
            # Sheet name
            # --------------------------------------------------
            if prefix:
                sheet_name = ( f"{prefix}_{analysis_type}_{result_type}" )
            else:
                sheet_name = ( f"{analysis_type}_{result_type}" )
            sheet_name = sheet_name[:31]
            # --------------------------------------------------
            # Export
            # --------------------------------------------------
            fl.export_to_excel_auto( result, filename=filename, sheet_name=sheet_name )
            exported_sheets.append(sheet_name)
        # ======================================================
        # 3. Export all candidate fitting results
        # ======================================================
        all_fits = contents.get("all_fits")
        if not isinstance(all_fits, dict):
            continue
        records = []
        # ------------------------------------------------------
        # Each parameter: Slope / Intercept
        # ------------------------------------------------------
        for parameter_name, model_fits in all_fits.items():
            if not isinstance(model_fits, dict):
                continue
            # --------------------------------------------------
            # Each candidate model
            # --------------------------------------------------
            for model_name, fit_result in model_fits.items():
                if not isinstance(fit_result, dict):
                    continue
                record = {
                    "Parameter": parameter_name,
                    "Model": model_name
                }
                # --------------------------------------------------
                # Add fitting results
                # --------------------------------------------------
                for key, value in fit_result.items():
                    # Skip arrays / Series
                    if isinstance( value, (list, tuple, np.ndarray, pd.Series) ):
                        continue
                    # Skip nested dictionaries
                    if isinstance(value, dict):
                        continue
                    # Function object -> save function name
                    if callable(value):
                        record[key] = value.__name__
                        continue
                    record[key] = value
                records.append(record)
        # ======================================================
        # 4. Export all_fits
        # ======================================================
        if records:
            df_all_fits = pd.DataFrame(records)
            if prefix:
                sheet_name = ( f"{prefix}_{analysis_type}_all" )
            else:
                sheet_name = ( f"{analysis_type}_all" )
            sheet_name = sheet_name[:31]
            fl.export_to_excel_auto( df_all_fits, filename=filename, sheet_name=sheet_name )
            exported_sheets.append(sheet_name)
    # ==========================================================
    # 5. Print information
    # ==========================================================
    print( f"\nParameter fitting results exported to: {filename}" )
    if exported_sheets:
        print("Sheets:")
        for sheet in exported_sheets:
            print(f"  - {sheet}")
def exportAnalysisToExcel( analysis, filename, prefix=None):
    """
    Export analysis DataFrames to Excel.
    Parameters
    ----------
    analysis : dict
        Nested analysis dictionary.
    filename : str
        Output Excel filename.
    prefix : str or None
        Optional prefix for sheet names.
    """
    exported_sheets = []
    for analysis_type, contents in analysis.items():
        if not isinstance(contents, dict):
            continue
        for result_type, df in contents.items():
            if not isinstance(df, pd.DataFrame):
                continue
            # ==================================================
            # Sheet name
            # ==================================================
            if prefix:
                sheet_name = ( f"{prefix}_{analysis_type}_{result_type}" )
            else:
                sheet_name = ( f"{analysis_type}_{result_type}" )
            # Excel sheet-name limitation
            sheet_name = sheet_name[:31]
            # ==================================================
            # Export
            # ==================================================
            fl.export_to_excel_auto( df, filename=filename, sheet_name=sheet_name )
            exported_sheets.append(sheet_name)
    print( f"\nAnalysis results exported to: {filename}" )
    if exported_sheets:
        print("Sheets:")
        for sheet in exported_sheets:
            print(f"  - {sheet}")
def exportPredictionMetrics( sample, current_metrics, out_path, sheet_name="Metrics"):
    if current_metrics is None:
        return
    # ==========================================
    # 1. Convert input to DataFrame
    # ==========================================
    if isinstance(current_metrics, pd.DataFrame):
        df_metrics = current_metrics.copy()
        if "Sample" in df_metrics.columns:
            df_metrics = df_metrics[ df_metrics["Sample"] == sample ].copy()
    elif isinstance(current_metrics, dict):
        df_metrics = pd.DataFrame([current_metrics])
    else:
        raise TypeError( "current_metrics must be a pandas DataFrame or dict." )
    if df_metrics.empty:
        return
    # ==========================================
    # 2. Remove Sample column
    # ==========================================
    if "Sample" in df_metrics.columns:
        df_metrics = df_metrics.drop( columns=["Sample"] )
    df_metrics = df_metrics.reset_index(drop=True)
    # ==========================================
    # 3. Convert every structure to Series
    # ==========================================
    metric_data = {}
    for col in df_metrics.columns:
        value = df_metrics.loc[0, col]
        if isinstance(value, pd.Series):
            metric_data[col] = value.reset_index(drop=True)
        elif isinstance(value, (np.ndarray, list, tuple)):
            metric_data[col] = pd.Series(value)
        else:
            metric_data[col] = pd.Series([value])
    # ==========================================
    # 4. Build DataFrame automatically
    # ==========================================
    df_export = pd.DataFrame(metric_data)
    # ==========================================
    # 5. Export
    # ==========================================
    fl.export_to_excel_auto( df_export, filename=out_path, sheet_name=sheet_name, )
def exportPsdWeightSampleDetail( sample, simulationDetails, expData, fitMetrics, out_path):
    detail = simulationDetails[sample]
    pressure = detail["pressures"]
    simPore = detail["simPore"]
    psdPore = detail["psdPore"]
    psdDV = detail["psdDV"]
    ignoredFraction = detail["ignoredFraction"]
    weight = detail["weight"]
    volume = detail["volume"]
    uptake = detail["uptake"]
    uptake_density = detail["uptake_density"]
    metrics = detail["metrics"]
    # ======================================================
    # 1. PSD
    # ======================================================
    df_psd = pd.DataFrame({
        "MergedPore Size": pd.Series(psdPore),
        "MergedPSD dV/dlogD": pd.Series(psdDV),
        "simPore": pd.Series(simPore),
        "simVolume": pd.Series(volume),
    })
    fl.export_to_excel_auto( df_psd, filename=out_path, sheet_name="PSD", )
    # ======================================================
    # 2. Contribution
    # ======================================================
    contribution = uptake["contribution"]
    contributionPercent = uptake["percent"]
    cumulative = uptake["cumulative"]
    df_contribution = pd.DataFrame({
        "simPore": simPore,
        "lower boundary": detail["boundary"]["lower"],
        "upper boundary": detail["boundary"]["upper"],
        "Volume": pd.Series(volume),
        "Weight": weight,
        "Contribution": contribution[:, -1],
        "Contribution (%)": contributionPercent[:, -1] * 100,
        "Cumulative (%)": cumulative[:, -1] * 100,
    })
    fl.export_to_excel_auto( df_contribution, filename=out_path, sheet_name="Contribution", )
    # ======================================================
    # 3. Isotherm
    # ======================================================
    df_isotherm = pd.DataFrame({
        "Pressure": pd.Series(pressure),
        "PSD-weighted Uptake": pd.Series(uptake["total"]),
        "expPressure": pd.Series(expData[sample]["pressure"]),
        "Experimental Uptake": pd.Series(expData[sample]["expUptake"]),
    })
    fl.export_to_excel_auto( df_isotherm, filename=out_path, sheet_name="Isotherm", )
    # ======================================================
    # 4. Heatmap
    # ======================================================
    heatmap = pd.DataFrame( contribution, index=simPore, 
                           columns=detail["pressures"])
    heatmap.index.name = "Pore Size"
    fl.export_to_excel_auto( heatmap, filename=out_path, sheet_name="Heatmap", index=True, )
    # ======================================================
    # 5. Fit metrics
    # ======================================================
    if sample in fitMetrics:
        current_metrics = fitMetrics[sample]
        exportPredictionMetrics( sample, current_metrics, out_path, sheet_name="Metrics-simVSexp" )
def exportPredictionSampleDetail( predictionResults, predictionAnalysis, file_path):
    df_data_pred = predictionAnalysis["data_pred"]
    df_data_all = predictionAnalysis["data_all"]
    df_fit_details = predictionAnalysis["fit_details"]
    df_fit_metrics = predictionAnalysis["fit_metrics"]
    for sample in predictionResults:
        # ==================================================
        # 1. Output path
        # ==================================================
        out_path = fl.get_expanded_name( file_path, sample, expand="PSD_weighted", expandPos=True, type="xlsx" )
        result = predictionResults[sample]
        # ==================================================
        # 2. Prediction
        # ==================================================
        df_prediction = pd.DataFrame({
            "Pressure": pd.Series(result["pressure"]),
            "PSD-weighted Uptake": pd.Series(result["simUptake"]),
            "Slope": pd.Series(result["slope"]),
            "Intercept": pd.Series(result["intercept"]),
            "Predicted Uptake": pd.Series(result["predicted"]),
        })
        fl.export_to_excel_auto( df_prediction, filename=out_path, sheet_name="Prediction", )
        # ==================================================
        # 3. Analysis
        # ==================================================
        current_pred = df_data_pred[ df_data_pred["Sample"] == sample ].copy()
        current_all = df_data_all[ df_data_all["Sample"] == sample ].copy()
        current_fit = pd.DataFrame()
        if not df_fit_details.empty:
            current_fit = df_fit_details[ df_fit_details["Sample"] == sample ].copy()
        analysis_data = {}
        # data_pred
        for col in current_pred.columns:
            if col != "Sample":
                analysis_data[f"Pred {col}"] = pd.Series( current_pred[col].values )
        # data_all
        for col in current_all.columns:
            if col != "Sample":
                analysis_data[f"All {col}"] = pd.Series( current_all[col].values )
        # fit details
        for col in current_fit.columns:
            if col != "Sample":
                analysis_data[f"Fit {col}"] = pd.Series( current_fit[col].values )
        df_analysis = pd.DataFrame(analysis_data)
        fl.export_to_excel_auto( df_analysis, filename=out_path, sheet_name="Analysis", )
        # ==================================================
        # 4. Metrics
        # ==================================================
        exportPredictionMetrics( sample, df_fit_metrics, out_path )
def exportPsdWeightedSimulationResult( file_path, sample, simulationDetails, expData, fitMetrics):
    out_path = fl.get_expanded_name( file_path, sample, expand="PSD_weighted", expandPos=True, type="xlsx" )
    exportPsdWeightSampleDetail( sample, simulationDetails, expData, fitMetrics, out_path )
def debugPsdWeightedSimulation( sample, simPore, pressures, uptake, weight=None, 
                               expData=None, uptake_density=None, flagUsingDensity=False):
    """
    Debug and visualize PSD-weighted GCMC simulation results.
    This function can be used for both:
        1. Original PSD-weighted simulation
        2. Pore-volume perturbation simulation
    Parameters
    ----------
    sample : str
        Sample name or perturbation case name.
    simPore : array-like
        Fixed GCMC pore-size windows.
    pressures : array-like
        Simulation pressure points.
    uptake : dict
        Uptake result from calculateUptakeByWeight() or
        calculateUptakeByDensity().
    weight : array-like, optional
        PSD/GCMC weighting values.
    expData : dict, optional
        Experimental data dictionary.
    uptake_density : dict, optional
        Density-based uptake result.
    flagUsingDensity : bool, default=False
        Whether density-based weighting is used as the main method.
    """
    print("-" * 70)
    print(f"Debugging PSD-weighted simulation: {sample}")
    print("-" * 70)
    # ======================================================
    # 1. Check required uptake data
    # ======================================================
    required_keys = [
        "total",
        "percent",
        "cumulative",
        "contribution",
    ]
    missing_keys = [ key for key in required_keys if key not in uptake ]
    if missing_keys:
        raise KeyError( f"Missing uptake keys for sample '{sample}': " f"{missing_keys}" )
    # ======================================================
    # 2. Basic information
    # ======================================================
    print(f"Number of GCMC pore windows : {len(simPore)}")
    print(f"Number of pressure points   : {len(pressures)}")
    if weight is not None:
        weight = np.asarray(weight, dtype=float)
        print(f"Number of weight values     : {len(weight)}")
        print(f"Weight sum                  : {np.sum(weight):.6g}")
    uptake_total = np.asarray( uptake["total"], dtype=float )
    print(f"Uptake array shape          : {uptake_total.shape}")
    # ======================================================
    # 3. Find contribution threshold
    # ======================================================
    findThreshold( simPore, uptake["cumulative"], uptake["percent"] )
    # ======================================================
    # 4. Contribution distribution
    # ======================================================
    plotContribution( simPore, uptake["percent"] )
    # ======================================================
    # 5. Cumulative contribution
    # ======================================================
    plotCumulative( simPore, uptake["cumulative"] )
    # ======================================================
    # 6. PSD-weighted contribution
    # ======================================================
    if weight is not None:
        plotPSDContribution( simPore, weight, uptake["percent"] )
    # ======================================================
    # 7. Contribution heatmap
    # ======================================================
    plotContributionHeatmap( simPore, pressures, uptake["contribution"] )
    # ======================================================
    # 8. Reconstructed isotherm
    # ======================================================
    plotReconstructedIsotherm( pressures, uptake["total"] )
    # ======================================================
    # 9. Experimental vs simulation
    # ======================================================
    if expData is not None and sample in expData:
        plt.ion()
        # --------------------------------------------------
        # 9.1 Experimental vs PSD-weighted simulation
        # --------------------------------------------------
        myPlt.plotCurve(
            data={ sample: ( expData[sample]["pressure"], expData[sample]["expUptake"] ), },
            fit={ sample: ( pressures, uptake["total"] ), },
            fit_label="PSD-weighted",
            marker=True,
            line=True,
            legendPosition="upper left"
        )
        # --------------------------------------------------
        # 9.2 Experimental vs density-based simulation
        # --------------------------------------------------
        if ( not flagUsingDensity and uptake_density is not None ):
            myPlt.plotCurve(
                data={ sample: ( expData[sample]["pressure"], expData[sample]["expUptake"] ), },
                fit={ sample: ( pressures, uptake_density["total"] ), },
                fit_label="PSD-weighted_density",
                marker=True,
                line=True,
                legendPosition="upper left"
            )
            # ------------------------------------------------
            # 9.3 Uptake vs density-based uptake
            # ------------------------------------------------
            myPlt.plotCurve(
                data={ sample: ( pressures, uptake["total"] ), },
                fit={ sample: ( pressures, uptake_density["total"] ), },
                fit_label="PSD-weighted_density",
                marker=True,
                line=False,
                legendPosition="upper left"
            )
    elif expData is not None:
        print( f"[Debug] Experimental data not found for sample: " f"{sample}" )
    print("-" * 70)
    print(f"Debug completed: {sample}")
    print("-" * 70)
def calculateWeightedSimulationFromVolume( volume, simData, ModelVolumeAcc, 
                                          flagUsingDensity=False, simData_density=None):
    weight = None
    uptake_density = None
    if flagUsingDensity:
        uptake = calculateUptakeByDensity( volume, simData )
    else:
        weight = calculateWeightByVolume( volume, ModelVolumeAcc )
        uptake = calculateUptakeByWeight( weight, simData )
        if simData_density is not None:
            uptake_density = calculateUptakeByDensity( volume, simData_density )
    return {
        "volume": volume,
        "weight": weight,
        "uptake": uptake,
        "uptake_density": uptake_density,
    }
def calculatePsdWeightedSimulation( file_path, sampleAll, expData, simPore, 
                                   pressures, simData, ModelVolumeAcc, flagUsingDensity=False, 
                                   simData_density=None, debug=False, export=False):
    boundary = gcmcBoundary(simPore)
    sampleResults = {}
    simulationDetails = {}
    fitMetrics = {}
    debug_done = False
    for sample in sampleAll:
        # ======================================================
        # 1. Prepare original PSD
        # ======================================================
        psdPore, psdDV = mergePSD(sampleAll[sample])
        ignoredFraction = calculateIgnoredFraction( psdPore.copy(), simPore )
        psdPore = extendPSDToSimulationRange( psdPore, simPore )
        # ======================================================
        # 2. PSD → pore-window volume
        # ======================================================
        volume = calculatePSDVolume( psdPore, psdDV, boundary )
        # ======================================================
        # 3. Common simulation calculation
        # ======================================================
        uptakeData = calculateWeightedSimulationFromVolume( volume=volume, simData=simData, ModelVolumeAcc=ModelVolumeAcc, 
                                                           flagUsingDensity=flagUsingDensity, simData_density=simData_density )
        weight = uptakeData["weight"]
        uptake = uptakeData["uptake"]
        uptake_density = uptakeData["uptake_density"]
        # ======================================================
        # 4. Build original results
        # ======================================================
        sampleResults[sample] = { "pressure": pressures, "simUptake": uptake["total"], }
        # ======================================================
        # 5. Fit metrics
        # ======================================================
        metrics = None
        if sample in expData:
            metrics = mf.computeFitMetrics( fit_x=pressures, fit_y=uptake["total"], 
                                           exp_x=expData[sample]["pressure"], exp_y=expData[sample]["expUptake"] )
            fitMetrics[sample] = metrics
        # ======================================================
        # 6. Store simulation details
        # ======================================================
        simulationDetails[sample] = {
            "pressures": pressures,
            "psdPore": psdPore,
            "psdDV": psdDV,
            "simPore": simPore,
            "boundary": boundary,
            "ignoredFraction": ignoredFraction,
            "weight": weight,
            "volume": volume,
            "uptake": uptake,
            "uptake_density": uptake_density,
            "metrics": metrics,
        }
        # ======================================================
        # 7. Debug
        # ======================================================
        if debug and not debug_done:
            debugPsdWeightedSimulation( sample=sample, simPore=simPore, pressures=pressures, 
                                       uptake=uptake, weight=weight, expData=expData, 
                                       uptake_density=uptake_density, flagUsingDensity=flagUsingDensity )
            debug_done = True
        # ======================================================
        # 8. Export
        # ======================================================
        if export:
            exportPsdWeightedSimulationResult( file_path=file_path, sample=sample, 
                                              simulationDetails=simulationDetails, expData=expData, 
                                              fitMetrics=fitMetrics )
    return sampleResults, simulationDetails, fitMetrics
def predictExperimentalFromSimulation( sampleResults, modelfit, validation_type="absolute"):
    """
    Predict experimental uptake from simulation uptake using
    the fitted pressure-dependent Slope and Intercept models.
    Model hierarchy
    ---------------
    Slope(P)     = f_slope(P)
    Intercept(P) = f_intercept(P)
    q_exp_pred(P) =
        Slope(P) * q_sim(P) + Intercept(P)
    Parameters
    ----------
    sampleResults : dict
        PSD-weighted simulation results.
        {
            sample: {
                "pressure": pressures,
                "uptake": simulation_uptake
            }
        }
    modelfit : dict
        Output from fitValidationParameters().
    validation_type : str
        Validation type, e.g. "absolute".
    Returns
    -------
    predictionResults : dict
        Predicted experimental results for each sample.
    """
    predictionResults = {}
    # ==========================================================
    # 1. Get complete fitting results
    # ==========================================================
    fit_results = modelfit[ validation_type ]["fits"]
    slope_fit = fit_results["Slope"]
    intercept_fit = fit_results["Intercept"]
    # ==========================================================
    # 2. Get fitted models
    # ==========================================================
    slope_model = slope_fit["FitFunc"]
    intercept_model = intercept_fit["FitFunc"]
    # ==========================================================
    # 3. Get fitted parameters
    # ==========================================================
    slope_params = slope_fit["Parameters"]
    intercept_params = intercept_fit["Parameters"]
    # ==========================================================
    # 4. Apply model to each sample
    # ==========================================================
    for sample, result in sampleResults.items():
        pressure = np.asarray( result["pressure"], dtype=float )
        q_sim = np.asarray( result["simUptake"], dtype=float )
        # ------------------------------------------------------
        # 4.1 Predict pressure-dependent Slope
        # ------------------------------------------------------
        slope = slope_model( pressure, **slope_params )
        # ------------------------------------------------------
        # 4.2 Predict pressure-dependent Intercept
        # ------------------------------------------------------
        intercept = intercept_model( pressure, **intercept_params )
        # ------------------------------------------------------
        # 4.3 Apply exp = f(sim)
        #
        # q_exp = Slope(P) * q_sim + Intercept(P)
        # ------------------------------------------------------
        q_pred = mf.linearModel( q_sim, slope, intercept )
        # ------------------------------------------------------
        # 4.4 Store
        # ------------------------------------------------------
        predictionResults[sample] = {
            "pressure": pressure,
            "simUptake": q_sim,
            "slope": slope,
            "intercept": intercept,
            "predicted": q_pred,
        }
    return predictionResults
def calculatePerturbationByAddition(
        sample,
        simulationDetails,
        predictionResults,
        simData,
        ModelVolumeAcc,
        flagUsingDensity=False,
        simData_density=None,
        modelfit=None,
        validation_type="absolute",
        delta_volume=0.0,
        plotIndex=None,
        debug=False):
    """
    Perform pore-volume addition perturbation analysis for one sample.
    For each GCMC pore-size window, a fixed pore volume increment
    delta_volume is added to that pore-size window.
    The function calculates:
        - modified PSD-weighted simulation uptake
        - change in simulated uptake
        - simulated uptake sensitivity
        - predicted experimental uptake
        - change in predicted uptake
        - predicted uptake sensitivity
    Parameters
    ----------
    sample : str
        Sample name.
    simulationDetails : dict
        Output containing PSD/GCMC weighting information for samples.
    predictionResults : dict
        Baseline prediction results obtained from
        predictExperimentalFromSimulation().
    simData : dict
        GCMC simulation data.
    ModelVolumeAcc : array-like
        Accumulated model pore volume information.
    flagUsingDensity : bool, default=False
        Whether density-weighted simulation is used.
    simData_density : dict, optional
        Density-based simulation data.
    modelfit : dict, optional
        Fitted simulation-to-experiment calibration model.
    validation_type : str, default="absolute"
        Validation model type used for prediction.
    delta_volume : float, default=0.0
        Pore volume added to each individual GCMC pore-size window.
    plotIndex : None, int, list, or "all", default=None
        Controls perturbation plotting.
        None      : no plotting
        int       : plot one pore index
        list/array: plot selected pore indices
        "all"     : plot all pore indices
    debug : bool, default=False
        Print detailed calculation information.
    Returns
    -------
    perturbationResults : dict
        Perturbation results for the specified sample.
    """
    # ==========================================================
    # 1. Basic validation
    # ==========================================================
    if delta_volume == 0:
        raise ValueError(
            "delta_volume must be non-zero."
        )
    if sample not in simulationDetails:
        raise KeyError(
            f"Sample '{sample}' not found in simulationDetails."
        )
    if sample not in predictionResults:
        raise KeyError(
            f"Sample '{sample}' not found in predictionResults."
        )
    # ==========================================================
    # 2. Get original simulation information
    # ==========================================================
    details = simulationDetails[sample]
    pressures = np.asarray(
        details["pressures"],
        dtype=float
    )
    volume = np.asarray(
        details["volume"],
        dtype=float
    ).copy()
    simPore = np.asarray(
        details["simPore"],
        dtype=float
    )
    # ==========================================================
    # 3. Get original prediction
    # ==========================================================
    originalPrediction = predictionResults[sample]
    originalSim = np.asarray(
        originalPrediction["simUptake"],
        dtype=float
    )
    originalPred = np.asarray(
        originalPrediction["predicted"],
        dtype=float
    )
    # ==========================================================
    # 4. Validate dimensions
    # ==========================================================
    if len(volume) != len(simPore):
        raise ValueError(
            f"Length mismatch for sample '{sample}': "
            f"volume={len(volume)}, simPore={len(simPore)}."
        )
    if len(originalSim) != len(pressures):
        raise ValueError(
            f"Length mismatch for sample '{sample}': "
            f"simUptake={len(originalSim)}, "
            f"pressures={len(pressures)}."
        )
    if len(originalPred) != len(pressures):
        raise ValueError(
            f"Length mismatch for sample '{sample}': "
            f"predicted={len(originalPred)}, "
            f"pressures={len(pressures)}."
        )
    # ==========================================================
    # 5. Process plotIndex
    # ==========================================================
    plotIndices = None
    if plotIndex is not None and plotIndex != "all":
        plotIndices = np.atleast_1d(
            plotIndex
        ).astype(int)
        invalidIndex = plotIndices[
            (plotIndices < 0) |
            (plotIndices >= len(simPore))
        ]
        if len(invalidIndex) > 0:
            raise IndexError(
                f"Invalid plotIndex for sample '{sample}': "
                f"{invalidIndex.tolist()}. "
                f"Valid range: 0-{len(simPore) - 1}."
            )
    # ==========================================================
    # 6. Initialize result
    # ==========================================================
    perturbationResults = {
        "sample": sample,
        # --------------------------
        # Original / baseline
        # --------------------------
        "pressure": pressures.copy(),
        "simPore": simPore.copy(),
        "original_volume": volume.copy(),
        "original_simUptake": originalSim.copy(),
        "original_predicted": originalPred.copy(),
        # --------------------------
        # Perturbation results
        # --------------------------
        "pore_results": []
    }
    # ==========================================================
    # 7. Loop over GCMC pore-size windows
    # ==========================================================
    for target_index in range(len(simPore)):
        target_pore = simPore[target_index]
        # ------------------------------------------------------
        # 7.1 Add pore volume
        # ------------------------------------------------------
        result = calculateSimByAddition(
            volume=volume,
            target_index=target_index,
            delta_volume=delta_volume,
            simPore=simPore,
            simData=simData,
            ModelVolumeAcc=ModelVolumeAcc,
            flagUsingDensity=flagUsingDensity,
            simData_density=simData_density
        )
        # ------------------------------------------------------
        # 7.2 Extract modified simulation
        # ------------------------------------------------------
        modified_volume = np.asarray(
            result["modified_volume"],
            dtype=float
        )
        modifiedSim = np.asarray(
            result["uptake"]["total"],
            dtype=float
        )
        # ------------------------------------------------------
        # 7.3 Prepare one-sample prediction input
        # ------------------------------------------------------
        modifiedSampleResults = {
            sample: {
                "pressure": pressures.copy(),
                "simUptake": modifiedSim.copy()
            }
        }
        # ------------------------------------------------------
        # 7.4 Predict experimental uptake
        # ------------------------------------------------------
        modifiedPredictionResults = predictExperimentalFromSimulation(
            sampleResults=modifiedSampleResults,
            modelfit=modelfit,
            validation_type=validation_type
        )
        modifiedPrediction = modifiedPredictionResults[sample]
        modifiedPred = np.asarray(
            modifiedPrediction["predicted"],
            dtype=float
        )
        # ------------------------------------------------------
        # 7.5 Calculate changes
        # ------------------------------------------------------
        deltaSim = (
            modifiedSim -
            originalSim
        )
        deltaPred = (
            modifiedPred -
            originalPred
        )
        # ------------------------------------------------------
        # 7.6 Calculate sensitivity
        # ------------------------------------------------------
        sensitivitySim = (
            deltaSim /
            delta_volume
        )
        sensitivityPred = (
            deltaPred /
            delta_volume
        )
        # ------------------------------------------------------
        # 7.7 Store pore-specific result
        # ------------------------------------------------------
        poreResult = {
            "target_index": target_index,
            "target_pore": target_pore,
            "delta_volume": delta_volume,
            "modified_volume": modified_volume.copy(),
            # Simulation
            "modified_simUptake": modifiedSim.copy(),
            "delta_simUptake": deltaSim.copy(),
            "sensitivity_sim": sensitivitySim.copy(),
            # Predicted experiment
            "modified_predicted": modifiedPred.copy(),
            "delta_predicted": deltaPred.copy(),
            "sensitivity_pred": sensitivityPred.copy(),
            # Calibration parameters
            "slope": modifiedPrediction["slope"],
            "intercept": modifiedPrediction["intercept"],
            "operation": "addition"
        }
        perturbationResults["pore_results"].append(
            poreResult
        )
        # ======================================================
        # 7.8 Optional plotting
        # ======================================================
        flagPlot = False
        if plotIndex == "all":
            flagPlot = True
        elif plotIndices is not None:
            flagPlot = (
                target_index in plotIndices
            )
        if flagPlot:
            plotPerturbationAddition( sample=sample, perturbationResults=perturbationResults, poreResult=poreResult )
        # ======================================================
        # 7.9 Debug output
        # ======================================================
        if debug:
            print(
                f"[Perturbation] "
                f"{sample} | "
                f"index={target_index} | "
                f"pore={target_pore:.4f} nm | "
                f"deltaV={delta_volume}"
            )
    plotSensitivityHeatmap( perturbationResults, sensitivity_type="pred" )
    return perturbationResults
def main(file_path=None):
    # ==== 输入参数 ====
    flagUsingDensity = False
    pressure_min = 10 #kPa
    ########### 读取数据 ###########
    file_path = select_folder()   # 改成你的文件路径
    HeliumFraction, ModelVolumeAcc = readVolumeAndHeliumVoidFraction(file_path,sheet_name="volume")
    if flagUsingDensity:
        simPore,pressures,density_simu_acc = readSimulationDensity(file_path,HeliumFraction,pressure_min=pressure_min)
        simData = density_simu_acc
    else:
        simPore,pressures,uptake_simu_acc = readSimulationUptake(file_path,pressure_min=pressure_min)
        simData = uptake_simu_acc
        _,_,density_simu_acc = readSimulationDensity(file_path,HeliumFraction,pressure_min=pressure_min)
        simData_density = density_simu_acc
    sampleAll = readPSD(file_path)
    expData =readExpUptake(file_path)
    ########### 读取数据 ###########
    sampleResults, simulationDetails, fitMetrics = \
        calculatePsdWeightedSimulation( file_path=file_path, sampleAll=sampleAll, 
                                       expData=expData, simPore=simPore, pressures=pressures, 
                                       simData=simData, ModelVolumeAcc=ModelVolumeAcc, flagUsingDensity=flagUsingDensity,
                                         simData_density=simData_density, debug=False, export=False )
    validation = validationSimAndExp( sampleResults=sampleResults, expData=expData, pressures=pressures, )
    out_path2 = fl.get_expanded_name(file_path, fileName="sim2exp", expandPos=True, type="xlsx")
    # exportAnalysisToExcel(validation,out_path2)
    modelfixed = "Exponential Saturation"
    modelfixed = None
    modelfit = fitValidationParameters(validation,plotflag=False,slope_models=modelfixed, pressure_min=pressure_min,
                                       intercept_models=modelfixed)
    # exportParameterFittingToExcel(modelfit,out_path2,prefix="fit")
    #############使用新的样品验证模型拟合是否合适###################
    sampleCheck = readPSD(file_path, sheet_name="checkPSD")
    expDataCheck = readExpUptake(file_path, sheet_name="checkExp")
    sampleResultsCheck, simulationDetailsCheck, _ = \
        calculatePsdWeightedSimulation( file_path=file_path, sampleAll=sampleCheck, 
                                       expData=expDataCheck, simPore=simPore, pressures=pressures, 
                                       simData=simData, ModelVolumeAcc=ModelVolumeAcc, 
                                       flagUsingDensity=flagUsingDensity, simData_density=simData_density, 
                                       debug=False, export=False )
    predictionResults = predictExperimentalFromSimulation( sampleResults=sampleResultsCheck, 
                                                          modelfit=modelfit, validation_type="absolute" )
    predictionAnalysis = preparePredictionAnalysisData( predictionResults=predictionResults, 
                                                       expData=expDataCheck, modelfit=modelfit, validation_type="absolute" )
    # plotPredictionAnalysisByPressure( predictionAnalysis, sample=None, pressure_type="merge",)
    plotPredictionAnalysisByPressure( predictionAnalysis, sample=None, pressure_type="sim",)
    # exportPredictionSampleDetail( predictionResults, predictionAnalysis, file_path)
    df_debug = debugPredictionResults( predictionResults, expDataCheck )
    #############使用新的样品验证模型拟合是否合适###################
    for sample in sampleCheck:
        perturbationResult = calculatePerturbationByAddition(
            sample=sample,
            simulationDetails=simulationDetailsCheck,
            predictionResults=predictionResults,
            simData=simData,
            ModelVolumeAcc=ModelVolumeAcc,
            flagUsingDensity=flagUsingDensity,
            simData_density=simData_density,
            modelfit=modelfit,
            validation_type="absolute",
            delta_volume=0.01,
            plotIndex=[5], #plotIndex="all", plotIndex=[0, 3, 7, 12]
            debug=True
            )
        # break
    plt.show(block=True)
if __name__ == "__main__": 
    f = sys.argv[1] if len(sys.argv) > 1 else None
    main(f)