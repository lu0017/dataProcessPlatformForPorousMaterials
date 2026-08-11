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
            "data": DataFrame,
            "plot_data": DataFrame,
            "fit_details": DataFrame
        }
    """
    data_all = []
    plot_data_all = []
    fit_details_all = []
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
        uptake_pred = np.asarray( predictionResults[sample]["predicted"], dtype=float )
        # Prediction interpolation
        _, _, pred_interp = mf.interpolateData( pressure_pred, uptake_pred, extrapolate=True )
        # ======================================================
        # 3. Common pressure grid
        # ======================================================
        pressure_plot = np.sort( np.unique( np.concatenate([ pressure_exp, pressure_pred ]) ) )
        # ======================================================
        # 4. Interpolated data
        # ======================================================
        uptake_exp_plot = exp_interp( pressure_plot )
        uptake_pred_plot = pred_interp( pressure_plot )
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
        for p, exp_q, pred_q, diff, rat in zip( pressure_pred, exp_at_pred, uptake_pred, difference, ratio):
            data_all.append({
                "Sample": sample,
                "Validation type": validation_type,
                "Pressure (kPa)": p,
                "Experimental": exp_q,
                "Predicted": pred_q,
                "Difference": diff,
                "Ratio": rat,
            })
        # ======================================================
        # 6. Plot data
        #
        #    Dense common grid
        # ======================================================
        for p, exp_q, pred_q in zip(
                pressure_plot,
                uptake_exp_plot,
                uptake_pred_plot):
            plot_data_all.append({
                "Sample": sample,
                "Pressure (kPa)": p,
                "Experimental": exp_q,
                "Predicted": pred_q,
            })
        # ======================================================
        # 7. Original prediction points
        #
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
    df_data = pd.DataFrame(data_all)
    df_plot_data = pd.DataFrame(plot_data_all)
    df_fit_details = pd.DataFrame(fit_details_all)
    return {
        "data": df_data,
        "plot_data": df_plot_data,
        "fit_details": df_fit_details,
    }
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
    required_columns = [
        pressure_col,
        parameter_col
    ]
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
    df_fit = df_parameter[
        required_columns
    ].copy()
    df_fit = df_fit.replace(
        [np.inf, -np.inf],
        np.nan
    )
    df_fit = df_fit.dropna(
        subset=required_columns
    )
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
        myPlt.plotSingleCorrelation(
            fit_result,
            xlabel=pressure_col,
            ylabel=parameter_name
        )
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
    required_columns = [
        pressure_col,
        parameter_col
    ]
    missing_columns = [
        col
        for col in required_columns
        if col not in df_parameter.columns
    ]
    if missing_columns:
        raise KeyError(
            f"Missing columns for {parameter_name}: "
            f"{missing_columns}"
        )
    df_fit = df_parameter[
        required_columns
    ].copy()
    df_fit = df_fit.replace(
        [np.inf, -np.inf],
        np.nan
    )
    df_fit = df_fit.dropna(
        subset=required_columns
    )
    # Sort by pressure
    df_fit = df_fit.sort_values(
        by=pressure_col
    ).reset_index(drop=True)
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
        selected_models = {
            models: mf.model_library[models]
        }
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
            selected_models[model_name] = (
                mf.model_library[model_name]
            )
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
            missing_config = [
                key
                for key in required_config
                if key not in config
            ]
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
                bounds=config.get(
                    "bounds",
                    (-np.inf, np.inf)
                ),
                param_names=config.get(
                    "param_names"
                ),
                plotflag=False
            )
            if fit_result is None:
                continue
            # --------------------------------------------------
            # Store complete fitting result
            # --------------------------------------------------
            fits[model_name] = fit_result
            # --------------------------------------------------
            # Convert fitting result to record
            # --------------------------------------------------
            record = mf.fittingResultToRecord(
                fit_result,
                parameter_name=parameter_name
            )
            # --------------------------------------------------
            # Number of fitted parameters
            # --------------------------------------------------
            if config["fit_method"] == "linear":
                k = 2
            else:
                param_names = config.get(
                    "param_names"
                )
                if param_names is None:
                    raise ValueError(
                        f"`param_names` is required for "
                        f"nonlinear model: {model_name}"
                    )
                k = len(param_names)
            # --------------------------------------------------
            # Calculate SSE
            # --------------------------------------------------
            rmse = record.get(
                "RMSE",
                np.nan
            )
            if not np.isfinite(rmse):
                sse = np.nan
            else:
                sse = n * rmse ** 2
            # --------------------------------------------------
            # Calculate AIC
            # --------------------------------------------------
            if (
                np.isfinite(sse)
                and sse > 0
            ):
                aic = (
                    n * np.log(sse / n)
                    + 2 * k
                )
            elif (
                np.isfinite(sse)
                and sse == 0
            ):
                aic = -np.inf
            else:
                aic = np.nan
            # --------------------------------------------------
            # Calculate AICc
            # --------------------------------------------------
            if (
                np.isfinite(aic)
                and n - k - 1 > 0
            ):
                aicc = (
                    aic
                    + 2 * k * (k + 1)
                    / (n - k - 1)
                )
            elif aic == -np.inf:
                aicc = -np.inf
            else:
                aicc = np.inf
            # --------------------------------------------------
            # Store model comparison record
            # --------------------------------------------------
            comparison_records.append({
                "Parameter": parameter_name,
                "Model": model_name,
                "N": n,
                "N_parameters": k,
                "R2": record.get(
                    "R2",
                    np.nan
                ),
                "RMSE": record.get(
                    "RMSE",
                    np.nan
                ),
                "AIC": aic,
                "AICc": aicc
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
    df_comparison = pd.DataFrame(
        comparison_records
    )
    if df_comparison.empty:
        print(
            f"No valid model was fitted for "
            f"{parameter_name}."
        )
        return None
    # ==========================================================
    # 6. Print model comparison
    # ==========================================================
    print(
        "\n================ MODEL COMPARISON ================"
    )
    print(
        df_comparison.to_string(
            index=False
        )
    )
    print(
        "==================================================="
    )
    # ==========================================================
    # 7. Find valid models
    # ==========================================================
    valid = df_comparison[
        df_comparison["AICc"].notna()
    ].copy()
    # Remove models whose fitting failed
    valid = valid[
        valid["Model"].isin(
            fits.keys()
        )
    ]
    if valid.empty:
        print(
            f"No valid fitted model for "
            f"{parameter_name}."
        )
        return None
    # ==========================================================
    # 8. Calculate Delta AICc
    # ==========================================================
    min_aicc = valid["AICc"].min()
    df_comparison["Delta_AICc"] = np.nan
    valid_indices = valid.index
    df_comparison.loc[
        valid_indices,
        "Delta_AICc"
    ] = (
        df_comparison.loc[
            valid_indices,
            "AICc"
        ]
        - min_aicc
    )
    # ==========================================================
    # 9. Select best model
    # ==========================================================
    best_idx = valid["AICc"].idxmin()
    best_model = df_comparison.loc[
        best_idx,
        "Model"
    ]
    best_fit = fits[
        best_model
    ]
    # ==========================================================
    # 10. Determine whether model was fixed
    # ==========================================================
    model_fixed = (
        isinstance(models, str)
        or (
            isinstance(models, (list, tuple))
            and len(models) == 1
        )
    )
    # Add selection information
    df_comparison["Selected"] = False
    df_comparison.loc[
        best_idx,
        "Selected"
    ] = True
    # ==========================================================
    # 11. Plot best model
    # ==========================================================
    if plotflag:
        myPlt.plotSingleCorrelation(
            best_fit,
            xlabel=pressure_col,
            ylabel=parameter_name
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
        df_metrics = validation[
            validation_type
        ]["metrics"]
        if df_metrics is None or df_metrics.empty:
            continue
        # ======================================================
        # 3. Required columns
        # ======================================================
        required_columns = [
            pressure_col,
            "Fit_slope",
            "Fit_intercept"
        ]
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
        df_parameter = df_metrics[
            required_columns
        ].copy()
        df_parameter = df_parameter.replace(
            [np.inf, -np.inf],
            np.nan
        )
        df_parameter = df_parameter.dropna(
            subset=required_columns
        )
        # ======================================================
        # 5. Pressure filter
        # ======================================================
        if pressure_min is not None:
            df_parameter = df_parameter[
                df_parameter[pressure_col]
                >= pressure_min
            ].copy()
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
        # ======================================================
        # 8. Build best-fit metrics
        # ======================================================
        slope_record = mf.fittingResultToRecord(
            slope_result["best_fit"],
            parameter_name="Slope"
        )
        intercept_record = mf.fittingResultToRecord(
            intercept_result["best_fit"],
            parameter_name="Intercept"
        )
        slope_record["Model"] = (
            slope_result["best_model"]
        )
        intercept_record["Model"] = (
            intercept_result["best_model"]
        )
        df_fit_metrics = pd.DataFrame([
            slope_record,
            intercept_record
        ])
        # ======================================================
        # 9. Store results
        # ======================================================
        parameter_fits[validation_type] = {
            "data": df_parameter,
            "metrics": df_fit_metrics,
            "fits": {
                "Slope":
                    slope_result["best_fit"],
                "Intercept":
                    intercept_result["best_fit"]
            },
            "model_comparison": {
                "Slope":
                    slope_result["model_comparison"],
                "Intercept":
                    intercept_result["model_comparison"]
            },
            "all_fits": {
                "Slope":
                    slope_result["fits"],
                "Intercept":
                    intercept_result["fits"]
            }
        }
    return parameter_fits
def analyzeNormalizedIsotherm(
        df,
        df_uptake_all,
        plotflag=False):
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
    #
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
        {
            "absolute": {
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
    plt.xlabel("Pressure (kPa)")
    plt.ylabel("Pore Size (nm)")
    plt.title("PSD-weighted Contribution")
    plt.tight_layout()
    plt.show()
def plotReconstructedIsotherm(pressures,totalUptake):
    plt.figure(figsize=(6,4))
    plt.plot(pressures,totalUptake,marker="o")
    plt.xscale("log")
    plt.xlabel("Pressure (kPa)")
    plt.ylabel("PSD-weighted Uptake (mol/kg)")
    plt.title("Reconstructed Isotherm")
    plt.tight_layout()
    plt.show()
def plotPredictionAnalysis( predictionAnalysis, sample=None, 
                           plot_original_prediction=False, fit_label="exp=f(sim)", legendPosition="upper left"):
    """
    Plot prediction analysis results.
    Parameters
    ----------
    predictionAnalysis : dict
        Output from preparePredictionAnalysisData().
    sample : str or None
        Sample to plot. If None, plot all samples.
    plot_original_prediction : bool
        Whether to overlay original prediction points.
    fit_label : str
        Label for prediction curve.
    """
    df_plot = predictionAnalysis["plot_data"]
    df_data = predictionAnalysis["data"]
    if sample is None:
        samples = df_plot["Sample"].unique()
    else:
        samples = [sample]
    for sample_name in samples:
        # ======================================================
        # Dense interpolated curves
        # ======================================================
        df_sample = df_plot[
            df_plot["Sample"] == sample_name
        ]
        myPlt.plotCurve(
            data={ sample_name: ( df_sample["Pressure (kPa)"].to_numpy(), df_sample["Experimental"].to_numpy() ) },
            fit={ sample_name: ( df_sample["Pressure (kPa)"].to_numpy(), df_sample["Predicted"].to_numpy() ) },
            fit_label=fit_label,
            marker=True,
            line=True,
            legendPosition=legendPosition
        )
        # ======================================================
        # Original prediction points
        # ======================================================
        if plot_original_prediction:
            df_original = df_data[ df_data["Sample"] == sample_name ]
            myPlt.plotCurve(
                data={
                    sample_name: (
                        df_original["Pressure (kPa)"].to_numpy(),
                        df_original["Predicted"].to_numpy()
                    )
                },
                marker=True,
                line=True,
                legendPosition=legendPosition
            )
    plt.show(block=False)
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
def exportPsdWeightDetail(simPore,boundary,pressures,weight,uptake,out_path):
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
    "Pressure (kPa)": pd.Series(pressures),
    "PSD-weighted Uptake (mol/kg)": pd.Series(totalUptake),
    })
    heatmap = pd.DataFrame(contribution)
    heatmap.index = simPore
    heatmap.columns = pressures
    with pd.ExcelWriter(out_path) as writer:
        result.to_excel(writer,sheet_name="Contribution",index=False)
        isotherm.to_excel(writer,sheet_name="Isotherm",index=False)
        heatmap.to_excel(writer,sheet_name="Heatmap")
def calculatePsdWeightedSimulation( file_path, sampleAll, expData, simPore, pressures, 
                                   simData, ModelVolumeAcc, flagUsingDensity=False, simData_density=None, 
                                   debug=False, export=False):
    """
    Calculate PSD-weighted simulation uptake for all samples,
    including diagnostic analysis, plotting, metrics calculation,
    and optional export of PSD-weight details.
    Parameters
    ----------
    file_path : str
        Input Excel file path.
    sampleAll : dict
        PSD information for all samples.
    expData : dict
        Experimental uptake data for all samples.
        Used for metrics calculation and comparison plots.
    simPore : array-like
        Simulation pore-size grid.
    pressures : array-like
        Simulation pressure grid.
        Unit: kPa.
    simData : array-like
        Simulation uptake/density data used for PSD weighting.
    ModelVolumeAcc : float
        Accessible simulation model volume.
    flagUsingDensity : bool
        If True:
            calculate uptake directly using pore volume and
            simulation density.
        If False:
            calculate PSD-volume weight using ModelVolumeAcc
            and calculate uptake using simulation uptake.
    simData_density : array-like or None
        Simulation density data.
        Used only when flagUsingDensity=False for comparison.
    debug : bool
        If True, perform threshold analysis and generate
        diagnostic plots.
    export : bool
        If True, export PSD-weighting details for each sample.
    Returns
    -------
    sampleResults : dict
        PSD-weighted simulation isotherms.
        {
            sample: {
                "pressure": pressures,
                "uptake": uptake["total"]
            }
        }
    simulationDetails : dict
        Detailed calculation results for each sample.
    fitMetrics : dict
        Fit metrics between experimental and PSD-weighted
        simulation uptake.
    """
    # ==========================================================
    # 1. Initialize
    # ==========================================================
    boundary = gcmcBoundary(simPore)
    sampleResults = {}
    simulationDetails = {}
    fitMetrics = {}
    # ==========================================================
    # 2. Calculate PSD-weighted simulation for each sample
    # ==========================================================
    for sample in sampleAll:
        print("=" * 70)
        print(f"Processing sample: {sample}")
        print("=" * 70)
        # ------------------------------------------------------
        # 2.1 Merge PSD
        # ------------------------------------------------------
        psdPore, psdDV = mergePSD( sampleAll[sample] )
        # ------------------------------------------------------
        # 2.2 Calculate ignored PSD fraction
        # ------------------------------------------------------
        ignoredFraction = calculateIgnoredFraction( psdPore.copy(), simPore )
        # ------------------------------------------------------
        # 2.3 Extend PSD to simulation range
        # ------------------------------------------------------
        psdPore = extendPSDToSimulationRange( psdPore, simPore )
        # ------------------------------------------------------
        # 2.4 Calculate PSD volume and weight
        # ------------------------------------------------------
        weight, volume = calculatePSDVolumeAndWeight( psdPore, psdDV, boundary )
        # ------------------------------------------------------
        # 2.5 Calculate PSD-weighted uptake
        # ------------------------------------------------------
        uptake_density = None
        if flagUsingDensity:
            uptake = calculateUptakeByDensity( volume, simData )
        else:
            # PSD volume → weighting factor
            volumeWweight = volume / ModelVolumeAcc
            weight = volumeWweight
            # PSD-weighted uptake
            uptake = calculateUptakeByWeight( weight, simData )
            # --------------------------------------------------
            # Optional density-based calculation
            # --------------------------------------------------
            if simData_density is not None:
                uptake_density = calculateUptakeByDensity( volume, simData_density )
        # ------------------------------------------------------
        # 2.6 Store final simulation result
        # ------------------------------------------------------
        sampleResults[sample] = { "pressure": pressures, "simUptake": uptake["total"], }
        # ------------------------------------------------------
        # 2.7 Calculate fit metrics against experiment
        # ------------------------------------------------------
        metrics = None
        if sample in expData:
            metrics = mf.computeFitMetrics( fit_x=pressures, fit_y=uptake["total"], 
                                        exp_x=expData[sample]["pressure"], exp_y=expData[sample]["expUptake"] )
            fitMetrics[sample] = metrics
        # ------------------------------------------------------
        # 2.8 Store detailed calculation results
        # ------------------------------------------------------
        simulationDetails[sample] = {
            "psdPore": psdPore,
            "psdDV": psdDV,
            "ignoredFraction": ignoredFraction,
            "weight": weight,
            "volume": volume,
            "uptake": uptake,
            "uptake_density": uptake_density,
            "metrics": metrics,
        }
        # ======================================================
        # 3. Debug / diagnostic analysis
        # ======================================================
        if debug:
            # --------------------------------------------------
            # 3.1 Threshold analysis
            # --------------------------------------------------
            findThreshold( simPore, uptake["cumulative"], uptake["percent"] )
            # --------------------------------------------------
            # 3.2 PSD contribution
            # --------------------------------------------------
            plotContribution( simPore, uptake["percent"] )
            # --------------------------------------------------
            # 3.3 Cumulative contribution
            # --------------------------------------------------
            plotCumulative( simPore, uptake["cumulative"] )
            # --------------------------------------------------
            # 3.4 PSD contribution weighted by pore volume
            # --------------------------------------------------
            plotPSDContribution( simPore, weight, uptake["percent"] )
            # --------------------------------------------------
            # 3.5 Contribution heatmap
            # --------------------------------------------------
            plotContributionHeatmap( simPore, pressures, uptake["contribution"] )
            # --------------------------------------------------
            # 3.6 Reconstructed isotherm
            # --------------------------------------------------
            plotReconstructedIsotherm( pressures, uptake["total"] )
            # --------------------------------------------------
            # 3.7 Experimental vs PSD-weighted simulation
            # --------------------------------------------------
            if sample in expData:
                plt.ion()
                myPlt.plotCurve(
                    data={ sample: ( expData[sample]["pressure"], expData[sample]["expUptake"] ), },
                    fit={ sample: ( pressures, uptake["total"] ), },
                    fit_label="PSD-weighted",
                    marker=True,
                    line=True,
                    legendPosition="upper left"
                )
                # ------------------------------------------------
                # 3.8 Experimental vs density-based simulation
                # ------------------------------------------------
                if ( not flagUsingDensity and uptake_density is not None ):
                    myPlt.plotCurve(
                        data={ sample: ( expData[sample]["pressure"], expData[sample]["expUptake"] ), },
                        fit={ sample: ( pressures, uptake_density["total"] ), },
                        fit_label="PSD-weighted_density",
                        marker=True,
                        line=True,
                        legendPosition="upper left"
                    )
                    # --------------------------------------------
                    # 3.9 Uptake vs density-based uptake
                    # --------------------------------------------
                    myPlt.plotCurve(
                        data={ sample: ( pressures, uptake["total"] ), },
                        fit={ sample: ( pressures, uptake_density["total"] ), },
                        fit_label="PSD-weighted_density",
                        marker=True,
                        line=False,
                        legendPosition="upper left"
                    )
        # ======================================================
        # 4. Export PSD weighting details
        # ======================================================
        if export:
            out_path = fl.get_expanded_name( file_path, sample, expand="PSD_weighted", expandPos=True, type="xlsx" )
            exportPsdWeightDetail( simPore, boundary, pressures, weight, uptake, out_path )
    return sampleResults, simulationDetails, fitMetrics
def predictExperimentalFromSimulation(
        sampleResults,
        modelfit,
        validation_type="absolute"):
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
def main(file_path=None):
    # ==== 输入参数 ====
    flagUsingDensity = False
    pressure_min = 0.0 #kPa
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
    out_path2 = fl.get_expanded_name(file_path, fileName="sim2expP1", expandPos=True, type="xlsx")
    exportAnalysisToExcel(validation,out_path2)
    modelfixed = "Exponential Saturation"
    modelfixed = None
    modelfit = fitValidationParameters(validation,plotflag=False,slope_models=modelfixed,intercept_models=modelfixed)
    exportAnalysisToExcel(modelfit,out_path2,prefix="fit")
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
    plotPredictionAnalysis( predictionAnalysis, sample=None )
    for sample in predictionResults:
        pressure_exp = expDataCheck[sample]["pressure"]
        _, _, exp_interp = mf.interpolateData( expDataCheck[sample]["pressure"],
                                                          expDataCheck[sample]["expUptake"], extrapolate=True)
        exp_uptake = exp_interp(pressures)
        uptake_exp = expDataCheck[sample]["expUptake"]
        pressure_pred = predictionResults[sample]["pressure"]
        uptake_pred = predictionResults[sample]["predicted"]
        myPlt.plotCurve(
            data={ sample: ( pressures, exp_uptake ), },
            fit={ sample: ( pressure_pred, uptake_pred ), },
            fit_label="exp=f(sim)",
            marker=True,
            line=False,
            legendPosition="upper left"
        )
        plt.show(block=False)
    df_debug = debugPredictionResults( predictionResults, expDataCheck )
    #############使用新的样品验证模型拟合是否合适###################
    plt.show(block=True)
if __name__ == "__main__": 
    f = sys.argv[1] if len(sys.argv) > 1 else None
    main(f)