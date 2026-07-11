##每一个子平台固定开头，用于找到依赖
import os
import sys
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from common import *

import constantsAndName as const

import fileSource.fileOperation as fl
import dataProcessSource.dataOperation as dop

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
def plotCompareIsotherm(sample, pressures, simData, expData):
    exp = expData[sample]
    simPressure = pressures / 1000
    simUptake = simData["total"]
    expPressure = exp["pressure"]
    expUptake = exp["expUptake"]
    plt.figure(figsize=(6,4))
    plt.plot( simPressure, simUptake, "-o", linewidth=2, markersize=5, label="PSD-weighted" )
    plt.scatter( expPressure, expUptake, marker="s", s=50, label="Experiment" )
    plt.xlabel("Pressure (kPa)")
    plt.ylabel("Uptake (mol/kg)")
    plt.title(sample)
    plt.legend()
    plt.tight_layout()
def plotCompareIsothermTwoMethod(sample, pressures, simData, simData_density):
    exp = simData_density
    simPressure = pressures / 1000
    simUptake = simData["total"]
    expPressure = simPressure
    expUptake = exp["total"]
    plt.figure(figsize=(6,4))
    plt.plot( simPressure, simUptake, "-o", linewidth=2, markersize=5, label="PSD-weighted" )
    plt.scatter( expPressure, expUptake, marker="s", s=50, label="PSD-weighted_density" )
    plt.xlabel("Pressure (kPa)")
    plt.ylabel("Uptake (mol/kg)")
    plt.title(sample)
    plt.legend()
    plt.tight_layout()

def exportResult(simPore,boundary,pressures,weight,uptake,out_path):
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
    isotherm = pd.DataFrame()
    isotherm["Pressure (kPa)"] = pressures / 1000
    isotherm["PSD-weighted Uptake (mol/kg)"] = totalUptake
    heatmap = pd.DataFrame(contribution)
    heatmap.index = simPore
    heatmap.columns = pressures
    with pd.ExcelWriter(out_path) as writer:
        result.to_excel(writer,sheet_name="Contribution",index=False)
        isotherm.to_excel(writer,sheet_name="Isotherm",index=False)
        heatmap.to_excel(writer,sheet_name="Heatmap")
def calculateUptakeByAdsorptionDensity():
    y = 0


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
    for sample in sampleAll:
        sample = "CC-Ca-800-2-1"
        print(sample)
        psdPore,psdDV = mergePSD(sampleAll[sample])
        ignoredFraction = calculateIgnoredFraction(psdPore.copy(),simPore)
        psdPore = extendPSDToSimulationRange(psdPore,simPore)
        
        weight, volume = calculatePSDVolumeAndWeight(psdPore,psdDV,boundary)
        if flagUsingDensity:
            uptake = calculateUptakeByDensity(volume, simData)
        else:
            # volumeWweight = volume/ModelVolumeAcc
            # weight = volumeWweight
            uptake = calculateUptakeByWeight(weight, simData)
            uptake_density = calculateUptakeByDensity(volume, simData_density)
            
        # print("ModelVolumeAcc= ", ModelVolumeAcc)
        print("sum(volume)= ", np.sum(volume))
        print("volume= ", volume)
        # print("volumeWweight= ", volume/ModelVolumeAcc)
        print("weight= ", weight)
        print("simPore= ",simPore)
        print("HeliumFraction= ",HeliumFraction)
        print("simData= ",np.max(simData,axis=1))
        findThreshold(simPore,uptake["cumulative"],uptake["percent"])
        # plotContribution(simPore,uptake["percent"])
        # plotCumulative(simPore,uptake["cumulative"])
        # plotPSDContribution(simPore,weight,uptake["percent"])
        # plotContributionHeatmap(simPore,pressures,uptake["contribution"])
        # plotReconstructedIsotherm(pressures,uptake["total"])
        plt.ion()
        plotCompareIsotherm(sample, pressures, uptake, expData)
        if not flagUsingDensity:
            plotCompareIsothermTwoMethod(sample, pressures, uptake, uptake_density)
        out_path = fl.get_expanded_name(file_path, sample, expand="PSD_weighted", expandPos=True, type="xlsx")
        # exportResult(simPore,boundary,pressures,weight,uptake,out_path)
        break
    plt.show(block=True)

if __name__ == "__main__":
    f = sys.argv[1] if len(sys.argv) > 1 else None
    main(f)