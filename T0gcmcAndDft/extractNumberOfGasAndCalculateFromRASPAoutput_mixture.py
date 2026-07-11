from collections import defaultdict
import re
import sys
import csv
import tkinter as tk
from tkinter import filedialog
import re
import numpy as np
import matplotlib.pyplot as plt

def select_folder():
    root = tk.Tk()
    root.withdraw()

    file_path = filedialog.askopenfilename(
        title="Select output file",
        filetypes=[("output files", "*.data"), ("All files", "*.*")]
    )

    if not file_path:
        raise ValueError("No file selected!")

    print("Selected file:", file_path)
    return file_path

# =========================
# 1. 读取 RASPA 输出
# =========================
def extract_temperature(file_path):
    with open(file_path, 'r') as f:
        for line in f:
            if "External temperature" in line:
                match = re.search(r'External temperature:\s*([\d\.]+)', line)
                if match:
                    return float(match.group(1))

    raise ValueError("Temperature not found in file.")

def extract_framework_mass(file_path):
    with open(file_path, 'r') as f:
        for line in f:
            if "Framework Mass:" in line:
                match = re.search(r'Framework Mass:\s*([\d\.]+)', line)
                if match:
                    return float(match.group(1))  # g/mol

    raise ValueError("Framework Mass not found in file.")

import re

def extract_components(file_path):
    """
    从RASPA output文件中提取所有气体组分名称
    返回:
        components : list
        例如:
        ["CO2", "N2"]
    """
    components = []
    with open(file_path, "r") as f:
        for line in f:
            # 新版RASPA
            match = re.search(
                r'Component\s+\d+\s+\[(.*?)\]',
                line
            )
            # 兼容旧版
            if match is None:
                match = re.search(
                    r'Component\s+\d+\s+\((.*?)\)',
                    line
                )
            if match:
                comp = match.group(1).strip()
                if comp not in components:
                    components.append(comp)
    if not components:
        raise ValueError("No components found in output file.")
    return components

def fill_component_properties(file_path, data_dict):
    current_comp = None
    with open(file_path, "r") as f:
        for line in f:
            match = re.search(r'Component\s+\d+\s+\[(.*?)\]', line)
            if match is None:
                match = re.search(r'Component\s+\d+\s+\((.*?)\)', line)
            if match:
                current_comp = match.group(1).strip()
                continue
            if current_comp is None:
                continue
            match = re.search(r'MolFraction:\s*([\d\.Ee+-]+)', line)
            if match:
                data_dict["components"][current_comp]["mol_fraction"] = float(match.group(1))
                continue
            match = re.search(r'Partial pressure:\s*([\d\.Ee+-]+)', line)
            if match:
                data_dict["components"][current_comp]["partial_pressure"] = float(match.group(1))
                continue
            match = re.search(r'Partial fugacity:\s*([\d\.Ee+-]+)', line)
            if match:
                data_dict["components"][current_comp]["fugacity"] = float(match.group(1))
                continue
    return data_dict

def initialize_component_container(file_path):
    components = extract_components(file_path)
    data_dict = {
        "cycle": [],
        "U": [],
        "components": {
            comp: {
                "mol_fraction": None,
                "partial_pressure": None,
                "fugacity": None,
                "uptake_FE": None,
                "uptake_RASPA": None,
                "uptake_RASPA_Error": None,
                "Qst_FE": None,
                "Qst_RASPA": None,
                "Qst_RASPA_Error": None,
                "N": []
            }
            for comp in components
        }
    }
    data_dict = fill_component_properties(file_path, data_dict)
    return data_dict

def extract_component_N(line):
    match = re.search(r'Component\s+\d+\s+\((.*?)\).*?molecules:\s*(\d+)/', line)
    if match:
        # print(match.group(1), match.group(2))
        return match.group(1).strip(), int(match.group(2))
    return None, None

def extract_multicomponent_data(file_path, data_dict):
    current_cycle = None
    is_init = False
    frame_data = None
    with open(file_path, 'r') as f:
        for line in f:
            if "Current cycle:" in line:
                match = re.search(r'Current cycle:\s*(\d+)', line)
                if match:
                    current_cycle = int(match.group(1))
                is_init = "[Init]" in line
                frame_data = {comp: None for comp in data_dict["components"]}
            comp, N = extract_component_N(line)
            if comp is not None:
                if comp in frame_data:
                    frame_data[comp] = N
            elif "Current total potential energy:" in line:
                match = re.search(r'Current total potential energy:\s*([-\d\.Ee+]+)', line)
                if match:
                    U = float(match.group(1))
                    if current_cycle is not None and not is_init:
                        data_dict["cycle"].append(current_cycle)
                        data_dict["U"].append(U)
                        for comp in data_dict["components"]:
                            data_dict["components"][comp]["N"].append(frame_data[comp])
    return data_dict

def compute_qst(T, N, U):

    """
    kB: 1.380649*E-23      unit J/K,                   Boltzman constant
    NA: alphogado
    R = kB * NA = 8.314:    J/(mol*K)
    T:              K,                          temperature unit
    NA: 6.02214076*E23, 阿伏伽德罗常数

    U*:             unit K (reduced energy),    RASPA output, 
    Uphys:          unit J (or J/mol depending on scaling), real physical energy
    Upar = kB * U*: unit J,                     single particle
    Umol = kB * NA * U* = R * U*:      unit J/mol,                 per mol

    fluctuation equation:

    RASPA reduced form:
    qst = cov(U*,N)/var(N) - T                 unit K

    Particle form:
    qst = cov(Upar,N)/var(N) - kB*T            unit J

    Molar form:
    qst = cov(Umol,N)/var(N) - R*T             unit J/mol

    cov(U,N)/var(N)=( <U*N> - <U><N>)/(<N*N> - <N><N>)
    """

    N = np.array(N)
    U = np.array(U)  # [K] from RASPA (E/kB), reduced energy

    if len(N) < 2:
        return np.nan

    N_mean = np.mean(N)
    U_mean = np.mean(U)

    cov_UN = np.mean(U * N) - U_mean * N_mean
    var_N  = np.mean(N**2) - N_mean**2

    if var_N < 1e-12:
        return np.nan

    # reduced-unit thermodynamic correction
    qst = (cov_UN / var_N - T)

    # convert K → kJ/mol
    R = 8.314
    qst = qst * R / 1000

    return qst

# =========================
# 3. 计算平均值
# =========================
def compute_average(N_values):
    return np.mean(N_values)


# =========================
# 4. 计算 uptake（mol/kg）
# =========================
def compute_uptake(N_avg, framework_mass_g_per_mol):
    """
    自动单位一致：
    输出：mmol/g
    """
    M_ads_kg = framework_mass_g_per_mol / 1000  # kg/mol
    uptake = N_avg / M_ads_kg  # mol/kg = mmol/g
    return uptake

def computeUptakeAndQst(T, data_dict, M_adsorbent):
# ===== 对每个组分计算 =====
    U = np.array(data_dict["U"])
    for comp, values in data_dict["components"].items():
        N = np.array(values["N"])

        print(comp)
        print("frames =", len(N))
        print("sum =", np.sum(N))
        print("mean =", np.mean(N))
        print("var =", np.var(N))

        N_avg = np.mean(N)
        uptake = compute_uptake(N_avg, M_adsorbent)
        qst = compute_qst(T, N, U)

        data_dict["components"][comp]["uptake_FE"] = uptake
        data_dict["components"][comp]["Qst_FE"] = qst

        print(f"--- {comp} ---")
        print(f"Average N = {N_avg:.5f}")
        print(f"Uptake = {uptake:.6f} mmol/g")
        print(f"Qst = {qst:.3f} kJ/mol")
        print()
    return data_dict

def extractQstAndUptakeFromOutput(file_path, data_dict):
    current_comp = None
    components = list(data_dict["components"].keys())
    with open(file_path, "r") as f:
        for line in f:
            match = re.search(r'Enthalpy of adsorption component\s+\d+\s+\[(.*?)\]', line)
            if match:
                current_comp = match.group(1).strip()
                continue
            if "Total enthalpy of adsorption" in line and len(components) == 1:
                current_comp = components[0]
                continue
            if current_comp is None:
                continue
            match = re.search(r'^\s*(-?\d+\.\d+)\s*\+/-\s*(\d+\.\d+)\s*\[KJ/MOL\]', line)
            if match:
                data_dict["components"][current_comp]["Qst_RASPA"] = abs(float(match.group(1)))
                data_dict["components"][current_comp]["Qst_RASPA_Error"] = float(match.group(2))
                print(f"--- {current_comp} ---")
                print(f"Qst_RASPA = {data_dict["components"][current_comp]["Qst_RASPA"]:.6f} kJ/mol")
                current_comp = None

    current_comp = None

    with open(file_path, "r") as f:
        for line in f:
            match = re.search(r'Component\s+\d+\s+\[(.*?)\]', line)
            if match:
                current_comp = match.group(1).strip()
                continue
            if current_comp is None:
                continue
            match = re.search(r'Average loading absolute \[mol/kg framework\]\s+([\d\.Ee+-]+)\s+\+/-\s+([\d\.Ee+-]+)', line)
            if match:
                data_dict["components"][current_comp]["uptake_RASPA"] = float(match.group(1))
                data_dict["components"][current_comp]["uptake_RASPA_Error"] = float(match.group(2))
                print(f"--- {current_comp} ---")
                print(f"uptake_RASPA = {data_dict["components"][current_comp]["uptake_RASPA"]:.6f} mmol/g")

    return data_dict
# =========================
# 5. 画收敛曲线
# =========================

def plot_components(data_dict):
    plt.figure()
    for comp, data in data_dict["components"].items():
        plt.plot(data_dict["cycle"], data["N"], label=f"{comp}-N", alpha=0.6)
    plt.xlabel("MC cycles")
    plt.ylabel("Number of molecules")
    plt.title("GCMC convergence")
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.show(block=True)


# =========================
# 6. 主程序
# 已兼容单组份和多组分提取和计算
# =========================
def main(file_path=None):
    # ==== 输入参数 ====
    file_path = select_folder()   # 改成你的文件路径

    M_adsorbent = extract_framework_mass(file_path)   # 示例值（必须修改！）
    # ==== 读取数据 ====
    data_dict0 = initialize_component_container(file_path)
    data_dict1 = extract_multicomponent_data(file_path, data_dict0)

    T = extract_temperature(file_path)
    print(f"Temperature = {T} K")

    data_dict2 = computeUptakeAndQst(T, data_dict1, M_adsorbent)
    data_dict3 = extractQstAndUptakeFromOutput(file_path, data_dict2)
    # ===== 对每个组分计算 =====

    # ===== 画图 =====
    plot_components(data_dict3)

    y =0

if __name__ == "__main__":
    f = sys.argv[1] if len(sys.argv) > 1 else None
    main(f)



 