from collections import defaultdict
import re
import sys
import os
import tkinter as tk
from tkinter import filedialog
import glob
import re
import numpy as np
import matplotlib.pyplot as plt

def select_folder():
    root = tk.Tk()
    root.withdraw()
    return filedialog.askdirectory(title="请选择包含 output 文件的文件夹")

def select_folder0():
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

def extract_multicomponent_data(file_path):
    cycles = []
    data_dict = defaultdict(lambda: {"N": [], "U": []})

    current_cycle = None
    is_init = False
    current_components = []
    current_U = None

    with open(file_path, 'r') as f:
        for line in f:

            # ===== cycle =====
            if "Current cycle:" in line:
                match = re.search(r'Current cycle:\s*(\d+)', line)
                if match:
                    current_cycle = int(match.group(1))
                is_init = "[Init]" in line
                current_components = []

            # ===== component =====
            elif "Component" in line and "(" in line:
                match = re.search(r'Component\s+\d+\s+\((.*?)\)', line)
                if match:
                    current_components.append(match.group(1))

            # ===== N =====
            elif "Number of Adsorbates:" in line:
                match = re.search(r'Number of Adsorbates:\s*(\d+)', line)
                if match and current_cycle is not None and not is_init:
                    N = int(match.group(1))

            # ===== U（关键新增）=====
            elif "Current total potential energy:" in line:
                match = re.search(r'Current total potential energy:\s*([-\d\.]+)', line)
                if match:
                    current_U = float(match.group(1))

                    # 记录（必须在同一个cycle里）
                    if current_cycle is not None and not is_init:
                        cycles.append(current_cycle)

                        if len(current_components) == 1:
                            comp = current_components[0]
                            data_dict[comp]["N"].append(N)
                            data_dict[comp]["U"].append(current_U)

    return np.array(cycles), data_dict

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

    if var_N == 0:
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

# =========================
# 5. 画收敛曲线
# =========================

def plot_components(cycles, data_dict):
    plt.figure()

    for comp, data in data_dict.items():
        plt.plot(cycles, data["N"], label=f"{comp}-N", alpha=0.6)
        # plt.plot(cycles, data["U"], label=f"{comp}-U", linestyle="--", alpha=0.6)

    plt.xlabel("MC cycles")
    plt.ylabel("Value")
    plt.title("GCMC convergence (N and U)")
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.show(block=True)

import pandas as pd
import extractDataFromRaspaOutput as EDF
EXTRACT_KEYS_VOLUME = {
    "unit_cell": {
        "keyword": "Unit cell size:",
        "type": "vector3",
        "scope": "global"
    },
    "pressure (Pa)": {
        "keyword": "Partial pressure",
        "type": "scalar",
        "scope": "molecule"
    }
}
# =========================
# 6. 主程序
# =========================
def main(file_path=None):
    all_results = []

# ==== 输入参数 ====
    folder = select_folder()   # 改成你的文件路径
    energy_results = EDF.find_all_outputs( folder, EXTRACT_KEYS_VOLUME )
    energy_dict = { item["path"]: item for item in energy_results } #先按照path建立索引

    for root, dirs, files in os.walk(folder):
        for file in files:
            if file.startswith("output_") and file.endswith(".data"):
                full_path = os.path.join(root, file)
                print(f"\nProcessing: {full_path}")
                # file_name = os.path.basename(file_path)
                M_adsorbent = extract_framework_mass(full_path)
                cycles, data_dict = extract_multicomponent_data(full_path)
        
    # ===== 对每个组分计算 =====
                T = extract_temperature(full_path)
                print(f"Temperature = {T} K")

                for comp, values in data_dict.items():  #依据adsorbate个数来循环
                    N = np.array(values["N"])
                    U = np.array(values["U"])

                    N_avg = np.mean(N)
                    uptake = compute_uptake(N_avg, M_adsorbent)

                    qst = compute_qst(T, N, U)

                    result = energy_dict[full_path].copy()

                    result.update({
                        "component": comp,
                        "avg_number": N_avg,
                        "T(K)": T,
                        "uptake(mmol/g)": uptake,
                        "Qst(kJ/mol)": qst
                    })

                    all_results.append(result)

                # 示例值（必须修改！）
                # ==== 读取数据 ====

                    # print(f"--- {comp} ---")
                    # print(f"Average N = {N_avg:.5f}")
                    # print(f"Uptake = {uptake:.6f} mmol/g")
                    # print(f"qst = {qst:.3f} kJ/mol\n")

                    # # ===== 画图 =====
                    # plot_components(cycles, data_dict)
    df = pd.DataFrame(all_results)
    file="summary_num.csv"
    full_path = os.path.join(folder, file)
    df.to_csv( full_path, index=False, encoding="utf-8-sig" )
    print(f"\nProcessed {len(all_results)} files.")

if __name__ == "__main__":
    f = sys.argv[1] if len(sys.argv) > 1 else None
    main(f)



 