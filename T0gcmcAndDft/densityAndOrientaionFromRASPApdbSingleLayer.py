import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import filedialog
import os
import pandas as pd

def get_safe_sin_theta(theta_centers_rad, theta_min_deg=5):
    """
    返回：
    - sin_theta_safe：避免除0（仅用于计算）
    - valid_theta：物理有效区域mask
    """
    sin_theta = np.sin(theta_centers_rad)

    theta_deg = np.degrees(theta_centers_rad)
    valid_theta = (theta_deg > theta_min_deg) & (theta_deg < 180 - theta_min_deg)

    # 只在计算时避免除0（但不改变物理意义）
    sin_theta_safe = sin_theta.copy()
    sin_theta_safe[~valid_theta] = np.nan  # 或者设为1也可以，但nan更安全

    return sin_theta_safe, valid_theta

def apply_rho_filter(rho, z_theta_counts, valid_theta, min_counts_z=10):
    """
    对 rho(z,θ) 进行统一过滤：
    - θ奇点过滤
    - z统计不足过滤
    """
    rho_clean = rho.copy()

    # --- θ过滤 ---
    rho_clean[:, ~valid_theta] = np.nan

    # --- z过滤 ---
    counts_z = z_theta_counts.sum(axis=1)
    valid_z = counts_z > min_counts_z
    rho_clean[~valid_z, :] = np.nan

    return rho_clean, valid_z
# =========================
# 1️⃣ 弹窗选择文件
# =========================
root = tk.Tk()
root.withdraw()

pdb_file = filedialog.askopenfilename(
    title="Select PDB file",
    filetypes=[("PDB files", "*.pdb"), ("All files", "*.*")]
)

if not pdb_file:
    raise ValueError("No file selected!")

print("Selected file:", pdb_file)

folder = os.path.dirname(pdb_file)
filename = os.path.basename(pdb_file)
print("filename:", filename)
if not filename.startswith("Movie_"):
    raise ValueError("Not a RASPA Movie file")

parts = filename.split('_')

# 解析
len0 = len(parts)
model_name = parts[len0-7]
unitcell = parts[len0-6]
temperature0 = parts[len0-5]
pressure0 = parts[len0-4]
temperature = f"{float(temperature0):.0f}K"
pressure = f"{float(pressure0)/1e5:.1f}bar"   # 如果原来是 Pa

# 提取分子名（component_CO2_0 → CO2）
if "component" in parts:
    comp_index = parts.index("component")
    molecule = parts[comp_index + 1]
else:
    raise ValueError("No component info found!")

print("模型:", model_name)
print("unitcell:", unitcell)
print("温度:", temperature)
print("压强:", pressure)
print("分子:", molecule)

parts = filename.split('_')
# unit部分（第3段）
structure_part = unitcell
# 提取 Z replication
rep_z = int(structure_part.split('.')[-1])
print("Z replication:", rep_z)

output_name = f"distributionDensityAndOrientation_{model_name}_{temperature}_{pressure}_{molecule}.xlsx"
output_file = os.path.join(os.path.dirname(pdb_file), output_name)

# =========================
# 📦 读取 PDB
# =========================
frames = []
box = None

with open(pdb_file, 'r') as f:
    atoms = []
    for line in f:
        if line.startswith("CRYST1"):
            parts = line.split()
            box = np.array([float(parts[1]), float(parts[2]), float(parts[3])])
        elif line.startswith("MODEL"):
            atoms = []
        elif line.startswith("ATOM"):
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            atoms.append([x, y, z])
        elif line.startswith("ENDMDL"):
            frames.append(np.array(atoms))

print(f"读取到 {len(frames)} 帧")

Lx, Ly, Lz = box

Lz_unit = Lz / rep_z

# =========================
# 🔁 PBC 最短距离函数
# =========================
def minimum_image(vec, box):
    return vec - box * np.round(vec / box)

# =========================
# 📊 数据容器
# =========================
z_all = []
cos_all = []
theta_all = []

# =========================
# 🔬 主循环
# =========================
for frame in frames:
    n_atoms = len(frame)

    for i in range(0, n_atoms, 3):
        O1 = frame[i]
        C  = frame[i+1]
        O2 = frame[i+2]

        COM = C.copy()

        # ===== orientation（考虑PBC）
        vec1 = minimum_image(O1 - C, box)
        vec2 = minimum_image(O2 - C, box)
        vec = vec2 - vec1

        norm = np.linalg.norm(vec)
        if norm == 0:
            continue
        
        cos_theta = vec[2] / norm

        # 数值安全（防止浮点误差 >1）
        cos_theta = np.clip(cos_theta, -1.0, 1.0)

        theta = np.degrees(np.arccos(cos_theta))  # 转成角度（°）

        z_folded = COM[2] % Lz_unit
        z_all.append(z_folded)
        cos_all.append(cos_theta)
        theta_all.append(theta)


z_all = np.array(z_all)
cos_all = np.array(cos_all)
theta_all = np.array(theta_all)

# =========================
# 📊 Density profile
# =========================
z_bins = 100   # z方向分箱数
z_hist, edges = np.histogram(z_all, bins=z_bins, range=(0, Lz_unit))
z_bin_centers = 0.5 * (edges[:-1] + edges[1:])

# 体积归一化积（xy面 × dz）
dz = Lz_unit / z_bins
A = Lx * Ly
volume_bin = A * dz
density = z_hist / (volume_bin * rep_z * len(frames))

N_recovered = np.sum(density * volume_bin)
print("Total COM:", len(z_all))
print("Frames:", len(frames))
print("Avg molecules/frame:", len(z_all)/len(frames))
print("每帧分子数:", N_recovered*rep_z)
# 计算直方图（归一化）
counts_normalized, bins_normalized = np.histogram(theta_all, bins=90, range=(0, 180), density=True)
# bin 中心
theta_bin_centers  = 0.5 * (bins_normalized[:-1] + bins_normalized[1:])

# =========================
# 📊 4️⃣ θ histogram（counts）
# =========================
theta_hist_counts, hist_bins = np.histogram(theta_all, bins=90, range=(0, 180))
theta_hist_centers = 0.5 * (hist_bins[:-1] + hist_bins[1:])

#  combined density and orientation distribution
z_theta_counts, z_edges, theta_edges = np.histogram2d(z_all, theta_all, bins=[z_bins, 90], range=[[0, Lz_unit], [0, 180]])

dtheta = np.radians(180 / 90)   # 必须用弧度！
theta_centers = 0.5 * (theta_edges[:-1] + theta_edges[1:])
theta_centers_rad = np.radians(theta_centers)
sin_theta = np.sin(theta_centers_rad)
sin_theta_safe, valid_theta = get_safe_sin_theta(theta_centers_rad,theta_min_deg=0)

rho_rad = z_theta_counts / (A * dz * dtheta * sin_theta_safe[np.newaxis, :] * rep_z  * len(frames))
rho_clean, valid_z = apply_rho_filter(rho_rad, z_theta_counts, valid_theta, min_counts_z=10)

rho_deg = rho_rad * (np.pi / 180.0)
denominator = np.sum(rho_rad * dtheta * sin_theta[np.newaxis, :], axis=1)
rho_norm = np.divide(rho_rad, denominator[:, None], out=np.zeros_like(rho_rad), where=denominator[:, None] > 1e-12)


def save_data_V2(x, y, x_name="", y_name=""):
    df = pd.DataFrame({
    x_name: x,
    y_name: y
    })  
    return df


# 📊 1️⃣ Density 数据
df_density = save_data_V2(z_bin_centers, density, x_name="z (Å)", y_name="Density (1/Å^3)")
# 📊 2️⃣ Orientation scatter（cosθ）
df_cos = save_data_V2(z_all, cos_all, x_name="z (Å)", y_name="cos(theta)")
# 📊 3️⃣ Orientation scatter（θ）
df_theta_z = save_data_V2(z_all, theta_all, x_name="z (Å)", y_name="theta (degree)")
df_theta_hist = save_data_V2(theta_hist_centers, theta_hist_counts, x_name="theta (degree)", y_name="Counts")
# =========================
# 📊 5️⃣ θ normalized distribution
# =========================
df_theta_norm = save_data_V2(theta_bin_centers, counts_normalized, x_name="theta (degree)", y_name="Normalized Frequency")

# 📊 7️⃣ ρ(z,θ) 长表
Z, THETA = np.meshgrid(z_bin_centers, theta_centers, indexing='ij')
df_rho_z_theta = pd.DataFrame({"z (Å)": Z.flatten(), "theta (degree)": THETA.flatten(), "rho (1/Å^3·degree)": rho_deg.flatten()})

# 💾 写入 Excel
with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
    df_density.to_excel(writer, sheet_name="Density", index=False)
    df_cos.to_excel(writer, sheet_name="cosTheta_vs_z", index=False)
    df_theta_z.to_excel(writer, sheet_name="theta_vs_z", index=False)
    df_theta_hist.to_excel(writer, sheet_name="theta_histogram", index=False)
    df_theta_norm.to_excel(writer, sheet_name="theta_normalized", index=False)
    df_rho_z_theta.to_excel(writer, sheet_name="rho_z_theta", index=False)
print(f"数据已保存到 {output_file}")

from mpl_toolkits.mplot3d import Axes3D

def plot_rho_surface(rho, z_vals, theta_vals):
    Z, TH = np.meshgrid(z_vals, theta_vals, indexing='ij')

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    surf = ax.plot_surface(
        Z, TH, rho,
        cmap='viridis',
        linewidth=0,
        antialiased=True
    )

    ax.set_xlabel("z (Å)")
    ax.set_ylabel("theta (degree)")
    ax.set_zlabel("ρ (1/Å³·degree)")
    ax.set_title("3D Surface of ρ(z,θ)")

    fig.colorbar(surf, shrink=0.5, aspect=10)

    plt.show()

def plot_rho_heatmap(rho, z_range, theta_range):
    plt.figure()
    plt.imshow(
        rho.T,
        extent=[z_range[0], z_range[1], theta_range[0], theta_range[1]],
        aspect='auto',
        origin='lower',
        interpolation='nearest'   # 👈 关键
    )
    plt.xlabel("z (Å)")
    plt.ylabel("theta (degree)")
    plt.title("Joint Density-Orientation Distribution ρ(z,θ)")
    plt.colorbar(label="ρ (1/Å³·degree)")

def plot_hist2d(x, y, bins, xlabel, ylabel, title, cbar_label="Counts"):
    plt.figure()
    plt.hist2d(x, y, bins=bins)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.colorbar(label=cbar_label)

def plot_hist(data, bins, xlabel, ylabel, title):
    plt.figure()
    plt.hist(data, bins=bins)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)

def plot_scatter(x, y, xlabel, ylabel, title):
    plt.figure()
    plt.scatter(x, y, s=2, alpha=0.3)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid()

def plot_line(x, y, xlabel, ylabel, title, xlim=None):
    plt.figure()
    plt.plot(x, y)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    if xlim:
        plt.xlim(xlim)
    plt.grid()
# =========================
# 📈 作图
# =========================
plt.ion()

# 🔹 Density
plot_line(z_bin_centers, density,
          "z (Å)", "Density (1/Å³)", "Density Profile")

# 🔹 Orientation scatter
plot_scatter(z_all, cos_all,
             "z (Å)", "cos(theta)", "Orientation Distribution")

plot_scatter(z_all, theta_all,
             "z (Å)", "theta (degree)", "Orientation Angle Distribution")

# 🔹 2D histogram
plot_hist2d(z_all, cos_all, [z_bins, 50],
            "z (Å)", "cos(theta)", "z vs Orientation")

plot_hist2d(z_all, theta_all, [100, 90],
            "z (Å)", "theta (degree)", "z vs theta")

# 🔹 θ histogram
plot_hist(theta_all, 90,
          "theta (degree)", "Counts", "Orientation Angle Histogram")

# 🔹 θ distribution（归一化）
plot_line(theta_bin_centers, counts_normalized,
          "Angle θ (degree)", "Normalized Frequency",
          "Orientation Angle Distribution", xlim=(0, 180))

# 🔹 ρ(z,θ)
plot_rho_heatmap(rho_deg, [0, Lz_unit], [0, 180])
plot_rho_surface(rho_norm, z_bin_centers, theta_centers)

plt.show(block=True)

