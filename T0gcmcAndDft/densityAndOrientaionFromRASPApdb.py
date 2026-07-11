# =========================================================
# Component configuration
# =========================================================
def get_component_config(component):
    configs = {
        "CO2": {
            "n_atoms": 3,
            "atom_order": ["O", "C", "O"],
            "com_index": 1,
            "orientation_atoms": (0, 2),
            "has_orientation": True
        },
        "N2": {
            "n_atoms": 2,
            "atom_order": ["N", "N"],
            "com_index": None,
            "orientation_atoms": (0, 1),
            "has_orientation": True
        },
        "CH4": {
            "n_atoms": 5,
            "atom_order": ["C", "H", "H", "H", "H"],
            "com_index": 0,
            "orientation_atoms": None,
            "has_orientation": False
        },
        "Ar": {
            "n_atoms": 1,
            "atom_order": ["Ar"],
            "com_index": 0,
            "orientation_atoms": None,
            "has_orientation": False
        }
    }
    if component not in configs:
        raise ValueError(f"Unsupported component: {component}")
    return configs[component]

# =========================================================
# Utility functions
# =========================================================
def minimum_image(vec, box):
    return vec - box * np.round(vec / box)

def unwrap_molecule(atoms, box):
    ref = atoms[0]
    unwrapped = [ref]
    for atom in atoms[1:]:
        vec = minimum_image(atom - ref, box)
        unwrapped.append(ref + vec)
    return np.array(unwrapped)

def get_safe_sin_theta(theta_centers_rad, eps=1e-6):
    sin_theta = np.sin(theta_centers_rad)
    sin_theta_safe = np.maximum(sin_theta, eps)
    valid_theta = np.ones_like(sin_theta, dtype=bool)
    return sin_theta_safe, valid_theta

def apply_rho_filter0(rho, z_theta_counts, valid_theta, min_percent=0.005):
    rho_clean = rho.copy()
    rho_clean[:, ~valid_theta] = np.nan
    total_counts = z_theta_counts.sum()
    threshold = max(2, int(np.ceil(total_counts * min_percent / 100.0)))
    low_count_mask = z_theta_counts < threshold
    removed = np.sum(low_count_mask)
    total = low_count_mask.size
    rho_clean[low_count_mask]  *= 0.1
    print(f"Total counts = {total_counts}")
    print(f"Threshold = {threshold:.2f}")
    print(f"Remove {removed}/{total} bins")
    # rho_clean = gaussian_filter(rho_clean,sigma=1.0)
    return rho_clean

def apply_rho_filter( rho, z_theta_counts, valid_theta,  min_percent=0.01):
    """
    对 rho(z,θ) 进行统一过滤：
    - θ奇点过滤
    - z统计不足过滤
    """
    rho_clean = rho.copy()
    rho_clean[:, ~valid_theta] = np.nan
    counts_z = z_theta_counts.sum(axis=1)
    total_counts = counts_z.sum()
    threshold = total_counts * min_percent / 100.0
    valid_z = counts_z > threshold
    rho_clean[~valid_z, :] = np.nan
    print("Min counts per z:", counts_z.min())
    print("Max counts per z:", counts_z.max())
    print("Mean counts per z:", counts_z.mean())
    print("valid z bins =", np.sum(counts_z > 100))

    for i in np.where(valid_z)[0]:
        print(
            "z =", i,
            "total =", z_theta_counts[i,:].sum(),
            "max_theta =", np.max(z_theta_counts[i,:]),
            "nonzero_theta_bins =", np.sum(z_theta_counts[i,:] > 0)
        )
    return rho_clean


def save_data_V2(x, y, x_name="", y_name=""):
    return pd.DataFrame({
        x_name: x,
        y_name: y
    })


# =========================================================
# Plot functions
# =========================================================
def plot_rho_surface(rho, z_vals, theta_vals):
    Z, TH = np.meshgrid( z_vals, theta_vals, indexing='ij')
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    surf = ax.plot_surface( Z, TH, rho, cmap='viridis', linewidth=0, antialiased=True )
    ax.set_xlabel("z (Å)")
    ax.set_ylabel("theta (degree)")
    ax.set_zlabel("ρ")
    ax.set_title("Orientation distribution conditioned on molecular position ρ(θ|z)")
    fig.colorbar(surf, shrink=0.5, aspect=10)
    plt.show()

def plot_rho_heatmap(rho, z_range, theta_range):
    plt.figure()
    plt.imshow(
        rho.T,
        extent=[
            z_range[0],
            z_range[1],
            theta_range[0],
            theta_range[1]
        ],
        aspect='auto',
        origin='lower',
        interpolation='nearest'
    )
    plt.xlabel("z (Å)")
    plt.ylabel("theta (degree)")
    plt.title("Joint Density-Orientation Distribution")
    plt.colorbar(label="ρ")

def plot_hist2d(  x, y, bins, xlabel, ylabel, title, cbar_label="Counts"):
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
    plt.scatter( x, y, s=2, alpha=0.3)
    plt.xlabel(xlabel, fontsize=24)
    plt.ylabel(ylabel, fontsize=24)
    plt.title(title, fontsize=24)
    plt.grid()


def plot_line(x, y, xlabel, ylabel, title, xlim=None):
    plt.figure()
    plt.plot(x, y)
    plt.xlabel(xlabel, fontsize=18)
    plt.ylabel(ylabel, fontsize=18)
    plt.title(title, fontsize=18)
    plt.tick_params(axis='both', labelsize=13)
    if xlim:
        plt.xlim(xlim)
    plt.grid()
    plt.tight_layout()

# =========================================================
# Select PDB file
# =========================================================
def loadFile():
    root = tk.Tk()
    root.withdraw()
    pdb_file = filedialog.askopenfilename(
        title="Select PDB file",
        filetypes=[
            ("PDB files", "*.pdb"),
            ("All files", "*.*")
        ]
    )
    if not pdb_file:
        raise ValueError("No file selected!")
    print("Selected file:", pdb_file)

# =========================================================
# Parse filename
# =========================================================
    folder = os.path.dirname(pdb_file)
    filename = os.path.basename(pdb_file)
    print("filename:", filename)
    parts = filename.replace(".pdb", "").split("_")

    # 1. 找 unitcell
    for i, item in enumerate(parts):
        if re.fullmatch(r"\d+\.\d+\.\d+", item):
            unitcell_idx = i
            break
    else:
        raise ValueError("Cannot find unitcell")

    # 2. model
    model_name = "_".join(parts[1:unitcell_idx])
    
    layerNumber = 1
    for part in parts:
        m = re.match(r"(\d+)[Ll]ayer", part)
        if m:
            layerNumber = int(m.group(1))
            break

    print("Layer number:", layerNumber)

    # 3. unitcell
    unitcell = parts[unitcell_idx]

    # 4. T、P（旧格式）
    temperature0 = parts[unitcell_idx + 1]
    pressure0 = parts[unitcell_idx + 2]
    temperature = f"{float(temperature0):.0f}K"
    pressure = f"{float(pressure0)/1e5:.3f}bar"

    # 5. molecule
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

    structure_part = unitcell
    rep_z = int(structure_part.split('.')[-1])
    print("Z replication:", rep_z)

    output_name = (
        f"distributionDensityAndOrientation_"
        f"{model_name}_{temperature}_{pressure}_{molecule}.xlsx"
    )

    output_file = os.path.join( os.path.dirname(pdb_file), output_name)
    return molecule, layerNumber, rep_z, pdb_file, output_file

# =========================================================
# Read PDB
# =========================================================
def extractData(molecule,rep_z,pdb_file,wallThickness):
    frames=[]
    box=None
    with open(pdb_file,'r') as f:
        atoms=[]
        for line in f:
            if line.startswith("CRYST1"):
                parts=line.split()
                box=np.array([
                    float(parts[1]),
                    float(parts[2]),
                    float(parts[3])
                ])
            elif line.startswith("MODEL"):
                atoms=[]
            elif line.startswith("ATOM"):
                # print(repr(line))
                x=float(line[30:38])
                y=float(line[38:46])
                z=float(line[46:54])
                atoms.append([x,y,z])
            elif line.startswith("ENDMDL"):
                frames.append(np.array(atoms))
    print(f"读取到 {len(frames)} 帧")

    box_size=box
    poreWidth = box_size[2] / rep_z - wallThickness
    Lz_unit=poreWidth + wallThickness
    if abs(Lz_unit-box_size[2]/rep_z)>1e-3:
        raise ValueError(
            f"Inconsistent geometry: poreWidth+wallThickness={Lz_unit:.4f}, "
            f"box_size[2]/rep_z={box_size[2]/rep_z:.4f}"
        )

    config=get_component_config(molecule)
    n_atoms_mol=config["n_atoms"]
    com_index=config["com_index"]
    orientation_atoms=config["orientation_atoms"]
    has_orientation=config["has_orientation"]

    z_all=[]
    z_orient_all=[]
    cos_all=[]
    theta_all=[]

    frams_num=len(frames)
    for frame in frames:
        n_atoms=len(frame)
        for i in range(0,n_atoms,n_atoms_mol):
            atoms=frame[i:i+n_atoms_mol]
            if len(atoms)<n_atoms_mol:
                continue
            atoms_unwrapped=unwrap_molecule(atoms,box)
            if com_index is not None:
                COM=atoms_unwrapped[com_index].copy()
            else:
                COM=atoms_unwrapped.mean(axis=0)

            z_pore=(COM[2]-wallThickness)%Lz_unit   #折叠到单周期空间
            if z_pore>poreWidth:
                continue
            z_all.append(z_pore)

            if has_orientation:
                a1,a2=orientation_atoms
                vec=minimum_image( atoms_unwrapped[a2]-atoms_unwrapped[a1], box )
                norm=np.linalg.norm(vec)
                if norm<1e-12:
                    continue

                cos_theta=vec[2]/norm
                cos_theta=np.clip(cos_theta,-1.0,1.0)
                theta=np.degrees(np.arccos(cos_theta))
                z_orient_all.append(z_pore)
                cos_all.append(cos_theta)
                theta_all.append(theta)

    z_all=np.array(z_all)
    z_orient_all=np.array(z_orient_all)
    cos_all=np.array(cos_all)
    theta_all=np.array(theta_all)
    return z_all,z_orient_all,cos_all,theta_all,box_size,frams_num, poreWidth

# =========================================================
# Density profile
# =========================================================
def densityFun(z_all, box_size, poreWidth, rep_z, frams_num, z_bins = 100):
    z_hist, edges = np.histogram(z_all, bins=z_bins, range=(0, poreWidth))
    z_bin_centers = 0.5 * (edges[:-1] + edges[1:])
    dz = poreWidth / z_bins
    A = box_size[0] * box_size[1]
    volume_bin = A * dz
    density_z = z_hist / ( volume_bin * rep_z * frams_num)
    N_recovered = np.sum( density_z * volume_bin)

    print("Total COM:", len(z_all))
    print("Frames:", frams_num)
    print("Avg molecules/frame:", len(z_all)/frams_num)
    print("每帧分子数:", N_recovered * rep_z)
    return z_bin_centers, density_z
# =========================================================
# Orientation statistics
# =========================================================
def orientationFun(has_orientation, theta_all, z_orient_all, frams_num, poreWidth, rep_z, box_size, z_bins = 100, angle_bins = 90):
    
    if has_orientation and len(theta_all) > 0:
        theta_normalized0, bins_normalized = np.histogram(theta_all, angle_bins, range=(0, 180), density=True)
        theta_normalized = gaussian_filter1d(theta_normalized0, sigma=1.5)
        theta_norm_bin_centers = 0.5 * (bins_normalized[:-1] + bins_normalized[1:])
        theta_histogram, hist_bins = np.histogram(theta_all, angle_bins, range=(0, 180))
        theta_hist_centers = 0.5 * (hist_bins[:-1] + hist_bins[1:])
    else:
        theta_normalized = None
        theta_norm_bin_centers = None
        theta_histogram = None
        theta_hist_centers = None

    # =========================================================
    # Joint distribution
    # =========================================================

    if has_orientation and len(theta_all) > 0:
        z_theta_counts, z_edges, theta_edges = np.histogram2d(z_orient_all, theta_all, bins=[z_bins, angle_bins], range=[[0, poreWidth], [0, 180]])
        dtheta = np.radians(180 / 90)
        theta_centers = 0.5 * (theta_edges[:-1] + theta_edges[1:])
        theta_centers_rad = np.radians(theta_centers)
        sin_theta = np.sin(theta_centers_rad)
        sin_theta_safe, valid_theta = get_safe_sin_theta(theta_centers_rad, eps=1e-6)
        dz = poreWidth / z_bins
        A = box_size[0] * box_size[1]
        volume_bin = A * dz
        rho_rad0 = z_theta_counts / (volume_bin * dtheta * sin_theta_safe[np.newaxis, :] * rep_z * frams_num)
        rho_rad = apply_rho_filter0(rho_rad0, z_theta_counts, valid_theta, min_percent=0.005)
        
        rho_zAndDeg = rho_rad * (np.pi / 180.0) #联合概率密度 P( theta, Z )
        denominator = np.sum(rho_rad * dtheta * sin_theta[np.newaxis, :], axis=1)
        rho_degAtZ = np.divide(rho_rad, denominator[:, None], out=np.zeros_like(rho_rad), where=denominator[:, None] > 1e-12) #条件概率密度 P( theta|Z )
    else:
        rho_rad = None
        rho_zAndDeg = None
        rho_degAtZ = None
        theta_centers = None
    return theta_normalized, theta_norm_bin_centers, theta_histogram, theta_hist_centers, theta_centers, rho_zAndDeg, rho_degAtZ

# =========================================================
# Save dataframes
# =========================================================
def saveData(output_file, z_bin_centers, density, has_orientation, theta_all, theta_centers, rho_zAndDeg, z_orient_all, cos_all, theta_hist_centers, theta_histogram, theta_norm_bin_centers, theta_normalized):
    df_density = save_data_V2( z_bin_centers, density, x_name="z (Å)", y_name="Density (1/Å^3)")

    if has_orientation and len(theta_all) > 0:
        df_cos = save_data_V2(z_orient_all, cos_all, x_name="z (Å)", y_name="cos(theta)")
        df_theta_z = save_data_V2(z_orient_all, theta_all, x_name="z (Å)", y_name="theta (degree)")
        df_theta_hist = save_data_V2(theta_hist_centers, theta_histogram, x_name="theta (degree)", y_name="Counts")
        df_theta_norm = save_data_V2(theta_norm_bin_centers, theta_normalized, x_name="theta (degree)", y_name="Normalized Frequency")
        Z, THETA = np.meshgrid(z_bin_centers, theta_centers, indexing='ij')
        df_rho_z_theta = pd.DataFrame({"z (Å)": Z.flatten(), "theta (degree)": THETA.flatten(), "rho (1/Å^3·degree)": rho_zAndDeg.flatten()})

    data_to_save = [
        (df_cos,         "cosTheta_vs_z"),
        (df_theta_z,     "theta_vs_z"),
        (df_theta_hist,  "theta_histogram"),
        (df_theta_norm,  "theta_normalized"),
        (df_rho_z_theta, "rho_z_theta"),
    ]
    # =========================================================
    # Save Excel
    # =========================================================
    fl.export_to_excel_auto(df_density, filename=output_file, sheet_name="Density")
    if has_orientation and len(theta_all) > 0:
        for df, sheet in data_to_save:
            fl.export_to_excel_auto( df, filename=output_file, sheet_name=sheet, )

    print(f"数据已保存到 {output_file}")


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

def main(file_path=None):
    #z_all : 分子坐标
    #z_orient_all：带方向的分子的坐标, CH4为无方向分子
    #density_z：局部密度分布，local density distribution
    #theta_normalized: 角度的概率分布密度函数
    #theta_norm_bin_centers：角度的概率分布密度函数对应的每一个小块的中心
    #theta_histogram：角度分布统计
    #theta_hist_centers：角度分布统计对应的每一个小块的中心
    #rho_zAndDeg：#联合概率密度 P( theta, Z )
    #rho_degAtZ：#条件概率密度 P( theta|Z )

    ##bins 采样分辨率
    z_bins = 200
    angle_bins = 90 
    thicknessOfInterlayers = 3.35
    
    molecule, layerNumber, rep_z, pdb_file, output_file = loadFile()
    wallThickness = thicknessOfInterlayers * (layerNumber - 1)
    config = get_component_config(molecule)
    has_orientation = config["has_orientation"]

    z_all,z_orient_all,cos_all,theta_all,box_size,frams_num, poreWidth = extractData(molecule,rep_z,pdb_file, wallThickness)
    z_bin_centers, density_z = densityFun(z_all, box_size, poreWidth, rep_z, frams_num, z_bins)
    theta_normalized, theta_norm_bin_centers, theta_histogram, theta_hist_centers, theta_centers, rho_zAndDeg, rho_degAtZ = orientationFun(has_orientation, theta_all, z_orient_all, frams_num, poreWidth, rep_z, box_size, z_bins, angle_bins)
    saveData(output_file, z_bin_centers, density_z, has_orientation, theta_all, theta_centers, rho_zAndDeg, z_orient_all, cos_all, theta_hist_centers, theta_histogram, theta_norm_bin_centers, theta_normalized)

    # =========================================================
    # Plot
    # =========================================================
    z_norm = z_bin_centers / poreWidth
    density_kmol = density_z * 1660.54
    plt.ion()
    plot_line( z_bin_centers, density_z, "z (Å)", "Density (molecule/Å³)", "Density Profile")
    plot_line( z_norm, density_kmol, "Normaralized Z (Å)", "Density (kmol/m³)", "Density Profile")
    if has_orientation and len(theta_all) > 0:
        # plot_scatter( z_orient_all, cos_all, "z (Å)", "cos(theta)", "Orientation Distribution" )
        plot_scatter( z_orient_all, theta_all, "z (Å)",  "theta (degree)", "Orientation Angle Distribution")
        # plot_hist2d( z_orient_all, cos_all, [z_bins, 50], "z (Å)", "cos(theta)", "z vs Orientation" )
        plot_hist2d( z_orient_all, theta_all, [z_bins, angle_bins], "z (Å)", "theta (degree)", "z vs theta" )
        # plot_hist( theta_all, angle_bins, "theta (degree)", "Counts", "Orientation Angle Histogram" )
        plot_line( theta_norm_bin_centers, theta_normalized, "Angle θ (degree)", "Normalized Frequency", "Orientation Angle Distribution", xlim=(0, 180) )
        rho_zAndDeg0 = rho_zAndDeg / np.max(rho_zAndDeg)
        plot_rho_heatmap( rho_zAndDeg0, [0, poreWidth], [0, 180] )
        plot_rho_heatmap( rho_degAtZ, [0, poreWidth], [0, 180] )
        # plot_rho_surface( rho_degAtZ, z_bin_centers, theta_centers )
    plt.show(block=True)

if __name__ == "__main__":
    f = sys.argv[1] if len(sys.argv) > 1 else None
    main(f)