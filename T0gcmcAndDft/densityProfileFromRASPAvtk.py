import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import filedialog

# =========================
# 1️⃣ 弹窗选择文件
# =========================
root = tk.Tk()
root.withdraw()

file_path = filedialog.askopenfilename(
    title="Select VTK file",
    filetypes=[("VTK files", "*.vtk"), ("All files", "*.*")]
)

if not file_path:
    raise ValueError("No file selected!")

print("Selected file:", file_path)

# =========================
# 2️⃣ 读取文件 + 自动解析结构
# =========================
nx = ny = nz = None
spacing = None
origin = None

data = []
start_reading = False

with open(file_path, 'r') as f:
    for line in f:
        line = line.strip()

        # --- 解析结构信息 ---
        if line.startswith("DIMENSIONS"):
            parts = line.split()
            nx, ny, nz = map(int, parts[1:4])

        elif line.startswith("SPACING"):
            spacing = list(map(float, line.split()[1:4]))

        elif line.startswith("ORIGIN"):
            origin = list(map(float, line.split()[1:4]))

        elif line.startswith("CELL_PARAMETERS"):
            cell = list(map(float, line.split()[1:7]))
            print("Cell parameters:", cell)

        # --- 找到数据开始 ---
        elif line.startswith("LOOKUP_TABLE"):
            start_reading = True
            continue

        # --- 读取数据 ---
        elif start_reading:
            if line == "":
                continue
            try:
                data.append(float(line))
            except:
                continue

# 转 numpy
data = np.array(data)

# =========================
# 3️⃣ 检查 & reshape
# =========================
print("\n--- Grid Info ---")
print(f"DIMENSIONS: {nx} x {ny} x {nz}")
print(f"Total expected: {nx*ny*nz}")
print(f"Actual data: {len(data)}")

if len(data) != nx * ny * nz:
    raise ValueError("❌ Data size mismatch! Check VTK file.")

grid = data.reshape((nz, ny, nx))   # 注意顺序
grid = np.transpose(grid, (2, 1, 0))  # 转成 (nx, ny, nz)
thicknessOfInterlayer = 3.35
layer = 2
thick = layer * thicknessOfInterlayer

print("✅ Grid loaded successfully")
print("Max density:", grid.max())
print("Min density:", grid.min())

# =========================
# 4️⃣ 2D slice（XY平面）
# =========================
z_slice = nz // 2

plt.figure()
plt.imshow(grid[:, :, z_slice])
plt.colorbar(label="Density")
plt.title(f"2D Slice (z={z_slice})")
plt.xlabel("X")
plt.ylabel("Y")
plt.tight_layout()
plt.show()

# =========================
# 5️⃣ 1D density profile（核心🔥）
# =========================
rho_z = grid.mean(axis=(0, 1))

# 用真实坐标（如果有spacing）
if spacing:
    z = np.arange(nz) * spacing[2] - thick
    xlabel = "Z (Å)"
else:
    z = np.arange(nz)
    xlabel = "Z index"

plt.figure()
plt.plot(z, rho_z)
plt.xlabel(xlabel)
plt.ylabel("Density")
plt.title("1D Density Profile")
plt.grid()
plt.tight_layout()
plt.show()

# =========================
# 6️⃣ log图（处理大量0）
# =========================
rho_z_safe = np.where(rho_z > 0, rho_z, 1e-12)

plt.figure()
plt.plot(z, np.log10(rho_z_safe))
plt.xlabel(xlabel)
plt.ylabel("log10(Density)")
plt.title("Log Density Profile")
plt.grid()
plt.tight_layout()
plt.show()

# =========================
# 7️⃣ 3D scatter（简单可视化）
# =========================
# 保留整数索引（用于取值）
threshold = np.percentile(grid, 99)
x_idx, y_idx, z_idx = np.where(grid > threshold)

# 转换成真实坐标（用于画图）
if spacing:
    x = x_idx * spacing[0]
    y = y_idx * spacing[1]
    z_real = z_idx * spacing[2]
else:
    x, y, z_real = x_idx, y_idx, z_idx

# 用整数索引取 density
values = grid[x_idx, y_idx, z_idx]

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

sc = ax.scatter(x, y, z_real, c=values, s=5)

plt.colorbar(sc, label="Density")

ax.set_xlabel("X (Å)")
ax.set_ylabel("Y (Å)")
ax.set_zlabel("Z (Å)")
ax.set_title("3D Density Scatter")

plt.show()