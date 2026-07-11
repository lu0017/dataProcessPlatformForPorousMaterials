from ase import Atoms
import numpy as np

# =========================
# Basic graphene geometry
# =========================
c2c = 1.42                         # C-C bond length (Å)
a = c2c * np.sqrt(3)              # lattice parameter

cell = [
    [a, 0, 0],
    [a/2, a*np.sqrt(3)/2, 0],
    [0, 0, 20.0]                  # temporary z cell
]

basis = [
    (0, 0, 0),
    (1/3, 1/3, 0)
]

# =========================
# Replication
# =========================
nx, ny = 1, 1

positions = []
symbols = []

for i in range(nx):
    for j in range(ny):
        for b in basis:
            r = (i + b[0]) * np.array(cell[0]) + (j + b[1]) * np.array(cell[1])
            positions.append(r)
            symbols.append("C")

# =========================
# Slit-pore settings
# =========================
d_slit = 10          # free pore width between inner surfaces (Å)
file = "3layer_LuSlitGraphene100.cif"

d_graphite = 3.35     # graphite interlayer spacing (Å)
# each wall has 3 graphene sheets:
# left wall  : z1 z2 z3
# right wall : z4 z5 z6

# total z height
c_vacuum = 0.0

# =========================
# Layer z positions
# =========================
z1 = c_vacuum
z2 = z1 + d_graphite
z3 = z2 + d_graphite
wall_thickness = 2 * d_graphite          # 3 layers = two spacings
# c_total = wall_thickness + d_slit + wall_thickness + 2*c_vacuum
c_total = wall_thickness + d_slit + c_vacuum

# z4 = z3 + d_slit
# z5 = z4 + d_graphite
# z6 = z5 + d_graphite

# z_list = [z1, z2, z3, z4, z5, z6]
z_list = [z1, z2, z3]

# =========================
# Build layers
# =========================
structure = Atoms()

for z in z_list:
    layer = Atoms(
        symbols=symbols,
        positions=positions,
        cell=cell,
        pbc=[True, True, False]
    )
    layer.translate([0, 0, z])
    structure += layer

# =========================
# Set final cell
# =========================
supercell = [
    np.array(cell[0]) * nx,
    np.array(cell[1]) * ny,
    np.array([0, 0, c_total])
]

structure.set_cell(supercell)
structure.set_pbc([True, True, True])

# =========================
# Write CIF
# =========================
structure.write(file)

# =========================
# Auto-fix CIF symmetry
# =========================
with open(file, "r") as f:
    lines = f.readlines()

new_lines = []
inserted = False

for line in lines:
    new_lines.append(line)
    if line.strip().startswith("data_") and not inserted:
        new_lines.append("_symmetry_space_group_name_H-M    'P1'\n")
        new_lines.append("_symmetry_Int_Tables_number       1\n")
        inserted = True

with open(file, "w") as f:
    f.writelines(new_lines)

print("Saved:", file)
print("Total carbon atoms:", len(structure))
print("Cell height (Å):", c_total)
print("Pore width (Å):", d_slit)
print("Each wall = 3 graphene layers")