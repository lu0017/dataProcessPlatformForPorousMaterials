from ase import Atoms
import numpy as np
import sys

def buildLattice(c2c):
    a = c2c * np.sqrt(3)

    cell = [
        [a, 0, 0],
        [a/2, a*np.sqrt(3)/2, 0],
        [0, 0, 1]
    ]

    basis = [
        (0, 0, 0),
        (1/3, 1/3, 0)
    ]
    return cell, basis

def buildBaseSupercell(cell, basis, nx=5, ny=5):

    positions = []
    symbols = []
    for i in range(nx):
        for j in range(ny):
            for b in basis:
                r = (i + b[0]) * np.array(cell[0]) + (j + b[1]) * np.array(cell[1])
                positions.append(r)
                symbols.append('C')

    atoms = Atoms(
    symbols=symbols,
    positions=positions
    )
    return atoms

def buildGrapheneSlitPore(d_slit, atoms, c_total):

    z_layer1 = (c_total - d_slit)/2

    layer1 = atoms.copy()
    layer1.translate([0, 0, z_layer1])

    layer2 = layer1.copy()
    layer2.translate([0, 0, d_slit])

    structure = layer1 + layer2

    return structure

def buildSimulationBox(cell, structure, c_total, nx=5, ny=5):
    box = structure.copy()
    # 修改整体 cell 高度
    supercell = [
        np.array(cell[0]) * nx,
        np.array(cell[1]) * ny,
        np.array([0, 0, c_total])
    ]
    box.set_cell(supercell)
    box.set_pbc([True, True, False])
    return box

def saveModel(file, box):
    box.write(file)

    # ===== 自动修复 CIF（关键）=====
    with open(file, 'r') as f:
        lines = f.readlines()

    # 在 data_ 后插入 space group
    new_lines = []
    inserted = False

    for line in lines:
        new_lines.append(line)
        if line.strip().startswith('data_') and not inserted:
            new_lines.append("_symmetry_space_group_name_H-M    'P1'\n")
            new_lines.append("_symmetry_Int_Tables_number       1\n")
            inserted = True

    with open(file, 'w') as f:
        f.writelines(new_lines)


def main(file_path=None):
    
    file = 'SG00635-QE.cif'
    c2c = 1.42
    nx, ny = 5, 5 #碳六环个数
    d_slit = 6.35       # slit pore 两 graphene 间距 A
    c_vaccum = 0  #真空层大小,对于GCMC，应该设置为0
    c_total = c_vaccum * 2 +  d_slit    # 整体 cell 高度，可调节

    cell, basis = buildLattice(c2c)
    atoms = buildBaseSupercell(cell, basis, nx, ny)
    structure = buildGrapheneSlitPore(d_slit, atoms, c_total)
    box = buildSimulationBox(cell, structure, c_total, nx=5, ny=5)
    saveModel(file, box)

if __name__ == "__main__":
    f = sys.argv[1] if len(sys.argv) > 1 else None
    main(f)