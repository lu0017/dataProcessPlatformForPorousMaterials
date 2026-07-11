from ase import Atoms
import numpy as np
import sys

def align_vectors(v1, v2):
    v1 = v1 / np.linalg.norm(v1)
    v2 = v2 / np.linalg.norm(v2)
    cross = np.cross(v1, v2)
    dot = np.dot(v1, v2)
    if dot > 0.999999:
        return np.eye(3)
    if dot < -0.999999:
        axis = np.array([1.0,0.0,0.0])
        if abs(v1[0]) > 0.9:
            axis = np.array([0.0,1.0,0.0])
        axis = np.cross(v1, axis)
        axis /= np.linalg.norm(axis)
        K = np.array([
            [0,-axis[2],axis[1]],
            [axis[2],0,-axis[0]],
            [-axis[1],axis[0],0]
        ])
        return np.eye(3) + 2*K@K
    K = np.array([
        [0,-cross[2],cross[1]],
        [cross[2],0,-cross[0]],
        [-cross[1],cross[0],0]
    ])
    R = np.eye(3)
    R += K
    R += K @ K * (1.0 / (1.0 + dot))
    return R

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

def buildGrapheneMonolayer(atoms, vacuum=20.0):
    structure = atoms.copy()
    structure.translate([0, 0, vacuum/2]) # 把graphene放在盒子中间
    return structure

def findCenterAtom(atoms):
    positions = atoms.get_positions()
    x_center = positions[:,0].mean()
    y_center = positions[:,1].mean()
    distance = (
        (positions[:,0]-x_center)**2 +
        (positions[:,1]-y_center)**2
    )
    center_id = np.argmin(distance)
    return center_id

def buildGrapheneSlitPore(d_slit, atoms, c_total):
    z_layer1 = (c_total - d_slit)/2
    layer1 = atoms.copy()
    layer1.translate([0, 0, z_layer1])
    layer2 = layer1.copy()
    layer2.translate([0, 0, d_slit])
    structure = layer1 + layer2
    return structure

import numpy as np
from ase import Atom


def terminate_edges(substrate, cutoff=1.75, bond_CH=1.09):
    """
    删除角点C并对边缘C进行H封端。
    规则：
    -------
    C邻居数=1:
        删除（角点）
    C邻居数=2:
        若第三个邻居不存在（没有O/N/H）
        则加H。
    Parameters
    ----------
    substrate : ASE Atoms
    cutoff : float
    bond_CH : float
        C-H键长（Å）
    Returns
    -------
    substrate
    """
    # =====================
    # 1. 删除角点
    # =====================
    corner_atoms = []
    for i, atom in enumerate(substrate):
        if atom.symbol != "C":
            continue
        neighbors = get_neighbors( substrate, i, cutoff )
        nC = sum(
            substrate[j].symbol == "C"
            for j in neighbors
        )
        if nC == 1:
            corner_atoms.append(i)
    # 必须倒序删除
    for idx in sorted(corner_atoms, reverse=True):
        del substrate[idx]
    # =====================
    # 2. 边缘H封端
    # =====================
    H_positions = []
    for i, atom in enumerate(substrate):
        if atom.symbol != "C":
            continue
        neighbors = get_neighbors( substrate, i, cutoff )
        carbon_neighbors = [
            j for j in neighbors
            if substrate[j].symbol == "C"
        ]
        # 必须只有两个碳邻居
        if len(carbon_neighbors) != 2:
            continue
        # 已经连接O/N/H则跳过
        if len(neighbors) != 2:
            continue
        p = substrate.positions[i]
        p1 = substrate.positions[carbon_neighbors[0]]
        p2 = substrate.positions[carbon_neighbors[1]]
        v1 = p - p1
        v2 = p - p2
        v1 /= np.linalg.norm(v1)
        v2 /= np.linalg.norm(v2)
        direction = v1 + v2
        norm = np.linalg.norm(direction)
        if norm < 1e-8:
            continue
        direction /= norm
        H_pos = p + bond_CH * direction
        H_positions.append(H_pos)
    # =====================
    # 3. 添加H
    # =====================
    for pos in H_positions:
        substrate.append( Atom("H", pos) )
    return substrate

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

import numpy as np


def get_neighbors(atoms, atom_id, cutoff=1.75):
    """
    返回 atom_id 的第一近邻。
    Parameters
    ----------
    atoms : ASE Atoms
    atom_id : int
    cutoff : float
        邻居截断距离，graphene推荐1.75 Å
    Returns
    -------
    list[int]
    """
    pos0 = atoms.positions[atom_id]
    neighbors = []
    for i, pos in enumerate(atoms.positions):
        if i == atom_id:
            continue
        d = np.linalg.norm(pos - pos0)
        if d < cutoff:
            neighbors.append((i, d))
    neighbors.sort(key=lambda x: x[1])
    return [i for i, d in neighbors]

from ase import Atoms, Atom
def dope_grafting( substrate, anchor_id, fg, fg_pos):
    start = 1 if fg.anchor_included else 0
    for elem, pos in zip( fg.atoms[start:], fg_pos[start:]):
        substrate.append( Atom(elem, pos) )
    return substrate

def dope_substitution(substrate, anchor_id, fg, fg_pos):
    substrate[anchor_id].symbol = fg.atoms[0]
    substrate[anchor_id].position = fg_pos[0]
    for elem, pos in zip(fg.atoms[1:], fg_pos[1:]):
        substrate.append( Atom(elem, position=pos) )
    return substrate

import numpy as np
from ase import Atom


def dope_vacancy(substrate, anchor_id, fg, fg_pos, anchor_context):
    anchor_pos = anchor_context["anchor_pos"]
    neighbor_pos = anchor_context["neighbor_pos"]
    # 默认选择前两个边缘碳
    p1 = neighbor_pos[0]
    p2 = neighbor_pos[1]
    midpoint = 0.5 * (p1 + p2)
    # vacancy中心指向边缘
    direction = anchor_pos - midpoint
    norm = np.linalg.norm(direction)
    if norm < 1e-8:
        direction = anchor_context["normal"]
    else:
        direction /= norm
    # pyridinic N键长
    CN = 1.35
    half_CC = np.linalg.norm(p1 - p2) / 2
    shift = np.sqrt(
        max(CN**2 - half_CC**2, 0.0)
    )
    N_pos = midpoint + shift * direction
    # 删除anchor
    del substrate[anchor_id]
    # 插入N
    substrate.append( Atom(fg.atoms[0], N_pos) )
    return substrate

def inject_anchor_id(fg, anchor_id):
    #把FunctionalGroup中anchor的-1编号替换为正常编号
    def rep(x):
        return anchor_id if x == -1 else x
    import copy
    fg = copy.deepcopy(fg)
    fg.bond_lengths = {
        tuple(rep(i) for i in k): v
        for k, v in fg.bond_lengths.items()
    }
    fg.bond_angles = {
        tuple(rep(i) for i in k): v
        for k, v in fg.bond_angles.items()
    }
    fg.dihedral_angles = {
        tuple(rep(i) for i in k): v
        for k, v in fg.dihedral_angles.items()
    }
    return fg

import numpy as np

def extract_anchor_context(substrate, anchor_id):
    # 提取待doping位置原子以及周围原子信息
    anchor_pos = substrate.positions[anchor_id]
    neighbors = get_neighbors(substrate, anchor_id)
    neighbor_pos = substrate.positions[neighbors]
    # ---------- 局部法向 ----------
    normal = np.array([0.0, 0.0, 1.0])
    if len(neighbor_pos) >= 2:
        v1 = neighbor_pos[0] - anchor_pos
        v2 = neighbor_pos[1] - anchor_pos
        normal = np.cross(v1, v2)
        norm = np.linalg.norm(normal)
        if norm > 1e-8:
            normal /= norm
        else:
            normal = np.array([0.0, 0.0, 1.0])
    # 保证法向朝 +z
    if normal[2] < 0:
        normal = -normal
    return {
        "anchor_id": anchor_id,
        "anchor_pos": anchor_pos,
        "neighbors": neighbors,
        "neighbor_pos": neighbor_pos,
        "normal": normal
    }

import buildFGS
def build_local_frame(anchor_context):
    anchor_pos = anchor_context["anchor_pos"]
    neighbors = anchor_context["neighbor_pos"]
    v1 = neighbors[0] - anchor_pos
    v2 = neighbors[1] - anchor_pos
    x = buildFGS.normalize(v1)
    z = buildFGS.normalize(np.cross(v1, v2))
    y = np.cross(z, x)
    return np.vstack([x, y, z])

def buildFGAtomsOnSubstrate(substrate, anchor_id, fg, anchor_context):
    # 1. 构建局部FG
    fg_local = buildFGS.buildFGAtoms(fg)
    fg_local.write(f"FGS_00.xyz")
    # --------------------------------
    # 默认方向：0→1
    # --------------------------------
    anchor_pos = anchor_context["anchor_pos"]
    normal = anchor_context["normal"]
    if len(fg_local) >= 2:
        fg_axis = ( fg_local.positions[1] - fg_local.positions[0] )
    else:
        fg_axis = np.array([1.0,0.0,0.0])
    fg_axis /= np.linalg.norm(fg_axis)
    R = align_vectors(fg_axis, normal)
    fg_global_positions = []
    for p in fg_local.positions:
        p_global = R @ p + anchor_pos
        fg_global_positions.append(p_global)
    # 4. topology分支处理（你要求放这里）
    if fg.topology == "substitution":
        return dope_substitution( substrate, anchor_id, fg, fg_global_positions )
    elif fg.topology == "vacancy":
        return dope_vacancy( substrate, anchor_id, fg, fg_global_positions, anchor_context )
    elif fg.topology == "grafting":
        return dope_grafting( substrate, anchor_id, fg, fg_global_positions )
    else:
        raise ValueError("Unknown topology")
    
import functionalGroups as FGS
import copy

import copy

def expand_anchor_for_doping(fg):
    """
    grafting专用预处理。
    原：
        -1
         |
         0
        / \
       1   2
    变：
         0(anchor)
            |
         1(connection atom)
        / \
       2   3
    FGS数据库本身完全不修改。
    """
    fg = copy.deepcopy(fg)
    # 已经展开过
    if fg.anchor_included:
        return fg
    # ------------------------------------------------
    # atoms
    # ------------------------------------------------
    fg.atoms = ["X"] + fg.atoms
    # ------------------------------------------------
    # bonds
    # ------------------------------------------------
    new_bonds = []
    for a, b in fg.bonds:
        new_bonds.append((a + 1, b + 1))
    # anchor → connection atom
    new_bonds.insert(0, (0, fg.connection_atom + 1))
    fg.bonds = new_bonds
    # ------------------------------------------------
    # bond lengths
    # ------------------------------------------------
    new_lengths = {}
    for (a, b), value in fg.bond_lengths.items():
        if a == -1:
            new_lengths[(0, b + 1)] = value
        else:
            new_lengths[(a + 1, b + 1)] = value
    fg.bond_lengths = new_lengths
    # ------------------------------------------------
    # bond angles
    # ------------------------------------------------
    new_angles = {}
    for (a, b, c), value in fg.bond_angles.items():
        aa = 0 if a == -1 else a + 1
        bb = b + 1
        cc = c + 1
        new_angles[(aa, bb, cc)] = value
    fg.bond_angles = new_angles
    # ------------------------------------------------
    # dihedrals
    # ------------------------------------------------
    new_dihedrals = {}
    for key, value in fg.dihedral_angles.items():
        a, b, c, d = key
        aa = 0 if a == -1 else a + 1
        bb = b + 1
        cc = c + 1
        dd = d + 1
        new_dihedrals[(aa, bb, cc, dd)] = value
    fg.dihedral_angles = new_dihedrals
    # ------------------------------------------------
    # bond order
    # ------------------------------------------------
    new_orders = {}
    for (a, b), value in fg.bond_orders.items():
        new_orders[(a + 1, b + 1)] = value
    fg.bond_orders = new_orders
    fg.anchor_included = True
    return fg

def dopeFunction(substrate, anchor_id, fg_name):
    if fg_name == None:
        return substrate
    fg = copy.deepcopy(FGS.FG_DATABASE[fg_name])
    if fg.topology == "grafting":
        fg = expand_anchor_for_doping(fg)
    # fg = inject_anchor_id(fg, anchor_id)
    anchor_context = extract_anchor_context(substrate, anchor_id)
    return buildFGAtomsOnSubstrate( substrate, anchor_id, fg, anchor_context )
    
def main(file_path=None):
    file = 'SG00635-QE.xyz'
    c2c = 1.42
    nx, ny = 5, 5 #碳六环个数
    d_slit = 6.35       # slit pore 两 graphene 间距 A
    c_vaccum = 0  #真空层大小,对于GCMC，应该设置为0
    c_total = c_vaccum * 2 +  d_slit    # 整体 cell 高度，可调节
    cell, basis = buildLattice(c2c)
    atoms = buildBaseSupercell(cell, basis, nx, ny)
    # structure = buildGrapheneSlitPore(d_slit, atoms, c_total) #GCMC 1layer用
    #QE优化前的model
    center_id = findCenterAtom(atoms)
    atoms = terminate_edges(atoms)
    FG = "COOH" #"NH2"
    atoms = dopeFunction(atoms, center_id, FG)
    atoms = terminate_edges(atoms)
    structure = buildGrapheneMonolayer(atoms, vacuum=20.0)
    box = buildSimulationBox(cell, structure, c_total, nx=5, ny=5)
    saveModel(file, box)

if __name__ == "__main__":
    f = sys.argv[1] if len(sys.argv) > 1 else None
    main(f)