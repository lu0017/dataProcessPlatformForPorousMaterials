import numpy as np
from ase import Atoms
# 
# buildFGAtoms.py
# ├─ FrameBuilder
# ├─ constraint_types
# ├─ constraint_collector
# ├─ solver_dispatcher
# ├─ geometry_solvers
# └─ buildFGAtoms
# ==========================================
# helper functions
# ==========================================
import numpy as np
from dataclasses import dataclass
from dataclasses import dataclass, field
@dataclass
class BondConstraint:
    atom: int
    length: float
@dataclass
class AngleConstraint:
    atom: int
    angle: float
@dataclass
class DihedralConstraint:
    atom1: int
    atom2: int
    dihedral: float
@dataclass
class ConstraintSet:
    parent: int
    bond_constraints: list[BondConstraint]
    angle_constraints: list[AngleConstraint]
    dihedral_constraints: list[DihedralConstraint]
    used_anchor_angle: bool = False
def collect_constraints(fg, atom, parent, positions, mode="built"):
    """
    mode:
        "built" -> 仅使用已构建原子
        "all"   -> 使用全部拓扑约束
    """
    d = get_bond_length(fg, parent, atom)
    angle_info = find_all_angles(fg, atom, parent)
    dihedral_info = find_all_dihedrals(fg, atom, parent)
    bond_constraints = [ BondConstraint( atom=parent, length=d ) ]
    angle_constraints = []
    dihedral_constraints = []
    # -----------------------------------
    # 判断参考原子是否可用
    # -----------------------------------
    def available(ref):
        if ref == -1:
            return True
        if mode == "all":
            return True
        return ref < atom
    # -----------------------------------
    # 统计真实约束数量
    # -----------------------------------
    real_angle_count = sum(
        1
        for ref, _ in angle_info
        if ref != -1 and available(ref)
    )
    real_dihedral_count = sum(
        1
        for a, b, _
        in dihedral_info
        if (
            a != -1 and
            b != -1 and
            available(a) and
            available(b)
        )
    )
    use_anchor = ( real_angle_count + real_dihedral_count < 2 )
    used_anchor_angle = False
    # =====================================
    # angle constraints
    # =====================================
    for ref, angle in angle_info:
        if not available(ref):
            continue
        if ref == -1:
            if not use_anchor:
                continue
            used_anchor_angle = True
        angle_constraints.append( AngleConstraint( atom=ref, angle=angle ) )
    # =====================================
    # dihedral constraints
    # =====================================
    for a, b, phi in dihedral_info:
        if not available(a):
            continue
        if not available(b):
            continue
        if a == -1 or b == -1:
            if not use_anchor:
                continue
        dihedral_constraints.append(
            DihedralConstraint( atom1=a, atom2=b, dihedral=phi )
        )
    return ConstraintSet(
        parent=parent,
        bond_constraints=bond_constraints,
        angle_constraints=angle_constraints,
        dihedral_constraints=dihedral_constraints,
        used_anchor_angle=used_anchor_angle
    )
from collections import deque

def build_parent_map(fg):
    graph = {i: [] for i in range(len(fg.atoms))}
    for a, b in fg.bonds:
        graph[a].append(b)
        graph[b].append(a)
    root = fg.connection_atom
    parent_map = {root: -1}
    queue = deque([root])
    while queue:
        current = queue.popleft()
        for nbr in graph[current]:
            if nbr in parent_map:
                continue
            parent_map[nbr] = current
            queue.append(nbr)
    return parent_map
def get_parent(atom, parent_map):
    return parent_map.get(atom, -1)
def get_grandparent(atom, parent_map):
    p = parent_map.get(atom, -1)
    if p == -1:
        return -1
    return parent_map.get(p, -1)
import numpy as np
def normalize(v):
    v = np.asarray(v, dtype=float)
    n = np.linalg.norm(v)
    if n < 1e-12:
        raise ValueError("Zero vector")
    return v / n
def get_bond_length(fg, i, j):
    if (i, j) in fg.bond_lengths:
        return fg.bond_lengths[(i, j)]
    if (j, i) in fg.bond_lengths:
        return fg.bond_lengths[(j, i)]
    raise ValueError(
        f"Missing bond length for ({i},{j})"
    )
def find_parent(fg, atom):
    candidates = []
    for a, b in fg.bonds:
        if b == atom and a < atom:
            candidates.append(a)
        elif a == atom and b < atom:
            candidates.append(b)
    if len(candidates) == 0:
        return None
    if len(candidates) > 1:
        raise ValueError(
            f"Multiple parents for atom {atom}: {candidates}"
        )
    return candidates[0]
def find_angle(fg, atom, parent):
    for key, value in fg.bond_angles.items():
        if len(key) != 3:
            continue
        a, b, c = key
        # parent - atom 关系可能在 b-c 或 a-b
        if b == parent and c == atom:
            return a, value
        if b == atom and a == parent:
            return c, value
    return None, 109.5
def find_dihedral(fg, atom, parent):
    for key, value in fg.dihedral_angles.items():
        if len(key) != 4:
            continue
        a, b, c, d = key
        # 标准方向
        if c == parent and d == atom:
            return a, value
        # 反向也允许
        if b == parent and a == atom:
            return d, value
    return None, 180.0
def find_specific_angles(fg, atom, parent, built_atoms):
    angle_atom = None
    angle = 109.0
    for key, value in fg.bond_angles.items():
        if len(key) != 3:
            continue
        a, b, c = key
        if b != parent:
            continue
        if c == atom:
            if a in built_atoms:
                angle_atom = a
                angle = value
        elif a == atom:
            if c in built_atoms:
                angle_atom = c
                angle = value
    return angle_atom, angle
def find_all_angles(fg, atom, parent):
    result = []
    for key,value in fg.bond_angles.items():
        if len(key) != 3:
            continue
        a,b,c = key
        if b != parent:
            continue
        if c == atom:
            result.append((a,value))
        elif a == atom:
            result.append((c,value))
    return result
def find_all_dihedrals(fg, atom, parent):
    result = []
    for key, value in fg.dihedral_angles.items():
        if len(key) != 4:
            continue
        a, b, c, d = key
        # 正向
        if c == parent and d == atom:
            result.append((a, b, value))
        # 反向（避免漏掉镜像定义）
        elif b == parent and a == atom:
            result.append((d, c, value))
    return result
def calc_angle(A, B, C):
    BA = normalize(A - B)
    BC = normalize(C - B)
    cos_theta = np.dot(BA, BC)
    cos_theta = np.clip( cos_theta, -1.0, 1.0 )
    return np.degrees( np.arccos(cos_theta) )
def get_refs(constraints, positions, parent_map):
    """
    A : dihedral reference
    B : angle reference
    C : parent
    """
    parent = constraints.parent
    C = positions[parent]
    if len(constraints.angle_constraints) == 0:
        raise RuntimeError(
            "Z-matrix requires at least one angle."
        )
    B_atom = constraints.angle_constraints[0].atom
    B = positions[B_atom]
    # 有真实二面角
    if len(constraints.dihedral_constraints) > 0:
        A_atom = constraints.dihedral_constraints[0].atom1
    # 无二面角，寻找B的父节点
    else:
        A_atom = parent_map.get(B_atom, B_atom)
        if A_atom == B_atom:
            # fallback
            return B, B, C
    A = positions[A_atom]
    return A, B, C
def build_local_frame_zmatrix(A, B, C):
    ex = normalize(C - B)
    n1 = np.cross(B - A, ex)
    if np.linalg.norm(n1) < 1e-6:
        test = np.array([0.0, 0.0, 1.0])
        if abs(np.dot(test, ex)) > 0.95:
            test = np.array([1.0, 0.0, 0.0])
        n1 = np.cross(test, ex)
    n1 = normalize(n1)
    n2 = normalize(np.cross(ex, n1))
    return ex, n1, n2
def calc_dihedral(A, B, C, D):
    b0 = A - B
    b1 = C - B
    b2 = D - C
    b1 = normalize(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    x = np.dot(v, w)
    y = np.dot( np.cross(b1, v), w )
    return np.degrees(
        np.arctan2(y, x)
    )

def solve_zmatrix_local(d, theta, phi):
    x = -np.cos(theta)
    y = np.sin(theta) * np.cos(phi)
    z = np.sin(theta) * np.sin(phi)
    return d * np.array([x, y, z])

def solve_two_angle_direction( constraints, positions):
    parent = constraints.parent
    parent_pos = positions[parent]
    atom1 = ( constraints.angle_constraints[0].atom )
    atom2 = ( constraints.angle_constraints[1].atom )
    P1 = positions[atom1]
    P2 = positions[atom2]
    theta1 = ( constraints.angle_constraints[0].angle )
    theta2 = ( constraints.angle_constraints[1].angle )
    d = ( constraints .bond_constraints[0].length )
    u1 = P1 - parent_pos
    u2 = P2 - parent_pos
    return solve_from_two_angles_local( u1, u2, d, theta1, theta2 )

def solve_from_two_angles_local(u1, u2, d, theta1, theta2, positive_z=True):
    u1 = normalize(u1)
    u2 = normalize(u2)
    cos12 = np.dot(u1, u2)
    c1 = np.cos(np.deg2rad(theta1))
    c2 = np.cos(np.deg2rad(theta2))
    denom = 1.0 - cos12**2
    if abs(denom) < 1e-8:
        raise ValueError("Two reference directions are nearly collinear")
    a = (c1 - c2 * cos12) / denom
    b = (c2 - c1 * cos12) / denom
    n = np.cross(u1, u2)
    if np.linalg.norm(n) < 1e-8:
        raise ValueError("Cannot define plane for double-angle construction")
    n = normalize(n)
    c_sq = 1.0 - a*a - b*b - 2*a*b*cos12
    c_sq = max(0.0, c_sq)
    c = np.sqrt(c_sq)
    if not positive_z:
        c = -c
    v = a*u1 + b*u2 + c*n
    return normalize(v)

from scipy.optimize import least_squares
def solve_constrained_position( constraints, positions, parent_map):
    parent = constraints.parent
    center = positions[parent]
    bond_length = ( constraints.bond_constraints[0].length )
    def residual(x):
        direction = normalize(x)
        pos = center + bond_length * direction
        errors = []
        # ------------------
        # angle constraints
        # ------------------
        for c in constraints.angle_constraints:
            ref = positions[c.atom]
            ref_v = normalize( ref - center )
            angle = np.degrees(
                np.arccos(
                    np.clip( np.dot( ref_v, direction ), -1.0, 1.0 )
                )
            )
            errors.append( angle - c.angle )
        # ------------------
        # dihedral constraints
        # ------------------
        for c in constraints.dihedral_constraints:
            A = positions[c.atom1]
            B = positions[c.atom2]
            phi = calc_dihedral( A, B, center, pos )
            diff = phi - c.dihedral
            while diff > 180:
                diff -= 360
            while diff < -180:
                diff += 360
            errors.append( 0.2 * diff )
        return errors
    # ------------------
    # initial guess
    # ------------------
    if len(constraints.angle_constraints) > 0:
        guess = np.zeros(3)
        for c in constraints.angle_constraints:
            ref = positions[c.atom]
            guess += ( center - ref )
        guess = normalize(guess)
    else:
        guess = np.array( [1.0, 0.0, 0.0] )
    result = least_squares( residual, guess, method="trf" )
    return normalize(result.x)

def solve_zmatrix( constraints, positions, parent_map):
    A, B, C = get_refs( constraints, positions, parent_map )
    ex, n1, n2 = build_local_frame_zmatrix( A, B, C )
    theta = np.deg2rad( constraints.angle_constraints[0].angle )
    if len(constraints.dihedral_constraints) > 0:
        phi = np.deg2rad( constraints.dihedral_constraints[0].dihedral )
    else:
        phi = 0.0
    d = constraints.bond_constraints[0].length
    p_local = solve_zmatrix_local( d, theta, phi )
    p_global = to_global( C, ex, n2, n1, p_local )
    return normalize( p_global - C )

def solve_two_angle(parent_pos, constraints):
    P1 = constraints.references[0]
    P2 = constraints.references[1]
    u1 = P1 - parent_pos
    u2 = P2 - parent_pos
    v = solve_from_two_angles_local(
        u1, u2,
        constraints.bond_length,
        constraints.target_angles[0],
        constraints.target_angles[1]
    )
    return parent_pos + constraints.bond_length * v

def solve_optimizer(parent_pos, constraints):
    refs = constraints.references - parent_pos  # 关键改动
    return solve_constrained_position(
        center=parent_pos,
        bond_length=constraints.bond_length,
        references=refs,
        target_angles=constraints.target_angles,
        dihedral_refs=constraints.dihedral_refs,
        target_dihedrals=constraints.target_dihedrals
    )

def solve_default_local( constraints, positions, parent_map):
    parent = constraints.parent
    # 尝试沿parent→grandparent反方向
    grandparent = parent_map.get(parent, -1)
    if grandparent != -1:
        v = ( positions[parent] - positions[grandparent] )
        if np.linalg.norm(v) > 1e-8:
            return normalize(v)
    return np.array([1.0,0.0,0.0])

def to_global(C, ex, n2, n1, p_local):
    return (
        C
        + p_local[0] * ex
        + p_local[1] * n2
        + p_local[2] * n1
    )

import numpy as np
def select_solver(constraints):
    """
    0个角度      -> default
    1个角度      -> zmatrix
    2个角度      -> two_angle
    >=3个角度    -> optimizer
    任意二面角   -> optimizer
    """
    n_angle = len(constraints.angle_constraints)
    n_dihedral = len(constraints.dihedral_constraints)
    if n_dihedral > 0:
        return "optimizer"
    if n_angle == 0:
        return "default"
    if n_angle == 1:
        return "zmatrix"
    if n_angle == 2 and not constraints.used_anchor_angle:
        return "two_angle"
    return "optimizer"

class LocalEngine:
    def __init__(self, fg):
        self.fg = fg
    def run( self, atom, parent, constraints, positions, parent_map):
        solver = select_solver(constraints)
        if solver == "default":
            direction = solve_default_local( constraints )
        elif solver == "zmatrix":
            direction = solve_zmatrix( constraints, positions, parent_map )
        elif solver == "two_angle": 
            direction = solve_two_angle_direction( constraints, positions )
        else:
            direction = solve_constrained_position( constraints, positions, parent_map )
        d = constraints.bond_constraints[0].length
        return positions[parent] + d * direction

def build_reference_frameOld(fg):
    n = len(fg.atoms)
    positions = np.zeros((n,3))
    # atom0
    positions[0] = [0.0,0.0,0.0]
    if n == 1:
        return positions
    # atom1
    parent = find_parent(fg,1)
    if parent is None:
        parent = 0
    d = get_bond_length( fg, parent, 1 )
    positions[1] = [d,0.0,0.0]
    if n == 2:
        return positions
    # atom2
    parent = find_parent(fg,2)
    if parent is None:
        parent = 0
    angle_atom, angle = find_angle(
        fg,
        2,
        parent
    )
    d = get_bond_length(
        fg,
        parent,
        2
    )
    theta = np.deg2rad(angle)
    positions[2] = np.array([ positions[parent][0] + d*np.cos(theta),
        positions[parent][1] + d*np.sin(theta), 0.0
    ])
    return positions
import numpy as np


def build_reference_frame(fg):
    """
    构建FGS内部局部参考系
    atom0 : 原点
    atom1 : X轴正方向
    atom2 : XY平面且Y>0
    返回
    -------
    positions : (n,3)
    """
    n = len(fg.atoms)
    positions = np.zeros((n, 3))
    # --------------------------
    # atom0
    # --------------------------
    positions[0] = [0.0, 0.0, 0.0]
    if n == 1:
        return positions
    # --------------------------
    # atom1
    # --------------------------
    parent = find_parent(fg, 1)
    if parent is None:
        parent = 0
    d = get_bond_length(fg, parent, 1)
    positions[1] = [d, 0.0, 0.0]
    if n == 2:
        return positions
    # --------------------------
    # atom2
    # --------------------------
    built_atoms = {0,1}
    parent = find_parent(fg, 2)
    if parent is None:
        parent = 1
    angle_atom, angle = find_specific_angles(fg, atom=2, parent=parent, built_atoms=built_atoms)
    if angle_atom is None:
        raise ValueError(
            "Cannot determine angle reference for atom 2."
        )
    d = get_bond_length( fg, parent, 2 )
    # P -> A
    u = positions[angle_atom] - positions[parent]
    norm = np.linalg.norm(u)
    if norm < 1e-10:
        raise ValueError( "Reference vector length is zero." )
    u /= norm
    # 法向方向
    nvec = np.array([ -u[1], u[0], 0.0 ])
    theta = np.deg2rad(angle)
    # 第一组解
    direction = ( np.cos(theta) * u + np.sin(theta) * nvec )
    candidate = positions[parent] + d * direction
    # 固定Y正方向
    if candidate[1] < positions[parent][1]:
        direction = ( np.cos(theta) * u - np.sin(theta) * nvec )
    positions[2] = ( positions[parent] + d * direction )
    return positions
def buildFGAtoms(fg):
    n = len(fg.atoms)
    if n <= 1 :
        return Atoms( symbols=fg.atoms, positions=[[0,0,0]] )
    if n == 2:
        d = get_bond_length(fg,0,1)
        positions = np.array([ [0,0,0], [d,0,0] ])
        return Atoms( symbols=fg.atoms, positions=positions )
    parent_map = build_parent_map(fg)
    positions = build_reference_frame(fg)
    engine = LocalEngine(fg)
    for atom in range(3,n):
        parent = find_parent(fg, atom)
        if parent is None:
            raise ValueError(f"No parent for atom {atom}")
        constraints = collect_constraints( fg, atom, parent, positions )
        positions[atom] = engine.run( atom, parent, constraints, positions, parent_map )
    return Atoms( symbols=fg.atoms, positions=positions )

def buildFGAtoms00(fg):
    if fg.topology in ["substitution", "vacancy"]:
        symbol = fg.dopant if fg.dopant else "C"
        return Atoms(
            symbols=[symbol],
            positions=[[0.0, 0.0, 0.0]]
        )
    parent_map = build_parent_map(fg)
    positions = build_reference_frame(fg)
    engine = LocalEngine(fg)
    n = len(fg.atoms)
    for atom in range(3, n):
        parent = find_parent(fg, atom)
        if parent is None:
            raise ValueError(f"No parent for atom {atom}")
        constraints = collect_constraints(fg, atom, parent, positions)
        positions[atom] = engine.run(
            atom=atom,
            parent=parent,
            constraints=constraints,
            positions=positions,
            parent_map=parent_map
        )
    return Atoms(
        symbols=fg.atoms,
        positions=positions
    )