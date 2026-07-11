from dataclasses import dataclass, field
@dataclass
class FunctionalGroup:
    # ========= 基本信息 =========
    name: str
    family: str
    topology: str
        # substitution
        # vacancy
        # grafting
    # ========= 结构拓扑 =========
    connection_atom: int = 0
    atoms: list[str] = field(default_factory=list)
    bonds: list[tuple] = field(default_factory=list) #按照atoms的原子顺序编码，从0开始
    # ========= 几何参数 =========
    bond_lengths: dict = field(default_factory=dict) #按照atoms的原子顺序编码，从0开始
    bond_orders: dict = field(default_factory=dict) #键类型，1 = single，2 = double 3 = triple 1.5 = aromatic
    bond_angles: dict = field(default_factory=dict) #键角，(1,0,2):107.0,按照atoms的原子顺序编码，从0开始，-1表示碳环基体
    dihedral_angles: dict = field(default_factory=dict) #二面角，(1,0,2,3):180；按照atoms的原子顺序编码，从0开始，-1表示碳环基体
    # ========= substitution =========
    dopant: str = None
    # ========= vacancy =========
    vacancy_size: int = None
    # ========= grafting =========
    anchor_included: bool = False
    # ========= 扩展 =========
    metadata: dict = field(default_factory=dict)
GRAPHITIC_N = FunctionalGroup(
    name="graphitic-N",
    family="N",
    topology="substitution",
    atoms=[
        "N"
    ],
    dopant="N"
)
PYRIDINIC_N = FunctionalGroup(
    name="pyridinic-N",
    family="N",
    topology="vacancy",
    dopant="N",
    atoms=[
        "N"
    ],
    vacancy_size=1
)
NH2 = FunctionalGroup(
    name="NH2",
    family="N",
    topology="grafting",
    connection_atom=0,
    atoms=[
        "N",  # 0
        "H",  # 1
        "H"   # 2
    ],
    bonds=[
        (0,1),
        (0,2)
    ],
    bond_lengths={
        (-1,0):1.47,
        (0,1):1.01,
        (0,2):1.01
    },
    bond_angles={
        (1,0,2):107.0,
        (-1,0,1):107.0,
        (-1,0,2):107.0
    },
    metadata={
        "reference":"https://doi.org/10.1021/jp065336v"
    }
)
OH = FunctionalGroup(
    name="OH",
    family="O",
    topology="grafting",
    connection_atom=0,
    atoms=[
        "O",  # 0
        "H"   # 1
    ],
    bonds=[
        (0,1)
    ],
    bond_lengths={
        (-1,0):1.36,
        (0,1):0.97
    },
    bond_angles={
        (-1,0,1):108.0
    },
    metadata={
        "reference":"doi:10.1016/j.colsurfa.2023.133113"
    }
)
C_O = FunctionalGroup(
    name="C=O",
    family="O",
    topology="grafting",
    connection_atom=0,
    atoms=[
        "O"  # 0
    ],
    bonds=[],
    bond_orders={
        (-1,0):2
    },
    bond_lengths={
        (-1,0):1.23
    },
    metadata={
        "reference":"doi:10.1016/j.colsurfa.2023.133113"
    }
)
COOH = FunctionalGroup(
    name="COOH",
    family="O",
    topology="grafting",
    connection_atom=0,
    atoms=[
        "C",  # 0
        "O",  # 1
        "O",  # 2
        "H"   # 3
    ],
    bonds=[
        (0,1),
        (0,2),
        (2,3)
    ],
        bond_orders={
        (0,1):2
    },
    bond_lengths={
        (-1,0):1.50,
        (0,1):1.23,
        (0,2):1.36,
        (2,3):0.97
    },
    bond_angles={
        (1,0,2):123.0,
        (0,2,3):108.0,
        (-1,0,1):111.0,
        (-1,0,2):120.0
    },
    metadata={
        "reference":"doi:10.1016/j.colsurfa.2023.133113"
    }
)
SO = FunctionalGroup(
    name="SO",
    family="S",
    topology="grafting",
    connection_atom=0,
    atoms=[
        "S",  # 0
        "O"   # 1
    ],
    bonds=[
        (0,1)
    ],
    bond_orders={
        (0,1):2
    },
    bond_lengths={
        (-1,0):1.78,
        (0,1):1.50
    },
    bond_angles={
        (-1,0,1):106.5
    },
    metadata={
        "reference":"doi:10.1107/S2056989017012464"
    }
)
SO2 = FunctionalGroup(
    name="SO2",
    family="S",
    topology="grafting",
    connection_atom=0,
    atoms=[
        "S",  # 0
        "O",  # 1
        "O"   # 2
    ],
    bonds=[
        (0,1),
        (0,2)
    ],
    bond_orders={
        (0,1):2,
        (0,2):2
    },
    bond_lengths={
        (-1,0):1.78,
        (0,1):1.45,
        (0,2):1.45
    },
    bond_angles={
        (1,0,2):118.0,
        (-1,0,1):106.0,
        (-1,0,2):106.0
    },
        metadata={
        "reference":"doi:https://doi.org/10.1007/s10870-009-9643-8"
    }
)
SO3H = FunctionalGroup(
    name="SO3H",
    family="S",
    topology="grafting",
    connection_atom=0,
    atoms=[
        "S",  # 0
        "O",  # 1
        "O",  # 2
        "O",  # 3
        "H"   # 4
    ],
    bonds=[
        (0,1),
        (0,2),
        (0,3),
        (3,4)
    ],
    bond_orders={
    (0,1):2,
    (0,2):2
    },
    bond_lengths={
        (-1,0):1.78,
        (0,1):1.45,
        (0,2):1.45,
        (0,3):1.60,
        (3,4):0.97
    },
    bond_angles={
        (1,0,2):120.0,
        (1,0,3):113.0,
        (2,0,3):113.0,
        (-1,0,2):106.0,
        (-1,0,3):106.0,
        (-1,0,1):106.0,
        (0,3,4):104.5
    },
    dihedral_angles={
        (1,0,3,4):180
    },
    metadata={
        "reference1":"doi: 10.1107/S2056989015015650",
        "reference1":"doi: https://doi.org/10.1107/S2056989019013367"
    }
)
GRAPHITIC_P = FunctionalGroup(
    name="graphitic-P",
    family="P",
    topology="substitution",
    atoms=[
        "P"
    ],
    dopant="P",
    metadata={
        "reference1":"doi: 10.1016/j.carbon.2019.10.018"
    }
)
C3_PO = FunctionalGroup(
    name="C3-P=O",
    family="P",
    topology="substitution",
    dopant="P",
    atoms=[
        "P",  # 0
        "O"   # 1
    ],
    bonds=[
        (0,1)
    ],
    bond_orders={
    (0,1):2
    },
    bond_lengths={
        (0,1):1.45
    },
    metadata={
        "oxide":True,
        "reference":"doi: 10.3390/molecules25122778"
    }
)
C_PO3H2 = FunctionalGroup(
    name="C-P(=O)(OH)2",
    family="P",
    topology="grafting",
    connection_atom=0,
    atoms=[
        "P",  # 0
        "O",  # 1 P=O
        "O",  # 2 P-OH
        "H",  # 3
        "O",  # 4 P-OH
        "H"   # 5
    ],
    bonds=[
        (0,1),
        (0,2),
        (2,3),
        (0,4),
        (4,5)
    ],
    bond_orders={
        (0,1):2
    },
    bond_lengths={
        (-1,0):1.78,  # C-P
        (0,1):1.50,         # P=O
        (0,2):1.54,         # P-OH
        (2,3):0.97,         # O-H
        (0,4):1.54,         # P-OH
        (4,5):0.97
    },
    bond_angles={
        (-1,0,1):110.0,
        (-1,0,2):109.5,
        (-1,0,4):109.5,
        (1,0,2):109.5,
        (1,0,4):109.5,
        (2,0,4):109.5,
        (0,2,3):108.0,
        (0,4,5):108.0
    },
    dihedral_angles={
        (-1,0,2,3):180.0,
        (-1,0,4,5):-60.0
    },
    metadata={
        "reference1":"doi: https://doi.org/10.3390/molecules25122778",
        "reference2":"doi: https://doi.org/10.1021/acsomega.8b03192",
    }
)
CO_PO3H2 = FunctionalGroup(
    name="C-O-P(=O)(OH)2",
    family="P",
    topology="grafting",
    connection_atom=0,
    atoms=[
        "O",  # 0 连接Anchor
        "P",  # 1
        "O",  # 2 P=O
        "O",  # 3 P-OH
        "H",  # 4
        "O",  # 5 P-OH
        "H"   # 6
    ],
    bonds=[
        (0,1),
        (1,2),
        (1,3),
        (3,4),
        (1,5),
        (5,6)
    ],
    bond_orders={
        (1,2):2
    },
    bond_lengths={
        (-1,0):1.43,  # C-O
        (0,1):1.62,         # O-P
        (1,2):1.50,         # P=O
        (1,3):1.54,         # P-OH
        (3,4):0.97,
        (1,5):1.54,
        (5,6):0.97
    },
    bond_angles={
        (-1,0,1):109.5,
        (0,1,2):109.5,
        (0,1,3):109.5,
        (0,1,5):109.5,
        (2,1,3):109.5,
        (2,1,5):109.5,
        (3,1,5):105.0,
        (1,3,4):108.0,
        (1,5,6):108.0
    },
    dihedral_angles={
        (-1,0,1,2):180.0,
        (-1,0,1,3):60.0,
        (-1,0,1,5):-60.0,
        # (0,1,3,4):-60.0
    },
    metadata={
        "oxide":True,
        "reference":"doi: https://doi.org/10.1021/acsomega.8b03192",
    }
)

@dataclass
class FunctionalGroupStructure:
    from ase import Atoms
    fg: FunctionalGroup
    atoms: Atoms = None
    optimized_atoms: Atoms = None
    charge: dict = None
    metadata: dict = None

# import old_buildFGS as buildFGsF
import buildFGS as buildFGsF
def createFGStructure(name):
    fg = FG_DATABASE[name]
    structure = FunctionalGroupStructure(fg)
    structure.atoms = buildFGsF.buildFGAtoms(fg)
    return structure

FG_DATABASE = {
    "GRAPHITIC_N":GRAPHITIC_N,
    "PYRIDINIC_N":PYRIDINIC_N,
    "NH2":NH2,
    "OH":OH,
    "C=O":C_O,
    "COOH":COOH,
    "SO":SO,
    "S=O2":SO2,
    "SO3H":SO3H,
    "graphitic-P":GRAPHITIC_P,
    "C3-P=O":C3_PO,
    "C_PO3H2":C_PO3H2,
    "C_O_PO3H2":CO_PO3H2
}
# fgs_name = "NH2"
# fgs = createFGStructure(fgs_name)
# fgs.atoms.write(f"FGS_{fgs_name}.xyz")