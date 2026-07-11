"""unit change factor"""
bar2kPa = 100
kPa2Pa = 1000
bar2Pa = 100000
atm2Pa = 101325
mol2mmol = 1000
C2K = 273.15
kJ2J = 1000

Pa2bar = 1 / bar2Pa
Pa2kPa = 1 / kPa2Pa
kPa2bar = 1 / bar2kPa
Pa2atm = 1 / atm2Pa
mmol2mol = 1 / mol2mmol
K2C = -C2K
J2kJ = 1 / kJ2J

"""constants"""
R = 8.314  # J/mol·K
NA = 6.02214076e23  # 阿伏伽德罗常数
M_CO2 = 44.01 # CO2摩尔质量 Molar mass of CO₂
Tb_C = -56.558 #CO2 三相界面温度，boiling point
Tb_K = Tb_C + C2K #CO2 三相界面温度，boiling point
rho_b = 1178.53 #CO2 三相界面液态密度，kg/m3 ## 来源Thermophysical Properties of Carbon Dioxide and CO2-Rich Mixtures
alph = 0.0025 #体积膨胀系数 isosteric coefficient of the expansion of the adsorbed volume


P0_kPa = { ##CO2在各温度下的饱和压强
    273.15:  # 0°C
        3485.1,  # 例子，单位 kPa（你填自己的）
    288.15:  # 15°C
        5087.1,
    298.15:  # 25°C
        6434.2,
    303.15:  # 30°C
        7213.7,
}
def get_P0(T_K):
    """unit is kPa"""
    return P0_kPa[T_K]

"""variables"""



"""folder"""
SELE_FOLDER_IAST_SOFT = "seleByiast++"
SELE_FOLDER_IAST_PY = "SelecByPyIast"
DSL_FOLDER = "dslFitting"
DA_FOLDER = "DAFitting"
QST_FOLDER = "Qst"

"""expand file name"""
SELE_FILE_IAST_SOFT = "SelecByIast++"
SELE_FILE_IAST_PY = "SelecByPyIast"
DSL_FILE = "DSL_fit"
DA_FILE = "DA_fit"
TSA_FILE = "TSA"
QST_DSL_FILE = "Qst_dsl"
SPLIT_FILE = "split"

"""sheet name"""
DSL_PARA_SHEET = "DSL_params"
DSL_LOW_PARA_SHEET = "DSL_lowpressure_params"
DSL_FITS_SHEET = "DSL_fits"
DSL_LOW_FITS_SHEET = "DSL_lowpressure_fits"
DA_PARA_SHEET = "DA_params"
DA_BYT_PARA_SHEET = "DA_params_byT"
DA_FITS_SHEET = "DA_fits"
DA_BYT_FITS_SHEET = "DA_fits_byT"
THERMO_SHEET = "thermo_df"
GIBBS_SHEET = "dg_df"
N2_FITS_SHEET = "N2_fit"
QST_EXP_SHEET = "Qst_exp"
QST_SIM_SHEET = "Qst_simul"
QST_SIM_CC_SHEET = "Qst_simul_CC"

