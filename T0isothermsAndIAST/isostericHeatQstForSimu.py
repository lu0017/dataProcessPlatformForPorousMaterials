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

# -----------------------------
# 定义 DSL 公式
# -----------------------------
def DSL_isotherm(P, qA, qB, bA, bB):
    """给定压力 P，计算 DSL 等温线吸附量"""
    return (qA * bA * P / (1 + bA * P) + qB * bB * P / (1 + bB * P))

def DSL_P(q_target, qA, qB, bA, bB):
    """数值求解 P: 根据 q_target 从 DSL 公式反算 P"""
    def func(P):
        return DSL_isotherm(P, qA, qB, bA, bB) - q_target
    try:
        P_root = brentq(func, 1e-12, 1e6)  # 搜索压力范围
    except ValueError:
        P_root = np.nan
    return P_root

def plotQst(Qst):
    # -----------------------------
    # 绘图
    # -----------------------------
    q = Qst["q_mmol_g"]
    qst = Qst["Qst_J_mol"]
    plt.figure()
    plt.plot(q, qst, marker="o", label="Qst")
    plt.xlabel("q (mmol/g)")
    plt.ylabel("Qst (kJ/mol)")
    plt.title("Isosteric Heat from DSL Data")
    plt.grid(True)
    plt.legend()
    plt.show(block=True)

def runCalculateQstByDSL(df):

    q_step = 0.1  # 吸附量步长
    q_point = 50

    # 构建每个温度的 DSL 参数字典
    qA = df[f"qA(mmol/g)"].iloc[0]
    qB = df[f"qB(mmol/g)"].iloc[0]

    T_list = df["T(K)"].values        # 温度数组
    bA_list = df["bA(1/kPa)"].values  # 对应温度的 bA
    bB_list = df["bB(1/kPa)"].values  # 对应温度的 bB

    q_min = 0.0
    q_min = 0.0381708

    # 计算各温度下的理论极限吸附量 q_max(T)（取 P=1e6 近似无穷大）
    q_max_list = []
    for T, bA, bB in zip(T_list, bA_list, bB_list):
        q_max_T = DSL_isotherm(1e6,  qA, qB, bA, bB)
        q_max_list.append(q_max_T)

    # q_max = min(q_max_list)  # 公共范围取最小值
    q_max = 3  # 公共范围取最小值
    q_range = np.arange(q_min, q_max, q_step)
    # q_range = np.linspace(q_min, q_max, q_point)

    print(f"全局吸附量范围: q_min={q_min:.3f}, q_max={q_max:.3f}")

    # -----------------------------
    # 计算 Qst
    # 目前只使用Clausius–Clapeyron来计算Qst
    # -----------------------------
    Qst_list = []
    valid_q = []
    results_full = []

    for q_target in q_range:
        P_at_q = []
        P_dict = {}   #用于导出excel
        skip = False
        for T, bA, bB in zip(T_list, bA_list, bB_list):
            P_q = DSL_P(q_target, qA, qB, bA, bB)
            if np.isnan(P_q):
                skip = True
                break
            P_at_q.append(P_q)
            col_name = f"P_{int(T)}K" if float(T).is_integer() else f"P_{T}K"
            P_dict[col_name] = P_q
        if skip:
            continue

        P_at_q = np.array(P_at_q)
        lnP = np.log(P_at_q)
        invT = 1 / np.array(T_list)
        slope, intercept, r_value, p_value, std_err = linregress(invT, lnP)
        Qst = -slope * const.R 
        Qst_list.append(Qst)
        valid_q.append(q_target)
        row = {
            "q_mmol_g": q_target,
            "Qst_J_mol": Qst,
            "Qst_kJ_mol": Qst/1000.0,
            "R2": r_value**2,
            "slope": slope,
            "intercept": intercept
        }
        row.update(P_dict)
        results_full.append(row)

    results = pd.DataFrame({"q_mmol_g": valid_q, "Qst_J_mol": Qst_list})
    # return results
    return pd.DataFrame(results_full)

def main(file_path=None):
    
    # 选择文件
    sheet_name = "DSL_params"
    file = fl.getFile()
    df_list, _, _ = fl.readFileBySheet(file, sheet_name, expand = const.QST_DSL_FILE)
    fileName = os.path.splitext(os.path.basename(file))[0]
    out_path = fl.get_expanded_name(file, fileName, expand=const.QST_DSL_FILE)
    df = pd.DataFrame(df_list)
    df_Qst = runCalculateQstByDSL(df)
    # 6. 保存结果
    fl.export_to_excel_auto(df_Qst, filename=out_path)
    plotQst(df_Qst)

if __name__ == "__main__":
    f = sys.argv[1] if len(sys.argv) > 1 else None
    main(f)