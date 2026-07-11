##每一个子平台固定开头，用于找到依赖
import os
import sys
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from common import *
import T1fileSource.fileOperation as fl
import constantsAndName as const


def calcuSelectivity(df):
     # ---- 动态提取列 ----
    pressure = df["Pressure"]
    uptake = df["Uptake"]

    # 找到所有 y(...) 列
    y_cols = [c for c in df.columns if c.startswith("y(")]
    # 找到所有 x(...) 列
    x_cols = [c for c in df.columns if c.startswith("x(")]

    if len(y_cols) >= 1:
        co2_y = df[y_cols[0]]  # 假设第一个 y 是 CO2
    else:
        co2_y = None

    if len(x_cols) >= 1:
        co2_x = df[x_cols[0]]  # 第一个 x 是 CO2
    if len(x_cols) >= 2:
        n2_x = df[x_cols[1]]  # 第二个 x 是 N2

    uptake_CO2 = co2_x * uptake
    uptake_N2 = n2_x * uptake
    selectivity = (uptake_CO2 / co2_y) / (uptake_N2 / (1-co2_y))

    # ---- 组织成新的 DataFrame ----
    df_out = pd.DataFrame({
        "Pressure": pressure,
        "y(CO2)": co2_y,
        "Uptake": uptake,
        "x(CO2)": co2_x,
        "x(N2)": n2_x,
        "CO2_uptake（mmol/g）": uptake_CO2,
        "N2_uptake（mmol/g）": uptake_N2,
        "CO2/N2 selectivity": selectivity
    })
    return df_out

def poltP2Sele(select):
    plt.figure(figsize=(6,4))
    plt.plot(select["Pressure"], select["CO2/N2 selectivity"], marker="o", label="Selectivity")

    plt.xlabel("Pressure (kPa)")
    plt.ylabel("Selectivity (-)")
    plt.title("P–Selectivity Curve")

    # 在图上写 y_co2
    y_value = select["y(CO2)"].iloc[0]   # 取第一个（如果是常数列）
    plt.text(
        x=0.7*select["Pressure"].max(),  # 横坐标：选在右上角
        y=0.9*select["CO2/N2 selectivity"].max(),  # 纵坐标：选在上方
        s=f"y_CO2 = {y_value:.2f}", 
        fontsize=10,
        bbox=dict(facecolor="white", alpha=0.6, edgecolor="gray")
    )
    plt.legend()
    plt.grid(True)  
    plt.show(block=True)

def main(file_path=None):

    file = fl.getFile(default_type="txt")
    df, out_path =fl.readFileByName(file, expand=const.SELE_FILE_IAST_SOFT)

    sele = calcuSelectivity(df)
    fl.export_to_excel_auto(sele, out_path)

    poltP2Sele(sele)

if __name__ == "__main__":
    f = sys.argv[1] if len(sys.argv) > 1 else None
    main(f)

# 举例提取不同的列
