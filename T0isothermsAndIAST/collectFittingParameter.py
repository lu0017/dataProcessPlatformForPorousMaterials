##平台固定开头，用于找到依赖
import os
import sys
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from common import *
import constantsAndName as const
import T1fileSource.fileOperation as fl
import T1dataProcessSource.dataOperation as dop
import T0isothermsAndIAST.modelAndFit as mf
import T0isothermsAndIAST.strongDSL as sl
import T0isothermsAndIAST.isostericHeatQstForRawData as qstByRawdata
# ===============================
# 文件夹
# ===============================
def extractModelParameter(data, model, prefix, out_path):
    results = []
    for name, df in data.items():
        sample = name.replace(prefix, "")
        if model == "DSL":
            df = df.set_index("T(K)")
            results.append({
                "sample": sample,
                "qA": df.iloc[0]["qA(mmol/g)"],
                "qB": df.iloc[0]["qB(mmol/g)"],
                "bA-T0": df.loc[273.15, "bA(1/kPa)"],
                "bA-T15": df.loc[288.15, "bA(1/kPa)"],
                "bA-T25": df.loc[298.15, "bA(1/kPa)"],
                "bB-T0": df.loc[273.15, "bB(1/kPa)"],
                "bB-T15": df.loc[288.15, "bB(1/kPa)"],
                "bB-T25": df.loc[298.15, "bB(1/kPa)"],
            })
        elif model == "D-A":
            results.append({
                "sample": sample,
                "W0 (cm3/g)": df.iloc[0]["W0 (cm3/g)"],
                "E (J/mol)": df.iloc[0]["E (J/mol)"],
                "n": df.iloc[0]["n"],
            })
        elif model == "Langmuir":
            # 找到 N2 的那一行
            row = df[df["gas"] == "N2"].iloc[0]
            results.append({
                "sample": sample,
                "q-N2 (mmol/g)": row["q(mmol/g)"],
                "b-N2 (1/kPa)": row["b(1/kPa)"],
            })
        elif model == "IAST":
            row = df[df["P (kPa)"] == 100].iloc[0]
            results.append({
                "sample": sample,
                "SELE pyIAST": row["CO2/N2 selectivity"],
            })
        elif model == "Henry":
            results.append({
                "sample": sample,
                "sele_henry_model": df.iloc[0]["sele_henry_model"],
            })
    summary = pd.DataFrame(results)
    summary = dop.naturalSortBy(summary, by="sample")
    fl.export_to_excel_auto( summary, filename=out_path, sheet_name=f"{model}_params" )
# ==========================================
# 主函数
# ==========================================
def main(file_path=None):
    models = {
        "DSL": {
            "folder": "dslFitting",
            "sheet": "DSL_params",
            "prefix": "DSL_fit_",
        },
        "D-A": {
            "folder": "DAFitting",
            "sheet": "DA_params",
            "prefix": "DA_fit_",
        },
        "Langmuir": {
                    "folder": "SelecByPyIast",
                    "sheet": "Model_params",
                    "prefix": "SelecByPyIast_",
                },
        "IAST": {
                    "folder": "SelecByPyIast",
                    "sheet": "P-S",
                    "prefix": "SelecByPyIast_",
                },
        "Henry": {
                    "folder": "SelecByPyIast",
                    "sheet": "S-HENRY",
                    "prefix": "SelecByPyIast_",
                },
    }
    # 选择文件夹
    root_folder = fl.select_folder()
    fileName = "kinetic"
    out_path = fl.get_expanded_name( root_folder, fileName=fileName, level=0 )
    for model, config in models.items():
        folder = os.path.join( root_folder, config["folder"] )
        data = fl.readFolderData( folder, extension=".xlsx", sheet_name=config["sheet"] )
        extractModelParameter( data=data, model=model, prefix=config["prefix"], out_path=out_path, )
if __name__ == "__main__":
    f = sys.argv[1] if len(sys.argv) > 1 else None
    main(f)