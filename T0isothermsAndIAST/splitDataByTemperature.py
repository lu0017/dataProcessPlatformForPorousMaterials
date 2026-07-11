import os
import sys
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from common import *
import constantsAndName as const

import T1fileSource.fileOperation as fl
import T1dataProcessSource.dataOperation as dop
from T0isothermsAndIAST import plotIsothermsByTem

def append_data_to_excel(listData, sample_name, out_path, firstSample=False):
    """
    将当前样品的四组温度数据写入到 out_path 的四个 sheet 中。
    - 若 sheet 不存在则新建；
    - 若样品已存在则覆盖；
    - 若不存在则追加；
    - 各列长度可不等。
    """
    T_Name= {}
    P_Name= {}
    q_Name= {}
    mapping = {}

    # # === 四组数据对应关系 ===

    # mapping = {
    #     '0C': ('0Absolute Pressure (kPa)', '0uptake (mmol/g)'),
    #     '15C': ('15Absolute Pressure (kPa)', '15uptake (mmol/g)'),
    #     '25C': ('25Absolute Pressure (kPa)', '25uptake (mmol/g)'),
    #     'N2-25C': ('N2_Absolute Pressure (kPa)', 'N2_25uptake (mmol/g)')
    # }
    
    _, _, T_list = dop.build_P_and_uptake_data(listData)
    T_list_C = T_list + const.K2C
    for T in T_list_C:
        T_Name[T] = dop.getName(T, unitChange = 0, prefix=None, suffix = "C")
        P_Name[T] = dop.getName(T, unitChange = 0, prefix=None, suffix = "Absolute Pressure (kPa)")
        q_Name[T] = dop.getName(T, unitChange = 0, prefix=None, suffix = "uptake (mmol/g)")
        # 更新 mapping，每个 T 对应 (P列名, q列名)
        mapping[T_Name[T]] = (P_Name[T], q_Name[T])

    N2_name = dop.getName(T_list_C.max(), unitChange = 0, prefix="N2-", suffix = "C")
    N2_Pname = dop.getName(T_list_C.max(), unitChange = 0, prefix="N2_", suffix = "Absolute Pressure (kPa)")
    N2_qname = dop.getName(T_list_C.max(), unitChange = 0, prefix="N2_", suffix = "uptake (mmol/g)")
    # 将 N2 添加到 mapping
    mapping[N2_name] = (N2_Pname, N2_qname)

    # 先删除旧文件（确保干净）
    if firstSample and os.path.exists(out_path):
        os.remove(out_path)

   # === 文件与工作簿准备 ===
    if not os.path.exists(out_path):
        # 🆕 创建新文件时保留默认 sheet，避免空工作簿错误
        wb = Workbook()
        default_sheet = wb.active
        default_sheet.title = "Init"
        wb.save(out_path)
        print(f"🆕 创建新文件: {out_path}")

    wb = load_workbook(out_path)

    # === 处理每个温度的 sheet ===
    for sheet_name, (p_col, u_col) in mapping.items():
        # 若当前温度数据不存在于 listData 中则跳过
        if p_col not in listData or u_col not in listData:
            continue

        pressure = listData[p_col]
        uptake = listData[u_col]

        # 获取或创建 sheet
        if sheet_name not in wb.sheetnames:
            ws = wb.create_sheet(sheet_name)

        else:
            ws = wb[sheet_name]

        # === 检查是否已存在该样品列 ===
        existing_samples = [ws.cell(row=1, column=c).value for c in range(2, ws.max_column + 1)]
        if sample_name in existing_samples:
            col_idx = existing_samples.index(sample_name) + 2  # +2 因为第一列是Sample标识列
            # 删除旧数据列
            ws.delete_cols(col_idx, 2)
            print(f"🧹 已删除旧列: {sheet_name} - {sample_name}")
        else:
            col_idx = ws.max_column + 1


        # === 写入样品名与数据 ===
        ws.cell(row=1, column=col_idx, value=sample_name)
        ws.cell(row=2, column=col_idx, value="Pressure (kPa)")
        ws.cell(row=2, column=col_idx + 1, value=f"uptake (mmol/g)")

        for i, val in enumerate(pressure, start=3):
            ws.cell(row=i, column=col_idx, value=val)
        for i, val in enumerate(uptake, start=3):
            ws.cell(row=i, column=col_idx + 1, value=val)

        print(f"✅ 写入完成: {sheet_name} - {sample_name}")

    wb.save(out_path)

def splitDataByTemperature(file_path):
    xls = pd.ExcelFile(file_path)
    fileName = os.path.splitext(os.path.basename(file_path))[0]
    out_path = fl.get_expanded_name(file_path, fileName, expand=const.SPLIT_FILE, expandPos=True, type="xlsx")
    firstSample = True

    for sheet in xls.sheet_names:
        if fl.should_skip(sheet):
            continue

        listData, _, _ = fl.readFileBySheet(file_path, sheet)

        if isinstance(listData, dict):
            listData0 = list(listData.values())[0]  # 取第一个 sheet 数据结构

        if not listData0 or len(listData0) <= 1:
            print(f"⚠️ {sheet}: 无有效数据，跳过。")
            continue

        # ✅ 每读取一个sheet，调用一次保存函数
        append_data_to_excel(listData, sheet, out_path, firstSample)
        firstSample = False
    return out_path


def main(file_path=None):

    file_path = fl.getFile()
    out_path = splitDataByTemperature(file_path)
    plotIsothermsByTem(out_path)

if __name__ == "__main__":
    f = sys.argv[1] if len(sys.argv) > 1 else None
    main(f)