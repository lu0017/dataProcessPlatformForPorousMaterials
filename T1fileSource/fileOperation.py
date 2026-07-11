from common import *
import constantsAndName as const
warnings.filterwarnings("ignore", category=RuntimeWarning)
# 默认列名映射（推荐写成最终想要在 Excel 看到的格式）
# -----------------------
# 默认列名映射（最终在 Excel 中显示的格式）
# -----------------------
DEFAULT_LABELS = {
    # Gibbs / Helmholtz
    "DeltaG_Jmol": "ΔG (J/mol)",
    "DeltaG_kJmol": "ΔG (kJ/mol)",
    "DeltaG_Jmol_byHS": "ΔG_byHS (J/mol)",
    "DeltaG_kJmol_byHS": "ΔG_byHS (kJ/mol)",
    # 热力学
    "dH_Jmol": "ΔH (J/mol)",
    "dH_kJmol": "ΔH (kJ/mol)",
    "dS_Jmol_K": "ΔS (J/mol·K)",
    # DSL 参数
    "q1": "qA (mmol/g)",
    "q2": "qB (mmol/g)",
    "b1": "bA (1/Pa)",
    "b2": "bB (1/Pa)",
    # D-A 参数
    "rho_mmol_cm3": "ρ (mmol/cm3)",
    "W0_cm3_g": "W0 (cm3/g)",
    "W_cm3_g": "W (cm3/g)",
    "E_J_mol": "E (J/mol)",
    "A_J_mol": "A (J/mol)",
    "P0_kPa": "P0 (kPa)",
    "n": "-",
    # 吸附量
    "q_mmol_g": "q (mmol/g)",
    "q_exp_mmol_g": "q_exp (mmol/g)",
    "q_fit_mmol_g": "q_fit (mmol/g)",
    # 压力 / 温度
    "P_Pa": "P (Pa)",
    "P_kPa": "P (kPa)",
    "T_K": "T (K)",
    "T_C": "T (°C)",
    # 选择性 /拟合相关
    "R2": "R²",
    "slope": "slope",
    "intercept": "intercept",
    "n": "n",
    "theta": "θ",
    # site 标识
    "site": "site",
    # Qst
    "Qst_J_mol": "Qst (J/mol)",
    "Qst_kJ_mol": "Qst (kJ/mol)",
    # 额外压力列表
    "P_list": "P_list (Pa)",
    "P_map": "P_map (Pa)",
    "P_list_str": "P_list_str (Pa)",
    # Raman
    "D_center": "D peak position (cm⁻¹)",
    "G_center": "G peak position (cm⁻¹)",
    "D_height": "D peak height (a.u.)",
    "G_height": "G peak height (a.u.)",
    "D_area": "D peak area (a.u.·cm⁻¹)",
    "G_area": "G peak area (a.u.·cm⁻¹)",
    "D_fwhm": "D FWHM (cm⁻¹)",
    "G_fwhm": "G FWHM (cm⁻¹)",
    "ID/IG": "ID/IG",
    "AD/AG": "AD/AG",
}
# -----------------------
# 默认单位字典（可用于拼接列名，如果 DEFAULT_LABELS 未指定显示名）
# -----------------------
DEFAULT_UNITS = {
    "b1": "1/Pa",
    "b2": "1/Pa",
    "DeltaG_Jmol": "J/mol",
    "DeltaG_kJmol": "kJ/mol",
    "DeltaG_Jmol_byHS": "J/mol",
    "DeltaG_kJmol_byHS": "kJ/mol",
    "dH_Jmol": "J/mol",
    "dH_kJmol": "kJ/mol",
    "dS_Jmol_K": "J/mol*K",
    "q1": "mmol/g",
    "q2": "mmol/g",
    "q_mmol_g": "mmol/g",
    "q_exp_mmol_g": "mmol/g",
    "q_fit_mmol_g": "mmol/g",
    "Qst_J_mol": "J/mol",
    "Qst_kJ_mol": "kJ/mol",
    "P_Pa": "Pa",
    "P_kPa": "kPa",
    "P_list": "Pa",
    "P_list_str": "Pa",
    "P_map": "Pa",
    "T_K": "K",
    "T_C": "°C",
    "R2": "-",
    "slope": "-",
    "intercept": "-",
    "n": "-",
    "theta": "-",
    "site": "-",
    "rho_mmol_cm3": "mmol/cm3",
    "W0_cm3_g": "cm3/g",
    "W_cm3_g": "cm3/g",
    "E_J_mol": "J/mol",
    "A_J_mol": "J/mol",
    "P0_kPa": "kPa",
}
# ==========================================
# 工具函数
# ==========================================
SkipSheets = [
    "总表",
    "说明",
    "参考数据",
    "SSA",
    "CC-BiCa-800-2-0.5",
    "CC-HyFe-800-2-1",
    "0C",
    "15C",
    "25C",
    "uio-66 Zn 40",
    "uio-66 Co 40"
]
splitSheet = [
    "0C",
    "15C",
    "25C",
    "N2-25C"
]
def validSheetforSplitFile(sheet_name: str) -> bool:
    return sheet_name in splitSheet
def should_skip(sheet_name: str) -> bool:
    return sheet_name in SkipSheets
def getFile(file_path, newName=True, sheet_name="Sheet1"):
    if not file_path:
        import tkinter as tk
    from tkinter import filedialog
    root=tk.Tk()
    root.withdraw()
    file_path=filedialog.askopenfilename(title="请选择数据文件",filetypes=[("Excel or CSV","*.xlsx *.xls *.csv")])
    if not file_path: return            
    print(f"已选择文件: {file_path}")
    # 读取
    ext=os.path.splitext(file_path)[1].lower()
    if ext in [".xlsx",".xls"]:
        df=pd.read_excel(file_path, sheet_name=sheet_name)
    elif ext==".csv": 
        df=pd.read_csv(file_path)
    else: 
        raise ValueError("只支持 Excel/CSV")
    output_dir = os.path.dirname(file_path)
    filename = os.path.splitext(os.path.basename(file_path))[0]
    if newName:
        out_path = os.path.join(output_dir, f"DSL_fit_{filename}.xlsx")
    else:
        out_path = file_path
    return df, out_path
def select_folder():
    root = tk.Tk()
    root.withdraw()
    return filedialog.askdirectory(title="请选择文件夹")
def getFile(file_path=None, default_type="excel"):
    """
    打开文件选择窗口，支持 Excel / CSV / TXT 文件。
    参数:
        file_path (str): 已有文件路径（可选）。如果提供，则直接读取该路径。
        default_type (str): 默认文件类型，可选值:
                            "excel" -> *.xlsx, *.xls
                            "csv"   -> *.csv
                            "txt"   -> *.txt
    返回:
        xls_or_df, file_path
        - ExcelFile 或 DataFrame 对象
        - 文件路径字符串
    """
    # 若未传入路径，则弹出对话框
    if not file_path:
        root = tk.Tk()
        root.withdraw()
        root.attributes("-topmost", True)  # ✅ 让窗口置顶
        # 设置文件过滤类型
        if default_type.lower() == "excel":
            filetypes = [("Excel 文件", "*.xlsx *.xls"), ("所有文件", "*.*")]
        elif default_type.lower() == "csv":
            filetypes = [("CSV 文件", "*.csv"), ("所有文件", "*.*")]
        elif default_type.lower() == "txt":
            filetypes = [("文本文件", "*.txt"), ("所有文件", "*.*")]
        else:
            filetypes = [("所有文件", "*.*")]
        file_path = filedialog.askopenfilename(
            title=f"请选择{default_type.upper()}文件",
            filetypes=filetypes
        )
        if not file_path:
            print("未选择文件。")
            return None, None
    # 打印选择结果
    print(f"已选择文件: {file_path}")  
    return file_path
def getPathOnly(file_path):
    if os.path.isdir(file_path):
        path = file_path  # 本身就是文件夹
    else:
        path = os.path.dirname(file_path)  # 文件 → 返回父目录
    return path
def check_or_create_folder(path, sub_folder_name=""):
    folder_path0 = getPathOnly(path)
    folder_path = os.path.join(folder_path0, sub_folder_name)
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    #     print(f"✅ 创建新文件夹: {folder_path}")
    # else:
    #     print(f"📁 文件夹已存在: {folder_path}")
    return folder_path
def get_expanded_name(file_path, fileName, expand="", expandPos=True, type="xlsx"):
    """用于拓展文件名"""
    output_dir = getPathOnly(file_path)
    if expand == "":
        out_path = os.path.join(output_dir, f"{fileName}.{type}")
    else:
        if expandPos:
            out_path = os.path.join(output_dir, f"{expand}_{fileName}.{type}")
        else:
            out_path = os.path.join(output_dir, f"{fileName}_{expand}.{type}")
    return  out_path
def readData( file, sheet_name=0, keyFlag=True, keyName=None, ncols=None, sep=None, ):
    """
    通用表格读取函数
    Parameters
    ----------
    file : str
        文件路径
    sheet_name : str/int, default=0
        Excel sheet名称或序号，仅xlsx有效
    keyFlag : bool, default=True
        True：第一行为列名
        False：无列名
    keyName : list or None
        指定列名。
        None 表示：
            keyFlag=True 使用原列名
            keyFlag=False 自动生成 Column1...
    ncols : int or None
        仅读取前 n 列
    sep : str or None
        txt/csv 分隔符。
        None 表示自动识别。
    Returns
    -------
    DataFrame
    """
    suffix = os.path.splitext(file)[1].lower()
    header = 0 if keyFlag else None
    usecols = None if ncols is None else range(ncols)
    # ---------------- Excel ----------------
    if suffix in [".xlsx", ".xls"]:
        df = pd.read_excel(
            file,
            sheet_name=sheet_name,
            header=header,
            usecols=usecols,
        )
    # ---------------- CSV ----------------
    elif suffix == ".csv":
        df = pd.read_csv(
            file,
            header=header,
            usecols=usecols,
            sep="," if sep is None else sep,
        )
    # ---------------- TXT ----------------
    elif suffix == ".txt":
        if sep is None:
            sep = r"\s+"
        df = pd.read_csv(
            file,
            header=header,
            usecols=usecols,
            sep=sep,
            engine="python",
        )
    else:
        raise ValueError(f"Unsupported file type: {suffix}")
    # 自动列名
    if keyName is not None:
        df.columns = keyName
    elif not keyFlag:
        df.columns = [f"Column{i+1}" for i in range(df.shape[1])]
    # 尝试全部转换成数值
    df = df.apply(pd.to_numeric, errors="ignore")
    return df
def readFolderData( folder, extension=".csv", keyFlag=True, keyName=None, ncols=None, sheet_name=0, sep=None, ):
    """
    批量读取文件夹中的数据。
    Returns
    -------
    dict
        {
            sample1 : DataFrame,
            sample2 : DataFrame,
            ...
        }
    """
    data = {}
    files = sorted( glob.glob( os.path.join(folder, f"**/*{extension}"), recursive=True, ) )
    for file in files:
        sample = os.path.splitext(os.path.basename(file))[0]
        df = readData( file=file, sheet_name=sheet_name, keyFlag=keyFlag, keyName=keyName, ncols=ncols, sep=sep, )
        data[sample] = df
    return data
def readFileBySheet(file_path, sheet_name, expand="", type="xlsx", header=0):
    """用sheet名，只允许excel, 返回dic格式"""
    df = []
    out_path = []
    validFile = True
    cols = {}
    xls = pd.ExcelFile(file_path)
    if sheet_name not in xls.sheet_names:
        print(f"Sheet '{sheet_name}' 不存在")
        validFile = False
    else:
        df = pd.read_excel(file_path, sheet_name=sheet_name, header= header)
        df.dropna(how="all", inplace=True)
        for c in df.columns:
            valid_vals = df[c].dropna().tolist()
            cols[c] = valid_vals
        output_dir = os.path.dirname(file_path) #without file name
        out_path = get_expanded_name(output_dir, sheet_name, expand=expand, type=type)
    return cols, out_path, validFile
def readTableBySheet(file_path, sheet_name, auto_header=True, auto_index=True):
    """
    Read a statistics table from an Excel sheet.
    Features
    --------
    ✓ Automatically remove empty rows/columns
    ✓ Automatically detect header row
    ✓ Automatically detect sample index column
    ✓ Automatically convert numeric columns
    ✓ Keep strings unchanged
    Parameters
    ----------
    file_path : str
        Excel file path.
    sheet_name : str
        Sheet name.
    auto_header : bool
        Automatically detect the header row.
    auto_index : bool
        Automatically use the first non-numeric column as index.
    Returns
    -------
    data : pandas.DataFrame
    meta : dict
        Information about the detected table.
    """
    # -----------------------------
    # Check file
    # -----------------------------
    if not os.path.exists(file_path):
        raise FileNotFoundError(file_path)
    xls = pd.ExcelFile(file_path)
    if sheet_name not in xls.sheet_names:
        raise ValueError(f"Sheet '{sheet_name}' does not exist.")
    # -----------------------------
    # Read raw sheet
    # -----------------------------
    raw = pd.read_excel(file_path,
                        sheet_name=sheet_name,
                        header=None)
    # remove empty rows/columns
    raw.dropna(axis=0, how="all", inplace=True)
    raw.dropna(axis=1, how="all", inplace=True)
    raw.reset_index(drop=True, inplace=True)
    # -----------------------------
    # Detect header
    # -----------------------------
    header_row = 0
    if auto_header:
        for i in range(len(raw)):
            row = raw.iloc[i]
            numeric = pd.to_numeric(row, errors="coerce")
            valid = numeric.notna().sum()
            # First row with few numbers is usually header
            if valid < len(row) * 0.5:
                # Next row should contain more numeric values
                if i + 1 < len(raw):
                    next_numeric = pd.to_numeric( raw.iloc[i + 1], errors="coerce" )
                    if next_numeric.notna().sum() >= len(row) * 0.5:
                        header_row = i
                        break
    # -----------------------------
    # Build dataframe
    # -----------------------------
    columns = raw.iloc[header_row].ffill().tolist()
    data = raw.iloc[header_row + 1:].copy()
    data.columns = columns
    data.reset_index(drop=True, inplace=True)
    # -----------------------------
    # Convert numeric columns
    # -----------------------------
    for col in data.columns:
        converted = pd.to_numeric(data[col], errors="coerce")
        # if most values are numeric -> convert
        if converted.notna().sum() >= len(data) * 0.8:
            data[col] = converted
    # -----------------------------
    # Detect sample column
    # -----------------------------
    index_name = None
    if auto_index:
        for col in data.columns:
            if not pd.api.types.is_numeric_dtype(data[col]):
                index_name = col
                data.set_index(col, inplace=True)
                break
    # -----------------------------
    # Meta
    # -----------------------------
    meta = {
        "header_row": header_row,
        "index": index_name,
        "shape": data.shape,
        "columns": list(data.columns)
    }
    return data, meta
def readFileBySheetWithMultiLevelHeader(file_path, sheet_name, expand="", type="xlsx", auto_header=True):
    """
    Read Excel sheet and automatically detect multi-level headers.
    Parameters
    ----------
    file_path : str
        Excel file path.
    sheet_name : str
        Sheet name.
    auto_header : bool
        Automatically detect header rows.
    Returns
    -------
    df : pandas.DataFrame
        DataFrame with MultiIndex columns (if multiple header rows exist).
    out_path : str
    validFile : bool
    meta : dict
        Information about header.
    """
    validFile = True
    if not os.path.exists(file_path):
        print("File does not exist.")
        return None, None, False, {}
    xls = pd.ExcelFile(file_path)
    if sheet_name not in xls.sheet_names:
        print(f"Sheet '{sheet_name}' 不存在")
        return None, None, False, {}
    # -----------------------------
    # Read without header
    # -----------------------------
    raw = pd.read_excel( file_path, sheet_name=sheet_name, header=None )
    # remove empty rows/columns
    raw.dropna(axis=0, how="all", inplace=True)
    raw.dropna(axis=1, how="all", inplace=True)
    raw.reset_index(drop=True, inplace=True)
    # -----------------------------
    # Auto detect header rows
    # -----------------------------
    header_rows = 1
    if auto_header:
        header_rows = None
        for i in range(len(raw)):
            row = raw.iloc[i]
            # 判断本行是否开始出现大量数字
            numeric = pd.to_numeric(row, errors="coerce")
            valid = numeric.notna().sum()
            # 超过50%的列都是数字
            if valid >= max(1, len(row) * 0.5):
                header_rows = i
                break
        if header_rows is None:
            header_rows = 1
    # -----------------------------
    # Build Header
    # -----------------------------
    header = raw.iloc[:header_rows].copy()
    # fill merged cells
    header = header.ffill(axis=1)
    header = header.ffill(axis=0)
    # build MultiIndex
    if header_rows == 1:
        columns = header.iloc[0].tolist()
    else:
        arrays = []
        for i in range(header_rows):
            arrays.append(header.iloc[i].tolist())
        columns = pd.MultiIndex.from_arrays(arrays)
    # -----------------------------
    # Data
    # -----------------------------
    data = raw.iloc[header_rows:].copy()
    data.columns = columns
    data.reset_index(drop=True, inplace=True)
    # convert numeric when possible
    data = data.apply(pd.to_numeric, errors="ignore")
    # -----------------------------
    # Output Path
    # -----------------------------
    output_dir = os.path.dirname(file_path)
    out_path = get_expanded_name( output_dir, sheet_name, expand=expand, type=type )
    # -----------------------------
    # Meta
    # -----------------------------
    meta = {
        "header_rows": header_rows,
        "column_levels": 1 if header_rows == 1 else header_rows,
        "shape": data.shape
    }
    return data, out_path, validFile, meta
def readFileByName(file_path, expand):
    """用文件名,读取多种类型"""
    if file_path.endswith((".xlsx", ".xls")):
        df = pd.read_excel(file_path)
    elif file_path.endswith(".csv"):
        df = pd.read_csv(file_path)
    elif file_path.endswith(".txt"):
        try:
            # 优先尝试制表符分隔
            df = pd.read_csv(file_path, sep="\t")
        except:
            # 如果失败，用空格分隔
            df = pd.read_csv(file_path, delim_whitespace=True)
    else:
        raise ValueError("不支持的文件格式: " + file_path)
    df.dropna(inplace=True)
    output_dir = os.path.dirname(file_path) #without file name
    filename = os.path.splitext(os.path.basename(file_path))[0]
    out_path = get_expanded_name(output_dir, filename, expand=expand)
    return df, out_path
def _excel_exists_and_valid(filename):
    """
    判断 Excel 文件是否存在且为合法 xlsx。
    """
    if not os.path.exists(filename):
        return False
    try:
        load_workbook(filename, read_only=True)
        return True
    except (BadZipFile, OSError, ValueError):
        return False
def export_to_excel_auto(df, filename="results.xlsx", sheet_name=None):
    """
    自动导出 DataFrame 到 Excel
    - 优先使用 DEFAULT_LABELS（含别名+单位）
    - 其次用 DEFAULT_UNITS 拼接
    - 否则保持原名
    - sheet_name None 时自动生成 Sheet1, Sheet2 ...
    - 若 sheet 已存在则覆盖
    """
    # 自动生成 sheet 名
    if sheet_name is None:
        base_name = "Sheet"
        idx = 1
        # 先判断Excel是否存在且可正常打开
        workbook = None
        if os.path.exists(filename):
            try:
                workbook = load_workbook(filename, read_only=True)
            except (FileNotFoundError, BadZipFile, OSError, ValueError):
                # 文件存在但不是合法Excel，后面会重新创建
                workbook = None
        while True:
            candidate = f"{base_name}{idx}"
            if workbook is None:
                sheet_name = candidate
                break
            if candidate not in workbook.sheetnames:
                sheet_name = candidate
                break
            idx += 1
    # 列名替换逻辑
    new_columns = {}
    for col in df.columns:
        if col in DEFAULT_LABELS:
            new_columns[col] = DEFAULT_LABELS[col]
        elif col in DEFAULT_UNITS:
            new_columns[col] = f"{col} ({DEFAULT_UNITS[col]})"
        else:
            new_columns[col] = col
    df_export = df.rename(columns=new_columns)
    # 写入 Excel
    excel_ok = _excel_exists_and_valid(filename)
    # 如果存在但损坏，直接删除重新创建
    if os.path.exists(filename) and not excel_ok:
        print(f"检测到损坏的 Excel，重新创建：{filename}")
        os.remove(filename)
    if excel_ok:
        with pd.ExcelWriter( filename, engine="openpyxl", mode="a", if_sheet_exists="replace", ) as writer:
            df_export.to_excel(writer, sheet_name=sheet_name, index=False)
    else:
        with pd.ExcelWriter( filename, engine="openpyxl", mode="w", ) as writer:
            df_export.to_excel(writer, sheet_name=sheet_name, index=False)
    print(f"已导出 sheet: {sheet_name} -> {filename}")
def saveToExcel(df, out_path, sheet_name=None):
    """
    保存 DataFrame 到 Excel
    - sheet_name=None -> 覆盖整个文件
    - sheet_name=str -> 保存到指定 sheet（存在则覆盖，不存在则创建）
    """
    if sheet_name is None:
        # 不指定 sheet -> 覆盖整个 Excel
        df.to_excel(out_path, index=False)
        print(f"文件已保存: {out_path}")
        return
    if df.empty:
        print(f"跳过空 sheet: {sheet_name}")
        return
    # 写入 Excel
    if os.path.exists(out_path):
        with pd.ExcelWriter(out_path, engine="openpyxl", mode="a", if_sheet_exists="replace") as writer:
            df.to_excel(writer, sheet_name=sheet_name, index=False)
    else:
        df.to_excel(out_path, sheet_name=sheet_name, index=False)
    # # 导出 DSL 拟合结果
    # export_to_excel_auto(fits_data, "results.xlsx", sheet_name="DSL_fits")
def saveDataframe2Excel( df, sample_name, out_path, sheet_name="Data", start_col=1, overwrite=True, ):
    """
    将 DataFrame 保存到 Excel。
    Excel 格式
    ----------
    Row1 : Sample Name
    Row2 : Column Name
    Row3~ : Data
    Parameters
    ----------
    df : pandas.DataFrame
    sample_name : str
    out_path : str
    sheet_name : str, default="Data"
    firstSample : bool, default=False
        True 时若文件已存在则删除重新创建。
    start_col : int, default=1
        起始写入列（Excel 从1开始）。
    """
    # ==================================================
    # 创建工作簿
    # ==================================================
    if not os.path.exists(out_path):
        wb = Workbook()
        ws = wb.active
        ws.title = sheet_name
        wb.save(out_path)
    wb = load_workbook(out_path)
    # ==================================================
    # Sheet
    # ==================================================
    if sheet_name in wb.sheetnames:
        ws = wb[sheet_name]
    else:
        ws = wb.create_sheet(sheet_name)
    ncols = df.shape[1]
    # ==================================================
    # 查找是否已有该 Sample
    # ==================================================
    col_idx = None
    for c in range(start_col, ws.max_column + 1):
        if ws.cell(row=1, column=c).value == sample_name:
            if overwrite:
                col_idx = c
                ws.delete_cols(col_idx, ncols)
                print(f"🧹 已覆盖样品：{sample_name}")
            else:
                print(f"⏩ 已存在样品，跳过：{sample_name}")
                wb.save(out_path)
                return
            break
    # ==================================================
    # 新样品
    # ==================================================
    if col_idx is None:
        if ws.max_column < start_col:
            col_idx = start_col
        else:
            col_idx = max(start_col, ws.max_column + 1)
    # ==================================================
    # Sample Name
    # ==================================================
    ws.cell(row=1, column=col_idx, value=sample_name)
    # ==================================================
    # Column Name
    # ==================================================
    for j, name in enumerate(df.columns):
        ws.cell( row=2, column=col_idx + j, value=name, )
    # ==================================================
    # Data
    # ==================================================
    for i, row in enumerate(df.itertuples(index=False), start=3):
        for j, value in enumerate(row):
            ws.cell( row=i, column=col_idx + j, value=value, )
    wb.save(out_path)
    print(f"✅ 保存完成：{sample_name}")
def exportCorrelationExcel( filename, X, Y, summary, X_matrix=None, Y_matrix=None, summary_matrix=None ):
    export_to_excel_auto(pd.DataFrame(summary), filename, "summary")
    if X_matrix is not None:
        export_to_excel_auto( pd.DataFrame(X_matrix, index=X, columns=X), filename, "X_matrix" )
    if Y_matrix is not None:
        export_to_excel_auto( pd.DataFrame(Y_matrix, index=Y, columns=Y), filename, "Y_matrix" )
    if summary_matrix is not None:
        export_to_excel_auto( pd.DataFrame(summary_matrix), filename, "summary_matrix" )