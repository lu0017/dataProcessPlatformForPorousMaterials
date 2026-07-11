import os
import re
import csv
import tkinter as tk
from tkinter import filedialog
# =========================
# 🔥 ENUM CONFIG（用户只改这里）
# =========================
EXTRACT_KEYS_VOLUME = {
    "Widom": {
        "keyword": "[helium] Average Widom Rosenbluth-weight",
        "type": "scalar",
        "scope": "global"
    },
    "BoxVolume (A^3)": {
        "keyword": "Average Volume:",
        "type": "scalar",
        "scope": "global"
    },
    "unit_cell": {
        "keyword": "Unit cell size:",
        "type": "vector3",
        "scope": "global"
    },
    # 🔥 Framework Density：单位选择策略（关键升级）
    "Framework Density (kg/m^3)" : {
        "keyword": "Framework Density",
        "type": "scalar",
        "scope": "global",
        "unit": "kg/m^3"
    }
}
EXTRACT_KEYS0 = {
    "ads_ads (K)": {
        "keyword": "Average Adsorbate-Adsorbate energy",
        "type": "scalar",
        "scope": "global"
    },
    "host_ads (K)": {
        "keyword": "Average Host-Adsorbate energy",
        "type": "scalar",
        "scope": "global"
    },
    "loading (mol/kg)": {
        "keyword": "Average loading absolute [mol/kg framework]",
        "type": "scalar",
        "scope": "global"
    },
    "density adsorption (cm^3 (STP)/cm^3)": {
        "keyword": "Average loading absolute [cm^3 (STP)/cm^3 framework]",
        "type": "scalar",
        "scope": "global"
    },
    "temperature (K)": {
        "keyword": "External temperature",
        "type": "scalar",
        "scope": "global"
    },
    "cycles": {
        "keyword": "Number of cycles",
        "type": "scalar",
        "scope": "global"
    },
    # 🔥 Qst：单位选择策略（关键升级）
    "Qst": {
        "keyword": "Enthalpy of adsorption",
        "type": "scalar",
        "scope": "global",
        "unit": "kJ/mol"
    },
    "MolFraction": {
        "keyword": "MolFraction",
        "type": "scalar",
        "scope": "molecule"
    },
    "pressure (Pa)": {
        "keyword": "Partial pressure",
        "type": "scalar",
        "scope": "molecule"
    },
    "unit_cell": {
        "keyword": "Unit cell size:",
        "type": "vector3",
        "scope": "global"
    },
    "cell_repeat": {
        "keyword": "Number of unitcells",
        "type": "vector3_int",
        "scope": "global"
    },
    "cell_width": {
        "keyword": "perpendicular cell widths:",
        "type": "vector3",
        "scope": "global"
    }
}
# =========================
# 🔥 BLOCK SPLITTER（核心升级🔥）
# =========================
def split_blocks(text):
    lines = text.splitlines()
    blocks = {}
    current_name = None
    current_block = []
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if i + 1 < len(lines):
            next_line = lines[i + 1]
            if re.match(r"=+", next_line):
                if current_name:
                    blocks[current_name] = "\n".join(current_block)
                current_name = line
                current_block = []
                i += 2
                continue
        if current_name:
            current_block.append(lines[i])
        i += 1
    if current_name:
        blocks[current_name] = "\n".join(current_block)
    return blocks
# =========================
# 🔧 COMPONENT SPLIT
# =========================
def split_components(text):
    pattern = r"(Component\s+\d+\s+\[.*?\].*?)(?=Component\s+\d+\s+\[|$)"
    return re.findall(pattern, text, re.S)

def split_molecules(text):
    """
    提取 MoleculeDefinitions 层（MolFraction / Pressure 所在层）
    """
    pattern = r"(MoleculeDefinitions:.*?)(?=MoleculeDefinitions:|$)"
    return re.findall(pattern, text, re.S)
# =========================
# 🔧 VECTOR
# =========================
def extract_vector3(text, keyword):
    pattern = rf"{re.escape(keyword)}\s*([\d\.\-E+]+)\s+([\d\.\-E+]+)\s+([\d\.\-E+]+)"
    m = re.search(pattern, text)
    return [float(m.group(i)) for i in range(1, 4)] if m else [None]*3


def extract_unitcells(text):
    a = re.search(r"Number of unitcells \[a\]:\s+(\d+)", text)
    b = re.search(r"Number of unitcells \[b\]:\s+(\d+)", text)
    c = re.search(r"Number of unitcells \[c\]:\s+(\d+)", text)
    return [
        int(a.group(1)) if a else None,
        int(b.group(1)) if b else None,
        int(c.group(1)) if c else None
    ]
# =========================
# 🔥 SMART SCALAR PARSER（核心🔥）
# =========================
def extract_scalar_from_block(block, meta=None):
    lines = block.splitlines()
    number = r"([-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?)"
    prefer_unit = meta.get("unit") if meta else None
    # =========================
    # 🔥 1. unit priority（最高优先级）
    # =========================
    if prefer_unit:
        for l in lines:
            if prefer_unit in l:
                m = re.search(number, l)
                if m:
                    return float(m.group(1))
    # =========================
    # 🔥 2. Average（但要防 unit override）
    # =========================
    avg_k = None
    avg_kj = None
    for l in lines:
        if l.strip().startswith("Average"):
            m = re.search(number, l)
            if m:
                avg_k = float(m.group(1))
        if "[KJ/MOL]" in l:
            m = re.search(number, l)
            if m:
                avg_kj = float(m.group(1))
    # 👉 如果用户要 kJ/mol
    if prefer_unit == "kJ/mol" and avg_kj is not None:
        return avg_kj
    # 👉 默认返回 K
    if avg_k is not None:
        return avg_k
    # =========================
    # fallback ±
    # =========================
    m = re.search(rf"{number}\s*\+/-\s*{number}", block)
    if m:
        return float(m.group(1))
    # fallback
    m = re.search(number, block)
    return float(m.group(1)) if m else None

def extract_molecule_scalar(block, keyword):
    lines = block.splitlines()
    number = r"([-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?)"
    for i, l in enumerate(lines):
        # case 1: MolFraction: 1.0
        if keyword in l:
            m = re.search(number, l)
            if m:
                return float(m.group(1))
            # case 2: MolFraction \n 1.0
            if i + 1 < len(lines):
                m2 = re.search(number, lines[i + 1])
                if m2:
                    return float(m2.group(1))
        # case 3: Partial pressure: 100000
        if keyword in l and ":" in l:
            m = re.search(number, l.split(":")[-1])
            if m:
                return float(m.group(1))
    return None

def extract_component_scalar(block, keyword):
    lines = block.splitlines()
    number = r"([-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?)"
    for i, l in enumerate(lines):
        # case 1: MolFraction 0.25
        if keyword in l:
            m = re.search(number, l)
            if m:
                return float(m.group(1))
            # case 2: MolFraction \n 0.25
            if i + 1 < len(lines):
                m2 = re.search(number, lines[i + 1])
                if m2:
                    return float(m2.group(1))
        # case 3: key : value
        if keyword in l and ":" in l:
            m = re.search(number, l.split(":")[-1])
            if m:
                return float(m.group(1))
    return None
# =========================
# 🔥 PARSER CORE（重构🔥）
# =========================
def parse_raspa_output(filepath, keys_dict):
    written_keys = set()
    with open(filepath, 'r') as f:
        text = f.read()
    results = {}
    blocks = split_blocks(text)
    # =====================
    # GLOBAL
    # =====================
    for user_key, cfg in keys_dict.items():
        if cfg["scope"] != "global":
            continue
        if cfg["type"] == "scalar":
            value = None
            for block_name, block_text in blocks.items():
                if cfg["keyword"] in block_name:
                    value = extract_scalar_from_block(
                        block_text,
                        meta=cfg
                    )
                    break
            if value is None:
                m = re.search(rf"{re.escape(cfg['keyword'])}.*?([-+]?\d*\.?\d+)", text, re.S)
                value = float(m.group(1)) if m else None
            key = user_key  # ⭐ global 不加后缀
            if key not in written_keys:
                results[key] = value
                written_keys.add(key)
        elif cfg["type"] == "vector3":
            vec = extract_vector3(text, cfg["keyword"])
            base = user_key
            for i, suf in enumerate(["x", "y", "z"]):
                k = f"{base}_{suf}"
                if k not in written_keys:
                    results[k] = vec[i]
                    written_keys.add(k)
        elif cfg["type"] == "vector3_int":
            vec = extract_unitcells(text)
            base = user_key
            for i, suf in enumerate(["a", "b", "c"]):
                k = f"{base}_{suf}"
                if k not in written_keys:
                    results[k] = vec[i]
                    written_keys.add(k)
    # =====================
    # COMPONENT
    # =====================
    comp_blocks = split_components(text)
    for block in comp_blocks:
        name_match = re.search(r"Component\s+\d+\s+\[(.*?)\]", block)
        gas = name_match.group(1) if name_match else "unknown"
        for user_key, cfg in keys_dict.items():
            if cfg["scope"] != "component":
                continue
            if cfg["type"] == "scalar":
                value = extract_component_scalar(block, cfg["keyword"])
                key = f"{user_key}_{gas}"  # ⭐ component统一后缀
                if key not in written_keys:
                    results[key] = value
                    written_keys.add(key)
    # =====================
    # MOLECULE
    # =====================
    mol_blocks = split_molecules(text)
    for block in mol_blocks:
        for user_key, cfg in keys_dict.items():
            if cfg["scope"] != "molecule":
                continue
            if cfg["type"] == "scalar":
                value = extract_molecule_scalar(block, cfg["keyword"])
                key = user_key  # ⭐ molecule 不再重复生成
                if value is not None and key not in written_keys:
                    results[key] = value
                    written_keys.add(key)
    return results
# =========================
# 🔧 FILE SCAN
# =========================
def find_all_outputs(base_dir, EXTRACT_KEYS):
    results = []
    for root, dirs, files in os.walk(base_dir):
        for file in files:
            if file.startswith("output_") and file.endswith(".data"):
                full_path = os.path.join(root, file)
                data = parse_raspa_output(full_path, EXTRACT_KEYS)
                ads_ads = data.get("ads_ads (K)")
                host_ads = data.get("host_ads (K)")
                if ads_ads is not None and host_ads not in [None, 0]:
                    R = abs(ads_ads) / abs(host_ads)
                else:
                    R = None
                data.update({
                    "R = |Ads/Host|": R,
                    "file": file,
                    "path": full_path
                })
                results.append(data)
    return results
# =========================
# 🔧 CSV
# =========================
def build_csv_keys(results, extract_keys):
    all_keys = set()
    for r in results:
        all_keys.update(r.keys())
    priority_keys = [
        "ads_ads (K)",
        "host_ads (K)",
        "R = |Ads/Host|",
        "loading (mol/kg)"
    ]
    # =====================
    # global keys（严格按 scope）
    # =====================
    global_keys = [
        k for k, v in extract_keys.items()
        if v["scope"] == "global"
    ]
    # =====================
    # molecule keys
    # =====================
    molecule_keys = [
        k for k, v in extract_keys.items()
        if v["scope"] == "molecule"
    ]
    # =====================
    # component keys（从结果反推）
    # =====================
    component_keys = sorted([
        k for k in all_keys
        if "_" in k and any(k.startswith(gk + "_") for gk in extract_keys)
    ])
    # =====================
    # structure keys
    # =====================
    structure_keys = sorted([
        k for k in all_keys
        if k.startswith("unit_cell") or k.startswith("cell_")
    ])
    keys = (
        [k for k in priority_keys if k in all_keys] +
        [k for k in global_keys if k in all_keys] +
        component_keys +
        structure_keys +
        [k for k in molecule_keys if k in all_keys] +
        ["file", "path"]
    )
    # ⭐ 强制去重（关键）
    keys = list(dict.fromkeys(keys))
    return keys

import os
import pandas as pd


def save_to_xlsx(results, base_dir, EXTRACT_KEYS):
    if not results:
        print("No results found.")
        return
    output_path = os.path.join(base_dir, "summary_data.xlsx")
    keys = build_csv_keys(results, EXTRACT_KEYS)
    # 数据部分
    df = pd.DataFrame(results)
    # 按 keys 排列列顺序
    df = df.reindex(columns=keys)
    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        # 第一行：文件名
        pd.DataFrame(
            [["FILE_LIST"] + [r["file"] for r in results]]
        ).to_excel( writer, sheet_name="Summary", index=False, header=False, startrow=0 )
        # 第二行：字段名
        pd.DataFrame([keys]).to_excel( writer, sheet_name="Summary", index=False, header=False, startrow=1 )
        # 第三行：关键词说明
        pd.DataFrame([[
            EXTRACT_KEYS[k]["keyword"] if k in EXTRACT_KEYS else ""
            for k in keys
        ]]).to_excel( writer, sheet_name="Summary", index=False, header=False, startrow=2 )
        # 数据部分
        df.to_excel( writer, sheet_name="Summary", index=False, startrow=3 )
    print(f"\nResults saved to: {output_path}")

def save_to_csv(results, base_dir, EXTRACT_KEYS):
    if not results:
        print("No results found.")
        return
    output_path = os.path.join(base_dir, "summary.csv")
    keys = build_csv_keys(results, EXTRACT_KEYS)
    with open(output_path, "w", newline='', encoding="utf-8-sig") as f:
        writer = csv.writer(f)
        writer.writerow(["FILE_LIST"] + [r["file"] for r in results])
        writer.writerow(keys)
        writer.writerow([EXTRACT_KEYS[k]["keyword"] if k in EXTRACT_KEYS else "" for k in keys])
        dict_writer = csv.DictWriter(f, fieldnames=keys)
        dict_writer.writerows(results)
    print(f"\nResults saved to: {output_path}")
# =========================
# 🔧 UI
# =========================
def select_folder():
    root = tk.Tk()
    root.withdraw()
    return filedialog.askdirectory(title="请选择包含 output 文件的文件夹")


def print_summary(results, EXTRACT_KEYS):
    for r in results:
        print(f"\nFile: {r.get('file','')}")
        for k in EXTRACT_KEYS.keys():
            print(f"  {k}: {r.get(k)}")
# =========================
# 🚀 MAIN
# =========================
if __name__ == "__main__":
    EXTRACT_KEYS = EXTRACT_KEYS0
    user_path = select_folder()
    if not user_path:
        print("未选择路径，退出")
        exit()
    print(f"选择路径: {user_path}")
    results = find_all_outputs(user_path, EXTRACT_KEYS)
    save_to_xlsx(results, user_path, EXTRACT_KEYS)
    print("\n=== Results ===")
    print_summary(results, EXTRACT_KEYS)
    print(f"\nProcessed {len(results)} files.")