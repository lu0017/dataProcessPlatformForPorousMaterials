import os
import glob
import re
from ase.io import read, write

def cif2xyz(file, filepath):

    # 找出所有 part_*.cif 文件
    file_list = glob.glob(os.path.join(filepath, "part_*.cif"))

    # 提取数字排序函数
    def extract_num(f):
        match = re.search(r"part_(\d+)", os.path.basename(f))
        return int(match.group(1)) if match else -1
    
    # 找出所有 part_*.cif 文件并排序
    file_list = glob.glob(os.path.join(filepath, "part_*.cif"))
    file_list_sorted = sorted(file_list, key=extract_num)

    # 将 .cif 转换为 .xyz，命名为 part_001.xyz 等
    for i, file_path in enumerate(file_list_sorted, 1):
        print(f"read: {file_path}")
        
        try:
            atoms = read(file_path)
            
            # 可选：保存为 .xyz 格式（用于OVITO或后处理）
            xyz_name = os.path.join(filepath, f"part_{i:03d}.xyz")
            write(xyz_name, atoms)
            
        except Exception as e:
            print(f"fail to read {file_path}: {e}")

def mergexyz(file, filepath,splitFlag=False):
    # 提取数字排序函数
    def extract_num(f):
        match = re.search(r"part_(\d+)", os.path.basename(f))
        return int(match.group(1)) if match else -1
    
    # 合并所有 .xyz 为 merged.xyz（多帧）
    xyz_files = glob.glob(os.path.join(filepath, "part_*.xyz"))
    xyz_sorted = sorted(xyz_files, key=extract_num)
    all_atom_lines = []
    if splitFlag:
        merged_path = os.path.join(filepath, file + ".xyz")
        with open(merged_path, "w") as outfile:
            for file_path in xyz_sorted:
                with open(file_path, "r") as infile:
                    lines = infile.readlines()
                    if len(lines) < 2:
                        continue
                    outfile.write(lines[0])  # 原子数
                    outfile.write(f"From: {os.path.basename(file_path)}\n")  # 注释行

                    for line in lines[2:]:
                        parts = line.strip().split()
                        if len(parts) >= 4:
                            # 仅保留前4列（元素 + x y z）
                            outfile.write(f"{parts[0]} {parts[1]} {parts[2]} {parts[3]}\n")
    else:
        for file_path in xyz_sorted:
            with open(file_path, "r") as infile:
                lines = infile.readlines()
                if len(lines) < 2:
                    continue
                # 处理每行原子数据，去掉多余列（只保留元素和三坐标）
                for line in lines[2:]:
                    parts = line.strip().split()
                    if len(parts) >= 4:
                        all_atom_lines.append(f"{parts[0]} {parts[1]} {parts[2]} {parts[3]}\n")

        merged_path = os.path.join(filepath, file + ".xyz")
        with open(merged_path, "w") as outfile:
            # 写入所有原子数（总数）
            outfile.write(f"{len(all_atom_lines)}\n")
            # 写注释行，可以自定义
            outfile.write(f"Merged from {len(xyz_sorted)} files\n")
            # 写入所有原子数据
            outfile.writelines(all_atom_lines)
    print(f"merge finish：{merged_path}")