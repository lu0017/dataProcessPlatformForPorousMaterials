import os
def split_cif_by_atoms(input_file, filepath,atoms_per_file):
  #  with open(input_file, 'r') as f:
 #       lines = f.readlines()
    lines = input_file
    header_lines = []
    atom_fields = []
    atom_data_lines = []
    in_atom_block = False
    collecting_fields = False
    find_C = False
    for i, line in enumerate(lines):
        stripped = line.strip()
        if stripped.startswith("C"):
            find_C = True
        # Detect start of atom loop
        if stripped == "loop_":
           # next_lines = lines[i+1:i+10]
        #    if any("_atom_site" in l for l in next_lines):
            if i+1 < len(lines) and lines[i+1].strip().startswith("_atom_site"):
                in_atom_block = True
                collecting_fields = True
                header_lines.append("loop_\n")
                continue
        if in_atom_block:
            if collecting_fields:
                if stripped.startswith("_atom_site"):
                    atom_fields.append(line)
                    header_lines.append(line)
                else:
                    collecting_fields = False  # fields done, now actual atom data starts
            if not collecting_fields and stripped:
                atom_data_lines.append(line)
        else:
            header_lines.append(line)
    # Now split atom data lines
    total_atoms = len(atom_data_lines)
    num_chunks = (total_atoms + atoms_per_file - 1) // atoms_per_file
    for idx in range(num_chunks):
        start = idx * atoms_per_file
        end = min(start + atoms_per_file, total_atoms)
        chunk = atom_data_lines[start:end]
        filename0 = f"part_{idx+1:03d}"
        filename = os.path.join(filepath, filename0 + ".cif")
        with open(filename, 'w') as fout:
            fout.writelines(header_lines)
            fout.writelines(chunk)
        print(f"写入 {filename}，包含原子数: {len(chunk)}")
# 使用方式
#split_cif_by_atoms("fixed.cif", atoms_per_file=2000)
def fix_cif_missing_labels(file, filepath, atoms_per_file):
    cif_path = os.path.join(filepath, file + ".cif")
    with open(cif_path, 'r') as f:
        lines = f.readlines()
    new_lines = []
    in_atom_block = False
    atom_index = 1
    for line in lines:
        stripped = line.strip()
        # 标识进入 _atom_site 数据字段区域
        if stripped.startswith("loop_"):
            in_atom_block = False
            new_lines.append(line)
            continue
        if stripped.startswith("_atom_site_label"):
            in_atom_block = True
            new_lines.append(line)
            continue
        if in_atom_block and stripped.startswith("_"):
            new_lines.append(line)
            continue
        # 如果是原子数据行，自动补标签
        if in_atom_block and stripped:
            parts = stripped.split()
            if len(parts) == 4:
                # 补全 label 和元素类型
                symbol = parts[0]
                x, y, z = parts[1:]
                label = f"{symbol}{atom_index}"
                atom_index += 1
                new_line = f"{label} {symbol} {x} {y} {z}\n"
                new_lines.append(new_line)
            else:
                new_lines.append(line)
        else:
            new_lines.append(line)
    split_cif_by_atoms(new_lines,filepath, atoms_per_file)
#    split_cif_file_by_lines(new_lines, lines_per_file=2000)
 #   with open(output_file, 'w') as f:
 #       f.writelines(new_lines)
    print(f"修复完成，已生成 {atom_index - 1} 个原子，结果已保存至：{filepath}")
