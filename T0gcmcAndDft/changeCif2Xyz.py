import os
import matplotlib.pyplot as plt
import addLabelAndSplitCIF
import cif2xyz
# 1 将大内存的cif文件分割为小cif,名称为part_x.cif
# 2 将每个小cif文件转换为xyz格式,名称为part_x.xyz
# 3 合并为一个xyz

atomNumberForOneFile = 2000 #每个小cif的原子数
bigCIF = "aCarbon-Marks-id006"  #要拆分的大内存cif文件

# 文件位置
if os.name == 'posix':
    output_dir = "/mnt/f/simulation/ASE/Tutorials"
else:
    output_dir = "F:/simulation/ASE/Tutorials"
os.makedirs(output_dir, exist_ok=True)

# addLabelAndSplitCIF.fix_cif_missing_labels(bigCIF, output_dir, atomNumberForOneFile)

# cif2xyz.cif2xyz(bigCIF,output_dir)
cif2xyz.mergexyz(bigCIF,output_dir)

