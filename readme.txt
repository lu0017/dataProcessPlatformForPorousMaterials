"""
##每一个子平台固定开头，用于找到依赖
import os
import sys
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from common import *

其它方法
方法一：需要在根目录luDataProcess下配置json文件，定义为工作文件夹，否则下层子文件无法启动 "cwd": "${workspaceFolder}",
json文件内容：
{
    "version": "0.2.0",
    "configurations": [
        {
            "name": "Python: Current File",
            "type": "debugpy",
            "request": "launch",
            "program": "${file}",
            "cwd": "${workspaceFolder}",
            "console": "integratedTerminal"
        }
    ]
}
方法二：ctrl+ shift +P : Python > Terminal: Execute In File Dir, 把它取消勾选（False）
方法三：直接在 settings.json 加：
{
    "python.terminal.executeInFileDir": false
}
"""

文件夹命名
T0：工程，应用，论文数据一级，如GCMC， DSL， TSA
T1：代码内部支撑文件