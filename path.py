"""Path"""
from common import *
import T1fileSource.fileOperation as fl

OUTPUT_DIR = None

def setOutputDir(file_path, folder_name=""):
    global OUTPUT_DIR
    OUTPUT_DIR = fl.check_or_create_folder(file_path, folder_name)
    return OUTPUT_DIR

def getOutputPath(fileName, expand="", expandPos=True, type="png", sub_folder="", parent_level=0):

    if OUTPUT_DIR is None:
        raise RuntimeError("OUTPUT_DIR has not been initialized.")

    output_dir = OUTPUT_DIR

    for _ in range(parent_level):
        output_dir = os.path.dirname(output_dir)

    if sub_folder:
        output_dir = fl.check_or_create_folder(output_dir, sub_folder)

    return fl.get_expanded_name( output_dir, fileName, expand=expand, expandPos=expandPos, type=type )