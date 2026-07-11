"""
common.py
整个工程共享的第三方依赖
自定义函数模块在_init_.py中
"""

from __future__ import annotations
# ======================
# Standard Library
# ======================
import os
import sys
import copy
import math
import json
import warnings
import pathlib
import re
import tkinter as tk
from tkinter import filedialog
import pyiast 
from copy import deepcopy

# ======================
# Scientific Computing
# ======================
import numpy as np
import pandas as pd

# ======================
# class / structure / enum
# ======================
from enum import Enum
from dataclasses import dataclass, field
from dataclasses import fields, is_dataclass


# ======================
# Plot
# ======================
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
from matplotlib.colors import (
    ListedColormap,
    to_hex,
    to_rgb,
)
import random
from itertools import cycle

# ======================
# Scipy
# ======================
from scipy.interpolate import interp1d
from scipy.stats import linregress
from scipy.optimize import brentq
from scipy.interpolate import PchipInterpolator
from scipy.interpolate import CubicSpline
from scipy.signal import savgol_filter
from scipy.ndimage import gaussian_filter
from scipy.ndimage import gaussian_filter1d
from scipy.stats import spearmanr

# ======================
# Fitting
# ======================
from lmfit import Parameters, minimize, report_fit
from pybaselines import Baseline
from lmfit.model import ModelResult
from lmfit.models import (
    LorentzianModel,
    GaussianModel,
    VoigtModel,
    PseudoVoigtModel,
)


# ======================
# Excel
# ======================
import glob
import csv
import openpyxl
from openpyxl import load_workbook, Workbook
from openpyxl.chart import ScatterChart, Reference, Series
from openpyxl.chart.axis import ChartLines
from openpyxl.chart.shapes import GraphicalProperties
from openpyxl.drawing.line import LineProperties

from zipfile import BadZipFile
from collections import defaultdict
