from common import *
from T1dataProcessSource import dataOperation as dop
LABEL_MAP = {
    "R2": r"$R^2$",
    "Pearson_r": r"$r$",
    # -------- Spectrum --------
    "Wavenumber": r"Wavenumber (cm$^{-1}$)",
    "Raman shift (cm-1)": r"Raman shift (cm$^{-1}$)",
    "Intensity (a.u.)": "Intensity (a.u.)",
    "Absorbance": "Absorbance",
    "Transmittance": "Transmittance (%)",
    # -------- Adsorption --------
    "Pressure": "Pressure (kPa)",
    "Relative Pressure": r"Relative pressure ($P/P_0$)",
    "Uptake": r"CO$_2$ uptake (mmol g$^{-1}$)",
    # -------- Pore --------
    "Pore Size": "Pore size (nm)",
    "Pore Volume": r"Pore volume (cm$^3$ g$^{-1}$)",
    "SSA": r"Surface area (m$^2$ g$^{-1}$)",
    # -------- Thermodynamics --------
    "Qst": r"$Q_{st}$ (kJ mol$^{-1}$)",
    "E": r"$E$ (J mol$^{-1}$)",
    "bA": r"Affinity of site $A$ (kPa$^{-1}$)",
    "DeltaG_Jmol": r"$\Delta G$ (J mol$^{-1}$)",
    "sele": r"Selectivity (-)",
    "SELE pyIAST": r"$S_{IAST}$",
    "sel_henry": r"$S_{Henry}$",
    "Sample": r"Sample",
    "Vultra (0.7)":   r"$V_{<0.7\,nm}$",
    "Vultra (1)":     r"$V_{<1\,nm}$",
    "Vmic(2)":        r"$V_{<2\,nm}$",
    "Vt":             r"$V_t$",
    "Vme":            r"$V_{meso}$",
    "Vmic(2)/Vt":     r"$V_{<2\,nm}/V_t$",
    "Vultra (1)/Vt":  r"$V_{<1\,nm}/V_t$",
    "BET":            r"$S_{BET}$",
    # "Vultra (0.7)": r"$V_{\mathrm{<0.7nm}}$",
    # "Vultra (1)": r"$V_{\mathrm{<1nm}}$",
    # "Vmic(2)": r"$V_{\mathrm{<2nm}}$",
    # "Vt": r"$V_{\mathrm{t}}$",
    # "Vme": r"$V_{\mathrm{meso}}$",
    # "Vmic(2)/Vt": r"$V_{\mathrm{<2}}/V_{\mathrm{t}}$",
    # "Vultra (1)/Vt": r"$V_{\mathrm{<1}}/V_{\mathrm{t}}$",
    # "BET": r"$S_{\mathrm{BET}}$",

    "Centroid":         r"$D_{\mathrm{c}}$ (nm)",
    "HighLowRatio":     r"HLR ($-$)",
    "CompetitionIndex": r"CI ($-$)",
    "Skewness":         r"Skewness ($-$)",

    "Centroid":         r"Centroid (nm)",
    "HighLowRatio":     r"HighLowRatio (-)",
    "CompetitionIndex": r"CompetitionIndex (-)",
    "Skewness":         r"Skewness (-)",

}
LABEL_MAP_MATRIX = {
    # -------- Spectrum --------
    "High": r"$V_{<0.5\,nm}$",
    "bA-T25": r"$b_{A}$",
    "E (J/mol)": "$E$",
    "ID/IG": r"$I_{D}/I_{G}$",
    "O-EDS": "O (Wt%)",
    "SELE pyIAST": r"$S_{IAST}$",
    "sel_henry": r"$S_{Henry}$",

    "GCMC_High": "$V_{<0.5\,nm}$",
    "GCMC_HighLowRatio": "HLR",
    "GCMC_Centroid": r"$D_{c}$",
    "GCMC_CompetitionIndex": r"$CI$",
    "GCMC_Skewness": "Skewness",

    "DFT_High": "$V_{<0.65\,nm}$",
    "DFT_HighLowRatio": "HLR",
    "DFT_Centroid": r"$D_{c}$",
    "DFT_CompetitionIndex": r"$CI$",
    "DFT_Skewness": "Skewness",
}
def unicodeUnit(text):
    """
    Convert common LaTeX unit exponents to Unicode.
    """
    return (
        text.replace("$^{-1}$", "⁻¹")
            .replace("$^{-2}$", "⁻²")
            .replace("$^{-3}$", "⁻³")
            .replace("$^{2}$", "²")
            .replace("$^{3}$", "³")
    )
def formatLabel(label,Unicode=False):
    """
    Convert common labels to publication format.
    """
    if label is None:
        return ""

    label = str(label)
    label = LABEL_MAP_MATRIX.get(label, label)
    # import matplotlib as mpl
    # print(mpl.rcParams["font.family"])
    # print(mpl.rcParams["font.sans-serif"])
    if Unicode:
        return (
            label.replace("$^{-1}$", "⁻¹")
                .replace("$^{-2}$", "⁻²")
                .replace("$^{-3}$", "⁻³")
                .replace("$^{2}$", "²")
                .replace("$^{3}$", "³")
        )
    else: 
        return label
def applyLabel(default, custom=None):
    """
    Return formatted label.

    Parameters
    ----------
    default : str
        Default variable name.
    custom : str or None
        User-defined label.

    Returns
    -------
    str
    """
    return formatLabel(default if custom is None else custom)
class ColorPalette:
    # -------- predefined palettes --------
    PAPER1 = [
        "#C00001",  # red
        "#CA2320",
        "#D3463F",
        "#DD695E",
        "#E68B7C",
        "#C77A88",  # rose
        "#AD6D8D",  # purple-red
        "#8E5FA8",  # violet
        "#744F9D",  # purple
        "#4F43B8",  # indigo
        "#0112BD",  # blue
    ]
    RED_PURPLE = [
        "#d62728",
        "#ff7f0e",
        "#f1c40f",
        "#2ca02c",
        "#17becf",
        "#1f77b4",
        "#9467bd",
    ]
    BLUE_RED = [
        "#1f77b4",
        "#17becf",
        "#2ca02c",
        "#f1c40f",
        "#ff7f0e",
        "#d62728",
    ]
    VIRIDIS = [
        "#440154",
        "#3b528b",
        "#21918c",
        "#5ec962",
        "#fde725",
    ]
    @staticmethod
    def get(colors, n, name="palette"):
        """
        Generate n colors.
        Parameters
        ----------
        colors : list or str
            Either
            - ColorPalette.PAPER1
            - "PAPER1"
        """
        if isinstance(colors, str):
            colors = getattr(ColorPalette, colors)
        cmap = LinearSegmentedColormap.from_list(name, colors)
        return cmap(np.linspace(0, 1, n))
    @staticmethod
    def cmap(name):
        """
        Return continuous colormap from ColorList.
        """
        colors = getattr(
            ColorPalette,
            name
        )
        return LinearSegmentedColormap.from_list(
            name,
            colors
        )
USING_COLOR = ColorPalette.PAPER1 
def sampleColors(colors, n):
    """
    Sample n colors uniformly from a color list.
    Parameters
    ----------
    colors : list
        Color palette.
    n : int
        Number of required colors.
    Returns
    -------
    list
    """
    if len(colors) <= n:
        return colors
    idx = np.linspace(
        0,
        len(colors) - 1,
        n
    ).round().astype(int)
    return [colors[i] for i in idx]
def splitXLim( xmin, xmax, xbreak, reverse_x=False, ):
    """
    Split x-axis limits for broken axis.
    Returns
    -------
    list
        [(xmin1, xmax1), (xmin2, xmax2)]
    """
    if xbreak is None:
        if reverse_x:
            return [(xmax, xmin)]
        else:
            return [(xmin, xmax)]
    b1, b2 = sorted(xbreak)
    if reverse_x:
        # 左图：高波数
        return [ (xmax, b2), (b1, xmin), ]
    else:
        # 左图：低波数
        return [ (xmin, b1), (b2, xmax), ]
def splitXLim( xmin, xmax, xbreak=None, reverse_x=False, ):
    """
    X limits corresponding to splitSpectrum().
    Returns
    -------
    list
        [(xmin1, xmax1), (xmin2, xmax2)]
    """
    if xbreak is None:
        if reverse_x:
            return [(xmax, xmin)]
        else:
            return [(xmin, xmax)]
    b1, b2 = sorted(xbreak)
    if reverse_x:
        # 左图：高波数
        return [ (xmax, b2), (b1, xmin), ]
    else:
        # 左图：低波数
        return [ (xmin, b1), (b2, xmax), ] 
def splitSpectrum( X, Y, xbreak=None, reverse_x=False, ):
    """
    Split spectrum for broken x-axis.
    Returns
    -------
    list
        [(X1, Y1), (X2, Y2)]
    """
    if xbreak is None:
        return [(X, Y)]
    b1, b2 = sorted(xbreak)
    low = X <= b1
    high = X >= b2
    if reverse_x:
        # 左图：高波数
        return [ (X[high], Y[high]), (X[low], Y[low]), ]
    else:
        # 左图：低波数
        return [ (X[low], Y[low]), (X[high], Y[high]), ]
def drawBreakMark(ax, side="right", size=4.5, linewidth=0.8):
    """
    Draw diagonal break mark with fixed visual size.
    Parameters
    ----------
    ax : matplotlib.axes.Axes
    side : {"left", "right"}
        Which side of axis to draw.
    size : float
        Length of diagonal mark in points.
    """
    # axes 坐标中心
    if side == "right":
        x = 1
    else:
        x = 0
    # 转换 points 到 axes fraction
    fig = ax.figure
    # axes 实际尺寸(pixel)
    bbox = ax.get_window_extent()
    # y方向比例修正
    dx = size / bbox.width
    dy = size / bbox.height
    kwargs = dict( color="black", clip_on=False, linewidth=linewidth, transform=ax.transAxes, )
    if side == "right":
        ax.plot( [x-dx, x+dx], [-dy, +dy], **kwargs, )
        ax.plot( [x-dx, x+dx], [1-dy, 1+dy], **kwargs, )
    else:
        ax.plot( [x-dx, x+dx], [-dy, +dy], **kwargs, )
        ax.plot( [x-dx, x+dx], [1-dy, 1+dy], **kwargs, )
def createAxes( ax=None, figsize=(5, 4), xbreak=None, xrange=None, reverse_x=False, ):
    """
    Create figure and axes.
    Parameters
    ----------
    ax : matplotlib.axes.Axes or None
    figsize : tuple
    xbreak : tuple or None
        (xmin, xmax)
    Returns
    -------
    fig
    axes
    """
    # ------------------------------
    # Existing Axes
    # ------------------------------
    if ax is not None:
        return ax.figure, [ax]
    # ------------------------------
    # Normal
    # ------------------------------
    if xbreak is None:
        fig, ax = plt.subplots(figsize=figsize)
        return fig, [ax]
    # ------------------------------
    # Broken X axis
    # ------------------------------
    if xrange is None:
        width_ratios = [1, 1]
    else:
        # 两段实际数据范围
        if reverse_x:
            width_ratios = [
                abs(xrange[1] - xbreak[1]),
                abs(xbreak[0] - xrange[0]),
            ]
        else:
            width_ratios = [
                abs(xbreak[0] - xrange[0]),
                abs(xrange[1] - xbreak[1]),
            ]
    fig, (ax1, ax2) = plt.subplots(
        1,
        2,
        figsize=figsize,
        sharey=True,
        gridspec_kw={
            "width_ratios": width_ratios,
            "wspace": 0.04,
        },
    )
    # 去掉中间边框
    ax1.spines["right"].set_visible(False)
    ax2.spines["left"].set_visible(False)
    # 右边不要Y轴
    ax2.yaxis.set_visible(False)
    # 去掉中间tick
    ax1.tick_params(right=False)
    ax2.tick_params(left=False)
    # ------------------------------
    # Draw //
    # ------------------------------
    drawBreakMark( ax1, side="right", )
    drawBreakMark( ax2, side="left", )
    return fig, [ax1, ax2]
def drawGradientBars( ax, bars, colors, bottom_color="#d9d9d9", steps=100, ):
    for bar, top_color in zip(bars, colors):
        height = bar.get_height()
        if height == 0:
            continue
        x = bar.get_x()
        width = bar.get_width()
        cmap = LinearSegmentedColormap.from_list( "bar_gradient", [ bottom_color, top_color, ] )
        y_start = 0 if height > 0 else height
        for i in range(steps):
            rect = Rectangle(
                ( x, y_start + abs(height)*i/steps ),
                width,
                abs(height)/steps,
                facecolor=cmap(i/steps),
                edgecolor="none",
                zorder=1,
            )
            ax.add_patch(rect)
        bar.set_facecolor("none")
        bar.set_edgecolor("black")
        bar.set_linewidth(0.8)
        bar.set_zorder(2)
    ax.relim()
    ax.autoscale_view()
def scatterSphere(
    ax,
    x,
    y,
    color="#d62728",
    edgecolor="black",
    s=60,
    linewidth=0.6,
    marker="o",
    label=None,
    zorder=2,
    highlight=True,
    highlight_size=0.28,
    highlight_alpha=0.55,
):
    """
    Draw sphere-like scatter markers.
    Parameters
    ----------
    ax : matplotlib.axes.Axes
    x, y : array-like
        Coordinates.
    color : color-like
        Marker face color.
    edgecolor : color-like
        Marker edge color.
    s : float
        Marker size.
    linewidth : float
        Marker edge width.
    marker : str
        Marker style, e.g. "o", "s", "^", "D".
    label : str, optional
        Legend label.
    zorder : int
        Drawing order.
    highlight : bool
        Whether to draw the specular highlight.
    highlight_size : float
        Relative size of the highlight.
    highlight_alpha : float
        Highlight transparency.
    """
    x = np.asarray(x)
    y = np.asarray(y)
    # ---------- Base marker ----------
    ax.scatter(
        x,
        y,
        s=s,
        c=color,
        edgecolors=edgecolor,
        linewidths=linewidth,
        marker=marker,
        label=label,
        zorder=zorder,
    )
    if not highlight:
        return
    # ---------- Highlight ----------
    if len(x) > 1:
        dx = (x.max() - x.min()) * 0.0035
    else:
        dx = 0
    if len(y) > 1:
        dy = (y.max() - y.min()) * 0.0035
    else:
        dy = 0
    ax.scatter(
        x - dx,
        y + dy,
        s=s * highlight_size,
        c="white",
        alpha=highlight_alpha,
        linewidths=0,
        marker=marker,
        zorder=zorder + 1,
    )
def buildCurveLegend(
    sample_names,
    color_list,
    marker=True,
    line=True,
    fit=False,
    fit_label="Fit",
    linewidth=2,
    linestyle="-",
    fit_linewidth=2,
    fit_linestyle="--",
    markersize=70,
    markeredgecolor="black",
):
    """
    Build legend handles for curve plots.
    Parameters
    ----------
    sample_names : list[str]
    color_list : list
        Colors returned by ColorPalette.get()
    marker : bool
        Whether experimental markers are shown.
    line : bool
        Whether experimental lines are shown.
    fit : bool
        Whether fitted curves are shown.
    fit_label : str
    Returns
    -------
    list
        Legend handles.
    """
    handles = []
    # scatter() 的 s 是面积，Line2D 的 markersize 是直径（pt）
    ms = np.sqrt(markersize)
    for name, color in zip(sample_names, color_list):
        # --------------------------
        # Experimental data
        # --------------------------
        handles.append(
            Line2D(
                [],
                [],
                color=color if line else "none",
                linestyle=linestyle if line else "None",
                linewidth=linewidth,
                marker="o" if marker else None,
                markersize=ms if marker else 0,
                markerfacecolor=color,
                markeredgecolor=markeredgecolor,
                label=name,
            )
        )
        # --------------------------
        # Fitted curve
        # --------------------------
        if fit:
            handles.append(
                Line2D(
                    [],
                    [],
                    color=color,
                    linestyle=fit_linestyle,
                    linewidth=fit_linewidth,
                    label=f"{name} {fit_label}",
                )
            )
    return handles
def removeDuplicateLegend(handles):
    """
    Remove duplicated legend handles according to labels.
    Parameters
    ----------
    handles : list
        Legend handles.
    Returns
    -------
    unique_handles : list
    labels : list[str]
    """
    unique_handles = []
    labels = []
    for h in handles:
        label = h.get_label()
        if label not in labels:
            labels.append(label)
            unique_handles.append(h)
    return unique_handles, labels
def annotateBars(
    ax,
    bars,
    values=None,
    fmt=".2f",
    offset=3,
    fontsize=12,
    color="black",
    ha="center",
    va="bottom",
    rotation=0,
    inside=False,
    position=0.5,
):
    """
    Add value labels to bar charts.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
    bars : BarContainer
    values : array-like, optional
        Values to display. If None, use bar heights.
    fmt : str, default=".2f"
        Number format.
    offset : float, default=3
        Vertical offset in points (used when inside=False).
    fontsize : int, default=12
    color : str, default="black"
    ha : str, default="center"
    va : str, default="bottom"
    rotation : float, default=0
        Text rotation angle.
    inside : bool, default=False
        Whether to place labels inside the bars.
    position : float, default=0.5
        Relative height inside the bar (0=bottom, 1=top).
    """
    if values is None:
        values = [bar.get_height() for bar in bars]

    for bar, value in zip(bars, values):

        if inside:
            xy = (
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() * position,
            )
            xytext = (0, 0)
            textcoords = "offset points"
            va = "center"
        else:
            xy = (
                bar.get_x() + bar.get_width() / 2,
                bar.get_height(),
            )
            xytext = (0, offset)
            textcoords = "offset points"

        text = ax.annotate(
            f"{value:{fmt}}",
            xy=xy,
            xytext=xytext,
            textcoords=textcoords,
            ha=ha,
            va=va,
            rotation=rotation,
            fontsize=fontsize,
            color=color,
        )

        # 标记为柱状图数值标签，避免 applyTextStyle 修改
        text.set_gid("barvalue")
def applyTextStyle( ax=None, fig=None, labelsize=12, ticksize=12, annotationsize=12, ):
    if ax is not None:
        ax.xaxis.label.set_size(labelsize)
        ax.yaxis.label.set_size(labelsize)
        ax.tick_params(axis="both", labelsize=ticksize)
        if annotationsize is not None:
            for text in ax.texts:
                # 保留柱状图数值标签原始字号
                if text.get_gid() == "barvalue":
                    continue
                text.set_fontsize(annotationsize)

    if fig is not None:
        if getattr(fig, "_supxlabel", None) is not None:
            fig._supxlabel.set_fontsize(labelsize)

        if getattr(fig, "_supylabel", None) is not None:
            fig._supylabel.set_fontsize(labelsize)
def applyLegend(
    fig,
    axes,
    show=True,
    handles=None,
    position="outside right",
    frameon=False,
    ncol=1,
    fontsize=None,
):
    """
    Apply legend automatically.
    Parameters
    ----------
    fig : matplotlib.figure.Figure
    axes : list[Axes]
        Axes list. Support normal and broken axis.
    show : bool
    handles : list, optional
        Custom legend handles. If None, legend handles will be collected
        automatically from axes.
    position : str
        Legend position:
        - "upper left"
        - "upper right"
        - "lower left"
        - "lower right"
        - "outside right"
        - "outside left"
    """
    if not show:
        return
    # ==================================================
    # Collect handles
    # ==================================================
    if handles is None:
        handles = []
        for ax in axes:
            h, _ = ax.get_legend_handles_labels()
            handles.extend(h)
    handles, labels = removeDuplicateLegend(handles)
    if not handles:
        return
    # ==================================================
    # Position configuration
    # ==================================================
    config = {
        # inside
        "upper left": {
            "loc": "upper left",
            "anchor": None,
            "side": "left",
        },
        "upper right": {
            "loc": "upper right",
            "anchor": None,
            "side": "right",
        },
        "lower left": {
            "loc": "lower left",
            "anchor": None,
            "side": "left",
        },
        "lower right": {
            "loc": "lower right",
            "anchor": None,
            "side": "right",
        },
        # outside
        "outside right": {
            "loc": "upper left",
            "anchor": (1, 1.03),
            "side": "right",
        },
        "outside left": {
            "loc": "center right",
            "anchor": (-0.02, 0.5),
            "side": "left",
        },
    }
    if position not in config:
        raise ValueError(f"Unknown legend position: {position}")
    cfg = config[position]
    # ==================================================
    # Choose host axis
    # ==================================================
    host = axes[0] if cfg["side"] == "left" else axes[-1]
    # ==================================================
    # Draw
    # ==================================================
    kwargs = dict(
        frameon=frameon,
        ncol=ncol,
        fontsize=fontsize,
    )
    if cfg["anchor"] is not None:
        kwargs["bbox_to_anchor"] = cfg["anchor"]
    host.legend(
        handles,
        labels,
        loc=cfg["loc"],
        **kwargs,
    )

def setMajorTicks( ax, axis="x", major_step=None, direction="in", length=5, width=1, ):
    """
    Set fixed major tick interval and style.
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axis.
    axis : str
        "x", "y", or "both"
    major_step : float or None
        Major tick interval.
    direction : str
        Tick direction: "in", "out", "inout"
    length : float
        Major tick length.
    width : float
        Major tick width.
    """
    if major_step is None:
        return
    locator = MultipleLocator(major_step)
    if axis in ("x", "both"):
        ax.xaxis.set_major_locator(locator)
    if axis in ("y", "both"):
        ax.yaxis.set_major_locator(locator)
    ax.tick_params( which="major", direction=direction, length=length, width=width, )
def setMinorTicks( ax, axis="x", minor_step=None, direction="in", length=3, width=0.8, ):
    """
    Set minor ticks with fixed interval.
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axis.
    axis : str
        Which axis to apply.
        Options:
        - "x"
        - "y"
        - "both"
    minor_step : float or None
        Fixed interval of minor ticks.
        If None, minor ticks are not modified.
    direction : str
        Tick direction.
    length : float
        Minor tick length.
    width : float
        Minor tick width.
    """
    if minor_step is None:
        return
    locator = MultipleLocator(minor_step)
    if axis in ("x", "both"):
        ax.xaxis.set_minor_locator(locator)
    if axis in ("y", "both"):
        ax.yaxis.set_minor_locator(locator)
    ax.tick_params( which="minor", direction=direction, length=length, width=width, )
def applyAxisStyle(ax):
    """Apply default axis style."""
    ax.tick_params(
        axis="both",
        which="both",
        direction="in",
        top=False,
        right=False
    )
def saveFigure( fig, savepath=None, dpi=300, transparent=False, ):
    """
    Save figure.
    """
    if savepath is None:
        return
    fig.savefig( savepath, dpi=dpi, bbox_inches="tight", transparent=transparent, )
def plotSpectrum(
    data,
    ax=None,
    figsize=(5, 4),
    x="Wavenumber",
    y="Transmittance (%)",
    xlabel=None,
    ylabel=None,
    show_yticks=True,
    reverse_x=True,
    offset=0,
    legend=True,
    colors="PAPER1",
    linewidth=1.5,
    linestyle="-",
    xbreak=None,
    textSize = 12,
    savepath=None,
    dpi=600,
    transparent=False,
):
    """
    Plot spectra (FTIR, Raman, XRD, UV-Vis, etc.).
    Parameters
    ----------
    data : dict
        {sample_name: DataFrame}
    xbreak : tuple, optional
        (xmin, xmax), remove the middle region.
    savepath : str or Path, optional
    """
    # ==================================================
    # X range
    # ==================================================
    xmin = np.inf
    xmax = -np.inf
    for sample, df in data.items():
        X = df[x].to_numpy()
        xmin = min(xmin, X.min())
        xmax = max(xmax, X.max())
    # ==================================================
    # Figure
    # ==================================================
    fig, axes = createAxes( ax=ax, figsize=figsize, xbreak=xbreak, xrange=(xmin, xmax),reverse_x=reverse_x)
    # ==================================================
    # Color
    # ==================================================
    if colors is None:
        colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    elif isinstance(colors, str):
        if hasattr(ColorPalette, colors):
            colors = sampleColors( getattr(ColorPalette, colors), len(data) )
        else:
            cmap = plt.get_cmap(colors)
            colors = [ cmap(i) for i in np.linspace(0, 1, len(data)) ]
    # ==================================================
    # Plot
    # ==================================================
    for i, (sample, df) in enumerate(data.items()):
        X = df[x].to_numpy()
        Y = df[y].to_numpy()
        if offset != 0:
            Y = Y + i * offset
        segments = splitSpectrum( X, Y, xbreak=xbreak, reverse_x=reverse_x, )
        for ax_i, (Xs, Ys) in zip(axes, segments):
            ax_i.plot(
                Xs,
                Ys,
                color=colors[i % len(colors)],
                linewidth=linewidth,
                linestyle=linestyle,
                label=sample,
            )
    # ==================================================
    # Axis Style
    # ==================================================
    for ax_i in axes:
        applyAxisStyle(ax_i)
        # applyTextStyle( ax_i,labelsize=textSize, ticksize=textSize, annotationsize=textSize,)
        setMajorTicks(ax_i,axis="x", major_step=200,)
        setMinorTicks(ax_i,axis="x", minor_step=100,)
        if not show_yticks:
            ax_i.set_yticks([])
    # --------------------------------------------------
    # Axis Label
    # --------------------------------------------------
    if len(axes) == 1:
        axes[0].set_xlabel( applyLabel(x, xlabel), )
        axes[0].set_ylabel( applyLabel(y, ylabel) )
    else:
        fig.supxlabel( applyLabel(x, xlabel) , y=-0.015)
        axes[0].set_ylabel( applyLabel(y, ylabel) )
    for ax_i in axes:
        applyTextStyle( ax_i,labelsize=textSize, ticksize=textSize, annotationsize=textSize,)
    applyTextStyle(fig=fig, labelsize=textSize)
    # --------------------------------------------------
    # X Limit
    # --------------------------------------------------
    xlims = splitXLim( xmin, xmax, xbreak=xbreak, reverse_x=reverse_x, )
    for ax_i, xlim in zip(axes, xlims):
        ax_i.set_xlim(*xlim)
    # ==================================================
    # Legend
    # ==================================================
    applyLegend( fig, axes, show=legend, fontsize=textSize)
    
    # ==================================================
    # Save
    # ==================================================
    if savepath is not None:
        saveFigure( fig, savepath=savepath, dpi=dpi)
    return axes



def plotCurve(
    data,
    ax=None,
    figsize=(5, 4),
    x="Pressure",
    y="Uptake",
    xlabel=None,
    ylabel=None,
    marker=True,
    line=True,
    fit=None,
    fit_label="Fit",
    colors="PAPER1",
    linewidth=2,
    linestyle="-",
    fit_linewidth=2,
    fit_linestyle="--",
    markersize=70,
    legend=True,
    legendPosition="outside right",
    savepath=None,
):
    """
    General XY curve plotting.
    Parameters
    ----------
    data : dict
        Experimental data.
        Supported formats:
        {
            "Sample1": DataFrame,
            "Sample2": DataFrame,
        }
        or
        {
            "Sample1": (x_array, y_array),
            "Sample2": (x_array, y_array),
        }
    fit : dict, optional
        Same format as data.
    Returns
    -------
    matplotlib.axes.Axes
    """
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    applyAxisStyle(ax)
    sample_names = list(data.keys())
    color_list = ColorPalette.get(colors, len(sample_names))
    def _extract_xy(curve):
        """Extract x and y from DataFrame or tuple/list."""
        if isinstance(curve, pd.DataFrame):
            xx = np.asarray(curve[x])
            yy = np.asarray(curve[y])
        elif isinstance(curve, (tuple, list)) and len(curve) == 2:
            xx = np.asarray(curve[0])
            yy = np.asarray(curve[1])
        else:
            raise TypeError(
                "Each dataset must be a pandas DataFrame or (x, y)."
            )
        idx = np.argsort(xx)
        return xx[idx], yy[idx]
    for color, name in zip(color_list, sample_names):
        xx, yy = _extract_xy(data[name])
        # Experimental line
        if line:
            ax.plot(
                xx,
                yy,
                color=color,
                linewidth=linewidth,
                linestyle=linestyle,
                label=name if fit is None else None,
                zorder=2,
            )
        # Experimental markers
        if marker:
            scatterSphere(
                ax,
                xx,
                yy,
                color=color,
                s=markersize,
                label=name if not line else None,
                zorder=3,
            )
        # Fit
        if fit is not None and name in fit:
            fx, fy = _extract_xy(fit[name])
            ax.plot(
                fx,
                fy,
                color=color,
                linewidth=fit_linewidth,
                linestyle=fit_linestyle,
                label=f"{name} {fit_label}",
                zorder=1,
            )
    ax.set_xlabel(xlabel or x)
    ax.set_ylabel(ylabel or y)

    if legend:

        handles = buildCurveLegend(
            sample_names=sample_names,
            color_list=color_list,
            marker=marker,
            line=line,
            fit=(fit is not None),
            fit_label=fit_label,
            linewidth=linewidth,
            linestyle=linestyle,
            fit_linewidth=fit_linewidth,
            fit_linestyle=fit_linestyle,
            markersize=markersize,
        )

        applyLegend(
            fig=ax.figure,
            axes=[ax],
            handles=handles,
            position=legendPosition
        )
    if savepath is not None:
        saveFigure( fig=ax.figure, savepath=savepath, dpi=300)

    return ax
def _plotCorrelation(ax, result,  xlabel=None, ylabel=None, text_position=(0.05, 0.95), figsize=(5, 4), 
                     pos="top", xscale="linear", yscale="linear"):
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure
    x = np.asarray(result["X"])
    y = np.asarray(result["Y"])
    y_fit = np.asarray(result["Y_fit"])
    idx = np.argsort(x)
    applyAxisStyle(ax)
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    ax.text(
        text_position[0],
        text_position[1],
        (
        f"{applyLabel('R2')} = {result['R2']:.3f}\n"
        f"{applyLabel('Pearson_r')} = {result['Pearson_r']:.3f}"
        ),
        transform=ax.transAxes,
        va=pos
    )
    ax.plot( x[idx], y_fit[idx], color="red", linewidth=2 )
    colors = ColorPalette.get(USING_COLOR, len(x))
    scatterSphere( ax, x, y, color=colors, s=70 )
    ax.set_xlabel( applyLabel(result["x_name"], xlabel) )
    ax.set_ylabel( applyLabel(result["y_name"], ylabel) )
    applyTextStyle(ax,labelsize=16, ticksize=14, annotationsize=16)

def plotSingleCorrelation(result,  xlabel=None, ylabel=None, text_position=(0.05, 0.95),figsize=(5, 4),
                           pos="top", xscale="linear", yscale="linear",savepath=None):
    fig, ax = plt.subplots(figsize=figsize)
    _plotCorrelation(ax, result, xlabel, ylabel, text_position, figsize, pos, xscale, yscale)
    plt.tight_layout()
    if savepath is not None:
        saveFigure( fig, savepath=savepath, dpi=600)
    plt.show(block=False)
    
def plotBatchCorrelation(results, topN=5, sort_by="R2"):
    """
    Batch plot top-N correlations.
    Parameters
    ----------
    results : dict
        Output from batchCorrelation().
    topN : int
        Number of figures to plot.
    sort_by : str
        "R2" or "Pearson_r".
    """
    # -----------------------------
    # flatten results
    # -----------------------------
    flat = []
    for xk in results:
        for yk in results[xk]:
            r = results[xk][yk]
            flat.append({
                "x": xk,
                "y": yk,
                "R2": r["R2"],
                "Pearson_r": r["Pearson_r"],
                "result": r
            })
    if len(flat) == 0:
        print("No correlation results.")
        return
    # -----------------------------
    # sort
    # -----------------------------
    flat = sorted(flat, key=lambda d: d[sort_by], reverse=True)
    # 防止 topN 超过总数
    topN = min(topN, len(flat))
    flat = flat[:topN]
    # -----------------------------
    # 自动计算布局
    # -----------------------------
    ncols = math.ceil(math.sqrt(topN))
    nrows = math.ceil(topN / ncols)
    fig, axes = plt.subplots( nrows, ncols, figsize=(5 * ncols, 4 * nrows) )
    # axes统一成一维数组
    axes = np.array(axes).reshape(-1)
    # -----------------------------
    # plot
    # -----------------------------
    for ax, item in zip(axes, flat):
        _plotCorrelation(ax, item["result"])
        ax.set_title( f"{item['x']} vs {item['y']}" )
    # 删除多余子图
    for ax in axes[topN:]:
        fig.delaxes(ax)
    plt.tight_layout()
    plt.show(block=False)
def plotBar( x, y, xlabel=None, ylabel=None, title=None, rotation=45, ax=None, figsize=(6, 4), 
            gradientFlag=True, showValue=False, valueFmt=".3f", valueOffset=3, valueFontsize=14, 
            valueRotation=0, valueInside=False, valuePosition=0.5, savepath=None, ):
    """
    Plot a single bar chart.
    Parameters
    ----------
    x : array-like
        X-axis labels.
    y : array-like
        Bar values.
    xlabel : str, optional
    ylabel : str, optional
    title : str, optional
    rotation : int, default=45
        Rotation angle of x tick labels.
    ax : matplotlib.axes.Axes, optional
        Existing axes. If None, create a new figure.
    Returns
    -------
    ax : matplotlib.axes.Axes
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    colors = ColorPalette.get(USING_COLOR, len(x))
    # gradientFlag = True
    if gradientFlag:
        bars = ax.bar(x, y, color="white", edgecolor="black", linewidth=0.8)
        # apply gradient
        drawGradientBars( ax, bars, colors, bottom_color="#d9d9d9" )
    else:
        ax.bar(x, y, color=colors)
    if xlabel is not None:
        ax.set_xlabel( applyLabel(x, xlabel) )
    if ylabel is not None:
        ax.set_ylabel( applyLabel(y, ylabel) )
    if title is not None:
        ax.set_title(title)
    ax.tick_params(axis="x", rotation=rotation)
    # 论文风格：刻度朝内，仅左、下有刻度
    if showValue:
        annotateBars(
            ax,
            bars,
            fmt=valueFmt,
            offset=valueOffset,
            fontsize=valueFontsize,
            color="white",
            rotation=valueRotation,
            inside=valueInside,
            position=valuePosition,
        )
    applyAxisStyle(ax)
    applyTextStyle( ax=ax, labelsize=16, ticksize=14, annotationsize=16 )
    if savepath is not None:
        saveFigure( fig, savepath=savepath, dpi=900)
    return ax
def plotBarByMetrics(metrics, columns=None, ylabel=None, title=None, ncols=1):
    """
    Plot one or multiple metric columns.
    Parameters
    ----------
    metrics : pandas.DataFrame
        Rows are samples and columns are metrics.
    columns : str or list[str], optional
        Columns to plot.
        Default is all columns.
    ylabel : str, optional
    title : str, optional
    ncols : int, default=1
        Number of subplot columns.
    """
    if columns is None:
        columns = list(metrics.columns)
    if isinstance(columns, str):
        columns = [columns]
    # 保证样品顺序
    metrics = dop.naturalSortData(metrics)
    n = len(columns)
    nrows = math.ceil(n / ncols)
    fig, axes = plt.subplots( nrows, ncols, figsize=(6 * ncols, 3.5 * nrows), squeeze=False )
    axes = axes.flatten()
    for ax, column in zip(axes, columns):
        plotBar( metrics.index, metrics[column], ylabel=ylabel, title=column, ax=ax, gradientFlag=True)
    # 删除多余子图
    for ax in axes[n:]:
        fig.delaxes(ax)
    if title is not None:
        fig.suptitle(title, fontsize=14)
    fig.tight_layout()
    plt.show(block=False)
def _buildCorrelationMatrices(
    results,
    order=None,
    circleSize="R2",
    circleColor="Pearson_r",
    textValue="R2",
    pKey="P_value",
    include=None,
    exclude=None,
):
    """
    Strict correlation matrix builder.

    Parameters
    ----------
    circleSize : str
        Variable controlling circle area.

    circleColor : str
        Variable controlling circle and text color.

    textValue : str
        Variable displayed in the upper triangle.
    """

    # =====================================================
    # Define order
    # =====================================================
    if order is None:
        order = list(results.keys())

    order = list(order)

    if include is not None:
        include_set = set(include)
        order = [v for v in order if v in include_set]

    if exclude is not None:
        exclude_set = set(exclude)
        order = [v for v in order if v not in exclude_set]

    n = len(order)

    # =====================================================
    # Initialize matrices
    # =====================================================
    R = np.full((n, n), np.nan)
    SIZE = np.full((n, n), np.nan)
    COLOR = np.full((n, n), np.nan)
    TEXT = np.full((n, n), np.nan)
    P = np.full((n, n), np.nan)

    np.fill_diagonal(R, 1.0)
    np.fill_diagonal(SIZE, 1.0)
    np.fill_diagonal(COLOR, 1.0)
    np.fill_diagonal(TEXT, 1.0)
    np.fill_diagonal(P, 0.0)

    # =====================================================
    # Fill matrices
    # =====================================================
    for i, x in enumerate(order):

        row = results.get(x)
        if row is None:
            continue

        for j, y in enumerate(order):

            item = row.get(y)
            if item is None:
                continue

            R[i, j] = item.get("Pearson_r", np.nan)
            SIZE[i, j] = item.get(circleSize, np.nan)
            COLOR[i, j] = item.get(circleColor, np.nan)
            TEXT[i, j] = item.get(textValue, np.nan)
            P[i, j] = item.get(pKey, np.nan)

    return order, R, SIZE, COLOR, TEXT, P
def _drawLowerTriangle(
    ax,
    SIZE,
    COLOR,
    cmap="RdBu_r",
    norm=None,
    circleScale=1200,
    sizeTransform="none",
    edgeColor="black",
    edgeWidth=0.5,
    alpha=1.0,
):
    """
    Draw the lower triangle of a correlogram.

    Circle
    ------
    Area :
        Controlled by SIZE and sizeTransform.

    Color :
        Controlled by COLOR.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axes.

    SIZE : ndarray
        Matrix controlling circle area.

    COLOR : ndarray
        Matrix controlling circle color.

    cmap : str or Colormap

    norm : matplotlib.colors.Normalize, optional

    circleScale : float

    sizeTransform : {"square", "abs", "none"}
        square : area = |SIZE|²
        abs    : area = |SIZE|
        none   : area = SIZE

    Returns
    -------
    scatter : PathCollection
    """
    if norm is None:
        norm = Normalize(vmin=-1, vmax=1)

    scatter = None
    n = SIZE.shape[0]

    for i in range(n):
        for j in range(i):

            size = SIZE[i, j]
            color = COLOR[i, j]

            if np.isnan(size) or np.isnan(color):
                continue

            if sizeTransform == "square":
                area = circleScale * (abs(size) ** 2)

            elif sizeTransform == "abs":
                area = circleScale * abs(size)

            elif sizeTransform == "none":
                area = circleScale * size

            else:
                raise ValueError(
                    f"Unknown sizeTransform: {sizeTransform}"
                )

            scatter = ax.scatter(
                j,
                i,
                s=area,
                c=[color],
                cmap=cmap,
                norm=norm,
                edgecolors=edgeColor,
                linewidths=edgeWidth,
                alpha=alpha,
                zorder=3,
            )

    return scatter
def _drawUpperTriangle(
    ax,
    COLOR,
    VALUE,
    P=None,
    decimals=3,
    cmap="RdBu_r",
    norm=None,
    showSignificance=True,
    significanceLevels=None,
    fontSize=10,
    fontWeight="normal",
):
    """
    Draw the upper triangle of a correlogram.

    Upper triangle
    --------------
    Text:
        Controlled by VALUE.

    Text color:
        Controlled by COLOR.

    Significance:
        Controlled by P.
    """

    if norm is None:
        norm = Normalize(vmin=-1, vmax=1)

    if significanceLevels is None:
        significanceLevels = {
            0.001: "***",
            0.01: "**",
            0.05: "*",
        }

    cmapObj = plt.get_cmap(cmap)

    n = VALUE.shape[0]

    for i in range(n):
        for j in range(i + 1, n):

            value = VALUE[i, j]
            colorValue = COLOR[i, j]

            if np.isnan(value):
                continue

            color = cmapObj(norm(colorValue))

            text = f"{value:.{decimals}f}"

            if showSignificance and P is not None:

                p = P[i, j]

                if not np.isnan(p):

                    for level, stars in significanceLevels.items():

                        if p <= level:
                            text += stars
                            break

            ax.text(
                j,
                i,
                text,
                ha="center",
                va="center",
                fontsize=fontSize,
                fontweight=fontWeight,
                color=color,
                zorder=4,
            )
def _drawDiagonal(ax, labels, mode="label", fontSize=11, fontWeight="bold"):
    """
    Draw diagonal cells.
    Parameters
    ----------
    ax : matplotlib.axes.Axes
    labels : list[str]
    mode : str
        "label"
            Display variable names.
        "blank"
            Leave diagonal blank.
    fontSize : int
    fontWeight : str
    """
    if mode == "blank":
        return
    n = len(labels)
    for i in range(n):
        ax.text(
            i,
            i,
            labels[i],
            ha="center",
            va="center",
            fontsize=fontSize,
            fontweight=fontWeight,
            zorder=5
        )
def _drawGrid(ax, n, color="lightgray", linewidth=0.8):
    """
    Draw matrix grid.
    Parameters
    ----------
    ax : matplotlib.axes.Axes
    n : int
    color : str
    linewidth : float
    """
    for k in range(n + 1):
        ax.axhline( k - 0.5, color=color, linewidth=linewidth, zorder=1 )
        ax.axvline( k - 0.5, color=color, linewidth=linewidth, zorder=1 )
def _drawColorbar( fig, ax, scatter, label="Pearson correlation", shrink=0.85, aspect=30, pad=0.02, n=None, labelSize=12, tickSize=12,):
    """
    Draw colorbar.
    Parameters
    ----------
    fig : matplotlib.figure.Figure
    ax : matplotlib.axes.Axes
    scatter : PathCollection
    label : str
    shrink : float
    aspect : float
    pad : float
    """
    if scatter is None:
        return None
    # =====================================================
    # Auto adjust based on matrix size
    # =====================================================
    if n is not None:
        # n 越大 → colorbar 越短
        shrink = max(0.5, min(0.95, 1.2 - 0.03 * n))
        # n 越大 → colorbar 越“细长”
        aspect = max(15, min(50, 35 - 0.3 * n))
    cbar = fig.colorbar( scatter, ax=ax, shrink=shrink, aspect=aspect, pad=pad )
    # Label
    cbar.set_label(label, fontsize=labelSize)

    # Tick labels
    cbar.ax.tick_params(labelsize=tickSize)
    return cbar
def _beautifyAxes(ax, labels, fontSize=10):
    """
    Beautify matrix axes.
    Parameters
    ----------
    ax : matplotlib.axes.Axes
    labels : list[str]
    fontSize : int
    """
    labels = [applyLabel(label) for label in labels]
    n = len(labels)
    ax.set_xlim(-0.5, n - 0.5)
    ax.set_ylim(n - 0.5, -0.5)
    ax.set_aspect("equal")
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels( labels, rotation=45, ha="right", fontsize=fontSize )
    ax.set_yticklabels( labels, fontsize=fontSize )
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
def _autoCircleScale(circleScale, n):
    """
    Automatically adjust circle scale according to
    the number of variables.
    """
    return circleScale * (10 / max(n, 10))
def _autoFontSize(n):
    return max(6, int(14 - 0.2*n))
def _autoTextColor(color):
    r,g,b = to_rgb(color)
    luminance = ( 0.299*r + 0.587*g + 0.114*b )
    return "white" if luminance < 0.5 else "black"
def _filterResults(results, include=None, exclude=None):
    """
    Filter raw results before building correlation matrices.
    This ensures all downstream steps are consistent.
    """
    if include is None and exclude is None:
        return results
    # 复制，避免污染原数据
    results_out = {}
    for k, v in results.items():
        # 如果 results 是 dict-of-dict / matrix结构
        # 假设 v 是 dict: {var: values} 或类似结构
        if isinstance(v, dict):
            v_out = {}
            for key, val in v.items():
                if include is not None and key not in include:
                    continue
                if exclude is not None and key in exclude:
                    continue
                v_out[key] = val
            results_out[k] = v_out
        else:
            # 如果是矩阵/数组类结构，直接保留（或你可扩展）
            results_out[k] = v
    return results_out
def plotCorrelogram(
        results,
        order=None, # ---------- 指定变量顺序 ----------
        include= None, # ---------- 排除或只留下指定变量 ----------
        exclude= None,
        # ---------- Visual Mapping ----------
        circleSize="R2",
        circleColor="Pearson_r",
        textValue="Pearson_r",
        ax=None,
        figsize=None,
        # ---------- Color ----------
        cmap="RdBu_r",
        # cmap="RdBu_r",
        norm=None,
        vmin=-1,
        vmax=1,
        # ---------- Circle ----------
        circleScale=1200,
        # ---------- Text ----------
        decimals=3,
        # ---------- Significance ----------
        showSignificance=True,
        significanceLevels=None,
        # ---------- Grid ----------
        showGrid=True,
        # ---------- Colorbar ----------
        showColorbar=True,
        colorbarLabel="Pearson correlation",
        colorbarKwargs=None,
        # ---------- Diagonal ----------
        diagonalMode="label",
        # ---------- Save ----------
        dpi=300,
        savePath=None
):
    """
    Plot a publication-quality correlogram.
    Returns
    -------
    fig : Figure
    ax : Axes
    """
    # results = _filterResults(results, include=include, exclude=exclude)
    # =====================================================
    # Build matrices
    # =====================================================
    order, R, SIZE, COLOR, TEXT, P = _buildCorrelationMatrices( results=results, order=order, circleSize=circleSize, 
                                                   circleColor=circleColor, textValue=textValue, include=include, exclude=exclude)
    n = len(order)
    # =====================================================
    # Automatic parameters
    # =====================================================
    if figsize is None:
        size = max(6, min(14, 0.6 * n + 2))
        figsize = (size, size)
    fontSize = _autoFontSize(n)
    circleScale = _autoCircleScale(circleScale, n)
    if norm is None:
        norm = Normalize(vmin=vmin, vmax=vmax)
    if isinstance(cmap, str):
        if hasattr(ColorPalette, cmap):
            cmap = ColorPalette.cmap(cmap)
        else:
            cmap = plt.get_cmap(cmap)
    if colorbarKwargs is None:
        colorbarKwargs = {}
    # =====================================================
    # Figure
    # =====================================================
    createdFigure = False
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
        createdFigure = True
    else:
        fig = ax.figure
    # =====================================================
    # Draw
    # =====================================================
    if circleSize in {"R2", "Adjusted_R2"}:
        sizeTransform = "none"
    else:
        sizeTransform = "square"
    scatter = _drawLowerTriangle( ax=ax, SIZE=SIZE, COLOR=COLOR, cmap=cmap, norm=norm, circleScale=circleScale, 
                                 sizeTransform=sizeTransform, )

    _drawUpperTriangle( ax=ax, COLOR=COLOR, VALUE=TEXT, P=P, decimals=decimals, cmap=cmap, norm=norm, 
                       showSignificance=showSignificance, significanceLevels=significanceLevels, fontSize=fontSize, )
    # _drawDiagonal( ax=ax, labels=order, mode=diagonalMode, fontSize=fontSize + 1 )
    if showGrid:
        _drawGrid( ax=ax, n=n )
    _beautifyAxes( ax=ax, labels=order, fontSize=fontSize + 2 )
    if showColorbar:
        kwargs = dict( label=colorbarLabel )
        kwargs.update(colorbarKwargs)
        _drawColorbar( fig=fig, ax=ax, scatter=scatter, n=n*1.02, labelSize=fontSize + 4, tickSize=fontSize + 2, **kwargs )
    # =====================================================
    # Save
    # =====================================================
    if savePath is not None:
        fig.savefig( savePath, dpi=dpi, bbox_inches="tight" )
    # =====================================================
    # Show
    # =====================================================
    if createdFigure:
        plt.tight_layout()
        plt.show(block=False)
    return fig, ax

