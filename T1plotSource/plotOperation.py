from common import *
from T1dataProcessSource import dataOperation as dop
LABEL_MAP = {
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
    "E": r"$E$ (kJ mol$^{-1}$)",
    "DeltaG_Jmol": r"$\Delta G$ (J mol$^{-1}$)",
}
def formatLabel(label):
    """
    Convert common labels to publication format.
    """
    return LABEL_MAP.get(label, label)
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
        Generate n colors from a custom palette.
        Parameters
        ----------
        colors : list
            List of color hex strings.
        n : int
            Number of colors.
        """
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
def applyLegend( fig, axes, show=True, position="outside right", frameon=False, ncol=1, fontsize=None, ):
    """
    Apply legend automatically.
    Parameters
    ----------
    fig : matplotlib.figure.Figure
    axes : list[Axes]
        Axes list. Support normal and broken axis.
    show : bool
    position : str
        Legend position:
        - "upper left"
        - "upper right"
        - "center left"
        - "center right"
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
    handles = []
    labels = []
    for ax in axes:
        h, l = ax.get_legend_handles_labels()
        for hh, ll in zip(h, l):
            if ll not in labels:
                handles.append(hh)
                labels.append(ll)
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
        raise ValueError( f"Unknown legend position: {position}" )
    cfg = config[position]
    # ==================================================
    # Choose host axis
    # ==================================================
    if cfg["side"] == "left":
        host = axes[0]
    else:
        host = axes[-1]
    # ==================================================
    # Draw
    # ==================================================
    kwargs = dict( frameon=frameon, ncol=ncol, fontsize=fontsize, )
    if cfg["anchor"] is not None:
        kwargs["bbox_to_anchor"] = cfg["anchor"]
    host.legend( handles, labels, loc=cfg["loc"], **kwargs, )
def saveFigure( fig, savepath=None, dpi=600, transparent=False, ):
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
        setMajorTicks(ax_i,axis="x", major_step=200,)
        setMinorTicks(ax_i,axis="x", minor_step=100,)
        if not show_yticks:
            ax_i.set_yticks([])
    # --------------------------------------------------
    # Axis Label
    # --------------------------------------------------
    if len(axes) == 1:
        axes[0].set_xlabel( formatLabel(x if xlabel is None else xlabel) )
        axes[0].set_ylabel( formatLabel(y if ylabel is None else ylabel) )
    else:
        fig.supxlabel( formatLabel(x if xlabel is None else xlabel) )
        axes[0].set_ylabel( formatLabel(y if ylabel is None else ylabel) )
    # --------------------------------------------------
    # X Limit
    # --------------------------------------------------
    xlims = splitXLim( xmin, xmax, xbreak=xbreak, reverse_x=reverse_x, )
    for ax_i, xlim in zip(axes, xlims):
        ax_i.set_xlim(*xlim)
    # ==================================================
    # Legend
    # ==================================================
    applyLegend( fig, axes, show=legend, )
    # ==================================================
    # Save
    # ==================================================
    if savepath is not None:
        saveFigure( fig, savepath=savepath, dpi=600)
    return axes
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
def scatterSphere(
    ax,
    x,
    y,
    color="#d62728",
    edgecolor="black",
    s=60,
    linewidth=0.6,
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
    color : str
        Marker color.
    edgecolor : str
        Marker edge color.
    s : float
        Marker size.
    linewidth : float
        Edge line width.
    highlight : bool
        Whether to draw highlight.
    highlight_size : float
        Fraction of marker size used for highlight.
    highlight_alpha : float
        Transparency of highlight.
    """
    # Base marker
    ax.scatter( x, y, s=s, c=color, edgecolors=edgecolor, linewidths=linewidth, zorder=2, )
    if highlight:
        x = np.asarray(x)
        y = np.asarray(y)
        # 根据数据范围自动计算偏移量
        dx = (x.max() - x.min()) * 0.008
        dy = (y.max() - y.min()) * 0.008
        ax.scatter( x - dx, y + dy, s=s * highlight_size, c="white", 
                   alpha=highlight_alpha, linewidths=0, zorder=3, )
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
def _plotCorrelation(ax, result):
    x = np.asarray(result["X"])
    y = np.asarray(result["Y"])
    y_fit = np.asarray(result["Y_fit"])
    idx = np.argsort(x)
    applyAxisStyle(ax)
    ax.text(
        0.05, 0.95,
        f"R² = {result['R2']:.3f}\n"
        f"r = {result['Pearson_r']:.3f}",
        transform=ax.transAxes,
        va="top"
    )
    ax.plot( x[idx], y_fit[idx], color="red", linewidth=2 )
    colors = ColorPalette.get(USING_COLOR, len(x))
    # ax.scatter(x, y,c=colors)
    scatterSphere( ax, x, y, color=colors, s=70 )
    ax.set_xlabel(result["x_name"])
    ax.set_ylabel(result["y_name"])
def plotSingleCorrelation(result):
    fig, ax = plt.subplots()
    _plotCorrelation(ax, result)
    plt.tight_layout()
    plt.show()
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
def plotBar(x, y, xlabel=None, ylabel=None, title=None, rotation=45, ax=None, gradientFlag=True):
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
        fig, ax = plt.subplots(figsize=(6, 4))
    colors = ColorPalette.get(USING_COLOR, len(x))
    # gradientFlag = True
    if gradientFlag:
        bars = ax.bar(x, y, color="white", edgecolor="black", linewidth=0.8)
        # apply gradient
        drawGradientBars( ax, bars, colors, bottom_color="#d9d9d9" )
    else:
        ax.bar(x, y, color=colors)
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if ylabel is not None:
        ax.set_ylabel(ylabel)
    if title is not None:
        ax.set_title(title)
    ax.tick_params(axis="x", rotation=rotation)
    # 论文风格：刻度朝内，仅左、下有刻度
    applyAxisStyle(ax)
    # print("axis:", id(ax))
    # print("colors:", colors.shape)
    # print(colors[0])
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
def _buildCorrelationMatrices( results, order=None, value="R2", pKey="P_value", include=None, exclude=None ):
    """
    Strict correlation matrix builder.
    Key principle:
    ----------------
    order defines matrix structure (NOT results).
    results is only used as a lookup table.
    """
    # =====================================================
    # 1. Define order (ONLY source of truth)
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
    index = {name: i for i, name in enumerate(order)}
    # =====================================================
    # 2. Initialize full matrix (strict grid)
    # =====================================================
    R = np.full((n, n), np.nan)
    VALUE = np.full((n, n), np.nan)
    P = np.full((n, n), np.nan)
    np.fill_diagonal(R, 1.0)
    np.fill_diagonal(VALUE, 1.0)
    np.fill_diagonal(P, 0.0)
    # =====================================================
    # 3. Fill matrix (order × order ONLY)
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
            VALUE[i, j] = item.get(value, np.nan)
            P[i, j] = item.get(pKey, np.nan)
    return order, R, VALUE, P
def _drawLowerTriangle(ax, R, cmap="RdBu_r", norm=None, circleScale=1200, edgeColor="black", edgeWidth=0.5, alpha=1.0):
    """
    Draw the lower triangle of a correlogram.
    Circle
    ------
    Color :
        Pearson correlation coefficient.
    Area :
        |Pearson_r|^2
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axes.
    R : ndarray
        Pearson correlation matrix.
    cmap : str or Colormap
        Colormap.
    norm : matplotlib.colors.Normalize, optional
        Color normalization.
        If None, Normalize(-1, 1) is used.
    circleScale : float
        Scale factor of circle area.
    edgeColor : str
        Circle edge color.
    edgeWidth : float
        Circle edge width.
    alpha : float
        Circle transparency.
    Returns
    -------
    scatter : PathCollection
        Scatter object for colorbar.
    """
    # ---------------------------------------
    # Default normalization
    # ---------------------------------------
    if norm is None:
        norm = Normalize(vmin=-1, vmax=1)
    scatter = None
    n = R.shape[0]
    # ---------------------------------------
    # Draw circles
    # ---------------------------------------
    for i in range(n):
        for j in range(i):
            r = R[i, j]
            if np.isnan(r):
                continue
            scatter = ax.scatter(
                j,
                i,
                s=circleScale * (abs(r) ** 2),   # area
                c=[r],
                cmap=cmap,
                norm=norm,
                edgecolors=edgeColor,
                linewidths=edgeWidth,
                alpha=alpha,
                zorder=3
            )
    return scatter
def _drawUpperTriangle(ax, R, VALUE, P=None, decimals=3, cmap="RdBu_r", 
                       norm=None, showSignificance=True, significanceLevels=None, fontSize=10, fontWeight="normal"):
    """
    Draw the upper triangle of a correlogram.
    Upper triangle
    --------------
    Display numerical statistics (e.g. R²).
    Text color
    ----------
    Pearson correlation coefficient.
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axes.
    R : ndarray
        Pearson correlation matrix.
    VALUE : ndarray
        Matrix displayed in the upper triangle.
    P : ndarray, optional
        P-value matrix.
    decimals : int
        Decimal places.
    cmap : str or Colormap
        Colormap.
    norm : Normalize, optional
        Color normalization.
    showSignificance : bool
        Whether to display significance stars.
    significanceLevels : dict, optional
        Example:
        {
            0.001: "***",
            0.01: "**",
            0.05: "*"
        }
    fontSize : int
    fontWeight : str
    """
    # ---------------------------------------
    # Default settings
    # ---------------------------------------
    if norm is None:
        norm = Normalize(vmin=-1, vmax=1)
    if significanceLevels is None:
        significanceLevels = {
            0.001: "***",
            0.01: "**",
            0.05: "*"
        }
    cmapObj = plt.get_cmap(cmap)
    n = VALUE.shape[0]
    # ---------------------------------------
    # Draw values
    # ---------------------------------------
    for i in range(n):
        for j in range(i + 1, n):
            if np.isnan(VALUE[i, j]):
                continue
            value = VALUE[i, j]
            r = R[i, j]
            color = cmapObj(norm(r))
            text = f"{value:.{decimals}f}"
            # -----------------------------
            # Significance stars
            # -----------------------------
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
                zorder=4
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
def _drawColorbar( fig, ax, scatter, label="Pearson correlation", shrink=0.85, aspect=30, pad=0.02, n=None ):
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
    cbar.set_label(label)
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
        value="R2",
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
    order, R, VALUE, P = _buildCorrelationMatrices( results=results, order=order, value=value,include=include, exclude=exclude)
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
    scatter = _drawLowerTriangle( ax=ax, R=R, cmap=cmap, norm=norm, circleScale=circleScale )
    _drawUpperTriangle( ax=ax, R=R, VALUE=VALUE, P=P, decimals=decimals, cmap=cmap, 
                       norm=norm, showSignificance=showSignificance, significanceLevels=significanceLevels, 
                       fontSize=fontSize )
    # _drawDiagonal( ax=ax, labels=order, mode=diagonalMode, fontSize=fontSize + 1 )
    if showGrid:
        _drawGrid( ax=ax, n=n )
    _beautifyAxes( ax=ax, labels=order, fontSize=fontSize )
    if showColorbar:
        kwargs = dict( label=colorbarLabel )
        kwargs.update(colorbarKwargs)
        _drawColorbar( fig=fig, ax=ax, scatter=scatter, n=n*1.02, **kwargs )
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
###########################
# 往下是标准化绘图
###########################
