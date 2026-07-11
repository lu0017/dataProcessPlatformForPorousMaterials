import os
import sys
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from common import *
import T1fileSource.fileOperation as fl

def add_xy_chart_to_excel(file_path, sheet_name="0C"):
    wb = load_workbook(file_path)
    ws = wb[sheet_name] if sheet_name else wb.active

    chart = ScatterChart()
    chart.title = "Adsorption Isotherms"
    chart.style = 2
    chart.x_axis.title = "Pressure (kPa)"
    chart.y_axis.title = "Uptake (mmol/g)"

    # ==== 自动找到数据起点 ====
    min_row, min_col = None, None
    max_row, max_col = 0, 0
    for row in ws.iter_rows():
        for cell in row:
            if cell.value not in (None, ""):
                if min_row is None or cell.row < min_row:
                    min_row = cell.row
                if min_col is None or cell.column < min_col:
                    min_col = cell.column
                if cell.row > max_row:
                    max_row = cell.row
                if cell.column > max_col:
                    max_col = cell.column

    if min_row is None:
        print("⚠️ 表格为空，未找到数据。")
        return

    print(f"🔍 数据范围: 行 {min_row}~{max_row}, 列 {min_col}~{max_col}")

    chart = ScatterChart()
    chart.title = "Adsorption Isotherms"
    chart.style = 2
    chart.x_axis.title = "Pressure (kPa)"
    chart.y_axis.title = "Uptake (mmol/g)"
    chart.smooth = True   # ✅ 平滑曲线（对应 Excel 的 "Smoothed lines"）
    chart.legend.position = 'r'  # 图例在右侧
    chart.legend.layout = None       # 保持默认布局
    chart.legend.overlay = False     # 确保图例在图表外部显示

    chart.x_axis.axId = 100
    chart.y_axis.axId = 101
    chart.x_axis.crossAx = chart.y_axis.axId
    chart.y_axis.crossAx = chart.x_axis.axId
    chart.x_axis.majorTickMark = "in"
    chart.y_axis.majorTickMark = "in"
    chart.x_axis.crosses = "min"        # X轴在底部交叉
    chart.y_axis.crosses = "min"        # Y轴在左侧交叉
    chart.x_axis.tickLblPos = "nextTo"  # 显示X轴刻度值
    chart.y_axis.tickLblPos = "nextTo"  # 显示Y轴刻度值
    chart.x_axis.majorGridlines = None
    chart.y_axis.majorGridlines = None
        # ✅ 关键：手动创建轴对象（强制启用 Axis）
    chart.x_axis.delete = False
    chart.y_axis.delete = False

    # ==== 每两个列为一组 (X, Y) ====
    for i in range(min_col, max_col + 1, 2):
        name_cell = ws.cell(row=min_row, column=i)
        if name_cell.value is None:
            continue

        # 防止最后一列孤立（没有Y列）
        if i + 1 > max_col:
            break

        # X-Y 数据区域
        xvalues = Reference(ws, min_col=i, min_row=min_row + 2, max_row=max_row)
        yvalues = Reference(ws, min_col=i + 1, min_row=min_row + 2, max_row=max_row)

        series_name = str(name_cell.value).strip()
        series = Series(yvalues, xvalues, title=series_name)
        series.graphicalProperties.line.width = 1500  # 单位EMU，默认12700，大约0.15pt
        chart.series.append(series)

    # ==== 自动放置图表位置 ====
    insert_col = max_col + 2  # 放在数据右侧
    insert_cell = ws.cell(row=min_row, column=insert_col)
    ws.add_chart(chart, insert_cell.coordinate)

    wb.save(file_path)
    print(f"✅ 图表已添加到：{sheet_name or ws.title}")

# 使用示例
def plotFunction(file_path):
    xls = pd.ExcelFile(file_path)
    for sheet in xls.sheet_names:
        if fl.validSheetforSplitFile(sheet): 
            add_xy_chart_to_excel(file_path, sheet_name=sheet)

def main(file_path=None):
    file_path = fl.getFile()
    plotFunction(file_path)

if __name__ == "__main__":
    f = sys.argv[1] if len(sys.argv) > 1 else None
    main(f)