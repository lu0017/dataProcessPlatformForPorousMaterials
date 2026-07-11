from common import *

def compute_fit_metrics(fit_x, fit_y, exp_x, exp_y, method='linear', use_interp=True,):
    from sklearn.metrics import r2_score
    """
    计算拟合曲线与实验数据的拟合指标
    包括：
    - 残差 residuals
    - 相对残差 rems
    - R² r2
    
    自动处理拟合数据和实验数据维度不一致，通过插值对齐。
    
    参数:
    fit_x, fit_y : 拟合曲线数据
    exp_x, exp_y : 实验数据
    method : 插值方法，默认线性插值
    
    返回:
    residuals, relative_residuals, r2, rmsd
    """
    # 拟合曲线插值到实验点
    if use_interp:
        interp_func = interp1d(fit_x, fit_y, kind=method, fill_value="extrapolate")
        fit_y_interp = interp_func(exp_x)
    else:   ##计算全局数据时，不能使用插值，必须在调用函数前保证数据点数和顺序的准确性
        fit_y_interp = fit_y
    
    # 残差
    residuals = exp_y - fit_y_interp
    
    # 相对残差
    with np.errstate(divide='ignore', invalid='ignore'):
        relative_residuals = np.where(exp_y != 0, residuals / exp_y, 0.0)
    
    # R²
    r2 = r2_score(exp_y, fit_y_interp)
    rmsd = np.sqrt(np.mean((residuals)**2))
    return residuals, relative_residuals, r2, rmsd

def compute_global_fit_metrics(exp_y, fit_y):
    """
    计算全局拟合评价指标（global metrics）

    Parameters
    ----------
    exp_y : array-like
        所有实验数据（可由多温度、多等温线拼接而成）
    fit_y : array-like
        与 exp_y 一一对应的拟合数据

    Returns
    -------
    metrics : dict
        {
            "residuals": ndarray,
            "relative_residuals": ndarray,
            "R2": float,
            "RMSD": float
        }
    """
    import numpy as np
    from sklearn.metrics import r2_score

    exp_y = np.asarray(exp_y)
    fit_y = np.asarray(fit_y)

    if exp_y.shape != fit_y.shape:
        raise ValueError("exp_y and fit_y must have the same shape.")

    # residuals
    residuals = exp_y - fit_y

    # relative residuals
    with np.errstate(divide='ignore', invalid='ignore'):
        relative_residuals = np.where(exp_y != 0, residuals / exp_y, 0.0)

    # R²
    r2 = r2_score(exp_y, fit_y)

    # RMSD
    rmsd = np.sqrt(np.mean(residuals ** 2))

    return residuals, relative_residuals, r2, rmsd