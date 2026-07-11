from common import *

def henry_model(P, KH):
    return KH * P

def fit_henry_model(df):
    """henry拟合"""
    from lmfit import Model
    model = Model(henry_model)
    P = df["P_Pa"].values
    q = df["q_mmol_g"].values
    params = model.make_params(KH=np.max(q) / np.max(P))
    # 限制参数范围（Henry常数应为正）
    params["KH"].min = 0
    params["KH"].max = 1e3  # 可根据量纲调节
    result = model.fit(q, params, P=P, method="leastsq")
    return result

def langmuir_model(P, q_sat, b):
    return q_sat * b * P / (1 + b * P)

def fit_single_langmuir(df):
    """单温度Langmuir拟合"""
    from lmfit import Model
    model = Model(langmuir_model)
    P = df["P_Pa"].values
    q = df["q_mmol_g"].values
    params = model.make_params(q_sat=np.max(q), b=1e-5)
    params["q_sat"].min = 0
    params["b"].min = 1e-10
    params["b"].max = 1e3
    result = model.fit(q, params, P=P, method="leastsq")
    return result

def get_initial_dsl_params_multitemp(T_list, grouped_data):
    """多温度Langmuir平均初值"""
    results = []
    for T in T_list:
        df = grouped_data[T]
        res = fit_single_langmuir(df)
        q_sat_L = res.params["q_sat"].value
        b_L = res.params["b"].value
        q_obs = df["q_mmol_g"].values
        q_fit = res.best_fit
        ss_res = np.sum((q_obs - q_fit)**2)
        ss_tot = np.sum((q_obs - np.mean(q_obs))**2)
        R2 = 1 - ss_res/ss_tot if ss_tot>0 else np.nan
        results.append({"T_K": T, "q_sat_L": q_sat_L, "b_L": b_L, "R2": R2})
    df_result = pd.DataFrame(results).sort_values("T_K")
    q_sat_mean = df_result["q_sat_L"].mean()
    b_mean = df_result["b_L"].mean()
    init_params = {
        "q1_sat": 0.7 * q_sat_mean,
        "q2_sat": 0.3 * q_sat_mean,
        "b1": 10 * b_mean,
        "b2": 0.1 * b_mean,
    }
    return init_params, df_result