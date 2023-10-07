# File: statistics
# Author: jiang tao
# Time: 2023/6/25 8:59
# Desc: 药代动力学计算
import numpy as np
import pandas as pd
from pathlib import Path

from matplotlib.font_manager import fontManager

from nca import NonCompart, RawDataError

import warnings

from matplotlib import pyplot as plt

warnings.filterwarnings(action='ignore')

TYPE = ("Intra-day", "Inter-day", "Recovery", "Matrix effect")
COL1 = ("Component Name", "Sample Name", "Calculated Concentration")
COL2 = ("nominal_concentration", "basic", "Condition")
COL3 = ("Nominal concentration\n(ng/mL)", "Measured concentration\n(mean ± SD, ng/mL)", "Accuracy\n(RE, %)",
        "Precision\n(RSD, %)", "Extraction recovery", "Matrix effect(%)", "RSD(%)", "Stability (%)")
COL4 = ("样本编号", "时间点")

MARKERS = ("o", "^", "s", "*", "d", "h", "H", "P", "D", "X")
TABLE = "src/tables"
FIG = "src/images"

CAPTION = ("表 4-2 基质样本中多组分精密度及准确度",
           "表 4-3 回收率及基质效应",
           "表 4-4 不同存储条件下多组分在基质样本中的稳定性",
           "表 4-5 样本中分析物浓度水平（ng_mL）",
           "图 4-5 分析物药时曲线",
           "表 4-6 药代动力学参数")


def calc_mean_std_rsd(df: pd.DataFrame):
    df[COL1[1]] = df[COL1[1]].apply(lambda x: x.split('-')[0])
    df.rename(columns={COL1[2]: "mean"}, inplace=True)
    df["std"] = df['mean']
    df = df.groupby(by=[COL1[0], COL1[1]], as_index=False).agg(
        {'mean': 'mean', 'std': 'std'})
    df['rsd'] = df.apply(lambda x: (x['std'] / x['mean']) * 100, axis=1)
    return df


def fill_nominal_basic(df: pd.DataFrame, nc: pd.DataFrame, basic=None):
    df[COL2[0]] = 1
    if isinstance(basic, pd.DataFrame):
        df[COL2[1]] = 0
    for i, row in df.iterrows():
        cn, sn = row[COL1[0]], row[COL1[1]]
        df.loc[i, COL2[0]] = nc.loc[
            (nc[COL1[0]] == cn) & (nc[COL1[1]] == sn), COL2[0]].values
        if isinstance(basic, pd.DataFrame):
            df.loc[i, COL2[1]] = basic.loc[(cn, "BS"), :].values
    if isinstance(basic, pd.DataFrame):
        df["accuracy"] = df.apply(lambda x: ((x['mean'] - x[COL2[1]]) / x[COL2[0]]) * 100, axis=1)
    else:
        df["accuracy"] = df.apply(lambda x: (x['mean'] / x[COL2[0]]) * 100, axis=1)
    df = df.applymap(lambda x: f"{x:.2f}" if not isinstance(x, str) else x)
    return df


def fill_nominal_basic_stable(df: pd.DataFrame, nc: pd.DataFrame, basic=None):
    df[COL2[0]] = 1
    if isinstance(basic, pd.DataFrame):
        df[COL2[1]] = 0
    for i, row in df.iterrows():
        cn, sn, cond = row[COL1[0]], row[COL1[1]], row[COL2[2]]
        df.loc[i, COL2[0]] = nc.loc[
            (nc[COL1[0]] == cn) & (nc[COL1[1]] == sn), COL2[0]].values
        if isinstance(basic, pd.DataFrame):
            df.loc[i, COL2[1]] = basic.loc[(cn, cond), :].values
    if isinstance(basic, pd.DataFrame):
        df["accuracy"] = df.apply(lambda x: ((x['mean'] - x[COL2[1]]) / x[COL2[0]]) * 100, axis=1)
    else:
        df["accuracy"] = df.apply(lambda x: (x['mean'] / x[COL2[0]]) * 100, axis=1)
    df = df.applymap(lambda x: f"{x:.2f}" if not isinstance(x, str) else x)
    return df


def rename_index(df: pd.DataFrame, _type: str, col_map: tuple = None, basic=False):
    # col_map=("Matrix effect(%)","RSD(%)", ) Extraction recovery
    # rename
    df[COL3[0]] = df[COL2[0]]
    df[COL3[1]] = df.apply(
        lambda x: x['mean'] + " ± " + x['std'], axis=1)
    if col_map:
        df[col_map[0]] = df["accuracy"]
        df[col_map[1]] = df["rsd"]
    else:
        df[COL3[2]] = df["accuracy"]
        df[COL3[3]] = df["rsd"]
    # clean
    df.drop(columns=["mean", "std", "rsd", "accuracy", COL2[0], COL1[1]], inplace=True)
    if basic:
        df.drop(columns=[COL2[1],], inplace=True)
    df.set_index(keys=[COL1[0], COL3[0]], drop=True, inplace=True)
    df.columns = pd.MultiIndex.from_tuples([(_type, col) for col in df.columns])
    df.sort_index(inplace=True)
    return df


def rename_index_stable(df: pd.DataFrame, col_map: tuple = None, basic=False):
    # col_map=("Stability (%)", "Precision\n(RSD, %)", )
    # rename
    df[COL3[0]] = df[COL2[0]]
    df[COL3[1]] = df.apply(
        lambda x: x['mean'] + " ± " + x['std'], axis=1)
    if col_map:
        df[col_map[0]] = df["accuracy"]
        df[col_map[1]] = df["rsd"]
    else:
        df[COL3[2]] = df["accuracy"]
        df[COL3[3]] = df["rsd"]
    # clean
    df.drop(columns=["mean", "std", "rsd", "accuracy", COL2[0], COL1[1]], inplace=True)
    if basic:
        df.drop(columns=[COL2[1],], inplace=True)
    df.set_index(keys=[COL1[0], COL2[2], COL3[0]], drop=True, inplace=True)
    df.sort_index(inplace=True)
    return df


def calc_precision_accuracy(df: pd.DataFrame, nc: pd.DataFrame, _type: str):
    """
        _type: [Intra-day, Inter-day]
    """
    return rename_index(fill_nominal_basic(calc_mean_std_rsd(df), nc), _type)


def calc_recovery_matrix_effect(df: pd.DataFrame, nc: pd.DataFrame, basic: pd.DataFrame, _type: str, col_map: tuple):
    return rename_index(fill_nominal_basic(calc_mean_std_rsd(df), nc, basic=basic), _type, col_map=col_map, basic=True)


def get_precision_accuracy(df: pd.DataFrame, nc: pd.DataFrame, f):
    """
        df: sheet_name='精密度及准确度'
        nc: sheet_name='nominal_concentration'
    """
    data = df[~df[COL1[1]].str.contains('S\d{1,2}|BLK', regex=True)]
    # intraday part using Day1
    data_intraday = data[data[COL1[1]].str.contains('-D1-')]
    data_intraday = calc_precision_accuracy(data_intraday, nc, TYPE[0])

    # inter-day part using Day1-5
    data_interday = data.copy(deep=True)
    data_interday = calc_precision_accuracy(data_interday, nc, TYPE[1])

    df_precision_accuracy = pd.concat([data_intraday, data_interday], axis=1)

    df_precision_accuracy.to_excel(f"{f}/{TABLE}/{CAPTION[0]}.xlsx")
    return df_precision_accuracy


def get_recovery_matrix_effect(df: pd.DataFrame, nc: pd.DataFrame, f):
    """
        df: sheet_name='提取回收率及基质效应'
        nc: sheet_name='nominal_concentration'
    """
    df = df[~df[COL1[1]].str.contains('S\d{1,2}|BLK', regex=True)]
    # 空白基质
    basic = df[df[COL1[1]].str.contains('BS-\d{1,2}', regex=True)]
    basic[COL1[1]] = basic[COL1[1]].apply(lambda x: x.split('-')[0])
    basic = basic.groupby(by=[COL1[0], COL1[1]]).agg("mean")

    # 提取回收率
    recovery = df[df[COL1[1]].str.contains('-R-')]
    recovery = calc_recovery_matrix_effect(recovery, nc, basic, TYPE[2], col_map=(COL3[4], COL3[6]))
    # 基质效应
    matrix_effect = df[df[COL1[1]].str.contains('-M-')]
    matrix_effect = calc_recovery_matrix_effect(matrix_effect, nc, basic, TYPE[3], col_map=(COL3[5], COL3[6]))

    df_recovery_matrix_effect = pd.concat([recovery, matrix_effect], axis=1)

    df_recovery_matrix_effect.to_excel(f"{f}/{TABLE}/{CAPTION[1]}.xlsx")
    return df_recovery_matrix_effect


def get_stability(df: pd.DataFrame, nc: pd.DataFrame, f):
    df = df[~df["Sample Name"].str.contains('S\d{1,2}|BLK', regex=True)]
    # basic
    basic = df[df[COL1[1]].str.contains('BS-', regex=True)]
    basic[COL1[1]] = basic[COL1[1]].apply(lambda x: x.rsplit('-', 1)[0].replace("BS-", ""))
    basic = basic.groupby(by=[COL1[0], COL1[1]]).aggregate("mean")

    # data
    data = df[~df[COL1[1]].str.contains('BS-', regex=True)]
    data[COL2[2]] = data[COL1[1]].str.split("-", n=1, expand=True)[1].apply(lambda x: x.rsplit("-", 1)[0])
    data[COL1[1]] = data[COL1[1]].str.split("-", n=1, expand=True)[0]

    # calc_mean_std_rsd
    data.rename(columns={COL1[2]: "mean"}, inplace=True)
    data["std"] = data['mean']
    data = data.groupby(by=[COL1[0], COL2[2], COL1[1]], as_index=False).agg(
        {'mean': 'mean', 'std': 'std'})
    data['rsd'] = data.apply(lambda x: (x['std'] / x['mean']) * 100, axis=1)

    # fill_nominal_basic
    data = fill_nominal_basic_stable(data, nc, basic=basic)
    # rename
    df_stability = rename_index_stable(data, col_map=(COL3[7], COL3[3]), basic=True)

    df_stability.to_excel(f"{f}/{TABLE}/{CAPTION[2]}.xlsx")
    return df_stability


def get_analyte_cst(df: pd.DataFrame, f):
    data = pd.pivot(df, index=COL1[1], columns=COL1[0], values=COL1[2])
    data = data.reset_index().rename(columns={COL1[1]: COL4[0]})
    data.to_excel(f"{f}/{TABLE}/{CAPTION[3]}.xlsx", index=False)
    return data


def get_curve_data(x: pd.DataFrame, y:pd.DataFrame):
    # y = pd.read_excel(p, sheet_name='药代动力学DEMO数据')[["Sample Name", "Component Name", "Calculated Concentration"]]
    # x = pd.read_excel(p, sheet_name='药代动力学数据说明')[["样本编号", "时间点"]]

    x = dict(zip(x[COL4[0]], x[COL4[1]]))
    y["x"] = y[COL1[1]].map(x)
    y[COL1[1]] = y[COL1[1]].apply(lambda x: x[0])
    y["y"] = y[COL1[2]]
    y["std"] = y[COL1[2]]
    data = y.groupby(by=[COL1[0], COL1[1]]).agg({"x": "mean", "y": "mean", "std": "std"})
    return data


def plot_drug_time_curve(data: pd.DataFrame, f):
    fontManager.addfont("/usr/share/fonts/ARIALUNI.TTF")
    plt.rcParams["font.sans-serif"] = ["Arial Unicode MS"]
    # [0.0625, 0.125, 0.25, 0.5, 1, 2, 4, 8, 16, 32]
    global xticks, xticklables
    fig, ax = plt.subplots(1, 1, figsize=(8, 5), dpi=200)
    for i, _ in enumerate(data.index.levels[0]):
        r = data.loc[_]
        xticks = np.linspace(0, r["x"].size - 1, r["x"].size)
        xticklables = r["x"].to_list()

        x_transformed = xticks
        ax.errorbar(x_transformed, r["y"].to_list(), yerr=r["std"].to_list(),
                    elinewidth=0.6, label=_, capsize=2, capthick=0.8, fmt=MARKERS[i] + "-", ms=3, lw=0.8)
        # plt.xscale('log')
        # plt.yscale('log')

    ax.set_xticks(xticks, labels=xticklables)
    ax.set_xlabel("Hours(h)")
    ax.set_ylabel("Concentration(ng/mL)")
    # ax.spines['right'].set_visible(False)
    # ax.spines['top'].set_visible(False)
    ax.legend()
    plt.savefig(f"{f}/{FIG}/{CAPTION[4]}.jpg")
    plt.savefig(f"{f}/{FIG}/{CAPTION[4]}.pdf")
    plt.close()


def calc_kinetics_params(desc, demo, f):
    """get kinetics params of Non-compartmental Analysis for Pharmacokinetic Data"""
    timemap = dict(zip(desc[COL4[0]], desc[COL4[1]]))
    dosemap = dict(map(lambda x: (x[0].replace('给药剂量（mg）', ''), x[1]), desc.iloc[0, 2:5].to_dict().items()))

    res = pd.DataFrame()
    for gn, grp in demo.groupby(by=COL1[0]):
        dose = dosemap[gn]
        grp["time"] = grp[COL1[1]].map(timemap)
        grp["group"] = grp[COL1[1]].apply(lambda x: x[1:])
        data = []
        idx = []
        for i, g in grp.groupby(by=["group"]):
            x = g["time"].to_numpy()
            y = g[COL1[2]].to_numpy()
            nca = NonCompart(x, y, dose=dose)
            row = nca.NCA()
            data.append(row)
            idx.append(i)
        df = pd.DataFrame(data=data,
                          columns=["Cmax(ng/mL)", "Tmax(hr)", "Ke(1/hr)", "AUC 0~24(hr*ng/mL)", "AUC 0~∞(hr*ng/mL)",
                                   "Vz/F(L/kg)", "CLz/F(L/hr/kg)", "MRT0~24(hr)", "T1/2z(hr)"], index=idx)
        df_filtered = df.loc[(df != 0).any(axis=1)]
        if df_filtered.shape[0] < 3:
            raise RawDataError(f"{gn}组：超过3组数据异常， 请检查原始数据！")
        df_describe = df_filtered.describe()
        res[gn] = df_describe.loc["mean"].apply(lambda x: f"{x:.2f}") + "±" + df_describe.loc["std"].apply(
            lambda x: f"{x:.2f}")
    res.to_excel(f"{f}/{TABLE}/{CAPTION[5]}.xlsx")


if __name__ == '__main__':
    ...
    # used_cols = ["Sample Name", "Component Name", "Calculated Concentration"]
    # save dir
    # f = "./"
    # test
    # p = Path(r"E:\研究生阶段数据\PycharmProjects\pythonBase\pharmacokinetics\data.xlsx")
    # nc = pd.read_excel(p, sheet_name='nominal_concentration')
    #
    # df1 = pd.read_excel(p, sheet_name='精密度及准确度')[used_cols]
    # get_precision_accuracy(df1, nc, f)
    #
    # df2 = pd.read_excel(p, sheet_name='提取回收率及基质效应')[used_cols]
    # get_recovery_matrix_effect(df2, nc, f)
    #
    # df3 = pd.read_excel(p, sheet_name='稳定性')[used_cols]
    # get_stability(df3, nc, f)
    #
    # df4 = pd.read_excel(p, sheet_name='药代动力学数据')[used_cols]
    # get_analyte_cst(df4)

    # df5 = pd.read_excel(p, sheet_name='药代动力学数据说明')[["样本编号", "时间点"]]
    # curve_data = get_curve_data(df5, df4)
    # plot_drug_time_curve(curve_data, f)
    #
    # desc = pd.read_excel(p, sheet_name="药代动力学数据说明")
    # demo = pd.read_excel(p,  sheet_name="药代动力学DEMO数据")[["Sample Name", "Component Name", "Calculated Concentration"]]
    # calc_kinetics_params(desc, demo)

