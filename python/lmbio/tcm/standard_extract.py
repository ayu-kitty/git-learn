#!/opt/conda/bin/python
import numpy as np
import pandas as pd
import sys
import pymzml
import datetime as dt
import os
import collections
from collections import OrderedDict
import glob
# from sqlalchemy import create_engine
from sortedcontainers import SortedDict
import argparse
import regex as re
from examine_filename import examin

'''

将质谱RAW文件转为mzml格式，并根据mz搜索RT和二级质谱
根据正负离子的RT，筛选合并得到化合物的RT。

Parameters
----------
file : sample info table

Comments
----------
tol 正离子 10ppm,负离子 15ppm
二级图谱搜索范围，前后共100张

'''

precursortype = ['M+K','M+Na','M+NH4','M+H','M+H-H2O','2M+H','M+2H','M+H-2H2O','M+',
                 'M+FA-H','M-H2O-H','M+Cl','M-H','2M-H','M-2H','M-']

precursortype_ = ["[M+K]+","[M+Na]+","[M+NH4]+","[M+H]+","[M+H-H2O]+","[2M+H]+","[M+2H]2+",
"[M+H-2H2O]+","[M]+","[M+FA-H]-","[M-H2O-H]-","[M+Cl]-","[M-H]-","[2M-H]-","[M-2H]2-","[M]-"]                 

atom_d = {"H": 1.00782503207, "He": 4.00260325415, "Li": 7.01600455, "Be": 9.0121822, "B": 11.0093054, "C": 12.0,
"N": 14.0030740048, "O": 15.99491461956, "F": 18.99840322, "Ne": 19.9924401754, "Na": 22.9897692809, "Mg": 23.9850417,
"Al": 26.98153863, "Si": 27.9769265325, "P": 30.97376163, "S": 31.972071, "Cl": 34.96885268, "Ar": 39.9623831225, 
"K": 38.96370668, "Ca": 39.96259098, "Sc": 44.9559119, "Ti": 47.9479463, "V": 50.9439595, "Cr": 51.9405075,
"Mn": 54.9380451, "Fe": 55.9349375, "Co": 58.933195, "Ni": 57.9353429, "Cu": 62.9295975, "Zn": 63.9291422,
"Ga": 68.9255736, "Ge": 73.9211778, "As": 74.9215965, "Se": 79.9165213, "Br": 78.9183371, "Kr": 83.911507, 
"Rb": 84.911789738, "Sr": 87.9056121, "Y": 88.9058483, "Zr": 89.9047044, "Nb": 92.9063781, "Mo": 97.9054082,
"Ru": 101.9043493, "Rh": 102.905504, "Pd": 105.903486, "Ag": 106.905097, "Cd": 113.9033585, "In": 114.903878, 
"Sn": 119.9021947, "Sb": 120.9038157, "Te": 129.9062244, "I": 126.904473, "Xe": 131.9041535, "Cs": 132.905451933,
"Ba": 137.9052472, "La": 138.9063533, "Ce": 139.9054387, "Pr": 140.9076528, "Nd": 141.9077233, "Sm": 151.9197324, 
"Eu": 152.9212303, "Gd": 157.9241039, "Tb": 158.9253468, "Dy": 163.9291748, "Ho": 164.9303221, "Er": 165.9302931, 
"Tm": 168.9342133, "Yb": 173.9388621, "Lu": 174.9407718, "Hf": 179.94655, "Ta": 180.9479958, "W": 183.9509312, 
"Re": 186.9557531, "Os": 191.9614807, "Ir": 192.9629264, "Pt": 194.9647911, "Au": 196.9665687, "Hg": 201.970643, 
"Tl": 204.9744275, "Pb": 207.9766521, "Bi": 208.9803987, "Th": 232.0380553, "Pa": 231.035884, "U": 238.0507882}

def mw(f):
    global atom_d
    f_list = re.findall(r'([A-Z][a-z]?)(\d*)',f)
    #转化为整数，且将空格替换为1
    for_list = [(atom,int(num) if num != "" else 1)for atom,num in f_list]
    mw_list = [atom_d[atom]*num for atom,num in for_list]
    return round(sum(mw_list),7)

def read_sample_info(sample_info):
    info = pd.read_excel(sample_info,header = 0, sheet_name=0, engine="openpyxl")
    if "M.Wt" not in info.columns:
        print("输入的表格无M.Wt,自动计算中！")
        info["M.Wt"] = info["Formula"].apply(mw)
    info.loc[:, "Formula"] = [_.strip() for _ in info["Formula"]]
    # info.set_index(pd.MultiIndex.from_arrays([info.Filename, info.Formula]), inplace = True)
    info.set_index(["Filename", "Formula"], inplace=True)
    return info

def convert_mzml():
    root_dir = os.getcwd()
    sep = os.path.sep
    file_name_pos = [os.path.basename(file).strip(r".raw") for file in glob.glob("pos/*.raw")]
    file_name_neg = [os.path.basename(file).strip(r".raw") for file in glob.glob("neg/*.raw")]
    if file_name_pos :
        file_name = file_name_pos
    elif file_name_neg:
        file_name = file_name_neg
    if len(file_name_pos) > 0 and len(file_name_neg) > 0:
        file_name = list(set(file_name_neg) | set(file_name_pos))
    elif len(file_name_pos) == 0 and len(file_name_neg) == 0:
        print("正负离子文件不存在！")
        sys.exit(1)
    file_name.sort()
    print(file_name)    
    for file in file_name:
        # 判断pymz文件是否存在，不存在则转换。
        if not os.path.exists("pos{sep}{file}.mzml".format(sep=sep, file=file)):
            cmd_pos = "msconvert pos{sep}{file}.raw --zlib -o pos".format(sep=sep, file=file)
            os.system(cmd_pos)
            print("{}.raw转换格式完成！".format(file))
        if not os.path.exists("neg{sep}{file}.mzml".format(sep=sep, file=file)):
            cmd_neg = "msconvert neg{sep}{file}.raw --zlib -o neg".format(sep=sep, file=file)
            os.system(cmd_neg)
            print("{}.raw转换格式完成！".format(file))
    print("所有RAW文件转换完成")
    pos_mzml = [root_dir + sep + _ for _ in glob.glob("pos/*.mzml")]
    neg_mzml = [root_dir + sep + _ for _ in glob.glob("neg/*.mzml")]
    return file_name, pos_mzml, neg_mzml


def read_from_db(sample_info,file_name):
    # 读取sample_info和pos&neg下的RAW文件，并按实际file_name对齐
    with engine.begin() as con:    
        info = read_sample_info(sample_info)
        info = info.loc[file_name]
        all_add = []

        for group in file_name:
            group_info = info.loc[(group,)]
            group_all = pd.DataFrame()
            for r in range(group_info.shape[0]):
                print(f"{r}" * 10)
                #以HMDB号开头或者纯数字形式<pubchemID>，按ID提取
                if re.match(r"HMDB|\d+",str(group_info.Identifiers.iloc[r])):
                    query = "SELECT * FROM `metamz-all` WHERE `id` = '{id}'".format(id=group_info.Identifiers.iloc[r])
                    print(query)
                    adduct = pd.read_sql_query(query, con)
                #其他ID号，以分子式提取  
                elif group_info.Formula.iloc[r]:
                    query = "SELECT * FROM `metamz-all` WHERE `Formula` = '{formula}'".format(
                        formula=group_info.Formula.iloc[r])
                    print(query)
                    adduct = pd.read_sql_query(query, con)
                adduct = adduct.drop_duplicates(subset="adduct", ignore_index=True).set_index("adduct",drop = False)
                adduct = adduct.reindex(precursortype,axis = 0).dropna(axis = 0,how = "all").reset_index(drop = True)
                adduct = adduct.sort_values("mode",axis = 0,ascending = False).drop(columns=["DatabaseFrom"])
                adduct.set_index(["Formula"], inplace=True)
                group_all = pd.concat([group_all, adduct])
            all_add.append(group_all)
        all_g = pd.concat(all_add, keys=file_name)
        all_g.index.names = ["Filename", "Formula"]
    return all_g

# 读取sample_info,并全部转添加adduct信息
def read_from_local(sample_info):
    #定义正负离子mz的计算数据
    pos_ion_d = {"[M+H]+":(1,1,mw("H")),"[M+NH4]+":(1,1,mw("NH4")),"[M+Na]+":(1,1,mw("Na")),"[M+H-H2O]+":(1,1,mw("H")-mw("H2O")),
                "[M+H-2H2O]+":(1,1,mw("H")-mw("H2O")*2),"[2M+H]+":(2,1,mw("H")),"[M+2H]2+":(1,2,mw("H")),"[M+K]+":(1,1,mw("K"))}
    neg_ion_d = {"[M-H]-":(1,1,-mw("H")),"[M-H2O-H]-":(1,1,-mw("H2O")-mw("H")),"[M+Cl]-":(1,1,mw("Cl")),"[M+FA-H]-":(1,1,mw("CH2O2")-mw("H")),
                "[2M-H]-":(2,1,-mw("H")),"[M-2H]2-":(1,2,-mw("H"))}

    ion_d = dict(pos_ion_d,**neg_ion_d)

    info_t = read_sample_info(sample_info)
    #按照sample_info逐项添加
    cols = info_t.columns.tolist() + ["adduct","charge","mode","mz"]
    
    all_g = pd.DataFrame()
    for index,s in info_t.iterrows():  
        df0 = pd.DataFrame(columns = cols)

        for key,(n,c,ms) in ion_d.items():
            s["adduct"] = key
            s["charge"] = c
            s["mode"] = "pos" if key.endswith("+") else "neg"
            s["mz"] = round(s["M.Wt"]*n/c+ms,7)
            df0 = pd.concat([df0,pd.DataFrame(s).T])
        #对于分子式以+结尾的，考虑M+离子    
        if s.name[1].endswith("+"):
            s["adduct"] = "M+"
            s["mode"] = "pos"
            s["mz"] = s["M.Wt"]
            s["charge"] = 1
            df0 = df0.append(s)
        elif s.name[1].endswith("-"):
            s["adduct"] = "[M]-"
            s["mode"] = "neg"
            s["mz"] = s["M.Wt"]
            s["charge"] = 1
            df0 = df0.append(s)
        df0.index = pd.MultiIndex.from_tuples(df0.index)
        all_g = all_g.append(df0)
        all_g.index.names = ["Filename", "Formula"]
    return all_g

def extract_peaks(mzml_file):
    run = pymzml.run.Reader(mzml_file)
    # level1,level2水平的质谱图
    level1_spectr = {x.ID: x for x in run if x.ms_level == 1}
    level2_spectr = {x.ID: x for x in run if x.ms_level == 2}

    def map_ll(level1, level2):
        ll = SortedDict()
        for i in range(len(level1) - 1):
            tmp = []
            for j in level2:
                if level1[i] < j and j < level1[i + 1]:
                    tmp.append(j)
            ll[level1[i]] = tmp
        return ll

    mapper = map_ll(list(level1_spectr.keys()), list(level2_spectr.keys()))
    # 提取centroid之后一级质谱峰
    level1_peaks_ctr = collections.defaultdict(list)
    for spectr in run:
        if spectr.ms_level == 1:
            level1_peaks_ctr[spectr.ID] = spectr.peaks("centroided")
    return level1_spectr, level2_spectr, mapper, level1_peaks_ctr


# 从保留时间列表中筛选间隔大于0.1的峰,返回ID-RT-peaks元组
def highest(rt, potential_peaks):
    df = pd.DataFrame([[abs(i - j) for i in rt.values()] for j in rt.values()])
    tmp = []
    for c, col in df.iteritems():
        if all(col[0:c] > 0.1):
            tmp.append(c)
    id_rt = [list(rt.items())[i] for i in tmp]

    return [[idx, rt, potential_peaks[idx].tolist()[0]] for idx, rt in id_rt]


def select_mz(adduct_mz,level1_spectr, level1_peaks_ctr, level2_spectr, mapper,tol=10e-6):
    # 设置加合物和误差度
    # tol = 1e-5
    adduct_max = adduct_mz * (1 + tol)
    adduct_min = adduct_mz * (1 - tol)
    # 提取范围内的加合物峰
    select_peaks = collections.defaultdict(list)
    select_peaks = {idx: peak[(adduct_min < peak[:, 0]) & (peak[:, 0] < adduct_max)] for idx, peak in
                    level1_peaks_ctr.items()}
    select_peaks = {key: value for key, value in select_peaks.items() if value.size > 0}

    # 返回噪音强度 和 离子峰强,筛选SN > 10的峰
    select_spectr = {idx: level1_spectr[idx] for idx in select_peaks.keys()}
    noise = {idx: select_spectr[idx].estimated_noise_level() for idx in select_spectr.keys()}
    intensity = {i: peak[0, 1] for i, peak in select_peaks.items()}
    potential_peak_i = {i: it for i, it in intensity.items() if it > noise[i] * 10}
    potential_peaks = {i: select_peaks[i] for i, it in potential_peak_i.items()}
    # 排序过的强度和保留时间
    inten0 = sorted(list(potential_peak_i.items()), key=lambda x: x[1], reverse=True)
    inten = collections.OrderedDict(inten0)
    RT0 = {i: level1_spectr[i].scan_time[0] for i in [x[0] for x in inten0]}
    RT = collections.OrderedDict(RT0)
    ####若是higest返回多组RT/peak
    h = highest(RT, potential_peaks)
    id_ = None
    rt0 = rt1 = rt2 = mz = i_ = np.nan
    if h:
        id_, rt0, (mz, i_) = h[0]
    if len(h) == 2:
        _, rt1, (_, _) = h[1]
    if len(h) == 3:
        _, rt2, (_, _) = h[2]
    rt = (rt0, rt1, rt2)

    def extract_ms2(mz=None, ID=None, mapper=None, level2_spectr=None):
        # 取当前一级ID以及前后五张一级ID对应的二级图谱
        if ID and ID in mapper.keys():
            ID_pos = list(mapper.keys()).index(ID)
            level2_list = []
            k = 1
            while len(level2_list) < 100:
                level1_list = mapper.keys()[ID_pos - k:ID_pos + k]
                level2_list = [j for i in level1_list for j in mapper[i]]
                k += 1
            pres = {id2: level2_spectr[id2].selected_precursors for id2 in level2_list}
            lv2_spectrs = [level2_spectr[id2] for id2, pre in pres.items() if
                           pre[0]["mz"] <= mz * (1 + tol) and pre[0]["mz"] >= mz * (1 - tol)]
            lv2_spectrs = [s for s in lv2_spectrs if s.has_peak(mz)]
            lv2_spectrs.sort(key=lambda spc: spc.has_peak(mz)[0][1], reverse=True)
            if len(lv2_spectrs) == 0:
                # 针对有1级碎片 无二级碎片
                return ["NO MS2","NO MS2","","","",""]
            elif len(lv2_spectrs) == 1:
                lv2_peak = lv2_spectrs[0].peaks("centroided")
                return (",".join([str(_) for _ in lv2_peak[:, 0].tolist()]),
                        ",".join([str(_) for _ in lv2_peak[:, 1].tolist()]),lv2_spectrs[0].ID,"","","")
            else:
                lv2_peak1 = lv2_spectrs[0].peaks("centroided")
                lv2_peak2 = lv2_spectrs[1].peaks("centroided")
                return (",".join([str(_) for _ in lv2_peak1[:, 0].tolist()]),
                        ",".join([str(_) for _ in lv2_peak1[:, 1].tolist()]),
                        lv2_spectrs[0].ID,
                        ",".join([str(_) for _ in lv2_peak2[:, 0].tolist()]),
                        ",".join([str(_) for _ in lv2_peak2[:, 1].tolist()]),
                        lv2_spectrs[1].ID)
        else:
            return [""] * 6

    level2_mz1, level2_i1,level2_id1,level2_mz2, level2_i2,level2_id2 = extract_ms2(mz, id_, mapper, level2_spectr)

    return rt, mz, i_, level2_mz1, level2_i1,level2_id1,level2_mz2, level2_i2,level2_id2


def extract(mode, pos_mzml=None, neg_mzml=None, file_name=None, info_df=None):
    partial_res = []
    keys = []
    if mode == "pos":
        for mzml in pos_mzml:
            #按照实际的pos_mzml对file_name进行筛选
            level1_spectr, level2_spectr, mapper, level1_peaks_ctr = extract_peaks(mzml)
            #处理G10文件时会被G1(file_name)替换掉。
            # file = [f for f in file_name if f in mzml][0]
            f_n = os.path.basename(mzml).split(".")[0]
            print(f_n)
            file = [f for f in file_name if f == f_n][0]
            #对每一组中的化合物进行迭代
            for idx in info_df.loc[file].index.drop_duplicates():
                pos_ion = []
                if isinstance(info_df.loc[file].loc[idx].adduct,str):
                    # 处理加合物只有M+这种情况,adduct只有一行,为str
                    adduct_mz = info_df.loc[file].loc[idx].mz
                    adduct_mode = info_df.loc[file].loc[idx].adduct
                    print(adduct_mz)
                    (rt1, rt2, rt3), mz, i_, level2_mz, level2_i,level2_id1,level2_mz2, level2_i2,level2_id2 = \
                        select_mz(adduct_mz, level1_spectr, level1_peaks_ctr,level2_spectr, mapper,tol = 10e-6)
                    id_rt_peaks = tuple([adduct_mode, rt1, rt2, rt3, mz, i_, level2_mz, level2_i,level2_id1,level2_mz2, level2_i2,level2_id2])
                    pos_ion.append(id_rt_peaks)
                
                else:
                    # 8种正离子在组内的index为0 - 8
                    #mz在第6列,adduct_mode在第三列
                    tmp_ = info_df.loc[file].loc[idx]
                    for idx,row in tmp_[tmp_["mode"] == "pos"].iterrows():
                        adduct_mz =row.mz
                        adduct_mode = row.adduct
                        print(adduct_mz)
                        (rt1, rt2, rt3), mz, i_, level2_mz, level2_i,level2_id1,level2_mz2, level2_i2,level2_id2= \
                            select_mz(adduct_mz, level1_spectr, level1_peaks_ctr,level2_spectr, mapper,tol = 10e-6)
                        id_rt_peaks = tuple([adduct_mode, rt1, rt2, rt3, mz, i_, level2_mz, level2_i,level2_id1,level2_mz2, level2_i2,level2_id2])
                        pos_ion.append(id_rt_peaks)
                tmp = pd.DataFrame(pos_ion, columns=['adduct', 'RT1', 'RT2', 'RT3', 'mz_measure', 'Intensity_measure',
                                                     'ms_level2_MZ1', 'ms_level2_Intensity1','level2_id1','ms_level2_MZ2', 'ms_level2_Intensity2','level2_id2'])
                merged = info_df.loc[file].loc[[idx]].merge(tmp, on="adduct")
                partial_res.append(merged)
                keys.append((file, idx))
        ret = pd.concat(partial_res, keys=keys)
        ret.index.names = ["Filename","Formula","ind"]
        ret.reset_index(inplace = True)
    if mode == "neg":
        for mzml in neg_mzml:
            level1_spectr, level2_spectr, mapper, level1_peaks_ctr = extract_peaks(mzml)
            # file = [f for f in file_name if f in mzml][0]
            f_n = os.path.basename(mzml).split(".")[0]
            print(f_n)
            file = [f for f in file_name if f == f_n][0]
            for idx in info_df.loc[file].index.drop_duplicates():
                neg_ion = []
                if isinstance(info_df.loc[file].loc[idx].adduct, str):
                    # 处理加合物只有M-这种情况,adduct只有一行,为str
                    adduct_mz = info_df.loc[file].loc[idx].mz
                    adduct_mode = info_df.loc[file].loc[idx].adduct
                    print(adduct_mz)
                    (rt1, rt2, rt3), mz, i_, level2_mz, level2_i,level2_id1,level2_mz2, level2_i2,level2_id2 = \
                        select_mz(adduct_mz, level1_spectr, level1_peaks_ctr,level2_spectr, mapper,tol = 15e-6)
                    id_rt_peaks = tuple([adduct_mode, rt1, rt2, rt3, mz, i_, level2_mz, level2_i,level2_id1,level2_mz2, level2_i2,level2_id2])
                    neg_ion.append(id_rt_peaks)
                else:
                    # 负离子在组内的index为8-14
                    tmp_ = info_df.loc[file].loc[idx]
                    for idx,row in tmp_[tmp_["mode"] == "neg"].iterrows():
                        adduct_mz =row.mz
                        adduct_mode = row.adduct
                        print(adduct_mz)
                        (rt1, rt2,
                         rt3), mz, i_, level2_mz, level2_i, level2_id1, level2_mz2, level2_i2, level2_id2 = select_mz(
                            adduct_mz, level1_spectr, level1_peaks_ctr,level2_spectr, mapper,tol = 15e-6)
                        id_rt_peaks = tuple(
                            [adduct_mode, rt1, rt2, rt3, mz, i_, level2_mz, level2_i, level2_id1, level2_mz2, level2_i2,level2_id2])
                        neg_ion.append(id_rt_peaks)
                tmp = pd.DataFrame(neg_ion, columns=['adduct', 'RT1', 'RT2', 'RT3', 'mz_measure', 'Intensity_measure',
                                                         'ms_level2_MZ1', 'ms_level2_Intensity1', 'level2_id1',
                                                         'ms_level2_MZ2', 'ms_level2_Intensity2', 'level2_id2'])
                partial_res.append(info_df.loc[file].loc[[idx]].merge(tmp, on="adduct"))
                keys.append((file, idx))
        ret = pd.concat(partial_res, keys=keys)
        ret.index.names = ["Filename","Formula","ind"]
        ret.reset_index(inplace = True)
    return ret

def filter_rt(ion,thresh = 0.1):
    #处理不止一个保留时间
    #提取RT1和RT2保留时间
    ionv = ion.loc[pd.notnull(ion.RT1),:]
    ion_d1 = OrderedDict({k:v for k,v in zip(ionv.adduct,ionv.RT1)})
    ionv2 = ion.loc[pd.notnull(ion.RT2),:]
    ion_d2 = OrderedDict({k:v for k,v in enumerate(ionv2.RT2) })
    ion_d1.update(ion_d2)
    
    def compare_rt(dic=None,thresh = None):
        #查找RT1之间或者RT1与RT2之间相差0.1min之内的保留时间
        df = pd.DataFrame([[abs(i-j) for i in dic.values()] for j in dic.values()])
        for c, col in df.iteritems():
            if any(col[c+1:] < thresh):
                return list(dic.items())[c]

    res = compare_rt(ion_d1,thresh)
    if res:
        ion.loc[ion["adduct"] == res[0],"RT3"] = res[1]
    else:
        #若无相近RT 则看单个峰强度 RT3需要剔除强度小于10^7峰
        if any(ion.Intensity_measure >= 10**7):
            print("found single strong peak!")
            ion.loc[ion.Intensity_measure > 10**7,"RT3"] = ion.loc[ion.Intensity_measure > 10**7,"RT1"]
    return ion

def compare_pos_neg(group):
    pn = group.sort_values(by = "Intensity_measure",ascending = False)
    if sum(pd.notnull(pn.RT1)) <= 1:
        if any(pn.Intensity_measure > 10**7):
            print("only one RT retained!")
            pn.RT3 = pn.RT1
    else:
        pn = filter_rt(pn,thresh = 0.1)
    return pn

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-f','--file',help = "Sample info file")
    args = parser.parse_args()

    # engine = create_engine("mysql+pymysql://lujw:lumingbio@192.168.10.200/meta")
    # engine =  create_engine("mysql+pymysql://xd:xd2021@localhost/msms")
    
    # 处理原始文件
    t0 = dt.datetime.now().replace(microsecond=0)
    print("开始处理时间: {}".format(t0.strftime("%Y-%m-%d %H:%M:%S")))
    try:
        examin()
    except AssertionError:
        print("--"*30)
        print("文件名有误！请先修改！") 
        sys.exit(-1)
    file_name, pos_mzml, neg_mzml = convert_mzml()
    
    # 信息表
    print("开始读入sample文件")
    info_table = args.file
    if os.path.exists("info_df.csv"):
        info_df = pd.read_csv("info_df.csv",index_col=["Filename","Formula"])
    else:
        # info_df = read_from_db(info_table,file_name)
        info_df = read_from_local(info_table)
        info_df.to_csv("info_df.csv")
    print("读入sample文件完成!用时{}".format(str(dt.datetime.now().replace(microsecond=0)-t0)))
    # sys.exit(0)
    t1 = dt.datetime.now().replace(microsecond=0)
    
    if pos_mzml:
        print("开始处理正离子")
        pos = extract("pos", pos_mzml=pos_mzml, file_name=file_name, info_df=info_df)
        print("正离子处理完！")
    else:
        pos = pd.DataFrame()
        print("无正离子！")
    if neg_mzml:
        print("开始处理负离子")
        neg = extract("neg", neg_mzml=neg_mzml, file_name=file_name, info_df=info_df)
        print("负离子处理完！")
    else:
        neg = pd.DataFrame()
        print("无负离子！")
    print("筛选mz完成!用时{}".format(str(dt.datetime.now().replace(microsecond=0)-t1)))
    #合并正负离子
    print("开始合并正负离子！")
    if pos.size > 0  or neg.size > 0:    
        total = pd.concat([pos, neg], axis=0)
        # total.rename({"level_0":"Filename","level_1":"Formula"},axis = 1,inplace = True)
        total.drop("ind",axis = 1,inplace = True)
        total.to_csv("pos_neg_merged.csv",index = False)
    else:
        total = pd.DataFrame()
    #按照样品信息表Compound Name对提取的结果中数据库中的name进行校正
    #同时起到校验HMDB号和分子式是否匹配的作用
    #如果是原表格中HMDB号或分子式填错，则查找会不成功报错；
    '''
    total = total.set_index(["Filename","Formula"]).sort_index(axis = 0)
    grouped_df = total.groupby(["Filename","Formula"],sort = False)
    df0 = read_sample_info(args.file)
    name = []

    for indx,row in grouped_df:
        # name = gdf.get_group(indx)["Name"].unique()[0]
        compdname = df0.loc[indx]["Compound Name"]
        name.extend([compdname for i in row.Name])
    total["Name"] = name
    '''
    #读取数据并分组
    total.loc[:,"RT3"] = np.nan
    cols = total.columns
    # 获取filename-formula对
    x = total.loc[:,["Filename","Formula"]].to_numpy()
    df = pd.DataFrame(None,columns = cols)
    for _,df0 in total.groupby(by=["Filename","Formula"]):
        df1 = compare_pos_neg(df0)
        df = pd.concat([df,df1],axis = 0)

    #保存同名的CSV文件
    fname = os.path.splitext(os.path.basename(args.file))[0]
    df.to_csv("{}.xls".format(fname),index = False,sep='\t')
    t2 = dt.datetime.now().replace(microsecond=0)
    print("{0}.csv文件保存完成,当前时间{1},总用时{2}".format(fname,t2.strftime("%Y-%m-%d %H:%M:%S"),str(t2-t0)))

