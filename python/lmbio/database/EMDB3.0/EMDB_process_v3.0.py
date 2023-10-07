import argparse
import pandas as pd
import os
import sys
import csv
import mysql.connector
import matplotlib.pyplot as plt
import get_SDF
import get_msp
import numpy as np
import re

# 建立与MySQL数据库的连接
def sql_cnx():
    cnx = mysql.connector.connect(user='tanfj', password='lumingbio123',
                                    host='192.168.10.200', port=3306,
                                    database='cosa')
    return cnx

def cid_map():
    #添加cid对应的name和真实保留时间等信息
    dir_path = os.path.dirname(os.path.abspath(__file__))
    cnx = sql_cnx()
    # 查询四个CID转换列表的表格，并将表头重新命名
    cid_hmdbid_df = pd.read_sql('SELECT hmdbid as "Compound ID", cid as cid FROM cid_hmdbid', con=cnx)
    cid_keggid_df = pd.read_sql('SELECT keggid as "Compound ID", cid as cid FROM cid_keggid', con=cnx)
    cid_lipidmapsid_df = pd.read_sql('SELECT lipidmapsid as "Compound ID", cid as cid FROM cid_lipidmapsid', con=cnx)
    cid_metlinid_df = pd.read_sql('SELECT metlinid as "Compound ID", cid as cid FROM cid_metlinid', con=cnx)
    # 将四个表格合并成一个ID_map数据框
    ID_map = pd.concat([cid_hmdbid_df, cid_keggid_df, cid_lipidmapsid_df, cid_metlinid_df], ignore_index=True)
    ID_map['cid'] = ID_map['cid'].astype(str)
    #添加cid 对应的Compound Name
    cid_identifier_df = pd.read_sql('SELECT name as "Compound Name", cid as cid FROM compound_identifier', con=cnx)
    cid_identifier_df['cid'] = cid_identifier_df['cid'].astype(str)
    ID_map = pd.merge(ID_map,cid_identifier_df,on='cid', how='left')
    #添加cid对应的inchikey,用于msp里的DATABASE_ID
    cid_structure = pd.read_sql('SELECT inchikey,cid FROM compound_structure', con=cnx)
    cid_structure['cid'] = cid_structure['cid'].astype(str)
    ID_map = pd.merge(ID_map,cid_structure,on='cid', how='left')
    # 关闭数据库连接
    cnx.close()
    #添加cid对应的预测保留时间和QC保留时间信息
    RT_df = pd.read_csv('/data/hstore4/database/EMDB/RT.txt',dtype={'cid': str, 'Predict': float,'QC_RetentionTime':float},sep='\t')
    ID_map = pd.merge(ID_map,RT_df,on='cid', how='left')
    return ID_map

def sample_map(info_file='./非靶项目信息.xlsx'):
    #读取xlsx里所有sheet的样本信息
    xlsx_file = pd.ExcelFile(info_file)
    df_list = []
    for sheet_name in xlsx_file.sheet_names:
        df = pd.read_excel(xlsx_file, sheet_name=sheet_name, usecols=['项目编号/订单号', '样本'])
        df_list.append(df)
    combined_df = pd.concat(df_list, ignore_index=True)
    combined_df.rename(columns={'项目编号/订单号':'Project_num','样本':'Sample'},inplace=True)
    return combined_df

def hist_plot(unique_df=None,out_path='./'):
    # 计算trusted_df2里的时间差异
    #删除QC_RetentionTime为NA的数据
    unique_df.dropna(subset=['QC_RetentionTime'], inplace=True)
    diff_time = abs(unique_df['QC_RetentionTime'] - unique_df['Retention time (min)'])
    # 绘制直方图
    bins=np.arange(diff_time.min(),diff_time.max(), 0.5)
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.hist(diff_time, bins=bins, edgecolor='black', linewidth=1.2)
    ax.set_xlabel('diff_time')
    ax.set_ylabel('frequency')
    ax.set_title('')
    plt.savefig(os.path.join(out_path,'时间差异分布图.png'), dpi=300, bbox_inches='tight')

def process_files(directory='./', out_path='./', fragmentation_score_threshold=70):
    #输出文件
    duplicates_file = os.path.join(out_path,'duplicated.csv')
    unique_file = os.path.join(out_path,'trusted_unique.xlsx')
    if os.path.exists(duplicates_file):
        os.remove(duplicates_file)
    if os.path.exists(unique_file):
        os.remove(unique_file)
    ID_map_df = cid_map()

#step1：单个文件处理后存储trusted_df1
    trusted_df1 = pd.DataFrame()
    suspicious_df = pd.DataFrame()
    for file_name in os.listdir(directory):
        if not file_name.endswith('ID.csv'):
            continue
        try:
            data = pd.read_csv(os.path.join(directory,file_name),encoding='ISO-8859-1')
        except pd.errors.ParserError:
            print(f"Error reading file: {file_name}. Skipping...")
            continue
        #第一行可能不是表头

        if "Compound" not in data.columns:
            data = data.iloc[0:]
            data.columns = data.iloc[0]
            data = data.drop(0)
        data['Project']= file_name
        #提取项目号,格式EMDB-DLM20223236-POS-ID.csv
        data['Project_num'] = re.sub(r'(EMDB-|-POS|-ID.csv|-NEG)', '', file_name)
        #与sample_map合并，提取样本类型
        sample_map_df = sample_map(info_file='/data/hstore4/database/EMDB/非靶项目信息.xlsx')
        data = pd.merge(data, sample_map_df, on='Project_num', how='left')
        
        if file_name.endswith('POS-ID.csv'):
            data['Charge']= 'Positive'
        else:
            data['Charge']= 'Negative'
        data = data[['Compound','Compound ID','Adducts','Formula','Fragmentation Score','m/z','Retention time (min)','Project','Sample','Charge']]
        # 将'Fragmentation Score'列的数据类型更改为float
        data['Fragmentation Score'] = data['Fragmentation Score'].astype(float)
        data['Retention time (min)'] = data['Retention time (min)'].astype(float)
        data = data[data['Fragmentation Score']>=float(fragmentation_score_threshold)]
        data['Adducts'] = data['Adducts'].map(lambda x:x.split(',')[0])
        #每个文件将Compound ID转换成cid
        data.dropna(subset=['Compound ID'], inplace=True)
        data = pd.merge(data,ID_map_df, on='Compound ID', how='left')
        data.dropna(subset=['cid'], inplace=True)

        #同一个Compound同一个cid取Fragmentation Score大的那个
        df_max_score = data.loc[data.groupby(['Compound', 'cid'])['Fragmentation Score'].idxmax()]

        #根据Compound进行分组
        grouped_data = df_max_score.groupby('Compound')
        # 遍历每个Compound的数据
        for name, group in grouped_data:
            # 判断Compound对应一个cid的情况
            if len(group) == 1:
                trusted_df1 = trusted_df1.append(group)
            # 判断Compound对应多个cid的情况
            else:
                max_score = group['Fragmentation Score'].max()
                second_max_score = group['Fragmentation Score'].nlargest(2).iloc[-1]
                if max_score - second_max_score >= 20:
                    trusted_group = group[group['Fragmentation Score'] == max_score]
                    trusted_df1 = trusted_df1.append(trusted_group)
                else:
                    suspicious_df = suspicious_df.append(group)

# Step 2: 将上述所有的可信矩阵合并,并根据相同cid筛选
    grouped_df = trusted_df1.groupby(['cid'])
    trusted_df2 = pd.DataFrame()
    for name, group in grouped_df:
        rt_min = group['Retention time (min)'].min()
        rt_max = group['Retention time (min)'].max()
        rt_diff = rt_max - rt_min
        # 判断相同Compound ID对应多个Retention time（min）的情况
        if rt_diff <= 0.2:
            trusted_df2 = trusted_df2.append(group)
        else:
            #统计在0.2min区间的比例
            counts = []
            for idx, rt in enumerate(group['Retention time (min)']):
                count = sum(abs(rt - other_rt) <= 0.2 for other_rt in group['Retention time (min)'])
                counts.append(count)
            # 计算满足条件的元素比例
            ratios = [count / len(group['Retention time (min)']) for count in counts]
            # 寻找比例最大的元素
            max_ratio_index = ratios.index(max(ratios))
            # 输出比例最大的元素对应的数据框
            if max(ratios) >= 0.8:
                trusted_df2 = trusted_df2.append(group.iloc[max_ratio_index])
            else:
            #如果不符合0.2区间，则选择预测保留时间与真实保留时间小于0.2的物质
                trusted_data = group[abs(group['Predict'] - group['Retention time (min)']) < 0.2]
                suspicious_data = group[abs(group['Predict'] - group['Retention time (min)']) >= 0.2]
                if not trusted_data.empty:
                    trusted_df2 = trusted_df2.append(trusted_data)
                if not suspicious_data.empty:
                    suspicious_df = suspicious_df.append(suspicious_data)
    trusted_df2.to_excel(os.path.join(out_path,"trusted.xlsx"), index=False)
    suspicious_df.to_csv(duplicates_file, header=True, index=False)
#step3:对可信矩阵trusted_df2中cid合并数据，保留时间取均值，Fragmentation Score取最大值
    final_df = trusted_df2.groupby('cid').agg({'Compound Name':lambda x: ','.join(x.astype(str).unique()),
                                                'Compound ID':lambda x: ','.join(x.unique()),
                                                'Compound': lambda x: ','.join(x.unique()),
                                                'Adducts' : lambda x: ','.join(x.unique()),
                                                'Retention time (min)': 'mean',
                                                'Predict': 'mean',
                                                'QC_RetentionTime': 'mean',
                                                'Fragmentation Score': 'max',
                                                'm/z': lambda x: ','.join(x.astype(str).unique()),
                                                'Formula': lambda x: ','.join(x.astype(str).unique()),
                                                'inchikey': lambda x: ','.join(x.astype(str).unique()),
                                                'Project': lambda x: ','.join(x.unique()),
                                                'Sample': lambda x: ','.join(x.fillna('').astype(str).unique()),
                                                'smiles': lambda x: ','.join(x.astype(str).unique()),
                                                'Charge': lambda x: ','.join(x.unique())}).reset_index()
    #输出数据
    final_df.to_excel(unique_file, index=False)
    final_df.to_csv(os.path.join(out_path,"trusted_unique.csv"), header=True, index=False)

    #输出数据库
    RT_database = final_df[['Compound Name','inchikey','Retention time (min)']]
    RT_database.rename(columns={'inchikey':'Compound ID'},inplace=True)
    RT_database.to_excel(os.path.join(out_path,"RT_database.xlsx"), index=False)
    #根据trusted_df2评估标准品保留时间与真实保留时间的差异情况
    hist_plot(unique_df=trusted_df2,out_path=out_path)
    #输出final_df对应的sdf文件
    get_SDF.save_SDF(infile=os.path.join(out_path,"trusted_unique.csv"),outdir=out_path)
    
    #输出trusted_df2 中Fragmentation Score最高情况对应的二级碎片对应的msp文件
    trusted_df2_max_score = trusted_df2.loc[trusted_df2.groupby('cid')['Fragmentation Score'].idxmax()]
    #trusted_df2_max_score.to_csv("trusted_df2.csv", header=True, index=False)
    trusted_df2_max_score_pos = trusted_df2_max_score[trusted_df2_max_score["Charge"]=="Positive"]
    #trusted_df2_max_score_pos.to_csv("trusted_df2_pos.csv", header=True, index=False)
    get_msp.save_msp(folder_path=directory,trusted_df=trusted_df2_max_score_pos,output_path=out_path, charge = "POS")
    trusted_df2_max_score_neg = trusted_df2_max_score[trusted_df2_max_score["Charge"]=="Negative"]
    #trusted_df2_max_score_neg.to_csv("trusted_df2_neg.csv", header=True, index=False)
    get_msp.save_msp(folder_path=directory,trusted_df=trusted_df2_max_score_neg,output_path=out_path, charge = "NEG")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory",default = "./",help = "数据目录")
    parser.add_argument("-f","--fragmentation_score_threshold",default = 70, help = "fragmentation_score过滤阈值")
    parser.add_argument("-o","--out_path",default = "./", help = "输出目录")    
    args = parser.parse_args()
    process_files(directory=args.directory, 
                out_path=args.out_path, 
                fragmentation_score_threshold=args.fragmentation_score_threshold)
