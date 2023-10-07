import os
import regex as re
import argparse
import pandas as pd
import csv
def save_msp(folder_path='./',trusted_df=None,output_path='./',charge = None):
    '''
    folder_path:msp文件路径
    unique_df:最终合并后的可信数据
    output_path：输出路径
    '''
    if os.path.exists(os.path.join(output_path,f'{charge}_trusted_df.msp')):
        os.remove(os.path.join(output_path,f'{charge}_trusted_df.msp'))
    # 读取所有msp文件
    msp_files = [f for f in os.listdir(folder_path) if f.endswith(".msp")]
    name_regex = re.compile(r"\((.*?)\)")
    # 创建一个空字典，用于存储msp数据
    msp_dict = {}
    # 循环遍历所有msp文件
    for msp_file in msp_files:
        with open(os.path.join(folder_path, msp_file), "r") as f:
            try:
                content = f.read()
            except UnicodeDecodeError:
                print(f"Error reading file: {msp_file}. Skipping...")
                continue
        # 按空行分割数据块
        data_blocks = content.split("\n\n")
        #最后一个为空，删除
        data_blocks.pop()
        # 循环遍历每个数据块
        for data_block in data_blocks:
            # 提取Name和value
            lines = data_block.split("\n")
            name_line = lines[0]
            name_match = name_regex.search(name_line)
            if name_match:
                name = "NAME:" + name_match.group(1)
                value = "\n".join(lines[1:])
                msp_dict[name] = value


    # 循环遍历trusted_df中的每一行,有重复的Compound
    trusted_df.drop_duplicates(subset='Compound', inplace=True)
    new_dict ={}
    for index, row in trusted_df.iterrows():
        # 提取对应的msp数据
        msp_key = "NAME:" + row["Compound"]
        ####
        if msp_key in msp_dict:
            try:
                msp_value = msp_dict[msp_key].split("\n")
            except AttributeError:
                print(f"Error : {msp_key}. Skipping...")
                print (msp_dict[msp_key])
                continue
            if any(line.startswith("Precursor_type:") for line in msp_value):
                Precursor_type = msp_value[0].split(':')[1]
            else:
                if row["Charge"] == "Positive":
                    Precursor_type = "["+row["Adducts"]+"]+"
                else:
                    Precursor_type = "["+row["Adducts"]+"]-"
            PrecursorMZ = msp_value[1].split(':')[1]
            msp_value = msp_value[3:]
            # 在msp_value数组中插入DATABASE_ID和Compound ID
            msp_value.insert(0, "NAME: "+str(row["Compound Name"]))
            msp_value.insert(1, "DATABASE_ID: "+row["inchikey"])
            msp_value.insert(2, f"PRECURSORMZ:{PrecursorMZ}")
            msp_value.insert(3, f"PRECURSORTYPE:{Precursor_type}")
            msp_value.insert(4, "FORMULA: "+row["Formula"])
            msp_value.insert(5, "RETENTIONTIME: "+str(row["Retention time (min)"]))
            msp_value.insert(6, "IONMODE: "+row["Charge"])
            msp_value.insert(7, "SMILES: "+str(row["smiles"]))
            new_dict[msp_key]="\n".join(msp_value)
            #to_csv输出含有转义字符
            #f_dict=pd.DataFrame(msp_dict[msp_key])
            #df_dict.to_csv(os.path.join(output_path,'trusted_df.msp'), mode='a', header=False,index=False,quoting=csv.QUOTE_NONE,escapechar='\\')
        else:
            print (f"{msp_key} not in msp")
            continue
    with open(os.path.join(output_path, f'{charge}_trusted_df.msp'), 'a') as f:
        for msp_key, msp_value in new_dict.items():
            f.write(msp_value+'\n\n\n')

