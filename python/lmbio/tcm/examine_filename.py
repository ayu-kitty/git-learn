#!/opt/conda/bin/python
import csv
import pandas as pd
import os
import glob
import sys

def extract_filename():
    file_name_pos = [os.path.basename(file).strip(r".raw") for file in glob.glob("pos/*.raw")]
    file_name_neg = [os.path.basename(file).strip(r".raw") for file in glob.glob("neg/*.raw")]
    print(f"正离子文件名：总计{len(file_name_pos)}个")
    # print(file_name_pos)
    print("--"*30)
    print(f"负离子文件名：总计{len(file_name_neg)}个")
    # print(file_name_neg)
    print("--"*30)

    if len(file_name_pos) > 0 or len(file_name_neg) > 0:
        pass
    else:
        print("正负离子文件不存在！")
        sys.exit(1)
    return file_name_pos,file_name_neg

def examin_filename(info_df,file_name):
    mismatch = []
    
    for fn in file_name:
        if fn in info_df["Filename"].values:
            continue
        else:
            mismatch.append(fn)
    if mismatch:
        # with open("mismatch.txt","w",newline = "") as f:
        print(f"找到{len(mismatch)}个文件名与上机表不匹配！")
        print(mismatch)
        # f.write("以下文件在sample_infomation表中未找到对应文件名:")
        # f.write("\n".join(mismatch))
    else:
        print("文件名无误！")
    return mismatch

def examin():
    print("请在pos和neg同级目录下运行此程序！")
    sample_info = "上机信息.xlsx"
    info = pd.read_excel(sample_info,header = 0, sheet_name=0, engine="openpyxl")
    file_name_pos,file_name_neg = extract_filename()
    print("正离子模式：")
    pos = examin_filename(info,file_name_pos)
    print("负离子模式：")
    neg = examin_filename(info,file_name_neg)
    assert len(pos+neg) == 0

if __name__ == "__main__":
    examin()
