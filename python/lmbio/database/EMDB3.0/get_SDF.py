import pandas as pd
from rdkit.Chem import Descriptors
import rdkit
from rdkit import Chem
from rdkit.Chem import  Draw
import numpy as np
import regex as re
import os
import glob
import mysql.connector

#def preprocess(csvfile):
#    try:
#        indf = pd.read_csv(csvfile)
#    except UnicodeDecodeError:
#        print("CSV文件包含非Unicode字符！")
#    indf = csvfile
    # indf = indf.apply(lambda x : x.replace(np.nan, ""),axis = 0)
#    if "smiles" not in indf.columns:
#        indf["smiles"] = indf["inchi"].apply(inchi2smiles)
#    return indf

def check_smiles(smi_list):
    flags = 0
    for smi in smi_list:
        try:
            mol = Chem.MolFromSmiles(smi)
            # Draw.MolToImage(mol)
        except ValueError as e:
            flags += 1
            print(smi,"is a invalid mol")

    return flags

def sdf(df,outdir = './'):
    props = ["_Name","COMMON_NAME","DATABASE_ID","SMILES"]
    # writer.SetProps(props)
    writer = Chem.SDWriter("tmp.sdf")
    for _,row in df.iterrows():
        # print(row)
        mol = Chem.MolFromSmiles(row.smiles)
        if mol:
            mol.SetProp(props[0],str(row["Compound Name"]))
            mol.SetProp(props[1],str(row["Compound Name"])) 
            mol.SetProp(props[2],str(row["inchikey"]))
            mol.SetProp(props[3],str(row["smiles"]))
            writer.write(mol)
    writer.close()

    with open("tmp.sdf","r") as f:
        with open(f"{outdir}/df_unique.sdf", "w") as fw:
            lines = f.readlines()
            for line in lines:
                # re.sub(r"^(>.*)( {2}\(\d\) )",lambda m:m.group(1),line)
                line = re.sub(r"(?<=^>.*)  \(\d+\) ","",line)
                fw.write(line)
    os.system("rm tmp.sdf")
    print(f"所有分子转换sdf完成！保存至{outdir}/df_unique.sdf文件")
    return 

#if __name__ == "__main__":
def save_SDF(infile = './trusted_unique.csv', outdir='./'):
    df_unique = pd.read_csv(infile)
    #去除smiles为空的值
    df_unique.dropna(subset=['smiles'], inplace=True)
    flags = check_smiles(df_unique["smiles"])
    if flags == 0:
        sdf(df_unique,outdir)