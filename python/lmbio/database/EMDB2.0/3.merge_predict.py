#!/opt/conda/bin/python
import pandas as pd
import numpy as np

if __name__ == "__main__":
    #将transferlearning 预测结果与HMDB_Database合并
    data = pd.read_csv("BP_2400个.csv",sep = "\t")
    hmdb = pd.read_excel("HMDB_DATABASE.xlsx")
    data1 = data.merge(hmdb,left_on = "Smiles",right_on = "SMILES",how = "left")
    data2 = data1.iloc[:,[2,3,4,6,5,1]].copy()
    data2.to_excel("hmdb_2400_pred.xlsx",index = False)
