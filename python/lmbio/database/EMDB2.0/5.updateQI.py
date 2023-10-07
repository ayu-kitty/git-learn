#!/opt/conda/bin/python
import pandas as pd
import numpy as np

if __name__ == "__main__":
    #准备标品1390文件
    bp = pd.read_csv("标品1390.csv")
    bp["Compound ID"].fillna(bp.Identifiers,inplace = True)

    bp["source"] = ''
    bp["Predict"] = np.nan
    cols = ["Compound ID","Formula","smiles",'inchikey','RT','Predict_RT','source']

    bp1 = bp.loc[:,["Compound ID","Formula_x","smiles","inchiKey","Retention time (min)","Predict","source"]].copy()
    bp1.columns = cols

    #准备Curated_HMDB文件
    hmdb = pd.read_excel("Curated_HMDB_1300.xlsx")
    hmdb =hmdb.loc[:,["Compound ID","Formula","SMILES","InChIKey","RT","Predict","Source"]].copy()
    hmdb.columns = cols

    hmdb.set_index("inchikey",inplace = True)
    bp1.set_index("inchikey",inplace = True)

    #合并BP和HMDB库
    mgd = hmdb.merge(bp1,left_index = True,right_index = True,suffixes = ("_hmdb","_bp"),how = "outer")
    mgd["Compound ID_hmdb"].fillna(mgd["Compound ID_bp"],inplace = True)
    mgd["Formula_hmdb"].fillna(mgd["Formula_bp"],inplace = True)
    mgd["smiles_hmdb"].fillna(mgd["smiles_bp"],inplace = True)
    # mgd.reset_index().to_excel("EMDB_8657.xlsx",index = False)

    #将标品对应的保留时间和source置空，并用BP填充
    mgd.loc[pd.notnull(mgd.RT_bp),"RT_hmdb"] = np.nan
    mgd.loc[pd.notnull(mgd.RT_bp),"source_hmdb"] = ""
    mgd.RT_hmdb.fillna(mgd.RT_bp,inplace = True)
    #选取前六列
    col = ["Compound ID","Formula","smiles",'RT','Predict_RT','source']
    mgd_hmdb = mgd.iloc[:,0:6]
    mgd_hmdb.columns = col

    mgd_hmdb.reset_index().to_excel("EMDB_database.xlsx",index = False)
