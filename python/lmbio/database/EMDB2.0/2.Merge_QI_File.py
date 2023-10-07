#!/opt/conda/bin/python
import numpy as np
import pandas as pd
import glob
import regex as re

def filter_score(qifile,thresh = 50):
	df = pd.read_csv(qifile)
    df["Source"] = qifile.strip(".csv")
    print(qifile.strip(".csv"))
    df.rename({"Retention time (min)":"RT"},axis = 1,inplace = True)
    res = df.loc[:,["Compound ID","Score","Fragmentation Score","RT","Source"]]

    res_hs = res[res.Score >= thresh]
    return res_hs

def concate(files):
	#合并各个文件的化合物
	cmpds = []
	for file in files:
	    cmpds.append(filter_score(file))

	all_cmpds = pd.concat(cmpds,axis = 0,ignore_index = True)

	#根据ID筛选出匹配到各个数据库中的物质
	hmdb = all_cmpds[all_cmpds["Compound ID"].apply(lambda x : True if re.search("^HMDB",x) else False)]
	# metlin = all_cmpds[all_cmpds["Compound ID"].apply(lambda x : True if re.search(r"^\d+$",x) else False)]
 	# lmaps = all_cmpds[all_cmpds["Compound ID"].apply(lambda x : True if re.search(r"^LM",x) else False)]
	# hml = pd.concat([hmdb,metlin,lmaps])
	hmdb_sorted = hmdb.sort_values(by = "Compound ID",axis = "index")
	# hmdb_sorted.set_index("Compound ID",inplace=True)
	hmdb_sorted.to_excel("Total_cmpd.xlsx",index = False)
	return hmdb_sorted

if __name__ == "__main__":
	files = glob.glob("*.csv")
	concate(files)
