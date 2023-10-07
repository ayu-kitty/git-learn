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

#合并后的文件与预测RT进行合并筛选
def merge_hm(q,h):	
	#保留RT与predict相差在0.5以内的行
	qh = q.merge(h,left_on = "Compound ID",right_on = "hmdbid",how = "left")
	qh["diff"] = qh.RT - qh.Predict
	qh_fl = qh[(-0.5 <= qh["diff"]) & (qh["diff"] <=0.5)]
	print("RT与prdeict 合并筛选完成")
	#根据CID分组，对同组的RT求均值
	data = []
	gpby = qh_fl.groupby("Compound ID")
	for _,d in gpby:
	    mrt = d.mean(numeric_only = True).RT
	    d["RT"] = mrt
	    data.append(d)

	res = pd.concat(data).drop_duplicates("Compound ID").loc[:,["Compound ID",'Formula',
		 'InChIKey', 'InChIIdentifier', 'SMILES','RT','Predict','Source']]
	res.to_excel("Curated_HMDB_2400.xlsx",index = False)
	print("已保存至Curated_HMDB文件")
	return res

def acc(hm,bpfile,thresh):
	bp = pd.read_csv(bpfile,sep = "\t")
	merged = hm.merge(bp,left_on = "InChIKey",right_on = "inchikey",how = "left")

	rts = merged.loc[:,["Compound ID","InChIKey","RT_x","RT_y"]]
	rt = rts[pd.notnull(rts.RT_y)]

	def percent(T,Y,thresh):
	    err = np.abs(T - Y)
	    pct = len(err[err < thresh]) / len(T)
	    return pct

	return percent(rt.RT_y,rt.RT_x,thresh)

if __name__ == "__main__":
	files = glob.glob("*.csv")
	bpfile = "标品total_2419_inchikey.txt"
	thresh = 0.5

	hmdb = pd.read_excel("hmdb_2400_pred.xlsx")
	# qi = concate(files)
	qi = pd.read_excel("Total_cmpd_uni.xlsx")
	hm = merge_hm(qi,hmdb)
	acc = acc(hm,bpfile,thresh)
	print(acc)
