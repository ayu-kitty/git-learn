import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import MinMaxScaler
import os
import cv2
from sklearn.cluster import KMeans
from sklearn.cluster import AgglomerativeClustering
import os
import sys

#剔除噪音离子(表达占比低，无模式)
def filter_peak(data,threhold):
    #对离子mz去重
    num = int(np.max(data[:, 0]))
    print(f"图像大小为{num}像素...",)
    mz_list = data[:, 1]
    mz_list = np.array(list(set(mz_list)))
    mz_list = np.array([float(i) for i in mz_list])
    mz_list = mz_list[np.argsort(mz_list)]

    return_mat = np.zeros((len(mz_list), num))
    
    #按照离子进行排列,离子的强度矩阵 -->unique_mz * pixel(34*5500)
    for k in range(num):
        mz = data[:, 1][np.where(data[:, 0] == (k + 1))[0]]
        intensity = data[:, 2][np.where(data[:, 0] == (k + 1))[0]]

        for mm in range(len(mz)):
            inn = np.where(mz_list == float(mz[mm]))[0]
            return_mat[inn, k] = intensity[mm]

    #筛选缺失值坐标占整个图像比例>5%的离子
    mask = np.where(return_mat == 0, 0, 1)
    index = np.sum(mask, axis=1)
    return_mat = return_mat[np.where(index / num > threhold)[0]]

    return return_mat,mz_list[np.where(index / num > threhold)[0]]

#画离子呈像图
def plot_ion(idata,peak,suf):
    data = idata.T.reshape(mm,nn,-1)
    print("num of pic:",data.shape[2])
    for i in range(data.shape[2]):
        plt.imshow(data[:,:,i])
        plt.colorbar(shrink=  0.5)
        plt.title("%.4f" %peak[i])
        plt.savefig(suf+"%.4f" %peak[i]+".jpg")
        plt.close()

#kmeans聚类区分背景和样品
def get_mask(data,mm,nn,n_clusters=2,method="kmeans"):
    #input:[pixel,ion]
    if method == "kmeans":
        labels = KMeans(n_clusters=n_clusters,random_state=0).fit_predict(data)
        
    elif method == "hirachical":
        linkage = "ward"  #OR "average", "complete", "single"):
        clustering = AgglomerativeClustering(linkage=linkage, n_clusters=n_clusters)
        clustering.fit(data)
        labels = clustering.labels_
    n = len(list(set(labels)))
    #画图
    for i in range(n):
        plt.subplot(1,n,i+1)
        to_img = np.where(labels == i, 1, 0)
        plt.imshow(to_img.reshape(mm,nn),cmap = "gray")
        plt.colorbar(shrink=0.3)
        plt.title(i)
    plt.show()
    
    tis_bgr = int(input("请选择背景(0)还是组织(1)："))
    sample_num = int(input("请根据图像输入样本对应编号："))
    if tis_bgr == 1:
        mask = np.where(labels == sample_num, 1, 0)
    elif  tis_bgr == 0:
        mask = np.where(labels == sample_num, 0, 1)
    return mask


def remvoe_noisy(data,peak,thre,n_clusters):
    #剔除背景 input: [pixel,ion],
    #mask: [pixel,1]
    mask = get_mask(data,mm,nn,n_clusters)
    mask2 = np.where(mask == 1, 0, 1)
    mask2 = mask2.reshape(-1,1)
    mask = mask.reshape(-1,1)
    m,n = data.shape

    returnlist = []
    for i in range(n):
        #分离样品组织和背景，组织区域非0数量/背景区域非0数量 > 1.5
        tissue = data[:,i].reshape(-1,1) * mask
        backgr = data[:,i].reshape(-1,1) * mask2
        one_tissue = len(np.where(tissue != 0)[0])/np.sum(mask)
        one_backgr = len(np.where(backgr != 0)[0])/np.sum(mask2)
        if one_tissue / one_backgr > thre:
            returnlist.append(i)

    returnlist = np.array(returnlist)

    return np.transpose(data)[returnlist],peak[returnlist]

def intensity_thre(data1, data2):
    data1 = np.where(data1 > np.percentile(data1, 95), 1, 0)
    data2 = np.where(data2 > np.percentile(data2, 95), 1, 0)

    data3 = data1 + data2
    score = len(np.where(data3 == 2)[0]) / len(np.where(data3 != 0)[0])
    return score

#根据离子强度相似性，对二级离子
def predict_ion(peak,data,mm,nn,mother_mzz,tolerate,score_threhold):
    #标准化 data:[ion, pixel]
    data = np.transpose(data)
    data = MinMaxScaler().fit_transform(data)

    #Gaussian模糊化
    data = data.reshape(mm, nn, -1)
    data_filter = np.zeros_like(data)
    for i in range(data.shape[2]):
        data_filter[:, :, i] = cv2.blur(data[:, :, i], (5, 5))

    data_filter = data_filter.reshape(mm * nn, -1)
    data = data.reshape(mm * nn, -1)

    #母离子范围 & 筛选
    mather_ions = peak[np.where((peak <= mother_mzz + tolerate) & (peak >= mother_mzz - tolerate))[0]]

    output_dict = {}

    for mother_ion in mather_ions:
        son_ion = []

        for i in peak:
            if i != mother_ion:
                score = intensity_thre(data[:, np.where(peak == i)[0]], data[:, np.where(peak == mother_ion)[0]])
                if score > score_threhold:
                    son_ion.append(i)

        output_dict['%s' % str(mother_ion)] = son_ion
    return output_dict

#将输出的二级离子列表等长
def equallen(l,ln):
    for i in range(ln-len(l)):
        l.append(np.nan)
    return l

if __name__ == '__main__':
    mother_mzz = 146.1649
    # mother_mzz = 203.2227
    # mother_mzz = 307.0433
    # mother_mzz = 798.539
    # mother_mzz = 826.5705
    mm, nn = 50, 111

    infile = './%.4f_bin.csv' % mother_mzz
    outfile = '%.4f_ms2.xlsx' % mother_mzz

    filter_threhold = 0.05
    noisy_threhold = 1.5
    motherion_threhold = 0.75
    score_threhold = 0.1
    n_clusters = 3

    data = pd.read_csv(infile, delimiter=',')
    data = data.values

    data_filter, peak_filter = filter_peak(data, filter_threhold)
    os.makedirs(f"{mother_mzz}_pic_flt",exist_ok = True)
    plot_ion(data_filter,peak_filter,suf = f"{mother_mzz}_pic_flt/")
    print("离子峰过滤完成...")

    data_remove, peak_remove = remvoe_noisy(data_filter.transpose(), peak_filter, noisy_threhold, n_clusters)
    os.makedirs(f"{mother_mzz}_pic_rm",exist_ok = True)
    plot_ion(data_filter,peak_filter,suf = f"{mother_mzz}_pic_rm/")
    print("去除背景离子完成...")
    
    out_dict = predict_ion(peak_remove, data_remove, mm, nn, mother_mzz, motherion_threhold, score_threhold)
    print("二级离子归属完成...")
    if len(out_dict ) > 0:
        maxlen = max([len(value) for value in out_dict.values()])
        out_df = pd.DataFrame({pk1:equallen(pk2,maxlen) for pk1,pk2 in out_dict.items()})
        out_df.to_excel(outfile,index =False)
    else:
        print("未找到二级离子...")
    print("文件已保存")
