import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.preprocessing import MinMaxScaler
from torchvision import models
import torch
import umap
import igraph as ig
import cv2

i = None
cu_peak = None
adj_matrix_weigh = None
adj_matrix = None
return_charge_signal = None
return_signal = None
return_score = None
peak_true = None

def convert_data(infile):
    data = pd.read_csv(infile,sep = "\t",index_col = "mz")
    coords = np.int32([coord.split("-") for coord in list(data.columns)])
    x = coords[:,0].max()
    y = coords[:,1].max()

    peak = np.array(data.index)
    data = data.values.T

    print(f"图片大小{x}*{y}")
    # np.savetxt("peak_list",peak)
    # np.savetxt("data_brain",data)
    # print(data.shape)
    return peak,data,y,x


def data_preprocessing(data,peak,mm,nn):
    # step_1 Scale

    data = MinMaxScaler().fit_transform(data)

    # step_2 Produce Filtered Data

    data = data.reshape(mm, nn, -1)
    data_filter = np.zeros_like(data)


    for i in range(data.shape[2]):
        data_filter[:, :, i] = cv2.blur(data[:, :, i], (5, 5))

    data_filter = data_filter.reshape(mm * nn, -1)
    data = data.reshape(mm * nn, -1)
    data = data[:, np.argsort(peak)]
    data_filter = data_filter[:, np.argsort(peak)]
    peak = peak[np.argsort(peak)]

    return data,data_filter,peak

def intensity_thre(data1, data2):
    data1 = np.where(data1 > np.percentile(data1, 95), 1, 0)
    data2 = np.where(data2 > np.percentile(data2, 95), 1, 0)

    data3 = data1 + data2
    score = len(np.where(data3 == 2)[0]) / len(np.where(data3 != 0)[0])
    return score

# 深度遍历
def graph_travel(G,x,marked):
    result=[]
    result.append(x)
    marked[x] = 1
    if sum(G[x,:])==1:
        return result
    else:
        for i in range(G.shape[1]):
            if G[x,i] != 0 and x!=i and marked[i]==0:
                result.extend(graph_travel(G,i,marked))
    return result


def detect_iso1(peak_plus,ion_mode, mode='ppm', tolerate=10, types = 1, charge=1, num_M = 1):
    global i,cu_peak,adj_matrix_weigh,adj_matrix,return_charge_signal,return_signal,peak_true
    #根据peak_true筛选各种加合离子模式，对应的家合离子
    numm = 0
    if ion_mode == 'positive':
        if mode == 'ppm':
            find_peak = (peak_true * num_M + peak_plus) / charge
            ind = np.where((peak <= find_peak * (1000000 + tolerate) / 1000000) & (
                    peak >= find_peak * (1000000 - tolerate) / 1000000))[0]

        if mode == 'Da':
            find_peak = (peak_true * num_M + peak_plus) / charge
            ind = np.where((peak <= find_peak + tolerate) & (peak >= find_peak - tolerate))[0]

    if ion_mode == 'negative':
        if mode == 'ppm':
            find_peak = (peak_true * num_M - peak_plus) / charge

            ind = np.where((peak <= find_peak * (1000000 + tolerate) / 1000000) & (
                    peak >= find_peak * (1000000 - tolerate) / 1000000))[0]

        if mode == 'Da':
            find_peak = (peak_true * num_M - peak_plus) / charge
            ind = np.where((peak <= find_peak + tolerate) & (peak >= find_peak - tolerate))[0]

    if len(ind) != 0:

        for j in ind:
            # print(peak[j])
            # print('type: '+ str(types))
            weight = intensity_thre(data_filter[:, i], data_filter[:, j])

            if peak[j] in cu_peak and weight >= 0.2 and j != i:
                numm = 1

    return numm

def detect_iso2(peak_plus,ion_mode, mode='ppm', tolerate=10, types = 1, charge=1, num_M = 1):
    global i,cu_peak,adj_matrix_weigh,adj_matrix,return_charge_signal,return_signal,peak_true,return_score
    #根据筛选出来的加合离子模式（peak_true)，挑选同组的加合离子
    numm = 0
    if ion_mode == 'positive':
        if mode == 'ppm':
            find_peak = (peak_true * num_M + peak_plus) / charge

            ind = np.where((peak <= find_peak * (1000000 + tolerate) / 1000000) & (
                    peak >= find_peak * (1000000 - tolerate) / 1000000))[0]
        if mode == 'Da':
            find_peak = (peak_true * num_M + peak_plus) / charge

            ind = np.where((peak <= find_peak + tolerate) & (peak >= find_peak - tolerate))[0]
    if ion_mode == 'negative':
        if mode == 'ppm':
            find_peak = (peak_true * num_M - peak_plus) / charge

            ind = np.where((peak <= find_peak * (1000000 + tolerate) / 1000000) & (
                    peak >= find_peak * (1000000 - tolerate) / 1000000))[0]
        if mode == 'Da':
            find_peak = (peak_true * num_M - peak_plus) / charge
            ind = np.where((peak <= find_peak + tolerate) & (peak >= find_peak - tolerate))[0]

    #从mz在误差范围内的离子中归为一组
    if len(ind) != 0:
        for j in ind:
            weight = intensity_thre(data_filter[:, i], data_filter[:, j])

            if peak[j] in cu_peak and weight >= 0.2 and j != i:
                adj_matrix_weigh[i, j] = weight
                adj_matrix[i,j] = 1
                numm = 1
                cu_peak.remove(peak[j])
                return_signal[j] = types
                return_score[j] = weight
    return numm


def return_adj(data, data_filter, peak,addunt_mz,ion_mode):

    global i,cu_peak,adj_matrix_weigh,adj_matrix,return_charge_signal,return_signal,peak_true,return_score

    charges = addunt_mz['charge'].values
    mass_shift = addunt_mz['m/z'].values
    num_Ms = addunt_mz['num'].values
    signal = addunt_mz['Adduct'].values

    adj_matrix = np.zeros((len(peak), len(peak)))
    adj_matrix_weigh = np.zeros((len(peak), len(peak)))

    cu_peak = list(peak)
    return_signal = np.empty(shape = peak.shape,dtype = object)
    return_score = np.zeros_like(peak)

    for i in range(len(peak)):

        if peak[i] in cu_peak:

            numss = []
            numss_id = []
            numss_id2 = []
            for ii in range(len(num_Ms)):
                if ion_mode == 'positive':
                    # peak_true = peak[i] - mass_shift[ii]
                    peak_true = (peak[i] * charges[ii] - mass_shift[ii]) / num_Ms[ii]
                if ion_mode == 'negative':
                    #正负离子模式都改为 -
                    peak_true = (peak[i] * charges[ii] + mass_shift[ii]) / num_Ms[ii]

                #根据peak_true对各种可能加合离子模式进行检测
                nums = 0
                for kk in range(len(num_Ms)):
                    if kk != ii:
                        numm = detect_iso1(peak_plus=mass_shift[kk],ion_mode = ion_mode, mode='ppm', tolerate=5, types=signal[kk],
                                           charge=charges[kk], num_M=num_Ms[kk])
                        nums += numm

                if nums != 0:
                    numss_id.append(ii)
                    numss.append(nums)

            if len(numss_id) != 0:
                #检测的加合离子数最多的离子模式(peak_true)作为结果
                sele = numss_id[np.argmax(np.array(numss))]
                if ion_mode == 'positive':
                    # peak_true = peak[i] - mass_shift[sele]
                    peak_true = (peak[i] * charges[sele] - mass_shift[sele]) / num_Ms[sele]

                if ion_mode == 'negative':
                    #正负离子模式都改为 -
                    peak_true = (peak[i] * charges[sele] + mass_shift[sele]) / num_Ms[sele]
                nums = 0
                for kk in range(len(num_Ms)):
                    if kk != sele:
                        numm = detect_iso2(peak_plus=mass_shift[kk],ion_mode = ion_mode, mode='ppm', tolerate=5, types=signal[kk],
                                           charge=charges[kk], num_M=num_Ms[kk])
                        nums += numm
                #adduct_type0的类型
                return_signal[i] = signal[sele]
                # return_score[i] = 

        if peak[i] in cu_peak:
            cu_peak.remove(peak[i])

    return adj_matrix, return_signal, return_charge_signal, adj_matrix_weigh, return_score

def output_max_ConnectMatrix(adj_matrix):

    for i in range(len(adj_matrix)):
        for j in range(len(adj_matrix)):
            if adj_matrix[i, j] == 1:
                adj_matrix[j, i] = 1
    
    le = adj_matrix.shape[0]
    for i in range(le):
        adj_matrix[i, i] = 1
    # print(G[0])

    marked = [0] * le
    # print(marked)

    results = []
    for i in range(le):
        if marked[i] == 0:
            result = graph_travel(adj_matrix, i,marked)
            # print(result)
            results.append(result)
    return results

def output_results(filename,results,peak,adj_matrix,return_signal,adj_matrix_weigh,return_score):

    adj_iso = adj_matrix
    peak_list = peak
    peak_list = peak_list[np.argsort(peak_list)]

    peak_list = np.around(peak_list, 4)
    max_num = np.max([len(i) for i in results])
    print("最大加合离子数为:",max_num)
    #df0...中存入peak
    for i in range(max_num):
        globals()['df%d' % i] = []

    for kk in results:
        for num in range(1, max_num + 1):
            if len(kk) == num:
                for i in range(len(kk)):
                    globals()['df%d' % i].append(peak_list[kk][i])
                for j in range(len(kk), max_num):
                    # print(peak_list[kk])
                    # print(peak_list[kk][j])
                    globals()['df%d' % j].append(0)

    data = np.empty(shape = (len(df0), max_num * 3-1),dtype = object)

    for i in range(max_num):
        data[:, i] = globals()['df%d' % i]

    #df_1*max_num...输出adduct类型
    for i in range(max_num, max_num + max_num):
        globals()['df%d' % i] = []

    for i in range(max_num):
        for j in range(len(df0)):

            if globals()['df%d' % (i)][j] != 0:
                signal = return_signal[np.where(peak_list == globals()['df%d' % (i)][j])[0]]

                globals()['df%d' % (i + max_num)].append(signal[0])
            else:

                globals()['df%d' % (i + max_num)].append("")

    #df_2*max_num...输出score
    for i in range(max_num*2, max_num*3-1):
        globals()['df%d' % i] = []

    for i in range(1,max_num):
        for j in range(len(df0)):

            if globals()['df%d' % (i)][j] != 0:
                score = return_score[np.where(peak_list == globals()['df%d' % (i)][j])[0]]

                globals()['df%d' % (i + max_num*2 -1)].append(score[0])
            else:

                globals()['df%d' % (i + max_num*2 -1)].append(0)


    for i in range(max_num, max_num *3 -1):
        data[:, i] = globals()['df%d' % i]

    columns = ['monoisotope']
    for i in range(1, max_num):
        columns.append('add_%d' % i)

    for i in range(max_num, max_num * 2):
        columns.append('add_type_%d' % (i - max_num))
    
    for i in range(max_num*2, max_num * 3 - 1):
        columns.append('score_%d' % (i - max_num*2+1))

    df = pd.DataFrame(data=data, columns=columns)
    df.to_csv(filename, index=False)

    return

def plot_iso(sub_graph,peak,data,x,y):
    to_pd = []
    peak = np.around(peak,4)
    for graphs in sub_graph:
        if len(graphs) > 1:
            graphs = [int(i) for i in graphs]
            # print(peak[graphs])
            to_pd.append(peak[graphs])
            num = 1
            for graph in graphs:
                plt.subplot(1,len(graphs),num)
                data2 = data[:, np.where(peak == peak[graph])[0]]
                plt.imshow(data2.reshape(y, x))
                plt.title(str(peak[graph]) + '\n' + str(graphs))
                num += 1
            plt.savefig(f"pic_/pic_{graphs}.jpg")
            plt.close()


if __name__ == '__main__':

    ion_mode = 'negative'
    infile = 'RatBrain-neg-process.txt'
    peak,data,y,x = convert_data(infile)

    raw_data = np.loadtxt(f'{infile.split(".")[0]}_data.csv',delimiter=',')
    peak = np.loadtxt(f'{infile.split(".")[0]}_peak.csv',delimiter=',')

    output_filename = f'{infile.split(".")[0]}_adduct.csv'
    addunt_mz = pd.read_excel(f'{ion_mode}.xlsx')

    data, data_filter, peak = data_preprocessing(raw_data, peak, y, x)
    print("读取单同位素峰完成....")
    adj_matrix, return_signal, return_charge_signal, adj_matrix_weigh,return_score = return_adj(data, data_filter, peak,addunt_mz,ion_mode)

    results = output_max_ConnectMatrix(adj_matrix)
    # print(results)
    output_results(output_filename, results, peak, adj_matrix, return_signal,  adj_matrix_weigh, return_score)
    print("加合离子鉴定完成!")
    # plot_iso(results,peak,raw_data,x,y)
    print("画图完成！")
