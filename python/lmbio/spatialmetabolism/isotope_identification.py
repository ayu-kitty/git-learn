import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import MinMaxScaler
import pandas as pd
import cv2

i = None
cu_peak = None
adj_matrix_weigh = None
adj_matrix = None
return_charge_signal = None
return_signal = None
return_score = None

#转换格式
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

#预处理--极大极小值标准化、高斯模糊、按mz排序
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


#相似度打分
def intensity_thre(data1, data2):
    data1 = np.where(data1 > np.percentile(data1, 95), 1, 0)
    data2 = np.where(data2 > np.percentile(data2, 95), 1, 0)

    data3 = data1 + data2
    score = len(np.where(data3 == 2)[0]) / len(np.where(data3 != 0)[0])
    return score

#鉴定单同位素
def detect_iso(peak_plus, mode='ppm', tolerate=10, type='C', charge=1):

    global i,cu_peak,adj_matrix_weigh,adj_matrix,return_charge_signal,return_signal,return_score
    num = 0
    if mode == 'ppm':
        ind = np.where((peak <= (peak[i] + peak_plus) * (1000000 + tolerate) / 1000000) & (
                        peak >= (peak[i] + peak_plus) * (1000000 - tolerate) / 1000000))[0]
    #
    if mode == 'Da':
        ind = np.where((peak <= (peak[i] + peak_plus) + tolerate) & (peak >= (peak[i] + peak_plus) - tolerate))[0]

    if len(ind) != 0:
        for j in ind:

            weight = intensity_thre(data_filter[:, i], data_filter[:, j])

            if peak[j] in cu_peak and np.sum(data[:, i]) > np.sum(data[:, j]) and weight >= 0.2:

                adj_matrix_weigh[i, j] = weight
                adj_matrix[i,j] = 1
                num = 1
                #从待筛选列表中剔除同位素峰
                cu_peak.remove(peak[j])

                if charge == 1:
                    return_charge_signal[j] = 1
                if charge == 2:
                    return_charge_signal[j] = 2

                if type == 'C':
                    return_signal[j] = 1
                if type == 'S':
                    return_signal[j] = 2
                if type == 'Cl':
                    return_signal[j] = 3
                if type == 'K':
                    return_signal[j] = 4
                if type == 'Br':
                    return_signal[j] = 5
                return_score[j] = weight

    return num


def return_adj(data,data_filter,peak,ion_mode):

    global i,cu_peak,adj_matrix_weigh,adj_matrix,return_charge_signal,return_signal,return_score

    adj_matrix = np.zeros((len(peak), len(peak)))
    adj_matrix_weigh = np.zeros((len(peak), len(peak)))

    cu_peak = list(peak)
    return_signal = np.zeros_like(peak)
    return_charge_signal = np.zeros_like(peak)
    return_score = np.zeros_like(peak)

    for i in range(len(adj_matrix)):

        for charge in [1, 2]:
            if charge == 1:
                if peak[i] in cu_peak:
                    #C鉴别 M+
                    num = detect_iso(peak_plus=1.00336, mode='Da', tolerate=0.0015, type='C', charge=charge)
                    if num == 1:
                        #C M+2
                        num2 = detect_iso(peak_plus=2.0067, mode='Da', tolerate=0.0025, type='C', charge=charge)
                        if num2 == 1:
                            #C M+3
                            num3 = detect_iso(peak_plus=3.01008, mode='Da', tolerate=0.01, type='C', charge=charge)
                            if num3 == 1:
                                #C M+4
                                num4 = detect_iso(peak_plus=4.01344, mode='Da', tolerate=0.01, type='C', charge=charge)
                    #S M+            
                    num = detect_iso(peak_plus=1.9958, mode='Da', tolerate=0.001, type='S', charge=charge)
                    if num == 1:
                        num2 = detect_iso(peak_plus=1.9958 * 2, mode='Da', tolerate=0.01, type='S', charge=charge)
                    
                    if ion_mode == 'positive':
                        num = detect_iso(peak_plus=1.99812, mode='Da', tolerate=0.001, type='K', charge=charge)
                        if num == 1:
                            num2 = detect_iso(peak_plus=1.99812 * 2, mode='Da', tolerate=0.01, type='K', charge=charge)
                    
                    if ion_mode == 'negative':
                        #Cl M+
                        num = detect_iso(peak_plus=1.99705, mode='Da', tolerate=0.001, type='Cl', charge=charge)
                        if num == 1:
                            num2 = detect_iso(peak_plus=1.99705 * 2, mode='Da', tolerate=0.01, type='Cl', charge=charge)
                        #Br M+
                        num = detect_iso(peak_plus=1.99795, mode='Da', tolerate=0.001, type='Br', charge=charge)
                        if num == 1:
                            num2 = detect_iso(peak_plus=1.99795 * 2, mode='Da', tolerate=0.01, type='Br', charge=charge)

            if charge == 2:

                if peak[i] in cu_peak:
                    num = detect_iso(peak_plus=0.50168, mode='Da', tolerate=0.0015, type='C', charge=charge)
                    if num == 1:
                        num2 = detect_iso(peak_plus=0.50168 * 2, mode='Da', tolerate=0.0025, type='C', charge=charge)

                        if num2 == 1:
                            num3 = detect_iso(peak_plus=0.50168 * 3, mode='Da', tolerate=0.01, type='C', charge=charge)

                            if num3 == 1:
                                num4 = detect_iso(peak_plus=0.50168 * 4, mode='Da', tolerate=0.01, type='C',
                                                  charge=charge)
        if peak[i] in cu_peak:
            cu_peak.remove(peak[i])

    return adj_matrix,return_signal,return_charge_signal,adj_matrix_weigh,return_score

# 深度遍历
def graph_travel(G,x,marked):
    result=[]
    result.append(x)
    marked[x] = 1
    if sum(G[x,:])==1:
        #如果行只有1，则无同位素峰
        return result
    else:
        for i in range(G.shape[1]):
            #否则找出所有的交叉点
            if G[x,i] != 0 and x!=i and marked[i]==0:
                result.extend(graph_travel(G,i,marked))
    return result

#对adj_matrix按单同位素离子组织
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

def output_results(filename,results,peak,adj_matrix,return_signal,return_charge_signal,return_score):
    adj_iso = adj_matrix

    peak_list = peak
    peak_list = peak_list[np.argsort(peak_list)]
    # peak_list = np.around(peak_list,4)
    max_num = np.max([len(i) for i in results])

    #创建df0...dfn等空列表
    for i in range(max_num):
        globals()['df%d' % i] = []

    #往df0...中传入peak值
    for kk in results:
        for num in range(1, max_num + 1):
            if len(kk) == num:
                for i in range(len(kk)):
                    globals()['df%d' % i].append(peak_list[kk][i])
                for j in range(len(kk), max_num):
                    globals()['df%d' % j].append(0)

    data = np.empty(shape = (len(df0), max_num + (max_num - 1) * 3),dtype=object)

    for i in range(max_num):
        data[:, i] = globals()['df%d' % i]

    #往df_1*max_num..中传入isotope_type
    for i in range(max_num, max_num + (max_num - 1) * 3):
        globals()['df%d' % i] = []

    for i in range(max_num - 1):
        for j in range(len(df0)):
            if globals()['df%d' % (i + 1)][j] != 0:
                signal = return_signal[np.where(peak_list == globals()['df%d' % (i + 1)][j])[0]]

                if signal == 1:
                    globals()['df%d'%(i + max_num)].append('C')
                if signal == 2:
                    globals()['df%d' % (i + max_num)].append('S')
                if signal == 3:
                    globals()['df%d' % (i + max_num)].append('C')
                if signal == 4:
                    globals()['df%d' % (i + max_num)].append('K')
                if signal == 5:
                    globals()['df%d' % (i + max_num)].append('Br')

                # globals()['df%d' % (i + max_num)].append(signal)
            else:

                globals()['df%d' % (i + max_num)].append("")
    
    #往df_2*max_num...中传入charge_signal
    for i in range(max_num - 1):
        for j in range(len(df0)):
            if globals()['df%d' % (i + 1)][j] != 0:
                signal = return_charge_signal[np.where(peak_list == globals()['df%d' % (i + 1)][j])[0]]

                globals()['df%d' % (i + max_num + max_num - 1)].append(signal[0])
            else:

                globals()['df%d' % (i + max_num + max_num - 1)].append(0)
    
    #往df_3*max_num...中传入score
    for i in range(max_num - 1):
        for j in range(len(df0)):
            if globals()['df%d' % (i + 1)][j] != 0:
                score = return_score[np.where(peak_list == globals()['df%d' % (i + 1)][j])[0]]
                #取数组的值;之前默认为单元素数组。
                globals()['df%d' % (i + max_num + (max_num - 1)*2)].append(score[0])
            else:

                globals()['df%d' % (i + max_num + (max_num - 1)*2)].append(0)
    


    #输出到data中，并添加列名
    for i in range(max_num, max_num + (max_num - 1) * 3):
        data[:, i] = globals()['df%d' % i]

    columns = ['monoisotope']
    for i in range(1, max_num):
        columns.append('isotope_%d' % i)

    for i in range(max_num, max_num * 2 - 1):
        columns.append('isotope_type_%d' % (i - max_num + 1))

    for i in range(max_num * 2 - 1, max_num * 3 - 2):
        columns.append('charge_%d' % (i - max_num * 2 + 1 + 1))
    
    for i in range(max_num * 3 - 2, max_num * 4 - 3):
        columns.append('weight_%d' % (i - max_num * 3 + 3))        

    df = pd.DataFrame(data=data, columns=columns)
    df.to_csv(filename, index=False)

def output_moin(results,data,peak,filename_peak,filename_data):
    #根据results 输出单同位素离子
    moin = np.array([k[0] for k in results])

    data_moin = data[:, moin]

    peak = peak[moin]

    np.savetxt(filename_data, data_moin, delimiter=',')
    np.savetxt(filename_peak, peak, delimiter=',')

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
            plt.savefig(f"pic/pic_{graphs}.jpg")
            plt.close()

if __name__ == '__main__':
    
    ion_mode = 'negative'
    infile = 'RatBrain-neg-process.txt'
    outfile = f'{infile.split(".")[0]}_iso.csv'
    filename_mono_data = f'{infile.split(".")[0]}_data.csv'
    filename_mono_peak = f'{infile.split(".")[0]}_peak.csv'
    
    peak,raw_data,y,x = convert_data(infile)
    print("转格式完成....")

    data, data_filter, peak = data_preprocessing(raw_data, peak, y, x)
    print("预处理完成....")
    adj_matrix, return_signal, return_charge_signal, adj_matrix_weigh,return_score = return_adj(data, data_filter, peak,ion_mode)

    results = output_max_ConnectMatrix(adj_matrix)
    # print(results)
    output_results(outfile, results, peak, adj_matrix, return_signal, return_charge_signal, return_score)
    print("鉴定单同位素完成!")
    output_moin(results,data,peak,filename_mono_peak,filename_mono_data)
    # plot_iso(results,peak,raw_data,x,y)
    print("画图完成！")

