#!/opt/conda/bin/python
import numpy as np
from numba import jit
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import argparse
import glob
import os
import warnings
from multiprocessing import Pool
import time
warnings.filterwarnings("ignore")

def forward_trans(img):
    f = np.fft.fft2(img)
    #将高频移到中心
    fshift = np.fft.fftshift(f)
    return fshift

#缩放到0~255
def scale(fshift):
    image = np.log1p(np.abs(fshift))
    max_ = image.max()
    min_ = image.min()
    image = np.uint8(255*(image-min_)/(max_-min_))
    return image

#根据中心区域表达量与整行均值比较，判断是否有横纹    
# def clf(amp,hrad = 2,vrad = 5,thresh=0.1):
#     center_h = amp.shape[0] // 2 
#     center_w = amp.shape[1] // 2 
#     high = amp[0:center_h-vrad,center_w-hrad:center_w+hrad].mean()
#     avg = amp[0:center_h-vrad,0:].mean()
#     # print("差异值:",(high-avg)/avg)
#     if (high-avg)/avg > thresh:
#         return 1
def clf(amp,hrad = 1,vrad = 5,step = 10,thresh1=0.2,thresh2 = 0.001):
    center_h = amp.shape[0] // 2 
    center_w = amp.shape[1] // 2 
    highs = []
    avgsv = []
    avgsw = []
    #纵向比较
    for v in np.arange(0,center_h-vrad-step):
        highs.append(amp[v:v+step,center_w-hrad:center_w+hrad].mean())
        avgsv.append(amp[v:v+step,:].mean())
    ratio = [(highs[i]-avgsv[i])/avgsv[i] for i in range(len(highs))]
    #横向比较
    avg = amp[0:center_h-vrad,center_w-hrad:center_w+hrad].mean()
    for w in np.arange(hrad,center_w):
        avgsw.append(amp[0:center_h-vrad,center_w-w:center_w+w].mean())
        
    if len(ratio) == 0:
      coef1 = 0
    else:
      coef1 = round(max(ratio),4)
    coef2 = round((max(avgsw)-avg)/avg,4)
    if coef1 > thresh1 and coef2 < thresh2:
        return 1,coef1,coef2
    else:
        return 0,coef1,coef2

#矩形低通滤波
def lpfilter(fshift,hrad=1,vrad=5,hs = True):
    #hrad,水平宽度;vrad,垂直高度
    f_fshift = fshift.copy()
    rows,cols = fshift.shape
    crow,ccol = rows//2,cols//2   
    if hs:
        #横纹
        f_fshift[:crow-vrad,ccol-hrad:ccol+hrad] = 0
        f_fshift[crow+vrad:,ccol-hrad:ccol+hrad] = 0
    else:
        #纵纹
        f_fshift[crow-vrad:crow+vrad,ccol+hrad:] = 0
        f_fshift[crow-vrad:crow+vrad,:ccol-hrad] = 0

    return f_fshift
  
def lpfilter1(fshift,hrad=1,vrad=1, alpha=0.3):
    
    f_fshift = fshift.copy()
    rows,cols = fshift.shape
    crow,ccol = rows//2,cols//2
    if crow>hrad and ccol>vrad:
        f_fshift[:crow-vrad,ccol-hrad:ccol+hrad] = alpha*f_fshift[:crow-vrad,ccol-hrad:ccol+hrad]
        f_fshift[crow+vrad:,ccol-hrad:ccol+hrad] = alpha*f_fshift[crow+vrad:,ccol-hrad:ccol+hrad]
        f_fshift[crow-vrad:crow+vrad,ccol+hrad:] = (1-alpha)*f_fshift[crow-vrad:crow+vrad,ccol+hrad:]
        f_fshift[crow-vrad:crow+vrad,:ccol-hrad] = (1-alpha)*f_fshift[crow-vrad:crow+vrad,:ccol-hrad]
    else:
        print("请减少hrad、vrad的值，使其小于x，y的一半")
        
    return f_fshift

@jit        
def gaussianFreqFilter(fshift, sigma, centX=None,centY=None):
    
    # 构造高斯核
    width,height = fshift.shape
    if centX is None and centY is None:
        centX = int(height/2)-1
        centY = int(width/2)-1
    Gauss_map1 = np.zeros((width,height))
    centX1 = centX
    centY1 = int(width/2)
    for i in range(width):
        for j in range(height):
            dis = np.sqrt((i - centY1) ** 2 + (j - centX1) ** 2)
            Gauss_map1[i, j] =1.0 -np.exp(-0.5 * dis / sigma)  # 1.0- 表明是高通滤波器

    Gauss_map2 =np.zeros((width,height))
    centX2 = height - centX1
    centY2 = int(width/2)
    for i in range(width):
        for j in range(height):
            dis =np.sqrt((i - centY2) ** 2 + (j - centX2) ** 2)
            Gauss_map2[i, j] =1.0 - np.exp(-0.5 * dis / sigma)

    Gauss_map3 = np.zeros((width,height))
    centX3 = int(height/2)
    centY3 = centY
    for i in range(width):
        for j in range(height):
            dis = np.sqrt((i - centY3) ** 2 + (j - centX3) ** 2)
            Gauss_map3[i, j] = 1.0 - np.exp(-0.5 * dis / sigma)

    Gauss_map4 = np.zeros((width,height))
    centX4 = int(height/2)
    centY4 = width - centY3

    for i in range(width):
        for j in range(height):
            dis = np.sqrt((i - centY4) ** 2 + (j - centX4) ** 2)
            Gauss_map4[i, j] = 1.0 - np.exp(-0.5 * dis / sigma)

    blur_fft = fshift*Gauss_map1*Gauss_map2*Gauss_map3*Gauss_map4

    return blur_fft
        
#移回
def backward_trans(fshift):
    ishift = np.fft.ifftshift(fshift)
    iimg = np.fft.ifft2(ishift)
    return np.abs(iimg)

#画图
def fixstrip_plot(image,mz,path,cmap, coef = False,coef1 = "",coef2 = "",ticks = True):
    if not os.path.exists(path):  # add
      os.makedirs(path)           # add
    fig = plt.figure(figsize = (5,5),dpi = 100)
    plt.imshow(image,aspect = "equal",cmap = cmap)
    if not ticks:        
        plt.xticks([])
        plt.yticks([])
    if coef:
        s = f'{coef1}\n{coef2}'
        plt.text(0.1,0.1,s,transform = fig.transFigure,fontsize = 10)
    plt.colorbar(shrink = 0.5)
    plt.title(mz)
    plt.savefig(f"{path}/{mz}.jpg")
    plt.close()


#调用低通滤波
def fixstrip_data_v1(data,pplot = False,plotpath = "",hrad = 1,vrad = 1, sigma=5, alpha=0.1, method=1, **kwargs):

    #根据colunm_name计算x,y值
    coords = np.int32([coord.split("-") for coord in list(data.columns)])
    x = coords[:,0].max()
    y = coords[:,1].max()

    print(f"读取图片大小为{x}*{y}！")
    print(f"共有{len(data)}个离子！")
    clist = ["black","blue","green","yellow","red"]
    cmp = LinearSegmentedColormap.from_list("sm",clist )
    #定义输出的dataFrame
    newdf = pd.DataFrame(np.zeros(data.shape,dtype =np.float64),index = data.index,columns = data.columns)
    fcs = []
    counter = 0
    for idx,row in data.iterrows():

        img = row.values.reshape(y,x)

        #傅里叶变换
        fshift = forward_trans(img)
        amp = scale(fshift)
        if pplot:
            fixstrip_plot(amp,mz = idx,path = plotpath+"freq",cmap = "gray",ticks = False)  #plot +str(sigma)+str(alpha)

        idf,coef1,coef2 = clf(amp,hrad = hrad,vrad = vrad)
        #显示原图
        if pplot:
            fixstrip_plot(img,mz = idx, coef = True, coef1 = coef1 ,coef2 = coef2,path = plotpath+"origin",cmap = cmp) # +str(sigma)+str(alpha)

        if idf and not np.isnan(coef1):
            #滤波;频率空间图-- 无ticks，灰度图。
            if method==1:
                print(f'using method: {method}(lpfilter)')
                f_fshift = lpfilter(fshift,hrad = hrad,vrad = vrad)
            elif method==2:
                print(f'using method: {method}(lpfilter1)')
                f_fshift = lpfilter1(fshift,hrad = hrad,vrad = vrad, alpha=alpha)
            else:
                print(f'using method: {method}(gaussianFreqFilter)')
                f_fshift = gaussianFreqFilter(fshift, sigma=sigma, centX=None,centY=None)
                
            counter += 1

            if pplot:
                fixstrip_plot(scale(f_fshift),mz = str(idx)+"_trim",coef = True,coef1 = coef1,coef2 = coef2,path = plotpath+"freq",cmap = "gray",ticks = False)
        else:
            f_fshift = fshift
            
        #转回空域,根据前后均值的变化倍数，回调转换后的值.
        iimg = backward_trans(f_fshift) if not np.isnan(coef1) else img
        fcs.append(iimg.mean()/img.mean())
        newdf.loc[idx,] = iimg.flatten()
        if pplot:            
            fixstrip_plot(iimg,mz = idx,coef = True,coef1 = coef1,coef2 = coef2,path = plotpath+"adjusted",cmap = cmp)
    fcs_arr = np.array(fcs)
    fcs_mean = np.round(fcs_arr[np.where(~np.isinf(fcs_arr) & ~(np.isnan(fcs_arr)))].mean(),4)
    
    # np.savetxt("fcs.txt",fcs_arr)
    print(f"转换完毕！共转换{counter}个离子！")
    print(f"转换后相比转换前倍数变化为{fcs_mean}")

    return newdf
  
def fixstrip_single(idx, row, x, y, method, hrad, vrad, cmp, pplot, plotpath, alpha, sigma, **kwargs):
        
        img = row.values.reshape(y,x)

        #傅里叶变换
        fshift = forward_trans(img)
        amp = scale(fshift)
        if pplot:
            fixstrip_plot(amp,mz = idx,path = plotpath+"freq",cmap = "gray",ticks = False)  #plot +str(sigma)+str(alpha)

        idf,coef1,coef2 = clf(amp,hrad = hrad,vrad = vrad)
        #显示原图
        if pplot:
            fixstrip_plot(img,mz = idx, coef = True, coef1 = coef1 ,coef2 = coef2,path = plotpath+"origin",cmap = cmp) # +str(sigma)+str(alpha)

        if idf and not np.isnan(coef1):
            #滤波;频率空间图-- 无ticks，灰度图。
            if method==1:
                f_fshift = lpfilter(fshift,hrad = hrad,vrad = vrad)
            elif method==2:
                f_fshift = lpfilter1(fshift,hrad = hrad,vrad = vrad, alpha=alpha)
            else:
                f_fshift = gaussianFreqFilter(fshift, sigma=sigma, centX=None,centY=None)
            count = 1   
            if pplot:
                fixstrip_plot(scale(f_fshift),mz = str(idx)+"_trim",coef = True,coef1 = coef1,coef2 = coef2,path = plotpath+"freq",cmap = "gray",ticks = False)
        else:
            count = 0
            f_fshift = fshift
            
        #转回空域,根据前后均值的变化倍数，回调转换后的值.
        iimg = backward_trans(f_fshift) if not np.isnan(coef1) else img
        if pplot:            
            fixstrip_plot(iimg,mz = idx,coef = True,coef1 = coef1,coef2 = coef2,path = plotpath+"adjusted",cmap = cmp)

        return (iimg.mean()/img.mean(), iimg.flatten(), count)
        
            
    

def fixstrip_data(data,pplot = False,plotpath = "",hrad = 1,vrad = 1, sigma=5, alpha=0.3, method=2, processes=10,**kwargs):

    #根据colunm_name计算x,y值
    coords = np.int32([coord.split("-") for coord in list(data.columns)])
    x = coords[:,0].max()
    y = coords[:,1].max()

    print(f"读取图片大小为{x}*{y}！")
    print(f"共有{len(data)}个离子！")
    print(f"运行用{processes}个Cpus")
    print(f'using fixstrip method: {method}')
    clist = ["black","blue","green","yellow","red"]
    cmp = LinearSegmentedColormap.from_list("sm",clist )
    #定义输出的dataFrame
    newdf = pd.DataFrame(np.zeros(data.shape,dtype =np.float64),index = data.index,columns = data.columns)
    
    processes = min(processes, data.shape[0])
    pool = Pool(processes=processes)
    fcs = []
    results = []
    for idx,row in data.iterrows():
        result = pool.apply_async(fixstrip_single, (idx, row, x,y, method, hrad, vrad, cmp, pplot, plotpath, alpha, sigma, ))
        results.append(result)
    pool.close()
    pool.join()
    
    counter = []
    for re,idx in zip(results, data.index):
        fcs.append(re.get()[0])
        newdf.loc[idx,] = re.get()[1]
        counter.append(re.get()[2])
    
    fcs_arr = np.array(fcs)
    fcs_mean = np.round(fcs_arr[np.where(~np.isinf(fcs_arr) & ~(np.isnan(fcs_arr)))].mean(),4)

    print(f"转换完毕！共转换{sum(counter)}个离子！")
    print(f"转换后相比转换前倍数变化为{fcs_mean}")
    
    return newdf

#调用低通滤波
def fixstrip(infile,outfile,**kwargs):
    start = time.time()
    # 数据读取
    data = pd.read_csv(infile,sep = "\t",index_col = "FeatureID")

    # 数据处理
    # newdf = fixstrip_data_v1(data = data,**kwargs)
    newdf = fixstrip_data(data = data,**kwargs)
    
    # 数据保存
    newdf.to_csv(outfile,sep = "\t")
    
    print(f'耗时：{time.time()-start}s')

if __name__ == "__main__":
    parser = argparse.ArgumentParser("添加输入输出路径参数！")

    parser.add_argument("-i","--input",type = str,required = True,help = "path for input text file")
    parser.add_argument("-o","--output",type = str, required = True, help = "path for output text file")
    parser.add_argument("-a","--alpha",type = float, required = False, default=0.3, help = "for frequency weakening, when method is 2")
    parser.add_argument("-s","--sigma",type = float, required = False, default=10, help = "for frequency weakening, when method is 1")
    parser.add_argument("-m","--method",type = int, required = False, default=2, help = "method choices for frequency weakening, method(1: lpfilter; 2: lpfilter1; othernum: gaussianFreqFilter)")
    parser.add_argument("-p","--processes",type = int, required = False, default=10, help = "Cores number for multi-process")
    args = parser.parse_args()

    fixstrip(infile=args.input, outfile=args.output, pplot = True, alpha=args.alpha, sigma=args.sigma, method=args.method, processes=args.processes)
    
    
