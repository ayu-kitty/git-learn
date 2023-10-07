# File: resnet
# Author: jiang tao
# Time: 2023/3/10 9:57
import concurrent.futures
import copy
import os.path
import time

import torch
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from torchvision.models import resnet50, ResNet50_Weights
from pathlib import Path

import argparse

# from PIL import Image
# from torchvision import models, transforms
torch.manual_seed(134)


def generate_pic(data):
    # 处理格式
    columns = pd.Series(data.columns)[:].apply(lambda x: (int(x.split("-")[1]), int(x.split("-")[0])))
    data.columns = columns
    # 根据colunm_name计算x,y最大值
    coords = tuple(zip(*columns.values))
    x_cor = max(coords[1])
    y_cor = max(coords[0])

    print(f"读取图片大小为{x_cor}*{y_cor}！")
    print(f"共有{len(data)}个离子！")

    l = len(data)
    bg = np.zeros((l, y_cor + 1, x_cor + 1))
    i = 0
    for _, row in data.iterrows():
        for coo in row.keys():
            bg[i][coo] = row[coo]
        i += 1
    return torch.from_numpy(bg).to(torch.float32)


def initialize_model(conv_channel=1, out_features=2048):
    # Step 1: Initialize model with the best available weights
    weights = ResNet50_Weights.DEFAULT
    model = resnet50(weights=weights)

    model.conv1 = torch.nn.Conv2d(conv_channel, 64, kernel_size=7, stride=2, padding=3, bias=False)
    model.fc = torch.nn.Linear(in_features=2048, out_features=out_features, bias=True)

    return model


def corr_in(data, thresh=0.6):
    # 样本内mz之间
    data0 = data.iloc[:-1, :]
    cormat = np.corrcoef(data0)
    cor_mat = np.unique(cormat[cormat < 0.999])
    # 百分比
    pct = len(cor_mat[np.abs(cor_mat) > thres]) / len(cor_mat)
    plt.hist(he_corr, bins=50, weights=np.ones_like(he_corr) / len(he_corr))

    plt.title("corrcoef of sample and HE")
    plt.show()
    return pct


def corr_bt(data, thresh=0.6, figpath="HE_corr.jpg"):
    # 样本与HE
    he = data.iloc[-1, :]
    data0 = data.iloc[:-1, :]
    he_corr = pd.Series()
    for i, row in data0.iterrows():
        he_corr[row.name] = np.corrcoef(row.values, he.values)[1, 0]
    hc = he_corr[(he_corr < 0.999) & (he_corr > thresh)]
    pct = len(hc) / len(he_corr)

    plt.hist(he_corr, bins=50, weights=np.ones_like(he_corr) / len(he_corr))
    plt.xlabel("Correlations")
    plt.ylabel("Percentage")
    plt.title("corrcoef of sample and HE")
    plt.savefig(figpath)
    return hc, pct

def predict(model, img):
    """model predict wrapper"""
    return model(img.unsqueeze(0).unsqueeze(0)).detach().numpy()


def exportResetFeatures(data: pd.DataFrame, out_path: str, out_features: int = 2048,
                        draw: bool = False, parallel: bool = True) -> pd.DataFrame:
    corr_coef1, corr_coef2 = None, None
    if draw:
        hc1, corr_coef1 = corr_bt(data, figpath=str(Path(out_path) / 'HE_corr_before.png'))
    print(f"处理前corrcoef > 0.6的比例:{corr_coef1}")
    imgs = generate_pic(data)
    model = initialize_model(out_features=out_features)
    model.eval()

    # multiprocess speed this step
    output = None
    # error in Linux
    if parallel:
        order = []
        features = []
        counter = 0
        total = len(imgs)

        with concurrent.futures.ThreadPoolExecutor() as executor:
            print("INFO - 线程池中开启{}个线程 - resnet模型异步处理开始...".format(executor._max_workers))
            future_to_index = {executor.submit(predict, model, img): index for index, img in zip(data.index, imgs)}
            START = time.time()
            for future in concurrent.futures.as_completed(future_to_index):
                counter += 1
                index = future_to_index[future]
                feature = future.result()
                if counter % 50 == 0:
                    TIME = time.time()
                    print(f"INFO - 当前处理进度 - {counter} / {total} - elapsed time: {TIME - START:.4f}s")
                features.append(feature)
                order.append(index)
            END = time.time()
            print(f"INFO - 当前处理进度 - {total} / {total} - elapsed time: {END - START:.4f}s")
        tmp = pd.DataFrame(np.squeeze(features, axis=1), index=order,
                           columns=[f"f{i}" for i in range(out_features)])
        # re-order
        output = tmp.loc[data.index, :]

        # print(f"INFO - resnet特征提取 - 启用 ThreadingBackend with 100 workers...")
        # # param require can speed this step, but not required!
        # res = Parallel(n_jobs=100, verbose=10, require="sharedmem")(delayed(model)(img.unsqueeze(0).unsqueeze(0)) for img in imgs)
        # out = [r.detach().numpy() for r in res]
        # output = pd.DataFrame(np.squeeze(out, axis=1), index=data.index, columns=[f"f{i}" for i in range(out_features)])

    #     order = []
    #     features = []
    #     counter = 0
    #     total = len(imgs)
    #
    #     PROCESSES = 10
    #     print('Creating pool with %d processes\n' % PROCESSES)
    #     import multiprocessing
    #     with multiprocessing.Pool(PROCESSES) as pool:
    #         tasks = [(model, img) for img in imgs]
    #         features = [pool.apply_async(predictstar, t) for t in tasks]
    #     with concurrent.futures.ProcessPoolExecutor(max_workers=2) as executor:
    #         print("开启{}个进程，resnet模型异步处理开始...".format(executor._max_workers))
    #         tasks = [(model, img) for img in imgs]
    #         features = executor.map(predictstar, tasks)
    #         # future_to_index = {executor.submit(predict, model, img): index for index, img in zip(data.index, imgs)}
    #         # for future in concurrent.futures.as_completed(future_to_index):
    #         #     index = future_to_index[future]
    #         #     feature = future.result()
    #         #     counter += 1
    #         #     print(f"当前处理进度{counter} / {total}...")
    #         #     order.append(index)
    #         #     features.append(feature)
    #     tmp = pd.DataFrame(np.squeeze(features, axis=1), index=data.index, columns=[f"f{i}" for i in range(out_features)])
    #     # re-order
    #     output = tmp.loc[data.index, :]
    else:
    # 对每一副image依次处理
        out = []
        counter = 0
        total = len(imgs)
        for img in imgs:
            counter += 1
            out.append(model(img.unsqueeze(0).unsqueeze(0)).detach().numpy())
            print(f"INFO - resnet特征提取 - 当前处理进度: {counter} / {total}")
        # squeeze to n_mz * embedding_size
        output = pd.DataFrame(np.squeeze(out, axis=1), index=data.index, columns=[f"f{i}" for i in range(out_features)])
    if draw:
        hc2, corr_coef2 = corr_bt(output, figpath=str(Path(out_path) / 'HE_corr_after.png'))

    print(f"处理后corrcoef > 0.6的比例:{corr_coef2}")
    # print("高相关性mz:",hc2.sort_values())
    # TODO: 直接return, 避免中间文件
    CSV = str(Path(out_path) / "resnet_data.csv")
    output.to_csv(CSV)
    print(f"文件已经保存至{CSV}")
    return output

if __name__ == "__main__":
    # print(torch.cuda.is_available())
    # print(torch.cuda.device_count())
    data = pd.read_csv(r"C:\Users\Administrator\Desktop\data_r.txt", header=0, index_col=0)
    print(data)
    exportResetFeatures(data, out_path="C:/Users/Administrator/Desktop/", draw=False)

    # parser = argparse.ArgumentParser(description=None)
    # parser.add_argument("-i", "--infile", type=str, default=None, help="输入文件")
    # parser.add_argument("-o", "--outfile", type=str, default=None, help="输出路径")
    # parser.add_argument("-d", "--draw", type=bool, default=False, help="是否输出图片")
    # args = parser.parse_args()
    #
    # out_features = 2048
    # infile = args.infile
    # outfile = args.outfile
    # draw = args.draw
    #
    # data = pd.read_csv(infile, header=0, index_col="mz")
    # # data = data.set_index("mz")
    # if draw:
    #     hc1, corr_coef1 = corr_bt(data, figname='HE_corr_before.jpg')
    # print(f"处理前corrcoef > 0.6的比例:{corr_coef1}")
    # imgs = generate_pic(data)
    #
    # model = initialize_model(out_features=out_features)
    # model.eval()
    # # 对每一副image依次处理
    # out = []
    # for img in imgs:
    #     out.append(model(img.unsqueeze(0).unsqueeze(0)).detach().numpy())
    # # squeeze to n_mz * embedding_size
    # output = pd.DataFrame(np.squeeze(out, axis=1), index=data.index, columns=[f"f{i}" for i in range(out_features)])
    # if draw:
    #     hc2, corr_coef2 = corr_bt(output, figname='HE_corr_after.jpg')
    #
    # print(f"处理后corrcoef > 0.6的比例:{corr_coef2}")
    # # print("高相关性mz:",hc2.sort_values())
    # output.to_csv(outfile)
    # print(f"文件已经保存至{outfile}")

