# File: picMatch
# Author: jiang tao
# Time: 2023/3/7 12:35
# !/opt/conda/bin/python

import cv2
import pandas as pd
import numpy as np
import argparse


# 读入呈像图,黑底白图
def rendering_pic(textfile, mode, out, write=True):
    data_rd = pd.read_csv(textfile, usecols=["x", "y"], header=0)
    h = np.max(data_rd.y) + 1
    w = np.max(data_rd.x) + 1
    img_rd = np.ones((h, w), dtype=np.uint8)
    for i, row in data_rd.iterrows():
        img_rd[row["y"], row["x"]] = 255
    if write:
        cv2.imwrite(f"{out}/rd_{mode}.png", img_rd)
    return img_rd, w, h


# rendering_pic修改
def load_msi(df: pd.DataFrame, mode, out, write=True):
    data_rd = df
    h = np.max(data_rd.y) + 1
    w = np.max(data_rd.x) + 1
    img_rd = np.ones((h, w), dtype=np.uint8)
    for i, row in data_rd.iterrows():
        img_rd[row["y"], row["x"]] = 255
    if write:
        cv2.imwrite(f"{out}/rd_{mode}.png", img_rd)
    return img_rd, w, h


# 处理染色图--要求为白底，转换为黑底白图
def HE_pic(he_pic, cutoff=200, write=False):
    # 灰度加载
    img_he = cv2.imread(he_pic, cv2.IMREAD_GRAYSCALE)
    # cutoff设为背景(白色)和图像的分割点
    thresh, img_he = cv2.threshold(img_he, cutoff, 255, cv2.THRESH_BINARY_INV)
    # 闭运算
    kernel1 = np.uint8((30, 30))
    img_he = cv2.morphologyEx(img_he, cv2.MORPH_CLOSE, kernel1, iterations=10)
    # img_he_dilate = cv2.dilate(img_he,kernel1,iterations= 10)
    if write:
        cv2.imwrite("HE.png", img_he)
    # 找轮廓
    # contours, hiearachy = cv2.findContours(img_he_close,cv2.RETR_EXTERNAL,cv2.CHAIN_APPROX_SIMPLE)
    # cnt = cv2.drawContours(img_he_close,contours,-1,(0,255,0))
    return img_he


def get_thum(img, size):
    # 缩放
    ret_img = cv2.resize(img, dsize=size)
    return ret_img


def rotate_nocrop(img, degree):
    h = img.shape[0]
    w = img.shape[1]
    M = np.zeros((2, 3), dtype=np.float32)
    alpha = np.cos(degree * np.pi / 180)
    beta = np.sin(degree * np.pi / 180)
    M[0, 0] = alpha
    M[1, 1] = alpha
    M[0, 1] = beta
    M[1, 0] = -beta
    cx = w / 2
    cy = h / 2
    tx = (1 - alpha) * cx - beta * cy
    ty = beta * cx + (1 - alpha) * cy
    M[0, 2] = tx
    M[1, 2] = ty
    # 计算图像边界
    bound_w = int(h * np.abs(beta) + w * np.abs(alpha))
    bound_h = int(h * np.abs(alpha) + w * np.abs(beta))
    # M = cv2.getRotationMatrix2D((w/2,h/2),degree,1)
    M[0, 2] += bound_w / 2 - cx
    M[1, 2] += bound_h / 2 - cy
    rt_img = cv2.warpAffine(img, M, dsize=(bound_w, bound_h))

    return rt_img


# 计算cos值来衡量相似度
def cos_value(img1, img2):
    images = [img1, img2]
    vectors = []
    norms = []
    for image in images:
        vector = []
        for pixel_tuple in image.ravel():
            vector.append(np.average(pixel_tuple))
            # vector.append(pixel_tuple)

        vectors.append(vector)
        norms.append(np.linalg.norm(vector, 2))
    a, b = vectors
    a_norm, b_norm = norms
    res = np.dot(a / a_norm, b / b_norm)
    return res


# 计算重合率来衡量相似度
def overlap(base_img, img):
    _, base_img = cv2.threshold(base_img, 1, 1, cv2.THRESH_BINARY)
    h, w = base_img.shape
    img = cv2.resize(img, (w, h))
    _, img = cv2.threshold(img, 1, 1, cv2.THRESH_BINARY)
    overlap_rate = 1 - np.sum(base_img ^ img) / (h * w)
    return overlap_rate


def match(base_img, img):
    # _,base_img=cv2.threshold(base_img,1,1,cv2.THRESH_BINARY)
    h, w = base_img.shape
    img = cv2.resize(img, (w, h))
    _, img = cv2.threshold(img, 1, 1, cv2.THRESH_BINARY)
    cor = cv2.matchTemplate(base_img, img, cv2.TM_CCORR_NORMED)[0, 0]
    # cor = cv2.matchTemplate(base_img,img,cv2.TM_CCOEFF_NORMED)[0,0]
    return cor


def find_match(he_img, rd_img, size):
    # 初始化：旋转角度0，相似度0
    max_similarity = [0, 0]
    similarities = []
    for degree in range(0, 360, 2):
        he_thum = rotate_nocrop(he_img, degree)
        he_thum = get_thum(he_thum, size)
        # cv2.imwrite("he_thum.jpg",he_thum)
        # 观测
        # time.sleep(0.1)

        similarity = match(he_thum, rd_img)
        similarities.append((degree, similarity))
        if similarity > max_similarity[1]:
            # print(degree,similarity)
            max_similarity[0] = degree
            max_similarity[1] = similarity

    return max_similarity, similarities


def write_first10(he_pic, similarities, flipcode, mode, outpath="", n=10):
    ss = sorted(similarities, key=lambda x: x[1], reverse=True)
    for i in range(n):
        print("Rotate degree:", ss[i][0], "----similarity:", ss[i][1])
        # cv2.imwrite(f"HE_Match_{ss[i][0]}_{mode}.png",rotate_nocrop(he_img,ss[i][0]))
        cropped_he = crop_img(he_pic, ss[i][0], flipcode)
        cv2.imwrite(f"{outpath}/HE_{ss[i][0]}_{mode}.png", cropped_he)
    return


def crop_img(he_pic, degree=0, flipcode=0):
    # 分别读入二值和彩色HE图
    img_he = cv2.imread(he_pic, cv2.IMREAD_GRAYSCALE)
    _, bin_img = cv2.threshold(img_he, 200, 255, cv2.THRESH_BINARY_INV)

    img = cv2.imread(he_pic)
    if flipcode:
        img = cv2.flip(img, 0)
        bin_img = cv2.flip(bin_img, 0)
    bin_img = rotate_nocrop(bin_img, degree)
    img = rotate_nocrop(img, degree)
    # 按照二值图进行裁剪
    row_mask = [True if np.any(row) else False for row in bin_img]
    col_mask = [True if np.any(x) else False for x in np.nditer(bin_img, order="F", flags=["external_loop"])]
    cropped_img = img[row_mask, :][:, col_mask]
    # 将黑色背景扣下来，替换为白色
    gray = cv2.cvtColor(cropped_img, cv2.COLOR_BGR2GRAY)
    _, mask = cv2.threshold(gray, 1, 255, cv2.THRESH_BINARY_INV)
    mask = np.stack([mask, mask, mask], axis=2)
    res = cv2.add(cropped_img, mask)
    res = cv2.cvtColor(res, cv2.COLOR_BGR2BGRA)
    # 提取白色背景
    k = 230
    for i in range(0, res.shape[0]):
        for j in range(0, res.shape[1]):
            if res[i, j, 0] > k and res[i, j, 1] > k and res[i, j, 2] > k:
                res[i, j, 3] = 0
    return res


def rotateMatch(he_pic, he_img, rd_img, size, mode=None, outpath="", cor=0.7):
    # 否考虑翻转
    # 有flip，则flipcode记为1；否则记为0
    no_flip_match, similarities1 = find_match(he_img, rd_img, size)

    flip_he_img = cv2.flip(he_img, 0)
    flip_match, similarities2 = find_match(flip_he_img, rd_img, size)
    if no_flip_match[1] > flip_match[1]:
        # cv2.imwrite(f"{outpath}/HE_bestMatch_{mode}.png",rotate_nocrop(he_img,no_flip_match[0]))
        no_flip_match.insert(0, 0)
        best_match = no_flip_match
        if best_match[2] < cor:
            write_first10(he_pic, similarities1, best_match[0], mode, outpath)

    else:
        # cv2.imwrite(f"{outpath}/HE_bestMatch_{mode}.png",rotate_nocrop(flip_he_img,flip_match[0]))
        flip_match.insert(0, 1)
        best_match = flip_match
        if best_match[2] < cor:
            write_first10(he_pic, similarities2, best_match[0], mode, outpath)

    cropped_he = crop_img(he_pic, best_match[1], best_match[0])
    cv2.imwrite(f"{outpath}/cropped_BestMatch_{mode}.png", cropped_he)
    print(f"Best Match for {mode} mode", best_match)
    return best_match


def rotateHEMatchMsi(df: pd.DataFrame, he_pic: str, mode: str, outpath:str):
    """
    auto match pipeline
    :param df: 成像图数据框
    :param he_pic: HE染色图绝对路径
    :param mode:    离子模式
    :param outpath: 输出文件绝对路径, 不要/
    :return:
    """
    rd_img, rw, rh = load_msi(df, mode, outpath)
    he_img = HE_pic(he_pic)
    size = (rw, rh)
    rotateMatch(he_pic, he_img, rd_img, size, mode, outpath)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=None)
    parser.add_argument("-o", "--out", type=str, default="./", help="输出路径,默认./")
    parser.add_argument("-t", "--textfile", type=str, default=None, help="成像图文件")
    parser.add_argument("-p", "--pic", type=str, default=None, help="HE染色图片")
    parser.add_argument("-m", "--mode", type=str, default=None, help="正负离子模式pos/neg")
    args = parser.parse_args()

    textfile = args.textfile
    he_pic = args.pic
    mode = args.mode
    outpath = args.out

    rd_img, rw, rh = rendering_pic(textfile, mode, outpath)
    he_img = HE_pic(he_pic)
    size = (rw, rh)
    # print(size)
    rotateMatch(he_pic, he_img, rd_img, size, mode, outpath)
