# File: network
# Author: jiang tao
# Time: 2023/5/23 10:42
import base64
import math

from io import BytesIO
from pathlib import Path

from typing import Union, List
from PIL import Image, ImageFont, ImageDraw
from sqlalchemy import create_engine,text
from .Constants import _engine
import attr
import numpy as np
from rdkit.Chem import MolToSmiles, MolFromInchi, MolFromSmiles, Draw, MolToInchiKey
import pandas as pd
import copy
import cv2
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np


engine = create_engine(_engine)
conn = engine.connect()

rules = pd.read_sql("herb_rules", conn)
# 支持反应类型换行
EN_REACTIONS = rules.set_index(["reaction"])["en"].apply(lambda x: x.replace("-", "\n")
                                                         .replace(" ", "\n")
                                                         .replace("_", "\n")).to_dict()
FLAGS = rules.set_index(["reaction"])["flag"].to_dict()
conn.close()

# 合并节点
ID_MAP = {}
ID_IMG = {}
LINKS = {}


@attr.s(auto_attribs=True)
class Node:
    id: str
    name: str
    category: str
    smiles: str
    phase: int
    rt: Union[float,None]=None
    reaction: Union[str,None]=None
    reactions: Union[str,None]=None
    clear: bool=True


@attr.s(auto_attribs=True)
class Link:
    source: str
    target: str
    value: str


def white_to_transparent(img):
    img = img.convert("RGBA")
    datas = img.getdata()
    newData = []
    for item in datas:
        if item[0] == 255 and item[1] == 255 and item[2] == 255:
            newData.append((255, 255, 255, 0))
        else:
            newData.append(item)
    img.putdata(newData)
    return img


def plot_nodes(row: pd.Series, border_in: int=1, border_out=150, bic: tuple=(0,0,0), size: tuple=(400,400), fmt='png', export_PIL = True):
    global FLAGS,ID_IMG
    # C:\Windows\Fonts\
    fontsize = 30
    font = ImageFont.truetype("arial.ttf", fontsize)

    id_plot = row["id_plot"]
    rt_plot = row["rt_plot"]

    _id = row.name
    rns, smi = ID_IMG.get(_id)

    img = np.array(Draw.MolToImage(MolFromSmiles(smi), size=size))

    # add border in
    img = cv2.copyMakeBorder(img, border_in, border_in, border_in, border_in, cv2.BORDER_CONSTANT, value=bic)
    # enlarge picture
    img = cv2.copyMakeBorder(img, border_out, border_out, border_out, border_out, cv2.BORDER_CONSTANT, value=(255,255,255))

    img = Image.fromarray(img)
    draw = ImageDraw.Draw(img)
    # put reactions 最大长度-OHOHOH
    if rns is not None:
        for i, r in enumerate(rns):
            f = FLAGS.get(r)
            draw.text((size[0]+1+border_out, size[0]/2+i*fontsize+border_out), f, font=font, fill=(0, 0, 0))
    # PUT rt
    if pd.notnull(rt_plot):
        n = len(rt_plot.split(","))
        #full_width = size[0] + 150
        if n <= 4:
            text_width = fontsize*4 + 4 + (fontsize*2+8)*n
            x = size[0]/2 - text_width/2+border_out
            y = size[0]+55+border_out
            rt = "RT: "+rt_plot
            draw.text((x,y), rt, font=font, fill=(0,0,0))
        else:
            text_width = fontsize*4 + 4 + (fontsize*2+8)*4
            x = size[0]/2 - text_width/2+border_out
            y = size[0]+55+border_out
            rt = "RT: "+",".join(rt_plot.split(",")[:4])+","
            draw.text((x, y), rt, font=font, fill=(0, 0, 0))
            # line2 超过四个换行处理
            rt = ",".join(rt_plot.split(",")[4:])
            text_width = fontsize * 4 + 4 + (fontsize * 2 + 8) * (n-4)
            x = size[0] / 2 - text_width / 2+border_out
            y = size[0] + 95+border_out
            draw.text((x, y), rt, font=font, fill=(0, 0, 0))

    # PUT id
    n = len(id_plot.split(","))
    text_width = (fontsize * 2 + 8) * n
    x = size[0] / 2 - text_width / 2+border_out
    y = size[0] + 15+border_out
    draw.text((x, y), id_plot, font=font, fill=(0,0,0))
    # alpha Image.fromarray(img.astype('uint8'))
    if export_PIL:
        img = white_to_transparent(img)
        return img
    img = white_to_transparent(img)
    #return img
    output_buffer = BytesIO()
    img.save(output_buffer, format=fmt,quality=95, subsampling=0)
    byte_data = output_buffer.getvalue()
    base64_str = base64.b64encode(byte_data).decode('utf-8')
    return f'image://data:image/{fmt};base64,' + base64_str


def coordinateGenerator(nodes: List[dict],step: int, interval: int=50, height: int=100):
    if step==0:
        for node in nodes:
            node["x"]=0
            node["y"]=3*height
    else:
        num = len(nodes)
        medium = math.floor(num/2)
        if num%2 == 0:
            for i, node in enumerate(nodes,1):
                if i <= medium:
                    node["x"]=-interval/2-interval*(medium-i)
                else:
                    node["x"]=interval/2+interval*(i-1-medium)
                node["y"]=(3-step)*height
        else:
            for i, node in enumerate(nodes,1):
                if i <= medium:
                    node["x"]=-interval*(medium-i+1)
                else:
                    node["x"]=interval*(i-medium-1)
                node["y"]=(3-step)*height
    return nodes


def nodeSorter(nodes_father: List[dict], nodes_son: List[dict])-> list:
    global LINKS
    sort_nodes = []
    keys = []
    for father in nodes_father: #layerUp
        fid = father["id"]
        for son in nodes_son: #layerNow
            sid = son["id"]
            if LINKS.get(sid) == fid:
                if sid not in keys:
                    keys.append(sid)
                    sort_nodes.append(son)
    return sort_nodes


def read_rawdata(rawdata: pd.DataFrame):
    if rawdata.shape[0] == 4 and "筛选" not in rawdata.columns:
        return
    rawdata["筛选"] = rawdata["筛选"].astype(str)
    rawdata["RT"] = rawdata["RT"].astype(float)
    # 双键还原 脱氢
    rawdata["Reaction Class"] = rawdata["Reaction Class"].apply(lambda x:x.replace('双键还原', '还原')
                                                                .replace('脱氢', '氧化')
                                                                .replace('氧化成醌', '氧化') if not pd.isna(x) else x)

    # add phase
    rawdata['phase'] = rawdata.apply(lambda x:x["Reaction Class"].count("+")+1 if pd.notnull(x["Reaction Class"]) else 0, axis=1)
    # 母体
    parent = rawdata.loc[rawdata["Reaction Class"].isnull(), :].reset_index(drop=True).loc[[0], :]
    # M0 RT保留两位小数
    parent["RT"] = parent["RT"].apply(lambda x: f"{x:.2f}" if pd.notnull(x) else x)

    rawdata = rawdata.query("筛选.str.contains('是') or 筛选.str.contains('候选')").reset_index(drop=True)
    if rawdata.empty:
        return
    rawdata["RT"] = rawdata["RT"].apply(lambda x: f"{x:.2f}" if pd.notnull(x) else x)

    # product
    product = rawdata.loc[rawdata["Reaction Class"].notnull(), :]\
        .dropna(subset=["筛选"])
    # 去除只画母核的情况
    if product.empty:
        return

    p1 = product.query("phase == 1").sort_values(["Reaction Class","RT", "筛选"])
    p2 = product.query("phase == 2").sort_values(["Reaction Class","RT", "筛选"])
    p3 = product.query("phase == 3").sort_values(["Reaction Class","RT", "筛选"])

    # merge
    rawdata = pd.concat([parent, p1, p2, p3]).reset_index(drop=True)

    # assign MID & IID
    counter = 0
    recorder = []
    flag = 0
    for i, row in rawdata.iterrows():
        rc = row["Reaction Class"].strip() if pd.notna(row["Reaction Class"]) else row["Reaction Class"]
        if i == 0:
            rawdata.loc[i, "Mid"] = f"M{counter}"
            # new
            recorder.append(0)
        else:
            metric = row["筛选"].strip()
            if metric.__contains__('是'):
                counter += 1
                rawdata.loc[i, "Mid"] = f"M{counter}"
                recorder.append(0)
            elif metric.__contains__('候选') and flag != rc:
                counter += 1
                recorder[-1] = 0
                recorder.append(recorder[-1]+1)
                rawdata.loc[i, "Mid"] = f"M{counter}-{recorder[-1]}"
            elif metric.__contains__('候选') and flag == rc:
                if recorder[-1] == 0:
                    counter += 1
                    recorder[-1] = 0
                recorder.append(recorder[-1]+1)
                rawdata.loc[i, "Mid"] = f"M{counter}-{recorder[-1]}"
            else:
                print("ERROR - 筛选列含有未知关键字...")
        flag = rc
    rawdata["I1id"] = rawdata["Mid"] + '-I1'
    rawdata["I2id"] = rawdata["Mid"] + '-I2'

    # add extra data
    rawdata["Intermediate 1 smiles"] = rawdata["Intermediate 1"].apply(lambda x:MolToSmiles(MolFromInchi(x)) if pd.notnull(x) else x)
    rawdata["Intermediate 1 inchikey"] = rawdata["Intermediate 1"].apply(lambda x:MolToInchiKey(MolFromInchi(x)) if pd.notnull(x) else x )
    rawdata["Intermediate 2 smiles"] = rawdata["Intermediate 2"].apply(lambda x:MolToSmiles(MolFromInchi(x)) if pd.notnull(x) else x)
    rawdata["Intermediate 2 inchikey"] = rawdata["Intermediate 2"].apply(lambda x:MolToInchiKey(MolFromInchi(x)) if pd.notnull(x) else x)

    return rawdata


def get_parent(rawdata):
    parent_id = rawdata.loc[0, "Mid"]
    parent_smi = rawdata.loc[0, "SMILES"]
    parent_rt = rawdata.loc[0, "RT"]
    parent_img = Draw.MolToImage(MolFromSmiles(parent_smi))

    return parent_id, parent_smi, parent_rt, parent_img


def get_nodes_links_df(rawdata, parent_id):
    global ID_MAP, ID_IMG, LINKS
    nodes, links = [], []
    # extract relationship
    for i, row in rawdata.iterrows():
        phase = row["phase"]
        reactions = row["Reaction Class"]
        flag = True if pd.notnull(row["筛选"]) and row["筛选"].__contains__('是') else False

        mid = row["Mid"]
        name = row["InChIKey"]
        smiles = row["SMILES"]
        rt = row["RT"]

        if phase == 0:
            nodes.append(Node(id=mid, name=name, category="Parent", smiles=smiles, phase=phase, rt=rt))
        elif phase == 1:
            reaction = reactions.strip()
            # fix bug jiangt@2023/5/24
            nodes.append(Node(id=mid, name=name, category="Product", smiles=smiles, reaction=reaction, clear=flag,
                              phase=phase, rt=rt, reactions=reactions))
            links.append(Link(source=parent_id, target=mid, value=reaction))
        elif phase == 2:
            _ = reactions.strip().split('+')
            reaction1 = _[0].strip()
            reaction2 = _[1].strip()

            i1_id = row["I1id"]
            i1_name = row["Intermediate 1 inchikey"]
            i1_smiles = row["Intermediate 1 smiles"]

            # Intermediate1
            nodes.append(Node(id=i1_id, name=i1_name, category="Intermediate", smiles=i1_smiles, reaction=reaction1, clear=False, phase=1,
                              reactions=reaction1))
            links.append(Link(source=parent_id, target=i1_id, value=reaction1))

            # Product2
            nodes.append(Node(id=mid, name=name, category="Product", smiles=smiles, reaction=reaction2, clear=flag, phase=phase,
                              rt=rt, reactions=reactions))
            links.append(Link(source=i1_id, target=mid, value=reaction2))
        elif phase == 3:
            _ = reactions.strip().split('+')
            reaction1 = _[0].strip()
            reaction2 = _[1].strip()
            reaction3 = _[2].strip()

            i1_id = row["I1id"]
            i1_name = row["Intermediate 1 inchikey"]
            i1_smiles = row["Intermediate 1 smiles"]

            i2_id = row["I2id"]
            i2_name = row["Intermediate 2 inchikey"]
            i2_smiles = row["Intermediate 2 smiles"]

            # Intermediate1
            nodes.append(Node(id=i1_id, name=i1_name, category="Intermediate", smiles=i1_smiles, reaction=reaction1, clear=False, phase=1,
                              reactions=reaction1))
            links.append(Link(source=parent_id, target=i1_id, value=reaction1))

            # Intermediate2
            nodes.append(Node(id=i2_id, name=i2_name, category="Intermediate", smiles=i2_smiles, reaction=reaction2, clear=False, phase=2,
                              reactions=f"{reaction1} + {reaction2}"))
            links.append(Link(source=i1_id, target=i2_id, value=reaction2))

            #Product3
            nodes.append(Node(id=mid, name=name, category="Product", smiles=smiles, reaction=reaction3, clear=flag, phase=phase,
                              rt=rt, reactions=reactions))
            links.append(Link(source=i2_id, target=mid, value=reaction3))

    nodes_df = pd.DataFrame(
        columns=["id", "name", "category", "smiles", "reaction", "clear", "phase", "rt", "reactions"])
    # arrange data into df
    for i, node in enumerate(nodes):
        nodes_df.loc[i, "id"] = node.id
        nodes_df.loc[i, "name"] = node.name
        nodes_df.loc[i, "category"] = node.category
        nodes_df.loc[i, "smiles"] = node.smiles
        nodes_df.loc[i, "reaction"] = node.reaction
        nodes_df.loc[i, "clear"] = node.clear
        nodes_df.loc[i, "phase"] = node.phase
        nodes_df.loc[i, "rt"] = node.rt
        nodes_df.loc[i, "reactions"] = node.reactions

    # 1.中间产物与产物合并， 保留产物
    # index saved -> nodes_df
    nodes_df_index_saved = []
    nodes_df_tmp = copy.deepcopy(nodes_df)
    nodes_df_tmp["reactions"] = nodes_df_tmp["reactions"].astype(str)
    uniq_ids = set(nodes_df_tmp["reactions"].to_list())
    for r in uniq_ids:
        # 保证 product 排在前面
        tmp = nodes_df_tmp.query("reactions == @r").sort_values(["category"], ascending=False)

        tmpIndex = tmp.index
        # 保留的第一个物质的信息
        id = tmp.loc[tmpIndex[0], "id"]
        _rt = tmp.loc[tmpIndex[0], "rt"]
        # 不需要合并
        if tmp.shape[0] == 1:
            nodes_df_index_saved.append(tmpIndex[0])
            if 'I' not in id:
                # new only for blood report -jiangt@2023/5/24 am
                # rm this @2023/5/24 pm
                # nodes_df.loc[tmpIndex[0], "clear"] = True
                if "-" in id:
                    nodes_df.loc[tmpIndex[0], "id_plot"] = id.split("-")[0]
                else:
                    nodes_df.loc[tmpIndex[0], "id_plot"] = id
            else:
                nodes_df.loc[tmpIndex[0], "id_plot"] = 'Intermediate'
            if _rt:
                nodes_df.loc[tmpIndex[0], "rt_plot"] = f"{_rt} min"
            else:
                nodes_df.loc[tmpIndex[0], "rt_plot"] = None
        # 需要合并
        else:
            # new
            nodes_df_index_saved.append(tmpIndex[0])
            if tmp.query("category == 'Product'").shape[0] == 1:
                #nodes_df.loc[tmpIndex[0], "clear"] = True
                ...
            else:
                nodes_df.loc[tmpIndex[0], "clear"] = False

            # 合并id 和 rt
            ids = []
            rts = []
            for i, row in tmp.iterrows():
                _ = row["id"]
                # 合并后的节点id map成第一个物质的id
                ID_MAP[_] = id
                # 多个产物合并时合并保留时间
                if row["rt"]:
                    rts.append(row["rt"])
                # 多个产物合并时合并Mid
                if 'I' not in _:
                    ids.append(_.split('-')[0])


            rt_plot = ",".join([str(t) for t in sorted(set([float(e) for e in rts]))]) + " min" if rts else None
            ids = np.asarray(list(set(ids)))
            arr = np.asarray([int(_[1:]) for _ in ids])
            ids = list(ids[np.argsort(arr)])
            id_plot = ",".join(ids) if ids else 'Intermediate'
            nodes_df.loc[tmpIndex[0], "id_plot"] = id_plot
            nodes_df.loc[tmpIndex[0], "rt_plot"] = rt_plot
    nodes_df_merged = nodes_df.loc[nodes_df_index_saved, :].set_index("id").sort_values(["phase"])

    # links
    links_df = pd.DataFrame(columns=["source", "target", "value"])
    # arrange data into df
    for i, link in enumerate(links):
        links_df.loc[i, "source"] = link.source
        links_df.loc[i, "target"] = link.target
        links_df.loc[i, "value"] = EN_REACTIONS.get(link.value)
    # ID_mapping
    for i, row in links_df.iterrows():
        source = row["source"]
        new_source = ID_MAP.get(source) if ID_MAP.get(source) else source
        links_df.loc[i, "source"] = new_source

        target = row["target"]
        new_target = ID_MAP.get(target) if ID_MAP.get(target) else target
        links_df.loc[i, "target"] = new_target
    links_df.drop_duplicates(inplace=True)

    LINKS = links_df.set_index(["target"])["source"].to_dict()

    # plot
    for id, row in nodes_df_merged.iterrows():
        if id == 'M0':
            # img_rns, img_smi
            ID_IMG[id] = (None, row["smiles"])
        else:
            clear = row["clear"]
            if clear:
                ID_IMG[id] = (None, row["smiles"])
            else:
                src = LINKS.get(id)
                reaction = row["reaction"]
                rns, smi = ID_IMG.get(src)
                rns = copy.deepcopy(rns)
                if rns is None:
                    ID_IMG[id] = ([reaction], smi)
                else:
                    rns.append(reaction)
                    ID_IMG[id] = (rns, smi)

    nodes_df_merged["symbol"] = nodes_df_merged.apply(lambda x: plot_nodes(x), axis=1)
    nodes_df_merged["symbolSize"] = 200
    nodes_df_merged.reset_index(drop=False, inplace=True)

    n_nodes_max = max(nodes_df_merged.query("phase == 1").shape[0],
        nodes_df_merged.query("phase == 2").shape[0],
        nodes_df_merged.query("phase == 3").shape[0])

    n_level = 1
    if nodes_df_merged.query("phase == 1").shape[0]:
        n_level = 2
    if nodes_df_merged.query("phase == 2").shape[0]:
        n_level = 3
    if nodes_df_merged.query("phase == 3").shape[0]:
        n_level = 4

    return nodes_df_merged, links_df, n_nodes_max, n_level


def prepare_plot_data(nodes_df_merged, links_df, interval=50, height=100):
    nodes0 = nodes_df_merged.loc[
        nodes_df_merged["phase"] == 0, ["id", "name", "category", "symbol", "symbolSize"]].to_dict('records')
    nodes1 = nodes_df_merged.loc[
        nodes_df_merged["phase"] == 1, ["id", "name", "category", "symbol", "symbolSize"]].to_dict('records')
    nodes2 = nodes_df_merged.loc[
        nodes_df_merged["phase"] == 2, ["id", "name", "category", "symbol", "symbolSize"]].to_dict('records')
    nodes3 = nodes_df_merged.loc[
        nodes_df_merged["phase"] == 3, ["id", "name", "category", "symbol", "symbolSize"]].to_dict('records')

    nodes0 = coordinateGenerator(nodes0, step=0, interval=interval, height=height)
    nodes1 = coordinateGenerator(nodes1, step=1, interval=interval, height=height)

    nodes2 = nodeSorter(nodes1, nodes2)
    nodes2 = coordinateGenerator(nodes2, step=2, interval=interval, height=height)

    nodes3 = nodeSorter(nodes2, nodes3)
    nodes3 = coordinateGenerator(nodes3, step=3, interval=interval, height=height)

    nodes = nodes0 + nodes1 + nodes2 + nodes3

    links = links_df[["source", "target", "value"]].to_dict('records')

    categories = [{"name": "Parent"}, {"name": "Product"}, {"name": "Intermediate"}]

    return nodes, links, categories


def run_one(data: pd.DataFrame):
    rawdata = read_rawdata(data)
    if rawdata is None:
        return
    # 导出重新排列的编号表
    # rawdata.set_index(["Mid"], drop=True)
    parent_id, parent_smi, parent_rt, parent_img = get_parent(rawdata)

    nodes_df_merged, links_df, n_nodes_max, n_level = get_nodes_links_df(rawdata, parent_id)

    nodes, links, categories = prepare_plot_data(nodes_df_merged, links_df, interval=20, height=100)
    return rawdata, n_nodes_max, nodes, links, categories


def plot_networkx(data: pd.DataFrame, save_path: str, cpd_name: str):
    global ID_MAP, ID_IMG, LINKS
    ID_MAP, ID_IMG, LINKS = {}, {}, {}
    # 1
    rawdata, n_nodes_max, nodes, links, categories = run_one(data)
    n_nodes_max = 5 if n_nodes_max <= 5 else n_nodes_max
    # 2
    G = nx.Graph()
    pos = {}
    for node in nodes:
        G.add_node(node['id'], **node)
        pos[node['id']] = np.array([node['x'], node['y']])
    for edge in links:
        G.add_edge(edge["source"], edge["target"], **edge)
    # 3
    # params: font_size=2.5, icon ~ 0.001, interval=20, height=75
    options = dict(arrows=True, arrowstyle="->", arrowsize=5, width=0.2,
                   min_source_margin=30, min_target_margin=26)

    fig, ax = plt.subplots(figsize=(n_nodes_max+1, 6), dpi=1500) # figsize=(9,5), dpi=100
    # networkx_edges
    nx.draw_networkx_edges(G, pos=pos, ax=ax, node_size=702, node_shape="o", **options)
    nx.draw_networkx_edge_labels(
        G, pos, edge_labels={(u, v): d["value"] for u, v, d in G.edges(data=True)}
        # {'center', 'top', 'bottom', 'baseline', 'center_baseline'}
        ,font_size=4, horizontalalignment="center", verticalalignment="center", rotate=True,label_pos=0.5,
        bbox=dict(boxstyle='round', ec=(1.0, 1.0, 1.0), fc=(1.0, 1.0, 1.0), alpha=1)
    )
    tr_figure = ax.transData.transform
    tr_axes = fig.transFigure.inverted().transform
    # (ax.get_xlim()[1] - ax.get_xlim()[0]) * 0.0005
    icon_size = 0.15
    icon_center = icon_size / 2
    for n in G.nodes:
        xf, yf = tr_figure(pos[n])
        xa, ya = tr_axes((xf, yf))
        # get overlapped axes and plot icon
        a = plt.axes([xa - icon_center, ya - icon_center, icon_size, icon_size])
        a.imshow(G.nodes[n]["symbol"])
        a.axis("off")
    fig.patch.set_alpha(1.)
    ax.axis('off')

    # save_path = str(Path(self.data) / project / "数据分析" / "4.网络图" / f"cpd——name")
    plt.savefig(f"{save_path}/png/{cpd_name}.png")
    plt.savefig(f"{save_path}/pdf/{cpd_name}.pdf")
    plt.close()

    return rawdata[["network_id", "Mid", "Reaction Class", "SMILES"]]
