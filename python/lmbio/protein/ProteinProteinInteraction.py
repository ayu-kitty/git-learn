#!/opt/conda/bin/python
import networkx as nx
from adjustText import adjust_text
import matplotlib.pyplot as plt
import matplotlib.colors as colors
#matplotlib.use('Agg')
import pandas as pd
import numpy as np
import os
import shutil
from collections import Counter
import argparse

def read(file):
    if file.endswith('.csv'):
        data = pd.read_csv(file,float_precision='round_trip')
    elif (file.endswith('.xlsx') or file.endswith('.xls')):
        try:
            data = pd.read_excel(file)
        except Exception:
            data = pd.read_csv(file,sep="\t",float_precision='round_trip')
    else:
        try:
            data = pd.read_csv(file,sep="\t",float_precision='round_trip')
        except ValueError as e:
            data = pd.read_csv(file,sep="\t",encoding='gb2312',float_precision='round_trip')
    return data

def _ppi(dproscore, epro, epath=None, dprofc=None, dprogene=None, figname="circular_layout",fontsize=None,font=None,
        layout="circular_layout", title=None, adjust=True,adjustsize=True,upcolor="red",dwcolor="green",shape=None):
    """
       Create a plot of ppi

    Parameters
    ----------
    dpath: dict, the dict of pathway and p-value

    dprofc: dict, Optional parameters,dict of protein and fc,default is None

    dproscore: dict,  the score of protein

    epath: list, connection of path and protein

    epro: list, connection of protein and protein

    figname: string, Optional parameters

    layout: string, layout of plot"""
    if dprofc is not None:
        max_fc = [x for x in dprofc.values() if x != np.inf]
        if len(max_fc) == 0:
            max_fc = [10]
        max_fc = max(max_fc)
        for i,j in dprofc.items():
            if j == np.inf:
                dprofc[i] = max_fc

    fig, ax = plt.subplots(figsize=(10, 10))

    # fontfamily=["Arial","Times","Verdana"]
    plt.rcParams['font.family'] = args.font
    # plt.rcParams['font.sans-serif'] = ['SimHei']
    #plt.rcParams['axes.unicode_minus'] = False

    G = nx.Graph()

    # if dpath is not None:
    #     dpath = sorted(dpath.items(), key=lambda x: len(x[0]), reverse=True)
    #     # 添加通路节点,和p值属性
    #     for path in dpath:
    #         G.add_node(path[0], pvalue=path[1])

    # 添加蛋白/基因节点,记录标准化后的得分---大小
    dpro = sorted(dproscore.items(), key=lambda x: x[1])
    s = np.array([x[1] for x in dpro])
    if len(s) == 0:
        return
    if max(s) == min(s):
        s = [300 for x in s]
    else:
        s = ((s - s.min()) / (s.max() - s.min())) * 380 + 100

    # 添加蛋白/基因节点
    for pro, score in zip(dpro, s):
        G.add_node(pro[0], score=score)

    if layout == "circular_layout":
        nodePos = nx.circular_layout(G)
    elif layout == "spring_layout":
        nodePos = nx.spring_layout(G)

    # 绘制通路节点
    # if dpath is not None:
        # nx.draw_networkx_nodes(G, nodePos, nodelist=[x[0] for x in dpath], node_shape="s",
                               # node_color=[x[1] for x in dpath], cmap=plt.get_cmap('cool'))

    colormap = {
    "red":[(1,0.27,0)],
    "orange":[(1,0.65,0)],
    "yellow":[(1,1,0)],
    "green":[(0,0.5,0)],
    "blue":[(0,0,1)],
    "purple":[(0.4,0.2,0.6)]

    }
    # 分别定义上下颜色映射
    color_up = colormap[str(upcolor)]
    cmup = colors.ListedColormap(color_up, 'indexed')

    color_dowm = colormap[str(dwcolor)]
    cmdown = colors.ListedColormap(color_dowm, 'indexed')

    # 绘制节点
    if dprofc is not None:
        uplist = [x for x in dprofc if float(dprofc[x]) > 1]
        upcolor_list = [dprofc[x] for x in uplist]
        s1 = np.array([G.nodes.get(x)["score"] for x in uplist])

        downlist = [x for x in dprofc if float(dprofc[x]) < 1]
        downcolor_list = [dprofc[x] for x in downlist]
        s2 = np.array([G.nodes.get(x)["score"] for x in downlist])
        #节点形状
        if adjustsize == False:
            
            s0 = 300
            nx.draw_networkx_nodes(G, nodePos, nodelist=uplist, node_shape="o",
                               node_color=upcolor_list,
                               cmap=plt.get_cmap(cmup), node_size=s0)
            nx.draw_networkx_nodes(G, nodePos, nodelist=downlist, node_shape="o",
                                   node_color=downcolor_list,
                                   cmap=plt.get_cmap(cmdown),node_size=s0)
        else:
            nx.draw_networkx_nodes(G, nodePos, nodelist=uplist, node_shape="o",
                                   node_color=upcolor_list,
                               cmap=plt.get_cmap(cmup), node_size=s1)
                                   
            nx.draw_networkx_nodes(G, nodePos, nodelist=downlist, node_shape="o",
                                    node_color=downcolor_list,
                               cmap=plt.get_cmap(cmdown), node_size=s2)
                                   
    else:
        # 没有上下调的颜色
        s0 = np.array([x[1] for x in dpro])
        s0 = ((s - s.min()) / (s.max() - s.min())) * 380 + 100
        if adjustsize == False:
            s0 = 380
            nx.draw_networkx_nodes(G, nodePos, nodelist=[x[0] for x in dpro], node_shape="o",
                               node_color="grey",
                               node_size=s0)
        else:      
            nx.draw_networkx_nodes(G, nodePos, nodelist=[x[0] for x in dpro], node_shape="o",
                               node_color="grey",
                               node_size=s0)
    if dprogene is None:
        if adjust is True:
            texts = [plt.text(1.06*j[0], 1.06*j[1], i, fontsize=fontsize) for i, j in nodePos.items()]
            adjust_text(texts)  
        else:
            for i, j in nodePos.items():
                # print(i,j,100)
                plt.text(1.06*j[0], 1.06*j[1], i, horizontalalignment='center')

    if dprogene is not None:
        if adjust is True:
            texts = [plt.text(j[0], j[1], dprogene[i], fontsize=fontsize) for i, j in nodePos.items()]
            adjust_text(texts)
        else:
            for i, j in nodePos.items():
                plt.text(1.06*j[0], 1.06*j[1], dprogene[i], horizontalalignment='center')
    # 绘制边
    # if epath is not None:
    #     nx.draw_networkx_edges(G, nodePos, edgelist=epath, alpha=0.2, style="--")
    nx.draw_networkx_edges(G, nodePos, edgelist=epro, alpha=0.2, style="-")
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    plt.title(title)
    plt.axis('off')
    plt.ylim(-1.2,1.2)
    plt.tight_layout()
    plt.savefig("ppi_result/{name}.png".format(name=figname))
    plt.savefig("ppi_result/{name}.pdf".format(name=figname))
    return

#From unionplot.py
def graphs(top=40,FC=False,upcolor="red",dwcolor="green",adjustsize=None,fontsize=None,font=None):
    def score(FC=FC,topl=None):
        # if term is True:
        #     data = pd.read_csv("ppi_selected_term.xls", sep="\t", index_col="term")
        #     term10 = data.index.tolist()

        # 根据 protein2protein_network背景 和 差异蛋白列表 准备 蛋白的连接 和 score

        dif = read(args.nodes)
        # 去重
        dif.drop_duplicates(subset=[dif.columns[0]],keep='first',inplace=True)

        ## 区分蛋白或基因绘图
        if "Accession" in dif.columns and "Gene Name" in dif.columns:
            #提取所有的基因,Accession:geneName的字典
            dif = dif.set_index("Accession")
            for i in dif.index:
                    if dif.loc[i,"Gene Name"] == "  ":
                        dif.loc[i,"Gene Name"] = i
                    if pd.isnull(dif.loc[i,"Gene Name"]):
                        dif.loc[i,"Gene Name"] = i
            dprogene = dif['Gene Name'].to_dict()
            w_plot = 'both'
        elif "Accession" in dif.columns and "Gene Name" not in dif.columns:
            dif = dif.set_index("Accession")
            dprogene = None
            w_plot = 'query'
        elif "Gene Name" in dif.columns and "Accession" not in dif.columns:
            dif = dif.set_index("Gene Name")
            dprogene = None
            w_plot = 'gene'

        # 根据原始 network 背景整理分析所需文件
        al = read(args.edges)
        edge = al[(al['node1'].isin(dif.index))&(al['node2'].isin(dif.index))]
        edge.to_csv("ppi_result/ppi_network.xls", sep="\t",index=0)

        # 控制 蛋白节点 输入
        if topl is None:
            #提取所有node1和node2的蛋白accesion（无重复）
            nodes_p = edge.iloc[:, 0].values.tolist() + edge.iloc[:,1].values.tolist()
            #nodes_p = [x for x in nodes_p if x not in term10]
        else:
            nodes_p = list(topl)

        # 计算权重得分
        dic_score = {}
        la = edge.loc[:, "node1"].values.tolist() + edge.loc[:, "node2"].values.tolist()
        #以蛋白名为键，以次数为值的字典
        lat = Counter(la).items()
        for i,j in lat:
            dic_score[i] = j
        #添加Degree
        d = Counter(la)
        dif["degree"] = 0
        for i in d:
            if i in dif.index:
                dif.loc[i,"degree"] = d[i]
        dif.to_csv("ppi_result/ppi_nodes.xls",sep="\t")
        #蛋白:次数
        dproscore = {x: dic_score[x] for x in nodes_p}
        #node1、node2的所有蛋白(有重复)
        bian = edge.iloc[:, [0, 1]].values.tolist()
        edgx = len(edge[edge["combined_score"].isnull()])
        epath = bian[0:edgx]
        # epath = [x for x in epath if x[0] in nodes_p and x[1] in nodes_p]
        epro = bian[edgx:]
        #嵌套列表，node1和node2
        epro = [x for x in epro if x[0] in nodes_p and x[1] in nodes_p]

        if FC is True:
            dprofc = dif["FC"][nodes_p].to_dict()
           # dprofc = dif["FC"].reindex(columns=nodes_p).to_dict()                
        else:
            dprofc = None



        return dproscore, epro, epath, dprofc, dprogene, w_plot

    dproscore, epro, epath, dprofc, dprogene , w_plot = score(FC=FC,topl=None)
    # print(epath,123)
    pro_path = [x for x in set([x[1] for x in epath])]
    # print(epro, 123)
    #
    ds = sorted(dproscore.items(), key=lambda x: x[1],reverse=True)

    if top is not None:
        ds = ds[0:25]
        topl = [x[0] for x in ds] + pro_path
    else:
        topl = [x[0] for x in ds]

    dproscore, epro, epath, dprofc, dprogene, w_plot = score(FC=FC,topl=topl)

    if w_plot == "both":
        _ppi(dproscore, epro, epath=epath, dprofc=dprofc, dprogene=None,upcolor=upcolor,dwcolor=dwcolor,figname="ppi_query",
                layout="circular_layout", adjust=True,title='',adjustsize=adjustsize,fontsize=fontsize,font=font)
        _ppi(dproscore, epro, epath=epath, dprofc=dprofc, dprogene=dprogene,upcolor=upcolor,dwcolor=dwcolor,figname="ppi_gene",
            layout="circular_layout", adjust=True,title='',adjustsize=adjustsize,fontsize=fontsize,font=font)
    else:
        _ppi(dproscore, epro, epath=epath, dprofc=dprofc, dprogene=None,upcolor=upcolor,dwcolor=dwcolor,figname="ppi_"+str(w_plot),
                layout="circular_layout", adjust=True,title=str(w_plot)+"_network",adjustsize=adjustsize,fontsize=fontsize,font=font)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Protein Protein interaction')
    parser.add_argument('-n','--nodes',help='input file for nodes',type=str)
    parser.add_argument('-e','--edges',help='input file for edges',type=str)
    parser.add_argument('-u','--upcolor',help='set up color',type=str,nargs="?",default='red')
    parser.add_argument('-d','--dwcolor',help='set down color',type=str,nargs="?",default='green')
    parser.add_argument('-a','--adjustsize',help='set scaling ratio',type=str,nargs='?',default='yes')
    parser.add_argument('-s','--fontsize',help='set scaling ratio',type=int,nargs="?",default=12)
    parser.add_argument('-f','--font',help='set font family',type=str,nargs="?",default="Arial")
    parser.add_argument('-c','--foldchange',help='include foldchange value',type=str,nargs='?',default='no')
    args = parser.parse_args()
    if args.foldchange == 'yes':
        FC = True
    else:
        FC = False
    if args.adjustsize == 'yes':
        adjustsize = True
    else:
        adjustsize = False
    path = os.getcwd()
    os.mkdir(path+"/ppi_result")
    graphs(top=40,FC=FC,upcolor=args.upcolor,dwcolor=args.dwcolor,adjustsize=
        adjustsize,fontsize=args.fontsize,font=args.font)
    zip_name = "ppi_result.zip"
    os.system(f'zip -r {zip_name} ppi_result > /dev/null')
    shutil.rmtree(path+"/ppi_result")
    #删除隐藏文件夹
    [shutil.rmtree(f) for f in os.listdir('.') if os.path.isdir(f) and f.startswith('.')]
    [os.remove(f) for f in os.listdir('.') if os.path.isfile(f) and f.startswith('.')]


