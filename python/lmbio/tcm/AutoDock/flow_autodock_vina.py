#!/opt/conda/bin/python
from pymol import cmd
import os
import pandas as pd
import shutil
import argparse
import dock_show
import PLIP_pymol_2D
import glob
import re

def autodock_visual(receptor = "./protein.pdbqt",
                    ligand = "./ligand.pdbqt",
                    alpha_color = "[0.0, 0.533, 0.682]",
                    beta_color = "yellow",
                    other_color = "[0.815, 0.325, 0.725]",
                    ligand_color = "green",
                    backcolor = "white",
                    labelcolor = "white",
                    out_path = "./"):
    '''
    蛋白质受体和配体对接，并对对接结果可视化
    '''
    os.makedirs(out_path, exist_ok=True)

    #1 蛋白质准备,pdb需要删掉配体以及水分子,加H原子
    #cmd1 = f"/root/ADFRsuite_x86_64Linux_1.0/bin/prepare_receptor \
    #            -r {receptor} \
    #            -A "hydrogens" \
    #            -o {out_path}/receptor.pdbqt"
    #os.system(cmd1)
    
    #2 配体准备:配体分子的3D结构,可以为mol/mol2/sdf格式,由于配体结构文件可能没有H原子，自动进行加氢处理
    #cmd2 = f"mk_prepare_ligand.py \
    #            -i {ligand} \
    #            -o {out_path}/ligand.pdbqt"
    #配体准备:配体分子的3D结构,可以为pdb or .mol2 or .pdbq format
    #cmd2 = f"/root/ADFRsuite_x86_64Linux_1.0/bin/prepare_ligand \
    #            -l {ligand} \
    #            -o {out_path}/ligand.pdbqt"
    #os.system(cmd2)
    
    #3 对接盒子(Grid Box)的设置: 通过AutoGrid4计算，自动生成
    os.chdir(os.path.abspath(out_path))
    cmd3 = f"python /opt/conda/lib/python3.11/site-packages/AutoDockTools/Utilities24/prepare_gpf4.py \
                -l {ligand} \
                -r {receptor} \
                -y"
    os.system(cmd3)
    
    #获取受体名称
    receptor_name = os.path.basename(receptor).rstrip(".pdbqt")

    cmd4 = f"/root/ADFRsuite_x86_64Linux_1.0/bin/autogrid4 \
                -p {receptor_name}.gpf \
                -l {receptor_name}.glg"
    os.system(cmd4)

    #4 如果采用AutoDock4力场进行分子对接，需要提供原子类型的映射表，即第3步中的map系列文件。并且在vina运行命令里面添加关键词--scoring ad4
    cmd5 = f"/root/vina_1.2.3_linux_x86_64 \
                --ligand ligand.pdbqt \
                --maps {receptor_name} \
                --scoring ad4 \
                --num_modes 1 \
                --exhaustiveness 32 \
                --out ligand_out.pdbqt"
    os.system(cmd5)
    
    #存储参数文件
    os.makedirs("config", exist_ok=True)
    os.system(f"mv {receptor_name}.*.map {receptor_name}.glg {receptor_name}.gpf {receptor_name}.maps.* config")

    #对接结果3D可视化
    dock_show.show_dockresult(receptor = receptor,
                    ligand = "ligand_out.pdbqt", 
                    alpha_color = alpha_color,
                    beta_color = beta_color,
                    other_color = other_color,
                    ligand_color = ligand_color,
                    backcolor = backcolor,
                    out_path = "3D_visual")
    #flow_autodock_visual -r 1iep_receptor.pdbqt -l 1iep_ligand_ad4_out.pdbqt -o visul
    #对接结果2D可视化
    PLIP_pymol_2D.pymol_plip_2D(receptor = "./protein.pdb",
                    ligand_out = "./ligand_out.pdbqt",
                    backcolor = backcolor,
                    labelcolor = labelcolor,
                    out_path = "./interaction_visual")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--receptor",default = "./protein.pdbqt",help = "受体结构pdbqt文件")
    parser.add_argument("-l", "--ligand",default = "./ligand.pdbqt", help = "配体结构文件pdbqt")
    parser.add_argument("-ac", "--alpha_color",default = "[0.0, 0.533, 0.682]",help = "受体α-螺旋颜色")
    parser.add_argument("-tc", "--beta_color",default = "yellow",help = "受体β-折叠颜色")
    parser.add_argument("-oc", "--other_color",default = "[0.815, 0.325, 0.725]",help = "除α-螺旋,β-折叠外其他结构颜色")
    parser.add_argument("-lc", "--ligand_color",default = "green",help = "配体颜色，默认是green")
    parser.add_argument("-bc","--backcolor",default = "white", help = "背景颜色")
    parser.add_argument("-lr","--labelcolor",default = "purple", help = "残基名称颜色")
    parser.add_argument("-o", "--out_path",default = "./", help = "输出目录")
    args = parser.parse_args()
    autodock_visual(receptor = args.receptor,
                    ligand = args.ligand,
                    alpha_color = args.alpha_color,
                    beta_color = args.beta_color,
                    other_color = args.other_color,
                    ligand_color = args.ligand_color,
                    backcolor = args.backcolor,
                    labelcolor = args.labelcolor,
                    out_path = args.out_path)
