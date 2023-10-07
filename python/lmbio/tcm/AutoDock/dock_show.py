#!/opt/conda/bin/python
from pymol import cmd
import os
import pandas as pd
import shutil
import argparse

def show_dockresult(receptor = "./protein.pdbqt",
                    ligand = "./ligand_out.pdbqt",
                    alpha_color = "[0.0, 0.533, 0.682]",
#                    alpha_color = "cyan",
                    beta_color = "yellow",
                    other_color = "[0.815, 0.325, 0.725]",
                    ligand_color = "green",
                    backcolor = "white",
                    out_path = "./out"):
    '''
    使用pymol对分子对接结果可视化
    '''
    if os.path.exists(out_path):
        shutil.rmtree(out_path)
        os.makedirs(out_path, exist_ok=True)
    else:
        os.makedirs(out_path, exist_ok=True)
    # 加载受体和配体分子结构文件
    cmd.load(ligand,'ligand')
    cmd.load(receptor,'receptor')
    cmd.show('cartoon', 'receptor')
    cmd.show('sticks', 'ligand')
    cmd.show('sticks', 'ligand within 5 of receptor')

    # 为α-螺旋着色
    cmd.set_color("helix_color", alpha_color)
    cmd.color("helix_color", "ss h")
    #cmd.color(alpha_color, "ss h")
    # 为β-折叠着色
    #cmd.set_color("h_color", beta_color)
    cmd.color("yellow", "ss s")
    #l+指的是loop及除beta sheet和helix以外的其他
    cmd.set_color("o_color", other_color)
    cmd.color("o_color", "ss l+")

    # 显示着色后的结构
    cmd.show("cartoon")
    cmd.hide("lines", "all")
    cmd.remove("hydrogens")
    cmd.bg_color(backcolor)
    # 设置颜色和透明度
    #cmd.color(receptor_color, 'receptor')
    cmd.set_color("l_color", "[0.60, 0.0, 0.0]")
    cmd.color("l_color", 'ligand')
    cmd.set('transparency', 0.6, 'receptor')
    cmd.set('transparency', 0.4, 'ligand')
    cmd.set("cartoon_fancy_helices",1)
    #cmd.set("cartoon_transparency", 0.8)
    cmd.set("valence",0)
    cmd.select("active_site", "byres all within 5 of ligand")
    cmd.show("lines", "active_site")
    cmd.set("line_width",1)

    cmd.color("grey", "active_site and (elem C*)")
    cmd.color(ligand_color, "ligand and (elem C*)")
    #cmd.set("cartoon_color", "red")
    cmd.distance("hbonds", "ligand", "active_site", 3.2, 2)
    cmd.hide("labels", "hbonds")
    cmd.set("dash_color", "black")
    cmd.set("dash_width",1.5)
    cmd.set("dash_length",0.3)
    # 添加阴影效果
    cmd.set('ray_trace_mode', 1)
    cmd.set('ambient', 0.3)
    cmd.set('specular', 0.6)
    cmd.set('shininess', 30)
    cmd.set("label_size",-0.8)
    cmd.set("label_font_id",7)
    cmd.set("depth_cue", 0)
    cmd.set("ray_trace_mode",0)
    cmd.set("ray_shadows",0)

    # 平移、旋转和缩放分子
    cmd.translate([10, 0, 0])
    #cmd.rotate('x', 90)
    cmd.zoom()
    # 创建表面表示
    #cmd.show('surface')

    # 添加标签或标注
    #cmd.label('atoms', '"text"')
    cmd.png(os.path.join(out_path,'docking_result1.png'),width=4000, height=4000, dpi=1500,ray=1)
    # 设置渲染参数
    cmd.rotate('x', 90)
    cmd.png(os.path.join(out_path,'docking_result2.png'),width=4000, height=4000, dpi=1500,ray=1)
    cmd.rotate('x', 90)
    cmd.png(os.path.join(out_path,'docking_result3.png'),width=4000, height=4000, dpi=1500,ray=1)
    cmd.rotate('x', 90)
    cmd.png(os.path.join(out_path,'docking_result4.png'),width=4000, height=4000, dpi=1500,ray=1)
    #cmd.quit()
    cmd.save(os.path.join(out_path,"3D_Modified.pse"))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--receptor",default = "./protein.pdbqt", help = "受体结构文件")
    parser.add_argument("-l","--ligand",default = "ligand_out.pdbqt", help = "分子对接结果文件")
    parser.add_argument("-ac", "--alpha_color",default = "[0.0, 0.533, 0.682]",help = "受体α-螺旋颜色")
    parser.add_argument("-tc", "--beta_color",default = "yellow",help = "受体β-折叠颜色")
    parser.add_argument("-oc", "--other_color",default = "[0.815, 0.325, 0.725]",help = "除α-螺旋,β-折叠外其他结构颜色")
    parser.add_argument("-lc", "--ligand_color",default = "green",help = "配体颜色，默认是green")
    parser.add_argument("-bc","--backcolor",default = "white", help = "背景颜色")
    parser.add_argument("-o","--out_path",default = "./out", help = "输出目录")
    args = parser.parse_args()
    show_dockresult(receptor = args.receptor,
                    ligand = args.ligand, 
                    alpha_color = args.alpha_color,
                    beta_color = args.beta_color,
                    other_color = args.other_color,
                    ligand_color = args.ligand_color,
                    backcolor = args.backcolor,
                    out_path = args.out_path)