#!/opt/conda/bin/python
from pymol import cmd
import re
import argparse
import os
import glob

def pymol_plip_2D(receptor = "./protein.pdb",
                    ligand_out = "./ligand_out.pdbqt",
                    backcolor = "white",
                    labelcolor = "purple",
                    out_path = "./interaction_visual"):
    '''
    使用pymolhe PLIP对分子对接结果复合物2D可视化
    '''
    #pymol.finish_launching()
    os.makedirs(out_path, exist_ok=True)
    
    ###生成蛋白配体复合物，现将ligand_out.pdbqt 转换成pdb
    os.system(f"/opt/conda/bin/obabel -ipdbqt {ligand_out} -opdb -O {out_path}/ligand_out.pdb")
    #分子对接结果生成复合物
    cmd.load(f"{out_path}/ligand_out.pdb",'ligand')
    cmd.load(receptor,'receptor')
    # 将蛋白质和配体合并到一个新的对象中
    cmd.create('complex', 'receptor or ligand')
    # 将配体移动到合适的位置
    cmd.align("ligand", "receptor")
    cmd.save(f"{out_path}/merged_complex.pdb", "complex")

    #PLIP分析
    cmd1 = f"/opt/conda/bin/plip -f {out_path}/merged_complex.pdb -pyxt -o {out_path}"
    # cmd1 = f"docker run --rm -v ${{PWD}}:/results -w /results -u $(id -u ${{USER}}):$(id -g ${{USER}}) pharmai/plip:latest -f {out_path}/merged_complex.pdb -pyxt -o {out_path}"
    # print(cmd1)
    os.system(cmd1)
    
    os.chdir(os.path.abspath(out_path))
    # 读取report.txt文件
    with open('report.txt', 'r') as file:
        report_text = file.read()
    # 定义正则表达式模式匹配链标识符（chain）、残基名称（resn）和原子名称（name）
    pattern = r"\|\s*(\d+)\s*\|\s*([A-Za-z]+)\s*\|\s*([A-Za-z])\s*\|"
    # 匹配并提取字段行
    matches = re.findall(pattern, report_text)
    ligand_residues = []
    for match in matches:
        resnr = match[0]
        restype = match[1]
        reschain = match[2]
        ligand_residues.append((resnr, restype, reschain))
    # 对数组进行去重
    unique_ligand_residues = list(set(ligand_residues))
    print (unique_ligand_residues)
    
    files = glob.glob("MERGED_COMPLEX_*.pse")
    for file in files:
        cmd.load(file, "ligand")
        cmd.color("green", "ligand and (elem C*)")
        #改作用残基颜色
        #cmd.select("active_site", "byres (ligand around 4)")
        #cmd.color("grey", "active_site and (elem C*)")
        # 遍历每个配体残基，显示名称
        for resnr,resnn,reschain in unique_ligand_residues:
            cmd.select("active_site", f"resi {resnr} and resn {resnn} and chain {reschain} and name N")
            #cmd.alter_state(1, "active_site", "(x, y, z) = (x, y, z + 1.0)")
            cmd.label("active_site", f"'{resnn}{resnr}{reschain}'")
        # 设置标签的颜色和其他属性
        cmd.set("label_color", labelcolor)
        cmd.set("label_bg_color","aquamarine")
        #cmd.set("label_font_id", 7)
        cmd.set("label_size", -1 )
        cmd.bg_color(backcolor)
        cmd.set("ray_trace_mode",0)
        cmd.set("ray_shadows",0)
        cmd.set("ray_trace_color","white")
        # 平移、旋转和缩放分子
        cmd.translate([0, 0, 0])
        # cmd.rotate('x', 90)
        # cmd.zoom()
        
        cmd.png("interaction_visualization1.png", width=1000, height=800, dpi=600,ray=1)
        cmd.rotate('y', 90)
        cmd.png("interaction_visualization2.png", width=1000, height=800, dpi=600,ray=1)
        cmd.rotate('y', 90)
        cmd.png("interaction_visualization3.png", width=1000, height=800, dpi=600,ray=1)
        cmd.rotate('y', 90)
        cmd.png("interaction_visualization4.png", width=1000, height=800, dpi=600,ray=1)
        cmd.save("2D_Modified.pse")
    #cmd.quit()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--receptor",default = "./protein.pdb",help = "受体结构文件,pdb格式")
    parser.add_argument("-lt", "--ligand_out",default = "./ligand_out.pdbqt",help = "分子对接结果,pdbqt格式")
    parser.add_argument("-bc","--backcolor",default = "white", help = "背景颜色")
    parser.add_argument("-lr","--labelcolor",default = "purple", help = "残基名称颜色")
    parser.add_argument("-o","--out_path",default = "./2D_visual", help = "输出目录")
    args = parser.parse_args()
    pymol_plip_2D(receptor = args.receptor,
                    ligand_out = args.ligand_out,
                    backcolor = args.backcolor,
                    labelcolor = args.labelcolor,
                    out_path = args.out_path)
