import os
import glob

#将脚本同级目录下的所有文件夹中的QI文件cp至TOTAL目录下
if __name__ == "__main__":
    for dirs in os.listdir():
    files = glob.glob(r"{}\*ID.csv".format(dirs))
    for file in files:
        os.system(f"cp {file} {targetdir}")
    targetdir = f"F:/EMDB/定性/TOTAL/"
    dirlist = os.listdir()
    dirlist.remove("TOTAL")
