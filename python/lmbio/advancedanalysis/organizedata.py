#!/opt/conda/bin/python
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr

def organizedata(rawfile = ["数据矩阵.xlsx"],
                 rawdatafile = [""],
                 rawclassfile = [""],
                 rawcomparefile = [""],
                 **kwargs):
  lmbior = importr("lmbio")
  
  path = lmbior.organizedata(rawfile = robjects.StrVector(rawfile),
                             rawdatafile = robjects.StrVector(rawdatafile),
                             rawclassfile = robjects.StrVector(rawclassfile),
                             rawcomparefile = robjects.StrVector(rawcomparefile),
                             **kwargs)

if __name__ == '__main__':
  organizedata()
