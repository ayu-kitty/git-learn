import rpy2.robjects as robjects
from rpy2.robjects.packages import importr

# 输出路径
def packagepath(path = ''):
  lmbior = importr("lmbio")
  return lmbior.packagepath(path = path)[0]
  
# 输出路径
def databasepath(database = "database/",
                 path = ''):
  lmbior = importr("lmbio")
  return lmbior.databasepath(database = database,path = path)[0]
