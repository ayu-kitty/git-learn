import rpy2.robjects as robjects
from rpy2.robjects.packages import importr

def selectfile(**kwargs):
  lmbior = importr("lmbio")
  
  path = lmbior.selectfile(**kwargs)[0]
  
  return path
  
