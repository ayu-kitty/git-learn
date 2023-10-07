#!/opt/conda/bin/python
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr

def metafetch(**kwargs):
  lmbior = importr("lmbio")
  
  data = lmbior.metafetch(**kwargs)
  
  return data

def profetch(**kwargs):
  lmbior = importr("lmbio")
  
  data = lmbior.profetch(**kwargs)
  
  return data[0]

def getpredealparams(**kwargs):
  lmbior = importr("lmbio")
  
  data = lmbior.getpredealparams(**kwargs)
  
  return data[0]
