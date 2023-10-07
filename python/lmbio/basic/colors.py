#!/opt/conda/bin/python
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr

def SelectColors(palette = "blindless",
                 n = robjects.NULL,
                 object = robjects.NULL,
                 subset = robjects.NULL,
                 outpalette = robjects.NULL):
  '''
  颜色选择函数
  palette: 颜色板名称
  n：颜色数量
  object：样本组别
  subset：组别名称
  outpalette：外部颜色板
  '''
  lmbior = importr("lmbio")
  
  if isinstance(subset,int):
    subset = [subset]
  elif isinstance(subset,tuple):
    subset = list(subset)
  
  if isinstance(subset,list):
    subset = robjects.StrVector(subset)
    
  if isinstance(object,int):
    object = [object]
  elif isinstance(subset,tuple):
    object = list(object)
  
  if isinstance(object,list):
    object = robjects.StrVector(object)
  
  col = lmbior.SelectColors(palette = palette,
                            n = n,
                            object = object ,
                            subset = subset,
                            outpalette = outpalette)
  
  return col

