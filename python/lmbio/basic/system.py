#!/opt/conda/bin/python

import os

def system(command):
  result = os.system(command)
  if result:
    raise ValueError(''.join(command)+"运行发生异常:")
  
