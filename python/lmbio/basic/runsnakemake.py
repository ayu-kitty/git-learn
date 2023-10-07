#!/opt/conda/bin/python
from lmbio.basic import path
import snakemake

def log_handler(msg):
    # print(msg)
    pass
    
def runsnakemake(snakefilepath,
                 config,
                 cores = 5,
                 configfiles = None,
                 targets = None,):
    snakefile = path.packagepath(path = snakefilepath)
    snakemake.snakemake(snakefile,
                        config = config,configfiles = configfiles,
                        cores = cores,targets = targets,
                        quiet = True,keepgoing = True,verbose = False,
                        log_handler = [])
