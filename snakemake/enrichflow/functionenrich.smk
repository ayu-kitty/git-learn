#!/opt/conda/bin/python
import glob
import os
import re
import pandas as pd

def collect_reads():
    files = glob.glob("*-diff-*.xls")
    files.sort()
    reads = [file.split('/')[-1].split('-diff-')[0] for file in files]
    other = files[0].split('/')[-1].split('-diff-')[1]
    return reads,other

groups = collect_reads()[0]
enclass = {"G":"GO","K":"KEGG","I":"InterPro","R":"Reactome","W":"Wikipathways"}
hz = collect_reads()[1]
kg = pd.read_csv(config['database']+"gene_kegg.backgroud.xls",sep="\t")
omic = config["omic"]
if omic == "P":
    if config['org'] in ["hsa","mmu","rno"]:
        if config['enrichtype']=="all":
            entypes = ["GO","KEGG","InterPro","Reactome","Wikipathways"]
            rule all:
                input:
                    expand(config['savepath']+"{entype}/{group}/",entype=entypes,group=groups)
        elif config['enrichtype']=="ppi":
            rule all:
                input:
                    expand(config['ppipath']+"{group}/",group=groups)
        else:
            entypes = [enclass[key] for key in config['enrichtype'].split(",")]
            rule all:
                input:
                    expand(config['savepath']+"{entype}/{group}/",entype=entypes,group=groups)
    else:
        if config['enrichtype']=="all":
            entypes = ["GO","KEGG","InterPro"]
            rule all:
                input:
                    expand(config['savepath']+"{entype}/{group}/",entype=entypes,group=groups)
        elif config['enrichtype']=="ppi":
            rule all:
                input:
                    expand(config['ppipath']+"{group}/",group=groups)
        else:
            entypes = [enclass[key] for key in config['enrichtype'].split(",")]
            rule all:
                input:
                    expand(config['savepath']+"{entype}/{group}/",entype=entypes,group=groups)
elif omic == "M":
    rule all:
        input:
            expand(config['savepath']+"KEGG/{group}/",group=groups),
            expand(config['savepath']+"KEGG_map/{group}/",group=groups)

rule go:
    input:
        diff = "{group}-diff-"+hz
    output:
        directory(config['savepath']+"GO/{group}/")
    threads:3
    params:
        db = config['database'],
        savedir = config['savepath'],
        inputdir = config['inputpath'],
        fontfamily = config['fontfamily']
    shell:
        "enrichpip -t G -if '{input.diff}' -db {params.db} -ip {params.inputdir} -s {params.savedir} -fa {params.fontfamily} "

rule kegg:
    input:
        diff = "{group}-diff-"+hz
    output:
        directory(config['savepath']+"KEGG/{group}/"),
        directory(config['savepath']+"KEGG_map/{group}/")
    threads:3
    params:
        org = config['org'],
        db = config['database'],
        kd = config['keggdb'],
        inputdir = config['inputpath'],
        savedir = config['savepath'],
        fontfamily = config['fontfamily']
    shell:
        "enrichpip -t K -if '{input.diff}' -db {params.db} -kd {params.kd} -o {params.org} -ip {params.inputdir} -s {params.savedir} -fa {params.fontfamily} "

rule interpro:
    input:
        diff = "{group}-diff-"+hz
    output:
        directory(config['savepath']+"InterPro/{group}/")
    params:
        db = config['database'],
        savedir = config['savepath'],
        inputdir = config['inputpath'],
        fontfamily = config['fontfamily']
    shell:
        "enrichpip -t I -if '{input.diff}' -db {params.db} -s {params.savedir} -ip {params.inputdir} -fa {params.fontfamily} "

rule reactome:
    input:
        diff = "{group}-diff-"+hz
    output:
        directory(config['savepath']+"Reactome/{group}/")
    params:
        db = config['database'],
        savedir = config['savepath'],
        inputdir = config['inputpath'],
        fontfamily = config['fontfamily']
    shell:
        "enrichpip -t R -if '{input.diff}' -db {params.db} -s {params.savedir} -ip {params.inputdir} -fa {params.fontfamily} "

rule wikipathway:
    input:
        diff = "{group}-diff-"+hz
    output:
        directory(config['savepath']+"Wikipathways/{group}/")
    params:
        db = config['database'],
        savedir = config['savepath'],
        inputdir = config['inputpath'],
        fontfamily = config['fontfamily']
    shell:
        "enrichpip -t W -if '{input.diff}' -db {params.db} -s {params.savedir} -ip {params.inputdir} -fa {params.fontfamily} "
    
rule ppi:
    input:
        diff = "{group}-diff-"+hz
    output:
        directory(config['ppipath']+"{group}/")
    params:
        db = config['database'],
        path = config['ppipath'],
        topn = config['number'],
        fontfamily = config['fontfamily']
    shell:
        "proppi -if '{input.diff}' -bp {params.db} -n {params.topn} -fa {params.fontfamily} -s {params.path} "

