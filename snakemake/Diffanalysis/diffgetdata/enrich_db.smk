
if "infofile" not in config["params"]["enrich"]:
  config["params"]["enrich"]["infofile"]="F"

rule enrich_db:
    output:
        directory(config["params"]["enrich"]["database"])
    threads: 1
    params:
        dbfrom=config["params"]["enrich"]["databasefrom"],
        org=config["params"]["enrich"]["org"],
        keggdbfrom=config["params"]["enrich"]["keggdb"],
        infofile=config["params"]["enrich"]["infofile"]
    shell:
        '''
        /opt/conda/bin/Rscript -e "lmbio::movebackground( \
        db = '{output}', \
        keggmapfrom = '{params.keggdbfrom}', \
        dbfrom = '{params.dbfrom}', \
        org = '{params.org}', \
        inputfile = '{params.infofile}' )"
        '''

rule get_oecloud_bg:
    input:
        bg = config["params"]["enrich"]["database"]
    output:
        zipbg = "{savepath}/OECloud_tools/background.oecloud"
    shell:
        '''
        tar -cf '{output.zipbg}' --directory='{input.bg}' .
        '''
