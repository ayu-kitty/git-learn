if config["params"]["omic"] == "M" or config["params"]["omic"] == "ML":
  rule get_plot_qcsampletree:
      input:
          config["path"]["datafile"]
      output:
          expand("{{savepath}}/SamplesTreeplot.{imagetype}",
                 imagetype=config["params"]["imagetype"])
      params:
          family = config["params"]["plot"]["family"],
          width = config["params"]["plot"]["width"],
          height = config["params"]["plot"]["height"],
          classfile = config["path"]["classtypefile"],
          imagetype = config["params"]["imagetype"]
      threads: 1
      shell:
          '''
          procre_corrtree \
          --inputfile '{input}' \
          --savepath '{wildcards.savepath}' \
          --imagetype {params.imagetype} \
          --fontfamily '{params.family}' \
          --width {params.width} \
          --height {params.height} \
          --classfile {params.classfile} \
          --scaledat T
          '''
else:
  rule get_plot_qcsampletree:
      input:
          config["path"]["datafile"]
      output:
          expand("{{savepath}}/SamplesTreeplot.{imagetype}",
                 imagetype=config["params"]["imagetype"])
      params:
          family = config["params"]["plot"]["family"],
          width = config["params"]["plot"]["width"],
          height = config["params"]["plot"]["height"],
          classfile = config["path"]["classtypefile"],
          imagetype = config["params"]["imagetype"]
      threads: 1
      shell:
          '''
          procre_corrtree \
          --inputfile '{input}' \
          --savepath '{wildcards.savepath}' \
          --imagetype {params.imagetype} \
          --fontfamily '{params.family}' \
          --width {params.width} \
          --height {params.height} \
          --classfile {params.classfile}
          '''

rule get_plot_qcboxplot:
    input:
        config["path"]["datafile"]
    output:
        expand("{{savepath}}/Boxplot.{imagetype}",
               imagetype=config["params"]["imagetype"])
    params:
        family = config["params"]["plot"]["family"],
        width = config["params"]["plot"]["width"],
        height = config["params"]["plot"]["height"],
        classfile = config["path"]["classtypefile"],
        imagetype = config["params"]["imagetype"]
    threads: 1
    shell:
        '''
        procre_boxplot \
        --inputfile '{input}' \
        --savepath '{wildcards.savepath}' \
        --imagetype {params.imagetype} \
        --fontfamily '{params.family}' \
        --width {params.width} \
        --height {params.height} \
        --classfile {params.classfile} \
        --log
        '''

if config["params"]["omic"] == "M" or config["params"]["omic"] == "ML":
  rule get_plot_qccor:
      input:
          config["path"]["datafile"]
      output:
          expand("{{savepath}}/QCcorrplot.{imagetype}",
                 imagetype=config["params"]["imagetype"])
      params:
          family = config["params"]["plot"]["family"],
          width = config["params"]["plot"]["width"],
          height = config["params"]["plot"]["height"],
          classfile = config["path"]["classtypefile"],
          imagetype = config["params"]["imagetype"]
      threads: 1
      shell:
          '''
          map_common_corrplot4 \
          --filename '{input}' \
          --savepath '{wildcards.savepath}' \
          --mapname 'QCcorrplot' \
          --imagetype {params.imagetype} \
          --family '{params.family}' \
          --width {params.width} \
          --height {params.height} \
          --classfile {params.classfile}
          '''
else:
  rule get_plot_qccor:
    input:
        config["path"]["datafile"]
    output:
        expand("{{savepath}}/QCcorrplot.{imagetype}",
               imagetype=config["params"]["imagetype"])
    params:
        family = config["params"]["plot"]["family"],
        width = config["params"]["plot"]["width"],
        height = config["params"]["plot"]["height"],
        classfile = config["path"]["classtypefile"],
        imagetype = config["params"]["imagetype"]
    threads: 1
    shell:
        '''
        procre_samcorr \
        --inputfile '{input}' \
        --savepath '{wildcards.savepath}' \
        --mapname 'QCcorrplot' \
        --imagetype {params.imagetype} \
        --fontfamily '{params.family}' \
        --width {params.width} \
        --height {params.height} \
        --classfile {params.classfile}
        '''


rule get_plot_samplecor:
  input:
      config["path"]["datafile"]
  output:
      expand("{{savepath}}/Samplecorrplot.{imagetype}",
             imagetype=config["params"]["imagetype"])
  params:
      family = config["params"]["plot"]["family"],
      width = config["params"]["plot"]["width"],
      height = config["params"]["plot"]["height"],
      classfile = config["path"]["classtypefile"],
      imagetype = config["params"]["imagetype"]
  threads: 1
  shell:
      '''
      procre_samcorr \
      --inputfile '{input}' \
      --savepath '{wildcards.savepath}' \
      --mapname 'Samplecorrplot' \
      --imagetype {params.imagetype} \
      --fontfamily '{params.family}' \
      --width {params.width} \
      --height {params.height} \
      --classfile {params.classfile}
      '''

rule get_plot_rsd:
  input:
      config["path"]["datafile"]
  output:
      expand("{{savepath}}/RSD-boxplot.{imagetype}",
             imagetype=config["params"]["imagetype"])
  params:
      family = config["params"]["plot"]["family"],
      width = config["params"]["plot"]["width"],
      height = config["params"]["plot"]["height"],
      classfile = config["path"]["classtypefile"],
      imagetype = config["params"]["imagetype"]
  threads: 1
  shell:
      '''
      procre_rsd \
      --inputfile '{input}' \
      --savepath '{wildcards.savepath}' \
      --imagetype {params.imagetype} \
      --fontfamily '{params.family}' \
      --width {params.width} \
      --height {params.height} \
      --classfile {params.classfile}
      '''

rule get_plot_density:
  input:
      config["path"]["datafile"]
  output:
      expand("{{savepath}}/Densityplot.{imagetype}",
             imagetype=config["params"]["imagetype"])
  params:
      family = config["params"]["plot"]["family"],
      width = config["params"]["plot"]["width"],
      height = config["params"]["plot"]["height"],
      classfile = config["path"]["classtypefile"],
      imagetype = config["params"]["imagetype"]
  threads: 1
  shell:
      '''
      procre_density \
      --inputfile '{input}' \
      --savepath '{wildcards.savepath}' \
      --imagetype {params.imagetype} \
      --fontfamily '{params.family}' \
      --width {params.width} \
      --height {params.height} \
      --classfile {params.classfile}
      '''
           