module mul_score:
    snakefile: "../scoremap/scoremap.smk"
    config: config
    skip_validation: False

use rule * from mul_score as mul_plot_*

module mul_score2:
    snakefile: "../scoremap2/scoremap2.smk"
    config: config
    skip_validation: False

use rule * from mul_score2 as mul_plot_*

module mul_3Dscore:
    snakefile: "../3Dscoremap/3Dscoremap.smk"
    config: config
    skip_validation: False

use rule * from mul_3Dscore as mul_plot_*

module mul_loading:
    snakefile: "../loadingmap/loadingmap.smk"
    config: config
    skip_validation: False

use rule * from mul_loading as mul_plot_*

module mul_loading2:
    snakefile: "../loadingmap2/loadingmap2.smk"
    config: config
    skip_validation: False

use rule * from mul_loading2 as mul_plot_*

module mul_splot:
    snakefile: "../splotmap/splotmap.smk"
    config: config
    skip_validation: False

use rule * from mul_splot as mul_plot_*

module mul_splot2:
    snakefile: "../splotmap2/splotmap2.smk"
    config: config
    skip_validation: False

use rule * from mul_splot2 as mul_plot_*

module mul_permutation:
    snakefile: "../permutation/permutation.smk"
    config: config
    skip_validation: False

use rule * from mul_permutation as mul_plot_*
