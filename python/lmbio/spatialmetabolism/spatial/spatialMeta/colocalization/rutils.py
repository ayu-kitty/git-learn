# File: rutils
# Author: jiang tao
# Time: 2023/3/10 13:17
from pathlib import Path

from rpy2.robjects.packages import importr
from rpy2.robjects import r

R_HOME = Path(__file__).parent.parent / 'R'

# r file
r.source(str(R_HOME / 'fig_deal.R'))
r.source(str(R_HOME / 'get_fig_data_for_resnet.R'))
r.source(str(R_HOME / 'spatial_rawdata_metabolite_correlation.R'))
r.source(str(R_HOME / 'spatial_resnetdata_figure_correlation.R'))

# r package
base = importr('base')
lmbio = importr('lmbio')
cardinal = importr('Cardinal')
ebi_image = importr('EBImage')

# r function
to_rdf = r("as.data.frame")
to_matrix = r("as.matrix")
fig_deal = r("fig_deal")
get_fig_data_for_resnet = r("get_fig_data_for_resnet")
spatial_rawdata_metabolite_correlation = r('spatial_rawdata_metacor')
spatial_resnetdata_figure_correlation = r('spatial_resnetdata_figcor')

# r var
NULL = r('as.null()')
