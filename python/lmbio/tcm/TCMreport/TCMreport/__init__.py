name = "TCMreport"
import warnings
from .__version__ import __version__
from .QI2report import Report
from .Blood2report import BloodReport
from .ProductResolve import ProductResolve
from .logger import _init_logger

_init_logger()
# silent pandas warnings
warnings.filterwarnings(action="ignore")

__author__ = "jiangtao"
__email__ = "jiangt@lumingbio.com"

__all__ = [
    "__version__",
    "Report",
    "BloodReport",
    "ProductResolve"
]
