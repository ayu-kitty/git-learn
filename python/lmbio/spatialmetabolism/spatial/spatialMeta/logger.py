# File: logger.py
# Author: jiang tao
# Time: 2023/3/10 14:10
import logging
import sys


def _init_logger(logger_name="SpatialMeta"):
    # 实例化logger对象
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.DEBUG)

    # 格式化器
    formatter = logging.Formatter(
        fmt='%(levelname)s - %(asctime)s - %(name)s - %(filename)-8s : line %(lineno)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

    # 将记录器产生的日志输出到屏幕 console
    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(logging.DEBUG)
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    #logger.info('Completed configuring SpatialMeta logger.')