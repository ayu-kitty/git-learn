# File: exception
# Author: jiang tao
# Time: 2023/7/12 9:48
class TcmException(Exception):
    pass


class DatabaseNotImplementError(TcmException):
    pass


class SaveOnlineRegisterTableError(TcmException):
    pass


class SqlExecuteError(TcmException):
    pass


class RawDataNotFoundError(TcmException):
    """no .raw format files found in rawdata folder"""
    pass


class ShellRunTimeError(TcmException):
    """error caused by sp.run """
    pass