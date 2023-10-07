#!/opt/conda/bin/python

class FilteredPrinter(object):
    def __init__(self, filtered_print, std):
        self._filtered_print = filtered_print
        self._std = std

    def _write(self, string):
        self._filtered_print(string, self._std)

    def __getattr__(self, attr):
        if attr == 'write':
            return self._write
        return getattr(self._std, attr)

def filtered_print(string, std):
    if string != "\n" and not string.__contains__("Select jobs"):
        std.write(string)
        pass

def filtered_print2(string, std):
    if string == "\n" or string.__contains__("oeWarning") or string.__contains__("oeInfo"):
        std.write(string)
        pass
