#!/opt/conda/bin/python
import requests
# from bs4 import BeautifulSoup
# import pandas as pd
# import numpy as np
import regex as re
import random
# from lxml import etree
import argparse

ua_list = [
    'Mozilla/4.0 (compatible; MSIE 7.0; Windows NT 5.1; Maxthon 2.0',
    'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_7_0) AppleWebKit/535.11 (KHTML, like Gecko) Chrome/17.0.963.56 Safari/535.11',
    'User-Agent:Opera/9.80 (Windows NT 6.1; U; en) Presto/2.8.131 Version/11.11',
    'Mozilla/5.0 (Windows NT 6.1; rv:2.0.1) Gecko/20100101 Firefox/4.0.1',
    'Mozilla/4.0 (compatible; MSIE 7.0; Windows NT 6.0)',
    'Mozilla/5.0 (Windows; U; Windows NT 6.1; en-us) AppleWebKit/534.50 (KHTML, like Gecko) Version/5.1 Safari/534.50',
    'Mozilla/5.0 (compatible; MSIE 9.0; Windows NT 6.1; Trident/5.0',
    'Mozilla/4.0 (compatible; MSIE 7.0; Windows NT 5.1',
    'Mozilla/4.0 (compatible; MSIE 6.0; Windows NT 5.1',
    'Mozilla/5.0 (Macintosh; Intel Mac OS X 10.6; rv:2.0.1) Gecko/20100101 Firefox/4.0.1',
]

def get_html(pid):
    base_url = f"https://alphafold.ebi.ac.uk/api/prediction/{pid}"
    res  = requests.get(base_url,headers ={"User-Agent":random.sample(ua_list,1)[0]})
    match = re.search(r'"pdbUrl":"(.*?)"',res.text)
    return match.group(1) if match else None


def get_pdb(prot,url):
    # with open(prot+".pdb","w") as f:
    with open("protein.pdb","w") as f:
        res = requests.get(url,headers ={"User-Agent":random.sample(ua_list,1)[0]})
        f.write(res.text)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-p","--pid",default = "O43252", help = "蛋白名称")
    args = parser.parse_args()
    
    pid = args.pid
    res = get_html(pid)
    if res:
        get_pdb(pid,res)
        print(pid,"下载完成！")
    else:
        print(pid,"无法下载！")
