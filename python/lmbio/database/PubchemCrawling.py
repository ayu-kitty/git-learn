#!/opt/conda/bin/python
import requests
# from bs4 import BeautifulSoup
import pandas as pd
import numpy as np
import regex as re
import random
from lxml import etree
from fake_useragent import UserAgent
import csv
from queue import Queue
from threading import Lock
from threading import Thread
import warnings

#将url压入queue中
#压入队列操作在出对列同时进行
#每个线程关联的函数会对同一文件各自写入内容

class pubchem_spyder:
    def __init__(self):
        self.proxypool_url = 'http://192.168.10.200:5555/random'
        self.base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/xml"
        self.header = ("cid","cas","cname","fml","mw","iupac","inchi","inchikey","smiles","lipidmapid")
        self.q = Queue()
        self.lock =Lock()
    
    #将1.6亿条拆100个区间，逐个区间下载
    def split_interval(self):
        j  = 1
        space = int(1.65e6)
        # space = int(10)
        start = [j+space*i for i in range(100)]
        stop = [j+space*i for i in range(1,100)]
        stop.append(163500000)
        # stop.append(31)
        intervals = list(zip(start,stop))

        return intervals

    def get_random_proxy(self):
        #get random proxy from proxypool return proxy_ip
        return requests.get(self.proxypool_url).text.strip()
    
    def get_html(self,url):
        # get_html文档
        headers = {"User-Agent":UserAgent().random}
        proxy = self.get_random_proxy()
        proxies = {'http': 'http://' + proxy,'https':'http://'+proxy}
        # print('get random proxy', proxy)
        try:
            xml = requests.get(url,headers = headers,proxies = proxies,timeout = 2,verify = False).text
        except Exception as e:
            xml = self.get_html(url)
            return xml
        else:
            return xml

    def crawl(self,interval,fname=None):
        i = 0
        while not self.q.empty():
            url = self.q.get()
            # print(url,i)
            xml = self.get_html(url)
            if re.search("no record found",xml,re.I):
                pass 
            else:
                xml = re.sub(r' encoding="UTF-8"',"",xml)
                print(url,"downloaded!")
                self.parse_html(xml,fname)
                # i += 1

    def url_in(self,start,stop):
        for cid in range(start,stop):
            url = self.base_url.format(cid = cid)
            self.q.put(url)


    def parse_html(self,xml,fname):
        html = etree.HTML(xml)
        #cid & title
        cid_n = html.xpath("//recordnumber/text()")
        cid = cid_n[0] if cid_n else ""
        cname_n = html.xpath("//recordtitle/text()")
        cname = cname_n[0] if cname_n else ""
       
        #depositor synonyms
        # syns_n = html.xpath('//section/tocheading[text()= "Depositor-Supplied Synonyms"]')
        # syns_l = syns_n[0].xpath('..//value//string/text()') if syns_n else []
        # syns = ";".join(syns_l)

        #iupac Name
        iupac_n = html.xpath('//section/tocheading[text()="IUPAC Name"]')
        iupac = iupac_n[0].xpath('..//string')[0].text if iupac_n else ""

        #INCHI
        inchi_n = html.xpath('//section/tocheading[text()="InChI"]')
        inchi = inchi_n[0].xpath('..//string')[0].text if inchi_n else ""

        #SMILES
        smiles_n = html.xpath('//section/tocheading[text()="Canonical SMILES"]')
        smiles = smiles_n[0].xpath('..//string')[0].text if smiles_n else ""

        #formula
        fml_n = html.xpath('//section/tocheading[text()="Molecular Formula"]')
        fml = fml_n[0].xpath('..//string')[0].text if fml_n else ""

        #CAS
        cas_n = html.xpath('//section/tocheading[text()="CAS"]')
        cas = cas_n[0].xpath('..//string')[0].text if cas_n else ""

        #CAS
        cas_n = html.xpath('//section/tocheading[text()="CAS"]')
        cas = cas_n[0].xpath('..//string')[0].text if cas_n else ""

        #MW
        mw_n = html.xpath('//section/tocheading[text()="Molecular Weight"]')
        mw = mw_n[0].xpath('..//string')[0].text if mw_n else ""

        #Inchikey
        inchikey_n = html.xpath('//section/tocheading[text()="InChI Key"]')
        inchikey = inchikey_n[0].xpath('..//string')[0].text if inchikey_n else ""

        #lipid_map_source
        lipidmap_n = html.xpath('//reference/sourcename[text()="LIPID MAPS"]')
        lipidmapid = lipidmap_n[0].xpath('..//sourceid')[0].text if lipidmap_n else ""

        res = (cid,cas,cname,fml,mw,iupac,inchi,inchikey,smiles,lipidmapid)
        # print()
        self.save(res,fname)

    def write_header(self,fname=None):
        with open(fname,'w',newline = '',encoding = "utf-8") as f:
            writer = csv.writer(f,dialect = "excel-tab")
            self.lock.acquire()
            writer.writerow(self.header)
            self.lock.release()
        return

    def save(self,row,fname=None):
        with open(fname,'a',newline = '',encoding = "utf-8") as f:
            writer = csv.writer(f,dialect = "excel-tab")
            self.lock.acquire()
            writer.writerow(row)
            self.lock.release()
        return

    def main(self,t,interval,fname):
        self.url_in(interval[0],interval[1])
        t_list = []
        for th in range(t):
            thr = Thread(target = self.crawl,args = (interval,fname))
            t_list.append(thr)
            thr.start()
        for th in t_list:
            th.join()
 

if __name__ == "__main__":
    warnings.filterwarnings('ignore')
    n_thd = 10
    spyder = pubchem_spyder()
    intervals = spyder.split_interval()

    for cid_interval in intervals:
        print(cid_interval)
        start,stop = cid_interval
        fname = f"{start}-{stop-1}.txt"
        spyder.write_header(fname)
        spyder.main(n_thd,cid_interval,fname)
