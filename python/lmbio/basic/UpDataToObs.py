#!/opt/conda/bin/python
from obs import ObsClient, PutObjectHeader, CompleteMultipartUploadRequest, CompletePart
import hashlib
import datetime
import tarfile
import base64
import traceback
import argparse
import math
import os
import re
import tempfile
import shutil
import sys

inner_print = print

def print(*arg):
    inner_print(*arg)
    inner_print(*arg, file=open("log.txt", "a", encoding="gbk"))


def dfs_showdir(path, 
                filename = "contents.txt",
                depth = 0,
                maxdepth = 2):
    if depth == 0:
        if os.path.exists(filename):
            os.remove(filename)
        inner_print("root:[" + path + "]", file=open(filename, "a", encoding="gbk"))

    for item in os.listdir(path):
        if '.git' not in item:
            inner_print("| " * depth + "|--" + item, file=open(filename, "a", encoding="gbk"))
            newitem = path +'/'+ item
            if depth < maxdepth-1:
                if os.path.isdir(newitem):
                    dfs_showdir(path = newitem, 
                                filename = filename,
                                depth = depth + 1,
                                maxdepth = maxdepth)


def make_targz(output_filename, source_dir,split = True,size = 10240):
    with tarfile.open(output_filename, "w") as tar:
        tar.add(source_dir, arcname=os.path.basename(source_dir))


def make_zip(output_filename, source_dir,split = True,size = "10240m"):
  if split:
    os.system("zip -1 -r -s '%s' '%s' '%s/' > /dev/null 2>&1" % (size,output_filename,source_dir))
  else:
    os.system("zip -1 -r '%s' '%s/' > /dev/null 2>&1" % (output_filename,source_dir))
    

def get_dir_size(dir):
    size = 0
    if not os.path.isfile(dir):
        for root, dirs, files in os.walk(dir):
            size += sum([os.path.getsize(os.path.join(root, name))
                         for name in files])
    else:
        size = os.path.getsize(dir)
    return size


def get_file_md5(file_name, isbyte=False):
    if not os.path.isfile(file_name):
        return
    m = hashlib.md5()  # 创建md5对象
    with open(file_name, 'rb') as fobj:
        while True:
            data = fobj.read(4096)
            if not data:
                break
            m.update(data)  # 更新md5对象
    if isbyte:
        return m.digest()  # 返回md5对象
    else:
        return m.hexdigest()  # 返回md5对象


def getFilemd5(path):
    if os.path.exists(path+"/md5.md5"):
        os.remove(path+"/md5.md5")
    for root, dirs, files in os.walk(path):
        for filename in files:
            if (filename == "md5.md5") | (filename == "upfiletoobs.exe"):
                continue
            filepath = root + os.path.sep + filename
            filepath = re.sub(r"\\", r"/", filepath)
            md5value = get_file_md5(filepath)
            newfilepath = re.sub("^"+path, ".", filepath)
            # print(md5value+" "+newfilepath)
            with open(path+"/md5.md5", "a", encoding="gbk") as f:
                f.write(md5value+" "+newfilepath+"\n")
    if os.listdir(path) == []:
        return False
    else:
        print(path+"目录md5运算完成")
        return True


def getsingleFilemd5(filename,
                     savename="md5.md5"):
    if os.path.exists(savename):
        os.remove(savename)
    md5value = get_file_md5(filename)
    with open(savename, "a", encoding="gbk") as f:
        f.write(md5value+" "+os.path.basename(filename)+"\n")
    print(filename+"目录md5运算完成")


i = 1


def callback(transferredAmount, totalAmount, totalSeconds):
    global i
    if i == 1000:
        # 获取上传平均速率(KB/S)
        uploadrate = transferredAmount * 1.0 / totalSeconds / 1024 / 1024
        uploadrate = str(round(uploadrate, 3))
        inner_print("上传速率:"+uploadrate+"MB/S")
        # 获取上传进度百分比
        uploadprocess = transferredAmount * 100.0 / totalAmount
        uploadprocess = str(round(uploadprocess, 1))+"%"
        inner_print("上传进度:"+uploadprocess)
        i = 1
    else:
        i = i+1


def putFiletoobs(filename,
                 obsClient,
                 bucketName,
                 savepath,
                 maxfilesize=1024*1024*1024):
    headers = PutObjectHeader()
    headers.contentType = 'binary/octet-stream'
    size = get_dir_size(filename)
    if size <= maxfilesize:
        print(filename+"文件进行上传中")
        MD5data = get_file_md5(filename, isbyte=True)
        MD5datab64 = base64.b64encode(MD5data).decode()
        headers.md5 = MD5datab64
        resp = obsClient.putFile(bucketName,
                                 savepath+os.path.basename(filename),
                                 filename,
                                 progressCallback=callback,
                                 headers=headers)
    else:
        print(filename+"文件大于最大上传限制，开启分段上传模式")
        resp = obsClient.initiateMultipartUpload(bucketName,
                                                 savepath +
                                                 os.path.basename(filename),
                                                 contentType='binary/octet-stream')
        if resp.status < 300:
            uploadId = resp.body.uploadId
        else:
            print('errorCode:', resp.errorCode)
            print('errorMessage:', resp.errorMessage)
            return resp
        # 上传段
        num = math.ceil(size/maxfilesize)
        print("将分为"+str(num)+"个分段进行上传:")
        part = []
        for i in range(1, num+1):
            n = 1
            while True:
                # [x] 优化报错
                try:
                    print("第"+str(i)+"个分段上传中")
                    resp = obsClient.uploadPart(bucketName,
                                                savepath+os.path.basename(filename),
                                                i,
                                                uploadId,
                                                filename,
                                                isFile=True,
                                                isAttachMd5=True,
                                                progressCallback=callback,
                                                offset=(i-1)*maxfilesize,
                                                partSize=maxfilesize)
                    if resp.status < 300:
                        part = part+[CompletePart(partNum=i,etag=resp.body.etag)]
                    else:
                        print('errorCode:', resp.errorCode)
                        print('errorMessage:', resp.errorMessage)
                        if n > 3:
                            print("第"+str(i)+"个分段上传错误3次以上,退出上传")
                            return resp
                        n = n+1
                        print("第"+str(i)+"个分段上传错误,重新上传中")
                        continue
                    break
                except:
                    if n > 3:
                        print("第"+str(i)+"个分段上传错误3次以上,退出上传")
                        resp.status = 400
                        return resp
                    n = n+1
                    print("第"+str(i)+"个分段上传错误,重新上传中")
                    continue

        completeMultipartUploadRequest = CompleteMultipartUploadRequest(parts=part)
        resp = obsClient.completeMultipartUpload(bucketName,
                                                 savepath +
                                                 os.path.basename(filename),
                                                 uploadId,
                                                 completeMultipartUploadRequest)
    if resp.status < 300:
        pass
    else:
        print("********"+filename+"上传存在以下问题:********")
        print('errorCode:', resp.errorCode)
        print('errorMessage:', resp.errorMessage)
        
    return resp


def updatatoobs(path,
                filename,
                bucketName="data-lm",
                savepath="111/"+datetime.datetime.now().strftime('%Y')+"/",
                upload=True,
                allmd5=False):
    # 创建ObsClient实例
    # 更新obs账号密码
    obsClient = ObsClient(
        #access_key_id='J1UHQL5QZKS1ZOR1ZGXM',
        access_key_id='0IUFFKNCTM4DQ72C9DV6',
        #secret_access_key='yhMGOwok94mY8gmfjFqQGtmXXWkWe5wvDFcALADH',
        secret_access_key='8GQzIqHetrCbbY8iCBO2Av3ccVNeb3pI45gi5mPV',
        server='obs.cn-east-2.myhuaweicloud.com',
        long_conn_mode=True)
    try:
        print(path+"目录上传开始，上传时间为" +
              datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))

        if os.path.isfile(path):
            size = get_dir_size(path)
            print('文件大小为: %.3f Gb' % (size/1024/1024/1024))
            tmpdirname = tempfile.mkdtemp(dir="/data/hstore1/updata")
            print('创建临时目录:'+tmpdirname)
            tarfilename = path
            md5file = tmpdirname+"/"+"md5.md5"
            getsingleFilemd5(filename=path,
                             savename=md5file)
            if upload:
              resp = putFiletoobs(filename=tarfilename,
                                  obsClient=obsClient,
                                  bucketName=bucketName,
                                  savepath=savepath+filename+"/")
        else:
            if allmd5:
                getFilemd5(path)
            size = get_dir_size(path)
            print('文件大小为: %.3f Gb' % (size/1024/1024/1024))
            tmpdirname = tempfile.mkdtemp(dir="/data/hstore1/updata")
            print('创建临时目录:'+tmpdirname)
            contentsfile = tmpdirname+"/contents.txt"
            print('形成目录文件:'+contentsfile)
            dfs_showdir(path = path, 
                        filename = contentsfile)
            # tarfilename = tmpdirname+"/"+filename+".tar"
            tarfilename = tmpdirname+"/"+filename+".zip"
            print(path+"目录正压缩为"+tarfilename)
            # make_targz(tarfilename, path)
            # size = get_dir_size(tarfilename)
            make_zip(tarfilename, path)
            size = get_dir_size(tmpdirname)
            print('压缩后文件大小为: %.3f Gb' % (size/1024/1024/1024))
            getFilemd5(tmpdirname)
            # md5file = tmpdirname+"/"+"md5.md5"
            # getsingleFilemd5(filename=tarfilename,
            #                  savename=md5file)
            # resp = putFiletoobs(filename=contentsfile,
            #                     obsClient=obsClient,
            #                     bucketName=bucketName,
            #                     savepath=savepath+filename+"/")
            
        if upload:
            print("上传位置:"+bucketName+"/"+savepath+filename)
            for file in os.listdir(tmpdirname):
              # print(tmpdirname+"/"+file)
              # resp = putFiletoobs(filename=md5file,
              #                     obsClient=obsClient,
              #                     bucketName=bucketName,
              #                     savepath=savepath+filename+"/")
              # resp = putFiletoobs(filename=tarfilename,
              #                     obsClient=obsClient,
              #                     bucketName=bucketName,
              #                     savepath=savepath+filename+"/")
              resp = putFiletoobs(filename=tmpdirname+"/"+file,
                                  obsClient=obsClient,
                                  bucketName=bucketName,
                                  savepath=savepath+filename+"/")
            print(path+"目录上传完成，结束时间为" +
                  datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
            shutil.rmtree(tmpdirname)
        else:
            print("数据不进行上传")
            print('数据在临时目录，请进行查看:'+tmpdirname)
    except:
        print(traceback.format_exc())
        print(path+"目录上传错误，结束时间为" +
              datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
        obsClient.close()
    # 关闭obsClient
    obsClient.close()


def md5andup(path,
             filename,
             bucketName="data-lm",
             savepath="111/"+datetime.datetime.now().strftime('%Y')+"/",
             upload=True,
             allmd5=False):
    print("当前运行目录："+path)
    # if getFilemd5(path):
    if re.match(r"\$", os.path.basename(path)):
        print(path+"为隐藏文件")
        return
    if re.match(r"\.", os.path.basename(path)):
        print(path+"为隐藏文件")
        return
    if os.path.isfile(path) or os.listdir(path) != []:
        updatatoobs(path=path,
                    filename=filename,
                    bucketName=bucketName,
                    savepath=savepath,
                    upload=upload,
                    allmd5=allmd5)
    else:
        print(path+"目录为空跳过运行")
        return


def updataprocess(path,
                  savepath,
                  filename,
                  bucketName="data-lm",
                  upload=True,
                  allmd5=False):
    savepath = re.sub("^/", "", savepath)
    print("-----------------")
    if savepath == "":
        savepath = "test/"
        print("路径输入为空，当前保存路径为："+savepath)
    else:
        savepath = savepath+"/"
        savepath = savepath.replace("//", "/")
        print("当前保存路径为："+savepath)
    md5andup(path=path,
             filename=filename,
             bucketName=bucketName,
             savepath=savepath,
             upload=upload,
             allmd5=allmd5)


def mulupdataprocess(savepath="",
                     bucketName="data-lm",
                     upload=True,
                     allmd5=False,
                     runpath="./"):
    if savepath == "":
        savepath = input("请输入保存路径:")

    if os.path.isdir(runpath):
        os.chdir(runpath)
        for path in os.listdir():
            if os.path.isdir(path):
                updataprocess(path=path,
                              filename=path,
                              bucketName=bucketName,
                              savepath=savepath,
                              upload=upload,
                              allmd5=allmd5)
            else:
                print("-----------------")
                print(path+"不是目录，跳过运行")
    else:
        print(runpath+"，未找到此目录")


def auto_updataprocess(args):
    if args.manual:
        updataprocess(path=args.path,
                      savepath=args.savepath,
                      filename=args.filename,
                      upload=args.upload,
                      allmd5=args.allmd5,
                      bucketName=args.bucketName)
    else:
        mulupdataprocess(allmd5=args.allmd5,
                         upload=args.upload,
                         savepath=args.savepath,
                         runpath=args.runpath,
                         bucketName=args.bucketName)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-u', '--upload', action='store_false',
                        default=True, help="是否上传文件")
    parser.add_argument('-a', '--allmd5', action='store_true',
                        default=False, help="是否针对所有文件进行md5运行后压缩")
    parser.add_argument('-s', '--savepath', type=str,
                        default="", help="指定上传的云端路径")
    parser.add_argument('-r', '--runpath', type=str,
                        default="./", help="指定工作目录")
    parser.add_argument('-m', '--manual', action='store_true',
                        default=False, help="使用后开启指定目录上传")
    parser.add_argument('-p', '--path', type=str,
                        default="test1", help="-m开启后，指定上传的本地路径")
    parser.add_argument('-f', '--filename', type=str,
                        default="test", help="-m开启后，指定上传后的名称")
    parser.add_argument('-b', '--bucketName', type=str,
                        default="data-lm", help="指定上传的桶") 
    args = parser.parse_args()
    auto_updataprocess(args)
