#!/bin/bash

# 获取当前工作目录路径
path=$(pwd)
SHELLFOLDER=$(dirname $(readlink -f "$0"))
imagename="harbor.oebiotech.com/oe-docker/rawtoimzml"

rand=$(sed "s/[^a-zA-Z0-9]//g" <<< $(openssl rand -base64 5))
echo '转格式容器编号:'$rand
echo '**********' >> runinfo.txt 
echo '转格式容器编号:'$rand >> runinfo.txt
echo '转格式镜像:'$imagename >> runinfo.txt
echo '转格式命令:RawToimzML.sh' >> runinfo.txt
echo '转格式用户:'$(whoami) >> runinfo.txt
echo '转格式开始时间:'$(TZ=UTC-8 date +%Y-%m-%d" "%H:%M:%S) >> runinfo.txt 

if [ -d "./raw" ]; then
docker run -i --rm \
-v "$path/raw":/mydata \
--cpuset-cpus="1-5" \
-m 6000M \
--memory-swap 6000M \
--name $rand \
$imagename \
RawToimzML.sh
else
echo "raw目录不存在,不进行转格式"
fi

echo '转格式结束时间:'$(TZ=UTC-8 date +%Y-%m-%d" "%H:%M:%S) >> runinfo.txt 
