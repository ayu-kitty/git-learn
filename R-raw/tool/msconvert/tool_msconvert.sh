#!/bin/bash

echo '**********' >> runinfo.txt
echo '开始时间:'$(TZ=UTC-8 date +%Y-%m-%d" "%H:%M:%S)  >> runinfo.txt
echo '运行命令:tool_msconvert '$*  >> runinfo.txt
echo '开始时间:'$(TZ=UTC-8 date +%Y-%m-%d" "%H:%M:%S) 
echo '运行命令:msconvert '$*

# 判断是否存在 harbor.oebiotech.com/oe-docker/msconvert 镜像
docker images | grep 'harbor.oebiotech.com/oe-docker/msconvert *latest' &> /dev/null
if [ $? -ne 0 ]
then
  echo "harbor.oebiotech.com/oe-docker/msconvert:latest is not existed,we will docker pull it!!!"
  docker pull harbor.oebiotech.com/oe-docker/msconvert:latest
fi

#判断是否存在 msconvertdeal 容器
docker ps | grep harbor.oebiotech.com/oe-docker/msconvert:latest | grep msconvert_run &> /dev/null

if [ $? -ne 0 ];then

echo "无harbor.oebiotech.com/oe-docker/msconvert:latest相关容器,开启新容器"

docker run -itd \
-v /luming:/luming \
-v /data:/data \
-v /public:/public \
-v /home:/home \
-e WINEDEBUG=-all \
--cpuset-cpus="21-30" \
--privileged=true \
--restart=always \
--name msconvert_run \
harbor.oebiotech.com/oe-docker/msconvert:latest \
bash

fi

docker exec -i -u root:root -w $(pwd) msconvert_run wine msconvert --mzML $*


echo '结束时间:'$(TZ=UTC-8 date +%Y-%m-%d" "%H:%M:%S)  >> runinfo.txt
echo '结束时间:'$(TZ=UTC-8 date +%Y-%m-%d" "%H:%M:%S)
