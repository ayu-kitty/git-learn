#!/opt/conda/bin/python
import paramiko

def windows_updataprocess():
    savepath = input("请输入保存路径:")
    datafromip = input("请输入NAS的ip地址:")
    runpath = input("请输入NAS上目录路径:")
    sshupdata(datafromip = datafromip,
              runpath = runpath,
              savepath = savepath)

def sshupdata(datafromip = "192.168.10.177",
              runpath = "",
              savepath = ""):

    if datafromip == "192.168.10.172":
        datapath = "/data/nas/172/"+runpath
    elif datafromip == "192.168.10.173":
        datapath = "/data/nas/173/"+runpath
    elif datafromip == "192.168.10.174":
        datapath = "/data/nas/174/"+runpath
    elif datafromip == "192.168.10.175":
        datapath = "/data/nas/175/"+runpath
    elif datafromip == "192.168.10.176":
        datapath = "/data/nas/176/"+runpath
    elif datafromip == "192.168.10.177":
        datapath = "/data/nas/177/"+runpath
    elif datafromip == "192.168.10.179":
        datapath = "/data/nas/179/"+runpath
    elif datafromip == "10.100.10.42":
        datapath = runpath
    else:
        print("NAS的ip地址不正确")
        savepath = input("请输入回车结束本次运行")
        return

    print("服务器上的目录为:"+datapath)

    splitdatapath = datapath.split("/")
    if len(splitdatapath) < 6:
        print("上传层次太低，怕误操作不进行上传")
        savepath = input("请输入回车结束本次运行")
        return

    ssh = paramiko.SSHClient() 
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect("10.100.10.42",
                "8011",
                "jiawei.lu",
                "Oebio@123")
    order = "tool_UpDataToObs " + "-r " + datapath + " -s "+savepath
    
    print("运行命令为:"+order)

    try:
        stdin,stdout,stderr = ssh.exec_command(order,get_pty=True)

        # 打印输出结果
        while not stdout.channel.exit_status_ready():
            result = stdout.readline()
            print(result,end="")
            # 由于在退出时，stdout还是会有一次输出，因此需要单独处理，处理完之后，就可以跳出了
            if stdout.channel.exit_status_ready():
                result = stdout.readlines()
                print(result)
                break

        # 打印命令执行错误信息
        for i in stderr.readlines():
            print(i,end="")
    except:
        ssh.close()

    ssh.close()

    savepath = input("请输入回车结束本次运行")

if __name__ == '__main__':
    windows_updataprocess()
