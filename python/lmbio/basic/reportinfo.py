import os

def oe_weclome(sys_logo = "LM"):
     if sys_logo == "LM":
         return '''

####  感谢您选择鹿明生物！

**关于网页报告使用有以下几点提示:**

1. **报告正文**中的**蓝色字体**均可以**点击**快速导引到感兴趣的章节, 便于您快速查看相应内容。
2. **报告正文**中的**图片**均可以点击后进行**放大查看**, 且图片点击后可以包含更多的细节信息, 比如左上角会显示具体的**图片数量**, 右上角有相关的**图片工具**, **键盘上的左右方向键/鼠标滚轮**可以对图片进行**切换**。
3. **报告正文**中**每个区域**中图片的数量**最多**只显示2张, 更多的图片可以在点击之后在弹出界面利用**键盘方向键/鼠标滚轮**查看。
4. 展示的表格均可以在各列**表头**进行**升序或者降序显示**, 当表格多于20行的时候, 表格会嵌入到网页之中, 可以利用滚动栏进行调整设置, 同时会在表格上方增加搜索设置, 便于快速查询表格信息。
5. 在报告中浏览, 需要返回顶部的时候, 可以使用网页**右下角的白色箭头标识**快速返回顶部。
6. 本提示可以点击右上方的不再显示, 则后续打开网页均不会显示本提示。
         '''
     elif sys_logo == "OE" or sys_logo == "QD":
         return '''

####  感谢您选择欧易生物！

**关于网页报告使用有以下几点提示:**

1. **报告正文**中的**蓝色字体**均可以**点击**快速导引到感兴趣的章节, 便于您快速查看相应内容。
2. **报告正文**中的**图片**均可以点击后进行**放大查看**, 且图片点击后可以包含更多的细节信息, 比如左上角会显示具体的**图片数量**, 右上角有相关的**图片工具**, **键盘上的左右方向键/鼠标滚轮**可以对图片进行**切换**。
3. **报告正文**中**每个区域**中图片的数量**最多**只显示2张, 更多的图片可以在点击之后在弹出界面利用**键盘方向键/鼠标滚轮**查看。
4. 展示的表格均可以在各列**表头**进行**升序或者降序显示**, 当表格多于20行的时候, 表格会嵌入到网页之中, 可以利用滚动栏进行调整设置, 同时会在表格上方增加搜索设置, 便于快速查询表格信息。
5. 在报告中浏览, 需要返回顶部的时候, 可以使用网页**右下角的白色箭头标识**快速返回顶部。
6. 本提示可以点击右上方的不再显示, 则后续打开网页均不会显示本提示。
         '''
     else:
         return '''

**关于网页报告使用有以下几点提示:**

1. **报告正文**中的**蓝色字体**均可以**点击**快速导引到感兴趣的章节, 便于您快速查看相应内容。
2. **报告正文**中的**图片**均可以点击后进行**放大查看**, 且图片点击后可以包含更多的细节信息, 比如左上角会显示具体的**图片数量**, 右上角有相关的**图片工具**, **键盘上的左右方向键/鼠标滚轮**可以对图片进行**切换**。
3. **报告正文**中**每个区域**中图片的数量**最多**只显示2张, 更多的图片可以在点击之后在弹出界面利用**键盘方向键/鼠标滚轮**查看。
4. 展示的表格均可以在各列**表头**进行**升序或者降序显示**, 当表格多于20行的时候, 表格会嵌入到网页之中, 可以利用滚动栏进行调整设置, 同时会在表格上方增加搜索设置, 便于快速查询表格信息。
5. 在报告中浏览, 需要返回顶部的时候, 可以使用网页**右下角的白色箭头标识**快速返回顶部。
6. 本提示可以点击右上方的不再显示, 则后续打开网页均不会显示本提示。
         '''

def companyinfo(report,
                sys_logo = "LM",
                yamlpath = "./"):
  
  if sys_logo.lower() == "lm":
      report.add_yaml_config(os.path.join(yamlpath,"yamlfiles/lm.yaml"))
      company = report.add_section('公司简介')
      company.add_plot('LM1*',caption='上海鹿明生物科技有限公司资质')
      company_instrument= company.add_section('公司仪器')
      company_instrument.add_plot('LM2.png',caption='上海鹿明生物科技有限公司仪器平台')
      company_vision= company.add_section('公司服务')
      company_vision.add_plot('LM3.png',caption='上海鹿明生物科技有限公司服务')
      declare = report.add_section('申明')
  elif sys_logo.lower() == "sg":
      report.add_yaml_config(os.path.join(yamlpath,"yamlfiles/sg.yaml"))
  elif sys_logo.lower() == "gy":
      report.add_yaml_config(os.path.join(yamlpath,"yamlfiles/gy.yaml"))
  elif sys_logo.lower() == "hy":
      report.add_yaml_config(os.path.join(yamlpath,"yamlfiles/hy.yaml"))
  elif sys_logo.lower() == "yz":
      report.add_yaml_config(os.path.join(yamlpath,"yamlfiles/yz.yaml"))
  elif sys_logo.lower() == "unlogo":
      report.add_yaml_config(os.path.join(yamlpath,"yamlfiles/unlogo.yaml"))
  elif sys_logo.lower() == "qds":
      report.add_yaml_config(os.path.join(yamlpath,"yamlfiles/qds.yaml"))
  elif sys_logo.lower() == "qd":
      report.add_yaml_config(os.path.join(yamlpath,"yamlfiles/qd.yaml"))
  elif sys_logo.lower() == "wnf":
      report.add_yaml_config(os.path.join(yamlpath,"yamlfiles/wnf.yaml"))
  else:
      report.add_yaml_config(os.path.join(yamlpath,"yamlfiles/oe.yaml"))
      declare = report.add_section('申明')

class getreportinfo:
    def __init__(self,sys_logo = "LM"):
        self.sys_logo = sys_logo
        self.oe_weclome = oe_weclome(sys_logo = sys_logo)
    def companyinfo(self,report,yamlpath = "/data/hstore1/database/report/2023-03-04~2023/src/"):
        return companyinfo(report = report,
                           sys_logo = self.sys_logo,
                           yamlpath = yamlpath)
