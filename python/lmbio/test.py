import re


def test():
  return 1

# import matplotlib.pyplot as plt
#
#
# fig, axes = plt.subplots(2, 1, figsize=(8, 7.4), dpi=100, gridspec_kw={'height_ratios': [1, 2.7]})
#
# ax1 = axes[0]
# ax2 = axes[1]
# print(ax2.get_subplotspec())
#
# print(ax1.get_subplotspec())
def filename_handler(filename: str, escape: bool = True) -> str:
  """
  :param filename: 文件名
  :param escape: 是否转义特殊字符
  :return: new name
  """
  _name = "".join(re.findall(r'[^\*"/:?\\|<>]', filename, re.S)) if escape else filename
  # 防止cpd_name过长而替换产物的xic图
  if _name.__contains__("_"):
    p = _name.split("_")[1]
    _name = _name.split("_")[0]
    if len(_name) > 40:
      # 13 ~ _Product_M2-1
      _name = _name[:30] + ".."
    _name = _name + "_" + p
  else:
    if len(_name) > 40:
      _name = _name[:37] + "..."
  return _name

print(filename_handler("Argininyl-fructosyl-glucose_qt"))