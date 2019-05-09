import re


def hex2color(hex_color, saturability=False):
    opt = re.findall(r'(.{2})', hex_color[1:])  # 将字符串两两分割
    return [int(i, 16) for i in opt] + ([saturability] if saturability else [])  # 用以存放最后结果


def rgb2hex(rgb):
    return "#" + ''.join([hex(i)[-2:] for j, i in enumerate(rgb) if j < 3]).replace("x", "0")
