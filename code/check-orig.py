#%%
from astropy.io import fits as pf
from astropy.time import Time
from datetime import timedelta

# 读取 FITS 文件
file = "../data_orig/3C138/20241207/3C138_multibeamcalibration-M01_W_0001.fits"
hdul = pf.open(file)
data = pf.getdata(file)
utobs = data['UTOBS']

# MJD 转换为 UTC
t = Time(utobs, format='mjd')
ut = t.utc.iso  # ISO 格式 UTC 时间

# 转换为北京时间（UTC+8）
bj_time = (t + timedelta(hours=8)).iso  # 直接加 8 小时

# 打印结果
print("MJD:", utobs)
print("UTC Time:", ut)
print("Beijing Time (UTC+8):", bj_time)
# %%
