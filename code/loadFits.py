#!/usr/bin/env python
# coding: utf-8

import os
import subprocess
import time
from glob import glob
from natsort import natsorted
import numpy as np
from astropy.io import fits as pf
from copy import deepcopy
from astropy.table import Table
chans_all = 65536

def load_fits(file_in, dir_out, avechans=64):
    """
    读取 FITS 文件，并对数据进行通道平均后保存。
    """
    source = file_in.split("/")[-3]
    obs_date = file_in.split("/")[-2]
    file_new = file_in.split("/")[-1].replace("_W_0001.fits", "_ave-64.fits")
    dir_new = f"{dir_out}/{source}_{obs_date}"
    os.makedirs(dir_new, exist_ok=True)
    realpath_file = f"{dir_new}/{file_new}"
    
    file_count = len(glob(file_in.replace("0001.fits", "*.fits")))
    obs_data_all = np.array([])
    obs_time_all = np.array([])
    print(f"Number of files: {file_count}")
    
    for num in range(1, file_count + 1):
        file_now = file_in.replace("0001.fits", f"000{num}.fits")
        print(f"Reading file: {file_now}")
        data = pf.getdata(file_now)

        if num ==  1:
            hdr = pf.getheader(file_now)
            hdr_new = deepcopy(hdr)
            hdr_new["FREQ0"] = data['FREQ'][0] + data['CHAN_BW'][0] * avechans / 2
            hdr_new["CHAN_BW"] = data['CHAN_BW'][0] * avechans
            
        with pf.open(file_now) as hdul:
            obs_data = hdul[1].data['DATA']
            obs_data = np.median(obs_data.reshape(obs_data.shape[0], chans_all // avechans, avechans, 4), axis=2)
            print(f"Now the shape is {obs_data.shape}")
            obs_time = hdul[1].data['UTOBS']
            obs_data_all = np.concatenate((obs_data_all, obs_data), axis=0) if obs_data_all.size else obs_data
            obs_time_all = np.concatenate((obs_time_all, obs_time), axis=0) if obs_time_all.size else obs_time
    # 创建一个空的hdr
    
    print(f"After combining: obs {obs_data_all.shape}, time {obs_time_all.shape}")
    table = Table()
    table['DATA'] = obs_data_all
    table['UTOBS'] = obs_time_all

    try:
        # 保存为fitstable,obs_data_all,obs_time_all
        hdu = pf.BinTableHDU(table, header=hdr_new) 
        hdu.writeto(realpath_file, overwrite=True)
        print(f"File saved: {realpath_file}")
    except Exception as e:
        print(f"Error: {e}")
    
    return realpath_file

def submit_sbatch(file_in, dir_out):
    """
    生成 sbatch 脚本，并提交到 SLURM。
    """
    file_base = file_in.split("/")[-2] + "_" + file_in.split("/")[-1].replace("_W_0001.fits", "")
    script_dir = "sbatch_scripts"
    log_dir = "logs"
    
    os.makedirs(script_dir, exist_ok=True)
    os.makedirs(log_dir, exist_ok=True)

    script_path = f"{script_dir}/{file_base}.sbatch"
    output_log = f"{log_dir}/{file_base}_%j.out"
    error_log = f"{log_dir}/{file_base}_%j.err"
    
    sbatch_script = f"""#!/bin/bash
#SBATCH --job-name={file_base}
#SBATCH --output={output_log}

python -c "from loadFits import load_fits; load_fits('{file_in}', '{dir_out}', avechans=64)" > {log_dir}/{file_base}.log"""

    with open(script_path, "w") as f:
        f.write(sbatch_script)

    print(f"Submitting job: {script_path}")
    subprocess.run(["sbatch", script_path], check=True)  # 阻塞式提交，确保提交成功
    time.sleep(0.5)  # 每次提交任务后休眠 1s，避免任务拥堵

def main():
    """
    主函数，遍历文件列表并提交任务。
    """
    dir_data_orig = "../data_orig"
    file_list = natsorted(glob(f"{dir_data_orig}/*/*/*_W_0001.fits"))

    file_list = file_list
    print(f"Total files to submit: {len(file_list)}")  # 确保只有 3 个任务

    dir_out = "../data_ave-64"
    os.makedirs(dir_out, exist_ok=True)

    for file_in in file_list:
        submit_sbatch(file_in, dir_out)

if __name__ == "__main__":
    main()