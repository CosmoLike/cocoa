"""
Utility script for compressing big text file to .npz. Run this in the root dir.
"""

import numpy as np

data_dir = "fortran/likes2020/actpollite_dr4_v4/data"
bp15 = np.genfromtxt(f"{data_dir}/coadd_bpwf_15mJy_191127_lmin2.txt", delimiter=None)
bp100 = np.genfromtxt(f"{data_dir}/coadd_bpwf_100mJy_191127_lmin2.txt", delimiter=None)

np.savez_compressed("data/coadd_bpwf_15mJy_191127_lmin2.npz", bpwf=bp15)
np.savez_compressed("data/coadd_bpwf_100mJy_191127_lmin2.npz", bpwf=bp100)
