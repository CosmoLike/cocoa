#!/bin/bash
# ---------------------------------------------------
# ---------------------------------------------------
# ---------------------------------------------------
PLANCK_CLIK='planck_clik.py'

sed --in-place --regexp-extended "s@data_path(\s{1,})=(\s{1,}).*@data_path = 'planck'@g" $PLANCK_CLIK
sed --in-place --regexp-extended "s@return os.path.join.*@return os.path.join\(os.getenv\('ROOTDIR'\),'.local/lib/python/site-packages/clik')@g" $PLANCK_CLIK

