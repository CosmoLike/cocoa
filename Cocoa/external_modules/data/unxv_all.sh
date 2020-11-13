#!/bin/bash

rm -rf ./sn_data/
tar xf sn_data.xz

rm -rf ./des_data/
tar xf des_data.xz

rm -rf ./bicep_keck_2015/
tar xf bicep_keck_2015.xz

rm -rf ./bao_data/
tar xf bao_data.xz

if [ -z "${THREAD_UNXZ}" ]; then
	cd ./planck
	echo 'DECOMPRESSING CAMSPEC (2018)'
	tar xf CamSpec2018.xz
	echo 'DECOMPRESSING SUPPLEMENTAL DATA AND COVARIANCES'
	tar xf planck_supp_data_and_covmats.xz
	echo 'DECOMPRESSING PLANCK-2015 (PLC-2.0)'
	tar xf plc_20.xz
	echo 'DECOMPRESSING PLANCK-2018 (PLC-3.0)'
	cd ./plc_3.0
	tar xf lensing.xz
	tar xf low_l.xz
	cd ./hi_l
	tar xf plik.xz
	tar xf plik_lite.xz
	cd ./camspec
	tar xf camspec_107HM_1400_TT_smallclik.xz
	tar xf camspec_10.7HM_1400_TTTEEEclik.xz
else
	cd ./planck
	tar xf CamSpec2018.xz &
	proc1=$!
	tar xf planck_supp_data_and_covmats.xz &
	proc2=$!
	tar xf plc_20.xz &
	proc3=$!
	cd ./plc_3.0
	tar xf lensing.xz &
	proc4=$!
	tar xf low_l.xz &
	proc5=$!
	cd ./hi_l
	tar xf plik.xz &
	proc6=$!
	tar xf plik_lite.xz	&
	proc7=$!
	cd ./camspec
	tar xf camspec_107HM_1400_TT_smallclik.xz &
	proc8=$!
	tar xf camspec_10.7HM_1400_TTTEEEclik.xz &
	proc9=$!

	echo 'DECOMPRESSING CAMSPEC (2018)'
	echo 'DECOMPRESSING SUPPLEMENTAL DATA AND COVARIANCES'
	echo 'DECOMPRESSING PLANCK-2015 (PLC-2.0)'
	echo 'DECOMPRESSING PLANCK-2018 (PLC-3.0)'
	echo 'DECOMPRESSION IS HAPPENING IN PARALLEL - WAITING ALL OF THEM TO FINISH'
	wait "$proc1" "$proc2" "$proc3" "$proc4" "$proc5" "$proc6" "$proc7" "$proc8" "$proc9"
fi
