#!/bin/bash

source $ROOTDIR/external_modules/data/clean_all.sh 
cd $ROOTDIR/external_modules/data

echo 'DECOMPRESSING SN DATA'
tar xf sn_data.xz
echo 'DECOMPRESSING BAO DATA'
tar xf bao_data.xz
echo 'DECOMPRESSING H0LICOW DATA'
tar xf h0licow_distance_chains.xz

if [ -z "${THREAD_UNXZ}" ]; then
	echo 'DECOMPRESSING SIMONS OBSERVATORY'
	tar xf simons_observatory.xz
	echo 'DECOMPRESSING SPT-3G HIELL DATA'
	tar xf spt_hiell_2020.xz
	echo 'DECOMPRESSING BICEP 2015 DATA'
	tar xf bicep_keck_2015.xz
	# ---------------------------------------------
	echo 'DECOMPRESSING SPT-3G Y1 DATA'
	tar xf spt_3g.xz
	# ---------------------------------------------
	echo 'DECOMPRESSING ACT-DR6 DATA'
	tar xf act.xz
	# ---------------------------------------------
	cd ./planck
	# ---------------------------------------------
	echo 'DECOMPRESSING SUPPLEMENTAL DATA AND COVARIANCES'
	tar xf planck_supp_data_and_covmats.xz
	# ---------------------------------------------
	echo 'DECOMPRESSING PLANCK-2015 (PLC-2.0) DATA'
	tar xf plc_20.xz
	# ---------------------------------------------
	echo 'DECOMPRESSING PLANCK-2018 (PLC-3.0) DATA'
	cd ./plc_3.0
	tar xf lensing.xz
	tar xf low_l.xz
	cd ./hi_l
	rm -rf ./plik/
	tar xf plik.xz
	rm -rf ./plik_lite/
	tar xf plik_lite.xz
	# ---------------------------------------------
else
	tar xf simons_observatory.xz &
	proc8=$!
	tar xf spt_hiell_2020.xz &
	proc9=$!
	tar xf bicep_keck_2015.xz &
	proc10=$!
	# ---------------------------------------------
	tar xf spt_3g.xz &
	proc1=$!
	# ---------------------------------------------
	cd ./planck
	# ---------------------------------------------
	tar xf planck_supp_data_and_covmats.xz &
	proc2=$!
	# ---------------------------------------------
	tar xf plc_20.xz &
	proc3=$!
	# ---------------------------------------------
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
	# ---------------------------------------------
	# ---------------------------------------------
	# ---------------------------------------------
	echo 'DECOMPRESSING SUPPLEMENTAL DATA AND COVARIANCES'
	echo 'DECOMPRESSING PLANCK-2015 (PLC-2.0) DATA'
	echo 'DECOMPRESSING PLANCK-2018 (PLC-3.0) DATA'
	echo 'DECOMPRESSING SPT-3G Y1 DATA'
	echo 'DECOMPRESSING SPT-3G HIELL DATA'
	echo 'DECOMPRESSING SIMONS OBSERVATORY'
	echo 'DECOMPRESSING BICEP 2015 DATA'
	echo 'DECOMPRESSION IS HAPPENING IN PARALLEL - WAITING ALL OF THEM TO FINISH'
	wait "$proc1" "$proc2" "$proc3" "$proc4" "$proc5" "$proc6" "$proc7" "$proc8" "$proc9" "$proc10" 
fi
