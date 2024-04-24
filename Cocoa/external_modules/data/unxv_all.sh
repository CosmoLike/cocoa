#!/bin/bash

if [ -z "${DEBUG_UNXV_CLEAN_ALL}" ]; then
  export OUTPUT_UNXV_ALL_1="/dev/null"
  export OUTPUT_UNXV_ALL_2="/dev/null"
else
  export OUTPUT_UNXV_ALL_1="/dev/tty"
  export OUTPUT_UNXV_ALL_2="/dev/tty"
fi

sh $ROOTDIR/external_modules/data/clean_all.sh 
cd $ROOTDIR/external_modules/data

echo 'DECOMPRESSING SN DATA' > ${OUTPUT_UNXV_ALL_1} 2> ${OUTPUT_UNXV_ALL_2}
tar xf sn_data.xz
echo 'DECOMPRESSING BAO DATA' > ${OUTPUT_UNXV_ALL_1} 2> ${OUTPUT_UNXV_ALL_2}
tar xf bao_data.xz
echo 'DECOMPRESSING H0LICOW DATA' > ${OUTPUT_UNXV_ALL_1} 2> ${OUTPUT_UNXV_ALL_2}
tar xf h0licow_distance_chains.xz

if [ -z "${THREAD_UNXZ}" ]; then
	echo 'DECOMPRESSING SIMONS OBSERVATORY' > ${OUTPUT_UNXV_ALL_1} 2> ${OUTPUT_UNXV_ALL_2}
	tar xf simons_observatory.xz
	echo 'DECOMPRESSING SPT-3G HIELL DATA' > ${OUTPUT_UNXV_ALL_1} 2> ${OUTPUT_UNXV_ALL_2}
	tar xf spt_hiell_2020.xz
	echo 'DECOMPRESSING BICEP 2015 DATA' > ${OUTPUT_UNXV_ALL_1} 2> ${OUTPUT_UNXV_ALL_2}
	tar xf bicep_keck_2015.xz
	# ---------------------------------------------
	echo 'DECOMPRESSING SPT-3G Y1 DATA' > ${OUTPUT_UNXV_ALL_1} 2> ${OUTPUT_UNXV_ALL_2}
	tar xf spt_3g.xz
	# ---------------------------------------------
	echo 'DECOMPRESSING ACT-DR6 DATA' > ${OUTPUT_UNXV_ALL_1} 2> ${OUTPUT_UNXV_ALL_2}
	tar xf act.xz
	# ---------------------------------------------
	cd ./planck
	# ---------------------------------------------
	echo 'DECOMPRESSING SUPPLEMENTAL DATA AND COVARIANCES' > ${OUTPUT_UNXV_ALL_1} 2> ${OUTPUT_UNXV_ALL_2}
	tar xf planck_supp_data_and_covmats.xz
	# ---------------------------------------------
	echo 'DECOMPRESSING PLANCK-2015 (PLC-2.0) DATA' > ${OUTPUT_UNXV_ALL_1} 2> ${OUTPUT_UNXV_ALL_2}
	tar xf plc_20.xz
	# ---------------------------------------------
	echo 'DECOMPRESSING PLANCK-2018 (PLC-3.0) DATA' > ${OUTPUT_UNXV_ALL_1} 2> ${OUTPUT_UNXV_ALL_2}
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
	tar xf act.xz &
	proc11=$!
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
	echo 'DECOMPRESSING SUPPLEMENTAL DATA AND COVARIANCES' > ${OUTPUT_UNXV_ALL_1} 2> ${OUTPUT_UNXV_ALL_2}
	echo 'DECOMPRESSING PLANCK-2015 (PLC-2.0) DATA' > ${OUTPUT_UNXV_ALL_1} 2> ${OUTPUT_UNXV_ALL_2}
	echo 'DECOMPRESSING PLANCK-2018 (PLC-3.0) DATA' > ${OUTPUT_UNXV_ALL_1} 2> ${OUTPUT_UNXV_ALL_2}
	echo 'DECOMPRESSING SPT-3G Y1 DATA' 			> ${OUTPUT_UNXV_ALL_1} 2> ${OUTPUT_UNXV_ALL_2}
	echo 'DECOMPRESSING ACT-DR6 DATA' 				> ${OUTPUT_UNXV_ALL_1} 2> ${OUTPUT_UNXV_ALL_2}
	echo 'DECOMPRESSING SIMONS OBSERVATORY'         > ${OUTPUT_UNXV_ALL_1} 2> ${OUTPUT_UNXV_ALL_2}
	echo 'DECOMPRESSING BICEP 2015 DATA'            > ${OUTPUT_UNXV_ALL_1} 2> ${OUTPUT_UNXV_ALL_2}
	echo 'DECOMPRESSION IS HAPPENING IN PARALLEL - WAITING ALL OF THEM TO FINISH' > ${OUTPUT_UNXV_ALL_1} 2> ${OUTPUT_UNXV_ALL_2}
	wait "$proc1" "$proc2" "$proc3" "$proc4" "$proc5" "$proc6" "$proc7" "$proc8" "$proc9" "$proc10" "$proc11"
fi
