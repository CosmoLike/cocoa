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

echo -e '\033[0;32m'"\t\t DECOMPRESSING SN DATA"'\033[0m'
tar xf sn_data.xz

echo -e '\033[0;32m'"\t\t DECOMPRESSING BAO DATA"'\033[0m'
tar xf bao_data.xz

echo -e '\033[0;32m'"\t\t DECOMPRESSING H0LICOW DATA"'\033[0m'
tar xf h0licow_distance_chains.xz

echo -e '\033[0;32m'"\t\t DECOMPRESSING ACT-DR6 DATA"'\033[0m'
mkdir -p act
cd $ROOTDIR/external_modules/data/act
mkdir -p lensing
cd $ROOTDIR/external_modules/data/act/lensing
wget https://lambda.gsfc.nasa.gov/data/suborbital/ACT/ACT_dr6/likelihood/data/ACT_dr6_likelihood_v1.2.tgz > ${OUTPUT_UNXV_ALL_1} 2> ${OUTPUT_UNXV_ALL_2}
tar -zxvf ACT_dr6_likelihood_v1.2.tgz > ${OUTPUT_UNXV_ALL_1} 2> ${OUTPUT_UNXV_ALL_2}
cd $ROOTDIR/external_modules/data

if [ -z "${THREAD_UNXZ}" ]; then
	echo -e '\033[0;32m'"\t\t DECOMPRESSING SIMONS OBSERVATORY"'\033[0m'
	tar xf simons_observatory.xz
	echo -e '\033[0;32m'"\t\t DECOMPRESSING BICEP 2015 DATA"'\033[0m'
	tar xf bicep_keck_2015.xz
	# ---------------------------------------------
	echo -e '\033[0;32m'"\t\t DECOMPRESSING SPT-3G Y1 DATA"'\033[0m'
	tar xf spt_3g.xz
	# ---------------------------------------------
	cd ./planck
	# ---------------------------------------------
	echo -e '\033[0;32m'"\t\t DECOMPRESSING SUPPLEMENTAL DATA AND COVARIANCES"'\033[0m'
	tar xf planck_supp_data_and_covmats.xz
	# ---------------------------------------------
	echo -e '\033[0;32m'"\t\t DECOMPRESSING PLANCK-2018 (PLC-3.0) DATA"'\033[0m'
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
	echo -e '\033[0;32m'"\t\t DECOMPRESSING THE REMAINING PACKAGES (SO, Planck, SPT...) IN PARALLEL"'\033[0m'
	tar xf simons_observatory.xz &
	proc8=$!
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
	wait "$proc1" "$proc2" "$proc4" "$proc5" "$proc6" "$proc7" "$proc8" "$proc10"
fi
