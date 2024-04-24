if [ -z "${ROOTDIR}" ]; then
    echo 'ERROR ROOTDIR not defined'
    return
fi

if [ -z "${DEBUG_UNXV_CLEAN_ALL}" ]; then
  export OUTPUT_CLEAN_ALL_1="/dev/null"
  export OUTPUT_CLEAN_ALL_2="/dev/null"
else
  export OUTPUT_CLEAN_ALL_1="/dev/tty"
  export OUTPUT_CLEAN_ALL_2="/dev/tty"
fi

echo 'CLEAN DATA FOLDERS - THAT MIGHT TAKE A WHILE' > ${OUTPUT_CLEAN_ALL_1} 2> ${OUTPUT_CLEAN_ALL_2}

cd $ROOTDIR/external_modules/data

rm -rf bao_data/
rm -rf bicep_keck_2015/
rm -rf sn_data/
rm -rf simons_observatory/
rm -rf spt_hiell_2020/
rm -rf spt_3g/
rm -rf act/
rm -rf planck/CamSpec2018
rm -rf planck/plc_2.0
rm -rf planck/plc_3.0/low_l
rm -rf planck/plc_3.0/lensing
rm -rf planck/spt3g_Y1_EETE.clik
rm -rf planck/planck_supp_data_and_covmats
rm -rf planck/plc_3.0/hi_l/plik
rm -rf planck/plc_3.0/hi_l/plik_lite
rm -rf planck/plc_3.0/hi_l/camspec/camspec_10.7HM_1400_TT_small.clik
rm -rf planck/plc_3.0/hi_l/camspec/camspec_10.7HM_1400_TTTEEE.clik
rm -rf h0licow_distance_chains

cd $ROOTDIR/

echo 'CLEAN DATA FOLDERS - DONE' > ${OUTPUT_CLEAN_ALL_1} 2> ${OUTPUT_CLEAN_ALL_2}