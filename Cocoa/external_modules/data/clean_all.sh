if [ -z "${ROOTDIR}" ]; then
    echo 'ERROR ROOTDIR not defined'
    return
fi

rm -rf bao_data
rm -rf bicep_keck_2015
rm -rf des_data
rm -rf sn_data
rm -rf simons_observatory
rm -rf spt_hiell_2020
rm -rf planck/CamSpec2018
rm -rf planck/plc_2.0
rm -rf planck/plc_3.0/low_l
rm -rf planck/plc_3.0/lensing
rm -rf planck/planck_supp_data_and_covmats
rm -rf planck/plc_3.0/hi_l/plik
rm -rf planck/plc_3.0/hi_l/plik_lite
rm -rf planck/plc_3.0/hi_l/camspec/camspec_10.7HM_1400_TT_small.clik
rm -rf planck/plc_3.0/hi_l/camspec/camspec_10.7HM_1400_TTTEEE.clik
rm -rf h0licow_distance_chains
rm -rf planck/spt3g_Y1_EETE.clik