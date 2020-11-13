#!/bin/bash

echo 'DECOMPRESSING PLANCK LIKELIHOOD ASSOCIATED FILES - THAT MIGHT TAKE A WHILE'

rm -f cib_1h_2h_100_353_Jsr-1_PS_2014_09.dat
xz --decompress -k cib_1h_2h_100_353_Jsr-1_PS_2014_09.xz
mv cib_1h_2h_100_353_Jsr-1_PS_2014_09 cib_1h_2h_100_353_Jsr-1_PS_2014_09.dat

rm -f cleak_eh_rd12rc3_v1.dat
xz --decompress -k cleak_eh_rd12rc3_v1.xz
mv cleak_eh_rd12rc3_v1 cleak_eh_rd12rc3_v1.dat

rm -f cnoise_e2e_v2.dat
xz --decompress -k cnoise_e2e_v2.xz
mv cnoise_e2e_v2 cnoise_e2e_v2.dat

rm -f cnoise_F100_143_217_353_v17.dat
xz --decompress -k cnoise_F100_143_217_353_v17.xz
mv cnoise_F100_143_217_353_v17 cnoise_F100_143_217_353_v17.dat

rm -f ksz_fromcamspec.dat
xz --decompress -k ksz_fromcamspec.xz
mv ksz_fromcamspec ksz_fromcamspec.dat

rm -f sbpx_tmpl_v4_hm.dat
xz --decompress -k sbpx_tmpl_v4_hm.xz
mv sbpx_tmpl_v4_hm sbpx_tmpl_v4_hm.dat

rm -f sky_template_v15_F100_143_217_353.dat
xz --decompress -k sky_template_v15_F100_143_217_353.xz
mv sky_template_v15_F100_143_217_353 sky_template_v15_F100_143_217_353.dat

rm -f sz_x_cib_template.dat
xz --decompress -k sz_x_cib_template.xz
mv sz_x_cib_template sz_x_cib_template.dat

rm -f tsz_143_eps0.50.dat
xz --decompress -k tsz_143_eps0.50.xz
mv tsz_143_eps0.50 tsz_143_eps0.50.dat