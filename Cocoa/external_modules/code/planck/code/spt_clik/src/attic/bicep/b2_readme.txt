### Readme.txt###
1/ Download the file b2_hl_likelihood.tgz in in a unix folder and untar. 
This will generate the followings files:
       -  bicep_r_wrapper.py, and bicep_util.py and b2_readme.txt
       -  B2_3yr_bpcm_no-sysuncer_20140310.txt, B2_3yr_likelihood_bandpowers_20140310.txt, B2_3yr_camb_planck_lensed_uK_20140226.txt, B2_3yr_camb_planck_withB_uK_20140226.txt
       -  windows/B2*
2/ The above data release files are identical to those provided on the BICEP2 data release page:  http://bicepkeck.org/bicep2_2014_release/
3/ If needed, make the top module bicep_r_wrapper.py executable (chmod +x bicep_r_wrapper.py)
4/ This package requires a working installation of the Numpy and Scipy packages.
5/ From the working directory, execute the top module: ./bicep_r_wrapper.py
6/ Use the help for instructions: ./bicep_r_wrapper.py -h
7/ As an example to generate the HL likelihood for BICEP2 versus r, run: ./bicep_r_wrapper.py -e bicep2
####