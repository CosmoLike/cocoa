# Baseline likelihoods release

This file contains 10 likelihood files that forms the baseline Planck likelihood data release. One is used for the low-&#8467; TT (`commander`), three provide the low-&#8467; polarization part (`simall`) to account for either the EE contribution (baseline) the BB one or both the EE and BB one (options). On the high-&#8467; side, two files are provided for the baseline TT and TTTEEE `plik` likelihood, and two others for the foreground and nuisance marginalized `plik_lite`. Finally, the lensing likelihood is provided in two versions, one model dependent, and the other marginalized over the CMB reconstructed in the `plik_lite` likelihood.

This file extracts to a directory hierarchy, containing the different data sets needed to compute different likelihoods.

## Low multipoles (&#8467; ≤ 29)

### TT only - `commander` 
This file allows for the computation of the CMB TT likelihood in the range &#8467;=2-29.
Using the optional tool in the code package, it can be modified to cover any multipole range 
up to &#8467;&lt;200. 

#### Production process
The likelihood is based on the results of the Commander approach, which implements a 
Bayesian component separation method in pixel space, sampling the posterior distribution 
of the parameters of a model that describes both the CMB and the foreground emissions
in a combination of the Planck maps. The samples
of this exploration are used to infer the foreground marginalized low-&#8467; likelihood 
for any <i>TT</i> CMB spectrum. 

#### Inputs:
* Planck 30- and 44-GHz frequency maps;
* Planck 70- to 857-GHz detector and detector-set maps;
* Commander  2018 confidence mask.

#### File name and usage
The commander likelihood is distributed in 
`plc_3.0/low_l/commander/commander_dx12_v3_2_29.clik`

When used with the library, this expects a vector of parameters consisting of the 
<i>TT</i> CMB power spectrum from &#8467;=0 to 29 (inclusive) and an extra nuisance parameter 
consisting of the overall Planck calibration.  Note that the vector really starts at &#8467;=0, although`
the first two entries are null.

### EE only - `simall`
This file allows for the computation of the EE likelihood 
in the range &#8467;=2-29. It should be used with the `commander` TT-only likelihood to form the baseline low-&#8467; likelihood.

#### Production process
The likelihood is estimated by comparing cross quasi maximum likelihood algorithm (QML) on the 100 and 143 GHz maps to high fidelity end-to-end simulation of the HFI instrument.
The galactic contamination (synchrotron and dust) in the input 100 and 143 GHz maps is mitigated by regressing the (LFI) 30 GHz and (HFI) 353 GHz maps. In order to avoid higly contaminated area of the sky, the regions near the galactic center and plane are masked and we retain only about 50% of the sky. The mask is built by applying a threshold on the 353 GHz map and combined with the `commander` confidence mask.
The end-to-end (E2E) simulations contain a realistic galactic contamination and our best knowledge of the HFI detector and analysis chain.
The E2E simulations only explore a single CMB realization, but at large scales, a CMB swapping procedure can be implemented to explore different cosmological parameter sets (within LambdaCDM) and different CMB realizations. The E2E realization are processed in the same way as the data (foreground mitigation, masking, QML) and using the distance between the data and the simulations, an approximation of the log likelihood is produced for each multipole. The final likelihood thus ignore multipole to multipole  correlations as well as TT to EE ones.
The likelihood files contains metadata allowing to compute the approximation for each multipole.

#### Inputs:
* Planck 30-, 100-, 143- and 353-GHz frequency maps;
* `Commander`  2018 confidence mask  and a specially designed galactic masks based on the 353 GHz map
* 300 single CMB and foreground realization end-to-end simulations. 300 000 large scale CMB only realization (300 map realizations for 1000 different random set of cosmological parameters.)

#### File name and usage
The low ell EE likelihood is distributed in 
`plc_3.0/low_l/simall/simall_100x143_offlike5_EE_Aplanck_B.clik`.

When used with the library, this expects a vector of parameters consisting of the 
<i>EE</i> CMB power spectrum from &#8467;=0 to 29 (inclusive) and by an extra nuisance parameter 
consisting of the overall Planck calibration.  Note that the entries really start at &#8467;=0, although the
&#8467;=0 and &#8467;=1 values will be null.

### BB only - simall
This file allows for the computation of the BB likelihood 
in the range &#8467;=2-29. It is not part of the baseline likelihood. The BB spectrum at large scale is compatible with zero.

#### Production process
The likelihood is estimated by comparing cross quasi maximum likelihood algorithm (QML) on the 100 and 143 GHz maps to high fidelity end-to-end simulation of the HFI instrument.
The galactic contamination (synchrotron and dust) in the input 100 and 143 GHz maps is mitigated by regressing the (LFI) 30 GHz and (HFI) 353 GHz maps. In order to avoid higly contaminated area of the sky, the regions near the galactic center and plane are masked and we retain only about 50% of the sky. The mask is built by applying a threshold on the 353 GHz map and combined with the commander confidence mask.
The end-to-end (E2E) simulations contain a realistic galactic contamination and our best knowledge of the HFI detector and analysis chain. 
The E2E simulations only explore a single CMB realization, but at large scales, a CMB swapping procedure can be implemented to explore different cosmological parameter sets (within LambdaCDM) and different CMB realizations. The E2E realization are processed in the same way as the data (foreground mitigation, masking, QML) and using the distance between the data and the simulations, an approximation of the log likelihood is produced for each multipole. The final likelihood thus ignore l-to-l correlations.
The likelihood files contains metadata allowing to compute the approximation for each multipole.

#### Inputs

* Planck 30-, 100-, 143- and 353-GHz frequency maps;
* Commander  2018 confidence mask and a specially designed galactic masks based on the 353 GHz map
* 300 single CMB and foreground realization end-to-end simulations. 300 000 large scale CMB only realization (300 map realizations for 1000 different random set of cosmological parameters.)

#### File name and usage
The low ell BB likelihood is distributed in 
`plc_3.0/low_l/simall/simall_100x143_offlike5_BB_Aplanck_B.clik`.

When used with the library, this expects a vector of parameters consisting of the 
<i>BB</i> CMB power spectrum from &#8467;=0 to 29 (inclusive) and by an extra nuisance parameter 
consisting of the overall Planck calibration.  Note that the entries really start at &#8467;=0, although the
&#8467;=0 and &#8467;=1 values will be null.

### EEBB - `simall`

This file allows for the computation of the EEBB likelihood 
in the range &#8467;=2-29. It is not part of the baseline likelihood. Since correlation between spectra and multipole are ignored in simall, using this file is equivalent (up to a possible normalization) to summing the log likelihood provided by the EE and BB only simall files.

#### Production process
The likelihood is estimated by comparing cross quasi maximum likelihood algorithm (QML) on the 100 and 143 GHz maps to high fidelity end-to-end simulation of the HFI instrument.
The galactic contamination (synchrotron and dust) in the input 100 and 143 GHz maps is mitigated by regressing the (LFI) 30 GHz and (HFI) 353 GHz maps. In order to avoid higly contaminated area of the sky, the regions near the galactic center and plane are masked and we retain only about 50% of the sky. The mask is built by applying a threshold on the 353 GHz map and combined with the commander confidence mask.
The end-to-end (E2E) simulations contain a realistic galactic contamination and our best knowledge of the HFI detector and analysis chain. 
The E2E simulations only explore a single CMB realization, but at large scales, a CMB swapping procedure can be implemented to explore different cosmological parameter sets (within LambdaCDM) and different CMB realizations. The E2E realization are processed in the same way as the data (foreground mitigation, masking, QML) and using the distance between the data and the simulations, an approximation of the log likelihood is produced for each multipole. The final likelihood thus ignore l-to-l correlations.
The likelihood files contains metadata allowing to compute the approximation for each multipole.

#### Inputs

* Planck 30-, 100-, 143- and 353-GHz frequency maps;
* Commander  2018 confidence mask ({{PlanckPapers|planck2016-l04}}) and a specially designed galactic masks based on the 353 GHz map
* 300 single CMB and foreground realization end-to-end simulations. 300 000 large scale CMB only realization (300 map realizations for 1000 different random set of cosmological parameters.)

#### File name and usage
The low ell BB likelihood is distributed in 
`plc_3.0/low_l/simall/simall_100x143_offlike5_EEBB_Aplanck_B.clik`.

When used with the library, this expects a vector of parameters consisting of the 
<i>EE</i> CMB power spectrum from &#8467;=0 to 29 (inclusive), followed by <i>BB</i> CMB power spectrum from &#8467;=0 to 29 (inclusive) and by an extra nuisance parameter 
consisting of the overall Planck calibration.  Note that the entries really start at &#8467;=0, although the
&#8467;=0 and &#8467;=1 values will be null.

## High-&#8467; likelihoods 

###TT only - Plik
This file allows for the computation of the CMB TT likelihood in the range &#8467;=30-2508.
Using the optional tools in the code package it can be modified to cover any multipole range within 29&lt;&#8467;&lt;2509, and to remove the contribution of any range of multipole from any of the cross spectra considered in the approximation.

#### Production process
The file contains the 100-GHz, 143-GHz, and 217-GHz binned half-mission <i>TT</i> cross-spectra. 
Only the 100&times;100, 143&times;143, 143&times;217, and 217&times;217 spectra are actually used. Masks and multipole 
ranges for each spectrum are different and described in {{PlanckPapers|planck2016-l05}}.
Masks are based on the CMB-cleaned 353-GHz map for the dust component, on the Planck catalogues for the point source part, and on the CO maps.
The file also contains
templates for the residual foreground contamination of each spectrum. The templates are needed
to allow for computation of the joint CMB and nuisance likelihood. The covariance matrix is 
computed using an analytical approximation, and corrected for the effect of point sources through Monte Carlo estimates.
The covariance matrix is computed for a fiducial 
cosmology and nuisance model that has been obtained using a first, less optimal estimate of the parameters.
The beam matrix computed for the specific masks and data cuts are applied to the 2015 TT LambdaCDM best fit spectra to predict leakage templates. 
Subpixel effect are predicted for the specific masks and data cuts.

#### Inputs

* Planck 100-, 143-, and 217-GHz half-mission <i>T</i> maps;
* CMB-cleaned (using the SMICA map) 353-GHz map, CO emission maps, and Planck catalogues for the masks;
* Planck 545-GHz maps for the dust residual contamination template;
* CIB, tSZ, kSZ, and CIB&times;SZ templates;
* Beam matrix (effective beams and leakages) and subpixel effect templates computed for the specific masks and sky fractions retained for the likelihood.

#### File name and usage
The high-&#8467;, Plik TT likelihood is distributed in 
`plc_3.0/hi_l/plik/plik_rd12_HM_v22_TT.clik`.

When used with the library, this expects a vector of parameters consisting of the 
<i>TT</i> CMB power spectrum from &#8467;=0 to 2508 (inclusive), followed by a vector of 20 nuisance parameters.
Those are, in order:

* `A_cib_217`, the CIB contamination at &#8467;=3000 in the 217-GHz Planck map;
* `cib_index`, the effective slope of the CIB spectrum, a parameter that should be set to -1.3;
* `xi_sz_cib`, the SZ&times;CIB cross-correlation;
* `A_sz`, the tSZ contamination at 143GHz;
* `ps_A_100_100`, the point source contribution in 100&times;100;
* `ps_A_143_143`, the point source contribution in 143&times;143;
* `ps_A_143_217`, the point source contribution in 143&times;217;
* `ps_A_217_217`, the point source contribution in 217&times;217;
* `ksz_norm`, the kSZ contamination;
* `gal545_A_100`, the dust residual contamination at &#8467;=200 in 100&times;100;
* `gal545_A_143`, the dust residual contamination at &#8467;=200 in 143&times;143;
* `gal545_A_143_217`, the dust residual contamination at &#8467;=200 in 143&times;217;
* `gal545_A_217`, the dust residual contamination at &#8467;=200 in 217&times;217;
* `A_sbpx_100_100_TT`, a rescaling amplitude for the subpixel effects at &#8467;=200 in 100&times;100 (1 is default)
* `A_sbpx_143_143_TT`, a rescaling amplitude for the subpixel effects at &#8467;=200 in 143&times;143 (1 is default)
* `A_sbpx_143_217_TT`, a rescaling amplitude for the subpixel effects at &#8467;=200 in 143&times;217 (1 is default)
* `A_sbpx_217_217_TT`, a rescaling amplitude for the subpixel effects at &#8467;=200 in 217&times;217 (1 is default)
* `calib_100T`, the relative calibration between the 100 and 143 spectra;
* `calib_217T`, the relative calibration between the 217 and 143 spectra;
* `A_planck`, the Planck absolute calibration.

Recommended priors can be found in the file `plc_3.0/hi_l/plik/plik_recommended_priors.txt`
For `cosmomc` users, a set of initialization files is available in `plc_3.0/cosmomc`

### TT+TE+EE - Plik
This file allows for the computation of the CMB joint TT, TE, and EE likelihood in the range &#8467;=30-2508 for TT and &#8467;=30-1996 for TE and EE.
Using the optional tools in the code package it can be modified to cover any multipole range within 29&lt;&#8467;&lt;2509, and to remove the contribution of any range of multipole from any of the cross spectra (in TT, TE EE and any cross-frequencies) considered in the approximation.

#### Production process
The file contains the 100-GHz, 143-GHz, and 217-GHz binned half-mission <i>T</i> and <i>E</i> cross-spectra. 
In temperature, only the 100&times;100, 143&times;143, 143&times;217, and 217&times;217 spectra are actually used, while
in <i>TE</i> and <i>EE</i> all of them are used. Masks and multipole 
ranges for each spectrum are different and described in {{PlanckPapers|planck2016-l05}}.
Masks are based on CMB-cleaned 353-GHz map for the dust component, on the Planck catalogues for the point source part, and on the CO maps.
The file also contains
templates for the residual foreground contamination of each spectrum. The templates are needed
to allow the computation of the joint CMB and nuisance likelihood. The covariance matrix is 
computed using an analytical approximation, and corrected for the effect of point sources through Monte Carlo estimates.
The covariance matrix is computed for a fiducial 
cosmology and nuisance model that has been obtained using a first, less optimal estimate of the parameters.
The beam matrix computed for the specific masks and data cuts are applied to the 2015 TT LambdaCDM best fit spectra to predict leakage templates. 
Subpixel effect are predicted for the specific masks and data cuts.
The 300 End-to-End HFI simulations from the FFP10 suite are used to estimate correlated noise residuals in the EE autospectra at 100GHz, 143GHz and 217GHz.

### Inputs 

* Planck 100-, 143-, and 217-GHz half-mission T+P maps;
* CMB-cleaned (using the SMICA map) 353-GHz map, CO emission maps, and Planck catalogues for the masks;
* Planck 545-GHz and 353-GHz maps for the dust residual contamination template;
* CIB, tSZ, kSZ, and CIB&times;SZ templates;
* Beam matrix (effective beams and leakages) and subpixel effect templates computed for the specific masks and sky fractions retained for the likelihood.
* 300 End-to-End HFI simulations and resulting smoothed correlated noise templates.

### File name and usage

The high-&#8467;, plik TTTEEE likelihood is distributed in 
`plc_3.0/hi_l/plik/plik_rd12_HM_v22b_TTTEEE.clik`.

This file should not be used with any other TT-only high-&#8467; file.

When used with the library, this expects a vector of parameters consisting of the 
<i>TT</i> CMB power spectrum from &#8467;=0 to 2508 (inclusive), followed the <i>EE</i> and <i>TE</i> spectra (same range) and by a vector of 47 nuisance parameters.
Those are, in g:

* `A_cib_217`, the CIB contamination at &#8467;=3000 in the 217-GHz Planck map;
* `cib_index`, the effective slope of the CIB spectrum, which should be set to -1.3;
* `xi_sz_cib`, the SZ&times;CIB cross-correlation;
* `A_sz`, the tSZ contamination at 143GHz;
* `ps_A_100_100`, the point source contribution in 100&times;100;
* `ps_A_143_143`, the point source contribution in 143&times;143;
* `ps_A_143_217`, the point source contribution in 143&times;217;
* `ps_A_217_217`, the point source contribution in 217&times;217;
* `ksz_norm`, the kSZ contamination;
* `gal545_A_100`, the dust residual contamination at &#8467;=200 in 100&times;100TT;
* `gal545_A_143`, the dust residual contamination at &#8467;=200 in 143&times;143TT;
* `gal545_A_143_217`, the dust residual contamination at &#8467;=200 in 143&times;217TT;
* `gal545_A_217`, the dust residual contamination at &#8467;=200 in 217&times;217TT;
* `galf_EE_A_100`, the dust residual contamination at &#8467;=500 in 100&times;100EE;
* `galf_EE_A_100_143`, the dust residual contamination at &#8467;=500 in 100&times;143EE;
* `galf_EE_A_100_217`, the dust residual contamination at &#8467;=500 in 100&times;217EE;
* `galf_EE_A_143`, the dust residual contamination at &#8467;=500 in 143&times;143EE;
* `galf_EE_A_143_217`, the dust residual contamination at &#8467;=500 in 143&times;217EE;
* `galf_EE_A_217`, the dust residual contamination at &#8467;=500 in 217&times;217EE;
* `galf_EE_index`, the dust EE template slope, which should be set to -2.4;
* `galf_TE_A_100`, the dust residual contamination at &#8467;=500 in 100&times;100TE;
* `galf_TE_A_100_143`, the dust residual contamination at &#8467;=500 in 100&times;143TE;
* `galf_TE_A_100_217`, the dust residual contamination at &#8467;=500 in 100&times;217TE;
* `galf_TE_A_143`, the dust residual contamination at &#8467;=500 in 143x143TE;
* `galf_TE_A_143_217`, the dust residual contamination at &#8467;=500 in 143&times;217TE;
* `galf_TE_A_217`, the dust residual contamination at &#8467;=500 in 217&times;217TE;
* `galf_TE_index`, the dust EE template slope, which should be set to -2.4;
* `A_cnoise_e2e_100_100_EE` Normalization for the end2end empirical correlated noise template at 100GHz (should be set to 1. Set to 0 to null effect of the template)
* `A_cnoise_e2e_143_143_EE` Normalization for the end2end empirical correlated noise template at 143GHz (should be set to 1. Set to 0 to null effect of the template)
* `A_cnoise_e2e_217_217_EE` Normalization for the end2end empirical correlated noise template at 217GHz (should be set to 1. Set to 0 to null effect of the template)
* `A_sbpx_100_100_TT` Normalization for the subpixel effect correction at 100&times;100 TT (should be set to 1)
* `A_sbpx_143_143_TT` Normalization for the subpixel effect correction at 143&times;143 TT (should be set to 1)
* `A_sbpx_143_217_TT` Normalization for the subpixel effect correction at 143&times;217 TT (should be set to 1)
* `A_sbpx_217_217_TT` Normalization for the subpixel effect correction at 217&times;217 TT (should be set to 1)
* `A_sbpx_100_100_EE` Normalization for the subpixel effect correction at 100&times;100 EE (should be set to 1)
* `A_sbpx_100_143_EE` Normalization for the subpixel effect correction at 100&times;143 EE (should be set to 1)
* `A_sbpx_100_217_EE` Normalization for the subpixel effect correction at 100&times;217 EE (should be set to 1)
* `A_sbpx_143_143_EE` Normalization for the subpixel effect correction at 143&times;143 EE (should be set to 1)
* `A_sbpx_143_217_EE` Normalization for the subpixel effect correction at 143&times;217 EE (should be set to 1)
* `A_sbpx_217_217_EE` Normalization for the subpixel effect correction at 217&times;217 EE (should be set to 1)
* `calib_100T`, the relative calibration between 100 and 143 TT spectra;
* `calib_217T`, the relative calibration between 217 and 143 TT spectra;
* `calib_100P`, the calibration of the 100 EE;
* `calib_143P`, the calibration of the 143 EE;
* `calib_217P`, the calibration of the 217 EE;
* `A_pol`, the calibration of the polarization relative to the temperature, which should be set to 1;
* `A_planck`, the Planck absolute calibration.


Recommended priors can be found in the file `plc_3.0/hi_l/plik/plik_recommended_priors.txt`
For `cosmomc` users, a set of initialization files is available in `plc_3.0/cosmomc`

### TT only - Plik lite
This file allows for the computation of the nuisance-marginalized CMB TT likelihood in the range &#8467;=30-2508. It should not be used with the regular high-&#8467; TT or TTTEEE files.


#### Production process
The Plik likelihood files described above have been explored using a Bayesian algorithm described in {{PlanckPapers|planck2015-l05}}. The joint posterior of the CMB <i>TT</i> spectrum, marginalized over the nuisance parameters, has been extracted from this analysis to build a high-&#8467; likelihood approximation for the 
CMB spectrum only. The nuisance parameters have been marginalized in the priors described above.

#### Inputs:

* Plik `plik_rd12_HM_v22b_TTTEEE` likelihood;
* dust residual, CIB, tSZ, kSZ, and CIB&times;SZ, leakage and subpixel templates.

#### File name and usage

The high-&#8467;, Plik TT, nuisance-marginalized likelihood is distributed in 
`plc_3.0/hi_l/plik_lite/plik_lite_v22_TT.clik`.

When used with the library, this expects a vector of parameters consisting of the 
<i>TT</i> CMB power spectrum from &#8467;=0 to 2508 (inclusive), followed by the Planck absolute calibration nuisance parameter.

### TT EE TE - Plik lite
This file allows for the computation of the nuisance marginalized CMB joint TT, TE, EE likelihood in the range &#8467;=30-2508 for temperature and &#8467;=30-1996 for TE and EE. It should not be used with the regular high-&#8467; TT or TTTEEE files.


#### Production process

The Plik likelihood file described above have been explored using a Bayesian algorithm described in the Planck 2018 likelihood paper. The joint posterior of the CMB <i>TT</i>, <i>TE</i>, and <i>EE</i> spectra, marginalized over the nuisance parameters has been extracted from this analysis to build a high-&#8467 likelihood approximation for the 
CMB spectrum only. The nuisance parameters have been marginalized in the priors described above.


#### Inputs

* Plik `plik_rd12_HM_v22b_TTTEEE` likelihood;
* dust residual, CIB, tSZ, kSZ, and CIB&times;SZ leakage, subpixel and correlated noise templates.

#### File name and usage

The high-&#8467;, Plik TTTEEE, nuisance-marginalized likelihood is distributed in 
`plc_3.0/hi_l/plik_lite/plik_lite_v22_TTTEEE.clik`.

When used with the library, this expects a vector of parameters consisting of the 
<i>TT</i> CMB power spectrum from &#8467;=0 to 2508 (inclusive), followed by the <i>EE</i>, and <i>TE</i> spectra (same range) and by the Planck absolute calibration nuisance parameter.

## Lensing likelihoods 

### CMB dependent likelihood
This file allows for the computation of the baseline lensing likelihood, using the SMICA temperature and polarization 
map-based lensing reconstruction. Only the lensing multipole range &#8467;=8-400
is used. The covariance matrix for the likelihood is based on Monte Carlos. The likelihood includes 
the computation of a model-dependent correction to the normalization and
the "N1 bias." To speed up those computations, both are performed using a linear approximation.

#### Production process
The SMICA <i>T</i> maps are filtered and correlated to reconstruct an optimal 
lensing reconstruction map. Biases are corrected using Monte Carlo simulations
for the dominant "mean field" and "N0" contributions, while the "N1" bias is computed using the CMB-based best-fit model. The power spectrum of this map is measured out of a large mask, and binned.
The covariance matrix of this binned spectrum is evaluated using numerical simulations.
The model-dependent corrections to the normalization and N1 biases are linearized and precomputed.

#### Inputs

* SMICA temperature and polarization CMB map;
* Galactic mask;
* best-fit Planck CMB model.

#### File name and usage

The lensing likelihood is distributed in 
`plc_3.0/hi_l/lensing/smicadx12_Dec5_ftl_mv2_ndclpp_p_teb_consext8.clik_lensing`.

When used with the library, this expects a vector of parameters consisting of the 
&phi;&phi; lensing spectrum for &#8467;=0 to 2500 (inclusive), followed by the 
followed by the 
<i>TT</i> CMB power spectrum for &#8467;=0 to 2500 (inclusive), the <i>EE</i> and <i>TE</i> power spectra (same range) and by the Planck absolute calibration nuisance parameter.

The <i>TT</i>, <i>EE</i>, and <i>TE</i> spectra are needed to compute the normalization and N1 corrections.

#### CMB marginalized lensing likelihood
This file allows for the computation of the baseline lensing likelihood, marginalized over the CMB power spectrum. The covariance is enlarged compared to the previous case to account for the change in N0 and N1 dues to possible variations of the CMB power spectrum. As in the previous case, only the lensing multipole in the range range &#8467;=8-400
are used. 

#### Production process
The covariance of the CMB dependent likelihood is enlarged, and the theory spectrum is shifted to account for th marginalization over the CMB spectrum, using a gaussian approximation. The `plik_lite` likelihood is used to provide the CMB spectra and covariances.


#### Inputs:  

* SMICA <i>T</i> and <i>E</i> CMB maps;
* Galactic mask;
* best-fit Planck CMB model.
* `plik_lite` likelihood

#### File name and usage

The CMB marginalized lensing likelihood is distributed in 
`plc_3.0/hi_l/lensing/smicadx12_Dec5_ftl_mv2_ndclpp_p_teb_consext8_CMBmarged.clik_lensing`.

When used with the library, this expects a vector of parameters consisting of the 
&phi;&phi; lensing spectrum for &#8467;=0 to 2500 (inclusive), followed by the Planck absolute calibration nuisance parameter.
