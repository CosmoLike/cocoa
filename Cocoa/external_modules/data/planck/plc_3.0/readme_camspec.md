# `camspec` extented release

This file contains 7 likelihood files, corresponding to the `camspec` likelihoods described in the 2018 likelihood paper and parameter papers. The `camspec` likelihood is a different implementation than the baseline `plik` one and explore a few different choices (noise model, **polarization efficiency corrections**, **polarization masks**, foreground templates). The likelihoods can be used instead of the  `plik` ones to form an alternative Planck CMB likelihood. They allow for the computation of the high-&#8467; likelihood using HFI data.

## TT only 
This file allows for the computation of the CMB TT likelihood in the range &#8467;=30-2500.


### Production process
The file contains the 100-GHz, 143-GHz, and 217-GHz half-mission <i>TT</i> cross-spectra. 
Only the 100&times;100, 143&times;143, 143&times;217, and 217&times;217 spectra are actually used. Masks and multipole 
ranges for each spectrum are different and described in the Planck 2018 likelihood paper.
Masks are based on the CMB-cleaned 353-GHz map for the dust component, on the Planck catalogues for the point source part, and on the CO maps.
The file also contains
templates for the residual foreground contamination of each spectrum. The templates are needed
to allow for computation of the joint CMB and nuisance likelihood. The covariance matrix is 
computed using an analytical approximation.
The covariance matrix is computed for a fiducial 
cosmology and nuisance model that has been obtained using a first, less optimal estimate of the parameters.


### Inputs
* Planck 100-, 143-, and 217-GHz half-mission <i>T</i> maps;
* CMB-cleaned (using the SMICA map) 353-GHz map, CO emission maps, and Planck catalogues for the masks;
* Planck 545-GHz maps for the dust residual contamination template;
* CIB, tSZ, kSZ, and CIB&times;SZ templates;
* Covariance matrix is restricted to the TT case.

### File name and usage
The camspec TT only likelihood is distributed in 
`plc_3.0/hi_l/camxpec/camspec_10.7HM_1400_TT_small.clik`.

When used with the library, this expects a vector of parameters consisting of the 
<i>TT</i> CMB power spectrum from &#8467;=0 to 2500 (inclusive), followed 23 nuisance parameters.

Nuisance parameters are, in order:

* `aps100`        : Point source contamination at 100x100           
* `aps143`        : Point source contamination at 143x143           
* `aps217`        : Point source contamination at 217x217           
* `acib143`       : CIB contamination at 143x143           (set to negative value for default)
* `acib217`       : CIB contamination at 217x217           
* `asz143`        : tSZ contamination at 143x143           
* `psr`           : Point source correlation coefficient between the 143- and 217-GHz         
* `cibr`          : CIB correlation coefficient between the 143- and 217-GHz         (set to negative value for default)
* `ncib143`       : CIB index at 143-GHz (set to negative value for default)
* `ncib`          : CIB index (set to 0 for default)
* `cibrun`        : CIB running (set to 0 for default)            
* `xi`            : CIB-SZ correlation
* `aksz`          : kSZ contamination at 143x143           
* `dust100`       : Galactic contamination normalization template (set to 1 for default)             
* `dust143`       : Galactic contamination normalization template (set to 1 for default)             
* `dust217`       : Galactic contamination normalization template (set to 1 for default)             
* `dust143x217`   : Galactic contamination normalization template (set to 1 for default)                 
* `calPlanck`     : Planck calibration
* `cal0`          : 100GHz TT calibration
* `cal1`          : 143GHz TT calibration
* `cal2`          : 217GHz TT calibration
* `calTE`         : TE spectrum recalibration (unused in TT only likelihood)
* `calEE`         : EE spectrum recalibration (unused in TT only likelihood)

The recommended values and priors are available in the file `plc_3.0/hi_l/camspec/camspec_recommendations.txt`

## EE only 
This file allows for the computation of the CMB EE likelihood in the range &#8467;=30-2000.

### Production process
The file contains the coadded 100-GHz, 143-GHz, and 217-GHz half-mission <i>EE</i> cross-spectra. 
Masks are based on the 353-GHz Q and U maps. 
Each of the different cross spectrum is corrected for dust contamination following an hybrid method using a map based correction up to  &#8467;=300 and a power spectrum template based one for higher multipoles. Details of teh dust correction, as well as multipole ranges retained for each cross spectra is described in the Planck 2018 likelihood paper.
Individual power spectra are corrected for polarization efficiency errors following the so called "spectrum-based" procedure also described in the paper.
Temperature to polarization leakages are corrected using the beam matrices computed for the specific mask and data cuts and using the temperature only 2015 best fit model.
The corrected cross spectra are coadded using the full EE covariance matrix. 

### Inputs
* Planck 100-, 143-, and 217-GHz half-mission <i>E</i> maps;
* Planck 353-GHz maps for the dust residual contamination correction and the mask;
* Covariance matrix is the full TTTEEE one restricted to the EE case at runtime.

### File name and usage
The camspec EE only likelihood is distributed in 
`plc_3.0/hi_l/camspec/camspec_10.7HM_1400_EE.clik`.

When used with the library, this expects a vector of parameters consisting of the 
<i>EE</i> CMB power spectrum from &#8467;=0 to 2500 (inclusive), followed 23 nuisance parameters.

The 2000<&#8467;<2501 entries are ignored.

Nuisance parameters are:

* `aps100`        : Point source contamination at 100x100          (unused in EE but required) 
* `aps143`        : Point source contamination at 143x143          (unused in EE but required)
* `aps217`        : Point source contamination at 217x217          (unused in EE but required)
* `acib143`       : CIB contamination at 143x143           (set to negative value for default) (unused in EE but required)
* `acib217`       : CIB contamination at 217x217           (unused in EE but required)
* `asz143`        : tSZ contamination at 143x143           (unused in EE but required)
* `psr`           : Point source correlation coefficient between the 143- and 217-GHz         (unused in EE but required)
* `cibr`          : CIB correlation coefficient between the 143- and 217-GHz         (set to negative value for default) (unused in EE but required)
* `ncib143`       : CIB index at 143-GHz (set to negative value for default) (unused in EE but required)
* `ncib`          : CIB index (set to 0 for default) (unused in EE but required)
* `cibrun`        : CIB running (set to 0 for default)            (unused in EE but required)
* `xi`            : CIB-SZ correlation (unused in EE but required)
* `aksz`          : kSZ contamination at 143x143            (unused in EE but required)
* `dust100`       : Galactic contamination  normalization template (set to 1 for default)             (unused in EE but required)
* `dust143`       : Galactic contamination normalization template (set to 1 for default)             (unused in EE but required)
* `dust217`       : Galactic contamination normalization template (set to 1 for default)             (unused in EE but required)
* `dust143x217`   : Galactic contamination normalization template (set to 1 for default)                 (unused in EE but required)
* `calPlanck`     : Planck calibration
* `cal0`          : 100GHz TT calibration (unused in EE but required)
* `cal1`          : 143GHz TT calibration (unused in EE but required)
* `cal2`          : 217GHz TT calibration (unused in EE but required)
* `calTE`         : TE spectrum recalibration  (unused in EE but required)
* `calEE`         : EE spectrum recalibration 
The recommended values and priors are available in the file `plc_3.0/hi_l/camspec/camspec_recommendations.txt`

## TE only 
This file allows for the computation of the CMB EE likelihood in the range &#8467;=30-2500.

### Production process
The file contains the coadded 100-GHz, 143-GHz, and 217-GHz half-mission <i>TE</i> cross-spectra. 
Masks are based on the 353-GHz Q and U maps. 
Each of the different cross spectrum is corrected for dust contamination following an hybrid method using a map based correction up to  &#8467;=300 and a power spectrum template based one for higher multipoles. Details of teh dust correction, as well as multipole ranges retained for each cross spectra is described in the Planck 2018 likelihood paper.
Individual power spectra are corrected for polarization efficiency errors following the so called "spectrum-based" procedure also described in the paper.
Temperature to polarization leakages are corrected using the beam matrices computed for the specific mask and data cuts and using the temperature only 2015 best fit model.
The corrected cross spectra are coadded using the full TE covariance matrix. 

### Inputs
* Planck 100-, 143-, and 217-GHz half-mission <i>T</i> <i>E</i> maps;
* Planck 353-GHz maps for the dust residual contamination correction and the mask;
* CMB-cleaned (using the SMICA map) 353-GHz map, CO emission maps, and Planck catalogues for the Temperature masks;
* Covariance matrix is the full TTTEEE one restricted to the TE case at runtime.

### File name and usage
The camspec TE only likelihood is distributed in 
`plc_3.0/hi_l/camspec/camspec_10.7HM_1400_TE.clik`.

When used with the library, this expects a vector of parameters consisting of the 
<i>TE</i> CMB power spectrum from &#8467;=0 to 2500 (inclusive), followed 23 nuisance parameter.


Nuisance parameters are:

* `aps100`        : Point source contamination at 100x100          (unused in TE but required) 
* `aps143`        : Point source contamination at 143x143          (unused in TE but required)
* `aps217`        : Point source contamination at 217x217          (unused in TE but required)
* `acib143`       : CIB contamination at 143x143           (set to negative value for default) (unused in TE but required)
* `acib217`       : CIB contamination at 217x217           (unused in TE but required)
* `asz143`        : tSZ contamination at 143x143           (unused in TE but required)
* `psr`           : Point source correlation coefficient between the 143- and 217-GHz         (unused in TE but required)
* `cibr`          : CIB correlation coefficient between the 143- and 217-GHz         (set to negative value for default) (unused in TE but required)
* `ncib143`       : CIB index at 143-GHz (set to negative value for default) (unused in TE but required)
* `ncib`          : CIB index (set to 0 for default) (unused in TE but required)
* `cibrun`        : CIB running (set to 0 for default)            (unused in TE but required)
* `xi`            : CIB-SZ correlation (unused in TE but required)
* `aksz`          : kSZ contamination at 143x143            (unused in TE but required)
* `dust100`       : Galactic contamination  normalization template (set to 1 for default)             (unused in TE but required)
* `dust143`       : Galactic contamination normalization template (set to 1 for default)             (unused in TE but required)
* `dust217`       : Galactic contamination normalization template (set to 1 for default)             (unused in TE but required)
* `dust143x217`   : Galactic contamination normalization template (set to 1 for default)                 (unused in TE but required)
* `calPlanck`     : Planck calibration
* `cal0`          : 100GHz TT calibration (unused in TE but required)
* `cal1`          : 143GHz TT calibration (unused in TE but required)
* `cal2`          : 217GHz TT calibration (unused in TE but required)
* `calTE`         : TE spectrum recalibration 
* `calEE`         : EE spectrum recalibration (unused in TE but required)

The recommended values and priors are available in the file `plc_3.0/hi_l/camspec/camspec_recommendations.txt`

## TTTEEE  
This file allows for the computation of the CMB TT TE and EE likelihood in the range &#8467;=30-2500 for TT and TE and &#8467;=30-2000 for EE.

### Production process
See the TT, TE and EE only likelihoods

### Inputs
See the TT, TE and EE only likelihoods

### File name and usage
The camspec TTTEEE only likelihood is distributed in 
`plc_3.0/hi_l/camspec/camspec_10.7HM_1400_TTTEEE.clik`.

When used with the library, this expects a vector of parameters consisting of the 
<i>TT</i> CMB power spectrum from &#8467;=0 to 2500 (inclusive), followed by the <i>EE</i> CMB power spectrum from &#8467;=0 to 2500 (inclusive), the <i>TE</i> CMB power spectrum from &#8467;=0 to 2500 (inclusive), and 23 nuisance parameters.

Nuisance parameters are:

* `aps100`        : Point source contamination at 100x100          
* `aps143`        : Point source contamination at 143x143          
* `aps217`        : Point source contamination at 217x217          
* `acib143`       : CIB contamination at 143x143           (set to negative value for default)
* `acib217`       : CIB contamination at 217x217           (unused in TE but required)
* `asz143`        : tSZ contamination at 143x143           
* `psr`           : Point source correlation coefficient between the 143- and 217-GHz         
* `cibr`          : CIB correlation coefficient between the 143- and 217-GHz         (set to negative value for default) 
* `ncib143`       : CIB index at 143-GHz (set to negative value for default) 
* `ncib`          : CIB index (set to 0 for default) 
* `cibrun`        : CIB running (set to 0 for default)            
* `xi`            : CIB-SZ correlation
* `aksz`          : kSZ contamination at 143x143            
* `dust100`       : Galactic contamination  normalization template (set to 1 for default)            
* `dust143`       : Galactic contamination normalization template (set to 1 for default)             
* `dust217`       : Galactic contamination normalization template (set to 1 for default)        
* `dust143x217`   : Galactic contamination normalization template (set to 1 for default)                 
* `calPlanck`     : Planck calibration
* `cal0`          : 100GHz TT calibration 
* `cal1`          : 143GHz TT calibration 
* `cal2`          : 217GHz TT calibration 
* `calTE`         : TE spectrum recalibration 
* `calEE`         : EE spectrum recalibration 

The recommended values and priors are available in the file `plc_3.0/hi_l/camspec/camspec_recommendations.txt`

## TTEE, TEEE and TTTE
3 extra files allow to compute the TTEE, TEEE and TTTE likelihoods.

  
### File names and usage
The camspec TTEE, TEEE and TTTE only likelihood are available in 
`plc_3.0/hi_l/camspec/camspec_10.7HM_1400_TTEE_b.clik`, 
`plc_3.0/hi_l/camspec/camspec_10.7HM_1400_TEEE_b.clik` and 
`plc_3.0/hi_l/camspec/camspec_10.7HM_1400_TTTE_b.clik`.

They contain the full covariance which is recomputed to the specific case at runtime.
