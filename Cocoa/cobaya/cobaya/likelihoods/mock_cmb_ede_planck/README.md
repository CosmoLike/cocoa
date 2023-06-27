# cobaya_mock_cmb
Mock CMB likelihood class for [Cobaya sampler](https://github.com/CobayaSampler/cobaya), and several specific experiment examples.

Introduced in [Rashkovetskyi et al 2021](https://arxiv.org/abs/2108.02747), so please cite it if you use the code.

## Python package contents (classes)

* MockCMBLikelihood - base mock CMB likelihood class. Ported [MontePython's 3.3 Likelihood_mock_cmb](https://github.com/brinckmann/montepython_public/blob/master/montepython/likelihood_class.py) with the help of [Cobaya's example of likelihood class](https://cobaya.readthedocs.io/en/latest/cosmo_external_likelihood_class.html), most loops replaced with Numpy vectorized operations. Unlensed Cl's supported, but note that Cobaya theory class wrappers didn't provide them [until Jan 2021](https://github.com/CobayaSampler/cobaya/commit/6af50e519a294b2bde82257d88bf636710c1ea64).
* MockSO - Simons Observatory (SO) model following [Sailer, Schaan and Ferraro 2020](https://arxiv.org/abs/2007.04325), based on MontePython config by [Julian Munoz](https://github.com/JulianBMunoz).
* MockSOBaseline - Simons Observatory (SO) with TT, EE deproj0 [noise curves](https://simonsobservatory.org/assets/supplements/20180822_SO_Noise_Public.tgz) (baseline sensitivity).
* MockSOGoal - Simons Observatory (SO) with TT, EE deproj0 [noise curves](https://simonsobservatory.org/assets/supplements/20180822_SO_Noise_Public.tgz) (goal sensitivity).
* MockCMBS4 - CMB-S4 model following [the science book](https://arxiv.org/abs/1907.04473), based on MontePython config by [Julian Munoz](https://github.com/JulianBMunoz).
* MockCMBS4sens0 - CMB-S4 with TT, EE deproj0 [noise curves](http://sns.ias.edu/~jch/S4_190604d_2LAT_Tpol_default_noisecurves.tgz).
* MockPlanck - Planck model following [Munoz et al 2016](https://arxiv.org/abs/1611.05883) with `f_sky=0.2` (fraction independent of SO and CMB-S4).

## External contents
* make_fiducial.py - example script to generate fiducial power spectra for the experiments
* mock_test.yaml - example config to run Cobaya with all the likelihoods (`cobaya-run mock_test.yaml` in this directory) to make sure they work. One can't simply combine all in real runs, because they cover the same sky.