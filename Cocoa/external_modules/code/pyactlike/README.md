# ACTPol DR4 CMB Power Spectrum Likelihood

![Python package](https://github.com/ACTCollaboration/pyactlike/workflows/Python%20package/badge.svg)

This is the **Data Release 4 (DR4)** CMB power spectrum likelihood measured by the Atacama Cosmology Telescope (ACT), from the 2013–2016 survey covering >15,000 sq. deg. This spectrum has already been marginalized over SZ and foreground emission. The polarization efficiency is the only nuisance parameter required to be sampled. It is based on the WMAP and ACT team's likelihood software.

**Cite:** Aiola et al. 2020, Choi et al. 2020. This package is based off of the Fortran implementation written by Erminia Calabrese and Jo Dunkley. Thanks to Tim Morton for helping interface with Cobaya, and special thanks to Leander Thiele for finding a bug in the Monte Python 3 wrapper.

<img src="https://act.princeton.edu/sites/act/files/styles/panopoly_image_original/public/media/angelapano.jpg" 
alt="panoramic image"/></a>

## Installation
To install, clone this repository and install it using pip.
```bash
git clone https://github.com/ACTCollaboration/pyactlike
cd pyactlike
pip install . --user
```

## Usage

This package can be called standalone or with [Cobaya][1]. For experimental MontePython support,
see the `montepython3/` directory in this repository. If you are using a Cobaya version &lt; 2.1.0 (the current stable branch), please see the example Jupyter notebook in
`notebooks/Example for Cobaya (stable).ipynb`.

If you are on Cobaya 2.1.0 
(currently the devel branch), using the likelihood is as easy as including it in your YAML 
or configuration dict. There's one nuisance parameter, the polarization efficiency `yp2`.

[1]: https://github.com/CobayaSampler/cobaya
[2]: https://github.com/brinckmann/montepython_public

```
likelihood:
    pyactlike.ACTPol_lite_DR4:

params:   
    yp2:
        prior:
            min: 0.5
            max: 1.5     
```

We provide additional likelihood configurations for common uses, 
* `pyactlike.ACTPol_lite_DR4_for_combining_with_Planck` for TT,TE,EE excluding the large scale ACT temperature
* `pyactlike.ACTPol_lite_DR4_onlyTT` for restricting only to TT
* `pyactlike.ACTPol_lite_DR4_onlyTE` for restricting only to TE
* `pyactlike.ACTPol_lite_DR4_onlyEE` for restricting only to EE


## Tests
Run `pytest` in the repository base directory run the tests.
