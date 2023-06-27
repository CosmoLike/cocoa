## Installation
To use this likelihood with Monte Python, first install the `pyactlike` package as usual.
Then copy contents of the `likelihoods/` directory to the MontePython
likelihoods folder. These represent the likelihood itself, and the suggested prior on the 
optical depth to reionization to be used with the likelihood. i.e.

```bash
cp -r likelihoods/* PATH_TO_MONTEPYTHON/montepython/likelihoods/
```

## Usage
An example param file is included, `ACTPol_lite_DR4.param`. If one is sampling with just
the ACT likelihood, we recommend the use of a prior on tau, i.e.

```python
data.experiments=['ACTPol_lite_DR4', 'tau_prior']
```

As in the Cobaya interface, there is one nuisance parameter, the polarization efficiency `yp2`,
which should go in the parameters section.

```python
data.parameters['yp2']         = [0.97,    0.9, 1.1,   0.015, 1,  'nuisance'] 
```

We recommend the use of halofit for the nonlinear matter power spectrum.
```python
data.cosmo_arguments['non linear'] = 'halofit'
```
