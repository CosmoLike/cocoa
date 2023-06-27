#!/usr/bin/env python3
# Example script to create fiducial values for mock CMB likelihoods
# with non-LCDM fiducial, namely with small-scale baryon clumping
from cobaya.model import get_model

# from best fit to Planck+SH0ES
# with fixed massless neutrinos and nuisance-marginalized high-l
fiducial_params = {
    # LambdaCDM parameters
    'H0': 70.91588,
    # '100*theta_s': 1.041771159e+00,
    'omega_b': 2.271421832e-02,
    'N_ur': 3.046,  # three massless neutrinos
    'omega_cdm': 1.199919146e-01,
    'A_s': 2.113151105e-09,
    # 'sigma8': 8.334149150e-01,
    'n_s': 9.655213736e-01,
    'tau_reio': 5.420878071e-02,
    # clumping space
    'delta_m': -9.549925860e-01,
    'delta_p': 1.319526438e+00,
    'f2V': 6.516101687e-01
}

fiducial_params_extra = {
    'recombination': 'recfast_3zones_lowlevel',
    'non linear': 'halofit'
}

fiducial_params_full = fiducial_params.copy()
fiducial_params_full.update(fiducial_params_extra)

info_fiducial = {
    'params': fiducial_params,
    'likelihood': {'cobaya_mock_cmb.MockSOClumping': {'python_path': '.'},
                   'cobaya_mock_cmb.MockCMBS4Clumping': {'python_path': '.'}},
    'theory': {'classy': {"extra_args": fiducial_params_extra}}}

model_fiducial = get_model(info_fiducial)

model_fiducial.logposterior({})

Cls = model_fiducial.provider.get_Cl(units="muK2")

for likelihood in model_fiducial.likelihood.values():
    likelihood.create_fid_values(Cls, fiducial_params_full, override=True)
