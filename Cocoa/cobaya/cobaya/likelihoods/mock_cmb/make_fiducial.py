#!/usr/bin/env python3
# Example script to create fiducial values for mock CMB likelihoods
from cobaya.model import get_model

# from best fit with fixed massless neutrinos and nuisance-marginalized high-l
fiducial_params = {
    # LambdaCDM parameters
    'H0': 68.1461,
    # '100*theta_s': 1.041920539e+00,
    'omega_b': 2.241606328e-02,
    'N_ur': 3.046,  # three massless neutrinos
    'omega_cdm': 1.194462584e-01,
    'A_s': 2.109371924e-09,
    # 'sigma8': 8.245006041e-01,
    'n_s': 9.660360599e-01,
    'tau_reio': 5.142494234e-02
}

fiducial_params_extra = {
    'recombination': 'recfast',
    'non linear': 'halofit'
}

fiducial_params_full = fiducial_params.copy()
fiducial_params_full.update(fiducial_params_extra)

info_fiducial = {
    'params': fiducial_params,
    'likelihood': {'cobaya_mock_cmb.MockSO': {'python_path': '.'},
                   'cobaya_mock_cmb.MockSOBaseline': {'python_path': '.'},
                   'cobaya_mock_cmb.MockSOGoal': {'python_path': '.'},
                   'cobaya_mock_cmb.MockCMBS4': {'python_path': '.'},
                   'cobaya_mock_cmb.MockCMBS4sens0': {'python_path': '.'},
                   'cobaya_mock_cmb.MockPlanck': {'python_path': '.'}},
    'theory': {'classy': {"extra_args": fiducial_params_extra}}}

model_fiducial = get_model(info_fiducial)

model_fiducial.logposterior({})

Cls = model_fiducial.provider.get_Cl(units="muK2")

for likelihood in model_fiducial.likelihood.values():
    likelihood.create_fid_values(Cls, fiducial_params_full, override=True)
