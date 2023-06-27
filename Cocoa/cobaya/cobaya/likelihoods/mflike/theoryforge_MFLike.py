import os
from typing import Union

import numpy as np


# Converts from cmb units to brightness. Numerical factors not included, it needs proper normalization when used.
def _cmb2bb(nu):
    # NB: numerical factors not included
    from scipy import constants

    T_CMB = 2.72548
    x = nu * constants.h * 1e9 / constants.k / T_CMB
    return np.exp(x) * (nu * x / np.expm1(x)) ** 2


# Provides the frequency value given the bandpass name. To be modified - it is ACT based!!!!!!
def _get_fr(array):
    a = array.split("_")[0]
    if a == "PA1" or a == "PA2":
        fr = 150
    if a == "PA3":
        fr = array.split("_")[3]
    if a == "PA4" or a == "PA5" or a == "PA6":
        fr = array.split("_")[3].replace("f", "")
    return np.int(fr), a


def _get_arrays_weights(arrays, polarized_arrays, freqs):
    array_weights = {}
    counter = []
    for array in arrays:
        fr, pa = _get_fr(array)
        if (pa in polarized_arrays) and (fr in freqs):
            if fr not in counter:
                counter.append(fr)
                array_weights[fr] = 1
            else:
                array_weights[fr] = array_weights[fr] + 1
    return array_weights


class TheoryForge_MFLike:
    def __init__(
        self,
        mflike=None,
    ):

        if mflike is None:
            self.data_folder = None
            self.freqs = np.array([145, 225, 93])
            self.foregrounds = {
                "normalisation": {"nu_0": 150.0, "ell_0": 3000, "T_CMB": 2.725},
                "components": {
                    "tt": ["kSZ", "tSZ_and_CIB", "cibp", "dust", "radio"],
                    "te": ["radio", "dust"],
                    "ee": ["radio", "dust"],
                },
            }
            self.l_bpws = np.arange(2, 6002)
            self.requested_cls = ["tt", "te", "ee"]
            self.bandint_nsteps = 1
            self.bandint_width = 0

            nuis_params = {
                "bandint_shift_93": 0,
                "bandint_shift_145": 0,
                "bandint_shift_225": 0,
            }

            self.bandint_freqs = self._bandpass_construction(**nuis_params)
        else:

            self.data_folder = mflike.data_folder
            self.freqs = mflike.freqs
            self.foregrounds = mflike.foregrounds
            self.l_bpws = mflike.l_bpws
            self.requested_cls = mflike.requested_cls
            self.expected_params_fg = mflike.expected_params_fg
            self.expected_params_nuis = mflike.expected_params_nuis
            self.spec_meta = mflike.spec_meta
            self.defaults_cuts = mflike.defaults

            # Initialize foreground model
            self._init_foreground_model()

            # Parameters for template from file
            self.systematics_template = mflike.systematics_template
            # Initialize template for marginalization, if needed
            if self.systematics_template["has_file"]:
                self._init_template_from_file()

            # Parameters for band integration
            self.bandint_nsteps = mflike.band_integration["nsteps"]
            self.bandint_width = mflike.band_integration["bandwidth"]
            self.bandint_external_bandpass = mflike.band_integration["external_bandpass"]
            self.polarized_arrays = mflike.band_integration["polarized_arrays"]

            # Bandpass construction for band integration
            if self.bandint_external_bandpass:
                path = os.path.normpath(os.path.join(self.data_folder, "external_bandpasses/"))
                arrays = os.listdir(path)
                self._init_external_bandpass_construction(path, arrays)
                self.array_weights = _get_arrays_weights(arrays, self.polarized_arrays, self.freqs)

    def get_modified_theory(self, Dls, **params):

        fg_params = {k: params[k] for k in self.expected_params_fg}
        nuis_params = {k: params[k] for k in self.expected_params_nuis}

        # Bandpass construction for band integration
        # Bandpasses are regenerated at each iteration, even if
        # "bandint_shift_{ch}" is not explored in the MCMC
        # this might be optimized in the future
        if self.bandint_external_bandpass:
            order, self.bandint_freqs = self._external_bandpass_construction(**nuis_params)
            fg_dict = self._get_foreground_model(freqs_order=order, **fg_params)

            for f1 in self.freqs:
                for f2 in self.freqs:
                    for s in self.requested_cls:
                        fg_dict[s, "all", f1, f2] = np.zeros(len(self.l_bpws))
                        for pa in self.polarized_arrays:
                            o1 = pa + "_" + str(f1)
                            o2 = pa + "_" + str(f2)
                            if (o1 in order) and (o2 in order):
                                fg_dict[s, "all", f1, f2] += (
                                    fg_dict[s, "all", o1, o2]
                                    / self.array_weights[f1]
                                    / self.array_weights[f2]
                                )

        else:
            self.bandint_freqs = self._bandpass_construction(**nuis_params)
            fg_dict = self._get_foreground_model(**fg_params)

        cmbfg_dict = {}
        # Sum CMB and FGs
        for f1 in self.freqs:
            for f2 in self.freqs:
                for s in self.requested_cls:
                    cmbfg_dict[s, f1, f2] = Dls[s] + fg_dict[s, "all", f1, f2]

        # Apply alm based calibration factors
        cmbfg_dict = self._get_calibrated_spectra(cmbfg_dict, **nuis_params)

        # Introduce spectra rotations
        cmbfg_dict = self._get_rotated_spectra(cmbfg_dict, **nuis_params)

        # Introduce templates of systematics from file, if needed
        if self.systematics_template["has_file"]:
            cmbfg_dict = self._get_template_from_file(cmbfg_dict, **nuis_params)

        # Built theory
        dls_dict = {}
        for m in self.spec_meta:
            p = m["pol"]
            if p in ["tt", "ee", "bb"]:
                dls_dict[p, m["nu1"], m["nu2"]] = cmbfg_dict[p, m["nu1"], m["nu2"]]
            else:  # ['te','tb','eb']
                if m["hasYX_xsp"]:  # not symmetrizing
                    dls_dict[p, m["nu1"], m["nu2"]] = cmbfg_dict[p, m["nu2"], m["nu1"]]
                else:
                    dls_dict[p, m["nu1"], m["nu2"]] = cmbfg_dict[p, m["nu1"], m["nu2"]]

                if self.defaults_cuts["symmetrize"]:  # we average TE and ET (as we do for data)
                    dls_dict[p, m["nu1"], m["nu2"]] += cmbfg_dict[p, m["nu2"], m["nu1"]]
                    dls_dict[p, m["nu1"], m["nu2"]] *= 0.5

        return dls_dict

    ###########################################################################
    ## This part deals with foreground construction and bandpass integration ##
    ###########################################################################

    # Initializes the foreground model. It sets the SED and reads the templates
    def _init_foreground_model(self):

        from fgspectra import cross as fgc
        from fgspectra import frequency as fgf
        from fgspectra import power as fgp

        template_path = os.path.join(os.path.dirname(os.path.abspath(fgp.__file__)), "data")
        cibc_file = os.path.join(template_path, "cl_cib_Choi2020.dat")

        # set pivot freq and multipole
        self.fg_nu_0 = self.foregrounds["normalisation"]["nu_0"]
        self.fg_ell_0 = self.foregrounds["normalisation"]["ell_0"]

        # We don't seem to be using this
        # cirrus = fgc.FactorizedCrossSpectrum(fgf.PowerLaw(), fgp.PowerLaw())
        self.ksz = fgc.FactorizedCrossSpectrum(fgf.ConstantSED(), fgp.kSZ_bat())
        self.cibp = fgc.FactorizedCrossSpectrum(fgf.ModifiedBlackBody(), fgp.PowerLaw())
        self.radio = fgc.FactorizedCrossSpectrum(fgf.PowerLaw(), fgp.PowerLaw())
        self.tsz = fgc.FactorizedCrossSpectrum(fgf.ThermalSZ(), fgp.tSZ_150_bat())
        self.cibc = fgc.FactorizedCrossSpectrum(fgf.CIB(), fgp.PowerSpectrumFromFile(cibc_file))
        self.dust = fgc.FactorizedCrossSpectrum(fgf.ModifiedBlackBody(), fgp.PowerLaw())
        self.tSZ_and_CIB = fgc.SZxCIB_Choi2020()

        components = self.foregrounds["components"]
        self.fg_component_list = {s: components[s] for s in self.requested_cls}

    # Gets the actual power spectrum of foregrounds given the passed parameters
    def _get_foreground_model(self, ell=None, freqs_order=None, **fg_params):
        # if ell = None, it uses the l_bpws, otherwise the ell array provided
        # useful to make tests at different l_max than the data
        if not hasattr(ell, "__len__"):
            ell = self.l_bpws
        ell_0 = self.fg_ell_0
        nu_0 = self.fg_nu_0

        # Normalisation of radio sources
        ell_clp = ell * (ell + 1.0)
        ell_0clp = ell_0 * (ell_0 + 1.0)

        model = {}
        model["tt", "kSZ"] = fg_params["a_kSZ"] * self.ksz(
            {"nu": self.bandint_freqs}, {"ell": ell, "ell_0": ell_0}
        )
        model["tt", "cibp"] = fg_params["a_p"] * self.cibp(
            {
                "nu": self.bandint_freqs,
                "nu_0": nu_0,
                "temp": fg_params["T_d"],
                "beta": fg_params["beta_p"],
            },
            {"ell": ell_clp, "ell_0": ell_0clp, "alpha": 1},
        )
        model["tt", "radio"] = fg_params["a_s"] * self.radio(
            {"nu": self.bandint_freqs, "nu_0": nu_0, "beta": -0.5 - 2.0},
            {"ell": ell_clp, "ell_0": ell_0clp, "alpha": 1},
        )
        model["tt", "tSZ"] = fg_params["a_tSZ"] * self.tsz(
            {"nu": self.bandint_freqs, "nu_0": nu_0},
            {"ell": ell, "ell_0": ell_0},
        )
        model["tt", "cibc"] = fg_params["a_c"] * self.cibc(
            {
                "nu": self.bandint_freqs,
                "nu_0": nu_0,
                "temp": fg_params["T_d"],
                "beta": fg_params["beta_c"],
            },
            {"ell": ell, "ell_0": ell_0},
        )
        model["tt", "dust"] = fg_params["a_gtt"] * self.dust(
            {"nu": self.bandint_freqs, "nu_0": nu_0, "temp": 19.6, "beta": 1.5},
            {"ell": ell, "ell_0": 500.0, "alpha": -0.6},
        )
        model["tt", "tSZ_and_CIB"] = self.tSZ_and_CIB(
            {
                "kwseq": (
                    {"nu": self.bandint_freqs, "nu_0": nu_0},
                    {
                        "nu": self.bandint_freqs,
                        "nu_0": nu_0,
                        "temp": fg_params["T_d"],
                        "beta": fg_params["beta_c"],
                    },
                )
            },
            {
                "kwseq": (
                    {"ell": ell, "ell_0": ell_0, "amp": fg_params["a_tSZ"]},
                    {"ell": ell, "ell_0": ell_0, "amp": fg_params["a_c"]},
                    {
                        "ell": ell,
                        "ell_0": ell_0,
                        "amp": -fg_params["xi"] * np.sqrt(fg_params["a_tSZ"] * fg_params["a_c"]),
                    },
                )
            },
        )

        model["ee", "radio"] = fg_params["a_psee"] * self.radio(
            {"nu": self.bandint_freqs, "nu_0": nu_0, "beta": -0.5 - 2.0},
            {"ell": ell_clp, "ell_0": ell_0clp, "alpha": 1},
        )
        model["ee", "dust"] = fg_params["a_gee"] * self.dust(
            {"nu": self.bandint_freqs, "nu_0": nu_0, "temp": 19.6, "beta": 1.5},
            {"ell": ell, "ell_0": 500.0, "alpha": -0.4},
        )

        model["te", "radio"] = fg_params["a_pste"] * self.radio(
            {"nu": self.bandint_freqs, "nu_0": nu_0, "beta": -0.5 - 2.0},
            {"ell": ell_clp, "ell_0": ell_0clp, "alpha": 1},
        )
        model["te", "dust"] = fg_params["a_gte"] * self.dust(
            {"nu": self.bandint_freqs, "nu_0": nu_0, "temp": 19.6, "beta": 1.5},
            {"ell": ell, "ell_0": 500.0, "alpha": -0.4},
        )

        fg_dict = {}
        if not hasattr(freqs_order, "__len__"):
            freqs = self.freqs
        else:
            freqs = freqs_order
        for c1, f1 in enumerate(freqs):
            for c2, f2 in enumerate(freqs):
                for s in self.requested_cls:
                    fg_dict[s, "all", f1, f2] = np.zeros(len(ell))
                    for comp in self.fg_component_list[s]:
                        if comp == "tSZ_and_CIB":
                            fg_dict[s, "tSZ", f1, f2] = model[s, "tSZ"][c1, c2]
                            fg_dict[s, "cibc", f1, f2] = model[s, "cibc"][c1, c2]
                            fg_dict[s, "tSZxCIB", f1, f2] = (
                                model[s, comp][c1, c2]
                                - model[s, "tSZ"][c1, c2]
                                - model[s, "cibc"][c1, c2]
                            )
                            fg_dict[s, "all", f1, f2] += model[s, comp][c1, c2]
                        else:
                            fg_dict[s, comp, f1, f2] = model[s, comp][c1, c2]
                            fg_dict[s, "all", f1, f2] += fg_dict[s, comp, f1, f2]

        return fg_dict

    # Takes care of the bandpass construction. It returns a list of nu-transmittance for each frequency or an array with the effective freqs.
    def _bandpass_construction(self, **params):

        if not hasattr(self.bandint_width, "__len__"):
            self.bandint_width = np.full_like(self.freqs, self.bandint_width, dtype=float)
        if np.any(np.array(self.bandint_width) > 0):
            assert self.bandint_nsteps > 1, "bandint_width and bandint_nsteps not coherent"
            assert np.all(
                np.array(self.bandint_width) > 0
            ), "one band has width = 0, set a positive width and run again"

            bandint_freqs = []
            for ifr, fr in enumerate(self.freqs):
                bandpar = "bandint_shift_" + str(fr)
                bandlow = fr * (1 - self.bandint_width[ifr] * 0.5)
                bandhigh = fr * (1 + self.bandint_width[ifr] * 0.5)
                nubtrue = np.linspace(bandlow, bandhigh, self.bandint_nsteps, dtype=float)
                nub = np.linspace(
                    bandlow + params[bandpar],
                    bandhigh + params[bandpar],
                    self.bandint_nsteps,
                    dtype=float,
                )
                tranb = _cmb2bb(nub)
                tranb_norm = np.trapz(_cmb2bb(nubtrue), nubtrue)
                bandint_freqs.append([nub, tranb / tranb_norm])
        else:
            bandint_freqs = np.empty_like(self.freqs, dtype=float)
            for ifr, fr in enumerate(self.freqs):
                bandpar = "bandint_shift_" + str(fr)
                bandint_freqs[ifr] = fr + params[bandpar]

        return bandint_freqs

    def _init_external_bandpass_construction(self, path, arrays):
        self.external_bandpass = []
        for array in arrays:
            fr, pa = _get_fr(array)
            if pa in self.polarized_arrays and fr in self.freqs:
                nu_ghz, bp = np.loadtxt(path + "/" + array, usecols=(0, 1), unpack=True)
                trans_norm = np.trapz(bp * _cmb2bb(nu_ghz), nu_ghz)
                self.external_bandpass.append([pa, fr, nu_ghz, bp / trans_norm])

    def _external_bandpass_construction(self, **params):
        bandint_freqs = []
        order = []
        for pa, fr, nu_ghz, bp in self.external_bandpass:
            bandpar = "bandint_shift_" + str(fr)
            nub = nu_ghz + params[bandpar]
            trans = bp * _cmb2bb(nub)
            bandint_freqs.append([nub, trans])
            order.append(pa + "_" + str(fr))

        return order, bandint_freqs

    ###########################################################################
    ## This part deals with calibration factors
    ## Here we implement an alm based calibration
    ## Each field {T,E,B}{freq1,freq2,...,freqn} gets an independent
    ## calibration factor, e.g. calT_145, calE_154, calT_225, etc..
    ## A global calibration factor calG_all is also considered.
    ###########################################################################

    def _get_calibrated_spectra(self, dls_dict, **nuis_params):

        from syslibrary import syslib_mflike as syl

        cal_pars = {}
        if "tt" in self.requested_cls or "te" in self.requested_cls:
            cal_pars["tt"] = nuis_params["calG_all"] * np.array(
                [nuis_params["calT_" + str(fr)] for fr in self.freqs]
            )

        if "ee" in self.requested_cls or "te" in self.requested_cls:
            cal_pars["ee"] = nuis_params["calG_all"] * np.array(
                [nuis_params["calE_" + str(fr)] for fr in self.freqs]
            )

        calib = syl.Calibration_alm(ell=self.l_bpws, spectra=dls_dict)

        return calib(cal1=cal_pars, cal2=cal_pars, nu=self.freqs)

    ###########################################################################
    ## This part deals with rotation of spectra
    ## Each freq {freq1,freq2,...,freqn} gets a rotation angle alpha_93, alpha_145, etc..
    ###########################################################################

    def _get_rotated_spectra(self, dls_dict, **nuis_params):

        from syslibrary import syslib_mflike as syl

        rot_pars = [nuis_params["alpha_" + str(fr)] for fr in self.freqs]

        rot = syl.Rotation_alm(ell=self.l_bpws, spectra=dls_dict, cls=self.requested_cls)

        return rot(rot_pars, nu=self.freqs)

    ###########################################################################
    ## This part deals with template marginalization
    ## A dictionary of template dls is read from yaml (likely to be not efficient)
    ## then rescaled and added to theory dls
    ###########################################################################

    # Initializes the systematics templates
    # This is slow, but should be done only once
    def _init_template_from_file(self):

        from syslibrary import syslib_mflike as syl

        # decide where to store systematics template.
        # Currently stored inside syslibrary package
        templ_from_file = syl.ReadTemplateFromFile(rootname=self.systematics_template["rootname"])
        self.dltempl_from_file = templ_from_file(ell=self.l_bpws)

    def _get_template_from_file(self, dls_dict, **nuis_params):

        # templ_pars=[nuis_params['templ_'+str(fr)] for fr in self.freqs]
        # templ_pars currently hard-coded
        # but ideally should be passed as input nuisance
        templ_pars = {
            cls: np.zeros((len(self.freqs), len(self.freqs))) for cls in self.requested_cls
        }

        for cls in self.requested_cls:
            for i1, f1 in enumerate(self.freqs):
                for i2, f2 in enumerate(self.freqs):
                    dls_dict[cls, f1, f2] += (
                        templ_pars[cls][i1][i2] * self.dltempl_from_file[cls, f1, f2]
                    )

        return dls_dict
