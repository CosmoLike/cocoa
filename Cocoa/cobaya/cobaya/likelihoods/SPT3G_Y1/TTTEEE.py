from cobaya.likelihoods.base_classes import planck_clik
from cobaya.log import LoggedError, get_logger
from cobaya.tools import are_different_params_lists
import numpy as np
from typing import Optional
import os



class TTTEEE(planck_clik.PlanckClik): 

    SPT3G_2018_TTTEEE_late_crop_msk  : Optional[str]
    SPT3G_2018_TTTEEE_beam_covariance_scale  : Optional[str]
    SPT3G_2018_TTTEEE_galdust_T  : Optional[str]
    SPT3G_2018_TTTEEE_CIB_T  : Optional[str]
    SPT3G_2018_TTTEEE_tSZ_cosmology_scaling  : Optional[str]
    SPT3G_2018_TTTEEE_kSZ_cosmology_scaling  : Optional[str]
    SPT3G_2018_TTTEEE_spectra_to_fit_bin_min  : Optional[str]
    SPT3G_2018_TTTEEE_spectra_to_fit_bin_max  : Optional[str]
    SPT3G_2018_TTTEEE_spectra_to_fit : Optional[str]
    
    def initialize(self):
        try:
            install_path = (
                lambda p: self.get_code_path(p) if p else None)(self.packages_path)
            # min_version here is checked inside get_clik_import_path, since it is
            # displayed in the folder name and cannot be retrieved from the module.
            clik = planck_clik.load_clik(
                "clik", path=self.path, install_path=install_path,
                get_import_path=lambda pth:os.path.join(planck_clik.get_clik_source_folder(pth), 'lib/python/site-packages'), logger=self.log,
                not_installed_level="debug")
        except planck_clik.VersionCheckError as excpt:
            raise planck_clik.VersionCheckError(
                str(excpt) + " Install new clik version following indications at https://github.com/SouthPoleTelescope/spt3g_y1_dist")
        except ComponentNotInstalledError as excpt:
            raise ComponentNotInstalledError(
                self.log, (f"Could not find clik: {excpt}. "
                           "To install follow indications at https://github.com/SouthPoleTelescope/spt3g_y1_dist"))
        if int(clik.version().split("_")[1].split(".")[0])<16:
            raise planck_clik.VersionCheckError("SPT3G likelihood requires clik v16+. See information here https://github.com/SouthPoleTelescope/spt3g_y1_dist")
        # Loading the likelihood data
        data_path = planck_clik.get_data_path(self.__class__.get_qualified_class_name())
        if not os.path.isabs(self.clik_file):
            self.path_data = getattr(self, "path_data", os.path.join(
                self.path or self.packages_path, "data", data_path))
            self.clik_file = os.path.join(self.path_data, self.clik_file)
        
        # get the options
        options_list = [v for v in dir(self) if v.startswith("SPT3G_2018_TTTEEE")]
        options = dict([(op,str(getattr(self,op))) for op in options_list if str(getattr(self,op))!=""])
        try:
            self.clik = clik.clik(self.clik_file,**options)
        except clik.lkl.CError:
            # Is it that the file was not found?
            if not os.path.exists(self.clik_file):
                raise ComponentNotInstalledError(
                    self.log, "The .clik file was not found where specified in the "
                              "'clik_file' field of the settings of this likelihood. "
                              "Maybe the 'path' given is not correct? The full path where"
                              " the .clik file was searched for is '%s'", self.clik_file)
            # Else: unknown clik error
            self.log.error("An unexpected error occurred in clik (possibly related to "
                           "multiple simultaneous initialization, or simultaneous "
                           "initialization of incompatible likelihoods; e.g. polarised "
                           "vs non-polarised 'lite' likelihoods. See error info below:")
            raise
        self.l_maxs = self.clik.get_lmax()
        # calculate requirements here so class can also be separately instantiated
        requested_cls = ["tt", "ee", "bb", "te", "tb", "eb"]
        has_cl = self.clik.get_has_cl()
        self.requested_cls = [cl for cl, i in zip(requested_cls, has_cl) if int(i)]
        self.l_maxs_cls = [lmax for lmax, i in zip(self.l_maxs, has_cl) if int(i)]
        self.expected_params = list(self.clik.extra_parameter_names)
        # Placeholder for vector passed to clik
        length = (len(self.clik.get_has_cl()))
        self.vector = np.zeros(np.sum(self.l_maxs) + length + len(self.expected_params))

    def log_likelihood(self, cl, **params_values):
        # fill with Cl's
        self.vector[:-len(self.expected_params)] = np.concatenate(
            [(cl[spectrum][:1 + lmax] if spectrum not in ["tb", "eb"]
              else np.zeros(1 + lmax))
             for spectrum, lmax in zip(self.requested_cls, self.l_maxs_cls)])
        # check for nan's: mey produce a segfault in clik
        # dot product is apparently the fastest way in threading-enabled numpy
        if np.isnan(np.dot(self.vector, self.vector)):
            return -np.inf
        # fill with likelihood parameters
        #first the nuisance
        self.vector[-len(self.expected_params):] = (
            [params_values[p] for p in self.expected_params])
        
        loglike = self.clik(self.vector)[0]
        # "zero" of clik, and sometimes nan's returned
        if np.allclose(loglike, -1e30) or np.isnan(loglike):
            loglike = -np.inf
        return loglike
    

_planck_get_data_path = planck_clik.get_data_path

def get_data_path(name):
    log = get_logger(name)
    if "spt" not in name.lower():
        return _planck_get_data_path(name)
    log.info("override default get_data_path from %s"%(_planck_get_data_path.__module__))
    return "spt_data"

planck_clik.get_data_path = get_data_path
