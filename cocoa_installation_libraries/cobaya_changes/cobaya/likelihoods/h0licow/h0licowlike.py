################################# LIKELIHOOD FROM H0LICOW #####################################
#################################  COBAYA IMPLEMENTATION  #####################################
# João Victor Silva Rebouças, May 2022

import numpy as np
import pandas as pd
from cobaya.likelihood import Likelihood
from .lensutils import *

################################# DEFINING log(likelihood) #####################################
# It uses AstroPy to calculate distances, should change to CAMB
# Must provide: angular diameter distance

def log_like_add(lens, self):
    """
    Computes the relevant angular diameter distance(s) of a given lens in a given cosmology,
    and evaluate its/their joint likelihood against the same modeled distances of the lens.

    param lens: either a GLEELens or LenstronomyLens instance.
    """
    # JVR: old implementation with AstroPy
    # dd = cosmo.angular_diameter_distance(z=lens.zlens).value
    # ds = cosmo.angular_diameter_distance(z=lens.zsource).value
    # dds = cosmo.angular_diameter_distance_z1z2(z1=lens.zlens, z2=lens.zsource).value
    # ddt = (1. + lens.zlens) * dd * ds / dds
    # New implementation with CAMB
    z_pair = (lens.zlens, lens.zsource)
    dd = self.provider.get_angular_diameter_distance(z=[lens.zlens])
    ds = self.provider.get_angular_diameter_distance(z=[lens.zsource])
    dds = self.provider.get_angular_diameter_distance_2([z_pair])
    ddt = (1. + lens.zlens) * dd * ds / dds # time-delay distance

    if isinstance(lens, GLEELens):
        if lens.loglikelihood_type in ["kde_full", "kde_hist"]:
            # because the newest sklearn kde want arrays, not numbers... 
            return lens.loglikelihood(np.array(ddt).reshape(1, -1))
        elif lens.loglikelihood_type in ["kde_hist_2d"]:
            return lens.loglikelihood(np.array([dd, ddt]).reshape(1, -1))
        elif lens.loglikelihood_type in ["sklogn_analytical_Dd"] :
            return lens.loglikelihood(ddt, dd)
        elif lens.loglikelihood_type in ["sklogn_analytical_Ddonly"] :
            return lens.loglikelihood(dd)
        else:
            return lens.loglikelihood(ddt)

    elif isinstance(lens, LenstronomyLens):
        return lens.loglikelihood(dd, ddt)

    else:
        sys.exit("I don't know what to do with %s, unknown instance" % lens)

################################# DEFINING LIKELIHOOD CLASS #####################################

class h0licowlike(Likelihood):

    def initialize(self):
        """
         Prepare any computation, importing any necessary code, files, etc.
        """
        declare_lens_systems(self)
        print("[H0liCOW] Using strong lenses:", [lens.name for lens in self.lenses])
    
    def get_requirements(self):
        """
         Return dictionary specifying quantities calculated by a theory code are needed
        """
        lens_redshifts = []
        source_redshifts = []
        for lens in self.lenses:
            lens_redshifts.append(lens.zlens)
            source_redshifts.append(lens.zsource)
        lens_redshifts = np.array(lens_redshifts)
        source_redshifts = np.array(source_redshifts)
        redshift_pairs = np.array(list(zip(lens_redshifts, source_redshifts)))
        all_redshifts = np.concatenate((lens_redshifts, source_redshifts))
        return {'angular_diameter_distance': {'z': all_redshifts},
                'angular_diameter_distance_2':{'z_pairs': redshift_pairs}
               }

    def logp(self, **params_values):
        """
        Taking a dictionary of (sampled) nuisance parameter values params_values
        and return a log-likelihood.
        """
        """
        Compute the likelihood of the given cosmological parameters against the
        modeled angular diameter distances of the lenses.

        param theta: list of loat, folded cosmological parameters.
        param lenses: list of lens objects (currently either GLEELens or LenstronomyLens).
        param cosmology: string, keyword indicating the choice of cosmology to work with.
        """
        # JVR - not using this now
        # lp = log_prior(theta, cosmology)
        logpost = 0
        for lens in self.lenses:
            logpost += log_like_add(lens=lens, self=self)
        return logpost
    
####################################### TECHNICAL SUBROUTINES #############################################

def declare_lens_systems(self):
        ################################# CREATING LENS OBJECTS #######################################
        """B1608"""
        # B1608 using Ddt only, first analysis from Suyu+2010
        B1608 = GLEELens(name="B1608", longname="B1608 (Suyu+2010)", zlens=0.6304, zsource=1.394,
                            mu=7.0531390, sigma=0.2282395, lam=4000.0,
                            loglikelihood_type="sklogn_analytical"
                           )
        # B1608 using Dd only, analysis from Jee+2019
        B1608_Ddonly = GLEELens(name="B1608Dd", longname="B1608 Dd (Jee+2019)", zlens=0.6304, zsource=1.394,
                                   mu=7.0531390, sigma=0.2282395, lam=4000.0, mu_Dd = 6.79671, sigma_Dd=0.1836, lam_Dd = 334.2,
                            loglikelihood_type="sklogn_analytical_Ddonly"
                           )
        # B1608 using both Ddt and Dd only, analysis from Suyu+2010, Jee+2019, used in Wong+2019
        B1608_DdDdt = GLEELens(name="B1608DdDdt", longname="B1608 (Suyu+2010, Jee+2019)", zlens=0.6304, zsource=1.394,
                            mu=7.0531390, sigma=0.2282395, lam=4000.0, mu_Dd = 6.79671, sigma_Dd=0.1836, lam_Dd = 334.2,
                            loglikelihood_type="sklogn_analytical_Dd"
                           )
        """J1206"""
        ddt_vs_dd_1206s = pd.read_csv("./external_modules/data/h0licow_distance_chains/J1206_final.csv")
        J1206 = LenstronomyLens(name="J1206", longname="J1206 (Birrer+2019)", zlens=0.745, zsource=1.789, 
                                   ddt_vs_dd_samples=ddt_vs_dd_1206s,
                                   loglikelihood_type="kde_full", kde_type="scipy_gaussian",
                                  )
        """WFI2033"""
        #preprocess the 2033 chains...
        ddt_2033s_bic = pd.read_csv("./external_modules/data/h0licow_distance_chains/wfi2033_dt_bic.dat")

        # remove the Ddt that are above 8000, as it makes it hard to have a decent kde fit.
        cutweights = [w for w, ddt in zip(ddt_2033s_bic["weight"], ddt_2033s_bic["Dt"]) if 0 < ddt < 8000]
        cutddts = [ddt for ddt in ddt_2033s_bic["Dt"] if 0 < ddt < 8000] 
        ddt_2033s_bic = pd.DataFrame.from_dict(data={"ddt": cutddts, "weight": cutweights})

        # create the lens object
        WFI2033 = GLEELens(name="WFI2033", longname="WFI2033 (Rusu+2019)", zlens=0.6575, zsource=1.662,
                              loglikelihood_type="kde_hist", kde_kernel="gaussian", ddt_samples=ddt_2033s_bic['ddt'],
                                  weights=ddt_2033s_bic['weight'],
                                  bandwidth=20, nbins_hist=400
                             )
        """HE0435"""
        #Using only HST data, analysis from Wong+2017
        HE0435_HST = GLEELens(name="HE0435_HST", longname="HE0435-HST (Wong+2017)", zlens=0.4546, zsource=1.693,
                             mu=7.57930024e+00, sigma=1.03124167e-01, lam=6.53901645e+02,
                             loglikelihood_type="sklogn_analytical"
                            )
        #Using HST + AO data, analysis from Chen+2019
        ddt_0435s_AO_HST = pd.read_csv("./external_modules/data/h0licow_distance_chains/HE0435_Ddt_AO+HST.dat", delimiter=" ", skiprows=1, names=("ddt",))
        ddt_0435s_AO_HST["weight"] = np.ones(len(ddt_0435s_AO_HST["ddt"]))

        HE0435_AO_HST = GLEELens(name="HE0435_AO_HST", longname="HE0435 (Wong+2017, Chen+2019)", zlens=0.4546, zsource=1.693,
                            loglikelihood_type="kde_hist", kde_kernel="gaussian", ddt_samples=ddt_0435s_AO_HST['ddt'],
                                weights = ddt_0435s_AO_HST["weight"], 
                              bandwidth=20, nbins_hist=400
                            )
        """RXJ1131"""
        #Using only HST data, analysis from Suyu+2014
        RXJ1131_HST = GLEELens(name="RXJ1131_HST", longname="RXJ1131-HST (Suyu+2014)", zlens=0.295, zsource=0.654, 
                              mu=6.4682, sigma=0.20560, lam=1388.8, 
                              loglikelihood_type="sklogn_analytical"
                             )
        #Using HST + AO data, analysis from Chen+2019
        dd_vs_ddt_1131s_AO_HST = pd.read_csv("./external_modules/data/h0licow_distance_chains/RXJ1131_AO+HST_Dd_Ddt.dat", 
                                          delimiter=" ", skiprows=1, names=("dd", "ddt"))
        RXJ1131_AO_HST = GLEELens(name="RXJ1131_AO_HST", longname="RXJ1131 (Suyu+2014, Chen+2019)", 
                                      zlens=0.295, zsource=0.654,
                                      loglikelihood_type="kde_hist_2d", kde_kernel="gaussian",
                                      bandwidth=20, nbins_hist=80,
                                      ddt_samples=dd_vs_ddt_1131s_AO_HST["ddt"], dd_samples=dd_vs_ddt_1131s_AO_HST["dd"]
                             )
        """PG1115"""
        #Using HST + AO data, analysis from Chen+2019
        dd_vs_ddt_1115s= pd.read_csv("./external_modules/data/h0licow_distance_chains/PG1115_AO+HST_Dd_Ddt.dat", 
                                          delimiter=" ", skiprows=1, names=("dd", "ddt"))
        PG1115 = GLEELens(name="PG1115", longname="PG1115 (Chen+2019)", 
                                      zlens=0.311, zsource=1.722,
                                      loglikelihood_type="kde_hist_2d", kde_kernel="gaussian",
                                      bandwidth=20, nbins_hist=80,
                                      ddt_samples=dd_vs_ddt_1115s["ddt"], dd_samples=dd_vs_ddt_1115s["dd"]
                             )
        # Which lenses set do you want to work on :
        # The full combination from H0LICOW XIII paper
        self.lenses = [B1608_DdDdt, RXJ1131_AO_HST, HE0435_AO_HST, J1206, WFI2033, PG1115]