###################################
# MOCK CMB TYPE LIKELIHOOD
# --> mock planck, cmbpol, etc.
# Based on MontePython 3.3 Likelihood_mock_cmb
# and Cobaya's example of likelihood class
###################################
from cobaya.likelihood import Likelihood
from cobaya.log import LoggedError
import numpy as np
import os
import sys


def fast_det(A, N):
    """
    determinant working faster than np.linalg.det for 2x2, 1x1 and 3x3
    assumes matrix indices N x N to be the first ones
    """

    if N == 2:  # the default case
        return A[1, 1]*A[0, 0]-A[1, 0]*A[0, 1]
    if N == 1:
        return A[0, 0]
    if N == 3:
        return A[0, 0]*A[1, 1]*A[2, 2]+A[0, 2]*A[1, 0]*A[2, 1] +\
               A[2, 0]*A[0, 1]*A[1, 2]-A[0, 2]*A[1, 1]*A[2, 0] -\
               A[0, 0]*A[1, 2]*A[2, 1]-A[2, 2]*A[1, 0]*A[0, 1]
    return np.linalg.det(A.T).T


class MockCMBLikelihood(Likelihood):

    def initialize(self):
        """
        Set default settings, call other functions to load/set up data
        """

        self.data_directory = getattr(self, 'data_directory', None)

        if not self.data_directory:
            path = sys.modules[self.__module__].__file__
            self.data_directory = os.path.dirname(os.path.realpath(path))

        # - ignore B modes by default:
        self.Bmodes = getattr(self, 'Bmodes', False)
        # - do not use delensing by default:
        self.delensing = getattr(self, 'delensing', False)
        # - do not include lensing extraction by default:
        self.LensingExtraction = getattr(self, 'LensingExtraction', False)
        # - neglect TD correlation by default:
        self.neglect_TD = getattr(self, 'neglect_TD', True)
        # - use the lensed TT, TE, EE by default:
        self.unlensed_clTTTEEE = getattr(self, 'unlensed_clTTTEEE', False)
        # - do not exclude TTEE by default:
        self.ExcludeTTTEEE = getattr(self, 'ExcludeTTTEEE', False)
        if self.ExcludeTTTEEE and not self.LensingExtraction:
            raise LoggedError(self.log,
                              "Mock CMB likelihoods where TTTEEE is not used "
                              "have only been implemented for the deflection "
                              "spectrum (i.e. not for B-modes), but you do not"
                              "seem to have lensing extraction enabled")
        # added by Siavash Yasini
        self.OnlyTT = getattr(self, 'OnlyTT', False)
        if self.OnlyTT and self.ExcludeTTTEEE:
            raise LoggedError(self.log, "OnlyTT and ExcludeTTTEEE cannot be "
                                        "used simultaneously.")

        self.init_noise()

        if self.delensing:
            self.init_delensing()

        self.load_fid_values()

        if self.fid_values_exist:
            self.precalc_lkl()

        # Explicitly display the flags
        # to be sure that likelihood does what you expect:
        self.log.info("Initialised MockCMBLikelihood with following options:")
        self.log.info("unlensed_clTTTEEE is %s" % str(self.unlensed_clTTTEEE))
        self.log.info("Bmodes is %s" % str(self.Bmodes))
        self.log.info("delensing is %s" % str(self.delensing))
        self.log.info("LensingExtraction is %s" % str(self.LensingExtraction))
        self.log.info("neglect_TD is %s" % str(self.neglect_TD))
        self.log.info("ExcludeTTTEEE is %s" % str(self.ExcludeTTTEEE))
        self.log.info("OnlyTT is %s" % str(self.OnlyTT))

    def init_noise(self):
        """
        Set up noise spectrum
        """

        self.noise_from_file = getattr(self, 'noise_from_file', False)

        self.noise_T = np.zeros(self.l_max+1, 'float64')
        self.noise_P = np.zeros(self.l_max+1, 'float64')

        if self.noise_from_file:

            self.noise_file = getattr(self, 'noise_file', None)
            if not self.noise_file:
                raise LoggedError(self.log,
                                  "For reading noise from file,"
                                  "you must provide noise_file")

            if self.LensingExtraction:
                self.Nldd = np.zeros(self.l_max+1, 'float64')

            noise_filename = os.path.join(self.data_directory, self.noise_file)
            if os.path.exists(noise_filename):
                noise_content = np.loadtxt(noise_filename).T
                ll = noise_content[0].astype(int)
                self.noise_T[ll] = noise_content[1]
                self.noise_P[ll] = noise_content[2]
                if self.LensingExtraction:
                    try:
                        # read noise for C_l^dd = l(l+1) C_l^pp
                        self.Nldd[ll] = noise_content[3]/(ll*(ll+1)/2./np.pi)
                    except IndexError:
                        raise LoggedError(self.log,
                                          "For reading lensing noise from "
                                          "file, you must provide one more "
                                          "column")
            else:
                raise LoggedError(self.log, "Could not find file %s" %
                                            noise_filename)

        else:

            # convert arcmin to radians
            self.theta_fwhm *= np.array([np.pi/60/180])
            self.sigma_T *= np.array([np.pi/60/180])
            self.sigma_P *= np.array([np.pi/60/180])

            # compute noise in muK**2
            ll = np.arange(self.l_min, self.l_max+1)
            for channel in range(self.num_channels):
                self.noise_T[ll] += self.sigma_T[channel]**-2 *\
                    np.exp(
                        -ll*(ll+1)*self.theta_fwhm[channel]**2/8/np.log(2))
                self.noise_P[ll] += self.sigma_P[channel]**-2 *\
                    np.exp(
                        -ll*(ll+1)*self.theta_fwhm[channel]**2/8/np.log(2))
            self.noise_T[ll] = 1/self.noise_T[ll]
            self.noise_P[ll] = 1/self.noise_P[ll]

        # trick to remove any information from polarisation for l<30
        self.no_small_l_pol = getattr(self, 'no_small_l_pol', False)

        if self.no_small_l_pol:
            ll = np.arange(self.l_min, 30)
            # plug a noise level of 100 muK**2
            # equivalent to no detection at all of polarisation
            self.noise_P[ll] = 100.

        # trick to remove any information from temperature above l_max_TT
        self.l_max_TT = getattr(self, 'l_max_TT', False)

        if self.l_max_TT:
            ll = np.arange(self.l_max_TT+1, self.l_max+1)
            # plug a noise level of 100 muK**2
            # equivalent to no detection at all of temperature
            self.noise_T[ll] = 100.

    def init_delensing(self):
        """
        Delensing noise: implemented by  S. Clesse
        """

        self.delensing_file = getattr(self, 'delensing_file', None)
        if not self.delensing_file:
            raise LoggedError(self.log,
                              "For delensing, you must provide delensing_file")

        self.noise_delensing = np.zeros(self.l_max+1)
        delensing_filename = os.path.join(self.data_directory,
                                          self.delensing_file)
        if os.path.exists(delensing_filename):
            delensing_cont = np.loadtxt(delensing_filename).T
            ll = delensing_cont[0].astype(int)
            self.noise_delensing[ll] = delensing_cont[2]/(ll*(ll+1)/2./np.pi)
            # change 3 to 4 in the above line for CMBxCIB delensing

        else:
            raise LoggedError(self.log, "Could not find file %s" %
                              delensing_filename)

    def load_fid_values(self):
        """
        Read fiducial power spectrum for TT, EE, TE,
        [eventually BB or phi-phi, phi-T]
        """

        # default:
        if not self.ExcludeTTTEEE:
            numCls = 3
        # default 0 if excluding TT EE
        else:
            numCls = 0

        # deal with BB:
        if self.Bmodes:
            self.index_B = numCls
            numCls += 1

        # deal with pp, pT (p = CMB lensing potential):
        if self.LensingExtraction:
            self.index_pp = numCls
            numCls += 1
            if not self.ExcludeTTTEEE:
                self.index_tp = numCls
                numCls += 1

            if not self.noise_from_file:
                # provide a file containing NlDD (noise for the extracted
                # deflection field spectrum) This option is temporary
                # because at some point this module will compute NlDD
                # itself, when logging the fiducial model spectrum.
                self.temporary_Nldd_file = getattr(self, 'temporary_Nldd_file',
                                                   None)
                if not self.temporary_Nldd_file:
                    raise LoggedError(self.log, "For lensing extraction, you "
                                      "must provide a temporary_Nldd_file")

                # read the NlDD file
                self.Nldd = np.zeros(self.l_max+1, 'float64')

                temporary_Nldd_fname = os.path.join(self.data_directory,
                                                    self.temporary_Nldd_file)
                if os.path.exists(temporary_Nldd_fname):
                    temporary_Nldd_content = np.loadtxt(temporary_Nldd_fname).T
                    ll = temporary_Nldd_content[0].astype(int)
                    # this line assumes that Nldd is stored in the 4th column
                    # (can be customised)
                    self.Nldd[ll] = temporary_Nldd_content[3]/(ll*(ll+1.) /
                                                               2./np.pi)
                else:
                    raise LoggedError(self.log, "Could not find file %s" %
                                      temporary_Nldd_fname)

        # If the file exists, initialize the fiducial values
        self.Cl_fid = np.zeros((numCls, self.l_max+1), 'float64')
        self.fid_values_exist = False
        fiducial_filename = os.path.join(
                self.data_directory, self.fiducial_file)
        if os.path.exists(fiducial_filename):
            fiducial_content = np.loadtxt(fiducial_filename).T
            ll = fiducial_content[0].astype(int)
            try:
                self.Cl_fid[:, ll] = fiducial_content[1:]
                self.fid_values_exist = True
            except IndexError:
                self.log.warning("Fiducial model file has wrong number of rows"
                                 " or columns, so not loaded")
        else:
            self.log.warning("Fiducial model not loaded")

        # Else the file should be created in the create_fid_values() function.

    def precalc_lkl(self):
        """
        Calculate things needed for likelihood:
        count number of modes, set up Cov_obs and det_obs
        """

        # count number of modes.
        # number of modes is different from number of spectra
        # modes = T,E,[B],[D=deflection]
        # spectra = TT,EE,TE,[BB],[DD,TD]
        # default:
        if not self.ExcludeTTTEEE:
            if self.OnlyTT:
                self.num_modes = 1
            else:
                self.num_modes = 2
        # default 0 if excluding TT EE
        else:
            self.num_modes = 0
        # add B mode:
        if self.Bmodes:
            self.num_modes += 1
        # add D mode:
        if self.LensingExtraction:
            self.num_modes += 1

        # set up likelihood information depending only on fiducials
        self.set_cov_obs()

    def get_requirements(self):
        """
        here we need C_L^{...} to l_max
        but which Cl is a question
        follows the logics of load/create_fid_values, logp and precalc_lkl
        """
        # inner tuple of requirements
        cl_req = (dict(), dict())
        # 0 (False) is lensed Cl and 1 (True) is unlensed
        keys = ("Cl", "unlensed_Cl")
        # fill the requirements
        if not self.ExcludeTTTEEE:
            if self.OnlyTT:
                cl_req[self.unlensed_clTTTEEE]['tt'] = self.l_max
            else:
                cl_req[self.unlensed_clTTTEEE].update({'tt': self.l_max,
                                                       'te': self.l_max,
                                                       'ee': self.l_max})
        if self.Bmodes:
            cl_req[self.delensing]['bb'] = self.l_max
        if self.LensingExtraction:
            cl_req[self.unlensed_clTTTEEE]['pp'] = self.l_max
            if not self.ExcludeTTTEEE:
                cl_req[self.unlensed_clTTTEEE]['tp'] = self.l_max
        # leave only not empty requirements for output
        cl_req_out = dict()
        for i, key in enumerate(keys):
            if len(cl_req[i]) > 0:
                cl_req_out[key] = cl_req[i]
        return cl_req_out

    def logp(self, **params_values):
        """
        Taking a dictionary of (sampled) nuisance parameter values
        and return a log-likelihood.

        here calculate chi^2  using cl's without any nuisance pars
        """
        # if fiducial values exist
        if self.fid_values_exist:
            if self.unlensed_clTTTEEE:
                # get unlensed Cl's from the cosmological code in muK**2
                cl = self.provider.get_unlensed_Cl(units='muK2')
                # exception: for non-delensed B modes
                # we need the lensed BB spectrum
                # (this case is usually not useful/relevant)
                if self.Bmodes and (not self.delensing):
                    cl_lensed = self.provider.get_Cl(units='muK2')
                    cl['bb'] = cl_lensed['bb']
            else:
                # get lensed Cl's from the cosmological code in muK**2
                # without l*(l+1)/(2*pi) factor (default)
                cl = self.provider.get_Cl(units='muK2')
                # exception: for delensed B modes we need the unlensed spectrum
                if self.Bmodes and self.delensing:
                    cl_unlensed = self.provider.get_unlensed_Cl(units='muK2')
                    cl['bb'] = cl_unlensed['bb']

            # get likelihood
            return self.compute_lkl(cl, params_values)

        # otherwise warn and return -inf
        self.log.warning("Fiducial model not present, set likelihood=0")
        return -np.inf

    def compute_lkl(self, cl, params_values):
        """
        compute likelihood
        """

        ll = np.arange(self.l_min, self.l_max+1)

        Cov_the = self.get_cov_the(ll, cl)
        # Cov_obs already computed in self.Cov_obs

        # get determinant of observational and theoretical covariance matrices
        # det_obs already computed in self.det_obs
        det_the = fast_det(Cov_the, self.num_modes)
        det_mix = np.zeros(self.l_max-self.l_min+1)

        # get determinant of mixed matrix (= sum of N theoretical
        # matrices with, in each of them, the nth column replaced
        # by that of the observational matrix)
        # here it's actually row instead of column, thanks to symmetry
        for i in range(self.num_modes):
            Cov_mix = np.copy(Cov_the)
            Cov_mix[i] = self.Cov_obs[i]
            det_mix += fast_det(Cov_mix, self.num_modes)

        chi2 = np.sum((2.*ll+1.)*self.f_sky *
                      (det_mix/det_the + np.log(det_the/self.det_obs) -
                       self.num_modes))

        return -chi2/2

    def get_cov_the(self, ll, cl):
        """
        Fill the theoretical covariance matrix
        """

        if self.Bmodes and self.LensingExtraction:
            raise LoggedError(self.log,
                              "We have implemented a version of the likelihood"
                              "with B modes, a version with lensing extraction"
                              ", but not yet a version with both at the same "
                              "time. You can implement it.")

        # case with B modes:
        elif self.Bmodes:
            # next 5 lines added by S. Clesse for delensing
            if self.delensing:
                Cov_the = np.array([
                    [cl['tt'][ll]+self.noise_T[ll], cl['te'][ll], 0],
                    [cl['te'][ll], cl['ee'][ll]+self.noise_P[ll], 0],
                    [0, 0, cl['bb'][ll]+self.noise_P[ll] +
                     self.noise_delensing[ll]]])
            else:
                Cov_the = np.array([
                    [cl['tt'][ll]+self.noise_T[ll], cl['te'][ll], 0],
                    [cl['te'][ll], cl['ee'][ll]+self.noise_P[ll], 0],
                    [0, 0, cl['bb'][ll]+self.noise_P[ll]]])

        # case with lensing
        # note that the likelihood is based on ClDD (deflection spectrum)
        # rather than Clpp (lensing potential spectrum)
        # But the Bolztmann code input is Clpp
        # So we make the conversion using ClDD = l*(l+1.)*Clpp
        # So we make the conversion using ClTD = sqrt(l*(l+1.))*Cltp

        # just DD, i.e. no TT or EE.
        elif self.LensingExtraction and self.ExcludeTTTEEE:
            cldd = ll*(ll+1.)*cl['pp'][ll]
            Cov_the = np.array([[cldd+self.Nldd[ll]]])

        # Usual TTTEEE plus DD and TD
        elif self.LensingExtraction:
            cldd = ll*(ll+1.)*cl['pp'][ll]
            if self.neglect_TD:
                cltd = 0.
            else:
                cltd = np.sqrt(ll*(ll+1.))*cl['tp'][ll]

            Cov_the = np.array([
                [cl['tt'][ll]+self.noise_T[ll], cl['te'][ll], cltd],
                [cl['te'][ll], cl['ee'][ll]+self.noise_P[ll], 0],
                [cltd, 0, cldd+self.Nldd[ll]]])

        # case with TT only (Added by Siavash Yasini)
        elif self.OnlyTT:
            Cov_the = np.array([[cl['tt'][ll]+self.noise_T[ll]]])

        # case without B modes nor lensing:
        else:
            Cov_the = np.array([[cl['tt'][ll]+self.noise_T[ll], cl['te'][ll]],
                                [cl['te'][ll], cl['ee'][ll]+self.noise_P[ll]]])

        return Cov_the

    def set_cov_obs(self):
        """
        Fill the observational covariance matrix and compute its determinant
        """

        ll = np.arange(self.l_min, self.l_max+1)

        if self.Bmodes and self.LensingExtraction:
            raise LoggedError(self.log,
                              "We have implemented a version of the likelihood"
                              "with B modes, a version with lensing extraction"
                              ", but not yet a version with both at the same "
                              "time. You can implement it.")

        # case with B modes:
        elif self.Bmodes:
            Cov_obs = np.array([
                [self.Cl_fid[0, ll], self.Cl_fid[2, ll], 0],
                [self.Cl_fid[2, ll], self.Cl_fid[1, ll], 0],
                [0, 0, self.Cl_fid[3, ll]]])

        # case with lensing
        # note that the likelihood is based on ClDD (deflection spectrum)
        # rather than Clpp (lensing potential spectrum)
        # But the Bolztmann code input is Clpp
        # So we make the conversion using ClDD = l*(l+1.)*Clpp
        # So we make the conversion using ClTD = sqrt(l*(l+1.))*Cltp

        # just DD, i.e. no TT or EE.
        elif self.LensingExtraction and self.ExcludeTTTEEE:
            cldd_fid = self.Cl_fid[self.index_pp, ll]
            Cov_obs = np.array([[cldd_fid]])

        # Usual TTTEEE plus DD and TD
        elif self.LensingExtraction:
            cldd_fid = self.Cl_fid[self.index_pp, ll]
            if self.neglect_TD:
                cltd_fid = 0.
            else:
                cltd_fid = self.Cl_fid[self.index_tp, ll]

            Cov_obs = np.array([
                [self.Cl_fid[0, ll], self.Cl_fid[2, ll], cltd_fid],
                [self.Cl_fid[2, ll], self.Cl_fid[1, ll], 0],
                [cltd_fid, 0, cldd_fid]])

        # case with TT only (Added by Siavash Yasini)
        elif self.OnlyTT:
            Cov_obs = np.array([[self.Cl_fid[0, ll]]])

        # case without B modes nor lensing:
        else:
            Cov_obs = np.array([[self.Cl_fid[0, ll], self.Cl_fid[2, ll]],
                                [self.Cl_fid[2, ll], self.Cl_fid[1, ll]]])

        # save as attribute
        self.Cov_obs = Cov_obs

        # compute determinant
        self.det_obs = fast_det(Cov_obs, self.num_modes)

    def create_fid_values(self, cl, params, override=False):
        # Write fiducial model spectra if needed
        # params should be the dictionary of cosmological parameters
        # corresponding to the cl provided
        fid_filename = os.path.join(self.data_directory,
                                    self.fiducial_file)
        if not self.fid_values_exist or override:
            # header string with parameter values
            header_str = 'Fiducial parameters: '
            for key, value in params.items():
                value = str(value)
                header_str += '%s = %s, ' % (key, value)
            header_str = header_str[:-2]
            # output arrays
            ll = np.arange(self.l_min, self.l_max+1)
            out_data = [ll]
            if not self.ExcludeTTTEEE:
                out_data.extend((cl['tt'][ll]+self.noise_T[ll],
                                cl['ee'][ll]+self.noise_P[ll],
                                cl['te'][ll]))
            if self.Bmodes:
                # next three lines added by S. Clesse for delensing
                if self.delensing:
                    out_data.append(cl['bb'][ll]+self.noise_P[ll] +
                                    self.noise_delensing[ll])
                else:
                    out_data.append(cl['bb'][ll]+self.noise_P[ll])
            if self.LensingExtraction:
                # we want to store clDD = l(l+1) clpp
                out_data.append(ll*(ll+1.)*cl['pp'][ll] + self.Nldd[ll])
                # and ClTD = sqrt(l(l+1)) Cltp
                if not self.ExcludeTTTEEE:
                    out_data.append(np.sqrt(ll*(ll+1.))*cl['tp'][ll])
            out_data = np.array(out_data).T
            # write data
            np.savetxt(fid_filename, out_data, "%.8g", header=header_str)
            self.log.info(
                "Writing fiducial model in %s" % fid_filename)
            # (re-)initialize fiducial power spectra
            # and prepare likelihood calculation
            self.load_fid_values()
            if self.fid_values_exist:
                self.precalc_lkl()
            else:
                self.log.warning("Fiducial model not read")
            return True
        self.log.warning(
            "Fiducial model in %s already exists" % fid_filename)
        return False
