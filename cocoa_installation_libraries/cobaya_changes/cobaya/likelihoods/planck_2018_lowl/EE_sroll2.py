import os
import numpy as np
from cobaya.likelihood import Likelihood


class EE_sroll2(Likelihood):
    type = "CMB"
    aliases = ["lowE"]

    _lmin = 2
    _lmax = 29
    _nstepsEE = 3000
    _stepEE = 0.0001
    _table_file_name = 'sroll2_prob_table.txt'

    @classmethod
    def get_bibtex(cls):
        return cls.get_associated_file_content('.bibtex')
    
    def initialize(self):
        self.path = os.getenv("ROOTDIR")
        self.path = self.path + "/external_modules/data/planck/"
        print(self.path)
        self.probEE = np.loadtxt(os.path.join(self.path, self._table_file_name))

    def get_can_support_params(self):
        return ['A_planck']

    def get_requirements(self):
        return {'Cl': {'ee': self._lmax}}

    def log_likelihood(self, cls_EE, calib=1):
        r"""
        Calculate log likelihood from CMB EE spectrum by using likelihood table

        :param cls_EE: L(L+1)C_L/2pi zero-based array in muK^2 units
        :param calib: optional calibration parameter
        :return: log likelihood
        """
        EE_index = (cls_EE[self._lmin:self._lmax + 1]
                    / (calib ** 2 * self._stepEE)).astype(int)
        try:
            return np.take_along_axis(self.probEE, EE_index[np.newaxis, :], 0).sum()
        except IndexError:
            self.log.warning("low EE multipole out of range, rejecting point")
            return -np.inf

    def logp(self, **params_values):
        cls = self.provider.get_Cl(ell_factor=True)['ee']
        return self.log_likelihood(cls, params_values.get('A_planck', 1))