# this file, along with actpols2_like_py.data, allows you
# to use this likelihood for Monte Python.
#

import os
import numpy as np
from montepython.likelihood_class import Likelihood

import pyactlike  # our likelihood


class ACTPol_lite_DR4(Likelihood):

    # initialization routine
    def __init__(self, path, data, command_line):

        Likelihood.__init__(self, path, data, command_line)

        self.need_cosmo_arguments(
            data,
            {
                "lensing": "yes",
                "output": "tCl lCl pCl",
                "l_max_scalars": 6000,
                "modes": "s",
            },
        )

        self.need_update = True
        self.use_nuisance = ["yp2"]
        self.nuisance = ["yp2"]

        self.act = pyactlike.ACTPowerSpectrumData()

        # \ell values 2, 3, ... 6000
        self.xx = np.array(range(2, 6001))

    # compute likelihood

    def loglkl(self, cosmo, data):

        # print "STARTING LIKELIHOOD------------ ", data.cosmo_arguments

        lkl = 0.0
        try:
            # call CLASS
            cl = self.get_cl(cosmo, 6000)

            # we follow the convention of operating with (l(l+1)/2pi) * C_l
            ee = cl["ee"][2:]
            te = cl["te"][2:]
            tt = cl["tt"][2:]
            tt = (self.xx) * (self.xx + 1) * tt / (2 * np.pi)
            te = (self.xx) * (self.xx + 1) * te / (2 * np.pi)
            ee = (self.xx) * (self.xx + 1) * ee / (2 * np.pi)

            yp = data.mcmc_parameters["yp2"]["current"]
            lkl = self.act.loglike(tt, te, ee, yp)

        except:
            lkl = -np.inf

        return lkl
