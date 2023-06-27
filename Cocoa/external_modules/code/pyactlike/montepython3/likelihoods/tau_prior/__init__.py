import os
from montepython.likelihood_class import Likelihood_prior


class tau_prior(Likelihood_prior):

    # initialisation of the class is done within the parent Likelihood_prior. For
    # this case, it does not differ, actually, from the __init__ method in
    # Likelihood class.
    def loglkl(self, cosmo, data):

        tau_reio = data.cosmo_arguments['tau_reio']
        loglkl = -0.5 * (tau_reio - self.mu) ** 2 / (self.sigma ** 2)
        return loglkl
