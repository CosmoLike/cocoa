"""
We start with a bunch of standard imports..
"""


import numpy as np
import pandas as pd
import math
# JVR - won't need emcee
# import emcee
# JVR - ending
import sys
from sklearn.neighbors import KernelDensity
import os

"""
Next, we import the KDELikelihood implemented in `lenstronomy` in order to properly sample lenstronomy's 2D likelihood (Dd versus Ddt, more on that below).

We also import cosmological models from `astropy`. These will be used to compute cosmological angular diameter distances from cosmological parameters in various cosmologies.
"""


from lenstronomy.Cosmo.kde_likelihood import KDELikelihood
from astropy.cosmology import FlatLambdaCDM, FlatwCDM, LambdaCDM, w0waCDM



"""
The output of the lens modeling are angular diameter distances; depending on the modeling code used, it is either the time-delay distance Ddt alone (using `GLEE`) or a joint time-delay distance Ddt versus observer-deflector angular diameter distance Dd (using `lenstronomy`).

In order to easily merge the output of these two softwares - and possibly combine with other softwares in the future - we create the classes GLEELens and LenstronomyLens that we initiate with the angular diameter distances posteriors predicted by the modeling codes and whose main goal is to evalute the likelihood of a chosen cosmological model with respect to the input posteriors of the lens.

    GLEE's time-delay distance posteriors are usually well fitted by a skewed log-normal likelihood. The analytical fit parameters are thus the input parameters of a GLEELens object. Alternatively, one can fit the angular diameter distance with a KDELikelihood if the skewed log-normal fit does not properly capture the distribution.

    Lenstronomy's time-delay distance and observer-deflector angular diameter distance joint posteriors are fitted with a KDELikelihood. They are the input parameters of a LenstronomyLens object.
"""

class StrongLensSystem(object):
    """
    This is a parent class, common to all lens modeling code outputs.
    It stores the "physical" parameters of the lens (name, redshifts, ...)
    """
    def __init__(self, name, zlens, zsource, longname=None):
        self.name = name
        self.zlens = zlens
        self.zsource = zsource
        self.longname = longname

    def __str__(self):
        return "%s\n\tzl = %f\n\tzs = %f" % (self.name, self.zlens, self.zsource)


class GLEELens(StrongLensSystem):
    """
    This class takes the output of GLEE (Ddt distribution) from which it evaluates the likelihood of a Ddt (in Mpc) predicted in a given cosmology.

    The default likelihood follows a skewed log-normal distribution. You can also opt for a normal distribution for testing purposes. In case no analytical form fits the Ddt distribution well, one can use a KDE log_likelihood instead - either fitted on the whole Ddt distribution (slow) or on a binned version of it (faster, no significant loss of precision).
    
    You can now also give mu,sigma and lambda parameter of the skewed likelihood for Dd. 
    """
    def __init__(self, name, zlens, zsource,
                 loglikelihood_type="normal_analytical",
                 mu=None, sigma=None, lam=None, explim=100.,
                 ddt_samples=None, dd_samples=None, weights = None, kde_kernel=None,
                 bandwidth=20, nbins_hist=200,
                 longname=None, mu_Dd=None, sigma_Dd=None, lam_Dd=None):

        StrongLensSystem.__init__(self, name=name, zlens=zlens, zsource=zsource, longname=longname)
        self.mu = mu
        self.sigma = sigma
        self.lam = lam
        self.explim = explim
        self.loglikelihood_type = loglikelihood_type
        self.ddt_samples = ddt_samples
        self.dd_samples = dd_samples
        if weights is None :
            if self.ddt_samples is not None :
                self.weights = np.ones(len(self.ddt_samples))
            else :
                self.weights = None
        else :
            self.weights = weights
        self.kde_kernel = kde_kernel
        self.bandwidth = bandwidth
        self.nbins_hist = nbins_hist
        self.mu_Dd = mu_Dd
        self.sigma_Dd = sigma_Dd
        self.lam_Dd = lam_Dd

        # do it only once at initialisation
        if loglikelihood_type == "hist_lin_interp":
            self.vals, self.bins = np.histogram(self.ddt_samples["ddt"], bins=self.nbins_hist, weights=self.weights)

        self.init_loglikelihood()

    def sklogn_analytical_likelihood(self, ddt):
        """
        Evaluates the likelihood of a time-delay distance ddt (in Mpc) against the model predictions, using a skewed log-normal distribution.
        """
        # skewed lognormal distribution with boundaries
        if (ddt <= self.lam) or ((-self.mu + math.log(ddt - self.lam)) ** 2 / (2. * self.sigma ** 2) > self.explim):
            return -np.inf
        else:
            llh = math.exp(-((-self.mu + math.log(ddt - self.lam)) ** 2 / (2. * self.sigma ** 2))) / (math.sqrt(2 * math.pi) * (ddt - self.lam) * self.sigma)

            if np.isnan(llh):
                return -np.inf
            else:
                return np.log(llh)

    def sklogn_analytical_likelihood_Dd(self, ddt, dd):
        """
        Evaluates the likelihood of a time-delay distance ddt and angular diameter distance Dd(in Mpc) against the model predictions, using a skewed log-normal distribution for both ddt and dd. The two distributions are asssumed independant and can be combined
        """
        # skewed lognormal distribution with boundaries
        if (ddt < self.lam) or ((-self.mu + math.log(ddt - self.lam)) ** 2 / (2. * self.sigma ** 2) > self.explim) or (dd < self.lam_Dd) or ((-self.mu_Dd + math.log(dd - self.lam_Dd)) ** 2 / (2. * self.sigma_Dd ** 2) > self.explim):
            return -np.inf
        else:
            llh = math.exp(-((-self.mu + math.log(ddt - self.lam)) ** 2 / (2. * self.sigma ** 2))) / (math.sqrt(2 * math.pi) * (ddt - self.lam) * self.sigma)

            llh_Dd = math.exp(-((-self.mu_Dd + math.log(dd - self.lam_Dd)) ** 2 / (2. * self.sigma_Dd ** 2))) / (math.sqrt(2 * math.pi) * (dd - self.lam_Dd) * self.sigma_Dd)

            if np.isnan(llh) or np.isnan(llh_Dd):
                return -np.inf
            else:
                return np.log(llh) + np.log(llh_Dd)

    def sklogn_analytical_likelihood_Ddonly(self,dd) :
        """
        Evaluates the likelihood of a angular diameter distance Dd(in Mpc) against the model predictions, using a skewed log-normal distribution for dd.
        """
        # skewed lognormal distribution with boundaries
        #if (dd < self.lam_Dd) or ((-self.mu_Dd + math.log(dd - self.lam_Dd)) ** 2 / (2. * self.sigma_Dd ** 2) > self.explim):
         #   return -np.inf
        #else:
        llh_Dd = math.exp(-((-self.mu_Dd + math.log(dd - self.lam_Dd)) ** 2 / (2. * self.sigma_Dd ** 2))) / (math.sqrt(2 * math.pi) * (dd - self.lam_Dd) * self.sigma_Dd)

        if np.isnan(llh_Dd):
            return -np.inf
        else:
            return np.log(llh_Dd)



    def normal_analytical_likelihood(self, ddt):
        """
        Evaluates the likelihood of a time-delay distance ddt (in Mpc) against the model predictions, using a normalised gaussian distribution.
        """
        # Normal distribution with boundaries

        if np.abs(ddt - self.mu) > 3*self.sigma:
            return -np.inf
        else:
            lh = math.exp(- (ddt - self.mu) **2 / (2. * self.sigma **2) ) / (math.sqrt(2 * math.pi) * self.sigma)
            if np.isnan(lh):
                return -np.inf
            else:
                return np.log(lh)


    def kdelikelihood_full(self, kde_kernel, bandwidth):
        """
        Evaluates the likelihood of a time-delay distance ddt (in Mpc) against the model predictions, using a loglikelihood sampled from a Kernel Density Estimator. the KDE is constructed using the full ddt samples.

        __ warning:: you should adjust bandwidth to the spacing of your samples chain!
        """
        kde = KernelDensity(kernel=kde_kernel, bandwidth=bandwidth).fit(self.ddt_samples, sample_weight=self.weights)
        return kde.score


    def kdelikelihood_hist(self, kde_kernel, bandwidth, nbins_hist):
        """
        Evaluates the likelihood of a time-delay distance ddt (in Mpc) against the model predictions, using a loglikelihood sampled from a Kernel Density Estimator. the KDE is constructed using a binned version of the full samples. Greatly improves speed at the cost of a (tiny) loss in precision

        __warning:: you should adjust bandwidth and nbins_hist to the spacing and size of your samples chain!

        """
        hist = np.histogram(self.ddt_samples, bins=nbins_hist, weights=self.weights)
        vals = hist[0]
        bins = [(h + hist[1][i+1])/2.0 for i, h in enumerate(hist[1][:-1])]

        # ignore potential zero weights, sklearn does not like them
        kde_bins = [(b,) for v, b in zip(vals, bins) if v>0]
        kde_weights = [v for v in vals if v>0]

        kde = KernelDensity(kernel=kde_kernel, bandwidth=bandwidth).fit(kde_bins, sample_weight=kde_weights)
        return kde.score


    def kdelikelihood_hist_2d(self, kde_kernel, bandwidth, nbins_hist):
        """
        Evaluates the likelihood of a angular diameter distance to the deflector Dd (in Mpc) versus its time-delay distance Ddt (in Mpc) against the model predictions, using a loglikelihood sampled from a Kernel Density Estimator. The KDE is constructed using a binned version of the full samples. Greatly improves speed at the cost of a (tiny) loss in precision

        __warning:: you should adjust bandwidth and nbins_hist to the spacing and size of your samples chain!

        __note:: nbins_hist refer to the number of bins per dimension. Hence, the final number of bins will be nbins_hist**2

        """
        hist, dd_edges, ddt_edges = np.histogram2d(x=self.dd_samples, y=self.ddt_samples, bins=nbins_hist, weights=self.weights)
        dd_vals = [(dd + dd_edges[i+1])/2.0 for i, dd in enumerate(dd_edges[:-1])]
        ddt_vals = [(ddt + ddt_edges[i+1])/2.0 for i, ddt in enumerate(ddt_edges[:-1])]

        # ugly but fast enough way to get the correct format for kde estimates
        dd_list, ddt_list, vals = [], [], []
        for idd, dd in enumerate(dd_vals):
            for iddt, ddt in enumerate(ddt_vals):
                dd_list.append(dd)
                ddt_list.append(ddt)
                vals.append(hist[idd, iddt])

        kde_bins = pd.DataFrame.from_dict({"dd": dd_list, "ddt": ddt_list})
        kde_weights = np.array(vals)

        # remove the zero weights values
        kde_bins = kde_bins[kde_weights > 0]
        kde_weights = kde_weights[kde_weights > 0]

        # fit the KDE
        kde = KernelDensity(kernel=kde_kernel, bandwidth=bandwidth).fit(kde_bins, sample_weight=kde_weights)
        return kde.score


    def hist_lin_interp_likelihood(self, ddt):
        """
        Evaluates the likelihood of a time-delay distance ddt (in Mpc) agains the model predictions, using linear interpolation from an histogram.

        __warning:: for testing purposes only - prefer kdelikelihood_hist, which gives similar results
        """
        if ddt <= self.bins[0] or ddt >= self.bins[-1]:
            return -np.inf
        else:
            indright = np.digitize(ddt, self.bins)
            #return np.log((self.vals[indright-1]+self.vals[indright])/2.0)
            return np.log(self.vals[indright-1])


    def init_loglikelihood(self):
        if self.loglikelihood_type == "sklogn_analytical":
            self.loglikelihood = self.sklogn_analytical_likelihood

        elif self.loglikelihood_type == "sklogn_analytical_Dd":
            self.loglikelihood = self.sklogn_analytical_likelihood_Dd

        elif self.loglikelihood_type == "sklogn_analytical_Ddonly":
            self.loglikelihood = self.sklogn_analytical_likelihood_Ddonly

        elif self.loglikelihood_type == "normal_analytical":
            self.loglikelihood = self.normal_analytical_likelihood

        elif self.loglikelihood_type == "kde_full":
            self.loglikelihood = self.kdelikelihood_full(kde_kernel=self.kde_kernel, bandwidth=self.bandwidth)

        elif self.loglikelihood_type == "kde_hist":
            self.loglikelihood = self.kdelikelihood_hist(kde_kernel=self.kde_kernel, bandwidth=self.bandwidth, nbins_hist=self.nbins_hist)

        elif self.loglikelihood_type == "kde_hist_2d":
            self.loglikelihood = self.kdelikelihood_hist_2d(kde_kernel=self.kde_kernel, bandwidth=self.bandwidth, nbins_hist=self.nbins_hist)

        elif self.loglikelihood_type == "hist_lin_interp":
            self.loglikelihood = self.hist_lin_interp_likelihood

        else:
            assert ValueError("unknown keyword: %s" % self.loglikelihood_type)
            # if you want to implement other likelihood estimators, do it here



class LenstronomyLens(StrongLensSystem):
    """
    This class takes the output of Lenstronomy (Dd versus Ddt distributions) from which it evaluates the likelihood of a Dd versus Ddt (in Mpc) predicted in a given cosmology.

    The default likelihood follows the KDE log-normal distribution implemented in Lenstronomy. You can change the type of kernel used. No other likelihoods have been implemented so far.
    """

    def __init__(self, name, zlens, zsource, ddt_vs_dd_samples, longname=None, loglikelihood_type="kde", kde_type="scipy_gaussian"):
        StrongLensSystem.__init__(self, name=name, zlens=zlens, zsource=zsource, longname=longname)

        self.ddt_vs_dd = ddt_vs_dd_samples     
        self.loglikelihood_type = loglikelihood_type
        self.kde_type = kde_type
        self.init_loglikelihood()

    def kdelikelihood(self, kde_type, bandwidth=20):
        """
        Evaluates the likelihood of a angular diameter distance to the deflector Dd (in Mpc) versus its time-delay distance Ddt (in Mpc) against the model predictions, using a loglikelihood sampled from a Kernel Density Estimator.
        """
        self.ddt = self.ddt_vs_dd["ddt"]
        self.dd = self.ddt_vs_dd["dd"]
        KDEl = KDELikelihood(self.dd.values, self.ddt.values, kde_type=kde_type, bandwidth=bandwidth)
        return KDEl.logLikelihood


    def init_loglikelihood(self):
        if self.loglikelihood_type == "kde_full":
            self.loglikelihood = self.kdelikelihood(kde_type=self.kde_type)
        else:
            assert ValueError("unknown keyword: %s" % self.loglikelihood_type)
            # if you want to implement other likelihood estimators, do it here


"""
The functions below evaluate the joint likelihood function of a set of cosmological parameters in a chosen cosmology against the modeled angular diameter distances of the lens systems.

The core of the process is the sample_params function, that works around the `emcee` MCMC sampler. Provided with uniform priors (log_prior), a list of lens objects and a suitable `astropy.cosmo` cosmology, it computes the angular diameter distances predicted by the cosmology (log_prob_ddt) and evaluate their log-likelihood against the angular diameter distances of the lens systems (log_like_add).

By sampling the cosmological parameters space, sample_params return chains of joint cosmological parameters that can be used to produce fancy posterior plots.
"""
        
                    
            
            



#     def log_like_add(lens, _self):
#     """
#     Computes the relevant angular diameter distance(s) of a given lens in a given cosmology,
#     and evaluate its/their joint likelihood against the same modeled distances of the lens.

#     param lens: either a GLEELens or LenstronomyLens instance.
#     param cosmo: an astropy cosmology object. 
#     """
#     # JVR: old implementation with AstroPy
#     # dd = cosmo.angular_diameter_distance(z=lens.zlens).value
#     # ds = cosmo.angular_diameter_distance(z=lens.zsource).value
#     # dds = cosmo.angular_diameter_distance_z1z2(z1=lens.zlens, z2=lens.zsource).value
#     # ddt = (1. + lens.zlens) * dd * ds / dds
#     # New implementation with CAMB
#     z_pair = (lens.zlens, lens.zsource)
#     dd = _self.provider.get_angular_diameter_distance(z=lens.zlens)
#     ds = _self.provider.get_angular_diameter_distance(z=lens.zsource)
#     dds = _self.provider.get_angular_diameter_distance_2(z_pair)
#     ddt = (1. + lens.zlens) * dd * ds / dds # time-delay distance

#     if isinstance(lens, GLEELens):
#         if lens.loglikelihood_type in ["kde_full", "kde_hist"]:
#             # because the newest sklearn kde want arrays, not numbers... 
#             return lens.loglikelihood(np.array(ddt).reshape(1, -1))
#         elif lens.loglikelihood_type in ["kde_hist_2d"]:
#             return lens.loglikelihood(np.array([dd, ddt]).reshape(1, -1))
#         elif lens.loglikelihood_type in ["sklogn_analytical_Dd"] :
#             return lens.loglikelihood(ddt, dd)
#         elif lens.loglikelihood_type in ["sklogn_analytical_Ddonly"] :
#             return lens.loglikelihood(dd)
#         else:
#             return lens.loglikelihood(ddt)

#     elif isinstance(lens, LenstronomyLens):
#         return lens.loglikelihood(dd, ddt)

#     else:
#         sys.exit("I don't know what to do with %s, unknown instance" % lens)

# JVR - Ported this function to main likelihood file h0licowlike.py
        
# def log_prob_ddt(theta, lenses, cosmology):
#     """
#     Compute the likelihood of the given cosmological parameters against the
#     modeled angular diameter distances of the lenses.

#     param theta: list of loat, folded cosmological parameters.
#     param lenses: list of lens objects (currently either GLEELens or LenstronomyLens).
#     param cosmology: string, keyword indicating the choice of cosmology to work with.
#     """

#     lp = log_prior(theta, cosmology)
#     if not np.isfinite(lp):
#         return -np.inf
#     else:
#         logprob = lp
#         if cosmology == "FLCDM":
#             h0, om = theta
#             cosmo = FlatLambdaCDM(H0=h0, Om0=om)
#         elif cosmology == "ULCDM":
#             h0 = theta[0]
#             cosmo = FlatLambdaCDM(H0=h0, Om0=0.3)
#         elif cosmology == "FwCDM":
#             h0, om, w = theta
#             cosmo = FlatwCDM(H0=h0, Om0=om, w0=w)
#         elif cosmology == "w0waCDM":
#             h0, om, w0, wa = theta
#             cosmo = w0waCDM(H0=h0, Om0=om, Ode0=1.0-om, w0=w0, wa=wa)
#         elif cosmology == "oLCDM":
#             h0, om, ok = theta
#             # assert we are not in a crazy cosmological situation that prevents computing the angular distance integral
#             if np.any([ok*(1.0+lens.zsource)**2 + om*(1.0+lens.zsource)**3 + (1.0-om-ok) <= 0 for lens in lenses]):
#                 return -np.inf
#             else:
#                 cosmo = LambdaCDM(H0=h0, Om0=om, Ode0=1.0-om-ok)
#         else:
#             raise ValueError("I don't know the cosmology %s" % cosmology)
#         for lens in lenses:
#             logprob += log_like_add(lens=lens, cosmo=cosmo)
#         return logprob

# JVR - From here, it's just a emcee sampler. We use cobaya, so won't need this part of the code, just logp.

# def sample_params(lenses, cosmology, nwalkers=32, nsamples=20000, save=True, filepath="temp.pkl", cluster = False):
#     """
#     High-level wrapper around the above functions. Explore the cosmological parameters space and
#     return their likelihood evaluated against the modeled angular diameter distances
#     of (multiple) lens system(s).

#     param lenses: list of lens objects (currently either GLEELens or LenstronomyLens).
#     param cosmology: string, keyword indicating the choice of cosmology to work with.
#     param nwalkers: int, number of emcee walkers used to sample the parameters space.
#     param nsamples: int, number of samples for an MCMC chain to converge. 
#         Make sure these are larger than the autocorrelation time!
#     param save: boolean, if True the combined, flattened chain is saved in filepath
#     param filepath: string, path of where the output chain is saved.
#     """


#     # Our starting point is a "decent" solution. You might want to check it does not impact the results, but unless you start from really crazy values, it shouldn't. We slightly randomize the starting point for each walker.
#     if cosmology == "FLCDM":
#         startpos = [72, 0.3]  # H0, Om
#     elif cosmology == "ULCDM":
#         startpos = [72] # H0
#     elif cosmology == "FwCDM":
#         startpos = [72, 0.3, -1]  # H0, Om, w
#     elif cosmology == "w0waCDM":
#         startpos = [72, 0.3, -1, 0] # H0, Om, w0, wa
#     elif cosmology == "oLCDM":
#         startpos = [72, 0.3, 0.0]  # H0, Om, Ok
#     else:
#         raise ValueError("I don't know the cosmology %s" % cosmology)

#     pos = startpos + 1e-4*np.random.randn(nwalkers, len(startpos))
#     nwalkers, ndim = pos.shape

#     sampler = emcee.EnsembleSampler(nwalkers, ndim, log_prob_ddt, args=([lenses, cosmology]))
#     if cluster :
#         sampler.run_mcmc(pos, nsamples)
#     else :
#         sampler.run_mcmc(pos, nsamples, progress =True)


#     tau = sampler.get_autocorr_time()
#     print("Autocorrelation time: ", tau)
#     flat_samples = sampler.get_chain(discard=1000, thin=4, flat=True)

#     if save:
#         pkl_file = open(filepath, 'wb')
#         pickle.dump(flat_samples, pkl_file, protocol=-1)
#         pkl_file.close()

#     return flat_samples


# def readpickle(filepath):
#     """
#     Small utility function to read pickle files written by the sample_params function above.

#     param filepath: path of the pickle file to load.
#     """
#     pkl_file = open(filepath, 'rb')
#     try:
#         obj = pickle.load(pkl_file, encoding='bytes')
#     except:
#         obj = pickle.load(pkl_file)
#     pkl_file.close()
#     return obj

# def confinterval(xs, weights=None, displayhist=False, conflevel=0.68, cumulative=False, bound='max'):
#     """
#     Hand-made weighted quantile/percentile computation tool, since there is no numpy equivalent

#     :param xs: list, distribution
#     :param weights: list, weights of the distribution
#     :param displayhist: bool. Show a nice plot
#     :param conflevel: float, quantile (between 0 and 1) of the conf interval
#     :param cumulative: bool. If set to True, return the "bound" (min or max) value of xs corresponding to the conflevel. For example, if bound=max, conflevel=0.68, xs=nnu, return the value such that nnu<=val at 68% confidence.
#     :param bound: 'max' or 'min', depending on the xs distribution. Valid only if cumulative is set to True
#     :return: return 3 floats: the median value, minimal value and maximal value matching your confidence level
    
#     __warning:: This function is not optimized, so it can be really slow if you apply it to very large data sets.
#                 However, for the samples we are dealing with here (around 50k), that's enough (fraction of seconds)
    
#     __note:: This function has been compared to numpy.percentiles for uniform weights 
#                 and give similar results (at 0.02%). For non-uniform weights, it is able to reproduce the
#                 results from the Planck collaboration with the same level of precision
    
#     """

#     hist, bin_edges = np.histogram(xs, weights=weights, bins=1000)
#     #meanval = np.average(xs, weights=weights)

#     if not cumulative:
#         frac = 0
#         for ind, val in enumerate(hist):
#             frac += val
#             if frac/sum(hist) >= 0.5:
#                 meanval = (bin_edges[ind]+bin_edges[ind-1])/2.0
#                 break

#         lowerhalf = 0
#         upperhalf = 0
#         for ind, val in enumerate(bin_edges[:-1]):
#             if val < meanval:
#                 lowerhalf += hist[ind]
#             else:
#                 upperhalf += hist[ind]


#         limfrac = (1.0 - conflevel) / 2.0
#         lowerfrac, upperfrac = 0, 0
#         for ind, val in enumerate(hist):
#             lowerfrac += val
#             if lowerfrac/sum(hist) > limfrac:
#                 bottomlimit = (bin_edges[ind]+bin_edges[ind-1])/2.0
#                 break

#         for ind, val in enumerate(hist):
#             upperfrac += val
#             if upperfrac/sum(hist) > 1.0 - limfrac:
#                 toplimit = (bin_edges[ind]+bin_edges[ind-1])/2.0
#                 break

#         return meanval, (meanval-bottomlimit), -1.0*(meanval-toplimit)

#     else:
#         frac = 0
#         for ind, val in enumerate(hist):
#             frac += val
#             if frac/sum(hist) >= conflevel:
#                 return (bin_edges[ind]+bin_edges[ind-1])/2.0
            
# def compute_percentile_from_dist(center_bins, value, percentile):
#     """
#     Function to return the percentile of a distribution providing the center of the bins and the number of sample in each bin
#     the value has to be normalized to one. 
    
#      :param center_bins: list, center of the bins of your distribution 
#      :param value: list, which have the same size as center_bins, 
#      :parma percentile : float between 0 and 1, desired percentile you want to obtain
#      :return: return 1 floats: the value of the percentile you asked for
#     """
    
#     if percentile > 1 or percentile <0 :
#         raise RuntimeError('Percentile should be between 0 and 1.')
#     if len(center_bins) != len(value) : 
#         raise RuntimeError('The length of the bins and the values are not the same.')

#     cumul = 0.0
#     ind = 0
#     for i in range(len(center_bins)):
#         cumul +=value[i]
#         if cumul >= percentile:
#             ind = i
#             break
#     return center_bins[ind]
