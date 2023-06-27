"""
.. module:: samplers.pocomc

Kunhao Zhong: A low-level implementation for the pocomc cobaya wrapper
TODO:   1. Does NOT support derived parameters; NOT blocking (fast/slow hierachy)
        2. @classmethod about installation
        3. _correct_unphysical_fraction() 
        4. Parallel: Not working it seems///schwimmbad (across nodes); ipyparallel missing; using dynesty's own pool only
"""
# Global
import os
import sys
import numpy as np
import logging
import inspect
from itertools import chain
from typing import Any, Callable, Optional
from tempfile import gettempdir
import re
import pocomc as pc


# Local
from cobaya.tools import read_dnumber, get_external_function, \
    find_with_regexp, NumberWithUnits, load_module, VersionCheckError
from cobaya.sampler import Sampler
from cobaya.mpi import is_main_process, share_mpi, sync_processes
from cobaya.collection import SampleCollection
from cobaya.log import LoggedError, get_logger
from cobaya.install import download_github_release, NotInstalledError
from cobaya.yaml import yaml_dump_file
from cobaya.conventions import derived_par_name_separator, packages_path_arg, Extension


class pocomc(Sampler):
    r"""
    TODO
    """

    _base_dir_suffix = ""

    # variables from yaml
    confidence_for_unbounded: float
    nlive: NumberWithUnits
    is_dynamial: bool

    def initialize(self):
        self.n_sampled = len(self.model.parameterization.sampled_params())
        self.n_derived = len(self.model.parameterization.derived_params())
        self.n_priors = len(self.model.prior)
        self.n_likes = len(self.model.likelihood)
        self.nDims = self.model.prior.d()
        if self.n_derived>0:
            raise LoggedError(
                self.log, "Does Not support derived parameters YET")

        self.prior_bounds = self.model.prior.bounds(
            confidence_for_unbounded=self.confidence_for_unbounded)

        self.prior_lower = self.prior_bounds[:,0]
        self.prior_upper = self.prior_bounds[:,1]

        # Import additional modules for parallel computing if requested
        self.pool = None

        # Prepare output folders and prefixes
        if self.output:
            self.file_root = self.output.prefix
            self.read_resume = self.output.is_resuming()
        else:
            output_prefix = share_mpi(hex(int(self._rng.random() * 16 ** 6))[2:]
                                      if is_main_process() else None)
            self.file_root = output_prefix
            # dummy output -- no resume!
            self.read_resume = False
        self.base_dir = self.get_base_dir(self.output)
        # self.output.create_folder(self.base_dir)
        self.mpi_info("Storing pocomc output to '%s'.", self.base_dir)

    # def log_prior(self, x):
    #     if np.any((x < self.prior_lower) | (x > self.upper_upper)):  # If any dimension is out of bounds, the log prior is -infinity
    #         return -np.inf
    #     else:
    #         return -const

    def run(self):
        """
        Prepares the likelihood and prior_transform function and calls ``Pocomc``'s ``run`` function.
        """
        self.mpi_info("Calling pocomc...")

        def log_likelihood(params_values):
            loglikes = np.zeros(len(params_values))
            #print("testing", np.shape(params_values))
            if params_values.ndim == 1:
                params_values = params_values.reshape(1,-1)
            for i in range(len(params_values)):
                #print("TESTING", i, params_values[i])
                result = self.model.logposterior(params_values[i])
                loglikes[i] = result.loglikes[0]
            return np.squeeze(loglikes)

        def log_prior(params_values):
            log_prior = np.zeros(len(params_values))
            if params_values.ndim == 1:
                params_values = params_values.reshape(1,-1)
            for i in range(len(params_values)):
                result = self.model.logposterior(params_values[i])
                log_prior[i] = result.logpriors[0]
            # should have extra argument if derived is present? 
            return np.squeeze(log_prior)

        # generate random prior_samples
        random_values = np.random.uniform(size=(self.n_walkers, self.prior_bounds.shape[0]))
        prior_samples = self.prior_lower + (self.prior_upper - self.prior_lower) * random_values

        if self.parallel.get("kind") == "multiprocessing":
            self.mpi_info("Warning: Check MPI")
            from concurrent.futures import ThreadPoolExecutor
            if self.parallel.get("args").get("threads") == -1:
                num_cores = int(os.environ.get("SLURM_CPUS_PER_TASK", 1))
            else:
                num_cores = self.parallel.get("args").get("threads")
            self.mpi_info("Running with %i cores", num_cores)
            #KZ: It seems without it it's also parallelized
            with ThreadPoolExecutor(max_workers=num_cores) as executor:
                pmc = pc.Sampler(
                    self.n_walkers,
                    self.nDims,
                    log_likelihood,
                    log_prior,
                    bounds=self.prior_bounds,
                    random_state=0,
                    pool = executor
                )
                pmc.run(prior_samples)
            results = pmc.results
            self.save_raw(samples=pmc.results['samples'], logz=pmc.results['logz'])
            self.save_sample(results, "1") # use cobaya SampleCollection
            self.mpi_info("pocomc finished")
        return

    def save_raw(self, samples, logz):
        if is_main_process():
            # How to get directly output directory?
            np.savetxt(self.base_dir + '.logz.txt',logz)
        return

    def save_sample(self, results, name):
        if is_main_process():
            collection = SampleCollection(self.model, self.output, name=str(name))
            sample_equal = results['samples'] # should be equal weighted?
            logpriors    = results['logprior']
            loglikes     = results['loglikelihood']

            for i in range(len(sample_equal)):
                collection.add(
                    sample_equal[i],
                    # derived=row[2 + self.n_sampled:2 + self.n_sampled + self.n_derived], # can't handel derived now
                    weight = 1, # equal weights
                    logpriors=[logpriors[i]],
                    loglikes=[loglikes[i]]
                    )
            # make sure that the points are written
            collection.out_update()
            return

    def products(self):
        """
        Auxiliary function to define what should be returned in a scripted call.

        Returns:
           The sample ``SampleCollection`` containing the sequentially
           discarded live points.
        """
        if is_main_process():
            products = {
                "sample": self.collection, "logZ": self.logZ, "logZstd": self.logZstd}
            return products
        else:
            return {}
    @property
    def raw_prefix(self):
        return os.path.join(
            self.pc_settings.base_dir, self.pc_settings.file_root)
    @classmethod
    def get_base_dir(cls, output):
        if output:
            return output.add_suffix(cls._base_dir_suffix, separator="")
        return os.path.join(gettempdir(), cls._base_dir_suffix)
