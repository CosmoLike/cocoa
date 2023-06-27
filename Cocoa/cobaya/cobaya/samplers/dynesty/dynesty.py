"""
.. module:: samplers.dynesty

Kunhao Zhong: A low-level implementation for the dynesty cobaya wrapper
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
import dynesty as dnt


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


class dynesty(Sampler):
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
        self.sampled_params_names = list(self.model.parameterization.sampled_params().keys())
        if self.n_derived>0:
            raise LoggedError(
                self.log, "Does Not support derived parameters YET")

        self.prior_bounds = self.model.prior.bounds(
            confidence_for_unbounded=self.confidence_for_unbounded)

        # Import additional modules for parallel computing if requested
        pool = None
        # Multithreading parallel using dynesty's own pool
        if self.parallel.get("kind") == "multiprocessing":
            self.mpi_info("Using dynesty's own pool")
            self.n_threads = self.parallel.get("args", dict(threads=1)).get("threads")

        self.rstate = np.random.default_rng(5647)

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
        # self.output.create_folder(self.base_dir) # Do not save to raw folder as PolyChord here
        self.mpi_info("Storing dynesty output to '%s'.", self.base_dir)


    def prior_transform(self, u):
        x = np.array(u) # copy u
        for i in range(len(x)):
            prior_min = self.prior_bounds[i][0]
            prior_max = self.prior_bounds[i][1]
            x[i] = u[i]*(prior_max-prior_min)+prior_min
        return x


    def run(self):
        """
        Prepares the likelihood and prior_transform function and calls ``Dynesty``'s ``run`` function.
        """
        self.mpi_info("Calling dynesty...")
        def loglikelihood(params_values):
            result = self.model.logposterior(params_values)
            loglikes = result.loglikes
            # should have extra argument if derived is present? 
            return np.squeeze(loglikes)

        sync_processes()
        if self.is_dynamical:
            self.mpi_info("Use Dynamical Nested Sampling")
            from dynesty import DynamicNestedSampler
            if self.parallel.get("kind") == "schwimmbad":
                #KZ: Not working as desired
                from schwimmbad import MPIPool
                #KZ: TODO
                pool = MPIPool()
                if not pool.is_master():
                    pool.wait()
                    sys.exit(0)
                sampler = DynamicNestedSampler(loglikelihood, self.prior_transform, self.nDims, bound='single', pool=pool)
                sampler.run_nested(nlive_init=self.nlive_init, nlive_batch=self.nlive_batch, use_stop=True, 
                        wt_kwargs={'pfrac': self.posterior_evidence_ratio}, print_progress=self.print_progress)
            elif self.parallel.get("kind") == "multiprocessing":
                self.mpi_info("Warning: Check MPI")
                # import dynesty.pool as dypool         
                # with dypool.Pool(self.n_threads, loglikelihood, self.prior_transform) as pool:
                #     # The important thing that we provide the loglikelihood/prior transform from
                #     # the pool
                #     sampler = dnt.DynamicNestedSampler(pool.loglike, pool.prior_transform, self.nDims, bound='single', pool=pool,
                #                                         rstate=self.rstate)
                #     sampler.run_nested(dlogz_init=0.05, nlive_init=self.nlive_init, nlive_batch=self.nlive_batch, use_stop=True, 
                #         wt_kwargs={'pfrac': self.posterior_evidence_ratio}, print_progress=self.print_progress)
                from concurrent.futures import ThreadPoolExecutor
                if self.parallel.get("args").get("threads") == -1:
                    num_cores = int(os.environ.get("SLURM_CPUS_PER_TASK", 1))
                else:
                    num_cores = self.parallel.get("args").get("threads")
                self.mpi_info("Running with %i cores", num_cores)
                with ThreadPoolExecutor(max_workers=num_cores) as executor:
                    sampler = dnt.DynamicNestedSampler(loglikelihood, self.prior_transform, self.nDims, bound='single', pool=executor,
                                                        rstate=self.rstate, queue_size=num_cores)
                    sampler.run_nested(dlogz_init=0.05, nlive_init=self.nlive_init, nlive_batch=self.nlive_batch, use_stop=True, 
                        wt_kwargs={'pfrac': self.posterior_evidence_ratio}, print_progress=self.print_progress)
            else:
                sampler = DynamicNestedSampler(loglikelihood, self.prior_transform, self.nDims, bound='single')
                sampler.run_nested(nlive_init=self.nlive_init, nlive_batch=self.nlive_batch, use_stop=True, 
                        wt_kwargs={'pfrac': self.posterior_evidence_ratio}, print_progress=self.print_progress)
        else:  
            self.mpi_info("Use Static Nested Sampling")
            if self.parallel.get("kind") == "multiprocessing":
                self.mpi_info("Warning: Check MPI")
                from concurrent.futures import ThreadPoolExecutor

                num_cores = int(os.environ.get("SLURM_CPUS_PER_TASK", 1))
                with ThreadPoolExecutor(max_workers=num_cores) as executor:
                    sampler = dnt.NestedSampler(loglikelihood, self.prior_transform, self.nDims, 
                                             nlive=self.nlive, sample='rslice',pool=executor,
                                                 rstate=self.rstate, queue_size=num_cores)
                    sampler.run_nested(dlogz=self.dlogz, print_progress=self.print_progress)

            else:
                sampler = dnt.NestedSampler(loglikelihood, self.prior_transform, self.nDims, nlive=self.nlive,
                                   rstate=self.rstate)
                sampler.run_nested(print_progress=self.print_progress)

        res = sampler.results
        # save equal weight samples; dynesty also output real samples with weights
        # logvol is log_prior
        self.save_raw(samples_equal=res.samples_equal(), logvol=res.logvol, logl=res.logl, logz=res.logz, logzerr=res.logzerr)
        self.save_sample(res, "1") # use cobaya SampleCollection
        self.mpi_info("dynesty finished")

    def save_raw(self, samples_equal, logvol, logl, logz, logzerr):
        if is_main_process():
            # move save samples to another function using cobaya sample collection
            # equal_weights  = np.ones(len(samples_equal))
            # minuslogpost   = -logvol-logl
            # chi2           = -2*logl
            # samples_cobaya = np.column_stack((equal_weights, minuslogpost, samples_equal, -logvol ,chi2))
            # headers = ["weight", "minuslogpost"] + self.sampled_params_names + ["minuslogprior","chi2"]
            # assert len(headers)==np.shape(samples_cobaya)[1]
            # headers = ['%10s' % h for h in headers] # format headers
            # header_str = ' '.join(headers)
            # data_format_str = ' '.join(['%10.7f']*samples_cobaya.shape[1]) # format the data
            # np.savetxt(self.base_dir + '.1.txt', samples_cobaya, header=header_str, comments='', fmt=data_format_str)
            np.savetxt(self.base_dir + '.logz.txt',logz)
            np.savetxt(self.base_dir + '.logzerr.txt',logzerr)
        return

    def save_sample(self, results, name):
        if is_main_process():
            collection = SampleCollection(self.model, self.output, name=str(name))
            sample_equal = results.samples_equal()
            logpriors    = results.logvol
            loglikes     = results.logl

            for i in range(len(results.samples_equal())):
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

    @classmethod
    def get_base_dir(cls, output):
        if output:
            return output.add_suffix(cls._base_dir_suffix, separator="")
        return os.path.join(gettempdir(), cls._base_dir_suffix)
