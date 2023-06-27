"""
.. module:: samplers.emcee

:Synopsis: Goodman & Weare's Affine Invariant Markov chain Monte Carlo Ensemble sampler
:Author: Daniel Foreman-Mackey, David W. Hogg, Dustin Lang, Jonathan Goodman (wrapped for cobaya by Xavier Garrido)
"""
import os
import re
import sys
from typing import Any, Union

import numpy as np
from cobaya.collection import SampleCollection
from cobaya.conventions import Extension
from cobaya.install import NotInstalledError, pip_install
from cobaya.log import LoggedError, get_logger, is_debug
from cobaya.sampler import Sampler
from cobaya.tools import NumberWithUnits, VersionCheckError, load_module

Extension.hdf5 = ".h5"


class EMCEE(Sampler):
    r"Goodman & Weare's Affine Invariant Markov chain Monte Carlo Ensemble sample \cite{emcee2013}"

    # Name of the emcee repo/folder and version to download
    _min_emcee_version = "v3.1.1"
    _emcee_repo_name = os.environ.get("EMCEE_REPO_NAME", "dfm/emcee")
    _emcee_repo_version = os.environ.get("EMCEE_REPO_VERSION", _min_emcee_version)

    file_base_name = "emcee"

    # instance variables from yaml
    path: str
    nwalkers: Union[int, NumberWithUnits]
    nsteps: Union[int, NumberWithUnits]
    thin_by: int
    moves: dict
    progress: bool
    output_format: str
    emcee_module: Any

    def set_instance_defaults(self):
        super().set_instance_defaults()
        self.parallel = dict(kind="")
        self.ndim = 0
        self.nblobs = 0

    def initialize(self):
        # Allow global import if no direct path specification
        allow_global = not self.path
        if not self.path and self.packages_path:
            self.path = self.get_path(self.packages_path)
        self.emcee_module = self.is_installed(
            path=self.path, allow_global=allow_global, check=False
        )
        if not self.emcee_module:
            raise NotInstalledError(self.log, "Could not find emcee. Check error message above.")
        self.mpi_info("Loaded 'emcee' version %s", self.get_version())

        self.ndim = self.model.prior.d()
        if self.ndim == 0:
            raise LoggedError(self.log, "No parameters being varied for sampler")
        self.log.debug("Number of parameters to sample: %i", self.ndim)

        for var in ["nwalkers", "nsteps"]:
            if not isinstance(getattr(self, var), int):
                nwu = NumberWithUnits(getattr(self, var), "d", dtype=int, scale=self.ndim)
                setattr(self, var, nwu.value)
            self.log.debug("Number of %s: %i", var, getattr(self, var))

        # Import additional modules for parallel computing if requested
        pool = None
        # Multithreading parallel computing via python-native multiprocessing module
        if self.parallel.get("kind") == "multiprocessing":
            from multiprocessing import Pool

            pool = Pool(
                self.parallel.get("args", dict(threads=1)).get("threads")
            )  # number of threads chosen by user
        # MPI parallel computing via external schwimmbad module
        elif self.parallel.get("kind") == "mpi":
            from schwimmbad import MPIPool

            pool = MPIPool()
            if not pool.is_master():  # Necessary bit for MPI
                pool.wait()
                sys.exit(0)
        else:
            self.log.warning("No parallel processing will be used!")

        # Prepare some inputs for the MCMC
        blobs_dtype = [("logprior", float)]
        blobs_dtype += [(f"loglike_{name}", float) for name in self.model.likelihood]
        blobs_dtype += [
            (f"{derived}", float) for derived in self.model.parameterization.derived_params().keys()
        ]
        self.nblobs = len(blobs_dtype)

        # Initialize output file
        backend = None
        if self.output:
            if self.output_format == "hdf5":
                backend = self.emcee_module.backends.HDFBackend(
                    os.path.join(self.output.folder, self.output.prefix + Extension.hdf5)
                )
                if not self.output.is_resuming():
                    backend.reset(self.nwalkers, self.ndim)
                else:
                    self.log.info("Loading %i previous steps...", backend.iteration)
            elif self.output_format == "txt":
                # One collection per walker
                self.collection = [
                    SampleCollection(
                        self.model,
                        self.output,
                        name=str(i),
                        resuming=self.output.is_resuming(),
                    )
                    for i in range(1, self.nwalkers + 1)
                ]
            else:
                raise LoggedError(
                    self.log,
                    f"Unkown '{self.output_format}' output format! "
                    "Must be either 'hdf5' or 'txt'.",
                )

        # Initialize the move schedule
        moves = getattr(self.emcee_module.moves, self.moves.get("kind", "StretchMove"))
        moves = moves(**self.moves.get("args", {}))
        self.sampler = self.emcee_module.EnsembleSampler(
            self.nwalkers,
            self.ndim,
            self._wrapped_logposterior,
            moves=moves,
            pool=pool,
            backend=backend,
            blobs_dtype=blobs_dtype,
        )
        self.mpi_info("Initialized!")

    def run(self):

        kwargs = dict(
            initial_state=self.get_initial_state(),
            iterations=self.nsteps,
            thin_by=self.thin_by,
            progress=self.progress,
        )

        for istep, result in enumerate(self.sampler.sample(**kwargs)):
            if istep == 0:
                nfinite = np.isfinite(result.log_prob).sum()
                if nfinite < (0.5 * self.nwalkers):
                    self.log.warning(
                        "Your chain will take time to converge: "
                        "only %i%% of your walkers are starting "
                        "at a finite value of the posterior. "
                        "Please check if your starting positions are correct, and/or use "
                        "debug mode to check your likelihoods.",
                        nfinite * 100 / self.nwalkers,
                    )

            if self.output_format == "txt":
                blobs_arr = result.blobs.view(dtype=np.float64).reshape(self.nwalkers, -1)
                nlkl = len(self.model.likelihood)
                for iwalk in range(self.nwalkers):
                    blobs = blobs_arr[iwalk]
                    self.collection[iwalk].add(
                        values=result.coords[iwalk],
                        logpost=result.log_prob[iwalk],
                        logpriors=blobs[:1],
                        loglikes=blobs[1 : 1 + nlkl],
                        derived=blobs[1 + nlkl :],
                    )
                    self.collection[iwalk].out_update()

        self.mpi_info("Reached maximum number of steps allowed (%s). Stopping.", self.nsteps)

    def get_initial_state(self):
        """
        Get initial/starting points either from paramters PDF or by loading the initial state from
        previous runs
        """
        if self.output_format == "hdf5" and self.output.is_resuming():
            return self.sampler._previous_state

        initial_state = np.empty((self.nwalkers, self.ndim))
        for i in range(self.nwalkers):
            if self.output_format == "txt" and self.output.is_resuming():
                last = len(self.collection[i]) - 1
                initial_point = (
                    self.collection[i][self.collection[i].sampled_params].iloc[last]
                ).to_numpy(dtype=np.float64, copy=True)
            else:
                initial_point, _results = self.model.get_valid_point(
                    max_tries=100 * self.ndim, random_state=self._rng
                )
                # initial_point = self.model.prior.reference()
            initial_state[i] = initial_point

        return initial_state

    def _wrapped_logposterior(self, param_values):
        results = self.model.logposterior(param_values)
        if results.logpost == -np.inf:
            results = [-np.inf] * (1 + self.nblobs)
        else:
            results = (
                [results.logpost] + results.logpriors + results.loglikes.tolist() + results.derived
            )
        return tuple(results)

    # Class methods
    @classmethod
    def output_files_regexps(cls, output, info=None, minimal=False):
        regexps = [output.collection_regexp(name=None)]
        # if minimal:
        #     return [(r, None) for r in regexps]
        regexps += [
            re.compile(output.prefix_regexp_str + re.escape(ext.lstrip(".")) + "$")
            for ext in [Extension.hdf5]
        ]
        return [(r, None) for r in regexps]

    def get_version(self):
        return getattr(self.emcee_module, "__version__", None)

    # Installation routines

    @classmethod
    def get_path(cls, path):
        return os.path.realpath(os.path.join(path, "code", cls.__name__))

    @classmethod
    def is_installed(cls, **kwargs):
        if not kwargs.get("code", True):
            return True
        log = get_logger(cls.__name__)
        check = kwargs.get("check", True)
        func = log.info if check else log.error
        path = kwargs["path"]
        if path is not None and path.lower() == "global":
            path = None
        if isinstance(path, str) and not kwargs.get("allow_global"):
            log.info("Importing *local* emcee from %s", path)
            if not os.path.exists(path):
                func("The given folder does not exist: '%s'", path)
                return False
            if not os.path.exists(os.path.join(path, "setup.py")):
                func(
                    "Downloading emcee may have failed and no 'setup.py' file has been found in %s",
                    path,
                )
                return False
        elif not path:
            log.info("Importing *global* emcee.")
            path = None
        else:
            log.info("Importing *auto-installed* emcee (but defaulting to *global*).")
        try:
            return load_module("emcee", path=path, min_version=cls._min_emcee_version)
        except ImportError:
            if path is not None and path.lower() != "global":
                func(
                    "Couldn't find the emcee module at '%s'. Are you sure it has been installed there?",
                    path,
                )
            elif not check:
                log.error(
                    "Could not import global emcee installation. "
                    "Specify a Cobaya or emcee installation path, "
                    "or install the 'emcee' Python package globally."
                )
            return False
        except VersionCheckError as err:
            log.error(str(err))
            return False

    @classmethod
    def install(cls, path=None, code=True, no_progress_bars=False, **_kwargs):
        log = get_logger(cls.__name__)
        if not code:
            log.info("Code not requested. Nothing to do.")
            return True
        log.info("Installing pre-requisites...")
        packages = ["setuptools_scm", "pep517", "h5py", "schwimmbad"]
        exit_status = pip_install(packages)
        if exit_status:
            log.error("Could not install pre-requisite: %s", packages)
            return False

        # Here we first clone and then checkout the release tag. emcee use setuptools_scm which
        # needs the .git directory and github release do not provide it (so download_github_release
        # function works but the compilation process fails).
        def git_process(cmd, cwd, err_msg):
            from subprocess import PIPE, Popen

            process_git = Popen(cmd, cwd=cwd, stdout=PIPE, stderr=PIPE)
            out, err = process_git.communicate()
            if process_git.returncode:
                log.info(out.decode())
                log.info(err.decode())
                log.error(err_msg)
                return False
            return True

        log.info("Downloading emcee...")
        cmd = ["git", "clone", r"https://github.com/" + cls._emcee_repo_name, cls.get_path(path)]
        if not git_process(
            cmd=cmd,
            cwd=None,
            err_msg=f"Download failed! "
            f"Maybe the {cls._emcee_repo_name} repository does not exist or you misspelled its name "
            f"Only github repositories are supported.",
        ):
            return False

        cmd = ["git", "checkout", "-b", cls._emcee_repo_version, cls._emcee_repo_version]
        if not git_process(
            cmd=cmd,
            cwd=cls.get_path(path),
            err_msg="Checkout of {0} release failed! Make sure the {0} release tag exists.".format(
                cls._emcee_repo_version
            ),
        ):
            return False

        log.info("Installing emcee...")
        exit_status = pip_install(cls.get_path(path))
        if exit_status:
            log.error("Could not install emcee")
            return False
        return True
