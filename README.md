# Table of contents
1. [Overview of the Cobaya-CosmoLike Joint Architecture (Cocoa)](#overview)
2. [Installation of Cocoa's required packages via Conda](#required_packages_conda)
3. [Installation of Cobaya base code](#cobaya_base_code)
4. [Running Cobaya Examples](#cobaya_base_code_examples)
5. [Running Cosmolike projects](#running_cosmolike_projects)
6. [Creating Cosmolike projects (external readme)](Cocoa/projects/)
7. [Appendix](#appendix)
    1. [Proper Credits](#appendix_proper_credits)
    2. [The whovian-cocoa docker container](#appendix_jupyter_whovian)
    3. [Miniconda Installation](#overview_miniconda)
    4. [Compiling Boltzmann, CosmoLike and Likelihood codes separatelly](#appendix_compile_separatelly)
    5. [Warning about Weak Lensing YAML files in Cobaya](#appendix_example_runs)
    6. [Manual Blocking of Cosmolike Parameters](#manual_blocking_cosmolike)
    7. [Adding a new modified CAMB/CLASS to Cocoa (external readme)](Cocoa/external_modules/code)
    8. [Fine-tunning CAMB Accuracy](#camb_accuracy)
    9. [Bash/C/C++ Notes](#lectnotes)
    10. [Installation of Cocoa's required packages via Cocoa's internal cache](#required_packages_cache)
    11. [Setting-up conda environment for Machine Learning emulators](#ml_emulators)

## Overview of the [Cobaya](https://github.com/CobayaSampler)-[CosmoLike](https://github.com/CosmoLike) Joint Architecture (Cocoa) <a name="overview"></a>

Cocoa allows users to run [CosmoLike](https://github.com/CosmoLike) routines inside the [Cobaya](https://github.com/CobayaSampler) framework. [CosmoLike](https://github.com/CosmoLike) can analyze data primarily from the [Dark Energy Survey](https://www.darkenergysurvey.org) and simulate future multi-probe analyses for LSST and Roman Space Telescope. Besides integrating [Cobaya](https://github.com/CobayaSampler) and [CosmoLike](https://github.com/CosmoLike), Cocoa introduces shell scripts and readme instructions that allow users to containerize [Cobaya](https://github.com/CobayaSampler). The container structure ensures that users will adopt the same compiler and libraries (including their versions), and that they will be able to use multiple [Cobaya](https://github.com/CobayaSampler) instances consistently. This readme file presents basic and advanced instructions for installing all [Cobaya](https://github.com/CobayaSampler) and [CosmoLike](https://github.com/CosmoLike) components.

## Installation of Cocoa's required packages via Conda <a name="required_packages_conda"></a>

There are two installation methods. Users must choose one of them:

1. Via Conda - easier, best overall.
2. Via Cocoa's internal cache - slow, not advisable. See Appendix [Installation of Cocoa's required packages via Cocoa's internal cache](#required_packages_cache) 

We also provide the docker image whovian-cocoa to facilitate the installation of Cocoa on Windows and MacOS. For further instructions, refer to the Appendix [whovian-cocoa docker container](#appendix_jupyter_whovian).

We assume here the user has previously installed either [Minicoda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://www.anaconda.com/products/individual). If this is not the case, refer to the Appendix [Miniconda Installation](#overview_miniconda) for further instructions.

Type the following commands to create the cocoa Conda environment.

        conda create --name cocoapy38 python=3.8 --quiet --yes \
           && conda install -n cocoapy38 --quiet --yes  \
           'conda-forge::libgcc-ng=12.3.0' \
           'conda-forge::libstdcxx-ng=12.3.0' \
           'conda-forge::libgfortran-ng=12.3.0' \
           'conda-forge::gxx_linux-64=12.3.0' \
           'conda-forge::gcc_linux-64=12.3.0' \
           'conda-forge::gfortran_linux-64=12.3.0' \
           'conda-forge::openmpi=4.1.5' \
           'conda-forge::sysroot_linux-64=2.17' \
           'conda-forge::git=2.40.0' \
           'conda-forge::git-lfs=3.3.0' \
           'conda-forge::fftw=3.3.10' \
           'conda-forge::cfitsio=4.0.0' \
           'conda-forge::hdf5=1.14.0' \
           'conda-forge::lapack=3.9.0' \
           'conda-forge::openblas=0.3.23' \
           'conda-forge::lapack=3.9.0' \
           'conda-forge::gsl=2.7' \
           'conda-forge::cmake=3.26.4' \
           'conda-forge::xz==5.2.6' \
           'conda-forge::armadillo=11.4.4' \
           'conda-forge::boost-cpp=1.81.0' \
           'conda-forge::expat=2.5.0' \
           'conda-forge::cython=0.29.35' \
           'conda-forge::scipy=1.10.1' \
           'conda-forge::pandas=1.5.3' \
           'conda-forge::numpy=1.23.5' \
           'conda-forge::matplotlib=3.7.1' \
           'conda-forge::mpi4py=3.1.4'
      
For those working on projects that utilize machine-learning-based emulators, the Appendix [Setting-up conda environment for Machine Learning emulators](#ml_emulators) provides additional commands for installing the necessary packages.

When adopting this installation method, users must activate the Conda environment whenever working with Cocoa, as shown below.

        $ conda activate cocoapy38
    
Furthermore, users must install GIT-LFS on the first loading of the Conda cocoa environment.

        $(cocoapy38) git-lfs install

Users can now proceed to the section [Installation of Cobaya base code](#cobaya_base_code).

## Installation of Cobaya base code <a name="cobaya_base_code"></a>

Assuming the user opted for the easier *Conda installation*, type:

        $ conda activate cocoapy38
    
        $(cocoapy38) git clone --depth 1 https://github.com/CosmoLike/cocoa.git

to clone the repository. 

Cocoa developers should drop the shallow clone option `--depth 1`; they should also authenticate to GitHub via ssh keys and use the command instead

        $(cocoapy38) git clone git@github.com:CosmoLike/cocoa.git

:warning: **Warning** :warning: We have a limited monthly quota in bandwidth for Git LFS files, and therefore we ask users to use good judgment in the number of times they clone files from Cocoa's main repository.
 
Cocoa is made aware of the chosen installation method of required packages via special environment keys located on the *Cocoa/set_installation_options* script, as shown below:

        [Extracted from set_installation_options script]
        # --------------------------------------------------------------------------------------
        # --------------------------------------------------------------------------------------
        # --------------------------------------------------------------------------------------
        # ----------------------- HOW COCOA SHOULD BE INSTALLED? -------------------------------
        # --------------------------------------------------------------------------------------
        # --------------------------------------------------------------------------------------
        # --------------------------------------------------------------------------------------
        export MINICONDA_INSTALLATION=1
        #export MANUAL_INSTALLATION=1
    
The user must uncomment the appropriate key (here, we assume `MINICONDA_INSTALLATION`), and then type the following command

        $(cocoapy38) cd ./Cocoa/
        $(cocoapy38) source setup_cocoa_installation_packages

The script `setup_cocoa_installation_packages` decompresses the data files, which only takes a few minutes, and installs any remaining necessary packages. Typical package installation time ranges, depending on the installation method, from a few minutes (installation via Conda) to ~1/2 hour (installation via Cocoa's internal cache). It is important to note that our scripts never install packages on `$HOME/.local`. All requirements for Cocoa are installed at

        Cocoa/.local/bin
        Cocoa/.local/include
        Cocoa/.local/lib
        Cocoa/.local/share

This behavior is critical to enable users to work on multiple instances of Cocoa simultaneously.

Finally, type

        $(cocoapy38) source compile_external_modules
    
to compile CAMB, CLASS, Planck and Polychord. If the user wants to compile only a subset of these packages, then refer to the appendix [Compiling Boltzmann, CosmoLike and Likelihood codes separatelly](#appendix_compile_separatelly).

:books: **expert** :books: The script *set_installation_options script* contains a few additional flags that may be useful if something goes wrong. They are shown below:

        [Extracted from set_installation_options script]
        # --------------------------------------------------------------------------------------
        # --------- IF TRUE, THEN COCOA USES CLIK FROM https://github.com/benabed/clik ---------
        # --------------------------------------------------------------------------------------
        export USE_SPT_CLIK_PLANCK=1

        # --------------------------------------------------------------------------------------
        # ----------------- CONTROL OVER THE COMPILATION OF EXTERNAL CODES ---------------------
        # --------------------------------------------------------------------------------------
        #export IGNORE_CAMB_COMPILATION=1
        export IGNORE_CLASS_COMPILATION=1
        #export IGNORE_COSMOLIKE_COMPILATION=1
        #export IGNORE_POLYCHORD_COMPILATION=1
        #export IGNORE_PLANCK_COMPILATION=1
        #export IGNORE_ACT_COMPILATION=1

        # --------------------------------------------------------------------------------------
        # IN CASE COMPILATION FAILS, THESE FLAGS WILL BE USEFUL. BY DEFAULT, THE COMPILATION'S -
        # OUTPUT IS NOT WRITTEN ON THE TERMINAL. THESE FLAGS ENABLE THAT OUTPUT ---------------- 
        # --------------------------------------------------------------------------------------
        #export DEBUG_PLANCK_OUTPUT=1
        #export DEBUG_CAMB_OUTPUT=1
        #export DEBUG_CLASS_OUTPUT=1
        #export DEBUG_POLY_OUTPUT=1
        #export DEBUG_ACT_OUTPUT=1

## Running Cobaya Examples <a name="cobaya_base_code_examples"></a>

Assuming the user opted for the easier *Conda installation* and located the terminal at the folder *where Cocoa was cloned*, this is how to run some example YAML files we provide (*no Cosmolike code involved*): 

:one: **Step 1 of 5**: activate the conda environment

        $ conda activate cocoapy38
     
:two: **Step 2 of 5**: go to the Cocoa main folder 

        $(cocoapy38) cd ./cocoa/Cocoa

3️⃣ **Step 3 of 5**: activate the private python environment

        $(cocoapy38) source start_cocoa

Users will see a terminal that looks like this: `$(cocoapy38)(.local)`. *This is a feature, not a bug*! 

Why did we choose to have two separate bash environments? Users should be able to manipulate multiple Cocoa instances seamlessly, which is particularly useful when running chains in one instance while experimenting with code development in another. Consistency of the environment across all Cocoa instances is crucial, and the start_cocoa/stop_cocoa scripts handle the loading and unloading of environmental path variables for each Cocoa. All of them, however, depends on many of the same prerequisites, so it is advantageous to maintain the basic packages inside the shared conda cocoa environment. 

:four: **Step 4 of 5**: select the number of OpenMP cores
    
        $(cocoapy38)(.local) export OMP_PROC_BIND=close; export OMP_NUM_THREADS=4

:five: **Step 5 of 5**: run `cobaya-run` on a the first example YAML files we provide.

One model evaluation:

        $(cocoapy38)(.local) mpirun -n 1 --oversubscribe --mca btl vader,tcp,self --bind-to core:overload-allowed --rank-by core --map-by numa:pe=${OMP_NUM_THREADS} cobaya-run  ./projects/example/EXAMPLE_EVALUATE1.yaml -f
        
PS: We offer the flag `COCOA_RUN_EVALUATE` as an alias (syntax-sugar) for `mpirun -n 1 --oversubscribe --mca btl vader,tcp,self --bind-to core:overload-allowed --rank-by core --map-by numa:pe=4 cobaya-run`. 

MCMC:

        $(cocoapy38)(.local) mpirun -n 4 --oversubscribe --mca btl vader,tcp,self --bind-to core:overload-allowed --rank-by core --map-by numa:pe=${OMP_NUM_THREADS} cobaya-run ./projects/example/EXAMPLE_MCMC1.yaml -f

PS: We offer the flag `COCOA_RUN_MCMC` as an alias (syntax-sugar) for `mpirun -n 4 --oversubscribe --mca btl vader,tcp,self --bind-to core:overload-allowed --rank-by core --map-by numa:pe=4 cobaya-run`. 

:books: **expert** :books: Why the `--mca btl vader,tcp,self` flag? Conda-forge developers don't [compile OpenMPI with Infiniband compatibility](https://github.com/conda-forge/openmpi-feedstock/issues/38).

:books: **expert** :books: Why the `--bind-to core:overload-allowed --map-by numa:pe=${OMP_NUM_THREADS}` flag? This flag enables efficient hybrid MPI + OpenMP runs on NUMA architecture.

Once the work is done, type:

        $(cocoapy38)(.local) source stop_cocoa
        $(cocoapy38) conda deactivate cocoa

## Running Cosmolike projects <a name="running_cosmolike_projects"></a> 

The *projects* folder was designed to include Cosmolike projects. Similar to the previous section, we assume the user opted for the more direct *Conda installation* method. We also presume the user's terminal is in the folder where Cocoa was cloned.

:one: **Step 1 of 5**: activate the Conda Cocoa environment
    
        $ conda activate cocoapy38

:two: **Step 2 of 5**: go to the project folder (`./cocoa/Cocoa/projects`) and clone a Cosmolike project, with fictitious name `XXX`:
    
        $(cocoapy38) cd ./cocoa/Cocoa/projects
        $(cocoapy38) $CONDA_PREFIX/bin/git clone git@github.com:CosmoLike/cocoa_XXX.git XXX

By convention, the Cosmolike Organization hosts a Cobaya-Cosmolike project named XXX at `CosmoLike/cocoa_XXX`. However, our provided scripts and template YAML files assume the removal of the `cocoa_` prefix when cloning the repository.

Example of cosmolike projects: [lsst_y1](https://github.com/CosmoLike/cocoa_lsst_y1).
 
:three: **Step 3 of 5**: go back to Cocoa main folder, and activate the private python environment
    
        $(cocoapy38) cd ../
        $(cocoapy38) source start_cocoa
 
:warning: (**warning**) :warning: Remember to run the start_cocoa script only after cloning the project repository. The script *start_cocoa* creates the necessary symbolic links and adds the *Cobaya-Cosmolike interface* of all projects to `LD_LIBRARY_PATH` and `PYTHONPATH` paths.

:four: **Step 4 of 5**: compile the project
 
        $(cocoapy38)(.local) source ./projects/XXX/scripts/compile_XXX
  
:five:  **Step 5 of 5**: select the number of OpenMP cores and run a template YAML file
    
        $(cocoapy38)(.local) export OMP_PROC_BIND=close; export OMP_NUM_THREADS=4
        $(cocoapy38)(.local) mpirun -n 1 --oversubscribe --mca btl vader,tcp,self --bind-to core --rank-by core --map-by numa:pe=${OMP_NUM_THREADS} cobaya-run ./projects/XXX/EXAMPLE_EVALUATE1.yaml -f

:warning: **Warning** :warning: Be careful when creating YAML for weak lensing projects in Cobaya using the $\Omega_m/\Omega_b$ parameterization. See Appendix [warning about weak lensing YAML files](#appendix_example_runs) for further details.

## Appendix <a name="appendix"></a>

### Proper Credits <a name="appendix_proper_credits"></a>

The following is not an exhaustive list of the codes we use

- [Cobaya](https://github.com/CobayaSampler) is a framework developed by Dr. Jesus Torrado and Prof. Anthony Lewis

- [Cosmolike](https://github.com/CosmoLike) is a framework developed by Prof. Elisabeth Krause and Prof. Tim Eifler

- [CAMB](https://github.com/cmbant/CAMB) is a Boltzmann code developed by Prof. Anthony Lewis

- [CLASS](https://github.com/lesgourg/class_public) is a Boltzmann code developed by Prof. Julien Lesgourgues and Dr. Thomas Tram

- [Polychord](https://github.com/PolyChord/PolyChordLite) is a sampler code developed by Dr. Will Handley, Prof. Lasenby, and Prof. M. Hobson

By no means, we want to discourage people from cloning code from their original repositories. We've included these codes as compressed [xz file format](https://tukaani.org/xz/format.html) in our repository for convenience in the initial development. The work of those authors is extraordinary, and they must be properly cited.

### The whovian-cocoa docker container <a name="appendix_jupyter_whovian"></a>

We provide the docker image [whovian-cocoa](https://hub.docker.com/r/vivianmiranda/whovian-cocoa) to facilitate the installation of Cocoa on Windows and MacOS. This appendix assumes the users already have the docker engine installed on their local PC. For instructions on installing the docker engine in specific operating systems, please refer to [Docker's official documentation](https://docs.docker.com/engine/install/). 

To download and run the container for the first time, type:

         docker run --platform linux/amd64 --hostname cocoa --name cocoa2023 -it -p 8080:8888 -v $(pwd):/home/whovian/host/ -v ~/.ssh:/home/whovian/.ssh:ro vivianmiranda/whovian-cocoa

Following the command above, users should see the following text on the screen terminal

<img width="1121" alt="Screen Shot 2023-06-09 at 2 48 23 PM" src="https://github.com/vivianmiranda/cocoa3_testing/assets/3210728/8ad62524-ed62-4560-9e8e-0187a8f8cb3b">

:books: **expert** :books: The flag `-v $(pwd):/home/whovian/host/` ensures that files on the host computer, where the user should install Cocoa to avoid losing work in case the docker image needs to be deleted, have been mounted to the directory `/home/whovian/host/`. Therefore, the user should see the host files of the directory where the whovian-cocoa container was initialized after typing `cd /home/whovian/host/; ls`

When running the container the first time, the user needs to init conda with `conda init bash` followed by `source ~/.bashrc`, as shown below.

        whovian@cocoa:~$ conda init bash
        whovian@cocoa:~$ source ~/.bashrc

The container already comes with conda Cocoa environment pre-installed:

        whovian@cocoa:~$ conda activate cocoa

When the user exits the container, how to restart it? Type 
    
        $ docker start -ai cocoa2023

How to run Jupyter Notebooks remotely when using Cocoa within the whovian-cocoa container? First, type the following command:

        whovian@cocoa:~$ jupyter notebook --no-browser --port=8080

The terminal will show a message similar to the following template:

        [... NotebookApp] Writing notebook server cookie secret to /home/whovian/.local/share/jupyter/runtime/notebook_cookie_secret
        [... NotebookApp] WARNING: The notebook server is listening on all IP addresses and not using encryption. This is not recommended.
        [... NotebookApp] Serving notebooks from local directory: /home/whovian/host
        [... NotebookApp] Jupyter Notebook 6.1.1 is running at:
        [... NotebookApp] http://f0a13949f6b5:8888/?token=XXX
        [... NotebookApp] or http://127.0.0.1:8888/?token=XXX
        [... NotebookApp] Use Control-C to stop this server and shut down all kernels (twice to skip confirmation).

Below, we assume the user runs the container in a server with the URL `your_sever.com`. We also presume the server can be accessed via ssh protocol. From a local PC, type:

        $ ssh your_username@your_sever.com -L 8080:localhost:8080

Finally, go to a browser and type `http://localhost:8080/?token=XXX`, where `XXX` is the previously saved token.

PS: The docker container also has the conda environment `cocoalite` that is useful in the rare case someone want to install Cocoa via the slow/not-advisable instructions on section [Installation of Cocoa's required packages via Cocoa's internal cache](#required_packages_cache)

### Miniconda Installation <a name="overview_miniconda"></a>

Download and run Miniconda installation script (please adapt `CONDA_DIR`):

        export CONDA_DIR=/gpfs/home/vinmirandabr/miniconda

        mkdir $CONDA_DIR

        wget https://repo.continuum.io/miniconda/Miniconda3-py38_4.12.0-Linux-x86_64.sh

        /bin/bash Miniconda3-py38_4.12.0-Linux-x86_64.sh -f -b -p $CONDA_DIR

After installation, users must source conda configuration file

        source $CONDA_DIR/etc/profile.d/conda.sh \
            && conda config --set auto_update_conda false \
            && conda config --set show_channel_urls true \
            && conda config --set auto_activate_base false \
            && conda config --prepend channels conda-forge \
            && conda config --set channel_priority strict \
            && conda init bash
    
### Compiling Boltzmann, CosmoLike and Likelihood codes separatelly <a name="appendix_compile_separatelly"></a>

To avoid excessive compilation times during development, users can use specialized scripts located at `Cocoa/installation_scripts/` that compile only a specific module. A few examples of these scripts are: 

        $(cocoapy38)(.local) source ./installation_scripts/compile_class
        $(cocoapy38)(.local) source ./installation_scripts/compile_camb
        $(cocoapy38)(.local) source ./installation_scripts/compile_planck
        $(cocoapy38)(.local) source ./installation_scripts/compile_act
        $(cocoapy38)(.local) source ./installation_scripts/setup_polychord
    
### :warning: Warning :warning: Weak Lensing YAML files in Cobaya <a name="appendix_example_runs"></a>

The CosmoLike pipeline takes $\Omega_m$ and $\Omega_b$, but the CAMB Boltzmann code only accepts $\Omega_c h^2$ and $\Omega_b h^2$ in Cobaya. Therefore, there are two ways of creating YAML compatible with CAMB and Cosmolike: 

1. CMB parameterization: $\big(\Omega_c h^2,\Omega_b h^2\big)$ as primary MCMC parameters and $\big(\Omega_m,\Omega_b\big)$ as derived quantities.

        omegabh2:
            prior:
                min: 0.01
                max: 0.04
            ref:
                dist: norm
                loc: 0.022383
                scale: 0.005
            proposal: 0.005
            latex: \Omega_\mathrm{b} h^2
        omegach2:
            prior:
                min: 0.06
                max: 0.2
            ref:
                dist: norm
                loc: 0.12011
                scale: 0.03
            proposal: 0.03
            latex: \Omega_\mathrm{c} h^2
        mnu:
            value: 0.06
            latex: m_\\nu
        omegam:
            latex: \Omega_\mathrm{m}
        omegamh2:
            derived: 'lambda omegam, H0: omegam*(H0/100)**2'
            latex: \Omega_\mathrm{m} h^2
        omegab:
            derived: 'lambda omegabh2, H0: omegabh2/((H0/100)**2)'
            latex: \Omega_\mathrm{b}
        omegac:
            derived: 'lambda omegach2, H0: omegach2/((H0/100)**2)'
            latex: \Omega_\mathrm{c}

2. Weak Lensing parameterization: $\big(\Omega_m,\Omega_b\big)$ as primary MCMC parameters and $\big(\Omega_c h^2, \Omega_b h^2\big)$ as derived quantities.

Adopting $\big(\Omega_m,\Omega_b\big)$ as main MCMC parameters can create a silent bug in Cobaya. The problem occurs when the option `drop: true` is absent in $\big(\Omega_m,\Omega_b\big)$ parameters, and there are no expressions that define the derived $\big(\Omega_c h^2, \Omega_b h^2\big)$ quantities. The bug is silent because the MCMC runs without any warnings, but the CAMB Boltzmann code does not update the cosmological parameters at every MCMC iteration. As a result, the resulting posteriors are flawed, but they may seem reasonable to those unfamiliar with the issue. It's important to be aware of this bug to avoid any potential inaccuracies in the results. 

The correct way to create YAML files with $\big(\Omega_m,\Omega_b\big)$ as primary MCMC parameters is exemplified below

        omegab:
            prior:
                min: 0.03
                max: 0.07
            ref:
                dist: norm
                loc: 0.0495
                scale: 0.004
            proposal: 0.004
            latex: \Omega_\mathrm{b}
            drop: true
        omegam:
            prior:
                min: 0.1
                max: 0.9
            ref:
                dist: norm
                loc: 0.316
                scale: 0.02
            proposal: 0.02
            latex: \Omega_\mathrm{m}
            drop: true
        mnu:
            value: 0.06
            latex: m_\\nu
        omegabh2:
            value: 'lambda omegab, H0: omegab*(H0/100)**2'
            latex: \Omega_\mathrm{b} h^2
        omegach2:
            value: 'lambda omegam, omegab, mnu, H0: (omegam-omegab)*(H0/100)**2-(mnu*(3.046/3)**0.75)/94.0708'
            latex: \Omega_\mathrm{c} h^2
        omegamh2:
            derived: 'lambda omegam, H0: omegam*(H0/100)**2'
            latex: \Omega_\mathrm{m} h^2

### Manual Blocking of Cosmolike Parameters <a name="manual_blocking_cosmolike"></a>

Cosmolike Weak Lensing pipeline contains parameters with different speed hierarchies. For example, Cosmolike execution time is reduced by approximately 50% when fixing the cosmological parameters. When varying only multiplicative shear calibration, Cosmolike execution time is reduced by two orders of magnitude. 

Cobaya can't automatically handle parameters associated with the same likelihood that have different speed hierarchies. Luckily, we can manually impose the speed hierarchy in Cobaya using the `blocking:` option. The only drawback of this method is that parameters of all adopted likelihoods need to be manually specified, not only the ones required by Cosmolike.

In addition to that, Cosmolike can't cache the intermediate products of the last two evaluations, which is necessary to exploit optimizations associated with dragging (`drag: True`). However, Cosmolike caches the intermediate products of the previous evaluation, thereby enabling the user to take advantage of the slow/fast decomposition of parameters in Cobaya's main MCMC sampler. 

Below we provide an example YAML configuration for an MCMC chain that with DES 3x2pt likelihood.

        likelihood: 
            des_y3.des_3x2pt:
            path: ./external_modules/data/des_y3
            data_file: DES_Y1.dataset
         
         (...)
         
        sampler:
            mcmc:
                covmat: "./projects/des_y3/EXAMPLE_MCMC22.covmat"
                covmat_params:
                # ---------------------------------------------------------------------
                # Proposal covariance matrix learning
                # ---------------------------------------------------------------------
                learn_proposal: True
                learn_proposal_Rminus1_min: 0.03
                learn_proposal_Rminus1_max: 100.
                learn_proposal_Rminus1_max_early: 200.
                # ---------------------------------------------------------------------
                # Convergence criteria
                # ---------------------------------------------------------------------
                # Maximum number of posterior evaluations
                max_samples: .inf
                # Gelman-Rubin R-1 on means
                Rminus1_stop: 0.02
                # Gelman-Rubin R-1 on std deviations
                Rminus1_cl_stop: 0.2
                Rminus1_cl_level: 0.95
                # ---------------------------------------------------------------------
                # Exploiting Cosmolike speed hierarchy
                # ---------------------------------------------------------------------
                measure_speeds: False
                drag: False
                oversample_power: 0.2
                oversample_thin: True
                blocking:
                - [1,
                    [
                        As_1e9, H0, omegab, omegam, ns
                    ]
                  ]
                - [4,
                    [
                        DES_DZ_S1, DES_DZ_S2, DES_DZ_S3, DES_DZ_S4, DES_A1_1, DES_A1_2,
                        DES_B1_1, DES_B1_2, DES_B1_3, DES_B1_4, DES_B1_5,
                        DES_DZ_L1, DES_DZ_L2, DES_DZ_L3, DES_DZ_L4, DES_DZ_L5
                    ]
                  ]
                - [25,
                    [
                        DES_M1, DES_M2, DES_M3, DES_M4, DES_PM1, DES_PM2, DES_PM3, DES_PM4, DES_PM5
                    ]
                  ]
                # ---------------------------------------------------------------------
                max_tries: 10000
                burn_in: 0
                Rminus1_single_split: 4
                


### Fine-tunning CAMB Accuracy <a name="camb_accuracy"></a>

The accurate computation of many CMB and large-scale-structure data vectors requires high `AccuracyBoost` values in CAMB. However, this parameter is particularly inefficient, causing an exponential increase in CAMB's runtime. This issue has been frequent enough that we provide below a simple but partial remedy. 

The underlying reason for `AccuracyBoost` inefficiency is that this flag raises the required accuracy of multiple modules in CAMB. The appropriate boost must be fine-tuned until the $\chi^2$ of the adopted experiments remain stable. However, we do not need to raise the boost factor in all CAMB modules by the same amount to achieve such stability. 

The Python function `set_accuracy`,  located in the file `$ROOTDIR/external_modules/code/CAMB/camb`, can be modified for a more fine-tuned change to CAMB accuracy. Below is an example of possible modifications: 

    def set_accuracy(self, AccuracyBoost=1., lSampleBoost=1., lAccuracyBoost=1., DoLateRadTruncation=True): 
        #COCOA: BEGINS
        self.Accuracy.AccuracyBoost = AccuracyBoost       
        
        #COCOA: below, actual boost = 1.0 + 2*(AccuracyBoost-1.0)
        self.Accuracy.BessIntBoost = (1.0 + 2*(AccuracyBoost-1.0))/AccuracyBoost
        self.Accuracy.IntkAccuracyBoost = (1.0 + 2*(AccuracyBoost-1.0))/AccuracyBoost
        self.Accuracy.TransferkBoost = (1.0 + 2*(AccuracyBoost-1.0))/AccuracyBoost
        self.Accuracy.KmaxBoost = (1.0 + 2*(AccuracyBoost-1.0))/AccuracyBoost
        self.Accuracy.TimeStepBoost = (1.0 + 2*(AccuracyBoost-1.0))/AccuracyBoost
        
        #COCOA: below, actual boost = 1.0 + 5*(AccuracyBoost-1.0)
        self.Accuracy.SourcekAccuracyBoost = (1.0 + 5*(AccuracyBoost-1.0))/AccuracyBoost
        self.Accuracy.BesselBoost = (1.0 + 5*(AccuracyBoost-1.0))/AccuracyBoost
        #COCOA: ENDS
        
        self.Accuracy.lSampleBoost = lSampleBoost
        self.Accuracy.lAccuracyBoost = lAccuracyBoost
        self.DoLateRadTruncation = DoLateRadTruncation
        return self
        
With the code above, the theoretical error in Simons Observatory $\chi^2$ seems to be under control (i.e., $\Delta \chi^2 =$ O(few)) with `AccuracyBoost: 1.06` and `lens_potential_accuracy: 4` even a bit away from the best-fit model so that chains can be later corrected via Importance Sampling. As a reminder, corrections based on Importance Sampling are much faster when compared to running MCMC chains with insane accuracy because they can be computed on thinned versions of converged chains and are trivially parallelizable. 

Out of caution, we have not implemented these changes in `$ROOTDIR/external_modules/code/CAMB/`.
        
### :book: Bash/C/C++ Notes :book: <a name="lectnotes"></a>

To effectively work with the Cobaya framework and Cosmolike codes at the developer level, a working knowledge of Python to understand Cobaya and Bash language to comprehend Cocoa's scripts is required. Proficiency in C and C++ is also needed to manipulate Cosmolike and the C++ Cobaya-Cosmolike C++ interface. Finally, users need to understand the Fortran-2003 language to modify CAMB.

Learning all these languages can be overwhelming, so to enable new users to do research that demands modifications on the inner workings of these codes, we include [here](cocoa_installation_libraries/LectNotes.pdf) a link to approximately 600 slides that provide an overview of Bash (slides 1-137), C (slides 138-371), and C++ (slides 372-599). In the future, we aim to add lectures about Python and Fortran. 

### Installation of Cocoa's required packages via Cocoa's internal cache <a name="required_packages_cache"></a>

This method is slow and not advisable :stop_sign::thumbsdown:. When Conda is unavailable, the user can still perform a local semi-autonomous installation on Linux based on a few scripts we implemented. We provide a local copy of almost all required packages on Cocoa's cache folder named *cocoa_installation_libraries*. We assume the pre-installation of the following packages:

   - [Bash](https://www.amazon.com/dp/B0043GXMSY/ref=cm_sw_em_r_mt_dp_x3UoFbDXSXRBT);
   - [Git](https://git-scm.com) v1.8+;
   - [Git LFS](https://git-lfs.github.com);
   - [gcc](https://gcc.gnu.org) v12.*;
   - [gfortran](https://gcc.gnu.org) v12.*;
   - [g++](https://gcc.gnu.org) v12.*;
   - [Python](https://www.python.org) v3.8.*;
   - [PIP package manager](https://pip.pypa.io/en/stable/installing/)
   - [Python Virtual Environment](https://www.geeksforgeeks.org/python-virtual-environment/)

The conda environment `cocoalitepy38` contains the minimum packages necessary for this installation method

        conda create --name cocoalitepy38 python=3.8 --quiet --yes \
        && conda install -n cocoalitepy38 --quiet --yes  \
            'conda-forge::libgcc-ng=12.3.0' \
            'conda-forge::libstdcxx-ng=12.3.0' \
            'conda-forge::libgfortran-ng=12.3.0' \
            'conda-forge::gxx_linux-64=12.3.0' \
            'conda-forge::gcc_linux-64=12.3.0' \
            'conda-forge::gfortran_linux-64=12.3.0' \
            'conda-forge::openmpi=4.1.5' \
            'conda-forge::sysroot_linux-64=2.17' \
            'conda-forge::git=2.40.0' \
            'conda-forge::git-lfs=3.3.0'
    
To perform the local semi-autonomous installation, users must modify flags written on the file *set_installation_options* because the default behavior corresponds to an installation via Conda. First, select the environmental key `MANUAL_INSTALLATION` as shown below:

    [Extracted from set_installation_options script] 
    
    # --------------------------------------------------------------------------------------
    # --------------------------------------------------------------------------------------
    # --------------------------------------------------------------------------------------
    # ----------------------- HOW COCOA SHOULD BE INSTALLED? -------------------------------
    # --------------------------------------------------------------------------------------
    # --------------------------------------------------------------------------------------
    # --------------------------------------------------------------------------------------
    #export MINICONDA_INSTALLATION=1
    export MANUAL_INSTALLATION=1
    
Finally, set the following environmental keys
 
    [Extracted from set_installation_options script]
  
    if [ -n "${MANUAL_INSTALLATION}" ]; then
        # --------------------------------------------------------------------------------------
        # IF SET, COCOA DOES NOT USE SYSTEM PIP PACKAGES (RELIES EXCLUSIVELY ON PIP CACHE FOLDER)
        # --------------------------------------------------------------------------------------
        export DONT_USE_SYSTEM_PIP_PACKAGES=1

        # --------------------------------------------------------------------------------------
        # IF SET, COCOA WILL NOT INSTALL TENSORFLOW, KERAS, PYTORCH, GPY
        # --------------------------------------------------------------------------------------
        export IGNORE_EMULATOR_CPU_PIP_PACKAGES=1
        export IGNORE_EMULATOR_GPU_PIP_PACKAGES=1

        # --------------------------------------------------------------------------------------
        # WE USE CONDA COLASLIM ENV WITH JUST PYTHON AND GCC TO TEST MANUAL INSTALLATION
        # --------------------------------------------------------------------------------------
        #conda create --name cocoalite python=3.8 --quiet --yes \
        #   && conda install -n cocoapy38 --quiet --yes  \
        #   'conda-forge::libgcc-ng=12.3.0' \
        #   'conda-forge::libstdcxx-ng=12.3.0' \
        #   'conda-forge::libgfortran-ng=12.3.0' \
        #   'conda-forge::gxx_linux-64=12.3.0' \
        #   'conda-forge::gcc_linux-64=12.3.0' \
        #   'conda-forge::gfortran_linux-64=12.3.0' \
        #   'conda-forge::openmpi=4.1.5' \
        #   'conda-forge::sysroot_linux-64=2.17' \
        #   'conda-forge::git=2.40.0' \
        #   'conda-forge::git-lfs=3.3.0'
        # --------------------------------------------------------------------------------------

        export GLOBAL_PACKAGES_LOCATION=$CONDA_PREFIX
        export GLOBALPYTHON3=$CONDA_PREFIX/bin/python${PYTHON_VERSION}
        export PYTHON_VERSION=3.8

        # --------------------------------------------------------------------------------------
        # COMPILER
        # --------------------------------------------------------------------------------------
        export C_COMPILER=$CONDA_PREFIX/bin/x86_64-conda-linux-gnu-cc
        export CXX_COMPILER=$CONDA_PREFIX/bin/x86_64-conda-linux-gnu-g++
        export FORTRAN_COMPILER=$CONDA_PREFIX/bin/x86_64-conda-linux-gnu-gfortran
        export MPI_CC_COMPILER=$CONDA_PREFIX/bin/mpicxx
        export MPI_CXX_COMPILER=$CONDA_PREFIX/bin/mpicc
        export MPI_FORTRAN_COMPILER=$CONDA_PREFIX/bin/mpif90

        # --------------------------------------------------------------------------------------
        # USER NEEDS TO SPECIFY THE FLAGS BELOW SO COCOA CAN FIND PYTHON / GCC
        # --------------------------------------------------------------------------------------
        export PATH=$CONDA_PREFIX/bin:$PATH

        export CFLAGS="${CFLAGS} -I$CONDA_PREFIX/include"

        export LDFLAGS="${LDFLAGS} -L$CONDA_PREFIX/lib"

        export C_INCLUDE_PATH=$CONDA_PREFIX/include:$C_INCLUDE_PATH
        export C_INCLUDE_PATH=$CONDA_PREFIX/include/python${PYTHON_VERSION}m/:$C_INCLUDE_PATH

        export CPLUS_INCLUDE_PATH=$CONDA_PREFIX/include:$CPLUS_INCLUDE_PATH
        export CPLUS_INCLUDE_PATH=$CONDA_PREFIX/include/python${PYTHON_VERSION}m/:$CPLUS_INCLUDE_PATH

        export PYTHONPATH=$CONDA_PREFIX/lib/python$PYTHON_VERSION/site-packages:$PYTHONPATH
        export PYTHONPATH=$CONDA_PREFIX/lib:$PYTHONPATH

        export LD_RUN_PATH=$CONDA_PREFIX/lib/python$PYTHON_VERSION/site-packages:$LD_RUN_PATH
        export LD_RUN_PATH=$CONDA_PREFIX/lib:$LD_RUN_PATH

        export LIBRARY_PATH=$CONDA_PREFIX/lib/python$PYTHON_VERSION/site-packages:$LIBRARY_PATH
        export LIBRARY_PATH=$CONDA_PREFIX/lib:$LIBRARY_PATH

        export CMAKE_INCLUDE_PATH=$CONDA_PREFIX/include/:$CMAKE_INCLUDE_PATH
        export CMAKE_INCLUDE_PATH=$CONDA_PREFIX/include/python${PYTHON_VERSION}m/:$CMAKE_INCLUDE_PATH    

        export CMAKE_LIBRARY_PATH=$CONDA_PREFIX/lib/python$PYTHON_VERSION/site-packages:$CMAKE_LIBRARY_PATH
        export CMAKE_LIBRARY_PATH=$CONDA_PREFIX/lib:$CMAKE_LIBRARY_PATH

        export INCLUDE_PATH=$CONDA_PREFIX/include/:$INCLUDE_PATH

        export INCLUDEPATH=$CONDA_PREFIX/include/:$INCLUDEPATH

        export INCLUDE=$CONDA_PREFIX/x86_64-conda-linux-gnu/include:$INCLUDE
        export INCLUDE=$CONDA_PREFIX/include/:$INCLUDE

        export CPATH=$CONDA_PREFIX/include/:$CPATH

        export OBJC_INCLUDE_PATH=$CONDA_PREFIX/include/:OBJC_INCLUDE_PATH

        export OBJC_PATH=$CONDA_PREFIX/include/:OBJC_PATH

        # --------------------------------------------------------------------------------------
        # DEBUG THE COMPILATION OF PREREQUISITES PACKAGES. BY DEFAULT, THE COMPILATION'S -------
        # OUTPUT IS NOT WRITTEN ON THE TERMINAL. THESE FLAGS ENABLE THAT OUTPUT ---------------- 
        # --------------------------------------------------------------------------------------
        #export DEBUG_CPP_PACKAGES=1
        #export DEBUG_C_PACKAGES=1
        #export DEBUG_FORTRAN_PACKAGES=1
        #export DEBUG_PIP_OUTPUT=1
        #export DEBUG_XZ_PACKAGE=1
        #export DEBUG_CMAKE_PACKAGE=1
        #export DEBUG_OPENBLAS_PACKAGE=1
        #export DEBUG_DISTUTILS_PACKAGE=1
        #export DEBUG_HDF5_PACKAGES=1
    
The fine-tunning over the use of system-wide packages instead of our local copies can be set via the environmental flags

        export IGNORE_XZ_INSTALLATION=1
        export IGNORE_HDF5_INSTALLATION=1
        export IGNORE_CMAKE_INSTALLATION=1
        export IGNORE_DISTUTILS_INSTALLATION=1
        export IGNORE_C_GSL_INSTALLATION=1
        export IGNORE_C_CFITSIO_INSTALLATION=1
        export IGNORE_C_FFTW_INSTALLATION=1
        export IGNORE_CPP_BOOST_INSTALLATION=1 
        export IGNORE_OPENBLAS_INSTALLATION=1
        export IGNORE_FORTRAN_LAPACK_INSTALLATION=1
        export IGNORE_CPP_ARMA_INSTALLATION=1
       
    
Users can now proceed to the section [Installation of Cobaya base code](#cobaya_base_code)

### Setting-up conda environment for Machine Learning emulators <a name="ml_emulators"></a>

If the user wants to add Tensorflow, Keras and Pytorch for an emulator-based project via Conda, then type

        $ conda activate cocoapy38 
      
        $(cocoapy38) $CONDA_PREFIX/bin/pip install --no-cache-dir \
            'tensorflow-cpu==2.12.0' \
            'keras==2.12.0' \
            'keras-preprocessing==1.1.2' \
            'torch==1.13.1+cpu' \
            'torchvision==0.14.1+cpu' \
            'torchaudio==0.13.1' --extra-index-url https://download.pytorch.org/whl/cpu

In case there are GPUs available, the following commands will install the GPU version of 
Tensorflow, Keras and Pytorch (assuming CUDA 11.6, click [here](https://pytorch.org/get-started/previous-versions/) for additional information).

        $(cocoapy38) $CONDA_PREFIX/bin/pip install --no-cache-dir \
            'tensorflow==2.12.0' \
            'keras==2.12.0' \
            'keras-preprocessing==1.1.2' \
            'torch==1.13.1+cu116' \
            'torchvision==0.14.1+cu116' \
            'torchaudio==0.13.1' --extra-index-url https://download.pytorch.org/whl/cu116

Based on our experience, we recommend utilizing the GPU versions to train the emulator while using the CPU versions to run the MCMCs. This is because our supercomputers possess a greater number of CPU-only nodes. It may be helpful to create two separate conda environments for this purpose. One could be named `cocoa` (CPU-only), while the other could be named `cocoaemu` and contain the GPU versions of the machine learning packages.

Commenting out the environmental flags shown below, located at *set_installation_options* script, will enable the installation of machine-learning-related libraries via pip.  

        # IF TRUE, THEN COCOA WON'T INSTALL TENSORFLOW, KERAS and PYTORCH
        #export IGNORE_EMULATOR_CPU_PIP_PACKAGES=1
        #export IGNORE_EMULATOR_GPU_PIP_PACKAGES=1

Unlike most installed pip prerequisites, which are cached at `cocoa_installation_libraries/pip_cache.xz`, the installation of the Machine Learning packages listed above requires an active internet connection.
