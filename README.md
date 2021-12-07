# Table of contents
1. [Overview of the Cobaya-CosmoLike Joint Architecture (Cocoa)](#overview)
2. [Installation of cocoa required packages](#required_packages)
    1. [Via Conda (best for Linux)](#required_packages_conda)
    2. [Via Docker (best for MacOS/Windows)](#required_packages_docker)
    3. [Via Homebrew](#required_packages_homebrew)
    4. [Via Cocoa's internal cache](#required_packages_cache)
3. [Installation of cocoa base code](#cocoa_base_code)
    1. [Running Examples](#cocoa_base_code_examples)
5. [Download/compiling/running specific cocoa projects](#cocoa_projects)
6. [Appendix](#appendix)
    1. [Compiling Boltzmann, CosmoLike and Likelihood codes separatelly](#appendix_compile_separatelly)
    2. [Running Jupyter Notebooks inside the Whovian-Cosmo docker container](#appendix_jupyter_whovian)
    3. [Summary Information about Cocoa's configuration files](#appendix_config_files)

## Overview of the [Cobaya](https://github.com/CobayaSampler)-[CosmoLike](https://github.com/CosmoLike) Joint Architecture (Cocoa) <a name="overview"></a>

Cocoa allows users to run [CosmoLike](https://github.com/CosmoLike) routines inside the [Cobaya](https://github.com/CobayaSampler) framework. Cosmolike is capable of analyzing data primarily from the [Dark Energy Survey](https://www.darkenergysurvey.org) (a.k.a DES) and simulating future multi-probe analyses for Rubin Observatory's Legacy Survey of Space and Time or the Roman Space Telescope. This readme file presents basic and advanced instructions for installing all Cocoa components, including the [Planck likelihood](https://wiki.cosmos.esa.int/planck-legacy-archive/index.php/Main_Page).

(**Proper Credits**): The following is not an exhaustive list of the codes we use

- [Cobaya](https://github.com/CobayaSampler) is a framework developed by Dr. Jesus Torrado and Prof. Anthony Lewis

- [Cosmolike](https://github.com/CosmoLike) is a framework developed by Prof. Elisabeth Krause and Prof. Tim Eifler

- [CAMB](https://github.com/cmbant/CAMB) is a Boltzmann code developed by Prof. Anthony Lewis

- [CLASS](https://github.com/lesgourg/class_public) is a Boltzmann code developed by Prof. Julien Lesgourgues and Dr. Thomas Tram

- [Polychord](https://github.com/PolyChord/PolyChordLite) is a sampler code developed by Dr. Will Handley, Prof. Lasenby, and Prof. M. Hobson

By no means, we want to discourage people from cloning code from their original repositories. We've included these codes as compressed [xz file format](https://tukaani.org/xz/format.html) in our repository for convenience in the initial development. The work of those authors is extraordinary, and they must be properly cited.

## Installation of cocoa required packages <a name="required_packages"></a>

Cosmolike and the interface between Cosmolike and Cocoa requires many packages to be installed, including [GSL](https://www.gnu.org/software/gsl/), [FFTW](https://www.fftw.org), [Armadillo](http://arma.sourceforge.net) and [Boost](https://www.boost.org). The plethora of GCC, Python, and package versions, each one with different bugs and regressions, can make the installation of any big code to be pure agony, especially given that CAMB, CLASS, Cosmolike, and Planck likelihood involves Fortran, C, C++, and Python languages. We try to simplify this process by offering a few installation options on Linux and macOS. 

### Via Conda (best for Linux) <a name="required_packages_conda"></a>

A simple way to install most prerequisites is via [Conda](https://github.com/conda/conda) environments. Cocoa's internal scripts will then install any remaining missing packages when compiling the base code. Assuming that the user had previously installed [Minicoda](https://docs.conda.io/en/latest/miniconda.html) (or [Anaconda](https://www.anaconda.com/products/individual)), the first step is to type the following commands to create the cocoa Conda environment.

    $ conda create --name cocoa python=3.7 --quiet --yes

and

    $ conda install -n cocoa --quiet --yes  \
      'conda-forge::libgcc-ng=10.3.0' \
      'conda-forge::libstdcxx-ng=10.3.0' \
      'conda-forge::libgfortran-ng=10.3.0' \
      'conda-forge::gxx_linux-64=10.3.0' \
      'conda-forge::gcc_linux-64=10.3.0' \
      'conda-forge::gfortran_linux-64=10.3.0' \
      'conda-forge::openmpi=4.1.1' \
      'conda-forge::sysroot_linux-64=2.17' \
      'conda-forge::git=2.33.1' \
      'conda-forge::git-lfs=3.0.2' \
      'conda-forge::hdf5=1.10.6' \
      'conda-forge::git-lfs=3.0.2' \
      'conda-forge::cmake=3.21.3' \
      'conda-forge::boost=1.76.0' \
      'conda-forge::gsl=2.7' \
      'conda-forge::fftw=3.3.10' \
      'conda-forge::cfitsio=4.0.0' \
      'conda-forge::openblas=0.3.18' \
      'conda-forge::lapack=3.9.0' \
      'conda-forge::armadillo=10.7.3'\
      'conda-forge::expat=2.4.1' \
      'conda-forge::cython=0.29.24' \
      'conda-forge::numpy=1.21.4' \
      'conda-forge::scipy=1.7.2' \
      'conda-forge::pandas=1.3.4' \
      'conda-forge::mpi4py=3.1.2' \
      'conda-forge::matplotlib=3.5.0' \
      'conda-forge::astropy=4.3.1' 
 
With this installation method, users must activate the Conda environment whenever working with Cocoa, as shown below 

    $ conda activate cocoa

Users can now proceed to the section [Installation of cocoa base code](#cocoa_base_code). 

### Via Docker (best for MacOS/Windows) <a name="required_packages_docker"></a>

Docker installation will allow users to run Cocoa inside an instantiation of the [Whovian-Cosmo](https://hub.docker.com/r/vivianmiranda/whovian-cosmo) docker image. Installation of the [docker engine](https://docs.docker.com/engine/) on local PCs is a straightforward process, but it does require `sudo` privileges (see Docker's [official documentation](https://docs.docker.com/engine/install/) for OS-specific instructions).

  On macOS, type:

    $ docker run -it -p 8080:8888 -v $(pwd):/home/whovian/host/ -v ~/.ssh:/home/whovian/.ssh:ro vivianmiranda/whovian-cosmo:version-1.0.3

Linux users must type the following command instead:

    $ docker run -it -p 8080:8888 --user $(id -u):$(id -g) -v $(pwd):/home/whovian/host/ -v ~/.ssh:/home/whovian/.ssh:ro vivianmiranda/whovian-cosmo:version-1.0.3

The bind `-v ~/.ssh:/home/whovian/.ssh:ro` allows users to pull, push and clone GitHub repositories from inside the container. Users must invoke the command on the parent directory of the path where access inside the docker container is sought. 

When running the `docker run (...)/whovian-cosmo:version-1.0.3` for the first time, the docker engine will automatically download the corresponding image. This step may take some time, as the [Whovian-Cosmo](https://hub.docker.com/r/vivianmiranda/whovian-cosmo) image has approximately 700 Megabytes.

  The last step is to access the folder `/home/whovian/host/` where the host files have been mounted:

    $ cd /home/whovian/host/

and proceed to the section [Installation of cocoa base code](#cocoa_base_code). 

(**Warning**) There isn't permanent storage outside `/home/whovian/host/`. Be aware of this fact to not lose any work

(**Warning**) Most HPC systems don't allow users to run docker containers via the standard [docker engine](https://docs.docker.com/engine/) for [security reasons](https://www.reddit.com/r/docker/comments/7y2yp2/why_is_singularity_used_as_opposed_to_docker_in/?utm_source=share&utm_medium=web2x&context=3). There is, however, an alternative engine called [Singularity](https://sylabs.io/guides/3.6/user-guide/index.html) that is in compliance with most HPC requirements. The [Singularity](https://sylabs.io/guides/3.6/user-guide/index.html) engine installation requires administrative privileges, but many HPC enviroments have already adopted it.

To run docker images with Singularity, go to the folder you want to store the image and type:

    $ singularity build whovian-cosmo docker://vivianmiranda/whovian-cosmo

This command will download the [Whovian-Cosmo](https://hub.docker.com/r/vivianmiranda/whovian-cosmo) image and convert it to a format that can be understood by Singularity (this might take a few minutes). To run the container interactively, type:

    $ singularity shell --no-home --bind /path/to/cocoa:/home/whovian/host --bind ~/.ssh:/home/whovian/.ssh:ro whovian-cosmo

### Via Homebrew <a name="required_packages_homebrew"></a>

This subsection assumes users adopt [Homebrew](https://brew.sh) as the macOS package manager and [BASH](https://www.howtogeek.com/444596/how-to-change-the-default-shell-to-bash-in-macos-catalina/) as the default shell. 

(**Warning**) The installation procedures for MacOS are always working in progress; there is no guarantee that our recommendations will work without additional tweaks given how fast [Homebrew](https://brew.sh) updates its packages and Apple changes the OS. *The Docker installation is preferred*.

(**Warning**) Currently, there is a regression in GCC-11, the default compiler on Homebrew, preventing us from adopting it. Therefore, users must build most brew numerical packages from their source via the command `--build-from-source`. Building from source is painfully slow and error prone and does not work in Apple M1 chips without Rosetta 2 emulator. We can do nothing but wait for the GCC developers to fix the bug. 

Below is a list of Homebrew and pip commands that install most dependencies

    brew install gcc@10
    alias brew='HOMEBREW_CC=gcc-10 HOMEBREW_CXX=g++-10 brew'
    brew install open-mpi --build-from-source
    brew install git
    brew install git-lfs
    git lfs install
    git lfs install --system
    brew install gmp
    brew install hdf5 --build-from-source
    brew install python@3.7 || (brew upgrade python@3.7 && brew cleanup python@3.7)
    brew install cmake
    brew install armadillo --build-from-source
    brew install boost --build-from-source
    brew install gsl --build-from-source
    brew install fftw --build-from-source
    brew install cfitsio --build-from-source
    brew install lapack --build-from-source
    brew install findutils
    brew install xz
    export PATH="/usr/local/opt/python@3.7/bin:$PATH"

    pip3.7 install virtualenv --upgrade
    pip3.7 install wget --upgrade
    pip3.7 install wheel --upgrade
    pip3.7 install setuptools --upgrade
    pip3.7 install six --upgrade
    pip3.7 install python-dateutil --upgrade
    pip3.7 install pytz --upgrade
    pip3.7 install mpmath --upgrade
    pip3.7 install PyYAML --upgrade
    pip3.7 install pyparsing --upgrade
    pip3.7 install fuzzywuzzy --upgrade
    pip3.7 install cycler --upgrade
    pip3.7 install kiwisolver --upgrade
    pip3.7 install enum34 --upgrade
    pip3.7 install pillow --upgrade
    pip3.7 install numpy --upgrade
    pip3.7 install scipy --upgrade
    pip3.7 install sympy --upgrade
    pip3.7 install cython  --upgrade
    pip3.7 install imageio --upgrade
    pip3.7 install pillow --upgrade
    pip3.7 install pandas --upgrade
    pip3.7 install py-bobyqa --upgrade
    pip3.7 install matplotlib --upgrade
    pip3.7 install pybind11 --upgrade
    pip3.7 install GetDist --upgrade
    pip3.7 install astropy --upgrade
    pip3.7 install pyfits --upgrade

Users must also update their `$PATH` so the OS can detect the [Homebrew](https://brew.sh) python installation, as shown below

    python3=/Users/XXX/Library/Python/3.7/bin
    export PATH=$python3:$PATH

Users can now proceed to the section [Installation of cocoa base code](#cocoa_base_code). 

### Via Cocoa's internal cache <a name="required_packages_cache"></a>

(**Warning**) This method is painfully slow, not advisable. It does, however, offer users that don't work with Conda the opportunity to encapsulate the installation of required packages. *We advise the adoption of Conda or Docker installation*. 

Whenever Conda or Docker installation procedures are unavailable, the user can still perform a local semi-autonomous installation on Linux based on a few scripts we implemented. We also provide a local copy of most required packages on Cocoa's cache folder [cocoa_installation_libraries](https://github.com/CosmoLike/cocoa/tree/main/cocoa_installation_libraries), as there are HPC machines where compute nodes that compile code don't have internet access (NASA Pleiades being one example). We, therefore, only assume the pre-installation of the following packages to perform the local setup:

   - [Bash](https://www.amazon.com/dp/B0043GXMSY/ref=cm_sw_em_r_mt_dp_x3UoFbDXSXRBT);
   - [Git](https://git-scm.com) v1.8+;
   - [Git LFS](https://git-lfs.github.com);
   - [gcc](https://gcc.gnu.org) v10.*;
   - [gfortran](https://gcc.gnu.org) v10.*;
   - [g++](https://gcc.gnu.org) v10.*;
   - [Python](https://www.python.org) v3.7.*;
   - [PIP package manager](https://pip.pypa.io/en/stable/installing/)
   - [Python Virtual Environment](https://www.geeksforgeeks.org/python-virtual-environment/)

To perform the local semi-autonomous installation, users should follow the procedures on section [Installation of cocoa base code](https://github.com/CosmoLike/cocoa#installation-of-cocoa-base-code), adding, however, the many additional configurations on [set_installation_options](https://github.com/CosmoLike/cocoa/blob/master/Cocoa/set_installation_options) script that are explained below.

The local installation via cocoa's internal cache is selected whenever the environmental key `MANUAL_INSTALLATION` is set:

    [Extracted from set_installation_options script] 
    
    #  ---------------------------------------------------------------------------
    # HOW COCOA BE INSTALLED? -------------------------------

    #export DOCKER_INSTALLATION=1
    #export MINICONDA_INSTALLATION=1
    #export MACOS_HOMEBREW_INSTALLATION=1
    export MANUAL_INSTALLATION=1
    
Users also need to set the following self-explanatory environmental keys on [set_installation_options](https://github.com/CosmoLike/cocoa/blob/master/Cocoa/set_installation_options):
 
    [Extracted from set_installation_options script]
  
    elif [ -n "${MANUAL_INSTALLATION}" ]; then

      export GLOBAL_PACKAGES_LOCATION=/usr/local
      export PYTHON_VERSION=3
      export FORTRAN_COMPILER=gfortran
    
      export C_COMPILER=gcc
      export CXX_COMPILER=g++
      export GLOBALPYTHON3=python3
      export MPI_FORTRAN_COMPILER=mpif90
      export MPI_CXX_COMPILER=mpicc
      export MPI_CC_COMPILER=mpicxx
    
      # In case global packages are available 
      #export IGNORE_DISTUTILS_INSTALLATION=1
      #export IGNORE_OPENBLAS_INSTALLATION=1
      #export IGNORE_XZ_INSTALLATION=1
      #export IGNORE_ALL_PIP_INSTALLATION=1
      #export IGNORE_CMAKE_INSTALLATION=1
      #export IGNORE_CPP_BOOST_INSTALLATION=1
      #export IGNORE_CPP_ARMA_INSTALLATION=1
      #export IGNORE_CPP_SPDLOG_INSTALLATION=1
      #export IGNORE_C_GSL_INSTALLATION=1
      #export IGNORE_C_CFITSIO_INSTALLATION=1
      #export IGNORE_C_FFTW_INSTALLATION=1 
      #export IGNORE_OPENBLAS_INSTALLATION=1
      #export IGNORE_FORTRAN_LAPACK_INSTALLATION=1
   
(**Warning**) Our scripts never install packages on `$HOME/.local`. Doing so could impose incompatibilities between Cobaya and different projects (or break the user's environment for other projects). All requirements for Cocoa are installed at

    Cocoa/.local/bin
    Cocoa/.local/include
    Cocoa/.local/lib
    Cocoa/.local/share

Users can now proceed to the section [Installation of cocoa base code](#cocoa_base_code). 

## Installation of cocoa base code <a name="cocoa_base_code"></a>

Type:

    $ git lfs clone git@github.com:CosmoLike/cocoa.git

to clone the repository. 

(**Warning**) We have a limited monthly quota in bandwidth for [Git LFS](https://git-lfs.github.com) files, and therefore we ask users to use good judgment in the number of times they clone repositories with large files (each clone will download around 5GB from Git LFS).

Cocoa chooses the preferred method of installation and compilation via special environment keys located on [set_installation_options](https://github.com/CosmoLike/cocoa/blob/master/Cocoa/set_installation_options) script. The key should reflect the user's choice in section [Installation of cocoa required packages](https://github.com/CosmoLike/cocoa#installation-of-cocoa-required-packages) to the installation method of cocoa required packages.

    [Extracted from set_installation_options script]
    #  ---------------------------------------------------------------------------
    # HOW COCOA BE INSTALLED? -------------------------------

    #export DOCKER_INSTALLATION=1
    #export MINICONDA_INSTALLATION=1
    #export MACOS_HOMEBREW_INSTALLATION=1
    #export MANUAL_INSTALLATION=1
    
Users must then type the following commands

    $ source setup_cocoa_installation_packages
    
    $ source compile_external_modules
    
`setup_cocoa_installation_packages` script will decompress the data files and install all packages that may have been left out in the Conda/Docker/Homebrew/Manual installation. File decompression should only take a few minutes, while package installation time may range from a few minutes (if installation via *Conda*, *Docker* or *Homebrew* was selected in the previous section) to more than one hour (if installation *via Cocoa's internal scripts and cache* was selected in the previous section).

(**warning**) The script [compile_external_modules](https://github.com/CosmoLike/cocoa/blob/master/Cocoa/compile_external_modules) will try to compile Camb, Planck, Class, and Polychord samplers. If users want to compile/recompile just one of these packages for any reason, we provide scripts for that as well, see appendix [Compiling Boltzmann, CosmoLike and Likelihood codes separatelly](#compile_separatelly) section.

### Running Examples <a name="cocoa_base_code_examples"></a>

After that, the Cobaya Framework should be ready, and the user can test a few examples following the commands below

    $ source start_cocoa
    
    $(.local) mpirun -n 1 cobaya-run ./projects/example/EXAMPLE_EVALUATE[1-4].yaml -f
    
    $(.local) mpirun -n 4 cobaya-run ./projects/example/EXAMPLE_MCMC1.yaml -f
    
These examples will evaluate various likelihoods at specific cosmologies. The `-f` ensures that the same YAML file can be run multiple times, overwriting output files from previous evaluations that are located at `./chains`.

(**warning**) In HPC environments, you may face this error when running `cobaya` within the Conda cocoa environment. 
    
    --------------------------------------------------------------------------
    It looks like MPI_INIT failed for some reason; your parallel process is
    likely to abort.  There are many reasons that a parallel process can
    fail during MPI_INIT; some of which are due to configuration or environment
    problems.  This failure appears to be an internal failure; here's some
    additional information (which may only be relevant to an Open MPI
    developer):

    mca_bml_base_open() failed
    --> Returned "Not found" (-13) instead of "Success" (0)
    --------------------------------------------------------------------------

The origin of this problem lies in the choice of conda-forge developers to not compile openmpi with Infiniband compatibility. The solution is to replace `mpirun -n XXX` with `mpirun -n XXX --mca btl tcp,self` which will enforce the TCP protocol for communication. Users outraged by the small overhead TCP will bring over Infiniband can perform the [installation via Cocoa's internal cache](required_packages_cache) that depends on the HPC module system to load the openmpi compiled by the system administrators.    

Once the work on the Cocoa environment is done, type:

    $(.local) source stop_cocoa

The script [start_cocoa](https://github.com/CosmoLike/cocoa/blob/master/Cocoa/start_cocoa) ensures that `LD_LIBRARY_PATH`, `CPATH`, `C_INCLUDE_PATH`, `CPLUS_INCLUDE_PATH` and `PATH` will point preferably to the local Cocoa installation. It also initializes the `$(.local)` private python environment (including `PYTHONPATH`). The script [stop_cocoa](https://github.com/CosmoLike/cocoa/blob/master/Cocoa/stop_cocoa), on the other hand, guarantees that the bash environment is clean after the work on Cocoa is completed. 

Users can, therefore, install multiple Cocoa instances. While one instantiation of Cocoa is running chains, for example, users can tweak Cosmolike in a second instance without affecting these MCMCs.

(**Warning**) Users that opt for the Conda installation will have a terminal that looks like this:

    $(Cocoa)(.local)

*This is a feature, not a bug*! The Conda environment can be the same for all Cocoa instances (projects rarely need local versions of GSL/FFTW/Boost...), with [start_cocoa](https://github.com/CosmoLike/cocoa/blob/master/Cocoa/start_cocoa)/[stop_cocoa](https://github.com/CosmoLike/cocoa/blob/master/Cocoa/stop_cocoa) loading/unloading the corresponding `(.local)` private python environments (with corresponding local Planck/CAMB/Cosmolike/CLASS/Cobaya...).

# Download/compiling/running specific cocoa projects <a name="cocoa_projects"></a> 

The `projects` folder was designed to include all the projects that are being developed by our group. Individual projects must be hosted on independent folders named `cocoa_XXX` where XXX is the project name. 

The majority of projects we are working on are not public (yet), and they are safeguarded on the private repositories listed on `project/clone_all.sh` (the backbone Cosmolike software, however, is publicly available at `external_modules/code`!). 

You can add your projects there, and the script `setup_cocoa_installation_packages` will try to clone all listed projects. Having inaccessible repositories listed at `project/clone_all.sh` will not cause any errors.

The `cocoa_XXX` folder that host the `XXX` project needs to have the more or less the following structure (taken from our private DES-Y3 project)

    +-- cocoa_des_y3
    |    +-- likelihood
    |    |   +-- _cosmolike_prototype_base.py
    |    |   +-- des_3x2pt.py
    |    |   +-- des_3x2pt.yaml
    |    |   +-- des_2x2pt.py
    |    |   +-- des_3x2pt.yaml
    |    |   +-- des_cosmic_shear.py
    |    |   +-- des_cosmic_shear.yaml
    |    +-- scripts
    |    |   +-- compile_des_y3
    |    |   +-- start_des_y3
    |    |   +-- stop_des_y3
    |    +-- data
    |    |   +-- DES.paramnames
    |    |   +-- DES_Y3.dataset
    |    |   +-- datavector.txt
    |    |   +-- covariance.txt
    |    |   +-- nzlens.txt
    |    |   +-- nzsource.txt
    |    |   +-- mask.mask
    |    +-- interface
    |    |   +-- MakefileCosmolike
    |    |   +-- cosmolike_des_y3_interface.py
    |    |   +-- interface.cpp
    |    |   +-- interface.hpp
    |    +-- chains
    |    |   +-- README

(**Warning**) Developers with access to UofA projects can download them by typing

    $(.local) source $ROOTDIR/projects/clone_all.sh

Here is a list of projects inside [clone_all](https://github.com/CosmoLike/cocoa/blob/main/Cocoa/projects/.clean_all.sh) script

    [Extracted from clone_all script]
    URLS="git@github.com:CosmoLike/cocoa_des_y1.git
      git@github.com:CosmoLike/cocoa_des_y3.git
      git@github.com:CosmoLike/cocoa_desxplanck.git
      git@github.com:CosmoLike/cocoa_desy1xplanck.git
      git@github.com:CosmoLike/cocoa_lsst_fourier.git
      git@github.com:CosmoLike/cocoa_kids.git
      git@github.com:CosmoLike/cocoa_des_4x2ptN.git"

## Appendix <a name="appendix"></a>

### Compiling Boltzmann, CosmoLike and Likelihood codes separatelly <a name="appendix_compile_separatelly"></a>

To avoid excessive compilation times during development, users can use specialized scripts that compile only the specific modules:

    $(.local) source ./installation_scripts/setup_class

    $(.local) source ./installation_scripts/setup_camb

    $(.local) source ./installation_scripts/setup_planck

    $(.local) source ./installation_scripts/setup_polychord

### Running Jupyter Notebooks inside the [Whovian-Cosmo](https://hub.docker.com/r/vivianmiranda/whovian-cosmo) docker container <a name="appendix_jupyter_whovian"></a>

[Cobaya](https://github.com/CobayaSampler), the framework that Cocoa heavily depends on has excellent integration with Jupyter notebooks. Below, some in-depth instructions to run notebooks inside the [Whovian-Cosmo](https://hub.docker.com/r/vivianmiranda/whovian-cosmo) docker container

To start a jupyter notebook, type the following command on the docker container:

    $ jupyter notebook --no-browser

The terminal will show a message similar to the following template:

    [... NotebookApp] Writing notebook server cookie secret to /home/whovian/.local/share/jupyter/runtime/notebook_cookie_secret
    [... NotebookApp] WARNING: The notebook server is listening on all IP addresses and not using encryption. This is not recommended.
    [... NotebookApp] Serving notebooks from local directory: /home/whovian/host
    [... NotebookApp] Jupyter Notebook 6.1.1 is running at:
    [... NotebookApp] http://f0a13949f6b5:8888/?token=XXX
    [... NotebookApp] or http://127.0.0.1:8888/?token=XXX
    [... NotebookApp] Use Control-C to stop this server and shut down all kernels (twice to skip confirmation).

where `XXX` in the line `[... NotebookApp] or http://127.0.0.1:8888/?token=XXX` is the token you need to save to access the notebook. We will assume you are running the docker container in a host server with URL `your_sever.com` that you are accessing via ssh. From your local PC type:

    $ ssh your_username@your_sever.com -L 8080:localhost:8080

   Finally, go to your browser and type `http://localhost:8080/?token=XXX`, where `XXX` is the previously saved token. For security, we do not allow password-based connections to the jupyter notebooks.
   
### Summary Information about Cocoa's configuration files <a name="appendix_config_files"></a>

The installation of Cocoa required packages, as well as Boltzmann and Likelihood codes, are managed via the following scripts located at `./Cocoa`.

 - [set_installation_options](https://github.com/CosmoLike/cocoa/blob/master/Cocoa/set_installation_options)

    This file contains environment variables that manage the installation process.

 - [setup_cocoa_installation_packages](https://github.com/CosmoLike/cocoa/blob/master/Cocoa/setup_cocoa_installation_packages)

    This file has instructions on how to install packages required by the Cocoa Framework.

 - [compile_external_modules](https://github.com/CosmoLike/cocoa/blob/master/Cocoa/compile_external_modules)

    This file has instructions on how to compile Boltzmann, Sampler and likelihood codes. 

 - [start_cocoa](https://github.com/CosmoLike/cocoa/blob/master/Cocoa/start_cocoa)

    This file has instructions on how to set up the Python virtual environment.

 - [stop_cocoa](https://github.com/CosmoLike/cocoa/blob/master/Cocoa/stop_cocoa)

    This file has instructions on how to unset the Python virtual environment - including recovering original `PYTHONPATH`, `LD_LIBRARY_PATH`, and `PATH`. 

 - [clean_all](https://github.com/CosmoLike/cocoa/blob/master/Cocoa/clean_all)

    This file has instructions on how to clean keys associated with the Python virtual environment and delete the compilation of the Boltzmann, Sampler, and likelihood codes, and local installation of the required packages installed by the [setup_cocoa_installation_packages].
