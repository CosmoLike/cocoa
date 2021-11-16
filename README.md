# Cocoa - The [Cobaya](https://github.com/CobayaSampler)-[CosmoLike](https://github.com/CosmoLike) Joint Architecture

## Overview

  Cocoa allows users to run, inside the [Cobaya](https://github.com/CobayaSampler) framework, [CosmoLike](https://github.com/CosmoLike) routines that 1) analyze data primarily from the [Dark Energy Survey](https://www.darkenergysurvey.org) (a.k.a DES), and 2) simulate future multi-probe analyses, e.g. for Rubin Observatory's Legacy Survey of Space and Time or the Roman Space Telescope. This readme file presents basic and advanced instructions for installing all Cocoa components, including the [Planck likelihood](https://wiki.cosmos.esa.int/planck-legacy-archive/index.php/Main_Page). **By no means, we want to discourage general users from cloning Cobaya, CAMB, CLASS, Polychord, and Planck data from their original repositories. Please check the appendix [Proper Credits](https://github.com/CosmoLike/cocoa#proper-credits)**. 

Cocoa developers explicitly copied all Boltzman codes, likelihoods, Cocoa Sampler, and all python/C/C++/Fortran requirements to create a very controllable environment that will ensure the reproducibility of our DES-Y3 and LSST-Y1 results. All the required python packages, including Numpy and Scipy, have been stored in [cocoa_installation_libraries](https://github.com/CosmoLike/cocoa/tree/main/cocoa_installation_libraries) folder (linux installation). Installation scripts enforce the use of the available local cache instead of internet connection or packages already installed on the host machine. 

PS: We are quite flexible on the installation procedures for MacOS (see appendix: [Prerequisites for MacOS](https://github.com/CosmoLike/cocoa#prerequisites-for-macos-system-installation)) because our group rarely runs production-ready code in such an environment.

# Installation

## Before Cloning the Repository

There are two options for installing the requirements to run Cocoa and clone the repository using [Git LFS](https://git-lfs.github.com). The first and more straightforward one is through Docker containers. The latter option is the manual installation of the requirements directly into the Linux/macOS (Windows is not supported).

If you want to use cocoa on your Mac OS system, we recommend using the Docker installation, see appendices [Docker Installation Part I: Further Information for PCs Environment](https://github.com/CosmoLike/cocoa/blob/master/README.md#docker-installation-part-i-further-information-for-pcs-environment) 

If you want to use Docker on HPCs please go here [Docker Installation Part II: Further Information for HPC Systems](https://github.com/CosmoLike/cocoa/blob/master/README.md#docker-installation-part-ii-further-information-for-hpc-systems).

Below we detail non-Docker, system specific installation. Please see appendices for prerequisites: [Prerequisites for MacOS (System Installation)](https://github.com/CosmoLike/cocoa#prerequisites-for-macos-system-installation) and [Prerequisites for Linux (System Installation)](https://github.com/CosmoLike/cocoa#prerequisites-for-linux-system-installation).

## OPTIONAL: MINICONDA SETUP
A simple way to install all prerequisites (plus most packages needed by Cocoa) is via [miniconda](https://docs.conda.io/en/latest/miniconda.html). This is almost as simple as using Docker

    conda create --name cocoa python=3.7 --quiet --yes
    
    conda install -n cocoa --quiet --yes  \
      'conda-forge::libgcc-ng=11.2.0' \
      'conda-forge::libstdcxx-ng=11.2.0' \
      'conda-forge::libgfortran-ng=11.2.0' \
      'conda-forge::gxx_linux-64=11.2.0' \
      'conda-forge::gcc_linux-64=11.2.0' \
      'conda-forge::gfortran_linux-64=11.2.0' \
      'conda-forge::sysroot_linux-64=2.17' \
      'conda-forge::git=2.33.1' \
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
      'conda-forge::pandas=1.3.4'
      
    source activate cocoa
        
## Cloning the Repository

Type

    $ git lfs clone https://github.com/CosmoLike/cocoa.git

(**Warning**) We have a monthly quota of only 150 GB in bandwidth for [Git LFS](https://git-lfs.github.com) files, and therefore we ask users to use good judgment in the number of times you clone repositories with large files (each clone will download around 5GB from Git LFS). See appendix [Git LFS](https://github.com/CosmoLike/cocoa/blob/master/README.md#git-lfs) for further information.

## Installation/Compilation of Numerical libraries

Cocoa chooses the preferred method of installation via special environment keys located on [set_installation_options](https://github.com/CosmoLike/cocoa/blob/master/Cocoa/set_installation_options) script (Docker is the default option)

    # ----------------------------------------------------------------------------
    # ----------------------------------------------------------------------------
    # -----------------------  OPTIONS FOR SYSTEM INSTALLATION -------------------
    # ----------------------------------------------------------------------------
    # ----------------------------------------------------------------------------

    export DOCKER_INSTALLATION=1
    #export MINICONDA_INSTALLATION=1
    (...)
    #export AMYPOND_SUPERCOMPUTER_INSTALLATION=1
    #export MACOS_HOMEBREW_INSTALLATION=1
    (...)
    
For a few HPC systems commonly used by the Cocoa developers, special keys with pre-defined configurations exist. See section [System Installation: Further Information for Linux]() for further instructions on `MANUAL_INSTALLATION` if the pre-defined configuration is not available for your environment. We do, however, strongly encourage the miniconda installation (way easier). For macOS system installation, the user needs to uncomment the line `export MACOS_HOMEBREW_INSTALLATION=1`, while making sure that all other special keys are unset


**Exporting, at the same, more than one of the special keys listed above produces undefined behavior!**

PS: On HPC systems, dealing with docker containers and [PBS](https://www.openpbs.org) and  [SLURM](https://slurm.schedmd.com/documentation.html) submission scripts can be a little annoying. Therefore, users are encouraged to learn how to do a manual installation on your HPC system when using Cocoa for running lots of MCMC chains.

## Sourcing the Configuration Files

The installation can be performed with the following commands:

    $ source setup_cocoa_installation_packages

**WARNING**: before running `source setup_cocoa_installation_packages`, make sure that you added the projects you have access (and want to work on) via ssh-keys on `projects/clone_all.sh` file. See [Project-Folder](https://github.com/CosmoLike/cocoa#the-projects-folder) section for further details.

    $ source compile_external_modules
    
    $ source start_cocoa

Sourcing [setup_cocoa_installation_packages](https://github.com/CosmoLike/cocoa/blob/master/Cocoa/setup_cocoa_installation_packages) may require a long time (~hours) depending on the supercomputer environment (System Option). It is possible to speed-up installation considerably by running `source setup_cocoa_installation_packages` on an interactive node w/ ~16 threads. For the Docker Installation, however, sourcing [setup_cocoa_installation_packages](https://github.com/CosmoLike/cocoa/blob/master/Cocoa/setup_cocoa_installation_packages) should take only a few minutes.

Once the work on the Cocoa environment is done, type:

    $ source stop_cocoa

PS: Users are highly encouraged to read appendix [Further Information about the configuration files](https://github.com/CosmoLike/cocoa/blob/master/README.md#further-information-about-the-configuration-files) to learn more about the configuration files.

## Running Examples

After that, the Cobaya Framework should be ready, and the user can test a few examples:

    mpirun -n 1 cobaya-run ./projects/example/EXAMPLE_EVALUATE[1-4].yaml -f
    mpirun -n 4 cobaya-run ./projects/example/EXAMPLE_MCMC1.yaml -f
    
These examples will evaluate various likelihoods at specific cosmologies. The `-f` ensures that the same YAML file can be run multiple times, overwriting output files from previous evaluations that are located at `./chains`.

# The Projects Folder

The `projects` folder includes all the projects that are being developed by our group. Individual projects must be hosted on independent folders named `cocoa_XXX` where XXX is the project name. The majority of projects we are working on are not public (yet), and they are safeguarded on the private repositories listed on `project/clone_all.sh` (the backbone Cosmolike software, however, is publicly available at `external_modules/code`!). You can add your projects there, and the script `setup_cocoa_installation_packages` will try to clone all listed projects. **Having inaccessible repositories listed at `project/clone_all.sh` will not cause any errors**. By default, cocoa will try to clone the private projects `cocoa_des_y3`, `cocoa_desxplanck` and `cocoa_des_y3u`.

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
      
# Appendix

## Prerequisites for Linux (System Installation)

We assume the user has the following packages installed:

   - [Bash](https://www.amazon.com/dp/B0043GXMSY/ref=cm_sw_em_r_mt_dp_x3UoFbDXSXRBT);
   - [Git](https://git-scm.com) v1.8+;
   - [Git LFS](https://git-lfs.github.com);
   - [gcc](https://gcc.gnu.org) v9.*;
   - [gfortran](https://gcc.gnu.org) v9.*;
   - [g++](https://gcc.gnu.org) v9.*;
   - [Python](https://www.python.org) v3.7.*;
   - [PIP package manager](https://pip.pypa.io/en/stable/installing/)
   - [Python Virtual Environment](https://www.geeksforgeeks.org/python-virtual-environment/)

We also assume all contributors have added their [ssh-keys](https://docs.github.com/en/enterprise/2.15/user/articles/adding-a-new-ssh-key-to-your-github-account) to Cocoa's private submodules. Note that

   - **zsh (the Z shell) is not a valid substitution for bash**,



## Prerequisites for MacOS (System Installation)

We assume the user adopts [Homebrew](https://brew.sh) as the OS package manager, and [BASH](https://www.howtogeek.com/444596/how-to-change-the-default-shell-to-bash-in-macos-catalina/) as the default shell. While alternatives, such as [MacPorts](https://www.macports.org) and [Fink](http://www.finkproject.org), may provide all the requirements, we don't  present instructions for them. We also don't provide documentation for installation on [zsh](https://www.howtogeek.com/444596/how-to-change-the-default-shell-to-bash-in-macos-catalina/) shell, the default shell on macOS Catalina. Finally, we also assume all contributors have added their ssh-keys to Cocoa's private submodules.

Here is a list of Homebrew and pip commands that install most dependencies

    brew install gcc@9
    alias brew='HOMEBREW_CC=gcc-9 HOMEBREW_CXX=g++-9 brew'
    brew install open-mpi || (brew upgrade open-mpi && brew cleanup open-mpi)
    brew install git || (brew upgrade git && brew cleanup git)
    brew install git-lfs || (brew upgrade git-lfs && brew cleanup git-lfs)
    git lfs install
    git lfs install --system
    brew install gmp --build-from-source
    brew install hdf5 --build-from-source
    brew install python@3.7 || (brew upgrade python@3.7 && brew cleanup python@3.7)
    brew install cmake || (brew upgrade cmake && brew cleanup cmake)
    brew install armadillo --build-from-source
    brew install boost --build-from-source
    brew install gsl || (brew upgrade gsl && brew cleanup gsl)
    brew install fftw || (brew upgrade fftw && brew cleanup fftw)
    brew install cfitsio || (brew upgrade cfitsio && brew cleanup cfitsio)
    brew install lapack || (brew upgrade lapack && brew cleanup lapack)
    brew install findutils || (brew upgrade findutils && brew cleanup findutils)
    brew install xz || (brew upgrade xz && brew cleanup xz)
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

**It is also necessary to update `$PATH`** so the OS can detect the [Homebrew](https://brew.sh) python installation. Adding to `~/.bash_profile` code similar to the lines below should be enough (don't forget to replace `XXX` by the user's home folder):

    python3=/Users/XXX/Library/Python/3.7/bin
    export PATH=$python3:$PATH
    
## Docker Installation Part I: Further Information for PCs Environment

We highly advise users to run Cocoa inside an instantiation of the [Whovian-Cosmo](https://hub.docker.com/r/vivianmiranda/whovian-cosmo) docker image. Installation of the [docker engine](https://docs.docker.com/engine/) on local PCs is a reasonably straightforward process - see the [official documentation](https://docs.docker.com/engine/install/) for OS-specific instructions.

  The only step left is the command on macOS:

    $ docker run -it -p 8080:8888 -v $(pwd):/home/whovian/host/ -v ~/.ssh:/home/whovian/.ssh:ro vivianmiranda/whovian-cosmo:version-1.0.2

For linux users, use the following command instead:

    $ docker run -it -p 8080:8888 --user $(id -u):$(id -g) -v $(pwd):/home/whovian/host/ -v ~/.ssh:/home/whovian/.ssh:ro vivianmiranda/whovian-cosmo:version-1.0.2

The bind of the `~/.ssh` folder in the commands above allows users to commit and clone from inside the container. Both commands need to be invoked on the parent folder so the docker container will have access to the listed folders on the host OS. When running the command `docker run` on a specific container for the first time, the docker engine will automatically download the corresponding [Docker Hub](https://hub.docker.com/) image.
This step may take some time, as the [Whovian-Cosmo](https://hub.docker.com/r/vivianmiranda/whovian-cosmo) image has approximately 700 Megabytes.

  The last step is to access the folder `/home/whovian/host/` where the host files have been mounted:

    $ cd /home/whovian/host/

and proceed to the section [Cloning the Repository](https://github.com/CosmoLike/cocoa#cloning-the-repository). **There isn't permanent storage outside `/home/whovian/host/`. Be aware of this fact to not lose any work**.

PS: If users want to use Jupiter notebooks from the container at any time during development (the [Cobaya](https://github.com/CobayaSampler) framework that Cocoa heavily depends on has excellent integration with Jupyter notebooks), then we advise them to read appendix [Running Jupyter Notebooks inside the Whovian-Cosmo container](https://github.com/CosmoLike/cocoa/blob/master/README.md#running-jupyter-notebooks-inside-the-whovian-cosmo-container) for in-depth instructions.


## Docker Installation Part II: Further Information for HPC Systems

HPC systems don't allow users to run docker containers using the standard [docker engine](https://docs.docker.com/engine/) for [security reasons](https://www.reddit.com/r/docker/comments/7y2yp2/why_is_singularity_used_as_opposed_to_docker_in/?utm_source=share&utm_medium=web2x&context=3). There is, however, an alternative engine called [Singularity](https://sylabs.io/guides/3.6/user-guide/index.html) that can run docker images in compliance with HPC security requirements. The [Singularity](https://sylabs.io/guides/3.6/user-guide/index.html) engine installation requires administrative privileges. However, many HPC systems have already adopted it.

To run docker images with Singularity, go to the folder you want to store the image and type:

    singularity build whovian-cosmo docker://vivianmiranda/whovian-cosmo

This command will download the [Whovian-Cosmo](https://hub.docker.com/r/vivianmiranda/whovian-cosmo) image and convert it to an image format that can be accessed via Singularity (Singularity takes a few minutes to do such conversion). To run the container interactively, type:

    singularity shell --no-home --bind /path/to/cocoa:/home/whovian/host --bind ~/.ssh:/home/whovian/.ssh:ro whovian-cosmo

after requesting interactive nodes successfully (never run jobs on login nodes!). The bind of the `~/.ssh` folder in the command above allows users to commit and clone from inside the container. The last step is to access the folder `/home/whovian/host/` where the host files have been mounted

    $ cd /home/whovian/host/

and proceed to the section [Cloning the Repository](https://github.com/CosmoLike/cocoa#cloning-the-repository). **There isn't permanent storage outside `/home/whovian/host/`. Be aware of this fact to not lose any work**

PS: If users want to use Jupiter notebooks from the container at any time during development (the [Cobaya](https://github.com/CobayaSampler) framework that Cocoa heavily depends on has excellent integration with Jupyter notebooks), then we advise them to read appendix [Running Jupyter Notebooks inside the Whovian-Cosmo container](https://github.com/CosmoLike/cocoa/blob/master/README.md#running-jupyter-notebooks-inside-the-whovian-cosmo-container) for in-depth instructions.

## System Installation: Further Information for Linux

If there is no preset special key for your particular enviroment, the user need to perform a manual installation (and select the `export MANUAL_INSTALLATION = 1`). Manually installing the required packages by advanced users has the additional advantage of avoiding repeated compilation of prerequisites packages that take a long time to build. Our scripts never install packages on `$HOME/.local` or other locations in the user's `$PATH` and `$LD_LIBRARY_PATH`. Doing so could impose incompatibilities between Cobaya and different projects. The decision to perform global upgrades on packages required by multiple projects is the user's sole responsibility. All requirements for Cocoa are installed at

    Cocoa/.local/bin
    Cocoa/.local/include
    Cocoa/.local/lib
    Cocoa/.local/share

The script `compile_external_modules` contains the lines:

    if [ -n "${DONT_USE_SYSTEM_PIP_PACKAGES}" ]; then
        $GLOBALPYTHON3 -m venv $ROOTDIR/.local/
    else
        $GLOBALPYTHON3 -m venv $ROOTDIR/.local/ --system-site-packages
    fi

The command `venv` activates an [isolated python environment](https://python-guide-kr.readthedocs.io/ko/latest/dev/virtualenvs.html) that utilizes system packages if available (and if the users allows it via the enviroment flag `DONT_USE_SYSTEM_PIP_PACKAGES`)  but does not affect them in any way in case upgrades are required. Finally, consistent use of the scripts [compile_external_modules](https://github.com/CosmoLike/cocoa/blob/master/Cocoa/compile_external_modules) and [stop_cocoa](https://github.com/CosmoLike/cocoa/blob/master/Cocoa/stop_cocoa) ensures that `PYTHONPATH`, `LD_LIBRARY_PATH`, and `PATH` are not affected once the user decides to switch projects, avoiding the potential use of Cobaya's required packages elsewhere.

Once the manual installation of all required packages is performed, the environment keys

    export IGNORE_ALL_PIP_INSTALLATION=1

    export IGNORE_XZ_INSTALLATION=1

    export IGNORE_CMAKE_INSTALLATION=1

    export IGNORE_CPP_INSTALLATION=1

    export IGNORE_C_INSTALLATION=1

    export IGNORE_FORTRAN_INSTALLATION=1

ensure that the script [setup_cocoa_installation_packages](https://github.com/CosmoLike/cocoa/blob/master/Cocoa/setup_cocoa_installation_packages) skips compilation of any requirements as long as the keys `[*]_SUPERCOMPUTER_INSTALLATION` are not set as they erase manual configurations. This is the correct behavior as these keys are intended to manage the installation for beginners. If all those keys are set, [setup_cocoa_installation_packages](https://github.com/CosmoLike/cocoa/blob/master/Cocoa/setup_cocoa_installation_packages) will only uncompress [xz files](https://tukaani.org/xz/format.html) on `external_modules/code` and `external_modules/data` as long as the following environment variables
`NO_UNXZ_EXTERNAL_MODULES_CODE` and `NO_UNXZ_EXTERNAL_MODULES_DATA` are not set. Finally, a more fine-tuned configuration can be arranged instead of the global `IGNORE_[C|CPP]_INSTALLATION` keys with the following flags

    export IGNORE_CPP_BOOST_INSTALLATION=1

    export IGNORE_CPP_ARMA_INSTALLATION=1

    export IGNORE_CPP_SPDLOG_INSTALLATION=1

    export IGNORE_C_CFITSIO_INSTALLATION=1

    export IGNORE_C_FFTW_INSTALLATION=1

    export IGNORE_C_GSL_INSTALLATION=1

Users that perform manual configurations also need to tell the script the compiler to be used, as shown below:

     export FORTRAN_COMPILER=gfortran

     export C_COMPILER=gcc

     export CXX_COMPILER=g++

If the system-wide CMAKE is old, but the user has a local [CMAKE](https://cmake.org) v3.11+, the keys

    export CMAKE_ROOT=local_installation

    export CMAKE=local_installation

can be set to point our scripts to the local [CMAKE](https://cmake.org) v3.11+ installation.

## Further Information about the configuration files

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
    
## Proper Credits
   The following is not an exhaustive list of the codes we use

- [Cobaya](https://github.com/CobayaSampler) is a framework developed by Dr. Jesus Torrado and Prof. Anthony Lewis

- [Cosmolike](https://github.com/CosmoLike) is a framework developed by Prof. Elisabeth Krause and Prof. Tim Eifler

- [CAMB](https://github.com/cmbant/CAMB) is a Boltzmann code developed by Prof. Anthony Lewis

- [CLASS](https://github.com/lesgourg/class_public) is a Boltzmann code developed by Prof. Julien Lesgourgues and Dr. Thomas Tram

- [Polychord](https://github.com/PolyChord/PolyChordLite) is a sampler code developed by Dr. Will Handley, Prof. Lasenby, and Prof. M. Hobson

By no means, we want to discourage people from cloning code from their original repositories. We've included these codes as compressed [xz file format](https://tukaani.org/xz/format.html) in our repository for convenience in the initial development. The work of those authors is extraordinary, and they must be properly cited. Once all Cocoa's submodules go public, we intend to replace these files by clones to their original repositories. The same applies to all the data products and auxiliary packages we use.

## Git LFS
  [Git LFS](https://git-lfs.github.com) requires Git v1.8+. We don't provide any script/file to aid the installation of the required Git package.

### Total bandwidth for Git LFS files

**We have a monthly quota of only 150 GB in bandwidth for [Git LFS](https://git-lfs.github.com) files, and therefore we ask users to use good judgment in the number of times you clone repositories with large files**. **Users can use `git clean` and `git reset` commands to reset the repository to its original state in case some catastrophic edit is made**.

## Compiling Boltzmann, CosmoLike and Likelihood codes

To avoid excessive compilation times during development, users can use specialized scripts that compile only the specific modules:

    source ./installation_scripts/setup_class

    source ./installation_scripts/setup_camb

    source ./installation_scripts/setup_planck

    source ./installation_scripts/setup_polychord


## Running Jupyter Notebooks inside the [Whovian-Cosmo](https://hub.docker.com/r/vivianmiranda/whovian-cosmo) container

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
