# Overview: The [Cobaya](https://github.com/CobayaSampler)-[CosmoLike](https://github.com/CosmoLike) Joint Architecture

Cocoa allows users to run [CosmoLike](https://github.com/CosmoLike) routines that analyze data primarily from the [Dark Energy Survey](https://www.darkenergysurvey.org) (a.k.a DES), and simulate future multi-probe analyses, e.g. for Rubin Observatory's Legacy Survey of Space and Time or the Roman Space Telescope, inside the [Cobaya](https://github.com/CobayaSampler) framework. This readme file presents basic and advanced instructions for installing all Cocoa components, including the [Planck likelihood](https://wiki.cosmos.esa.int/planck-legacy-archive/index.php/Main_Page).

(**Proper Credits**): The following is not an exhaustive list of the codes we use

- [Cobaya](https://github.com/CobayaSampler) is a framework developed by Dr. Jesus Torrado and Prof. Anthony Lewis

- [Cosmolike](https://github.com/CosmoLike) is a framework developed by Prof. Elisabeth Krause and Prof. Tim Eifler

- [CAMB](https://github.com/cmbant/CAMB) is a Boltzmann code developed by Prof. Anthony Lewis

- [CLASS](https://github.com/lesgourg/class_public) is a Boltzmann code developed by Prof. Julien Lesgourgues and Dr. Thomas Tram

- [Polychord](https://github.com/PolyChord/PolyChordLite) is a sampler code developed by Dr. Will Handley, Prof. Lasenby, and Prof. M. Hobson

By no means, we want to discourage people from cloning code from their original repositories. We've included these codes as compressed [xz file format](https://tukaani.org/xz/format.html) in our repository for convenience in the initial development. The work of those authors is extraordinary, and they must be properly cited. Once all Cocoa's submodules go public, we intend to replace these files by clones to their original repositories. The same applies to all the data products and auxiliary packages we use.

# Installation of cocoa required packages

### Via Conda (best for Linux)

A simple way to install most prerequisites is via Conda environment. Assuming that the user had previously installed Minicoda/Anaconda, the first step is to type the following commands to create the cocoa Conda environment.

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
 
Then, every time the user wants to work on Cocoa, they must first activate the Cocoa Conda environment via

    $ conda activate cocoa

### Via Docker (best for MacOS/Windows)

Docker installation will allow users to run Cocoa inside an instantiation of the [Whovian-Cosmo](https://hub.docker.com/r/vivianmiranda/whovian-cosmo) docker image. Installation of the [docker engine](https://docs.docker.com/engine/) on local PCs is a reasonably straightforward process, but it does require `sudo` privileges (see Docker's [official documentation](https://docs.docker.com/engine/install/) for OS-specific instructions).

  On macOS, type:

    $ docker run -it -p 8080:8888 -v $(pwd):/home/whovian/host/ -v ~/.ssh:/home/whovian/.ssh:ro vivianmiranda/whovian-cosmo:version-1.0.2

For linux users, use the following command instead:

    $ docker run -it -p 8080:8888 --user $(id -u):$(id -g) -v $(pwd):/home/whovian/host/ -v ~/.ssh:/home/whovian/.ssh:ro vivianmiranda/whovian-cosmo:version-1.0.2

The bind of the `~/.ssh` folder in the commands above allows users to commit and clone from inside the container. Both commands need to be invoked on the parent folder so the docker container can access to the listed folders on the host OS. When running the command `docker run` on a specific container for the first time, the docker engine will automatically download the corresponding [Docker Hub](https://hub.docker.com/) image.
This step may take some time, as the [Whovian-Cosmo](https://hub.docker.com/r/vivianmiranda/whovian-cosmo) image has approximately 700 Megabytes.

  The last step is to access the folder `/home/whovian/host/` where the host files have been mounted:

    $ cd /home/whovian/host/

and proceed to the section [Cloning the Repository](https://github.com/CosmoLike/cocoa#cloning-the-repository). 

(**Warning**) HPC systems don't allow users to run docker containers using the standard [docker engine](https://docs.docker.com/engine/) for [security reasons](https://www.reddit.com/r/docker/comments/7y2yp2/why_is_singularity_used_as_opposed_to_docker_in/?utm_source=share&utm_medium=web2x&context=3). There is, however, an alternative engine called [Singularity](https://sylabs.io/guides/3.6/user-guide/index.html) that can run docker images in compliance with HPC security requirements. The [Singularity](https://sylabs.io/guides/3.6/user-guide/index.html) engine installation requires administrative privileges. However, many HPC systems have already adopted it.

To run docker images with Singularity, go to the folder you want to store the image and type:

    $ singularity build whovian-cosmo docker://vivianmiranda/whovian-cosmo

This command will download the [Whovian-Cosmo](https://hub.docker.com/r/vivianmiranda/whovian-cosmo) image and convert it to an image format that can be accessed via Singularity (Singularity takes a few minutes to do such conversion). To run the container interactively, type:

    $ singularity shell --no-home --bind /path/to/cocoa:/home/whovian/host --bind ~/.ssh:/home/whovian/.ssh:ro whovian-cosmo

after requesting interactive nodes successfully (never run jobs on login nodes!). The bind of the `~/.ssh` folder in the command above allows users to commit and clone from inside the container. The last step is to access the folder `/home/whovian/host/` where the host files have been mounted

    $ cd /home/whovian/host/

(**Warning**) There isn't permanent storage outside `/home/whovian/host/`. Be aware of this fact to not lose any work

### Via homebrew (need tweaks - not advisable)

We assume the user adopts [Homebrew](https://brew.sh) as the OS package manager and [BASH](https://www.howtogeek.com/444596/how-to-change-the-default-shell-to-bash-in-macos-catalina/) as the default shell. The installation procedures for MacOS are always working in progress; therefore, this section only offers a general guideline. Given how frequently [Homebrew](https://brew.sh) updates its packages, there is no guarantee that our recommendations will work without additional tweaks.  Below is a list of Homebrew and pip commands that install most dependencies

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

(**Warning**) It is also necessary to update `$PATH` so the OS can detect the [Homebrew](https://brew.sh) python installation. Adding to `~/.bash_profile` code similar to the lines below should be enough (don't forget to replace `XXX` by the user's home folder):

    python3=/Users/XXX/Library/Python/3.7/bin
    export PATH=$python3:$PATH
    
### Via Cocoa's internal cache (cumbersome - not advisable)

Whenever Conda or Docker installation procedures are unavailable, the user can still perform a local semi-autonomous installation on Linux based on scripts we provide. We provide a local copy of most required packages on Cocoa's cache folder [cocoa_installation_libraries](https://github.com/CosmoLike/cocoa/tree/main/cocoa_installation_libraries), as there are HPC machines where compute nodes that compile code don't have internet access, NASA Pleiades being one example. Therefore, we only assume the pre-installation of the following packages to perform the local setup:

   - [Bash](https://www.amazon.com/dp/B0043GXMSY/ref=cm_sw_em_r_mt_dp_x3UoFbDXSXRBT);
   - [Git](https://git-scm.com) v1.8+;
   - [Git LFS](https://git-lfs.github.com);
   - [gcc](https://gcc.gnu.org) v10.*;
   - [gfortran](https://gcc.gnu.org) v10.*;
   - [g++](https://gcc.gnu.org) v10.*;
   - [Python](https://www.python.org) v3.7.*;
   - [PIP package manager](https://pip.pypa.io/en/stable/installing/)
   - [Python Virtual Environment](https://www.geeksforgeeks.org/python-virtual-environment/)

To perform the local semi-autonomous installation, users should just follow the procedures on section [Installation of cocoa base code](https://github.com/CosmoLike/cocoa#installation-of-cocoa-base-code), adding, however, extra configurations on [set_installation_options](https://github.com/CosmoLike/cocoa/blob/master/Cocoa/set_installation_options) that are explained below.

The local installation via cocoa's internal cache is selected when the environment key `MANUAL_INSTALLATION` is set as shown below

    [Extracted from set_installation_options script] 
    
    #  ---------------------------------------------------------------------------
    # HOW COCOA BE INSTALLED? -------------------------------

    #export DOCKER_INSTALLATION=1
    #export MINICONDA_INSTALLATION=1
    #export MACOS_HOMEBREW_INSTALLATION=1
    export MANUAL_INSTALLATION=1
    
Users also need to set the following self-explanatory keys on [set_installation_options](https://github.com/CosmoLike/cocoa/blob/master/Cocoa/set_installation_options):
 
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
   
(**Warning**) Our scripts never install packages on `$HOME/.local` or other locations in the user's `$PATH` and `$LD_LIBRARY_PATH`. Doing so could impose incompatibilities between Cobaya and different projects. The decision to opt for global upgrades on packages required by multiple projects is the user's sole responsibility. All requirements for Cocoa are installed at

    Cocoa/.local/bin
    Cocoa/.local/include
    Cocoa/.local/lib
    Cocoa/.local/share

# Installation of cocoa base code

Type:

    $ git lfs clone git@github.com:CosmoLike/cocoa.git

(**Warning**) We have a monthly quota of only 150 GB in bandwidth for [Git LFS](https://git-lfs.github.com) files, and therefore we ask users to use good judgment in the number of times they clone repositories with large files (each clone will download around 5GB from Git LFS).

Cocoa chooses the preferred method of installation/compilation via special environment keys located on [set_installation_options](https://github.com/CosmoLike/cocoa/blob/master/Cocoa/set_installation_options) script (Docker is the default option). 
Choose accordingly to the used installation method of cocoa required packages  

    [Extracted from set_installation_options script]
    #  ---------------------------------------------------------------------------
    # HOW COCOA BE INSTALLED? -------------------------------

    #export DOCKER_INSTALLATION=1
    #export MINICONDA_INSTALLATION=1
    #export MACOS_HOMEBREW_INSTALLATION=1
    #export MANUAL_INSTALLATION=1
    
After the setup of the proper environment key, the installation/compilation of cocoa base code can be performed with the following commands:

    $ source setup_cocoa_installation_packages
    
    $ source compile_external_modules
    
    $ source start_cocoa

`setup_cocoa_installation_packages` this will decompress data files (from Planck for example) and install required packages that may have been left out in the Conda/Docker/homebrew procedure. File decompression should take a few minutes, while package installation time may range from a few seconds (conda/docker) to hours (manual installation).  

### Running Examples

After that, the Cobaya Framework should be ready, and the user can test a few examples:

    mpirun -n 1 cobaya-run ./projects/example/EXAMPLE_EVALUATE[1-4].yaml -f
    mpirun -n 4 cobaya-run ./projects/example/EXAMPLE_MCMC1.yaml -f
    
These examples will evaluate various likelihoods at specific cosmologies. The `-f` ensures that the same YAML file can be run multiple times, overwriting output files from previous evaluations that are located at `./chains`.

Once the work on the Cocoa environment is done, type:

    $ source stop_cocoa

PS: Users are encouraged to read appendix [Further Information about the configuration files](https://github.com/CosmoLike/cocoa/blob/master/README.md#further-information-about-the-configuration-files) to learn more about the configuration files.

# Download/compiling/running specific cocoa projects

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

### Compiling Boltzmann, CosmoLike and Likelihood codes separatelly

To avoid excessive compilation times during development, users can use specialized scripts that compile only the specific modules:

    source ./installation_scripts/setup_class

    source ./installation_scripts/setup_camb

    source ./installation_scripts/setup_planck

    source ./installation_scripts/setup_polychord

### Running Jupyter Notebooks inside the [Whovian-Cosmo](https://hub.docker.com/r/vivianmiranda/whovian-cosmo) container

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
   
### Summary Information about the configuration files

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
