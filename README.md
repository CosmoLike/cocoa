# Table of contents
1. [Overview of the Cobaya-CosmoLike Joint Architecture (Cocoa)](#overview)
2. [Installation of core packages via Conda](#required_packages_conda)
3. [Installation and Compilation of base code (Cobaya and external modules)](#cobaya_base_code)
4. [Running Examples (not involving Cosmolike)](#cobaya_base_code_examples)
5. [Running Cosmolike projects](#running_cosmolike_projects)
6. [Creating Cosmolike projects (external readme)](Cocoa/projects/)
7. [Appendix](#appendix)
    1. [Credits](#appendix_proper_credits)
    2. [Additional Installation Notes For Experts and Developers](#additional_notes)
    3. [FAQ: What if installation or compilation goes wrong?](#running_wrong)
    4. [FAQ: How to compile the Boltzmann, CosmoLike, and Likelihood codes separately](#appendix_compile_separately)
    5. [FAQ: How to run cocoa on a laptop? The docker image named *whovian-cocoa*](#appendix_jupyter_whovian)
    6. [FAQ: How to use an available Anaconda module on HPC?](#overview_anaconda)
    7. [FAQ: What if there is no Conda? Miniconda installation](#overview_miniconda)
    8. [FAQ: How to the Slow/Fast decomposition on MCMC chains with Cosmolike? Manual Blocking](#manual_blocking_cosmolike)
    9. [FAQ: How to switch Cocoa's adopted CAMB/CLASS/Polychord? (external readme)](Cocoa/external_modules/code)
    10. [FAQ: How to download modern CMB data? (external readme)](Cocoa/external_modules/data)
    11. [FAQ: How to set the environment for projects involving Machine Learning emulators?](#ml_emulators)
    12. [FAQ: How can users improve their Bash/C/C++ knowledge to develop Cocoa/Cosmolike?](#lectnotes)
    13. [Warning about Weak Lensing YAML files in Cobaya](#appendix_example_runs)
    14. [FAQ: How to install Cocoa without conda](#required_packages_cache)
    15. [FAQ: How to push changes to the Cocoa main branch?](#push_main)

## Overview of the [Cobaya](https://github.com/CobayaSampler)-[CosmoLike](https://github.com/CosmoLike) Joint Architecture (Cocoa) <a name="overview"></a>

Cocoa allows users to run [CosmoLike](https://github.com/CosmoLike) routines inside the [Cobaya](https://github.com/CobayaSampler) framework. [CosmoLike](https://github.com/CosmoLike) can analyze data primarily from the [Dark Energy Survey](https://www.darkenergysurvey.org) and simulate future multi-probe analyses for LSST and Roman Space Telescope. 

Besides integrating [Cobaya](https://github.com/CobayaSampler) and [CosmoLike](https://github.com/CosmoLike), Cocoa introduces shell scripts that allow users to containerize [Cobaya](https://github.com/CobayaSampler), the Boltzmann codes and multiple likelihoods. The container structure ensures that users can adopt consistent versions for the Fortran/C/C++ compilers and libraries across multiple machines; this greatly simplifies debugging. 

Our scripts never install packages, including Python modules, on `$HOME/.local` as that would make them global to the user. This behavior enables users to work on multiple instances of Cocoa simultaneously, similar to what was possible with [CosmoMC](https://github.com/cmbant/CosmoMC). 

This readme file presents basic and advanced instructions for installing all [Cobaya](https://github.com/CobayaSampler) and [CosmoLike](https://github.com/CosmoLike) components.

## Installation of core packages via Conda <a name="required_packages_conda"></a>

**Step :one:**: Download the file `cocoapy39.yml` yml file, create the cocoa environment, activate it, and create symbolic links that will give better names for the GNU compiler installed by Conda.

    conda env create --name cocoa --file=cocoapy39.yml
    conda activate cocoa
    ln -s "${CONDA_PREFIX}"/bin/x86_64-conda_cos6-linux-gnu-gcc "${CONDA_PREFIX}"/bin/gcc
    ln -s "${CONDA_PREFIX}"/bin/x86_64-conda_cos6-linux-gnu-g++ "${CONDA_PREFIX}"/bin/g++
    ln -s "${CONDA_PREFIX}"/bin/x86_64-conda_cos6-linux-gnu-gfortran "${CONDA_PREFIX}"/bin/gfortran
    ln -s "${CONDA_PREFIX}"/bin/x86_64-conda-linux-gnu-gcc-ar "${CONDA_PREFIX}"/bin/gcc-ar
    ln -s "${CONDA_PREFIX}"/bin/x86_64-conda-linux-gnu-gcc-ranlib "${CONDA_PREFIX}"/bin/gcc-ranlib

:interrobang: What if the user wants to install the Cocoa environment on a supercomputer? Many HPC environments offer Anaconda as an external module. If this is the case, check the Appendix [FAQ: How to use an available Anaconda module on HPC?](#overview_anaconda).

:interrobang: What if the user does not have conda installed? If the user is not working on an HPC environment that offers Anaconda, check the Appendix [FAQ: What if there is no Conda? Miniconda installation](#overview_miniconda).

**Step :two:**: Install `git-lfs` when loading the Conda cocoa environment for the first time.

    git-lfs install

## Installation and Compilation of base code (Cobaya and external modules) <a name="cobaya_base_code"></a>

**Step :one:**: We assume you are still in the Conda cocoa environment from the previous `conda activate cocoa` command. Now, clone the repository and go to the `cocoa` main folder,

    "${CONDA_PREFIX}"/bin/git clone --depth 1 https://github.com/CosmoLike/cocoa.git --branch v4.0-beta1 cocoa
    cd ./cocoa/Cocoa

*This will download the latest release, not the latest commit (more stable!)*.

:interrobang: If the user is a developer, then type the following instead *(at your own risk!)*

    "${CONDA_PREFIX}"/bin/git clone git@github.com:CosmoLike/cocoa.git cocoa
    cd ./cocoa/Cocoa

**Step :two:**: Run the script `setup_cocoa.sh` via
        
    source setup_cocoa.sh

The script `setup_cocoa.sh` decompresses the data files and installs a few necessary packages that have not been installed via conda.

**Step :three:**: Run the script `compile_cocoa.sh` by typing 

    source compile_cocoa.sh
    
This command compiles CAMB and Class Boltzmann codes, Planck likelihood, and Polychord sampler. 

:interrobang: Cocoa ignores a few external modules (code and likelihoods) by default, but users may find them helpful. In this case, check the available options on the script `set_installation_options.sh` and restart steps 2-3. 

## Running Examples (not involving Cosmolike)  <a name="cobaya_base_code_examples"></a>

We assume that you are still in the Conda cocoa environment from the previous `conda activate cocoa` command and that you are in the cocoa main folder `cocoa/Cocoa`, 

 **Step :one:**: Activate the private Python environment by sourcing the script `start_cocoa.sh`

    source start_cocoa.sh

Users will see a terminal like this: `$(cocoa)(.local)`. *This is a feature, not a bug*! 

 **Step :two:**: Select the number of OpenMP cores. Below, we set it to 4 for concreteness. 
    
    export OMP_PROC_BIND=close; export OMP_NUM_THREADS=4

The OpenMP affinity flag `OMP_PROC_BIND` is set to close place cores as closely as possible, which is important as current architectures severely limit communication bandwidth between distant cores.

 **Step :three:**: Run the `cobaya-run` command on the first example of the YAML files we provide.

One model evaluation:

    mpirun -n 1 --oversubscribe --mca btl vader,tcp,self --bind-to core:overload-allowed --rank-by core --map-by numa:pe=${OMP_NUM_THREADS} cobaya-run  ./projects/example/EXAMPLE_EVALUATE1.yaml -f
        
MCMC:

    mpirun -n 4 --oversubscribe --mca btl vader,tcp,self --bind-to core:overload-allowed --rank-by core --map-by numa:pe=${OMP_NUM_THREADS} cobaya-run ./projects/example/EXAMPLE_MCMC1.yaml -f

Once the work is done, clean your environment via :

    source stop_cocoa.sh
    
and

    conda deactivate cocoa

## Running Cosmolike projects <a name="running_cosmolike_projects"></a> 

The *projects* folder was designed to include Cosmolike projects. We assume that you are still in the Conda cocoa environment from the previous `conda activate cocoa` command and that you are in the cocoa main folder `cocoa/Cocoa`, 

**Step :one:**: Go to the project folder (`./cocoa/Cocoa/projects`) and clone a Cosmolike project with the fictitious name `XXX`:
    
    cd ./cocoa/Cocoa/projects
    "${CONDA_PREFIX}/bin/git" clone git@github.com:CosmoLike/cocoa_XXX.git XXX

By convention, the Cosmolike Organization hosts a Cobaya-Cosmolike project named XXX at `CosmoLike/cocoa_XXX`. When cloning the repository, the `cocoa_` prefix must be dropped. In the above `git clone` command, we eliminate the `cocoa_` prefix with the second `XXX` argument.

Example of cosmolike projects: [lsst_y1](https://github.com/CosmoLike/cocoa_lsst_y1).
 
**Step :two:**: Go back to the Cocoa main folder and activate the private Python environment
    
    cd ../
    source start_cocoa.sh
 
:warning: :warning: The `start_cocoa.sh` script must be run after cloning the project repository. 

Users will see a terminal like this: `$(cocoa)(.local)`. *This is a feature, not a bug*!

**Step :three:**: Compile the project, as shown below
 
    source "${ROOTDIR:?}"/projects/XXX/scripts/compile_XXX
  
**Step :four:**: Select the number of OpenMP cores and run a template YAML file
    
    export OMP_PROC_BIND=close; export OMP_NUM_THREADS=4
    mpirun -n 1 --oversubscribe --mca btl vader,tcp,self --bind-to core --rank-by core --map-by numa:pe=${OMP_NUM_THREADS} cobaya-run ./projects/XXX/EXAMPLE_EVALUATE1.yaml -f

## Appendix <a name="appendix"></a>

### Credits <a name="appendix_proper_credits"></a>

The following is not an exhaustive list of the codes we use/download/adopt

- [Cobaya](https://github.com/CobayaSampler) is a framework developed by Dr. Jesus Torrado and Prof. Anthony Lewis

- [Cosmolike](https://github.com/CosmoLike) is a framework developed by Prof. Elisabeth Krause and Prof. Tim Eifler

- [CAMB](https://github.com/cmbant/CAMB) is a Boltzmann code developed by Prof. Anthony Lewis

- [CLASS](https://github.com/lesgourg/class_public) is a Boltzmann code developed by Prof. Julien Lesgourgues and Dr. Thomas Tram

- [Polychord](https://github.com/PolyChord/PolyChordLite) is a sampler code developed by Dr. Will Handley, Prof. Lasenby, and Prof. M. Hobson

- [CLIK](https://github.com/benabed/clik) is the likelihood code used to analyze Planck and SPT data, maintained by Prof. Karim Benabed

- [SPT](https://github.com/SouthPoleTelescope/spt3g_y1_dist) is the official likelihood of the South Pole Telescope 3G Year 1

- [MFLike](https://github.com/simonsobs/LAT_MFLike) is the official likelihood of the Simons Observatory

- [ACTLensing](https://github.com/ACTCollaboration/act_dr6_lenslike) is the official lensing likelihood of the ACT collaboration developed by Prof. Mathew Madhavacheril

- [HiLLiPoP CMB likelihood](https://github.com/planck-npipe/hillipop.git) is a multifrequency CMB likelihood for Planck data.

- [Lollipop CMB likelihood](https://github.com/planck-npipe/lollipop.git) is a Planck low-l polarization likelihood.
  
Following best practices, Cocoa scripts download most external modules from their original repositories, including Cobaya, CAMB, Class, Polychord, ACT-DR6, HiLLiPoP, and Lollipop. We do not want to discourage people from cloning code from their original repositories. Our repository has included a few likelihoods as compressed [xz file format](https://tukaani.org/xz/format.html). The work of those authors is extraordinary, and users **must cite them** appropriately.

### Additional Installation Notes for experts and developers <a name="additional_notes"></a>

:books::books: *Installation of Cocoa's core packages via Conda* :books::books:
 
- For those working on projects that utilize machine-learning-based emulators, the Appendix [Setting-up conda environment for Machine Learning emulators](#ml_emulators) provides additional commands for installing the necessary packages.

- We provide a docker image named *whovian-cocoa* that facilitates cocoa installation on Windows and MacOS. For further instructions, refer to the Appendix [FAQ: How do you run cocoa on your laptop? The docker container is named *whovian-cocoa*](#appendix_jupyter_whovian).

We assume here that the user has previously installed either [Minicoda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://www.anaconda.com/products/individual) so that conda environments can be created. If this is not the case, refer to the Appendix [What should you do if you do not have Miniconda installed? Installation Guide](#overview_miniconda) for further instructions.

The conda installation method should be chosen in the overwhelming majority of cases. In the rare instances in which the user cannot work with Conda, refer to the Appendix [Installation of Cocoa's core packages without Conda](#required_packages_cache), as it contains instructions for a much slower (and prone to errors) but conda-independent installation method.

:books::books: *Installation of Cobaya base code* :books::books:

- If the user wants to compile only a subset of these packages, refer to the appendix [Compiling Boltzmann, CosmoLike, and Likelihood codes separately](#appendix_compile_separately).
        
- Cocoa developers should drop the shallow clone option `--depth 1`; they should also authenticate to GitHub via SSH keys and use the command instead.

      "${CONDA_PREFIX}"/bin/git clone git@github.com:CosmoLike/cocoa.git cocoa
  
- Our scripts never install packages on `$HOME/.local` as that would make them global to the user. All requirements for Cocoa are installed at

      Cocoa/.local/bin
      Cocoa/.local/include
      Cocoa/.local/lib
      Cocoa/.local/share

This behavior enables users to work on multiple instances of Cocoa simultaneously, similar to what was possible with [CosmoMC](https://github.com/cmbant/CosmoMC).

:books::books: *Running Cobaya Examples* :books::books:

- We offer the flag `COCOA_RUN_EVALUATE` as an alias (syntax-sugar) for `mpirun -n 1 --oversubscribe --mca btl vader,tcp,self --bind-to core:overload-allowed --rank-by core --map-by numa:pe=4 cobaya-run`.

- We offer the flag `COCOA_RUN_MCMC` as an alias (syntax-sugar) for `mpirun -n 4 --oversubscribe --mca btl vader,tcp,self --bind-to core:overload-allowed --rank-by core --map-by numa:pe=4 cobaya-run`. 

- Additional explanations about our `mpirun` flags: Why the `--mca btl vader,tcp,self` flag? Conda-forge developers don't [compile OpenMPI with Infiniband compatibility](https://github.com/conda-forge/openmpi-feedstock/issues/38).

- Additional explanations about our `mpirun` flags: Why the `--bind-to core:overload-allowed --map-by numa:pe=${OMP_NUM_THREADS}` flag? This flag enables efficient hybrid MPI + OpenMP runs on NUMA architecture.

- Additional explanations about the functioning `start_cocoa.sh`/`stop_cocoa.sh` scripts, starting from why we created two separate shell environments, `(cocoa)` and `(.local)`? Users should be able to manipulate multiple Cocoa instances seamlessly, which is particularly useful when running chains in one instance while experimenting with code development in another. Consistency of the environment across all Cocoa instances is crucial, and the `start_cocoa.sh`/`stop_cocoa.sh` scripts handle the loading and unloading of environmental path variables for each Cocoa.

### :interrobang: FAQ: What if installation or compilation goes wrong? <a name="running_wrong"></a>

- The script *set_installation_options script* contains a few additional flags that may be useful. Some of these flags are shown below:

      [Adapted from Cocoa/set_installation_options.sh shell script] 
      
      # ------------------------------------------------------------------------------
      # VERBOSE AS DEBUG TOOL --------------------------------------------------------
      # ------------------------------------------------------------------------------
      #export COCOA_OUTPUT_VERBOSE=1

      # ------------------------------------------------------------------------------
      # If set, COSMOLIKE will compile with DEBUG flags ------------------------------
      # ------------------------------------------------------------------------------
      #export COSMOLIKE_DEBUG_MODE=1
  
      # ------------------------------------------------------------------------------
      # The flags below allow users to skip downloading specific datasets ------------
      # ------------------------------------------------------------------------------
      # export IGNORE_BICEP_CMB_DATA=1
      # export IGNORE_HOLICOW_STRONG_LENSING_DATA=1
      # export IGNORE_SN_DATA=1
      # export IGNORE_BAO_DATA=1
      # export IGNORE_SPT_CMB_DATA=1
      export IGNORE_SIMONS_OBSERVATORY_CMB_DATA=1
      # export IGNORE_PLANCK_CMB_DATA=1
      export IGNORE_CAMSPEC_CMB_DATA=1
      export IGNORE_LIPOP_CMB_DATA=1

      (...)

      # ------------------------------------------------------------------------------
      # The keys below control which packages will be installed and compiled 
      # ------------------------------------------------------------------------------
      #export IGNORE_CAMB_CODE=1
      #export IGNORE_CLASS_CODE=1
      #export IGNORE_COSMOLIKE_CODE=1
      #export IGNORE_POLYCHORD_SAMPLER_CODE=1
      #export IGNORE_PLANCK_LIKELIHOOD_CODE=1
      #export IGNORE_ACTDR4_CODE=1
      #export IGNORE_ACTDR6_CODE=1
      #export IGNORE_CPP_CUBA_INSTALLATION=1
      #export IGNORE_VELOCILEPTORS_CODE=1
      #export IGNORE_SIMONS_OBSERVATORY_LIKELIHOOD_CODE=1
      #export IGNORE_CAMSPEC_LIKELIHOOD_CODE=1
      #export IGNORE_LIPOP_LIKELIHOOD_CODE=1

      (...)
 
Steps to debug Cocoa

- The first step is to define the `COCOA_OUTPUT_VERBOSE` and `COSMOLIKE_DEBUG_MODE` flags to obtain a more detailed output. To accomplish that, we advise users to uncomment the lines below that are part of the `set_installation_options.sh` script and then restart the cocoa private environment by running `source stop_cocoa.sh; source start_cocoa.sh`

      [Adapted from Cocoa/set_installation_options.sh shell script] 

      # ------------------------------------------------------------------------------
      # VERBOSE AS DEBUG TOOL --------------------------------------------------------
      # ------------------------------------------------------------------------------
      export COCOA_OUTPUT_VERBOSE=1

      # ------------------------------------------------------------------------------
      # If set, COSMOLIKE will compile with DEBUG flags ------------------------------
      # ------------------------------------------------------------------------------
      export COSMOLIKE_DEBUG_MODE=1

      (....)

- The second step consists of rerunning the particular script that failed with the verbose output set. The scripts `setup_cocoa.sh` and `compile_cocoa.sh` run many shell scripts. Users may find it advantageous to run only the routine that failed. For further information on how to do that, see the appendix [FAQ: How to compile the Boltzmann, CosmoLike, and Likelihood codes separately](#appendix_compile_separately).

After fixing a particular issue, users should rerun the shell scripts `setup_cocoa.sh` and `compile_cocoa.sh` to ensure all packages are installed and compiled correctly.

### :interrobang: FAQ: How to compile the Boltzmann, CosmoLike, and Likelihood codes separately <a name="appendix_compile_separately"></a>

To avoid excessive compilation or download times during development, users can use specialized scripts located at `Cocoa/installation_scripts/` that compile only a specific module or download only a particular dataset. A few examples of these scripts are: 

     $(cocoa)(.local) cd "${ROOTDIR:?}"
     $(cocoa)(.local) source "${ROOTDIR:?}"/installation_scripts/compile_act_dr4.sh
     $(cocoa)(.local) source "${ROOTDIR:?}"/installation_scripts/compile_camb.sh
     $(cocoa)(.local) source "${ROOTDIR:?}"/installation_scripts/compile_class.sh
     $(cocoa)(.local) source "${ROOTDIR:?}"/installation_scripts/compile_planck.sh
     $(cocoa)(.local) source "${ROOTDIR:?}"/installation_scripts/compile_polychord.sh

Above and below, the `$(cocoa)(.local)` emphasizes they should run after activating the cocoa environments. The shell subroutines that download external modules from their original Git repositories are shown below.

     $(cocoa)(.local) cd "${ROOTDIR:?}"
     $(cocoa)(.local) source "${ROOTDIR:?}"/installation_scripts/setup_act_dr4.sh
     $(cocoa)(.local) source "${ROOTDIR:?}"/installation_scripts/setup_camb.sh
     $(cocoa)(.local) source "${ROOTDIR:?}"/installation_scripts/setup_class.sh
     $(cocoa)(.local) source "${ROOTDIR:?}"/installation_scripts/setup_polychord.sh
     $(cocoa)(.local) source "${ROOTDIR:?}"/installation_scripts/setup_cobaya.sh

To ensure these scripts can download and install these packages, users must be sure that the environment keys below are *NOT* set. These keys are shown on `set_installation_options.sh`. The command `unset -v` unset them. 
      
     unset -v IGNORE_CAMB_CODE
     unset -v IGNORE_CLASS_CODE
     unset -v IGNORE_POLYCHORD_SAMPLER_CODE
     unset -v IGNORE_PLANCK_LIKELIHOOD_CODE
     unset -v IGNORE_ACTDR4_CODE
     unset -v IGNORE_COBAYA_CODE

Below, we show the shell subroutines that download and unpack data from multiple experiments. 

     $(cocoa)(.local) cd "${ROOTDIR:?}"
     $(cocoa)(.local) source "${ROOTDIR:?}"/installation_scripts/unxv_act_dr6.sh
     $(cocoa)(.local) source "${ROOTDIR:?}"/installation_scripts/unxv_bao.sh
     $(cocoa)(.local) source "${ROOTDIR:?}"/installation_scripts/unxv_bicep.sh
     $(cocoa)(.local) source "${ROOTDIR:?}"/installation_scripts/unxv_camspec.sh
     $(cocoa)(.local) source "${ROOTDIR:?}"/installation_scripts/unxv_h0licow.sh
     $(cocoa)(.local) source "${ROOTDIR:?}"/installation_scripts/unxv_lipop.sh
     $(cocoa)(.local) source "${ROOTDIR:?}"/installation_scripts/unxv_planck2018_basic.sh
     $(cocoa)(.local) source "${ROOTDIR:?}"/installation_scripts/unxv_simons_observatory.sh
     $(cocoa)(.local) source "${ROOTDIR:?}"/installation_scripts/unxv_sn.sh
     $(cocoa)(.local) source "${ROOTDIR:?}"/installation_scripts/unxv_spt.sh

To ensure these scripts can download these datasets, users must be sure that the environment keys below are *NOT* set. These keys are shown on `set_installation_options.sh`. The command `unset -v` unset them. 

     unset -v IGNORE_ACTDR6_DATA
     unset -v IGNORE_BAO_DATA
     unset -v IGNORE_BICEP_CMB_DATA
     unset -v IGNORE_CAMSPEC_CMB_DATA
     unset -v IGNORE_HOLICOW_STRONG_LENSING_DATA
     unset -v IGNORE_LIPOP_CMB_DATA
     unset -v IGNORE_PLANCK_CMB_DATA
     unset -v IGNORE_SIMONS_OBSERVATORY_CMB_DATA
     unset -v IGNORE_SN_DATA
     unset -v IGNORE_SPT_CMB_DATA

### :interrobang: FAQ: How to run cocoa on a laptop? The docker image named *whovian-cocoa* <a name="appendix_jupyter_whovian"></a>

We provide the docker image [whovian-cocoa](https://hub.docker.com/r/vivianmiranda/whovian-cocoa) to facilitate the installation of Cocoa on Windows and MacOS. This appendix assumes the users already have the docker engine installed on their local PC. For instructions on installing the docker engine in specific operating systems, please refer to [Docker's official documentation](https://docs.docker.com/engine/install/). 

To download the docker image *whovian-cocoa*, name the associated container `cocoa2023`, and run the container for the first time, type:

      docker run --platform linux/amd64 --hostname cocoa --name cocoa2023 -it -p 8080:8888 -v $(pwd):/home/whovian/host/ -v ~/.ssh:/home/whovian/.ssh:ro vivianmiranda/whovian-cocoa

Following the command above, users should see the following text on the screen terminal:

<img width="872" alt="Screenshot 2023-12-19 at 1 26 50 PM" src="https://github.com/CosmoLike/cocoa/assets/3210728/eb1fe7ec-e463-48a6-90d2-2d84e5b61aa1">

The user needs to init Conda when running the container the first time, as shown below.

      conda init bash
      source ~/.bashrc

Now, proceed with the standard cocoa installation. 

Once installation is complete, the user must learn how to start, use, and exit the container. Below, we answer a few common questions about using/managing Docker containers.  

- :interrobang: FAQ: How do users restart the container when they exit?

    Assuming the user maintained the container name `cocoa2023` via the flag `--name cocoa2023` on the `docker run` command, type:
    
      docker start -ai cocoa2023

- :interrobang: FAQ: How do I run Jupyter Notebooks remotely when using Cocoa within the *whovian-cocoa* container?

    First, type the following command:

      jupyter notebook --no-browser --port=8080

    The terminal will show a message similar to the following template:

      [... NotebookApp] Writing notebook server cookie secret to /home/whovian/.local/share/jupyter/runtime/notebook_cookie_secret
      [... NotebookApp] WARNING: The notebook server is listening on all IP addresses and not using encryption. This is not recommended.
      [... NotebookApp] Serving notebooks from local directory: /home/whovian/host
      [... NotebookApp] Jupyter Notebook 6.1.1 is running at:
      [... NotebookApp] http://f0a13949f6b5:8888/?token=XXX
      [... NotebookApp] or http://127.0.0.1:8888/?token=XXX
      [... NotebookApp] Use Control-C to stop this server and shut down all kernels (twice to skip confirmation).

    If you are running the Docker container on your laptop, there is only one remaining step. The flag `-p 8080:8888` in the `docker run` command maps the container port `8888` to the host (your laptop) port `8080`. Therefore, open a browser and enter `http://localhost:8080/?token=XXX`, where `XXX` is the token displayed in the output line `[... NotebookApp] or http://127.0.0.1:8888/?token=XXX`, to access the Jupyter Notebook. 

    If you need to use a different port than `8080`, adjust the flag `-p 8080:8888` in the `docker run` command accordingly.

- :interrobang: **FAQ: How to manipulate files on the host computer from within the Docker container?**

    The flag `-v $(pwd):/home/whovian/host/` in the `docker run` command ensures that files on the host computer have been mounted to the directory `/home/whovian/host/`. Files within the folder where the Docker container was initialized are accessible in the `/home/Whovian/host/` folder. When the user accesses the container, they should see the host's directory where the Docker container was initiated after typing: 

      cd /home/whovian/host/
      ls 

    Users should work inside the `/home/whovian/host/` directory to avoid losing work in case the docker image needs to be deleted,

- :interrobang: **FAQ: How to run the docker container on a remote server?**

    Below, we assume the user runs the container in a server with the URL `your_sever.com`. We also presume the server can be accessed via SSH protocol. On your local PC/laptop, type:

      ssh your_username@your_sever.com -L 8080:localhost:8080

    This will bind the port `8080` on the server to the local. Then, go to a browser and type `http://localhost:8080/?token=XXX`, where `XXX` is the previously saved token displayed in the line `[... NotebookApp] or http://127.0.0.1:8888/?token=XXX`.  

### :interrobang: FAQ: How to use an available Anaconda module on HPC? <a name="overview_anaconda"></a>

Below, we list users' most common issues when installing Cocoa conda environments in a supercomputer environment using a globally defined Anaconda module. 

- :interrobang: **Conda command not found**.
  
Anaconda is not usually set by default on HPC environments but may be available as a module. For example, on the Midway HPC cluster, it can be loaded using the command below.

    module load Anaconda3/2022.10

If you are not sure about the name of the available Anaconda module, type the command.

    module avail An

to show all modules with names that start with `An`. The output should resemble the following.

<img width="700" alt="Anaconda" src="https://github.com/user-attachments/assets/09326f5f-49e0-45b5-a157-25fe2b09918e">

- :interrobang: **Installation seems to take forever**.

There are various reasons why installing the Cocoa conda environment may take a long time. Here is a checklist of good practices to overcome this problem.

:one: *Never install conda environments using the login node*. 

Instead, request an interactive job with a few cores. However, users must know that **some supercomputers do not provide internet access on computing nodes** (e.g., the Midway HPC at the University of Chicago). Ask the HPC staff for a **queue dedicated to installing and compiling code** in this case; they should exist in a well-designed HPC environment. For example, the `build partition` on the Midway supercomputer can be accessed with the command.

    sinteractive --nodes=1 --ntasks=1 --cpus-per-task=5 --time=3:00:00 --account=pi-XXX --partition=build

:two: *Set conda verbosity to DEBUG* in case the installation still takes too much time after step :one:.

    conda config --set verbosity 2

The DEBUG mode will ensure conda outputs many more intermediate installation steps, which helps the user check whether Conda is stuck on some intermediate step. Users can later reset the verbosity level with the command `conda config --set verbosity 0`.

- **Conda installation is interrupted due to quota limitations**.

Supercomputers usually enforce strict quota limits on home folders. These limits apply to the total file size and the number of files. By default, Anaconda modules install new environments at `$HOME/.conda/envs`. Anaconda also stores Gigabytes of downloaded packages in the `$HOME/.conda/pkgs` folder; `pkgs` is used by Anaconda as a package cache folder. Therefore, reasonable and widely applied quota limitations to the home folder significantly hinder the installation of new environments without the proposed changes below. 

:one: Create an Anaconda folder on a project folder outside `$HOME` with significantly more tolerant quota restrictions. For instance, on the Midway supercomputer I create the Anaconda folder at KICP projects folder with the command.

    mkdir /project2/kicp/XXX/anaconda/

:two: Set the pkgs package cache folder to `anaconda/pkgs`.

    conda config --add pkgs_dirs /project2/kicp/XXX/anaconda/pkgs

:tree: Set the env folder to `anaconda/envs/` 

    conda config --add envs_dirs /project2/kicp/XXX/anaconda/envs

:four: Set the following flags for a proper conda environment installation

    conda config --set auto_update_conda false
    conda config --set show_channel_urls true
    conda config --set auto_activate_base false
    conda config --prepend channels conda-forge
    conda config --set channel_priority strict
    conda init bash

:five: Source the `.bashrc` file so the new Anaconda settings are loaded. 

    source ~/.bashrc

:six: After completing steps :one:-:four:, check that the `$HOME/.condarc` file with the command

    more $HOME/.condarc

to make sure it resembles the one below.

    auto_update_conda: false
    show_channel_urls: true
    auto_activate_base: false
    channels:
      - conda-forge
      - defaults
    channel_priority: strict
    verbosity: 0
    pkgs_dirs:
      - /project2/kicp/XXX/anaconda/pkgs
    envs_dirs:
      - /project2/kicp/XXX/anaconda/envs/

### :interrobang: FAQ: What if there is no Miniconda? Miniconda installation <a name="overview_miniconda"></a>

Download and run the Miniconda installation script. 

    export CONDA_DIR="/gpfs/home/XXX/miniconda"
    
    mkdir "${CONDA_DIR:?}"
    
    wget https://repo.continuum.io/miniconda/Miniconda3-py38_23.9.0-0-Linux-x86_64.sh
    
    /bin/bash Miniconda3-py38_23.9.0-0-Linux-x86_64.sh -f -b -p "${CONDA_DIR:?}"

Please don't forget to adapt the path assigned to `CONDA_DIR` in the command above:

After installation, users must source the conda configuration file, as shown below:

    source $CONDA_DIR/etc/profile.d/conda.sh \
          && conda config --set auto_update_conda false \
          && conda config --set show_channel_urls true \
          && conda config --set auto_activate_base false \
          && conda config --prepend channels conda-forge \
          && conda config --set channel_priority strict \
          && conda init bash

### :interrobang: FAQ: How to set the Slow/Fast decomposition on MCMC chains with Cosmolike? Manual Blocking <a name="manual_blocking_cosmolike"></a>

The Cosmolike Weak Lensing pipelines contain parameters with different speed hierarchies. For example, Cosmolike execution time is reduced by approximately 50% when fixing the cosmological parameters. When varying only multiplicative shear calibration, Cosmolike execution time is reduced by two orders of magnitude. 

Cobaya cannot automatically handle parameters associated with the same likelihood that have different speed hierarchies. Luckily, we can manually impose the speed hierarchy in Cobaya using the `blocking:` option. The only drawback of this method is that parameters of all adopted likelihoods, not only the ones required by Cosmolike, must be manually specified.

In addition to that, Cosmolike can't cache the intermediate products of the last two evaluations, which is necessary to exploit optimizations associated with dragging (`drag: True`). However, Cosmolike caches the intermediate products of the previous evaluation, thereby enabling the user to take advantage of the slow/fast decomposition of parameters in Cobaya's main MCMC sampler. 

Below, we provide an example YAML configuration for an MCMC chain with DES 3x2pt likelihood.

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
                Rminus1_stop: 0.015
                # Gelman-Rubin R-1 on std deviations
                Rminus1_cl_stop: 0.17
                Rminus1_cl_level: 0.95
                # ---------------------------------------------------------------------
                # Exploiting Cosmolike speed hierarchy
                # ---------------------------------------------------------------------
                measure_speeds: False # We provide the approximate speeds in the blocking
                # drag = false. The drag sampler requires the intermediate products of the last
                # two evaluations to be cached. Cosmolike can only cache the last evaluation.
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
                max_tries: 100000
                burn_in: 0
                Rminus1_single_split: 4

### :interrobang: FAQ: How do users set the environment for projects involving Machine Learning emulators? <a name="ml_emulators"></a>

Commenting out the environmental flags shown below, located at *set_installation_options* script, will enable the installation of machine-learning-related libraries via pip.  

    [Adapted from Cocoa/set_installation_options.sh shell script] 
     
    # ------------------------------------------------------------------------------
    # If not set, Cocoa/installation_scripts/setup_pip_core_packages.sh will install several ML packages
    # ------------------------------------------------------------------------------
    #export IGNORE_EMULATOR_CPU_PIP_PACKAGES=1
    #export IGNORE_EMULATOR_GPU_PIP_PACKAGES=1
        
We recommend using the GPU versions to train the emulator while using the CPU versions to run the MCMCs. This can be achieved by creating two separate environments, i.e., two independent CoCoA copies. Alternatively, users can manually incorporate the pip ML installation commands located on the shell script `Cocoa/installation_scripts/setup_pip_core_packages.sh` in their conda environments (two conda environments, one for the GPU and one for the CPU runs).

    [Adapted from Cocoa/installation_scripts/setup_pip_core_packages.sh script - the command may not be entirely up to date] 
     
    #CPU VERSION
    env CXX="${CXX_COMPILER:?}" CC="${C_COMPILER:?}" ${PIP3:?} install \
        'tensorflow-cpu==2.12.0' \
        'tensorflow_probability-0.21.0' \
        'keras==2.12.0' \
        'keras-preprocessing==1.1.2' \
        'torch==1.13.1+cpu' \
        'torchvision==0.14.1+cpu' \
        'torchaudio==0.13.1' \
        'tensiometer==0.1.2' \
      --extra-index-url "https://download.pytorch.org/whl/cpu" \
      --prefix="${ROOTDIR:?}/.local" 

    (...)

    #GPU VERSION
    env CXX="${CXX_COMPILER:?}" CC="${C_COMPILER:?}" ${PIP3:?} install \
        'tensorflow==2.12.0' \
        'tensorflow_probability-0.21.0' \
        'keras==2.12.0' \
        'keras-preprocessing==1.1.2' \
        'torch==1.13.1+cu116' \
        'torchvision==0.14.1+cu116' \
        'torchaudio==0.13.1' \
        'tensiometer==0.1.2' \
      --extra-index-url "https://download.pytorch.org/whl/cu116" \
      --prefix="${ROOTDIR:?}/.local"
      
### :interrobang: FAQ: How can users improve their Bash/C/C++ knowledge to develop Cocoa/Cosmolike? :book::book: <a name="lectnotes"></a>

A working knowledge of Python is required to understand the Cobaya framework at the developer level. Users must also know the Bash language to understand Cocoa's scripts. Proficiency in C and C++ is also needed to manipulate Cosmolike and the C++ Cobaya-Cosmolike C++ interface. Finally, users need to understand the Fortran-2003 language to modify CAMB.

Learning all these languages can be overwhelming, so to enable new users to do research that demands modifications on the inner workings of these codes, we include [here](cocoa_installation_libraries/LectNotes.pdf) a link to approximately 600 slides that provide an overview of Bash (slides ~1-137), C (slides ~138-371), and C++ (slides ~372-599). In the future, we aim to add lectures about Python and Fortran. 

### :warning::warning: Warning about Weak Lensing YAML files in Cobaya <a name="appendix_example_runs"></a>

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

Adopting $\big(\Omega_m,\Omega_b\big)$ as main MCMC parameters can create a silent bug in Cobaya; *we are unsure if this problem persists in newer Cobaya versions, so this report should be a warning*. The problem occurs when the option `drop: true` is absent in $\big(\Omega_m,\Omega_b\big)$ parameters, and there are no expressions that define the derived $\big(\Omega_c h^2, \Omega_b h^2\big)$ quantities. *The bug is silent* because the MCMC runs without any warnings, but the CAMB Boltzmann code does not update the cosmological parameters at every MCMC iteration. As a result, the resulting posteriors are flawed, but they may seem reasonable to those unfamiliar with the issue. Please be aware of this bug to avoid any potential inaccuracies in the results. 

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
            
### 💀 ☠️ :stop_sign::thumbsdown: FAQ: How to install Cocoa without conda (not recommended) <a name="required_packages_cache"></a>

This method is slow and not advisable :stop_sign::thumbsdown:. When Conda is unavailable, the user can still perform a local semi-autonomous installation on Linux based on a few scripts we implemented. We require the pre-installation of the following packages:

   - [Bash](https://www.amazon.com/dp/B0043GXMSY/ref=cm_sw_em_r_mt_dp_x3UoFbDXSXRBT);
   - [Git](https://git-scm.com) v1.8+;
   - [Git LFS](https://git-lfs.github.com);
   - [gcc](https://gcc.gnu.org) v12.*+;
   - [gfortran](https://gcc.gnu.org) v12.*+;
   - [g++](https://gcc.gnu.org) v12.*+;
   - [Python](https://www.python.org) v3.8.*;
   - [wget](https://www.gnu.org/software/wget/) v1.16+;
   - [curl](https://curl.se)
   - [PIP package manager](https://pip.pypa.io/en/stable/installing/);
   - [Python Virtual Environment](https://www.geeksforgeeks.org/python-virtual-environment/);
    
To perform the local semi-autonomous installation, users must modify flags written on the shell scripts *set_installation_options.sh* and `installation_scripts/flags_manual_installation.sh` because the default behavior corresponds to an installation via Conda. First, select the environmental key `MANUAL_INSTALLATION` as shown below:

    [Adapted from Cocoa/set_installation_options.sh script] 

    # ------------------------------------------------------------------------------------
    # HOW COCOA SHOULD BE INSTALLED? -----------------------------------------------------
    # ------------------------------------------------------------------------------------
    #export MINICONDA_INSTALLATION=1
    export MANUAL_INSTALLATION=1
    
Finally, set the following environmental keys:
 
    [Adapted from Cocoa/installation_scripts/flags_manual_installation.sh shell script]
    # ------------------------------------------------------------------------------------
    # IF SET, COCOA DOES NOT USE SYSTEM PIP PACKAGES -------------------------------------
    # ------------------------------------------------------------------------------------
    export DONT_USE_SYSTEM_PIP_PACKAGES=1

    (...)

    # ------------------------------------------------------------------------------------
    # USER NEEDS TO SPECIFY THE FLAGS BELOW SO COCOA CAN FIND PYTHON / GCC ---------------
    # ------------------------------------------------------------------------------------
	
    export PYTHON_VERSION=XXX
    export GLOBALPYTHON3=XXX
    export GLOBAL_PACKAGES_LOCATION=XXX
    export GLOBALPIP3=XXX
    export GIT=XXX
    export WGET=XXX
    export CURL=XXX

    # ------------------------------------------------------------------------------------

    export PATH=XXX:$PATH
    export CFLAGS="${CFLAGS} -IXXX"
    export LDFLAGS="${LDFLAGS} -LXXX"
    export C_INCLUDE_PATH=XXX:$C_INCLUDE_PATH
    export CPLUS_INCLUDE_PATH=XXX:$CPLUS_INCLUDE_PATH
    export PYTHONPATH=XXX:$PYTHONPATH
    export LD_RUN_PATH=XXX:$LD_RUN_PATH
    export LIBRARY_PATH=XXX:$LIBRARY_PATH
    export CMAKE_INCLUDE_PATH=XXX:$CMAKE_INCLUDE_PATH
    export CMAKE_LIBRARY_PATH=XXX:$CMAKE_LIBRARY_PATH
    export INCLUDE_PATH=XXX:$INCLUDE_PATH
    export INCLUDE=XXX:$INCLUDE
    export CPATH=XXX:$CPATH
    export OBJC_INCLUDE_PATH=XXX:$OBJC_INCLUDE_PATH
    export OBJC_PATH=XXX:$OBJC_PATH
                
    # ------------------------------------------------------------------------------------
    # COMPILER ---------------------------------------------------------------------------
    # ------------------------------------------------------------------------------------
    export C_COMPILER=
    export CXX_COMPILER=
    export FORTRAN_COMPILER=
    export MPI_CC_COMPILER=
    export MPI_CXX_COMPILER=    
    export MPI_FORTRAN_COMPILER=

    # ------------------------------------------------------------------------------------
    # FINE-TUNNING OVER THE USE OF SYSTEM-WIDE PACKAGES ---------------------------------
    # ------------------------------------------------------------------------------------
    #export IGNORE_XZ_INSTALLATION=1
    #export IGNORE_DISTUTILS_INSTALLATION=1
    #export IGNORE_C_GSL_INSTALLATION=1
    #export IGNORE_C_CFITSIO_INSTALLATION=1
    #export IGNORE_C_FFTW_INSTALLATION=1
    #export IGNORE_CPP_BOOST_INSTALLATION=1
    #export IGNORE_CMAKE_INSTALLATION=1
    #export IGNORE_OPENBLAS_INSTALLATION=1
    #export IGNORE_FORTRAN_LAPACK_INSTALLATION=1
    #export IGNORE_CPP_ARMA_INSTALLATION=1
    #export IGNORE_HDF5_INSTALLATION=1
    #export IGNORE_EXPAT_CORE_PACKAGE=1
    #export IGNORE_PIP_CORE_PACKAGES=1

### :interrobang: FAQ: How to push changes to the Cocoa main branch? <a name="push_main"></a>

Until recently, Cocoa development was unstructured, and developers were allowed to push directly to the `main` branch. Small commits were also not discouraged. This loose development rules will soon change. In particular, we will soon protect the `main` branch by requiring every push to be reviewed by Cocoa's leading developers. We also want to reduce the number of commits from now on. 

We will implement the philosophy that **a commit in the main branch should contain an atomic change that takes code from a working state to another working state with improvements that are meaningful and well tested.** So, developers should propose changes to the `main` branch in larger chunks, as shown below. 

Important note: we **strongly** advise developers to use the up-to-date `git` given by Cocoa conda environment. 

**Step :one:**: create a development branch from the `main` branch. Do not call the development branch `dev`, it is reserved for work done by the primary Cocoa developers. Let's call this new branch `xyzdev` for concreteness (tip: the use of developers' initials helps make the branch name unique)

    # (main) = developers must input the command below on the main branch.
    (main) git checkout -b xyzdev

In the case where the `xyzdev` already exists, use `git switch` instead of `git checkout -b`. Developers also have freedom to push their develop branches to server via the command

    # (xyzdev) = developers must input the command below on the xyzdev branch.
    (xyzdev) git push -u origin xyzdev


**Step :two:**: develop the proposed changes. We advise developers to commit frequently. In your branch, a commit does not need to be atomic, changing the code from one working state to another well-tested meaningful working state. In your branch, you have absolute freedom.

**Step :three:**: Once there is an atomic, meaningful, and well-tested improvement to Cocoa, the developer needs first to merge any subsequent changes made in `main` while the developer has been working on the `xyzdev` branch.

    # (xyzdev) = developers must input the command below on the xyzdev branch.
    (xyzdev) git merge main

This step may create conflicts that must be addressed before step four. 

**Step :four:**: Once you have merge recent changes made on the `main` branch, we will push to the main branch the changes made on the `xyzdev` branch by first **squashing all your changes into a single commit** as shown below

    # (xyzdev) = developers must input the command below on the xyzdev branch.
    (xyzdev) git switch main
    
    # (main) = developers must input the command below on the main branch.
    (main) git merge --squash xyzdev

    (main) git commit -m "squash merge - xyzdev branch: added development on abc features"

    (main) git push origin main

Important note: **never** revert the branch ordering on squash merging by squashing the `main` changes to the `xyzdev` branch.
