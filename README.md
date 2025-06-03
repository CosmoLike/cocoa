
<img width="780" alt="Screenshot 2024-09-16 at 10 55 14 PM" src="https://github.com/user-attachments/assets/f327e3bf-22c7-46e4-9bb5-14c8e7afc4c1">

# Table of contents
1. [Overview of the Cobaya-CosmoLike Joint Architecture (Cocoa)](#overview) 
2. [Installation of core packages via Conda](#required_packages_conda)
3. [Installation and Compilation of external modules](#cobaya_base_code)
4. [Running Examples](#cobaya_base_code_examples)
5. [Running ML emulators](#cobaya_base_code_examples_emul)
6. [Creating Cosmolike projects (external readme)](Cocoa/projects/)
7. [Appendix](#appendix)
    1. [Credits](#appendix_proper_credits)
    2. [FAQ: What if installation or compilation goes wrong?](#running_wrong)
    3. [FAQ: How do we compile the Boltzmann, CosmoLike, and Likelihood codes separately](#appendix_compile_separately)
    4. [FAQ: How do we run Cocoa on a laptop? The docker image named *whovian-cocoa*](#appendix_jupyter_whovian)
    5. [FAQ: How do we use an available Anaconda module on HPC?](#overview_anaconda)
    6. [FAQ: What if there is no Conda? Miniconda installation](#overview_miniconda)
    7. [FAQ: How do we make the Slow/Fast decomposition on MCMC chains with Cosmolike? Manual Blocking](#manual_blocking_cosmolike)
    8. [FAQ: How do we switch Cocoa's adopted CAMB/CLASS/Polychord? (external readme)](Cocoa/external_modules/code)
    9. [FAQ: How do we download modern CMB data? (external readme)](Cocoa/external_modules/data)
    10. [FAQ: How do we set the environment for Machine Learning projects?](#ml_emulators)
    11. [FAQ: How can users improve their Bash/C/C++ knowledge to develop Cocoa/Cosmolike?](#lectnotes)
    12. [Warning about Weak Lensing YAML files in Cobaya](#appendix_example_runs)
    13. [FAQ: How do we push changes to the Cocoa main branch? A few git hacks](#push_main)
    14. [FAQ: How do we develop from a Git tag? A few more Git hacks](#dev_from_tag)
    15. [FAQ: How do we download and run Cosmolike projects?](#running_cosmolike_projects)
    16. [FAQ: How do we add a new package to Cocoa? The Dark Emulator Example](#add_package_v1)
    17. [FAQ: How do we add a new package to Cocoa? The MGCAMB Example](#add_package_v2)
    
## Overview of the [Cobaya](https://github.com/CobayaSampler)-[CosmoLike](https://github.com/CosmoLike) Joint Architecture (Cocoa) <a name="overview"></a>

Cocoa allows users to run [CosmoLike](https://github.com/CosmoLike) routines inside the [Cobaya](https://github.com/CobayaSampler) framework. [CosmoLike](https://github.com/CosmoLike) can analyze data primarily from the [Dark Energy Survey](https://www.darkenergysurvey.org) and simulate future multi-probe analyses for LSST and Roman Space Telescope. 

Besides integrating [Cobaya](https://github.com/CobayaSampler) and [CosmoLike](https://github.com/CosmoLike), Cocoa introduces shell scripts that allow users to containerize [Cobaya](https://github.com/CobayaSampler), the Boltzmann codes, and multiple likelihoods. The container structure of Cocoa ensures that users can adopt consistent versions for the Fortran/C/C++ compilers and libraries across multiple machines. Such a systematic approach greatly simplifies the debugging process. 

Our scripts never install packages or Python modules on a globlal folder such as `$HOME/.local`. Here, `$HOME` denotes a shell environment variable that points to the user's home folder. Doing so would force cocoa packages to be global to the user, possibly breaking environments. Our scripts enables users to work on multiple Cocoa instances simultaneously, similar to what was possible with [CosmoMC](https://github.com/cmbant/CosmoMC). 

This readme file presents basic and advanced instructions for installing all [Cobaya](https://github.com/CobayaSampler) and [CosmoLike](https://github.com/CosmoLike) components.

## Installation of core packages via Conda <a name="required_packages_conda"></a>

Core packages include compilers and numerical libraries that users never modify. We install most of these core packages via Conda, as shown below.

**Step :one:**: Download the file `cocoapy310.yml` yml file

    wget https://raw.githubusercontent.com/CosmoLike/cocoa/refs/heads/main/cocoapy310.yml

create the cocoa environment,

    conda env create --name cocoa --file=cocoapy310.yml

activate it

    conda activate cocoa

**Step :two:**: When and only when loading the conda cocoa environment for the first time, create symbolic links that will give better names for the GNU compilers

    ln -s "${CONDA_PREFIX}"/bin/x86_64-conda_cos6-linux-gnu-gcc "${CONDA_PREFIX}"/bin/gcc
    ln -s "${CONDA_PREFIX}"/bin/x86_64-conda_cos6-linux-gnu-g++ "${CONDA_PREFIX}"/bin/g++
    ln -s "${CONDA_PREFIX}"/bin/x86_64-conda_cos6-linux-gnu-gfortran "${CONDA_PREFIX}"/bin/gfortran
    ln -s "${CONDA_PREFIX}"/bin/x86_64-conda-linux-gnu-gcc-ar "${CONDA_PREFIX}"/bin/gcc-ar
    ln -s "${CONDA_PREFIX}"/bin/x86_64-conda-linux-gnu-gcc-ranlib "${CONDA_PREFIX}"/bin/gcc-ranlib
    ln -s "${CONDA_PREFIX}"/bin/x86_64-conda-linux-gnu-ld "${CONDA_PREFIX}"/bin/ld

and install `git-lfs`

    git-lfs install

Users can now proceed to the **next section**.

> [!TIP]
> To install the `cocoa` conda environment on a supercomputer, users may take advantage of the fact that several HPC environments provide the [Anaconda installer](https://www.anaconda.com) as an external module. If this applies to you, then check the appendix [FAQ: How do we use an available Anaconda module on HPC?](#overview_anaconda).

> [!TIP]
> Users working on an HPC environment that does not offer Anaconda or Miniconda may want to check the appendix [FAQ: What if there is no Conda? Miniconda installation](#overview_miniconda).

> [!TIP]
> We provide a Docker image named *whovian-cocoa* that provides cocoa pre-installed and pre-compiled. For further instructions, refer to the appendix [FAQ: How do you run Cocoa on your laptop? The Docker container named *whovian-cocoa*](#appendix_jupyter_whovian).

## Installation and Compilation of external modules <a name="cobaya_base_code"></a>

In this section, we assume users have previously activated the Cocoa conda environment.

**Step :one:**: Download Cocoa's latest release and go to its main folder (`cocoa/Cocoa`),

    git clone https://github.com/CosmoLike/cocoa.git --branch v4.0-beta26 cocoa

and

    cd ./cocoa/Cocoa

**Step :two:**: Run the script `setup_cocoa.sh` via
        
    source setup_cocoa.sh

This script downloads and decompresses external modules, requiring internet access. Therefore, users cannot run this script on an interactive compute node in a HPC environment where only the cluster login node can access the web.

**Step :three:**: Run the script `compile_cocoa.sh` by typing 

    source compile_cocoa.sh
    
This script compiles external modules selected for installation on `set_installation_options.sh` (e.g., CAMB) and does not require internet access. Code compilation is a CPU-intensive operation; therefore, running  `compile_cocoa.sh` on a cluster login node can be against HPC policy. Users should then run `setup_cocoa.sh` in a login node and `compile_cocoa.sh` on an interactive compute node.

Users can now proceed to **the next section**.

> [!TIP]
> If you want to work from the latest commit, then clone the repository with the following command 
>
> (SSH)
> 
>     git clone git@github.com:CosmoLike/cocoa.git cocoa
> 
> (HTTP)
> 
>     git clone https://github.com/CosmoLike/cocoa.git cocoa
>
>
> Users who want to develop from a release version (e.g., `v4.0-beta20`) should read the appendix [FAQ: How do we push changes to the cocoa main branch? A few git hacks](#push_main)

> [!TIP]
> Cocoa does not install all the available external modules by default. If the user requires additional packages, refer to the appendix [Compiling Boltzmann, CosmoLike, and Likelihood codes separately](#appendix_compile_separately).

> [!NOTE]
> In case users need to run `setup_cocoa.sh` more than once, Cocoa will not download previously installed packages, cosmolike projects, or large datasets, unless the following keys are set on `set_installation_options.sh`
>
>     [Adapted from Cocoa/set_installation_options.sh shell script]
>
>     # ------------------------------------------------------------------------------
>     # OVERWRITE_EXISTING_XXX_CODE=1 -> setup_cocoa overwrites existing PACKAGES ----
>     # overwrite: delete the existing PACKAGE folder and install it again -----------
>     # redownload: delete the compressed file and download data again ---------------
>     # These keys are only relevant if you run setup_cocoa multiple times -----------
>     # ------------------------------------------------------------------------------
>     (...)
>     export OVERWRITE_EXISTING_ALL_PACKAGES=1    # except cosmolike projects
>     #export OVERWRITE_EXISTING_COSMOLIKE_CODE=1 # dangerous (possible lost of uncommit work)
>                                                 # if unset, users must manually delete cosmolike projects
>     #export REDOWNLOAD_EXISTING_ALL_DATA=1      # warning: some data are many GB

## Running Examples  <a name="cobaya_base_code_examples"></a>

We assume that you are still in the Conda cocoa environment from the previous `conda activate cocoa` command and that you are in the cocoa main folder `cocoa/Cocoa`, 

 **Step :one:**: Activate the private Python environment by sourcing the script `start_cocoa.sh`

    source start_cocoa.sh

Users will see a terminal like this: `$(cocoa)(.local)`. *This is a feature, not a bug*!

 **Step :two:**: Select the number of OpenMP cores (below, we set it to 8).
    
    export OMP_PROC_BIND=close; export OMP_NUM_THREADS=8

### Examples not involving Cosmolike

 **Step :three:**: The folder `projects/example` contains a few examples involving different likelihoods. So, run the `cobaya-run` on the first example following the commands below.

One model evaluation:

    mpirun -n 1 --oversubscribe --mca pml ^ucx --mca btl vader,tcp,self --bind-to core:overload-allowed --rank-by slot --map-by numa:pe=${OMP_NUM_THREADS} cobaya-run  ./projects/example/EXAMPLE_EVALUATE1.yaml -f
        
MCMC (we run MCMCs with 32 cores):

    mpirun -n 4 --oversubscribe --mca pml ^ucx --mca btl vader,tcp,self --bind-to core:overload-allowed --rank-by slot --map-by numa:pe=${OMP_NUM_THREADS} cobaya-run ./projects/example/EXAMPLE_MCMC1.yaml -f

### Examples involving Cosmolike

 **Step :three:**: The folder `projects/lsst_y1` contains a dozen examples involving different combinations of two-point correlation functions. So, run the `cobaya-run` on the first example following the commands below.

One model evaluation:

    mpirun -n 1 --oversubscribe --mca pml ^ucx --mca btl vader,tcp,self --bind-to core:overload-allowed --rank-by slot --map-by numa:pe=${OMP_NUM_THREADS} cobaya-run ./projects/lsst_y1/EXAMPLE_EVALUATE1.yaml -f


MCMC (we run MCMCs with 32 cores):

    mpirun -n 4 --oversubscribe --mca pml ^ucx --mca btl vader,tcp,self --bind-to core:overload-allowed --rank-by slot --map-by numa:pe=${OMP_NUM_THREADS} cobaya-run ./projects/lsst_y1/EXAMPLE_MCMC1.yaml -f

Profile:

    cd ./projects/lsst_y1

and

    export NMPI=4

and

    mpirun -n ${NMPI} --oversubscribe --mca pml ^ucx --mca btl vader,tcp,self --bind-to core:overload-allowed --rank-by slot --map-by numa:pe=${OMP_NUM_THREADS} python -m mpi4py.futures EXAMPLE_PROFILE1.py --mpi $((${NMPI}-1)) --profile 1 --tol 0.05 --AB 1.0 --outroot 'profile' --minmethod 5 --maxiter 1 --maxfeval 250 

**End of basic instructions**.

> [!Tip]
> Cocoa provides several cosmolike projects, not all of which are installed by default. To activate them, refer to the appendix [Compiling Boltzmann, CosmoLike, and Likelihood codes separately](#appendix_compile_separately).

> [!TIP]
> Assuming Cocoa is installed on a local (not remote!) machine, type the command below after step 2️⃣ to run Jupyter Notebooks.
>
>     jupyter notebook --no-browser --port=8888
>
> The terminal will then show a message similar to the following template:
>
>     (...)
>     [... NotebookApp] Jupyter Notebook 6.1.1 is running at:
>     [... NotebookApp] http://f0a13949f6b5:8888/?token=XXX
>     [... NotebookApp] or http://127.0.0.1:8888/?token=XXX
>     [... NotebookApp] Use Control-C to stop this server and shut down all kernels (twice to skip confirmation).
>
> Now go to the local internet browser and type `http://127.0.0.1:8888/?token=XXX`, where XXX is the previously saved token displayed on the line
> 
>     [... NotebookApp] or http://127.0.0.1:8888/?token=XXX
>
> The project lsst-y1 contains jupyter notebook examples located at `projects/lsst_y1`.

> [!NOTE]
> Why did we choose to work with two separate shell environments, `(cocoa)` and `(.local)`? Our scripts enable users to work on multiple Cocoa instances, similar to what was possible with [CosmoMC](https://github.com/cmbant/CosmoMC). On each instance, our scripts install packages at
>
>      Cocoa/.local/bin
>      Cocoa/.local/include
>      Cocoa/.local/lib
>      Cocoa/.local/share

## Running ML emulators <a name="cobaya_base_code_examples_emul"></a>

Cocoa contains a few transformer-based neural network emulators capable of simulating CMB, cosmolike, matter power spectrum, and distances. We provide a few scripts that exemplify their API. To run them, we assume users have commented out the following lines on `set_installation_options.sh` prior to running the `setup_cocoa.sh` and `compile_cocoa.sh` installation scripts.

      [Adapted from Cocoa/set_installation_options.sh shell script] 
      
      # inset # symbol in the lines below (i.e., unset these environmental keys)
      #export IGNORE_EMULTRF_CODE=1  #SaraivanovZhongZhu (SZZ) transformer-based emul
      #export IGNORE_EMULTRF_DATA=1  #SaraivanovZhongZhu (SZZ) transformer-based emul
      
      #export IGNORE_SIMONS_OBSERVATORY_LIKELIHOOD_CODE=1    # to run EXAMPLE_EVALUATE26.yaml
      #export IGNORE_SIMONS_OBSERVATORY_CMB_DATA=1           # to run EXAMPLE_EVALUATE26.yaml
 
Now, users must follow all the steps below.

 **Step :one:**: Activate the private Python environment by sourcing the script `start_cocoa.sh`

    source start_cocoa.sh

 **Step :two:**: Select the number of OpenMP cores (no need to thread via OpenMP; however, `OMP_NUM_THREADS=3` did reduce the TRF emulator runtime by a factor of two on our test machine).
    
    export OMP_PROC_BIND=close; export OMP_NUM_THREADS=1

 **Step :three:** So, run `cobaya-run` on the first emulator example following the commands below.

One model evaluation:

    mpirun -n 1 --oversubscribe --mca pml ^ucx --mca btl vader,tcp,self --bind-to core:overload-allowed --rank-by slot --map-by numa:pe=${OMP_NUM_THREADS} cobaya-run ./projects/example/EXAMPLE_EVALUATE22.yaml -f
        
MCMC:

    mpirun -n 4 --oversubscribe --mca btl vader,tcp,self --bind-to core:overload-allowed --rank-by slot --map-by numa:pe=${OMP_NUM_THREADS} cobaya-run ./projects/example/EXAMPLE_MCMC22.yaml -f

PolyChord:

    mpirun -n 8 --oversubscribe --mca pml ^ucx --mca btl vader,tcp,self --bind-to core:overload-allowed --rank-by slot --map-by numa:pe=${OMP_NUM_THREADS} cobaya-run ./projects/example/EXAMPLE_POLY22.yaml -f

> [!NOTE]
> What should users do if they don't unset ML-related environment keys before running `setup_cocoa.sh` and `compile_cocoa.sh`, as rerunning these scripts can require a long time? Instead, run the following commands.
>
>      source start_cocoa.sh      # even if (.local) is already active, users must run start_cocoa.sh again to update key values
> 
> and
>
>      source ./installation_scripts/setup_pip_core_packages.sh     # install pip packages required by ML emulators
>      source ./installation_scripts/setup_emultrf.sh               # download emulator codes
>      source ./installation_scripts/unxv_emultrf.sh                # download pre-trained emulators
>
> and
>
>      source ./installation_scripts/setup_simons_observatory.sh       # download Simons Observatory Likelihood code (to run EXAMPLE_EVALUATE26.yaml)
>      source ./installation_scripts/unxv_simons_observatory.sh        # download Simons Observatory Likelihood data  (to run EXAMPLE_EVALUATE26.yaml)
>      source ./installation_scripts/compile_simons_observatory.sh     # compile Simons Observatory Likelihood code     (to run EXAMPLE_EVALUATE26.yaml)
> 
> Finally, rerun all the steps presented in this section, including step one. Users must reload the `(.local)` environment by rerunning `start_cocoa.sh` so Cocoa can create appropriate symlinks that expose the emulators to Cobaya.

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
  
Following best practices, Cocoa scripts download most external modules from their original repositories, including Cobaya, CAMB, Class, Polychord, ACT-DR6, HiLLiPoP, and Lollipop. Although our repository includes a few likelihoods in compressed [xz file format](https://tukaani.org/xz/format.html), we do not want to discourage users from cloning their code/data from their original repositories.  The work of those authors is extraordinary, and users **must cite them** appropriately.

### :interrobang: FAQ: What if installation or compilation goes wrong? <a name="running_wrong"></a>

Below, we present a few suggested steps to debug Cocoa.

**Step :one:**: define the `COCOA_OUTPUT_VERBOSE` and `COSMOLIKE_DEBUG_MODE` flags on `set_installation_options.sh` to obtain a more detailed output, as shown below
  
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

**Step :two:**: restart the Cocoa private environment by rerunning `source start_cocoa.sh` (every time users edit `set_installation_options.sh`, they must reload the `(.local)` environment by rerunning `start_cocoa.sh`).

**Step :three:**: compile the failed package by following the instructions in the appendix [FAQ: How do we compile the Boltzmann, CosmoLike, and Likelihood codes separately](#appendix_compile_separately).

**Step 4️⃣**: rerun `setup_cocoa.sh` and `compile_cocoa.sh` to ensure all packages are installed and compiled correctly.

### :interrobang: FAQ: How do we compile the Boltzmann, CosmoLike, and Likelihood codes separately <a name="appendix_compile_separately"></a>

To avoid excessive compilation or download times during development, users can use scripts located at `Cocoa/installation_scripts/` to download and compile only specific modules (or datasets). To take full advantage of them, users must first unset the appropriate keys on `set_installation_options.sh`, as exemplified below 

     [Adapted from Cocoa/set_installation_options.sh shell script]

     # ------------------------------------------------------------------------------
     # The flags below allow users to skip downloading specific datasets ------------
     # ------------------------------------------------------------------------------
     #export IGNORE_ACTDR6_DATA=1                  # ACT-DR6 likelihood data
     (...)
     #export IGNORE_SIMONS_OBSERVATORY_CMB_DATA=1  # SO likelihood data

     (...)

     # ------------------------------------------------------------------------------
     # The keys below control which packages will be installed and compiled 
     # ------------------------------------------------------------------------------
     #export IGNORE_SIMONS_OBSERVATORY_LIKELIHOOD_CODE=1     # SO likelihood code
     (...)
     #export IGNORE_ACTDR6_CODE=1                            # ACT-DR6 likelihood code

Everytime users edit `set_installation_options.sh`, they need to reload `(.local)` by rerunning `start_cocoa.sh`. Then, users should run the commands below

     cd ./cocoa/Cocoa
     source start_cocoa.sh      # even if (.local) is already active, users must run start_cocoa.sh again to update key values

 and
 
     source ./installation_scripts/setup_act_dr6.sh                # download likelihood code
     source ./installation_scripts/setup_simons_observatory.sh     # download likelihood code

     source ./installation_scripts/compile_act_dr6.sh              # compile likelihood code
     source ./installation_scripts/compile_simons_observatory.sh   # compile likelihood code
     
     source ./installation_scripts/unxv_act_dr6.sh                 # download and unpack likelihood data
     source ./installation_scripts/unxv_simons_observatory.sh      # download and unpack likelihood data

Finally, cocoa also provide several cosmolike projects. To activate them, manipulate the following lines on `set_installation_options.sh` 

     [Adapted from Cocoa/set_installation_options.sh shell script]

     # ------------------------------------------------------------------------------
     # The keys below control which cosmolike projects will be installed and compiled
     # ------------------------------------------------------------------------------
     export IGNORE_COSMOLIKE_LSSTY1_CODE=1
     export IGNORE_COSMOLIKE_DES_Y3_CODE=1
     #export IGNORE_COSMOLIKE_ROMAN_FOURIER_CODE=1
     #export IGNORE_COSMOLIKE_ROMAN_REAL_CODE=1

     (...)

     # ------------------------------------------------------------------------------
     # OVERWRITE_EXISTING_XXX_CODE=1 -> setup_cocoa overwrites existing PACKAGES ----
     # overwrite: delete the existing PACKAGE folder and install it again -----------
     # redownload: delete the compressed file and download data again ---------------
     # These keys are only relevant if you run setup_cocoa multiple times -----------
     # ------------------------------------------------------------------------------
     (...)
     export OVERWRITE_EXISTING_COSMOLIKE_CODE=1 # dangerous (possible lost of uncommit work)
                                                # if unset, users must manually delete cosmolike projects

     (...)
 
     # ------------------------------------------------------------------------------
     # Cosmolike projects below -------------------------------------------
     # ------------------------------------------------------------------------------
     (...)
     export ROMAN_REAL_URL="https://github.com/CosmoLike/cocoa_roman_real.git"
     export ROMAN_REAL_NAME="roman_real"
     #BRANCH: if unset, load the latest commit on the specified branch
     #export ROMAN_REAL_BRANCH="dev"
     #COMMIT: if unset, load the specified commit
     export ROMAN_REAL_COMMIT="a5cf62ffcec7b862dda5bf343bf6bb19124bb5d0"
     #TAG: if unset, load the specified TAG
     #export ROMAN_REAL_TAG="v4.0-beta17"
 
Everytime users edit `set_installation_options.sh`, they need to reload `(.local)` by rerunning `start_cocoa.sh`. Then, run the following commands:

      cd ./cocoa/Cocoa
      source start_cocoa.sh # even if (.local) is already active, users must run start_cocoa.sh again to update key values

 and
 
      source ./installation_scripts/setup_cosmolike_projects.sh   # download all cosmolike projects  
      source ./installation_scripts/compile_all_projects.sh       # compile  all cosmolike project

In case users only want to compile a single cosmolike project (let's say the `roman_real` project)

      source ./projects/roman_real/scripts/compile_roman_real.sh
     
### :interrobang: FAQ: How do we run cocoa on a laptop? The docker image named *whovian-cocoa* <a name="appendix_jupyter_whovian"></a>

We provide the Docker image [whovian-cocoa](https://hub.docker.com/r/vivianmiranda/whovian-cocoa) to facilitate the installation of Cocoa on Windows and macOS. This appendix assumes that users have already installed the Docker Engine on their local PC. For instructions on installing the Docker engine on specific operating systems, refer to [Docker's official documentation](https://docs.docker.com/engine/install/). 

 **Step :one:**: Create a folder and go to the location on the host computer where you want to provide access to the Docker container, as shown below. 

     mkdir -p cocoa_docker
     cd ./cocoa_docker

 **Step :two:**: Download the docker image *whovian-cocoa*, name the associated container `cocoa2023`, and run the container for the first time, type:

    docker run --platform linux/amd64 --hostname cocoa --name cocoa2025 -it -p 8888:8888 -v $(pwd):/home/whovian/host/ -v ~/.ssh:/home/whovian/.ssh:ro vivianmiranda/whovian-cocoa

This is a large image with a size of approximately 16GB, as it already contains the cocoa `v4.0beta26` installed and pre-compiled. 

 **Step :three:**: As shown in the picture below, users can follow the instructions provided in Section [Running Examples](#cobaya_base_code_examples) to run a few non-cosmolike-based examples, as well as examples within `LSST-Y1`, `ROMAN_REAL`, and `ROMAN_FOURIER` projects. 

![screenshot_2025-06-02_at_7 23 48_pm](https://github.com/user-attachments/assets/2b15ce75-3d43-4a65-ab7b-e70077492b32)

> [!TIP]
> Once installation is complete, the user must learn how to **start** and **exit** the Docker container. Assuming the user maintained the container name `cocoa2025` set on the flag `--name cocoa2025`, type:
>    
>      docker start -ai cocoa2025
>
>  to restart the container.

> [!TIP]
> To run Jupyter Notebooks within the *whovian-cocoa* docker container installed on a local machine, type the following command:
>
>      jupyter notebook --no-browser --port=8888
>
> The terminal will show a message similar to the following template:
>
>      [... NotebookApp] Writing notebook server cookie secret to /home/whovian/.local/share/jupyter/runtime/notebook_cookie_secret
>      [... NotebookApp] WARNING: The notebook server is listening on all IP addresses and not using encryption. This is not recommended.
>      [... NotebookApp] Serving notebooks from local directory: /home/whovian/host
>      [... NotebookApp] Jupyter Notebook 6.1.1 is running at:
>      [... NotebookApp] http://f0a13949f6b5:8888/?token=XXX
>      [... NotebookApp] or http://127.0.0.1:8888/?token=XXX
>      [... NotebookApp] Use Control-C to stop this server and shut down all kernels (twice to skip confirmation).
>
> Now go to the local internet browser and type `http://127.0.0.1:8888/?token=XXX`, where XXX is the previously saved token displayed on the line
> 
>     [... NotebookApp] or http://127.0.0.1:8888/?token=XXX  
      
> [!TIP]
> To run the Jupyter Notebook on the *whovian-cocoa* docker container installed in a remote server, adjust the command below
>
>     ssh your_username@your_sever.com -L 8888:localhost:8888
>
> before typing `http://127.0.0.1:8888/?token=XXX` on the local desktop/laptop. This will bind the server port `8888` to the local port `8888`.
> 
>     [... NotebookApp] or http://127.0.0.1:8888/?token=XXX  

> [!TIP]
> To delete a particular container, assuming the container name `cocoa2025`, type
>
>      docker rm -f cocoa2025
>

> [!Tip]
> The flag `-v $(pwd):/home/whovian/host/` in the `docker run` command ensures that files in the host computer located within the folder where the Docker container was initialized are accessible inside the container.

> [!Warning]
> Do not allow the Docker container to have system-wide access to your files. Accidents happen, especially when dealing with dangerous bash commands such as `rm` (deletion).

### :interrobang: FAQ: How do we use an available Anaconda module on HPC? <a name="overview_anaconda"></a>

Below, we list the most common issues users encounter when installing Cocoa conda environments in a supercomputer environment using a globally defined Anaconda module. 

- :interrobang: **Conda command not found**.
  
Anaconda is not usually set by default on HPC environments, but may be available as a module. For example, on the Midway HPC cluster, it can be loaded using the command below.

    module load Anaconda3/2022.10

If you are not sure about the name of the available Anaconda module, type the command

    module avail An

to show all modules with names that start with `An`. The output should resemble the following.

<img width="700" alt="Anaconda" src="https://github.com/user-attachments/assets/09326f5f-49e0-45b5-a157-25fe2b09918e">

- :interrobang: **Installation seems to take forever**.

There are several reasons why installing the Cocoa conda environment may take a long time. Here is a checklist of best practices for troubleshooting installations.

:one: *Never install conda environments using the login node*. 

Instead, request an interactive job with a few cores. However, users should be aware that some supercomputers do not provide internet access on computing nodes. Ask the HPC staff for a queue dedicated to installing and compiling code; it should be part of a well-designed HPC environment.

- :interrobang: **Conda installation is interrupted due to quota limitations**.

Supercomputers usually enforce strict quota limits on home folders. These limits apply to the total file size and the number of files. By default, Anaconda modules install new environments at `$HOME/.conda/envs`. Anaconda also stores Gigabytes of downloaded packages in the `$HOME/.conda/pkgs` folder; the `pkgs` folder is used by Anaconda as a package cache.

:one: Create an Anaconda folder in a partition containing a more tolerant quota restriction than `$HOME`. For instance, we used the command below on the Midway supercomputer to create the Anaconda folder in the KICP projects partition.

    mkdir /project2/kicp/XXX/anaconda/

:two: Set the `pkgs` package cache folder to `anaconda/pkgs`.

    conda config --add pkgs_dirs /project2/kicp/XXX/anaconda/pkgs

3️⃣: Set the env folder to `anaconda/envs/` 

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

:six: After completing steps :one:-:five:, check the `$HOME/.condarc` file with the command

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

### :interrobang: FAQ: What if there is no Conda? Miniconda installation <a name="overview_miniconda"></a>

**Step :one:**: Download and run the Miniconda installation script. 

    export CONDA_DIR="/gpfs/home/XXX/miniconda"
    
    mkdir "${CONDA_DIR:?}"
    
    wget https://repo.anaconda.com/miniconda/Miniconda3-py39_25.3.1-1-Linux-x86_64.sh
    
    /bin/bash Miniconda3-py39_25.3.1-1-Linux-x86_64.sh -f -b -p "${CONDA_DIR:?}"

    /bin/bash 

Please don't forget to adapt the path assigned to `CONDA_DIR` in the command above:

**Step :two:**: After installation, users must source the conda configuration file, as shown below:

    source $CONDA_DIR/etc/profile.d/conda.sh \
          && conda config --set auto_update_conda false \
          && conda config --set show_channel_urls true \
          && conda config --set auto_activate_base false \
          && conda config --prepend channels conda-forge \
          && conda config --set channel_priority strict \
          && conda init bash

**Step :three:**: After running this command, you will see a message in the terminal that ends with the statement *For changes to take effect, close and re-open your current shell*. Then, type

    source ~/.bashrc

After that, the `conda` command will be available.

### :interrobang: FAQ: How do we set the Slow/Fast decomposition on MCMC chains with Cosmolike? Manual Blocking <a name="manual_blocking_cosmolike"></a>

The Cosmolike Weak Lensing pipelines contain parameters with different speed hierarchies. For example, Cosmolike execution time is reduced by approximately 50% when fixing the cosmological parameters. When varying only the multiplicative shear calibration, Cosmolike execution time is reduced by two orders of magnitude. 

Cobaya cannot automatically handle parameters associated with the same likelihood with different speed hierarchies. Luckily, we can manually impose the speed hierarchy in Cobaya using the `blocking:` option. The only drawback of this method is that parameters of all adopted likelihoods, not only the ones required by Cosmolike, must be manually specified.

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

### :interrobang: FAQ: How do users set the environment for Machine Learning projects? <a name="ml_emulators"></a>

Commenting out the environmental flags shown below, located at *set_installation_options* script, will enable the installation of machine-learning-related libraries via pip.  

    [Adapted from Cocoa/set_installation_options.sh shell script] 
     
    # ------------------------------------------------------------------------------
    # If not set, Cocoa/installation_scripts/setup_pip_core_packages.sh will install several ML packages
    # ------------------------------------------------------------------------------
    #export IGNORE_EMULATOR_GPU_PIP_PACKAGES=1
              
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

### :interrobang: FAQ: How do we push changes to the Cocoa main branch? A few git hacks <a name="push_main"></a>

Until recently, Cocoa development was a bit unstructured. Developers could push directly to the `main` branch, and small commits were not discouraged. Such flexible development rules will soon change when `v4.0` leaves the beta phase. We will protect the `main` branch by requiring every push to be reviewed by Cocoa's leading developers. Our new philosophy establishes that *a commit in the main branch should contain an atomic change that takes code from one working state to another working state with meaningful and well-tested improvements*. Therefore, developers should propose changes to the `main` branch in larger chunks (*via squash commits*), as shown below.

- :interrobang: **How to apply squash commits?**
  
**Step :one:**: create a development branch. Do not call the development branch `dev`, as `dev` is reserved for work done by the leading Cocoa developers. For concreteness, let's name the new branch `xyzdev`

    # Developers must input the command below on the main branch.
    git switch -c xyzdev

> [!TIP]
> The use of developers' initials followed by `dev` helps make the branch easily identifiable.

> [!TIP]
> If the branch `xyzdev` already exists, use the `git switch` command without the `-c` flag. 

**Step :two:**: develop the proposed changes. We advise developers to commit frequently. In your branch, a commit does not need to be atomic, changing the code from one working state to another well-tested, meaningful working state. Developers can push to the server via the command.

    git push -u origin xyzdev   # run on the xyzdev branch

**Step :three:**: Once the developers created an atomic, meaningful, and well-tested improvement to Cocoa, the developer needs to merge any subsequent changes made in `main`.

    git merge main              # run on the xyzdev branch

This step may create conflicts that must be addressed before step four. 

**Step :four:**: Once the developers have merged recent changes made on the `main` branch, they must push to the main branch the modifications made on the `xyzdev` branch by first **squashing all your changes into a single commit**, as shown below

    git switch main             # run on the xyzdev branch

and

    git merge --squash xyzdev   # run on the main branch

and

    git commit -m "merge xyzdev branch"  # run on the main branch

and

    git push origin main # run on the main branch

### :interrobang: FAQ: How do we develop from a git tag? A few more git hacks <a name="dev_from_tag"></a>

A useful git hack is related to developing Cocoa from a git tag. We reserve Git tags to set milestones in our development, so they serve as good starting points for coding localized new features (e.g., changes to a file that other developers have not recently modified) or bug fixes.

**Step :one: (optional)** If the developer has cloned the repository using the `https` URL address, we change the URL to the SSH-key-based address

    git remote set-url origin git@github.com:CosmoLike/cocoa.git

**Step :two:** Move the detached state to a new local branch via the command

    git switch -c xyzlocdev

Now, all commits will be associated with this local branch. 

> [!TIP]
> The use of developers' initials followed by `dev` helps make the branch easily identifiable.

**Step :three:** The developer has two options at the end of development. They can **either** create a new remote branch

    git push origin xyzlocdev     # run on the xyzlocdev branch

**or** they can fetch and download the remote `xyzdev` branch, which will later absorb the changes made on `xyzlocdev`

    git switch -c xyzdev origin/xyzdev  # run on the xyzlocdev branch

Finally, the developer needs to merge the changes made on `xyzlocdev`.

    git merge --squash xyzlocdev  # run on the xyzdev branch

If this merge does not create any merge conflicts, type

    git push origin xyzdev        # run on the xyzdev branch

### :interrobang: FAQ: How do we download and run Cosmolike projects? <a name="running_cosmolike_projects"></a> 

In the instructions below, we assume that you are still in the Conda cocoa environment from the previous `conda activate cocoa` command and that you are in the cocoa main folder `cocoa/Cocoa`, 

**Step :one:**: Go to the projects folder (`./projects`) and clone a Cosmolike project with the fictitious name `XXX`:
    
    cd ./projects

and 

    "${CONDA_PREFIX}/bin/git" clone https://github.com/CosmoLike/cocoa_lsst_XXX.git XXX

By convention, the Cosmolike Organization hosts a Cobaya-Cosmolike project named XXX at `CosmoLike/cocoa_XXX`. The `cocoa_` prefix must be dropped when cloning the repository.
 
**Step :two:**: Go back to the Cocoa main folder and activate the private Python environment
    
    cd ../

and

    source start_cocoa.sh
 
> [!Warning]
> Users must run `start_cocoa.sh` after cloning the project repository, so Cocoa can create appropriate soft-links.

**Step :three:**: Compile the project, as shown below (two possibilities)
 
    source "${ROOTDIR:?}"/projects/XXX/scripts/compile_XXX

or

    # this will compile all cosmolike courses
    source "${ROOTDIR:?}"/installation_scripts/compile_all_projects.sh
 
**Step :four:**: Select the number of OpenMP cores and run a template YAML file
    
    export OMP_PROC_BIND=close; export OMP_NUM_THREADS=8

and 

    mpirun -n 1 --oversubscribe --mca pml ^ucx --mca btl vader,tcp,self --bind-to core:overload-allowed --rank-by slot --map-by numa:pe=${OMP_NUM_THREADS} cobaya-run ./projects/XXX/EXAMPLE_EVALUATE1.yaml -f

If users want to make a particular cosmolike project widely available in Cocoa, implement the following steps:

**Step :one:**: Add the following env keys on `set_installation_options.sh`

    [adapted from ${ROOTDIR:?}/set_installation_options.sh script]
    
    #the flag below allows users to skip the downloading project XXX
    #export IGNORE_COSMOLIKE_XXX_CODE=1

    (...)
   
    export XXX_URL="https://github.com/.../cocoa_lsst_XXX.git"
    export XXX_NAME="XXX"
    #Key XXX_COMMIT is optional, but we strongly recommend its inclusion 
    export XXX_COMMIT="a5cf62ffcec7b862dda5b1234f6bb19124bb5d0"

**Step :two:**: Adapt and add the keys below to `flags_impl_unset_keys.sh` 

    [adapted from ${ROOTDIR:?}/installation_scripts/flags_impl_unset_keys.sh]

    unset -v XXX_URL XXX_NAME XXX_COMMIT IGNORE_COSMOLIKE_XXX_CODE

This will ensure that `stop_cocoa.sh` unsets them before exiting Cocoa.

**Step :three:** Add and adapt the following block to `setup_cosmolike_projects.sh`.

     [adapted from ${ROOTDIR:?}/installation_scripts/setup_cosmolike_projects.sh script]
     
     if [ -z "${IGNORE_COSMOLIKE_XXX_CODE}" ]; then 
       ptop "GETTING XXX" || return 1

       if [ -n "${XXX_COMMIT}" ]; then
         gitact2 "${XXX_NAME:?}" "${XXX_URL:?}" "${XXX_COMMIT:?}"  || return 1
       fi

       pbottom "GETTING XXX" || return 1
     fi

> [!NOTE]
> *projects* was designed to include all Cosmolike projects, and Cocoa contains two scripts
>
>      "${ROOTDIR:?}"/installation_scripts/setup_cosmolike_projects.sh
>      "${ROOTDIR:?}"/installation_scripts/compile_all_projects.sh
> 
> that `setup` and `compile` all projects defined there. As a standard, we defined the project name `XXX` to be stored on a GitHub repository with the name `cocoa_XXX`. The prefix `cocoa_` helps developers in the Cosmolike organization differentiate legacy Cosmolike projects from matching ones designed for Cocoa. 

> [!Warning]
> Never delete a folder from `projects` without first running `stop_cocoa.sh`; otherwise, Cocoa will have ill-defined links.

### :interrobang: FAQ: How do we add a new package to Cocoa? The Dark Emulator Example <a name="add_package_v1"></a> 

Cocoa provides a verbose but methodical method for adding packages to its environment. Science packages that are still in development by the authors should be cloned in development mode with commit flags clearly delimited. This includes pure Python packages, which are typically installed using the pip command. Every package should have a `setup` and a `compile` script, and the' compile' script should never depend on an internet connection. As an example, we list below the steps we implemented to add the Dark Emulator package to Cocoa. 

