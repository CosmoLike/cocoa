
<img width="780" alt="Screenshot 2024-09-16 at 10 55 14 PM" src="https://github.com/user-attachments/assets/f327e3bf-22c7-46e4-9bb5-14c8e7afc4c1">

# Table of contents
1. [Overview of the Cobaya-CosmoLike Joint Architecture (Cocoa)](#overview) 
2. [Installation of core packages via Conda](#required_packages_conda)
3. [Installation and Compilation of external modules](#cobaya_base_code)
4. [Running Examples](#cobaya_base_code_examples)
5. [Running Examples based on Machine Learning emulators](#cobaya_base_code_examples_emul)
6. [Creating Cosmolike projects (external readme)](Cocoa/projects/)
7. [Appendix](#appendix)
    1. [Credits](#appendix_proper_credits)
    2. [Additional Installation Notes](#additional_notes)
    3. [FAQ: What if installation or compilation goes wrong?](#running_wrong)
    4. [FAQ: How do we compile the Boltzmann, CosmoLike, and Likelihood codes separately](#appendix_compile_separately)
    5. [FAQ: How do we run Cocoa on a laptop? The docker image named *whovian-cocoa*](#appendix_jupyter_whovian)
    6. [FAQ: How do we use an available Anaconda module on HPC?](#overview_anaconda)
    7. [FAQ: What if there is no Conda? Miniconda installation](#overview_miniconda)
    8. [FAQ: How do we make the Slow/Fast decomposition on MCMC chains with Cosmolike? Manual Blocking](#manual_blocking_cosmolike)
    9. [FAQ: How do we switch Cocoa's adopted CAMB/CLASS/Polychord? (external readme)](Cocoa/external_modules/code)
    10. [FAQ: How do we download modern CMB data? (external readme)](Cocoa/external_modules/data)
    11. [FAQ: How do we set the environment for Machine Learning projects?](#ml_emulators)
    12. [FAQ: How can users improve their Bash/C/C++ knowledge to develop Cocoa/Cosmolike?](#lectnotes)
    13. [Warning about Weak Lensing YAML files in Cobaya](#appendix_example_runs)
    14. [FAQ: How do we install Cocoa without conda?](#required_packages_cache)
    15. [FAQ: How do we push changes to the Cocoa main branch? A few git hacks](#push_main)
    16. [FAQ: How do we develop from a git tag? A few more git hacks](#dev_from_tag)
    17. [FAQ: How do we download and run Cosmolike projects?](#running_cosmolike_projects)
    18. [FAQ: How do we add a new package to Cocoa? The Dark Emulator Example](#add_package_v1)
    19. [FAQ: How do we add a new package to Cocoa? The MGCAMB Example](#add_package_v2)
    
## Overview of the [Cobaya](https://github.com/CobayaSampler)-[CosmoLike](https://github.com/CosmoLike) Joint Architecture (Cocoa) <a name="overview"></a>

Cocoa allows users to run [CosmoLike](https://github.com/CosmoLike) routines inside the [Cobaya](https://github.com/CobayaSampler) framework. [CosmoLike](https://github.com/CosmoLike) can analyze data primarily from the [Dark Energy Survey](https://www.darkenergysurvey.org) and simulate future multi-probe analyses for LSST and Roman Space Telescope. 

Besides integrating [Cobaya](https://github.com/CobayaSampler) and [CosmoLike](https://github.com/CosmoLike), Cocoa introduces shell scripts that allow users to containerize [Cobaya](https://github.com/CobayaSampler), the Boltzmann codes, and multiple likelihoods. The container structure of Cocoa ensures that users can adopt consistent versions for the Fortran/C/C++ compilers and libraries across multiple machines. Such a systematic approach greatly simplifies the debugging process. 

Our scripts never install packages and Python modules on `$HOME/.local`, where `$HOME` is a shell environment variable that points to the user's home folder, as that would make them global to the user. Such behavior would interfere with all previously installed packages, potentially creating dependency incompatibilities between different projects the user works on. This behavior enables users to work on multiple instances of Cocoa simultaneously, similar to what was possible with [CosmoMC](https://github.com/cmbant/CosmoMC). 

This readme file presents basic and advanced instructions for installing all [Cobaya](https://github.com/CobayaSampler) and [CosmoLike](https://github.com/CosmoLike) components.

## Installation of core packages via Conda <a name="required_packages_conda"></a>

Core packages include compilers and numerical libraries (e.g., GSL and FFTW), which users typically do not modify. We install most of these core packages via Conda, as shown below.

**Step :one:**: Download the file `cocoapy310.yml` yml file

    wget https://raw.githubusercontent.com/CosmoLike/cocoa/refs/heads/main/cocoapy310.yml

create the cocoa environment,

    conda env create --name cocoa --file=cocoapy310.yml

activate it

    conda activate cocoa

and create symbolic links that will give better names for the GNU compilers

    ln -s "${CONDA_PREFIX}"/bin/x86_64-conda_cos6-linux-gnu-gcc "${CONDA_PREFIX}"/bin/gcc
    ln -s "${CONDA_PREFIX}"/bin/x86_64-conda_cos6-linux-gnu-g++ "${CONDA_PREFIX}"/bin/g++
    ln -s "${CONDA_PREFIX}"/bin/x86_64-conda_cos6-linux-gnu-gfortran "${CONDA_PREFIX}"/bin/gfortran
    ln -s "${CONDA_PREFIX}"/bin/x86_64-conda-linux-gnu-gcc-ar "${CONDA_PREFIX}"/bin/gcc-ar
    ln -s "${CONDA_PREFIX}"/bin/x86_64-conda-linux-gnu-gcc-ranlib "${CONDA_PREFIX}"/bin/gcc-ranlib
    ln -s "${CONDA_PREFIX}"/bin/x86_64-conda-linux-gnu-ld "${CONDA_PREFIX}"/bin/ld

Users can now proceed to **step :two:**. 

> [!TIP]
> To install the `cocoa` conda environment on a supercomputer, users may take advantage of the fact that many HPC environments provide the [Anaconda installer](https://www.anaconda.com) as an external module. Check the appendix [FAQ: How do we use an available Anaconda module on HPC?](#overview_anaconda). Conversely, users working on an HPC environment that does not offer Anaconda or Miniconda may want to check the appendix [FAQ: What if there is no Conda? Miniconda installation](#overview_miniconda).

**Step :two:**: Install `git-lfs` when loading the conda cocoa environment for the first time.

    git-lfs install

## Installation and Compilation of external modules <a name="cobaya_base_code"></a>

> [!WARNING]
> In this section, we assume users have previously activated the Cocoa conda environment
>
> 
**Step :one:**: Download Cocoa's latest release and go to its main folder (`cocoa/Cocoa`),

    "${CONDA_PREFIX}"/bin/git clone --depth 1 https://github.com/CosmoLike/cocoa.git --branch v4.0-beta25 cocoa

and

    cd ./cocoa/Cocoa

Users can now proceed to **step :two:**.

> [!TIP]
> If you want to work from the latest commit, then clone the repository with the following command 
>
> (SSH)
> 
>     "${CONDA_PREFIX}"/bin/git clone git@github.com:CosmoLike/cocoa.git cocoa
> 
> (HTTP)
> 
>     "${CONDA_PREFIX}"/bin/git clone https://github.com/CosmoLike/cocoa.git cocoa
>
>
> Users who want to develop from a release version (e.g., `v4.0-beta20`) may want to read the appendix [FAQ: How do we push changes to the cocoa main branch? A few git hacks](#push_main)

**Step :two:**: Run the script `setup_cocoa.sh` via
        
    source setup_cocoa.sh

This script downloads and decompresses external modules, requiring internet access to run successfully.

> [!NOTE]
> If you run `setup_cocoa.sh` multiple times, Cocoa will not download previously installed packages. To overwrite this behavior, export the key `OVERWRITE_EXISTING_ALL_PACKAGES` on `set_installation_options.sh`. Even with this optimization disabled, Cocoa will not download large datasets repeatedly unless the key `REDOWNLOAD_EXISTING_ALL_DATA` is also set.

**Step :three:**: Run the script `compile_cocoa.sh` by typing 

    source compile_cocoa.sh
    
This script compiles external modules selected for installation on `set_installation_options.sh` (e.g., CAMB and Class). Cocoa does not install many external modules by default, but users may find them helpful in a particular project. In this case, check the many available options on the `set_installation_options.sh` shell script. Then, rerun steps :two: and :three:. 

> [!NOTE]
> In some HPC environments, the compute nodes cannot access the web. So, by design, the script `compile_cocoa.sh` does not require internet access to run successfully. Code compilation is a CPU-intensive operation, so running  `compile_cocoa.sh` on a cluster login node can be against HPC policy. Users should then run `setup_cocoa.sh` in a login node and `compile_cocoa.sh` in a compute node.

## Running Examples  <a name="cobaya_base_code_examples"></a>

We assume that you are still in the Conda cocoa environment from the previous `conda activate cocoa` command and that you are in the cocoa main folder `cocoa/Cocoa`, 

 **Step :one:**: Activate the private Python environment by sourcing the script `start_cocoa.sh`

    source start_cocoa.sh

Users will see a terminal like this: `$(cocoa)(.local)`. *This is a feature, not a bug*! 

 **Step :two:**: Select the number of OpenMP cores (below, we set it to 8).
    
    export OMP_PROC_BIND=close; export OMP_NUM_THREADS=8

### Examples not involving Cosmolike

 **Step :three:**: The folder `projects/example` contains a dozen examples involving different likelihoods. So, run the `cobaya-run` on the first example following the commands below.

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

> [!Note]
> We provide several cosmolike projects that can be loaded and compiled using `setup_cocoa.sh` and `compile_cocoa.sh` scripts. To activate them, comment the following lines on `set_installation_options.sh` 
> 
>     [Adapted from Cocoa/set_installation_options.sh shell script]
>     (...)
>
>     # ------------------------------------------------------------------------------
>     # The keys below control which cosmolike projects will be installed and compiled
>     # ------------------------------------------------------------------------------
>     export IGNORE_COSMOLIKE_LSSTY1_CODE=1
>     export IGNORE_COSMOLIKE_DES_Y3_CODE=1
>     #export IGNORE_COSMOLIKE_ROMAN_FOURIER_CODE=1
>     #export IGNORE_COSMOLIKE_ROMAN_REAL_CODE=1
>
>     (...)
>
>     # ------------------------------------------------------------------------------
>     # OVERWRITE_EXISTING_XXX_CODE=1 -> setup_cocoa overwrites existing PACKAGES ----
>     # overwrite: delete the existing PACKAGE folder and install it again -----------
>     # redownload: delete the compressed file and download data again ---------------
>     # These keys are only relevant if you run setup_cocoa multiple times -----------
>     # ------------------------------------------------------------------------------
>     (...)
>     export OVERWRITE_EXISTING_COSMOLIKE_CODE=1 # dangerous (possible lost of uncommit work)
>                                                # if unset, users must manually delete
>                                                # project if wants setup_cocoa to reclone it
>
>     (...)
> 
>     # ------------------------------------------------------------------------------
>     # Cosmolike projects below -------------------------------------------
>     # ------------------------------------------------------------------------------
>     (...)
>     export ROMAN_REAL_URL="https://github.com/CosmoLike/cocoa_roman_real.git"
>     export ROMAN_REAL_NAME="roman_real"
>     #BRANCH: if unset, load the latest commit on the specified branch
>     #export ROMAN_REAL_BRANCH="dev"
>     #COMMIT: if unset, load the specified commit
>     export ROMAN_REAL_COMMIT="a5cf62ffcec7b862dda5bf343bf6bb19124bb5d0"
>     #TAG: if unset, load the specified TAG
>     #export ROMAN_REAL_TAG="v4.0-beta17"
> 
> If users comment these lines (unsetting the corresponding IGNORE keys) after running `setup_cocoa.sh` and `compile_cocoa.sh`, there is no need to rerun these general scripts, which would reinstall many packages (slow). Instead, run the following three commands:
>
>      source start_cocoa.sh
>
> and
> 
>      source ./installation_scripts/setup_cosmolike_projects.sh
>
> and
> 
>       source ./installation_scripts/compile_all_projects.sh
>
> or in case users just want to compile a single project (let's say the `roman_real` project)
>
>       source ./projects/roman_real/scripts/compile_roman_real.sh
  
> [!TIP]
> To run Jupyter Notebook, assuming Cocoa is installed on a local machine, type, after step 2️⃣, the command 
> 
>     jupyter notebook --no-browser --port=8888
>
> The terminal will show a message similar to the following template:
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

## Running Examples based on Machine Learning emulators <a name="cobaya_base_code_examples_emul"></a>

Cocoa contains a few transformer-based neural network emulators capable of simulating CMB, cosmolike, matter power spectrum, and distances. We provide a few scripts that exemplify their API. To run them, we assume users have commented out the following lines on `set_installation_options.sh` prior to running the `setup_cocoa.sh` and `compile_cocoa.sh` installation scripts.

      [Adapted from Cocoa/set_installation_options.sh shell script] 
      
      # inset # symbol in the lines below (i.e., unset these environmental keys)
      #export IGNORE_EMULTRF_CODE=1  #SaraivanovZhongZhu (SZZ) transformer-based emul
      #export IGNORE_EMULTRF_DATA=1  #SaraivanovZhongZhu (SZZ) transformer-based emul

 **Step :one:**: Activate the private Python environment by sourcing the script `start_cocoa.sh`

    source start_cocoa.sh

 **Step :two:**: Select the number of OpenMP cores (no need to thread via OpenMP).
    
    export OMP_PROC_BIND=close; export OMP_NUM_THREADS=1

 **Step :three:** So, run `cobaya-run` on the first emulator example following the commands below.

One model evaluation:

    mpirun -n 1 --oversubscribe --mca pml ^ucx --mca btl vader,tcp,self --bind-to core:overload-allowed --rank-by slot --map-by numa:pe=${OMP_NUM_THREADS} cobaya-run ./projects/example/EXAMPLE_EVALUATE22.yaml -f
        

MCMC (we run MCMCs with 12 cores):

    mpirun -n 4 --oversubscribe --mca btl vader,tcp,self --bind-to core:overload-allowed --rank-by slot --map-by numa:pe=${OMP_NUM_THREADS} cobaya-run ./projects/example/EXAMPLE_MCMC22.yaml -f


PolyChord (we run Nested Sampling with 24 cores):

    mpirun -n 8 --oversubscribe --mca pml ^ucx --mca btl vader,tcp,self --bind-to core:overload-allowed --rank-by slot --map-by numa:pe=${OMP_NUM_THREADS} cobaya-run ./projects/example/EXAMPLE_POLY22.yaml -f
    
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

### Additional Installation Notes <a name="additional_notes"></a>

- *Installation of core packages via Conda* 
 
> [!TIP]
> For those working on projects that utilize machine-learning-based emulators, the appendix [Setting-up conda environment for Machine Learning emulators](#ml_emulators) provides additional commands for installing the necessary packages.

> [!TIP]
> We provide a Docker image named *whovian-cocoa* that facilitates cocoa installation on Windows and MacOS. For further instructions, refer to the appendix [FAQ: How do you run Cocoa on your laptop? The Docker container is named *whovian-cocoa*](#appendix_jupyter_whovian).

> [!NOTE]
> The conda installation method should be chosen in most cases. In the rare instances in which users cannot work with Conda, refer to the appendix [Installation of Cocoa's core packages without Conda](#required_packages_cache), as it contains instructions for a much slower (and prone to errors) but Conda-independent installation method.

- *Installation and Compilation of external modules* 

> [!TIP]
> If users want to compile only a subset of these packages, refer to the appendix [Compiling Boltzmann, CosmoLike, and Likelihood codes separately](#appendix_compile_separately).
          
- *Running Examples*

> [!NOTE]
> `OMP_PROC_BIND=close` bound OpenMP threads to physically close cores (within the same chiplet on chiplet-based architectures).

> [!NOTE]
> Additional explanations about our `mpirun` flags: Why the `--bind-to core:overload-allowed --map-by numa:pe=${OMP_NUM_THREADS}` flag? This flag enables efficient hybrid MPI+OpenMP runs on NUMA architecture.

> [!NOTE]
> Additional explanations about the functioning `start_cocoa.sh`/`stop_cocoa.sh` scripts: why two separate shell environments, `(cocoa)` and `(.local)`? Users should be able to manipulate multiple Cocoa instances seamlessly, which is particularly useful when running chains in one instance while experimenting with code development in another. Consistency of the environment across all Cocoa instances is crucial, and the `start_cocoa.sh`/`stop_cocoa.sh` scripts handle the loading and unloading of environmental path variables. Our scripts never install packages on `$HOME/.local` as that would make them global to the user. Instead, on each instance, they are installed at
>
>      Cocoa/.local/bin
>      Cocoa/.local/include
>      Cocoa/.local/lib
>      Cocoa/.local/share
>
> This behavior enables users to work on multiple instances of Cocoa simultaneously, similar to what was possible with [CosmoMC](https://github.com/cmbant/CosmoMC).

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
      export IGNORE_ACTDR6_DATA=1
      # export IGNORE_BAO_DATA=1
      export IGNORE_BICEP_CMB_DATA=1
      # export IGNORE_HOLICOW_STRONG_LENSING_DATA=1
      # export IGNORE_SN_DATA=1
      export IGNORE_SPT_CMB_DATA=1
      export IGNORE_SIMONS_OBSERVATORY_CMB_DATA=1
      #export IGNORE_PLANCK_CMB_DATA=1
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
      export IGNORE_HYREC_CODE=1
      export IGNORE_COSMOREC_CODE=1
      #export IGNORE_EUCLID_EMULATOR_V2_CODE=1
      #export IGNORE_COSMOLIKE_LSSTY1_CODE=1

      (...)

      # ------------------------------------------------------------------------------
      # OVERWRITE_EXISTING_XXX_CODE=1 -> setup_cocoa overwrites existing PACKAGES ----
      # overwrite means: delete the existing PACKAGE folder and install it again -----------
      # redownload: delete the compressed file and download the data again ---------------
      # These keys are only relevant if you run setup_cocoa multiple times -----------
      # ------------------------------------------------------------------------------
       export OVERWRITE_EXISTING_ALL_PACKAGES=1      # except cosmolike projects
       #export OVERWRITE_EXISTING_COSMOLIKE_CODE=1   # dangerous (possible lost of uncommit work)
       #export REDOWNLOAD_EXISTING_ALL_DATA=1        # warning: some data are many GB

Steps to debug Cocoa

- The first step is to define the `COCOA_OUTPUT_VERBOSE` and `COSMOLIKE_DEBUG_MODE` flags to obtain a more detailed output. To accomplish that, we advise users to uncomment the lines below that are part of the `set_installation_options.sh` script and then restart the Cocoa private environment by running `source stop_cocoa.sh; source start_cocoa.sh`

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

- The second step involves rerunning the failed script with the verbose output set. The scripts `setup_cocoa.sh` and `compile_cocoa.sh` run many shell scripts, so users may find it advantageous to run only the routine that failed. For further information on how to do that, see the appendix [FAQ: How do we compile the Boltzmann, CosmoLike, and Likelihood codes separately](#appendix_compile_separately).

After fixing a particular issue, users should rerun the shell scripts `setup_cocoa.sh` and `compile_cocoa.sh` to ensure all packages are installed and compiled correctly.

### :interrobang: FAQ: How do we compile the Boltzmann, CosmoLike, and Likelihood codes separately <a name="appendix_compile_separately"></a>

To avoid excessive compilation or download times during development, users can use scripts located at `Cocoa/installation_scripts/` that compile only a specific module or download only a particular dataset. A few examples of these scripts are: 

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

To ensure these scripts can download and install these packages, users must be sure that the environment keys below are *NOT* set. These keys are shown on `set_installation_options.sh`. The command `unset -v` unsets them. 
      
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

### :interrobang: FAQ: How do we run cocoa on a laptop? The docker image named *whovian-cocoa* <a name="appendix_jupyter_whovian"></a>

We provide the Docker image [whovian-cocoa](https://hub.docker.com/r/vivianmiranda/whovian-cocoa) to facilitate the installation of Cocoa on Windows and macOS. This appendix assumes that users already have the Docker Engine installed on their local PC. For instructions on installing the Docker engine on specific operating systems, please refer to [Docker's official documentation](https://docs.docker.com/engine/install/). 

 **Step :one:**: Create a folder and go to the location on the host computer where you want to provide access to the Docker container, as shown below. 

     mkdir -p cocoa_docker
     cd ./cocoa_docker

> [!NOTE]
> The flag `-v $(pwd):/home/whovian/host/` in the `docker run` command ensures that files on the host computer have been mounted to the directory `/home/whovian/host/`. Files within the folder where the Docker container was initialized are accessible at `/home/Whovian/host/`. Users should work inside this directory to avoid losing work if the Docker image needs to be deleted.

> [!WARNING]
>  Do not run the Docker container on a general folder (like the host's home directory); this would provide too much access to the Docker container. Accidents happen, especially when dealing with dangerous bash commands such as `rm` (deletion).

 **Step :two:**: Download the docker image *whovian-cocoa*, name the associated container `cocoa2023`, and run the container for the first time, type:

    docker run --platform linux/amd64 --hostname cocoa --name cocoa2023 -it -p 8888:8888 -v $(pwd):/home/whovian/host/ -v ~/.ssh:/home/whovian/.ssh:ro vivianmiranda/whovian-cocoa

Following the command above, users should see the following text on the screen terminal:

<img width="872" alt="Screenshot 2023-12-19 at 1 26 50 PM" src="https://github.com/CosmoLike/cocoa/assets/3210728/eb1fe7ec-e463-48a6-90d2-2d84e5b61aa1">

 **Step :three:**: Init Conda when running the container the first time, as shown below.

    conda init bash
    source ~/.bashrc
    conda activate cocoa
    git-lfs install

 **Step :four:**: Access the host computer on `/home/whovian/host/`
 
    cd ~/host/

Now proceed with the standard cocoa installation in section [Installation and Compilation of external modules](#cobaya_base_code)
 
Once installation is complete, the user must learn how to **start**, use, and **exit** the container. Below, we address a few common questions about using and managing Docker containers.  

> [!TIP]
> Assuming the user maintained the container name `cocoa2023` via the flag `--name cocoa2023` on the `docker run` command, type:
>    
>      docker start -ai cocoa2023
>
>  to restart the container after the first exit.

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

### :interrobang: FAQ: How do we use an available Anaconda module on HPC? <a name="overview_anaconda"></a>

Below, we list users' most common issues when installing Cocoa conda environments in a supercomputer environment using a globally defined Anaconda module. 

- :interrobang: **Conda command not found**.
  
Anaconda is not usually set by default on HPC environments, but may be available as a module. For example, on the Midway HPC cluster, it can be loaded using the command below.

    module load Anaconda3/2022.10

If you are not sure about the name of the available Anaconda module, type the command

    module avail An

to show all modules with names that start with `An`. The output should resemble the following.

<img width="700" alt="Anaconda" src="https://github.com/user-attachments/assets/09326f5f-49e0-45b5-a157-25fe2b09918e">

- :interrobang: **Installation seems to take forever**.

There are various reasons why installing the Cocoa conda environment may take a long time. Here is a checklist of best practices to troubleshoot installation.

:one: *Never install conda environments using the login node*. 

Instead, request an interactive job with a few cores. However, users should be aware that **some supercomputers do not provide internet access on computing nodes**. Ask the HPC staff for a **queue dedicated to installing and compiling code**; it should exist in a well-designed HPC environment.

- :interrobang: **Conda installation is interrupted due to quota limitations**.

Supercomputers usually enforce strict quota limits on home folders. These limits apply to the total file size and the number of files. By default, Anaconda modules install new environments at `$HOME/.conda/envs`. Anaconda also stores Gigabytes of downloaded packages in the `$HOME/.conda/pkgs` folder; the `pkgs` folder is used by Anaconda as a package cache.

:one: Create an Anaconda folder in a project folder outside `$HOME` with significantly more tolerant quota restrictions. For instance, we used the command below on the Midway supercomputer to create an Anaconda folder in the KICP projects partition.

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

The Cosmolike Weak Lensing pipelines contain parameters with different speed hierarchies. For example, Cosmolike execution time is reduced by approximately 50% when fixing the cosmological parameters. When varying only multiplicative shear calibration, Cosmolike execution time is reduced by two orders of magnitude. 

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
            
### 💀 ☠️ :stop_sign::thumbsdown: FAQ: How do we install Cocoa without conda (not recommended) <a name="required_packages_cache"></a>

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

### :interrobang: FAQ: How do we push changes to the Cocoa main branch? A few git hacks <a name="push_main"></a>

Until recently, Cocoa development was unstructured, and developers could push directly to the `main` branch. Small commits were not discouraged, but such flexible development rules will soon change. We will protect the `main` branch by requiring every push to be reviewed by Cocoa's leading developers. We aim to reduce the number of commits in the future. Our new philosophy establishes that **a commit in the main branch should contain an atomic change that takes code from one working state to another working state with meaningful and well-tested improvements.** So, developers should propose changes to the `main` branch in larger chunks (*via squash commits*), as shown below.

Important note: We strongly advise developers to use the up-to-date `git` provided by the Cocoa conda environment. 

- :interrobang: **How to apply squash commits?**
  
**Step :one:**: create a development branch from the `main` branch. Do not call the development branch `dev`, it is reserved for work done by the primary Cocoa developers. Let's call this new branch `xyzdev` for concreteness (tip: the use of developers' initials helps make the branch name unique)

    # Developers must input the command below on the main branch.
    "${CONDA_PREFIX}"/bin/git switch -c xyzdev

If the `xyzdev` already exists, use `git switch` without the `-c` flag. Developers can also push their development branches to the server via the command.

    # Developers must input the command below on the xyzdev branch.
    "${CONDA_PREFIX}"/bin/git push -u origin xyzdev

**Step :two:**: develop the proposed changes. We advise developers to commit frequently. In your branch, a commit does not need to be atomic, changing the code from one working state to another well-tested, meaningful working state. In your branch, you have absolute freedom.

**Step :three:**: Once the developers created an atomic, meaningful, and well-tested improvement to Cocoa, the developer needs to merge any subsequent changes made in `main` while the developer has been working on the `xyzdev` branch.

    # Developers must input the command below on the xyzdev branch.
    "${CONDA_PREFIX}"/bin/git merge main

This step may create conflicts that must be addressed before step four. 

**Step :four:**: Once the developer has merged recent changes made on the `main` branch, the developer must push to the main branch the modifications made on the `xyzdev` branch by first **squashing all your changes into a single commit** as shown below

    # Developers must input the command below on the xyzdev branch.
    "${CONDA_PREFIX}"/bin/git switch main

and

    # Developers must input the command below on the main branch.
    "${CONDA_PREFIX}"/bin/git merge --squash xyzdev

and

    # Developers must input the command below on the main branch.
    "${CONDA_PREFIX}"/bin/git commit -m "squash merge - xyzdev branch: added development on abc features"

and

    # Developers must input the command below on the main branch.
    "${CONDA_PREFIX}"/bin/git push origin main

Important note: **never** revert the branch ordering on squash merging by squashing the `main` changes to the `xyzdev` branch.

### :interrobang: FAQ: How do we develop from a git tag? A few more git hacks <a name="dev_from_tag"></a>

A useful git hack is related to developing Cocoa from a git tag. We reserve git tags to set milestones in our development, so they are good starting points for coding localized new features (e.g., changes to a file that other developers have not recently modified) or bug fixes.

**Step :one: (optional)** If the developer has cloned the repository using the `https` URL address, we change the URL to the SSH-key-based address

    "${CONDA_PREFIX}"/bin/git remote set-url origin git@github.com:CosmoLike/cocoa.git

**Step :two:** Move the detached state to a new local branch via the command

    "${CONDA_PREFIX}"/bin/git switch -c mylocalbranch

Now, all commits will be associated with this local branch. The developer can then make a series of git commits to implement a new feature or fix a bug.

**Step :three:** The developer has two options at the end of development. They can **either** create a new remote branch named `mylocalbranch`

    # Developers must input the command below on the mylocalbranch branch.
    "${CONDA_PREFIX}"/bin/git push origin mylocalbranch

**or** they can fetch and download the remote branch, named `myremotebranch`, that will hold the changes made on `mylocalbranch`

    # Developers must input the command below on the mylocalbranch branch.
    "${CONDA_PREFIX}"/bin/git switch -c myremotebranch origin/myremotebranch

This will switch the repository to the `myremotebranch`. Now, the developer needs to merge the changes made on `mylocalbranch`. If `mylocalbranch` is **NOT the main branch**, type

    # Developers must input the command below on the myremotebranch != main branch.
    "${CONDA_PREFIX}"/bin/git merge mylocalbranch

If the developer wants to merge their changes to the `mylocalbranch = main` branch instead, do a `squash` merge

    # Developers must input the command below on the main branch.
    "${CONDA_PREFIX}"/bin/git merge --squash mylocalbranch

If this does not create any merge conflicts, type

    # Developers must input the command below on the myremotebranch branch.
    "${CONDA_PREFIX}"/bin/git push origin myremotebranch

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
 
> [!Warning] Users must run `start_cocoa.sh` after cloning the project repository, so Cocoa can create appropriate soft-links.
>
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

