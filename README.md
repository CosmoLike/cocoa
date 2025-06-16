
<img width="780" alt="Screenshot 2024-09-16 at 10 55 14 PM" src="https://github.com/user-attachments/assets/f327e3bf-22c7-46e4-9bb5-14c8e7afc4c1">

# Table of contents
1. [Overview](#overview) 
2. [Installation of core packages](#required_packages_conda)
3. [Installation and Compilation of external modules](#cobaya_base_code)
4. [Running Examples](#cobaya_base_code_examples)
5. [Running ML emulators](#cobaya_base_code_examples_emul)
6. [Creating Cosmolike projects (external readme)](Cocoa/projects/)
7. [Appendix](#appendix)
    1. [FAQ: How do we debug Cocoa? Suggested steps](#running_wrong)
    2. [FAQ: How do we compile external modules?](#appendix_compile_separately)
    3. [FAQ: How do we run Cocoa with Docker?](#appendix_jupyter_whovian)
    4. [FAQ: How do we use an available Anaconda module on HPC?](#overview_anaconda)
    5. [FAQ: How do we install Conda?](#overview_miniconda)
    6. [FAQ: How do we set the environment for Machine Learning projects?](#ml_emulators)
    7. [FAQ: How do we push changes to the Cocoa main branch?](#push_main)
    8. [FAQ: How do we develop from a Git tag? A few more Git hacks](#dev_from_tag)
    9. [FAQ: How do we download additional likelihood data? (external readme)](Cocoa/external_modules/data)
   10. [FAQ: Where do we find common FAQs about external modules? (external readme)](Cocoa/external_modules/code)
   11. [FAQ: Where do we find common FAQs about Cosmolike? (external readme)](Cocoa/projects/)
   12. [FAQ: How can we improve our Bash/C/C++ knowledge?](#lectnotes)
   13. [Credits](#appendix_proper_credits)

# Overview <a name="overview"></a>

Cocoa allows users to run [CosmoLike](https://github.com/CosmoLike) routines inside the [Cobaya](https://github.com/CobayaSampler) framework. [CosmoLike](https://github.com/CosmoLike) can analyze data primarily from the [Dark Energy Survey](https://www.darkenergysurvey.org) and simulate future multi-probe analyses for LSST and Roman Space Telescope. 

Besides integrating [Cobaya](https://github.com/CobayaSampler) and [CosmoLike](https://github.com/CosmoLike), Cocoa introduces shell scripts that allow users to containerize [Cobaya](https://github.com/CobayaSampler), the Boltzmann codes, and multiple likelihoods. The container structure of Cocoa ensures that users can adopt consistent versions for the Fortran/C/C++ compilers and libraries across multiple machines. Such a systematic approach greatly simplifies the debugging process. 

Our scripts never install packages or Python modules in a global folder such as `$HOME/.local`. Here, `$HOME` denotes a shell environment variable that points to the user's home folder. Doing so would force cocoa packages to be global to the user, possibly breaking environments. Our scripts enable users to work on multiple Cocoa instances simultaneously, similar to what was possible with [CosmoMC](https://github.com/cmbant/CosmoMC). 

This readme file presents basic and advanced instructions for installing all [Cobaya](https://github.com/CobayaSampler) and [CosmoLike](https://github.com/CosmoLike) components.

# Installation of core packages <a name="required_packages_conda"></a>

Core packages include compilers and numerical libraries that users typically do not modify. We install most of these core packages via Conda, as shown below.

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
> Users working in an HPC environment that does not offer Anaconda or  may want to check the appendix [FAQ: How do we install Conda?](#overview_conda).

> [!TIP]
> We provide a Docker image named *whovian-cocoa* with Cocoa pre-installed and pre-compiled. For further instructions, refer to the appendix [FAQ: How do we run Cocoa with Docker?](#appendix_jupyter_whovian).

# Installation and Compilation of external modules <a name="cobaya_base_code"></a>

In this section, we assume users have previously activated the Cocoa conda environment.

**Step :one:**: Download Cocoa's latest release and go to its main folder (`cocoa/Cocoa`),

    git clone https://github.com/CosmoLike/cocoa.git --branch v4.0-beta27 cocoa

and

    cd ./cocoa/Cocoa

**Step :two:**: Run the script `setup_cocoa.sh` via
        
    source setup_cocoa.sh

This script downloads and decompresses external modules, requiring internet access. Therefore, users cannot run this script on an interactive compute node in an HPC environment where only the cluster login node can access the web.

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
> Cocoa does not install all the available external modules by default. If the user requires additional packages, refer to the appendix [FAQ: How do we compile external modules?](#appendix_compile_separately).

> [!NOTE]
> In case users need to run `setup_cocoa.sh` more than once, Cocoa will not download previously installed packages, cosmolike projects, or large datasets, unless the following keys are set on `set_installation_options.sh`
>
>     [Adapted from Cocoa/set_installation_options.sh shell script]
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

# Running Examples  <a name="cobaya_base_code_examples"></a>

We assume that you are still in the Conda cocoa environment from the previous `conda activate cocoa` command and that you are in the cocoa main folder `cocoa/Cocoa`, 

 **Step :one:**: Activate the private Python environment by sourcing the script `start_cocoa.sh`

    source start_cocoa.sh

Users will see a terminal like this: `$(cocoa)(.local)`. *This is a feature, not a bug*!

 **Step :two:**: Select the number of OpenMP cores (below, we set it to 8).
    
    export OMP_PROC_BIND=close; export OMP_NUM_THREADS=8

## Examples not involving Cosmolike

 **Step :three:**: The folder `projects/example` contains a few examples involving different likelihoods. So, run the `cobaya-run` on the first example following the commands below.

One model evaluation:

    mpirun -n 1 --oversubscribe --mca pml ^ucx --mca btl vader,tcp,self --bind-to core:overload-allowed --rank-by slot --map-by numa:pe=${OMP_NUM_THREADS} cobaya-run  ./projects/example/EXAMPLE_EVALUATE1.yaml -f
        
MCMC (we run MCMCs with 32 cores):

    mpirun -n 4 --oversubscribe --mca pml ^ucx --mca btl vader,tcp,self --bind-to core:overload-allowed --rank-by slot --map-by numa:pe=${OMP_NUM_THREADS} cobaya-run ./projects/example/EXAMPLE_MCMC1.yaml -f

## Examples involving Cosmolike

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

> [!Tip]
> Cocoa provides several cosmolike projects, not all of which are installed by default. To activate them, refer to the appendix [FAQ: How do we compile external modules?](#appendix_compile_separately).

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
> Why did we choose to work with two distinct shell environments, `(cocoa)` and `(.local)`? Our scripts enable users to work on multiple Cocoa instances, similar to what was possible with [CosmoMC](https://github.com/cmbant/CosmoMC). On each instance, our scripts install packages at
>
>      Cocoa/.local/bin
>      Cocoa/.local/include
>      Cocoa/.local/lib
>      Cocoa/.local/share
>
> Another reason 

# Running ML emulators <a name="cobaya_base_code_examples_emul"></a>

Cocoa contains a few transformer-based neural network emulators capable of simulating CMB, cosmolike, matter power spectrum, and distances. We provide a few scripts that exemplify their API. To run them, we assume users have commented out the following lines on `set_installation_options.sh` prior to running the `setup_cocoa.sh` and `compile_cocoa.sh` installation scripts.

      [Adapted from Cocoa/set_installation_options.sh shell script] 
      
      # inset # symbol in the lines below (i.e., unset these environmental keys)
      #export IGNORE_EMULTRF_CODE=1  #SaraivanovZhongZhu (SZZ) transformer-based emul
      #export IGNORE_EMULTRF_DATA=1  #SaraivanovZhongZhu (SZZ) transformer-based emul

      
      #export IGNORE_FGSPECTRA_CODE=1                        # to run EXAMPLE_EVALUATE26.yaml
      #export IGNORE_SIMONS_OBSERVATORY_LIKELIHOOD_CODE=1    # to run EXAMPLE_EVALUATE26.yaml
      #export IGNORE_SIMONS_OBSERVATORY_CMB_DATA=1           # to run EXAMPLE_EVALUATE26.yaml
 
Now, users must follow all the steps below.

 **Step :one:**: Activate the private Python environment by sourcing the script `start_cocoa.sh`

    source start_cocoa.sh

 **Step :two:**: Select the number of OpenMP cores.
    
    # No need to thread via OpenMP. However, 2-3 threads can reduce the TRF emulator runtime from ~0.2s to 0.1s.
    export OMP_PROC_BIND=close; export OMP_NUM_THREADS=1

 **Step :three:** Run `cobaya-run` on the first emulator example following the commands below.

One model evaluation:

    mpirun -n 1 --oversubscribe --mca pml ^ucx --mca btl vader,tcp,self --bind-to core:overload-allowed --rank-by slot --map-by numa:pe=${OMP_NUM_THREADS} cobaya-run ./projects/example/EXAMPLE_EVALUATE22.yaml -f
        
MCMC:

    mpirun -n 4 --oversubscribe --mca btl vader,tcp,self --bind-to core:overload-allowed --rank-by slot --map-by numa:pe=${OMP_NUM_THREADS} cobaya-run ./projects/example/EXAMPLE_MCMC22.yaml -f

PolyChord:

    mpirun -n 8 --oversubscribe --mca pml ^ucx --mca btl vader,tcp,self --bind-to core:overload-allowed --rank-by slot --map-by numa:pe=${OMP_NUM_THREADS} cobaya-run ./projects/example/EXAMPLE_POLY22.yaml -f

> [!NOTE]
> What should users do if they have not configured ML-related keys before running `setup_cocoa.sh` and `compile_cocoa.sh`, as rerunning these scripts can require a long time? Instead, run the following commands.
>
>      source start_cocoa.sh # even if (.local) is already active, users must run start_cocoa.sh again to update bash environment values
> 
> and
>
>      source ./installation_scripts/setup_pip_core_packages.sh     # install pip packages required by ML emulators
>      source ./installation_scripts/setup_emultrf.sh               # download emulator codes
>      source ./installation_scripts/unxv_emultrf.sh                # download pre-trained emulators
>
> and
>
>      source  ./installation_scripts/setup_fgspectra.sh               # download Simons Observatory foreground library
>      source  ./installation_scripts/compile_fgspectra.sh             # download Simons Observatory foreground library
>      source ./installation_scripts/setup_simons_observatory.sh       # download Simons Observatory Likelihood code (to run EXAMPLE_EVALUATE26.yaml)
>      source ./installation_scripts/unxv_simons_observatory.sh        # download Simons Observatory Likelihood data  (to run EXAMPLE_EVALUATE26.yaml)
>      source ./installation_scripts/compile_simons_observatory.sh     # compile Simons Observatory Likelihood code     (to run EXAMPLE_EVALUATE26.yaml)
> 
> Finally, rerun all the steps presented in this section, including step one. Users must reload the `(.local)` environment by rerunning `start_cocoa.sh` so Cocoa can create appropriate symlinks that expose the emulators to Cobaya.

# Appendix <a name="appendix"></a>

## :interrobang: FAQ: How do we debug Cocoa? Suggested steps <a name="running_wrong"></a>

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
      (...)

**Step :two:**: restart the Cocoa private environment by rerunning `source start_cocoa.sh` (every time users edit `set_installation_options.sh`, they must reload the `(.local)` environment by rerunning `start_cocoa.sh`).

**Step :three:**: compile the failed package by following the instructions in the appendix [FAQ: How do we compile external modules?](#appendix_compile_separately).

**Step 4️⃣**: rerun `setup_cocoa.sh` and `compile_cocoa.sh` to ensure all packages are installed and compiled correctly.

## :interrobang: FAQ: How do we compile external modules? <a name="appendix_compile_separately"></a>

To avoid excessive compilation or download times during development, users can use scripts located at `Cocoa/installation_scripts/` to download and compile only specific modules (or datasets). To take full advantage of them, users must first unset the appropriate keys on `set_installation_options.sh`, as exemplified below.

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
     source start_cocoa.sh # even if (.local) is already active, users must run start_cocoa.sh again to update bash environment values

 and
 
     source ./installation_scripts/setup_act_dr6.sh                # download likelihood code
     source ./installation_scripts/setup_simons_observatory.sh     # download likelihood code
     source ./installation_scripts/compile_act_dr6.sh              # compile likelihood code
     source ./installation_scripts/compile_simons_observatory.sh   # compile likelihood code
     source ./installation_scripts/unxv_act_dr6.sh                 # download and unpack likelihood data
     source ./installation_scripts/unxv_simons_observatory.sh      # download and unpack likelihood data

Finally, cocoa's `set_installation_options.sh` master script includes instructions to install several cosmolike projects. To activate them, manipulate the following lines on `set_installation_options.sh` 

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
     export OVERWRITE_EXISTING_COSMOLIKE_CODE=1 # dangerous (possible loss of uncommitted work)
                                                # If unset, users must manually delete cosmolike projects
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
      source start_cocoa.sh # even if (.local) is already active, users must run start_cocoa.sh again to update bash environment values
      
 and
 
      source ./installation_scripts/setup_cosmolike_projects.sh   # download all cosmolike projects  
      source ./installation_scripts/compile_all_projects.sh           # compile  all cosmolike project

In case users only want to compile a single cosmolike project (let's say the `roman_real` project)

      source ./projects/roman_real/scripts/compile_roman_real.sh
     
## :interrobang: FAQ: How do we run Cocoa with Docker? <a name="appendix_jupyter_whovian"></a>

We provide the Docker image [whovian-cocoa](https://hub.docker.com/r/vivianmiranda/whovian-cocoa) to facilitate the installation of Cocoa on Windows and macOS. This appendix assumes that users have already installed the Docker Engine on their local PC. For instructions on installing the Docker engine on specific operating systems, refer to [Docker's official documentation](https://docs.docker.com/engine/install/). 

 **Step :one:**: Create a folder and go to the location on the host computer where you want to provide access to the Docker container, as shown below. 

     mkdir -p cocoa_docker
     cd ./cocoa_docker

 **Step :two:**: Download the docker image *whovian-cocoa*, name the associated container `cocoa2023`, and run the container for the first time, type:

    docker run --platform linux/amd64 --hostname cocoa --name cocoa2025 -it -p 8888:8888 -v $(pwd):/home/whovian/host/ -v ~/.ssh:/home/whovian/.ssh:ro vivianmiranda/whovian-cocoa

This is a large image with a size of approximately 16GB, as it already contains the Cocoa version `v4.0beta26` installed and pre-compiled. 

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
> To run the Jupyter Notebook on the *whovian-cocoa* docker container installed on a remote server, adjust the command below
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

## :interrobang: FAQ: How do we use an available Anaconda module on HPC? <a name="overview_anaconda"></a>

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

:six: After completing steps :one:-:five:, check the internals of the `$HOME/.condarc` file with the command

    more "${HOME:?}"/.condarc

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
      - /project2/kicp/XXX/anaconda/pkgs  # adjust this directory path
    envs_dirs:
      - /project2/kicp/XXX/anaconda/envs/ # adjust this directory path

## :interrobang: FAQ: How do we install Conda? <a name="overview_miniconda"></a>

**Step :one:**: Download and run the Miniconda installation script. 

    export CONDA_DIR="/gpfs/home/XXX/miniconda" # replace this string!
    
    mkdir "${CONDA_DIR:?}"
    
    wget https://repo.anaconda.com/miniconda/Miniconda3-py310_25.3.1-1-Linux-x86_64.sh
    
    /bin/bash Miniconda3-py310_25.3.1-1-Linux-x86_64.sh -f -b -p "${CONDA_DIR:?}"

    /bin/bash 

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

## :interrobang: FAQ: How do users set the environment for Machine Learning projects? <a name="ml_emulators"></a>

Commenting out the environmental flags below *before running* `setup_cocoa.sh` will enable the installation of machine-learning-related libraries via pip.  

    [Adapted from Cocoa/set_installation_options.sh shell script] 
    # ------------------------------------------------------------------------------
    # If not set, pip_core_packages.sh will install several ML package
    # ------------------------------------------------------------------------------
    #export IGNORE_EMULATOR_GPU_PIP_PACKAGES=1
    (...)

In case users have already run `setup_cocoa.sh`, then run the command.

    source start_cocoa.sh # even if (.local) is already active, users must run start_cocoa.sh again to update bash environment values
   
and

    source ./installation_scripts/setup_pip_core_packages.sh
              
## :interrobang: FAQ: How do we push changes to the Cocoa main branch? <a name="push_main"></a>

Until recently, Cocoa development was a bit unstructured. Developers could push directly to the `main` branch, and small commits were not discouraged. Such flexible development rules will soon change when `v4.0` leaves the beta phase. We will protect the `main` branch by requiring every push to be reviewed by Cocoa's leading developers. Our new philosophy establishes that *a commit in the main branch should contain an atomic change that takes code from one working state to another working state with meaningful and well-tested improvements*. Therefore, developers should propose changes to the `main` branch in larger chunks (*via squash commits*), as shown below.

- :interrobang: **How to apply squash commits?**
  
**Step :one:**: create a development branch. Do not call the development branch `dev`, as `dev` is reserved for work done by the leading Cocoa developers. For concreteness, let's name the new branch `xyzdev`

   
    git switch -c xyzdev  # run on the main branch

> [!TIP]
> The use of developers' initials followed by `dev` helps make the branch easily identifiable.

> [!TIP]
> If the branch `xyzdev` already exists, use the `git switch` command without the `-c` flag. 

**Step :two:**: develop the proposed changes. We advise developers to commit frequently. In your branch, a commit does not need to be atomic, changing the code from one working state to another well-tested, meaningful working state. Developers can push to the server via the command.

    git push -u origin xyzdev   # run on the xyzdev branch

**Step :three:**: Once the developers created an atomic, meaningful, and well-tested improvement to Cocoa, the developer needs to merge any subsequent changes made in `main`.

    git merge main   # run on the xyzdev branch

This step may create conflicts that must be addressed before step four. 

**Step :four:**: Once the developers have merged recent changes made on the `main` branch, they must push to the main branch the modifications made on the `xyzdev` branch by first **squashing all your changes into a single commit**, as shown below

    git switch main  # run on the xyzdev branch

and

    git merge --squash xyzdev   # run on the main branch

and

    git commit -m "merge xyzdev branch"  # run on the main branch

and

    git push origin main  # run on the main branch

## :interrobang: FAQ: How do we develop from a git tag? A few more git hacks <a name="dev_from_tag"></a>

A useful git hack is related to developing Cocoa from a git tag. We reserve Git tags to set milestones in our development, so they serve as good starting points for coding localized new features (e.g., changes to a file that other developers have not recently modified) or bug fixes.

**Step :one: (optional)** If the developer has cloned the repository using the `https` URL address, then change the URL to the SSH-key-based address (if the developer has previously uploaded a public key to their GitHub account)

    git remote set-url origin git@github.com:CosmoLike/cocoa.git # that would allow users to push without typing a password

**Step :two:** Move the detached state to a new local branch via the command

    git switch -c xyzlocdev

Now, all commits will be associated with this local branch. 

> [!TIP]
> The use of developers' initials followed by `dev` helps make the branch easily identifiable.

**Step :three:** The developer has two options at the end of development. They can **either** create a new remote branch

    git push origin xyzlocdev # run on the xyzlocdev branch

**or** they can fetch and download the remote `xyzdev` branch, which will later absorb the changes made on `xyzlocdev`

    git switch -c xyzdev origin/xyzdev # run on the xyzlocdev branch

Finally, the developer needs to merge the changes made on `xyzlocdev`.

    git merge --squash xyzlocdev # run on the xyzdev branch

If this merge does not create any merge conflicts, type

    git push origin xyzdev # run on the xyzdev branch

## :interrobang: FAQ: How can we improve our Bash/C/C++ knowledge? <a name="lectnotes"></a>

A working knowledge of Python is required to understand the Cobaya framework at the developer level. Users must also be familiar with the Bash language to understand Cocoa's scripts. Proficiency in C and C++ is also needed to manipulate Cosmolike and the C++ Cobaya-Cosmolike C++ interface. Finally, users need to understand the Fortran-2003 language to modify CAMB.

Learning all these languages can be overwhelming, so to enable new users to do research that demands modifications on the inner workings of these codes, we include [here](cocoa_installation_libraries/LectNotes.pdf) a link to approximately 600 slides that provide an overview of Bash (slides ~1-137), C (slides ~138-371), and C++ (slides ~372-599). In the future, we aim to add lectures about Python and Fortran. 

## Credits <a name="appendix_proper_credits"></a>

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
