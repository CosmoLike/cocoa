<p align="center">
<img width="780" alt="Screenshot 2024-09-16 at 10 55 14 PM" src="https://github.com/user-attachments/assets/f327e3bf-22c7-46e4-9bb5-14c8e7afc4c1">
</p>

# Table of contents
1. [Overview](#overview) 
2. [Installation of core packages](#required_packages_conda)
3. [Installation and Compilation of external modules](#cobaya_base_code)
4. [Running Examples](#cobaya_base_code_examples)
5. [Running ML emulators](#cobaya_base_code_examples_emul)
6. [Running Hybrid Cosmolike-ML emulators](#cobaya_base_code_examples_emul2)
7. [Creating Cosmolike projects (external readme)](Cocoa/projects/)
8. [Credits](#appendix_proper_credits)
9. [Appendix](#appendix)
    1. [FAQ: How can users debug Cocoa? Suggested steps](#running_wrong)
    2. [FAQ: How can users compile external modules (not involving Cosmolike)?](#appendix_compile_separately)
    3. [FAQ: How can users install Cosmolike projects?](#appendix_compile_cosmolike_separately)
    4. [FAQ: How can users run Cocoa with Docker?](#appendix_jupyter_whovian)
    5. [FAQ: How can users run Cocoa on Google Colab?](#overview_google_colab)
    6. [FAQ: How can users install Conda?](#overview_miniforge)
    7. [FAQ: How can users set the appropriate environment for ML?](#ml_emulators)
    8. [FAQ: How can developers push changes to the Cocoa main branch?](#push_main)
    9. [FAQ: How can developers develop from a Git tag?](#dev_from_tag)
   10. [FAQ: How can users fork Cocoa correctly?](#fork_cocoa_tag)
   11. [FAQ: How can users download additional likelihood data? (external readme)](Cocoa/external_modules/data)
   12. [FAQ: Where do users find common FAQs about external modules? (external readme)](Cocoa/external_modules/code)
   13. [FAQ: Where do users find common FAQs about Cosmolike? (external readme)](Cocoa/projects/)
   14. [FAQ: Note on Planck-2018 low-ell SimAll EE and Gibbs TT likelihoods?](#planck2018lowell)
   15. [FAQ: How can users improve our Bash/C/C++ knowledge?](#lectnotes)

# Overview <a name="overview"></a>

Cocoa allows users to run [CosmoLike](https://github.com/CosmoLike) routines inside the [Cobaya](https://github.com/CobayaSampler) framework. [CosmoLike](https://github.com/CosmoLike) can analyze data from the [Dark Energy Survey](https://www.darkenergysurvey.org) and simulate future multi-probe analyses for LSST and Roman Space Telescope. 

Besides integrating [Cobaya](https://github.com/CobayaSampler) and [CosmoLike](https://github.com/CosmoLike), Cocoa also introduces shell scripts that allow users to containerize [Cobaya](https://github.com/CobayaSampler), the Boltzmann codes, and multiple likelihoods. The container structure of Cocoa ensures that users can adopt consistent versions for the Fortran/C/C++ compilers and libraries across multiple machines. Such a systematic approach greatly simplifies the debugging process. 

Our scripts never install packages or Python modules in a global folder such as `$HOME/.local`. Here, `$HOME` denotes a shell environment variable that points to the user's home folder. Doing so would force cocoa packages to be global to the user, possibly breaking environments. Our scripts enable users to work on multiple Cocoa instances simultaneously, similar to what was possible with [CosmoMC](https://github.com/cmbant/CosmoMC). 

This Readme file presents basic and advanced instructions for installing all [Cobaya](https://github.com/CobayaSampler) and [CosmoLike](https://github.com/CosmoLike) components **on linux** or **macOS-arm**.

We provide the Docker image [whovian-cocoa](https://hub.docker.com/r/vivianmiranda/whovian-cocoa) to facilitate the installation of Cocoa on Windows. 

# Installation of core packages <a name="required_packages_conda"></a>

Core packages include compilers and numerical libraries that users typically do not modify.

**Step :one:**: Download the appropriate `Python-3.10` compatible `yml` file

> [!Note]
> We recommend but do not yet impose the use of the package `conda-lock` and its tailored YML files to install the Cocoa conda environment, *see the [Tip] section below*

  - Linux
    
         wget https://raw.githubusercontent.com/CosmoLike/cocoa/refs/heads/dev/cocoapy310.yml

  - macOS (arm)
    
    Users on macOS may not have `wget` installed. If that is the case, run the following.
  
        conda activate; conda install wget

    This will activate the conda base environment (the prefix `(base)`) and install `wget`. Then, type.

        wget https://raw.githubusercontent.com/CosmoLike/cocoa/refs/heads/dev/cocoapy310-osxarm-base.yml
   
**Step :two:**: Create the Cocoa environment,

  - Linux
  
        conda env create --name cocoa --file=cocoapy310.yml

  - macOS (arm)

        conda env create --name cocoa --file=cocoapy310-osxarm-base.yml
    
and activate it

    conda activate cocoa

**Step :three:**: When and only when loading the conda cocoa environment for the first time, create the following symbolic links

  - Linux
    
        ln -s "${CONDA_PREFIX}"/bin/x86_64-conda_cos6-linux-gnu-gcc "${CONDA_PREFIX}"/bin/gcc
        ln -s "${CONDA_PREFIX}"/bin/x86_64-conda_cos6-linux-gnu-g++ "${CONDA_PREFIX}"/bin/g++
        ln -s "${CONDA_PREFIX}"/bin/x86_64-conda_cos6-linux-gnu-gfortran "${CONDA_PREFIX}"/bin/gfortran
        ln -s "${CONDA_PREFIX}"/bin/x86_64-conda-linux-gnu-gcc-ar "${CONDA_PREFIX}"/bin/gcc-ar
        ln -s "${CONDA_PREFIX}"/bin/x86_64-conda-linux-gnu-gcc-ranlib "${CONDA_PREFIX}"/bin/gcc-ranlib
        ln -s "${CONDA_PREFIX}"/bin/x86_64-conda-linux-gnu-ld "${CONDA_PREFIX}"/bin/ld

  - macOS (arm)

        ln -s "${CONDA_PREFIX}"/bin/clang "${CONDA_PREFIX}"/bin/gcc
        ln -s "${CONDA_PREFIX}"/bin/clang++ "${CONDA_PREFIX}"/bin/g++
  
**Step :four:**: When and only when loading the conda cocoa environment for the first time, install `git-lfs`

    git-lfs install

Users can now proceed to the **next section**.

> [!Warning]
> We advise users to stay away from all repositories managed by `Anaconda` due to draconian license restrictions that can affect researchers outside universities (e.g., NASA-JPL). See the Appendix [FAQ: How can we install Conda?](#overview_miniforge) for instructions on how to install `Miniforge`, which is a minimal installer of conda that downloads default packages from the `conda-forge` community-driven channel.

> [!Tip]
> We advise users to maintain *exact reproducibility* (across time) of the Cocoa conda environment by installing it via `conda-lock`,
> following the slightly more convoluted instructions below.
> 
> **Step :one:** Install the package `conda-lock` in a private conda environment to avoid conflicts.
> 
>     conda create -n lockenv -c conda-forge python=3.10 conda-lock=2.* wget
>
> and
> 
>     conda activate lockenv
>
> **Step :two:** Download the file appropriate conda-lock compatible `yml` file.
>   - Linux
>     
>         wget https://raw.githubusercontent.com/CosmoLike/cocoa/refs/heads/dev/cocoapy310-linux.yml
>
>   - macOS (arm)
>     
>         wget https://raw.githubusercontent.com/CosmoLike/cocoa/refs/heads/dev/cocoapy310-osxarm.yml
>
> **Step :three:** Create the conda environment
>   - Linux
>     
>         conda-lock install -n cocoa cocoapy310-linux.yml
>
>   - macOS (arm)
>
>         conda-lock install -n cocoa cocoapy310-osxarm.yml   
>
>  and activate it
>
>     conda activate cocoa
>
> Finally, proceed to **step :three:** in the general installation instructions. 

# Installation and Compilation of external modules <a name="cobaya_base_code"></a>

In this section, we assume users have previously activated the Cocoa conda environment.

**Step :one:**: Download Cocoa's latest release and go to its main folder (`cocoa/Cocoa`),

  - Stable Version
    
        git clone https://github.com/CosmoLike/cocoa.git --branch v4.05 cocoa

  - Testing beta release
    
        git clone https://github.com/CosmoLike/cocoa.git --branch v5.0beta1 cocoa

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
> Users who want to develop from a release version (e.g., `v4.0-beta20`) should read the appendix [FAQ: How can we push changes to the Cocoa main branch?](#push_main)

> [!TIP]
> Cocoa does not install all the available external modules by default. If the user needs additional packages, please refer to the appendix [FAQ: How can we compile external modules?](#appendix_compile_separately).

> [!NOTE]
> In case users need to rerun `setup_cocoa.sh`, Cocoa will not download previously installed packages, cosmolike projects, or large datasets, unless the following keys are set on `set_installation_options.sh`
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
>     #export OVERWRITE_EXISTING_COSMOLIKE_CODE=1 # dangerous (possible loss of uncommitted work)
>                                                 # if unset, users must manually delete cosmolike projects
>     #export REDOWNLOAD_EXISTING_ALL_DATA=1      # warning: some data is many GB

# Running Examples  <a name="cobaya_base_code_examples"></a>

We assume that you are still in the Conda cocoa environment from the previous `conda activate cocoa` command and that you are in the cocoa main folder `cocoa/Cocoa`, 

 **Step :one:**: Activate the private Python environment by sourcing the script `start_cocoa.sh`

    source start_cocoa.sh

Users will see a terminal like this: `$(cocoa)(.local)`. *This is a feature, not a bug*!

 **Step :two:**: Select the number of OpenMP cores (below, we set it to 8).

  - Linux
    
        export OMP_NUM_THREADS=8; export OMP_PROC_BIND=close; export OMP_PLACES=cores; export OMP_DYNAMIC=FALSE

  - macOS (arm)
    
        export OMP_NUM_THREADS=8; export OMP_PROC_BIND=disabled; export OMP_PLACES=cores; export OMP_DYNAMIC=FALSE
    
## Examples not involving Cosmolike

 **Step :three:**: The folder `projects/example` contains a few examples involving different likelihoods. So, run the `cobaya-run` on the first example following the commands below.

- **One model evaluation**:

  - Linux
  
        mpirun -n 1 --oversubscribe --mca pml ^ucx --mca btl vader,tcp,self \
           --bind-to core:overload-allowed --mca mpi_yield_when_idle 1 --report-bindings  \
           --rank-by slot --map-by numa:pe=${OMP_NUM_THREADS} \
           cobaya-run ./projects/example/EXAMPLE_EVALUATE1.yaml -f

  - macOS (arm)
  
        mpirun -n 1 --oversubscribe cobaya-run ./projects/example/EXAMPLE_EVALUATE1.yaml -f
    
- **MCMC (Metropolis-Hastings Algorithm)**:

    - Linux
  
          mpirun -n 4 --oversubscribe --mca pml ^ucx --mca btl vader,tcp,self \
             --bind-to core:overload-allowed --mca mpi_yield_when_idle 1 --report-bindings  \
             --rank-by slot --map-by numa:pe=${OMP_NUM_THREADS} \
             cobaya-run ./projects/example/EXAMPLE_MCMC1.yaml -f

  - macOS (arm)

        mpirun -n 4 --oversubscribe cobaya-run ./projects/example/EXAMPLE_MCMC1.yaml -f

## Examples involving Cosmolike

 **Step :three:**: The folder `projects/lsst_y1` contains a dozen examples involving different combinations of two-point correlation functions. So, run the `cobaya-run` on the first example following the commands below.

- **One model evaluation**:

  - Linux

        mpirun -n 1 --oversubscribe --mca pml ^ucx --mca btl vader,tcp,self \
           --bind-to core:overload-allowed --mca mpi_yield_when_idle 1 --report-bindings  \
           --rank-by slot --map-by numa:pe=${OMP_NUM_THREADS} \
           cobaya-run ./projects/lsst_y1/EXAMPLE_EVALUATE1.yaml -f

  - macOS (arm)

        mpirun -n 1 --oversubscribe  cobaya-run ./projects/lsst_y1/EXAMPLE_EVALUATE1.yaml -f

- **MCMC (Metropolis-Hastings Algorithm)**:

  - Linux
    
        mpirun -n 4 --oversubscribe --mca pml ^ucx --mca btl vader,tcp,self \
           --bind-to core:overload-allowed --mca mpi_yield_when_idle 1 --report-bindings  \
           --rank-by slot --map-by numa:pe=${OMP_NUM_THREADS} \
           cobaya-run ./projects/lsst_y1/EXAMPLE_MCMC1.yaml -f

  - macOS (arm)

        mpirun -n 4 --oversubscribe cobaya-run ./projects/lsst_y1/EXAMPLE_MCMC1.yaml -f
     
> [!Tip]
> Cocoa provides several Cosmolike projects, not all of which are installed by default. To activate them, please take a look at the appendix [FAQ: How can we compile external modules?](#appendix_compile_separately).

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
> The project `lsst-y1` contains Jupyter notebook examples located at `projects/lsst_y1`.

> [!NOTE]
> Why did we choose to work with two distinct shell environments, `(cocoa)` and `(.local)`? Our scripts enable users to work on multiple Cocoa instances, similar to what was possible with [CosmoMC](https://github.com/cmbant/CosmoMC). In each instance, our scripts install packages at
>
>      Cocoa/.local/bin
>      Cocoa/.local/include
>      Cocoa/.local/lib
>      Cocoa/.local/share
>
> Consistency of the environment across all Cocoa instances is crucial, and the `start_cocoa.sh`/`stop_cocoa.sh` scripts handle the loading and unloading of environmental path variables. 

# Running ML emulators <a name="cobaya_base_code_examples_emul"></a>

Cocoa contains a few transformer- and CNN-based neural network emulators capable of simulating the CMB, cosmolike outputs, matter power spectrum, and distances. We provide a few scripts that exemplify their API. To run them, users ensure the following lines are commented out in `set_installation_options.sh` before running the `setup_cocoa.sh` and `compile_cocoa.sh`. By default, these lines should be commented out, but it is worth checking.

      [Adapted from Cocoa/set_installation_options.sh shell script] 
      # insert the # symbol (i.e., unset these environmental keys  on `set_installation_options.sh`)
      #export IGNORE_EMULTRF_CODE=1              #SaraivanovZhongZhu (SZZ) transformer/CNN-based emulators
      #export IGNORE_EMULTRF_DATA=1            
      #export IGNORE_LIPOP_LIKELIHOOD_CODE=1     # to run EXAMPLE_EMUL_(EVALUATE/MCMC/NAUTILUS/EMCEE1).yaml
      #export IGNORE_LIPOP_CMB_DATA=1           
      #export IGNORE_ACTDR6_CODE=1               # to run EXAMPLE_EMUL_(EVALUATE/MCMC/NAUTILUS/EMCEE1).yaml
      #export IGNORE_ACTDR6_DATA=1         
      #export IGNORE_NAUTILUS_SAMPLER_CODE=1     # to run PROJECTS/EXAMPLE/EXAMPLE_EMUL_NAUTILUS1.py
      #export IGNORE_POLYCHORD_SAMPLER_CODE=1    # to run PROJECTS/EXAMPLE/EXAMPLE_EMUL_POLY1.yaml
      #export IGNORE_GETDIST_CODE=1              # to run EXAMPLE_TENSION_METRICS.ipynb
      #export IGNORE_TENSIOMETER_CODE=1          # to run EXAMPLE_TENSION_METRICS.ipynb
          

Now, users must follow all the steps below.

> [!Note]
> We provide SLURM job script examples in the `projects/example/script` folder, which allow users to run the examples below in an HPC environment.

> [!Note]
> What if users have not configured ML-related keys before sourcing `setup_cocoa.sh`?
> 
> Answer: Comment the keys below before rerunning `setup_cocoa.sh`.
> 
>     [Adapted from Cocoa/set_installation_options.sh shell script]
>     # These keys are only relevant if you run setup_cocoa multiple times
>     #export OVERWRITE_EXISTING_ALL_PACKAGES=1    
>     #export OVERWRITE_EXISTING_COSMOLIKE_CODE=1 
>     #export REDOWNLOAD_EXISTING_ALL_DATA=1

 **Step :one:**: Activate the private Python environment by sourcing the script `start_cocoa.sh`

    source start_cocoa.sh

 **Step :two:**: Ensure OpenMP is **OFF**.
    
    export OMP_NUM_THREADS=1
    
 **Step :three:** Run `cobaya-run` on the first emulator example following the commands below.

- **One model evaluation**:

  - Linux
    
        mpirun -n 1 --oversubscribe --mca pml ^ucx --mca btl vader,tcp,self --rank-by slot \
          --bind-to core:overload-allowed --map-by slot --mca mpi_yield_when_idle 1 \
          cobaya-run ./projects/example/EXAMPLE_EMUL_EVALUATE1.yaml -f

  - macOS (arm)

        mpirun -n 1 --oversubscribe cobaya-run ./projects/example/EXAMPLE_EMUL_EVALUATE1.yaml -f
    
- **MCMC (Metropolis-Hastings Algorithm)**:

  - Linux
    
        mpirun -n 4 --oversubscribe --mca pml ^ucx --mca btl vader,tcp,self --rank-by slot \
            --bind-to core:overload-allowed --map-by slot --mca mpi_yield_when_idle 1 \
            cobaya-run ./projects/example/EXAMPLE_EMUL_MCMC1.yaml -r

  - macOS (arm)

        mpirun -n 4 --oversubscribe cobaya-run ./projects/example/EXAMPLE_EMUL_MCMC1.yaml -r

> [!Note]
> The examples below may require a large number of MPI workers. Before running them, it may be necessary to increase 
> the limit of threads that can be created (at *UofA/SBU HPC* type `ulimit -u 1000000`), otherwise users 
> may encounter the error `libgomp: Thread creation failed`

- **PolyChord**:

  - Linux
    
        mpirun -n 90 --oversubscribe --mca pml ^ucx --mca btl vader,tcp,self --rank-by slot \
            --bind-to core:overload-allowed --map-by slot --mca mpi_yield_when_idle 1 \
            cobaya-run ./projects/example/EXAMPLE_EMUL_POLY1.yaml -r

  - macOS (arm)
 
        mpirun -n 12 --oversubscribe cobaya-run ./projects/example/EXAMPLE_EMUL_POLY1.yaml -r
    
> [!Note]
> When running `PolyChord` or any of our scripts in more than one node, replace `--mca btl vader,tcp,self` by `--mca btl tcp,self`.

The `Nautilis`, `Minimizer`, `Profile`, and `Emcee` scripts below contain an internally defined `yaml_string` that specifies priors, 
likelihoods, and the theory code, all following Cobaya Conventions.

- **Nautilus**:

  - Linux
    
        mpirun -n 90 --oversubscribe --mca pml ^ucx --mca btl vader,tcp,self --rank-by slot \
            --bind-to core:overload-allowed --map-by slot --mca mpi_yield_when_idle 1 \
            python -m mpi4py.futures ./projects/example/EXAMPLE_EMUL_NAUTILUS1.py \
                --root ./projects/example/ --outroot "EXAMPLE_EMUL_NAUTILUS1"  \
                --maxfeval 450000 --nlive 2048 --neff 15000 --flive 0.01 --nnetworks 5

  - macOS (arm)

        mpirun -n 12 --oversubscribe python -m mpi4py.futures ./projects/example/EXAMPLE_EMUL_NAUTILUS1.py \
                --root ./projects/example/ --outroot "EXAMPLE_EMUL_NAUTILUS1"  \
                --maxfeval 450000 --nlive 2048 --neff 15000 --flive 0.01 --nnetworks 5

> [!NOTE]
> What if the user runs an `Nautilus` chain with `maxeval` insufficient for producing `neff` samples? `Nautilus` creates a checkpoint at `chains/outroot_checkpoint.hdf5`.

- **Emcee**:

    - Linux
      
          mpirun -n 21 --oversubscribe --mca pml ^ucx --mca btl vader,tcp,self --rank-by slot \
              --bind-to core:overload-allowed --map-by slot --mca mpi_yield_when_idle 1 \
              python ./projects/example/EXAMPLE_EMUL_EMCEE1.py --root ./projects/example/ \
                  --outroot "EXAMPLE_EMUL_EMCEE1" --maxfeval 80000

    - macOS (arm)

          mpirun -n 12 --oversubscribe python ./projects/example/EXAMPLE_EMUL_EMCEE1.py \
                --root ./projects/example/ --outroot "EXAMPLE_EMUL_EMCEE1" --maxfeval 80000

  The number of steps per MPI worker is $n_{\mathrm{sw}} =  \mathrm{maxfeval}/n_{\mathrm{w}}$,
  with the number of walkers being $n_{\mathrm{w}}=\mathrm{max}(3n_{\mathrm{params}},n_{\mathrm{MPI}})$.
  For proper convergence, each walker should traverse 50 times its autocorrelation length ($\tau$),
  which is provided in the header of the output chain file. A reasonable rule of thumb is to assume
  $\tau > 200$ and therefore set $\mathrm{maxfeval} > 10,000 \times n_{\mathrm{w}}$.

  With these numbers, users may ask when `Emcee` is preferable to `Metropolis-Hastings`?
  Here are a few numbers based on a `Planck CMB (l < 396) + SN + BAO + LSST-Y1` chain with 38 parameters in total.
  1) `MH` achieves convergence with $n_{\mathrm{sw}} \sim 150,000$, but only requires four walkers.
  2) `Emcee` has $\tau \sim 300$, so it requires $n_{\mathrm{sw}} \sim 15,000$ when running with $n_{\mathrm{w}}=114$.
  
  Conclusion: `Emcee` requires $\sim 3$ more evaluations in this case, but the number of evaluations per MPI worker (assuming one MPI worker per walker) is reduced by $\sim 10$.
  Therefore, `Emcee` seems well-suited for cases where the evaluation of a single cosmology is time-consuming (and there is no slow/fast decomposition).

> [!NOTE]
> What if the user runs an `Emcee` chain with `maxeval` insufficient for convergence? `Emcee` creates a checkpoint at `chains/outroot.h5`.

- **Sampler Comparison**

    The script that generated the plot below is provided at `projects/example/scripts/EXAMPLE_PLOT_COMPARE_CHAINS.py`. The Google Colab notebook [Example Sampler Comparison](https://github.com/CosmoLike/CoCoAGoogleColabExamples/blob/main/Cocoa_Example_(Sampler_Comparison).ipynb) can also reconstruct a similar version of this figure.

    <p align="center">
    <img width="750" height="750" alt="projects_example_sampler_comparison" src="https://github.com/user-attachments/assets/d3639673-36ea-4fd9-9c91-1f5b97845fe0" />
    </p>
  
- **Global Minimizer**:

  Our minimizer is a reimplementation of `Procoli`, developed by Karwal et al ([arXiv:2401.14225](https://arxiv.org/abs/2401.14225)) 

    - Linux
      
          mpirun -n 21 --oversubscribe --mca pml ^ucx --mca btl vader,tcp,self --rank-by slot \
              --bind-to core:overload-allowed --map-by slot --mca mpi_yield_when_idle 1 \
              python ./projects/example/EXAMPLE_EMUL_MINIMIZE1.py --root ./projects/example/ \
                  --outroot "EXAMPLE_EMUL_MIN1" --nstw 200
  
    - macOS (arm)

          mpirun -n 12 --oversubscribe python ./projects/example/EXAMPLE_EMUL_MINIMIZE1.py \
                --root ./projects/example/ --outroot "EXAMPLE_EMUL_MIN1" --nstw 200
       
  The number of steps per Emcee walker per temperature is $n_{\mathrm{stw}}$,
  and the number of walkers is $n_{\mathrm{w}}=\mathrm{max}(3n_{\mathrm{params}},n_{\mathrm{MPI}})$. The minimum number of total evaluations is then 
  $3n_{\mathrm{params}} \times n_{\mathrm{T}} \times n_{\mathrm{stw}}$, which can be distributed among $n_{\mathrm{MPI}} = 3n_{\mathrm{params}}$ MPI processes for faster results.
  Do maintain $n_{\mathrm{stw}} > 200$ for reliable convergence in LCDM (see plot below).
  The same rule applies to *Profile* and *Scan* codes, as they are all based on the same minimization strategy.

  The script that generated the plot below is provided at `projects/example/scripts/EXAMPLE_MIN_COMPARE_CONV.py`. The Google Colab notebook [Test Minimizer Convergence](https://github.com/CosmoLike/CoCoAGoogleColabExamples/blob/main/Cocoa_Example_Test_Minimizer_Convergence.ipynb) can also reconstruct a similar version of this figure. 

  <p align="center">
  <img width="700" height="470" alt="Screenshot 2025-08-04 at 7 05 53 AM" src="https://github.com/user-attachments/assets/a48b267a-beba-4e53-9dbf-e3c5a24daff1" />
  </p>

  Below we show a case with $n_{\mathrm{param}} = 38$ that illustrates the need for performing convergence tests on a case-by-case basis.
  In this example, the total number of evaluations for a reliable minimum is approximately $319,200$ ($n_{\mathrm{stw}} \sim 700$), distributed among $n_{\mathrm{MPI}} = 114$ processes for faster results.

  <p align="center">
  <img width="700" height="470" alt="Screenshot 2025-08-13 at 5 29 59 PM" src="https://github.com/user-attachments/assets/c43b8eea-ee2e-443d-a497-cb9b2dae2fc3" />
  </p>

- **Profile**: 

    - Linux
      
          mpirun -n 21 --oversubscribe --mca pml ^ucx --mca btl vader,tcp,self --rank-by slot \
              --bind-to core:overload-allowed --map-by slot --mca mpi_yield_when_idle 1 \
              python ./projects/example/EXAMPLE_EMUL_PROFILE1.py \
                  --root ./projects/example/ --cov 'chains/EXAMPLE_EMUL_MCMC1.covmat' \
                  --outroot "EXAMPLE_EMUL_PROFILE1" --factor 3 --nstw 200 --numpts 10 \
                  --profile 1 --minfile "./projects/example/chains/EXAMPLE_EMUL_MIN1.txt"

    - macOS (arm)

          mpirun -n 12 --oversubscribe python ./projects/example/EXAMPLE_EMUL_PROFILE1.py \
                  --root ./projects/example/ --cov 'chains/EXAMPLE_EMUL_MCMC1.covmat' \
                  --outroot "EXAMPLE_EMUL_PROFILE1" --factor 3 --nstw 200 --numpts 10 \
                  --profile 1 --minfile "./projects/example/chains/EXAMPLE_EMUL_MIN1.txt"
      
    Profile provides the optional argument `minfile`, as it is significantly faster to run the profile script with a previously provided global minimum. 
    The profile also provides the optional argument `cov`. Again, it is considerably more efficient to employ a covariance matrix from a converged chain. 

    The argument `factor` specifies the start and end of the parameter being profiled:

          start value ~ mininum value - factor*np.sqrt(np.diag(cov))
          end   value ~ mininum value + factor*np.sqrt(np.diag(cov))

    We advise $\mathrm{factor} \sim 3$ for parameters that are well constrained by the data when a covariance matrix is provided.
    If `cov` is not supplied, the code estimates one internally from the prior.
    If a parameter is poorly constrained or `cov` is not given, we recommend $\mathrm{factor} \ll 1$.

    The script that generated the plot below is provided at `projects/example/scripts/EXAMPLE_PLOT_PROFILE1.py`. The Google Colab notebook [Example Profile Likelihood](https://github.com/CosmoLike/CoCoAGoogleColabExamples/blob/main/Cocoa_Example_Profile_Likelihoods.ipynb) can also reconstruct a similar version of this figure. 
  
    <p align="center">
    <img width="832.32" height="617.76" alt="Screenshot 2025-08-02 at 8 42 41 PM" src="https://github.com/user-attachments/assets/22182688-2865-4b15-a80b-783ddd21f715" />
    </p>

- **Profile method 2**:

    If the dimensionality of the problem is not large, and the spacing between values of the parameter
    being profiled is small, it can be considerably faster to use a simple scipy `Nelder-Mead`
    to calculate the profile. Here, the `minfile` and `cov` options are mandatory.

    - Linux
      
          mpirun -n 1 --mca pml ^ucx --mca btl vader,tcp,self --rank-by slot \
              --bind-to core:overload-allowed --map-by slot --mca mpi_yield_when_idle 1 \
              python ./projects/example/EXAMPLE_EMUL_PROFILE_SCIPY1.py \
                  --root ./projects/example/ --cov 'chains/EXAMPLE_EMUL_MCMC1.covmat' \
                  --outroot "EXAMPLE_EMUL_PROFILE1M2" --factor 3 --maxfeval 5000 --numpts 10 \
                  --profile 1 --minfile "./projects/example/chains/EXAMPLE_EMUL_MIN1.txt"

    - macOS (arm)

          mpirun -n 1 python ./projects/example/EXAMPLE_EMUL_PROFILE_SCIPY1.py \
                --root ./projects/example/ --cov 'chains/EXAMPLE_EMUL_MCMC1.covmat' \
                --outroot "EXAMPLE_EMUL_PROFILE1M2" --factor 3 --maxfeval 5000 --numpts 10 \
                --profile 1 --minfile "./projects/example/chains/EXAMPLE_EMUL_MIN1.txt"
      
    The script that generated the plot below is provided at `projects/example/scripts/EXAMPLE_PLOT_PROFILE1_COMP.py`
  
    <p align="center">
    <img width="832.32" height="617.76" alt="example_profile_comp" src="https://github.com/user-attachments/assets/ba0c0629-bd3b-4274-9f24-9db5929dc35c" />
    </p>

- **Scan**: 

  This profile code has a different MPI strategy. It scans one parameter on the entire prior,
  with each MPI being assigned to one minimization (not Emcee walker!). This is a strategy when probing 
  beyond-LCDM parameters with oscilatory behavior (e.g., Monodromic Dark Energy).

    - Linux
      
          mpirun -n 90 --oversubscribe --mca pml ^ucx --mca btl vader,tcp,self \
              --bind-to core:overload-allowed --map-by slot --mca mpi_yield_when_idle 1 \
              python -m mpi4py.futures ./projects/example/EXAMPLE_EMUL_SCAN1.py \
                  --root ./projects/example/ --outroot "EXAMPLE_EMUL_SCAN1" --nstw 200 --profile 1

    - macOS (arm)

          mpirun -n 12 --oversubscribe python -m mpi4py.futures ./projects/example/EXAMPLE_EMUL_SCAN1.py \
              --root ./projects/example/ --outroot "EXAMPLE_EMUL_SCAN1" --nstw 200 --profile 1
      
- **Tension Metrics**

    We provide the Jupyter Notebook and the SLURM script, located at

      projects/example/EXAMPLE_EMUL_MCMC_TENSION_METRICS/EXAMPLE_TENSION_METRICS.ipynb
      projects/example/scripts/EXAMPLE_EMUL_TM.sbatch

    that exemplifies how emulators can be used to quickly assess the tension between data sets.
    All Tension metrics were computed using the package `Tensiometer`.

    The plot below, taken from `EXAMPLE_TENSION_METRICS.ipynb`, exemplifies the tension between CMB+BAO vs. Type Ia SN.

<p align="center">
<img width="790" height="333" alt="Unknown" src="https://github.com/user-attachments/assets/b5ef808e-0efa-4bb9-852a-854f6531e3ed" />
</p>

> [!NOTE]
> What about [CosmoPower](https://alessiospuriomancini.github.io/cosmopower/)? CosmoPower is a suite of popular emulators developed by Prof. Alessio Mancini
> and collaborators. Our Cocoa port relies on the work done by Hidde Jense @HTJense (kudos to his/their work) as well as the Simons Observatory SOLikeT code.
> Even though they are not the suite of AI-powered emulators adopted by Cocoa developers, we do provide limited support for running them.
> To set up and compile Cosmopower and also the corresponding Cobaya Wrapper, comment the following keys before running `setup_cocoa.sh` and `compile_cocoa.sh`.
> 
>     [Adapted from Cocoa/set_installation_options.sh shell script]
>     #export IGNORE_COSMOPOWER_DATA=1  
>     #export IGNORE_COSMOPOWER_CODE=1
>     (...)
>     # Edit below if you have your own set of CosmoPower Emulators
>     export COSMOPOWER_URL_DATA="https://github.com/cosmopower-organization/jense_2024_emulators.git"
>     export COSMOPOWER_URL_DATA_COMMIT="4317635eed70289ee1ec6b3df828027173071e36"
> 
> We also provide a few examples using CosmoPower (they all end with `CP` suffix)
>
> - **One model evaluation**:
>
>  - Linux
>    
>        mpirun -n 1 --oversubscribe --mca pml ^ucx --mca btl vader,tcp,self --rank-by slot \
>          --bind-to core:overload-allowed --map-by slot --mca mpi_yield_when_idle 1 \
>          cobaya-run ./projects/example/EXAMPLE_EMUL_EVALUATE2_CP.yaml -f
>
>  - macOS (arm)
>
>        mpirun -n 1 --oversubscribe cobaya-run ./projects/example/EXAMPLE_EMUL_EVALUATE2_CP.yaml -f
> 

# Running Hybrid Cosmolike-ML emulators <a name="cobaya_base_code_examples_emul2"></a>

> [!Warning]
> The code and examples associated with this section are still in alpha stage

Our main line of research involves emulators that simulate the entire Cosmolike data vectors, and each project (LSST, Roman, DES) contains its own README with emulator examples. The speed of such emulators is incredible, especially when GPUs are available, and our emulators do take advantage of the CPU-GPU integration on Apple MX chips. For example, the average timing of lsst-y1 cosmic shear data vector emulation is around 0.005s ($\sim$ 200828 evaluations in $\sim$ 850.5 seconds) on a macOS M2 Pro.

While the data vector emulators are incredibly fast, there is an intermediate approach that emulates only the Boltzmann outputs (comoving distance, linear and nonlinear matter power spectrum). This hybrid-ML case can offer greater flexibility, especially in the initial phases of a research project, as changes to the modeling of nuisance parameters or to the assumed galaxy distributions do not require retraining of the network. 

Examples in the hybrid case all have the prefix **EXAMPLE_EMUL2** (note the `2`). The required flags on `set_installation_options.sh` are similar to what we showed in the previous emulator section.

Now, users must follow all the steps below.

 **Step :one:**: Activate the private Python environment by sourcing the script `start_cocoa.sh`

    source start_cocoa.sh

 **Step :two:**: Select the number of OpenMP cores. Below, we set it to 4, the ideal setting for hybrid examples.

  - Linux
    
        export OMP_NUM_THREADS=4; export OMP_PROC_BIND=close; export OMP_PLACES=cores; export OMP_DYNAMIC=FALSE

  - macOS (arm)
    
        export OMP_NUM_THREADS=4; export OMP_PROC_BIND=disabled; export OMP_PLACES=cores; export OMP_DYNAMIC=FALSE
    
 **Step :three:** Run `cobaya-run` on the first emulator example, following the commands below (here we only provide lsst-y1 examples).

- **One model evaluation**:

  - Linux
    
        mpirun -n 1 --oversubscribe --mca pml ^ucx --mca btl vader,tcp,self \
           --bind-to core:overload-allowed --mca mpi_yield_when_idle 1 --report-bindings  \
           --rank-by slot --map-by numa:pe=${OMP_NUM_THREADS} \
           cobaya-run ./projects/lsst_y1/EXAMPLE_EMUL2_EVALUATE1.yaml -f

  - macOS (arm)
    
        mpirun -n 1 --oversubscribe  cobaya-run ./projects/lsst_y1/EXAMPLE_EMUL2_EVALUATE1.yaml -f
    
- **MCMC (Metropolis-Hastings Algorithm)**:

  - Linux
    
        mpirun -n 4 --oversubscribe --mca pml ^ucx --mca btl vader,tcp,self --rank-by slot \
            --bind-to core:overload-allowed --map-by slot --mca mpi_yield_when_idle 1 \
            cobaya-run ./projects/lsst_y1/EXAMPLE_EMUL_EMUL2_MCMC1.yaml -r

  - macOS (arm)

        mpirun -n 4 --oversubscribe cobaya-run ./projects/lsst_y1/EXAMPLE_EMUL2_MCMC1.yaml -r
    
Details on the matter power spectrum emulator designs will be presented in the [emulator_code](https://github.com/CosmoLike/emulators_code) repository. Basically, we apply standard neural network techniques to generalize the *syren-new* Eq. 6 of [arXiv:2410.14623](https://arxiv.org/abs/2410.14623) formula for the linear power spectrum (w0waCDM with a fixed neutrino mass of $0.06$ eV) to new models, extended ranges, or higher precision. Similarly, we use networks to generalize the *syren-Halofit* LCDM nonlinear boost fit (Eq. 11 of [arXiv:2402.17492](https://arxiv.org/abs/2402.17492)).

> [!NOTE] 
> Users can decide not to correct the *syren-new* formula for the linear power spectrum (flag in the yaml). Although we have not conducted extensive studies of the caveats of the syren-new approximation, it appears sufficient for w0waCDM forecasts when combined with the Euclid Emulator to compute the nonlinear boost.
>
> For back-of-the-envelope LCDM calculations (e.g., to test cosmolike features), users can also choose not to correct the *syren-Halofit* formula for the LCDM nonlinear boost (see figure below). In this case, the overhead on top of cosmolike computations is minimum, at the order of $0.01$ seconds on a macOS M2Pro laptop. 
>
> <p align="center">
>  <span style="display:flex; justify-content:center; gap:16px; flex-wrap:wrap;">
>    <img width="450" alt="compare_emul_hemul" src="https://github.com/user-attachments/assets/7ada9b3b-db01-499f-8170-a2db5ef90636" />
>    <img width="450" alt="Screenshot 2026-01-06 at 2 13 06 AM" src="https://github.com/user-attachments/assets/d19f0900-cebc-41ca-9028-1d4e0bd40cc9" />
>  </span>
> </p>




## Credits <a name="appendix_proper_credits"></a>

- **Acknowledgments**:

  The entire CoCoA team is deeply grateful to everyone who has contributed to the development of our code.

  **Special Thanks to:**

  - Profs. Tim Eifler and Elisabeth Krause for their support of this idea since its inception in 2018 and all Cosmolike-related development.
  - Profs. Antony Lewis and Jesús Torrado for helping me understand Cobaya since its early days in 2019.
  - Jonathan Gordon, Joshua Kable, João Rebouças, Evan Saraivanov, Diogo Souza, Jiachuan Xu, Yijie Zhu, and KunHao Zhong for working on CoCoA on many fruitful projects at Stony Brook Univ. and the Univ. of Arizona.
  - Evan Saraivanov, Yijie Zhu, and KunHao Zhong for developing the emulator interface within the CoCoA framework (and training emulators).
  - Victoria Loyd for developing the hybrid-emulator interface within the CoCoA framework (and training emulators).
  - Haley Bowden, Kali Cao, Nihar Dalal, Yu-Hsiu Huang, Niko, Junzhou Zhang, and members of the Roman HLIS Cosmology PIT for all Roman-specific development and testing.

The following is not an exhaustive list of the codes we use/download/adopt

- [Cobaya](https://github.com/CobayaSampler) is a framework developed by Dr. Jesus Torrado and Prof. Anthony Lewis
- [Cosmolike](https://github.com/CosmoLike) is a framework developed by Prof. Elisabeth Krause and Prof. Tim Eifler
- [CAMB](https://github.com/cmbant/CAMB) is a Boltzmann code developed by Prof. Anthony Lewis
- [CLASS](https://github.com/lesgourg/class_public) is a Boltzmann code developed by Prof. Julien Lesgourgues and Dr. Thomas Tram
- [Polychord](https://github.com/PolyChord/PolyChordLite) is a sampler code developed by Dr. Will Handley, Prof. Lasenby, and Prof. M. Hobson
- [CLIK](https://github.com/benabed/clik) is the likelihood code used to analyze Planck and SPT data, maintained by Prof. Karim Benabed
- [SPT](https://github.com/SouthPoleTelescope/spt3g_y1_dist) is the official likelihood of the South Pole Telescope 3G Year 1
- [MFLike](https://github.com/simonsobs/LAT_MFLike) is the official likelihood of the Simons Observatory
- [ACTLensing](https://github.com/ACTCollaboration/act_dr6_lenslike) is the official lensing likelihood of the ACT collaboration, developed by Prof. Mathew Madhavacheril
- [HiLLiPoP CMB likelihood](https://github.com/planck-npipe/hillipop.git) is a multifrequency CMB likelihood for Planck data.
- [Lollipop CMB likelihood](https://github.com/planck-npipe/lollipop.git) is a Planck low-l polarization likelihood.
  
Following best practices, Cocoa scripts download most external modules from their original repositories. Although our repository includes a few likelihoods in compressed xz file format, we do not want to discourage users from cloning code and data from their original repositories.  The work of those authors is extraordinary, and users **must cite them** appropriately.

# Appendix <a name="appendix"></a>

## :interrobang: FAQ: How can users debug Cocoa? Suggested steps <a name="running_wrong"></a>

**Step :one:**: define the `COCOA_OUTPUT_VERBOSE` and `COSMOLIKE_DEBUG_MODE` flags on `set_installation_options.sh` to obtain a more detailed output, as shown below (for even more output, define `COCOA_OUTPUT_DEBUG` as well)
  
      [Adapted from Cocoa/set_installation_options.sh shell script] 
      # ------------------------------------------------------------------------------
      # VERBOSE AS DEBUG TOOL --------------------------------------------------------
      # ------------------------------------------------------------------------------
      export COCOA_OUTPUT_VERBOSE=1
      #export COCOA_OUTPUT_DEBUG=1 # turn on bash strict mode (set -exo pipefail) on  
                                   # instalation_scripts/setup/compile_x.sh scripts 
      # ------------------------------------------------------------------------------
      # If set, COSMOLIKE will compile with DEBUG flags ------------------------------
      # ------------------------------------------------------------------------------
      export COSMOLIKE_DEBUG_MODE=1
      (...)

**Step :two:**: restart the Cocoa private environment by rerunning `source start_cocoa.sh` (every time users edit `set_installation_options.sh`, they must reload the `(.local)` environment by rerunning `start_cocoa.sh`).

**Step :three:**: compile the failed package by following the instructions in the appendix [FAQ: How can we compile external modules?](#appendix_compile_separately).

**Step 4️⃣**: rerun `setup_cocoa.sh` and `compile_cocoa.sh` to ensure all packages are installed and compiled correctly.

## :interrobang: FAQ: How can users compile external modules (not involving Cosmolike)? <a name="appendix_compile_separately"></a>

To avoid excessive compilation or download times during development, users may run scripts located at `Cocoa/installation_scripts/` directly to download and compile only specific modules (or datasets). To take full advantage of them, users must first unset the appropriate keys on `set_installation_options.sh`, as exemplified below.

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
9(
Whenever users edit `set_installation_options.sh`, they must reload the Cocoa private environment (.local) by sourcing `start_cocoa.sh` (even if (.local) is already active). To do that, follow the commands below.

     cd ./cocoa/Cocoa

and

     source start_cocoa.sh # even if (.local) is already active, users must run start_cocoa.sh again to update bash environment values

Similar to `setup_cocoa.sh` and `compile_cocoa.sh`, each code package contains a `setup` scripts that download it from the internet and a `compile` scripts that compile it. E.g., the `ACT-DR6` likelihood code can be downloaded and installed via the bash commands
 
     source ./installation_scripts/setup_act_dr6.sh                # download likelihood code

and

     source ./installation_scripts/compile_act_dr6.sh              # compile likelihood code

In addition to `setup` and `compile` scripts, Cocoa contains `unxv` scripts that download data sets. For instance, to download the ACT-DR6 data, users must source the script

     source ./installation_scripts/unxv_act_dr6.sh                 # download and unpack likelihood data

## :interrobang: FAQ: How can users install Cosmolike projects?  <a name="appendix_compile_cosmolike_separately"></a>

The script `set_installation_options.sh` includes instructions for installing several Cosmo-like-based projects. To activate them, modify the following lines in `set_installation_options.sh` by inserting the symbol `#` before the name of the module users want to be installed and compiled.

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
 
Once more, anytime `set_installation_options.sh` is modified, we need to reload `(.local)` by rerunning `start_cocoa.sh`. Then, run the following commands:

      cd ./cocoa/Cocoa

and

      source start_cocoa.sh # even if (.local) is already active, users must run start_cocoa.sh again to update bash environment values
      
Now, to download and compile all Cosmolike projects, type
 
      source ./installation_scripts/setup_cosmolike_projects.sh   # download all cosmolike projects  

and

      source ./installation_scripts/compile_all_projects.sh       # compile  all cosmolike project

> [!NOTE]
> In case users need to rerun `setup_cocoa.sh` (or `setup_cosmolike_projects.sh`) , Cocoa will not download previously installed cosmolike projects (this avoids loss of uncommitted work), unless the following key is set on `set_installation_options.sh`
>
>     [Adapted from Cocoa/set_installation_options.sh shell script]
>     #export OVERWRITE_EXISTING_COSMOLIKE_CODE=1 # dangerous (possible loss of uncommitted work)
>                                                 # if unset, users must manually delete cosmolike projects

In case users only want to compile a single Cosmolike project (let's say the `roman_real` project)

      source ./projects/roman_real/scripts/compile_roman_real.sh
     
## :interrobang: FAQ: How can users run Cocoa with Docker? <a name="appendix_jupyter_whovian"></a>

We provide the Docker image [whovian-cocoa](https://hub.docker.com/r/vivianmiranda/whovian-cocoa) to facilitate installing Cocoa on Windows and macOS. This appendix assumes that users have already installed the Docker Engine on their local PC. For instructions on installing the Docker engine on specific operating systems, refer to [Docker's official documentation](https://docs.docker.com/engine/install/). 

 **Step :one:**: Create a folder and go to the location on the host computer where you want to provide access to the Docker container, as shown below. 

     mkdir -p cocoa_docker
     cd ./cocoa_docker

 **Step :two:**: Download the Docker image *whovian-cocoa*, name the associated container `cocoa2025` (flag `--name cocoa2025` in the command below), and run the container for the first time, type:

    docker run --platform linux/amd64 --hostname cocoa --name cocoa2025 -it -p 8888:8888 -v $(pwd):/home/whovian/host/ -v ~/.ssh:/home/whovian/.ssh:ro vivianmiranda/whovian-cocoa:thin

This is a large image, approximately 13GB, as it already contains the conda cocoa environment. Users can now proceed to the [Installation and Compilation of external modules](#cobaya_base_code) section to continue the installation. 

> [!TIP]
> Once installation is complete, the user must learn how to **start** and **exit** the Docker container. Assuming the user maintained the container name `cocoa2025` set on the flag `--name cocoa2025`, type:
>    
>      docker start -ai cocoa2025
>
>  to restart the container.

> [!TIP]
> To run Jupyter Notebooks within the *whovian-cocoa* Docker container installed on a local machine, type the following command:
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

## :interrobang: FAQ: How can users run Cocoa on Google Colab? <a name="overview_google_colab"></a>

[Google Colab](https://colab.research.google.com/) provides a convenient platform for users to run MCMCs, likelihood minimizations, and profiles, as long as Machine-Learning Emulators are used to compute the data vectors. In the repository [CoCoAGoogleColabExamples](https://github.com/CosmoLike/CoCoAGoogleColabExamples), we provide a few examples along with explanatory notes. 

Installing Cocoa requires time and also strains our limited Git-LFS quota, which is especially relevant given that **the entire local drive is wiped when a Colab notebook is disconnected**. To prevent this problem, we provide instructions for saving and loading Cocoa immediately after the initial installation. 

There are a few differences users should be aware of when running Cocoa on Google Colab.

  - Running Collab Notebook for the first time

    - **Cell :one:**: Connect the notebook to your Google Drive account (will be important later)

          from google.colab import drive
          drive.mount('/content/drive')

    Below, we provide instructions for installing Cocoa. *Google Colab does provide direct terminal access if users prefer to follow the standard installation procedure*
    
    - **Cell 2️⃣**: Install Miniforge (Similar to our documentation in section [FAQ: How can users install Conda?](#overview_miniforge))

          %%bash
          export CONDA_DIR="/content/conda"
          mkdir "${CONDA_DIR:?}"
          curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
          /bin/bash Miniforge3-$(uname)-$(uname -m).sh -f -b -p "${CONDA_DIR:?}"
          /bin/bash
          source $CONDA_DIR/etc/profile.d/conda.sh \
                && conda config --set auto_update_conda false \
                && conda config --set show_channel_urls true \
                && conda config --set auto_activate_base false \
                && conda config --prepend channels conda-forge \
                && conda config --add allowlist_channels conda-forge \
                && conda config --set channel_priority strict \
                && conda init bash
          source ~/.bashrc

    - **Cell 3️⃣**: Install Conda cocoa env (similar to our documentation in section [Installation of core packages](#required_packages_conda))

          %%bash
          source "/content/conda/etc/profile.d/conda.sh"
          conda activate lockenv
          conda-lock install -n cocoa cocoapy310-linux.yml
          conda activate cocoa 
          ln -s "${CONDA_PREFIX}"/bin/x86_64-conda_cos6-linux-gnu-gcc "${CONDA_PREFIX}"/bin/gcc
          ln -s "${CONDA_PREFIX}"/bin/x86_64-conda_cos6-linux-gnu-g++ "${CONDA_PREFIX}"/bin/g++
          ln -s "${CONDA_PREFIX}"/bin/x86_64-conda_cos6-linux-gnu-gfortran "${CONDA_PREFIX}"/bin/gfortran
          ln -s "${CONDA_PREFIX}"/bin/x86_64-conda-linux-gnu-gcc-ar "${CONDA_PREFIX}"/bin/gcc-ar
          ln -s "${CONDA_PREFIX}"/bin/x86_64-conda-linux-gnu-gcc-ranlib "${CONDA_PREFIX}"/bin/gcc-ranlib
          git-lfs install

    - **Cell 4️⃣**: Clone CoCoA (similar to our documentation in section [Installation and Compilation of external modules](#cobaya_base_code))

          %%bash
          source "/content/conda/etc/profile.d/conda.sh"
          conda activate cocoa                                  
          git clone https://github.com/CosmoLike/cocoa.git --branch v4.0 cocoa # users can adjust this line

    - **Cell 5️⃣**: run `setup_cocoa.sh`

           %%bash
           source "/content/conda/etc/profile.d/conda.sh" 
           conda activate cocoa
           cd ./cocoa/Cocoa/
           source setup_cocoa.sh

    - **Cell 6️⃣**: run `compile_cocoa.sh`

          %%bash
          source "/content/conda/etc/profile.d/conda.sh"
          conda activate cocoa
          cd ./cocoa/Cocoa/
          source compile_cocoa.sh

    - **Cell 7️⃣**: Save Cocoa on Drive (does not work with local runtime)

          %%bash
          DEST="/content/drive/MyDrive/ColabBackups"
          ARCHIVE="$DEST/colab_basic_cocoa.tar.gz"
          if [[ -f "$ARCHIVE" ]]; then
            echo "Backup already exists: $ARCHIVE — skipping."
            exit 0
          fi
          mkdir -p "$DEST"
          tar -czf "$DEST/colab_basic_cocoa.tar.gz" \
            --exclude='/content/drive' \
            --exclude='**/__pycache__' \
            --exclude='**/.ipynb_checkpoints' \
            /content
          echo "Created: $ARCHIVE"

  - Running Collab Notebook with Cocoa pre-installed, loaded from Drive (does not work with local runtime)

    - **Cell :one:**: Connect the notebook to your Google Drive account

          from google.colab import drive
          drive.mount('/content/drive')

    - **Cell 2️⃣**: Load Cocoa from Drive

          %%bash
          DEST="/content/drive/MyDrive/ColabBackups"
          ARCHIVE="$DEST/colab_basic_cocoa.tar.gz"
          SENTINEL="/content/conda/etc/profile.d/conda.sh"  # exists when your env is restored
          if [[ -e "$SENTINEL" ]]; then
            echo "Found $SENTINEL — environment already restored. Skipping untar."
            exit 0
          fi
          test -f "$ARCHIVE"
          ARCHIVE="/content/drive/MyDrive/ColabBackups/colab_basic_cocoa.tar.gz"
          tar -xzf "$ARCHIVE" -C /

> [!Note]
> From now on, users must start every subsequent shell with 
>
>        %%bash
>        source "/content/conda/etc/profile.d/conda.sh"
>        conda activate cocoa`
>        cd ./cocoa/Cocoa/
>        source start_cocoa.sh


- Saving/Loading checkpoints

  Not reserving time to copy the `/content` folder to the user's Google Drive —an expensive operation —can result in up to 24 hours of lost computation. To prevent such a catastrophe, the code below creates and loads *checkpoints* that users can add after computationally intensive cells.

  *This solution is not valid when running Colab with local runtime* (see [Google documentation](https://research.google.com/colaboratory/local-runtimes.html) for additional information on how to link notebooks to local resources). The good news here is that local storage is persistent, so there is no need to create backups on Google Drive.
  
  - Saving checkpoints: compress and copy the `/content` folder from the local disk to the user's Drive
    
        %%bash
        ROOT="colab_name_notebook"
        DEST="/content/drive/MyDrive/ColabBackups"
        mkdir -p "$DEST"
        ARCHIVE="$DEST/$ROOT_$(date +%F_%H-%M).tar.gz"
        tar -czf "$ARCHIVE" \
            --exclude='/content/drive' \
            --exclude='**/__pycache__' \
            --exclude='**/.ipynb_checkpoints' \
            /content
        echo "Created: $ARCHIVE"

  - Loading checkpoints: decompress and copy the `/content` folder from the user's Drive to the local disk
    
        %%bash
        SENTINEL="/content/conda/etc/profile.d/conda.sh"  # exists when your env is restored
        if [[ -e "$SENTINEL" ]]; then
          echo "Found $SENTINEL — environment already restored. Skipping untar."
          exit 0
        fi
        ARCHIVE="CHECKPOINT_FILE"
        test -f "$ARCHIVE"
        tar -xzf "$ARCHIVE" -C /

## :interrobang: FAQ: How can users install Conda? <a name="overview_miniforge"></a>

**Step :one:**: Download and run the Miniforge installation script. 

    export CONDA_DIR="/gpfs/home/XXX/miniforge" # replace this string!

and

    mkdir "${CONDA_DIR:?}"

and

    curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"

and

    /bin/bash Miniforge3-$(uname)-$(uname -m).sh -f -b -p "${CONDA_DIR:?}"

and

    /bin/bash 

**Step :two:**: After installation, users must source the conda configuration file, as shown below:

    source $CONDA_DIR/etc/profile.d/conda.sh \
          && conda config --set auto_update_conda false \
          && conda config --set show_channel_urls true \
          && conda config --set auto_activate_base false \
          && conda config --prepend channels conda-forge \
          && conda config --add allowlist_channels conda-forge \
          && conda config --set channel_priority strict \
          && conda init bash

> [!Note]
> The strict use of Miniconda and conda-forge packages appears to mitigate the significant licensing issue involving Anaconda packages and academic/research institutions.

**Step :three:**: After running this command, you will see a message in the terminal that ends with the statement *For changes to take effect, close and re-open your current shell*. Then, type

    source ~/.bashrc

After that, the `conda` command will be available.

## :interrobang: FAQ: How can users set the appropriate environment for ML? <a name="ml_emulators"></a>

Commenting out the environmental flags below *before running* `setup_cocoa.sh` will enable the installation of machine-learning-related libraries via pip.  

    [Adapted from Cocoa/set_installation_options.sh shell script] 
    # ------------------------------------------------------------------------------
    # If not set, pip_core_packages.sh will install several ML package
    # ------------------------------------------------------------------------------
    #export IGNORE_EMULATOR_GPU_PIP_PACKAGES=1
    (...)

If users have already run `setup_cocoa.sh`, run the command.

    source start_cocoa.sh # even if (.local) is already active, users must run start_cocoa.sh again to update bash environment values
   
and

    source ./installation_scripts/setup_pip_core_packages.sh
              
## :interrobang: FAQ: How can developers push changes to the Cocoa main branch? <a name="push_main"></a>

Until recently, Cocoa development was a bit unstructured. Developers could push directly to the `main` branch, and small commits to the main were not discouraged. Such flexible development rules will soon change when `v4.0` leaves the beta phase. We will protect the `main` branch by requiring every push to be reviewed by Cocoa's leading developers. 

Our new philosophy establishes that **a commit in the main branch should contain an atomic change that moves code from one working state to another working state with meaningful and well-tested improvements**. 

Therefore, developers should propose changes to the `main` branch in larger chunks (*via squash commits*), as shown below.

> [!NOTE]
> In a developer branch, users are encouraged to make small commits so that work can be tracked easily. Our policy regarding squash atomic changes only applies to the main branch.

- :interrobang: **How to apply squash commits?**
  
**Step :one:**: create a development branch. Do not call the development branch `dev`, as `dev` is reserved for work done by the leading Cocoa developers. For concreteness, let's name the new branch `xyzdev`

   
    git switch -c xyzdev  # run on the main branch

> [!TIP]
The use of developers' initials followed by `dev` makes the branch easily identifiable.

> [!TIP]
> If the branch `xyzdev` already exists, use the `git switch` command without the `-c` flag. 

**Step :two:**: develop the proposed changes. We advise developers to commit frequently. In your branch, a commit does not need to be *atomic*, meaning the commit is not require to change the code from one well-tested working state to another well-tested state. Developers can push to the server via the command.

    git push -u origin xyzdev   # run on the xyzdev branch

**Step :three:**: Once the developers have created an atomic, meaningful, and well-tested improvement to Cocoa, the developer needs to merge any subsequent changes made in `main`.

    git merge main   # run on the xyzdev branch

This step may create conflicts that must be addressed before step four. If there are conflicts, the best way to identify the files that need to be edited is via the command

    git diff --name-only --diff-filter=U # shows only files with conflict

**Step :four:**: Once the developers have merged recent changes made on the `main` branch, they must push to the main branch the modifications made on the `xyzdev` branch by first **squashing all your changes into a single commit**, as shown below

    git switch main  # run on the xyzdev branch

and

    git merge --squash xyzdev   # run on the main branch

Again, if there are conflicts, use the command. 

    git diff --name-only --diff-filter=U
   
to identify files that need to be edited to solve the conflict state. Finally, make the merger commit (please solve all conflicts before commit)
 
    git commit -m "merge xyzdev branch"  # run on the main branch

and

    git push origin main  # run on the main branch

## :interrobang: FAQ: How can developers develop from a Git tag? <a name="dev_from_tag"></a>

A useful Git hack is related to developing Cocoa from a Git tag. We reserve Git tags to set milestones in our development, so they serve as good starting points for coding localized new features (e.g., changes to a file that other developers have not recently modified) or bug fixes.

**Step :one: (optional)** If the developer has cloned the repository using the `https` URL address, then change the URL to the SSH-key-based address (if the developer has previously uploaded a public key to their GitHub account)

    git remote set-url origin git@github.com:CosmoLike/cocoa.git # that would allow users to push without typing a password

**Step :two:** Move the detached state to a new local branch via the command

    git switch -c xyzlocdev

Now, all commits will be associated with this local branch. 

> [!TIP]
The use of developers' initials followed by `dev` makes the branch easily identifiable.

**Step :three:** The developer has two options at the end of development. They can **either** create a new remote branch

    git push origin xyzlocdev # run on the xyzlocdev branch

**or** they can fetch and download the remote `xyzdev` branch, which will later absorb the changes made on `xyzlocdev`

    git switch -c xyzdev origin/xyzdev # run on the xyzlocdev branch

Finally, the developer needs to merge the changes made on `xyzlocdev`.

    git merge --squash xyzlocdev # run on the xyzdev branch

If this merge does not create any merge conflicts, type

    git push origin xyzdev # run on the xyzdev branch

## :interrobang: FAQ: How can users fork Cocoa correctly? <a name="fork_cocoa_tag"></a>

Cocoa manages a few core repositories that must all be forked to provide users with developer-editing capabilities across all our core projects. There are a few well-known Cocoa forks in addition to the main code hosted on the [Cosmolike organization](https://github.com/CosmoLike). They are

* Forks hosted by [SBU-Cosmolike organization](https://github.com/SBU-COSMOLIKE)
	- Here is where Prof. Miranda and her group at Stony Brook are developing Cocoa, while the University of Arizona leads development at Cosmolike org. 
    	+ Code is more unstable at SBU-Cosmolike, as the *atomic requirements* described in the appendix [FAQ: How can developers push changes to the Cocoa main branch?](#push_main) are not so strictly observed.
     	+ There are differences in commit counts between the SBU-Cosmolike and Cosmolike organizations, as changes to the Cosmolike org are typically pushed via squash commits. 
	- Some cutting-edge features are developed within SBU-Cosmolike, especially the ones related to machine-learning emulators (users should use forks to develop cutting-edge code!).
    - We do not provide any support for people outside Stony Brook when using code from the SBU-Cosmolike org, but they are still open-source.
    - Documentation within the SBU-Cosmolike org assumes users can clone repositories using SSH key authentication. So, non-members of the SBU-Cosmolike org must adapt the URL addresses to install Cocoa successfully.  
	
* Forks hosted by [Roman HLIS Cosmology PIT](https://github.com/Roman-HLIS-Cosmology-PIT)
	- Here is where we port the code that is intended to be branded as official releases by the [HLIS Cosmology Project Infrastructure Team](https://roman-hlis-cosmology.caltech.edu). 
	
Assuming, for concreteness, that users want to fork relevant repositories from SBU-Cosmolike (**not advised**, better to clone from Cosmolike org), we list below the settings that must be adjusted in `set_installation_options.sh`. 

    [Adapted from Cocoa/set_installation_options.sh shell script]    
    export COSMOLIKE_URL="git@github.com:SBU-COSMOLIKE/cocoa-cosmolike-core.git"  # main repository with cosmolike code
    
    (...)   
    
    # Below is the code that implements ML-learning emulators
    export EMULTRF_URL="git@github.com:SBU-COSMOLIKE/emulators_code.git"
    # Below are the files that implement ML-learning emulators
    export EMULTRF_DATA_URL="https://github.com/SBU-COSMOLIKE/emulators_data_lcdm.git" 
    
    (...)
    
    # Cosmolike projects below
    export LSST_Y1_URL="git@github.com:SBU-COSMOLIKE/cocoa_lsst_y1.git"
    (...)
    export DES_Y3_URL="git@github.com:SBU-COSMOLIKE/cocoa_des_y3.git" # This is DES-Y3 redmagic
	(...) 
    export ROMAN_FOURIER_URL="git@github.com:SBU-COSMOLIKE/cocoa_roman_fourier.git" 
    (...)
    export ROMAN_REAL_URL="git@github.com:SBU-COSMOLIKE/cocoa_roman_real.git"
    (...)
    export DESXPLANCK_URL="git@github.com:SBU-COSMOLIKE/cocoa_desy1xplanck.git"  # This is DES-Y3 maglim


## :interrobang: FAQ: Note on Planck-2018 low-ell SimAll EE and Gibbs TT likelihoods <a name="planck2018lowell"></a>

Cocoa does not adopt Cobaya's Python reimplementation of SimAll EE and Gibbs TT likelihoods. Therefore, we patch the files. 

    cobaya/cobaya/likelihoods/planck_2018_lowl/TT.py
    cobaya/cobaya/likelihoods/planck_2018_lowl/EE.py

to call the original Planck-2018 Clik likelihoods, as shown below.

<img width="1146" height="498" alt="Screenshot 2025-10-20 at 8 38 56 PM" src="https://github.com/user-attachments/assets/4ab8beff-bda1-4bf8-b2c6-56f5501e3773" />
<img width="1236" height="514" alt="Screenshot 2025-10-20 at 8 39 04 PM" src="https://github.com/user-attachments/assets/22b21a5d-c716-45dc-949c-59f31be91def" />

This ensures backward consistency in our code, as `TT.py` and `EE.py` used to point to Planck-2018 Clik before Cobaya's authors moved them to `EE_clik.py` and `TT_clik.py`.

## :interrobang: FAQ: How can users improve our Bash/C/C++ knowledge? <a name="lectnotes"></a>

A working knowledge of Python is required to understand the Cobaya framework at the developer level. Users must also be familiar with the Bash language to understand Cocoa's scripts. Proficiency in C and C++ is also needed to manipulate Cosmolike and the C++ Cobaya-Cosmolike C++ interface. Finally, users need to understand the Fortran-2003 language to modify CAMB.

Learning all these languages can be overwhelming, so to enable new users to do research that demands modifications on the inner workings of these codes, we include [here](cocoa_installation_libraries/LectNotes.pdf) a link to approximately 600 slides that provide an overview of Bash (slides ~1-137), C (slides ~138-371), and C++ (slides ~372-599). In the future, we aim to add lectures about Python and Fortran. 


