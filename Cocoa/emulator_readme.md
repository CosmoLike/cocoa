Additional packages to install:
    - `pytorch`
    - `pyDOE`
    - `tqdm`
    - `george`

Always use: 

`export OMP_NUM_THREADS=1`. Otherwise interference with pytorch. Refer to Vivian and Evan.

To run the code:

```
mpirun -n 16 --mca btl tcp,self $PYTHON3 train_emulator.py ./projects/lsst_y1/template_yamls/emulator/EXAMPLE.yaml 
```

## Configure using the yaml file

In order to train the emulator, you can set up the requisite config file by adding a `emulator` section in yout config file. A typical `emulator` configuration will look like the following: 

```
emulator:
  io:
      savedir: ./projects/des_y3/emulator_output/baseline_3x2      
      save_train_data: true             
      save_intermediate_model: true
  training:
      dv_fid: ./projects/des_y3/data/emu_files/dv_fid.txt
      chi_sq_cut: 8e+4
      n_train_iter: 6
      n_lhs: 4800
      n_resample: 9600
      emu_type: nn
      batch_size: 32
      n_epochs: 100
  shear_calib:
      mask: ./projects/des_y3/data/emu_files/shear_calib_mask.npy
      prior_std: 0.005, 0.005, 0.005, 0.005
  baryons:
      n_pcas_baryon: 0
      prior_Q1: -3., 12.
  sampling:
      scalecut_mask: ./projects/des_y3/data/3x2pt_baseline.mask
      n_mcmc: 10000
      n_emcee_walkers: 120
      n_burn_in: 5000
      n_thin: 5
      temper0: 0.05
      temper_increment: 0.1
```

### I/O

The io section points where to save the emulator models and whether or not to save the intermediate models / training data.
```
  io:
      savedir: WHERE TO SAVE THE MODELS AND TRAINING DATA
      save_train_data: WHETHER TO SAVE THE TRAINING DATA.
      save_intermediate_model: WHETHER TO SAVE THE INTERMEDIATE EMULATOR MODELS. IF SET TO FALSE, ONLY THE FINAL MODEL IS SAVED.
```    

### Training

The `training` section of the yaml file set the hyperparameters required for training the emulators.  

```  
training:
  dv_fid:       PATH OF THE FIDUCIAL DATA VECTOR REQUIRED FOR OUR EMULATOR TRAINING. 
  chi_sq_cut:   INCLUDE TRAINING SAMPLES WHICH chi_sq < chi_sq_cut. Default value is 1e+5.
  n_train_iter: NUMBER OF ITERATIONS REQUIRED FOR TRAINING THE EMULATOR.
  n_lhs:        NUMBER OF INITIAL LATIN HYPERCUBE SAMPLES.
  n_resample:   NUMBER OF TRAINING DATA ADDED IN EACH ITERATION
  emu_type:     TYPE OF EMULATOR. EITHER 'nn' or 'gp'
  batch_size:   BATCH SIZE USED FOR OPTIMIZATION OF THE NEURAL NETWORK
  n_epochs:     NUMBER OF TRAINING EPOCH FOR NN TRAINING
```        


### Configuring the fast parameters

In our emulator, the impact of shear calibration bias and baryons are added exactly a-posteriori. That is, they are not used for the emulation.

#### Shear calibration

In order to make sure that the shear calibration parameters are not included in the emulator, make sure that they are set to 0 in the `params` section in the yaml file. E.g, 

```
params:
  ....
  ....
  LSST_M1:
    value: 0.
    latex: m_\mathrm{LSST}^1
  LSST_M2:
    value: 0. 
    latex: m_\mathrm{LSST}^2
  LSST_M3:
    value: 0.
    latex: m_\mathrm{LSST}^3
  LSST_M4:
    value: 0.
    latex: m_\mathrm{LSST}^4
  LSST_M5:
    value: 0.
    latex: m_\mathrm{LSST}^5
```

The shear calibration is set to zero while training. The impact of shear calibration bias is to rescale the data vectors as,   
$$\xi^{ij}_{\pm} \rightarrow (1 + m_i) (1 + m_j) \xi^{ij}_{\pm}$$ and $$\gamma^{ij}_t \rightarrow (1 + m_i) \gamma^{ij}_t$$. 

A script to find the elements of the data vector that scale accordingly is provided in the `dv_grad.py` script.

Once calculated, you set the scaling with the shear calibration `mask`. 

While sampling, we impose a prior on the shear calibration parameters. It is assumed to be a Gaussian prior with some standard deviation. The standard deviation of the shear calibration prior is set by the `prior_std` field. 

```
shear_calib:
  mask: PATH to the shear calibration mask
  prior_std: 
```

#### Baryons

Similar to the shear calibration bias, the impact of baryons is added a-posteriori. The Baryon PCAs are precalculated and set in the `.dataset` file. While sampling with the emulator, we set the number of PCA modes we want to add. This is set using the `n_pcas_baryon` field. Set it to 0 if you do not want to include baryons.

The priors on the PCA modes are set by using the `prior_Q1`. We assume a flat prior between some minimum and maximum value of the parameters. E.g, if we set, `prior_Q1: -3., 12.`, we impose a prior between -3 and 12 for the first PC mode.

```
baryons:
  n_pcas_baryon: 1
  prior_Q1: -3., 12.
```


### Sampling

The sampling field can be used for two different purposes: i) While training the emulator, we sample from the tempered posterior to get new training data. ii) To sample from the posterior using the trained emulator.

You can use the `sample_emulator.py` script for sampling using the trained emulator. You can sample using the trained emulator using the following command:

```
$PYTHON3 sample_emulator.py ./projects/lsst_y1/template_yamls/emulator/EXAMPLE.yaml 
```

```
  sampling:
      scalecut_mask: ./projects/des_y3/data/3x2pt_baseline.mask
      n_mcmc: 10000
      n_emcee_walkers: 120
      n_burn_in: 5000
      n_thin: 5
      temper0: 0.05
      temper_increment: 0.1
```

#### Scale cut

Our emulator learns the full data vector without any mask. However, while sampling from the tempered posterior, we need to impose a scale cut. To ensure this, we always need to set the mask file in the `.dataset` file to `ones.mask` where there is no scale cut. 

To set the mask file while sampling, you can set the requisite mask file in the `.yaml` file under:


#### Temper value

In order to train the emulator, we sample the training samples from the tempered posterior. The value of tempering is set to some initial value, `alpha = temper0` and then increased by `temper_increment` at each iteration. So at a given `iteration`, its tempering value is 

```
    alpha = temper0 + (iteration - 1) * temper_increment
```    