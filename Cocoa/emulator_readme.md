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

In order to train the 

```
emulator:
  io:
      # Directory where to save models, data and chains
      savedir: ./projects/des_y3/emulator_output/baseline_3x2      
      # whether to save the trining data
      save_train_data: true             
      # whether to save the intermediate emulator models
      save_intermediate_model: true
  shear_calib:
      # whether to save the intermediate emulator models
      mask: ./projects/des_y3/data/emu_files/shear_calib_mask.npy
      prior_std: 0.005, 0.005, 0.005, 0.005
  baryons:
      n_pcas_baryon: 0
      prior_Q1: -3., 12.
  training:
      dv_fid: ./projects/des_y3/data/emu_files/dv_fid.txt
      n_train_iter: 6
      n_lhs: 4800
      n_resample: 9600
      emu_type: nn
      batch_size: 8
      n_epochs: 100
  sampling:
      scalecut_mask: ./projects/des_y3/data/3x2pt_baseline.mask
      n_mcmc: 10000
      n_emcee_walkers: 120
      n_burn_in: 5000
      n_thin: 5
      temper0: 0.05
      temper_increment: 0.1
```

### Scale cut

Our emulator learns the full data vector without any mask. However, while sampling from the tempered posterior, we need to impose a scale cut. To ensure this, we always need to set the mask file in the `.dataset` file to `ones.mask` where there is no scale cut. 

To set the mask file while sampling, you can set the requisite mask file in the `.yaml` file under:

```
emulator:
  sampling:
      scalecut_mask: ./projects/des_y3/data/3x2pt_baseline.mask
      ....
```

### Fast parameters

#### Shear calibration

#### Baryons

`.dataset`

