# Subroutine for changing l_max of Planck Likelihood

The environment is much easiler in `Cocoa`, so I’ll assume we use cocoa. The first thing to do is `conda activate cocoa` and `source start_cocoa`, just in the normal way as we are going to use cobaya.

The tools are for python2, so we need to make a little modification. Go to `Cocoa/external_modules/code/planck/code/plc_3.0/plc-3.01/src/python/tools` and change line 121 of `clik_change_lrange` from `siginv.shape=(siginv.shape[0]**.5),-1)` to 

```
siginv.shape=(int(siginv.shape[0]**.5),-1)
```

Then we use the script to compile planck again, so `source installation_scripts/compile_plank`.

The usage is following the python file, for example we can do : 

```
clik_change_lrange Cocoa/external_modules/data/planck/plc_3.0/hi_l/plik/plik_rd12_HM_v22_TT.clik -1 2400 Cocoa/external_modules/data/planck/plc_3.0/hi_l/plik/plik_rd12_HM_v22_TT_lmax_2400.clik
```

which will create a new likelihood with the new l range. Here -1 means keep the l_min the same, and change l_max to 2400. 

Then we can use `clik_print` to check the new l range. 

Please note that the l range passed to this script is not exactly the same as the outcome because the likelihood is binned here, but the real range will be printed once you call the script.

