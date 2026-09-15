# The `example` project

This folder contains the example YAML files used throughout [Cocoa's main README](../../../README.md),
which also provides the appropriate `mpirun` commands to run them. The table below lists, for each
example, the datasets (likelihoods) it runs and the theory code that computes the observables:
`CAMB`, `CLASS`, `EMUL` (Cocoa's neural-network emulators), or `EMUL2` (hybrid emulator, which
emulates only the Boltzmann outputs). This folder currently has no `EMUL2` examples; for those,
see the Cosmolike project READMEs (e.g., `projects/lsst_y1` and `projects/des_y3`).

| Example | Datasets | Theory |
| ------- | -------- | ------ |
| [EXAMPLE_EVALUATE1.yaml](EXAMPLE_EVALUATE1.yaml) | Planck 2018 TTTEEE (Plik and Plik-lite); Planck 2018 low-ℓ TT and EE; Pantheon, Pantheon+, and DES-Y5 SNe; DESI DR1 and DR2 BAO; ACT DR6 lensing | CAMB |
| [EXAMPLE_EVALUATE2.yaml](EXAMPLE_EVALUATE2.yaml) | Planck 2018 TTTEEE (Plik-lite); Planck 2018 low-ℓ TT and EE; Pantheon, DES-Y5, and Roman (conservative and optimistic) SNe; DESI DR2, 6dFGS, SDSS DR7 MGS, and SDSS DR12 BAO; ACT DR6 lensing; SO LAT TTTEEE (simulated, MFLike) | CAMB |
| [EXAMPLE_EVALUATE3.yaml](EXAMPLE_EVALUATE3.yaml) | SPT-3G Y1 TTTEEE; ACT DR6 lensing | CAMB |
| [EXAMPLE_EVALUATE4.yaml](EXAMPLE_EVALUATE4.yaml) | ACT DR4 TTTEEE (lite) | CAMB |
| [EXAMPLE_EVALUATE5.yaml](EXAMPLE_EVALUATE5.yaml) | Planck 2018 TTTEEE (Plik) | CLASS |
| [EXAMPLE_EVALUATE6.yaml](EXAMPLE_EVALUATE6.yaml) | SO LAT TTTEEE (simulated, MFLike) | CAMB |
| [EXAMPLE_EVALUATE7.yaml](EXAMPLE_EVALUATE7.yaml) | ACT DR6 TTTEEE (CMB-only); ACT DR6 lensing; Planck 2018 low-ℓ EE (SRoll2) | CAMB |
| [EXAMPLE_EVALUATE8.yaml](EXAMPLE_EVALUATE8.yaml) | ACT DR6 TTTEEE (MFLike); ACT DR6 lensing | CAMB |
| [EXAMPLE_EVALUATE9.yaml](EXAMPLE_EVALUATE9.yaml) | Planck 2018 TTTEEE (Plik); Planck 2018 low-ℓ EE | CAMB |
| [EXAMPLE_EVALUATE10.yaml](EXAMPLE_EVALUATE10.yaml) | Planck 2018 TTTEEE (CamSpec 2021); Planck 2018 low-ℓ TT and EE | CAMB |
| [EXAMPLE_EVALUATE11.yaml](EXAMPLE_EVALUATE11.yaml) | Planck 2020 TTTEEE (HiLLiPoP); Planck 2018 low-ℓ TT; Planck 2020 low-ℓ EE (LoLLiPoP) | CAMB |
| [EXAMPLE_MCMC1.yaml](EXAMPLE_MCMC1.yaml) | Pantheon SNe; Planck 2018 TT (Plik) | CAMB |
| [EXAMPLE_MCMC2.yaml](EXAMPLE_MCMC2.yaml) | Roman (conservative) SNe; Planck 2018 TT (Plik) | CAMB |
| [EXAMPLE_POLY1.yaml](EXAMPLE_POLY1.yaml) | Pantheon SNe; Planck 2018 TT (Plik) | CAMB |
| [EXAMPLE_POST1.yaml](EXAMPLE_POST1.yaml) | ACT DR4 TTTEEE (lite) | CAMB |
| [EXAMPLE_EMUL_EVALUATE1.yaml](EXAMPLE_EMUL_EVALUATE1.yaml) | Planck 2018 TTTEEE (Plik-lite); Planck 2018 low-ℓ TT and EE; Pantheon, DES-Y5, and Roman (conservative and optimistic) SNe; DESI DR2, 6dFGS, SDSS DR7 MGS, and SDSS DR12 BAO; ACT DR6 lensing; SO LAT TTTEEE (simulated, MFLike) | EMUL |
| [EXAMPLE_EMUL_EVALUATE2.yaml](EXAMPLE_EMUL_EVALUATE2.yaml) | Planck 2018 TTTEEE (Plik-lite); Planck 2018 low-ℓ TT and EE; Pantheon, Pantheon+, and DES-Y5 SNe; DESI DR1 and DR2 BAO; ACT DR6 lensing | EMUL |
| [EXAMPLE_EMUL_EVALUATE2_CP.yaml](EXAMPLE_EMUL_EVALUATE2_CP.yaml) | Planck 2018 TTTEEE (Plik-lite); Planck 2018 low-ℓ TT and EE; ACT DR6 lensing | EMUL (CosmoPower) |
| [EXAMPLE_EMUL_EVALUATE3.yaml](EXAMPLE_EMUL_EVALUATE3.yaml) | DES-Y5 SNe; DESI DR2 BAO; ACT DR6 lensing; SO LAT TTTEEE (simulated, MFLike) | EMUL |
| [EXAMPLE_EMUL_MCMC1.yaml](EXAMPLE_EMUL_MCMC1.yaml) | Planck 2018 TTTEEE (Plik-lite); Planck 2018 low-ℓ TT and EE (SRoll2); DES-Y5 SNe; DESI DR2 BAO; ACT DR6 lensing | EMUL |
| [EXAMPLE_EMUL_MCMC2.yaml](EXAMPLE_EMUL_MCMC2.yaml) | Planck 2018 TTTEEE (Plik-lite); Planck 2018 low-ℓ TT and EE (SRoll2) | EMUL |
| [EXAMPLE_EMUL_MCMC3.yaml](EXAMPLE_EMUL_MCMC3.yaml) | Planck 2018 TTTEEE (Plik-lite); Planck 2018 low-ℓ TT and EE (SRoll2); DES-Y5 SNe; DESI DR2 BAO; ACT DR6 lensing | EMUL |
| [EXAMPLE_EMUL_POLY1.yaml](EXAMPLE_EMUL_POLY1.yaml) | Planck 2018 TTTEEE (Plik-lite); Planck 2018 low-ℓ TT and EE (SRoll2); DES-Y5 SNe; DESI DR2 BAO; ACT DR6 lensing | EMUL |

The Python scripts in this folder (`EXAMPLE_MINIMIZE1.py`, `EXAMPLE_PROFILE1.py`, and the
`EXAMPLE_EMUL_*.py` Emcee/Nautilus/Minimize/Profile/Scan examples) are documented in the
main README; the plain scripts use `CAMB`, and the `EMUL_` scripts use Cocoa's emulators.
