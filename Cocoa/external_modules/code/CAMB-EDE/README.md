# CAMB-EDE
CAMB code with the addition of an Early Dark Energy component. EDE can be modelled as an effective fluid or as a scalar field, which can have either the AxionEDE or the monomial potentials.
## Installation
 - In your terminal and in the desired directory, clone this repository with `git clone git@github.com:CoupleDE-UNESP/CAMB-EDE.git`
 - Go to the CAMB-EDE folder with `cd CAMB-EDE` and install it using `pip install -e .`
## Command-Line Usage
In your terminal, use `camb <inifile.ini>` to run the cosmology set in the inifile. In the `inifiles/` directory, you will see a template for dark energy `ede.ini`.
## Python Usage
In the `docs/` directory there is a Jupyter notebook `cambEDEtest.ipynb` showing how to configure Early Dark Energy. You can set the specific model, scalar field potential parameters, use the <img src="https://render.githubusercontent.com/render/math?math=z_c, f_{de}"> parametrization,  There, you will see the parameters needed to set in MCMC chains. The Early Dark Energy model is specified within the `camb.set_params()` function. Its arguments are:
 - `dark_energy_model = "EarlyDarkEnergy"` chooses scalae field Early Dark Energy + a cosmological constant as your model, whereas `dark_energy_model = "AxionEffectiveFluid"` chooses the effective fluid parametrization.
 - For the effective fluid EDE, you need to set: `zc`, the critical transition redshift; `fde_zc`, the fractional amount of EDE at the critical redshift; and `w_n`, the EoS parameter EDE will transition to.
 - For the scalar field EDE, `which_potential` chooses the scalar field potential to use. `which_potential = 1` uses a monomial potential, while `which_potential = 2` uses the AxionEDE potential.
 - For the scalar field EDE, `use_zc` specifies whether you want to parametrize your model with the potential parameters or with the peak energy contribution `fde_zc` and redshift `zc`. If `use_zc = False`, you will need to specify the potential parameters. If `use_zc = True`, two of the potential parameters will be substituted for `zc` and `fde_zc`.
 - The monomial potential is parametrized as `V = V0 * phi^(2n)`. If `use_zc =  False`, you will need to insert `V0` (in reduced Planck mass, Mpl^4), `n` and `initial_phi`. If `use_zc = True`, you have to insert `zc`, `fde_zc` and `n`.
 - The AxionEDE potential is parametrized as `V = m^2 * f2 * (1 - cos(phi/f))^n`. If `use_zc =  False`, you will need to insert `m` (in reduced Planck mass, Mpl), `f` (also in Mpl), `n` and `theta_i`, the initial field displacement from the bottom. If `use_zc = True`, you have to insert `zc`, `fde_zc`, `n` and `theta_i`.
## Configuring an inifile
In the inifiles directory, there are some templates. To use Early Dark Energy, set `dark_energy_model = earlydarkenergy`
