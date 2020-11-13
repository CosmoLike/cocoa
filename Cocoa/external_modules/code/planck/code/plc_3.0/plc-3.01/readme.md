# plc 3.01 
August 2019

``plc`` is the public Planck Likelihood Code.  It provides C and Fortran 
libraries that allow users to compute the log likelihoods of the temperature, 
polarization, and lensing maps.  Optionally, it also provides a python version
of this library, as well as tools to modify the predetermined options for some 
likelihoods (e.g. changing the high-ell and low-ell lmin and lmax values of the
temperature).
 
Note that all of the old likelihood files (from plc 1 and plc 2) should still 
work correctely with the new version.

## Content
This release contains the code package `plc-3.0.tar.bz2`, as well as the low l 
(in `low_l/`), high l CMB likelihood (in `hi_l/`) and the lensing likelihood 
files (in `lensing/`). `cosmomc` files for the high l `plik` likelihood are 
available in the `cosmomc` directory.
All likelihood files need to be unzipped and untarred before use.
The Planck Legacy Archive distributes in separate packages the code and data 
files.

In more details:

- the low l (l=2 to 29) TT likelihood is provided by the file 
`low_l/commander/commander_dx12_v3_2_29.clik`

- the low l (l=2 to 29) EE likelihood is provided by 
`low_l/simall/simall_100x143_offlike5_EE_Aplanck_B.clik`

- For the high l (l=30 to 2508 in TT, 1997 in EE and TE), the joint TTTEEE 
likelihood is computed by the `hi_l/plik/plik_rd12_HM_v22b_TTTEEE.clik` file. 
(*Alternatively, one can use the TT only version of the high l likelihood 
provided by `hi_l/plik/plik_rd12_HM_v22_TT.clik`*). Note that the file 
`hi_l/plik/plik_recommended_priors.txt` contains the list of priors and settings
for the high l likelihood, translated in cosmomc form in the `cosmomc/` 
directory.

- The foreground marginalized likelihoods 
`hi_l/plik_lite/plik_lite_v22_TTTEEE.clik` and TT only 
`hi_l/plik_lite/plik_lite_v22_TT.clik` are also available.

- the lensing likelihood uses the file 
`lensing/smicadx12_Dec5_ftl_mv2_ndclpp_p_teb_consext8.clik_lensing` which uses 
a new conservative range of multipole between L=8 to 400. A CMB marginalized file,
`lensing/smicadx12_Dec5_ftl_mv2_ndclpp_p_teb_consext8_CMBmarged.clik_lensing` 
is also available. 

To compute the Planck TTTEEE+lowE+lensing likelihood, one should add the output 
from the two low-l files, one high-l and one lensing file. Planck baseline is to
use the  `hi_l/plik/plik_rd12_HM_v22b_TTTEEE.clik` high-l and 
`lensing/smicadx12_Dec5_ftl_mv2_ndclpp_p_teb_consext8.clik_lensing` lensing files.

Extra product are available (in dedicated files).

- At low-ell, the updated lowTEB likelihood (`bflike`) can be found in
`low_l/bflike/lowl_SMW_70_dx12_R1.8x_Ap.clik`

- At high-ell, the `hi_l/camspec` directory contains all the different camspec 
cases (TT, TE, EE, TTTE, TTEE, TEEE and TTTEEE), as well as 
`camspec_recommendations.txt`, the priors

- At high-ell, the `hi_l/plik_lite` directory also contains the TE and EE only 
files

- At high-ell, the `hi_l/plik` directory also contains EE and TE only 
likelihoods, as well as the unbinned cases for TT, TE, EE, and TTTEEE 
(`bin1` files)

- The agressive range (L<2000) files  
`lensing/smicadx12_Dec5_ftl_mv2_ndclpp_p_teb_agr2.clik_lensing` and CMB 
marginalized agressive range files  
`lensing/smicadx12_Dec5_ftl_mv2_ndclpp_p_teb_agr2_CMBmarged.clik_lensing` 
are also available.

All the extended cases are available in different PLA files.

## code 
The code is at version 3.01

###changes sinces v3.00

- Improve parsing of intel fortran dryrun output to fix an issue with the building of the intel fortran link line which was broken for version 2019
- Change an openmp scheduling option in src/plik/smica.c to go aroung an API incompatibility between intel 2019 openmp and gcc9 openmp which was preventing the code to correctely compile in this particular case.


## Installing

In the current version ``plc`` can only be installed using the waf configuration
and installation tool or using make and the associated `Makefile`. 
Note that in the latter case, the installation procedure will
not test the availability of the tools and libraries needed by plc.

### Prerequisites

**Mandatory**

- C compiler, either gcc, clang or icc. 
- fortran compiler, either ifort of gfortran
- ``blas/lapack`` library 
- ``cfitsio`` library 

**Optional**

- Python v>=2.7.8, including the header and libraries.
- numpy 
- pyfits or astropy
- cython v>=0.28


There are some incompatibilities between different versions of the C and fortran
compilers. In particular, gcc v4.9 and ifort v<14.1 are not compatible.

When installing manually the mandatory library, make sure that both the 
``blas/lapack`` and ``cfitsio`` are  built as shared libraries
(option ``make --shared`` of the ``cfitsio`` install).

The optional prerequisites are only needed for the optional ``plc`` tools 
(i.e. python tools). 
They are not required for basic uses of the library.

### installing with waf

Waf (http://waf.io) is a tool based on python. Installing with waf requires 
python v>2.7.8. The tool first need to be configured with

	waf configure [OPTIONS]

A complete list of the install options can be obtained by doing

	waf --help


The building and installing of the code is done using

	waf install


#### Simple configuration recipes

Detailed instructions on installing are found below. Here are a few simple and
classical cases.

##### MACOS

Use the following instructions to install on macos, letting the tool install all 
of the optional prerequisites and using the stock blas/lapack library. The C 
compiler will be either gcc, icc or clang (tested in that order). The fortran 
compiler will be eiher ifort or gfortran (tested in that order). Installation 
will be performed in the current plc directory.

	waf configure --install_all_deps 
	waf install


On macos, the installer will test whether the optional openmp clang support is 
installed. Please follow the isntruction at 
[https://iscinumpy.gitlab.io/post/omp-on-high-sierra/]()
to install openmp support on macos.

##### linux with mkl

Use the following instructions to install on a linux machine with ifort and mkl, 
letting the tool install all the optional prerequisites. The mkl library is 
available on the computer at the path ``$MKLROOT`` and its version is >=10.3 
(probably the case if your mkl library is post 2011). The C compiler will be 
either gcc or clang (tested in that order). The C compiler will be either gcc, 
icc or clang (tested in that order). The fortran compiler will be eiher ifort or
gfortran (tested in that order). 
Installation will be performed in the current plc directory.

	waf configure --install_all_deps --lapack_mkl=$MKLROOT
	waf install


#### Detailed configuration instructions

##### Install automatically the prerequisites


	waf configure --install_all_deps [OPTIONS]
 
asks the tool to *try* to install all of the absent prerequisites.

##### Changing the install path

By default, plc is installed in the source directory. To select a different
directory, use

	waf configure --prefix=/some/other/path [OPTIONS]


##### Selecting C and fortran compiler.

By default, the tool will check availability of gcc, icc and clang, in that 
order and select the first available compiler. To force gcc, icc or clang use one
of ``--gcc``, ``--icc`` or ``--clang``.

By default, the tool will check availability of ifort and gfortran, in that 
order and select the first available compiler. To force ifort or gfortran, use 
one of ``--ifort`` or ``--gfortran``.

##### Selecting a particular version of python 

By default, the library python wrapper and tools are compiled using the default python. To impose a particular version of python, run the configuration step by calling this particular version

	/path/to/some/pythonbin waf configure [OPTIONS]

and then compile and install wil

	/path/to/some/pythonbin waf install [OPTIONS]
  
##### ``blas/lapack``

By default, on macos, the stock ``blas/lapack`` will be used. On linux, the waf 
tool  needs to be pointed to a particular library
To use an alternative library on macos, or to point a particular library on 
linux one can use the following options.

###### ``mkl``

The simplest option is to use the intel ``mkl`` library. Assuming that the library
is installed at ``$MKLROOT`` and that the version of the library is 10.3 or 
better (which should be the case is the library was installed post 2011), use

	waf configure --lapack_mkl=$MKLROOT [OPTIONS]

To use an earlier version of the library (10.0, 10.1 or 10.2), for example 
version 10.2 use

	waf configure --lapack_mkl=$MKLROOT --lapack_mkl_version=10.2 [OPTIONS]


##### Other ``blas/lapack``

One can provide the location of the blas/lapack install using either
- ``--lapack_prefix=/some/path``, the include and library path will be 
``/some/path/lib`` and ``/some/path/include``, and the installer will try to 
link to ``libblas.so`` and ``liblapack.so``
- ``--lapack_lib=/some/path/lib --lapack_include=/some/path/include``, allow to 
define the library and include path. The installer will try to 
link to ``libblas.so`` and ``liblapack.so``
-- ``--lapack_link='-I/some/path/include -L/some/path/lib -lmyblas -lmylapack',
allow to define the full link line for the lapack/blas install.

##### Installing a basic blas/lapack automatically

Using either ```--install_all_deps``` or ```--install_lapack``` will ask the 
tool to install a plain vanilla, non-optimized blas/lapack library.


#### ``cfitsio``

``plc`` requires a ``cftisio`` library compiled as a shared library (option 
``--shared`` when building cfitsio). The tool will look for the include and 
library in the usual locations. If it fails, one can point it to a particular 
version of the ``cfitsio`` library using either
- ``--cfitsio_prefix=/some/path``, the include and library will be 
``/some/path/include`` and ``/some/path/lib``
- ``--cfitsio_lib=/some/path/lib --cfitsio_include=/some/path/include``, allow to 
define the library and include path. 

Alternatively, the ``cfitsio`` library can be installed by the tool using either
```--install_all_deps``` or ```--install_cfitsio```.

### Installing with make

All prerequisites must be installed before running make.  The ``Makefile`` file 
must be adapted to the particular computer configuration.  The C and Fortran
compiler must be selected, and in the case of the ifort compiler also the list
of runtime libraries (needed to link with C).  Please look at the contents of
the file for details.

Once the ``Makefile`` file is correct, build and install with 
```
make install
```
The optional python tools can be built using
```
make python_install
```

## Testing that the code is correctly built and that the data files are correct

First source the ``clik_profile.[c]sh`` file to import the environment variables.
Then, running ``clik_example_C`` on any likelihood file should perform an
automatic test.
For example

  	$> clik_example_C plik_rd12_HM_v22b_TTTEEE.clik/
  	----
  	clik version 5fde30f43612
      smica
  	Checking likelihood '../../src/baseline/plc_3.0/hi_l/plik/plik_rd12_HM_v22b_TTTEEE.clik/' on test data. got -1172.47 expected -1172.47 (diff -4.34053e-07)
  	----
  	...


## Basic use

The base installation of plc consists in a library (libclik.so) along with a C 
include and a fortran module, which allows to
- initialize a likelihood
- ask a likelihood for its characteristics (list of spectra, lmin, lmax and 
list of nuisance parameters)
- compute a log likelihood
- deallocate a likelihood

### Setting up the environment

To set the environment path and variables, source the ``clik_profile.sh`` (or 
``clik_profile.csh``) file. This can be added to your ``.login`` or other 
shell initialization script.
Once the script is executed, the variable ``$CLIK_PATH`` will point toward the
root directory of your plc install.

### Compiling with plc

The compilation and link options to compile a C code against the clik library
can be obtained using the script ``clik-config``.
Similarly, ``clik-config_f90`` provides the same information when compiling
a fortran code.

### Clik file

The data for the library are stored in directories which will be named 
*clik files* in the following. In the directories, one can find ascii metadata
files (``_mdb``) and data in multiple binary format (fits for numerical data).

### Using the C Library

Include ``clik.h`` to use the plc library. 

#### Initialization

The function 
```
clik_object* clik_init(char* filepath, error **err);
``` 
initializes a CMB likelihood object of type ``clik_object`` from a clik file at 
location ``filepath``. The error variable can be set to ``NULL``, in which case
the function will exit in case of error. Details on the optional error reporting
can be found in ``errorlist.h``.
In most cases, more than one likelihood can be initialized in a code. 

Similarly,
```
clik_object* clik_lensing_init(char* filepath, error **err);
``` 
initializes a lensing likelihood.

During the initialization, an automatic test is performed and displayed on stdout.

#### Obtaining the list of spectra, lrange, nuisance parameters

##### CMB likelihood

The function
```
void clik_get_lmax(clik_object *clikid, int lmax[6],error **err);
```
fills a array of 6 integer ``lmax``. Each element of the array correspond to the
lmax of one of the CMB spectra needed by the CMB likelihood referred by ``clikid``.
Ordering is **TT EE BB TE TB EB**. The array elements are set to -1 for each spectra
unused by the current likelihood.
As above, error can be set to ``NULL`` to ignore error management.
For example, when used with the high-ell temperature likelihood, ``lmax`` will be 
set to `` 2508 -1 -1 -1 -1 -1``, meaning that only the TT spectrum is needed by 
the likelihood, up to lmax = 2508.

The function 
```
int clik_get_extra_parameter_names(clik_object* clikid, parname **names, error **err);
``` 
returns the number of nuisance parameters. If the ``names`` argument is not ``NULL``
it will fill it with a pointer to a list of names of the nuisance parameters.
After usage, the memory pointed by ``names`` must be freed by the user.
As above, error can be set to ``NULL`` to ignore error management.

### Lensing likelihood

The function 
```
int clik_lensing_get_lmaxs(clik_lensing_object *lclik, int lmax[7], error **err);
```
works like the CMB function ```clik_get_lmax```.
Note that the lmax arrays is longer than in the CMB case. Order in this case is
**phiphi TT EE BB TE TB EB**. In the future, this array can be extended to 
include phi-CMB correlations. 

The function 
```
int clik_lensing_get_extra_parameter_names(clik_object* clikid, parname **names, error **err);
``` 
works like the CMB one.

#### Computing a log likelihood

The function
```
double clik_compute(clik_object* clikid, double* cl_and_pars,error **err);
```
computes the log likelihood for the CMB likelihood pointed by ``clikid`` and for
the cl and nuisance parameters contained in the array ``cl_and_pars``. As above, 
error can be set to ``NULL`` to ignore error management.
The ``cl_and_pars`` array is an array of double that contains the CMB cls and 
nuisance parameter. The order in the array is as follows :

- first the Cls, from l=0 to lmax (included) in the order **TT EE BB TE TB EB**. 
Only the spectra whose ``lmax`` (from the function ``clik_get_lmax``) is >-1 are
present. 

- then the nuisance parameters in the order defined by 
``clik_get_extra_parameter_names``.

For example, for the low ell TEB likelihood whose lmax = ``29 29 29 29 -1 -1`` 
and which has a unique nuisance parameter, the planck absolute calibration, the
``cl_and_pars`` array should be
```
c0_TT
c1_TT
c2_TT
...
c28_TT
c29_TT
c0_EE
c1_EE
...
c29_EE
c0_BB
...
c29_BB
c0_TE
...
c29_TE
planck_calibration
```

**Please note that the function expects Cl and not l(l+1)Cl/2pi.**
**Also note that the function returns a log likelihood and not a chi2.**

The function
```
double clik_lensing_compute(clik_object* clikid, double* cl_and_pars,error **err);
```
provides the log likelihood for a lensing, and works similarly to the CMB 
function. The only difference is in the ordering of the ``cl_and_pars`` array, 
which must include the **phiphi** spectrum, before the CMB spectra. Note that
CMB spectra are needed to compute the lensing normalization and biases.

#### Cleaning up

A CMB likelihood can be deallocated using 
```
void clik_cleanup(clik_object** pclikid);
``` 
and a lensing one with 
```
void clik_lensing_cleanup(clik_lensing_object **plclik);
```.

#### Example code

The file in ``src/clik_example_C.c`` demonstrates the use of the C library.
This code compiles to an executable which reads the path of a clik file as its
first argument and optionally Cl+nuisance ascii files as it next ones.
The likelihood is initialized, some information (lrange, nuisance parameters) is
displayed on screen. Then each optional Cl+nuisance ascii file is read and the 
log likelihood for each one is computed and displayed.

The Cl+nuisance files must be ascii files, containing one value on each line, 
forming a vector of values for the Cl (not Dl!) and nuisance parameters in the 
order expected for the compute functions. 


### Using the fortran Library

Use the ```clik``` module in your code.

#### Initialization

The subroutine 
```
subroutine clik_init(clikid,filepath)
``` 
initializes a CMB likelihood with the data at ``filepath`` and fills ``clikid`` 
with an integer referring to the likelihood. 

Similarily 
```
subroutine clik_lensing_init(clikid,filepath)
``` 
initializes a lensing likelihood.

During the initialization, an automatic test is performed and displayed on stdout.

#### Obtaining the list of spectra, lrange, nuisance parameters

##### CMB likelihood

The subroutine
```
subroutine clik_get_lmax(clikid,lmax)
```
fills a array of 6 integer ``lmax``. Each element of the array correspond to the
lmax of one of the CMB spectra needed by the CMB likelihood refered by ``clikid``.
Ordering is **TT EE BB TE TB EB**. The array element are set to -1 for each spectra
unused by the current likelihood.
For example, when used with the high-ell temperature likelihood, ``lmax`` will be 
set to `` 2508 -1 -1 -1 -1 -1``, meaning that only the TT spectrum is needed by 
the likelihood, up to lmax = 2508.

The function 
```
integer(kind=4) function clik_get_extra_parameter_names(clikid,names)
``` 
returns the number of nuisance parameters. ``names`` is filled with a list of 
names of the nuisance parameters.
After usage, the memory pointed by ``names`` must be freed by the user.

### Lensing likelihood

The subroutine
```
subroutine clik_lensing_get_lmax(clikid,lmax)
```
works like the CMB subroutine ```clik_get_lmax```.
Note that the lmax arrays is longer than in the CMB case. Order in this case is
**phiphi TT EE BB TE TB EB**. In the future, this array can be extended to 
include phi-CMB correlations. 

The function 
```
integer(kind=4) function clik_lensing_get_extra_parameter_names(clikid,names)
``` 
works like the CMB one.

#### Computing a log likelihood

The function
```
real(kind=8) function clik_compute(clikid,cl_and_pars)
```
computes the log likelihood for the CMB likelihood refered by ``clikid`` and for
the cl and nuisance parameters contained in the array ``cl_and_pars``.
The ``cl_and_pars`` array is an array of double that contains the CMB cls and 
nuisance parameter. The order in the array is as follows :

- first the Cls, from l=0 to lmax (included) in the order **TT EE BB TE TB EB**. 
Only the spectra whose ``lmax`` (from the function ``clik_get_lmax``) is >-1 are
present. 

- then the nuisance parameters in the order defined by 
``clik_get_extra_parameter_names``.

For example, for the low ell TEB likelihood whose lmax = ``29 29 29 29 -1 -1`` 
and which has a unique nuisance parameter, the planck absolute calibration, the
``cl_and_pars`` array should be
```
c0_TT
c1_TT
c2_TT
...
c28_TT
c29_TT
c0_EE
c1_EE
...
c29_EE
c0_BB
...
c29_BB
c0_TE
...
c29_TE
planck_calibration
```

**Please note that the function expects Cl and not l(l+1)Cl/2pi.**
**Also note that the function returns a log likelihood and not a chi2.**

The function
```
real(kind=8) function clik_lensing_compute(clikid,cl_and_pars)
```
provides the log likelihood for a lensing map, and works similarly to the CMB 
function. The only difference is in the ordering of the ``cl_and_pars`` array, 
which must include the **phiphi** spectrum, before the CMB spectra. Note that
CMB spectra are needed to compute the lensing normalization and biases.

#### Cleaning up

A CMB likelihood can be deallocated using 
```
subroutine clik_cleanup(clikid)
``` 
and a lensing one with 
```
subroutine clik_lensing_cleanup(clikid)
```

#### Example code

The file in ``src/clik_example_F90.f90`` demonstrates the use of the C library.
This code compiles to an executable which reads the path of a clik file as its
first argument and optionally Cl+nuisance ascii files as it next ones.
The likelihood is initialized, some information (lrange, nuisance parameters) is
displayed on screen. Then each optional Cl+nuisance ascii file is read and the 
log likelihood for each one is computed and displayed.

The Cl+nuisance files must be ascii files, containing one value on each line, 
forming a vector of values for the Cl (not Dl!) and nuisance parameters in the 
order expected for the compute functions. 

### Test tools

The executable 
```
clik_example_C /some/path/to/likelihood [Cl+nuisance.txt [Cl+nuisance.file [...]]]
``` 
and 
```
clik_example_F90 /some/path/to/likelihood [Cl+nuisance.txt [Cl+nuisance.file [...]]]
``` 
performs tests on a given likelihood file, display some information and allow to
compute the log likelihood of Cl+nuisance vectors. 
The Cl+nuisance files must be ascii files, containing one value on each line, 
forming a vector of values for the Cl (not Dl!) and nuisance parameters in the 
order expected for the compute functions. 

## Optional tools

The optional tools are only available if the optional requirements are met.

### Using the python library

The library can be called from python by importing the ``clik`` python package.

#### Initialization

CMB Likelihoods are represented by instance of the ``clik`` objects. 
Lensing likelihoods by instance of ``clik_lensing`` objects.
Both are initialized from a clik file.

```
import clik
CMBlkl = clik.clik("/some/path/to/clikfile")
lenslkl = clik.clik_lensing("/some/path/to/cliklensingfile")
```

#### Obtaining the list of spectra, lrange, nuisance parameters

##### CMB likelihood

The ``get_lmax`` method of the ``clik`` object returns a 6 element tuple 
containing the maximum multipole for each of the CMB spectra. Please see the C
description above for explanation.

The ``get_extra_parameter_names`` method of the ``clik`` object returns a tuple
containing the names of the nuisance parameters

##### Lensing likelihood

The ``get_lmax`` method of the ``clik_lensing`` object returns a 6 element tuple 
containing the maximum multipole for each of the phi and CMB spectra. Please see 
the C description above for explanation.

The ``get_extra_parameter_names`` method of the ``clik_lensing`` object returns 
a tuple containing the names of the nuisance parameters

The ``get_clpp_fid`` method of the ``clik_lensing`` object returns the fiducial
lensing spectrum, used to perform the automatic test, and to compute the N0 and 
N1 biases. Correction to those biases are computed as perturbations around this
spectrum. It returns an array which contains the phiphi lensing spectra from 0 
up to an lmax (included) given by the first element of the tuple return by the 
``get_lmax``.

The ``get_clcmb_fid`` method of the ``clik_lensing`` object returns the fiducial
CMB spectra, used to perform the automatic test, and to compute the N0 and 
N1 biases. Correction to those biases are computed as perturbations around those
spectra. it returns an array containing all of the CMB spectra from 0 to the 
lmax given by the values at index >=1 in the tuple returned by the ``get_lmax`` 
method.

#### Computing a log likelihood

The ``clik`` and ``clik_lensing`` object are callable. They expect a single 
argument cl_and_pars which must be either a one dimensional or a two dimensional 
array of float. The shape of the array must be either ``(ntot,)`` or ``(i,ntot)``.
If the array is one dimensional, it will be reshaped internally to ``(1,ntot)``.
``ntot`` is the number of element of the vector expected by the
``compute`` function described in the C and F90 version of the library.  The
function will compute the log likelihood for each vector of spectra and nuisance
parameters and return an array of log likelihood values. 

#### Error handling

In the case where the library encounter an error, a ``clik.lkl.CError`` error
object will be raised.

### Example code

The script ``clik_example_py`` demonstrate how to use the python library.
Its usage is similar to the ``clik_example_C`` and ``clik_exmaple_f90`` tools
described above

### Displaying some information with ``clik_print``

The script ``clik_print`` display on the screen some more information than the
``clik_example_XXX`` tools. It can also be used to check the correct building
and installation of the library.
It expects a single argument on its command line giving the path to a likelihood
file.
```
clik_print /where/is/my/clikfile
```

### Retrieving the selfcheck vector

The script `` clik_get_selfchek`` allows to retrieve the data used in the 
automatic test performed at each initialization. The scripts expects 2 arguments 
on his command line, the path to a clik file, and the name of a file in which 
the cl and parameters vector will be saved as an ascii file, in a format suitable
for the ``clik_example_XXX`` tools. Indeed, calling any of the ``clik_example_XXX``
with the same likelihood and using the data vector obtained this way should 
produce the same log likelihood value displayed during the automatic test.
```
clik_get_selfcheck /some/clikfile /where/to/save/datavector
```

### Changing the lmin and lmax of some of the likelihoods

As in version 2.0, the ``clik_change_lrange`` allows to change the lmin and lmax of the plik and 
commander likelihoods. This allow to reproduce the tests performed on the 
hybridization multipole, and on the high ell lmax described in the paper.
```
clik_change_lrange input_clik lmin lmax output_clik
```
The script expects an input clik file, the new lmin and lmax and an output clik 
file. The input file will not be modified. A value of -1 for lmin (resp. lmax)
asks the script not to modify the value.
The lmin and lmax are targets. The output file will have lmin and lmax as close 
as possible to those value, given the available data and the binning scheme 
(for plik). They will be displayed.
The script will fail if the path ``output_clik`` already exists.
This tool will only work on commander and TT plik files.

A new (more complex) tool, ``clik_change_lrange_bychannel_pol`` allows to cut in the plik likelihood, frequency by frequency, in Temperature and Polarization. The script expects a parameter files which provides the name of the input likelihood, the name of the output file and a matrix of lmin and lmax to define the new lmin and lmax for each cross spectrum. Here follows an example in the case a TTTEEE plik likelihood is reduced to the TEEE only part for the 143x143 and the 100x217, ensuring that the new lmin for each crosspectrum NUxMU that is kept is max(lminNUxMU, 800) and new lmax is min(lmaxNUxMU,1500)

    input_clik = some_TTTEEE.clik

    output_clik = cut_TEEE.clik

    #set to -1 to use default value
    #set to anything else to change lmin
    #      T100 T143 T217 E100 E143 E217
    lmin =   000  000  000  000  000  800 & 
             000  000  000  000  800  000 &
             000  000  000  800  000  000 &
             000  000  800  000  000  800 &
             000  800  000  000  800  000 &
             800  000  000  800  000  000 
 
    #set to -1 to use default value
    #set to anything else to change lmax
    #set any non negative value inferior than lmin 
    # (0 is a good example) to remove a cross spectra

    #      T100 T143 T217 E100 E143 E217
    lmax =   000  000  000  000  000 1500 & 
             000  000  000  000 1500  000 &
             000  000  000 1500  000  000 &
             000  000 1500  000  000 1500 &
             000 1500  000  000 1500  000 &
            1500  000  000 1500  000  000 


### Exploring the content of clik files

The scripts ``clik_cldf_ls`` and ``clik_cldf_dump`` allow to explore the content
of a clik file. ``clik_cldf_ls`` allow to browse the data tree of 
the file, while ``clik_cldf_dump`` will display the content of a given leaf. 
Please note that when the leaf contain long arrays, only a subset of the first
and last part of the vector will be displayed.
The file can also be browsed in python using the ``clik.cldf`` module.


    import clik.cldf

    # open file
    clf = clik.cldf.File("plik_dx11dr2_HM_v18_TT.clik")

    # get the content of the root
    print clf.keys()

    # is the path at clik a subtree or a leaf with data ?
    print isinstance(clf["clik"],clik.cldf.File)

    # get the list of content at path clik/lkl_0
    print clf["clik/lkl_0"].keys()

    # get the number of T channels in the likelihood
    print  clf["clik/lkl_0/m_channel_T"]










