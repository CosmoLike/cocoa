Installing
==========

The package can be installed using two different tools, `waf <http://waf.googlecode.com>`_ or `make <http://http://www.gnu.org/software/make/>`_.
Using waf, the installer will test for the different dependencies and try to install them if they are missing. Using make, one would have to
modify the ``Makefile`` file for your particular computer and install the dependencies. Besides, when using make, no test of the availability of the requisite will be made.

The package has a set of core utilities, consisting in the ``clik`` library with a C and F90 API, as well as the programs :program:`clik_example_C` and :program:`clik_example_F90` which allows to compute a log likelihood for a Cl and nuisance parameter file, and also doubles as example of how to interface with the library in this two languages.

The package also have a set of optional utilities, consisting in a python wrapper of the library, along with a few python scripts allowing to explore the content of likelihood files, 
and to manipulate them. Those tools will have a few more requirements to be built.

Finally, when using waf, and provided that the optional python utilities can be installed, a wrapper to the wmap9 likelihood will also be installed.

Requisites
----------

Mandatory requisites 
^^^^^^^^^^^^^^^^^^^^

**modern c and fortran compilers** (icc, gcc, ifort anf gfortran are ok) are absolute requisites. 
Gcc version must be >= 4.2 and gfortran version must be >=4.3.

To use the waf tool, one also need `Python <http://python.org>`_ (>=2.5).

Mandatory requisites that can be installed when using waf
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Those requisites will not be installed automatically when using make, and must be manually installed. 
`cfitsio <http://heasarc.gsfc.nasa.gov/fitsio/>`_,  as well as blas and lapack distribution (preferably intel MKL on linux machine. Macos computers already have a reasonnable parallelized blas/lapack that clik can use) are needed for the core functionalities of clik. If absent (or not available in the correct flavor), they can be installed automatically using waf and the options described below.
Beware that the shared version of the cfitsio library must be available. This can be obtained using ``make shared`` when building the cfisio library.



Optional requisites 
^^^^^^^^^^^^^^^^^^^

A python distribution including the header and library (if absent, this will not be installed automatically by waf), and the `pyfits <http://http://www.stsci.edu/institute/software_hardware/pyfits/>`_ (2.4<=version),  `numpy <http://numpy.scipy.org/>`_ (version>1.1) and `cython <http://cython.org/>`_ python package (version>1.12) are needed to provide the (optional) clik python wrapper and tools. The three python package will only be installed automatically when using waf if the python header and library are available on the system.

Install with waf
----------------

The library and executables must be installed with the 'waf' tool. It is distributed in the package. Please have a look at `the waf webpage <http://waf.googlecode.com>`_.

The installer must first be configured using::

    $> ./waf configure

This will test for the presence of all the required dependencies (as described above). Command line options are available to help the configuration by setting the location of some of those dependency, and or chose between different compiler options. An option can also be used to cause waf to try to download, compile and install for the user the required dependencies that are not found by the automatic discovery system. For a complete list of options do::

    $> ./waf configure --help


After this (possibly lengthy) configuration step, clik per-se can be compiled and installed with::

    $> ./waf install

Simplest case
^^^^^^^^^^^^^

In the simplest case, we will assume that you want all the dependencies absent from the usual locations to be installed automatically. Note that this translate into a rather slow clik library, since the lapack library will probably be compiled from a simple, non-parallel, version of the lib. That being said , the simplest configuration (if slow) line would be::

    $> ./waf configure --install_all_deps

followed by::

    $> ./waf install

Note that the option ``--install_all_deps`` is more powerful than simply installing all the dependency, this will be described in a latter section.

Simplest case with mkl
^^^^^^^^^^^^^^^^^^^^^^

Using the automatically installed lapack library results in a slow clik library. One can solve easily this problem by letting the configuration know about an existing mkl library. To do so the configuration line has to be changed into::

    $> ./waf configure --install_all_deps \
          --lapack_mkl=/PATH/TO/YOUR/MKL/ --lapack_mkl_version=MKL.VERSION

The ``/PATH/TO/YOUR/MKL/`` can be replaced by ``${MKLROOT}`` for most recent mkl installs. Otherwise, it has to be a directory containing a ``lib`` and an ``include`` subdirectories. The ``MKL.VERSION`` can be any of 10.0, 10.1, 10.2, 10.3.

Simplest case with apple lapack
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

MacOS X has a standard parallelized lapack distribution. One can use it to accelerate clik.
To do so the configuration line has to be changed into::

    $> ./waf configure --install_all_deps --lapack_apple

(planck insider) I am installing on magique3
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Use this::

    $> /softs/python/2.7.2/bin/python waf configure \
          --install_all_deps --cfitsio_prefix=/softs/cfitsio/3.24/
          --lapack_mkl=/softs/intel/mkl/10.2.6.038 --lapack_mkl_version=10.2 

and then::

    $> ./waf install

See below for the other otions used here.

(planck insider) I am installing on ccin2p3
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Use this::

    $> ./waf configure --install_all_deps \
             --lapack_mkl=/usr/local/intel/mkl/10.3.8/ --lapack_mkl_version=10.3

and then::

    $> ./waf install




ADVANCED: Installing with a particular Python executable
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It is possible to install clik with a python install different from the default one. For example if the default python installation does not contains the required header and libraries. To do so, call waf this way::

    $> /path/to/special/python waf configure 

and then::

    $> /path/to/special/python waf install 


ADVANCED: Bypassing the default compilers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To bypass the c compiler detection, set the ``CC`` environment variable. 
To bypass the fortran compiler detection, set the ``FC`` environment variable. Beware, you can only set the ``FC`` environment variable to either an intel fortran compiler or a gfortran compiler. 

Shortcuts for some classical cases are provided:

    * ``--icc`` causes the installer to use icc as c compiler.
    * ``--ifort`` causes the installer to use ifort as fortran compiler.
    * ``--gcc`` causes the installer to use gcc as c compiler.
    * ``--gfortran`` causes the installer to use gfortran as fortran compiler.


ADVANCED: Setting the architecture
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The architecture (32 or 64bits) can be set using the ``--m32`` or ``--m64`` flags. 64bits is the default.

ADVANCED: Setting installation path
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The installation path can be set using the ``--prefix=SOMEPATH`` option. Default is to install in the current directory.


ADVANCED: More on the automatic installation of dependencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

There are three levels of automatic installation. If one wants to *always* install the dependencies, one can use the ``--force_install_all_deps``::

    $> ./waf configure --forceinstall_all_deps

If one wants to install only the dependencies that are not present in the usual location (or that are present, but not compiled in a way suitable for clik), one can use the ``--install_all_deps`` option, already described above. Since this option first tests for the presence of each library, it can be used to upgrade a clik install, avoiding to reinstall everything.

Finally, each dependency can be installed on a dependency by dependency basis, using the ``--XXX_install`` or ``--XXX_installifneeded`` options where ``XXX`` is the name of the dependency. The former install all the time the dependency, the latter install it only if it is not found in the usual locations. In that sense, ``--forceinstall_all_deps`` works as if all possible ``--XXX_install`` options has been set, and ``--install_all_deps`` as if all ``--XXX_installifneeded`` options have been set.

One should also note that ``--forceinstall_all_deps`` and ``--install_all_deps`` are also unactivated on a dependency by dependency basis if any of the ``--XXX_prefix``, ``--XXX_lib``, ``--XXX_include``, or other dependency specific options are present. In that case, the the ``XXX`` dependency, the configuration script will look in the locations described by those option and if the package is not found will report an error.


ADVANCED: Setting the location of a library
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The location of the library dependencies (gsl, hdf5, healpix, blas/lapack) must be known to the installer. By default, it will look for them in the classical system 
locations:  ``/usr/lib``, ``/usr/lib64``, ``/usr/local/lib``, ``/usr/local/lib64`` for the library, ``/usr/include`` and ``/usr/local/include`` for the include files. One can 
change the lookup path on a library by library basis. If a given dependency, ``XXX``, is installed on the system such that its lib are in ``SOMEPREFIXPATH/lib`` and its 
include files in ``SOMEPREFIXPATH/include``, setting the command line option ``--XXX_prefix=SOMEPREFIXPATH``  will allow the clik install system. If ``SOMEPREFIXPATH`` is identical to the the install path of clik, this option can be replaced by ``-XXX_islocal``.

If the library are at 
``SOMEWEIRDPATH`` and the includes at ``SOMEDIFFERENTPATH``, then setting the two options  ``--XXX_lib=SOMEWEIRDPATH --XXX_include=SOMEDIFFERENTPATH`` will allow the clik 
install system to find them.

Finally, if the name of the library files differs from the usual ones one can set the option ``--XXX_link=THELINKLINE``.

Using these options allow to point the installer to a pmclib install in order to allow the linking of clik with pmclib.


ADVANCED: Special case: the mkl library
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This option is only for advanced users.
The blas/lapack distribution installed automatically is a very inefficient one. To improve the performance of clik (especially the low-l pixel based likelihood), one is advised to use the MKL library, which is fully supported and allow the use of shared memory computer architectures.

A special option is present to simplify the install using the intel MKL library: setting the option ``--lapack_mkl=PATH_OF_THE_MKL_INSTALL`` together with ``--lapack_mkl_version=SOMEVERSION`` will allow clik to pick the correct set of libraries for the particular version of the mkl package (version 10.0, 10.1, 10.2 and 10.3 only).
Setting this option will cancel the ``--install_all_deps`` option for the lapack dependency only.

On a MacOS X computer, one can use Apple provided lapack by setting ``--lapack_apple``.


ADVANCED: Special case: WMAP likelihood
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Clik can provide a wrapper to the wmap9 likelihood. It need to now where the sources of the likelihood are located to compile against them. One must set the option ``--wmap_src=WMAP7SRCPATH`` or let the install system download it for you by setting the option ``--wmap_install``. Note that to actually use this likelihood, one must also download the data files and prepare clik likelihood files from them. Look at :ref:`WMAP`. The ``--install_all_deps`` and ``--forceinstall_all_deps`` options will automatically download the sources, as if ``-wmap_install`` was set.


ADVANCED: Putting it all together
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following command::

    $> ./waf configure --install_all_deps

will tell the clik install system to install all the possible external dependency in the current directory. 

The following command::

    $> ./waf configure --lapack_mkl=/opt/intel/mkl \ --lapack_mkl_version=10.2
       --cfistio_prefix=/usr/local/cfitsio --cython_install

will tell the clik install system to install cython. The cfitsio library will be looked for in the unusual dir ``/usr/local/cfitsio``. /All the other dependency will be looked up in the classical locations. The blas/lapack library 
will be the one from an mkl install located at --lapack_mkl=/opt/intel/mkl. Clik will be compiled in 64bit and installed in the current directory.

 
ADVANCED: Best advanced choice 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Use a mkl lapack install and let the other dependencies on auto install::

    $> ./waf configure --install_all_deps  \
          --lapack_mkl=/opt/intel/mkl --lapack_mkl_version=10.2 

This will use your mkl libraries from ``/opt/intel/mkl``, test if numpy, cython and gsl are installed on your computer (often the case) if not install them, 
and finally install all the other requirements (helpaix, hdf5 and its python wrapper).

Installing with make
--------------------

The first lines of the ``Makfile`` file must be checked and modified before compiling.
In particular, one must set, 
 * the location of the cfitsio library
 * the location and version of the mkl library (or other lapack library)
 * the location and list of library needed to link c with fortran
 
To build and install the core utilities, use the following command::

  $> make install

To build and install the optional utilities, use the following command::

  $> make install_python

Environment variables
---------------------

Depending of your shell, a configuration file named ``clik_profile.sh`` of ``clik_profile.csh`` will be installed in the ``bin`` directory at the install location of clik. One can source it on the command line, or include it in its startup configuration file to set the environment variable needed by clik. This tool is installed both bty waf and make.


