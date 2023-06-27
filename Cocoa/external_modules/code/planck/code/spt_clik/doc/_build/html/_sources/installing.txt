Installing
==========

The library has many external dependencies. We are striving at reducing the number of dependency, but many of them will remain given the current set of likelihood codes we are including.

Only a few packages or tools must be installed before installing clik. All the other can be downloaded and installed automatically by the clik installer.

Furthermore some of the dependency can be absent, and will only reduce the number of extra functionalities provided by clik (like merging likelihood, or simulating them).

Requisites
----------

Mandatory requisites that cannot be installed automatically
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

`Python <http://python.org>`_ (>=2.5) to run the installer (waf), **modern c and fortran compilers** (icc, gcc, ifort anf gfortran are ok) are absolute requisites. 
Gcc version must be >= 4.2 and gfortran version must be >=4.3.
Having a full python installation (i.e. including the library and header, often called *python-devel* or *python-dev* in package managers) is very strongly advised.

Mandatory requisites that will be installed automatically if absent
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The `hdf5 <http://www.hdfgroup.org/HDF5>`_ (version>1.8), `cfitsio`_, `healpix c and f90 <http://healpix.jpl.nasa.gov/>`_ (version>2.20a shared) libraries with their dependency (i.e. cfitsio) as well as blas and lapack distribution (preferably intel MKL) are needed for the core functionalities of clik. If absent (or not available in the correct flavor), they can be installed automatically using the options described below.

Optional requisites 
^^^^^^^^^^^^^^^^^^^

A python distribution including the header and library (if absent, this will not be installed automatically), and the `h5py <http://alfven.org/wp/hdf5-for-python/>`_ (1.3<=version<2),  `numpy <http://numpy.scipy.org/>`_ (version>1.1) and `cython <http://cython.org/>`_ python package (version>1.12) are needed to provide the (optional) clik python wrapper and the simulation, merging , splitting and printing tools. The three python package will only be installed automatically if the python header and library are available on the system.

Special dependency : pmclib
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Clik uses some of the facilities available in pmclib. Instead of requiring pmclib, those facilities are included in the clik distribution. This prevents a stock clik install to be linked with pmclib. An option is available to point the installer to a pmclib install, resolving this issue.


Install tool
------------

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

I am installing on magique3
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Use this::

    $> /softs/python/2.7.2/bin/python waf configure \
          --install_all_deps --gsl_prefix=/softs/gsl/1.15/ \
          --lapack_mkl=/softs/intel/mkl/10.2.6.038 --lapack_mkl_version=10.2 

and then::

    $> ./waf install

See below for the other otions used here.

I am installing on ccin2p3
^^^^^^^^^^^^^^^^^^^^^^^^^^

Use this::

    $> ./waf configure --install_all_deps \
             --lapack_mkl=/usr/local/intel/mkl/10.3.8/ --lapack_mkl_version=10.3

and then::

    $> ./waf install



Advanced configuration options
------------------------------

Installing with a particular Python executable
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It is possible to install clik with a python install different from the default one. For example if the default python installation does not contains the required header and libraries. To do so, call waf this way::

    $> /path/to/special/python waf configure 

and then::

    $> /path/to/special/python waf install 


Bypassing the default compilers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To bypass the c compiler detection, set the ``CC`` environment variable. 
To bypass the fortran compiler detection, set the ``FC`` environment variable. Beware, you can only set the ``FC`` environment variable to either an intel fortran compiler or a gfortran compiler. 

Shortcuts for some classical cases are provided:

    * ``--icc`` causes the installer to use icc as c compiler.
    * ``--ifort`` causes the installer to use ifort as fortran compiler.
    * ``--gcc`` causes the installer to use gcc as c compiler.
    * ``--gfortran`` causes the installer to use gfortran as fortran compiler.


Setting the architecture
^^^^^^^^^^^^^^^^^^^^^^^^

The architecture (32 or 64bits) can be set using the ``--m32`` or ``--m64`` flags. 64bits is the default.

Setting installation path
^^^^^^^^^^^^^^^^^^^^^^^^^

The installation path can be set using the ``--prefix=SOMEPATH`` option. Default is to install in the current directory.


More on the automatic installation of dependencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

There are three levels of automatic installation. If one wants to *always* install the dependencies, one can use the ``--force_install_all_deps``::

    $> ./waf configure --forceinstall_all_deps

If one wants to install only the dependencies that are not present in the usual location (or that are present, but not compiled in a way suitable for clik), one can use the ``--install_all_deps`` option, already described above. Since this option first tests for the presence of each library, it can be used to upgrade a clik install, avoiding to reinstall everything.

Finally, each dependency can be installed on a dependency by dependency basis, using the ``--XXX_install`` or ``--XXX_installifneeded`` options where ``XXX`` is the name of the dependency. The former install all the time the dependency, the latter install it only if it is not found in the usual locations. In that sense, ``--forceinstall_all_deps`` works as if all possible ``--XXX_install`` options has been set, and ``--install_all_deps`` as if all ``--XXX_installifneeded`` options have been set.

One should also note that ``--forceinstall_all_deps`` and ``--install_all_deps`` are also unactivated on a dependency by dependency basis if any of the ``--XXX_prefix``, ``--XXX_lib``, ``--XXX_include``, or other dependency specific options are present. In that case, the the ``XXX`` dependency, the configuration script will look in the locations described by those option and if the package is not found will report an error.


Setting the location of a library
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


Special case: the mkl library
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This option is only for advanced users.
The blas/lapack distribution installed automatically is a very inefficient one. To improve the performance of clik (especially the low-l pixel based likelihood), one is advised to use the MKL library, which is fully supported and allow the use of shared memory computer architectures.

A special option is present to simplify the install using the intel MKL library: setting the option ``--lapack_mkl=PATH_OF_THE_MKL_INSTALL`` together with ``--lapack_mkl_version=SOMEVERSION`` will allow clik to pick the correct set of libraries for the particular version of the mkl package (version 10.0, 10.1, 10.2 and 10.3 only).
Setting this option will cancel the ``--install_all_deps`` option for the lapack dependency only.

On a MacOS X computer, one can use Apple provided lapack by setting ``--lapack_apple``.


Special case: WMAP likelihood
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Clik can provide a wrapper to the wmap7 likelihood. It need to now where the sources of the likelihood are located to compile against them. One must set the option ``--wmap_src=WMAP7SRCPATH`` or let the install system download it for you by setting the option ``--wmap_install``. Note that to actually use this likelihood, one must also download the data files and prepare clik likelihood files from them. Look at :ref:`WMAP`. The ``--install_all_deps`` and ``--forceinstall_all_deps`` options will automatically download the sources, as if ``-wmap_install`` was set.


Special case: Healpix
^^^^^^^^^^^^^^^^^^^^^

Clik requires a specialy build healpix library. Namely, it must link with a repositionable (or better shared) version of the healpix library. 
This option is currently not available for the fortran version of the lib (as of version 2.20a). The configuration script knows how to produce this special version
of healpix for you. Thus except if you really know what you are doing, and even if you already have healpix installed on your system, 
using the option ``--healpix_install`` is very strongly recommanded.

Putting it all together
^^^^^^^^^^^^^^^^^^^^^^^

The following command::

    $> ./waf configure --install_all_deps

will tell the clik install system to install all the possible external dependency in the current directory. 

The following command::

    $> ./waf configure --lapack_mkl=/opt/intel/mkl \ --lapack_mkl_version=10.2
       --healpix_install --hdf5_install --h5py_install --gsl_prefix=/usr/local/gsl

will tell the clik install system to install healpix, hdf5 and h5py. The gsl library will be looked for in the unusual dir ``/usr/loca/gsl``. /All the other dependency will be looked up in the classical locations. The blas/lapack library 
will be the one from an mkl install located at --lapack_mkl=/opt/intel/mkl. Clik will be compiled in 64bit and installed in the current directory.

 
Best advanced choice 
^^^^^^^^^^^^^^^^^^^^

Use a mkl lapack install and let the other dependencies on auto install::

    $> ./waf configure --install_all_deps  \
          --lapack_mkl=/opt/intel/mkl --lapack_mkl_version=10.2 

This will use your mkl libraries from ``/opt/intel/mkl``, test if numpy, cython and gsl are installed on your computer (often the case) if not install them, 
and finally install all the other requirements (helpaix, hdf5 and its python wrapper).

Environment variables
---------------------

Depending of your shell, a configuration file named ``clik_profile.sh`` of ``clik_profile.csh`` will be installed in the ``bin`` directory at the install location of clik. One can source it on the command line, or include it in its startup configuration file to set the environment variable needed by clik.


