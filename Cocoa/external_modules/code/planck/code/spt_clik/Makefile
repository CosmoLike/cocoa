# to build 
# $> make

# to install the core utilities
# $> make install

# to install the optional python utilities
# $> make python_install

# cleanup
# $> make clean

# here are the thing that you should modify 
################################################################################################

# set your prefix to where you want to install clik.
# default is to let it in the current directory
PREFIX := $(shell pwd)

# set the path of the cfitsio lib. 
# bewared that cfitsio must have been compiled with the option "make shared"
CFITSIOPATH := /usr/local
#CFITSIOPATH := /softs/cfitsio/3.24
# you have a CFITSIO lib in a weird location, also set those
CFITSIO_INCPATH := $(CFITSIOPATH)/include
CFITSIO_LIBPATH := $(CFITSIOPATH)/lib

#define your compilers and stuff
CC = gcc
FC = ifort
#FC = gfortran

# ifort
# if you are using ifort set here where its lib are installed
# and check the runtime libs
# PLEASE note that gcc 4.9 and ifort <14.0.4 have an imcompativbility
# (see https://software.intel.com/en-us/articles/gcc-49-openmp-code-cannot-be-linked-with-intel-openmp-runtime)

# on my mac I got
IFORTLIBPATH = /usr/bin/ifort-2011-base/compiler/lib
IFORTRUNTIME = -L$(IFORTLIBPATH) -lintlc -limf -lsvml -liomp5 -lifportmt -lifcoremt -lpthread

# on a linux machine, ifort 11.1
#IFORTLIBPATH = /softs/intel/fce/11.1.075/lib/intel64
#IFORTRUNTIME = -L$(IFORTLIBPATH) -lintlc -limf -lsvml -liomp5 -lifport -lifcoremt -lpthread

# on a linux machine, ifort 14.0
#IFORTLIBPATH = /softs/intel/xe/composer_xe_2013_sp1.1.106/compiler/lib/intel64
#IFORTRUNTIME = -L$(IFORTLIBPATH) -lintlc -limf -lsvml -liomp5 -lifport -lifcoremt -lirc -lpthread

# gfortran
# if you are using gfortran set here where the lib are installed
# and check the runtime libs
GFORTRANLIBPATH = /usr/lib
GFORTRANRUNTIME = -L$(GFORTRANLIBPATH) -lgfortran -lgomp

# if you are on linux and using mkl, you need to set this 
MKLROOT = 
LAPACKLIBPATHMKL = 
#some example
#MKLROOT = /softs/intel/mkl/10.2.6.038/
# on mkl 10.3
#LAPACKLIBPATHMKL = -L$(MKLROOT)/lib/intel64
# on mkl 10.2
#LAPACKLIBPATHMKL = -L$(MKLROOT)/lib/em64t

#if you want to point to your own version of lapack set the following variables
#LAPACK = -L/some/path -lsomefortranlapack -lsomedependencyforyourlapack
#LAPACKLIBPATH = /some/path


# pretty colors (comment to remove pretty colors or try to change echo to echo -e)
COLORS = 1

#set echo to echo -e to have colourized output on some shell
ECHO = echo 


# what is the openmp option for your C compiler (leave empty to compile without openmp)
#COPENMP =
COPENMP = -fopenmp
# what is the openmp option for your F90 compiler (leave empty to cmpile without openmp)
FOPENMP = -openmp

# what is the 32/64 bit option for your C compiler (leave empty if you don't want to know)
CM64 =
#CM64 = -arch x86_64 #macos
#CM64 = -m64 #linux

# what is the 32/64 bit option for your F90 compiler (leave empty if you don't want to know)
FM64 = 
#FM64 = -arch x86_64 #macos
#FM64 = -m64

# set the variable to the python cli to compile and install the python tools
PYTHON = python


################################################################################################

# you should not need to modify anything below


#temporary dirs
BDIR := $(shell pwd)/buildir
ODIR := $(shell pwd)/buildir/tmp

# tools
LD = gcc
INSTALL = install

# get the os
UNAME := $(shell uname -s)

ifeq ($(UNAME),Darwin)
OS = macos
else
OS = linux
endif

#defines for macos
SOMACOS = dylib
LIBPATHNAMEMACOS = DYLD_LIBRARY_PATH
#defines for linux
SOLINUX = so
LIBPATHNAMELINUX = LD_LIBRARY_PATH

ifeq ($(OS),macos)
SO = $(SOMACOS)
LIBPATHNAME = $(LIBPATHNAMEMACOS)
else
SO = $(SOLINUX)
LIBPATHNAME = $(LIBPATHNAMELINUX)
endif

#ifort
IFORTMODULEPATH = -module

#gfortran
GFORTRANMODULEPATH = -J

# this picks either ifort or gfortran, change those lines to set FRUNTIME and FMODULEPATH for your special case
ifeq ($(FC),ifort)
FLIBPATH = $(IFORTLIBPATH)
FRUNTIME = $(IFORTRUNTIME)
FMODULEPATH = $(IFORTMODULEPATH)
FFLAGS =
else
FLIBPATH = $(GFORTRANLIBPATH)
FRUNTIME = $(GFORTRANRUNTIME)
FMODULEPATH = $(GFORTRANMODULEPATH)
FFLAGS = -ffree-line-length-0
endif


# some defines (shared, relocatable openmp, etc)
CFPIC = -fPIC
FFPIC = -fPIC

# check here that the SHARED variable contain the correct invocation for your CC
ifeq ($(OS),macos)
SHARED = -dynamiclib
else
SHARED = -shared -Bdynamic
endif

# get version of the code from the svn version
git describe --abbrev=12 --always > svnversion
VERSION = $(strip $(shell cat svnversion)) MAKEFILE
#VERSION = MAKEFILE

# some more defines
#macos
DEFINESMACOS = -D HAS_RTLD_DEFAULT
#linux
DEFINESLINUX = 

DEFINESCOMMON = -D HAS_LAPACK -D LAPACK_CLIK -D NOHEALPIX -D CLIK_LENSING -D 'CLIKSVNVERSION="$(VERSION)"' -D CAMSPEC_V1


ifeq ($(OS),macos)
DEFINES = $(DEFINESMACOS) $(DEFINESCOMMON)
ifndef CM64
CM64 = -arch x86_64
endif
ifndef FM64
FM64 = -arch x86_64
endif
else
DEFINES = $(DEFINESLINUX) $(DEFINESCOMMON)
ifndef CM64
CM64 = -m64
endif
ifndef FM64
FM64 = -m64
endif
endif

INCLUDES = -I$(CFITSIO_INCPATH)

# final CFLAG and FFLAGS
CFLAGS = $(CM64) $(COPENMP) $(CFPIC) $(DEFINES) -I src -I src/cldf -I src/minipmc -I src/lenslike/plenslike -I src/plik $(INCLUDES)
FFLAGS += $(FM64) $(FOPENMP) $(FFPIC) $(DEFINES) $(FMODULEPATH) $(ODIR)


# Lapack section

#macos I advise you to use the builtin blas lapack that are reasonnably efficient
LAPACKLIBPATHMACOS = /System/Library/Frameworks/Accelerate.framework/Versions/Current/Frameworks/vecLib.framework/Versions/Current
LAPACKMACOS = -L$(LAPACKLIBPATHMACOS) -lBLAS -lLAPACK

# mkl I am assuming that the env variable MKLROOT contains the MKL root path
# if not define it here

LAPACKMKLCORELIB = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core
LAPACKMKL = -L$(LAPACKLIBPATHMKL) $(LAPACKMKLCORELIB)  -liomp5 -lpthread -lm

#LAPACK_FUNC := dtrsv  dpotrf  dpotrs  dpotri  dtrtri  dtrmm  dtrmv  dgeqrf  dormqr  dsyev  dgesvd  dsymv  dgemv  dgemm  dsyrk  dsyr2k  daxpy  dtrsm  dsymm  dsyr  ddot
LAPACK_FUNC = $(shell python -c"print(\" \".join(open(\"waf_tools/lapack_funcs.txt\").read().strip().split()))")
MKL_TO_INCLUDE := $(addprefix -u ,$(addsuffix _,$(LAPACK_FUNC)))
MKL_LIB_FULLPATH := $(filter $(addsuffix .a,$(addprefix %/lib,$(subst -l,,$(filter -l%,$(LAPACKMKLCORELIB))))),$(wildcard $(subst -L,,$(filter -L%,$(LAPACKMKL)))/lib*.a))

# pick lapack version
LAPACKDEP =
LAPACK_INSTALL =

ifndef LAPACK
ifeq ($(OS),macos)
#macos lapack
LAPACK = $(LAPACKMACOS)
LAPACKLIBPATH = $(LAPACKLIBPATHMACOS)
LAPACKDEP =
endif

ifdef MKLROOT
#mkl !
LAPACK = -L$(BDIR) -llapack_clik
LAPACKLIBPATH = $(LAPACKLIBPATHMKL)
LAPACKDEP = $(BDIR)/liblapack_clik.$(SO)
LAPACK_INSTALL = -L$(PREFIX)/lib -llapack_clik
endif
endif

#if you want to point to your own version of lapack set the following variables
#LAPACK = -L/some/path -lsomefortranlapack -lsomedependencyforyourlapack
#LAPACKLIBPATH = /some/path
# leave this one empty
#LAPACKDEP = 

#CFITSIO link
CFITSIO =  -L$(CFITSIO_LIBPATH) -lcfitsio

#final LDFLAG
LDFLAG = $(CM64) $(CFITSIO) $(LAPACK) $(FRUNTIME) -ldl -lm -lpthread

# define some path to find the codes
SRCPATHLIST := src src/minipmc src/cldf src/camspec src/bflike src/simall src/lenslike/plenslike src/cmbonly src/gibbs src/actspt src/spt3g src/simall src/lowlike src/plik src/plik/component_plugin/rel2015 
vpath %.c $(SRCPATHLIST)
vpath %.f90  $(SRCPATHLIST)
vpath  %.F90 $(SRCPATHLIST)

# define color output if needed
ifeq ($(COLORS),1)
NO_COLOR=\x1b[0m
GREEN_COLOR=\x1b[32;11m
RED_COLOR=\x1b[31;01m
BLUE_COLOR=\x1b[36;11m
PINK_COLOR=\x1b[35m
endif

# all the code
TOOLS := $(addprefix $(ODIR)/,errorlist.o io.o distribution.o cldf.o clik_dic.o)
CLIKMAIN := $(addprefix $(ODIR)/,clik.o lklbs.o lowly_common.o clik_helper.o)

LENSLKL := $(addprefix $(ODIR)/,plenslike_dat_mono.o plenslike_dat_quad.o plenslike_dat_qecl.o plenslike_dat_full.o qest.o wignerd.o clik_lensing.o)
ACTSPTLKL := $(addprefix $(ODIR)/,Highell_options.f90.o Highell_subroutines.f90.o  Foregrounds_loading.f90.o ACT_equa_likelihood.f90.o SPT_reichardt_likelihood.f90.o ACT_south_likelihood.f90.o  SPT_keisler_likelihood.f90.o  Highell_likelihood.f90.o clik_actspt.f90.o clik_actspt.o)
CAMSPECLKL := $(addprefix $(ODIR)/,CAMspec.f90.o clik_CAMspec.f90.o clik_CAMspec.o)
LOWLIKELKL := $(addprefix $(ODIR)/,healpix_types.f90.o read_archive_map.f90.o read_fits.f90.o br_mod_dist.f90.o Planck_options.f90.o  Planck_teeebb_pixlike.f90.o  Planck_likelihood.f90.o clik_lowlike.f90.o clik_lowlike.o)
GIBBSLKL := $(addprefix $(ODIR)/,comm_br_mod.f90.o comm_gauss_br_mod.f90.o comm_gauss_br_mod_v3.f90.o comm_lowl_mod_dist.f90.o clik_gibbs.f90.o clik_gibbs.o)
BFLIKELKL := $(addprefix $(ODIR)/,long_intrinsic_smw.f90.o fitstools_smw.f90.o bflike_QUonly.f90.o bflike.f90.o bflike_smw.f90.o clik_bflike.f90.o clik_bflike.o)
PLIKLITELKL := $(addprefix $(ODIR)/,plik_cmbonly.f90.o clik_cmbonly.f90.o clik_cmbonly.o)
PLIKLKL := $(addprefix $(ODIR)/, smica.o clik_hfipack.o clik_parametric.o clik_parametric_addon.o fg2015.o corrnoise.o leakage.o)
SIMALLLKL := $(addprefix $(ODIR)/, clik_simall.o)
SPTLKL := $(addprefix $(ODIR)/, clik_spt3g.o clik_spt3g.f90.o)

CMBLKL:= $(ACTSPTLKL) $(CAMSPECLKL) $(GIBBSLKL) $(LOWLIKELKL) $(BFLIKELKL) $(SIMALLLKL) $(SPTLKL) $(PLIKLITELKL) $(PLIKLKL)
CLIKLIB := $(TOOLS) $(CLIKMAIN) $(CMBLKL) $(LENSLKL) $(LAPACKDEP)


all: $(BDIR)/libclik.$(SO) $(BDIR)/libclik_f90.$(SO) $(BDIR)/clik_example_C $(BDIR)/clik_example_f90

install_dir: 
	@mkdir -p $(PREFIX)/bin
	@mkdir -p $(PREFIX)/lib
	@mkdir -p $(PREFIX)/include
	@mkdir -p $(PREFIX)/share/clik

DATAPLIK := $(addprefix src/plik/component_plugin/rel2015/,tsz_143_eps0.50.dat sz_x_cib_template.dat ksz_fromcamspec.dat cib_1h_2h_100_353_Jsr-1_PS_2014_09.dat sky_template_v15_F100_143_217_353.dat cnoise_F100_143_217_353_v17.dat cleak_eh_rd12rc3_v1.dat cnoise_e2e_v2.dat sbpx_tmpl_v4_hm.dat)

install_data: | install_dir
	@mkdir -p $(PREFIX)/share/clik/rel2015
	@$(ECHO) "install template data $(BLUE_COLOR) $(DATAPLIK) $(NO_COLOR) in $(BLUE_COLOR)$(PREFIX)/share/clik/rel2015 $(NO_COLOR)"
	$(INSTALL) $(DATAPLIK) $(PREFIX)/share/clik/rel2015

install: $(BDIR)/libclik.$(SO) $(BDIR)/libclik_f90.$(SO) $(BDIR)/clik_example_C $(BDIR)/clik_example_f90 $(LAPACKDEP) $(BDIR)/clik_profile.sh $(BDIR)/clik_profile.csh $(BDIR)/clik-config $(BDIR)/clik-config_f90 install_data | install_dir
	@$(ECHO) "install libs $(BLUE_COLOR)libclik.$(SO) libclik_f90.$(SO)$(NO_COLOR) in $(BLUE_COLOR)$(PREFIX)/lib $(NO_COLOR)"
	@$(INSTALL)  $(BDIR)/libclik.$(SO) $(BDIR)/libclik_f90.$(SO) $(LAPACKDEP) $(PREFIX)/lib
	@$(ECHO) "install includes $(BLUE_COLOR)clik.h clik.mod$(NO_COLOR) in $(BLUE_COLOR)$(PREFIX)/include $(NO_COLOR)"
	@$(INSTALL)  src/clik.h src/minipmc/maths_base.h src/minipmc/errorlist.h src/minipmc/io.h src/lapack_clik.h src/minipmc/pmc.h $(ODIR)/clik.mod $(PREFIX)/include
	@$(ECHO) "install clik_profile & clik-config$(BLUE_COLOR)clik_profile.sh clik_profile.csh clik-config clik-config_f90$(NO_COLOR) in $(BLUE_COLOR)$(PREFIX)/bin $(NO_COLOR)"
	@$(INSTALL)  $(BDIR)/clik_profile.sh $(BDIR)/clik_profile.csh $(BDIR)/clik-config $(BDIR)/clik-config_f90 $(PREFIX)/bin
	@$(ECHO) "install exec tools $(BLUE_COLOR)clik_example_C clik_example_f90$(NO_COLOR) in $(BLUE_COLOR)$(PREFIX)/bin $(NO_COLOR)"
	@$(INSTALL)  $(BDIR)/clik_example_C $(BDIR)/clik_example_f90 $(PREFIX)/bin
	@$(ECHO) "\n$(PINK_COLOR)*----------------------------------------------------*"
	@$(ECHO) "$(PINK_COLOR)|$(NO_COLOR)                                                    $(PINK_COLOR)|"
	@$(ECHO) "$(PINK_COLOR)|$(NO_COLOR)   Source clik_profile.sh (or clik_profile.csh)     $(PINK_COLOR)|"
	@$(ECHO) "$(PINK_COLOR)|$(NO_COLOR)   to set the environment variables needed by clik  $(PINK_COLOR)|"
	@$(ECHO) "$(PINK_COLOR)|$(NO_COLOR)                                                    $(PINK_COLOR)|"
	@$(ECHO) "$(PINK_COLOR)*----------------------------------------------------*\n$(NO_COLOR)"

ifdef PYTHON
PYTHONPATH = $(PREFIX)/lib/`$(PYTHON) -c"import sys;print('python%s/site-packages'%sys.version[0:3])"`
PYTHONEXE := `which $(PYTHON)`
else
PYTHONPATH := 
endif

$(BDIR)/clik_profile.sh: src/clik_profile.sh.template |$(BDIR)
	@sed "s!PREFIX!$(PREFIX)!g;s/DYLD_LIBRARY_PATH/$(LIBPATHNAME)/g;s@CFITSIOLIBPATH@$(CFITSIO_LIBPATH)@g;s!FORTRANLIBPATH!$(FLIBPATH)!g;s!LAPACKLIBPATH!$(LAPACKLIBPATH)!g;s!MPYTHONPATH!$(PYTHONPATH)!g" <$< >$@

$(BDIR)/clik_profile.csh: src/clik_profile.csh.template |$(BDIR)
	@sed "s!PREFIX!$(PREFIX)!g;s/DYLD_LIBRARY_PATH/$(LIBPATHNAME)/g;s@CFITSIOLIBPATH@$(CFITSIO_LIBPATH)@g;s!FORTRANLIBPATH!$(FLIBPATH)!g;s!LAPACKLIBPATH!$(LAPACKLIBPATH)!g;s!MPYTHONPATH!$(PYTHONPATH)!g" <$< >$@

$(BDIR):
	@mkdir $(BDIR)

$(ODIR): | $(BDIR)
	@mkdir $(ODIR)

$(CLIKLIB): | $(ODIR) $(ODIR)/.print_info $(ODIR)/.test_cfitsio

$(BDIR)/libclik.$(SO): $(CLIKLIB) 
	@$(ECHO) "build $(BLUE_COLOR)$(@) $(NO_COLOR)"
	@$(LD)  $(SHARED)  $(LAPACK) $(LDFLAG) $^ -o $@

$(BDIR)/libclik_f90.$(SO): $(BDIR)/libclik.$(SO) $(addprefix $(ODIR)/,clik_fortran.o clik.f90.o)
	@$(ECHO) "build $(BLUE_COLOR)$(@) $(NO_COLOR)"
	@$(LD) $(SHARED)  $(LDFLAG) $(LAPACK) -L$(BDIR) -lclik $^ -o $@

$(BDIR)/clik_example_C: $(ODIR)/clik_example_c.o $(BDIR)/libclik.$(SO)
	@$(ECHO) "build $(BLUE_COLOR)$(@) $(NO_COLOR)"
	@$(CC) $(LDFLAG) $(LAPACK) -L$(BDIR) -lclik $< -o $@

$(BDIR)/clik_example_f90: $(ODIR)/clik_example_f90.f90.o $(BDIR)/libclik_f90.$(SO)
	@$(ECHO) "build $(BLUE_COLOR)$(@) $(NO_COLOR)"
	@$(FC) $(LDFLAG) $(LAPACK)  -L$(BDIR) -lclik_f90 -lclik $< -o $@

$(BDIR)/liblapack_clik.$(SO): |$(BDIR)
ifndef MKL_LIB_FULLPATH
	@$(ECHO) "$(RED_COLOR)I suspect an error with your MKLROOT, or MKL_LIB_FULLPATH, please check$(NO_COLOR)"
endif
	@$(ECHO) "build $(BLUE_COLOR)$(@) $(NO_COLOR),"
	@$(ECHO) "(see chapter 5 in http://software.intel.com/sites/products/documentation/hpc/mkl/lin/)"
	@$(ECHO) "using the following command line:"
	gcc $(SHARED)  $(MKL_TO_INCLUDE) -Wl,--start-group $(MKL_LIB_FULLPATH) -Wl,--end-group -L$(IFORTLIBPATH) -L/lib -L/lib64 -liomp5 -lpthread -lm -o $@

svnversion:
	@git describe --abbrev=12 --always > svnversion

$(ODIR)/%.o : %.c | svnversion
	@$(ECHO) "$(GREEN_COLOR)$< $(NO_COLOR) -> $(GREEN_COLOR) $(@) $(NO_COLOR)"
	@$(CC) -c $(CFLAGS) $< -o$(@)

$(ODIR)/%.f90.o : %.f90 | svnversion
	@$(ECHO) "$(GREEN_COLOR)$< $(NO_COLOR) -> $(GREEN_COLOR) $(@) $(NO_COLOR)"
	@$(FC) -c $(FFLAGS) $< -o$(@)

$(ODIR)/%.f90.o : %.F90 | svnversion
	@$(ECHO) "$(GREEN_COLOR)$< $(NO_COLOR) -> $(GREEN_COLOR) $(@) $(NO_COLOR)"
	@$(FC) -c $(FFLAGS) $< -o$(@)

$(ODIR)/%.py: src/python/tools/%.py 
	@sed "s@PYTHONEXE@$(PYTHONEXE)@g;s@REPLACEPATH@$(PYTHONPATH)@g" <$< >$@
	@$(INSTALL) $@ $(PREFIX)/bin/$(subst .py,,$(@F))
	@$(ECHO) "install python tools $(BLUE_COLOR)  $(subst .py,,$(@F)) $(NO_COLOR) in $(BLUE_COLOR)$(PREFIX)/bin $(NO_COLOR)"

$(ODIR)/.print_info: |$(ODIR)
	@$(ECHO) "\n$(BLUE_COLOR)Compile$(NO_COLOR) clik $(VERSION) "
	@$(ECHO) "$(BLUE_COLOR)Using $(NO_COLOR) CC = $(CC)"
	@$(ECHO) "$(BLUE_COLOR)Using $(NO_COLOR) FC = $(FC)"
	@$(ECHO) "$(BLUE_COLOR)Using $(NO_COLOR) CFLAGS = $(CFLAGS)"
	@$(ECHO) "$(BLUE_COLOR)Using $(NO_COLOR) FFLAGS = $(FFLAGS)"
	@$(ECHO) "$(BLUE_COLOR)Using the following lapack link line:$(NO_COLOR) $(LAPACK)"
	@$(ECHO) "$(BLUE_COLOR)Using the following cfitsio link line:$(NO_COLOR) $(CFITSIO)"
	@$(ECHO) "$(BLUE_COLOR)Using the following fortran runtime link line:$(NO_COLOR) $(FRUNTIME)"
	@$(ECHO) "$(BLUE_COLOR)Build dir:$(NO_COLOR) $(BDIR)"
	@$(ECHO)
	@touch $(@)

PYTOOLS := $(shell cd src/python/tools/;ls *.py;cd ../../../)

install_python: install $(addprefix $(ODIR)/, $(PYTOOLS)) |$(ODIR)
	@LINK_CLIK="$(LDFLAG) $(LAPACK) -L$(PREFIX)/lib -lclik " $(PYTHON) setup.py build --build-base=$(ODIR) install --install-lib=$(PYTHONPATH)
	@$(ECHO) "\n$(PINK_COLOR)*----------------------------------------------------*"
	@$(ECHO) "$(PINK_COLOR)|$(NO_COLOR)                                                    $(PINK_COLOR)|"
	@$(ECHO) "$(PINK_COLOR)|$(NO_COLOR)   Source clik_profile.sh (or clik_profile.csh)     $(PINK_COLOR)|"
	@$(ECHO) "$(PINK_COLOR)|$(NO_COLOR)   to set the environment variables needed by clik  $(PINK_COLOR)|"
	@$(ECHO) "$(PINK_COLOR)|$(NO_COLOR)                                                    $(PINK_COLOR)|"
	@$(ECHO) "$(PINK_COLOR)*----------------------------------------------------*\n$(NO_COLOR)"

HAS_CFITSIO_INC := $(shell [ -f $(CFITSIO_INCPATH)/fitsio.h ] && echo OK)
HAS_CFITSIO_LIB := $(shell [ -f $(CFITSIO_LIBPATH)/libcfitsio.$(SO) ] && echo OK)

$(ODIR)/.test_cfitsio: |$(ODIR)
ifneq ($(HAS_CFITSIO_INC),OK)
	@$(ECHO) "\n$(RED_COLOR)Cannot find cfisio includes ($(CFITSIO_INCPATH)/fitsio.h)$(NO_COLOR)"
	@false
endif
ifneq ($(HAS_CFITSIO_LIB),OK)
	@$(ECHO) "\n$(RED_COLOR)Cannot find cfisio lib ($(CFITSIO_LIBPATH)/libcfitsio.$(SO))$(NO_COLOR)"
	@false
endif
	@touch $(@)


INSTALL_CFLAG = $(subst ",\",$(subst ',\',$(CM64) $(COPENMP) $(CFPIC) $(DEFINES) -I $(PREFIX)/include $(INCLUDES)))
INSTALL_FFLAG = $(subst ",\",$(subst ',\',$(FM64) $(FOPENMP) $(FFPIC) $(DEFINES) $(FMODULEPATH) $(PREFIX)/include))
INSTALL_CLIB = $(CM64) $(CFITSIO) $(FRUNTIME) -ldl -lm -lpthread $(CFITSIO) $(FRUNTIME) 
INSTALL_FLIB = $(CM64) $(CFITSIO) $(FRUNTIME) -ldl -lm -lpthread $(CFITSIO) $(FRUNTIME) 
ifdef LAPACK_INSTALL
INSTALL_CLIB += $(LAPACK_INSTALL)
INSTALL_FLIB += $(LAPACK_INSTALL)
else
INSTALL_CLIB += $(LAPACK)
INSTALL_FLIB += $(LAPACK)
endif
INSTALL_CLIB += -L$(PREFIX)/lib -lclik
INSTALL_FLIB += -L$(PREFIX)/lib -lclik -lclik_f90

$(BDIR)/clik-config: src/clik-config.template |$(BDIR)
	@sed "s@CFLAG@$(INSTALL_CFLAG)@g;s@LIB@$(INSTALL_CLIB)@g" <$< >$@

$(BDIR)/clik-config_f90: src/clik-config.template |$(BDIR)
	@sed "s@CFLAG@$(INSTALL_FFLAG)@g;s@LIB@$(INSTALL_FLIB)@g" <$< >$@


clean:
	@$(ECHO) "$(BLUE_COLOR)Removing all in $(BDIR)$(NO_COLOR)"
	@rm -rf $(BDIR)

.PHONY :clean  LAPACK_PRINT LAPACK_DEP
