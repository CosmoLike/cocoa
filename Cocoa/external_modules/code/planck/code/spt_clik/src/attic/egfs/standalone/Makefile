default: libegfs.so

FFLAGS = -m64 -fPIC -shared -g -fall-intrinsics -ffree-line-length-none 
libegfs.so: keysvalues.f90 egfs.f90 clik_egfs.f90
	gfortran $(FFLAGS) -c keysvalues.f90 -o keysvalues_f90.o
	gfortran $(FFLAGS) -c egfs.f90 -o egfs_f90.o
	gfortran $(FFLAGS) -c clik_egfs.f90 -o clik_egfs_f90.o
	gfortran -m64 -shared -g keysvalues_f90.o egfs_f90.o clik_egfs_f90.o  -o libegfs.so
	
clean: 
	-rm -f keysvalues_f90.o egfs_f90.o clik_egfs_f90.o *.mod libegfs.so
