/* emulator.cxx
*  =========
*  This file is part of EuclidEmulator2
*  Copyright (c) 2020 Mischa Knabenhans
*
*  EuclidEmulator2 is free software: you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation, either version 3 of the License, or
*  (at your option) any later version.
*
*  EuclidEmulator2 is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <iostream>
#include <fstream>
#include <string>
#include <sys/stat.h> // struct stat and fstat() function
#include <sys/mman.h> // mmap() function
#include <fcntl.h>    // declaration of O_RDONLY
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline2d.h>
#include <gsl/gsl_sf_legendre.h>
#include <math.h>
#include <assert.h>
#include "emulator.h"

#include <errno.h>
#include <string.h>

#include <stdio.h>
#include <limits.h>
#include <unistd.h>

using namespace std;

/* CONSTRUCTOR */
EuclidEmulator::EuclidEmulator():
	lmax(16),
	nz(101),
	nk(613),
	n_coeffs{53, 53, 117, 117, 53, \
			117, 117, 117, 117, 521, \
			117, 1539, 173, 457}
	{

	read_in_ee2_data_file();

	pc_2d_interp();

}

/* DESTRUCTOR */
EuclidEmulator::~EuclidEmulator(){
	for(int i=0; i<15; i++) {
		gsl_spline2d_free(logklogz2pc_spline[i]);
		gsl_interp_accel_free(logk2pc_acc[i]);
		gsl_interp_accel_free(logz2pc_acc[i]);
	}
	for (int ipar=0; ipar < 8; ipar++){
		delete [] univ_legendre[ipar];
	}
	// Free data stored when file is opened. Size is taken from having previously printed it out.
	munmap(pc[0], 7698560);

}

/* FUNCTION TO READ IN THE DATA FILE */
void EuclidEmulator::read_in_ee2_data_file(){
	/// VARIABLE DECLARATIONS ///
	off_t size;
    struct stat s;
	double *data;
	double *kptr;
	int i, ik, iz, idx = 0;

	for(iz = 0; iz < nz; iz++){
		for(ik = 0; ik < nk; ik++){
			this->Bvec[iz][ik] = 0.0; // initialize the resulting nlc vector
		}
	}

	// ==== LOAD EUCLIDEMULATOR2 DATA FILE ==== //
	int fp = open(PATH_TO_EE2_DATA_FILE1, O_RDONLY);

	if(fp==-1) {
		fp = open(PATH_TO_EE2_DATA_FILE2, O_RDONLY);
	}
  if(fp==-1) {
	  fp = open("./ee2_bindata.dat", O_RDONLY);
  }
	if(fp==-1) {
		cerr << "Unable to open ee2_bindata.dat\n";
		cerr << "Error is: " << strerror(errno) << endl;
    exit(1);
	}

	// Get the size of the file. //
    int status = fstat(fp, & s);
    size = s.st_size;

	// Map the file into memory //
	data = (double *) mmap (0, size, PROT_READ, MAP_PRIVATE, fp, 0);

	// Reading in principal components //
	for (i=0;i<15;i++) {
    	this->pc[i] = &data[idx];  // pc[0] = PCA mean
    	idx += nk*nz;
    }
	// Reading in PCE coefficients //
	for (i=0;i<14;i++) {
    	this->pce_coeffs[i] = &data[idx];
    	idx += n_coeffs[i];
    }

	// Reading in PCE multi-indices //
	for (i=0;i<14;i++) {
    	this->pce_multiindex[i] = &data[idx];
    	idx += 8*n_coeffs[i];
    }
	// vector of k modes
    kptr = &data[idx];

    for (i=0;i<nk;i++) {
		this->kvec[i] = kptr[i];
    }
    idx += nk;
	// Check if all data has been read in properly
	assert(idx == size/sizeof(double));

	close(fp);

}

/* 2D INTERPOLATION OF PRINCIPAL COMPONENTS */
void EuclidEmulator::pc_2d_interp(){
	double logk[nk];
	double stp[nz];
	int i;

	for (i=0; i<nk; i++) logk[i] = log(this->kvec[i]);
	for (i=nz-1; i>=0; i--) stp[i] = i;
	for (int i=0; i<15; i++){
		logk2pc_acc[i] = gsl_interp_accel_alloc();
        logz2pc_acc[i] = gsl_interp_accel_alloc();
    	logklogz2pc_spline[i] = gsl_spline2d_alloc(gsl_interp2d_bicubic, nk, nz);
    	gsl_spline2d_init(logklogz2pc_spline[i], logk, stp, pc[i], nk, nz);
	}
	for (int ipar=0; ipar < 8; ipar++){
		univ_legendre[ipar] = new double[lmax+1];
	}
}

/* COMPUTE NLC */
//void EuclidEmulator::compute_nlc(Cosmology csm, double* redshift, int n_redshift, double* kmodes, int n_kmodes){
void EuclidEmulator::compute_nlc(Cosmology* csm, vector<double> redshift, int n_redshift){

	double pc_weight;
	double basisfunc;
	double stp_no[n_redshift];

	//printf("There are %d redshifts to compute.\n", n_redshift);

	// Convert all redshifts into step numbers
	// As we are looping through all redshifts anyway, we can just
	// as well use the same loop to declare Bvec[iz]
	for(int iz=0; iz<n_redshift; iz++) {
		if(redshift.at(iz) > 10.0 || redshift.at(iz) < 0.0){
			std::cout << "ERROR: EuclidEmulator2 accepts only redshifts in the interval [0.0, 10.0]\n" \
					  << "The current redshift z = " << redshift.at(iz) << " is therefore ignored." << std::endl;
			continue;
		}
		stp_no[iz] = csm->compute_step_number(redshift.at(iz));
        //printf("nStep(%.2f) = %.4f\n", redshift.at(iz), stp_no[iz]);
	}
	//printf("Redshifts mapped to nStep\n");

	// Pre-compute all Legendre polynomials up to order lmax
	for (int ipar=0; ipar < 8; ipar++){
		//univ_legendre[ipar] = new double[lmax+1];
		gsl_sf_legendre_Pl_array(lmax, csm->cosmo_tf[ipar], univ_legendre[ipar]);
		for (int l=0; l<=lmax; l++){
			univ_legendre[ipar][l] *= sqrt(2.0*l + 1.0); //normalization
		}
	}
	//printf("Legendre polynomials computed\n");

	// Initialize with PCA mean
	for(int iz=0; iz<n_redshift; iz++){
		for(int ik=0; ik<nk; ik++){
			Bvec[iz][ik] = gsl_spline2d_eval(logklogz2pc_spline[0], log(this->kvec[ik]), stp_no[iz], logk2pc_acc[0], logz2pc_acc[0]);
		}
		//printf("B(k_max, z=%.2f) = %.2f\n", redshift.at(iz), Bvec[iz][nk-1]);
	}
	//printf("PCA initialized\n");

	// Loop over principal components
	for(int ipc=1; ipc<15; ipc++){
		pc_weight = 0.0;
		// assemble PCE to get the PCA weight according
        // to inner sum of eq. 27 in EE2 paper
        for(int ic=0; ic<n_coeffs[ipc-1]; ic++){
        	basisfunc = 1.0;
            for(int ipar=0; ipar<8 ; ipar++){
                basisfunc *= univ_legendre[ipar][int(pce_multiindex[ipc-1][ic*8 + ipar])];
            }
            pc_weight += pce_coeffs[ipc-1][ic]*basisfunc;
        }
		// assemble PCA to get the final NLC according
        // to outer sum of eq. 27 in EE2 paper
		for(int iz=0; iz<n_redshift; iz++){
			//printf("Treating z[%d] = %.3f ==> nStep = %.3f\n", iz, redshift.at(iz), stp_no[iz]);
			for(int ik=0; ik<nk; ik++){
				Bvec[iz][ik] += (pc_weight*gsl_spline2d_eval(logklogz2pc_spline[ipc], log(this->kvec[ik]), stp_no[iz], logk2pc_acc[ipc], logz2pc_acc[ipc]));
			}
		}
	}

	//printf("PCA assembled\n");
}

/* WRITE NLC TO FILE */
void EuclidEmulator::write_nlc2file(const string& filename, vector<double> zvec, int n_redshift){
	ofstream fp_out (filename);
	// Writing an informative header line
	string header = "#k [h/Mpc]";
	for(int iz=0; iz<n_redshift; iz++){
		header.append("\tB(k,z="+to_string(zvec[iz])+")");
	}
	fp_out << header << endl;

	// Writing the data
	if (fp_out.is_open()){
		for(int ik=0; ik<nk; ik++){
			string line = to_string(this->kvec[ik]);
			for(int iz=0; iz<n_redshift; iz++){
				line.append("\t");
  				line.append(to_string(pow(10.0,Bvec[iz][ik])));
			}
			fp_out << line << endl;
		}
	}
	else cout << "Unable to open file" << endl;
}

/* PRINT INFO */
void EuclidEmulator::print_info(){
	int i, j, ip, ic;

	for (i=1;i<14;++i) {
    	fprintf(stderr,"%1d ",i);
    	for (j=0;j<2*n_coeffs[i]-1;++j) fprintf(stderr,"-");
    	fprintf(stderr,"\n");
    	for (ic=0;ic<n_coeffs[i];++ic) {
        	fprintf(stderr,"%.3g ",pce_coeffs[i][ic]);
        }
    	fprintf(stderr,"\n");
    	for (j=0;j<2*n_coeffs[i]-1;++j) fprintf(stderr,"-");
    	fprintf(stderr,"\n");
    	for (ip=0;ip<8;++ip) {
      		for (ic=0;ic<n_coeffs[i];++ic) {
        		fprintf(stderr,"%1d ",(int)pce_multiindex[i][ic*8 + ip]);
      		}
      		fprintf(stderr,"\n");
    	}
    	fprintf(stderr,"\n");
    }
}
