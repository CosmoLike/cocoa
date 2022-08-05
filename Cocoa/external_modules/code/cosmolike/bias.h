#ifndef __COSMOLIKE_BIAS_H
#define __COSMOLIKE_BIAS_H
#ifdef __cplusplus
extern "C" {
#endif

// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// bias evolution with redshift
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------

// model b1 using one non-evoling parameter per redshift bin
double b1_per_bin(double z, int ni);

// model b1 using one parameter per redshift
// bin, power-law evolution within each bin
double b1_per_bin_evolv(double z, int ni);


// model b1 using one parameter per redshift bin + passive evolution
double b1_per_bin_pass_evolv(double z, int ni);

// model b1 assuming b1(z) = b_{1,0}*G(z)
double b1_growth_scaling(double z, int ni);

// bias evolution within redshift bin, used by
// clustering/G-G-lensing routines without HOD modeling
double bgal_z(double z, int ni);

double b2_from_b1(double b1); // fitting function for b_2(b_1)

double bs2_from_b1(double b1); // theory prediction for b_s2(b_1)

double b3nl_from_b1(double b1); // theory prediction b3nl(b1) = b1-1

#ifdef __cplusplus
}
#endif
#endif // HEADER GUARD