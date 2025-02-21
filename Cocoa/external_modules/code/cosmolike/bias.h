#ifndef __COSMOLIKE_BIAS_H
#define __COSMOLIKE_BIAS_H
#ifdef __cplusplus
extern "C" {
#endif

#define B1_PER_BIN 0
#define B1_PER_BIN_EVOLV 1
#define B1_PER_BIN_PASS_EVOLV 2
#define B1_GROWTH_SCALING 3
#define B1_POWER_LAW 4

#define B2_PER_BIN 0
#define B2_FROM_B1 1

#define BS2_PER_BIN 0
#define BS2_FROM_B1 1

#define B3_PER_BIN 0
#define B3_FROM_B1 1

#define BMAG_PER_BIN 0

// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// galaxy bias evolution with redshift
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------

double gb1(const double z, const int ni); // g = galaxy, b1 = linear galaxy bias

double gb2(const double z, const int ni);

double gbs2(const double z, const int ni);

double gb3(const double z, const int ni);

double gbmag(const double z, const int ni);

#ifdef __cplusplus
}
#endif
#endif // HEADER GUARD