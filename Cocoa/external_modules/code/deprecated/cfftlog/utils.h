#ifndef __CFFTLOG_UTILS_H
#define __CFFTLOG_UTILS_H
#ifdef __cplusplus
extern "C" {
#endif

void extrap_log_linear_cfft(double* fk, int N_origin, int N_extra, double* large_fk);

double gamma_lanczos_real_cfft(double z);

double lngamma_lanczos_real_cfft(double z);

#ifdef __cplusplus
}
#endif
#endif // HEADER GUARD