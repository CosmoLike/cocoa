#include "clik_parametric.h"
#include "clik_parametric_addon.h"

#define PRM_NU0 143.
void powerlaw_triangle_norm_derivative(parametric * egl, int iv, double *Rq, double *dRq, error **err);
void powerlaw_tanh_norm_derivative(parametric * egl, int iv, double *Rq, double *dRq, error **err);
double powerlaw_tensor_step_index(parametric* egl, int m1,int m2);