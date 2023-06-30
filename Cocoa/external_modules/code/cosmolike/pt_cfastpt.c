#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cfastpt/cfastpt.h"
#include "basics.h"
#include "cosmo3D.h"
#include "pt_cfastpt.h"
#include "recompute.h"
#include "structs.h"

#include "log.c/src/log.h"

static double tab_d1d3[800] = {
 3.910971e-18,  4.182876e-18,  4.466778e-18,  4.789119e-18,
 5.133766e-18,  5.508468e-18,  5.890354e-18,  6.301033e-18,  6.736776e-18,  7.223413e-18,
 7.741152e-18,  8.297983e-18,  8.873722e-18,  9.495336e-18,  1.016005e-17,  1.089240e-17, 
 1.166949e-17,  1.250109e-17,  1.337112e-17,  1.431206e-17,  1.532126e-17,  1.642177e-17,
 1.758856e-17,  1.883528e-17,  2.015129e-17,  2.157442e-17,  2.310156e-17,  2.475463e-17,
 2.650821e-17,  2.838234e-17,  3.037282e-17,  3.252296e-17,  3.482887e-17,  3.731253e-17,
 3.995090e-17,  4.277283e-17,  4.578202e-17,  4.902737e-17,  5.250472e-17,  5.623852e-17,
 6.021148e-17,  6.446465e-17,  6.901086e-17,  7.390553e-17,  7.914593e-17,  8.476280e-17,
 9.074936e-17,  9.716243e-17,  1.040262e-16,  1.114045e-16,  1.192997e-16,  1.277547e-16,
 1.367788e-16,  1.464498e-16,  1.568067e-16,  1.679257e-16,  1.798203e-16,  1.925535e-16,
 2.061593e-16,  2.207428e-16,  2.363632e-16,  2.531167e-16,  2.710377e-16,  2.902211e-16,
 3.107366e-16,  3.327257e-16,  3.562774e-16,  3.815202e-16,  4.085236e-16,  4.374315e-16, 
 4.683643e-16,  5.015149e-16,  5.370178e-16,  5.750523e-16,  6.157460e-16,  6.593147e-16,
 7.059517e-16,  7.559237e-16,  8.094360e-16,  8.667483e-16,  9.280793e-16,  9.937481e-16,
 1.064056e-15,  1.139381e-15,  1.220034e-15,  1.306400e-15,  1.398836e-15,  1.497819e-15,
 1.603806e-15,  1.717337e-15,  1.838892e-15,  1.969048e-15,  2.108370e-15,  2.257566e-15,
 2.417324e-15,  2.588429e-15,  2.771623e-15,  2.967772e-15,  3.177757e-15,  3.402623e-15,
 3.643414e-15,  3.901278e-15,  4.177355e-15,  4.472956e-15,  4.789437e-15,  5.128343e-15,
 5.491237e-15,  5.879840e-15,  6.295897e-15,  6.741376e-15,  7.218341e-15,  7.729100e-15,
 8.276001e-15,  8.861614e-15,  9.488599e-15,  1.015993e-14,  1.087872e-14,  1.164841e-14,  
 1.247254e-14,  1.335498e-14,  1.429977e-14,  1.531137e-14,  1.639451e-14,  1.755432e-14,  
 1.879613e-14,  2.012575e-14,  2.154932e-14,  2.307357e-14,  2.470557e-14,  2.645301e-14,  
 2.832397e-14,  3.032717e-14,  3.247190e-14,  3.476827e-14,  3.722695e-14,  3.985946e-14,
 4.267797e-14,  4.569567e-14,  4.892653e-14,  5.238568e-14,  5.608923e-14,  6.005448e-14,
 6.429970e-14,  6.884470e-14,  7.371067e-14,  7.892032e-14,  8.449782e-14,  9.046915e-14,
 9.686203e-14,  1.037061e-13,  1.110331e-13,  1.188773e-13,  1.272749e-13,  1.362649e-13,
 1.458890e-13,  1.561921e-13,  1.672216e-13,  1.790287e-13,  1.916679e-13,  2.051984e-13,
 2.196820e-13,  2.351860e-13,  2.517825e-13,  2.695480e-13,  2.885645e-13,  3.089191e-13,
 3.307068e-13,  3.540279e-13,  3.789886e-13,  4.057053e-13,  4.343016e-13,  4.649073e-13,
 4.976628e-13,  5.327211e-13,  5.702419e-13,  6.103966e-13,  6.533695e-13,  6.993587e-13,
 7.485736e-13,  8.012380e-13,  8.575955e-13,  9.179030e-13,  9.824331e-13,  1.051480e-12,
 1.125362e-12,  1.204414e-12,  1.288989e-12,  1.379476e-12,  1.476289e-12,  1.579862e-12,
 1.690659e-12,  1.809186e-12,  1.935983e-12,  2.071613e-12,  2.216680e-12,  2.371855e-12,
 2.537830e-12,  2.715330e-12,  2.905159e-12,  3.108183e-12,  3.325293e-12,  3.557438e-12,
 3.805677e-12,  4.071119e-12,  4.354910e-12,  4.658304e-12,  4.982675e-12,  5.329436e-12,
 5.700068e-12,  6.096228e-12,  6.519693e-12,  6.972257e-12,  7.455855e-12,  7.972665e-12,
 8.524937e-12,  9.114957e-12,  9.745293e-12,  1.041877e-11,  1.113821e-11,  1.190655e-11,
 1.272719e-11,  1.360374e-11,  1.453974e-11,  1.553907e-11,  1.660614e-11,  1.774543e-11,
 1.896143e-11,  2.025929e-11,  2.164471e-11,  2.312318e-11,  2.470046e-11,  2.638335e-11,
 2.817901e-11,  3.009426e-11,  3.213666e-11,  3.431506e-11,  3.663812e-11,  3.911438e-11,
 4.175395e-11,  4.456807e-11,  4.756713e-11,  5.076206e-11,  5.416619e-11,  5.779340e-11,
 6.165619e-11,  6.576908e-11,  7.014935e-11,  7.481305e-11,  7.977574e-11,  8.505704e-11,
 9.067847e-11,  9.665871e-11,  1.030175e-10,  1.097809e-10,  1.169744e-10,  1.246197e-10,
 1.327440e-10,  1.413797e-10,  1.505555e-10,  1.602981e-10,  1.706442e-10,  1.816334e-10,
 1.932977e-10,  2.056710e-10,  2.188017e-10,  2.327340e-10,  2.475040e-10,  2.631596e-10,
 2.797599e-10,  2.973530e-10,  3.159795e-10,  3.357058e-10,  3.566016e-10,  3.787166e-10,
 4.021044e-10,  4.268516e-10,  4.530301e-10,  4.806911e-10,  5.099134e-10,  5.408041e-10,
 5.734293e-10,  6.078450e-10,  6.441635e-10,  6.825040e-10,  7.229214e-10,  7.654973e-10,
 8.103779e-10,  8.576636e-10,  9.074091e-10,  9.597406e-10,  1.014830e-09,  1.072749e-09,
 1.133551e-09,  1.197424e-09,  1.264528e-09,  1.334908e-09,  1.408653e-09,  1.485997e-09,
 1.567055e-09,  1.651829e-09,  1.740492e-09,  1.833311e-09,  1.930300e-09,  2.031443e-09,
 2.137024e-09,  2.247260e-09,  2.362022e-09,  2.481426e-09,  2.605786e-09,  2.735158e-09,
 2.869325e-09,  3.008539e-09,  3.153150e-09,  3.302949e-09,  3.457718e-09,  3.617862e-09,
 3.783529e-09,  3.954278e-09,  4.130068e-09,  4.311383e-09,  4.497999e-09,  4.689217e-09,
 4.885375e-09,  5.086866e-09,  5.292994e-09,  5.502998e-09,  5.717604e-09,  5.936659e-09,
 6.159071e-09,  6.384533e-09,  6.613920e-09,  6.846239e-09,  7.080067e-09,  7.315722e-09,
 7.554118e-09,  7.793413e-09,  8.032312e-09,  8.272199e-09,  8.512558e-09,  8.751289e-09,
 8.987998e-09,  9.224358e-09,  9.458389e-09,  9.688348e-09,  9.915256e-09,  1.014032e-08,
 1.036112e-08,  1.057713e-08,  1.079017e-08,  1.099941e-08,  1.120347e-08,  1.140262e-08,
 1.159826e-08,  1.178967e-08,  1.197683e-08,  1.216077e-08,  1.234352e-08,  1.252460e-08,
 1.270536e-08,  1.288644e-08,  1.306892e-08,  1.325413e-08,  1.344348e-08,  1.363640e-08,
 1.383587e-08,  1.404387e-08,  1.426042e-08,  1.448590e-08,  1.472414e-08,  1.497596e-08,
 1.523894e-08,  1.551527e-08,  1.580494e-08,  1.610922e-08,  1.642317e-08,  1.674735e-08,
 1.708034e-08,  1.741793e-08,  1.775793e-08,  1.809333e-08,  1.842083e-08,  1.873667e-08,
 1.903524e-08,  1.930225e-08,  1.953807e-08,  1.974093e-08,  1.990135e-08,  2.000492e-08,
 2.006105e-08,  2.007143e-08,  2.002413e-08,  1.992219e-08,  1.978235e-08,  1.960302e-08,
 1.939254e-08,  1.916213e-08,  1.892781e-08,  1.869676e-08,  1.849077e-08,  1.831458e-08,
 1.817593e-08,  1.809311e-08,  1.807251e-08,  1.809947e-08,  1.818275e-08,  1.832320e-08,
 1.849831e-08,  1.869062e-08,  1.888601e-08,  1.905009e-08,  1.917073e-08,  1.922758e-08,
 1.918755e-08,  1.904472e-08,  1.882046e-08,  1.850958e-08,  1.811549e-08,  1.768764e-08,
 1.725675e-08,  1.685279e-08,  1.651102e-08,  1.624339e-08,  1.606020e-08,  1.597476e-08,
 1.596559e-08,  1.600036e-08,  1.604524e-08,  1.605254e-08,  1.599359e-08,  1.584568e-08,
 1.559096e-08,  1.524024e-08,  1.484203e-08,  1.443508e-08,  1.405558e-08,  1.374559e-08,
 1.350801e-08,  1.335093e-08,  1.325448e-08,  1.317625e-08,  1.307030e-08,  1.290136e-08,
 1.265488e-08,  1.235120e-08,  1.202121e-08,  1.169558e-08,  1.141953e-08,  1.119538e-08,
 1.102557e-08,  1.087283e-08,  1.070563e-08,  1.050542e-08,  1.025642e-08,  9.995836e-09,
 9.731042e-09,  9.494515e-09,  9.290400e-09,  9.113016e-09,  8.922902e-09,  8.726490e-09,
 8.510911e-09,  8.290442e-09,  8.075002e-09,  7.877449e-09,  7.693646e-09,  7.512469e-09,
 7.325777e-09,  7.130147e-09,  6.934499e-09,  6.755936e-09,  6.578155e-09,  6.409595e-09,
 6.236333e-09,  6.064748e-09,  5.894874e-09,  5.731432e-09,  5.572715e-09,  5.416408e-09,
 5.260924e-09,  5.109414e-09,  4.964454e-09,  4.822387e-09,  4.676328e-09,  4.537761e-09,
 4.397341e-09,  4.275379e-09,  4.143258e-09,  4.020034e-09,  3.890930e-09,  3.777839e-09,
 3.654600e-09,  3.534674e-09,  3.427619e-09,  3.316918e-09,  3.212418e-09,  3.108418e-09,  
 3.004245e-09,  2.915067e-09,  2.808451e-09,  2.720520e-09,  2.629657e-09,  2.548750e-09,
 2.455967e-09,  2.374318e-09,  2.294497e-09,  2.211696e-09,  2.131324e-09,  2.055306e-09,
 1.983690e-09,  1.912214e-09,  1.842391e-09,  1.785656e-09,  1.725352e-09,  1.666089e-09,
 1.604502e-09,  1.558114e-09,  1.472218e-09,  1.449952e-09,  1.385192e-09,  1.338616e-09,
 1.278733e-09,  1.244693e-09,  1.186379e-09,  1.146996e-09,  1.095423e-09,  1.060161e-09,
 1.021840e-09,  9.849930e-10,  9.507326e-10,  9.127424e-10,  8.747899e-10,  8.477043e-10,
 8.056461e-10,  7.744497e-10,  7.290331e-10,  7.121668e-10,  6.903934e-10,  6.674064e-10,
 6.281759e-10,  6.010071e-10,  5.749332e-10,  5.554663e-10,  5.298597e-10,  5.189734e-10,
 4.964601e-10,  4.809650e-10,  4.507484e-10,  4.472166e-10,  3.958664e-10,  4.131065e-10,  
 3.739812e-10,  3.731632e-10,  3.615851e-10,  3.489243e-10,  3.149575e-10,  3.094438e-10,
 2.887533e-10,  2.757441e-10,  2.716195e-10,  2.644490e-10,  2.485680e-10,  2.458788e-10,
 2.097910e-10,  2.333494e-10,  2.033126e-10,  1.858978e-10,  1.741653e-10,  1.855059e-10,
 1.690043e-10,  1.727941e-10,  1.485467e-10,  1.370789e-10,  1.481281e-10,  1.398454e-10,
 1.037015e-10,  1.195089e-10,  1.328648e-10,  1.346805e-10,  9.696735e-11,  1.133417e-10,
 7.546955e-11,  1.126585e-10,  7.378848e-11,  9.953503e-11,  6.532070e-11,  9.209956e-11,
 6.487498e-11,  6.532770e-11,  4.733612e-11,  8.576345e-11,  5.489123e-11,  9.260568e-11,
 7.600372e-11,  4.889013e-11,  5.337598e-11,  5.098376e-11,  2.873587e-11,  3.305308e-11,
 3.763090e-11,  4.290017e-11,  2.498463e-11,  4.902413e-11,  1.021268e-11,  2.354661e-11,
 2.696824e-11,  3.225043e-11, -1.558926e-12,  3.863233e-11,  4.937074e-11,  4.689884e-11,
 1.216311e-11,  3.139190e-11,  1.029372e-11,  7.390549e-12, -2.375731e-11,  4.697778e-11,
 2.998759e-12,  8.898919e-12,  3.034634e-12, -7.404769e-13,  2.543173e-11, -1.259460e-12,
-1.860863e-11,  2.071894e-11,  2.333235e-11,  6.307092e-13, -1.339409e-11, -1.071323e-11, 
 2.772548e-11, -1.337781e-11,  1.499155e-12,  8.061537e-13, -6.163683e-12, -1.508189e-11,
-3.111544e-11, -1.480802e-11, -1.843287e-11,  5.841283e-12, -2.366329e-11, -4.346045e-12,
 2.233708e-11, -3.300497e-11, -1.337509e-11,  3.595905e-11,  1.032686e-11, -3.353857e-11,
-2.810819e-11, -3.542621e-11, -4.214006e-11,  2.188873e-11,  4.306275e-12,  2.296783e-11,
 8.032307e-12, -3.207177e-11,  3.874715e-11,  3.834196e-11,  4.142161e-11, -1.285456e-11,
 6.229334e-12,  2.723868e-11,  3.837333e-11, -1.746537e-11,  1.543089e-11,  3.763768e-11,
 4.253352e-11,  5.613461e-11,  2.585776e-11,  6.655403e-12,  1.645950e-11,  4.789137e-11,
 5.834712e-11,  1.040917e-11,  4.302331e-12,  5.299197e-11, -2.310800e-12,  5.641187e-11,
 4.248075e-11,  1.065697e-11, -1.199625e-11,  8.731526e-11,  2.885130e-11,  1.500047e-11,
-1.012494e-11, -6.677577e-11, -1.255079e-11, -9.465239e-12,  5.511950e-11,  5.732885e-11,
 7.122139e-11, -4.896827e-11,  3.573271e-11,  2.349440e-11,  3.891886e-11, -2.483904e-11,
 5.467823e-11,  5.348753e-11, -4.485613e-11,  7.102118e-11,  3.034649e-11, -6.745813e-12,
 6.752224e-12,  6.130347e-11, -6.302750e-11,  1.372323e-11, -1.137712e-11, -5.648667e-11,
-3.766979e-11,  6.993598e-11, -3.987549e-11,  3.700542e-12,  3.789813e-11,  9.736969e-11,
 5.728901e-12,  1.074817e-10,  2.648058e-12,  2.012776e-11,  6.875324e-11, -5.272685e-11,
-2.073128e-11,  3.198252e-11,  1.871731e-11,  7.345423e-11,  4.530856e-11, -7.082568e-11,
 1.207592e-11, -3.980394e-11,  5.667161e-13,  5.317470e-11,  7.651207e-11,  1.348856e-10,
-4.019966e-11, -6.497377e-12,  7.451400e-12,  8.908720e-11, -8.415314e-11,  1.818953e-11,
 4.746137e-11, -7.546239e-11, -3.478833e-11,  1.413114e-11, -6.129904e-12, -5.290320e-11,
-7.740161e-11,  3.485255e-11,  1.126023e-10,  2.201863e-11,  1.361764e-10, -7.434679e-11,
 4.632411e-11,  2.665100e-11, -1.277972e-10, -4.906905e-11, -1.204729e-10,  4.313078e-11,
 1.016915e-10,  1.498392e-10, -7.396595e-11, -2.384477e-11, -1.684383e-11,  8.601377e-13,
-1.705284e-11, -7.631001e-12, -2.219474e-11, -6.215454e-11, -8.256433e-12,  6.812224e-11,
 6.907224e-11,  3.283077e-11, -2.583649e-11, -4.747581e-11,  2.512696e-11,  5.003908e-11,
-1.049408e-10,  6.498571e-11, -7.710137e-12,  1.002388e-10, -5.816910e-11,  1.330833e-10,
 7.744255e-11, -7.337186e-11, -5.724791e-12,  7.066245e-11,  8.812115e-11, -2.133001e-11,
 9.036031e-12, -1.468138e-10,  1.120680e-10, -9.316044e-12,  4.694973e-11,  5.859625e-11,
 1.375245e-10,  1.570844e-10,  6.863663e-11, -4.041614e-11,  8.088534e-11,  4.535038e-11,
-7.937927e-11,  3.631226e-11,  2.014300e-10,  2.837585e-11};

double K_CH0(double k_mpch) { 
  return k_mpch * cosmology.coverH0; 
}

void FPT_input(double k[FPT.N], double P[FPT.N]) 
{
  if ((int) log10(FPT.k_max / FPT.k_min) * FPT.N_per_dec != FPT.N) 
  {
    log_fatal("inconsistent k-range and number of bins for FPT");
    log_fatal("FPT.k_min=%e, FPT.k_max=%e, FPT.N_per_dec=%d; FPT.N=%d", FPT.k_min, FPT.k_max, 
      FPT.N_per_dec, FPT.N);
    exit(1);
  }
  if (FPT.k_min < limits.k_min_cH0 || FPT.k_max > limits.k_max_cH0) 
  {
    log_fatal("k_min/k_max out of range FPT.k_min = %e, FPT.k_max =%e, "
      "limits.k_min = %e, and limits.k_max = %e", FPT.k_min, FPT.k_max,
      limits.k_min_cH0, limits.k_max_cH0);
    exit(1);
  }

  const double dlgk = log(10.) / (double)FPT.N_per_dec;
  const double tmp = log(FPT.k_min);
  {
    const int i = 0;
    k[i] = exp(tmp  + i*dlgk);
    P[i] = p_lin(k[i], 1.0);
  }
  #pragma omp parallel for
  for (int i=1; i<FPT.N; i++) 
  {
    k[i] = exp(tmp  + i*dlgk);
    P[i] = p_lin(k[i], 1.0);
  }
}

void get_FPT_bias(void) 
{
  static cosmopara C;
  
  if (recompute_cosmo3D(C)) 
  {
    update_cosmopara(&C);
    
    if (FPT.tab_AB == 0) 
    { // if table doesn't exit yet, create it
      FPT.tab_AB = (double**) malloc(sizeof(double*)*FPT.N_AB);
      for(int i=0; i<FPT.N_AB; i++)
      {
        FPT.tab_AB[i] = (double*) malloc(sizeof(double)*FPT.N);
      }

      if (FPT.k_min < limits.k_min_cH0) 
      {
        FPT.k_min = K_CH0(FPT.k_min);
        FPT.k_max = K_CH0(FPT.k_max);
      }
    }
    
    double k[FPT.N], Pin[FPT.N], Pout[FPT.N];
    
    FPT_input(k, Pin);
    
    Pd1d2(k, Pin, FPT.N, Pout);
    for (int i = 0; i < FPT.N; i++) 
    {
      FPT.tab_AB[0][i] = Pout[i]; // Pd1d2
    }
    
    Pd2d2(k, Pin, FPT.N, Pout);
    for (int i = 0; i < FPT.N; i++) 
    {
      FPT.tab_AB[1][i] = Pout[i]; // Pd2d2
    }
    
    Pd1s2(k, Pin, FPT.N, Pout);
    for (int i = 0; i < FPT.N; i++) 
    {
      FPT.tab_AB[2][i] = Pout[i]; // Pd1s2
    }
    
    Pd2s2(k, Pin, FPT.N, Pout);
    for (int i = 0; i < FPT.N; i++) 
    {
      FPT.tab_AB[3][i] = Pout[i]; // Pd2s2
    }
    
    Ps2s2(k, Pin, FPT.N, Pout);
    for (int i = 0; i < FPT.N; i++) 
    {
      FPT.tab_AB[4][i] = Pout[i]; // Pd2s2
    }
  }
}

double PT_d1d2(double k_coverH0) 
{ // interpolate FPT.tab_AB[0] - Pd1d2
  static cosmopara C;
  static double logkmin = 0., logkmax = 0., dlgk = 0.;
  
  if (recompute_cosmo3D(C)) 
  {
    get_FPT_bias();  // only call FASTPT if cosmology changed since last call
    update_cosmopara(&C);
    logkmin = log(FPT.k_min);
    logkmax = log(FPT.k_max);
    dlgk = log(10.) / (double) FPT.N_per_dec;
  }
  
  double lgk = log(k_coverH0);
  if (lgk < logkmin || lgk >= logkmax) 
  {
    return 0.;
  }

  return interpol(FPT.tab_AB[0], FPT.N, logkmin, logkmax, dlgk, lgk, 0., 0.);
}

double PT_d1d3(double k_coverH0) { // interpolate FPT.tab_AB[0] - Pd1d3
  static cosmopara C;
  static double logkmin = 0., logkmax = 0., dlgk = 0.;
 
  if (recompute_cosmo3D(C)) 
  {
    get_FPT_bias(); // only call FASTPT if cosmology changed since last call
    update_cosmopara(&C);
    logkmin = log(FPT.k_min);
    logkmax = log(FPT.k_max);
    dlgk = log(10.)/(double) FPT.N_per_dec;
  }
  
  double lgk = log(k_coverH0);
  if (lgk < logkmin || lgk >= logkmax) 
  {
    return 0.;
  }
  return interpol(tab_d1d3, FPT.N, logkmin, logkmax, dlgk,lgk, 0.,0.);
}

double PT_d2d2(double k_coverH0) 
{ // interpolate FPT.tab_AB[1]
  static cosmopara C;
  static double logkmin = 0., logkmax = 0., dlgk = 0.;
  
  if (recompute_cosmo3D(C)) 
  {
    get_FPT_bias(); // only call FASTPT if cosmology changed since last call
    update_cosmopara(&C);
    logkmin = log(FPT.k_min);
    logkmax = log(FPT.k_max);
    dlgk = log(10.) / (double)FPT.N_per_dec;
  }
  
  double lgk = log(k_coverH0);
  return interpol(FPT.tab_AB[1], FPT.N, logkmin, logkmax, dlgk, lgk, 0., 0.);
}

double PT_d1s2(double k_coverH0) 
{ // interpolate FPT.tab_AB[2]
  static cosmopara C;
  static double logkmin = 0., logkmax = 0., dlgk = 0.;
  
  if (recompute_cosmo3D(C)) 
  {
    get_FPT_bias(); // only call FASTPT if cosmology changed since last call
    update_cosmopara(&C);
    logkmin = log(FPT.k_min);
    logkmax = log(FPT.k_max);
    dlgk = log(10.) / (double)FPT.N_per_dec;
  }
  
  double lgk = log(k_coverH0);
  if (lgk < logkmin || lgk >= logkmax) 
  {
    return 0.;
  }
  return interpol(FPT.tab_AB[2], FPT.N, logkmin, logkmax, dlgk, lgk, 0., 0.);
}

double PT_d2s2(double k_coverH0) 
{ // interpolate FPT.tab_AB[3]
  static cosmopara C;
  static double logkmin = 0., logkmax = 0., dlgk = 0.;
  
  if (recompute_cosmo3D(C))
  {
    get_FPT_bias(); // only call FASTPT if cosmology changed since last call
    update_cosmopara(&C);
    logkmin = log(FPT.k_min);
    logkmax = log(FPT.k_max);
    dlgk = log(10.) / (double)FPT.N_per_dec;
  }
  
  double lgk = log(k_coverH0);
  if (lgk < logkmin || lgk >= logkmax) 
  {
    return 0.;
  }
  return interpol(FPT.tab_AB[3], FPT.N, logkmin, logkmax, dlgk, lgk, 0., 0.);
}

double PT_s2s2(double k_coverH0) 
{ // interpolate FPT.tab_AB[4]
  static cosmopara C;
  static double logkmin = 0., logkmax = 0., dlgk = 0.;
 
  if (recompute_cosmo3D(C)) 
  {
    get_FPT_bias();  // only call FASTPT if cosmology changed since last call
    update_cosmopara(&C);
    logkmin = log(FPT.k_min);
    logkmax = log(FPT.k_max);
    dlgk = log(10.) / (double)FPT.N_per_dec;
  }
  
  double lgk = log(k_coverH0);
  if (lgk < logkmin || lgk >= logkmax) 
  {
    return 0.;
  }
  return interpol(FPT.tab_AB[4], FPT.N, logkmin, logkmax, dlgk, lgk, 0., 0.);
}

double PT_sigma4(double k_coverH0 __attribute__((unused))) 
{
  static cosmopara C;
  
  if (recompute_cosmo3D(C)) 
  {
    get_FPT_bias(); // only call FASTPT if cosmology changed since last call
    update_cosmopara(&C);
  }
  
  return 0.0;
}

void get_FPT_IA(void) 
{
  static cosmopara C;
  if (recompute_cosmo3D(C)) 
  {
    update_cosmopara(&C);

    if (FPT.tab_IA == 0) 
    { // if table doesn't exit yet, create it
      FPT.tab_IA = (double**) malloc(sizeof(double*)*FPT.N_IA);
      for(int i=0; i<FPT.N_IA; i++)
      {
        FPT.tab_IA[i] = (double*) malloc(sizeof(double)*FPT.N);
      }

      if (FPT.k_min < limits.k_min_cH0) 
      {
        FPT.k_min = K_CH0(FPT.k_min);
        FPT.k_max = K_CH0(FPT.k_max);
      }
    }

    double k[FPT.N], Pin[FPT.N];
    double IA_tt_EE[FPT.N], IA_tt_BB[FPT.N];
    double IA_ta_dE1[FPT.N], IA_ta_dE2[FPT.N], IA_ta_0E0E[FPT.N],
        IA_ta_0B0B[FPT.N];
    double IA_mix_A[FPT.N], IA_mix_B[FPT.N], IA_mix_DEE[FPT.N],
        IA_mix_DBB[FPT.N];

    FPT_input(k, Pin);

    IA_tt(k, Pin, FPT.N, IA_tt_EE, IA_tt_BB);

    IA_ta(k, Pin, FPT.N, IA_ta_dE1, IA_ta_dE2, IA_ta_0E0E, IA_ta_0B0B);

    IA_mix(k, Pin, FPT.N, IA_mix_A, IA_mix_B, IA_mix_DEE, IA_mix_DBB);

    #pragma omp parallel for
    for (int i = 0; i < FPT.N; i++) 
    {
      FPT.tab_IA[0][i] = IA_tt_EE[i]; // tt_EE
      FPT.tab_IA[1][i] = IA_tt_BB[i]; // tt_BB
      FPT.tab_IA[2][i] = IA_ta_dE1[i];  // ta_dE1
      FPT.tab_IA[3][i] = IA_ta_dE2[i];  // ta_dE2
      FPT.tab_IA[4][i] = IA_ta_0E0E[i]; // ta_EE
      FPT.tab_IA[5][i] = IA_ta_0B0B[i]; // ta_BB
      FPT.tab_IA[6][i] = IA_mix_A[i];      // mix_A
      FPT.tab_IA[7][i] = IA_mix_B[i] * 4.; // mix_B
      FPT.tab_IA[8][i] = IA_mix_DEE[i];    // mix_D_EE
      FPT.tab_IA[9][i] = IA_mix_DBB[i];    // mix_D_BB
    }
  }
}

double TATT_II_EE(double k_coverH0, double a __attribute__((unused)), double C1,
double C2, double b_ta, double C1_2, double C2_2, double b_ta_2,
double growfac_a, double pdelta_ak) 
{
  double ta_EE = 0., tt_EE = 0., mix_EE = 0.;
  const double nla_EE = C1 * C1_2 * pdelta_ak;
  
  static cosmopara C;
  static double logkmin = 0., logkmax = 0., dlgk = 0.;

  // if TATT is specified, i.e. C2 !=0 (TT), or b_ta != 0 (TA)
  if (C2 != 0 || b_ta != 0)
  {
    // only call FASTPT if cosmology changed since last call
    if (recompute_cosmo3D(C))
    {
      get_FPT_IA();
      update_cosmopara(&C);
      logkmin = log(FPT.k_min);
      logkmax = log(FPT.k_max);
      dlgk = log(10.) / (double)FPT.N_per_dec;
    }
    double lgk = log(k_coverH0);
    // outside FASTPT interpolation range; return NLA contribution only
    if (lgk <= logkmin || lgk >= logkmax) {
      return nla_EE;
    }
    double g4 = growfac_a*growfac_a*growfac_a*growfac_a;

    const double P_ta_EE = g4 * interpol(FPT.tab_IA[4], FPT.N, logkmin, logkmax, dlgk, lgk,
                            0., 0.);
    const double P_ta_dE1 = g4 * interpol(FPT.tab_IA[2], FPT.N, logkmin, logkmax, dlgk, lgk,
                             0., 0.);
    const double P_ta_dE2 = g4 * interpol(FPT.tab_IA[3], FPT.N, logkmin, logkmax, dlgk, lgk,
                             0., 0.);

    ta_EE = C1 * C1_2 *
            (b_ta * b_ta_2 * P_ta_EE + (b_ta + b_ta_2) * (P_ta_dE1 + P_ta_dE2));

    if (C2) 
    {
      const double P_mix_EE = g4 * interpol(FPT.tab_IA[8], FPT.N, logkmin, logkmax, dlgk,
                               lgk, 0., 0.);
      const double P_mix_A = g4 * interpol(FPT.tab_IA[6], FPT.N, logkmin, logkmax, dlgk, lgk,
                              0., 0.);
      const double P_mix_B = g4 * interpol(FPT.tab_IA[7], FPT.N, logkmin, logkmax, dlgk, lgk,
                              0., 0.);

      mix_EE = (C1 * C2_2 + C1_2 * C2) * (P_mix_A + P_mix_B) +
               (C1 * b_ta * C2_2 + C1_2 * b_ta_2 * C2) * P_mix_EE;

      tt_EE = C2 * C2_2 * g4 *
          interpol(FPT.tab_IA[0], FPT.N, logkmin, logkmax, dlgk, lgk, 0., 0.);
    }
  }

  return nla_EE + ta_EE + mix_EE + tt_EE;
}

double TATT_II_BB(double k_coverH0, double a __attribute__((unused)), double C1,
double C2, double b_ta, double C1_2, double C2_2, double b_ta_2, double growfac_a) 
{
  double ta_BB = 0., tt_BB = 0., mix_BB = 0.;
  static cosmopara C;
  static double logkmin = 0., logkmax = 0., dlgk = 0.;
  
  if (C2 != 0 || b_ta != 0) 
  {
    if (recompute_cosmo3D(C)) 
    {
      get_FPT_IA(); // only call FASTPT if cosmology changed since last call
      update_cosmopara(&C);
      logkmin = log(FPT.k_min);
      logkmax = log(FPT.k_max);
      dlgk = log(10.) / (double)FPT.N_per_dec;
    }
    
    double lgk = log(k_coverH0);
    if (lgk < logkmin || lgk >= logkmax) 
    {
      return 0.0;
    }

    const double g4 = growfac_a*growfac_a*growfac_a*growfac_a;

    ta_BB = C1 * C1_2 * b_ta * b_ta_2 * g4 *
            interpol(FPT.tab_IA[5], FPT.N, logkmin, logkmax, dlgk, lgk, 0., 0.);

    if (C2) {
      mix_BB = (C1 * b_ta * C2_2 + C1_2 * b_ta_2 * C2) * g4 *
          interpol(FPT.tab_IA[9], FPT.N, logkmin, logkmax, dlgk, lgk, 0., 0.);

      tt_BB = C2 * C2_2 * g4 *
          interpol(FPT.tab_IA[1], FPT.N, logkmin, logkmax, dlgk, lgk, 0., 0.);
    }
  }

  return ta_BB + mix_BB + tt_BB;
}

double TATT_GI_E(double k_coverH0, double a __attribute__((unused)), double C1,
double C2, double b_ta, double growfac_a, double pdelta_ak) 
{
  double nla_GI = 0., ta_GI = 0., tt_GI = 0.;
  nla_GI = C1 * pdelta_ak;
  
  static cosmopara C;
  static double logkmin = 0., logkmax = 0., dlgk = 0.;
  
  if (C2 != 0 || b_ta != 0) 
  { // if TATT is specified, i.e. C2 !=0 (TT), or b_ta != 0 (TA)
    if (recompute_cosmo3D(C)) 
    {
      get_FPT_IA();
      update_cosmopara(&C); // only call FASTPT if cosmology changed since last call
      logkmin = log(FPT.k_min);
      logkmax = log(FPT.k_max);
      dlgk = log(10.) / (double)FPT.N_per_dec;
    }
    
    double lgk = log(k_coverH0);
    if (lgk < logkmin || lgk >= logkmax) 
    {
      return nla_GI;
    }
    
    const double g4 = growfac_a*growfac_a*growfac_a*growfac_a;

    double P_ta_dE1, P_ta_dE2;
    P_ta_dE1 = g4 * interpol(FPT.tab_IA[2], FPT.N, logkmin, logkmax, dlgk, lgk,
                             0., 0.);
    P_ta_dE2 = g4 * interpol(FPT.tab_IA[3], FPT.N, logkmin, logkmax, dlgk, lgk,
                             0., 0.);

    ta_GI = C1 * b_ta * (P_ta_dE1 + P_ta_dE2);

    if (C2) 
    {
      double P_mix_A, P_mix_B;
      P_mix_A = g4 * interpol(FPT.tab_IA[6], FPT.N, logkmin, logkmax, dlgk, lgk,
                              0., 0.);
      P_mix_B = g4 * interpol(FPT.tab_IA[7], FPT.N, logkmin, logkmax, dlgk, lgk,
                              0., 0.);

      tt_GI = C2 * (P_mix_A + P_mix_B);
    }
  }

  return nla_GI + ta_GI + tt_GI;
}
