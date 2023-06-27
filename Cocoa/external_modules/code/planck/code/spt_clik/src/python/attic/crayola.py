# color correction yeah
import numpy as nm
import pyfits as pf

class const:
  k = 1.3806580e-23
  c =  2.9979246e+08
  h =  6.6260755e-34

def dBdT(nu, nu0):
  nu = nm.array(nu)
  x0 = nu0/56.78
  ex0 = nm.exp(x0)
  x = nu/56.78
  ex = nm.exp(x)
  res0 = (x0**4)*ex0/((ex0-1.)**2)
  res = (x**4)*ex/((ex-1.)**2)
  return res/res0

def dust_spectrum(nu, T_dust, beta_dust, nu0):
  nu = nm.array(nu)
  T_dust = nm.array(T_dust)
  beta_dust = nm.array(beta_dust)
  h_over_kT = 1./(T_dust * 20.836) # Assumes frequencies in GHz
  x0 = nu0 * h_over_kT
  ex0 = nm.exp(x0)
  res0 = (nu0**(3.0+beta_dust))/(ex0 - 1.) # Overall normalization will go away
  x = nu * h_over_kT
  ex = nm.exp(x);
  res = (nu**(3.0+beta_dust))/(ex - 1.)
  # Return dust emissivity normalized to one at nu0, in dT (CMB)
  #return res/res0/dBdT(nu,nu0)
  return res/res0

def dust_spectrum_function(T_dust,beta_dust,nu0):
  def fnc(nu):
    return dust_spectrum(nu,T_dust,beta_dust,nu0)
  return fnc

def non_thermal_spectrum(nu, alpha_non_thermal, nu0): 
  nu = nm.array(nu)
  alpha_non_thermal = nm.array(alpha_non_thermal)
  
  return nm.exp(alpha_non_thermal*nm.log(nu/nu0))/dBdT(nu,nu0)


def rimo_bandpass(path,cutoff=1e-10):
  rimo = pf.open(path)
  rep = {}
  for rr in rimo[1:]:
    extname = rr.header["EXTNAME"]
    if not extname.startswith("BANDPASS_"):
      continue
    name = extname[len("BANDPASS_"):]
    fq = rr.data.field("WAVENUMBER")*100*const.c/1e9
    tr = rr.data.field("TRANSMISSION")
    fq = fq[tr>cutoff]
    tr = tr[tr>cutoff]
    rep[name] = (fq[1:],tr[1:])
  return rep

def trapz(xs,ys):
  return nm.sum((xs[1:]-xs[:-1])*(ys[1:]+ys[:-1])/2.)

def color_correction(bandpass,spectrum_function,nu0):
  sp = spectrum_function(bandpass[0])
  s0 = spectrum_function(nu0)
  sp = sp/s0
  dbdt= dBdT(bandpass[0],nu0)
  conv = trapz(bandpass[0],bandpass[1]*dbdt)
  iras = trapz(bandpass[0],(bandpass[1]*(nu0/bandpass[0])))
  bd = trapz(bandpass[0],(bandpass[1]*sp))
  return bd/iras,conv/iras


def planck_blackbody(nu,T):
  nu = nu*1e9
  return 2*const.h*nu**3/const.c**2 * 1./(nm.exp(const.h*nu/(const.k*T))-1)

def dust_radiance(nu,T,beta,nu0):
  bd = planck_blackbody(nu,T)
  return bd * (nu/nu0)**beta

def dust_radiance_function(T,beta,nu0):
  def f(nu):
    return dust_radiance(nu,T,beta,nu0)
  return f