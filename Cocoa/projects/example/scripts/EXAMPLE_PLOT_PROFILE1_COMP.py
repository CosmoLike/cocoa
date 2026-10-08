#!/usr/bin/env python3

"""Compare the profiles of the two profile methods of the example project.

Each panel overlays two profiles of one parameter, both computed with the same
grid, covariance and minimum (scripts/EXAMPLE_EMUL_PROFILE1.sbatch and
scripts/EXAMPLE_EMUL_PROFILE1M2.sbatch):
  chains/EXAMPLE_EMUL_PROFILE1.<name>.txt    annealed emcee minimization at each
      grid point (EXAMPLE_EMUL_PROFILE1.py; labeled "Procoli", the method it
      reimplements), blue diamonds and a solid parabola;
  chains/EXAMPLE_EMUL_PROFILE1M2.<name>.txt  scipy's Nelder-Mead minimizer
      (EXAMPLE_EMUL_PROFILE_SCIPY1.py), red triangles and a dotted parabola.
Column 0 of a file holds the fixed value of the parameter and column 1 the chi2
= -2 log posterior minimized over the other parameters. Each profile is plotted
as Delta chi2 = chi2 - min(chi2) with a parabola fitted to its points; where the
two curves agree, the two minimizers found the same conditional minima.

Panels: logA, ns, 100 omegabh2 and 10 omegach2 on the top row; 100 theta_*, tau
and A_planck on the bottom row. omegabh2 and omegach2 are multiplied by 100 and
10 so that their tick labels stay short. The axis ranges of each panel follow
the Nelder-Mead profile, the one drawn last.

Input: the fourteen profile files, found through the ROOTDIR environment
variable (the Cocoa folder, exported by start_cocoa.sh).
Output: chains/EXAMPLE_PLOT_PROFILE1_COMP.pdf in the project folder (savefig
adds the .pdf extension of the savefig.format setting).

Run from the Cocoa folder after "source start_cocoa.sh":

    python ./projects/example/scripts/EXAMPLE_PLOT_PROFILE1_COMP.py
"""

import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import math
# Figure style: STIX fonts for text and mathematics, a light grid, and PDF output
# cropped to the drawn area.
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
matplotlib.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
matplotlib.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
matplotlib.rcParams['xtick.bottom'] = True
matplotlib.rcParams['xtick.top'] = False
matplotlib.rcParams['ytick.right'] = False
matplotlib.rcParams['axes.edgecolor'] = 'black'
matplotlib.rcParams['axes.linewidth'] = '1.0'
matplotlib.rcParams['axes.labelsize'] = 'medium'
matplotlib.rcParams['axes.grid'] = True
matplotlib.rcParams['grid.linewidth'] = '0.0'
matplotlib.rcParams['grid.alpha'] = '0.18'
matplotlib.rcParams['grid.color'] = 'lightgray'
matplotlib.rcParams['legend.labelspacing'] = 0.77
matplotlib.rcParams['savefig.bbox'] = 'tight'
matplotlib.rcParams['savefig.format'] = 'pdf'
# ------------------------------------------------------------------------------
fig = plt.figure(figsize=(15.1, 12.1))
# ------------------------------------------------------------------------------
# Layout: a master grid of two rows, four panels on top and three below, one per
# profiled parameter. ax is the list of the seven axes; each list comprehension
# adds one axis per column of its row.
master = gridspec.GridSpec(2, 
                           1, 
                           height_ratios=[1,1.1], 
                           hspace=0.225) # Master grid
top = gridspec.GridSpecFromSubplotSpec(1, 
                                       4,
                                       subplot_spec=master[0],
                                       wspace=0.275)
ax = [fig.add_subplot(top[0,i]) for i in range(4)]
bottom = gridspec.GridSpecFromSubplotSpec(1,
                                          3,
                                          subplot_spec=master[1],
                                          wspace=0.35)
ax += [fig.add_subplot(bottom[0,i]) for i in range(3)]
# ------------------------------------------------------------------------------
# Profile files are <root><name>.txt (annealed emcee) and <rootM2><name>.txt
# (Nelder-Mead), with name taken from params, in panel order.
root   = os.environ['ROOTDIR'] + "/projects/example/chains/EXAMPLE_EMUL_PROFILE1."
rootM2 = os.environ['ROOTDIR'] + "/projects/example/chains/EXAMPLE_EMUL_PROFILE1M2."
params = ['logA', 'ns', 'omegabh2', 'omegach2', 'thetastar', 'tau', 'A_planck' ]
latex  = ["$\\log(10^{10} A_\\mathrm{s})$", "$n_\\mathrm{s}$", 
          "$100\\Omega_\\mathrm{b} h^2$", "$10\\Omega_\\mathrm{c} h^2$", 
          "$100\\theta_*$", "$\\tau_\\mathrm{reio}$", "$A_{\\rm Planck}$" ]
# ------------------------------------------------------------------------------
for i in range(7):
# ------------------------------------------------------------------------------
    # Annealed emcee profile. Column 0: fixed parameter value; column 1: minimized
    # chi2. omegabh2 (i = 2) and omegach2 (i = 3) are rescaled by 100 and 10 to
    # match their labels.
    data = np.loadtxt(root + params[i] + '.txt', comments="#",)
    if i == 2:
        data[:,0] = 100*data[:,0]
    if i == 3:
        data[:,0] = 10*data[:,0]
    x  = data[:, 0]
    # Delta chi2 above the smallest chi2 of this profile's grid.
    y  = data[:, 1]-np.min(data[:,1])

    ax[i].plot(x, y, 
               marker='D',c='blue', linestyle='None', markersize=4,
               alpha=1.0,lw=1.0,
               label=None)
    
    # Least-squares parabola Delta chi2 = a x^2 + b x + c, the shape of a profile
    # when the posterior is Gaussian in this parameter.
    coeffs = np.polyfit(x, y, deg=2)
    xfit = np.linspace(np.min(x), np.max(x), 300)
    yfit = np.polyval(coeffs, xfit)
    ax[i].plot(xfit, 
               yfit,
               linestyle='solid', 
               color='blue', 
               lw=1.0, 
               alpha=0.8, 
               label='Procoli')
# ------------------------------------------------------------------------------
    # Nelder-Mead profile, with the same columns and rescaling. data, x and y now
    # hold this profile, so the axis limits below follow it.
    data = np.loadtxt(rootM2 + params[i] + '.txt', comments="#",)
    if i == 2:
        data[:,0] = 100*data[:,0]
    if i == 3:
        data[:,0] = 10*data[:,0]
    x  = data[:, 0]
    y  = data[:, 1]-np.min(data[:,1])

    ax[i].plot(x, y, 
               marker='v',c='red', linestyle='None', markersize=4,
               alpha=1.0,lw=1.0,
               label=None)
    
    # Parabola fit of the Nelder-Mead profile.
    coeffs = np.polyfit(x, y, deg=2)
    xfit = np.linspace(np.min(x), np.max(x), 300)
    yfit = np.polyval(coeffs, xfit)
    ax[i].plot(xfit, 
               yfit, 
               color='red', 
               linestyle=':', 
               lw=3.0, 
               alpha=0.8,
               label='SciPy (Nelder-Mead)')
# ------------------------------------------------------------------------------
    if i==4:
        ax[i].legend(fontsize=12, frameon=False)
# ------------------------------------------------------------------------------
    ax[i].grid(True)
    ax[i].grid(True, 
               which='minor', 
               color='black',
               linestyle='--', 
               linewidth=0.25, 
               alpha=0.1)
    ax[i].minorticks_on()
    ax[i].tick_params(axis='both', 
                      which='major', 
                      labelsize=15)
    ax[i].tick_params(axis='both', 
                      which='minor', 
                      labelsize=15)
    ax[i].set_xlabel(latex[i],fontsize = 19)
    if i == 0 or i==4:
        ax[i].set_ylabel('$\\Delta \\chi^2$',fontsize = 19)
    ax[i].set_ylim(np.min(y),np.max(y))
    ax[i].set_xlim(data[0,0]-0.075*(data[-1,0]-data[0,0]),
                   x[-1]+0.075*(x[-1]-x[0]))
# ------------------------------------------------------------------------------
plt.subplots_adjust(bottom=0.25, left = 0.2)
plt.savefig(os.environ['ROOTDIR'] + "/projects/example/chains/EXAMPLE_PLOT_PROFILE1_COMP")