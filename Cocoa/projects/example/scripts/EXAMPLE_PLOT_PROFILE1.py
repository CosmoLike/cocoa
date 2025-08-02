#!/usr/bin/env python3

import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# GENERAL PLOT OPTIONS
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

#fig, ax = plt.subplots(2, 4, sharex=False)
#ax = ax.flat
#fig.set_figheight(12.1)
#fig.set_figwidth(12.1)
#fig.subplots_adjust(hspace=0.225, wspace=0.275)

fig = plt.figure(figsize=(12.1, 12.1))

# Master grid: 2 rows (top and bottom blocks)
outer_gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1.1], hspace=0.225)

# Top row: 4 columns
gs_top = gridspec.GridSpecFromSubplotSpec(1,4,subplot_spec=outer_gs[0],wspace=0.275)
ax = [fig.add_subplot(gs_top[0, i]) for i in range(4)]

# Bottom row: 3 centered columns in full width
gs_bottom = gridspec.GridSpecFromSubplotSpec(1,3,subplot_spec=outer_gs[1], wspace=0.35)
ax += [fig.add_subplot(gs_bottom[0, i]) for i in range(3)]

# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------
# read minimum
minmodel=np.loadtxt(os.environ['ROOTDIR'] + "/projects/example/chains/EXAMPLE_EMUL_MIN1.txt")
# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------

root = os.environ['ROOTDIR'] + "/projects/example/chains/EXAMPLE_EMUL_PROFILE1."
params = ['logA', 'ns', 'omegabh2', 'omegach2', 'thetastar', 'tau', 'A_planck' ]
latex  = ["$\\log(10^{10} A_\\mathrm{s})$", "$n_\\mathrm{s}$", 
          "$100\\Omega_\\mathrm{b} h^2$", "$10\\Omega_\\mathrm{c} h^2$", 
          "$100\\theta_*$", "$\\tau_\\mathrm{reio}$", "$A_{\\rm Planck}$" ]
for i in range(7):
    data = np.loadtxt(root + params[i] + '.txt', comments="#",)
    if i == 2:
        data[:,0] = 100*data[:,0]
    if i == 3:
        data[:,0] = 10*data[:,0]
    ax[i].plot(data[:,0], data[:,1]-minmodel[-1], 
               marker='D',c='black', linestyle='solid', markersize=4,
               alpha=1.0,lw=1.0,
               label=params[i])
    ax[i].grid(True)
    ax[i].grid(True, which='minor', color='black',
           linestyle='--', linewidth=0.25, alpha=0.1)
    ax[i].minorticks_on()
    ax[i].tick_params(axis='both', which='major', labelsize=15)
    ax[i].tick_params(axis='both', which='minor', labelsize=15)
    ax[i].set_xlabel(latex[i],fontsize = 19)
    if i == 0 or i==4:
        ax[i].set_ylabel('$\\Delta \\chi^2$',fontsize = 19)
    ax[i].set_ylim(0.99*np.min(data[:,1]-minmodel[-1]),
                   1.01*np.max(data[:,1]-minmodel[-1])) 
    ax[i].set_xlim(data[0,0]-0.075*(data[-1,0]-data[0,0]),
                   data[-1,0]+0.075*(data[-1,0]-data[0,0]))

plt.subplots_adjust(bottom=0.25, left = 0.2)
plt.savefig(os.environ['ROOTDIR'] + "/projects/example/chains/EXAMPLE_PLOT_PROFILE1")