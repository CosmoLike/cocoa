"""Overlay the posteriors that four samplers found for the example emulator model.

The four chains, read from the chains/ folder of the project:
  EXAMPLE_EMUL_MCMC1      Cobaya's Metropolis-Hastings MCMC
                          (cobaya-run ./projects/example/EXAMPLE_EMUL_MCMC1.yaml);
  EXAMPLE_EMUL_NAUTILUS1  Nautilus nested sampling (EXAMPLE_EMUL_NAUTILUS1.py);
  EXAMPLE_EMUL_EMCEE1     emcee ensemble sampling (EXAMPLE_EMUL_EMCEE1.py);
  EXAMPLE_EMUL_POLY1      PolyChord nested sampling
                          (cobaya-run ./projects/example/EXAMPLE_EMUL_POLY1.yaml).
GetDist, the chain-analysis package, draws a triangle plot: the marginalized
posterior of each parameter on the diagonal and the 68% and 95% contours of each
pair of parameters below it. Agreement of the four sets of contours checks the
samplers against each other.

The plot also shows chi2v2 = -2 log posterior, built to mean the same in every
chain. Cobaya's chains store chi2 = -2 log likelihood and minuslogprior =
-log prior in separate columns, so chi2v2 = chi2 + 2*minuslogprior there; the
Nautilus and emcee scripts already wrote chi2 = -2 log posterior. The chains with
the new column are saved under hidden roots, chains/.VM_P1_TMP1 to
chains/.VM_P1_TMP4, which the plot reads; the original chain files do not change.

The four runs are not one configuration: the ACT DR6 lensing likelihood has
lens_only: False in EXAMPLE_EMUL_MCMC1.yaml and in the Nautilus script, and
lens_only: True in the emcee script and in EXAMPLE_EMUL_POLY1.yaml.

Output: chains/example_compare_chains.pdf in the project folder.

Run from the Cocoa folder after "source start_cocoa.sh" (the chain folder is
found through ROOTDIR, the Cocoa folder):

    python ./projects/example/scripts/EXAMPLE_PLOT_COMPARE_CHAINS.py
"""

import getdist.plots as gplot
from getdist import MCSamples
from getdist import loadMCSamples
import os
import matplotlib
import subprocess
import matplotlib.pyplot as plt
import numpy as np

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

# Parameters of the triangle plot; chi2v2 is the column added below.
parameter = [u'omegach2', u'logA', u'ns', u'omegabh2', u'tau', u'chi2v2']
chaindir  = os.environ['ROOTDIR'] + "/projects/example/chains/"

# GetDist settings. smooth_scale 0.25: Gaussian smoothing kernels a quarter of a
# standard deviation wide, in one and two dimensions. ignore_rows: the fraction of
# each chain dropped from the start as burn-in, 30% for the Metropolis-Hastings
# chain and 0 for the others (nested sampling has no burn-in, and the emcee
# script already dropped it). range_confidence: the tail probability GetDist uses
# to choose each axis range.
analysissettings={'smooth_scale_1D':0.25, 
                  'smooth_scale_2D':0.25,
                  'ignore_rows': u'0.3',
                  'range_confidence' : u'0.005'}

analysissettings2={'smooth_scale_1D':0.25,
                   'smooth_scale_2D':0.25,
                   'ignore_rows': u'0.0',
                   'range_confidence' : u'0.005'}

# Chain roots in chains/: Metropolis-Hastings, Nautilus, emcee, PolyChord.
root_chains = (
  'EXAMPLE_EMUL_MCMC1',
  'EXAMPLE_EMUL_NAUTILUS1',
  'EXAMPLE_EMUL_EMCEE1',
  'EXAMPLE_EMUL_POLY1',
)

# --------------------------------------------------------------------------------
# Metropolis-Hastings (Cobaya): chi2v2 = chi2 + 2*minuslogprior = -2 log
# posterior. getParams() gives access to the columns by name (p.chi2), and
# saveAsText writes the samples, new column included, under a hidden root.
samples=loadMCSamples(chaindir + root_chains[0],settings=analysissettings)
p = samples.getParams()
samples.addDerived(p.chi2+2*p.minuslogprior,name='chi2v2',label='{\\chi^2}')
samples.saveAsText(chaindir + '/.VM_P1_TMP1')
# --------------------------------------------------------------------------------
# Nautilus: its chi2 column is already -2 log posterior.
samples=loadMCSamples(chaindir+ root_chains[1], settings=analysissettings2)
p = samples.getParams()
samples.addDerived(p.chi2,name='chi2v2',label='{\\chi^2}')
samples.saveAsText(chaindir + '/.VM_P1_TMP2')
# --------------------------------------------------------------------------------
# emcee: its chi2 column is already -2 log posterior.
samples=loadMCSamples(chaindir+ root_chains[2],settings=analysissettings2)
p = samples.getParams()
samples.addDerived(p.chi2,name='chi2v2',label='{\\chi^2}')
samples.saveAsText(chaindir + '/.VM_P1_TMP3')
# --------------------------------------------------------------------------------
# PolyChord (Cobaya): the same columns as the Metropolis-Hastings chain.
samples=loadMCSamples(chaindir+ root_chains[3],settings=analysissettings2)
p = samples.getParams()
samples.addDerived(p.chi2+2*p.minuslogprior,name='chi2v2',label='{\\chi^2}')
samples.saveAsText(chaindir + '/.VM_P1_TMP4')
# --------------------------------------------------------------------------------

# GetDist plotter: analysis_settings applies to the chains it loads by root
# name. ignore_rows is 0 for all four hidden roots: the Metropolis-Hastings copy
# was saved after its burn-in rows were dropped. The g.settings lines set fonts,
# line widths and the legend frame.
g=gplot.getSubplotPlotter(chain_dir=chaindir,
                          analysis_settings=analysissettings2,
                          width_inch=10.5)
g.settings.axis_tick_x_rotation=65
g.settings.lw_contour=1.0
g.settings.legend_rect_border = False
g.settings.figure_legend_frame = False
g.settings.axes_fontsize = 15.0
g.settings.legend_fontsize = 15.5
g.settings.alpha_filled_add = 0.85
g.settings.lab_fontsize=15.5
g.legend_labels=False

# The lists give one entry per root, in the order of roots: filled contours for
# the Metropolis-Hastings and PolyChord chains, lines for the other two.
g.triangle_plot(
  params=parameter,
  roots=[chaindir + '/.VM_P1_TMP1',
         chaindir + '/.VM_P1_TMP2',
         chaindir + '/.VM_P1_TMP3',
         chaindir + '/.VM_P1_TMP4'],
  plot_3d_with_param=None,
  line_args=[{'lw': 1.0,'ls': 'solid', 'color':'lightcoral'},
              {'lw': 1.2,'ls': '--', 'color':'black'},
              {'lw': 2.1,'ls': 'dotted', 'color': 'maroon'},
              {'lw': 1.6,'ls': '-.', 'color': 'indigo'}
            ],
  contour_colors=['lightcoral','black','maroon', 'indigo'],
  contour_ls=['solid','--','dotted','-.'], 
  contour_lws=[1.0,1.2,2.1,1.6],
  filled=[True,False,False,True],
  shaded=False,
  # The numbers in the labels are typed in by hand from one set of runs; update
  # them after new runs. The two log(Z) values are not comparable: the Nautilus
  # script counts the prior density twice, which raises its log Z by -ln V = 15.0
  # for its prior box (EXAMPLE_EMUL_NAUTILUS1.py explains), and the two values
  # below differ by 14.98.
  legend_labels=[
    'MH, 4-walkers, $(R-1)_{\\rm median}$=0.015, $(R-1)_{\\rm std dev}$ = 0.18, burn-in=0.3',
    'Nautilus, $n_{\\rm live}=1024$, $\\log(Z)=-1172.64$, $n_{\\rm eval} \\sim 64,000$',
    'EMCEE $n_{\\rm walkers}=21$, $n_{\\rm eval} \\sim 80,000$',
    'PolyChord $n_{\\rm live}=512$, $n_{\\rm repeat}=3D$, $\\log(Z)=-1187.62 \\pm 0.20$',
  ],
  legend_loc=(0.3, 0.85))
g.export(os.path.join(chaindir,"example_compare_chains.pdf"))