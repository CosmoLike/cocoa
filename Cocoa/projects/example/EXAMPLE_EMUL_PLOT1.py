import getdist.plots as gplot
from getdist import MCSamples
from getdist import loadMCSamples
import os
import matplotlib
import subprocess
import matplotlib.pyplot as plt
import numpy as np

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

parameter = [u'omegach2', u'logA', u'ns', u'omegabh2', u'tau']
chaindir=os.getcwd()

analysissettings={'smooth_scale_1D':0.3, 'smooth_scale_2D':0.3,'ignore_rows': u'0.3',
'range_confidence' : u'0.005'}

analysissettings2={'smooth_scale_1D':0.3,'smooth_scale_2D':0.3,'ignore_rows': u'0.0',
'range_confidence' : u'0.005'}

root_chains = (
  '/chains/EXAMPLE_EMUL_MCMC1',
  '/chains/EXAMPLE_EMUL_NAUTILUS1',
  '/chains/EXAMPLE_EMUL_EMCEE1',
)

# --------------------------------------------------------------------------------
samples=loadMCSamples(chaindir + root_chains[0],settings=analysissettings)
p = samples.getParams()
samples.saveAsText(chaindir + '/.VM_P1_TMP1')
# --------------------------------------------------------------------------------
samples=loadMCSamples(chaindir+ root_chains[1],settings=analysissettings2)
p = samples.getParams()
samples.saveAsText(chaindir + '/.VM_P1_TMP2')
# --------------------------------------------------------------------------------
samples=loadMCSamples(chaindir+ root_chains[1],settings=analysissettings2)
p = samples.getParams()
samples.saveAsText(chaindir + '/.VM_P1_TMP3')
# --------------------------------------------------------------------------------

#GET DIST PLOT SETUP
g=gplot.getSubplotPlotter(chain_dir=chaindir,
  analysis_settings=analysissettings2,width_inch=10.5)
g.settings.axis_tick_x_rotation=65
g.settings.lw_contour=1.0
g.settings.legend_rect_border = False
g.settings.figure_legend_frame = False
g.settings.axes_fontsize = 15.0
g.settings.legend_fontsize = 15.5
g.settings.alpha_filled_add = 0.85
g.settings.lab_fontsize=15.5
g.legend_labels=False

param_3d = None
g.triangle_plot([chaindir + '/.VM_P1_TMP1',
                 chaindir + '/.VM_P1_TMP2',
                 chaindir + '/.VM_P1_TMP3'],
parameter,
plot_3d_with_param=param_3d,line_args=[
{'lw': 1.0,'ls': 'solid', 'color':'lightcoral'},
{'lw': 1.1,'ls': 'solid', 'color':'black'},
{'lw': 1.6,'ls': 'solid', 'color': 'indigo'},
],
contour_colors=['lightcoral','black','indigo'],
contour_ls=['solid','solid','-.'], 
contour_lws=[1.0,1.1,1.6],
filled=[True,False,False],
shaded=False,
legend_labels=[
'MH, 4-walkers, $(R-1)_{\\rm median}$=0.015, $(R-1)_{\\rm std dev}$ = 0.18, burn-in=0.3',
'Nautilus, $n_{\\rm eval} \sim 64k$, $n_{\\rm live}=1024$',
'EMCEE $n_{\\rm eval} = 200k$, $n_{\\rm walkers}=28$, burn-in=0.3',
],
legend_loc=(0.3, 0.85))
g.export()