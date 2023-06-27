import getdist.plots as gplot
from getdist import MCSamples
from getdist import loadMCSamples
import os
import matplotlib
import subprocess
import matplotlib.pyplot as plt
import numpy as np

#DELETE TMP FILES
subprocess.Popen("rm .VM_P4_TMP[0-9].*", shell=True, cwd=".")

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

# ADDING NEFF AS DERIVED PARAMETER
parameter = [u'w',u'wa',u'omegam']
chaindir=r'.'
root_chains = ('EXAMPLE_MCMC2','EXAMPLE_POLY2')

analysissettings={'smooth_scale_1D':0.45,'smooth_scale_2D':0.45,'ignore_rows': u'0.5',
'range_confidence' : u'0.025'}

analysissettings2={'smooth_scale_1D':0.45,'smooth_scale_2D':0.45,'ignore_rows': u'0.0',
'range_confidence' : u'0.025'}

# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------
samples=loadMCSamples('./' + root_chains[1],settings=analysissettings2) 
p = samples.getParams()
samples.saveAsText('.VM_P4_TMP2')
# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------
samples=loadMCSamples('./' + root_chains[0],settings=analysissettings)
p = samples.getParams()
# --------------------------------------------------------------------------------
samples.addDerived(np.exp(p.logA)/10.0,name='1e9As',label='{10^9\\rm A_s}',
range=[np.exp(samples.getLower('logA'))/10.0,np.exp(samples.getUpper('logA'))/10.0])
# --------------------------------------------------------------------------------
samples.addDerived(p.omegach2*10.0,name='10omegach2',label='{10 \Omega_c h^2}',
range=[samples.getLower('omegach2')*10.0,samples.getUpper('omegach2')*10.0])
# --------------------------------------------------------------------------------
samples.addDerived(p.omegabh2*100.0,name='100omegabh2',label='{100 \Omega_b h^2}',
range=[samples.getLower('omegabh2')*100.0,samples.getUpper('omegabh2')*100.0])
# --------------------------------------------------------------------------------
samples.addDerived(p.ns*10.0,name='10ns',label='{10 n_s}',
range=[samples.getLower('ns')*10.0,samples.getUpper('ns')*10.0])
# --------------------------------------------------------------------------------
samples.addDerived(p.w0pwa-p.w,name='wa',label='{w_a}')
# --------------------------------------------------------------------------------
samples.saveAsText('.VM_P4_TMP1')


#GET DIST PLOT SETUP
g=gplot.getSubplotPlotter(chain_dir=chaindir,
analysis_settings=analysissettings2,width_inch=5.9)
g.settings.lw_contour = 1.2
g.settings.legend_rect_border = False
g.settings.figure_legend_frame = False
g.settings.axes_fontsize = 13.5
g.settings.legend_fontsize = 15.5
g.settings.alpha_filled_add = 0.7
g.settings.lab_fontsize=15
g.legend_labels=False

param_3d = None
g.triangle_plot(['.VM_P4_TMP1','.VM_P4_TMP2'],parameter,
plot_3d_with_param=param_3d,line_args=[
{'lw': 1.0,'ls': 'solid', 'color':'royalblue'},
{'lw': 1.6,'ls': 'solid', 'color':'firebrick'},
{'lw': 1.9,'ls': 'dashed', 'color':'black'}],
contour_colors=['royalblue','firebrick','black'],
filled=True,shaded=False,
legend_labels=['Roman C - MCMC','Roman C - Poly'],legend_loc=(0.58,0.765),
imax_shaded=0)
g.export()

#DELETE TMP FILES
subprocess.Popen("rm .VM_P4_TMP[0-9].*", shell=True, cwd=".")