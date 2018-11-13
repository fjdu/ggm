import sys
import datetime
import math
import os
import shutil
import stat
import string
import numpy
import matplotlib
import matplotlib.cm as cm
from matplotlib.ticker import MultipleLocator
#import matplotlib.ticker as tk
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages

time_start = datetime.datetime.now()
print 'Started at ', time_start


#matplotlib.use('cairo.ps')
#matplotlib.rc('text', usetex=True)
import matplotlib.pyplot as plt

# ++++ Initialization begin ++++
pathAbsScript = os.path.dirname(os.path.abspath(sys.argv[0]))

execfile(os.path.join(pathAbsScript, 'config.py'))
execfile(os.path.join(pathAbsScript, 'multiRunPreProCommon.py'))
execfile(os.path.join(pathAbsScript, 'MySubsoutines.py'))

data_files = [\
{'Temperature': 10.0, 'density': 1.0E+03, 'path': "run_0.01/run_0.01_10.0_1.0E+03_1.0E-07_2.0_1.0E-02/"},
{'Temperature': 10.0, 'density': 1.0E+04, 'path': "run_0.01/run_0.01_10.0_1.0E+04_1.0E-07_15.0_1.0E-02/"},
{'Temperature': 10.0, 'density': 1.0E+05, 'path': "run_0.01/run_0.01_10.0_1.0E+05_1.0E-07_15.0_1.0E-02/"},
{'Temperature': 10.0, 'density': 1.0E+06, 'path': "run_0.01/run_0.01_10.0_1.0E+06_1.0E-07_15.0_1.0E-02/"},
#{'Temperature': 15.0, 'density': 1.0E+03, 'path': "run_0.01/run_0.01_15.0_1.0E+03_1.0E-07_2.0_1.0E-02/"},
#{'Temperature': 15.0, 'density': 1.0E+04, 'path': "run_0.01/run_0.01_15.0_1.0E+04_1.0E-07_15.0_1.0E-02/"},
#{'Temperature': 15.0, 'density': 1.0E+05, 'path': "run_0.01/run_0.01_15.0_1.0E+05_1.0E-07_15.0_1.0E-02/"},
#{'Temperature': 15.0, 'density': 1.0E+06, 'path': "run_0.01/run_0.01_15.0_1.0E+06_1.0E-07_15.0_1.0E-02/"},
#{'Temperature': 20.0, 'density': 1.0E+03, 'path': "run_0.01/run_0.01_20.0_1.0E+03_1.0E-07_2.0_1.0E-02/"},
#{'Temperature': 20.0, 'density': 1.0E+04, 'path': "run_0.01/run_0.01_20.0_1.0E+04_1.0E-07_15.0_1.0E-02/"},
#{'Temperature': 20.0, 'density': 1.0E+05, 'path': "run_0.01/run_0.01_20.0_1.0E+05_1.0E-07_15.0_1.0E-02/"},
#{'Temperature': 20.0, 'density': 1.0E+06, 'path': "run_0.01/run_0.01_20.0_1.0E+06_1.0E-07_15.0_1.0E-02/"},
#{'Temperature': 30.0, 'density': 1.0E+03, 'path': "run_0.01/run_0.01_30.0_1.0E+03_1.0E-07_2.0_1.0E-02/"},
#{'Temperature': 30.0, 'density': 1.0E+04, 'path': "run_0.01/run_0.01_30.0_1.0E+04_1.0E-07_15.0_1.0E-02/"},
#{'Temperature': 30.0, 'density': 1.0E+05, 'path': "run_0.01/run_0.01_30.0_1.0E+05_1.0E-07_15.0_1.0E-02/"},
#{'Temperature': 30.0, 'density': 1.0E+06, 'path': "run_0.01/run_0.01_30.0_1.0E+06_1.0E-07_15.0_1.0E-02/"},
]
path_save_fig = 'run_0.01'
filename_bin = "evolution_moment__bin.dat"
filename_ascii = "evolution_moment__ascii.dat"
get_data(useNumPerGrain=True)
divideIntoLayers()

speciesPlot_s = [
  [
     {'name': 'mH2O'},
     {'name': 'mCO'},
     {'name': 'mCO2'},
     {'name': 'mCH3OH'},
     {'name': 'mH2CO'},
     {'name': 'mCH4'},
     {'name': 'mNH3'},
     {'name': 'mN2'},
     {'name': 'mO2'},
     {'name': 'mO3'},
  ],
              ]
xRange = (0, 100)
yRange = (6e-3, 2.0)
ncol_legend = 3
color_range = 0.99

#speciesPlot_s = [
#  [
#     {'name': 'mH2O'},
#     {'name': 'mCO'},
#     {'name': 'mCO2'},
#     {'name': 'mCH3OH'},
#     {'name': 'mH2CO'},
#  ],
#              ]
#xRange = (0, 100)
#yRange = (6e-3, 2.0)
#ncol_legend = 4

#speciesPlot_s = [
#  [
#     {'name': 'mCH4'},
#     {'name': 'mNH3'},
#     {'name': 'mN2'},
#     {'name': 'mO2'},
#     {'name': 'mO3'},
#  ],
#              ]
#xRange = (0, 100)
#yRange = (6e-3, 2.0)
#ncol_legend = 4

tmp_list = [speciesPlot_s[0][i]['name'] for i in xrange(len(speciesPlot_s[0]))]
fig_filename = 'new_grid_layers_' + '_'.join(tmp_list) + '.eps'

# ---- Initialization end ----

nMarkPoints = 10
linewidth = 2.0
nMarkEvery = 10
markersize = 10.0
markeredgewidth = 1.5
markeredgecolor = 'black'
plotWidth = 0.86
plotheight = 0.80
figsize = (6,8)
x_begin = 0.11
y_begin = 0.06
#colormap = cm.jet
execfile(os.path.join(pathAbsScript, 'make_colors_auto_my.py'))
colormap = my_custom_colormap

data_files.sort(key = lambda item: item['density'])
data_files.sort(key = lambda item: item['Temperature'])

n_files = len(data_files)

densities = sorted(list_unique(list(set([data_files[i]['density'] for i in range(n_files)]))))
temperatures = sorted(list_unique(list(set([data_files[i]['Temperature'] for i in range(n_files)]))))
idx_den_dic = dict([(densities[i], i) for i in range(len(densities))])
idx_T_dic = dict([(temperatures[i], i) for i in range(len(temperatures))])

nx = len(temperatures)
ny = len(densities)
pan_x = plotWidth / nx
pan_y = plotheight / ny

linestyle_density = dict([(linestyle[l], linestyle_dens[l]) for l in xrange(len(linestyle))])
color_temperature = dict([(c, RGB2UV(c)) for c in colors])
if len(idx_den_dic) > 1:
  linestyle.sort(key=lambda x: linestyle_density[x])
if len(idx_T_dic) > 1:
  colors.sort(key=lambda x: color_temperature[x][0])
if len(idx_den_dic) < len(linestyle):
  linestyle = linestyle[len(linestyle)-len(idx_den_dic):]
if len(idx_T_dic) < len(colors):
  colors = colors[len(colors)-len(idx_T_dic):]


for speciesPlot in speciesPlot_s:
  nSpeciesPlot = len(speciesPlot)
  for i in xrange(0,nSpeciesPlot):
    speciesPlot[i].update({'idx':[[] for j in xrange(n_files)]})
    for j in xrange(n_files):
      speciesPlot[i]['idx'][j] = \
        get_column(speciesPlot[i]['name'], data_files[j])

  fig = plt.figure(figsize=figsize)
  for j in xrange(n_files):
    print j, speciesPlot[i]['name']
    dt = data_files[j]
    T = dt['Temperature']
    den = dt['density']
    ix = temperatures.index(T)
    iy = densities.index(den)
    ixy = iy * nx + ix + 1
    this_position = [x_begin+ix*pan_x, y_begin+(ny-iy-1)*pan_y, pan_x, pan_y]
    idx_den = idx_den_dic[dt['density']]
    idx_T = idx_T_dic[dt['Temperature']]
    xlabel = 'Layer number' if iy==ny-1 else ''
    ylabel = 'Fraction' if ix==0 else ''
    ax = fig.add_axes(this_position,
      autoscalex_on=False, autoscaley_on=False, 
      xscale='linear', yscale='log', xlim=xRange, ylim=yRange,
      xlabel=xlabel, ylabel=ylabel)
    ax.xaxis.set_minor_locator(MultipleLocator(5))
    if iy < (ny-1):
      ax.set_xticklabels([])
    if ix > 0:
      ax.set_yticklabels([])
    for i in xrange(0,nSpeciesPlot):
      if speciesPlot[i]['idx'][j] == None:
        continue
      x_layers = dt['numL']
      y_layers = numpy.zeros(dt['nLayers'])
     #y_layers[0] = dt['bin'][dt['div'][0], speciesPlot[i]['idx'][j]]/dt['N_S']
     ##for ii in xrange(dt['nLayers']):
     ##  print ii, x_layers[ii], y_layers[ii]
     ##print len(numpy.diff(dt['bin'][dt['div'], speciesPlot[i]['idx'][j]])/(dt['N_S']))
     #y_layers[1:] = numpy.diff(dt['bin'][dt['div'], speciesPlot[i]['idx'][j]])/(dt['N_S'])
     #y_layers = smooth(y_layers)
      y_layers[0] = dt['bin'][dt['div'][0], speciesPlot[i]['idx'][j]]/dt['N_in'][0]
      y_layers[1:] = numpy.diff(dt['bin'][dt['div'], speciesPlot[i]['idx'][j]])/(dt['N_in'][1:])
      #thiscolor = colormap((i)*color_range/float(nSpeciesPlot), 1)
      thiscolor = colormap(i, 1)
      ax.plot(x_layers, y_layers,
        linestyle='-',
        color=thiscolor, linewidth=linewidth,
        marker=markers[(i)%nmarkers], markersize=markersize,
        markeredgewidth=markeredgewidth, markeredgecolor=markeredgecolor,
        markerfacecolor='None', markevery=nMarkEvery)
      if j == 0:
        ax.plot([],[],
          linestyle='-',
          color=thiscolor, linewidth=linewidth,
          marker=markers[(i)%nmarkers], markersize=markersize,
          markeredgewidth=markeredgewidth, markeredgecolor=markeredgecolor,
          markerfacecolor='None', markevery=nMarkEvery,
          label=getNameForPlot(speciesPlot[i]['name']))

    props = dict(boxstyle='round', facecolor='white', alpha=0.5)
    textstr = '(${0:2.0f}$ K, '.format(T) + \
              scientific2latex('{0:7.1E}'.format(den)) + ' cm$^{{-3}}$)'
    ax.text(0.05, 0.95, textstr, transform=ax.transAxes,
            verticalalignment='top', bbox=props)
    if j == 0:
      lgd = ax.legend(loc='lower left', bbox_to_anchor=(0.0, 1.02), fancybox=False, shadow=False, ncol=ncol_legend)
      lgd.legendPatch.set_alpha(0.5)
      lgd.get_frame().set_facecolor('white')
      lgd.get_frame().set_boxstyle('round')

fig.savefig(os.path.join('.', path_save_fig, fig_filename))
plt.close()
print 'gv ' + os.path.join('.', path_save_fig, fig_filename) + ' &'
