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
from matplotlib.patches import Polygon
from matplotlib.ticker import AutoMinorLocator
matplotlib.use('Agg')

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
{'Temperature': 15.0, 'density': 1.0E+03, 'path': "run_0.01/run_0.01_15.0_1.0E+03_1.0E-07_2.0_1.0E-02/"},
{'Temperature': 15.0, 'density': 1.0E+04, 'path': "run_0.01/run_0.01_15.0_1.0E+04_1.0E-07_15.0_1.0E-02/"},
{'Temperature': 15.0, 'density': 1.0E+05, 'path': "run_0.01/run_0.01_15.0_1.0E+05_1.0E-07_15.0_1.0E-02/"},
{'Temperature': 15.0, 'density': 1.0E+06, 'path': "run_0.01/run_0.01_15.0_1.0E+06_1.0E-07_15.0_1.0E-02/"},
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
     {'ntot':  7, 'i0': -0.2,  'cm': cm.Blues,  'name': 'mH2O'},
     {'ntot':  7, 'i0': -0.2,  'cm': cm.Blues,  'name': 'mHDO'},
     {'ntot':  7, 'i0': -0.2,  'cm': cm.Blues,  'name': 'mH2O2'},
     {'ntot':  7, 'i0': -0.2,  'cm': cm.Blues,  'name': 'mHDO2'},
     {'ntot':  7, 'i0': -0.2,  'cm': cm.Blues,  'name': 'mO2H'},
     {'ntot':  7, 'i0': -0.2,  'cm': cm.Blues,  'name': 'mOH'},
     {'ntot':  7, 'i0': -0.2,  'cm': cm.Blues,  'name': 'mO'},
     {'ntot':  7, 'i0': -0.2,  'cm': cm.Blues,  'name': 'mO2'},
     {'ntot':  7, 'i0': -0.2,  'cm': cm.Blues,  'name': 'mO3'},
     {'ntot': 11, 'i0': 9,  'cm': cm.OrRd ,  'name': 'mCO'},
     {'ntot': 11, 'i0': 9,  'cm': cm.OrRd ,  'name': 'mCO2'},
     {'ntot': 11, 'i0': 9,  'cm': cm.OrRd ,  'name': 'mH2CO'},
     {'ntot': 11, 'i0': 9,  'cm': cm.OrRd ,  'name': 'mHDCO'},
     {'ntot': 11, 'i0': 9,  'cm': cm.OrRd ,  'name': 'mD2CO'},
     {'ntot': 11, 'i0': 9,  'cm': cm.OrRd ,  'name': 'mCH3OH'},
     {'ntot': 11, 'i0': 9,  'cm': cm.OrRd ,  'name': 'mCH2DOH'},
     {'ntot': 11, 'i0': 9,  'cm': cm.OrRd ,  'name': 'mCHD2OH'},
     {'ntot': 11, 'i0': 9,  'cm': cm.OrRd ,  'name': 'mCD3OH'},
#    {'ntot': 11, 'i0': 9,  'cm': cm.OrRd ,  'name': 'mCH3OD'},
     {'ntot': 11, 'i0': 9,  'cm': cm.OrRd ,  'name': 'mHCOOH'},
#    {'ntot': 11, 'i0': 9,  'cm': cm.OrRd ,  'name': 'mDCOOH'},
     {'ntot': 11, 'i0': 9,  'cm': cm.OrRd ,  'name': 'mCH4'},
#    {'ntot': 11, 'i0': 9,  'cm': cm.OrRd ,  'name': 'mCH3D'},
     {'ntot':  2, 'i0':19,  'cm': cm.YlGnBu ,  'name': 'mN2'},
     {'ntot':  6, 'i0':20,  'cm': cm.cool ,  'name': 'mNO2'},
     {'ntot':  6, 'i0':20,  'cm': cm.cool ,  'name': 'mNH3'},
     {'ntot':  6, 'i0':20,  'cm': cm.cool ,  'name': 'mNH2D'},
#    {'ntot':  6, 'i0':20,  'cm': cm.cool ,  'name': 'mN2H2'},
#    {'ntot':  6, 'i0':20,  'cm': cm.cool ,  'name': 'mN2HD'},
     {'ntot':  6, 'i0':20,  'cm': cm.cool ,  'name': 'mHCN'},
     {'ntot':  6, 'i0':20,  'cm': cm.cool ,  'name': 'mHNC'},
     {'ntot':  6, 'i0':20,  'cm': cm.cool ,  'name': 'mHNO'},
#    {'ntot':  6, 'i0':20,  'cm': cm.cool ,  'name': 'mDNO'},
     {'ntot':  6, 'i0':20,  'cm': cm.cool ,  'name': 'mOCN'},
  ],
              ]
xRange = (0.01, 1.0)
yRange = (0, 100)
ncol_legend = 6
color_range = 0.9

tmp_list = [speciesPlot_s[0][i]['name'] for i in xrange(len(speciesPlot_s[0]))]
fig_filename = 'grid_layers_goodLooking_' + '_'.join(tmp_list) + '_New3.eps'

# ---- Initialization end ----

nMarkPoints = 10
linewidth = 2.0
nMarkEvery = 10
markersize = 10.0
markeredgewidth = 1.5
plotWidth = 0.9
plotheight = 0.8
figsize = (8,8.5)
x_begin = 0.08
y_begin = 0.05

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

nLayers_max = 150

for speciesPlot in speciesPlot_s:
  nSpeciesPlot = len(speciesPlot)
  for i in xrange(0,nSpeciesPlot):
    speciesPlot[i].update({'idx':[[] for j in xrange(n_files)]})
    for j in xrange(n_files):
      speciesPlot[i]['idx'][j] = \
        get_column(speciesPlot[i]['name'], data_files[j])

  fig = plt.figure(figsize=figsize)
  for j in xrange(n_files):
    dt = data_files[j]
    T = dt['Temperature']
    den = dt['density']
    ix = temperatures.index(T)
    iy = densities.index(den)
    ixy = iy * nx + ix + 1
    this_position = [x_begin+ix*pan_x, y_begin+(ny-iy-1)*pan_y, pan_x, pan_y]
    idx_den = idx_den_dic[dt['density']]
    idx_T = idx_T_dic[dt['Temperature']]
    xlabel = 'Fraction' if iy==ny-1 else ''
    ylabel = 'Layers' if ix==0 else ''
    ax = fig.add_axes(this_position,
      autoscalex_on=False, autoscaley_on=False, 
      xscale='linear', yscale='linear', xlim=xRange, ylim=yRange,
      xlabel=xlabel, ylabel=ylabel)
    if iy < (ny-1):
      ax.set_xticklabels([])
    if ix > 0:
      ax.set_yticklabels([])
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    x_layers = dt['numL']
    layer_acc = numpy.zeros(nLayers_max)
    for i in xrange(0,nSpeciesPlot):
      if speciesPlot[i]['idx'][j] == None:
        continue
      print j, i, speciesPlot[i]['name']
      y_layers = numpy.zeros(dt['nLayers'])
     #y_layers[0] = dt['bin'][dt['div'][0], speciesPlot[i]['idx'][j]]/dt['N_S']
     #y_layers[1:] = numpy.diff(dt['bin'][dt['div'], speciesPlot[i]['idx'][j]])/(dt['N_S'])
      y_layers[0] = dt['bin'][dt['div'][0], speciesPlot[i]['idx'][j]]/dt['N_in'][0]
      y_layers[1:] = numpy.diff(dt['bin'][dt['div'], speciesPlot[i]['idx'][j]])/(dt['N_in'][1:])
     #y_layers = smooth(y_layers)
      layer_prev = numpy.copy(layer_acc)
     #print y_layers, layer_acc[0:dt['nLayers']], len(y_layers), len(layer_acc[0:dt['nLayers']]), dt['nLayers']
      layer_acc[0:dt['nLayers']] += y_layers
      pxy = make_polygon_tran(x_layers, layer_prev[0:dt['nLayers']], layer_acc[0:dt['nLayers']])
     #print layer_prev[0:dt['nLayers']], layer_acc[0:dt['nLayers']]
      colormap = speciesPlot[i]['cm']
      thiscolor = colormap((i-speciesPlot[i]['i0'])*color_range/float(speciesPlot[i]['ntot']), 1)
      poly = Polygon(pxy, closed=True, facecolor=thiscolor, edgecolor='black', linewidth=0.01)
      ax.add_patch(poly)
      ax.plot([],[],marker='s', markerfacecolor=thiscolor, markeredgecolor='black', markersize=15, markeredgewidth=0.01, linestyle='None',
          label=getNameForPlot(speciesPlot[i]['name']))

    props = dict(boxstyle='round', facecolor='None')
    textstr = '(T = ${0:2.0f}$ K, '.format(T) + 'n$_{\\rm H}$ = ' + \
              scientific2latex('{0:7.1E}'.format(den)) + ' cm$^{{-3}}$)'
    #ax.text(0.05, 0.976, textstr, transform=ax.transAxes,
    #        verticalalignment='top', bbox=props)
    #if j == 0:
    #  lgd = ax.legend(loc='lower left', bbox_to_anchor=(0.0, 1.005), borderpad=0.01, fancybox=False, shadow=False, ncol=ncol_legend, numpoints=1)
    # #lgd.legendPatch.set_alpha(0.5)
    # #lgd.get_frame().set_facecolor('white')
    #  lgd.get_frame().set_boxstyle('round')
    # #lgd.get_frame().set_visible(False)

fig.savefig(os.path.join('.', path_save_fig, fig_filename))
plt.close()
print 'gv ' + os.path.join('.', path_save_fig, fig_filename) + ' &'
