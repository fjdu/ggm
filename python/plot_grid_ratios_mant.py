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

#speciesPlot_s = [
#  [
#     {'name': 'H2'},
#     {'name': 'mCH3OH'},
#     {'name': 'H2'},
#     {'name': 'mCH2DOH'},
#     {'name': 'H2'},
#     {'name': 'mCHD2OH'},
#     {'name': 'H2'},
#     {'name': 'mCD3OH'},
#     {'name': 'H2'},
#     {'name': 'mCH3OD'},
#  ],
#              ]
#xRange = (9.1e3, 9.9e6)
#ncol_legend = 2

#speciesPlot_s = [
#  [
#     {'name': 'mCH3OH'},
#     {'name': 'mCH2DOH'},
#     {'name': 'mCH3OH'},
#     {'name': 'mCHD2OH'},
#     {'name': 'mCH3OH'},
#     {'name': 'mCD3OH'},
#     {'name': 'mCH3OH'},
#     {'name': 'mCH3OD'},
#  ],
#              ]
#xRange = (9.1e3, 9.9e6)
#ncol_legend = 2

#speciesPlot_s = [
#  [
#     {'name': 'mH2O'},
#     {'name': 'mHDO'},
#     {'name': 'mH2CO'},
#     {'name': 'mHDCO'},
#     {'name': 'mCH3OH'},
#     {'name': 'mCH2DOH'},
#  ],
#              ]
#xRange = (9.1e3, 9.9e6)
#ncol_legend = 2

#speciesPlot_s = [
#  [
#     {'name': 'gH'},
#     {'name': 'gD'},
#     {'name': 'H'},
#     {'name': 'D'},
#  ],
#              ]
#xRange = (9.1e3, 9.9e6)
#ncol_legend = 2

#speciesPlot_s = [
#  [
#     {'name': 'H2'},
#     {'name': 'H'},
#     {'name': 'H2'},
#     {'name': 'D'},
#  ],
#              ]
#xRange = (1.1e-3, 9.9e6)
#yRange = (1e-8, 3e-1)
#ncol_legend = 2
#color_range = 0.7

#speciesPlot_s = [
#  [
#     {'name': 'H2'},
#     {'name': 'gH'},
#     {'name': 'H2'},
#     {'name': 'gD'},
#  ],
#              ]
#xRange = (1.1e-4, 9.9e6)
#yRange = (1e-27, 3e-14)
#ncol_legend = 2

#speciesPlot_s = [
#  [
#     {'name': 'H2'},
#     {'name': 'gO'},
#     {'name': 'H2'},
#     {'name': 'gO2'},
#     {'name': 'H2'},
#     {'name': 'gO3'},
#  ],
#              ]
#xRange = (1.1e-4, 9.9e6)
#yRange = (1e-15, 1e-7)
#ncol_legend = 2

#speciesPlot_s = [
#  [
#     {'name': 'H2'},
#     {'name': 'O'},
#     {'name': 'H2'},
#     {'name': 'O2'},
#     {'name': 'H2'},
#     {'name': 'CO'},
#  ],
#              ]
#xRange = (1.1e-4, 9.9e6)
#yRange = (1e-6, 6e-4)
#ncol_legend = 2

speciesPlot_s = [
  [
     {'name': 'H2'},
     {'name': 'mH2O'},
     {'name': 'H2'},
     {'name': 'mCO'},
     {'name': 'H2'},
     {'name': 'mCO2'},
     {'name': 'H2'},
     {'name': 'mCH3OH'},
     {'name': 'H2'},
     {'name': 'mH2CO'},
  ],
              ]
xRange = (3E2, 1E7)
yRange = (1e-7, 1E-3)
ncol_legend = 3
color_range = 0.99

#speciesPlot_s = [
#  [
#     {'name': 'mH2O'},
#     {'name': 'mCO'},
#     {'name': 'mH2O'},
#     {'name': 'mCO2'},
#     {'name': 'mH2O'},
#     {'name': 'mCH3OH'},
#     {'name': 'mH2O'},
#     {'name': 'mH2CO'},
#  ],
#              ]
#xRange = (3E2, 1E7)
#yRange = (1e-4, 3E1)
#ncol_legend = 2
#color_range = 0.99

speciesPlot_s = [
  [
     {'name': 'totalMant'},
  ],
              ]
xRange = (3E2, 1E7)
yRange = (0, 110)
ncol_legend = 2
color_range = 0.9

tmp_list = [speciesPlot_s[0][i]['name'] for i in xrange(len(speciesPlot_s[0]))]
fig_filename = 'grid_' + '_'.join(tmp_list) + '.eps'

# ---- Initialization end ----

ymin_, ymax_ = 1E-16, 4E2
nMarkPoints = 10
linewidth = 2.0
markersize = 10.0
markeredgewidth = 1.5
markeredgecolor = 'black'
plotWidth = 0.87
plotheight = 0.80
figsize = (8,9)
x_begin = 0.1
y_begin = 0.06
colormap = cm.jet

nMarkEvery = numpy.sum(numpy.logical_and\
  (xRange[0]<=data_files[0]['bin'][:,0],\
   xRange[1]>=data_files[0]['bin'][:,0])) \
  / nMarkPoints

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
  for i in xrange(0,nSpeciesPlot,1):
    speciesPlot[i].update({'idx':[[] for j in xrange(n_files)]})
    speciesPlot[i].update({'ratio':[[] for j in xrange(n_files)]})
    for j in xrange(n_files):
      speciesPlot[i]['idx'][j] = \
        get_column(speciesPlot[i]['name'], data_files[j])
      if speciesPlot[i]['idx'][j] != None:
        speciesPlot[i]['ratio'][j] = \
          data_files[j]['bin'][:,speciesPlot[i]['idx'][j]] / \
          (4.0*math.pi*((GrainRadius_s[0]*1e2)**2) * SitesDensity)

  maxAb_s = numpy.array([numpy.nanmax(data_files[j]['bin']\
    [numpy.where\
       (numpy.logical_and(xRange[0]<=data_files[j]['bin'][:,0],
                          xRange[1]>=data_files[j]['bin'][:,0])),
     speciesPlot[i]['idx'][j]]) for i in xrange(0,nSpeciesPlot-1,2) for j in xrange(n_files)])
  maxRa_s_y = numpy.array([numpy.nanmax(speciesPlot[i]['ratio'][j][numpy.where\
       (numpy.logical_and(\
        numpy.logical_and(xRange[0]<=data_files[j]['bin'][:,0],
                          xRange[1]>=data_files[j]['bin'][:,0]), \
        data_files[j]['bin'][:,speciesPlot[i]['idx'][j]] > min(maxAb_s)/1E4))])
     for j in xrange(n_files) for i in xrange(0,nSpeciesPlot-1,2)])
  minRa_s_y = numpy.array([numpy.nanmin(speciesPlot[i]['ratio'][j][numpy.where\
       (numpy.logical_and(\
        numpy.logical_and(xRange[0]<=data_files[j]['bin'][:,0],
                          xRange[1]>=data_files[j]['bin'][:,0]), \
        data_files[j]['bin'][:,speciesPlot[i]['idx'][j]] > min(maxAb_s)/1E4))])
     for j in xrange(n_files) for i in xrange(0,nSpeciesPlot-1,2)])
  if 'yRange' not in locals():
    yRange = (max(ymin_,min(minRa_s_y)/2, min(maxRa_s_y)/5), min(ymax_,max(maxRa_s_y)*2.0))

  fig = plt.figure(figsize=figsize)
  for j in xrange(n_files):
    print j, speciesPlot[0]['name']
    dt = data_files[j]
    T = dt['Temperature']
    den = dt['density']
    ix = temperatures.index(T)
    iy = densities.index(den)
    ixy = iy * nx + ix + 1
    this_position = [x_begin+ix*pan_x, y_begin+(ny-iy-1)*pan_y, pan_x, pan_y]
    idx_den = idx_den_dic[dt['density']]
    idx_T = idx_T_dic[dt['Temperature']]
    xlabel = 'Time (yr)' if iy==ny-1 else ''
    ylabel = 'Abundance ratio' if ix==0 else ''
    ax = fig.add_axes(this_position,
      autoscalex_on=False, autoscaley_on=False, 
      xscale='log', yscale='linear', xlim=xRange, ylim=yRange,
      xlabel=xlabel, ylabel=ylabel)
    if iy < (ny-1):
      ax.set_xticklabels([])
    if ix > 0:
      ax.set_yticklabels([])
    for i in xrange(0,nSpeciesPlot,1):
      if speciesPlot[i]['idx'][j] == None:
        continue
      thiscolor = colormap((i+1)*color_range/float(nSpeciesPlot), 1)
      ax.plot(dt['bin'][:,0], speciesPlot[i]['ratio'][j],
        linestyle='-',
        color=thiscolor, linewidth=linewidth,
        marker=markers[(i/2)%nmarkers], markersize=markersize,
        markeredgewidth=markeredgewidth, markeredgecolor=markeredgecolor,
        markerfacecolor='None', markevery=nMarkEvery)
      if j == 0:
        ax.plot([],[],
          linestyle='-',
          color=thiscolor, linewidth=linewidth,
          marker=markers[(i/2)%nmarkers], markersize=markersize,
          markeredgewidth=markeredgewidth, markeredgecolor=markeredgecolor,
          markerfacecolor='None', markevery=nMarkEvery,
          label='Number of mantle layers')
      
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
