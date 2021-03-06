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

def get_data():
  for dt in data_files:
    ShapeEvolMo = getNumOfColRowInFile(os.path.join('.',
                  dt['path'], filename_ascii))
    ShapeEvolMo = (ShapeEvolMo[1]-1, ShapeEvolMo[0])

    fEvolMoBin = open(os.path.join('.',
        dt['path'], filename_bin), 'rb')
    dt.update({'bin': numpy.fromfile(file=fEvolMoBin,
        dtype=numpy.dtype('f8')).reshape(ShapeEvolMo)})
    dt['bin'][:,1:] *= funcRatioGrainToH()
    fEvolMoBin.close()
    f_tmp = open(os.path.join('.', dt['path'], filename_ascii))
    dt.update({'nameSpecies': f_tmp.readline().split()})
    f_tmp.close()

def get_column(speName, dt):
  for i in xrange(len(dt['nameSpecies'])):
    if dt['nameSpecies'][i] == speName:
      return i
      break
  return None

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
get_data()

speciesPlot_s = [
  [
     {'name': 'H2'},
     {'name': 'gCH3OH'},
     {'name': 'H2'},
     {'name': 'gCH2DOH'},
     {'name': 'H2'},
     {'name': 'gCHD2OH'},
     {'name': 'H2'},
     {'name': 'gCD3OH'},
     {'name': 'H2'},
     {'name': 'gCH3OD'},
  ],
              ]
fig_filename = 'vs_n_gCH3OHgCH2DOHgCHD2OHgCD3OHgCH3OD.eps'

# ---- Initialization end ----

t_check = [1E6]
xRange = (9.1E2, 1.2E6)
yRange = (1E-12, 2E-6)
ymin_, ymax_ = 1E-16, 4E2
nMarkPoints = 10
linewidth = 2.0
markersize = 10.0
markeredgewidth = 1.5
markeredgecolor = 'black'
ncol_legend = 3

data_files.sort(key = lambda item: item['density'])
data_files.sort(key = lambda item: item['Temperature'])

n_files = len(data_files)

T_fixed_s = [10.0,15]
densities = sorted(list_unique(list(set([data_files[i]['density'] for i in range(n_files)]))))
temperatures = sorted(list_unique(list(set([data_files[i]['Temperature'] for i in range(n_files)]))))
idx_den_dic = dict([(densities[i], i) for i in range(len(densities))])
idx_T_dic = dict([(temperatures[i], i) for i in range(len(temperatures))])

nx = len(T_fixed_s)
ny = 1
plotWidth = 0.85
plotheight = 0.70
x_begin = 0.1
y_begin = 0.1
x_separation = 0.03
pan_x = plotWidth / nx
pan_y = plotheight / ny
figsize = (8,6)

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
  for i in xrange(0,nSpeciesPlot-1,2):
    speciesPlot[i].update({'idx':[[] for j in xrange(n_files)]})
    speciesPlot[i+1].update({'idx':[[] for j in xrange(n_files)]})
    speciesPlot[i].update({'ratio':[[] for j in xrange(n_files)]})
    for j in xrange(n_files):
      speciesPlot[i]['idx'][j] = \
        get_column(speciesPlot[i]['name'], data_files[j])
      speciesPlot[i+1]['idx'][j] = \
        get_column(speciesPlot[i+1]['name'], data_files[j])
      if speciesPlot[i]['idx'][j] != None and speciesPlot[i+1]['idx'][j] != None:
        speciesPlot[i]['ratio'][j] = \
          data_files[j]['bin'][:,speciesPlot[i+1]['idx'][j]] / \
          data_files[j]['bin'][:,speciesPlot[i]['idx'][j]]

fig = plt.figure(figsize=figsize)
for i_T in xrange(len(T_fixed_s)):
  T_fixed = T_fixed_s[i_T]
  den_plot = []
  idx_check = []
  idx_file = []
  for j in xrange(n_files):
    print j, speciesPlot[0]['name']
    dt = data_files[j]
    T = dt['Temperature']
    if T != T_fixed:
      continue
    den = dt['density']
    den_plot.append(den)
    idx_check.append(getApproxIdxAsc_min(dt['bin'][:,0], t_check[0]))
    idx_file.append(j)

  ix = i_T
  iy = 0
  ixy = iy * nx + ix + 1
  idx_T = idx_T_dic[T_fixed]
  if ix==0:
    this_position = [x_begin+ix*pan_x, y_begin+(ny-iy-1)*pan_y, pan_x, pan_y]
  else:
    this_position = [x_begin+ix*pan_x+x_separation, y_begin+(ny-iy-1)*pan_y, pan_x, pan_y]
  xlabel = 'n$_{\\rm H}$ (cm$^{-3}$)'
  ylabel = 'Abundance ratio' if i_T==0 else ''
  ax = fig.add_axes(this_position,
    autoscalex_on=False, autoscaley_on=False, 
    xscale='log', yscale='log', xlim=xRange, ylim=yRange,
    xlabel=xlabel, ylabel=ylabel)
  if iy < (ny-1):
    ax.set_xticklabels([])
  if ix > 0:
    ax.set_yticklabels([])
  for i in xrange(0,nSpeciesPlot-1,2):
    if speciesPlot[i]['idx'][j] == None:
      continue
    ratios_plot = [speciesPlot[i]['ratio'][idx_file[iCh]][idx_check[iCh]] for iCh in xrange(len(idx_check))]
    ax.plot(den_plot, ratios_plot,
      linestyle='-',
      color=cm.hot((idx_T+1)*0.6/float(nx), 1), linewidth=linewidth,
      marker=markers[(i/2)%nmarkers], markersize=markersize,
      markeredgewidth=markeredgewidth, markeredgecolor=markeredgecolor,
      markerfacecolor='None', markevery=1)
    if i_T == 0:
      ax.plot([],[], linestyle='None',
        marker=markers[(i/2)%nmarkers], markersize=markersize,
        markeredgewidth=markeredgewidth, markeredgecolor=markeredgecolor,
        markerfacecolor='None', markevery=1,
        label='[' + getNameForPlot(speciesPlot[i+1]['name']) + \
          '/' + getNameForPlot(speciesPlot[i]['name']) + ']')
    
  props = dict(boxstyle='round', facecolor='white', alpha=0.5)
  textstr = 'T = ${0:2.0f}$ K'.format(T_fixed)
  ax.text(0.05, 0.95, textstr, transform=ax.transAxes,
          verticalalignment='top', bbox=props)
  if i_T == 0:
    lgd = ax.legend(loc='lower left', bbox_to_anchor=(0.0, 1.02), fancybox=False, shadow=False, ncol=ncol_legend)
    lgd.legendPatch.set_alpha(0.5)
    lgd.get_frame().set_facecolor('white')
    lgd.get_frame().set_boxstyle('round')

  fig.savefig(os.path.join('.', path_save_fig, fig_filename))
  plt.close()


time_end = datetime.datetime.now()
print 'Finished at ', time_end
print 'A time span of ', time_end - time_start, ' has elapsed.'
