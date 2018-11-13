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

data_files = [\
{'Temperature':  7.5, 'density': 1.0E+04, 'path': "run_0.01/run_0.01_7.5_1.0E+04_1.0E-07_15.0_1.0E-02/"},
{'Temperature': 10.0, 'density': 1.0E+04, 'path': "run_0.01/run_0.01_10.0_1.0E+04_1.0E-07_15.0_1.0E-02/"},
{'Temperature': 12.5, 'density': 1.0E+04, 'path': "run_0.01/run_0.01_12.5_1.0E+04_1.0E-07_15.0_1.0E-02/"},
{'Temperature': 15.0, 'density': 1.0E+04, 'path': "run_0.01/run_0.01_15.0_1.0E+04_1.0E-07_15.0_1.0E-02/"},
{'Temperature': 17.5, 'density': 1.0E+04, 'path': "run_0.01/run_0.01_17.5_1.0E+04_1.0E-07_15.0_1.0E-02/"},
#{'Temperature': 19.5, 'density': 1.0E+04, 'path': "run_0.01/run_0.01_19.5_1.0E+04_1.0E-07_15.0_1.0E-02/"},
{'Temperature':  7.5, 'density': 1.0E+05, 'path': "run_0.01/run_0.01_7.5_1.0E+05_1.0E-07_15.0_1.0E-02/"},
{'Temperature': 10.0, 'density': 1.0E+05, 'path': "run_0.01/run_0.01_10.0_1.0E+05_1.0E-07_15.0_1.0E-02/"},
{'Temperature': 12.5, 'density': 1.0E+05, 'path': "run_0.01/run_0.01_12.5_1.0E+05_1.0E-07_15.0_1.0E-02/"},
{'Temperature': 15.0, 'density': 1.0E+05, 'path': "run_0.01/run_0.01_15.0_1.0E+05_1.0E-07_15.0_1.0E-02/"},
{'Temperature': 17.5, 'density': 1.0E+05, 'path': "run_0.01/run_0.01_17.5_1.0E+05_1.0E-07_15.0_1.0E-02/"},
#{'Temperature': 19.5, 'density': 1.0E+05, 'path': "run_0.01/run_0.01_19.5_1.0E+05_1.0E-07_15.0_1.0E-02/"},
]

#speciesPlot_s = [
#  [
#     {'name': 'H'},
#     {'name': 'D'},
#     {'name': 'gH'},
#     {'name': 'gD'},
#  ],
#              ]
#t_check = [5E6]
#yRange = (8E-3, 2.5)
#ncol_legend = 2

speciesPlot_s = [
  [
     {'name': 'H2'},
     {'name': 'CO'},
     {'name': 'H2'},
     {'name': 'H'},
     {'name': 'H2'},
     {'name': 'D'},
  ],
              ]
yRange = (1E-7, 1E-1)
ncol_legend = 4
t_check = [1E1]
color_range = 0.9


nH_fixed_s = [1E4, 1E5]
xRange = (7.0, 20.0)

# ++++ Initialization begin ++++
pathAbsScript = os.path.dirname(os.path.abspath(sys.argv[0]))

execfile(os.path.join(pathAbsScript, 'config.py'))
execfile(os.path.join(pathAbsScript, 'multiRunPreProCommon.py'))
execfile(os.path.join(pathAbsScript, 'MySubsoutines.py'))

path_save_fig = 'run_0.01'
filename_bin = "evolution_moment__bin.dat"
filename_ascii = "evolution_moment__ascii.dat"

get_data()

# ---- Initialization end ----

linewidth = 2.0
markersize = 10.0
markeredgewidth = 1.5
markeredgecolor = 'black'
colormap = cm.prism

tmp_list = [speciesPlot_s[0][i]['name'] for i in xrange(len(speciesPlot_s[0]))]
fig_filename = 'vs_T_' + '_'.join(tmp_list) + '_at_{0:7.1E}'.format(t_check[0]).replace(' ','') + '.eps'

data_files.sort(key = lambda item: item['density'])
data_files.sort(key = lambda item: item['Temperature'])

n_files = len(data_files)

densities = sorted(list_unique(list(set([data_files[i]['density'] for i in range(n_files)]))))
temperatures = sorted(list_unique(list(set([data_files[i]['Temperature'] for i in range(n_files)]))))
idx_den_dic = dict([(densities[i], i) for i in range(len(densities))])
idx_T_dic = dict([(temperatures[i], i) for i in range(len(temperatures))])

nx = len(nH_fixed_s)
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
for i_nH in xrange(len(nH_fixed_s)):
  nH_fixed = nH_fixed_s[i_nH]
  T_plot = []
  idx_check = []
  idx_file = []
  for j in xrange(n_files):
    dt = data_files[j]
    nH = dt['density']
    if nH != nH_fixed:
      continue
    T = dt['Temperature']
    T_plot.append(T)
    idx_approx = getApproxIdxAsc_min(dt['bin'][:,0], t_check[0])
    t_ratio = t_check[0]/dt['bin'][idx_approx,0]
    print j, speciesPlot[0]['name']
    if t_ratio > 1.3 or t_ratio<0.8:
      print dt['path'], t_ratio
      raise Warning("The plotted curve might be wrong!")
    idx_check.append(idx_approx)
    idx_file.append(j)

  ix = i_nH
  iy = 0
  ixy = iy * nx + ix + 1
  idx_nH = idx_den_dic[nH_fixed]
  if ix==0:
    this_position = [x_begin+ix*pan_x, y_begin+(ny-iy-1)*pan_y, pan_x, pan_y]
  else:
    this_position = [x_begin+ix*pan_x+x_separation, y_begin+(ny-iy-1)*pan_y, pan_x, pan_y]
  xlabel = 'T (K)'
  ylabel = 'Abundance ratio' if i_nH==0 else ''
  ax = fig.add_axes(this_position,
    autoscalex_on=False, autoscaley_on=False, 
    xscale='linear', yscale='log', xlim=xRange, ylim=yRange,
    xlabel=xlabel, ylabel=ylabel)
  if iy < (ny-1):
    ax.set_xticklabels([])
  if ix > 0:
    ax.set_yticklabels([])
  for i in xrange(0,nSpeciesPlot-1,2):
    if speciesPlot[i]['idx'][j] == None:
      continue
    ratios_plot = [speciesPlot[i]['ratio'][idx_file[iCh]][idx_check[iCh]] for iCh in xrange(len(idx_check))]
    thiscolor = colormap((i+1)*color_range/float(nSpeciesPlot), 1)
    ax.plot(T_plot, ratios_plot,
      linestyle='-',
      color=thiscolor, linewidth=linewidth,
      marker=markers[(i/2)%nmarkers], markersize=markersize,
      markeredgewidth=markeredgewidth, markeredgecolor=markeredgecolor,
      markerfacecolor='None', markevery=1)
    if i_nH == 0:
      ax.plot([],[],
        linestyle='-',
        color=thiscolor,
        marker=markers[(i/2)%nmarkers], markersize=markersize,
        markeredgewidth=markeredgewidth, markeredgecolor=markeredgecolor,
        markerfacecolor='None', markevery=1,
        label='[' + getNameForPlot(speciesPlot[i+1]['name']) + \
          '/' + getNameForPlot(speciesPlot[i]['name']) + ']')
    
  props = dict(boxstyle='round', facecolor='white', alpha=0.5)
  textstr = 'n$_{\\rm H}$ = ' + \
              scientific2latex('{0:7.1E}'.format(nH_fixed)) + ' cm$^{{-3}}$'
  ax.text(0.05, 0.95, textstr, transform=ax.transAxes,
          verticalalignment='top', bbox=props)
  if i_nH == 0:
    lgd = ax.legend(loc='lower left', bbox_to_anchor=(0.0, 1.02), fancybox=False, shadow=False, ncol=ncol_legend)
    lgd.legendPatch.set_alpha(0.5)
    lgd.get_frame().set_facecolor('white')
    lgd.get_frame().set_boxstyle('round')

fig.savefig(os.path.join('.', path_save_fig, fig_filename))
plt.close()
print 'gv ' + os.path.join('.', path_save_fig, fig_filename) + ' &'


time_end = datetime.datetime.now()
print 'Finished at ', time_end
print 'A time span of ', time_end - time_start, ' has elapsed.'
