import sys
import datetime
import math
import os
import shutil
import stat
import string
import numpy
import matplotlib
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
execfile(os.path.join(pathAbsScript, 'data_files_plot.py'))

get_data()

pp = PdfPages(os.path.join('.', path_save_fig, 'Abundances_plot.pdf'))

# ---- Initialization end ----

xRange = (1e1, 1e7)
#yRange = (1e-4, 5e1)
nMarkPoints = 10
linewidth = 3.0
markersize = 10.0
markeredgewidth = 1.5
markeredgecolor = 'black'
figsize = (8, 6)
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
  for i in xrange(nSpeciesPlot):
    speciesPlot[i].update({'idx':[[] for j in xrange(n_files)]})
    speciesPlot[i].update({'ratio':[[] for j in xrange(n_files)]})
    for j in xrange(n_files):
      speciesPlot[i]['idx'][j] = \
        get_column(speciesPlot[i]['name'], data_files[j])

  if speciesPlot[0]['idx'][0] == None:
    continue

  maxAb_s = numpy.array([numpy.nanmax(data_files[j]['bin']\
    [:,speciesPlot[i]['idx'][j]]) for j in xrange(n_files) for i in xrange(nSpeciesPlot)])

  yRange = (min(min(maxAb_s)/1e2, max(maxAb_s)/1e4), max(maxAb_s)*2.0)

  print speciesPlot[0]['name'], yRange
  fig = plt.figure(figsize=figsize)
  ax = fig.add_axes([0.1,0.1,0.85,0.85],
    autoscalex_on=False, autoscaley_on=False, 
    xscale='log', yscale='log', xlim=xRange, ylim=yRange,
    xlabel='Time (yr)', ylabel='Abundance')
  for j in xrange(n_files):
    dt = data_files[j]
    idx_den = idx_den_dic[dt['density']]
    idx_T = idx_T_dic[dt['Temperature']]
    for i in xrange(0,nSpeciesPlot):
      if speciesPlot[i]['idx'][j] == None:
        continue
      ax.plot(dt['bin'][:,0], dt['bin'][:,speciesPlot[i]['idx'][j]],
        linestyle=linestyle[idx_den%nlinestyle],
        color=colors[idx_T%ncolors], linewidth=linewidth,
        marker=markers[(i)%nmarkers], markersize=markersize,
        markeredgewidth=markeredgewidth, markeredgecolor=markeredgecolor,
        markerfacecolor='None', markevery=nMarkEvery)
      if j == 0:
        ax.plot([],[], linestyle='None',
          marker=markers[(i)%nmarkers], markersize=markersize,
          markeredgewidth=markeredgewidth, markeredgecolor=markeredgecolor,
          markerfacecolor='None', markevery=nMarkEvery,
          label=getNameForPlot(speciesPlot[i]['name']))
  lgd = ax.legend(loc='upper center', fancybox=False, shadow=True)
  lgd.legendPatch.set_alpha(0.5)

  plot_1 =[ax.plot([], [],
    linewidth=5, linestyle='-',
    color=colors[idx_T_dic[T]%ncolors]) for T in temperatures]
  legend_1 = plt.legend(plot_1,
             ['T={0:2.0f} K'.format(T) for T in temperatures],
             loc='upper left', fancybox=False, shadow=True, title='Colour code')
  legend_1.legendPatch.set_alpha(0.5)
  plot_2 = [ax.plot([], [],
    linewidth=3, color=colors[-1],
    linestyle=linestyle[idx_den_dic[den]]) for den in densities]
  legend_2 = plt.legend(plot_2,
             ['n=' + scientific2latex('{0:7.1E}'.format(den)) + ' cm$^{{-3}}$' for den in densities],
             loc='upper right', fancybox=False, shadow=True, title='Line style')
  legend_2.legendPatch.set_alpha(0.5)
  ax.add_artist(legend_1)
  ax.add_artist(lgd)

  fig.savefig(pp, format='pdf')
  plt.close()

pp.close()

time_end = datetime.datetime.now()
print 'Finished at ', time_end
print 'A time span of ', time_end - time_start, ' has elapsed.'
