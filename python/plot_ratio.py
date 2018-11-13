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

pp = PdfPages(os.path.join('.', path_save_fig, 'ratio_plot.pdf'))

# ---- Initialization end ----

xRange = (1e1, 1e7)
#yRange = (1e-4, 5e1)
nMarkPoints = 10
linewidth = 3.0
markersize = 10.0
markeredgewidth = 1.5
markeredgecolor = 'black'
figsize = (8, 12)
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
      if i > 0 and speciesPlot[i]['idx'][j] != None:
        speciesPlot[i]['ratio'][j] = \
          data_files[j]['bin'][:,speciesPlot[i]['idx'][j]] / \
          data_files[j]['bin'][:,speciesPlot[0]['idx'][j]]

  if speciesPlot[0]['idx'][0] == None:
    continue

  maxAb_s = numpy.array([numpy.nanmax(data_files[j]['bin']\
    [:,speciesPlot[0]['idx'][j]]) for j in xrange(n_files)])
  maxRa_s = numpy.array([numpy.nanmax(speciesPlot[i]['ratio'][j][numpy.where(data_files[j]['bin'][:,speciesPlot[0]['idx'][j]] > min(maxAb_s)/1E5)])
     for j in xrange(n_files) for i in xrange(1,nSpeciesPlot)])

  yRange = (min(min(maxAb_s)/1e2, max(maxAb_s)/1e4), max(maxAb_s)*2.0)
  fig = plt.figure(figsize=figsize)
  ax = fig.add_axes([0.1,0.6,0.85,0.35],
    autoscalex_on=False, autoscaley_on=False,
    xscale='log', yscale='log', xlim=xRange, ylim=yRange,
    xticklabels=[],
    xlabel='', ylabel='[' + \
      getNameForPlot(speciesPlot[0]['name'] + ']'))
  for j in xrange(n_files):
    dt = data_files[j]
    idx_den = idx_den_dic[dt['density']]
    idx_T = idx_T_dic[dt['Temperature']]
    ax.plot(dt['bin'][:,0], dt['bin'][:,speciesPlot[0]['idx'][j]],
      linestyle=linestyle[idx_den%nlinestyle],
      color=colors[idx_T%ncolors], linewidth=3.0,
      label='T={0:4.1f} K, '.format(dt['Temperature']) + \
      'n=' + scientific2latex('{0:7.1e}'.format(dt['density'])) + ' cm$^{-3}$')
  lgd = ax.legend(loc='best', fancybox=True, shadow=True)
  lgd.legendPatch.set_alpha(0.5)

  yRange = (min(maxRa_s[numpy.where(numpy.isfinite(maxRa_s))])/1e2, 
            max(maxRa_s[numpy.where(numpy.isfinite(maxRa_s))])*2.0)
  print speciesPlot[0]['name'], yRange
  ax = fig.add_axes([0.1,0.1,0.85,0.5],
    autoscalex_on=False, autoscaley_on=False, 
    xscale='log', yscale='log', xlim=xRange, ylim=yRange,
    xlabel='Time (yr)', ylabel='Relative to [' + \
      getNameForPlot(speciesPlot[0]['name'] + ']'))
  for j in xrange(n_files):
    dt = data_files[j]
    idx_den = idx_den_dic[dt['density']]
    idx_T = idx_T_dic[dt['Temperature']]
    for i in xrange(1,nSpeciesPlot):
      if speciesPlot[i]['idx'][j] == None:
        continue
      ax.plot(dt['bin'][:,0], speciesPlot[i]['ratio'][j],
        linestyle=linestyle[idx_den%nlinestyle],
        color=colors[idx_T%ncolors], linewidth=linewidth,
        marker=markers[(i-1)%nmarkers], markersize=markersize,
        markeredgewidth=markeredgewidth, markeredgecolor=markeredgecolor,
        markerfacecolor='None', markevery=nMarkEvery)
      if j == 0:
        ax.plot([],[], linestyle='None',
          marker=markers[(i-1)%nmarkers], markersize=markersize,
          markeredgewidth=markeredgewidth, markeredgecolor=markeredgecolor,
          markerfacecolor='None', markevery=nMarkEvery,
          label=getNameForPlot(speciesPlot[i]['name']))
  lgd = ax.legend(loc='best', fancybox=True, shadow=True)
  lgd.legendPatch.set_alpha(0.5)

 #fig.savefig(os.path.join('.', path_save_fig,
 #    'Relative_to_' + speciesPlot[0]['name'] + '.png'))
  fig.savefig(pp, format='pdf')
  plt.close()

pp.close()

time_end = datetime.datetime.now()
print 'Finished at ', time_end
print 'A time span of ', time_end - time_start, ' has elapsed.'
