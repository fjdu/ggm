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
execfile(os.path.join(pathAbsScript, 'data_files_plot.py'))

get_data()

pp = PdfPages(os.path.join('.', path_save_fig, 'Ratios_grid.pdf'))

# ---- Initialization end ----

xRange = (1.1e0, 1e7)
ymin_, ymax_ = 1E-16, 4E2
nMarkPoints = 20
linewidth = 2.0
markersize = 10.0
markeredgewidth = 1.5
markeredgecolor = 'black'

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
plotWidth = 0.89
plotheight = 0.89
x_begin = 0.07
y_begin = 0.07
pan_x = plotWidth / nx
pan_y = plotheight / ny
figsize = (8*(1.2+(nx-1)/3.0), 6*(1.2+(ny-1)/3.0))
figsize = (8*(1.5), 6*(2.5))

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
      xscale='log', yscale='log', xlim=xRange, ylim=yRange,
      xlabel=xlabel, ylabel=ylabel)
    if iy < (ny-1):
      ax.set_xticklabels([])
    if ix > 0:
      ax.set_yticklabels([])
    for i in xrange(0,nSpeciesPlot-1,2):
      if speciesPlot[i]['idx'][j] == None:
        continue
      ax.plot(dt['bin'][:,0], speciesPlot[i]['ratio'][j],
        linestyle='-',
        color=cm.hot((idx_T+1)*0.6/float(nx), 1), linewidth=linewidth,
        marker=markers[(i/2)%nmarkers], markersize=markersize,
        markeredgewidth=markeredgewidth, markeredgecolor=markeredgecolor,
        markerfacecolor='None', markevery=nMarkEvery)
      if j == 0:
        ax.plot([],[], linestyle='None',
          marker=markers[(i/2)%nmarkers], markersize=markersize,
          markeredgewidth=markeredgewidth, markeredgecolor=markeredgecolor,
          markerfacecolor='None', markevery=nMarkEvery,
          label='[' + getNameForPlot(speciesPlot[i+1]['name']) + \
            ']/[' + getNameForPlot(speciesPlot[i]['name']) + ']')
    if nSpeciesPlot%2 == 1:
      if j==0:
        speciesPlot[nSpeciesPlot-1].update({'idx':[[] for jj in xrange(n_files)]})
      speciesPlot[nSpeciesPlot-1]['idx'][j] =\
        get_column(speciesPlot[nSpeciesPlot-1]['name'], data_files[j])
      idx_tmp = numpy.where\
       (numpy.logical_and(xRange[0]<=data_files[j]['bin'][:,0],
                          xRange[1]>=data_files[j]['bin'][:,0]))
      max_tmp = numpy.nanmax(dt['bin'][idx_tmp, speciesPlot[nSpeciesPlot-1]['idx'][j]])
      min_tmp = numpy.nanmin(dt['bin'][idx_tmp, speciesPlot[nSpeciesPlot-1]['idx'][j]])
      scale_fac = 0.9 * (numpy.sqrt(yRange[1]*yRange[0])-yRange[0]) / (max_tmp-min_tmp)
      ax.plot(dt['bin'][:,0], (dt['bin'][:,speciesPlot[nSpeciesPlot-1]['idx'][j]] - min_tmp)*scale_fac + yRange[0] * 2.0,
        linestyle=':',
        color=cm.hot((idx_T+1)*0.6/float(nx), 1), linewidth=linewidth,
        marker=markers[(nSpeciesPlot-1)%nmarkers], markersize=markersize,
        markeredgewidth=markeredgewidth, markeredgecolor=markeredgecolor,
        markerfacecolor='None', markevery=nMarkEvery)
      ax.plot([],[], linestyle='None',
        marker=markers[(nSpeciesPlot-1)%nmarkers], markersize=markersize,
        markeredgewidth=markeredgewidth, markeredgecolor=markeredgecolor,
        markerfacecolor='None', markevery=nMarkEvery,
        label='[' + getNameForPlot(speciesPlot[nSpeciesPlot-1]['name']) + ']')
      
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    textstr = '(${0:2.0f}$ K, '.format(T) + \
              scientific2latex('{0:7.1E}'.format(den)) + ' cm$^{{-3}}$)'
    ax.text(0.05, 0.95, textstr, transform=ax.transAxes,
            verticalalignment='top', bbox=props)
    if j == 0:
      lgd = ax.legend(loc='lower left', fancybox=False, shadow=False)
      lgd.legendPatch.set_alpha(0.5)
      lgd.get_frame().set_facecolor('wheat')
      lgd.get_frame().set_boxstyle('round')

  fig.savefig(pp, format='pdf')
  plt.close()

pp.close()

time_end = datetime.datetime.now()
print 'Finished at ', time_end
print 'A time span of ', time_end - time_start, ' has elapsed.'
