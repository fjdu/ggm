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

T_plot = 20.0
T_plot_str = '20.0'

folder_type1 = 'HME'
folder_type2 = 'HME'

data_files = [\
{'Temperature': T_plot, 'density': 1.0E+03,
  'path':        folder_type1 + "_0.01/"    + T_plot_str + "_1.0E+03_1.0E-07_2.0_1.0E-02/",
  'comparewith': folder_type2 + "_0.01/" + T_plot_str + "_1.0E+03_1.0E-07_2.0_1.0E-02/"},
{'Temperature': T_plot, 'density': 1.0E+04,
  'path':        folder_type1 + "_0.01/"    + T_plot_str + "_1.0E+04_1.0E-07_15.0_1.0E-02/",
  'comparewith': folder_type2 + "_0.01/" + T_plot_str + "_1.0E+04_1.0E-07_15.0_1.0E-02/"},
{'Temperature': T_plot, 'density': 1.0E+05,
  'path':        folder_type1 + "_0.01/"    + T_plot_str + "_1.0E+05_1.0E-07_15.0_1.0E-02/",
  'comparewith': folder_type2 + "_0.01/" + T_plot_str + "_1.0E+05_1.0E-07_15.0_1.0E-02/"},
{'Temperature': T_plot, 'density': 1.0E+06,
  'path':        folder_type1 + "_0.01/"    + T_plot_str + "_1.0E+06_1.0E-07_15.0_1.0E-02/",
  'comparewith': folder_type2 + "_0.01/" + T_plot_str + "_1.0E+06_1.0E-07_15.0_1.0E-02/"},
]
path_save_fig = ''
filename_bin = "evolution_moment__bin.dat"
filename_ascii = "evolution_moment__ascii.dat"
get_data()
get_data_cmp()
figname_prefix = 'T_' + T_plot_str + '_' + folder_type1 + '_'

#data_files = [\
#{'Temperature': 10.0, 'density': 1.0E+03,
#  'path':        "HME_0.01/"    + "10.0_1.0E+03_1.0E-07_2.0_1.0E-02/",
#  'comparewith': "HME_0.01/" + "10.0_1.0E+03_1.0E-07_2.0_1.0E-02/"},
#{'Temperature': 10.0, 'density': 1.0E+04,
#  'path':        "HME_0.01/"    + "10.0_1.0E+04_1.0E-07_15.0_1.0E-02/",
#  'comparewith': "HME_0.01/" + "10.0_1.0E+04_1.0E-07_15.0_1.0E-02/"},
#{'Temperature': 10.0, 'density': 1.0E+05,
#  'path':        "HME_0.01/"    + "10.0_1.0E+05_1.0E-07_15.0_1.0E-02/",
#  'comparewith': "HME_0.01/" + "10.0_1.0E+05_1.0E-07_15.0_1.0E-02/"},
#{'Temperature': 10.0, 'density': 1.0E+06,
#  'path':        "HME_0.01/"    + "10.0_1.0E+06_1.0E-07_15.0_1.0E-02/",
#  'comparewith': "HME_0.01/" + "10.0_1.0E+06_1.0E-07_15.0_1.0E-02/"},
#]
#path_save_fig = ''
#filename_bin = "evolution_moment__bin.dat"
#filename_ascii = "evolution_moment__ascii.dat"
#get_data()
#get_data_cmp()
#figname_prefix = 'HME_all_n_'

#data_files = [\
#{'Temperature': 10.0, 'density': 1.0E+03,
#  'path':        "RE_0.01/"    + "10.0_1.0E+03_1.0E-07_2.0_1.0E-02/",
#  'comparewith': "RE_0.01/" + "10.0_1.0E+03_1.0E-07_2.0_1.0E-02/"},
#{'Temperature': 10.0, 'density': 1.0E+04,
#  'path':        "RE_0.01/"    + "10.0_1.0E+04_1.0E-07_15.0_1.0E-02/",
#  'comparewith': "RE_0.01/" + "10.0_1.0E+04_1.0E-07_15.0_1.0E-02/"},
#{'Temperature': 10.0, 'density': 1.0E+05,
#  'path':        "RE_0.01/"    + "10.0_1.0E+05_1.0E-07_15.0_1.0E-02/",
#  'comparewith': "RE_0.01/" + "10.0_1.0E+05_1.0E-07_15.0_1.0E-02/"},
#{'Temperature': 10.0, 'density': 1.0E+06,
#  'path':        "RE_0.01/"    + "10.0_1.0E+06_1.0E-07_15.0_1.0E-02/",
#  'comparewith': "RE_0.01/" + "10.0_1.0E+06_1.0E-07_15.0_1.0E-02/"},
#]
#path_save_fig = ''
#filename_bin = "evolution_moment__bin.dat"
#filename_ascii = "evolution_moment__ascii.dat"
#get_data()
#get_data_cmp()
#figname_prefix = 'RE_all_n_'

#data_files = [\
#{'Temperature': 10.0, 'density': 1.0E+03,
#  'path':        "highHbind_0.01/"    + "10.0_1.0E+03_1.0E-07_2.0_1.0E-02/",
#  'comparewith': "highHbind_0.01/" + "10.0_1.0E+03_1.0E-07_2.0_1.0E-02/"},
#{'Temperature': 10.0, 'density': 1.0E+04,
#  'path':        "highHbind_0.01/"    + "10.0_1.0E+04_1.0E-07_15.0_1.0E-02/",
#  'comparewith': "highHbind_0.01/" + "10.0_1.0E+04_1.0E-07_15.0_1.0E-02/"},
#{'Temperature': 10.0, 'density': 1.0E+05,
#  'path':        "highHbind_0.01/"    + "10.0_1.0E+05_1.0E-07_15.0_1.0E-02/",
#  'comparewith': "highHbind_0.01/" + "10.0_1.0E+05_1.0E-07_15.0_1.0E-02/"},
#{'Temperature': 10.0, 'density': 1.0E+06,
#  'path':        "highHbind_0.01/"    + "10.0_1.0E+06_1.0E-07_15.0_1.0E-02/",
#  'comparewith': "highHbind_0.01/" + "10.0_1.0E+06_1.0E-07_15.0_1.0E-02/"},
#]
#path_save_fig = ''
#filename_bin = "evolution_moment__bin.dat"
#filename_ascii = "evolution_moment__ascii.dat"
#get_data()
#get_data_cmp()
#figname_prefix = 'highHbind_all_n_'

speciesPlot_s = [
  [
     {'name': 'mH2O'},
     {'name': 'mHDO'},
     {'name': 'mHDO'},
     {'name': 'mD2O'},
  ],
              ]
xRange = (1e2, 1e7)
yRange = (1e-4, 1e0)
ncol_legend = 1
nMarkPoints = 10
linewidth = 1.0
markersize = 10.0
markeredgewidth = 1.5
markeredgecolor = 'black'
plotWidth = 0.87
plotheight = 0.87
figsize = (6,6)
x_begin = 0.12
y_begin = 0.1
execfile(os.path.join(pathAbsScript, 'make_colors_auto_my.py'))
colormap = my_colors

# ---- Config end ----

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
#ny = len(densities)
ny = 1
pan_x = plotWidth / nx
pan_y = plotheight / ny

linestyle_density = dict([(linestyle[l], linestyle_dens[l]) for l in xrange(len(linestyle))])
color_temperature = dict([(c, RGB2UV(c)) for c in colors])
#if len(idx_den_dic) > 1:
#  linestyle.sort(key=lambda x: linestyle_density[x])
#if len(idx_T_dic) > 1:
#  colors.sort(key=lambda x: color_temperature[x][0])
#if len(idx_den_dic) < len(linestyle):
#  linestyle = linestyle[len(linestyle)-len(idx_den_dic):]
#if len(idx_T_dic) < len(colors):
#  colors = colors[len(colors)-len(idx_T_dic):]

for speciesPlot in speciesPlot_s:
  nSpeciesPlot = len(speciesPlot)
  for i in xrange(0,nSpeciesPlot-1,2):
    speciesPlot[i].update({'idx':[[] for j in xrange(n_files)]})
    speciesPlot[i+1].update({'idx':[[] for j in xrange(n_files)]})
    speciesPlot[i].update({'ratio':[[] for j in xrange(n_files)]})
    speciesPlot[i].update({'finalRatio':[[] for j in xrange(n_files)]})
    for j in xrange(n_files):
      speciesPlot[i]['idx'][j] = \
        get_column(speciesPlot[i]['name'], data_files[j])
      speciesPlot[i+1]['idx'][j] = \
        get_column(speciesPlot[i+1]['name'], data_files[j])
      if speciesPlot[i]['idx'][j] != None and speciesPlot[i+1]['idx'][j] != None:
        speciesPlot[i]['ratio'][j] = \
          data_files[j]['bin'][:,speciesPlot[i+1]['idx'][j]] / \
          data_files[j]['bin'][:,speciesPlot[i]['idx'][j]]
        if 'multiply' in speciesPlot[i]:
          speciesPlot[i]['ratio'][j] *= speciesPlot[i]['multiply']
    speciesPlot[i].update({'idx_cmp':[[] for j in xrange(n_files)]})
    speciesPlot[i+1].update({'idx_cmp':[[] for j in xrange(n_files)]})
    speciesPlot[i].update({'ratio_cmp':[[] for j in xrange(n_files)]})
    for j in xrange(n_files):
      speciesPlot[i]['idx_cmp'][j] = \
        data_files[j]['nameSpecies_cmp'].index(speciesPlot[i]['name'])
      speciesPlot[i+1]['idx_cmp'][j] = \
        data_files[j]['nameSpecies_cmp'].index(speciesPlot[i+1]['name'])
      if speciesPlot[i]['idx_cmp'][j] != None and speciesPlot[i+1]['idx'][j] != None:
        speciesPlot[i]['ratio_cmp'][j] = \
          data_files[j]['bin_cmp'][:,speciesPlot[i+1]['idx'][j]] / \
          data_files[j]['bin_cmp'][:,speciesPlot[i]['idx'][j]]
        if 'multiply' in speciesPlot[i]:
          speciesPlot[i]['ratio_cmp'][j] *= speciesPlot[i]['multiply']

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
  this_position = [x_begin, y_begin, pan_x, pan_y]
  xlabel = 'Time (yr)'
  ylabel = 'Abundance ratio'
  ax = fig.add_axes(this_position,
    autoscalex_on=False, autoscaley_on=False, 
    xscale='log', yscale='log', xlim=xRange, ylim=yRange,
    xlabel=xlabel, ylabel=ylabel)
  #if iy < (ny-1):
  #  ax.set_xticklabels([])
  #if ix > 0:
  #  ax.set_yticklabels([])
  for j in xrange(n_files):
    dt = data_files[j]
    T = dt['Temperature']
    den = dt['density']
    iy = densities.index(den)
    thiscolor = colormap[iy]
    for i in xrange(0,nSpeciesPlot-1,2):
      if len(speciesPlot[i]['ratio'][j]) > 2000:
        speciesPlot[i]['finalRatio'][j] = \
          speciesPlot[i]['ratio'][j][2000] / speciesPlot[i]['ratio_cmp'][j][2000]
        print '{0:.1f}, {1:.1e}, {2:16s}, {3:.1e}, {4:.1e}, {5:.1e}'.\
          format(T, den, speciesPlot[i+1]['name'] + '/' + speciesPlot[i]['name'], \
          speciesPlot[i]['ratio'][j][2000], speciesPlot[i]['ratio_cmp'][j][2000], \
          speciesPlot[i]['ratio'][j][2000] / speciesPlot[i]['ratio_cmp'][j][2000])
      if speciesPlot[i]['idx'][j] == None:
        continue
      thislinestyle = linestyle[(i/2)%nlinestyle]
      ax.plot(dt['bin'][:,0], speciesPlot[i]['ratio'][j],
        linestyle=thislinestyle,
        color=thiscolor, linewidth=linewidth)
        #marker=markers[(i/2)%nmarkers], markersize=markersize,
        #markeredgewidth=markeredgewidth, markeredgecolor=markeredgecolor,
        #markerfacecolor='None', markevery=nMarkEvery)
      #ax.plot(dt['bin_cmp'][:,0], speciesPlot[i]['ratio_cmp'][j],
      #  linestyle=':',
      #  color=thiscolor, linewidth=linewidth*2,
      #  marker=markers[(i/2)%nmarkers], markersize=markersize,
      #  markeredgewidth=markeredgewidth, markeredgecolor=markeredgecolor,
      #  markerfacecolor='None', markevery=nMarkEvery)
      if 'drawVLineAt' in speciesPlot[i]:
        maxVal = max(speciesPlot[i]['ratio'][j])
        idx_drawVL = getApproxIdxAsc_min(speciesPlot[i]['ratio'][j], maxVal*speciesPlot[i]['drawVLineAt'])
        ax.axvline(x=dt['bin'][idx_drawVL,0], ymin=yRange[0], ymax=yRange[1],
          linestyle='--',
          color=thiscolor, linewidth=linewidth/2)
      if j == 0:
        label_str = getNameForPlot(speciesPlot[i+1]['name']) + \
            '/' + getNameForPlot(speciesPlot[i]['name'])
            #+ '  n$_{\\rm H}$ = ' +\
            #scientific2latex('{0:7.1E}'.format(den)) + ' cm$^{-3}$'
        ax.plot([],[],
          linestyle=linestyle[(i/2)%nlinestyle],
          color=thiscolor, linewidth=linewidth,
          label=label_str)
      
    if j == 0:
      lgd = ax.legend(loc='upper left', bbox_to_anchor=(0.01, 0.990), #borderpad=0.2,
                      fancybox=False, shadow=False, ncol=ncol_legend, numpoints=2)
      #lgd.get_frame().set_boxstyle('round')

    plot_2 = [ax.plot([], [],
      linewidth=3, color=colormap[densities.index(den)],
      linestyle='-') for den in densities]
    legend_2 = ax.legend(plot_2, ['n=' + scientific2latex('{0:7.1E}'.format(den)) + ' cm$^{{-3}}$' for den in densities],
               bbox_to_anchor=(0.01, 0.85), loc='upper left', fancybox=False, shadow=False)
    #legend_2.legendPatch.set_alpha(0.5)
    ax.add_artist(lgd)

  tmp_list = [speciesPlot[i]['name'] for i in xrange(len(speciesPlot))]
  fig_filename = figname_prefix + '_'.join(tmp_list) + '.eps'
  fig.savefig(os.path.join('.', path_save_fig, fig_filename))
  plt.close()
  print 'gv ' + os.path.join('.', path_save_fig, fig_filename) + ' &'
