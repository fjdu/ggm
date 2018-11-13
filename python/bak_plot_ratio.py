import sys
import datetime
import math
import os
import shutil
import stat
import string
import numpy
import matplotlib
import matplotlib.ticker as tk

time_start = datetime.datetime.now()
print 'Started at ', time_start

matplotlib.use('cairo.ps')
#matplotlib.rc('text', usetex=True)
import matplotlib.pyplot as plt

# ++++ Initialization begin ++++
pathAbsScript = os.path.dirname(os.path.abspath(sys.argv[0]))

execfile(os.path.join(pathAbsScript, 'config.py'))
execfile(os.path.join(pathAbsScript, 'multiRunPreProCommon.py'))
execfile(os.path.join(pathAbsScript, 'MySubsoutines.py'))
execfile(os.path.join(pathAbsScript, 'make_colors_auto.py'))
execfile(os.path.join(pathAbsScript, 'data_files_plot.py'))

# ---- Initialization end ----

xRange = (1e3, 1e7)
yRange = (1e-4, 5e1)
nSplit = 50
nMarkEvery = 5

markers    = ['+', '*', '.', '1', '2', '3', '4', '<', '>',
              'D', 'H', '^', 'd', 'h', 'o', 'p', 's', 'v', 'x', '|',
              '+', '*', '.', '1', '2', '3', '4', '<', '>',
              'D', 'H', '^', 'd', 'h', 'o', 'p', 's', 'v', 'x', '|']
markersize = [10 ,  6 ,  6 ,  8 ,  8 ,  8 ,  8 ,  4 ,  4 ,
               4 ,  6 ,  4 ,  4 ,  6 ,  4 ,  4 ,  4 ,  4 ,  4 ,  6 ,
              10 ,  6 ,  6 ,  8 ,  8 ,  8 ,  8 ,  4 ,  4 ,
               4 ,  6 ,  4 ,  4 ,  6 ,  4 ,  4 ,  4 ,  4 ,  4 ,  6 ]
linestyle = ['-', '--', ':', '-.']
colors = [(1.0,   0,   0),
          (  0, 0.7,   0),
          (0.0, 0.0, 1.0),
          (  0, 0.8, 1.0),
          (0.8,   0, 1.0),
          (1.0, 1.0,   0),
          (0.0, 0.0,   0),
         ]
ncolors = len(colors)
nlinestyle = len(linestyle)

data_files.sort(key = lambda item: item['density'])
data_files.sort(key = lambda item: item['Temperature'])

n_files = len(data_files)

densities = list(set([data_files[i]['density'] for i in range(n_files)]))
n_densities = len(densities)
n_T = n_files/n_densities


for speciesPlot in speciesPlot_s:
  nSpeciesPlot = len(speciesPlot)
  for i in xrange(nSpeciesPlot):
      speciesPlot[i].update({'color':colors[i]})
      speciesPlot[i].update({'marker':markers[i]})
      speciesPlot[i].update({'markersize':markersize[i]})
      speciesPlot[i].update({'data':[[] for j in range(n_T)]})
      speciesPlot[i].update({'time':[[] for j in range(n_T)]})
      speciesPlot[i].update({'T':[[] for j in range(n_T)]})
  #    speciesPlot[i].update({'data_all':range(nSplit*n_files)})
  
  fig = plt.figure(figsize=(8, 6))
  ax = fig.add_subplot(1, 1, 1, autoscalex_on=False, autoscaley_on=False, 
     xscale='log', yscale='log', xlim=xRange, ylim=yRange,
  #   title=r'n$_{{\rm H}}$ = {0:7.1E} cm$^{{-3}}$'.format(density),
     xlabel='Time (yr)', ylabel='[' + getNameForPlot(speciesPlot[i]['name']) + ']/[' + \
     getNameForPlot(speciesPlot[0]['name'] + ']'))
  # ++++ Loop over all the directories ++++
  for i_density in xrange(len(densities)):
    density = densities[i_density]
    i_T = -1
    for i_file in xrange(n_files):
      dt = data_files[i_file]
      if dt['density'] != density:
        continue
      else:
        i_T += 1
      ShapeEvolMo = getNumOfColRowInFile(os.path.join('.',
                      dt['path'], filename_ascii))
      ShapeEvolMo = (ShapeEvolMo[1]-1, ShapeEvolMo[0])
  
      fEvolMoBin = open(os.path.join('.',
            dt['path'], filename_bin), 'rb')
      read_data = numpy.fromfile(file=fEvolMoBin,
            dtype=numpy.dtype('f8')).reshape(ShapeEvolMo)
      fEvolMoBin.close()
  
      f_tmp = open(os.path.join('.', dt['path'], filename_ascii))
      nameSpecies = f_tmp.readline().split()
      f_tmp.close()
  
      for i1 in xrange(nSpeciesPlot):
        flag_found = False
        for i in xrange(ShapeEvolMo[1]):
          if nameSpecies[i] == speciesPlot[i1]['name']:
            idx = i
            flag_found = True
            break
        if not flag_found:
          print 'Cannot find ', speciesPlot[i1]['name']
        else:
          speciesPlot[i1]['data'][i_T] = read_data[:,idx] * funcRatioGrainToH(GrainRadius)
          speciesPlot[i1]['time'][i_T] = read_data[:,0]
          speciesPlot[i1]['T'][i_T] = dt['Temperature']
  
    for i in xrange(1, nSpeciesPlot):
      for i_T in xrange(n_T):
        ratio = speciesPlot[0]['time'][i_T].copy()
        for j in xrange(len(ratio)):
          ratio[j] = speciesPlot[i]['data'][i_T][j] / speciesPlot[0]['data'][i_T][j]
        ax.plot(speciesPlot[i]['time'][i_T], ratio, linestyle=linestyle[i_density%nlinestyle],
          color=colors[i_T%ncolors], marker='None', markevery=1, linewidth=3.0,
          label='T={0:7.1f}  n={1:7.1e}'.format(speciesPlot[i]['T'][i_T], densities[i_density]))
  
  lgd = ax.legend(loc='best', fancybox=True, shadow=True)
  lgd.legendPatch.set_alpha(0.5)
  fig.savefig(os.path.join('.', path_save_fig, speciesPlot[i]['name'] + \
      '_divided_by_' + speciesPlot[0]['name'] + '.pdf'))
  plt.close()

time_end = datetime.datetime.now()
print 'Finished at ', time_end
print 'A time span of ', time_end - time_start, ' has elapsed.'
