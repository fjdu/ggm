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

nSplit = 80
nMarkEvery = 5
figure_size = (8, 6)

markers    = ['+', '*', '.', '1', '2', '3', '4', '<', '>',
              'D', 'H', '^', 'd', 'h', 'o', 'p', 's', 'v', 'x', '|',
              '+', '*', '.', '1', '2', '3', '4', '<', '>',
              'D', 'H', '^', 'd', 'h', 'o', 'p', 's', 'v', 'x', '|']
markersize = [10 ,  6 ,  6 ,  8 ,  8 ,  8 ,  8 ,  4 ,  4 ,
               4 ,  6 ,  4 ,  4 ,  6 ,  4 ,  4 ,  4 ,  4 ,  4 ,  6 ,
              10 ,  6 ,  6 ,  8 ,  8 ,  8 ,  8 ,  4 ,  4 ,
               4 ,  6 ,  4 ,  4 ,  6 ,  4 ,  4 ,  4 ,  4 ,  4 ,  6 ]
linestyle = [':', '--', '-.', '-']

data_files.sort(key = lambda item: item['density'])
data_files.sort(key = lambda item: item['Temperature'])
#Temperature_s = numpy.array([tmp['Temperature'] for tmp in data_files])
#n_H_s = numpy.array([tmp['density'] for tmp in data_files])

n_files = len(data_files)

for speciesPlot in speciesPlot_s:
  nSpeciesPlot = len(speciesPlot)
  for i in xrange(nSpeciesPlot):
    speciesPlot[i].update({'color':colors[i]})
    speciesPlot[i].update({'marker':markers[i]})
    speciesPlot[i].update({'markersize':markersize[i]})
    speciesPlot[i].update({'time':[range(nSplit) for i1 in range(n_files)]})
    speciesPlot[i].update({'data':[range(nSplit) for i1 in range(n_files)]})
  #  speciesPlot[i].update({'data_all':range(nSplit*n_files)})
  
  # ++++ Loop over all the directories ++++
  for i_file in xrange(n_files):
    dt = data_files[i_file]
    #print dt['path']
    ShapeEvolMo = getNumOfColRowInFile(os.path.join('.',
                    dt['path'], filename_ascii))
    ShapeEvolMo = (ShapeEvolMo[1]-1, ShapeEvolMo[0])
  
    fEvolMoBin = open(os.path.join('.',
          dt['path'], filename_bin), 'rb')
    read_data = numpy.fromfile(file=fEvolMoBin,
          dtype=numpy.dtype('f8')).reshape(ShapeEvolMo)
    fEvolMoBin.close()
    # The abundance is relative to the number of H nuclei.
    read_data[:,1:] *= funcRatioGrainToH(GrainRadius)
  
    f_tmp = open(os.path.join('.', dt['path'], filename_ascii))
    nameSpecies = f_tmp.readline().split()
    f_tmp.close()
  
    ratio_step = 1.1
    #t_range = (read_data[50,0], read_data[-5,0])
    t_range = (1E4, read_data[-5,0])
    step_0 = (t_range[1] - t_range[0]) / ((ratio_step**nSplit-1.0) / (ratio_step-1.0))
    t_range_s = [(0,0) for i in range(nSplit)]
    t_range_s[0]=(t_range[0], t_range[0]+step_0)
    for i in range(1,nSplit):
      step_0 = step_0 * ratio_step
      t_range_s[i] = (t_range_s[i-1][1], t_range_s[i-1][1]+step_0)
  
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
        for i_split in xrange(nSplit):
          iStart = getApproxIdxAsc(read_data[:,0], t_range_s[i_split][0])
          iEnd = getApproxIdxAsc(read_data[:,0], t_range_s[i_split][1])
          #print iStart, iEnd
          speciesPlot[i1]['time'][i_file][i_split] = \
            numpy.average(read_data[iStart:iEnd,0],
              weights=read_data[iStart:iEnd,0]-read_data[iStart-1:iEnd-1,0])
          speciesPlot[i1]['data'][i_file][i_split] = \
            numpy.average(read_data[iStart:iEnd,idx],
              weights=read_data[iStart:iEnd,0]-read_data[iStart-1:iEnd-1,0])
  
  #for i in xrange(nSpeciesPlot):
  #    speciesPlot[i]['data_all'] = [x for xx in speciesPlot[i]['data'] for x in xx]
  
  for i in xrange(1, nSpeciesPlot):
    fig = plt.figure(figsize=figure_size)
    ax = fig.add_subplot(1, 1, 1, autoscale_on=True, 
         xscale='log', yscale='log',
         xlabel='Abundance of '+getNameForPlot(speciesPlot[0]['name']), ylabel='Abundance of '+getNameForPlot(speciesPlot[i]['name']))
    for i_file in xrange(n_files):
      ax.plot([],[], linestyle=linestyle[i_file%4],
          color=colors[i_file], marker='None', markevery=1,
          label='{0:2d}: '.format(i_file+1) + \
                'T=' + '{0:4.1F}'.format(data_files[i_file]['Temperature']) + ' K, ' + \
                'n=' + scientific2latex('{0:7.1E}'.format(data_files[i_file]['density'])) + ' cm$^{{-3}}$')
      ax.plot(speciesPlot[0]['data'][i_file],
          speciesPlot[i]['data'][i_file], linestyle=linestyle[i_file%4],
          color=colors[i_file], marker='None', markevery=1)
      ax.plot(speciesPlot[0]['data'][i_file][0],
          speciesPlot[i]['data'][i_file][0], 
          color=colors[i_file], marker='+')
      ax.plot(speciesPlot[0]['data'][i_file][-1],
          speciesPlot[i]['data'][i_file][-1], 
          color=colors[i_file], marker='|')
     #plt.legend(loc='upper left')
      lgd = ax.legend(loc='best', fancybox=True, shadow=True)
      lgd.legendPatch.set_alpha(0.5)
     #plt.text(speciesPlot[0]['data'][i_file][0], speciesPlot[i]['data'][i_file][0], 
     #  '{0:2d}: '.format(i_file+1) + scientific2latex('{0:9.1E}'.format(speciesPlot[i]['time'][i_file][0])) + ' yr')
     #plt.text(speciesPlot[0]['data'][i_file][-1], speciesPlot[i]['data'][i_file][-1], 
     #  '{0:2d}: '.format(i_file+1) + scientific2latex('{0:9.1E}'.format(speciesPlot[i]['time'][i_file][-1])) + ' yr')
    fig.savefig(os.path.join('.', path_save_fig, speciesPlot[i]['name'] + '_vs_' + speciesPlot[0]['name'] + '.pdf'))
    plt.close()

time_end = datetime.datetime.now()
print 'Finished at ', time_end
print 'A time span of ', time_end - time_start, ' has elapsed.'
