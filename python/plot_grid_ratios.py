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

#data_files = [\
#{'Temperature': 10.0, 'density': 1.0E+03,
#  'path':        "RE_0.01/"    + "10.0_1.0E+03_1.0E-07_2.0_1.0E-02/"},
#{'Temperature': 10.0, 'density': 1.0E+04,
#  'path':        "RE_0.01/"    + "10.0_1.0E+04_1.0E-07_15.0_1.0E-02/"},
#{'Temperature': 10.0, 'density': 1.0E+05,
#  'path':        "RE_0.01/"    + "10.0_1.0E+05_1.0E-07_15.0_1.0E-02/"},
#{'Temperature': 10.0, 'density': 1.0E+06,
#  'path':        "RE_0.01/"    + "10.0_1.0E+06_1.0E-07_15.0_1.0E-02/"},
#]
#path_save_fig = ''
#filename_bin = "evolution_moment__bin.dat"
#filename_ascii = "evolution_moment__ascii.dat"
#get_data()
#figname_prefix = 'RE_ratios_'

#data_files = [\
#{'Temperature': 10.0, 'density': 1.0E+03,
#  'path':        "RE_0.01/"    + "10.0_1.0E+03_1.0E-07_2.0_1.0E-02/"},
#{'Temperature': 10.0, 'density': 1.0E+04,
#  'path':        "RE_0.01/"    + "10.0_1.0E+04_1.0E-07_15.0_1.0E-02/"},
#{'Temperature': 10.0, 'density': 1.0E+05,
#  'path':        "RE_0.01/"    + "10.0_1.0E+05_1.0E-07_15.0_1.0E-02/"},
#{'Temperature': 10.0, 'density': 1.0E+06,
#  'path':        "RE_0.01/"    + "10.0_1.0E+06_1.0E-07_15.0_1.0E-02/"},
#{'Temperature': 15.0, 'density': 1.0E+03,
#  'path':        "RE_0.01/"    + "15.0_1.0E+03_1.0E-07_2.0_1.0E-02/"},
#{'Temperature': 15.0, 'density': 1.0E+04,
#  'path':        "RE_0.01/"    + "15.0_1.0E+04_1.0E-07_15.0_1.0E-02/"},
#{'Temperature': 15.0, 'density': 1.0E+05,
#  'path':        "RE_0.01/"    + "15.0_1.0E+05_1.0E-07_15.0_1.0E-02/"},
#{'Temperature': 15.0, 'density': 1.0E+06,
#  'path':        "RE_0.01/"    + "15.0_1.0E+06_1.0E-07_15.0_1.0E-02/"},
#{'Temperature': 20.0, 'density': 1.0E+03,
#  'path':        "RE_0.01/"    + "20.0_1.0E+03_1.0E-07_2.0_1.0E-02/"},
#{'Temperature': 20.0, 'density': 1.0E+04,
#  'path':        "RE_0.01/"    + "20.0_1.0E+04_1.0E-07_15.0_1.0E-02/"},
#{'Temperature': 20.0, 'density': 1.0E+05,
#  'path':        "RE_0.01/"    + "20.0_1.0E+05_1.0E-07_15.0_1.0E-02/"},
#{'Temperature': 20.0, 'density': 1.0E+06,
#  'path':        "RE_0.01/"    + "20.0_1.0E+06_1.0E-07_15.0_1.0E-02/"},
#]
#path_save_fig = ''
#filename_bin = "evolution_moment__bin.dat"
#filename_ascii = "evolution_moment__ascii.dat"
#get_data()
#figname_prefix = 'RE_ratios_'

data_files = [\
{'Temperature': 10.0, 'density': 1.0E+03,
  'path':        "HME_0.01/"    + "10.0_1.0E+03_1.0E-07_15.0_1.0E-02/"},
{'Temperature': 10.0, 'density': 1.0E+04,
  'path':        "HME_0.01/"    + "10.0_1.0E+04_1.0E-07_15.0_1.0E-02/"},
{'Temperature': 10.0, 'density': 1.0E+05,
  'path':        "HME_0.01/"    + "10.0_1.0E+05_1.0E-07_15.0_1.0E-02/"},
{'Temperature': 10.0, 'density': 1.0E+06,
  'path':        "HME_0.01/"    + "10.0_1.0E+06_1.0E-07_15.0_1.0E-02/"},
{'Temperature': 15.0, 'density': 1.0E+03,
  'path':        "HME_0.01/"    + "15.0_1.0E+03_1.0E-07_2.0_1.0E-02/"},
{'Temperature': 15.0, 'density': 1.0E+04,
  'path':        "HME_0.01/"    + "15.0_1.0E+04_1.0E-07_15.0_1.0E-02/"},
{'Temperature': 15.0, 'density': 1.0E+05,
  'path':        "HME_0.01/"    + "15.0_1.0E+05_1.0E-07_15.0_1.0E-02/"},
{'Temperature': 15.0, 'density': 1.0E+06,
  'path':        "HME_0.01/"    + "15.0_1.0E+06_1.0E-07_15.0_1.0E-02/"},
{'Temperature': 20.0, 'density': 1.0E+03,
  'path':        "HME_0.01/"    + "20.0_1.0E+03_1.0E-07_2.0_1.0E-02/"},
{'Temperature': 20.0, 'density': 1.0E+04,
  'path':        "HME_0.01/"    + "20.0_1.0E+04_1.0E-07_15.0_1.0E-02/"},
{'Temperature': 20.0, 'density': 1.0E+05,
  'path':        "HME_0.01/"    + "20.0_1.0E+05_1.0E-07_15.0_1.0E-02/"},
{'Temperature': 20.0, 'density': 1.0E+06,
  'path':        "HME_0.01/"    + "20.0_1.0E+06_1.0E-07_15.0_1.0E-02/"},
{'Temperature': 30.0, 'density': 1.0E+03,
  'path':        "HME_0.01/"    + "30.0_1.0E+03_1.0E-07_2.0_1.0E-02/"},
{'Temperature': 30.0, 'density': 1.0E+04,
  'path':        "HME_0.01/"    + "30.0_1.0E+04_1.0E-07_15.0_1.0E-02/"},
{'Temperature': 30.0, 'density': 1.0E+05,
  'path':        "HME_0.01/"    + "30.0_1.0E+05_1.0E-07_15.0_1.0E-02/"},
{'Temperature': 30.0, 'density': 1.0E+06,
  'path':        "HME_0.01/"    + "30.0_1.0E+06_1.0E-07_15.0_1.0E-02/"},
]
path_save_fig = ''
filename_bin = "evolution_moment__bin.dat"
filename_ascii = "evolution_moment__ascii.dat"
get_data()
figname_prefix = 'HME_ratios_'

#speciesPlot_s = [
#  [
#     {'name': 'mH2O'},
#     {'name': 'mHDO'},
#     {'name': 'mH2CO'},
#     {'name': 'mHDCO'},
#     {'name': 'mCH3OH'},
#     {'name': 'mCH2DOH'},
#     {'name': 'mCH3OH'},
#     {'name': 'mCH3OD'},
#  ],
#              ]
#xRange = (1e2, 1e6)
#yRange = (1e-3, 4.0)
#ncol_legend = 2
#nMarkPoints = 10
#linewidth = 1.0
#markersize = 10.0
#markeredgewidth = 1.5
#markeredgecolor = 'black'
#plotWidth = 0.87
#plotheight = 0.7
#figsize = (8,10)
#x_begin = 0.1
#y_begin = 0.07
#execfile(os.path.join(pathAbsScript, 'make_colors_auto_my.py'))
#colormap = my_colors

#speciesPlot_s = [
#  [
#     {'name': 'H2', 'drawVLineAt': 0.9},
#     {'name': 'mH2O'},
#     {'name': 'H2', 'drawVLineAt': 0.9},
#     {'name': 'mHDO'},
#     {'name': 'H2', 'drawVLineAt': 0.9},
#     {'name': 'mD2O'},
#     {'name': 'gH', 'multiply': 1e-7},
#     {'name': 'gD'},
#  ],
#              ]
#xRange = (1e2, 9e6)
#yRange = (1e-10, 5e-4)
#ncol_legend = 2
#nMarkPoints = 6
#linewidth = 1.0
#markersize = 10.0
#markeredgewidth = 1.5
#markeredgecolor = 'black'
#plotWidth = 0.87
#plotheight = 0.77
#figsize = (8,8)
#x_begin = 0.1
#y_begin = 0.07
#execfile(os.path.join(pathAbsScript, 'make_colors_auto_my.py'))
#colormap = my_colors

speciesPlot_s = [
  [
     {'name': 'H2', 'drawHLineAt': 6e-9},
     {'name': 'NH'},
     {'name': 'H2', 'drawHLineAt': 3e-9},
     {'name': 'NH2'},
     {'name': 'H2', 'drawHLineAt': 3e-9},
     {'name': 'NH3'},
     #{'name': 'H2', 'drawHLineAt': 1e-10},
     #{'name': 'NH+'},
  ],
              ]
xRange = (1e2, 9e6)
yRange = (1e-13, 5e-7)
ncol_legend = 2
nMarkPoints = 6
linewidth = 1.0
markersize = 10.0
markeredgewidth = 1.5
markeredgecolor = 'black'
plotWidth = 0.87
plotheight = 0.77
figsize = (10,10)
x_begin = 0.1
y_begin = 0.07
execfile(os.path.join(pathAbsScript, 'make_colors_auto_my.py'))
colormap = my_colors

#speciesPlot_s = [
#  [
#     {'name': 'H2'},
#     {'name': 'H'},
#     {'name': 'H2'},
#     {'name': 'D'},
#     {'name': 'H2'},
#     {'name': 'CO'},
#  ],
#              ]
#xRange = (1e-4, 1e7)
#yRange = (1e-8, 1.0)
#ncol_legend = 2
#nMarkPoints = 10
#linewidth = 1.0
#markersize = 10.0
#markeredgewidth = 1.5
#markeredgecolor = 'black'
#plotWidth = 0.87
#plotheight = 0.7
#figsize = (8,10)
#x_begin = 0.1
#y_begin = 0.07
#execfile(os.path.join(pathAbsScript, 'make_colors_auto_my.py'))
#colormap = my_colors

#speciesPlot_s = [
#  [
#     {'name': 'H'},
#     {'name': 'D'},
#     {'name': 'gH'},
#     {'name': 'gD'},
#  ],
#              ]
#xRange = (1e-4, 1e7)
#yRange = (1e-3, 1e4)
#ncol_legend = 2
#nMarkPoints = 10
#linewidth = 1.0
#markersize = 10.0
#markeredgewidth = 1.5
#markeredgecolor = 'black'
#plotWidth = 0.87
#plotheight = 0.7
#figsize = (8,10)
#x_begin = 0.1
#y_begin = 0.07
#execfile(os.path.join(pathAbsScript, 'make_colors_auto_my.py'))
#colormap = my_colors

#speciesPlot_s = [
#  [
#     {'name': 'H2'},
#     {'name': 'CO'},
#     {'name': 'H2'},
#     {'name': 'O2'},
#  ],
#              ]
#xRange = (1e2, 1e7)
#yRange = (1e-7, 5e-4)
#ncol_legend = 2
#nMarkPoints = 10
#linewidth = 1.0
#markersize = 10.0
#markeredgewidth = 1.5
#markeredgecolor = 'black'
#plotWidth = 0.87
#plotheight = 0.7
#figsize = (8,10)
#x_begin = 0.1
#y_begin = 0.07
#execfile(os.path.join(pathAbsScript, 'make_colors_auto_my.py'))
#colormap = my_colors

#speciesPlot_s = [
#  [
#     {'name': 'H2'},
#     {'name': 'H2O2'},
#     {'name': 'H2'},
#     {'name': 'O2H'},
#     {'name': 'H2'},
#     {'name': 'H2O'},
#     {'name': 'H2'},
#     {'name': 'O2'},
#  ],
#              ]
#xRange = (1e3, 1e7)
##yRange = (1e-6, 60)
#yRange = (1e-11, 4e-4)
#ncol_legend = 2
#color_range = 1.8

#speciesPlot_s = [
#  [
#     {'name': 'mH2CO'},
#     {'name': 'mHDCO'},
#     {'name': 'mH2CO'},
#     {'name': 'mD2CO'},
#     {'name': 'H2', 'multiply':1e5},
#     {'name': 'mH2CO'},
#  ],
#              ]
#xRange = (1e3, 1e7)
##yRange = (1e-6, 60)
#yRange = (1e-2, 3)
#ncol_legend = 2
#color_range = 1.8

#speciesPlot_s = [
#  [
#     {'name': 'H2CO'},
#     {'name': 'HDCO'},
#     {'name': 'H2CO'},
#     {'name': 'D2CO'},
#     {'name': 'H2', 'multiply':1e9},
#     {'name': 'H2CO'},
#  ],
#              ]
#xRange = (1e3, 1e7)
#yRange = (1e-6, 60)
#ncol_legend = 1
#color_range = 1.8

#speciesPlot_s = [
#  [
#     {'name': 'H'},
#     {'name': 'D'},
#     {'name': 'gH'},
#     {'name': 'gD'},
#     {'name': 'H2', 'multiply':1e0},
#     {'name': 'H'},
#     {'name': 'H2', 'multiply':1e15},
#     {'name': 'gH'},
#     {'name': 'H2'},
#     {'name': 'D'},
#     {'name': 'H2'},
#     {'name': 'CO'},
#     {'name': 'H2'},
#     {'name': 'mCO2'},
#     {'name': 'H2'},
#     {'name': 'O'},
#     {'name': 'H2'},
#     {'name': 'O2'},
#     {'name': 'H2'},
#     {'name': 'N2'},
#     {'name': 'H2'},
#     {'name': 'mN2'},
#  ],
#              ]
#xRange = (1e3, 1e7)
#yRange = (1e-8, 1)
#ncol_legend = 2
#color_range = 0.9

#speciesPlot_s = [
#  [
#     {'name': 'mCH3OH'},
#     {'name': 'mCH2DOH'},
#     {'name': 'mCH3OH'},
#     {'name': 'mCH3OD'},
#     {'name': 'mCH3OH'},
#     {'name': 'mCHD2OH'},
#     {'name': 'mCH3OH'},
#     {'name': 'mCD3OH'},
#     {'name': 'H2', 'multiply':1e7},
#     {'name': 'mCH3OH'},
#  ],
#              ]
#xRange = (1e3, 1e7)
#yRange = (1e-1, 20)
#ncol_legend = 2
#color_range = 1.8

#speciesPlot_s = [
#  [
#     {'name': 'OH'},
#     {'name': 'OD'},
#     {'name': 'HDO'},
#     {'name': 'OD'},
#     {'name': 'H2O'},
#     {'name': 'OH'},
#     {'name': 'H2', 'multiply':1e7},
#     {'name': 'HDO'},
#     {'name': 'H2', 'multiply':1e7},
#     {'name': 'OD'},
#  ],
#              ]
#xRange = (3e2, 9.9e6)
#yRange = (3e-4, 1e2)
#ncol_legend = 2
#color_range = 0.99

#speciesPlot_s = [
#  [
#     {'name': 'H2O'},
#     {'name': 'HDO'},
#     {'name': 'OH'},
#     {'name': 'OD'},
#     {'name': 'HDO'},
#     {'name': 'OD'},
#     {'name': 'H2O'},
#     {'name': 'OH'},
#     {'name': 'H2', 'multiply':1e7},
#     {'name': 'HDO'},
#     {'name': 'H2', 'multiply':1e7},
#     {'name': 'OH'},
#     {'name': 'H2', 'multiply':1e7},
#     {'name': 'OD'},
#     {'name': 'H2', 'multiply':1e4},
#     {'name': 'D'},
#  ],
#              ]
#xRange = (3e2, 9.9e6)
#yRange = (3e-4, 1e2)
#ncol_legend = 2
#color_range = 0.99

#speciesPlot_s = [
#  [
#     {'name': 'mH2O'},
#     {'name': 'mHDO'},
#     {'name': 'mH2CO'},
#     {'name': 'mHDCO'},
#     {'name': 'mCH3OH'},
#     {'name': 'mCH2DOH'},
#     {'name': 'gH'},
#     {'name': 'gD'},
#     {'name': 'gOH'},
#     {'name': 'gOD'},
#     {'name': 'gHCO'},
#     {'name': 'gDCO'},
#  ],
#              ]
#xRange = (3e2, 9.9e6)
#yRange = (1e-4, 2)
#ncol_legend = 2
#color_range = 0.99

#speciesPlot_s = [
#  [
#     {'name': 'gH2O'},
#     {'name': 'gHDO'},
#     {'name': 'gH2CO'},
#     {'name': 'gHDCO'},
#     {'name': 'gCH3OH'},
#     {'name': 'gCH2DOH'},
#     {'name': 'gH'},
#     {'name': 'gD'},
#     {'name': 'gOH'},
#     {'name': 'gOD'},
#     {'name': 'gHCO'},
#     {'name': 'gDCO'},
#  ],
#              ]
#xRange = (3e2, 9.9e6)
#yRange = (1e-4, 2)
#ncol_legend = 2
#color_range = 0.99

#speciesPlot_s = [
#  [
#     {'name': 'H2', 'multiply':1e-0},
#     {'name': 'mCH3OH'},
#     {'name': 'H2'},
#     {'name': 'mCH2DOH'},
#     {'name': 'H2', 'multiply':1e-0},
#     {'name': 'mH2CO'},
#     {'name': 'H2'},
#     {'name': 'mHDCO'},
#     {'name': 'H', 'multiply':1e-6},
#     {'name': 'D'},
#     {'name': 'gH', 'multiply':1e-6},
#     {'name': 'gD'},
#  ],
#              ]
#xRange = (3e2, 9.9e6)
#yRange = (1e-9, 3e-4)
#ncol_legend = 2
#color_range = 0.99

#speciesPlot_s = [
#  [
#     {'name': 'H2', 'multiply':1e-0},
#     {'name': 'mH2O'},
#     {'name': 'H2', 'multiply':1e-0},
#     {'name': 'mHDO'},
#     {'name': 'H2', 'multiply':1e-0},
#     {'name': 'mCH3OH'},
#     {'name': 'H2'},
#     {'name': 'mCH2DOH'},
#     {'name': 'H2', 'multiply':1e-0},
#     {'name': 'mH2CO'},
#     {'name': 'H2'},
#     {'name': 'mHDCO'},
#     {'name': 'H2', 'multiply':1e0},
#     {'name': 'H'},
#     {'name': 'H2', 'multiply':1e0},
#     {'name': 'D'},
#     {'name': 'H', 'multiply':1e-3},
#     {'name': 'D'},
#  ],
#              ]
#xRange = (3e2, 9.9e6)
#yRange = (5e-9, 2e-2)
#ncol_legend = 3
#color_range = 0.99

#speciesPlot_s = [
#  [
#     {'name': 'H2', 'multiply':1e-0},
#     {'name': 'mHDO'},
#     {'name': 'H2'},
#     {'name': 'mCH2DOH'},
#     {'name': 'H2', 'multiply':1e-0},
#     {'name': 'mHDCO'},
#     {'name': 'H2', 'multiply':1e0},
#     {'name': 'H'},
#     {'name': 'H2', 'multiply':1e0},
#     {'name': 'D'},
#     {'name': 'H', 'multiply':1e-3},
#     {'name': 'D'},
#  ],
#              ]
#xRange = (3e2, 9.9e6)
#yRange = (5e-9, 2e-2)
#ncol_legend = 2
#color_range = 0.99

#speciesPlot_s = [
#  [
#     {'name': 'H2', 'multiply':1e-0},
#     {'name': 'CO'},
#     {'name': 'H2', 'multiply':1e-0},
#     {'name': 'E-'},
#     #{'name': 'H2'},
#     #{'name': 'mCH2DOH'},
#     #{'name': 'H2', 'multiply':1e-0},
#     #{'name': 'mHDCO'},
#     {'name': 'H2', 'multiply':1e6},
#     {'name': 'H2D+'},
#     {'name': 'H2', 'multiply':1e6},
#     {'name': 'HD2+'},
#     {'name': 'H2', 'multiply':1e6},
#     {'name': 'D3+'},
#     {'name': 'H2', 'multiply':1e0},
#     {'name': 'D'},
#     {'name': 'H', 'multiply':1e-2},
#     {'name': 'D'},
#     {'name': 'gH', 'multiply':1e-4},
#     {'name': 'gD'},
#  ],
#              ]
#xRange = (3e2, 9.9e6)
#yRange = (1e-9, 2e-2)
#ncol_legend = 2
#color_range = 0.99

#speciesPlot_s = [
#  [
#     {'name': 'mH2O'},
#     {'name': 'mD2O'},
#     {'name': 'mH2CO'},
#     {'name': 'mD2CO'},
#     {'name': 'mCH3OH'},
#     {'name': 'mCHD2OH'},
#  ],
#              ]
#xRange = (3e2, 1e7)
#yRange = (7e-7, 6.0)
#ncol_legend = 2
##color_range = 1.9

#speciesPlot_s = [
#  [
#     {'name': 'H2'},
#     {'name': 'mH2O'},
#     {'name': 'H2'},
#     {'name': 'mHDO'},
#     {'name': 'H2'},
#     {'name': 'mD2O'},
#     {'name': 'H2', 'multiply':1e9},
#     {'name': 'gOD'},
#     {'name': 'H2'},
#     {'name': 'D'},
#  ],
#              ]
#xRange = (3e2, 9.9e6)
#yRange = (1e-10, 5e-2)
#ncol_legend = 2
#color_range = 1.9

#speciesPlot_s = [
#  [
#     {'name': 'gH'},
#     {'name': 'gD'},
#     {'name': 'H'},
#     {'name': 'D'},
#     #{'name': 'H3+'},
#     #{'name': 'H2D+'},
#     #{'name': 'H3+'},
#     #{'name': 'HD2+'},
#     #{'name': 'H3+'},
#     #{'name': 'D3+'},
#     {'name': 'H2', 'multiply':1e4},
#     {'name': 'CO'},
#     #{'name': 'H2', 'multiply':1e0},
#     #{'name': 'mCO'},
#     #{'name': 'H2', 'multiply':1e0},
#     #{'name': 'mCO2'},
#     #{'name': 'H2', 'multiply':1e0},
#     #{'name': 'O2'},
#     #{'name': 'H2', 'multiply':1e0},
#     #{'name': 'H'},
#     #{'name': 'H2', 'multiply':1e0},
#     #{'name': 'D'},
#     #{'name': 'H2', 'multiply':1e5},
#     #{'name': 'H3+'},
#     #{'name': 'H2', 'multiply':1e5},
#     #{'name': 'H2D+'},
#     #{'name': 'H2', 'multiply':1e5},
#     #{'name': 'E-'},
#     #{'name': 'H2', 'multiply':1e5},
#     #{'name': 'mN2HD'},
#     #{'name': 'H2', 'multiply':1e5},
#     #{'name': 'gN2HD'},
#     #{'name': 'H2', 'multiply':1e8},
#     #{'name': 'gN2D'},
#  ],
#              ]
#xRange = (3e2, 1e7)
#yRange = (1e-4, 4.0e1)
#ncol_legend = 3
#color_range = 0.7

#speciesPlot_s = [
#  [
#     {'name': 'mH2O'},
#     {'name': 'mHDO'},
#     {'name': 'gH2O'},
#     {'name': 'gHDO'},
#     {'name': 'H'},
#     {'name': 'D'},
#     {'name': 'gH'},
#     {'name': 'gD'},
#     {'name': 'gOH'},
#     {'name': 'gOD'},
#  ],
#              ]
#xRange = (3e2, 9.9e6)
#yRange = (5e-4, 30.0)
#ncol_legend = 2
#color_range = 1.0

#speciesPlot_s = [
#  [
#     {'name': 'H'},
#     {'name': 'D'},
#     {'name': 'gH'},
#     {'name': 'gD'},
#     {'name': 'gOH'},
#     {'name': 'gOD'},
#     {'name': 'gH2O'},
#     {'name': 'gHDO'},
#     {'name': 'mH2O'},
#     {'name': 'mHDO'},
#     {'name': 'H2', 'multiply':1e4, 'drawVLineAt':0.9},
#     {'name': 'mH2O'},
#     {'name': 'H2', 'multiply':1e4},
#     {'name': 'mHDO'},
#     {'name': 'H2', 'multiply':1e6},
#     {'name': 'gH2O'},
#     {'name': 'H2', 'multiply':1e6},
#     {'name': 'gHDO'},
#     #{'name': 'H2', 'multiply':1e4, 'drawVLineAt':0.1},
#     #{'name': 'CO'},
#     #{'name': 'H2', 'multiply':1e4},
#     #{'name': 'O2'},
#     #{'name': 'H2', 'multiply':1e4},
#     #{'name': 'N2'},
#     #{'name': 'H2', 'multiply':1e6},
#     #{'name': 'gH2O'},
#     #{'name': 'H2', 'multiply':1e6},
#     #{'name': 'gHDO'},
#  ],
#              ]
#xRange = (3e2, 9.9e6)
#yRange = (5e-4, 30.0)
#ncol_legend = 2
#color_range = 1.0

#speciesPlot_s = [
#  [
#     {'name': 'H2', 'multiply':1e0},
#     {'name': 'gH2O'},
#     {'name': 'H2', 'multiply':1e0},
#     {'name': 'gHDO'},
#     {'name': 'H2', 'multiply':1e0},
#     {'name': 'mH2O'},
#     {'name': 'H2', 'multiply':1e0},
#     {'name': 'mHDO'},
#  ],
#              ]
#xRange = (3e2, 9.9e6)
#yRange = (5e-9, 30.0e-4)
#ncol_legend = 2
#color_range = 0.8

#speciesPlot_s = [
#  [
#     {'name': 'gH2O'},
#     {'name': 'gHDO'},
#     {'name': 'mH2O'},
#     {'name': 'mHDO'},
#  ],
#              ]
#xRange = (3e2, 9.9e6)
#yRange = (6e-4, 3)
#ncol_legend = 2
#color_range = 0.9

#speciesPlot_s = [
#  [
#     {'name': 'H2'},
#     {'name': 'mH2O'},
#     {'name': 'H2'},
#     {'name': 'mHDO'},
#     {'name': 'mH2O'},
#     {'name': 'mHDO'},
#     {'name': 'gH2O'},
#     {'name': 'gHDO'},
#     {'name': 'gH'},
#     {'name': 'gD'},
#  ],
#              ]
#xRange = (3e2, 9.9e6)
#yRange = (6e-8, 2)
#ncol_legend = 2
#color_range = 0.9

#speciesPlot_s = [
#  [
#     {'name': 'H2', 'multiply':1e5},
#     {'name': 'mCH3OH'},
#     {'name': 'H2', 'multiply':1e5},
#     {'name': 'mCH2DOH'},
#     {'name': 'gCH3OH'},
#     {'name': 'gCH2DOH'},
#     {'name': 'gH2CO'},
#     {'name': 'gHDCO'},
#     {'name': 'gH'},
#     {'name': 'gD'},
#  ],
#              ]
#xRange = (3e2, 9.9e6)
#yRange = (1e-4, 30)
#ncol_legend = 2
#color_range = 0.9

#speciesPlot_s = [
#  [
#     {'name': 'gCH3O'},
#     {'name': 'gCH2DO'},
#     {'name': 'gH2CO'},
#     {'name': 'gHDCO'},
#     {'name': 'gHCO'},
#     {'name': 'gDCO'},
#     {'name': 'gOH'},
#     {'name': 'gOD'},
#     {'name': 'gH'},
#     {'name': 'gD'},
#     {'name': 'H'},
#     {'name': 'D'},
#  ],
#              ]
#xRange = (3e2, 9.9e6)
#yRange = (1e-4, 30)
#ncol_legend = 2
#color_range = 0.9

#speciesPlot_s = [
#  [
#     {'name': 'H2'},
#     {'name': 'H'},
#     {'name': 'H2'},
#     {'name': 'D'},
#  ],
#              ]
#xRange = (1.1, 1e7)
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

#speciesPlot_s = [
#  [
#     {'name': 'H2'},
#     {'name': 'mH2O'},
#     {'name': 'H2'},
#     {'name': 'mCO'},
#     {'name': 'H2'},
#     {'name': 'mCO2'},
#     {'name': 'H2'},
#     {'name': 'mCH3OH'},
#     {'name': 'H2'},
#     {'name': 'mH2CO'},
#  ],
#              ]
#xRange = (3E2, 1E7)
#yRange = (1e-7, 1E-3)
#ncol_legend = 3
#color_range = 0.99

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
#     {'name': 'mH2O'},
#     {'name': 'mCH4'},
#     {'name': 'mH2O'},
#     {'name': 'mNH3'},
#     {'name': 'mH2O'},
#     {'name': 'mO2'},
#     {'name': 'mH2O'},
#     {'name': 'mO3'},
#  ],
#              ]
#xRange = (3E2, 1E7)
#yRange = (1e-4, 5E1)
#ncol_legend = 2
#color_range = 0.99

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
#     {'name': 'mH2O'},
#     {'name': 'mCH4'},
#     {'name': 'mH2O'},
#     {'name': 'mNH3'},
#  ],
#              ]
#xRange = (3E2, 1E7)
#yRange = (1e-4, 3E1)
#ncol_legend = 2
#color_range = 0.99

#speciesPlot_s = [
#  [
#     {'name': 'H2'},
#     {'name': 'mH2O'},
#     {'name': 'H2'},
#     {'name': 'mCO'},
#     {'name': 'H2'},
#     {'name': 'mCO2'},
#     {'name': 'H2'},
#     {'name': 'mCH3OH'},
#     {'name': 'H2'},
#     {'name': 'mH2CO'},
#     {'name': 'H2'},
#     {'name': 'mCH4'},
#     {'name': 'H2'},
#     {'name': 'mNH3'},
#  ],
#              ]
#xRange = (3E2, 1E7)
#yRange = (1e-7, 1E-3)
#ncol_legend = 3
#color_range = 0.99

#speciesPlot_s = [
#  [
#     {'name': 'totalSurf'},
#     {'name': 'totalMant'},
#  ],
#              ]
#xRange = (3E2, 1E7)
#yRange = (0, 110)
#ncol_legend = 2
#color_range = 0.9

# ---- Initialization end ----

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
      xscale='log', yscale='log', xlim=xRange, ylim=yRange,
      xlabel=xlabel, ylabel=ylabel)
    if iy < (ny-1):
      ax.set_xticklabels([])
    if ix > 0:
      ax.set_yticklabels([])
    for i in xrange(0,nSpeciesPlot-1,2):
      if speciesPlot[i]['idx'][j] == None:
        continue
      thiscolor = colormap[i/2]
      thismarker = markers[(i/2)%nmarkers]
      ax.plot(dt['bin'][:,0], speciesPlot[i]['ratio'][j],
        linestyle='-',
        color=thiscolor, linewidth=linewidth,
        marker=thismarker, markersize=markersize,
        markeredgewidth=markeredgewidth, markeredgecolor=markeredgecolor,
        markerfacecolor='None', markevery=nMarkEvery)
      if 'drawVLineAt' in speciesPlot[i]:
        maxVal = max(speciesPlot[i]['ratio'][j])
        idx_drawVL = getApproxIdxAsc_min(speciesPlot[i]['ratio'][j], maxVal*speciesPlot[i]['drawVLineAt'])
        ax.axvline(x=dt['bin'][idx_drawVL,0],# ymin=yRange[0], ymax=yRange[1],
          linestyle='--',# alpha=1.0, visible=True,
          color=thiscolor, linewidth=linewidth/2.0)
      if 'drawHLineAt' in speciesPlot[i]:
        ax.axhline(y=speciesPlot[i]['drawHLineAt'],
          linestyle= '--' if i/2%2 else ':',# alpha=1.0, visible=True,
          color=thiscolor, linewidth=linewidth/2.0)
        #ax.plot([dt['bin'][idx_drawVL,0], dt['bin'][idx_drawVL,0]], [yRange[0], yRange[1]], linestyle='--', color=thiscolor, linewidth=linewidth/2.0)
        #print dt['bin'][idx_drawVL,0], yRange[0], yRange[1]
      if j == 0:
        label_str = '[' + getNameForPlot(speciesPlot[i+1]['name']) + \
            '/' + getNameForPlot(speciesPlot[i]['name']) + ']'
        if 'multiply' in speciesPlot[i]:
          if speciesPlot[i]['multiply'] != 1e0:
            label_str = scientific2latex('{0:7.1E}'.format(speciesPlot[i]['multiply'])) + '$\\times$ ' + label_str
        ax.plot([],[],
          linestyle='-',
          #linestyle=linestyle[(i/2)%nlinestyle],
          color=thiscolor, linewidth=linewidth,
          marker=thismarker, markersize=markersize,
          markeredgewidth=markeredgewidth, markeredgecolor=markeredgecolor,
          markerfacecolor='None', markevery=nMarkEvery,
          label=label_str)
      
    props = dict(boxstyle='round', facecolor='None')
    textstr = '(${0:2.0f}$ K, '.format(T) + \
              scientific2latex('{0:7.1E}'.format(den)) + ' cm$^{{-3}}$)'
    ax.text(0.05, 0.95, textstr, transform=ax.transAxes,
            verticalalignment='top', bbox=props)
    if j == 0:
      lgd = ax.legend(loc='lower left', bbox_to_anchor=(0.0, 1.005), borderpad=0.01, fancybox=False, shadow=False, ncol=ncol_legend, numpoints=2)
      lgd.get_frame().set_boxstyle('round')

  tmp_list = [speciesPlot[i]['name'] for i in xrange(len(speciesPlot))]
  fig_filename = figname_prefix + '_'.join(tmp_list) + '.eps'
  fig.savefig(os.path.join('.', path_save_fig, fig_filename))
  plt.close()
  print 'gv ' + os.path.join('.', path_save_fig, fig_filename) + ' &'
