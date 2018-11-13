my_colors = \
[
  (0.0, 0.0, 0.0), #black
  (1.0, 0.0, 0.0), #r
  (0.0, 0.8, 0.0), #g
  (0.0, 0.0, 1.0), #b

  (0.3, 0.3, 0.3), #black
  (1.0, 0.0, 1.0),
  (0.0, 0.7, 0.7),
  (0.7, 0.7, 0.0),

  (0.6, 0.6, 0.6), #black
  (0.8, 0.0, 0.0), #r
  (0.0, 0.5, 0.0), #g
  (0.0, 0.0, 0.6), #b

  (1.0, 0.0, 0.5),
  (0.0, 0.3, 0.7),
  (0.7, 0.3, 0.0),

  (1.0, 0.4, 0.5),
  (0.3, 0.3, 0.7),
  (0.7, 0.3, 0.4),

  (0.8, 0.6, 0.2),
  (0.6, 0.8, 0.2),
  (0.6, 0.2, 0.8),

  (1.0, 1.0, 1.0),
]

#print len(my_colors)
#
#from pylab import *
#
#clf()
#for i in xrange(len(my_colors)):
#  plot([],[], color=my_colors[i], label=str(i), marker='s', markersize=15)
#
#legend(handlelength=10, ncol=2)
#ion()
#show()

def my_custom_colormap(i,j):
  return my_colors[i]
