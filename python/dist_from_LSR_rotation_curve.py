# All angles are in RADIANS, unless keyword degree is set.

# Given the Galactic longitude and latitude, and the
# LSR velocity, calculate the distance to the Sun.

# Reference:
# Fich, M., Blitz, L., Stark, A., 1989, ApJ, 342, 272
#w0 = v0 / R0
#V_II = 4.2 ; km/s
#R = sqrt(R0^2 + d^2 -2*R0*d*cos(l))
#Vr = (rotation_curve(R) - 1) * V0 * sin(l) *cos(b)

# Author: DFJ

import math
import numpy
import pylab

def dist_from_LSR_rotation_curve(l, b, vLSR, 
    OutOne = True, noNaN = True,
    use_degree = True, noV_II = True):
    

  R0 = 8.5
  V_0 = 220.0
  V_II = 0.0 #4.2 * (~keyword_set(noV_II))
  
  if use_degree:
    l_tmp = l * math.pi / 180.0
    b_tmp = b * math.pi / 180.0
  else:
    l_tmp = l
    b_tmp = b
  
  sinL = math.sin(l_tmp)
  cosL = math.cos(l_tmp)
  cosB = math.cos(b_tmp)
  
  R = 1.00746 * R0 / \
     (vLSR / (V_0 * sinL * cosB) - (V_II/V_0)*(cosL/sinL) \
           + 1.017112)
  
  R_to_GalCen = R
  
  tmp1 = R0 * cosL
  
  if noNaN:
    tmp2 = math.sqrt(max(R*R - (R0*sinL)*(R0*sinL), 0.0))
  else:
    tmp2 = math.sqrt(R*R  - (R0*sinL)*(R0*sinL))
  
  if OutOne:
    return (tmp1+tmp2) if tmp1 < tmp2 else (tmp1-tmp2)
  else:
    return [tmp1-tmp2, tmp1+tmp2]

n = 50
lon_s = numpy.linspace(5.0, 175.0, n)
lat = 0.0
vLSR = 50.0
dV = 0.1
ddv = numpy.empty_like(lon_s)

for i in xrange(len(lon_s)):
  lon = lon_s[i]
  d1 = dist_from_LSR_rotation_curve(lon, lat, vLSR)
  d2 = dist_from_LSR_rotation_curve(lon, lat, vLSR+dV)
  ddv[i] = abs((d1-d2)/dV)
  print '{0:7.2f}, {1:7.2f}, {2:7.5f}'.format(lon, d1, ddv[i])

pylab.plot(lon_s, ddv)
pylab.show()
