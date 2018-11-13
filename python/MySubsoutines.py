markers    = ['+', '*', '1', '2', '3', '4', '<', '>', '.',
              'D', 'H', '^', 'd', 'h', 'o', 'p', 's', 'v', 'x', '|']
#markersize = [10 ,  6 ,  6 ,  8 ,  8 ,  8 ,  8 ,  4 ,  4 ,
#               4 ,  6 ,  4 ,  4 ,  6 ,  4 ,  4 ,  4 ,  4 ,  4 ,  6 ,
#              10 ,  6 ,  6 ,  8 ,  8 ,  8 ,  8 ,  4 ,  4 ,
#               4 ,  6 ,  4 ,  4 ,  6 ,  4 ,  4 ,  4 ,  4 ,  4 ,  6 ]
linestyle = ['-', '--', ':', '-.']
linestyle_dens = \
            [1.0, 0.5, 0.3, 0.4]
colors = [(1.0,   0,   0),
          (  0, 0.7,   0),
          (0.0, 0.0, 1.0),
          (  0, 0.8, 1.0),
          (0.8,   0, 1.0),
          (0.8, 0.8,   0),
          (0.0, 0.0,   0),
         ]
nmarkers = len(markers)
nlinestyle = len(linestyle)
ncolors = len(colors)


def RGB2XYZ(rgb):
  M = [[0.5767309, 0.1855540, 0.1881852],
       [0.2973769, 0.6273491, 0.0752741],
       [0.0270343, 0.0706872, 0.9911085],
      ]
  xyz = [0.0, 0.0, 0.0]
  for i in xrange(3):
    for j in xrange(3):
      xyz[i] += M[i][j] * rgb[j]
  return xyz

def RGB2UV(rgb):
  xyz = RGB2XYZ(rgb)
  u = 4 * xyz[0] / (-2*xyz[0] + 12*xyz[1] + 3)
  v = 6 * xyz[1] / (-2*xyz[0] + 12*xyz[1] + 3)
  return [u, v]

def getNumOfColRowInFile(filename):
    # Ref: http://stackoverflow.com/questions/845058/how-to-get-line-count-cheaply-in-python
    import mmap, string
    f = open(filename, "r+")
    buf = mmap.mmap(f.fileno(), 0)
    rows = 0
    readline = buf.readline
    while readline():
        rows += 1
    f.seek(0,0)
    cols = len(string.split(f.readline()))
    f.close()
    return cols, rows


def List2Dict(L):
    d = {}
    for i in L:
        d.setdefault(i)
    return d


def getApproxIdxAsc(L, key):
    if key < L[0]:
      print key, L[0]
      raise NameError('The key is too small!')
    for i in xrange(len(L)-1):
      if L[i] <= key and L[i+1] > key:
        return i
    if L[-1] != key:
      print key, L[-1]
      raise NameError('The key is too big!')
    else:
      return len(L)-1

def getApproxIdxAsc_min(L, key):
  import numpy
  return numpy.argmin(numpy.abs(L-key))

def species2latex(strSpecies):
    import re as RegularExpression
    str_tmp = RegularExpression.sub(r'(\d+)', r'$_{\1}$', strSpecies)
    str_tmp = RegularExpression.sub(r'([+-]+)$', r'$^{\1}$', str_tmp)
    return str_tmp


def getNameForPlot(strSpecies):
    str_tmp = species2latex(strSpecies)
    if str_tmp[0] == 'g':
      str_tmp = 'g' + str_tmp[1:]
    if str_tmp[0] == 'm':
      str_tmp = 'm' + str_tmp[1:]
    return str_tmp


def scientific2latex(sstr):
    import re as RegularExpression
    if sstr.startswith('1.0'):
      str_tmp = RegularExpression.sub(r'([0-9.]+)[eEdD]{1}([+-]{0,1})(0*)([0-9.]+)', r'$10^{\2\4}$', sstr)
    else:
      str_tmp = RegularExpression.sub(r'([0-9.]+)[eEdD]{1}([+-]{0,1})(0*)([0-9.]+)', r'$\1{\\times}10^{\2\4}$', sstr)
    str_tmp = RegularExpression.sub(r'([+-]|[{]+)\.', r'\1 0.', str_tmp)
    str_tmp = RegularExpression.sub(r' ', r'', str_tmp)
    str_tmp = RegularExpression.sub(r'[+]{1}(\d)', r'\1', str_tmp)
    return str_tmp


def list_unique(seq, idfun=None): 
# http://www.peterbe.com/plog/uniqifiers-benchmark
   # order preserving
   if idfun is None:
       def idfun(x): return x
   seen = {}
   result = []
   for item in seq:
       marker = idfun(item)
       # in old Python versions:
       # if seen.has_key(marker)
       # but in new ones:
       if marker in seen: continue
       seen[marker] = 1
       result.append(item)
   return result


def get_data(useNumPerGrain=False):
  for dt in data_files:
    ShapeEvolMo = getNumOfColRowInFile(os.path.join('.',
                  dt['path'], filename_ascii))
    ShapeEvolMo = (ShapeEvolMo[1]-1, ShapeEvolMo[0])
    fEvolMoBin = open(os.path.join('.',
        dt['path'], filename_bin), 'rb')
    dt.update({'bin': numpy.fromfile(file=fEvolMoBin,
        dtype=numpy.dtype('f8')).reshape(ShapeEvolMo)})
    if not useNumPerGrain:
      dt['bin'][:,1:] *= funcRatioGrainToH()
    fEvolMoBin.close()
    f_tmp = open(os.path.join('.', dt['path'], filename_ascii))
    dt.update({'nameSpecies': f_tmp.readline().split()})
    f_tmp.close()



def get_data_cmp(useNumPerGrain=False):
  for dt in data_files:
    ShapeEvolMo = getNumOfColRowInFile(os.path.join('.',
                  dt['comparewith'], filename_ascii))
    ShapeEvolMo = (ShapeEvolMo[1]-1, ShapeEvolMo[0])
    fEvolMoBin = open(os.path.join('.',
        dt['comparewith'], filename_bin), 'rb')
    dt.update({'bin_cmp': numpy.fromfile(file=fEvolMoBin,
        dtype=numpy.dtype('f8')).reshape(ShapeEvolMo)})
    if not useNumPerGrain:
      dt['bin_cmp'][:,1:] *= funcRatioGrainToH()
    fEvolMoBin.close()
    f_tmp = open(os.path.join('.', dt['comparewith'], filename_ascii))
    dt.update({'nameSpecies_cmp': f_tmp.readline().split()})
    f_tmp.close()



def divideIntoLayers():
  N_sites_per_grain = SitesDensity * 4.0*math.pi * (GrainRadius_s[0]*1e2)**2
  strMant = 'totalMant'
  for dt in data_files:
    icol_Mant = dt['nameSpecies'].index(strMant)
    nL_s = numpy.int_(numpy.floor(dt['bin'][:, icol_Mant] / N_sites_per_grain))
    idx_div = numpy.where(numpy.diff(nL_s) >= 1)[0]
    dt.update({'div': idx_div})
    dt.update({'numL': 1+nL_s[idx_div]})
    dt.update({'nLayers': len(idx_div)})
    dt.update({'N_S':N_sites_per_grain})
    idx_div_ = numpy.empty(len(idx_div)+1, dtype=numpy.int)
    idx_div_[0] = 0
    idx_div_[1:] = idx_div
    dt.update({'N_in': numpy.diff(dt['bin'][idx_div_, icol_Mant])})

def get_column(speName, dt):
  #for i in xrange(len(dt['nameSpecies'])):
  #  if dt['nameSpecies'][i] == speName:
  #    return i
  #    break
  #return None
  return dt['nameSpecies'].index(speName)


def smooth(x, nwin=10):
  win = [1.0/nwin]*nwin
  return numpy.convolve(x, win, 'same')


def make_polygon(x, y1, y2):
  print x.shape, y1.shape, y2.shape
  nx = len(x)
  n = 2*nx + 1
  p = numpy.zeros((n,2))
  p[0:nx, 0] = x
  p[0:nx, 1] = y1
  p[nx:(n-1), 0] = x[-1::-1]
  p[nx:(n-1), 1] = y2[-1::-1]
  p[n-1,:] = p[0,:]
  return p


def make_polygon_tran(x, y1, y2):
  print x.shape, y1.shape, y2.shape
  nx = len(x)
  n = 2*nx + 1
  p = numpy.zeros((n,2))
  p[0:nx, 0] = y1
  p[0:nx, 1] = x
  p[nx:(n-1), 0] = y2[-1::-1]
  p[nx:(n-1), 1] = x[-1::-1]
  p[n-1,:] = p[0,:]
  return p

def copy_by_cmp(src, dest, chmodexe=True):
  import os, filecmp
  if not os.access(dest, os.F_OK):
    shutil.copyfile(src, dest)
    if chmodexe:
      os.chmod(dest, stat.S_IRWXU | stat.S_IWOTH)
  if not filecmp.cmp(src, dest):
    shutil.copyfile(src, dest)
    if chmodexe:
      os.chmod(dest, stat.S_IRWXU | stat.S_IWOTH)
