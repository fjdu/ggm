import pylab as pl
import scipy.constants as C

def get_k11(nu, Ediff, N_s):
  return nu * pl.exp(-2*a/C.hbar * pl.sqrt(2*C.m_p*C.k*Ediff)) / N_s

def get_ke(nu, Eb, T):
  return nu * pl.exp(-Eb/T)

def get_gH_HME(H2, k_cos, k11, ke):
  b = k11 / (k11 + ke)
  c = ke / (H2 * k_cos)
  return 2.0 / (b + pl.sqrt(b*b + 4.0*b*c))

def get_ka(nH, r_grain, T, m):
  v = pl.sqrt(8*C.k*T / (C.pi * m)) * 1e2
  return C.pi * r_grain*r_grain * v * nH

def get_gH_RE(H2, k_cos, k11):
  return pl.sqrt(H2 * k_cos / (2.0*k11))

def get_H_RE(H2, k_cos, ke, gH, ka):
  return (H2 * k_cos + ke * gH) / ka

T = 10.0
H2 = 1.7887E+011
k_cos = 1.3e-17

N_sites = 1e6
nu = 2.4e12
a = 1.0e-10

Eb = pl.linspace(100.0, 900.0, num=100)
Ediff_to_Eb = 0.5
nH = 1e5
r_grain = 0.1e-4

k11 = [get_k11(nu, Ediff_to_Eb * Eb_, N_sites) for Eb_ in Eb]
ke = [get_ke(nu, Eb_, T) for Eb_ in Eb]
ka = get_ka(nH, r_grain, T, C.m_p)

gH_HME = [get_gH_HME(H2, k_cos, k11[i], ke[i]) for i in range(len(k11))]
gH_RE =  [get_gH_RE(H2, k_cos, k11[i]) for i in range(len(k11))]
H_HME =  [(k_cos*H2 + ke[i]*gH_HME[i])/ka * nH for i in range(len(k11))]
H_RE =   [(k_cos*H2 + ke[i]*gH_RE[i] )/ka * nH for i in range(len(k11))]

#for i in xrange(len(k11)):
#  print i, Eb[i], H_HME[i], H_RE[i]

pl.loglog(Eb, gH_HME, label='HME')
pl.loglog(Eb, gH_RE, label='RE')
#pl.loglog(Eb, H_HME, label='HME')
#pl.loglog(Eb, H_RE, label='RE')
#pl.ylim((1e-6, 2.0))
pl.legend()
pl.show()
