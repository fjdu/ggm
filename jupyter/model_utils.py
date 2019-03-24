import numpy as np


def load_time_evol_data(fname):
    with open(fname, 'r') as f:
        headerline = f.readline()
    header = headerline.split()
    arr = np.loadtxt(fname, skiprows=1)    
    return {header[i]: arr[:, i] for i in range(len(header))}


def dust2gas(r_dust_micron=0.1,
             dust2gas_mass=0.01, mean_mol_weight=1.4, dust_density=2.0):
    m_proton = 1.67e-24
    m_dust = 4*np.pi/3 * (r_dust_micron*1e-4)**3 * dust_density
    return dust2gas_mass * (m_proton * mean_mol_weight) / m_dust


def getMonolayers(numMantle, r_dust_micron=0.1, N_Site=1e15):
    surfAread = 4*np.pi * (r_dust_micron*1e-4)**2
    return numMantle / (surfAread * N_Site)


def getNSites(r_dust_micron=0.1, N_Site=1e15):
    surfAread = 4*np.pi * (r_dust_micron*1e-4)**2
    return surfAread * N_Site


def make_initial_condition_CO_H2O_ice(d2g, Nsites, C_ab=1.4e-4, O_ab=3.2e-4):
    abun = {
      'H2': 5.00E-1,
      'HD': 2.00E-5,
      'He': 0.09E-0,
      'N':  7.50E-5,
      'S':  8.00E-8,
      'Si': 8.00E-9,
      'Na': 2.00E-8,
      'Mg': 7.00E-9,
      'Fe': 3.00E-9,
      'P':  3.00E-9,
      'F':  2.00E-8,
      'Cl': 4.00E-9,
    }
    CO_ab = C_ab
    H2O_ab = max(0.0, O_ab - CO_ab)
    nCO_perdust = CO_ab / d2g
    nH2O_perdust = H2O_ab / d2g
    nCO_surf = min(Nsites, nCO_perdust)
    nH2O_surf = max(Nsites - nCO_surf, 0)
    nCO_mantle = max(nCO_perdust - nCO_surf, 0)
    nH2O_mantle = max(nH2O_perdust - nH2O_surf, 0)
    abun['gCO'] = nCO_surf * d2g
    abun['mCO'] = nCO_mantle * d2g
    abun['gH2O'] = nH2O_surf * d2g
    abun['mH2O'] = nH2O_mantle * d2g
    return abun


def write_initial_condition(fname, d_ab):
    with open(fname, 'w') as f:
        for k in d_ab:
            f.write('{k:12s}{v:7.2e}\n'.format(k=k, v=d_ab[k]))
    return
