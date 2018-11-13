import os
import string

rate06_exe = 'rate06_for_ini'
config_for_init = 'config_ini.dat'
reaction_file_for_ini = 'rate06_dipole_reduced_20120426_50K_com_red_com.dat'
init_file_for_gas = 'initial_condition_Garrod08.dat'
file_ini_out = 'initial_condition_steady.dat'
dir_ini_out = os.path.join(CommonFilesDir,  'initial_condition/')
out_absorb_file = 'tmp.dat'

fU = open(os.path.join(dir_ini_out, config_for_init), 'w')
fU.writelines(string.join([\
      '&PhysicalParameters',
      'Temperature = {0:8.2F}'.format(T_for_ini),
      'n_H = {0:15.3E}'.format(n_H_for_ini),
      'Av = 15.0D0',
      'omega_Albedo = 0.5D0',
      'rateCosIon = 3.0D-17',
      'rateHHH2 = 3.003953e-18',
      'ratioDustGasMass = 1.3D-2',
      'GrainDensity = 2.0D0',
      'GrainRadius = 1.0D-1',
      'GrainRadiusBase = 0.02',
      'aGrainMin = 1.0D-5',
      'aGrainMax = 1.0D-5',
      'timeUnit = 1.0D2',
      'nSpecies_Est = 512',
      'rateThreshold = 0.0D0',
      '/',
      '&ODEParameters',
      'nIteration = 750',
      'RTOL = 1.0D-6',
      'ATOL = 1.0D-90',
      '/',
      '&PlotParameters',
      'yPGmax = -3.0',
      'yPGmin = -13.5',
      '/',
      '&Paths',
      'path = "{0:s}"'.format(dir_ini_out),
      'fReactions = "{0:s}"'.format(reaction_file_for_ini),
      'fInitialCondition = "{0:s}"'.format(init_file_for_gas),
      'fDissRecom = "Pagani2009.dat"',
      'fPlotList = "saveHistoryList_2.txt"',
      'fReactionsSave = "reactions_saved.dat"',
      'fSpeciesSave = "species_saved.dat"',
      'fRatesSave = "rates_saved.dat"',
      'fSaveFinalResult = "{0:s}"'.format(file_ini_out),
      'fSaveALLHistory = "evolveHistory.dat"',
      'fSaveElementalAb = "elemental_abundance.txt"',
      'fPGout = "ChemEvol.ps/CPS"',
      'fAnalyseResult = "AnalyseResult.txt"',
      '/'], '\n'))
fU.close()

command_run = os.path.join(dir_ini_out, rate06_exe) + ' ' \
            + os.path.join(dir_ini_out, config_for_init \
            + ' > ' \
            + os.path.join(dir_ini_out, out_absorb_file))
print 'Making the initial condition...'
os.system(command_run)
print 'Initial condition is ready.'
