#<timestamp>2011-07-14 Thu 16:34:39</timestamp>

GrainRadius_s = [0.1E-6]
Temperature_s = [15.0]
n_H_s = [1e3, 1e4, 1e5, 1e6]

ChemDesorpPreFactor = 0.01

CommonFilesDir = '~/projects/DeuteriumGasGrain/streamline/'

output_prefix = 'run_{0:6G}'.format(ChemDesorpPreFactor).replace(' ', '')

reaction_file = 'rate06_dipole_reduced_20120522_50K_com_red_com_noH2OH.dat'
enthalpy_file = 'Species_enthalpy.dat'
initial_condition = 'initial_condition_molecular.dat'
fSurfRadical = 'surface_radicals_RE.dat'

flag_use_old_moment = 0
only_config = 0
nNumMMEst = 3000

allowExposureByReaction = '.TRUE.'
allowAlwaysExpose = '.FALSE.'
useCosPhotoDiss = '.FALSE.'
useZeroPointEnergy = '.TRUE.'
use_gas_steady_state_as_ini = True

useThreeBody = '.TRUE.'
useThreeBodyPrefac = '.TRUE.'

file_exe_moment = 'RE'

if use_gas_steady_state_as_ini:
  reaction_file_for_ini = 'rate06_dipole_reduced_20120513_50K_com_red_com_for_ini.dat'
  init_file_for_gas = 'initial_condition_molecular.dat'
  file_ini_out = 'initial_condition_steady.dat'
  config_file_ini = 'config_for_ini.dat'
