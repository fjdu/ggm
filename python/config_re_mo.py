#<timestamp>2011-07-14 Thu 16:34:39</timestamp>

GrainRadius_s = [0.1E-6]
Temperature_s = [10.0]
n_H_s = [1e3, 1e4, 1e5, 1e6]

ChemDesorpPreFactor = 0.01

CommonFilesDir = '~/projects/DeuteriumGasGrain/streamline/'

output_prefix = 'HandD_{0:6G}'.format(ChemDesorpPreFactor).replace(' ', '')

reaction_file = 'rate06_dipole_reduced_20120513_50K_com_red_com.dat'
enthalpy_file = 'Species_enthalpy.dat'
initial_condition = 'initial_condition_molecular.dat'
#fSurfRadical = 'surface_radicals_HME.dat'
#fSurfRadical = 'surface_radicals_RE.dat'
fSurfRadical = 'surface_radicals_HME_HandD.dat'

flag_use_old_moment = 0
only_config = 0
nNumMMEst = 20000
#numGrainLayerThreshold = 1.0

allowExposureByReaction = '.TRUE.'
allowAlwaysExpose = '.FALSE.'
useCosPhotoDiss = '.FALSE.'
useZeroPointEnergy = '.TRUE.'
use_gas_steady_state_as_ini = True

useThreeBody = '.FALSE.'
useThreeBodyPrefac = '.TRUE.'

file_exe_moment = 'MO'
file_exe_moment_ini = 'RE'

if use_gas_steady_state_as_ini:
  reaction_file_for_ini = 'rate06_dipole_reduced_20120513_50K_com_red_com_for_ini.dat'
  init_file_for_gas = 'initial_condition_molecular.dat'
  file_ini_out = 'initial_condition_steady.dat'
  config_file_ini = 'config_for_ini.dat'
