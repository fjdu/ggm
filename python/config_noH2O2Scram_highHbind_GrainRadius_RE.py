#<timestamp>2011-07-14 Thu 16:34:39</timestamp>

#GrainRadius_s = [0.01e-6, 1e-6]
GrainRadius_s = [0.1e-6] #[0.001e-6, 10e-6]
Temperature_s = [10.0]
n_H_s = [1e5]

ChemDesorpPreFactor = 0.01

CommonFilesDir = '~/projects/DeuteriumGasGrain/streamline/'

output_prefix = 'RE_{0:6G}'.format(ChemDesorpPreFactor).replace(' ', '')

reaction_file = 'rate06_dipole_red_com_deut_20120711_noH2O2Scram_highHbind.dat'
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

useThreeBody = '.FALSE.'
useThreeBodyPrefac = '.TRUE.'

file_exe_moment = 'MO'
file_exe_moment_ini = 'RE'

if use_gas_steady_state_as_ini:
  reaction_file_for_ini = 'rate06_dipole_red_com_deut_20120711_highHbind_for_ini.dat'
  init_file_for_gas = 'initial_condition_molecular.dat'
  file_ini_out = 'initial_condition_steady.dat'
  config_file_ini = 'config_for_ini.dat'
