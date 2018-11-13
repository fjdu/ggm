#<timestamp>2011-07-14 Thu 16:34:39</timestamp>

GrainRadius_s = [0.1E-6]
Temperature_s = [10.0, 15.0, 20.0, 30.0]
n_H_s = [1.0E4, 1.0E5, 1.0E6]

ChemDesorpPreFactor = 0.01

CommonFilesDir = '~/projects/DeuteriumGasGrain/streamline/'

output_prefix = 'run_{0:6G}'.format(ChemDesorpPreFactor).replace(' ', '')

reaction_file = 'rate06_dipole_red_50K_deut_mant_20120424_no_negative.dat'
enthalpy_file = 'Species_enthalpy.dat'
initial_condition = 'initial_condition.dat'
fSurfRadical = 'surface_radicals_HME.dat'

flag_use_old_moment = 0
only_config = 0
nNumMMEst = 20000

useCosPhotoDiss = '.FALSE.'
useZeroPointEnergy = '.TRUE.'
use_gas_steady_state_as_ini = True

file_exe_moment = 'mo'
