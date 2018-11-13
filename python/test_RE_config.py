#<timestamp>2011-07-14 Thu 16:34:39</timestamp>

GrainRadius_s = [0.1E-6]
Temperature_s = [10.0, 20.0, 30.0]
n_H_s = [1.0E5, 1.0E6]

ChemDesorpPreFactor = 0.1

CommonFilesDir = '~/projects/DeuteriumGasGrain/streamline/'

output_prefix = 'run_{0:6G}'.format(ChemDesorpPreFactor).replace(' ', '')

reaction_file = 'network_test_HME.dat'
enthalpy_file = 'Species_enthalpy.dat'
initial_condition = 'initial_condition_test_HME.dat'
fSurfRadical = 'surface_radicals_test_RE.dat'

flag_use_old_moment = 0
only_config = 0
nNumMMEst = 1024

file_exe_moment = 'mo'
