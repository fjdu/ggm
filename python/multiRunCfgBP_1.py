#<timestamp>2011-07-14 Thu 16:34:39</timestamp>

GrainRadius_s = [0.1E-6]
Temperature_s = [15.0, 30.0]
n_H_s = [1.0E4]

ChemDesorpPreFactor = 0.1

CommonFilesDir = '/homes/fjdu/CommonFiles/H2O2_model_20111010_D/'

output_prefix = 'run_{0:6G}'.format(ChemDesorpPreFactor).replace(' ', '')

reaction_file = 'rate06_dipole_reformated_reduced_20111010_D.dat'
enthalpy_file = 'Species_enthalpy.dat'
initial_condition = 'initial_condition_Garrod08.dat'

flag_use_old_moment = 0
only_config = 0

file_exe_moment = 'mo'
