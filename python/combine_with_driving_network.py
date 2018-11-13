# This file is used to compare the gas phase driving network and an
# automatically generated gas phase network, and combine them together.
# file_deuterated = '../deuterate_network/rate06_dipole_reduced_20120404_50K_Deuterated.dat'
# file_driving = '../driving_network/gas_phase_deuteration_Roberts2000_2004.dat'
# file_combined = '../combine_network/rate06_dipole_reduced_20120404_50K_combined.dat'
#file_deuterated = '../reduce_network/rate06_dipole_reduced_20120426_50K_combined_red.dat'
#file_driving = '../reduce_network/rate06_dipole_reduced_20120426_50K_combined_red_1.dat'
#file_combined = '../reduce_network/rate06_dipole_reduced_20120426_50K_combined_red_com.dat'
#file_deuterated = '../reduce_network/rate06_dipole_reduced_20120426_50K_very.dat'
#file_driving = '../reduce_network/rate06_dipole_reduced_20120407_50K.dat'
#file_combined = '../reduce_network/rate06_dipole_reduced_20120407_50K_com.dat'
#file_deuterated = '../deuterate_network/rate06_dipole_reduced_20120407_50K_com_Deuterated.dat'
#file_driving = '../driving_network/gas_phase_deuteration_Roberts2000_2004.dat'
#file_combined = '../reduce_network/rate06_dipole_reduced_20120407_50K_com_deu.dat'
file_deuterated = '../reduce_network/rate06_dipole_reduced_20120426_50K_com_red_very.dat'
file_driving = '../reduce_network/rate06_dipole_reduced_20120426_50K_com_red_0.dat'
file_combined = '../reduce_network/rate06_dipole_reduced_20120426_50K_com_red_com.dat'

f_driving = open(file_driving, 'r')
f_deuterated = open(file_deuterated, 'r')
f_combined = open(file_combined, 'w')

str_driving = f_driving.readlines()
str_deuterated = f_deuterated.readlines()

for str_de in str_deuterated:
  flag = True
  if str_de[0] != '!':
    for str_dr in str_driving:
      if str_dr[0] == '!':
        continue
      s_dr_reac = sorted(str_dr[0:25].split())
      s_de_reac = sorted(str_de[0:25].split())
      s_dr_prod = sorted(str_dr[25:85].split())
      s_de_prod = sorted(str_de[25:85].split())
      if s_dr_reac == s_de_reac and s_dr_prod == s_de_prod:
        flag = False
       #print str_dr, str_de
        break
  if flag:
    f_combined.write(str_de)

for str_dr in str_driving:
  f_combined.write(str_dr)

f_deuterated.close()
f_driving.close()
f_combined.close()
