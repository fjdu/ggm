# 2011-05-05 Thu 19:44:41

def funcAv(n_H):
    #return 5.0
    return 15.0 # if n_H > 2E3 else 2.0

def funcRatioGrainToH(GrainRadius=1E-7):
    return ratioDustGasMass * (mProton * MeanMolWeight) \
      / (4.0E0*math.pi/3.0E0 * (GrainRadius)**3 *
          GrainDensity)

rateCosIon = 1.36E-17
omega_Albedo = 0.6 # Woodall et al. 2007 said the typical value should be 0.6.
ratioDustGasMass = 1.0E-2
mProton = 1.67262158E-27
MeanMolWeight = 1.4E0
GrainDensity = 2.0E3
SitesDensity = 1.0E15
ReactBarrierWidth = 1E-10
DiffBarrierWidth = 1E-10

CosmicDesorpGrainT = 70.0
CosmicDesorpPreFactor = 3.16E-19 # Hasegawa93, Eq 15
#CosmicDesorpPreFactor = 0.0E0
Diff2DesorRatio = 0.5

SimpleChemDesorption = '.TRUE.'

nSpecies_Est = 1024

#RTOL = 1.0E-10
RTOL = 1.0E-8
ATOL = 1.0E-35

tFinal = 1.0E7
nIteration = 2020
sto_threshold = 1.0
nOrderLim = 2

output_postfix = ''

running_dir = '.'
src_exe_dir = 'fortran'
src_dat_dir = 'network_for_use'

file_evolution_MC = 'evolution_MC_'
file_evolution_MC_bin = 'evolution_MC_bin_'
file_evolution_moment_ascii = 'evolution_moment_' + output_postfix + '_ascii.dat'
file_evolution_moment_bin = 'evolution_moment_' + output_postfix + '_bin.dat'

config_file_moment = 'config_moment_' + output_postfix + '.dat'

n_T_s = Temperature_s.__len__()
n_n_H_s = n_H_s.__len__()
n_GrainRadius_s = GrainRadius_s.__len__()

nCount_all = n_T_s * n_n_H_s * n_GrainRadius_s

outputInfo_s = [[[{'basepath':os.path.join('.', output_prefix),
                'path':os.path.join('.', output_prefix, 
                #output_prefix + '_' + \
                ('{0:6.1F}_{1:10.1E}_{2:10.1E}_{3:6.1F}_{4:10.1E}/'.
                 format(Temperature, n_H, GrainRadius, funcAv(n_H), ChemDesorpPreFactor)).
               replace(' ', '')), 'nH':n_H, 'T_gas':Temperature, 'T_dust':Temperature,
               'GrainRadius':GrainRadius}
            for n_H in n_H_s]
            for Temperature in Temperature_s]
            for GrainRadius in GrainRadius_s]
