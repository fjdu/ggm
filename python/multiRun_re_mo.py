#<timestamp>2012-03-26 Mon 18:29:37</timestamp>
import sys
import datetime
import os
import shutil
import stat
import string
import filecmp

time_start = datetime.datetime.now()
print '\tStarted at ', time_start

# ++++ Initialization begin ++++
pathAbsScript = os.path.dirname(os.path.abspath(sys.argv[0]))

if len(sys.argv) == 1:
    print '\tWarning: no configuration file is provided.'
    exit()
else:
    print '\tLoading config file ', sys.argv[1]
    #execfile(sys.argv[1])
    execfile(os.path.join(pathAbsScript, sys.argv[1]))

execfile(os.path.join(pathAbsScript, 'multiRunPreProCommon.py'))
execfile(os.path.join(pathAbsScript, 'MySubsoutines.py'))

# ---- Initialization end ----

# ++++ Loop over all the directories ++++
iCount_all = 0

CommonFilesDir = os.path.expanduser(CommonFilesDir)

for info2 in  outputInfo_s:
 for info1 in info2:
  for info in info1:
    #print '\t', info
    iCount_all += 1
    print '\t{0:5d} of {1:5d}'.format(iCount_all, nCount_all)
    if flag_use_old_moment == 0:
        outputdir_base = info['basepath']
        outputdir = info['path']
        n_H = info['nH']
        T_gas = info['T_gas']
        T_dust = info['T_dust']
        GrainRadius = info['GrainRadius']
        Av = funcAv(n_H)

        if not os.access(outputdir_base, os.F_OK):
            os.mkdir(outputdir_base)
        if not os.access(outputdir, os.F_OK):
            os.mkdir(outputdir)

        print '\tOutput files are in ' + outputdir

        exe_src  = os.path.join(CommonFilesDir, src_exe_dir, file_exe_moment)
        exe_dest = os.path.join(running_dir, file_exe_moment)
        if not os.access(exe_dest, os.F_OK):
          shutil.copyfile(exe_src, exe_dest)
          os.chmod(exe_dest, stat.S_IRWXU | stat.S_IWOTH)
        if not filecmp.cmp(exe_src, exe_dest):
          shutil.copyfile(exe_src, exe_dest)
          os.chmod(exe_dest, stat.S_IRWXU | stat.S_IWOTH)

        shutil.copyfile(os.path.join(CommonFilesDir, src_dat_dir, reaction_file),
            os.path.join(outputdir, reaction_file))

        shutil.copyfile(os.path.join(CommonFilesDir, src_dat_dir, enthalpy_file),
            os.path.join(outputdir, enthalpy_file))

        shutil.copyfile(os.path.join(CommonFilesDir, src_dat_dir, fSurfRadical),
            os.path.join(outputdir, fSurfRadical))

        if use_gas_steady_state_as_ini:
            T_for_ini = T_gas
            n_H_for_ini = n_H
            execfile(os.path.join(pathAbsScript, 'make_ini.py'))
        else:
            shutil.copyfile(os.path.join(CommonFilesDir, src_dat_dir, initial_condition),
              os.path.join(outputdir, initial_condition))

        fU = open(os.path.join(outputdir, config_file_moment), 'w')
        fU.writelines(string.join([\
              '&PhysicalParameters', 
              'T_gas = {0:8.2F}'.format(T_gas), 
              'T_dust = {0:8.2F}'.format(T_dust), 
              'n_H = {0:15.3E}'.format(n_H), 
              'rateCosIon = {0:15.6E}'.format(rateCosIon), 
              'Av = {0:15.3E}'.format(Av), 
              'omega_Albedo = {0:15.3E}'.format(omega_Albedo), 
              'ReactBarrierWidth = {0:15.3E}'.format(ReactBarrierWidth), 
              'DiffBarrierWidth = {0:15.3E}'.format(DiffBarrierWidth), 
              'MeanMolWeight = {0:15.6E}'.format(MeanMolWeight), 
              'ratioDustGasMass = {0:15.6E}'.format(ratioDustGasMass), 
              'GrainDensity = {0:15.6E}'.format(GrainDensity), 
              'GrainRadius = {0:15.6E}'.format(GrainRadius), 
              'SitesDensity = {0:15.6E}'.format(SitesDensity), 
              'CosmicDesorpGrainT = {0:15.3E}'.format(CosmicDesorpGrainT), 
              'CosmicDesorpPreFactor = {0:15.3E}'.format(CosmicDesorpPreFactor), 
              'ChemDesorpPreFactor = {0:15.3E}'.format(ChemDesorpPreFactor), 
              'SimpleChemDesorption = {0:8s}'.format(SimpleChemDesorption), 
              'Diff2DesorRatio = {0:15.3E}'.format(Diff2DesorRatio), 
              'nSpecies_Est = {0:12d}'.format(nSpecies_Est), 
              'sto_threshold = {0:15.6E}'.format(sto_threshold), 
              'rateThreshold = 0.0E0', 
              'nOrderLim = {0:6d}'.format(nOrderLim), 
              'tFinal = {0:15.3E}'.format(tFinal), 
              'nNumMMEst = {0:8d}'.format(nNumMMEst), 
              'allowExposureByReaction = {0:s}'.format(allowExposureByReaction), 
              'allowAlwaysExpose = {0:s}'.format(allowAlwaysExpose), 
              'useCosPhotoDiss = {0:s}'.format(useCosPhotoDiss), 
              'useZeroPointEnergy = {0:s}'.format(useZeroPointEnergy), 
              'useThreeBody = {0:s}'.format(useThreeBody), 
              'useThreeBodyPrefac = {0:s}'.format(useThreeBodyPrefac), 
              'numGrainLayerThreshold = {0:8.2F}'.format(numGrainLayerThreshold), 
              '/', 
              '&ODEParameters', 
              'nIteration = {0:6d}'.format(nIteration), 
              'RTOL = {0:15.6E}'.format(RTOL), 
              'ATOL = {0:15.6E}'.format(ATOL), 
              '/', 
              '&Paths', 
              'path = "' + outputdir + '"', 
              'fReactions = "' + reaction_file + '"', 
              'fSpeciesEnthalpy = "' + enthalpy_file + '"', 
              'fInitialCondition = "' + initial_condition + '"', 
              'fSurfRadical = "' + fSurfRadical + '"', 
              'fSaveALLHistory = "' +file_evolution_moment_bin + '"', 
              'fSaveALLHistory_ASCII = "' +file_evolution_moment_ascii + '"', 
              'fSaveElementalAb = "elemental_abundance_moment.dat"', 
              'fSpeciesSave = "species_saved.dat"', 
              'fRatesSave = "rates_saved.dat"', 
              'fNameMoments = "NameAllMoments.dat"', 
              'fOutLog = "Log.dat"', 
              'fPhyParHistory = "PhyParHistory.dat"', 
              'fAnalyseResult = "AnalyseResult.dat"', 
              '/'], '\n'))
        fU.close()

        command_moment = os.path.join(running_dir, file_exe_moment) + ' ' + \
                         os.path.join(outputdir, config_file_moment)
        print '\tRunning command:  ', command_moment
        if only_config == 0:
            os.system(command_moment)

time_end = datetime.datetime.now()
print '\tFinished at ', time_end
print '\tA time span of ', time_end - time_start, ' has elapsed.'
