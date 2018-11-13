#<timestamp>2012-05-14 Mon 00:56:32</timestamp>

# This modifies the value set by the original config file.
initial_condition = file_ini_out

RTOL_make_ini = 1E-11

time_start = datetime.datetime.now()
print '\tStarted at ', time_start

print '\tOutput files are in ' + outputdir

shutil.copyfile(os.path.join(CommonFilesDir, src_dat_dir, reaction_file_for_ini),
    os.path.join(outputdir, reaction_file_for_ini))

shutil.copyfile(os.path.join(CommonFilesDir, src_dat_dir, init_file_for_gas),
  os.path.join(outputdir, init_file_for_gas))

fU = open(os.path.join(outputdir, config_file_ini), 'w')
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
      '/',
      '&ODEParameters',
      'nIteration = {0:6d}'.format(nIteration),
      'RTOL = {0:15.6E}'.format(RTOL_make_ini),
      'ATOL = {0:15.6E}'.format(ATOL),
      '/',
      '&Paths',
      'path = "' + outputdir + '"',
      'fReactions = "' + reaction_file_for_ini + '"',
      'fSpeciesEnthalpy = "' + enthalpy_file + '"',
      'fInitialCondition = "' + init_file_for_gas + '"',
      'fSurfRadical = "' + fSurfRadical + '"',
      'fSaveFinalResult = "' + file_ini_out + '"',
      '/'], '\n'))
fU.close()

command_moment = os.path.join(running_dir, file_exe_moment_ini) + ' ' + \
                 os.path.join(outputdir, config_file_ini)
print 'Making the initial condition...'
print '\tRunning command:  ', command_moment
os.system(command_moment)
print 'Initial condition is ready.'

time_end = datetime.datetime.now()
print '\tFinished at ', time_end
print '\tA time span of ', time_end - time_start, ' has elapsed.'

