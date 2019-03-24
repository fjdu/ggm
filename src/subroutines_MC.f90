      SUBROUTINE init_random_seed()
        INTEGER :: i, n, clock
        INTEGER, DIMENSION(:), ALLOCATABLE :: seed

        CALL RANDOM_SEED(size = n)
        ALLOCATE(seed(n))

        !write (*, '(A, I6)') 'Seed size: ', n

        CALL SYSTEM_CLOCK(COUNT=clock)

        seed = clock + 37 * (/ (i - 1, i = 1, n) /)
        CALL RANDOM_SEED(PUT = seed)

        DEALLOCATE(seed)
      END SUBROUTINE



      subroutine initialize (fU, fileInitialize)
      ! Read all the configuration parameters
      use CMDT
      implicit none
      character(LEN=128) fileInitialize
      integer fU, ios

      namelist /PhysicalParameters/ &
        T_gas, T_dust, n_H, &
        Av, omega_Albedo, rateCosIon, MeanMolWeight, &
        BarrierWidth, ratioGrainToH, ratioDustGasMass, &
        GrainRadius, GrainDensity, SitesDensity, &
        nSpecies_Est, nIteration, nMCRepeat, &
        CosmicDesorpGrainT, CosmicDesorpPreFactor, &
        ChemDesorpPreFactor, &
        Diff2DesorRatio, tFinal

      namelist /Paths/ &
        path, fReactions, fDissRecom, &
        fInitialCondition, &
        fReactionsSave, fSpeciesSave, fRatesSave, &
        fSaveFinalResult, fSaveElementalAb, fSaveALLHistory, &
        fSaveALLHistory_ASCII, fAnalyseResult, fPhyParHistory, &
        fSpeciesEnthalpy

      CALL openFileSequentialRead (fU, fileInitialize, 999)
      read (UNIT=fU, IOSTAT=ios, NML=PhysicalParameters)
      read (UNIT=fU, IOSTAT=ios, NML=Paths)
      close (UNIT=fU, IOSTAT=ios, STATUS='KEEP')

      end subroutine initialize
