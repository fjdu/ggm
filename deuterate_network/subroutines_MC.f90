      subroutine SaveMiscBefore (fU)
      use CMDT
      implicit none
      logical IsWordChar

      integer, parameter :: recordLen=99999
      integer i, fU, ios
      character FMTstr*128

      if (IsWordChar(fReactionsSave(1:1))) then
        write (*,'(/A)') 'Saving reactions...'
        CALL openFileSequentialWrite &
          (fU, trim(path)//trim(fReactionsSave), recordLen)
        write (FMTstr, FMT='("(", I2, "I8)")') nParticipants+2
        do i=1, nReactions
          write (UNIT=fU, FMT=FMTstr, IOSTAT=ios) &
            reactions(:, i), nRealReactants(i), nRealProducts(i)
        end do
        close (UNIT=fU, IOSTAT=ios, STATUS='KEEP')
      end if

      if (IsWordChar(fSpeciesSave(1:1))) then
        write (*,'(/A)') 'Saving species...'
        CALL openFileSequentialWrite &
          (fU, trim(path)//trim(fSpeciesSave), recordLen)
        write (FMTstr, FMT='("(A", I2, ", ", I2, "I4, ES14.5)")') &
          constLenNameSpecies, nElement
        do i=1, nSpecies
          write (UNIT=fU, FMT=FMTstr, IOSTAT=ios) &
            nameSpecies(i), SpeciesElements(:, i), massSpecies(i)
        end do
        close (UNIT=fU, IOSTAT=ios, STATUS='KEEP')
      end if

      if (IsWordChar(fRatesSave(1:1))) then
        write (*,'(/A)') 'Saving rates...'
        CALL openFileSequentialWrite &
          (fU, trim(path)//trim(fRatesSave), recordLen)
        write (FMTstr, FMT='("(I4, 2X, ", I2, "A", I2, &
          & ", ES15.4, ES11.2, F11.2, F11.2, I4)")') &
          nParticipants, constLenNameSpecies
        do i=1, nReactions
          write (UNIT=fU, FMT=FMTstr, &
            IOSTAT=ios) i, strReactants(1:nReactants, i), &
            strProducts(1:nProducts, i), &
            rates(i), dblABC(:,i), typeReac(i)
        end do
        close (UNIT=fU, IOSTAT=ios, STATUS='KEEP')
      end if
      
      end subroutine SaveMiscBefore




      subroutine initialize (fU, fileInitialize)
      ! Read all the configuration parameters
      use CMDT
      implicit none
      character(LEN=128) fileInitialize
      integer fU, ios

      namelist /PhysicalParameters/ &
        Temperature, n_H, rateCosIon, &
        Av, omega_Albedo, MeanMolWeight, &
        BarrierWidth, ratioGrainToH, ratioDustGasMass, &
        GrainRadius, GrainDensity, SitesDensity, &
        nSpecies_Est, nMCSteps, nMCRepeat, nRecordsSave, &
        GrainRadiusBase, stickCoeffChargeGrain

      namelist /Paths/ &
        path, fReactions, fMobility, fDissRecom, &
        fInitialCondition, &
        fReactionsSave, fSpeciesSave, fRatesSave, &
        fSaveFinalResult, fSaveElementalAb, fSaveALLHistory, &
        fSaveALLHistory_ASCII, fAnalyseResult, fPhyParHistory

      CALL openFileSequentialRead (fU, fileInitialize, 999)
      read (UNIT=fU, IOSTAT=ios, NML=PhysicalParameters)
      read (UNIT=fU, IOSTAT=ios, NML=Paths)
      close (UNIT=fU, IOSTAT=ios, STATUS='KEEP')

      if (nRecordsSave .GT. nMCSteps) nRecordsSave = nMCSteps

      end subroutine initialize
