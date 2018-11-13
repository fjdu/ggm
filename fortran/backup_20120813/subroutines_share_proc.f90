      subroutine CalcRates
      use CMDT
      implicit none
      integer i
      double precision dblTmp, T_gas_300

      T_gas_300 = T_gas / 300.0D0
      ! T_gas_Reduced = kBoltzmann * T_gas &
      !   / (elementaryCharge**2 * coulombConst / (GrainRadius*1D-6))
      ! JNegaPosi = (1D0 + 1D0/T_gasReduced) &
      !           * (1D0 + SQRT(2D0/(2D0+T_gasReduced)))
      ! JChargeNeut = (1D0 + SQRT(DPI/2D0/T_gasReduced))

      ! Set the default value of several parameters.
      if (ratioGrainToH .LT. tiny(0D0)) &
        ratioGrainToH = & ! ratioGrainToH := N_Grain / N_H
          ratioDustGasMass * (mProton * MeanMolWeight) &
          / (4.0D0*DPI/3.0D0 * (GrainRadius)**3 * GrainDensity)
      if (Diff2DesorRatio .LT. tiny(0D0)) &
          Diff2DesorRatio = 0.3D0
      if (ChemDesorpPreFactor .LT. tiny(0D0)) &  ! The default
          ChemDesorpPreFactor = 0.0D0 ! is to disallow chem-desorption.

      SitesPerGrain = 4D0*DPI * GrainRadius**2 * (SitesDensity*1D4)
      VolumnPerGrain = 1D0 / (ratioGrainToH * n_H) ! in cm3

      CALL CalcMobilities
      CALL CalcBranchRatioOnGrain

      do i=1, nReactions
        if ((nRealReactants(i) .LE. 0) .OR. &
            (nRealProducts(i) .LE. 0)) then
          rates(i) = 0.0D0
          cycle
        end if
        if (.NOT. useCosPhotoDiss) then
          if ((typeReac(i) .GE. 71) .AND. (typeReac(i) .LE. 73)) then
            rates(i) = 0D0
            cycle
          end if
        end if
        select case (typeReac(i))
        case (5, 53) ! Two body
          rates(i) = dblABC(1, i) * (T_gas_300**dblABC(2, i)) &
              * exp(-dblABC(3, i)/T_gas) &
              / VolumnPerGrain
        ! <timestamp>2011-05-25 Wed 23:31:21</timestamp>
        case (1, 71) ! One body; cosmic ray ionization
          rates(i) = dblABC(1, i) * (rateCosIon / rateCosIon_base)
        case (2, 72) ! One body; cosmic ray induced ionization
          rates(i) = dblABC(1, i) * (T_gas_300**dblABC(2, i)) &
              * dblABC(3, i) / (1-omega_Albedo) &
              * (rateCosIon / rateCosIon_base)
        case (3, 73) ! One body; photo-ionization
          rates(i) = dblABC(1, i) * exp(-dblABC(3, i) * Av)
        case (61) ! Adsorption
          ! rates(i) * Population =  number of i molecules 
          !   accreted per grain per unit time
          ! Pi * r**2 * V * n
          ! Also take into account the possible effect of
          !   stick coefficient and other temperature dependence 
          !   (e.g., coloumb focus).
          rates(i) = dblABC(1,i) * T_gas**dblABC(2,i) &
            * exp(-dblABC(3,i)/T_gas) &
            * DPI * GrainRadius**2 &
            * sqrt(8D0/DPI*kBoltzmann*T_gas &
            / (massSpecies(reactions(1, i)) * mProton)) &
            / (VolumnPerGrain*1D-6) ! VolumnPerGrain is in cm3.
        case (62) ! Desorption
          ! <timestamp>2011-06-10 Fri 18:07:51</timestamp>
          !     A serious typo is corrected.
          rates(i) = vibFreqSpecies(reactions(1,i)) * dblABC(1,i) * ( &
              T_dust**dblABC(2,i) * exp(-dblABC(3,i)/T_dust) &
            ! Cosmic ray desorption rate from Hasegawa1993.
            + CosmicDesorpPreFactor * CosmicDesorpGrainT**dblABC(2,i) &
                * exp(-dblABC(3,i)/CosmicDesorpGrainT))
            !  Quantum tunneling not included (people never talk about the
            !  possibility of desorption through chemical tunneling).
            !  The value of BarrierWidth is uncertain.
            !  Notice the difference bwtween hbarPlanck and hPlanck.
        case (63) ! A + A -> xxx
          rates(i) = mobilitySpecies(reactions(1,i)) &
            / SitesPerGrain * BranchingRatios(i)
        case (64) ! A + B -> xxx
          rates(i) = &
            (mobilitySpecies(reactions(1,i)) &
             + mobilitySpecies(reactions(2,i))) &
            / SitesPerGrain * BranchingRatios(i)
        case (65, 66)
          cycle
        case default
          rates(i) = 0D0
        end select
      end do
      end subroutine CalcRates



      subroutine CalcBranchRatioOnGrain
      use CMDT
      use CMTP
      implicit none
      integer i, j, i1
      double precision sumWeights
      double precision CalcWeightSlave
      ! BranchingRatios is defined for all the reactions in the network.
      ! But it is not needed for in the case of gas phase reactions,
      ! so a default value of 1 is assigned to it.
      BranchingRatios = 1D0
      do i=1, nGrReactions
        if (BranchInfoOnGrain(i)%nSlave .LT. 0) then
          ! The branching ratio will be set to the default (0.1), or will be
          ! set by the master reaction.
          cycle
        else
          sumWeights = 0D0
          do j=1, BranchInfoOnGrain(i)%nSlave
            ! i1 is the absolute index (i.e. its index in the whole network)
            ! of a slave reaction.
            i1 = BranchInfoOnGrain(i)%idxSlave(j)
            BranchInfoOnGrain(i)%WeightSlave(j) = CalcWeightSlave(i1)
            sumWeights = sumWeights + &
                BranchInfoOnGrain(i)%WeightSlave(j)
          end do
          ! For example, if the two branching ratios are 0.5 and 1.5, then they
          ! are modified into 0.25 and 0.75.
          ! However, if the two branching ratios are 0.4 and 0.5, then they
          ! will not be modified.
          if (sumWeights .LT. 1D0) sumWeights = 1D0
          do j=1, BranchInfoOnGrain(i)%nSlave
            i1 = BranchInfoOnGrain(i)%idxSlave(j)
            BranchingRatios(i1) = &
              BranchInfoOnGrain(i)%WeightSlave(j) / sumWeights
          end do
        end if
      end do
      end subroutine CalcBranchRatioOnGrain




      subroutine MakeBranchInfo
      use CMDT
      use CMTP
      implicit none
      integer i, j, ii, ij
      integer, dimension(8) :: idxTmp
      allocate (BranchInfoOnGrain(nGrReactions))
      do i=1, nGrReactions
        BranchInfoOnGrain(i)%nSlave = 1
      end do
      do i=1, nGrReactions
        if (BranchInfoOnGrain(i)%nSlave .EQ. -1) cycle
        ii = idxGrReactions(i)
        idxTmp(1) = ii
        do j=i+1, nGrReactions
          ij = idxGrReactions(j)
          ! <timestamp>2012-04-11 Wed 01:02:06</timestamp>
          ! Correct a bug: the branching ratio should only be calculated for
          ! those with the sam reaction type.
          if ((sum(abs(reactions(1:2,ii) - reactions(1:2,ij))) .EQ. 0) &
              .AND. (typeReac(ii) .EQ. typeReac(ij))) then
            BranchInfoOnGrain(i)%nSlave = &
              BranchInfoOnGrain(i)%nSlave + 1
            idxTmp(BranchInfoOnGrain(i)%nSlave) = ij
            BranchInfoOnGrain(j)%nSlave = -1
          end if
        end do
        if (BranchInfoOnGrain(i)%nSlave .GT. 0) then
          allocate (&
            BranchInfoOnGrain(i)%idxSlave( &
              BranchInfoOnGrain(i)%nSlave), &
            BranchInfoOnGrain(i)%WeightSlave( &
              BranchInfoOnGrain(i)%nSlave))
          BranchInfoOnGrain(i)%idxSlave = &
            idxTmp(1:BranchInfoOnGrain(i)%nSlave)
        end if
      end do
      end subroutine MakeBranchInfo



      double precision function CalcWeightSlave (idxReac)
      use CMDT
      use CMTP
      implicit none
      integer, intent(in) :: idxReac
      integer i, DoVibF
      double precision CalcWeightSlave, E_reac, E_desorb, &
        dblTmp
      integer getDoVibFSpecies
      select case (typeReac(idxReac))
      case (63, 64, 67)
        ! Take the maximum between the thermal hopping and quantum tunneling
        ! rates.
        if (dblABC(3, idxReac) .NE. 0D0) then
          CalcWeightSlave = dblABC(1, idxReac) * exp(max( &
            -dblABC(3, idxReac) / T_dust, &
            -2D0 * dblABC(2, idxReac) * 1D-10 / hbarPlanck * &
              sqrt(2D0 * tunneling_mass(idxReac) * mProton &
                * kBoltzmann * dblABC(3, idxReac))))
        else
          CalcWeightSlave = dblABC(1, idxReac)
        end if
        if (IsGas(reactions(1+nReactants, idxReac))) then
          ! Chemical desorption
          if (.NOT. SimpleChemDesorption) then
            E_desorb = 0D0
            DoVibF = 0
            E_reac = enthalpySpecies(reactions(1, idxReac)) + &
              enthalpySpecies(reactions(2, idxReac))
            do i=1, nRealProducts(idxReac)
              E_reac = E_reac - enthalpySpecies(reactions(i+nReactants, &
                idxReac))
              E_desorb = E_desorb + E_desorp_species( &
                reactions(i+nReactants, idxReac))
              DoVibF = DoVibF + &
                getDoVibFSpecies(reactions(i+nReactants, idxReac))
            end do
            if (IsNaN(E_reac) .OR. IsNaN(E_desorb)) then
              CalcWeightSlave = 0D0
              return
            end if
            if (E_desorb .GT. E_reac) then
              write (*,*) "In CalcWeightSlave:"
              write (*,*) "Something might be wrong!"
              write (*,*) idxReac, strReactants(1:2, idxReac), &
                strProducts(1:2, idxReac), &
                E_desorb, E_reac, &
                E_desorp_species(reactions(1+nReactants, idxReac))
              stop
            end if
            dblTmp = 1D0 - E_desorb/E_reac
            do i=1, DoVibF-1
              CalcWeightSlave = CalcWeightSlave * dblTmp
            end do
            ! write (*,*) idxReac, strReactants(1:2, idxReac), dblTmp, &
            !   DoVibF, dblTmp**(DoVibF-1)
          end if
          CalcWeightSlave = CalcWeightSlave * ChemDesorpPreFactor
        end if
      case default
        CalcWeightSlave = 1D0
      end select
      return
      end function CalcWeightSlave




      integer function getDoVibFSpecies(idxSpecies)
      ! Get degree of vibrational freedom.
      use CMDT
      implicit none
      integer, intent(in) :: idxSpecies
      integer getDoVibFSpecies
      getDoVibFSpecies = sum(SpeciesElements(4:nElement, idxSpecies))
      ! For number of atoms =
      !     1: DoVibF = 1
      !     2: DoVibF = 2
      !   n>2: DoVibF = 3*n-5
      if (getDoVibFSpecies .GT. 2) then
        getDoVibFSpecies = 3*getDoVibFSpecies - 5
      end if
      return
      end function getDoVibFSpecies




      subroutine ImportSpeciesEnthalpy
      use CMDT
      implicit none
      integer i, j, i1, nLineAll, fU, ios
      double precision dblTmp, dblNaN
      character(Len=128) FMTstr, strTMP
      character commentChar
      character(LEN=constLenNameSpecies) nameSpecies_tmp
      logical IsWordChar, getFileUnit
      ! The output enthalpies are in K.
      enthalpySpecies = dblNaN()
      i1 = 0
      commentChar = '!'
      if (IsWordChar(fSpeciesEnthalpy(1:1))) then
        if (.NOT. getFileUnit(fU)) then
          write (*,*) 'In subroutine ImportSpeciesEnthalpy:'
          write (*,*) 'Cannot allocate an output file unit!'
          stop
        end if
        CALL openFileSequentialRead &
            (fU, trim(path)//trim(fSpeciesEnthalpy), 99999)
        CALL GetNLineF (fU, nLineAll, nMobilityList, commentChar)
        rewind (UNIT=fU, IOSTAT=ios)
        write (FMTstr, FMT=&
          '("(", "A", I2, ", F", I1, ".0)")') &
          constLenNameSpecies, 9
        do i=1, nLineAll
          read (UNIT=fU, FMT='(A128)', IOSTAT=ios) strTMP
          if ((strTMP(1:1) .EQ. commentChar) .OR. &
              (strTMP(1:1) .EQ. ' ')) cycle
          read (strTMP, FMT=FMTstr, IOSTAT=ios) &
            nameSpecies_tmp, dblTmp
          if (ios .NE. 0) then
            write (*, *) 'Error in importing enthalpies: ios = ', ios
            stop
          end if
          do j=1, nSpecies
            if (trim(nameSpecies(j)) .EQ. &
              trim(nameSpecies_tmp)) then
              ! Convert from kJ/mol to K.
              enthalpySpecies(j) = dblTmp * 1D3 / IdealGasConst
              i1 = i1 + 1
              exit
            end if
          end do
        end do
        close (UNIT=fU, IOSTAT=ios, STATUS='KEEP')
        j = 0
        do i=1, nGrSpecies
          if (IsNan(enthalpySpecies(idxGrSpecies(i)))) then
            j = j + 1
           !write (*,*) 'In ImportSpeciesEnthalpy:  '
           !write (*,*) nameSpecies(idxGrSpecies(i)), &
           !    ' does not have an enthalpy!!!!!!'
            !stop
          end if
        end do
        if (j .GT. 0) then
          write (*,*) j, ' species do not have enthalpies!'
        end if
      end if
      end subroutine ImportSpeciesEnthalpy




      subroutine CalcMobilities
      use CMDT
      implicit none
      integer i, i1, i2
      ! Species without a diffusion barrier provided is considered to be
      ! immobile.
      mobilitySpecies = 0D0
      do i=1, nGrSpecies
        i1 = idxGrSpecies(i)
        if (IsNaN(E_desorp_species(i1))) cycle
        ! The greater of the thermal hopping rates and quantum tunneling
        ! rates are selected.
        ! Note that according to Hasegawa et al. 1992, the energy appears in
        ! the following formula is the desorption energy barrier, not the
        ! diffusion energy barrier.
        ! Only H, H2, and D are allowed to hop through quantum tunneling.
        !if (massSpecies(i1) .GT. 2D0) then
        !  mobilitySpecies(i1) = &
        !    vibFreqSpecies(i1) * exp( &
        !      -E_desorp_species(i1) * Diff2DesorRatio / T_dust)
        !else
          mobilitySpecies(i1) = &
            vibFreqSpecies(i1) * exp(max( &
              -E_desorp_species(i1) * Diff2DesorRatio / T_dust, &
              -2D0 * DiffBarrierWidth / hbarPlanck * &
              sqrt(2D0 * massSpecies(i1) * mProton &
                * kBoltzmann * E_desorp_species(i1) * Diff2DesorRatio)))
        !end if
      end do
      end subroutine CalcMobilities




      subroutine getDesorEnergy
      use CMDT
      implicit none
      integer i, i1
      double precision dblNaN
      ! Set the energy barriers to a high value by default.
      E_desorp_species = dblNaN()
      do i=1, nGrReactions
        i1 = idxGrReactions(i)
        if (typeReac(i1) .EQ. 62) then
          if (useZeroPointEnergy) then
            dblABC(3, i1) = max(0D0, dblABC(3, i1) - &
              hPlanck/kBoltzmann/2D0 * vibFreqSpecies(Reactions(1, i1)))
          end if
          E_desorp_species(Reactions(1, i1)) = dblABC(3, i1)
          E_desorp_species(Reactions(nReactants+1, i1)) &
            = dblABC(3, i1)
        end if
      end do
      end subroutine getDesorEnergy



      subroutine getVibFreq
      use CMDT
      implicit none
      integer i, i1
      vibFreqSpecies = 0D0
      do i=1, nGrReactions
        i1 = idxGrReactions(i)
        if (typeReac(i1) .EQ. 62) then
          vibFreqSpecies(reactions(1,i1)) = &
            sqrt(2D0 * (SitesDensity * 1D4) * &
               kBoltzmann * dblABC(3,i1) / (DPI*DPI) / &
               (mProton*massSpecies(reactions(1,i1))))
        end if
      end do
      end subroutine getVibFreq



      subroutine getSurfGainLossReac
      use CMDT
      implicit none
      integer i, i1, j, nMantReac, nMantProd
      logical flag
      n_surf_mant_inc = 0
      n_surf_mant_dec = 0
      do i=1, nGrReactions
        i1 = idxGrReactions(i)
        if ((typeReac(i1) .EQ. 65) .OR. &
            (typeReac(i1) .EQ. 66) .OR. &
            (typeReac(i1) .EQ. 61) .OR. &
            (typeReac(i1) .EQ. 62)) cycle
        flag = .TRUE.
        do j=1, nRealReactants(i1)
          if (IsSurfWithMan(reactions(j, i1))) then
            flag = .FALSE.
            exit
          end if
        end do
        do j=1, nRealProducts(i1)
          if (IsSurfWithMan(reactions(j+nReactants, i1))) then
            flag = .FALSE.
            exit
          end if
        end do
        if (flag) cycle
        call get_num_surfMant(i1, nMantReac, nMantProd) 
        if (nMantReac .EQ. nMantProd) then
          cycle
        else if (nMantReac .LT. nMantProd) then
          n_surf_mant_inc = n_surf_mant_inc + 1
          surf_mant_inc(n_surf_mant_inc)%idx = i1
          surf_mant_inc(n_surf_mant_inc)%gain_loss = &
            dble(nMantProd - nMantReac)
        else
          n_surf_mant_dec = n_surf_mant_dec + 1
          surf_mant_dec(n_surf_mant_dec)%idx = i1
          surf_mant_dec(n_surf_mant_dec)%gain_loss = &
            dble(nMantReac - nMantProd)
        end if
      end do
      end subroutine getSurfGainLossReac



      subroutine get_num_surfMant (iReac, nMantReac, nMantProd)
      use CMDT
      implicit none
      integer i, j, iReac, nMantReac, nMantProd
      nMantReac = 0
      nMantProd = 0
      do i=1, nRealReactants(iReac)
        if (IsSurfWithMan(reactions(i, iReac))) &
          nMantReac = nMantReac + 1
      end do
      do i=1, nRealProducts(iReac)
        if (IsSurfWithMan(reactions(i+nReactants, iReac))) &
          nMantProd = nMantProd + 1
      end do
      end subroutine get_num_surfMant



      subroutine ReadReactions (fU, nLineAll, commentChar)
      ! Read the reaction file.
      use CMDT
      implicit none
      integer i, j, k, fU, ios, nLineAll
      character(LEN=256) strtmp
      character(LEN=64) FMTstring
      character commentChar

      nRealReactants = 0
      nRealProducts = 0

      write (FMTstring, FMT=&
        '("(", I1, "A", I2, ", ", I1, "E", I1, ".", I1, ", F", I1, &
        &".1, 6X, ", "I", I1, ")")') &
        nParticipants, constLenNameSpecies, 3, 9, 2, 6, 3

      rewind (UNIT=fU, IOSTAT=ios)

      i = 1

      do
        read (UNIT=fU, FMT='(A256)', IOSTAT=ios) strtmp
        if (ios .LT. 0) exit
        ! Commented or empty lines are treated the same.
        if ((strtmp(1:1) .EQ. commentChar) .OR. &
            (strtmp(1:1) .EQ. ' ')) &
          cycle
        do j=1, constLenNameSpecies*nParticipants-1
          if ((strtmp(j:j) .EQ. ' ') .AND. &
              (strtmp(j+1:j+1) .NE. ' ')) then
            if (mod(j, constLenNameSpecies) .NE. 0) then
              write(*,*) 'Warning: a typo in the reaction file?'
              write(*,*) i, strtmp
            end if
          end if
        end do
        read (UNIT=strtmp, FMT=FMTstring, IOSTAT=ios) &
          strReactants(:,i), strProducts(:,i), dblABC(:,i), &
          tunneling_mass(i),  typeReac(i)
        do j=1, nReactants
          do k=1, constLenNameSpecies
            if (strReactants(j, i)(k:k) .NE. ' ') then
              nRealReactants(i) = nRealReactants(i) + 1
              exit
            end if
          end do
        end do
        do j=1, nProducts
          do k=1, constLenNameSpecies
            if (strProducts(j, i)(k:k) .NE. ' ') then
              nRealProducts(i) = nRealProducts(i) + 1
              exit
            end if
          end do
        end do
        i = i + 1
      end do
      end subroutine ReadReactions




      subroutine MakeReactions
      use CMDT
      implicit none

      integer i, j, k
      logical flag

      reactions = 0

      nameSpecies(1) = strReactants(1, 1)
      nSpecies = 1
      do i=1, nReactions
        do k=1, nRealReactants(i)
          flag = .TRUE.
          do j=1, nSpecies
            if (trim(nameSpecies(j)) .EQ. trim(strReactants(k, i))) then
              flag = .FALSE.
              reactions(k, i) = j
              exit
            end if
          end do
          if (flag) then
            nSpecies = nSpecies + 1
            nameSpecies(nSpecies) = strReactants(k, i)
            reactions(k, i) = nSpecies
          end if
        end do
        do k=1, nRealProducts(i)
          flag = .TRUE.
          do j=1, nSpecies
            if (trim(nameSpecies(j)) .EQ. trim(strProducts(k, i))) then
              flag = .FALSE.
              reactions(k+nReactants, i) = j
              exit
            end if
          end do
          if (flag) then
            nSpecies = nSpecies + 1
            nameSpecies(nSpecies) = strProducts(k, i)
            reactions(k+nReactants, i) = nSpecies
          end if
        end do
      end do
      end subroutine MakeReactions



      logical function getPhyPar (t)
      ! Input
      !         t:  time. The unit should be in year.
      ! Output
      !         getPhyPar:
      !              .TRUE.: at least one physical parameter has changed
      !             .FALSE.: no physical parameter has changed
      ! If there are some parameters don't change with time, then
      !   just ignore them.
      use CMDT
      implicit none
      logical getPhyPar
      double precision t, T_gas_Tmp
      double precision, parameter :: t0 = 3.0D5, t1 = 3.001D5
      double precision, parameter :: &
        T_gas0 = 1.0D1, T_gas1 = 6.0D1, ratioChange=1D-1
      T_gas_Tmp = T_gas
      !
      if (t .LT. t0) then
        T_gas = T_gas0
      else if (t .LT. t1) then
        T_gas = T_gas0 + &
          (t-t0)/(t1-t0) * (T_gas1 - T_gas0)
      else
        T_gas = T_gas1
      end if
      !
      if (abs(T_gas_Tmp-T_gas) .GE. &
        T_gas_Tmp*ratioChange) then
        getPhyPar = .TRUE.
      else
        T_gas = T_gas_Tmp
        getPhyPar = .FALSE.
      end if
      return
      end function getPhyPar


      logical function ExternalAction (t, y, NEQ)
      ! Some external actions which modify the system status
      ! Input
      !         t:  time. The unit should be in year.
      ! Output
      !         ExternalAction:
      !              .TRUE.: something important has changed
      !             .FALSE.: nothing important has changed
      use CMDT
      implicit none
      logical ExternalAction
      double precision, parameter :: t0 = 3.0D5, t1 = 3.03D5
      integer i, i1, NEQ
      double precision t, y(NEQ)
      if ((t .GE. t0) .AND. (t .LT. t1)) then
        do i=1, nGrReactions
          i1 = idxGrReactions(i)
          if (typeReac(i1) .EQ. 62) then
            y(reactions(1+nReactants, i1)) = &
              y(reactions(1+nReactants, i1)) + y(reactions(1, i1))
            y(reactions(1, i1)) = 0D0
          end if
        end do
        do i=nSpecies+1, NEQ
          y(i) = 0D0
        end do
        ExternalAction = .TRUE.
      else
        ExternalAction = .FALSE.
      end if
      return
      end function ExternalAction



      !  Save miscellaneous preparatory informations.
      subroutine SaveMiscBefore
      use CMDT
      implicit none
      logical IsWordChar, getFileUnit
      integer, parameter :: recordLen=99999
      integer i, fU, ios
      character FMTstr*128

      if (.NOT. getFileUnit(fU)) then
        write (*,*) 'Cannot allocate an output file unit!'
        stop
      end if
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
        write (FMTstr, '("(A", I2, ", ", I2, "I4, ES14.5, 3ES11.2)")') &
          constLenNameSpecies, nElement
        do i=1, nSpecies
          write (UNIT=fU, FMT=FMTstr, IOSTAT=ios) &
            nameSpecies(i), SpeciesElements(:, i), massSpecies(i), &
            mobilitySpecies(i), vibFreqSpecies(i), E_desorp_species(i)
        end do
        close (UNIT=fU, IOSTAT=ios, STATUS='KEEP')
      end if

      if (IsWordChar(fRatesSave(1:1))) then
        write (*,'(/A)') 'Saving rates...'
        CALL openFileSequentialWrite &
          (fU, trim(path)//trim(fRatesSave), recordLen)
        write (FMTstr, FMT='("(I4, 2X, ", I2, "A", I2, &
          & ", ES15.4, ES11.2, ES11.2, 3F11.2, I4)")') &
          nParticipants, constLenNameSpecies
        do i=1, nReactions
          write (UNIT=fU, FMT=FMTstr, &
            IOSTAT=ios) i, strReactants(1:nReactants, i), &
            strProducts(1:nProducts, i), &
            rates(i), BranchingRatios(i), dblABC(:,i), &
            tunneling_mass(i), typeReac(i)
        end do
        close (UNIT=fU, IOSTAT=ios, STATUS='KEEP')
      end if

      end subroutine SaveMiscBefore


      subroutine getSurfRadical()
      use CMDT
      use CMTP
      implicit none
      integer i, j, i1, fU, nLineAll, nLineData, ios
      character, parameter :: commentChar = '!'
      character(len=12) nameSpecies_tmp
      logical getFileUnit
      if (.NOT. getFileUnit(fU)) then
        write (*,*) 'Cannot allocate an input file unit!'
        stop
      end if
      CALL openFileSequentialRead &
          (fU, trim(path)//trim(fSurfRadical), 9999)
      CALL GetNLineF (fU, nLineAll, nLineData, commentChar)
      nSurfRadical = nLineData
      rewind (UNIT=fU, IOSTAT=ios)
      IsSurfRadical = .FALSE.
      idxSurfRadical = 0
      i1 = 1
      do i=1, nLineAll
        read (UNIT=fU, FMT='(A12)', IOSTAT=ios) nameSpecies_tmp
        if (ios .NE. 0) then
          write (*,*) 'ios = ', ios; stop
        end if
        if ((nameSpecies_tmp(1:1) .EQ. commentChar) .OR. &
            (nameSpecies_tmp(1:1) .EQ. ' ')) cycle
        do j=1, nSpecies
          if (trim(nameSpecies(j)) .EQ. &
              trim(nameSpecies_tmp)) then
            IsSurfRadical(j) = .TRUE.
            idxSurfRadical(i1) = j
            i1 = i1 + 1
            exit
          end if
        end do
      end do
      close (UNIT=fU, IOSTAT=ios, STATUS='KEEP')
      end subroutine getSurfRadical



      subroutine getSpeInfo
      use CMDT
      use CMTP
      implicit none
      integer i
      logical flag
      nGrSpecies = 0
      nMantleSpecies = 0
      idxToMant = 0
      IsGas = .FALSE.
      IsMan = .FALSE.
      do i=1, nSpecies
        CALL getElements(nameSpecies(i), nameElements, nElement, &
          SpeciesElements(:, i))
        massSpecies(i) = sum(SpeciesElements(:, i) * ElementMassNumber)
        if (nameSpecies(i)(1:1) .EQ. 'g') then
          nGrSpecies = nGrSpecies + 1
          idxGrSpecies(nGrSpecies) = i
        else if (nameSpecies(i)(1:1) .EQ. 'm') then
          IsMan(i) = .TRUE.
          nMantleSpecies = nMantleSpecies + 1
          idxToMant(i) = nMantleSpecies
        else
          IsGas(i) = .TRUE.
        end if
      end do
      end subroutine getSpeInfo



      subroutine getSurfMantInfo
      use CMDT
      use CMTP
      implicit none
      integer i, i1, j
      logical flag
      i1 = 0
      do i=1, nMantleSpecies
        mantleSpe(i)%idxGas = 0
        mantleSpe(i)%idxSurf = 0
        mantleSpe(i)%idx_acc = 0
        mantleSpe(i)%idx_eva = 0
      end do
      do i=1, nSpecies
        if (IsMan(i)) then
          i1 = i1 + 1
          mantleSpe(i1)%idxMant = i
        end if
      end do
      IsSurfWithMan = .FALSE.
      do i=1, nMantleSpecies
        flag = .TRUE.
        do j= nSpecies, 1, -1
          if (nameSpecies(j) .EQ. &
              nameSpecies(mantleSpe(i)%idxMant)(2: &
                constLenNameSpecies)) then
            mantleSpe(i)%idxGas = j
            flag = .FALSE.
            exit
          end if
        end do
        if (flag) then
          write (*,*) nameSpecies(mantleSpe(i)%idxMant), &
            " does not have a gas phase counterpart!"
        end if
        flag = .TRUE.
        do j= nSpecies, 1, -1
          if ((nameSpecies(j)(1:1) .EQ. 'g') .AND. &
              (nameSpecies(j)(2:constLenNameSpecies) .EQ. &
               nameSpecies(mantleSpe(i)%idxMant)(2: &
                 constLenNameSpecies))) then
            mantleSpe(i)%idxSurf = j
            IsSurfWithMan(j) = .TRUE.
            flag = .FALSE.
            exit
          end if
        end do
        if (flag) then
          write (*,*) nameSpecies(mantleSpe(i)%idxMant), &
            " does not have a surface counterpart!"
        end if
        do j=nReactions, 1, -1
          select case(typeReac(j))
          case(61) ! Accretion
            if (reactions(1, j) .EQ. mantleSpe(i)%idxGas) then
              mantleSpe(i)%idx_acc = j
            end if
          case(62) ! Evaporation
            if (reactions(1, j) .EQ. mantleSpe(i)%idxSurf) then
              mantleSpe(i)%idx_eva = j
            end if
          case(65) ! Mantle species gets exposed.
            if (reactions(1, j) .EQ. mantleSpe(i)%idxMant) then
              mantleSpe(i)%idx_M2S = j
            end if
          case(66) ! Surface species gets covered.
            if (reactions(1, j) .EQ. mantleSpe(i)%idxSurf) then
              mantleSpe(i)%idx_S2M = j
            end if
          end select
        end do
      end do
      end subroutine getSurfMantInfo




      subroutine Neutralize(y, NEQ)
      use CMDT
      use CMTP
      implicit none
      double precision, dimension(NEQ) :: y
      integer NEQ
      double precision totalCharge
      double precision, parameter :: chargeThreshold = 1D-13
      totalCharge = dot_product(y(1:nSpecies), &
        SpeciesElements(idxCharge, :))
      if (abs(totalCharge) .GT. chargeThreshold/ratioGrainToH) then
        write (*,*) 'The initial condition is not neutral!'
        write (*,*) 'chargeThreshold/ratioGrainToH = ', &
          chargeThreshold/ratioGrainToH
        write (*,*) 'Will try to neutralize the system.'
        !end if
        !if (totalCharge .GE. 0D0) then
        y(idxElectron) = y(idxElectron) + totalCharge
        !else
        !  y(idxElectron) = y(idxElectron) - (-totalCharge)
        !end if
        if (y(idxElectron) .LT. 0D0) then
          write (*,*) 'I cannot neutralized the system!!!'
          write (*,*) 'totalCharge = ', totalCharge
          write (*,*) 'y(idxElectron) = ', y(idxElectron)
          write (*,*) 'Will run any way.'
        end if
      end if
        
      end subroutine Neutralize
