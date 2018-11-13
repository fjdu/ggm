! Runs fluently on 2010-08-09
!
! 2010-09-06:
!   Adapted for gas-grain.
!   Corrected a severe error:
!     !polyReac(i)%vecMon(2)%vecExp(1:nRealReactants(i)) = 1
!     ->  polyReac(i)%vecMon(2)%vecExp(1:nRealProducts(i)) = 1
!   This error does not show up previously, because previously all
!     the reactions have relatively small number of products.
!   Treat the gas species the same way as grain species.
!
! XXX Start final revision again on 2010-10-26
!
! 2011-01-17 Mon 17:44:28
!
! Its previous filename is moment_rate06_Drop.f90, copied from
!   chemical_modeling/moment_eq_20101026
! Modification for changing physical conditions.
!
! 2011-01-20 Thu 22:17:47
! This one is only used for analysis.
!
! 2011-04-13 Wed 11:23:02
! The only modification is to change the string format, namely, the space padding.
!
!gfortran common_MC_MO.o common_MC_MO_TP.o subroutines_share_triv.o mo_analyse.o subroutines_share_proc.o subroutines_mo.o opkd*.o -fbounds-check -lpgplot -lX11 -o mo_analyse
!gfortran common_MC_MO.f90 common_MC_MO_TP.f90 subroutines_share_triv.f90 mo_analyse.f90 subroutines_share_proc.f90 subroutines_mo.f90 opkd*.o -fbounds-check -lpgplot -lX11 -o mo_analyse
!cp /homes/fjdu/CommonFiles/H2O2_model_20110714_B/fortran/mo_analyse ./
!./mo_analyse run_0.1/run_0.1_21.0_6.0E+05_1.0E-07_15.0_1.0E-01/config_moment_.dat species_ana.dat
      program moment_auto
      use CMDT
      use CMTP
      implicit none

      external f, jac, getPhyPar

      integer fU, fUOutLog, fUyHistory_ascii, fUyHistory_bin, fUPhyPar
      character, parameter :: commentChar = '!'
      integer ios, statALLOC
      integer nLineAll, nLineData
      character strTMP*128, strTMP1*128
      character FMTstr*128, FMTstryHistory*128, FMTstrPhyPar*128
      character(LEN=constLenNameSpecies) nameSpecies_tmp
      logical IsWordChar

      integer i, j, k, h, i1, i2, i3, nNumMMTmp
      double precision, dimension(3) :: dblTmpVec
      double precision initialAbundance_tmp
      integer, dimension(8) :: intTmpvec
      logical flag, flag1, flag2, flagChange

      real ProgStartTime, ProgCurrentTime

      double precision, dimension(:), allocatable :: &
        initialElementalAb, finalElementalAb
      integer, dimension(:), allocatable :: &
        nElementsReac, nElementsProd

      integer nIterationBase, nTargetDecades, nerr, nChange
      double precision ratioTout, tStep
      double precision t, tout
      double precision, dimension(:), allocatable :: y, yBak, RWORK
      double precision, dimension(:), allocatable :: ydotTMP
      integer IOPT, ISTATE, ITASK, ITOL, LIW, LRW, NEQ, MF, NNZ
      integer, dimension(:), allocatable :: IWORK

      character(len=256), dimension(dimvecMon) ::  StrPolynomial
      character(len=16384) StrEqn

      integer, dimension(:), allocatable :: iWhenDropped

      logical getPhyPar, ExternalAction, getFileUnit, FileUnitOpened

      character fileSpeciesAna*128, fileElementsReservoir*128, &
        fileWhatsHappening*128
      integer getIdxSpecies, idxSpeciesAna
      integer nSpeciesAna
      character(len=constLenNameSpecies), dimension(:), allocatable :: &
        nameSpeciesAna
      double precision, dimension(:,:), allocatable :: elementsReservoir
      integer, dimension(:), allocatable :: idxSpeEleSorted
     !integer, parameter :: nSpeciesAna = 28
     !character(len=constLenNameSpecies), dimension(nSpeciesAna) :: &
     !  nameSpeciesAna = &
     !    (/'gH2O2   ', 'H2O2    ', 'gCH3OH  ', 'CH3OH   ', &
     !      'gH2O    ', 'H2O     ', 'gO2     ', 'O2      ', &
     !      'gO2H    ', 'O2H     ', 'gCO2    ', 'CO2     ', &
     !      'gO      ', 'O       ', 'gH      ', 'H       ', &
     !      'gO3     ', 'O3      ', 'gH2CO   ', 'H2CO    ', &
     !      'gOH     ', 'OH      ', 'HCN     ', 'HNC     ', &
     !      'C2H     ', 'C2H2    ', 'gNH3    ', 'gCH4    '/)

      integer nReactionContrib
      double precision, dimension(:,:), allocatable :: reactionContrib
      double precision, dimension(:,:), allocatable :: reactionAccum
      double precision, dimension(:), allocatable :: reactionContribAve
      double precision, dimension(:), allocatable :: reactionAccumAve
      double precision, dimension(:), allocatable :: allReacRates
      integer, dimension(:), allocatable :: idxYdotOneReactionSorted
      integer, dimension(:), allocatable :: idxAccumSorted
      integer, dimension(:), allocatable :: idxWhenJump
      double precision, dimension(:,:), allocatable :: ratesSpecies
      double precision dblNaN
      double precision, dimension(2), parameter :: &
        timerange_an=(/1D4, 5D6/)
      integer, dimension(2) :: idxrange_an
      integer imaxloc
      double precision maxAb

      integer nPGpoints, PGBEG
      real, dimension(:), allocatable :: xPGPLOT, yPGPLOT
      real xPGmin, xPGmax, yPGmax, yPGmin, xPGstep, yPGstep
      character(len=128) getFilePreName
      character(len=256) combineStrArr
      real, dimension(18) :: RGB_R, RGB_G, RGB_B
      INTERFACE
        SUBROUTINE setEnvPG (ColorIdx, LineWidth, LineStyle, &
          ifClip, CharSize, BGColor)
          INTEGER, INTENT (IN), OPTIONAL :: ColorIdx
          INTEGER, INTENT(IN), OPTIONAL :: LineWidth
          INTEGER, INTENT(IN), OPTIONAL :: LineStyle
          INTEGER, INTENT(IN), OPTIONAL :: ifClip
          REAL, INTENT(IN), OPTIONAL :: CharSize
          INTEGER, INTENT(IN), OPTIONAL :: BGColor
        END SUBROUTINE setEnvPG
      END INTERFACE
      RGB_R = (/1.0, 0.0, 0.0, 0.8, 0.0, 1.0, 0.7, 0.1, 0.3, 0.5, 0.1, 0.8, 1.0, 0.5, 0.7, 1.0, 0.7, 0.5/)
      RGB_G = (/0.0, 0.8, 0.0, 0.8, 0.8, 0.0, 0.2, 0.5, 0.3, 0.5, 0.5, 0.1, 0.7, 0.8, 0.5, 0.5, 0.8, 0.5/)
      RGB_B = (/0.0, 0.0, 1.0, 0.0, 0.8, 1.0, 0.2, 0.1, 0.7, 0.1, 0.5, 0.8, 0.5, 0.7, 1.0, 0.8, 0.5, 1.0/)

      fileElementsReservoir = 'Elements_reservoir.dat'
      fileWhatsHappening = 'whatsHappening.dat'

      CALL CPU_TIME(ProgStartTime)
      call date_and_time(DATE=strTMP, TIME=strTMP1)

      write (*,'(A, A60, A)') &
        char(27)//'[47m'//char(27)//'[34m'//&
        'Program start', '', char(27)//'[0m'
      write (*,'("Current time: ", A8, 2X, A10)') strTMP, strTMP1

      CALL GET_COMMAND_ARGUMENT(1, strTMP, i, ios)
      if ((ios .NE. 0) .OR. (i .EQ. 0)) strTMP = 'config.dat'
      CALL GET_COMMAND_ARGUMENT(2, fileSpeciesAna, i, ios)
      if ((ios .NE. 0) .OR. (i .EQ. 0)) &
        fileSpeciesAna = 'SpeciesToAna.dat'

      ! This initialization imports, and only imports all the configuration
      ! parameters.

      write (*, '(/A)') 'Initializing...'
      CALL initialize (strTMP)
      write (*,*) 'useZeroPointEnergy = ', useZeroPointEnergy

      if (.NOT. IsWordChar(fAnalyseResult(1:1))) then
        write (*,*) 'fAnalyseResult is not provided!'
        stop
      end if

      if (.NOT. getFileUnit(fU)) then
        write (*,*) 'Cannot allocate an output file unit!'
        stop
      end if
      CALL openFileSequentialRead(fU, fileSpeciesAna, 999)
      CALL GetNLineF (fU, nLineAll, nLineData, commentChar)
      rewind (UNIT=fU, IOSTAT=ios)
      nSpeciesAna = nLineData
      allocate(nameSpeciesAna(nSpeciesAna))
      i = 1
      do j=1, nLineAll
        read (UNIT=fU, FMT='(A128)', IOSTAT=ios) strTMP
        if ((strTMP(1:1) .NE. commentChar) .AND. &
            (strTMP(1:1) .NE. ' ')) then
          nameSpeciesAna(i) = strTMP(1:constLenNameSpecies)
          i = i+1
        end if
      end do
      close (fU)

      ! Find a file unit for message output.
      !if (IsWordChar(fOutLog(1:1))) then
      !  if (.NOT. getFileUnit(fUOutLog)) then
      !    write (*,*) 'Cannot allocate an output file unit!'
      !    stop
      !  end if
      !  CALL openFileSequentialWrite &
      !    (fUOutLog, trim(path)//trim(fOutLog), 9999999)
      !  CALL XSETUN(fUOutLog)
      !else
      !  do i=16, 2, -1
      !    if (ISATTY(unit=i)) then
      !      INQUIRE(UNIT=i, action=strTMP)
      !      if (trim(strTMP) .EQ. 'WRITE') then
      !        fUOutLog = i; exit
      !      end if
      !    end if
      !  end do
      !end if

! Import all the reactions

      write (*, '(/A)') 'Importing reactions...'

      if (.NOT. getFileUnit(fU)) then
        write (*,*) 'Cannot allocate an output file unit!'
        stop
      end if
      CALL openFileSequentialRead &
        (fU, trim(path)//trim(fReactions), 999)
      CALL GetNLineF (fU, nLineAll, nLineData, commentChar)

      nReactions = nLineData

      if (nSpecies_Est .LE. 0) & ! This estimation would be surely
        nSpecies_Est = nReactions * nParticipants ! enough.

      allocate &
        (strReactants(nReactants, nReactions), &
         strProducts(nProducts, nReactions), &
         nameSpecies(nSpecies_Est), &
         nRealReactants(nReactions), &
         nRealProducts(nReactions), &
         reactions(nParticipants, nReactions), &
         dblABC(3, nReactions), &
         tunneling_mass(nReactions), &
         typeReac(nReactions), &
         rates(nReactions), &
         BranchingRatios(nReactions), STAT=statALLOC)

      CALL ReadReactions (fU, nLineAll, commentChar)
      close (UNIT=fU, IOSTAT=ios, ERR=999, STATUS='KEEP')

      CALL MakeReactions

      allocate &
        (SpeciesElements(nElement, nSpecies), &
         massSpecies(nSpecies), &
         mobilitySpecies(nSpecies), &
         initialElementalAb(nElement), &
         finalElementalAb(nElement), &
         nElementsReac(nElement), &
         nElementsProd(nElement), &
         nReac_AsReactants(nSpecies), &
         nReac_AsProducts(nSpecies), &
         idxReac_AsReactants(nReactions, nSpecies), &
         idxReac_AsProducts(nReactions, nSpecies), &
         idxGrSpecies(nSpecies), &
         idxGrReactions(nReactions), &
         nChangeSpe(nSpecies), &
         idxChangeSpe(16, nSpecies), &
         IsStoSpecies(nSpecies), &
         IsGas(nSpecies), &
         IsMan(nSpecies), &
         IsSurfWithMan(nSpecies), &
         IsEndProd(nSpecies), &
         enthalpySpecies(nSpecies), &
         E_desorp_species(nSpecies), &
         vibFreqSpecies(nSpecies), &
         idxSurfRadical(nSpecies), &
         IsSurfRadical(nSpecies), &
         idxToMant(nSpecies), &
         idxReacMM(nReactions), &
         STAT=statALLOC)

      call getSpeInfo
      allocate(mantleSpe(nMantleSpecies), STAT=statALLOC)
      call getSurfMantInfo

      nGrReactions = 0
      nReac_AsReactants = 0
      nReac_AsProducts = 0
      do i=1, nReactions
        if (nRealReactants(i) .EQ. 2) then ! Never exceed two.
        ! Put the species in a reaction in order
          if (reactions(1, i) .GT. reactions(2, i)) then
            CALL SwapInt(reactions(1, i), reactions(2, i))
          end if
        end if
        if (nRealProducts(i) .GE. 2) then
          CALL SORT_Asc_idx_Int(nRealProducts(i), &
            reactions(nReactants+1:nReactants+nRealProducts(i), i), &
            intTmpvec(1:nRealProducts(i)))
          reactions(nReactants+1:nReactants+nRealProducts(i), i) = &
            reactions(nReactants+intTmpvec(1:nRealProducts(i)), i)
        end if

        do j=1, nRealReactants(i)
          flag = .TRUE.
          do k=1, j-1
            if (reactions(k, i) .EQ. reactions(j, i)) then
              flag = .FALSE.
              exit
            end if
          end do
          if (flag) then
            nReac_AsReactants(reactions(j, i)) = &
              nReac_AsReactants(reactions(j, i)) + 1
            idxReac_AsReactants(nReac_AsReactants(reactions(j, i)), &
              reactions(j, i)) = i
          end if
        end do

        do j=1+nReactants, nRealProducts(i)+nReactants
          flag = .TRUE.
          do k=1+nReactants, j-1
            if (reactions(k, i) .EQ. reactions(j, i)) then
              flag = .FALSE.
              exit
            end if
          end do
          if (flag) then
            nReac_AsProducts(reactions(j, i)) = &
              nReac_AsProducts(reactions(j, i)) + 1
            idxReac_AsProducts(nReac_AsProducts(reactions(j, i)), &
              reactions(j, i)) = i
          end if
        end do

        if (typeReac(i) .GE. 60) then
          nGrReactions = nGrReactions + 1
          idxGrReactions(nGrReactions) = i
        end if

        nElementsReac = 0
        nElementsProd = 0
        do j=1, nRealReactants(i)
          nElementsReac = nElementsReac + &
            SpeciesElements(:, reactions(j, i))
        end do
        do j=1, nRealProducts(i)
          nElementsProd = nElementsProd + &
            SpeciesElements(:, reactions(nReactants+j, i))
        end do
        if ((abs(nElementsReac(1) - nElementsProd(1)) + &
          sum(abs(nElementsReac(3:nElement) - &
          nElementsProd(3:nElement)))) .NE. 0) then
          write (*, '(2A, I6, A, 2X, 2A12, " -> ", 5A12)') &
            'Elements not conserved [discarded]: ', &
            char(27)//'[41m', i, char(27)//'[0m', &
            strReactants(1:nReactants, i), strProducts(1:nProducts, i)
          nRealReactants(i) = 0
          nRealProducts(i) = 0
        end if
        do j=1, i-1
          if ((typeReac(j) .EQ. typeReac(i)) .AND. &
              (sum(abs(reactions(:, j)-reactions(:, i))) .EQ. 0)) then
            write (*,'(2A, 2I6, A)') &
              'Duplicate reaction pair: ', &
              char(27)//'[45m', i, j, char(27)//'[0m'
            write (fUOutLog,'(2A, 2I6, A)') &
              'Duplicate reaction pair: ', &
              char(27)//'[45m', i, j, char(27)//'[0m'
          end if
        end do
      end do
      deallocate (nElementsReac, nElementsProd, STAT=statALLOC)

      allocate(surf_mant_inc(nGrReactions), surf_mant_dec(nGrReactions))
      call getSurfGainLossReac

      call getSurfRadical

      if (useThreeBody) then
        nThreeBodySurfReacs = 0
        do i=1, nGrReactions
          if (typeReac(idxGrReactions(i)) .EQ. 67) then
            nThreeBodySurfReacs = nThreeBodySurfReacs + 1
          end if
        end do
        allocate(ThreeBodySurfReacs(nThreeBodySurfReacs), &
          rateThreeBodySurf(nThreeBodySurfReacs))
        i1 = 1
        do i=1, nGrReactions
          if (typeReac(idxGrReactions(i)) .EQ. 67) then
            ThreeBodySurfReacs(i1) = idxGrReactions(i)
            i1 = i1 + 1
          end if
        end do
      end if

      write (*, FMT='(7(/, A32, I16))') &
        'Number of species: ', nSpecies, &
        'Number of grain species: ', nGrSpecies, &
        'Number of reactions: ', nReactions, &
        'Number of surface reactions: ', nGrReactions, &
        'Max number of reactants: ', maxval(nRealReactants), &
        'Max number of products: ', maxval(nRealProducts), &
        'Number of surface radicals: ', nSurfRadical

      write (*, '(/A)') 'Calculating rates...'

      ! These three do not depend on the physical conditions.
      CALL getVibFreq
      CALL getDesorEnergy ! From the reaction file.
      CALL ImportSpeciesEnthalpy
      CALL MakeBranchInfo

      ! These rates does depend on the physical conditions.
      CALL CalcRates

      rates = rates * SecondsPerYear

      write (*, '(3(/A32, ES16.6))') &
        'Grain to H number ratio: ', ratioGrainToH, &
        'Sites per grain: ', SitesPerGrain, &
        'Volumn containing one grain: ', VolumnPerGrain

      where (rates .LT. rateThreshold)
        rates = 0D0
        nRealReactants = 0
        nRealProducts = 0
      end where

! Now I need to prepare all the possibly needed moments.

      write (*, '(A)') 'Making moments...'
      ! 128 is for a small network; 4 is for bigger one.
      !nNumMMEst = nSpecies*16 + nReactions * 128
      !nNumMMEst = nSpecies*16 + nReactions * 4
      allocate (&
        polyReac(nReactions), &
        MM(dimvecDeriv, nNumMMEst), &
        nnzMM(nNumMMEst), &
        dtMM(3, nReactions*2, nNumMMEst), & ! The number of terms might
        ndtMM(nNumMMEst), &  ! be larger than the number of reactions.
        indMM(nNumMMEst, dimvecDeriv), &
        dauMM(nIterMoment+2, nNumMMEst), &
        ndauMM(nNumMMEst))

      do i=1, nReactions ! Make the generating function of each reaction
        call initpoly(polyReac(i)) ! Initialization is important.
        polyReac(i)%nnonzero = 2
        polyReac(i)%vecMon(1)%coeff = -1
        polyReac(i)%vecMon(1)%vecVar(1:nRealReactants(i)) = &
          reactions(1:nRealReactants(i), i)
        polyReac(i)%vecMon(1)%vecExp(1:nRealReactants(i)) = 1
        polyReac(i)%vecMon(1)%nvecDeriv = nRealReactants(i)
        polyReac(i)%vecMon(1)%vecDeriv(1:nRealReactants(i)) = &
          reactions(1:nRealReactants(i), i)
        polyReac(i)%vecMon(2)%coeff = 1
        polyReac(i)%vecMon(2)%vecVar(1:nRealProducts(i)) = &
          reactions(1+nReactants:nRealProducts(i)+nReactants, i)
        polyReac(i)%vecMon(2)%vecExp(1:nRealProducts(i)) = 1
        polyReac(i)%vecMon(2)%nvecDeriv = nRealReactants(i)
        polyReac(i)%vecMon(2)%vecDeriv(1:nRealReactants(i)) = &
          reactions(1:nRealReactants(i), i)
      end do

      nnzMM = 0
      ndtMM = 0
      ndauMM = 0
      nAllMM = 0
      indMM = 0

      nAllMM(1) = nSpecies
      do i=1, nSpecies ! All the first order moments.
        MM(1, i) = i
        nnzMM(i) = 1
        indMM(i, 1) = i
      end do
      nNumMM = nSpecies
      nNewMM = (/1, nSpecies/)

      flag2 = .FALSE.
      i3 = 0
      do i=1, 1!nIterMoment
        if (flag2) exit
        flag1 = .TRUE.
        flag2 = .TRUE.
        i3 = i3 + 1
        write (*,*) 'Loop ', i
        do j=nNewMM(1), nNewMM(2)
          if (nnzMM(j) .GT. nOrderLim) cycle
          if (nnzMM(j) .GT. 1) then
            if ((getFstGasSpe(j) .NE. 0) .OR. &
                (getFstNoneRadical(j) .NE. 0)) cycle
          end if
          do k=1, nReactions
            if (rates(i) .LT. rateThreshold) cycle
            call derivpolymulti &
              (MM(1:nnzMM(j), j), nnzMM(j), polyReac(k), pderivTmp)
            pmomTmp = polyMon2Moment(pderivTmp)
            do h=1, pmomTmp%nnonzero
              flag = .TRUE.
              orderTmp = pmomTmp%vecMoment(h)%nvecDeriv
              do i2=1, nNumMM
                if (nnzMM(i2) .NE. orderTmp) cycle
                if (IsEqualVec(pmomTmp%vecMoment(h)%&
                  vecDeriv(1:orderTmp), &
                  MM(1:orderTmp, i2), orderTmp)) then
                  flag = .FALSE.
                  ndtMM(j) = ndtMM(j) + 1
                  dtMM(:, ndtMM(j), j) = &
                    (/i2, k, pmomTmp%vecMoment(h)%coeff/)
                  exit
                end if
              end do
              if (flag) then
                flag2 = .FALSE.
                nNumMM = nNumMM + 1
                if (nNumMM .GT. nNumMMEst) then
                  write (*, *) 'nNumMM .GT. nNumMMEst'
                  write (*, *) nNumMM, ' .GT. ', nNumMMEst
                  write (*, *) '  i = ', i
                  stop
                end if
                MM(1:orderTmp, nNumMM) = &
                  pmomTmp%vecMoment(h)%vecDeriv(1:orderTmp)
                nnzMM(nNumMM) = orderTmp
                ndtMM(j) = ndtMM(j) + 1
                dtMM(:, ndtMM(j), j) = &
                  (/nNumMM, k, pmomTmp%vecMoment(h)%coeff/)
                nAllMM(orderTmp) = nAllMM(orderTmp) + 1
                indMM(nAllMM(orderTmp), orderTmp) = nNumMM
                if (flag1) then
                  nNewMM(1) = nNumMM
                  flag1 = .FALSE.
                else
                  nNewMM(2)  = nNumMM
                end if
                if (orderTmp .GT. 1) then
                  nNumMMTmp = nNumMM
                  call makedau &
                    (MM(1:orderTmp, nNumMM), orderTmp, nNumMMTmp, flag1)
                end if
              end if
            end do
          end do
        end do
      end do
      write (*, *) &
        'Number of iteration for making the moments: ', i3
      if (i3 .EQ. (nIterMoment)) write (*, *) &
        '*** Maybe the iteration is not enough! ***'
      CALL CPU_TIME(ProgCurrentTime)
      write (*, '(/A, F10.3/)') 'Seconds elapsed: ', &
        ProgCurrentTime - ProgStartTime

      do i=1, nReactions
        idxReacMM(i) = getReactionMM(i)
      end do

      IsEndProd = .TRUE.
      do i=1, nGrReactions
        i1 = idxGrReactions(i)
        if ((typeReac(i1) .LT. 63) .OR. &
            (typeReac(i1) .GT. 64)) cycle
        do j=1, nRealReactants(i1)
          IsEndProd(reactions(j, i1)) = .FALSE.
        end do
      end do

      deallocate(polyReac, indMM)

      allocate(IsY(nNumMM), countMMNotUsed(nNumMM))
      IsY(1:nSpecies) = .TRUE.
      do i=nSpecies+1, nNumMM
        if (nnzMM(i) .GT. nOrderLim) then
          IsY(i) = .FALSE.
        else
          if ((getFstGasSpe(i) .GT. 0) .OR. &
              (getFstNoneRadical(i) .GT. 0)) then
            IsY(i) = .FALSE.
          else
            IsY(i) = .TRUE.
          end if
        end if
      end do

      NEQ = count(IsY)
      write (*,*) 'NEQ = ', NEQ
      allocate(indY(NEQ), invIndY(nNumMM), y(NEQ), yBak(nNumMM), &
        ydotTMP(nNumMM))
      i1 = 0
      do i=1, nNumMM
        if (IsY(i)) then
          i1 = i1 + 1
          indY(i1) = i
          invIndY(i) = i1
        else
          invIndY(i) = 0
        end if
      end do

      write (*,'(/A32, I8)') 'Number of moments: ', nNumMM
      do i=1, dimvecDeriv
        if (nAllMM(i) .EQ. 0) exit
        write (*, '("Number of ", I1, "-moments: ", I8)') &
          i, nAllMM(i)
      end do

      CALL CPU_TIME(ProgCurrentTime)
      write (*, '(/A, F10.3/)') 'Seconds elapsed: ', &
        ProgCurrentTime - ProgStartTime

      CALL openFileSequentialRead &
        (fUyHistory_ascii, trim(path)//trim(fSaveALLHistory_ASCII), 999)
      CALL GetNLineF (fUyHistory_ascii, &
        nLineAll, nLineData, commentChar)
      close (fUyHistory_ascii)

      nHistoryLen = nLineAll-1
      allocate(yHistory(nSpecies+nAuxVar, nHistoryLen), &
               ydotHistory(nSpecies, nHistoryLen), &
               PhyParHistory(nPhyPar, nHistoryLen), &
               touts(nHistoryLen), &
               elementsReservoir(nSpecies, nElement), &
               idxSpeEleSorted(nSpecies), &
               ratesSpecies(nHistoryLen, nSpecies), STAT=statALLOC)

      if (IsWordChar(fSaveALLHistory(1:1))) then
        if (.NOT. getFileUnit(fUyHistory_bin)) then
          write (*,*) 'Cannot allocate an output file unit!'
          stop
        end if
        open (UNIT=fUyHistory_bin, &
             FILE=trim(path)//trim(fSaveALLHistory), &
             IOSTAT=ios, RECL=kind(1D0)*(nSpecies+1+nAuxVar), &
             ACCESS='DIRECT', FORM='UNFORMATTED', ACTION='READ')
        if (ios .LT. 0) then
          write (*, *) 'OPEN FILE ERROR: IOSTAT=', ios
          stop
        end if
        do i=1, nHistoryLen
          read (fUyHistory_bin, rec=i) touts(i), yHistory(:, i)
        end do
        close (UNIT=fUyHistory_bin, IOSTAT=ios, ERR=999, STATUS='KEEP')
      end if

      IsStoSpecies = .FALSE.
      do i=1, nGrSpecies ! Only Surface species can be stochastic.
        i1 = idxGrSpecies(i)
        if ((.NOT. IsGas(i1)) .AND. &
            (yHistory(i1,1) .LT. sto_threshold)) then
          IsStoSpecies(i1) = .TRUE.
        end if
      end do

      do i=1, nElement
        initialElementalAb(i) = dot_product(SpeciesElements(i, :), &
          yHistory(1:nSpecies, 1))
      end do

      write (*,*) 'Elemental reservoir...'
      if (.NOT. getFileUnit(fU)) then
        write (*,*) 'Cannot allocate an output file unit!'
        stop
      end if
      call openFileSequentialWrite(fU, &
        trim(path)//fileElementsReservoir, 999)
      do i=1, nHistoryLen, 10
        write (*,*) i, ' ele.'
        elementsReservoir = 0D0
        do j=1, nSpecies
          elementsReservoir(j, :) = &
            yHistory(j, i) * SpeciesElements(:, j)
        end do
        write (fU, '("At ", ES12.4)') touts(i)
        do j=1, nElement
          CALL SORT_Asc_idx(nSpecies, &
            -abs(elementsReservoir(:, j)), idxSpeEleSorted)
          write (fU, '(2X, "Element: ", A7)') nameElements(j)
          do k=1, min(10, nSpecies)
            i1 = idxSpeEleSorted(k)
            if (elementsReservoir(i1, j) .LE. 1D-30) exit
            write (fU, '(4X, I2, 2X, A12, 2X, 2(ES12.4, 2X), &
              & F6.3, "%")') &
              k, nameSpecies(i1), yHistory(i1, i), &
              elementsReservoir(i1, j), &
              elementsReservoir(i1, j)/initialElementalAb(j)*1D2
          end do
        end do
      end do
      close(fU)

      if (.NOT. getFileUnit(fU)) then
        write (*,*) 'Cannot allocate an output file unit!'
        stop
      end if
      allocate (idxYdotOneReactionSorted(nReactions), &
        allReacRates(nReactions))
      call openFileSequentialWrite(fU, &
        trim(path)//fileWhatsHappening, 999)
      do i=1, nHistoryLen, 10
        write (*,*) i, ' what.'
        allReacRates = 0D0
        do k=1, nReactions
          if (typeReac(k) .LT. 60) then
            if (nRealReactants(k) .EQ. 1) then
              allReacRates(k) = rates(k) * yHistory(reactions(1, k), i)
            else
              allReacRates(k) = rates(k) * &
                yHistory(reactions(1, k), i) * &
                yHistory(reactions(2, k), i)
            end if
          end if
        end do
        CALL SORT_Asc_idx(nReactions, &
          -abs(allReacRates), idxYdotOneReactionSorted)
        write (fU, '("[", I0, "] ", "Time: ", ES12.4)') i, touts(i)
        write (fU, '(X, "--Gas Phase Reactions--")')
        do k=1, min(nReactions, 30)
          i1 = idxYdotOneReactionSorted(k)
          if (allReacRates(i1) .LE. allReacRates(1)*1E-10) exit
          write (fU, &
            '(2X, I3, ": " 7A12, 2X, 3ES13.4, 2X, ES8.2, 2F8.2)') &
            k, strReactants(1:2, i1), ' -> ', strProducts(1:4, i1), &
            allReacRates(i1), yHistory(reactions(1, i1), i), &
            rates(i1), dblABC(:, i1)
        end do
        allReacRates = 0D0
        call updateMantSurfRates(NEQ, yHistory(1:nSpecies,i))
        do j=1, nGrReactions
          k = idxGrReactions(j)
          if (typeReac(k) .LT. 63) cycle
          if ((typeReac(k) .EQ. 67) .AND. (useThreeBody)) then
            call getThreeBodySurfRate(nSpecies, yHistory(1:nSpecies, i))
            do i2=1, nThreeBodySurfReacs
              i1 = ThreeBodySurfReacs(i2)
              if (k .EQ. i1) then
                allReacRates(k) = rateThreeBodySurf(i2)
                exit
              end if
            end do
          else
            if (nRealReactants(k) .EQ. 1) then
              allReacRates(k) = rates(k) * yHistory(reactions(1, k), i)
            else
              allReacRates(k) = rates(k) * &
                yHistory(reactions(1, k), i) * &
                yHistory(reactions(2, k), i)
            end if
          end if
        end do
        CALL SORT_Asc_idx(nReactions, &
          -abs(allReacRates), idxYdotOneReactionSorted)
        write (fU, '(X, "--Surface Reactions--")')
        do k=1, min(nReactions, 30)
          i1 = idxYdotOneReactionSorted(k)
          if (allReacRates(i1) .LE. allReacRates(1)*1E-10) exit
          write (fU, &
            '(2X, I3, ": " 7A12, 2X, 3ES13.4, 2X, ES8.2, 2F8.2)') &
            k, strReactants(1:2, i1), ' -> ', strProducts(1:4, i1), &
            allReacRates(i1), yHistory(reactions(1, i1), i), &
            rates(i1), dblABC(:, i1)
        end do
        allReacRates = 0D0
        do j=1, nGrReactions
          k = idxGrReactions(j)
          if (typeReac(k) .EQ. 61) then
            allReacRates(k) = rates(k) * yHistory(reactions(1, k), i)
          end if
        end do
        CALL SORT_Asc_idx(nReactions, &
          -abs(allReacRates), idxYdotOneReactionSorted)
        write (fU, '(X, "--Accretion--")')
        do k=1, min(nReactions, 20)
          i1 = idxYdotOneReactionSorted(k)
          if (allReacRates(i1) .LE. allReacRates(1)*1E-10) exit
          write (fU, &
            '(2X, I3, ": " 7A12, 2X, 3ES13.4, 2X, ES8.2, 2F8.2)') &
            k, strReactants(1:2, i1), ' -> ', strProducts(1:4, i1), &
            allReacRates(i1), yHistory(reactions(1, i1), i), &
            rates(i1), dblABC(:, i1)
        end do
        allReacRates = 0D0
        do j=1, nGrReactions
          k = idxGrReactions(j)
          if (typeReac(k) .EQ. 62) then
            allReacRates(k) = rates(k) * yHistory(reactions(1, k), i)
          end if
        end do
        CALL SORT_Asc_idx(nReactions, &
          -abs(allReacRates), idxYdotOneReactionSorted)
        write (fU, '(X, "--Evaporation--")')
        do k=1, min(nReactions, 20)
          i1 = idxYdotOneReactionSorted(k)
          if (allReacRates(i1) .LE. allReacRates(1)*1E-10) exit
          write (fU, &
            '(2X, I3, ": " 7A12, 2X, 3ES13.4, 2X, ES8.2, 2F8.2)') &
            k, strReactants(1:2, i1), ' -> ', strProducts(1:4, i1), &
            allReacRates(i1), yHistory(reactions(1, i1), i), &
            rates(i1), dblABC(:, i1)
        end do
        write (fU, '(X, "--tot_man_acc_rate(reac) = ", 2ES15.4)') &
          tot_man_acc_rate, tot_man_acc_rate_reac
        write (fU, '(X, "--tot_man_eva_rate(reac) = ", 2ES15.4)') &
          tot_man_eva_rate, tot_man_eva_rate_reac
      end do
      close(fU)

      if (.NOT. getFileUnit(fU)) then
        write (*,*) 'Cannot allocate an output file unit!'
        stop
      end if

      do i=1, nSpecies
        ratesSpecies(1,i) = dblNaN()
        do j=2, nHistoryLen
          ratesSpecies(j,i) = (yHistory(i,j)-yHistory(i,j-1)) &
            / (touts(j)-touts(j-1))
        end do
      end do
      CALL openFileSequentialWrite &
        (fU, trim(path)//trim(getFilePreName(fAnalyseResult))//&
          "_ALLSpecies.dat", 999)
      do i=1, nSpecies
        write(fU, '(A12)') nameSpecies(i)
        do j=2, nHistoryLen-1, 10
          if (((yHistory(i,j) .GT. yHistory(i,j-1)) .AND. &
               (yHistory(i,j) .GT. yHistory(i,j+1))) .OR. &
              ((yHistory(i,j) .LT. yHistory(i,j-1)) .AND. &
               (yHistory(i,j) .LT. yHistory(i,j+1)))) then
            write(fU, '(2X, "Critical: ", I5, ES12.4)') j, touts(j)
          end if
          if (abs(log(yHistory(i,j+1) / yHistory(i,j))) .GT. &
              abs(log(yHistory(i,j) / yHistory(i,j-1))) + 2D0) then
            write(fU, '(2X, "Jump:     ", I5, ES12.4)') j, touts(j)
          end if
        end do
      end do
      close (fU)

      CALL openFileSequentialWrite &
        (fU, trim(path)//trim(getFilePreName(fAnalyseResult))//&
          "_ALLTime.dat", 999)
      do i=2, nHistoryLen-1, 10
        flag = .TRUE.
        do i1=1, nSpeciesAna
          j = getIdxSpecies(nameSpeciesAna(i1))
          if ((j .LT. 1) .OR. (j .GT. nSpecies)) then
            if (i .LT. 3) write(*,*) nameSpeciesAna(i1), &
              ' is probably not in the network!'
            cycle
          end if
          if (((yHistory(j,i) .GT. yHistory(j,i-1)) .AND. &
               (yHistory(j,i) .GT. yHistory(j,i+1))) .OR. &
              ((yHistory(j,i) .LT. yHistory(j,i-1)) .AND. &
               (yHistory(j,i) .LT. yHistory(j,i+1)))) then
            if (flag) then
              write(fU, '("At ", I4, ES12.4)') i, touts(i)
              flag = .FALSE.
            end if
            write(fU, '(2X, "Critical: ", A12)') nameSpecies(j)
          end if
          if (abs(log(yHistory(j,i+1) / yHistory(j,i))) .GT. &
              abs(log(yHistory(j,i) / yHistory(j,i-1))) + 2D0) then
            if (flag) then
              write(fU, '("At ", I4, ES12.4)') i, touts(i)
              flag = .FALSE.
            end if
            write(fU, '(2X, "Jump:     ", A12)') nameSpecies(j)
          end if
        end do
      end do
      close (fU)

      allocate (idxAccumSorted(nReactions))
      do i=1, nSpeciesAna
        idxSpeciesAna = getIdxSpecies(nameSpeciesAna(i))
        if ((idxSpeciesAna .LT. 1) .OR. &
            (idxSpeciesAna .GT. nSpecies)) then
         !write(*,*) nameSpeciesAna(i), &
         !  ' is probably not in the network!'
          cycle
        end if
        nReactionContrib = ndtMM(idxSpeciesAna)
        if (allocated(reactionContrib)) deallocate(reactionContrib)
        if (allocated(reactionAccum)) deallocate(reactionAccum)
        if (allocated(reactionContribAve)) &
          deallocate(reactionContribAve)
        if (allocated(reactionAccumAve)) &
          deallocate(reactionAccumAve)
        allocate(reactionContrib(nReactionContrib, nHistoryLen))
        allocate(reactionAccum(nReactionContrib, nHistoryLen))
        reactionAccum = 0D0
        allocate(reactionContribAve(nReactionContrib))
        allocate(reactionAccumAve(nReactionContrib))
        CALL openFileSequentialWrite &
          (fU, trim(path)//trim(getFilePreName(fAnalyseResult))//&
            "_"//trim(nameSpeciesAna(i))//".dat", 999)
        write (fU, '(A, 4X, "#", I4, 4X, "nReactionContrib: ", I6)') &
          nameSpeciesAna(i), idxSpeciesAna, nReactionContrib
        intTmpvec(1:1) = maxloc(yHistory(idxSpeciesAna, :))
        maxAb = yHistory(idxSpeciesAna, intTmpvec(1))
        write (fU, '(2X, A, ES12.4, A, ES10.2, 2(A,ES12.4)/)') &
          '!Max: ', maxAb, " At ", touts(intTmpvec(1)), &
          '  Init Ab: ', yHistory(idxSpeciesAna, 1), &
          '  Final Ab: ', yHistory(idxSpeciesAna, nHistoryLen)
        write (*,*) 'Looking for whats happenning...'
        do j=2, nHistoryLen
          !if (getPhyPar(touts(j))) then
          !  CALL CalcMobilities
          !  CALL CalcRates
          !  rates = rates * SecondsPerYear
          !  do i1=1, nNumMM
          !    if (ndtMM(i1) .GT. 0) then
          !      do i2=1, ndtMM(i1)
          !        dtMMMo(i1)%coeff(i2) = rates(dtMM(2, i2, i1)) * &
          !          DBLE(dtMM(3, i2, i1))
          !      end do
          !    end if
          !  end do
          !end if
          call f (NEQ, touts(j), yHistory(1:nSpecies, j), &
            ydotTMP(1:NEQ))
          if (IsMan(idxSpeciesAna) .OR. IsSurfWithMan(idxSpeciesAna)) then
            call updateMantSurfRates(NEQ, yHistory(1:nSpecies,j))
          end if
          do k=1, nReactionContrib
            if ((typeReac(dtMM(2, k, idxSpeciesAna)) .EQ. 67) .AND. &
                (useThreeBody)) then
              call getThreeBodySurfRate(nSpecies, &
                yHistory(1:nSpecies,j))
              do i2=1, nThreeBodySurfReacs
                i1 = ThreeBodySurfReacs(i2)
                if (dtMM(2, k, idxSpeciesAna) .EQ. i1) then
                  if (dtMM(3, k, idxSpeciesAna) .LT. 0) then
                    reactionContrib(k, j) = -rateThreeBodySurf(i2)
                  else
                    reactionContrib(k, j) = rateThreeBodySurf(i2)
                  end if
                  exit
                end if
              end do
            else
              reactionContrib(k, j) = &
                rates(dtMM(2, k, idxSpeciesAna)) * &
                DBLE(dtMM(3, k, idxSpeciesAna)) * &
                getMMProd(dtMM(1, k, idxSpeciesAna), &
                  yHistory(1:nSpecies,j), nSpecies)
            end if
            reactionAccum(k, j) = reactionAccum(k, j-1) + &
                reactionContrib(k, j) * (touts(j) - touts(j-1))
          end do
          CALL SORT_Asc_idx(nReactionContrib, &
            -abs(reactionContrib(:, j)), idxYdotOneReactionSorted)
          write (fU, '("[", I0, "] ", "Time: ", ES12.4, 2X, &
            "Abundance: ", ES12.4, "(", F4.2, ")", 2X, &
            "Total rate: ", ES12.4, 2X, &
            "vs ", ES12.4, 2X, "y(j)-y(j-1): ", ES12.4)') &
            j, touts(j), yHistory(idxSpeciesAna, j), &
            yHistory(idxSpeciesAna, j)/maxAb, &
            sum(reactionContrib(:, j)), &
            ydotTMP(idxSpeciesAna), &
            yHistory(idxSpeciesAna, j)-yHistory(idxSpeciesAna, j-1)
          do k=1, min(nReactionContrib, 50)
            i1 = idxYdotOneReactionSorted(k)
            i2 = dtMM(2, i1, idxSpeciesAna)
            write (fU, &
              '(2X, I3, ES11.2, 2X, ES11.2, ": " &
              & 7A8, 2X, ES15.4, 4X, ES8.2, 2F8.2)') &
              k, reactionContrib(i1, j), reactionAccum(i1, j), &
              strReactants(1:2, i2), ' -> ', strProducts(1:4, i2), &
              rates(dtMM(2, i1, idxSpeciesAna)) * &
              DBLE(dtMM(3, i1, idxSpeciesAna)), dblABC(:, i2)
          end do
          if (j .GT. 0) then
           if ((yHistory(idxSpeciesAna, j-1) .LT. maxAb * 0.1) .AND. &
               (yHistory(idxSpeciesAna, j) .GE. maxAb * 0.1)) then
             write (fU, '(A)') '  !Reached 10%!'
           end if
           if ((yHistory(idxSpeciesAna, j-1) .LT. maxAb * 0.5) .AND. &
               (yHistory(idxSpeciesAna, j) .GE. maxAb * 0.5)) then
             write (fU, '(A)') '  !Reached 50%!'
           end if
           if ((yHistory(idxSpeciesAna, j-1) .LT. maxAb * 0.9) .AND. &
               (yHistory(idxSpeciesAna, j) .GE. maxAb * 0.9)) then
             write (fU, '(A)') '  !Reached 90%!'
           end if
           if ((yHistory(idxSpeciesAna, j-1) .GT. maxAb * 0.9) .AND. &
               (yHistory(idxSpeciesAna, j) .LE. maxAb * 0.9)) then
             write (fU, '(A)') '  !Dropped below 90%!'
           end if
           if ((yHistory(idxSpeciesAna, j-1) .GT. maxAb * 0.5) .AND. &
               (yHistory(idxSpeciesAna, j) .LE. maxAb * 0.5)) then
             write (fU, '(A)') '  !Dropped below 50%!'
           end if
           if ((yHistory(idxSpeciesAna, j-1) .GT. maxAb * 0.1) .AND. &
               (yHistory(idxSpeciesAna, j) .LE. maxAb * 0.1)) then
             write (fU, '(A)') '  !Dropped below 10%!'
           end if
          end if
        end do
        close (UNIT=fU, IOSTAT=ios, ERR=999, STATUS='KEEP')

        do i1=1, nHistoryLen
          if (touts(i1) .GE. timerange_an(1)) then
            idxrange_an(1) = i1
            exit
          end if
        end do
        do i1=nHistoryLen, 1, -1
          if (touts(i1) .LE. timerange_an(2)) then
            idxrange_an(2) = i1
            exit
          end if
        end do
        reactionContribAve = 0D0
        reactionAccumAve = 0D0
        do i1=idxrange_an(1), idxrange_an(2)
          reactionContribAve = reactionContribAve + &
            reactionContrib(:, i1) * (touts(i1+1) - touts(i1))
          reactionAccumAve = reactionAccumAve + &
            reactionAccum(:, i1)
        end do
        CALL SORT_Asc_idx(nReactionContrib, &
            -abs(reactionContribAve), idxYdotOneReactionSorted)
        CALL SORT_Asc_idx(nReactionContrib, &
            -abs(reactionAccumAve), idxAccumSorted)

        write (*,'(/A, I6, 2X, A)') 'Plotting...', i, nameSpeciesAna(i)
        j = PGBEG(0, trim(path)//&
          trim(getFilePreName(fAnalyseResult))//"_"// &
          trim(nameSpeciesAna(i))//'.ps/CPS', -1, 4)
        IF (j .NE. 1) STOP
        CALL PGPAP(15.0, 0.7)
        nPGpoints = nHistoryLen
        allocate(xPGPLOT(nPGpoints), yPGPLOT(nPGpoints), STAT=statALLOC)
        xPGPLOT = log10(touts + 1D-10)
        xPGmin = -2.0 !xPGPLOT(2)
        xPGmax = xPGPLOT(nPGpoints)
        do i1=1, nPGpoints
          if (xPGPLOT(i1+1) .GT. xPGmin) then
            intTmpvec(1) = i1
            exit
          end if
        end do
        do i1=nPGpoints, 1, -1
          if (xPGPLOT(i1-1) .LT. xPGmax) then
            intTmpvec(2) = i1
            exit
          end if
        end do
        yPGmax = log10(maxval(abs(reactionContrib(:, &
                            intTmpvec(1):intTmpvec(2))))) + 0.2
        yPGmin = max(log10(minval(abs(reactionContrib(:, &
                            intTmpvec(1):intTmpvec(2))))), yPGmax-10.0)
        xPGstep = (xPGmax - xPGmin) / 10.0
        yPGstep = (yPGmax - yPGmin) / 25.0
        CALL setEnvPG (ColorIdx=1, LineWidth=3, LineStyle=1)
        CALL PGPAGE ! Creat a plot page
        CALL PGSVP (0.1, 0.65, 0.00, 0.80) ! page range
        CALL PGSWIN(xPGmin, xPGmax, yPGmin, yPGmax)
        CALL PGBOX ('BCSTL', 0.0, 0, 'BCNSTVL', 0.0, 0)
        CALL PGMTXT('L', 4.0, 0.5, 0.5, &
          'Absolute Reaction Rates (yr\u-1\d)')
       !CALL PGMTXT('B', 2.5, 0.5, 0.5, 'Time (year)')
        do k=1, min(nReactionContrib, 18)
          j = idxYdotOneReactionSorted(k)
          if (dtMM(3, j, idxSpeciesAna) .GT. 0) then
            yPGPLOT = log10(reactionContrib(j, :))
          else
            yPGPLOT = log10(-reactionContrib(j, :))
          end if
          !if (maxval(abs(yPGPLOT)) .LT. &
          !  min(yPGmin+1.0, yPGmax-1.0)) then
          !  cycle
          !end if
          i1 = mod(k-1,18) + 1
          CALL PGSCR (16, RGB_R(i1), RGB_G(i1), RGB_B(i1))
          !!CALL PGSHLS (16, mod(k*17, 360)*1.0, &
          !!    mod(k,3)/10.0+0.4, mod(j,2)/4.0+0.7)
          !CALL PGMOVE(xPGmax+xPGstep*2.1, yPGmax - k*yPGstep)
          !CALL PGDRAW(xPGmax+xPGstep*3.2, yPGmax - k*yPGstep)
          i1 = dtMM(2, j, idxSpeciesAna)
          strTMP = trim(combineStrArr(strReactants(1:2, i1), &
                2, "\(2284)"))//" \(2261) "//&
              trim(combineStrArr(strProducts(1:4, i1), 4, "\(2284)"))
          if (dtMM(3, j, idxSpeciesAna) .GT. 0) then
            CALL setEnvPG (ColorIdx=16, LineWidth=5, LineStyle=1, &
                    ifClip=1)
            strTMP = "F: "//strTMP
          else
            CALL setEnvPG (ColorIdx=16, LineWidth=5, LineStyle=2, &
                    ifClip=1)
            strTMP = "D: "//strTMP
          end if
          !CALL PGMOVE(xPGmax+xPGstep*2.1, yPGmax - k*yPGstep)
          !CALL PGDRAW(xPGmax+xPGstep*3.2, yPGmax - k*yPGstep)
          CALL PGLINE(nPGpoints, xPGPLOT, yPGPLOT)
          CALL setEnvPG (ColorIdx=16, LineWidth=2, LineStyle=1, &
                         CharSize=1.2, BGColor=0)
          CALL PGPTXT (xPGmax+xPGstep*0.1, &
              yPGmax - (k+0.25)*yPGstep, &
              0.0, 0.0, trim(strTMP))
        end do

        CALL setEnvPG (ColorIdx=1, LineWidth=2, LineStyle=1, &
                       CharSize=1.0, BGColor=0)
        CALL PGPTXT (xPGmax+xPGstep*0.1, &
              yPGmax - (k+1.5)*yPGstep, &
              0.0, 0.0, "For the formation")
        CALL PGPTXT (xPGmax+xPGstep*0.1, &
              yPGmax - (k+2.5)*yPGstep, &
              0.0, 0.0, "and destruction of")
        CALL PGPTXT (xPGmax+xPGstep*0.1, &
              yPGmax - (k+3.5)*yPGstep, &
              0.0, 0.0, trim(nameSpeciesAna(i)))


        yPGmax = log10(maxval(abs(reactionAccum(:, &
                            intTmpvec(1):intTmpvec(2))))) + 0.2
        yPGmin = max(log10(minval(abs(reactionAccum(:, &
                            intTmpvec(1):intTmpvec(2))))), yPGmax-10.0)
        xPGstep = (xPGmax - xPGmin) / 10.0
        yPGstep = (yPGmax - yPGmin) / 25.0
        CALL setEnvPG (ColorIdx=1, LineWidth=3, LineStyle=1)
        CALL PGPAGE ! Creat a plot page
        CALL PGSVP (0.1, 0.65, 0.20, 1.0) ! page range
        CALL PGSWIN(xPGmin, xPGmax, yPGmin, yPGmax)
        CALL PGBOX ('BCSTL', 0.0, 0, 'BCNSTVL', 0.0, 0)
        CALL PGMTXT('L', 4.0, 0.5, 0.5, &
          'Absolute Reaction Rates (yr\u-1\d)')
       !CALL PGMTXT('B', 2.5, 0.5, 0.5, 'Time (year)')
        do k=1, min(nReactionContrib, 18)
          j = idxAccumSorted(k)
          if (dtMM(3, j, idxSpeciesAna) .GE. 0) then
            yPGPLOT = log10(reactionAccum(j, :))
          else
            yPGPLOT = log10(-reactionAccum(j, :))
          end if
          i1 = mod(k-1,18) + 1
          CALL PGSCR (16, RGB_R(i1), RGB_G(i1), RGB_B(i1))
          i1 = dtMM(2, j, idxSpeciesAna)
          strTMP = trim(combineStrArr(strReactants(1:2, i1), &
                2, "\(2284)"))//" \(2261) "//&
              trim(combineStrArr(strProducts(1:4, i1), 4, "\(2284)"))
          if (dtMM(3, j, idxSpeciesAna) .GE. 0) then
            CALL setEnvPG (ColorIdx=16, LineWidth=5, LineStyle=1, &
                    ifClip=1)
            strTMP = "F: "//strTMP
          else
            CALL setEnvPG (ColorIdx=16, LineWidth=5, LineStyle=2, &
                    ifClip=1)
            strTMP = "D: "//strTMP
          end if
          CALL PGLINE(nPGpoints, xPGPLOT, yPGPLOT)
          CALL setEnvPG (ColorIdx=16, LineWidth=2, LineStyle=1, &
                         CharSize=1.2, BGColor=0)
          CALL PGPTXT (xPGmax+xPGstep*0.1, &
              yPGmax - (k+0.25)*yPGstep, &
              0.0, 0.0, trim(strTMP))
        end do

        CALL setEnvPG (ColorIdx=1, LineWidth=2, LineStyle=1, &
                       CharSize=1.0, BGColor=0)
        CALL PGPTXT (xPGmax+xPGstep*0.1, &
              yPGmax - (k+1.5)*yPGstep, &
              0.0, 0.0, "For the formation")
        CALL PGPTXT (xPGmax+xPGstep*0.1, &
              yPGmax - (k+2.5)*yPGstep, &
              0.0, 0.0, "and destruction of")
        CALL PGPTXT (xPGmax+xPGstep*0.1, &
              yPGmax - (k+3.5)*yPGstep, &
              0.0, 0.0, trim(nameSpeciesAna(i)))

        yPGmax = log10(maxval(yHistory(idxSpeciesAna, &
                        intTmpvec(1):intTmpvec(2)))) + 0.2
        yPGmin = max(log10(minval(yHistory(idxSpeciesAna, &
                        intTmpvec(1):intTmpvec(2)))+1D-30), &
                yPGmax-10.0)
        CALL PGPAGE ! Creat a plot page
        CALL PGSVP (0.1, 0.65, 0.70, 1.20)
        CALL PGSWIN(xPGmin, xPGmax, yPGmin, yPGmax)
        CALL PGBOX ('BCSTL', 0.0, 0, 'BCNSTVL', 0.0, 0)
        CALL PGMTXT('L', 4.0, 0.5, 0.5, &
          'Number per Grain')
        yPGPLOT = log10(yHistory(idxSpeciesAna, :))
        CALL setEnvPG (ColorIdx=1, LineWidth=5, LineStyle=1, &
                         ifClip=1)
        CALL PGLINE(nPGpoints, xPGPLOT, yPGPLOT)
        CALL setEnvPG (ColorIdx=1, LineWidth=2, LineStyle=1, &
                       CharSize=1.0, BGColor=0)
        yPGstep = (yPGmax - yPGmin) / 10.0
        CALL PGPTXT (xPGmax+xPGstep*0.1, &
              yPGmax - 1.0*yPGstep, &
              0.0, 0.0, "The evolution of ")
        CALL PGPTXT (xPGmax+xPGstep*0.1, &
              yPGmax - 2.0*yPGstep, &
              0.0, 0.0, trim(nameSpeciesAna(i)))

        yPGPLOT = log10(abs(ratesSpecies(:,idxSpeciesAna)))
        yPGmax = maxval(yPGPLOT)
        yPGmin = max(minval(yPGPLOT), yPGmax-10.0)
        CALL PGPAGE ! Creat a plot page
        CALL PGSVP (0.1, 0.65, 1.20, 1.70)
        CALL PGSWIN(xPGmin, xPGmax, yPGmin, yPGmax)
        CALL PGBOX ('BCNSTL', 0.0, 0, 'BCNSTVL', 0.0, 0)
        CALL PGMTXT('L', 4.0, 0.5, 0.5, &
          'Changing rate (per year)')
        CALL PGMTXT('B', 2.5, 0.5, 0.5, 'Time (year)')
        where (ratesSpecies(:,idxSpeciesAna) .LT. 0D0)
          yPGPLOT = dblNaN()
        end where
        CALL setEnvPG (ColorIdx=1, LineWidth=5, LineStyle=1, &
                ifClip=1)
        CALL PGLINE(nPGpoints, xPGPLOT, yPGPLOT)
        yPGPLOT = log10(abs(ratesSpecies(:,idxSpeciesAna)))
        where (ratesSpecies(:,idxSpeciesAna) .GT. 0D0)
          yPGPLOT = dblNaN()
        end where
        CALL setEnvPG (ColorIdx=2, LineWidth=5, LineStyle=1, &
                ifClip=1)
        CALL PGLINE(nPGpoints, xPGPLOT, yPGPLOT)
        CALL PGEND
      end do

998   CALL CPU_TIME(ProgCurrentTime)
      write (*, '(/A, F10.3)') 'Seconds elapsed: ', &
        ProgCurrentTime - ProgStartTime

     !if (FileUnitOpened(fUOutLog)) then
     !  close (UNIT=fUOutLog, IOSTAT=ios, ERR=999, STATUS='KEEP')
     !end if

999   write (*,'(A, A60, A)') &
        char(27)//'[47m'//char(27)//'[34m'//&
        'Program end', '', char(27)//'[0m'

      end program



      subroutine setEnvPG (ColorIdx, LineWidth, LineStyle, &
        ifClip, CharSize, BGColor)
        INTEGER, INTENT (IN), OPTIONAL :: ColorIdx
        INTEGER, INTENT(IN), OPTIONAL :: LineWidth
        INTEGER, INTENT(IN), OPTIONAL :: LineStyle
        INTEGER, INTENT(IN), OPTIONAL :: ifClip
        REAL, INTENT(IN), OPTIONAL :: CharSize
        INTEGER, INTENT(IN), OPTIONAL :: BGColor
        if (present(ColorIdx)) &
          CALL PGSCI(ColorIdx) ! Color index
        if (present(LineWidth)) &
          CALL PGSLW(LineWidth) ! Line width
        if (present(LineStyle)) &
          CALL PGSLS(LineStyle) ! Line style
        if (present(ifClip)) &
          CALL PGSCLP(ifClip) ! Clip
        if (present(CharSize)) &
          CALL PGSCH(CharSize) ! Character size
        if (present(BGColor)) &
          CALL PGSTBG(BGColor) ! Background color
      end subroutine setEnvPG
