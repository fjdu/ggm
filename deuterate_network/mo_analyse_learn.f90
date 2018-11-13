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
!gfortran common_MC_MO.o common_MC_MO_TP.o subroutines_share_triv.o mo_analyse.o subroutines_share_proc.o subroutines_mo.f90 opkd*.o -fbounds-check -lpgplot -lX11 -o mo_analyse
!cp ../src/mo_analyse ./
!./mo_analyse proj_Ext_10.0_2.0E+05_1.0E-07_15.0/config_moment_withTime.dat
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
        nElementReac, nElementProd

      integer nIterationBase, nTargetDecades, nerr, nChange
      double precision ratioTout, tStep
      double precision t, tout, tFinal
      double precision, dimension(:), allocatable :: y, yBak, RWORK
      double precision, dimension(:), allocatable :: ydotTMP
      integer IOPT, ISTATE, ITASK, ITOL, LIW, LRW, NEQ, MF, NNZ
      integer, dimension(:), allocatable :: IWORK

      character(len=256), dimension(dimvecMon) ::  StrPolynomial
      character(len=16384) StrEqn

      integer, dimension(:), allocatable :: iWhenDropped

      logical getPhyPar, ExternalAction, getFileUnit, FileUnitOpened
      integer, parameter :: nPhyPar = 8

      integer getIdxSpecies, idxSpeciesAna
      integer, parameter :: nSpeciesAna = 10
      character(len=constLenNameSpecies), dimension(nSpeciesAna) :: &
        nameSpeciesAna = &
          (/'H3O+', 'H2O+', 'H3+', 'E-', &
            'HCO+', 'CO', 'H2O', 'HCN', &
            'HO2', 'gHO2'/)
      integer nReactionContrib
      double precision, dimension(:,:), allocatable :: reactionContrib
      double precision, dimension(:), allocatable :: reactionContribAve
      integer, dimension(:), allocatable :: idxYdotOneReactionSorted
      
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
      RGB_R = (/1.0, 0.0, 0.0, 0.8, 0.0, 1.0, 0.8, 0.0, 0.0, 0.5, 0.0, 0.8, 1.0, 0.5, 0.7, 1.0, 0.7, 0.5/)
      RGB_G = (/0.0, 0.8, 0.0, 0.8, 0.8, 0.0, 0.0, 0.5, 0.0, 0.5, 0.5, 0.0, 0.7, 0.8, 0.5, 0.5, 0.8, 0.5/)
      RGB_B = (/0.0, 0.0, 1.0, 0.0, 0.8, 1.0, 0.0, 0.0, 0.8, 0.0, 0.5, 0.8, 0.5, 0.7, 1.0, 0.8, 0.5, 1.0/)

      CALL CPU_TIME(ProgStartTime)
      call date_and_time(DATE=strTMP, TIME=strTMP1)

      write (*,'(A, A60, A)') &
        char(27)//'[47m'//char(27)//'[34m'//&
        'Program start', '', char(27)//'[0m'
      write (*,'("Current time: ", A8, 2X, A10)') strTMP, strTMP1

      CALL GET_COMMAND_ARGUMENT(1, strTMP, i, ios)
      if (i .EQ. 0) strTMP = 'config.dat'

      ! This initialization imports, and only imports all the configuration
      ! parameters.
      write (*, '(/A)') 'Initializing...'
      CALL initialize (strTMP)

      ! Find a file unit for message output.
      if (IsWordChar(fOutLog(1:1))) then
        if (.NOT. getFileUnit(fUOutLog)) then
          write (*,*) 'Cannot allocate an output file unit!'
          stop
        end if
        CALL openFileSequentialWrite &
          (fUOutLog, trim(path)//trim(fOutLog), 9999999)
        CALL XSETUN(fUOutLog)
      else
        do i=16, 2, -1
          if (ISATTY(unit=i)) then
            INQUIRE(UNIT=i, action=strTMP)
            if (trim(strTMP) .EQ. 'WRITE') then
              fUOutLog = i; exit
            end if
          end if
        end do
      end if

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
         typeReac(nReactions), &
         rates(nReactions), STAT=statALLOC)

      CALL ReadReactions (fU, nLineAll, commentChar)
      close (UNIT=fU, IOSTAT=ios, ERR=999, STATUS='KEEP')

      CALL MakeReactions

      allocate &
        (SpeciesElements(nElement, nSpecies), &
         massSpecies(nSpecies), &
         mobilitySpecies(nSpecies), &
         initialElementalAb(nElement), &
         finalElementalAb(nElement), &
         nElementReac(nElement), &
         nElementProd(nElement), &
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
         IsEndProd(nSpecies), STAT=statALLOC)

      nGrSpecies = 0
      do i=1, nSpecies
        CALL getElements(nameSpecies(i), nameElements, nElement, &
          SpeciesElements(:, i))
        massSpecies(i) = sum(SpeciesElements(:, i) * ElementMassNumber)
        if (nameSpecies(i)(1:1) .EQ. 'g') then
          nGrSpecies = nGrSpecies + 1
          idxGrSpecies(nGrSpecies) = i
          IsGas(i) = .FALSE.
        else
          IsGas(i) = .TRUE.
        end if
      end do

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

        nElementReac = 0
        nElementProd = 0
        do j=1, nRealReactants(i)
          nElementReac = nElementReac + &
            SpeciesElements(:, reactions(j, i))
        end do
        do j=1, nRealProducts(i)
          nElementProd = nElementProd + &
            SpeciesElements(:, reactions(nReactants+j, i))
        end do
        if ((abs(nElementReac(1) - nElementProd(1)) + &
          sum(abs(nElementReac(3:nElement) - &
          nElementProd(3:nElement)))) .NE. 0) then
          write (*, '(2A, I6, A, 2X, 2A12, " -> ", 5A12)') &
            'Elements not conserved [discarded]: ', &
            char(27)//'[41m', i, char(27)//'[0m', &
            strReactants(1:nReactants, i), strProducts(1:nProducts, i)
          nRealReactants(i) = 0
          nRealProducts(i) = 0
        end if
        do j=1, i-1
          if ((sum(abs(reactions(:, j)-reactions(:, i))) .EQ. 0) &
              .AND. (typeReac(j) .EQ. typeReac(i))) then
            write (*,'(2A, 2I4, A)') &
              'Duplicate reaction pair: ', &
              char(27)//'[45m', i, j, char(27)//'[0m'
          end if
        end do
      end do
      deallocate (nElementReac, nElementProd, STAT=statALLOC)

      write (*, FMT='(6(/, A32, I16))') &
        'Number of species: ', nSpecies, &
        'Number of gas species: ', count(IsGas), &
        'Number of reactions: ', nReactions, &
        'Number of surface reactions: ', nGrReactions, &
        'Max number of reactants: ', maxval(nRealReactants), &
        'Max number of products: ', maxval(nRealProducts)

! Read the mobilities.

      write (*, '(/A)') 'Importing mobilities...'

      CALL CalcMobilities

      write (*, '(/A)') 'Calculating rates...'

      CALL CalcRates
      rates = rates * SecondsPerYear

      write (*, '(3(/A32, ES16.6))') &
        'Grain to H number ratio: ', ratioGrainToH, &
        'Sites per grain: ', SitesPerGrain, &
        'Volumn containing one grain: ', VolumnPerGrain

      do i=1, nReactions
        if (rates(i) .LT. rateThreshold) then
          rates(i) = 0D0
          nRealReactants(i) = 0
          nRealProducts(i) = 0
        end if
      end do

! Now I need to prepare all the possibly needed moments.

      write (*, '(A)') 'Making moments...'
      ! 128 is for a small network; 4 is for bigger one.
      !nNumMMEst = nSpecies*16 + nReactions * 128
      nNumMMEst = nSpecies*16 + nReactions * 4
      allocate (&
        polyReac(nReactions), &
        MM(dimvecDeriv, nNumMMEst), &
        nnzMM(nNumMMEst), &
        dtMM(3, nReactions*4, nNumMMEst), & ! The number of terms might
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
            if (getFstGasSpe(j) .NE. 0) cycle
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

      IsEndProd = .TRUE.
      do i=1, nGrSpecies
        i1 = idxGrSpecies(i)
        if (IsEndProd(i1)) then
          do j=nSpecies+1, nNumMM
            if (IsEndProd(i1)) then
              do k=1, nnzMM(j)
                if (MM(k, j) .EQ. i1) then
                  IsEndProd(i1) = .FALSE.
                  exit
                end if
              end do
            end if
          end do
        end if
      end do

      deallocate(polyReac, indMM)

      allocate(IsY(nNumMM), countMMNotUsed(nNumMM))
      IsY(1:nSpecies) = .TRUE.
      do i=nSpecies+1, nNumMM
        if (nnzMM(i) .GT. nOrderLim) then
          IsY(i) = .FALSE.
        else
          if (getFstGasSpe(i) .GT. 0) then
            IsY(i) = .FALSE.
          else
            IsY(i) = .TRUE.
          end if
        end if
      end do

      NEQ = count(IsY)
      allocate(indY(NEQ), invIndY(nNumMM), y(NEQ), yBak(nNumMM), &
        dtMMMo(nNumMM), ydotTMP(nNumMM))
      i1 = 0
      do i=1, nNumMM
        if (IsY(i)) then
          i1 = i1 + 1
          indY(i1) = i
          invIndY(i) = i1
        else
          invIndY(i) = 0
        end if
        if (ndtMM(i) .GT. 0) then
          allocate(dtMMMo(i)%coeff(ndtMM(i)), STAT=statALLOC)
          do j=1, ndtMM(i)
            dtMMMo(i)%coeff(j) = rates(dtMM(2, j, i)) * &
              DBLE(dtMM(3, j, i))
          end do
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

      if (nIteration .LT.  0) then
844     write (*, FMT='(/A)') 'Number of iteration (please input):'
        READ (*, '(I8)', ERR=845, IOSTAT=ios) nIteration
845     if (ios .NE. 0) then
          write (*,'(A)') &
            char(27)//'[43m'//char(27)//'[31m'//&
            'Invalid Input!'//char(27)//'[0m'
          go to 844
        end if
      end if

      write (*, '(/A, I8)') 'Number of iteration: ', nIteration
      if (nIteration .LT. 1) go to 999

      nHistoryLen = nIteration + 1
      allocate(yHistory(nSpecies, nHistoryLen), &
               ydotHistory(nSpecies, nHistoryLen), &
               PhyParHistory(nPhyPar, nHistoryLen), &
               touts(nHistoryLen), STAT=statALLOC)

      if (IsWordChar(fSaveALLHistory(1:1))) then
        if (.NOT. getFileUnit(fUyHistory_bin)) then
          write (*,*) 'Cannot allocate an output file unit!'
          stop
        end if
        open (UNIT=fUyHistory_bin, &
             FILE=trim(path)//trim(fSaveALLHistory), &
             IOSTAT=ios, RECL=kind(1D0)*(nSpecies+1), &
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

      if (.NOT. getFileUnit(fU)) then
        write (*,*) 'Cannot allocate an output file unit!'
        stop
      end if
      if (.NOT. IsWordChar(fAnalyseResult(1:1))) then
        write (*,*) fAnalyseResult, ' is not provided!'
        stop
      end if
      allocate (idxYdotOneReactionSorted(nReactions))
      do i=1, nSpeciesAna
        idxSpeciesAna = getIdxSpecies(nameSpeciesAna(i))
        nReactionContrib = ndtMM(idxSpeciesAna)
        if (allocated(reactionContrib)) deallocate(reactionContrib)
        if (allocated(reactionContribAve)) &
          deallocate(reactionContribAve)
        allocate(reactionContrib(nReactionContrib, nHistoryLen))
        allocate(reactionContribAve(nReactionContrib))
        CALL openFileSequentialWrite &
          (fU, trim(path)//trim(getFilePreName(fAnalyseResult))//&
            "_"//trim(nameSpeciesAna(i))//".dat", 999)
        write (fU, '(A, 4X, "#", I6, 4X, "nReaction: ", I6)') &
          nameSpeciesAna(i), idxSpeciesAna, nReactionContrib
        do j=1, nHistoryLen
          if (getPhyPar(touts(j))) then
            CALL CalcMobilities
            CALL CalcRates
            rates = rates * SecondsPerYear
            do i1=1, nNumMM
              if (ndtMM(i1) .GT. 0) then
                do i2=1, ndtMM(i1)
                  dtMMMo(i1)%coeff(i2) = rates(dtMM(2, i2, i1)) * &
                    DBLE(dtMM(3, i2, i1))
                end do
              end if
            end do
          end if
          do k=1, nReactionContrib
            reactionContrib(k, j) = dtMMMo(idxSpeciesAna)%coeff(k) * &
              getMMProd(dtMM(1, k, idxSpeciesAna), &
                yHistory(:,j), nSpecies)
          end do
          if (j .EQ. 1) cycle
          CALL SORT_Asc_idx(nReactionContrib, &
            -abs(reactionContrib(:, j)), idxYdotOneReactionSorted)
          write (fU, '(2X, "[", I0, "] ", "At time: ", ES10.2, 2X, &
            "Abundance: ", ES10.2, 2X, "Total rate: ", Es10.2)') &
            j, touts(j), yHistory(idxSpeciesAna, j), &
            sum(reactionContrib(:, j))
          do k=1, min(nReactionContrib, 50)
            i1 = idxYdotOneReactionSorted(k)
            i2 = dtMM(2, i1, idxSpeciesAna)
            write (fU, &
              '(4X, I3, ES11.2, 2X, I5, ": " &
              & 7A8, 2X, ES15.4, 4X, ES8.2, 2F8.2)') &
              k, reactionContrib(i1, j), i2, &
              strReactants(1:2, i2), ' -> ', strProducts(1:4, i2), &
              dtMMMo(idxSpeciesAna)%coeff(i1), dblABC(:, i2)
          end do
        end do
        close (UNIT=fU, IOSTAT=ios, ERR=999, STATUS='KEEP')

        reactionContribAve = sum(reactionContrib, 2)
        CALL SORT_Asc_idx(nReactionContrib, &
            -abs(reactionContribAve), idxYdotOneReactionSorted)

        write (*,'(/A, I6, 2X, A)') 'Plotting...', i, nameSpeciesAna(i)
        j = PGBEG(0, trim(path)//&
          trim(getFilePreName(fAnalyseResult))//"_"// &
          trim(nameSpeciesAna(i))//'.ps/CPS', 1, 1)
        IF (j .NE. 1) STOP
        nPGpoints = nHistoryLen
        allocate(xPGPLOT(nPGpoints), yPGPLOT(nPGpoints), STAT=statALLOC)
        xPGPLOT = log10(touts + 1D-10)
        xPGmin = 2.0 !xPGPLOT(2)
        xPGmax = xPGPLOT(nPGpoints)
        yPGmax = log10(maxval(abs(reactionContrib)))
        yPGmin = max(log10(minval(abs(reactionContrib))), yPGmax-7.0)
        xPGstep = (xPGmax - xPGmin) / 10.0
        yPGstep = (yPGmax - yPGmin) / 25.0
        CALL setEnvPG (ColorIdx=1, LineWidth=3, LineStyle=1)
        CALL PGPAGE ! Creat a plot page
        CALL PGSVP (0.15, 0.7, 0.15, 0.80) ! page range
        CALL PGSWIN(xPGmin, xPGmax, yPGmin, yPGmax)
        CALL PGBOX ('BCNSTL', 0.0, 0, 'BCNSTVL', 0.0, 0)
        CALL PGMTXT('L', 4.0, 0.5, 0.5, &
          'Absolute Reaction Rates (yr\u-1\d)')
        CALL PGMTXT('B', 2.5, 0.5, 0.5, 'Time (year)')
        do k=1, min(nReactionContrib, 18)
          j = idxYdotOneReactionSorted(k)
          if (dtMMMo(idxSpeciesAna)%coeff(j) .GE. 0D0) then
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
          CALL setEnvPG (ColorIdx=15, LineWidth=16, LineStyle=1, &
                         ifClip=0)
          !CALL PGMOVE(xPGmax+xPGstep*2.1, yPGmax - k*yPGstep)
          !CALL PGDRAW(xPGmax+xPGstep*3.2, yPGmax - k*yPGstep)
          i1 = dtMM(2, j, idxSpeciesAna)
          strTMP = trim(combineStrArr(strReactants(1:2, i1), &
                2, "\(2284)"))//" \(2261) "//&
              trim(combineStrArr(strProducts(1:4, i1), 4, "\(2284)"))
          if (dtMMMo(idxSpeciesAna)%coeff(j) .GE. 0D0) then
            CALL setEnvPG (ColorIdx=16, LineWidth=5, LineStyle=1)
            strTMP = "F: "//strTMP
          else
            CALL setEnvPG (ColorIdx=16, LineWidth=5, LineStyle=2)
            strTMP = "D: "//strTMP
          end if
          !CALL PGMOVE(xPGmax+xPGstep*2.1, yPGmax - k*yPGstep)
          !CALL PGDRAW(xPGmax+xPGstep*3.2, yPGmax - k*yPGstep)
          CALL setEnvPG (ifClip=1)
          CALL PGLINE(nPGpoints, xPGPLOT, yPGPLOT)
          CALL setEnvPG (ColorIdx=16, LineWidth=3, LineStyle=1, &
                         CharSize=0.6, BGColor=0)
          CALL PGPTXT (xPGmax+xPGstep*0.1, &
              yPGmax - (k+0.25)*yPGstep, &
              0.0, 0.0, trim(strTMP))
        end do
        CALL setEnvPG (ColorIdx=1, LineWidth=3, LineStyle=1, &
                       CharSize=0.7, BGColor=0)
        CALL PGPTXT (xPGmax+xPGstep*0.1, &
              yPGmax - (k+1.5)*yPGstep, &
              0.0, 0.0, "For the formation")
        CALL PGPTXT (xPGmax+xPGstep*0.1, &
              yPGmax - (k+2.5)*yPGstep, &
              0.0, 0.0, "and destruction of")
        CALL PGPTXT (xPGmax+xPGstep*0.1, &
              yPGmax - (k+3.5)*yPGstep, &
              0.0, 0.0, trim(nameSpeciesAna(i)))
        CALL PGEND
      end do

998   CALL CPU_TIME(ProgCurrentTime)
      write (*, '(/A, F10.3)') 'Seconds elapsed: ', &
        ProgCurrentTime - ProgStartTime

      if (FileUnitOpened(fUOutLog)) then
        close (UNIT=fUOutLog, IOSTAT=ios, ERR=999, STATUS='KEEP')
      end if

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
