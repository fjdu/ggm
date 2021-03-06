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
!gfortran common_MC_MO.o common_MC_MO_TP.o subroutines_share_triv.o mo.o subroutines_share_proc.o subroutines_mo.o opkd*.o -fbounds-check -o mo
!
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

      CALL SaveMiscBefore

      do i=1, nReactions
        if (rates(i) .LT. rateThreshold) then
          rates(i) = 0D0
          nRealReactants(i) = 0
          nRealProducts(i) = 0
        end if
      end do

      if (IsWordChar(fSaveALLHistory(1:1))) then
        if (.NOT. getFileUnit(fUyHistory_bin)) then
          write (*,*) 'Cannot allocate an output file unit!'
          stop
        end if
        open (UNIT=fUyHistory_bin, &
             FILE=trim(path)//trim(fSaveALLHistory), &
             IOSTAT=ios, RECL=kind(1D0)*(nSpecies+1), &
             STATUS='REPLACE', ACCESS='DIRECT', &
             FORM='UNFORMATTED', ACTION='WRITE')
        if (ios .LT. 0) then
          write (*, *) 'OPEN FILE ERROR: IOSTAT=', ios
          stop
        end if
      end if
      if (IsWordChar(fSaveALLHistory_ASCII(1:1))) then
        if (.NOT. getFileUnit(fUyHistory_ascii)) then
          write (*,*) 'Cannot allocate an output file unit!'
          stop
        end if
        CALL openFileSequentialWrite (fUyHistory_ascii, &
          trim(path)//trim(fSaveALLHistory_ASCII), &
          kind('a')*(nSpecies+1)*constLenNameSpecies*2)
        write (FMTstryHistory, &
          FMT='("(A12, 3X, ", I4, "(A", I2, ", 2X))")') &
          nSpecies, constLenNameSpecies
        write (UNIT=fUyHistory_ascii, FMT=FMTstryHistory, IOSTAT=ios) &
          'Time_(yr)', nameSpecies(1:nSpecies)
        write (FMTstryHistory, &
          FMT='("(ES12.3, 2X, ", I4, "(ES", I2, ".4E3, 2X))")') &
          nSpecies, constLenNameSpecies
      end if
      if (IsWordChar(fPhyParHistory(1:1))) then
        if (.NOT. getFileUnit(fUPhyPar)) then
          write (*,*) 'Cannot allocate an output file unit!'
          stop
        end if
        CALL openFileSequentialWrite &
          (fUPhyPar, trim(path)//trim(fPhyParHistory), &
          kind('a')*(nPhyPar+1)*constLenNameSpecies*2)
        write (FMTstrPhyPar, &
          FMT='("(A12, 2X, ", I4, "(A", I2, ", 6X))")') &
          nPhyPar, constLenNameSpecies
        write (UNIT=fUPhyPar, FMT=FMTstrPhyPar, IOSTAT=ios) &
          "Time_(yr)", "Temperature", "n_H", "GrainRadius", "Av", &
            "rateCosIon", "omega_Albedo", "GrainDensity", "SitesDensity"
        write (FMTstrPhyPar, &
          FMT='("(ES12.3, 2X, ", I4, "(ES", I2, ".4E3, 6X))")') &
          nPhyPar, constLenNameSpecies
      end if

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
      do i=1, nIterMoment
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

      ! Save the moments.
      allocate (nameMoments(nNumMM))
      if (IsWordChar(fNameMoments(1:1))) then
        if (.NOT. getFileUnit(fU)) then
          write (*,*) 'Cannot allocate an output file unit!'
          stop
        end if
        CALL openFileSequentialWrite &
          (fU, trim(path)//trim(fNameMoments), 9999999)
        do i=1, nNumMM
          call vec2str(MM(1:nnzMM(i), i), nnzMM(i), &
            nameSpecies, strTMP)
          nameMoments(i) = strTMP
          write (fU, '(I6, A32, 3I6)') &
            i, trim(strTMP), nnzMM(i), ndtMM(i), ndauMM(i)
        end do
        close (UNIT=fU, IOSTAT=ios, STATUS='KEEP')
      end if

      ! Save the equations.
      i2 = 0
      if (IsWordChar(fEqMoments(1:1))) then
        if (.NOT. getFileUnit(fU)) then
          write (*,*) 'Cannot allocate an output file unit!'
          stop
        end if
        CALL openFileSequentialWrite &
          (fU, trim(path)//trim(fEqMoments), 9999999)
        do i=1, nNumMM
          if (ndtMM(i) .EQ. 0) cycle
          i2 = i2 + 1
          call vec2str(MM(1:nnzMM(i), i), nnzMM(i), &
            nameSpecies, strTMP)
          StrEqn = trim(strTMP) // '_t = '
          do j=1, ndtMM(i)
            i1 = dtMM(1, j, i)
            call vec2str(MM(1:nnzMM(i1), i1), nnzMM(i1), &
              nameSpecies, strTMP)
            if (dtMM(3, j, i) .EQ. 1) then
              write (strTMP1, '("+k(", I3.2, ")")') dtMM(2, j, i)
            else if (dtMM(3, j, i) .EQ. -1) then
              write (strTMP1, '("-k(", I3.2, ")")') dtMM(2, j, i)
            else
              write (strTMP1, '(I3, "k(", I3.2, ")")') &
                dtMM(3, j, i), dtMM(2, j, i)
            end if
            StrEqn = trim(StrEqn) // trim(strTMP1) // trim(strTMP)
          end do
          write (fU, '(2I6, 2X, A)') i2, i, trim(StrEqn)
        end do
        close (UNIT=fU, IOSTAT=ios, STATUS='KEEP')
      end if

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
      write (*, '(A32, I8)') 'Number of initial variables: ', NEQ

! Read initial conditions

      y = 0D0

      if (.NOT. getFileUnit(fU)) then
        write (*,*) 'Cannot allocate an output file unit!'
        stop
      end if
      CALL openFileSequentialRead &
          (fU, trim(path)//trim(fInitialCondition), 9999)
      CALL GetNLineF (fU, nLineAll, nLineData, commentChar)
      nInitialSpecies = nLineData
      write (*, '(/A32, I16)') 'Number of initial species: ', &
        nInitialSpecies
      rewind (UNIT=fU, IOSTAT=ios)
      do i=1, nLineAll
        read (UNIT=fU, FMT='(A128)', IOSTAT=ios) strTMP
        if ((strTMP(1:1) .EQ. commentChar) .OR. &
            (strTMP(1:1) .EQ. ' ')) cycle
        if (ios .NE. 0) then
          write (*,*) 'ios = ', ios; stop
        end if
        read (strTMP, FMT='(A12, ES10.2)', IOSTAT=ios) &
          nameSpecies_tmp, initialAbundance_tmp
        do j=1, nSpecies
          if (trim(nameSpecies(j)) .EQ. &
              trim(nameSpecies_tmp)) then
            if (nameSpecies(j)(1:1) .EQ. 'g') then
              y(j) = initialAbundance_tmp
            else
              !y(j) = VolumnPerGrain * (n_H * 0.5D0) &
              y(j) = VolumnPerGrain * n_H &
                * initialAbundance_tmp
              ! The abundance from the rate06 code is relative to H2,
              !   while the gas-to-dust ratio is based on H nuclei.
            end if
            exit
          end if
        end do
      end do
      close (UNIT=fU, IOSTAT=ios, ERR=999, STATUS='KEEP')

      do i=nSpecies+1, NEQ
        y(i) = getMMProd(indY(i), y, NEQ)
      end do

      IsStoSpecies = .FALSE.
      do i=1, nGrSpecies ! Only Surface species can be stochastic.
        i1 = idxGrSpecies(i)
        if ((.NOT. IsGas(i1)) .AND. (y(i1) .LT. sto_threshold)) then
          IsStoSpecies(i1) = .TRUE.
        end if
      end do

      do i=1, nElement
        initialElementalAb(i) = &
          dot_product(y(1:nSpecies), SpeciesElements(i, :))
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
      nIterationBase = 2000
      nTargetDecades = 9
      ratioTout = &
        exp(real(nTargetDecades) / real(nIterationBase) * log(10D0))

      t = 0.0D0 ! Start time
      tStep = 1.0D-4 ! First time step
      tout = t + tStep ! First output time
      touts(1) = t
      tFinal = tout + tStep*(ratioTout**nIteration-1D0)/(ratioTout-1D0)
      yHistory(:, 1) = y(1:nSpecies)
      ydotHistory(:, 1) = 0D0
      PhyParHistory(:, 1) = &
        (/Temperature, n_H, GrainRadius, Av, &
          rateCosIon, omega_Albedo, GrainDensity, SitesDensity/)
      if (FileUnitOpened(fUyHistory_bin)) &
        write (fUyHistory_bin, rec=1) touts(1), yHistory(:, 1)
      if (FileUnitOpened(fUyHistory_ascii)) &
        write (fUyHistory_ascii, FMT=FMTstryHistory) &
          touts(1), yHistory(:, 1)
      if (FileUnitOpened(fUPhyPar)) &
        write (fUPhyPar, FMT=FMTstrPhyPar) &
          touts(1), PhyParHistory(:, 1)

      IOPT = 1 ! 1: allow optional input; 0: disallow
      ITOL = 1 ! scalar control of tolerance
      ITASK = 1 ! normal, allow overshoot
      ISTATE = 1 ! first call
      MF = 22 !021 ! 221 also works. 021: use Jac; 022: not.

      nerr = 0 ! Count the number of errors during iteration.
      nChange = 0
      nChangeSpe = 0
      idxChangeSpe = 0
      flagChange = .TRUE.
      countMMNotUsed = 0
      flagMMNotUsed = .FALSE.

      allocate(iWhenDropped(nNumMM))
      iWhenDropped = 0

      write (*, '(/A, ES10.2, A, ES10.2)') &
        'Time span (years): ', t*timeUnit, ' -> ', tFinal*timeUnit

      write (*, '(/A//)') 'Now start to evolve...'

      do i=1, nIteration
        call date_and_time(DATE=strTMP, TIME=strTMP1)
        write (*, '(A, I6, " (", F5.1, "%)  tStep = ", &
          & ES10.2, " t = ", ES10.2, " @ ", A8, 2X, A10)') &
          CHAR(27)//'[AIterating... ', &
          i, real(i*100)/nIteration, tStep, tout,  strTMP, strTMP1
        ! Check whether the physical parameter has changed
        if (getPhyPar(t)) then
          write (*,*) 'Updating the rates...'
          CALL CalcMobilities
          CALL CalcRates
          rates = rates * SecondsPerYear
          do i1=1, nNumMM
            if (ndtMM(i1) .GT. 0) then
              do j=1, ndtMM(i1)
                dtMMMo(i1)%coeff(j) = rates(dtMM(2, j, i1)) * &
                  DBLE(dtMM(3, j, i1))
              end do
            end if
          end do
          ISTATE = 1
        end if
        if (flagChange) then
          ISTATE = 1
          intTmpvec(1:3) = (/NEQ, LRW, LIW/)
          do j=1, nNumMM
            if (IsY(j)) then
              if (nnzMM(j) .GT. 1) then
                if (getFstDetSpe(j, y, NEQ) .GT. 0) then
                  IsY(j) = .FALSE.
                  iWhenDropped(j) = i
                end if
              end if
            end if
            if (IsY(j)) yBak(j) = y(invIndY(j))
            !+ This part is new....
            if (.NOT. IsY(j)) then
              if ((getFstDetSpe(j, y, NEQ) .EQ. 0) .AND. &
                  (nnzMM(j) .LE. nOrderLim)) then
                if ((i - iWhenDropped(j)) .GT. 50) then
                  IsY(j) = .TRUE.
                  yBak(j) = getMMProd(j, y, NEQ)
                end if
              end if
            end if
            !- This part is new....
          end do
          NEQ = count(IsY)
          if (NEQ .GT. intTmpvec(1)) then
            deallocate(indY, y)
            allocate(indY(NEQ), y(NEQ))
            if (allocated(sparseMaskJac)) deallocate(sparseMaskJac)
            allocate(sparseMaskJac(NEQ, NEQ))
          end if
          if (.NOT. allocated(sparseMaskJac)) &
            allocate(sparseMaskJac(NEQ, NEQ))
          i1 = 0
          do j=1, nNumMM
            if (IsY(j)) then
              i1 = i1 + 1
              indY(i1) = j
              invIndY(j) = i1
              y(i1) = yBak(j)
            else
              invIndY(j) = 0
            end if
          end do
          call makeSparse(y, NEQ)
          NNZ = COUNT(sparseMaskJac)
          LRW = 20 + NEQ * (12 + 1) + 3 * NEQ + &
            4 * NNZ + 2 * NEQ + (NNZ + 10 * NEQ)
          LRW = LRW * 2
          LIW = 31 + NEQ + NNZ

          if (LRW .GT. intTmpvec(2)) then
            if (allocated(RWORK)) deallocate(RWORK)
            allocate(RWORK(LRW))
          end if
          if (LIW .GT. intTmpvec(3)) then
            if (allocated(IWORK)) deallocate(IWORK)
            allocate(IWORK(LIW))
          end if
          if (.NOT. allocated(RWORK)) allocate(RWORK(LRW))
          if (.NOT. allocated(IWORK)) allocate(IWORK(LIW))

          RWORK(5:10) = 0D0
          IWORK(5:10) = 0
          IWORK(5) = 12
          IWORK(6) = 5000
          IWORK(31) = 1
          k = 1
          do i1=1, NEQ
            do j=1, NEQ
              if (sparseMaskJac(j, i1)) then
                IWORK(31 + NEQ + k) = j
                k = k + 1
              end if
            end do
            IWORK(31+i1) = k
          end do
          write (*, &
            FMT='(/A32, I16, 2X, "(", F6.2, "%)", A16, I6/)') &
            'Number of nonzeros in Jac: ', NNZ, &
            real(NNZ*100)/(NEQ*NEQ), 'NEQ: ', NEQ
          write (*, *) "Waiting for the solver..."

          deallocate(sparseMaskJac, STAT=statALLOC)
          flagChange = .FALSE.
          flagMMNotUsed = .FALSE.
        end if

        call DLSODES &
             (f, NEQ, y, t, tout, &
              ITOL, RTOL, ATOL, &
              ITASK, ISTATE, IOPT, &
              RWORK, LRW, IWORK, LIW, &
              jac, MF)
        if (ExternalAction(t, y, NEQ)) then
          write (*,*) 'External action...'
          ISTATE = 1
        end if

        CALL DINTDY (T, 1, RWORK(IWORK(22)), NEQ, ydotTMP(1:NEQ), &
          intTmpvec(1))
        ydotHistory(:, i+1) = ydotTMP(1:nSpecies)

        ! Saving history
        touts(i+1) = t
        yHistory(:, i+1) = y(1:nSpecies)
        PhyParHistory(:, i+1) = &
          (/Temperature, n_H, GrainRadius, Av, &
            rateCosIon, omega_Albedo, GrainDensity, SitesDensity/)
        if (FileUnitOpened(fUyHistory_bin)) &
          write (fUyHistory_bin, rec=(i+1)) touts(i+1), yHistory(:, i+1)
        if (FileUnitOpened(fUyHistory_ascii)) &
          write (fUyHistory_ascii, FMT=FMTstryHistory) &
            touts(i+1), yHistory(:, i+1)
        if (FileUnitOpened(fUPhyPar)) &
          write (fUPhyPar, FMT=FMTstrPhyPar) &
            touts(i+1), PhyParHistory(:, i+1)
        FLUSH (fUyHistory_bin)
        FLUSH (fUyHistory_ascii)
        FLUSH (fUPhyPar)

        if (ISTATE .LT. 0) then
          nerr = nerr + 1
          write (*, '(3(A, I4, 4X)/)') &
            'Step: ', i, "Error number: ", &
            nerr, 'Error code = ', ISTATE
          ISTATE = 3
        end if

        if (mod(i, 50) .EQ. 0) then
          do j=nSpecies+1, NEQ
            i1 = indY(j)
            dblTmpVec(1) = getMMProd(i1, y, NEQ)
            if (((y(j) .GT. 2D0*dblTmpVec(1)) .AND. &
                 (getMaxYMM(i1, y, NEQ) .GT. 0.5D0)) &
                .OR. (y(j) .LT. 0D0) &
                .OR. (abs(y(j)-dblTmpVec(1)) .LT. 1D-1*y(j))) then
              IsY(i1) = .FALSE.
              iWhenDropped(i1) = i
              flagChange = .TRUE.
              write (fUOutLog, '(A8, I8, 2A32, 2ES13.4)') 'Step: ', i, &
                'Ruling out ', nameMoments(i1), y(j), dblTmpVec(1)
            end if
          end do
        end if

        do j=1, nGrSpecies
          i1 = idxGrSpecies(j)
          if (IsEndProd(i1)) cycle ! End species is never a component.
          if ((yHistory(i1, i) .LT. sto_threshold) .AND. &
              (y(i1) .GE. sto_threshold)) then
            if (nChangeSpe(i1) .GT. 0) then
              if ((i - idxChangeSpe(nChangeSpe(i1), i1)) .LT. 5) then
                cycle
              end if
            end if
            IsStoSpecies(i1) = .FALSE.
            nChange = nChange + 1
            nChangeSpe(i1) = nChangeSpe(i1) + 1
            idxChangeSpe(nChangeSpe(i1), i1) = i
            flagChange = .TRUE.
            write (*, '(/, I6, A25, A12, F8.2, A, F8.2, /)') &
              i, ' Changing upwards: ', nameSpecies(i1), &
              yHistory(i1, i), " -> ", y(i1)
            cycle
          end if
          if ((yHistory(i1, i) .GE. sto_threshold*1D0) .AND. &
              (y(i1) .LT. sto_threshold*1D0)) then
            if (nChangeSpe(i1) .GT. 0) then
              if ((i - idxChangeSpe(nChangeSpe(i1), i1)) .LT. 5) then
                cycle
              end if
            end if
            IsStoSpecies(i1) = .TRUE.
            nChange = nChange + 1
            nChangeSpe(i1) = nChangeSpe(i1) + 1
            idxChangeSpe(nChangeSpe(i1), i1) = i
            flagChange = .TRUE.
            write (*, '(/, I6, A25, A12, F8.2, A, F8.2, /)') i, &
              ' Changing downwards: ', nameSpecies(i1), &
              yHistory(i1, i), " -> ", y(i1)
            cycle
          end if
        end do

        tStep = tStep * ratioTout
        tout = t + tStep
      end do

      write (*, '(/A32, I6, "(", F6.2, "%)")') &
        'Number of errors: ', &
        nerr, real(nerr*100)/real(nIteration)
      write (*, '(A32, I6)') 'Number of changes: ', nChange

      !touts = touts * timeUnit ! useless

      if (FileUnitOpened(fUyHistory_bin)) &
        close (UNIT=fUyHistory_bin, IOSTAT=ios, ERR=999, STATUS='KEEP')
      if (FileUnitOpened(fUyHistory_ascii)) &
        close (UNIT=fUyHistory_ascii, IOSTAT=ios, ERR=999, &
          STATUS='KEEP')
      if (FileUnitOpened(fUPhyPar)) &
        close (UNIT=fUPhyPar, IOSTAT=ios, ERR=999, STATUS='KEEP')

      if (IsWordChar(fFinalJac(1:1))) then
        if (.NOT. getFileUnit(fU)) then
          write (*,*) 'Cannot allocate an output file unit!'
          stop
        end if
        CALL openFileSequentialWrite &
          (fU, trim(path)//trim(fFinalJac), 9999999)
        call printjac(fU, NEQ, t, y)
        close (UNIT=fU, IOSTAT=ios, STATUS='KEEP')
      end if

      if (IsWordChar(fSaveElementalAb(1:1))) then
        if (.NOT. getFileUnit(fU)) then
          write (*,*) 'Cannot allocate an output file unit!'
          stop
        end if
        do i=1, nElement
          finalElementalAb(i) = &
            dot_product(y(1:nSpecies), SpeciesElements(i, :))
        end do
        CALL openFileSequentialWrite &
          (fU, trim(path)//trim(fSaveElementalAb), 999)
        write (fU, '(A)') '! Elemental conservation:'
        write (fU, '(A12, A17, A17, A17)') &
            '!    Element', 'Initial', 'Final', '(F-I)/I'
        do i=1, nElement
          write (fU, '(A12, ES17.8E2, ES17.8E2, ES17.8E2)') &
              adjustr(nameElements(i)), &
              initialElementalAb(i), finalElementalAb(i), &
              (finalElementalAb(i) - initialElementalAb(i)) &
              / (abs(initialElementalAb(i)) + 1D-30)
        end do
        close (UNIT=fU, IOSTAT=ios, ERR=999, STATUS='KEEP')
      end if

      if (IsWordChar(fAnalyseResult(1:1))) then
        CALL AnalyseReaction
      end if

      if (IsWordChar(fSaveFinalResult(1:1))) then
        if (.NOT. getFileUnit(fU)) then
          write (*,*) 'Cannot allocate an output file unit!'
          stop
        end if
        CALL openFileSequentialWrite &
          (fU, trim(path)//trim(fSaveFinalResult), 999)
        do i=1, nSpecies
          write (UNIT=fU, FMT='(A12, 2ES16.6E3)', IOSTAT=ios) &
            nameSpecies(i), yHistory(i, nHistoryLen), &
            ydotHistory(i, nHistoryLen)
        end do
        close (UNIT=fU, IOSTAT=ios, ERR=999, STATUS='KEEP')
      end if

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
