! Note added on 2010-09-01: I remember this one is faster than MC_acc.f90.
! For larger networks, both of them should be tried.
!
! 2010-09-06:
!   Adapted for gas grain.
! 2010-09-07:
!   Be cause some very abundant species like H2 and He are included,
!     I choose to use kind=8 integers to hold the populations.
!     So for a 64 bit integer, it is able to hold numbers up to
!       2^63 ~ 9E18, more than enough.
!   Droped this line:
!       t_Save = t_Save/SecondsPerYear
!     This is because previous;y the unit used for the rates are in seconds.
!     Now the MC and moment share the same code for calculating the rates,
!       so the rates are always expressed in years.
! 2010-11-16 Tue 21:57:41 
!   Save the final state into a binary file for potential later use.
!
! gfortran common_MC_MO.o subroutines_MC.o sub_pre_share.o MC.o subroutines_share.o common_MC_MO_TP.o -fbounds-check -o MC

      program MC
      use CMDT
      use CMTP
      implicit none

      integer fU, ios, statALLOC, lenArgCMD
      integer nLineAll, nLineData
      character commentChar, strArgCMD*128, strTMP*128, FMTstr*128
      character(LEN=constLenNameSpecies)  nameInitialSpecies_tmp
      logical IsWordChar, flag

      integer(kind=8) i, j, k, h, i1, i2, iReactionNext
      double precision RandProbTmp, dblTmp
      double precision, dimension(2) :: rand_tmp_vec
      double precision, dimension(3) :: dblTmp_3
      double precision TotalPropensity, AccumPropensity
      double precision initialAbundance_tmp
      integer, dimension(16) :: intTmpvec

      integer ClockCount0, ClockCount1, CountRate

      double precision, dimension(:), allocatable :: &
        initialElementalAb, finalElementalAb
      integer, dimension(:), allocatable :: &
        nElementReac, nElementProd

      double precision tNext, t_Accum
      double precision, dimension(:), allocatable :: PopSpeciesAccum
      integer(kind=8), dimension(:), allocatable :: PopSpeciesInit
      integer(kind=8) nStepPerSave

      CALL SYSTEM_CLOCK (ClockCount0, CountRate)
      call date_and_time(DATE=strTMP, TIME=FMTstr)
      
      write (*,'(A, A60, A)') &
        char(27)//'[47m'//char(27)//'[34m'//&
        'Program start', '', char(27)//'[0m'
      write (*,'("Current time: ", A8, 2X, A10)') strTMP, FMTstr

      fU = 1

      commentChar = '!' ! A blank space also means comment.

      CALL GET_COMMAND_ARGUMENT(1, strArgCMD, lenArgCMD, ios)
      if (lenArgCMD .EQ. 0) &
        strArgCMD = 'config.dat'

      write (*,'(/A)') 'Initializing...'
      CALL initialize (fU, strArgCMD)

! Import all the reactions

      write (*,'(/A)') 'Importing reactions...'

      CALL openFileSequentialRead &
        (fU, trim(path)//trim(fReactions), 999)

      CALL GetNLineF (fU, nLineAll, nLineData, commentChar)

      nReactions = nLineData

      allocate &
        (strReactants(nReactants, nReactions), &
         strProducts(nProducts, nReactions), &
         nameSpecies(nSpecies_Est), &
         nRealReactants(nReactions), nRealProducts(nReactions), &
         reactions(nParticipants, nReactions), &
         dblABC(3, nReactions), typeReac(nReactions), &
         rates(nReactions), propensities(nReactions), STAT=statALLOC)

      CALL ReadReactions (fU, nLineAll, commentChar)

      close (UNIT=fU, IOSTAT=ios, ERR=999, STATUS='KEEP')

      CALL MakeReactions

      allocate &
        (SpeciesElements(nElement, nSpecies), &
         massSpecies(nSpecies), mobilitySpecies(nSpecies), &
         initialElementalAb(nElement), finalElementalAb(nElement), &
         nElementReac(nElement), nElementProd(nElement), &
         PopSpecies(nSpecies), PopSpeciesInit(nSpecies), &
         PopSpeciesSave(nSpecies, nRecordsSave), &
         IsReactionChanged(nReactions), idxSpeUpdated(8), &
         idxReac_AsReactants(nReactions, nSpecies), &
         idxReac_AsProducts(nReactions, nSpecies), &
         nReac_AsReactants(nSpecies), nReac_AsProducts(nSpecies), &
         t_Save(nRecordsSave), PopSpeciesAccum(nSpecies), &
         STAT=statALLOC)

      do i=1, nSpecies
        CALL getElements(nameSpecies(i), nameElements, nElement, &
           SpeciesElements(:, i))
        massSpecies(i) = sum(SpeciesElements(:, i) * ElementMassNumber)
      end do

      nReac_AsReactants = 0
      nReac_AsProducts = 0

      do i=1, nReactions
        if (nRealReactants(i) .EQ. 2) then
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

      write (*, FMT='(4(/, A32, I16))') &
        'Number of species: ', nSpecies, &
        'Number of reactions: ', nReactions, &
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

      CALL SaveMiscBefore (fU)

! Read initial conditions

      PopSpeciesInit = 0

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
          nameInitialSpecies_tmp, initialAbundance_tmp
        if (ios .NE. 0) then
          write (*,*) 'ios = ', ios; stop
        end if
        if ((nameInitialSpecies_tmp(1:1) .EQ. commentChar) .OR. &
            (nameInitialSpecies_tmp(1:1) .EQ. ' ')) cycle
        do j=1, nSpecies
          if (trim(nameSpecies(j)) .EQ. &
              trim(nameInitialSpecies_tmp)) then
            if (nameSpecies(j)(1:1) .EQ. 'g') then
              PopSpeciesInit(j) = int(initialAbundance_tmp, 8)
            else if (nameSpecies(j)(1:5) .EQ. 'Grain') then
              PopSpeciesInit(j) = int(initialAbundance_tmp, 8)
            else
              PopSpeciesInit(j) = int(VolumnPerGrain * n_H &
                * initialAbundance_tmp, 8)
            end if
            exit
          end if
        end do
      end do
      close (UNIT=fU, IOSTAT=ios, ERR=999, STATUS='KEEP')

      do i=1, nElement
        initialElementalAb(i) = &
          dot_product(dble(PopSpeciesInit(:)), SpeciesElements(i, :))
      end do

      CALL SYSTEM_CLOCK (ClockCount1, CountRate)

      write (*, '(/A, F10.3)') 'Seconds elapsed: ', &
          (ClockCount1 - ClockCount0) / real(CountRate)

! Now start the Monte Carlo

      nStepPerSave = (nMCSteps-1) / (nRecordsSave-1)
      nMCSteps = nStepPerSave * (nRecordsSave-1)

      write (*,'(/A)') 'Monte Carlo start...'
      write (*, '(4(/A32, I16))') &
        'Number of repeats: ', nMCRepeat, &
        'Number of steps of each repeat: ', nMCSteps, &
        'Number of steps per save: ', nStepPerSave, &
        'Number of records in total: ', nRecordsSave

      do i1=1, nMCRepeat

        write (*, '(/A16, I16, /)') 'Repeat: ', i1

        CALL init_random_seed()

        PopSpecies = PopSpeciesInit

        propensities = 0D0

        k = 1
        t_Save(1) = 0D0
        PopSpeciesSave(:, 1) = PopSpecies
        PopSpeciesAccum = 0D0
        t_Accum = 0D0
        IsReactionChanged = .TRUE.

        do i=1, nMCSteps
          TotalPropensity = 0D0
          do j=1, nReactions
            if (IsReactionChanged(j)) then
              select case (typeReac(j))
              case (64)
                propensities(j) = rates(j) &
                  * DBLE(PopSpecies(reactions(1,j))) &
                  * DBLE(PopSpecies(reactions(2,j)))
              case (63)
                propensities(j) = rates(j) &
                  * DBLE(PopSpecies(reactions(1,j))) &
                  * DBLE(PopSpecies(reactions(1,j)) - 1)
              case (61, 62)
                propensities(j) = rates(j) &
                  * DBLE(PopSpecies(reactions(1,j)))
              case (5, 53)
                propensities(j) = rates(j) &
                  * DBLE(PopSpecies(reactions(1,j))) &
                  * DBLE(PopSpecies(reactions(2,j)))
              case (1, 2, 3)
                propensities(j) = rates(j) &
                  * DBLE(PopSpecies(reactions(1,j)))
              end select
            end if
            TotalPropensity = TotalPropensity + propensities(j)
          end do

          if (TotalPropensity .EQ. 0D0) then
            write (*,*) 'No reaction can happen any more.'
            write (*,*) 'At step ', i
            k = k + 1
            t_Save(k:nRecordsSave) = t_Save(k-1) + t_Accum
            do j=k, nRecordsSave
              PopSpeciesSave(:, j) = DBLE(PopSpecies)
            end do
            exit
          end if

          CALL RANDOM_NUMBER(rand_tmp_vec)
          tNext = -DLOG(rand_tmp_vec(1)) / TotalPropensity
          RandProbTmp = TotalPropensity * rand_tmp_vec(2)

          PopSpeciesAccum = PopSpeciesAccum + DBLE(PopSpecies) * tNext
          t_Accum = t_Accum + tNext

          if (MOD(i, nStepPerSave) .EQ. 0) then
            CALL init_random_seed()
            write (*, '(A32, I16, " (", F6.2, "%)")') &
              CHAR(27)//'[ARunning: ', i, real(i)*100.0/nMCSteps
            k = k + 1
            t_Save(k) = t_Save(k-1) + t_Accum
            PopSpeciesSave(:, k) = PopSpeciesAccum / t_Accum
            PopSpeciesAccum = 0D0
            t_Accum = 0D0
          end if

          ! A pre-sort might make this faster.
          AccumPropensity = 0D0

          do iReactionNext=1, nReactions
            AccumPropensity = AccumPropensity &
              + propensities(iReactionNext)
            if (AccumPropensity .GT. RandProbTmp) exit
          end do

          nSpeUpdated = 0
          do j=1, nRealReactants(iReactionNext)
            PopSpecies(reactions(j, iReactionNext)) = &
              PopSpecies(reactions(j, iReactionNext)) - 1
            nSpeUpdated = nSpeUpdated + 1
            idxSpeUpdated(nSpeUpdated) = reactions(j, iReactionNext)
          end do
          do j=1+nReactants, nRealProducts(iReactionNext)+nReactants
            PopSpecies(reactions(j, iReactionNext)) = &
              PopSpecies(reactions(j, iReactionNext)) + 1
            nSpeUpdated = nSpeUpdated + 1
            idxSpeUpdated(nSpeUpdated) = reactions(j, iReactionNext)
          end do
          IsReactionChanged = .FALSE.
          do j=1, nSpeUpdated
            i2 = idxSpeUpdated(j)
            flag = .FALSE.
            do h=1, j-1
              if (idxSpeUpdated(h) .EQ. i2) then
                flag = .TRUE.
                exit
              end if
            end do
            if (flag) cycle
            do h=1, nReac_AsReactants(i2)
              IsReactionChanged(idxReac_AsReactants(h, i2)) = .TRUE.
            end do
          end do
        end do

        write (*, '(A16, I16)') 'Finish Repeat: ', i1

        CALL SYSTEM_CLOCK (ClockCount1, CountRate)
        write (*, '(/A, F10.3)') 'Seconds elapsed: ', &
            (ClockCount1 - ClockCount0) / real(CountRate)

        write (*, '(/A)') 'Saving results ...'

        write (strTMP, '(I6.6)') i1

        if (IsWordChar(fSaveALLHistory(1:1))) then
          open (UNIT=fU, FILE=trim(path)//trim(fSaveALLHistory)// &
            trim(adjustl(strTMP))//'.dat', IOSTAT=ios, &
            RECL=kind(1D0)*(nSpecies+1), STATUS='REPLACE', &
            ACCESS='DIRECT', FORM='UNFORMATTED', ACTION='WRITE')
          if (ios .LT. 0) then
            write (*,*) 'OPEN FILE ERROR: IOSTAT=', ios
            stop
          end if
          do i=1, nRecordsSave
            write (fU, rec=i)  t_Save(i), &
              PopSpeciesSave(:, i)
          end do
          close (UNIT=fU, IOSTAT=ios, ERR=999, STATUS='KEEP')
          open (UNIT=fU, FILE=trim(path)//trim(fSaveALLHistory)// &
            'FinalState_'//trim(adjustl(strTMP))//'.dat', IOSTAT=ios, &
            RECL=(kind(1D0) + 8*nSpecies), STATUS='REPLACE', &
            ACCESS='DIRECT', FORM='UNFORMATTED', ACTION='WRITE')
          if (ios .LT. 0) then
            write (*,*) 'OPEN FILE ERROR: IOSTAT=', ios
            stop
          end if
          write (fU, rec=1)  t_Save(nRecordsSave), PopSpecies
          close (UNIT=fU, IOSTAT=ios, ERR=999, STATUS='KEEP')
        end if

        if (IsWordChar(fSaveALLHistory_ASCII(1:1))) then
          CALL openFileSequentialWrite &
            (fU, trim(path)//trim(fSaveALLHistory_ASCII)// &
              trim(adjustl(strTMP))//'.dat', 99999)
          write (FMTstr, FMT='("(A12, 3X, ", I4, "(A", I2, ", 2X))")') &
            nSpecies, constLenNameSpecies
          write (UNIT=fU, FMT=FMTstr, IOSTAT=ios) &
            'Time_(yr)', nameSpecies(1:nSpecies)
          write (FMTstr, &
            FMT='("(ES12.3, 2X, ", I4, "(ES", I2, ".4E3, 2X))")') &
            nSpecies, constLenNameSpecies
          do i=1, nRecordsSave
            write (UNIT=fU, FMT=FMTstr) t_Save(i), &
              PopSpeciesSave(:, i)
          end do
          close (UNIT=fU, IOSTAT=ios, ERR=999, STATUS='KEEP')
        end if

        if (IsWordChar(fSaveElementalAb(1:1))) then
          do j=1, nElement
            finalElementalAb(j) = &
              dot_product(PopSpeciesSave(:, nRecordsSave), &
              SpeciesElements(j, :))
          end do
          CALL openFileSequentialWrite &
            (fU, trim(path)//trim(fSaveElementalAb)// &
              trim(adjustl(strTMP))//'.dat', 999)
          write (fU, '(A)') '! Elemental conservation:'
          write (fU, '(A12, A17, A17, A17)') &
              '!    Element', 'Initial', 'Final', '(F-I)/I'
          do j=1, nElement
            write (fU, '(A12, ES17.8E2, ES17.8E2, ES17.8E2)') &
                adjustr(nameElements(j)), &
                initialElementalAb(j), finalElementalAb(j), &
                (finalElementalAb(j) - initialElementalAb(j)) &
                / (abs(initialElementalAb(j)) + 1D-30)
          end do
          close (UNIT=fU, IOSTAT=ios, ERR=999, STATUS='KEEP')
        end if

        CALL SYSTEM_CLOCK (ClockCount1, CountRate)
        write (*, '(/A, F10.3)') 'Seconds elapsed: ', &
            (ClockCount1 - ClockCount0) / real(CountRate)
      end do

      write (*,'(/A)') 'Monte Carlo finish!'

998   CALL SYSTEM_CLOCK (ClockCount1, CountRate)
      write (*, '(/A, F10.3)') 'Seconds elapsed: ', &
          (ClockCount1 - ClockCount0) / real(CountRate)

999   write (*,'(A, A60, A)') &
        char(27)//'[47m'//char(27)//'[34m'//&
        'Program end', '', char(27)//'[0m'

      end program




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
