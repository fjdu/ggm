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
! <timestamp>2011-07-01 Fri 20:10:59</timestamp>
!
! gfortran common_MC_MO.o subroutines_MC.o MC.o common_MC_MO_TP.o subroutines_share_triv.o subroutines_share_proc.o -fbounds-check -o MC

      program MC
      use CMDT
      use CMTP
      implicit none

      integer fU, fUyHistory_ascii, fUyHistory_bin
      integer ios, statALLOC, lenArgCMD
      integer nLineAll, nLineData
      character commentChar, strArgCMD*128, strTMP*128, FMTStr*128
      character(LEN=constLenNameSpecies) nameInitialSpecies_tmp
      logical flag

      integer(kind=LongIntKind) i, j, k, h, i1, i2, iReactionNext
      double precision RandProbTmp, dblTmp
      double precision, dimension(2) :: rand_tmp_vec
      double precision TotalPropensity, AccumPropensity
      double precision initialAbundance_tmp
      integer, dimension(8) :: intTmpvec

      real ProgStartTime, ProgCurrentTime

      double precision, dimension(:), allocatable :: &
        initialElementalAb, finalElementalAb
      integer, dimension(:), allocatable :: &
        nElementReac, nElementProd

      double precision tNext, t_accum
      double precision, dimension(:), allocatable :: PopSpeciesAccum
      integer(kind=LongIntKind), dimension(:), allocatable :: &
        PopSpeciesInit
      double precision, dimension(:), allocatable :: steps_preset
      integer(kind=LongIntKind), parameter :: &
        nStepsLimit = HUGE(0_LongIntKind)
      double precision t_save
      double precision, dimension(:), allocatable :: PopSpeciesSave

      logical IsWordChar, getPhyPar, ExternalAction, getFileUnit, &
        FileUnitOpened

      write (*,'(A, A60, A)') &
        char(27)//'[47m'//char(27)//'[34m'//&
        'Program start', '', char(27)//'[0m'

      CALL CPU_TIME(ProgStartTime)
      call date_and_time(DATE=strTMP, TIME=FMTStr)
      write (*,'("Current time: ", A8, 2X, A10)') strTMP, FMTStr

      if (.NOT. getFileUnit(fU)) then
        write (*,*) 'Cannot allocate an output file unit!'
        stop
      end if

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
         rates(nReactions), propensities(nReactions), &
         BranchingRatios(nReactions), STAT=statALLOC)

      CALL ReadReactions (fU, nLineAll, commentChar)
      close (UNIT=fU, IOSTAT=ios, ERR=999, STATUS='KEEP')

      CALL MakeReactions

      allocate &
        (SpeciesElements(nElement, nSpecies), &
         massSpecies(nSpecies), mobilitySpecies(nSpecies), &
         initialElementalAb(nElement), finalElementalAb(nElement), &
         nElementReac(nElement), nElementProd(nElement), &
         PopSpecies(nSpecies), PopSpeciesInit(nSpecies), &
         PopSpeciesSave(nSpecies), &
         IsReactionChanged(nReactions), idxSpeUpdated(8), &
         idxReac_AsReactants(nReactions, nSpecies), &
         idxReac_AsProducts(nReactions, nSpecies), &
         nReac_AsReactants(nSpecies), nReac_AsProducts(nSpecies), &
         PopSpeciesAccum(nSpecies), &
         idxGrSpecies(nSpecies), &
         idxGrReactions(nReactions), &
         IsGas(nSpecies), &
         enthalpySpecies(nSpecies), &
         E_desorp_species(nSpecies), &
         vibFreqSpecies(nSpecies), &
         STAT=statALLOC)

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
            write (*,'(2A, 2I6, A)') &
              'Duplicate reaction pair: ', &
              char(27)//'[45m', i, j, char(27)//'[0m'
          end if
        end do
      end do
      deallocate (nElementReac, nElementProd, STAT=statALLOC)

      write (*, FMT='(4(/, A32, I16))') &
        'Number of species: ', nSpecies, &
        'Number of grain species: ', nGrSpecies, &
        'Number of reactions: ', nReactions, &
        'Number of surface reactions: ', nGrReactions, &
        'Max number of reactants: ', maxval(nRealReactants), &
        'Max number of products: ', maxval(nRealProducts)

      write (*, '(/A)') 'Calculating rates...'

      ! These three do not depend on the physical conditions.
      CALL getDesorEnergy ! From the reaction file.
      CALL getVibFreq
      CALL ImportSpeciesEnthalpy
      CALL MakeBranchInfo

      ! These rates does depend on the physical conditions.
      CALL CalcRates

      rates = rates * SecondsPerYear

      write (*, '(3(/A32, ES16.6))') &
        'Grain to H number ratio: ', ratioGrainToH, &
        'Sites per grain: ', SitesPerGrain, &
        'Volumn containing one grain: ', VolumnPerGrain

      CALL SaveMiscBefore

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
        read (strTMP, FMT='(A12, ES16.6)', IOSTAT=ios) &
          nameInitialSpecies_tmp, initialAbundance_tmp
        if (ios .NE. 0) then
          write (*,*) 'ios = ', ios; stop
        end if
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
          dot_product(DBLE(PopSpeciesInit(:)), SpeciesElements(i, :))
      end do

      CALL CPU_TIME(ProgCurrentTime)
      write (*, '(/A, F10.3/)') 'Seconds elapsed: ', &
        ProgCurrentTime - ProgStartTime

      if (.NOT. IsWordChar(fSaveALLHistory(1:1))) then
        write (*,*) 'fSaveALLHistory is not provided!'
        stop
      end if
      if (.NOT. IsWordChar(fSaveALLHistory_ASCII(1:1))) then
        write (*,*) 'fSaveALLHistory_ASCII is not provided!'
        stop
      end if

! Now start the Monte Carlo

      allocate(steps_preset(nIteration), STAT=statALLOC)

      steps_preset(1) = &
        tFinal * eta_step / (exp(log(1+eta_step)*(nIteration-1)) - 1)
      dblTmp = 1D0 + eta_step
      do i=2, nIteration
        steps_preset(i) = steps_preset(i-1) * dblTmp
      end do

      write (*,'(/A)') 'Monte Carlo start...'
      write (*,*) 'Number of repeats: ', nMCRepeat
      write (*,*) 'Initial time step for smooth = ', steps_preset(1)
      write (*,*) 'Final time to be reached = ', sum(steps_preset)

      do i1=1, nMCRepeat

        write (*, '(/A16, I16, /)') 'Repeat: ', i1

        write (strTMP, '(I6.6)') i1

        if (.NOT. getFileUnit(fUyHistory_bin)) then
          write (*,*) 'Cannot allocate an output file unit!'
          stop
        end if
        open (UNIT=fUyHistory_bin, &
          FILE=trim(path)//trim(fSaveALLHistory)// &
            trim(adjustl(strTMP))//'.dat', &
          IOSTAT=ios, RECL=kind(1D0)*(nSpecies+1), &
          STATUS='REPLACE', ACCESS='DIRECT', &
          FORM='UNFORMATTED', ACTION='WRITE')
        if (ios .LT. 0) then
          write (*, *) 'OPEN FILE ERROR: IOSTAT=', ios
          stop
        end if
        if (.NOT. getFileUnit(fUyHistory_ascii)) then
          write (*,*) 'Cannot allocate an output file unit!'
          stop
        end if
        CALL openFileSequentialWrite (fUyHistory_ascii, &
          trim(path)//trim(fSaveALLHistory_ASCII)// &
            trim(adjustl(strTMP))//'.dat', &
          kind('a')*(nSpecies+1)*constLenNameSpecies*2)
        write (FMTStr, &
          FMT='("(A12, 3X, ", I4, "(A", I2, ", 2X))")') &
          nSpecies, constLenNameSpecies
        write (UNIT=fUyHistory_ascii, FMT=FMTStr, &
          IOSTAT=ios) 'Time_(yr)', nameSpecies(1:nSpecies)
        write (FMTStr, &
          FMT='("(ES12.3, 2X, ", I4, "(ES", I2, ".4E3, 2X))")') &
          nSpecies, constLenNameSpecies
        write (UNIT=fUyHistory_ascii, FMT=FMTStr, &
          IOSTAT=ios) 0D0, DBLE(PopSpeciesInit)

        CALL init_random_seed()

        PopSpecies = PopSpeciesInit

        propensities = 0D0

        PopSpeciesAccum = 0D0
        t_accum = 0D0
        IsReactionChanged = .TRUE.

        k = 1
        do i=1, nStepsLimit
          TotalPropensity = 0D0
          do j=1, nReactions
            if (IsReactionChanged(j)) then
              if (nRealReactants(j) .EQ. 1) then
                propensities(j) = rates(j) &
                  * DBLE(PopSpecies(reactions(1,j)))
              else
                if (reactions(1,j) .NE. reactions(2,j)) then
                  propensities(j) = rates(j) &
                    * DBLE(PopSpecies(reactions(1,j))) &
                    * DBLE(PopSpecies(reactions(2,j)))
                else
                  propensities(j) = rates(j) &
                    * DBLE(PopSpecies(reactions(1,j))) &
                    * DBLE(PopSpecies(reactions(1,j)) - 1)
                end if
              end if
            end if
            TotalPropensity = TotalPropensity + propensities(j)
          end do

          if (TotalPropensity .EQ. 0D0) then
            write (*,*) 'No reaction can happen any more.'
            write (*,*) 'At step ', i
            t_save = t_save + t_accum
            PopSpeciesSave = PopSpeciesAccum / t_accum
            write(fUyHistory_bin, rec=k) t_save, PopSpeciesSave
            write(fUyHistory_ascii, FMT=FMTStr) t_save, PopSpeciesSave
            exit
          end if

          CALL RANDOM_NUMBER(rand_tmp_vec)
          tNext = -DLOG(rand_tmp_vec(1)) / TotalPropensity
          RandProbTmp = TotalPropensity * rand_tmp_vec(2)

          PopSpeciesAccum = PopSpeciesAccum + DBLE(PopSpecies) * tNext
          t_accum = t_accum + tNext

          if (t_accum .GE. steps_preset(k)) then
            CALL init_random_seed()
            t_save = t_save + t_accum
            PopSpeciesSave = PopSpeciesAccum / t_accum
            write (*, '(A8, "Step/Record/Time/TimeStep: ", &
              & I16, " / ", I8, 2(" / ", ES12.3))') &
              CHAR(27)//'[A', i, k, t_save, t_accum
            write(fUyHistory_bin, rec=k) t_save, PopSpeciesSave
            write(fUyHistory_ascii, FMT=FMTStr) t_save, PopSpeciesSave
            FLUSH (fUyHistory_bin)
            FLUSH (fUyHistory_ascii)
            if (t_save .GT. tFinal) then
              write (*,*) 'Reach tFinal at step ', i
              write (*,*) 'Total records saved: ', k
              exit
            end if
            PopSpeciesAccum = 0D0
            t_accum = 0D0
            k = k + 1
          end if

          AccumPropensity = 0D0

          do iReactionNext=1, nReactions
            AccumPropensity = AccumPropensity &
              + propensities(iReactionNext)
            if (AccumPropensity .GE. RandProbTmp) exit
          end do
          if (iReactionNext .GT. nReactions) then
            write (*,*) 'Warning: iReactionNext > nReactions!'
            write (*,*) AccumPropensity, RandProbTmp
            write (*,*) '  ... continue anyway.'
            write (*,*) '  Set iReactionNext = nReactions'
            write (*,*)
            iReactionNext = nReactions
          end if

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

        if (FileUnitOpened(fUyHistory_bin)) &
          close (UNIT=fUyHistory_bin, IOSTAT=ios, STATUS='KEEP')
        if (FileUnitOpened(fUyHistory_ascii)) &
          close (UNIT=fUyHistory_ascii, IOSTAT=ios, STATUS='KEEP')

        if (IsWordChar(fSaveElementalAb(1:1))) then
          do j=1, nElement
            finalElementalAb(j) = &
              dot_product(PopSpeciesSave, SpeciesElements(j, :))
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

        CALL CPU_TIME(ProgCurrentTime)
        write (*, '(/A, F10.3/)') 'Seconds elapsed: ', &
          ProgCurrentTime - ProgStartTime
      end do

      write (*,'(/A)') 'Monte Carlo finish!'

999   write (*,'(A, A60, A)') &
        char(27)//'[47m'//char(27)//'[34m'//&
        'Program end', '', char(27)//'[0m'

      end program
