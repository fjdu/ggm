! Prototype finished on Oct 15th, 2009.
! Run fluently on Oct 23rd, 2009.
! Modification for RATE06, Nov 3rd, 2009
! A minor erros corrected.
!   2011-03-31 Thu 19:00:55
! <timestamp>2012-03-17 Sat 16:36:15</timestamp>
! Delete the analysis and plotting part.
! gfortran rate06_for_ini.f90 common.f90 opkd*.o -fbounds-check -o rate06_for_ini

      program rate06_new
      use CMDT
      implicit none

      external f, jac

      integer fU, ios, statALLOC, lenArgCMD
      integer nLineAll, nLineData, nLineHeader
      character commentChar, strArgCMD*128

      integer i, j, k, i1
      real, dimension(4) :: tmpVecReal
      integer, dimension(1) :: tmpVecInt
      integer, dimension(8) :: tmpVecInt8
      logical flag
      character(LEN=constLenNameSpecies) strTmp
      double precision dblTmp, dblTmp1

      integer ClockCount0, ClockCount1, CountRate

      double precision, dimension(:), allocatable :: &
          initialElementalAb, finalElementalAb
      integer, dimension(:), allocatable :: nElementReac, nElementProd

      integer  nIterationBase, nTargetDecades, nerr, nHistoryLen
      double precision ratioTout, tStep
      double precision, dimension(:), allocatable :: &
        initialAbundances, touts
      double precision, dimension(:, :), allocatable :: yHistory
      double precision, dimension(:, :), allocatable :: ydotHistory

      integer IOPT, ISTATE, ITASK, ITOL, LIW, LRW, NEQ, MF, NNZ
      double precision t, tout, tFinal
      double precision, dimension(:), allocatable :: y, RWORK
      integer, dimension(:), allocatable :: IWORK
      logical, dimension(:, :), allocatable :: sparseMaskJac

      integer idxSpeciesCheck, nNonzero, iCheck
      character(LEN=constLenNameSpecies) nameSpeciesCheck
      double precision rThreshold, tCheck
      integer, dimension(:), allocatable :: &
        idxReactionsNonzero, idxYdotSpeciesSorted, &
        idxYdotOneSpeciesSorted
      double precision, dimension(:), allocatable :: &
        ydotSpecies, ydotOneSpecies

      write (*,'(A)') &
        char(27)//'[47m'//char(27)//'[34m'//&
        'Program start'//char(27)//'[0m'

      CALL SYSTEM_CLOCK (ClockCount0, CountRate)

      fU = 1

      commentChar = '!'

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

      nLineHeader = nLineAll - nLineData
      nReactions = nLineData

      allocate &
        (strReactants(nReactants, nReactions), &
         strProducts(nProducts, nReactions), &
         nameSpecies(0:nSpecies_Est), dblABC(3, nReactions), &
         nRealReactants(nReactions), nRealProducts(nReactions), &
         reactions(nParticipants, nReactions), &
         iType(nReactions), &
         strType(nReactions), strQuality(nReactions), &
         dblTLTU(2, nReactions), rates(nReactions), &
         STAT=statALLOC)

      CALL ReadRate06_New (fU, nReactions, nLineHeader)

      close (UNIT=fU, IOSTAT=ios, ERR=999, STATUS='KEEP')

!      do i=1,nreactions
!      write (*, '(8A8,3E9.2,3I2)') &
!        strReactants(:,i), strProducts(:,i), dblABC(:,i), iType(i), &
!        nRealReactants(i), nRealProducts(i)
!      end do

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

      allocate &
        (SpeciesElements(nElement, nSpecies), &
         initialElementalAb(nElement), finalElementalAb(nElement), &
         nElementReac(nElement), nElementProd(nElement), &
         STAT=statALLOC)

      do i=1, nSpecies
        CALL getElements(nameSpecies(i), nameElements, nElement, &
          SpeciesElements(:, i))
      end do

      write (*, FMT=25) nSpecies, nReactions, &
        maxval(nRealReactants), maxval(nRealProducts)
25    FORMAT(/I8, ' species', /I8, ' reactions', &
             //'Max number of reactants: ', I4, &
             /'Max number of products:  ', I4)


! Calculate the rates.

      write (*,'(/A)') 'Calculating rates...'

      do i=1, nReactions
        if ((dblABC(3, i) .LT. -100.0) .AND. &
            (minval(dblTLTU(:,i)/Temperature) .GE. 10.0)) then
          rates(i) = 0D0
          write (*, '(2A, I6, 9A, F0.1)') &
            'Negative barrier reaction [discarded]: ', &
            char(27)//'[43m', i, char(27)//'[0m', &
            " ", strReactants(1:2,i), " -> ", strProducts(1:4,i), &
            dblABC(3, i)
          cycle
        endif

        select case (iType(i))
          case (5, 53)
            rates(i) = dblABC(1, i) * (T300**dblABC(2, i)) &
                * exp(-dblABC(3, i)/Temperature)
          case (1)
            rates(i) = dblABC(1, i)
          case (2)
            rates(i) = dblABC(1, i) * (T300**dblABC(2, i)) &
                * dblABC(3, i) / (1-omega_Albedo)
          case (3)
            rates(i) = dblABC(1, i) * exp(-dblABC(3, i) * Av)
          case (0)
            rates(i) = dblABC(1, i) * rateHHH2 / 2D0 &
                * sqrt(Temperature)
        end select

        do j=1, i-1
          if ((strType(i) .EQ. strType(j)) .AND. &
              (iType(i) .EQ. iType(j)) .AND. &
            (sum(abs(reactions(:, j)-reactions(:, i))) .EQ. 0)) then
            tmpVecReal(1:2) = abs(dblTLTU(:,j) - Temperature)
            tmpVecReal(3:4) = abs(dblTLTU(:,i) - Temperature)
            tmpVecInt = minloc(tmpVecReal)
            i1 = tmpVecInt(1)
            if ((i1 .EQ. 1) .OR. (i1 .EQ. 2)) then
              write (*,'(A, 2I6, 6X, 2A, I6, A)') &
                'Duplicate reactions: ', i, j, &
                'Use: ', char(27)//'[45m', j, char(27)//'[0m'
              rates(i) = 0D0
              exit
            end if
            if ((i1 .EQ. 3) .OR. (i1 .EQ. 4)) then
              write (*,'(A, 2I6, 6X, 2A, I6, A)') &
                'Duplicate reactions: ', i, j, &
                'Use: ', char(27)//'[45m', i, char(27)//'[0m'
              rates(j) = 0D0
              cycle
            end if
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
          rates(i) = 0D0
          write (*, '(2A, I6, A)') &
            'Elements not conserved [discarded]: ', &
            char(27)//'[41m', i, char(27)//'[0m'
          cycle
        end if
      end do

      deallocate(nElementReac, nElementProd, STAT=statALLOC)

! Write the retrieved reactions and rates into files.
! Saved rates are expressed in time unit of one second.
! Rates are saved in binary form to keep precision.
      write (*,'(/A)') 'Saving reactions...'
      CALL openFileSequentialWrite &
        (fU, trim(path)//trim(fReactionsSave), 999)
      do i=1, nReactions
        write (UNIT=fU, FMT='(9(I8))', IOSTAT=ios) &
          reactions(:, i), &
          nRealReactants(i), nRealProducts(i)
      end do
      close (UNIT=fU, IOSTAT=ios, ERR=999, STATUS='KEEP')

      write (*,'(/A)') 'Saving species...'
      CALL openFileSequentialWrite &
        (fU, trim(path)//trim(fSpeciesSave), 999)
      do i=1, nSpecies
        write (UNIT=fU, FMT='(A, 17I4)', IOSTAT=ios) &
          nameSpecies(i), SpeciesElements(:, i)
      end do
      close (UNIT=fU, IOSTAT=ios, ERR=999, STATUS='KEEP')

!      open (UNIT=fU, FILE=trim(path)//trim(fRatesSave), &
!           IOSTAT=ios, STATUS='REPLACE', ACCESS='DIRECT', &
!           RECL=KIND(1D0)*nReactions, FORM='UNFORMATTED', &
!           ACTION='WRITE')
!      if (ios .LT. 0) then
!          write (*,*) 'OPEN FILE ERROR: IOSTAT=', ios
!          STOP
!      end if
!      write (UNIT=fU, REC=1, IOSTAT=ios) rates
!      close (UNIT=fU, IOSTAT=ios, ERR=999, STATUS='KEEP')

      write (*,'(/A)') 'Saving rates...'
      CALL openFileSequentialWrite &
        (fU, trim(path)//trim(fRatesSave), 999)
      do i=1, nReactions
        write (UNIT=fU, FMT=89, &
          IOSTAT=ios) i, strReactants(1:2, i), strProducts(1:4,i), &
          rates(i), dblABC(:,i), iType(i), strQuality(i), strType(i)
      end do
89    FORMAT(I4, 2X, 6A, ES13.4E4, 4X, ES9.2E2, 2X, F8.2, &
        2X, F10.1, I4, X, A, X, A)
      close (UNIT=fU, IOSTAT=ios, ERR=999, STATUS='KEEP')

      CALL SYSTEM_CLOCK (ClockCount1, CountRate)
      write (*, '(/A, F10.3)') 'Seconds elapsed: ', &
          (ClockCount1 - ClockCount0) / real(CountRate)

! Change the time unit.
! Thus dt = 1 means a time step of 1 tiemUnit years.
      rates = rates * yearSeconds * timeUnit
      where ((rates .LE. rateThreshold) .OR. (iType .EQ. 53))
        rates = 0D0
        nRealReactants = 0
        nRealProducts = 0
      end where
      where (nRealReactants .EQ. 2)
        rates = rates * n_H
      end where

      NEQ = nSpecies

! Make sparse structure

      write (*,'(/A)') 'Making sparse structure...'

      allocate (sparseMaskJac(NEQ, NEQ), STAT=statALLOC)

      sparseMaskJac = .FALSE.

      do i=1, nReactions
        do j=1, nRealReactants(i)
          do k=1, nRealReactants(i)
            sparseMaskJac &
              (reactions(k, i), reactions(j, i)) = .TRUE.
          end do
          do k=1+nReactants, nRealProducts(i)+nReactants
            sparseMaskJac &
              (reactions(k, i), reactions(j, i)) = .TRUE.
          end do
        end do
      end do

      NNZ = COUNT (sparseMaskJac)

      write (*, FMT=830) NNZ, real(NNZ*100)/(NEQ*NEQ)
830   FORMAT (/'Number of nonzeros in Jac: ', I8, ' (', F0.2, '%)')

! Initialize variables of ODE solver.

      IOPT = 1 ! 1: allow optional input; 0: disallow
      LRW = 20 + NEQ * (12 + 1) + 3 * NEQ + & ! 500000
            4 * NNZ + 2 * NEQ + (NNZ + 10 * NEQ) / 1
      LIW = 31 + NEQ + NNZ
      MF = 021

      allocate (y(NEQ), RWORK(LRW), IWORK(LIW), STAT=statALLOC)

      y = 0D0
      RWORK(5:10) = 0D0
      IWORK(5:10) = 0
      IWORK(6) = 5000 ! Maximum number of steps

      IWORK(31) = 1

      k = 1
      do i=1, NEQ
        do j=1, NEQ
          if (sparseMaskJac(j, i)) then
            IWORK(31 + NEQ + k) = j
            k = k + 1
          end if
        end do
        IWORK(31+i) = k
      end do

      deallocate(sparseMaskJac, STAT=statALLOC)

      CALL SYSTEM_CLOCK (ClockCount1, CountRate)
      write (*, '(/A, F10.3)') 'Seconds elapsed: ', &
          (ClockCount1 - ClockCount0) / real(CountRate)

      write (*, '(/A)') 'Reading initial condition...'

! Read initial conditions

      CALL openFileSequentialRead &
          (fU, trim(path)//trim(fInitialCondition), 999)
      CALL GetNLineF (fU, nLineAll, nLineData, '')

      nInitialSpecies = nLineData
      allocate (nameInitialSpecies(nInitialSpecies), &
          initialAbundances(nInitialSpecies), &
          idxInitialSpecies(nInitialSpecies), STAT=statALLOC)

      write (*,'(/A, I8)') 'Number of initial species: ', &
          nInitialSpecies

      rewind (UNIT=fU, IOSTAT=ios)

      i1 = 0
      do i=1, nInitialSpecies
        read (UNIT=fU, FMT='(A12, F16.0)', IOSTAT=ios) &
          nameInitialSpecies(i), initialAbundances(i)

        if (ios .NE. 0) then
          write (*,*) 'ios = ', ios
          stop
        end if

        flag = .FALSE.
        do j=1, nSpecies
          if (trim(nameSpecies(j)) .EQ. &
            trim(nameInitialSpecies(i))) then
            idxInitialSpecies(i) = j
            y(j) = initialAbundances(i)
            if (SpeciesElements(3,j) .NE. 0) &
              y(j) = y(j) / ratioHDust
            i1 = i1 + 1
            flag = .TRUE.
            exit
          end if
        end do
        if (.NOT. flag) write (*, '(I4, 2X, 2A10)') &
          i, nameInitialSpecies(i), " not used!"
      end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      y = y / 2.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      write (*,'(/A, I8)') 'Actually used initial species: ', i1

      close (UNIT=fU, IOSTAT=ios, ERR=999, STATUS='KEEP')

      do i=1, nElement
        initialElementalAb(i) = &
          dot_product(y, SpeciesElements(i, :))
      end do

      CALL SYSTEM_CLOCK (ClockCount1, CountRate)
      write (*, '(/A, F10.3)') 'Seconds elapsed: ', &
          (ClockCount1 - ClockCount0) / real(CountRate)

      if (nIteration .LT. -1) then
844     write (*,FMT='(/A)') 'Number of iteration (please input):'
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

! Chemical evolution

      nHistoryLen = nIteration + 1
      allocate(yHistory(nSpecies, nHistoryLen), &
               ydotHistory(nSpecies, nHistoryLen), &
               touts(nHistoryLen), STAT=statALLOC)

      nIterationBase = 1000
      nTargetDecades = 9
      ratioTout = exp(real(nTargetDecades) / real(nIterationBase) &
                  * log(10D0))

      t = 0.0D0 ! Start time
      tStep = 1.0D-2 ! First time step
      tout = t + tStep ! First output time

      ITOL = 1 ! scalar control of tolerance
      ITASK = 1 ! normal, allow overshoot
      ISTATE = 1 ! first call

      yHistory(:, 1) = y
      touts(1) = t

      nerr = 0 ! for counting number of errors in iteration

      tFinal = tout + tStep*(ratioTout**nIteration-1D0)/(ratioTout-1D0)

      write (*, '(/A/)') 'Time span:'
      write (*, '(E10.2, A, ES10.2, A)') &
        t*timeUnit, ' to ', tFinal*timeUnit, ' years'

      write (*, '(/A)') 'Program timing reset to zero.'
      CALL SYSTEM_CLOCK (ClockCount0, CountRate)

      write (*, '(/A//)') 'Now start to evolve...'

      do i=1, nIteration
        write (*, FMT=850) &
          CHAR(27)//'[A', i, real(i*100)/nIteration, tStep
850     FORMAT (A, 'Iterating... ', I8, &
          ' (', F5.1, '%)', '  tStep = ', E10.2)

        call DLSODES &
             (f, NEQ, y, t, tout, &
              ITOL, RTOL, ATOL, &
              ITASK, ISTATE, IOPT, &
              RWORK, LRW, IWORK, LIW, &
              jac, MF)

        yHistory(:, i+1) = y
        touts(i+1) = t

        tStep = tStep * ratioTout
        tout = t + tStep

        if (ISTATE .LT. 0) then
          nerr = nerr + 1
          ISTATE = 3
        end if
      end do

      touts = touts * timeUnit

      write (*, FMT=860) nerr
860   FORMAT (/'Iteration completed!', &
              /'Number of errors in iteration: ', I8)

      CALL SYSTEM_CLOCK (ClockCount1, CountRate)
      write (*, '(/A, F10.3)') 'Seconds elapsed: ', &
          (ClockCount1 - ClockCount0) / real(CountRate)

      deallocate (RWORK, IWORK)

! Save evolve history in binary form.

      write (*, '(A)') 'Saving history...'
      open (UNIT=fU, FILE=trim(path)//trim(fSaveALLHistory), &
           IOSTAT=ios, RECL=kind(1D0)*(NEQ+1), STATUS='REPLACE', &
           ACCESS='DIRECT', FORM='UNFORMATTED', ACTION='WRITE')

      if (ios .LT. 0) then
          write (*,*) 'OPEN FILE ERROR: IOSTAT=', ios
          stop
      end if

      do i=1, nHistoryLen
        write (fU, rec=i) touts(i), yHistory(:, i)
      end do

      close (UNIT=fU, IOSTAT=ios, ERR=999, STATUS='KEEP')
      
      do i=1, nElement
        finalElementalAb(i) = &
          dot_product(y, SpeciesElements(i, :))
      end do

      CALL openFileSequentialWrite &
        (fU, trim(path)//trim(fSaveElementalAb), 999)
      write (fU, '(A)') '! Elemental conservation:'
      write (fU, '(A12, A12, A12, A12)') &
          '!    Element', 'Initial', 'Final', '(F-I)/I'
      do i=1, nElement
        write (fU, '(A12, ES12.3E2, ES12.3E2, ES12.3E2)') &
            adjustr(nameElements(i)), &
            initialElementalAb(i), finalElementalAb(i), &
            (finalElementalAb(i) - initialElementalAb(i)) &
            / (abs(initialElementalAb(i)) + 1D-30)
      end do

      close (UNIT=fU, IOSTAT=ios, ERR=999, STATUS='KEEP')

! Save the final abundance (the last y) and ydot into a file.

      CALL openFileSequentialWrite &
        (fU, trim(path)//trim(fSaveFinalResult), 999)

      ! 2011-05-16 Mon 19:12:56
      !     The output abundance is relative to n_H, not n(H2).
      ! yHistory = 2.0D0 * yHistory

      do i=1, nSpecies
        write (UNIT=fU, FMT=910, IOSTAT=ios) &
          nameSpecies(i), &
          yHistory(i, nHistoryLen), &
          yHistory(i, 1), &
          ydotHistory(i, nHistoryLen) &
            / (abs(yHistory(i, nHistoryLen))+1D-50), &
          (yHistory(i, nHistoryLen)-yHistory(i, 1)) &
            /(yHistory(i, 1)+1D-30)
910     FORMAT (A12, 2ES16.6E3, 4X, ES12.2E2, F8.2)
      end do

      close (UNIT=fU, IOSTAT=ios, ERR=999, STATUS='KEEP')


      CALL SYSTEM_CLOCK (ClockCount1, CountRate)
      write (*, '(/A, F10.3)') 'Seconds elapsed: ', &
          (ClockCount1 - ClockCount0) / real(CountRate)

999   write (*,'(A)') &
        char(27)//'[47m'//char(27)//'[34m'//&
        'Program end'//char(27)//'[0m'

      end program




      subroutine f (NEQ, t, y, ydot)
      use CMDT
      implicit none
      integer NEQ, i, j
      double precision t, y(NEQ), ydot(NEQ), rtmp

      ydot = 0.0D0

      do i=1, nReactions
        select case (iType(i))
          case (5)
            rtmp = rates(i) * y(reactions(1, i)) * y(reactions(2, i))
            ydot(reactions(1, i)) = ydot(reactions(1, i)) - rtmp
            ydot(reactions(2, i)) = ydot(reactions(2, i)) - rtmp
            do j=1+nReactants, nRealProducts(i)+nReactants
                ydot(reactions(j, i)) = &
                  ydot(reactions(j, i)) + rtmp
            end do
          case (1, 2, 3)
            rtmp = rates(i) * y(reactions(1, i))
            ydot(reactions(1, i)) = ydot(reactions(1, i)) - rtmp
            do j=1+nReactants, nRealProducts(i)+nReactants
              ydot(reactions(j, i)) = &
                ydot(reactions(j, i)) + rtmp
            end do
          case (0)
            rtmp = rates(i) * y(reactions(2, i))
            ydot(reactions(1, i)) = ydot(reactions(1, i)) - rtmp
            ydot(reactions(2, i)) = ydot(reactions(2, i)) - rtmp
            do j=1+nReactants, nRealProducts(i)+nReactants
              ydot(reactions(j, i)) = &
                ydot(reactions(j, i)) + rtmp
            end do
          case (53)
            if (nRealReactants(i) .EQ. 1) then
              rtmp = rates(i) * y(reactions(1, i))
              ydot(reactions(1, i)) = ydot(reactions(1, i)) - rtmp
              do j=1+nReactants, nRealProducts(i)+nReactants
                ydot(reactions(j, i)) = &
                  ydot(reactions(j, i)) + rtmp
              end do
            end if
            if (nRealReactants(i) .EQ. 2) then
              rtmp = rates(i) * y(reactions(1, i)) * y(reactions(2, i))
              ydot(reactions(1, i)) = ydot(reactions(1, i)) - rtmp
              ydot(reactions(2, i)) = ydot(reactions(2, i)) - rtmp
              do j=1+nReactants, nRealProducts(i)+nReactants
                  ydot(reactions(j, i)) = &
                    ydot(reactions(j, i)) + rtmp
              end do
            end if
        end select
      end do
      end subroutine f




      subroutine jac (NEQ, t, y, j, ian, jan, pdj)
      use CMDT
      double precision t, y, ian(*), jan(*), pdj, rtmp
      dimension y(NEQ), pdj(NEQ)
      integer NEQ, i, j, k

      do i=1, nReactions
        if ((j .NE. reactions(1, i)) .AND. &
            (j .NE. reactions(2, i))) cycle

        select case (iType(i))
          case (5)
            if (j .EQ. reactions(1, i)) &
              rtmp = rates(i) * y(reactions(2, i))
            if (j .EQ. reactions(2, i)) &
              rtmp = rates(i) * y(reactions(1, i))
            pdj(reactions(1, i)) = pdj(reactions(1, i)) - rtmp
            pdj(reactions(2, i)) = pdj(reactions(2, i)) - rtmp
            do k=1+nReactants, nRealProducts(i)+nReactants
              pdj(reactions(k, i)) = &
                pdj(reactions(k, i)) + rtmp
            end do
          case (1, 2, 3)
            pdj(reactions(1, i)) = pdj(reactions(1, i)) - rates(i)
            do k=1+nReactants, nRealProducts(i)+nReactants
              pdj(reactions(k, i)) = &
                pdj(reactions(k, i)) + rates(i)
            end do
          case (0)
            if (j .EQ. reactions(2, i)) then
              pdj(reactions(1, i)) = pdj(reactions(1, i)) - rates(i)
              pdj(reactions(2, i)) = pdj(reactions(2, i)) - rates(i)
              do k=1+nReactants, nRealProducts(i)+nReactants
                pdj(reactions(k, i)) = &
                  pdj(reactions(k, i)) + rates(i)
              end do
            end if
          case (53)
            if (nRealReactants(i) .EQ. 1) then
              pdj(reactions(1, i)) = pdj(reactions(1, i)) - rates(i)
              do k=1+nReactants, nRealProducts(i)+nReactants
                pdj(reactions(k, i)) = &
                  pdj(reactions(k, i)) + rates(i)
              end do
            end if
            if (nRealReactants(i) .EQ. 2) then
              if (j .EQ. reactions(1, i)) &
                rtmp = rates(i) * y(reactions(2, i))
              if (j .EQ. reactions(2, i)) &
                rtmp = rates(i) * y(reactions(1, i))
              pdj(reactions(1, i)) = pdj(reactions(1, i)) - rtmp
              pdj(reactions(2, i)) = pdj(reactions(2, i)) - rtmp
              do k=1+nReactants, nRealProducts(i)+nReactants
                pdj(reactions(k, i)) = &
                  pdj(reactions(k, i)) + rtmp
              end do
            end if
        end select
      end do
      end subroutine jac




      subroutine fOneSpecies (NEQ, y, idxSpecies, rThreshold, &
        ydot, idxReactionsNonzero, nNonzero)
      use CMDT
      implicit none
      integer i, j, idx1, idx2, idxSpecies, NEQ, nNonzero
      double precision rThreshold, y(NEQ), ydot(nReactions)
      integer idxReactionsNonzero(nReactions)
      logical flag

      ydot = 0.0D0
      nNonzero = 0

      do i=1, nReactions
        flag = .FALSE.
        do j=1, nRealReactants(i) !2
          if (reactions(j,i) .EQ. idxSpecies) then
            if (.NOT. flag) then
              nNonzero = nNonzero + 1
              idxReactionsNonzero(nNonzero) = i
              flag = .TRUE.
            end if
            select case (iType(i))
            case (5)
              ydot(nNonzero) = ydot(nNonzero) &
                -rates(i) * y(reactions(1, i)) * y(reactions(2, i))
            case (1, 2, 3)
              ydot(nNonzero) = ydot(nNonzero) &
                -rates(i) * y(reactions(1, i))
            case (0)
              ydot(nNonzero) = ydot(nNonzero) &
                -rates(i) * y(reactions(2, i))
            case (53)
              if (nRealReactants(i) .EQ. 1) then
                ydot(nNonzero) = ydot(nNonzero) &
                  -rates(i) * y(reactions(1, i))
              end if
              if (nRealReactants(i) .EQ. 2) then
                ydot(nNonzero) = ydot(nNonzero) &
                  -rates(i) * y(reactions(1, i)) * y(reactions(2, i))
              end if
            end select
!             if (abs(ydot(nNonzero)) .LT. rThreshold) then
!               ydot(nNonzero) = 0.0D0
!               nNonzero = nNonzero - 1
!             end if
          end if
        end do
        do j=1, nRealProducts(i) !4
          if (reactions(j+nReactants,i) .EQ. idxSpecies) then
            if (.NOT. flag) then
              nNonzero = nNonzero + 1
              idxReactionsNonzero(nNonzero) = i
              flag = .TRUE.
            end if
            select case(iType(i))
            case (5)
              ydot(nNonzero) = ydot(nNonzero) + &
                rates(i) * y(reactions(1, i)) * y(reactions(2, i))
            case (1, 2, 3)
              ydot(nNonzero) = ydot(nNonzero) + &
                rates(i) * y(reactions(1, i))
            case (0)
              ydot(nNonzero) = ydot(nNonzero) + &
                rates(i) * y(reactions(2, i))
            case (53)
              if (nRealReactants(i) .EQ. 1) then
                ydot(nNonzero) = ydot(nNonzero) + &
                  rates(i) * y(reactions(1, i))
              end if
              if (nRealReactants(i) .EQ. 2) then
                ydot(nNonzero) = ydot(nNonzero) + &
                  rates(i) * y(reactions(1, i)) * y(reactions(2, i))
              end if
            end select
!             if (abs(ydot(nNonzero)) .LT. rThreshold) then
!               ydot(nNonzero) = 0.0D0
!               nNonzero = nNonzero - 1
!             end if
          end if
        end do
      end do
      end subroutine fOneSpecies




      subroutine GetNLineF (fU, nLineAll, nLineData, commentChar)
      integer fU, ios, nLineAll, nLineData
      character commentChar, chartmp, strtmp*16

      rewind (UNIT=fU, IOSTAT=ios)

      nLineAll = 0
      nLineData = 0
      do
        read (UNIT=fU, FMT='(A16)', IOSTAT=ios) strtmp
        if (ios .LT. 0) exit
        nLineAll = nLineAll + 1
        if (strtmp(1:1) .EQ. ' ') strtmp = adjustl(strtmp)
        chartmp = strtmp(1:1)
        if (chartmp .NE. commentChar) nLineData = nLineData + 1
      end do
      end subroutine GetNLineF




      subroutine ReadRate06_New (fU, nLineData, nLineHeader)
      USE CMDT
      implicit none
      integer i, j, k, fU, ios, nLineData, nLineHeader
      character strtmp

      nRealReactants = 0
      nRealProducts = 0

      rewind (UNIT=fU, IOSTAT=ios)

      do i=1, nLineHeader
        read (UNIT=fU, FMT='(A1)', IOSTAT=ios) strtmp
      end do

      do i=1, nLineData
        read (UNIT=fU, FMT=&
          '(7(A12), 3F9.0, 2F6.0, I3, X, A1, X, A2)', &
          IOSTAT=ios) &
          strReactants(:,i), strProducts(:,i), dblABC(:,i), &
          dblTLTU(:, i), iType(i), strQuality(i), strType(i)
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
      end do

      do i=1, nLineData
        do j=1, nRealReactants(i)
          if (trim(strReactants(j, i)) .EQ. 'PHOTON') then
            nRealReactants(i) = nRealReactants(i) - 1
            exit
          end if
          if (trim(strReactants(j, i)) .EQ. 'CRPHOT') then
            nRealReactants(i) = nRealReactants(i) - 1
            exit
          end if
          if (trim(strReactants(j, i)) .EQ. 'CRP') then
            nRealReactants(i) = nRealReactants(i) - 1
            exit
          end if
        end do
        do j=1, nRealProducts(i)
          if (trim(strProducts(j, i)) .EQ. 'PHOTON') then
            nRealProducts(i) = nRealProducts(i) - 1
            exit
          end if
        end do
      end do

      end subroutine ReadRate06_New




      subroutine openFileSequentialRead (fU, filename, maxRecLen)
      integer fU, maxRecLen, ios
      character(len=*) filename

      open (UNIT=fU, FILE=trim(filename), IOSTAT=ios, &
           STATUS='OLD', ACCESS='SEQUENTIAL', FORM='FORMATTED', &
           RECL=maxRecLen, BLANK='NULL', POSITION='REWIND', &
           ACTION='READ', DELIM='NONE', PAD='YES')

      if (ios .NE. 0) then
        write (*, '(/A, I8, /A, A)') 'OPEN FILE ERROR: IOSTAT=', ios, &
          'Filename: ', filename
        stop
      end if
      end subroutine openFileSequentialRead



      subroutine openFileSequentialWrite (fU, filename, maxRecLen)
      integer fU, maxRecLen, ios
      character(len=*) filename

      open (UNIT=fU, FILE=trim(filename), IOSTAT=ios, &
           STATUS='REPLACE', ACCESS='SEQUENTIAL', FORM='FORMATTED', &
           RECL=maxRecLen, BLANK='NULL', POSITION='REWIND', &
           ACTION='WRITE', DELIM='NONE', PAD='YES')
      if (ios .NE. 0) then
        write (*, '(/A, I8, /A, A)') 'OPEN FILE ERROR: IOSTAT=', ios, &
          'Filename: ', filename
        stop
      end if
      end subroutine openFileSequentialWrite




      subroutine getElements &
        (nameSpec, listElements, nElements, arrNElements)
      character(len=*) nameSpec, listElements(nElements)
      integer, dimension(nElements) :: arrNElements
      integer i, j, k, ntmp, lenName, lenEle
      integer, dimension(32) :: belongto
      logical, dimension(32) :: used
      logical flagReplace
      integer, parameter :: chargePos = 1

      arrNElements = 0

      lenName = len(trim(nameSpec))
      belongto = 0
      used = .FALSE.

      do i=1, nElements
        lenEle = len(trim(listElements(i)))
        do j=1, lenName-lenEle+1
          if (nameSpec(j:(j+lenEle-1)) .EQ. &
              listElements(i)(1:(lenEle))) then
            flagReplace = .TRUE.
            do k=j, (j+lenEle-1)
              if (used(k)) then
                if (len(trim(listElements(belongto(k)))) .GE. &
                    len(trim(listElements(i)))) then
                  flagReplace = .FALSE.
                  exit
                else
                  arrNElements(belongto(k)) = &
                    arrNElements(belongto(k)) - 1
                end if
              end if
            end do
            if (flagReplace) then
              belongto(j:(j+lenEle-1)) = i
              used(j:(j+lenEle-1)) = .TRUE.
              arrNElements(i) = arrNElements(i) + 1
            end if
          end if
        end do
      end do

      do i=1, lenName
        if (.NOT. used(i)) then
          do j=1, (i-1)
            if (used(i-j)) then
              belongto(i) = belongto(i-j)
              exit
            end if
          end do
          if (((nameSpec(i-1:i-1) .GT. '9') .OR. &
            (nameSpec(i-1:i-1) .LT. '0')) .AND. &
            (nameSpec(i:i) .LE. '9') .AND. &
            (nameSpec(i:i) .GE. '0')) then
            if ((nameSpec(i+1:i+1) .LE. '9') .AND. &
              (nameSpec(i+1:i+1) .GE. '0')) then
              read (nameSpec(i:i+1), '(I2)') ntmp
              if (ntmp .EQ. 0) cycle
              arrNElements(belongto(i)) = &
                arrNElements(belongto(i)) + ntmp - 1
            else
              read (nameSpec(i:i), '(I1)') ntmp
              if (ntmp .EQ. 0) cycle
              arrNElements(belongto(i)) = &
                arrNElements(belongto(i)) + ntmp - 1
            end if
          else if (nameSpec(i:i) .EQ. '+') then
            arrNElements(chargePos) = 1
          else if (nameSpec(i:i) .EQ. '-') then
            arrNElements(chargePos) = -1
          end if
        end if
      end do
      end subroutine getElements




      subroutine initialize (fU, fileInitialize)
      use CMDT
      character(LEN=128) fileInitialize
      integer fU, ios

      namelist /PhysicalParameters/ &
        Temperature, n_H, &
        Av, omega_Albedo, ratioHDust, &
        rateHHH2, rateCosIon, &
        ratioDustGasMass, stickCoeffHH, stickCoeffChargeGrain, &
        GrainDensity, GrainRadius, GrainRadiusBase, &
        aGrainMin, aGrainMax, &
        rateThreshold, timeUnit, nSpecies_Est

      namelist /ODEParameters/ &
        ATOL, RTOL, nIteration
        
      namelist /PlotParameters/ &
        xPGmin, xPGmax, yPGmin, yPGmax

      namelist /Paths/ path, &
        fReactions, fDissRecom, &
        fInitialCondition, fPlotList, &
        fReactionsSave, fSpeciesSave, fRatesSave, &
        fSaveFinalResult, fSaveElementalAb, fSaveALLHistory, &
        fPGout, fAnalyseResult

      CALL openFileSequentialRead (fU, fileInitialize, 999)
      read (UNIT=fU, IOSTAT=ios, NML=PhysicalParameters)
      read (UNIT=fU, IOSTAT=ios, NML=ODEParameters)
      read (UNIT=fU, IOSTAT=ios, NML=PlotParameters)
      read (UNIT=fU, IOSTAT=ios, NML=Paths)
      close (UNIT=fU, IOSTAT=ios, STATUS='KEEP')

      T300 = Temperature / 300.0D0
      TemperatureReduced = kBoltzmann * Temperature / &
        (elementaryCharge**2 * coulombConst / (GrainRadius*1D-6))
      JNegaPosi = (1D0 + 1D0/TemperatureReduced) * &
                  (1D0 + SQRT(2D0/(2D0+TemperatureReduced)))
      JChargeNeut = (1D0 + SQRT(DPI/2D0/TemperatureReduced))

      if (ratioHDust .LT. 1D-80) &
        ratioHDust = 1D0 / & ! N_H / N_Grain
          (ratioDustGasMass / (1D0-0.25D0) * mProton / & 
          (4D0*DPI/3D0 * (GrainRadius*1D-4)**3 * GrainDensity))
      if (rateHHH2 .LT. 1D-80) &
        rateHHH2 = 0.5D0 * stickCoeffHH &
          * 3D0 / 4D0 * 1.4D0 * mProton &
          * ratioDustGasMass / GrainDensity &
          / sqrt(aGrainMin * aGrainMax) &
          * sqrt(8D0/DPI*kBoltzmann/(mProton*1D-3)) * 1D2
      end subroutine initialize




      subroutine SORT_Asc_idx (nDim, Y, idxSorted)
      implicit none
      integer nDim, i, j, itmp
      integer idxSorted(nDim)
      double precision Y(nDim)

      do i=1, nDim
        idxSorted(i) = i
      end do

      do i=1, nDim
        do j=i+1, nDim
          if (Y(idxSorted(j)) .LT. Y(idxSorted(i))) then
            itmp = idxSorted(i)
            idxSorted(i) = idxSorted(j)
            idxSorted(j) = itmp
          end if
        end do
      end do

      end subroutine SORT_Asc_idx
