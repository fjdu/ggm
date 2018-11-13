      subroutine f (NEQ, t, y, ydot)
      use CMDT
      use CMTP
      implicit none
      integer NEQ, i, j, i1
      double precision t, y(NEQ), ydot(NEQ)

      ydot = 0D0

      do i=1, nSpecies
        do j=1, ndtMM(i)
          ydot(i) = ydot(i) + dtMMMo(i)%coeff(j) * &
            getMomValue(dtMM(1, j, i), i, y, NEQ)
        end do
      end do
      do i=nSpecies+1, NEQ
        i1 = indY(i)
        do j=1, ndtMM(i1)
          ydot(i) = ydot(i) + dtMMMo(i1)%coeff(j) * &
            getMomValue(dtMM(1, j, i1), i, y, NEQ)
        end do
      end do
      end subroutine f




      subroutine jac (NEQ, t, y, jj, ian, jan, pdj)
!     Return: \partial Ydot / \partial Y(jj)
      use CMDT
      use CMTP
      implicit none
      double precision t, y, ian(*), jan(*), pdj
      dimension y(NEQ), pdj(NEQ)
      integer NEQ, i, j, jj, i1, j1
      end subroutine jac




      subroutine makeSparse (y, NEQ)
      use CMDT
      use CMTP
      implicit none
      integer i, j, NEQ
      double precision, dimension(NEQ) :: y
      sparseMaskJac = .FALSE.
      do i=1, NEQ
        do j=1, NEQ
          sparseMaskJac(i, j) = getSparse(indY(i), indY(j), y, NEQ)
        end do
      end do
      end subroutine makeSparse



      subroutine printjac (fU, NEQ, t, y)
!     Return: \partial Ydot / \partial Y(jj)
      use CMDT
      implicit none
      integer fU, NEQ, jj
      double precision t
      double precision, dimension(NEQ) :: y, pdj
      double precision, dimension(1) :: ian, jan
      character(LEN=128) FMTstr
      write (FMTstr, '("(", I5, "A12)")') NEQ
      write (UNIT=fU, FMT=FMTstr) &
        nameSpecies(1:nSpecies), nameMoments
      write (FMTstr, '("(", I5, "ES12.2E3)")') NEQ
      do jj = 1, NEQ
        pdj = 0D0
        CALL jac(NEQ, t, y, jj, ian, jan, pdj)
        write (UNIT=fU, FMT=FMTstr) pdj
      end do
      end subroutine printjac




      subroutine fOneSpecies (NEQ, y, iSpecies, ydot, &
        idxReactionsNonzero, nNonzero)
      use CMDT
      use CMTP
      implicit none
      integer NEQ, nNonzero, iSpecies, j
      double precision t, y(NEQ), ydot(nReactions)
      integer idxReactionsNonzero(nReactions)
      nNonzero = ndtMM(iSpecies)
      do j=1, ndtMM(iSpecies)
        ydot(j) = dtMMMo(iSpecies)%coeff(j) * &
          getMMProd(dtMM(1, j, iSpecies), y, NEQ)
        idxReactionsNonzero(j) = dtMM(2, j, iSpecies)
      end do
      end subroutine fOneSpecies




      subroutine AnalyseReaction
      use CMDT
      use CMTP
      implicit none
      integer fU, i, i1, j, j1, k, iCheck, nNonzero
      integer, dimension(8) :: tmpVecInt8
      double precision, dimension(:), allocatable :: ydotOneSpecies
      integer, dimension(:), allocatable :: idxReactionsNonzero
      integer, dimension(:), allocatable :: idxYdotOneReactionSorted
      logical flagTmp, getFileUnit, getPhyPar
      double precision dblTmp, dblTmp1

      allocate (ydotOneSpecies(nReactions), &
        idxReactionsNonzero(nReactions), &
        idxYdotOneReactionSorted(nReactions))

      if (.NOT. getFileUnit(fU)) then
        write (*,*) 'Cannot allocate an output file unit!'
        stop
      end if
      CALL openFileSequentialWrite &
        (fU, trim(path)//trim(fAnalyseResult), 999)

      write (fU, '("! ", A/)') '----Analysis----'
      write (fU, '(1("! ", A/))') &
        'For each species find out at what time it changes fastest'

      do i=1, nSpecies
        write (fU, '(A, "  #", I0)') nameSpecies(i), i
        tmpVecInt8(1) = 1
        tmpVecInt8(2) = nHistoryLen
        CALL SORT_Asc_idx(nHistoryLen, ydotHistory(i, :), &
          idxYdotOneReactionSorted(1:nHistoryLen))
        tmpVecInt8(3) = idxYdotOneReactionSorted(nHistoryLen)
        tmpVecInt8(4) = idxYdotOneReactionSorted(1)
        CALL SORT_Asc_idx(nHistoryLen, &
          ydotHistory(i, :) * touts / (yHistory(i, :) + 1D-30), &
          idxYdotOneReactionSorted(1:nHistoryLen))
        tmpVecInt8(5) = idxYdotOneReactionSorted(nHistoryLen)
        tmpVecInt8(6) = idxYdotOneReactionSorted(1)
        do j=1, 6
          iCheck = tmpVecInt8(j)
          !++++++++++++++
          ! 2011-01-20 Thu 23:58:54
          ! In analysing the rates have to be re-calculated.
          if (getPhyPar(touts(iCheck))) then
            dblTmp = 0D0 ! dumb statement
          end if
          !-------------
          CALL CalcMobilities
          CALL CalcRates
          rates = rates * SecondsPerYear
          do i1=1, nNumMM
            if (ndtMM(i1) .GT. 0) then
              do j1=1, ndtMM(i1)
                dtMMMo(i1)%coeff(j1) = rates(dtMM(2, j1, i1)) * &
                  DBLE(dtMM(3, j1, i1))
              end do
            end if
          end do
          CALL fOneSpecies(nSpecies, yHistory(:, iCheck), i, &
            ydotOneSpecies, idxReactionsNonzero, nNonzero)
          CALL SORT_Asc_idx(nNonzero, &
            -abs(ydotOneSpecies(1:nNonzero)), idxYdotOneReactionSorted)
          dblTmp = ydotHistory(i, iCheck)
          dblTmp1 = 0D0
          write (fU, '(2X, "[", I0, "] ", "At time: ", ES10.2, 2X, &
            "Abundance: ", ES10.2, 2X, "Rate: ", &
            ES10.2, 2X, "Number of Effective Reactions: ", I4)') &
            j, touts(iCheck), yHistory(i, iCheck), dblTmp, nNonzero
          do k=1, min(nNonzero,50)
            i1 = idxYdotOneReactionSorted(k)
            if (abs(ydotOneSpecies(i1)) .LE. abs(1D-3 * dblTmp)) &
              exit
            dblTmp1 = dblTmp1 + ydotOneSpecies(i1)
            write (fU, &
              '(4X, I3, ES11.2, 2ES10.2, 2X, I5, ": " &
              & 7A8, 2X, ES8.2, 2F8.2)') &
              k, ydotOneSpecies(i1), &
              ydotOneSpecies(i1)/dblTmp, dblTmp1/dblTmp, &
              idxReactionsNonzero(i1), &
              strReactants(1:2, idxReactionsNonzero(i1)), ' -> ', &
              strProducts(1:4, idxReactionsNonzero(i1)), &
              dblABC(:, idxReactionsNonzero(i1))
          end do
        end do
      end do

      close (UNIT=fU, IOSTAT=i, STATUS='KEEP')
      end subroutine AnalyseReaction




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
        write (*,'(A/)') 'Saving reactions...'
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
        write (*,'(A/)') 'Saving species...'
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
        write (*,'(A/)') 'Saving rates...'
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




      subroutine initialize (fileInitialize)
      ! Read all the configuration parameters
      use CMDT
      implicit none
      character(LEN=128) fileInitialize
      integer fU, ios
      logical getFileUnit

      namelist /PhysicalParameters/ &
        Temperature, n_H, &
        Av, omega_Albedo, rateCosIon, MeanMolWeight, &
        BarrierWidth, ratioGrainToH, ratioDustGasMass, &
        GrainRadius, GrainDensity, SitesDensity, &
        nSpecies_Est, nIteration, timeUnit, sto_threshold, &
        rateThreshold, nOrderLim

      namelist /ODEParameters/ &
        ATOL, RTOL, nIteration

      namelist /Paths/ &
        path, fReactions, fMobility, fDissRecom, &
        fInitialCondition, &
        fReactionsSave, fSpeciesSave, fRatesSave, &
        fSaveFinalResult, fSaveElementalAb, fSaveALLHistory, &
        fAnalyseResult, fSaveALLHistory_ASCII, fFinalJac, &
        fNameMoments, fEqMoments, fOutLog, fPhyParHistory

      if (.NOT. getFileUnit(fU)) then
        write (*,*) 'Cannot allocate an output file unit!'
        stop
      end if

      CALL openFileSequentialRead (fU, fileInitialize, 999)
      read (UNIT=fU, IOSTAT=ios, NML=PhysicalParameters)
      read (UNIT=fU, IOSTAT=ios, NML=ODEParameters)
      read (UNIT=fU, IOSTAT=ios, NML=Paths)
      close (UNIT=fU, IOSTAT=ios, STATUS='KEEP')

      end subroutine initialize
