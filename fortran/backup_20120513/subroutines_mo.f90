      subroutine f (NEQ, t, y, ydot)
      use CMDT
      use CMTP
      implicit none
      integer NEQ, i, j, i1
      double precision t, y(NEQ), ydot(NEQ)

      ydot = 0D0

      call updateMantSurfRates(NEQ, y)

      do i=1, nSpecies
        do j=1, ndtMM(i)
          ydot(i) = ydot(i) + rates(dtMM(2, j, i)) * &
            DBLE(dtMM(3, j, i)) * &
            getMomValue(dtMM(1, j, i), i, y, NEQ)
          !     getMMProd(dtMM(1, j, i), y, NEQ)
        end do
      end do
      do i=nSpecies+1, NEQ
        i1 = indY(i)
        do j=1, ndtMM(i1)
          ydot(i) = ydot(i) + rates(dtMM(2, j, i1)) * &
            DBLE(dtMM(3, j, i1)) * &
            getMomValue(dtMM(1, j, i1), i, y, NEQ)
        end do
      end do
      if (useThreeBody) then
        call getThreeBodySurfRate(NEQ, y)
        do i=1, nThreeBodySurfReacs
          i1 = ThreeBodySurfReacs(i)
          ydot(reactions(1:nRealReactants(i1), i1)) = &
            ydot(reactions(1:nRealReactants(i1), i1)) &
              - rateThreeBodySurf(i)
          ydot(reactions(1+nReactants: &
            nReactants+nRealProducts(i1), i1)) = &
            ydot(reactions(1+nReactants: &
              nReactants+nRealProducts(i1), i1)) &
            + rateThreeBodySurf(i)
        end do
      end if
      end subroutine f




      subroutine updateMantSurfRates (NEQ, y)
      use CMDT
      use CMTP
      implicit none
      integer NEQ, i, i1
      double precision, dimension(NEQ) :: y
      tot_man_mant_ab = sum(y(mantleSpe%idxMant))
      tot_man_surf_ab = sum(y(mantleSpe%idxSurf))
      tot_man_acc_rate = dot_product(y(mantleAccReacs%iS), &
                                 rates(mantleAccReacs%iR))
      tot_man_eva_rate = dot_product(y(mantleEvaReacs%iS), &
                                 rates(mantleEvaReacs%iR))
      tot_man_acc_rate_reac = 0D0
      tot_man_eva_rate_reac = 0D0
      do i=1, n_surf_mant_inc
        i1 = surf_mant_inc(i)%idx
        tot_man_acc_rate_reac = tot_man_acc_rate_reac + &
          rates(i1) * surf_mant_inc(i)%gain_loss * &
            getMomValue(idxReacMM(i1), &
            reactions(1, i1),  y, NEQ)
      end do
      do i=1, n_surf_mant_dec
        i1 = surf_mant_dec(i)%idx
        tot_man_eva_rate_reac = tot_man_eva_rate_reac + &
          rates(i1) * surf_mant_dec(i)%gain_loss * &
            getMomValue(idxReacMM(i1), &
            reactions(1, i1),  y, NEQ)
      end do
      tot_man_acc_rate = tot_man_acc_rate + tot_man_acc_rate_reac
      tot_man_eva_rate = tot_man_eva_rate + tot_man_eva_rate_reac
      if (.NOT. allowAlwaysExpose) then
        if (tot_man_acc_rate .GT. tot_man_eva_rate) then
          tot_man_acc_rate = tot_man_acc_rate - tot_man_eva_rate
          tot_man_eva_rate = 0D0
        else
          tot_man_eva_rate = tot_man_eva_rate - tot_man_acc_rate
          tot_man_acc_rate = 0D0
        end if
      end if
      rates(mantleSpe%idx_S2M) = tot_man_acc_rate / SitesPerGrain
      rates(mantleSpe%idx_M2S) = tot_man_eva_rate / &
        max(tot_man_surf_ab, tot_man_mant_ab, tiny(1D0))
      end subroutine updateMantSurfRates


      subroutine updateMantSurfRates_obsolete (NEQ, y)
      use CMDT
      use CMTP
      implicit none
      integer NEQ, i, i1
      double precision, dimension(NEQ) :: y
      tot_man_mant_ab = sum(y(mantleSpe%idxMant))
      tot_man_surf_ab = sum(y(mantleSpe%idxSurf))
      tot_man_acc_rate = sum(y(mantleSpe%idxGas) * &
                             rates(mantleSpe%idx_acc))
      tot_man_eva_rate = sum(y(mantleSpe%idxSurf) * &
                             rates(mantleSpe%idx_eva))
      tot_man_acc_rate_reac = 0D0
      tot_man_eva_rate_reac = 0D0
      if (allowExposureByReaction) then
        do i=1, n_surf_mant_inc
          i1 = surf_mant_inc(i)%idx
          tot_man_acc_rate_reac = tot_man_acc_rate_reac + &
            rates(i1) * surf_mant_inc(i)%gain_loss * &
              getMomValue(idxReacMM(i1), &
              reactions(1, i1),  y, NEQ)
        end do
        do i=1, n_surf_mant_dec
          i1 = surf_mant_dec(i)%idx
          tot_man_eva_rate_reac = tot_man_eva_rate_reac + &
            rates(i1) * surf_mant_dec(i)%gain_loss * &
              getMomValue(idxReacMM(i1), &
              reactions(1, i1),  y, NEQ)
        end do
      else
        do i=1, n_surf_mant_inc
          i1 = surf_mant_inc(i)%idx
          if (IsGas(reactions(1+nReactants, i1))) &
            tot_man_acc_rate_reac = tot_man_acc_rate_reac + &
              rates(i1) * surf_mant_inc(i)%gain_loss * &
                getMomValue(idxReacMM(i1), &
                reactions(1, i1),  y, NEQ)
        end do
        do i=1, n_surf_mant_dec
          i1 = surf_mant_dec(i)%idx
          if (IsGas(reactions(1+nReactants, i1))) &
            tot_man_eva_rate_reac = tot_man_eva_rate_reac + &
              rates(i1) * surf_mant_dec(i)%gain_loss * &
                getMomValue(idxReacMM(i1), &
                reactions(1, i1),  y, NEQ)
        end do
      end if
      tot_man_acc_rate = tot_man_acc_rate + tot_man_acc_rate_reac
      tot_man_eva_rate = tot_man_eva_rate + tot_man_eva_rate_reac
      if (.NOT. allowAlwaysExpose) then
        if (tot_man_acc_rate .GT. tot_man_eva_rate) then
          tot_man_acc_rate = tot_man_acc_rate - tot_man_eva_rate
          tot_man_eva_rate = 0D0
        else
          tot_man_eva_rate = tot_man_eva_rate - tot_man_acc_rate
          tot_man_acc_rate = 0D0
        end if
      end if
      rates(mantleSpe%idx_S2M) = tot_man_acc_rate / SitesPerGrain
      rates(mantleSpe%idx_M2S) = tot_man_eva_rate / &
        max(tot_man_surf_ab, tot_man_mant_ab, tiny(1D0))
      end subroutine updateMantSurfRates_obsolete



      subroutine getThreeBodySurfRate(NEQ, y)
      use CMDT
      use CMTP
      implicit none
      integer NEQ, iReac, i, j, k, L, i1, i2
      double precision, dimension(NEQ) :: y
      double precision prefac, theta, rProd, rtmp
      rateThreeBodySurf = 0D0
      do j=1, nThreeBodySurfReacs
        iReac = ThreeBodySurfReacs(j)
        do i=1, nRealReactants(iReac)
          theta = min(1D0, y(reactions(nRealReactants(iReac)+1-i, &
            iReac)) / SitesPerGrain)
          i1 = reactions(i, iReac)
          if (useThreeBodyPrefac) then
            prefac = &
              BranchingRatios(iReac)*vibFreqSpecies(i1) / &
              (BranchingRatios(iReac)*vibFreqSpecies(i1) &
               + vibFreqSpecies(i1) * ( &
                   exp(-E_desorp_species(i1)/T_dust) &
                   + CosmicDesorpPreFactor * &
                     exp(-E_desorp_species(i1) / CosmicDesorpGrainT)) &
               + mobilitySpecies(i1))
          else
            prefac = 1D0
          end if
          rProd = 0D0
          do k=1, nReac_AsProducts(i1)
            i2 = idxReac_AsProducts(k, i1)
            rtmp = 1D0
            do L=1, nRealReactants(i2)
              rtmp = rtmp * y(reactions(L, i2))
            end do
            rProd = rProd + rates(i2) * rtmp
          end do
          rateThreeBodySurf(j) = rateThreeBodySurf(j) &
            + prefac * rProd * theta
        end do
      end do

      end subroutine getThreeBodySurfRate



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
      if (IsMan(iSpecies) .OR. IsSurfWithMan(iSpecies)) then
        call updateMantSurfRates(NEQ, y)
      end if

      do j=1, ndtMM(iSpecies)
        ydot(j) = rates(dtMM(2, j, iSpecies)) * &
          DBLE(dtMM(3, j, iSpecies)) * &
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
        idxYdotOneReactionSorted(max(nReactions, nHistoryLen)))

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
         !do i1=1, nNumMM
         !  if (ndtMM(i1) .GT. 0) then
         !    do j1=1, ndtMM(i1)
         !      dtMMMo(i1)%coeff(j1) = rates(dtMM(2, j1, i1)) * &
         !        DBLE(dtMM(3, j1, i1))
         !    end do
         !  end if
         !end do
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



      subroutine initialize (fileInitialize)
      ! Read all the configuration parameters
      use CMDT
      use CMTP
      implicit none
      character(LEN=128) fileInitialize
      integer fU, ios
      logical getFileUnit

      namelist /PhysicalParameters/ &
        T_gas, T_dust, n_H, &
        Av, omega_Albedo, rateCosIon, MeanMolWeight, &
        ReactBarrierWidth, DiffBarrierWidth, &
        ratioGrainToH, ratioDustGasMass, &
        GrainRadius, GrainDensity, SitesDensity, &
        nSpecies_Est, sto_threshold, &
        rateThreshold, nOrderLim, &
        CosmicDesorpGrainT, CosmicDesorpPreFactor, &
        ChemDesorpPreFactor, SimpleChemDesorption, &
        Diff2DesorRatio, tFinal, nNumMMEst, &
        useCosPhotoDiss, useZeroPointEnergy, &
        useThreeBody, useThreeBodyPrefac, &
        allowExposureByReaction, allowAlwaysExpose

      namelist /ODEParameters/ &
        ATOL, RTOL, nIteration

      namelist /Paths/ &
        path, fReactions, fDissRecom, &
        fInitialCondition, &
        fReactionsSave, fSpeciesSave, fRatesSave, &
        fSaveFinalResult, fSaveElementalAb, fSaveALLHistory, &
        fAnalyseResult, fSaveALLHistory_ASCII, fFinalJac, &
        fNameMoments, fEqMoments, fOutLog, fPhyParHistory, &
        fSpeciesEnthalpy, fSurfRadical

      if (.NOT. getFileUnit(fU)) then
        write (*,*) 'In initialize:'
        write (*,*) 'Cannot allocate an output file unit!'
        stop
      end if

      CALL openFileSequentialRead (fU, fileInitialize, 999)
      read (UNIT=fU, IOSTAT=ios, NML=PhysicalParameters)
      read (UNIT=fU, IOSTAT=ios, NML=ODEParameters)
      read (UNIT=fU, IOSTAT=ios, NML=Paths)
      close (UNIT=fU, IOSTAT=ios, STATUS='KEEP')

      end subroutine initialize
