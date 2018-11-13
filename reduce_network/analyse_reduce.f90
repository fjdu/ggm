      module analyse_reduce
      use CMDT
      implicit none
      integer n_time_care
      double precision, dimension(32) :: time_care
      integer, dimension(:), allocatable :: n_cons_reacs, n_prod_reacs
      integer, dimension(:,:), allocatable :: cons_reacs, prod_reacs
      double precision, dimension(:,:), allocatable :: &
        cons_rates, prod_rates
      integer n_imp_reacs
      integer, dimension(:), allocatable :: imp_reacs
      integer n_imp_species, n_imp_species_init
      integer, dimension(:), allocatable :: imp_species
      character(len=constLenNameSpecies), dimension(:), allocatable :: &
        name_imp_species
      character(len=128) file_imp_species_in, file_imp_reacs_out
      integer, parameter :: n_recursion = 64
      double precision, dimension(n_recursion) :: ratio_tobe_imp
      integer, dimension(:), allocatable :: idxSorted
      logical allow_dead_end
      double precision abundMaxHeavy, abundMaxSecondary
      integer, parameter :: nPhyParMax = 32
      integer nPhyPar
      double precision, dimension(nPhyParMax) :: Temperature_s, n_H_s
      integer, dimension(11), parameter :: &
        idx_HeavySpecies = (/ &
        7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17/)
      integer, parameter :: idxD=5, idxC = 7
      double precision tScaleMax, abundMaxC
      integer idx_record, i_rec, i_care, nHeavyEleMax,  nMaxD, nMaxC
      logical use_acum_threshold
      logical, dimension(:), allocatable :: &
        flag_imp_reacs, flag_care_spe, flag_reac_saved
      double precision, dimension(:), allocatable :: &
        priority_spe, priority_reac
      double precision prio_reac_threshold, prio_spe_threshold
      character(len=128) fLog
      integer fULog


      contains

      subroutine get_reacs_of_species()
      use CMDT
      implicit none
      integer i, j
      logical flag
      n_cons_reacs = 0
      n_prod_reacs = 0
      do i=1, nReactions
        do j=1, nRealReactants(i)
          flag = .FALSE.
          if (n_cons_reacs(reactions(j,i)) .EQ. 0) then
            flag = .TRUE.
          else if (cons_reacs(n_cons_reacs(reactions(j,i)), &
                   reactions(j,i)) .NE. i) then
            flag = .TRUE.
          end if
          if (flag) then
            n_cons_reacs(reactions(j, i)) = &
              n_cons_reacs(reactions(j, i)) + 1
            cons_reacs(n_cons_reacs(reactions(j, i)), &
                reactions(j, i)) = i
          end if
        end do
        do j=nReactants+1, nReactants+nRealProducts(i)
          flag = .FALSE.
          if (n_prod_reacs(reactions(j,i)) .EQ. 0) then
            flag = .TRUE.
          else if (prod_reacs(n_prod_reacs(reactions(j,i)), &
                   reactions(j,i)) .NE. i) then
            flag = .TRUE.
          end if
          if (flag) then
            n_prod_reacs(reactions(j, i)) = &
              n_prod_reacs(reactions(j, i)) + 1
            prod_reacs(n_prod_reacs(reactions(j, i)), &
                reactions(j, i)) = i
          end if
        end do
      end do
      end subroutine get_reacs_of_species



      subroutine get_cons_rates(y, ny)
      use CMDT
      implicit none
      integer i,j, ny
      double precision, dimension(ny) :: y
      double precision tmp
      do i=1, nSpecies
        tmp = y(i)
        y(i) = 1D0
        do j=1, n_cons_reacs(i)
          cons_rates(j,i) = get_a_reac_rate(cons_reacs(j,i), y, ny)
        end do
        y(i) = tmp
      end do
      end subroutine get_cons_rates


      subroutine get_prod_rates(y, ny)
      use CMDT
      implicit none
      integer i,j, ny
      double precision, dimension(ny) :: y
      do i=1, nSpecies
        do j=1, n_prod_reacs(i)
          prod_rates(j,i) = get_a_reac_rate(prod_reacs(j,i), y, ny)
        end do
      end do
      end subroutine get_prod_rates


      double precision function get_a_reac_rate(i, y, ny)
      use CMDT
      implicit none
      integer i, ny
      double precision get_a_reac_rate
      double precision, dimension(ny) :: y

      select case (iType(i))
      case (5)
        get_a_reac_rate = rates(i) * y(reactions(1, i)) &
            * y(reactions(2, i))
      case (1, 2, 3)
        get_a_reac_rate = rates(i) * y(reactions(1, i))
      case (0)
        get_a_reac_rate = rates(i) * y(reactions(2, i))
      case (53)
        if (nRealReactants(i) .EQ. 1) then
          get_a_reac_rate = rates(i) * y(reactions(1, i))
        end if
        if (nRealReactants(i) .EQ. 2) then
          get_a_reac_rate = rates(i) * y(reactions(1, i)) &
            * y(reactions(2, i))
        end if
      end select
      end function get_a_reac_rate


      subroutine get_imp_reacs()
      use CMDT
      implicit none
      integer, dimension(:), allocatable :: idx_time_care
      integer idx_bg_spe, idx_ed_spe, idx_bg_reac, idx_ed_reac
      integer n_tmp
      logical IsWordChar
      if (.NOT. allocated(idxSorted)) &
        allocate(idxSorted(nReactions))
      do i_care=1, n_time_care
        idx_record = get_closest_record(time_care(i_care), &
          touts, nHistoryLen)
        call get_cons_rates(yHistory(:, idx_record), nSpecies)
        call get_prod_rates(yHistory(:, idx_record), nSpecies)
        idx_bg_spe = 1
        idx_ed_spe = n_imp_species
        if (IsWordChar(fLog(1:1))) &
          write (fULog,'("Time = ", ES8.2)') touts(idx_record)
        do i_rec=1, n_recursion
          n_tmp = n_imp_species
          if (IsWordChar(fLog(1:1))) &
            write (fULog,'("i_rec = ", I3)') i_rec
          call add_imp_reacs_spe &
            (idx_bg_spe, idx_ed_spe, idx_bg_reac, idx_ed_reac)
          call add_imp_spe_reac &
            (idx_bg_reac, idx_ed_reac, idx_bg_spe, idx_ed_spe)
          if (n_tmp .EQ. n_imp_species) then
            write (*,'(2X, "At time=", ES8.2, 2X, "i_rec=", I2, &
              & 2X, "n_imp_reacs=", I5, 2X, "n_imp_species=", I4)') &
              touts(idx_record), i_rec, n_imp_reacs, n_imp_species
            exit
          end if
        end do
      end do
      end subroutine get_imp_reacs


      subroutine add_imp_reacs_spe &
        (i_bg_spe, i_ed_spe, i_bg_reac, i_ed_reac)
      use CMDT
      implicit none
      integer i_bg_spe, i_ed_spe, i_bg_reac, i_ed_reac
      integer i, j, i_spe, idx_ed
      double precision rate_sum, threshold, &
        priority_reac_here
      logical IsWordChar
      i_bg_reac = n_imp_reacs + 1
      do i=i_bg_spe, i_ed_spe
        i_spe = imp_species(i)
        if (IsWordChar(fLog(1:1))) &
          write (fULog, '(A12, 2X, F8.3)') &
            nameSpecies(i_spe), priority_spe(i_spe)
        call SORT_Dec_idx(n_cons_reacs(i_spe), &
          cons_rates(1:n_cons_reacs(i_spe), i_spe), idxSorted)
        if (use_acum_threshold) then
          rate_sum = sum(cons_rates(1:n_cons_reacs(i_spe), i_spe))
          threshold = rate_sum * ratio_tobe_imp(i_rec)
          call get_idx_range_for_Acc &
            (cons_rates(1:n_cons_reacs(i_spe), i_spe), &
             idxSorted(1:n_cons_reacs(i_spe)), &
             n_cons_reacs(i_spe), threshold, idx_ed)
        else
          threshold = ratio_tobe_imp(i_rec) * &
            cons_rates(idxSorted(1), i_spe)
          call get_idx_range_max &
            (cons_rates(1:n_cons_reacs(i_spe), i_spe), &
             idxSorted(1:n_cons_reacs(i_spe)), &
             n_cons_reacs(i_spe), threshold, idx_ed)
        end if
        if (IsWordChar(fLog(1:1))) write (fULog, '(2X,"=>")')
        do j=1, idx_ed
          if (1D0 / cons_rates(idxSorted(j), i_spe) &
              .GT. tScaleMax) cycle
          priority_reac_here = priority_spe(i_spe) * &
            cons_rates(idxSorted(j), i_spe) / &
            (cons_rates(idxSorted(1), i_spe) + 1D-100)
          call add_to_imp_reac(cons_reacs(idxSorted(j), i_spe), &
              priority_reac_here)
          if (IsWordChar(fLog(1:1))) &
            write (fULog, '(4X, F8.3, 3X, 7A12)') priority_reac_here, &
              strReactants(:, cons_reacs(idxSorted(j), i_spe)), &
              strProducts(:, cons_reacs(idxSorted(j), i_spe))
        end do
        call SORT_Dec_idx(n_prod_reacs(i_spe), &
          prod_rates(1:n_prod_reacs(i_spe), i_spe), idxSorted)
        if (use_acum_threshold) then
          rate_sum = sum(prod_rates(1:n_prod_reacs(i_spe), i_spe))
          threshold = rate_sum * ratio_tobe_imp(i_rec)
          call get_idx_range_for_Acc &
            (prod_rates(1:n_prod_reacs(i_spe), i_spe), &
             idxSorted(1:n_prod_reacs(i_spe)), &
             n_prod_reacs(i_spe), threshold, idx_ed)
        else
          threshold = ratio_tobe_imp(i_rec) * &
            prod_rates(idxSorted(1), i_spe)
          call get_idx_range_max &
            (prod_rates(1:n_prod_reacs(i_spe), i_spe), &
             idxSorted(1:n_prod_reacs(i_spe)), &
             n_prod_reacs(i_spe), threshold, idx_ed)
        end if
        if (IsWordChar(fLog(1:1))) write (fULog, '(2X,"<=")')
        do j=1, idx_ed
          priority_reac_here = priority_spe(i_spe) * &
            prod_rates(idxSorted(j), i_spe) / &
            (prod_rates(idxSorted(1), i_spe) + 1D-100)
          call add_to_imp_reac(prod_reacs(idxSorted(j), i_spe), &
            priority_reac_here)
          if (IsWordChar(fLog(1:1))) &
            write (fULog, '(4X, F8.3, 3X, 7A12)') priority_reac_here, &
              strReactants(:, prod_reacs(idxSorted(j), i_spe)), &
              strProducts(:, prod_reacs(idxSorted(j), i_spe))
        end do
      end do
      i_ed_reac = n_imp_reacs
      end subroutine add_imp_reacs_spe


      subroutine add_imp_spe_reac &
            (i_bg_reac, i_ed_reac, i_bg_spe, i_ed_spe)
      use CMDT
      implicit none
      integer i_bg_spe, i_ed_spe, i_bg_reac, i_ed_reac
      integer i, j, i_reac
      double precision priority_spe_here
      i_bg_spe = n_imp_species + 1
      do i=i_bg_reac, i_ed_reac
        i_reac = imp_reacs(i)
        priority_spe_here = priority_reac(i_reac)
        do j=1, nRealReactants(i_reac)
          call add_to_imp_spe(reactions(j, i_reac), priority_spe_here)
        end do
        if (.NOT. allow_dead_end) then
        do j=nReactants+1, nReactants+nRealProducts(i_reac)
          call add_to_imp_spe(reactions(j, i_reac), priority_spe_here)
        end do
        end if
      end do
      i_ed_spe = n_imp_species
      end subroutine add_imp_spe_reac


      subroutine add_to_imp_reac(idx_reac, priority_reac_here)
      use CMDT
      implicit none
      integer i, idx_reac
      logical flag
      double precision priority_reac_here
      if (priority_reac_here .LT. prio_reac_threshold) return
      flag = .TRUE.
      do i=1, n_imp_reacs
        if (imp_reacs(i) .EQ. idx_reac) then
          flag = .FALSE.
          exit
        end if
      end do
      if ((rates(idx_reac) .GE. 0D0) .AND. &
          (nRealReactants(idx_reac) .GT. 0) .AND. &
          (nRealProducts(idx_reac) .GT. 0) .AND. &
          (flag)) then
        n_imp_reacs = n_imp_reacs + 1
        imp_reacs(n_imp_reacs) = idx_reac
        priority_reac(idx_reac) = max(priority_reac(idx_reac), &
          priority_reac_here)
      end if
      end subroutine add_to_imp_reac


      subroutine add_to_imp_spe(idx_spe, priority_spe_here)
      use CMDT
      implicit none
      integer i, idx_spe
      logical flag
      double precision priority_spe_here
      if (priority_spe_here .LT. prio_spe_threshold) return
      flag = .TRUE.
      do i=1, n_imp_species
        if (imp_species(i) .EQ. idx_spe) then
          flag = .FALSE.
          exit
        end if
      end do
      if (flag) then
        ! Caution must be taken.
        if ((sum(SpeciesElements(idx_HeavySpecies, idx_spe)) &
            .GT. nHeavyEleMax) .AND. &
            (yHistory(idx_spe, idx_record) .LT. abundMaxHeavy)) &
          return
        if (SpeciesElements(idxC, idx_spe) .GT. nMaxC) then
          if (yHistory(idx_spe, idx_record) .LT. abundMaxC) &
            return
        end if
        if (SpeciesElements(idxD, idx_spe) .GT. nMaxD) return
        if (i_rec .GT. 1) then
          if (yHistory(idx_spe, idx_record) .LT. &
            abundMaxSecondary) return
        end if
        n_imp_species = n_imp_species + 1
        imp_species(n_imp_species) = idx_spe
        priority_spe(idx_spe) = max(priority_spe(idx_spe), &
          priority_spe_here)
      end if
      end subroutine add_to_imp_spe




      subroutine import_reac_care()
      use CMDT
      implicit none
      integer fU, i, j, nLineAll, nLineData, ios
      character(len=constLenNameSpecies) strTmp
      fU = 1
      CALL openFileSequentialRead &
          (fU, trim(path)//trim(file_imp_species_in), 999)
      CALL GetNLineF (fU, nLineAll, nLineData, '!')
      allocate(name_imp_species(nLineData))
      rewind (UNIT=fU, IOSTAT=ios)
      n_imp_species = 0
      flag_care_spe = .FALSE.
      do i=1, nLineAll
        read (UNIT=fU, FMT='(A12)', IOSTAT=ios) strTmp
        if (ios .NE. 0) then
          write (*,*) 'ios = ', ios
          stop
        end if
        do j=1, nSpecies
          if (trim(nameSpecies(j)) .EQ. trim(strTmp)) then
            n_imp_species = n_imp_species + 1
            imp_species(n_imp_species) = j
            name_imp_species(n_imp_species) = strTmp
            flag_care_spe(j) = .TRUE.
            exit
          end if
        end do
      end do
      write (*,'(/A, I8)') 'Number of species you most care about: ', &
        nLineData
      write (*,'(A, I8)') 'Actually used care species: ', &
        n_imp_species
      n_imp_species_init = n_imp_species
      close (UNIT=fU, IOSTAT=ios, STATUS='KEEP')
      end subroutine import_reac_care



      subroutine save_reacs(fU, i_reac)
      use CMDT
      implicit none
      integer fU, i_reac
      if (.NOT. flag_imp_reacs(i_reac)) return
      if (flag_reac_saved(i_reac)) return
      if ((nRealReactants(i_reac) .EQ. 0) .OR. &
          (nRealProducts(i_reac) .EQ. 0)) then
        write (*,*) strReactants(:, i_reac), strProducts(:, i_reac)
        return
      else
        strReactants(nRealReactants(i_reac)+1:nReactants, i_reac) = ''
        strProducts(nRealProducts(i_reac)+1:nProducts, i_reac) = ''
        write (fU, '(7(A12), ES8.2E2, F9.2, F9.1, 2F6.0, I3, &
          &X, A1, X, A2)') &
          strReactants(:,i_reac), strProducts(:,i_reac), &
          dblABC(:,i_reac), dblTLTU(:, i_reac), iType(i_reac), &
          strQuality(i_reac), strType(i_reac)
        flag_reac_saved(i_reac) = .TRUE.
      end if
      end subroutine save_reacs



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

      subroutine SORT_Dec_idx (nDim, Y, idxSorted)
      implicit none
      integer nDim, i, j, itmp
      integer idxSorted(nDim)
      double precision Y(nDim)
      do i=1, nDim
        idxSorted(i) = i
      end do
      do i=1, nDim
        do j=i+1, nDim
          if (Y(idxSorted(j)) .GT. Y(idxSorted(i))) then
            itmp = idxSorted(i)
            idxSorted(i) = idxSorted(j)
            idxSorted(j) = itmp
          end if
        end do
      end do
      end subroutine SORT_Dec_idx


      integer function get_closest_record(key, vec, Len)
      implicit none
      double precision key
      integer Len
      double precision, dimension(Len) :: vec
      integer get_closest_record
      integer, dimension(2) :: tmp_vec
!     tmp_vec = minloc(abs(vec-key))
      get_closest_record = minloc(abs(vec-key), 1)!tmp_vec(1)
      end function  get_closest_record


      subroutine get_idx_range_max(vec, idx_sorted, &
        nLen, threshold, idx_ed)
      implicit none
      integer nLen, idx_ed, i
      double precision, dimension(nLen) :: vec
      integer, dimension(nLen) :: idx_sorted
      double precision threshold
      idx_ed = nLen
      do i=1, nLen
        if (vec(idx_sorted(i)) .LE. threshold) then
          idx_ed = i
          exit
        end if
      end do
      end subroutine get_idx_range_max


      subroutine get_idx_range_for_Acc(vec, idx_sorted, &
        nLen, threshold, idx_ed)
      implicit none
      integer nLen, idx_ed, i
      double precision, dimension(nLen) :: vec
      integer, dimension(nLen) :: idx_sorted
      double precision threshold, acc
      acc = 0D0
      idx_ed = nLen
      do i=1, nLen
        acc = acc + vec(idx_sorted(i))
        if (acc .GE. threshold) then
          idx_ed = i
          exit
        end if
      end do
      end subroutine get_idx_range_for_Acc

      end module analyse_reduce
