      module analyse_reduce
      use CMDT
      implicit none
      integer n_time_care
      double precision, dimension(8) :: time_care
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
      integer n_recursion
      double precision ratio_tobe_imp
      integer, dimension(:), allocatable :: idxSorted


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
      do i=1, nSpecies
        do j=1, n_cons_reacs(i)
          cons_rates(j,i) = get_a_reac_rate(cons_reacs(j,i), y, ny)
        end do
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
      integer idx_record, i_rec, i_care
      integer n_tmp
      n_imp_reacs = 0
      allocate(idxSorted(nReactions))
      do i_care=1, n_time_care
        idx_record = get_closest_record(time_care(i_care), &
          touts, nHistoryLen)
        call get_cons_rates(yHistory(:, idx_record), nSpecies)
        call get_prod_rates(yHistory(:, idx_record), nSpecies)
        idx_bg_spe = 1
        idx_ed_spe = n_imp_species
        do i_rec=1, n_recursion
          n_tmp = n_imp_species
          call add_imp_reacs_spe &
            (idx_bg_spe, idx_ed_spe, idx_bg_reac, idx_ed_reac)
          call add_imp_spe_reac &
            (idx_bg_reac, idx_ed_reac, idx_bg_spe, idx_ed_spe)
          if (n_tmp .EQ. n_imp_species) then
            write (*,'(A16, ES12.2)') 'For time = ', touts(idx_record)
            write (*,'(20X, "i_rec = ", I6)') i_rec
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
      double precision rate_sum, threshold
      i_bg_reac = n_imp_reacs + 1
      do i=i_bg_spe, i_ed_spe
        i_spe = imp_species(i)
        call SORT_Dec_idx(n_cons_reacs(i_spe), &
          cons_rates(1:n_cons_reacs(i_spe), i_spe), idxSorted)
        rate_sum = sum(cons_rates(1:n_cons_reacs(i_spe), i_spe))
        threshold = rate_sum * ratio_tobe_imp
        call get_idx_range_for_Acc &
          (cons_rates(1:n_cons_reacs(i_spe), i_spe), &
           idxSorted(1:n_cons_reacs(i_spe)), &
           n_cons_reacs(i_spe), &
           threshold, idx_ed)
        do j=1, idx_ed
          call add_to_imp_reac(cons_reacs(idxSorted(j), i_spe))
        end do
        call SORT_Dec_idx(n_prod_reacs(i_spe), &
          prod_rates(1:n_prod_reacs(i_spe), i_spe), idxSorted)
        rate_sum = sum(prod_rates(1:n_prod_reacs(i_spe), i_spe))
        threshold = rate_sum * ratio_tobe_imp
        call get_idx_range_for_Acc &
          (prod_rates(1:n_prod_reacs(i_spe), i_spe), &
           idxSorted(1:n_prod_reacs(i_spe)), &
           n_prod_reacs(i_spe), &
           threshold, idx_ed)
        do j=1, idx_ed
          call add_to_imp_reac(prod_reacs(idxSorted(j), i_spe))
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
      i_bg_spe = n_imp_species + 1
      do i=i_bg_reac, i_ed_reac
        i_reac = imp_reacs(i)
        do j=1, nRealReactants(i_reac)
          call add_to_imp_spe(reactions(j, i_reac))
        end do
        do j=nReactants+1, nReactants+nRealProducts(i_reac)
          call add_to_imp_spe(reactions(j, i_reac))
        end do
      end do
      i_ed_spe = n_imp_species
      end subroutine add_imp_spe_reac


      subroutine add_to_imp_reac(idx_reac)
      use CMDT
      implicit none
      integer i, idx_reac
      logical flag
      flag = .TRUE.
      do i=1, n_imp_reacs
        if (imp_reacs(i) .EQ. idx_reac) then
          flag = .FALSE.
          exit
        end if
      end do
      if (flag) then
        n_imp_reacs = n_imp_reacs + 1
        imp_reacs(n_imp_reacs) = idx_reac
      end if
      end subroutine add_to_imp_reac


      subroutine add_to_imp_spe(idx_spe)
      use CMDT
      implicit none
      integer i, idx_spe
      logical flag
      flag = .TRUE.
      do i=1, n_imp_species
        if (imp_species(i) .EQ. idx_spe) then
          flag = .FALSE.
          exit
        end if
      end do
      if (flag) then
        n_imp_species = n_imp_species + 1
        imp_species(n_imp_species) = idx_spe
      end if
      end subroutine add_to_imp_spe




      subroutine import_reac_care()
      use CMDT
      implicit none
      integer fU, i, i1, j, nLineAll, nLineData, ios
      logical flag
      fU = 1
      CALL openFileSequentialRead &
          (fU, trim(path)//trim(file_imp_species_in), 999)
      CALL GetNLineF (fU, nLineAll, nLineData, '!')
      n_imp_species = nLineData
      allocate(name_imp_species(n_imp_species))
      write (*,'(/A, I8)') 'Number of species you most care about: ', &
          n_imp_species
      n_imp_species_init = n_imp_species
      rewind (UNIT=fU, IOSTAT=ios)
      i1 = 0
      do i=1, n_imp_species
        read (UNIT=fU, FMT='(A12)', IOSTAT=ios) name_imp_species(i)
        if (ios .NE. 0) then
          write (*,*) 'ios = ', ios
          stop
        end if
        flag = .FALSE.
        do j=1, nSpecies
          if (trim(nameSpecies(j)) .EQ. &
            trim(name_imp_species(i))) then
            imp_species(i) = j
            i1 = i1 + 1
            flag = .TRUE.
            exit
          end if
        end do
        if (.NOT. flag) write (*, '(I4, 2X, 2A10)') &
          i, name_imp_species(i), " is not in the network!"
      end do
      write (*,'(/A, I8)') 'Actually used initial species: ', i1
      close (UNIT=fU, IOSTAT=ios, STATUS='KEEP')
      end subroutine import_reac_care



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


      subroutine get_idx_range_for_Acc(vec, idx_sorted, &
        nLen, threshold, idx_ed)
      implicit none
      integer nLen, idx_ed, i
      double precision, dimension(nLen) :: vec
      integer, dimension(nLen) :: idx_sorted
      double precision threshold, acc
      acc = 0D0
      idx_ed = 0
      do i=1, nLen
        acc = acc + vec(idx_sorted(i))
        if (acc .GE. threshold) then
          idx_ed = i
          exit
        end if
      end do
      end subroutine get_idx_range_for_Acc

      end module analyse_reduce
