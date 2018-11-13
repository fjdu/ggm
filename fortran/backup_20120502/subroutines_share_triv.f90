!  Open a file for sequential read.
      subroutine openFileSequentialRead (fU, filename, maxRecLen)
      implicit none
      integer fU, maxRecLen, ios
      character(len=*) filename
      open (UNIT=fU, FILE=trim(filename), IOSTAT=ios, &
           STATUS='OLD', ACCESS='SEQUENTIAL', FORM='FORMATTED', &
           RECL=maxRecLen, BLANK='NULL', POSITION='REWIND', &
           ACTION='READ', DELIM='NONE', PAD='YES')
      if (ios .NE. 0) then
        write (*, '(A)') 'In openFileSequentialRead:'
        write (*, '(/A, I8, /A, A)') 'Open File Error: IOSTAT=', ios, &
          'Filename: ', filename
        ! Some people say that a subroutine should not 
        ! terminate the whole program. But, ...
        stop
      end if
      end subroutine openFileSequentialRead



!  Open a file for sequential write.
      subroutine openFileSequentialWrite (fU, filename, maxRecLen)
      implicit none
      integer fU, maxRecLen, ios
      character(len=*) filename
      open (UNIT=fU, FILE=trim(filename), IOSTAT=ios, &
           STATUS='REPLACE', ACCESS='SEQUENTIAL', FORM='FORMATTED', &
           RECL=maxRecLen, BLANK='NULL', POSITION='REWIND', &
           ACTION='WRITE', DELIM='NONE', PAD='YES')
      if (ios .NE. 0) then
        write (*, '(A)') 'In openFileSequentialWrite'
        write (*, '(/A, I8, /A, A)') 'Open File Error: IOSTAT=', ios, &
          'Filename: ', filename
        stop
      end if
      end subroutine openFileSequentialWrite



!  Get a file unit to output.
      logical function getFileUnit (fU)
      implicit none
      integer fU
      logical uExist, uOpened, getFileUnit
      ! 0: standard error
      ! 5: standard input
      ! 6: standard output
      ! 101:
      ! 102:
      do fU=10, 99
        inquire (unit=fU, exist=uExist, opened=uOpened)
        if (uExist .AND. .NOT. (uOpened)) then
          getFileUnit = .TRUE.
          return
        end if
      end do
      getFileUnit = .FALSE.
      fU = -1
      return
      end function getFileUnit



!  Check if a file unit is opened.
      logical function FileUnitOpened (fU)
      implicit none
      integer fU
      logical FileUnitOpened
      inquire (unit=fU, opened=FileUnitOpened)
      return
      end function FileUnitOpened



      integer function getIdxSpecies (SpeciesName)
      use CMDT
      implicit none
      integer getIdxSpecies
      character(*) SpeciesName
      do getIdxSpecies=1, nSpecies
        if (nameSpecies(getIdxSpecies) .EQ. SpeciesName) then
          return
        end if
      end do
      getIdxSpecies = -1
      return
      end function getIdxSpecies



      !  
      function getFilePreName (strFileName)
      implicit none
      integer ntrim, i
      character(len=*) strFileName
      character(len=128) getFilePreName
      ntrim = len_trim(strFileName)
      do i=ntrim, 1, -1
        if (strFileName(i:i) .EQ. '.') exit
      end do
      if (i .EQ. 1) then
        getFilePreName = strFileName(1:ntrim)
      else
        getFilePreName = strFileName(1:(i-1))
      end if
      return
      end function getFilePreName



      !  
      function combineStrArr (strArr, ndim, strSep)
      implicit none
      integer ndim, i
      character(len=*), dimension(ndim) :: strArr
      character(len=*) strSep
      character(len=256) combineStrArr
      combineStrArr = strArr(1)
      do i=2, ndim
        if (len_trim(strArr(i)) .EQ. 0) cycle
        combineStrArr = trim(combineStrArr)//strSep//strArr(i)
      end do
      return
      end function combineStrArr




!  Get the number of lines and data lines in a file.
      subroutine GetNLineF (fU, nLineAll, nLineData, commentChar)
      implicit none
      integer fU, ios, nLineAll, nLineData
      character commentChar, strtmp
      rewind (UNIT=fU, IOSTAT=ios)
      nLineAll = 0
      nLineData = 0
      do
        read (UNIT=fU, FMT='(A1)', IOSTAT=ios) strtmp
        if (ios .LT. 0) exit
        nLineAll = nLineAll + 1
        ! Any lines beginning with the commentChar (a single letter)
        ! or a blanck space are treated as commented.
        ! DON'T use Tab; use space.
        ! Added on 2010-04-05.
        if ((strtmp(1:1) .EQ. commentChar) .OR. &
            (strtmp(1:1) .EQ. ' ')) &
          cycle
        nLineData = nLineData + 1
      end do
      end subroutine GetNLineF



!  Get the number of lines in a file.
      subroutine GetFileLen (fU, FileName, nFileLen)
      implicit none
      integer fU, nFileLen, ios
      character(len=*) FileName
      character strtmp
      CALL openFileSequentialRead &
        (fU, FileName, 999)
      nFileLen = 0
      do
        read (UNIT=fU, FMT='(A)', IOSTAT=ios) strtmp
        if (ios .LT. 0) exit
        nFileLen = nFileLen + 1
      end do
      close (UNIT=fU, IOSTAT=ios, STATUS='KEEP')
      end subroutine GetFileLen



!  Return TRUE if the character in in the range [0-9A-Za-z_]
      logical function IsWordChar (ch)
      implicit none
      character ch
      logical IsWordChar
      IsWordChar = &
        (LGE(ch, '0') .AND. LLE(ch, '9')) .OR. &
        (LGE(ch, 'A') .AND. LLE(ch, 'Z')) .OR. &
        (LGE(ch, 'a') .AND. LLE(ch, 'z')) .OR. &
        (ch .EQ. '_')
      return
      end function IsWordChar



      subroutine waittostop
      implicit none
      character charTMP
      write (*,'(/A)') 'Press any key then <Enter> to exit:'
      read (*,*) charTMP
      stop
      end subroutine



      subroutine waithere
      implicit none
      character charTMP
      write (*,'(/A)') 'Press any key then <Enter> to continue:'
      read (*,*) charTMP
      end subroutine



! Count the number of occurrences of Ch in Str
      FUNCTION CharCountsInStr(Str, Ch)
      implicit none
      character(len=*) Str
      character Ch
      integer :: CharCountsInStr, nTmp, i
      CharCountsInStr = 0
      nTmp = len(Str)
      do i=1, nTmp
        if (Str(i:i) .EQ. Ch) &
          CharCountsInStr = CharCountsInStr + 1
      end do
      END FUNCTION CharCountsInStr



! Sort in ascending order
! The vector to be sorted is not changed, only the index
! of the sorted array is returned.
! To be optimized.
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



! Sort in ascending order
! The vector to be sorted is not changed, only the index
! of the sorted array is returned.
! To be optimized.
      subroutine SORT_Asc_idx_Int (nDim, Y, idxSorted)
      implicit none
      integer nDim, i, j, itmp
      integer idxSorted(nDim)
      integer Y(nDim)
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
      end subroutine SORT_Asc_idx_Int


      subroutine SORT_Asc_in_place (nDim, Y)
      implicit none
      integer nDim, i, j, itmp
      integer Y(nDim)
      do i=1, nDim
        do j=i+1, nDim
          if (Y(j) .LT. Y(i)) then
            itmp = Y(i)
            Y(i) = Y(j)
            Y(j) = itmp
          end if
        end do
      end do
      end subroutine SORT_Asc_in_place



! Get the elemental composition of each molecule.
      subroutine getElements &
        (nameSpec, listElements, nElements, arrNElements)
      character(len=*) nameSpec, listElements(nElements)
      integer, dimension(nElements) :: arrNElements
      integer i, j, k, ntmp, lenName, lenEle, nElements
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

      ! <timestamp>2011-06-10 Fri 23:16:38</timestamp>
      !     Correct a bug which has been hidden for a long time.
      !     Change the starting value of i from 1 into 2.
      do i=2, lenName
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



!  Get a double precision NaN.
      double precision function dblNaN ()
      implicit none
      double precision dblNaN, a
      a = 0D0
      dblNaN = a/a
      return
      end function dblNaN
