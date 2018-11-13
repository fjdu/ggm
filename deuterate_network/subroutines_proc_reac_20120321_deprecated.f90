! Created on 2011-02-01 Tue 18:17:19

MODULE CMReac

implicit none

integer, parameter :: nItemInASpecies = 32, nMaxGrp = 64, &
  idxH = 6, idxD = 7, nDeutSideMax = 16, lenStrSideMax = 64, &
  idxDumb = 20, nHDeutMax=5, nOtherDeutMax=5, nDeutDegree=2

TYPE :: EleSpecies
  integer nItem
  integer, dimension(nItemInASpecies) :: vecIdEle, vecCtEle
END TYPE EleSpecies

TYPE :: EleGroup
  integer nItem, nEle, nWeight
  integer, dimension(:), allocatable :: EleList
  integer, dimension(:,:), allocatable :: ArrCt
END TYPE EleGroup

type (EleSpecies), dimension(:), allocatable :: SpeciesEleAll

CONTAINS

subroutine Name2Ele (SpeciesName, SpeciesEle)
  use CMDT
  implicit none
  character(len=*) SpeciesName
  type (EleSpecies) SpeciesEle
  integer i, i1, i2, nlen
  nlen = len(SpeciesName)
  SpeciesEle%nItem=0
  i=1
  do
    ! Start from position i, try to match the longest possible element
    ! name, return the index of the matched element in i1, and the end
    ! position of the matched element in i2.
    if (i .GT. nlen) exit
    i1 = MatchEle(SpeciesName, nlen, i, i2)
    if (i1 .GT. 0) then
      SpeciesEle%nItem = SpeciesEle%nItem + 1
      SpeciesEle%vecIdEle(SpeciesEle%nItem) = i1
      SpeciesEle%vecCtEle(SpeciesEle%nItem) = &
        getThisEleCount(SpeciesName, nlen, i2+1, i1)
      i = i1 + 1
    else
      i = i + 1
    end if
  end do
  SpeciesEle%vecIdEle(SpeciesEle%nItem+1 : nItemInASpecies) = 0
  SpeciesEle%vecCtEle(SpeciesEle%nItem+1 : nItemInASpecies) = 0
end subroutine Name2Ele


integer function MatchEle(SpeciesName, nlen, i1, i2)
  use CMDT
  implicit none
  character(len=*) SpeciesName
  integer i, i1, i2, nlen, MatchEle, lenMatch, lenMatchTmp
  MatchEle = 0
  lenMatch = 0
  do i=1, nElement
    if (MatchStr(SpeciesName(i1:nlen), trim(nameElements(i)), &
      lenMatchTmp)) then
      if (lenMatchTmp .GT. lenMatch) then
        MatchEle = i
        lenMatch = lenMatchTmp
      end if
    end if
  end do
  if (MatchEle .GT. 0) then
    i2 = i1 + lenMatch - 1
  else
    i2 = i1
  end if
end function MatchEle


integer function getThisEleCount(SpeciesName, nlen, i1, i2)
  use CMDT
  implicit none
  character(len=*) SpeciesName
  integer i, i1, i2, nlen, getThisEleCount
  logical flagFound
  flagFound = .FALSE.
  do i=i1, nlen
    if (IsDigit(SpeciesName(i:i))) then
      i2 = i
      flagFound = .TRUE.
    else
      exit
    end if
  end do
  if (flagFound) then
    read (SpeciesName(i1:i2), '(I16)') getThisEleCount
  else
    getThisEleCount = 1
    i2 = i1 - 1
  end if
end function getThisEleCount


logical function MatchStr(Str1, Str2, len2)
! Determine if Str2 is a substring of Str1.
! The match point must be the leading position.
  implicit none
  character(len=*) Str1, Str2
  integer len1, len2, i
  logical MatchStr
  len1 = len_trim(Str1)
  len2 = len_trim(Str2)
  MatchStr = .TRUE.
  if (len1 .GE. len2) then
    do i=1, len2
      if (Str1(i:i) .NE. Str2(i:i)) then
        MatchStr = .FALSE.
        exit
      end if
    end do
  else
    MatchStr = .FALSE.
  end if
end function MatchStr


logical function IsDigit(ch)
  logical IsDigit
  character ch
  if (LGE(ch, '0') .AND. LLE(ch, '9')) then
    IsDigit = .TRUE.
  else
    IsDigit = .FALSE.
  end if
end function IsDigit


subroutine Ele2Elements(SpeciesEle, SpeElements)
  use CMDT
  implicit none
  type (EleSpecies) SpeciesEle
  integer, dimension(nElement) :: SpeElements
  integer i
  SpeElements = 0
  do i=1, SpeciesEle%nItem
    SpeElements(SpeciesEle%vecIdEle(i)) = &
      SpeElements(SpeciesEle%vecIdEle(i)) + &
      SpeciesEle%vecCtEle(i)
  end do
end subroutine Ele2Elements


logical function IsEquiv(SpeciesEle1, SpeciesEle2)
  use CMDT
  implicit none
  integer i, nLenTmp
  logical IsEquiv
  type (EleSpecies) SpeciesEle1, SpeciesEle2
  IsEquiv = .TRUE.
  if (SpeciesEle1%nItem .EQ. SpeciesEle2%nItem) then
    if (SpeciesEle1%nItem .EQ. 1) then
      if ((SpeciesEle1%vecIdEle(1) .NE. SpeciesEle2%vecIdEle(1)) .OR. &
          (SpeciesEle1%vecCtEle(1) .NE. SpeciesEle2%vecCtEle(1))) then
        IsEquiv = .FALSE.
      end if
    else
      ! H3COH will be considered to be equivalent to CH3OH.
      if (((SpeciesEle1%vecIdEle(1) .NE. SpeciesEle2%vecIdEle(1)) .OR. &
           (SpeciesEle1%vecCtEle(1) .NE. SpeciesEle2%vecCtEle(1)) .OR. &
           (SpeciesEle1%vecIdEle(2) .NE. SpeciesEle2%vecIdEle(2)) .OR. &
           (SpeciesEle1%vecCtEle(2) .NE. SpeciesEle2%vecCtEle(2))) &
           .AND. &
          ((SpeciesEle1%vecIdEle(1) .NE. SpeciesEle2%vecIdEle(2)) .OR. &
           (SpeciesEle1%vecCtEle(1) .NE. SpeciesEle2%vecCtEle(2)) .OR. &
           (SpeciesEle1%vecIdEle(2) .NE. SpeciesEle2%vecIdEle(1)) .OR. &
           (SpeciesEle1%vecCtEle(2) .NE. SpeciesEle2%vecCtEle(1)))) then
        IsEquiv = .FALSE.
      else
        do i=3, SpeciesEle1%nItem
          if ((SpeciesEle1%vecIdEle(i) .NE. SpeciesEle2%vecIdEle(i)) .OR. &
              (SpeciesEle1%vecCtEle(i) .NE. SpeciesEle2%vecCtEle(i))) then
            IsEquiv = .FALSE.
            exit
          end if
        end do
      end if
      ! With the following part, HCN and NCH will be considered equivalent.
      ! Use this part with caution!
      if (.NOT. IsEquiv) then
        IsEquiv = .TRUE.
        if ((SpeciesEle1%vecIdEle(SpeciesEle1%nItem) .EQ. 1) .OR. &
            (SpeciesEle1%vecIdEle(SpeciesEle1%nItem) .EQ. 2)) then
          nLenTmp = SpeciesEle1%nItem - 1
        else
          nLenTmp = SpeciesEle1%nItem
        end if
        ! Assume the charge symbol is always put at the last position.
        if ((SpeciesEle1%vecIdEle(SpeciesEle1%nItem) .NE. SpeciesEle2%vecIdEle(SpeciesEle1%nItem)) .OR. &
            (SpeciesEle1%vecCtEle(SpeciesEle1%nItem) .NE. SpeciesEle2%vecCtEle(SpeciesEle1%nItem))) then
          IsEquiv = .FALSE.
        else
          do i=1, nLenTmp
            if ((SpeciesEle1%vecIdEle(i) .NE. SpeciesEle2%vecIdEle(nLenTmp+1-i)) .OR. &
                (SpeciesEle1%vecCtEle(i) .NE. SpeciesEle2%vecCtEle(nLenTmp+1-i))) then
              IsEquiv = .FALSE.
              exit
            end if
          end do
        end if
      end if
      !!!!!!!!!
    end if
  else
    IsEquiv = .FALSE.
  end if
end function IsEquiv


subroutine DeutGroups (HydrGroup, nDeu, vecDeutGroup, nDeutGroup)
  use CMDT
  implicit none
  integer i, j, k
  integer nDeu, nDeutGroup, nTotalH, nGamb, binRepre, binRepreBak
  integer, dimension(:), allocatable :: ArrCt
  type (EleGroup), dimension(*) :: vecDeutGroup
  type (EleGroup) HydrGroup
  logical flagMatch
  nDeutGroup = 0
  nTotalH = sum(HydrGroup%ArrCt(1:HydrGroup%nItem, 1))
  if (nTotalH .LT. nDeu) then
    !write (*,*) 'nTotalH .LT. nDeu'
    return
  end if
  if (nTotalH .GT. 32) then
    write (*,*) 'nTotalH .GT. 32'
    return
  end if
  nGamb = binomcoeff(nTotalH, nDeu)
  allocate(ArrCt(HydrGroup%nItem))
  binRepre = ISHFT(1, nDeu) - 1
  do i=1, nGamb
    binRepreBak = binRepre
    ArrCt = 0
    do j=1, HydrGroup%nItem
      do k=1, HydrGroup%ArrCt(j,1)
        ArrCt(j) = ArrCt(j) + IAND(binRepre, 1)
        binRepre = ISHFT(binRepre, -1)
      end do
    end do
    binRepre = getNextColex(binRepreBak)
    flagMatch = .FALSE.
    do j=1, nDeutGroup
      flagMatch = .TRUE.
      do k=1, HydrGroup%nItem
        if (vecDeutGroup(j)%ArrCt(k,2) .NE. ArrCt(k)) then
          flagMatch = .FALSE.
          exit
        end if
      end do
      if (flagMatch) then
        vecDeutGroup(j)%nWeight = vecDeutGroup(j)%nWeight + 1
        exit
      end if
    end do
    if (.NOT. flagMatch) then
      nDeutGroup = nDeutGroup + 1
      vecDeutGroup(nDeutGroup)%nItem = HydrGroup%nItem
      vecDeutGroup(nDeutGroup)%nEle = 2
      vecDeutGroup(nDeutGroup)%nWeight = 1
      vecDeutGroup(nDeutGroup)%EleList(1) = idxH
      vecDeutGroup(nDeutGroup)%EleList(2) = idxD
      do j=1, HydrGroup%nItem
        vecDeutGroup(nDeutGroup)%ArrCt(j, 2) = ArrCt(j)
        vecDeutGroup(nDeutGroup)%ArrCt(j, 1) = &
          HydrGroup%ArrCt(j, 1) - ArrCt(j)
      end do
    end if
  end do
  deallocate(ArrCt)
end subroutine DeutGroups


subroutine DeutCleanGroups (vecDeutGroup, nDeutedSide, nDeutedSideCleaned)
  implicit none
  type (EleGroup), dimension(:) :: vecDeutGroup
  integer i, j, nDeutedSide, nDeutedSideCleaned
  nDeutedSideCleaned = nDeutedSide
  do i= nDeutedSide, 1, -1
    do j=1, i-1
      if (IsEquivSide(vecDeutGroup(i), vecDeutGroup(j))) then
        nDeutedSideCleaned = nDeutedSideCleaned - 1
        vecDeutGroup(j)%nWeight = vecDeutGroup(j)%nWeight + &
            vecDeutGroup(i)%nWeight
        vecDeutGroup(i)%nWeight = 0
        exit
      end if
    end do
  end do
end subroutine DeutCleanGroups



logical function IsEquivSide(EleGroup1, EleGroup2)
  implicit none
  type (EleGroup) EleGroup1, EleGroup2
  logical IsEquivSide
  integer, dimension(:), allocatable :: Ele
  if (EleGroup1%nItem .NE. EleGroup2%nItem) then
    IsEquivSide = .FALSE.
  else
    sort()
  end if
end function IsEquivSide



subroutine SortEleGroup(EleGroup1)
  implicit none
  type (EleGroup) EleGroup1
  integer i
  do i=1, EleGroup1%nItem
    do j=i+1, EleGroup1%nItem
      if (EleGroup1LT2(EleGroup1%ArrCt(j,:), EleGroup1%ArrCt(i,:))) then
        swap
      end if
    end do
  end do
end subroutine SortEleGroup





subroutine DeutIns (DeutGroup, SpeciesEle, SpeciesEleD)
  implicit none
  type (EleSpecies) SpeciesEle, SpeciesEleD
  type (EleGroup) DeutGroup
  integer i, i1, i2
  i1 = 0
  i2 = 0
  do i=1, SpeciesEle%nItem
    if (SpeciesEle%vecIdEle(i) .NE. idxH) then
      i1 = i1 + 1
      SpeciesEleD%vecIdEle(i1) = SpeciesEle%vecIdEle(i)
      SpeciesEleD%vecCtEle(i1) = SpeciesEle%vecCtEle(i)
    else
      i2 = i2 + 1
      i1 = i1 + 1
      SpeciesEleD%vecIdEle(i1) = idxH
      SpeciesEleD%vecCtEle(i1) = DeutGroup%ArrCt(i2,1)
      i1 = i1 + 1
      SpeciesEleD%vecIdEle(i1) = idxD
      SpeciesEleD%vecCtEle(i1) = DeutGroup%ArrCt(i2,2)
    end if
  end do
  SpeciesEleD%nItem = i1
end subroutine DeutIns


subroutine DeutReac (iReac, nDeut, fU)
  use CMDT
  implicit none
  integer i, j, iReac, nDeut, nDeutThis, nDeutedLeft, nDeutedRight, fU
  integer TotalWeight, nSplittedLeft, nSplittedRight
  character(len=lenStrSideMax), dimension(nDeutSideMax) :: &
    StrDeutedLeft, StrDeutedRight
  integer, dimension(nDeutSideMax) :: WeightsLeft, WeightsRight
  character(len=constLenNameSpecies), dimension(nReactants) :: &
    StrSplittedLeft
  character(len=constLenNameSpecies), dimension(nProducts) :: &
    StrSplittedRight

  do nDeutThis=1, nDeut
    call DeutSide(Reactions(1:nRealReactants(iReac), iReac), &
           nRealReactants(iReac), nDeutThis, &
           StrDeutedLeft, WeightsLeft, nDeutedLeft)
    call DeutSide( &
           Reactions(nReactants+1:nReactants+nRealProducts(iReac), iReac), &
           nRealProducts(iReac), nDeutThis, &
           StrDeutedRight, WeightsRight, nDeutedRight)
    TotalWeight = sum(WeightsRight(1:nDeutedRight))
    do i=1, nDeutedLeft
      StrSplittedLeft=""
      call SplitSideStr(StrDeutedLeft(i), StrSplittedLeft, nSplittedLeft)
      do j=1, nDeutedRight
        StrSplittedRight=""
        call SplitSideStr(StrDeutedRight(j), StrSplittedRight, nSplittedRight)
        write (fU, '(7A12, ES9.2, F9.2, F9.1, 12X, I3, 7X, "!", 4I3)') &
          StrSplittedLeft, StrSplittedRight, &
          dblABC(1, iReac) * dble(WeightsRight(j))/dble(TotalWeight), &
          dblABC(2:3, iReac), typeReac(iReac), &
          WeightsRight(j), j, i, nDeutThis
      end do
    end do
  end do
end subroutine DeutReac


subroutine DeutSide (SpeciesGroup, &
           nGroupSize, nDeutThis, &
           StrDeutedSide, WeightsSide, nDeutedSide)
  implicit none
  type (EleSpecies) SpeciesEleD, SpeciesEleBig
  type (EleGroup) HydrGroup
  type (EleGroup), dimension(nMaxGrp) :: vecDeutGroup
  integer i, nDeutedSide, nGroupSize, nDeutThis
  integer, dimension(:) :: WeightsSide, SpeciesGroup
  character(len=*), dimension(:) :: StrDeutedSide
  call makeBigSpecies (SpeciesGroup, nGroupSize, SpeciesEleBig)
  HydrGroup%nItem = getHAprCt(SpeciesEleBig)
  HydrGroup%nEle = 1
  if (allocated(HydrGroup%EleList)) deallocate(HydrGroup%EleList)
  if (allocated(HydrGroup%ArrCt)) deallocate(HydrGroup%ArrCt)
  allocate(HydrGroup%EleList(1))
  allocate(HydrGroup%ArrCt(HydrGroup%nItem,1))
  do i=1, nMaxGrp
    if (.NOT. allocated(vecDeutGroup(i)%EleList)) &
      allocate(vecDeutGroup(i)%EleList(2))
    if (allocated(vecDeutGroup(i)%ArrCt)) &
      deallocate(vecDeutGroup(i)%ArrCt)
    allocate(vecDeutGroup(i)%ArrCt(HydrGroup%nItem, 2))
  end do
  HydrGroup%EleList(1) = idxH
  call getHGroup (SpeciesEleBig, HydrGroup)
  call DeutGroups (HydrGroup, nDeutThis, vecDeutGroup, nDeutedSide)
  call DeutCleanGroups (vecDeutGroup, nDeutedSide, nDeutedSideCleaned)
  do i=1, nDeutedSide
    call DeutIns(vecDeutGroup(i), SpeciesEleBig, SpeciesEleD)
    StrDeutedSide(i) = ele2str(SpeciesEleD)
    WeightsSide(i) = vecDeutGroup(i)%nWeight
  end do
end subroutine DeutSide



subroutine makeBigSpecies (SpeciesGroup, nGroupSize, SpeciesEleBig)
  implicit none
  integer, dimension(:) :: SpeciesGroup
  integer nGroupSize, i
  type (EleSpecies) :: SpeciesEleBig
  SpeciesEleBig = SpeciesEleAll(SpeciesGroup(1))
  !Combine the EleSpecies structures of all the species in a group together.
  do i=2, nGroupSize
    SpeciesEleBig%nItem = SpeciesEleBig%nItem + 1
    SpeciesEleBig%vecIdEle(SpeciesEleBig%nItem) = idxDumb
    SpeciesEleBig%vecCtEle(SpeciesEleBig%nItem) = 1
    SpeciesEleBig%vecIdEle(SpeciesEleBig%nItem+1 : &
      SpeciesEleBig%nItem+SpeciesEleAll(SpeciesGroup(i))%nItem) = &
      SpeciesEleAll(SpeciesGroup(i))%vecIdEle(1:&
        SpeciesEleAll(SpeciesGroup(i))%nItem)
    SpeciesEleBig%vecCtEle(SpeciesEleBig%nItem+1 : &
      SpeciesEleBig%nItem+SpeciesEleAll(SpeciesGroup(i))%nItem) = &
      SpeciesEleAll(SpeciesGroup(i))%vecCtEle(1:&
        SpeciesEleAll(SpeciesGroup(i))%nItem)
    SpeciesEleBig%nItem = SpeciesEleBig%nItem + &
      SpeciesEleAll(SpeciesGroup(i))%nItem
  end do
end subroutine makeBigSpecies


subroutine SplitSideStr (SideStr, StrSplitted, nSpeSplitted)
  use CMDT
  implicit none
  character(len=lenStrSideMax) SideStr
  character(len=*), dimension(:) :: StrSplitted
  integer nSpeSplitted, i, j
  nSpeSplitted = 1
  j = 1
  do i = 1, len_trim(SideStr)
    if (SideStr(i:i) .NE. nameElements(idxDumb)) then
      StrSplitted(nSpeSplitted)(j:j) = SideStr(i:i)
      j = j + 1
    else
      nSpeSplitted = nSpeSplitted + 1
      j = 1
    end if
  end do
end subroutine SplitSideStr


function ele2str (EleSpe)
  use CMDT
  implicit none
  character(len=constLenNameSpecies) ele2str
  character(len=2) ctmp
  TYPE (EleSpecies) EleSpe
  integer i
  ele2str=""
  do i=1, EleSpe%nItem
    if (EleSpe%vecCtEle(i) .GT. 1) then
      write(ctmp, '(I2)') , EleSpe%vecCtEle(i)
      ele2str = trim(ele2str)//trim(nameElements(EleSpe%vecIdEle(i)))//&
        trim(ADJUSTL(ctmp))
    end if
    if (EleSpe%vecCtEle(i) .EQ. 1) then
      ele2str = trim(ele2str)//trim(nameElements(EleSpe%vecIdEle(i)))
    end if
  end do
end function ele2str

integer function getHAprCt (EleSpeciesA)
  implicit none
  integer i, getHAprCt
  TYPE(EleSpecies) EleSpeciesA
  getHAprCt = 0
  do i=1, EleSpeciesA%nItem
    if (EleSpeciesA%vecIdEle(i) .EQ. idxH) then
      getHAprCt = getHAprCt + 1
    end if
  end do
end function getHAprCt


subroutine getHGroup (EleSpeciesA, HydrGroup)
  use CMDT
  implicit none
  integer i, i1
  TYPE(EleSpecies) EleSpeciesA
  type (EleGroup) HydrGroup
  i1 = 0
  do i=1, EleSpeciesA%nItem
    if (EleSpeciesA%vecIdEle(i) .EQ. idxH) then
      i1 = i1 + 1
      HydrGroup%ArrCt(i1,1)  =  EleSpeciesA%vecCtEle(i)
    end if
  end do
end subroutine getHGroup


integer function binomcoeff (n, m)
  implicit none
  integer binomcoeff, n, m, i
  binomcoeff = 1
  do i=max(m,n-m)+1, n
    binomcoeff = binomcoeff * i
  end do
  do i=2, min(m,n-m)
    binomcoeff = binomcoeff / i
  end do
end function binomcoeff


integer function getNextColex (x)
! Gosper's hack.
  implicit none
  integer getNextColex, x, s, r
  s = IAND(x, -x)
  r = s + x
  getNextColex = IOR(r, ISHFT(IEOR(x,r), -2)/s)
end function getNextColex


subroutine intExchange(x, y)
! Or a even easier one: x=x+y; y=x-y; x=x-y
  implicit none
  integer x, y
  x = ieor(x, y)
  y = ieor(x, y)
  x = ieor(x, y)
end subroutine intExchange


integer function getAltAB(x, A, B)
  implicit none
  integer getAltAB, x, A, B
  getAltAB = A - x + B
  !getAltAB = IEOR(A, IEOR(B, x))
end function getAltAB



END MODULE
