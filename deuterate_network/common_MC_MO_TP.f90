! Modification for third order
! 2010-10-26

MODULE CMTP
implicit none

! dimvecVar: max number of (different reactants/products)
! dimvecDeriv: max order of moments
! dimvecMon: max number of terms
! The long versions are only used as temporary variables.

integer, parameter :: dimvecVar=4, dimvecDeriv=8, dimvecMon=32, & !16
  dimvecMonlong=32, nameLen=1, nIterMoment=256
TYPE :: monomial
  integer coeff, nvecDeriv
  integer, dimension(dimvecVar) :: vecVar, vecExp
  integer, dimension(dimvecDeriv) :: vecDeriv
END TYPE monomial

TYPE :: moment
  integer coeff, nvecDeriv
  integer, dimension(dimvecDeriv) :: vecDeriv
END TYPE moment

TYPE :: polynomial
  integer nnonzero
  type (monomial), dimension(dimvecMon) :: vecMon
END TYPE polynomial

TYPE :: polynomiallong
  integer nnonzero
  type (monomial), dimension(dimvecMonlong) :: vecMon
END TYPE polynomiallong

TYPE :: polymoment
  integer nnonzero
  type (moment), dimension(dimvecMon) :: vecMoment
END TYPE polymoment

type :: dtMo
  double precision, dimension(:), allocatable :: coeff
end type dtMo

integer, dimension(:, :), allocatable :: MM
integer, dimension(:), allocatable :: nnzMM
integer, dimension(:, :, :), allocatable :: dtMM
integer, dimension(:), allocatable :: ndtMM
integer, dimension(dimvecDeriv) :: nAllMM
integer, dimension(:, :), allocatable :: indMM
integer, dimension(:, :), allocatable :: dauMM
integer, dimension(:), allocatable :: ndauMM
integer nNumMMEst, nNumMM, orderTmp
integer, dimension(2) :: nNewMM
logical, dimension(:), allocatable :: IsY
integer, dimension(:), allocatable :: indY, invIndY
logical, dimension(:), allocatable :: IsGas

type (polynomial), dimension(:), allocatable :: polyReac
type (polynomial) pderivTmp
TYPE (polymoment) pmomTmp
type (dtMo), dimension(:), allocatable :: dtMMMo

double precision, dimension(:), allocatable :: MomVals
integer, dimension(:), allocatable :: nTmMomVals
logical, dimension(:), allocatable :: IsYPre
double precision, dimension(:), allocatable :: TScale

logical flagMMNotUsed
integer, dimension(:), allocatable :: countMMNotUsed

CONTAINS

SUBROUTINE derivpolymulti (x, n, p, pderiv)
! p must be sorted.
  integer i
  integer, intent(in) :: n
  integer, dimension(n), intent(in) :: x
  type (polynomial), intent(in) :: p
  type (polynomial), intent(out) :: pderiv
  type (polynomial) ptmp
  ptmp = p
  do i=1, n
    call derivpoly(x(i), ptmp, pderiv)
    ptmp = pderiv
  end do
END SUBROUTINE derivpolymulti


SUBROUTINE derivpoly (x, p, pderiv)
  implicit none
  integer x, i, j, i1
  type (polynomial) p, pderiv
  type (polynomiallong) ptmp
  call initpoly(pderiv)
  call initpolylong(ptmp)
  i1 = 0
  do i=1, p%nnonzero
    do j=1, dimvecVar
      if ((p%vecMon(i)%vecVar(j) .EQ. x) .AND. &
        (p%vecMon(i)%vecExp(j) .NE. 0)) then
        i1 = i1 + 1
        ptmp%vecMon(i1) = p%vecMon(i)
        ptmp%vecMon(i1)%coeff = ptmp%vecMon(i1)%coeff * &
          ptmp%vecMon(i1)%vecExp(j)
        ptmp%vecMon(i1)%vecExp(j) = ptmp%vecMon(i1)%vecExp(j) - 1
        if (ptmp%vecMon(i1)%vecExp(j) .EQ. 0) then
          ptmp%vecMon(i1)%vecVar(j) = 0
        end if
      end if
    end do
    i1 = i1 + 1
    ptmp%vecMon(i1) = p%vecMon(i)
    ptmp%vecMon(i1)%nvecDeriv = ptmp%vecMon(i1)%nvecDeriv + 1
    ptmp%vecMon(i1)%vecDeriv(ptmp%vecMon(i1)%nvecDeriv) = x
  end do
  ptmp%nnonzero = i1
  call CombineLongpoly (ptmp)
  pderiv%nnonzero = ptmp%nnonzero
  pderiv%vecMon(1:ptmp%nnonzero) = ptmp%vecMon(1:ptmp%nnonzero)
END SUBROUTINE derivpoly



SUBROUTINE CombineLongpoly (p)
  implicit none
  integer i, j, k
  type (polynomiallong) p
  do i=1, p%nnonzero
    call sortmon(p%vecMon(i))
  end do
  do i=1, p%nnonzero
    !if (p%vecMon(i)%coeff .EQ. 0) cycle
    do j=i+1, p%nnonzero
      if (p%vecMon(j)%coeff .EQ. 0) cycle
      if (IsEqualMon(p%vecMon(i), p%vecMon(j))) then
        p%vecMon(i)%coeff = p%vecMon(i)%coeff + p%vecMon(j)%coeff
        p%vecMon(j)%coeff = 0
        p%vecMon(j)%vecVar = 0
        p%vecMon(j)%vecExp = 0
      end if
    end do
  end do
  call sortpolylong(p)
  call nnonzeropolylong(p)
END SUBROUTINE CombineLongpoly


logical function IsEqualMon (m1, m2)
! Assuming m1 and m2 are already sorted.
  implicit none
  integer i
  type (monomial) m1, m2
  logical :: IsEqualMon
  IsEqualMon = .TRUE.
  if (m1%nvecDeriv .NE. m2%nvecDeriv) then
    IsEqualMon = .FALSE.
    return
  end if
  do i=1, dimvecVar
    if ((m1%vecVar(i) .NE. m2%vecVar(i)) .OR. &
        (m1%vecExp(i) .NE. m2%vecExp(i))) then
      IsEqualMon = .FALSE.
      return
    end if
  end do
  do i=1, dimvecDeriv
    if (m1%vecDeriv(i) .NE. m2%vecDeriv(i)) then
      IsEqualMon = .FALSE.
      return
    end if
  end do
  return
END function IsEqualMon


SUBROUTINE sortpoly (p) ! Sort in desceding order; zeros are last
  implicit none
  integer i, j
  type (polynomial) p
  do i=1, dimvecMon
    do j=i+1, dimvecMon
      if (((p%vecMon(j)%coeff .GT. p%vecMon(i)%coeff) .AND. &
          (p%vecMon(j)%coeff .NE. 0)) .OR. &
          ((p%vecMon(j)%coeff .LT. p%vecMon(i)%coeff) .AND. &
          (p%vecMon(i)%coeff .EQ. 0))) then
        call SwapMon(p%vecMon(i), p%vecMon(j))
      end if
    end do
  end do
END SUBROUTINE sortpoly


SUBROUTINE sortpolylong (p) ! Sort in desceding order; zeros are last
! Only supposed to be used in CombineLongpoly because of the
!   p%nnonzero thing requires the nonzero monomials of p are all
!   in the first p%nnonzero positions.
  implicit none
  integer i, j
  type (polynomiallong) p
  do i=1, p%nnonzero
    do j=i+1, p%nnonzero
      if (((p%vecMon(j)%coeff .GT. p%vecMon(i)%coeff) .AND. &
          (p%vecMon(j)%coeff .NE. 0)) .OR. &
          ((p%vecMon(j)%coeff .LT. p%vecMon(i)%coeff) .AND. &
          (p%vecMon(i)%coeff .EQ. 0))) then
        call SwapMon(p%vecMon(i), p%vecMon(j))
      end if
    end do
  end do
END SUBROUTINE sortpolylong


SUBROUTINE sortpolymon (p) ! Sort in desceding order; zeros are last
  implicit none
  integer i, j
  type (polynomial) p
  do i=1, dimvecMon
    call sortmon(p%vecMon(i))
  end do
  do i=1, dimvecMon
    do j=i+1, dimvecMon
      if (((p%vecMon(j)%coeff .GT. p%vecMon(i)%coeff) .AND. &
          (p%vecMon(j)%coeff .NE. 0)) .OR. &
          ((p%vecMon(j)%coeff .LT. p%vecMon(i)%coeff) .AND. &
          (p%vecMon(i)%coeff .EQ. 0))) then
        call SwapMon(p%vecMon(i), p%vecMon(j))
      end if
    end do
  end do
END SUBROUTINE sortpolymon


SUBROUTINE sortmon (m)
  implicit none
  integer i, j
  type (monomial) m
  do i=1, dimvecVar
    do j=i+1, dimvecVar
      if (m%vecVar(j) .GT. m%vecVar(i)) then
        call SwapInt (m%vecVar(i), m%vecVar(j))
        call SwapInt (m%vecExp(i), m%vecExp(j))
      end if
    end do
  end do
  do i=1, dimvecDeriv
    do j=i+1, dimvecDeriv
      if (m%vecDeriv(j) .GT. m%vecDeriv(i)) then
        call SwapInt (m%vecDeriv(i), m%vecDeriv(j))
      end if
    end do
  end do
END SUBROUTINE sortmon


SUBROUTINE sortpolymoment (p) ! Sort in desceding order; zeros are last
  implicit none
  integer i, j
  type (polymoment) p
  do i=1, p%nnonzero
    do j=i+1, p%nnonzero
      if (((p%vecMoment(j)%coeff .GT. p%vecMoment(i)%coeff) .AND. &
          (p%vecMoment(j)%coeff .NE. 0)) .OR. &
          ((p%vecMoment(j)%coeff .LT. p%vecMoment(i)%coeff) .AND. &
          (p%vecMoment(i)%coeff .EQ. 0))) then
        call SwapMoment(p%vecMoment(i), p%vecMoment(j))
      end if
    end do
  end do
END SUBROUTINE sortpolymoment


SUBROUTINE nnonzeropoly (p) ! Assuming p is sorted.
  implicit none
  integer i
  type (polynomial) p
  do i=1, dimvecMon
    if (p%vecMon(i)%coeff .EQ. 0) exit
  end do
  p%nnonzero = i - 1
END SUBROUTINE nnonzeropoly


SUBROUTINE nnonzeropolymoment (p) ! Assuming p is sorted.
  implicit none
  integer i
  type (polymoment) p
  do i=1, dimvecMon
    if (p%vecMoment(i)%coeff .EQ. 0) exit
  end do
  p%nnonzero = i - 1
END SUBROUTINE nnonzeropolymoment


SUBROUTINE nnonzeropolylong (p) ! Assuming p is sorted.
  implicit none
  integer i
  type (polynomiallong) p
  do i=1, dimvecMonlong
    if (p%vecMon(i)%coeff .EQ. 0) exit
  end do
  p%nnonzero = i - 1
END SUBROUTINE nnonzeropolylong


SUBROUTINE initpoly (p)
  implicit none
  integer i
  type (polynomial) p
  p%nnonzero = 0
  do i=1, dimvecMon
    p%vecMon(i)%coeff = 0
    p%vecMon(i)%nvecDeriv = 0
    p%vecMon(i)%vecVar = 0
    p%vecMon(i)%vecExp = 0
    p%vecMon(i)%vecDeriv = 0
  end do
END SUBROUTINE initpoly


SUBROUTINE initpolylong (p)
  implicit none
  integer i
  type (polynomiallong) p
  p%nnonzero = 0
  do i=1, dimvecMonlong
    p%vecMon(i)%coeff = 0
    p%vecMon(i)%nvecDeriv = 0
    p%vecMon(i)%vecVar = 0
    p%vecMon(i)%vecExp = 0
    p%vecMon(i)%vecDeriv = 0
  end do
END SUBROUTINE initpolylong


SUBROUTINE initpolymoment (p)
  implicit none
  integer i
  type (polymoment) p
  p%nnonzero = 0
  do i=1, dimvecMon
    p%vecMoment(i)%coeff = 0
    p%vecMoment(i)%nvecDeriv = 0
    p%vecMoment(i)%vecDeriv = 0
  end do
END SUBROUTINE initpolymoment


SUBROUTINE SwapInt (i, j)
  implicit none
  integer i, j, k
  k = i
  i = j
  j = k
END SUBROUTINE SwapInt


SUBROUTINE SwapMon (m1, m2)
  implicit none
  type (monomial) m1, m2
  integer i
  call SwapInt(m1%coeff, m2%coeff)
  call SwapInt(m1%nvecDeriv, m2%nvecDeriv)
  do i=1, dimvecVar
    call SwapInt(m1%vecVar(i), m2%vecVar(i))
    call SwapInt(m1%vecExp(i), m2%vecExp(i))
  end do
  do i=1, dimvecDeriv
    call SwapInt(m1%vecDeriv(i), m2%vecDeriv(i))
  end do
END SUBROUTINE SwapMon


SUBROUTINE SwapMoment (m1, m2)
  implicit none
  type (moment) m1, m2
  integer i
  call SwapInt(m1%coeff, m2%coeff)
  call SwapInt(m1%nvecDeriv, m2%nvecDeriv)
  do i=1, dimvecDeriv
    call SwapInt(m1%vecDeriv(i), m2%vecDeriv(i))
  end do
END SUBROUTINE SwapMoment


logical function IsEqualMonDeriv (m1, m2)
  implicit none
  integer i
  type (monomial) m1, m2
  logical :: IsEqualMonDeriv
  IsEqualMonDeriv = .TRUE.
  if (m1%nvecDeriv .NE. m2%nvecDeriv) then
    IsEqualMonDeriv = .FALSE.
    return
  end if
  do i=1, dimvecDeriv
    if (m1%vecDeriv(i) .NE. m2%vecDeriv(i)) then
      IsEqualMonDeriv = .FALSE.
      return
    end if
  end do
  return
END function IsEqualMonDeriv


logical function IsEqualVec (v1, v2, n)
  implicit none
  integer i, n
  integer, dimension(:) :: v1, v2
  logical :: IsEqualVec
  IsEqualVec = .TRUE.
  do i=1, n
    if (v1(i) .NE. v2(i)) then
      IsEqualVec = .FALSE.
      return
    end if
  end do
  return
END function IsEqualVec


logical function IsEqualDeriv (D1, D2)
  implicit none
  integer, dimension(dimvecDeriv) :: D1, D2
  logical :: IsEqualDeriv
  IsEqualDeriv = IsEqualVec(D1, D2, dimvecDeriv)
  return
END function IsEqualDeriv


TYPE (polymoment) function polyMon2Moment (p)
  implicit none
  integer i, j, i1
  type (polynomial) p
  TYPE (polymoment) polyMon2Moment
  logical flag
  call initpolymoment(polyMon2Moment)
  if (p%nnonzero .LT. 1) then
    return
  end if
  i1 = 1
  polyMon2Moment%vecMoment(1)%coeff = p%vecMon(1)%coeff
  polyMon2Moment%vecMoment(1)%nvecDeriv = p%vecMon(1)%nvecDeriv
  polyMon2Moment%vecMoment(1)%vecDeriv = p%vecMon(1)%vecDeriv
  do i=2, p%nnonzero
    flag = .TRUE.
    do j=1, i1
      if (IsEqualDeriv(p%vecMon(i)%vecDeriv, &
        polyMon2Moment%vecMoment(j)%vecDeriv)) then
        polyMon2Moment%vecMoment(j)%coeff = &
          polyMon2Moment%vecMoment(j)%coeff + p%vecMon(i)%coeff
        flag = .FALSE.
        exit
      end if
    end do
    if (flag) then
      i1 = i1 + 1
      polyMon2Moment%vecMoment(i1)%coeff = p%vecMon(i)%coeff
      polyMon2Moment%vecMoment(i1)%nvecDeriv = p%vecMon(i)%nvecDeriv
      polyMon2Moment%vecMoment(i1)%vecDeriv = p%vecMon(i)%vecDeriv
    end if
  end do
  polyMon2Moment%nnonzero = i1
  call sortpolymoment(polyMon2Moment)
  call nnonzeropolymoment(polyMon2Moment)
  return
end function polyMon2Moment


SUBROUTINE monomial2str (m, StrNameVec, StrMonomial)
  implicit none
  integer i
  type (monomial) m
  character(len=*), dimension(:) :: StrNameVec
  character(len=*) StrMonomial
  character(len=1) CharTmp
  StrMonomial = ''
  write (StrMonomial, '(I3, "*(")') m%coeff
  StrMonomial = ADJUSTL(StrMonomial)
  do i=1, dimvecVar
    if (m%vecVar(i) .EQ. 0) exit
    write (CharTmp, '(I1)') m%vecExp(i)
    StrMonomial = trim(StrMonomial) // &
      ((StrNameVec(m%vecVar(i)))) // '^' // &
      CharTmp
  end do
  StrMonomial = trim(StrMonomial) // ')_'
  do i=1, m%nvecDeriv
    StrMonomial = trim(StrMonomial) // &
      trim(ADJUSTL(StrNameVec(m%vecDeriv(i))))
  end do
END SUBROUTINE monomial2str


SUBROUTINE polynomial2str (p, StrNameVec, StrPolynomial)
  implicit none
  integer i
  type (polynomial) p
  character(len=*), dimension(:) :: StrNameVec
  character(len=*), dimension(:) :: StrPolynomial
  do i=1, p%nnonzero
    call monomial2str (p%vecMon(i), StrNameVec, StrPolynomial(i))
  end do
END SUBROUTINE polynomial2str


SUBROUTINE moment2str (m, StrNameVec, StrMoment)
  implicit none
  integer i, i1
  type (moment) m
  character(len=*), dimension(:) :: StrNameVec
  character(len=*) StrMoment
  character(len=((nameLen+5)*dimvecDeriv+6)) StrTmp2
  StrMoment = ''
  StrTmp2 = ''
  if (m%coeff .EQ. 1) then
    StrMoment = '<'
  else if (m%coeff .EQ. -1) then
    StrMoment = '-<'
  else
    write (StrMoment, '(I3, "*<")') m%coeff
  end if
  StrMoment = ADJUSTL(StrMoment)
  StrMoment = trim(StrMoment) // &
    trim(ADJUSTL(StrNameVec(m%vecDeriv(1))))
  i1 = 0
  do i=2, m%nvecDeriv
    if (m%vecDeriv(i) .EQ. m%vecDeriv(i-1)) then
      i1 = i1 - 1
      write (StrTmp2, '(I3)') i1
      StrTmp2 = '(' // trim(adjustl(StrNameVec(m%vecDeriv(i)))) // &
        trim(adjustl(StrTmp2)) // ')'
    else
      i1 = 0
      StrTmp2 = trim(adjustl(StrNameVec(m%vecDeriv(i))))
    end if
    StrMoment = trim(StrMoment) // StrTmp2
  end do
  StrMoment = trim(StrMoment) // '>'
END SUBROUTINE moment2str


SUBROUTINE polymoment2str (p, StrNameVec, StrPolymoment)
  implicit none
  integer i
  type (polymoment) p
  character(len=*), dimension(:) :: StrNameVec, StrPolymoment
  do i=1, p%nnonzero
    call moment2str (p%vecMoment(i), StrNameVec, StrPolymoment(i))
  end do
END SUBROUTINE polymoment2str


SUBROUTINE vec2str (v, n, StrNameVec, StrOut)
  implicit none
  integer i, i1, n
  integer, dimension(n) :: v
  character(len=*), dimension(:) :: StrNameVec
  character(len=*) StrOut
  character(len=128) StrTmp2
  StrOut = ''
  if (n .LT. 1) return
  StrTmp2 = ''
  StrOut = '<' // trim(StrNameVec(v(1)))
  i1 = 0
  do i=2, n
    if (v(i) .EQ. v(i-1)) then
      i1 = i1 - 1
      write (StrTmp2, '(I3)') i1
      StrTmp2 = '(' // trim(adjustl(StrNameVec(v(i)))) // &
        trim(adjustl(StrTmp2)) // ')'
    else
      i1 = 0
      StrTmp2 = trim(adjustl(StrNameVec(v(i))))
    end if
    StrOut = trim(StrOut) // '.' // StrTmp2
  end do
  StrOut = trim(StrOut) // '>'
end SUBROUTINE vec2str


SUBROUTINE vec2strbasic (v, n, StrNameVec, StrOut)
  implicit none
  integer i, n
  integer, dimension(n) :: v
  character(len=*), dimension(:) :: StrNameVec
  character(len=*) StrOut
  StrOut = ''
  if (n .LT. 1) return
  StrOut = StrNameVec(v(1))
  do i=2, n
    StrOut = trim(StrOut) // '.' // StrNameVec(v(i))
  end do
end SUBROUTINE vec2strbasic


integer function getFstGasSpe (n)
  implicit none
  integer n, getFstGasSpe
  do getFstGasSpe=1, nnzMM(n)
    if (IsGas(MM(getFstGasSpe, n))) then
      return
    end if
  end do
  getFstGasSpe = 0
  return
end function getFstGasSpe


integer function getFstDetSpe1 (n)
  use CMDT
  implicit none
  integer n, getFstDetSpe1
  do getFstDetSpe1=1, nnzMM(n)
    if (IsGas(MM(getFstDetSpe1, n)) .OR. &
        (.NOT. IsStoSpecies(MM(getFstDetSpe1, n)))) then
      return
    end if
  end do
  getFstDetSpe1 = 0
  return
end function getFstDetSpe1


integer function getFstDetSpe2 (n, y, NEQ)
  use CMDT
  implicit none
  integer n, getFstDetSpe2
  integer, intent(in) :: NEQ
  double precision ,dimension(NEQ), intent(in) :: y
  do getFstDetSpe2=1, nnzMM(n)
    if (IsGas(MM(getFstDetSpe2, n)) .OR. &
        (.NOT. IsStoSpecies(MM(getFstDetSpe2, n))) .OR. &
        (y(MM(getFstDetSpe2, n)) .GE. 2D0*sto_threshold)) then
      return
    end if
  end do
  getFstDetSpe2 = 0
  return
end function getFstDetSpe2



integer function getFstDetSpe (n, y, NEQ)
  use CMDT
  implicit none
  integer n, getFstDetSpe
  integer, intent(in) :: NEQ
  double precision ,dimension(NEQ), intent(in) :: y
  do getFstDetSpe=1, nnzMM(n)
    if (IsGas(MM(getFstDetSpe, n)) .OR. &
        (.NOT. IsStoSpecies(MM(getFstDetSpe, n)))) then
        !(y(MM(getFstDetSpe, n)) .GE. sto_threshold)) then
      return
    end if
  end do
  getFstDetSpe = 0
  return
end function getFstDetSpe



double precision function GetMinComp (n, y, NEQ)
  use CMDT
  implicit none
  integer n, i
  double precision GetMinComp
  integer, intent(in) :: NEQ
  double precision ,dimension(NEQ), intent(in) :: y
  GetMinComp = y(MM(1, n))
  do i=2, nnzMM(n)
    if (y(MM(i, n)) .LT. GetMinComp) then
      GetMinComp = y(MM(i, n))
    end if
  end do
  return
end function GetMinComp



function getdau(v, n, ndrop)
  implicit none
  integer n, ndrop
  integer, dimension(n) :: v
  integer, dimension(n-1) :: getdau
  getdau(1:(ndrop-1)) = v(1:(ndrop-1))
  if (ndrop .LE. (n-1)) getdau(ndrop:(n-1)) = v((ndrop+1):n)
  return
end function getdau


integer function getdauDrop(n, ndrop)
  implicit none
  integer n, ndrop, n1, i
  integer, dimension(dimvecDeriv) :: v
  integer getdauDrop
  n1 = nnzMM(n)
  v(1:(n1-1)) = getdau(MM(1:n1, n), n1, ndrop)
  do i=1, ndauMM(n)
    getdauDrop = dauMM(i, n)
    if (IsEqualVec(v, MM(1:(n1-1), getdauDrop), n1-1)) then
      return
    end if
  end do
  getdauDrop = 0
  return
end function getdauDrop


double precision function getMMProd (n, y, NEQ)
  implicit none
  integer i, n, NEQ
  double precision getMMProd
  double precision, dimension(NEQ) :: y
  getMMProd = y(MM(1, n))
  do i=2, nnzMM(n)
    getMMProd = getMMProd * y(MM(i, n))
  end do
end function getMMProd


recursive subroutine makedau (v, n, ind, flagReNew)
  implicit none
  integer n, n1, i, j, i1, ind, nNumMMTmp
  integer, dimension(n) :: v
  integer, dimension(n-1) :: tmpM
  logical flag, flagReNew
  if (n .LE. 1) return
  n1 = n - 1 ! A 1-moment does not have a dau.
  do i=1, n
    tmpM = getdau(v, n, i)
    flag = .TRUE.
    do j=1, ndauMM(ind)
      i1 = dauMM(j, ind)
      if (nnzMM(i1) .EQ. n1) then
        if (IsEqualVec(tmpM, MM(1:n1, i1), n1)) then
          flag = .FALSE. ! Not a new dau.
          exit
        end if
      end if
    end do
    if (flag) then
      do j=1, nNumMM
        if (nnzMM(j) .EQ. n1) then
          if (IsEqualVec(tmpM, MM(1:n1, j), n1)) then
            flag = .FALSE. ! A dau already exist.
            ndauMM(ind) = ndauMM(ind) + 1
            dauMM(ndauMM(ind), ind) = j
            exit
          end if
        end if
      end do
      if (flag) then ! A new-born dau.
        nNumMM = nNumMM + 1
        MM(1:n1, nNumMM) = tmpM
        nnzMM(nNumMM) = n1
        nAllMM(n1) = nAllMM(n1) + 1
        indMM(nAllMM(n1), n1) = nNumMM
        ndauMM(ind) = ndauMM(ind) + 1
        dauMM(ndauMM(ind), ind) = nNumMM
        if (flagReNew) then
          nNewMM(1) = nNumMM
          flagReNew = .FALSE.
        else
          nNewMM(2)  = nNumMM
        end if
        if (n1 .GT. 1) then
          nNumMMTmp = nNumMM
          call makedau (tmpM, n1, nNumMMTmp, flagReNew)
        end if
      end if
    end if
  end do
end subroutine makedau


recursive function getMomValue(n, m, y, NEQ) result (MomValue)
  use CMDT
  IMPLICIT NONE
  double precision MomValue
  integer, intent(in) :: n, m, NEQ
  integer i, i1
  double precision, dimension(NEQ), intent(in) :: y
  if (nnzMM(n) .EQ. 1) then
    MomValue = y(n)
    return
  end if
  if (IsY(n)) then
    if (n .EQ. indY(m)) then
      MomValue = y(m)
    else
      !MomValue = min(y(invIndY(n)), GetMinComp(n, y, NEQ))
      MomValue = min(y(invIndY(n)), getMMProd(n, y, NEQ))
      !if (MomValue .LT. y(invIndY(n))) then
      !  flagMMNotUsed = .TRUE.
      !  countMMNotUsed(n) = countMMNotUsed(n) + 1
      !else
      !  MomValue = y(invIndY(n))
      !end if
    end if
  else
    i = getFstDetSpe(n, y, NEQ)
    if (i .GT. 0) then
      i1 = getdauDrop(n, i)
      MomValue = y(MM(i, n)) * getMomValue(i1, m, y, NEQ)
    else
      if (nnzMM(n) .GT. nOrderLim) then
        MomValue = 0D0
      else
        MomValue = getMMProd(n, y, NEQ)
      end if
    end if
  end if
  return
end function getMomValue



double precision function derivDet (ni, nj, y, NEQ) ! nj should be
  use CMDT    ! single.
  IMPLICIT NONE
  integer i, i1, ni, nj, NEQ
  double precision, dimension(NEQ) :: y
  double precision derivDet
  derivDet = 0D0
  i1 = nnzMM(ni)
  do i=1, i1
    if (MM(i, ni) .EQ. nj) then
      derivDet = derivDet + &
        getVecProd(getdau(MM(1:i1, ni), i1, i), i1-1, y, NEQ)
    end if
  end do
  return
end function derivDet


double precision function getVecProd (v, n, y, NEQ)
  use CMDT ! The elements of v should be no larger  than
  IMPLICIT NONE   ! nSpecies.
  integer i, i1, n, NEQ
  double precision, dimension(NEQ) :: y
  integer, dimension(n) :: v
  double precision getVecProd
  if (n .EQ. 0) then
    getVecProd = 0D0
    return
  end if
  getVecProd = y(v(1))
  do i=2, n
    getVecProd = getVecProd * y(v(i))
  end do
  return
end function getVecProd


logical function getSparse(ni, nj, y, NEQ)
  use CMDT
  IMPLICIT NONE
  logical getSparse
  integer i, i1, ni, nj, NEQ
  double precision, dimension(NEQ), intent(in) :: y
  double precision rtmp
  do i=1, ndtMM(ni)
    i1 = dtMM(1, i, ni)
    if (IfContain(i1, nj)) then
      getSparse = .TRUE.
      return
    end if
  end do
  getSparse = .FALSE.
  return
end function getSparse



recursive logical function IfContain(ni, nj) result (Contain)
  use CMDT
  IMPLICIT NONE
  logical Contain
  integer i, i1, ni, nj
  if (nnzMM(ni) .LT. nnzMM(nj)) then
    Contain = .FALSE.
    return
  end if
  if (nnzMM(ni) .EQ. nnzMM(nj)) then
    if (ni .EQ. nj) then
      Contain = .TRUE.
      return
    else
      Contain = .FALSE.
      return
    end if
  end if
  Contain = .FALSE.
  do i=1, ndauMM(ni)
    Contain = (Contain .OR. IfContain(dauMM(i, ni), nj))
    if (Contain) return
  end do
  return
end function IfContain


double precision function getMaxYMM (n, y, NEQ)
  use CMDT
  IMPLICIT NONE
  integer, intent(in) :: n, NEQ
  double precision, dimension(NEQ), intent(in) :: y
  double precision getMaxYMM
  integer i
  getMaxYMM = y(MM(1, n))
  do i=2, nnzMM(n)
    if (getMaxYMM .LT. y(MM(i, n))) getMaxYMM = y(MM(i, n))
  end do
  return
end function getMaxYMM

END MODULE CMTP
