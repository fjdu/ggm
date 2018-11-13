! Created on 2011-02-01 Tue 18:17:19


MODULE CMReac
implicit none

TYPE :: FunGroup
! Example:
!   For (CH2)3,
!     nSize = 2, nCopy=3
!     vecEle=[7, 4], vecCount=[1, 2]
  integer nSize, nCopy
  integer, dimension(:), allocatable :: vecEle, vecCount
END TYPE FunGroup


TYPE :: ChemSpecies
! Example:
!   For CH3(CH2)3OH,
!     nGroup = 5
!     vecGroup = [C, H3, (CH2)3, O, H]
  integer nGroup
  type (FunGroup), dimension(:), allocatable :: vecGroup
END TYPE ChemSpecies

!subroutine EquivSpecies (SpeciesA, SpeciesB)
!  use CMDT
!  implicit none
!  integer i
!  double precision dblTmp
!end subroutine EquivSpecies

END MODULE
