      module CMDT
      implicit none

      double precision, parameter :: yearSeconds = 3600D0 * 24D0 * 365D0
      double precision, parameter :: elementaryCharge = 1.602176487D-19
      double precision, parameter :: mProton = 1.67262158D-24 ! gram
      double precision, parameter :: coulombConst = 8.9875517873681764D9
      double precision, parameter :: kBoltzmann = 1.3806503D-23
      double precision, parameter :: DPI = 3.1415926535897932384626433D0

      integer, parameter :: nElement = 17
      character(LEN=8), dimension(nElement) :: &
        nameElements = &
          (/'+-      ', 'E       ', 'Grain   ', 'H       ', 'D       ', &
            'He      ', 'C       ', 'N       ', 'O       ', 'Si      ', &
            'S       ', 'Fe      ', 'Na      ', 'Mg      ', 'Cl      ', &
            'P       ', 'F       '/)
      ! Reformatted
      ! 2011-03-31 Thu 19:01:10

      integer, parameter :: constLenNameSpecies = 12, &
        nReactants = 3, nProducts = 4
      integer, parameter :: nParticipants = nReactants + nProducts

      double precision Temperature, n_H, &
        Av, omega_Albedo, ratioHDust, &
        rateHHH2, rateCosIon, &
        ratioDustGasMass, stickCoeffHH, stickCoeffChargeGrain, &
        GrainRadius, GrainDensity, GrainRadiusBase, &
        aGrainMin, aGrainMax
      double precision timeUnit, rateThreshold

      double precision T300, TemperatureReduced, nHatom, NHtotal, &
        JNegaPosi, JChargeNeut

      character(LEN=128) path, &
        fReactions, &
        fInitialCondition, &
        fReactionsSave, fSpeciesSave, fRatesSave, &
        fSaveFinalResult, fSaveElementalAb, fSaveALLHistory

      double precision ATOL, RTOL

      integer nReactions, nSpecies, nSpecies_Est, nInitialSpecies, &
        nDissRecom, nIteration

      character(LEN=constLenNameSpecies), dimension(:,:), &
        allocatable :: strReactants, strProducts
      character(LEN=constLenNameSpecies), dimension(:), &
        allocatable :: nameSpecies, nameInitialSpecies
      character(LEN=2), dimension(:), allocatable :: strType
      character(LEN=1), dimension(:), allocatable :: strQuality

      double precision, dimension(:, :), allocatable :: &
        dblABC, dblTLTU, ratesDissRecom
      double precision, dimension(:), allocatable :: rates

      integer, dimension(:, :), allocatable :: &
        reactions, SpeciesElements
      integer, dimension(:), allocatable :: iType, nRealReactants, &
        nRealProducts, idxInitialSpecies

      integer  nIterationBase, nTargetDecades, nerr, nHistoryLen
      double precision ratioTout, tStep
      double precision, dimension(:), allocatable :: &
        initialAbundances, touts
      double precision, dimension(:, :), allocatable :: yHistory
      double precision, dimension(:, :), allocatable :: ydotHistory

      INTERFACE
        SUBROUTINE setEnvPG (ColorIdx, LineWidth, LineStyle, &
          ifClip, CharSize, BGColor)
          INTEGER, INTENT (IN), OPTIONAL :: ColorIdx
          INTEGER, INTENT(IN), OPTIONAL :: LineWidth
          INTEGER, INTENT(IN), OPTIONAL :: LineStyle
          INTEGER, INTENT(IN), OPTIONAL :: ifClip
          REAL, INTENT(IN), OPTIONAL :: CharSize
          INTEGER, INTENT(IN), OPTIONAL :: BGColor
        END SUBROUTINE setEnvPG
      END INTERFACE

      end module CMDT
