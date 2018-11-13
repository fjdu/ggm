      MODULE CMDT
      implicit none
      ! Common data
      
      ! Physical constants. Which unit system to use is an issue.
      double precision, parameter :: SecondsPerYear = 3600D0*24D0*365D0
      double precision, parameter :: elementaryCharge = 1.602176487D-19
      double precision, parameter :: mProton = 1.67262158D-27 ! kg
      double precision, parameter :: CoulombConst = 8.9875517873681764D9
      double precision, parameter :: kBoltzmann = 1.3806503D-23
      double precision, parameter :: DPI = 3.1415926535897932384626433D0
      double precision, parameter :: hPlanck = 6.62606896D-34
      double precision, parameter :: hbarPlanck = 1.054571628D-34
      double precision, parameter :: GravitationConst = 6.67428D-11
      double precision, parameter :: SpeedOfLight = 299792458D0

      ! 2011-02-03 Thu 12:26:49
      ! Split +- into + and -.
      ! Possible issues of backward compatibility.
      integer, parameter :: nElement = 20
      character(LEN=8), dimension(nElement) :: &
        nameElements = &
          (/'+       ', '-       ', 'E       ', 'g       ', &
            'Grain   ', 'H       ', 'D       ', 'He      ', &
            'C       ', 'N       ', 'O       ', 'Si      ', &
            'S       ', 'Fe      ', 'Na      ', 'Mg      ', &
            'Cl      ', 'P       ', 'F       ', '@       '/)
      double precision, dimension(nElement), parameter :: &
        ElementMassNumber = &
          (/0D0       , 0D0       , 5.45D-4   , 0D0       , &
            0D0       , 1D0       , 2D0       , 4D0       , &
            12D0      , 14D0      , 16D0      , 28D0      , &
            32D0      , 56D0      , 23D0      , 24D0      , &
            35.5D0    , 31D0      , 19D0      , 0D0/)
      double precision, dimension(nElement), parameter :: &
        ElementTypicalAbundance = &
          (/0D0       , 0D0       , 5D-8      , 0D0       , &
            3D-12     , 1D0       , 2D-5      , 1.4D-1    , &
            7.3D-5    , 2.14D-5   , 1.76D-4   , 3.0D-9    , &
            2.0D-8    , 3.0D-9    , 3.0D-9    , 3.0D-9    , &
            3.0D-9    , 3.0D-9    , 2.0D-8    , 0D0/)


      integer, parameter :: constLenNameSpecies = 12, &
        nReactants = 3, nProducts = 4
      integer, parameter :: nParticipants = nReactants + nProducts

      double precision Temperature, n_H, &
        rateCosIon, Av, omega_Albedo, &
        MeanMolWeight, ratioGrainToH, ratioDustGasMass, &
        GrainRadius, GrainDensity, SitesDensity, &
        SitesPerGrain, VolumnPerGrain, BarrierWidth, &
        GrainRadiusBase, sto_threshold, rateThreshold, timeUnit, &
        stickCoeffChargeGrain
      integer nOrderLim

      double precision T300, TemperatureReduced, nHatom, NHtotal, &
        JNegaPosi, JChargeNeut

      character(LEN=128) path, &
        fReactions, fMobility, fDissRecom, &
        fInitialCondition, &
        fReactionsSave, fSpeciesSave, fRatesSave, &
        fSaveFinalResult, fSaveElementalAb, fSaveALLHistory, &
        fSaveALLHistory_ASCII, fAnalyseResult, fFinalJac, &
        fNameMoments, fEqMoments, fOutLog, &
        fPhyParHistory

      double precision ATOL, RTOL

      integer nReactions, nSpecies, nSpecies_Est, nInitialSpecies, &
        nMobilityList, nDissRecom, nIteration, nRecordsSave
      integer(kind=8) nMCSteps, nMCRepeat

      character(LEN=constLenNameSpecies), dimension(:,:), &
        allocatable :: strReactants, strProducts
      character(LEN=constLenNameSpecies), dimension(:), &
        allocatable :: nameSpecies
      character(LEN=128), dimension(:), &
        allocatable :: nameMoments
      
      double precision, dimension(:, :), allocatable :: &
        dblABC, ratesDissRecom
      double precision, dimension(:), allocatable :: &
        massSpecies, rates, propensities, mobilitySpecies
      double precision, dimension(:), allocatable :: t_Save
      
      logical, dimension(:), allocatable :: IsReactionChanged
      integer, dimension(:), allocatable :: idxSpeUpdated
      integer nSpeUpdated

      integer, dimension(:, :), allocatable :: &
        reactions, SpeciesElements, moments, &
        invIdx, invIdxM, &
        invIdxSpecies, invIdxSpeciesM, &
        idxReac_AsReactants, idxReac_AsProducts, &
        idxChangeSpe
      integer, dimension(:), allocatable :: &
        momentsType, n_invIdxSpecies, &
        nReac_AsReactants, nReac_AsProducts, &
        nChangeSpe
      integer, dimension(:), allocatable :: typeReac, nRealReactants, &
        nRealProducts
      double precision, dimension(:, :), allocatable :: PopSpeciesSave
      !double precision, dimension(:), allocatable :: PopSpecies
      integer(kind=8), dimension(:), allocatable :: PopSpecies
      logical, dimension(:, :), allocatable :: &
        sparseMaskJac, sparseMaskJac1
      integer nGrReactions, nGrSpecies
      integer, dimension(:), allocatable :: &
        idxGrReactions, idxGrSpecies
      logical, dimension(:), allocatable :: &
        IsStoSpecies, IsStoMoments, IsStoReac, IsEndProd, IsWatch

      integer nHistoryLen
      double precision, dimension(:), allocatable :: touts
      double precision, dimension(:,:), allocatable :: yHistory
      double precision, dimension(:,:), allocatable :: ydotHistory
      double precision, dimension(:,:), allocatable :: PhyParHistory
      end module CMDT
