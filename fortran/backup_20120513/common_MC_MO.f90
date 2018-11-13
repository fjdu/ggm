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
      double precision, parameter :: IdealGasConst = 8.314472D0
      double precision, parameter :: AvogadroConst = 6.02214179D23

      integer, parameter :: nElement = 17
      character(LEN=8), dimension(nElement) :: &
        nameElements = &
          (/'+-      ', 'E       ', 'Grain   ', 'H       ', &
            'D       ', 'He      ', 'C       ', 'N       ', &
            'O       ', 'Si      ', 'S       ', 'Fe      ', &
            'Na      ', 'Mg      ', 'Cl      ', 'P       ', &
            'F       '/)
      double precision, dimension(nElement), parameter :: &
        ElementMassNumber = &
          (/0D0,        5.45D-4,    0D0,        1D0,        &
            2D0,        4D0,        12D0,       14D0,       &
            16D0,       28D0,       32D0,       56D0,       &
            23D0,       24D0,       35.5D0,     31D0,       &
            19D0/)

      integer, parameter :: constLenNameSpecies = 12, &
        nReactants = 3, nProducts = 4
      integer, parameter :: nParticipants = nReactants + nProducts

      ! <timestamp>2011-05-25 Wed 23:32:21</timestamp>
      ! From Woodall et al. 2007
      double precision, parameter :: rateCosIon_base = 1.36D-17

      integer, parameter :: LongIntKind = selected_int_kind(13)

      double precision, parameter :: eta_step = 1D-2

      double precision T_gas, T_dust, n_H, &
        rateCosIon, Av, omega_Albedo, &
        MeanMolWeight, ratioGrainToH, ratioDustGasMass, &
        GrainRadius, GrainDensity, SitesDensity, &
        SitesPerGrain, VolumnPerGrain, ReactBarrierWidth, &
        DiffBarrierWidth, &
        GrainRadiusBase, sto_threshold, rateThreshold, &
        stickCoeffChargeGrain, &
        CosmicDesorpGrainT, CosmicDesorpPreFactor, &
        ChemDesorpPreFactor, &
        Diff2DesorRatio, tFinal
      integer nOrderLim
      logical useCosPhotoDiss, useZeroPointEnergy

      ! double precision T_gas_Reduced, nHatom, NHtotal, &
      !   JNegaPosi, JChargeNeut

      character(LEN=128) path, &
        fReactions, fDissRecom, &
        fInitialCondition, &
        fReactionsSave, fSpeciesSave, fRatesSave, &
        fSaveFinalResult, fSaveElementalAb, fSaveALLHistory, &
        fSaveALLHistory_ASCII, fAnalyseResult, fFinalJac, &
        fNameMoments, fEqMoments, fOutLog, &
        fPhyParHistory, fSpeciesEnthalpy, fSurfRadical

      double precision ATOL, RTOL

      integer nReactions, nSpecies, nSpecies_Est, nInitialSpecies, &
        nMobilityList, nDissRecom, nIteration, nRecordsSave

      character(LEN=constLenNameSpecies), dimension(:,:), &
        allocatable :: strReactants, strProducts
      character(LEN=constLenNameSpecies), dimension(:), &
        allocatable :: nameSpecies
      character(LEN=128), dimension(:), &
        allocatable :: nameMoments
      
      double precision, dimension(:, :), allocatable :: &
        dblABC, ratesDissRecom
      double precision, dimension(:), allocatable :: &
        massSpecies, rates, propensities, mobilitySpecies, &
        enthalpySpecies, BranchingRatios, E_desorp_species, &
        vibFreqSpecies, &
        tunneling_mass
      logical SimpleChemDesorption
      integer, parameter :: nPhyPar = 8, nAuxVar = 2
      
      logical, dimension(:), allocatable :: IsReactionChanged
      integer, dimension(:), allocatable :: idxSpeUpdated
      integer nMCRepeat, nSpeUpdated

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
      !double precision, dimension(:), allocatable :: PopSpecies
      integer(kind=8), dimension(:), allocatable :: PopSpecies
      logical, dimension(:, :), allocatable :: &
        sparseMaskJac
      integer nGrReactions, nGrSpecies, nSurfRadical
      integer, dimension(:), allocatable :: &
        idxGrReactions, idxGrSpecies, idxSurfRadical
      logical, dimension(:), allocatable :: &
        IsStoSpecies, IsStoMoments, IsStoReac, IsEndProd, IsWatch

      type :: ManSpe
        integer idxSurf, idxMant, idxGas, &
          idx_acc, idx_eva, idx_S2M, idx_M2S
      end type ManSpe
      integer nMantleSpecies
      integer, dimension(:), allocatable :: idxToMant
      logical, dimension(:), allocatable :: IsMan, IsSurfWithMan
      type (ManSpe), dimension(:), allocatable :: mantleSpe
      double precision totalMantleSpe, totalSurfSpe
      type :: ThreePhaseReac
        integer idx
        double precision gain_loss
      end type ThreePhaseReac
      type (ThreePhaseReac), dimension(:), allocatable :: &
        surf_mant_inc, surf_mant_dec
      integer n_surf_mant_inc, n_surf_mant_dec
      logical allowExposureByReaction, allowAlwaysExpose
      double precision tot_man_surf_ab, tot_man_mant_ab, &
        tot_man_acc_rate, tot_man_eva_rate, &
        tot_man_acc_rate_reac, tot_man_eva_rate_reac

      integer nThreeBodySurfReacs
      integer, dimension(:), allocatable :: ThreeBodySurfReacs
      double precision, dimension(:), allocatable :: rateThreeBodySurf
      logical useThreeBodyPrefac, useThreeBody

      integer, dimension(:), allocatable :: idxReacMM

      integer nHistoryLen
      double precision, dimension(:), allocatable :: touts
      double precision, dimension(:,:), allocatable :: yHistory
      double precision, dimension(:,:), allocatable :: ydotHistory
      double precision, dimension(:,:), allocatable :: PhyParHistory
      end module CMDT
