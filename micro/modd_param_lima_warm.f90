!MNH_LIC Copyright 2013-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ###########################
      MODULE MODD_PARAM_LIMA_WARM
!     ###########################
!
!!****  *MODD_PARAM_LIMA_WARM* - declaration of some descriptive parameters and
!!                               microphysical factors extensively used in
!!                               the LIMA warm scheme.
!!    AUTHOR
!!    ------
!!      J.-P. Pinty      * Laboratoire d'Aerologie*
!!      S.    Berthet    * Laboratoire d'Aerologie*
!!      B.    Vié        * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original             ??/??/13
!!
!-------------------------------------------------------------------------------
USE MODD_PARAMETERS, ONLY: JPSVNAMELGTMAX
!
IMPLICIT NONE
!
!*       1.   DESCRIPTIVE PARAMETERS
!             ----------------------
!
TYPE PARAM_LIMA_WARM_t
REAL      ::  XLBC, XLBEXC,          & ! shape parameters of the cloud droplets
              XLBR, XLBEXR, XNR        ! shape parameters of the raindrops
!
REAL      :: XAR,XBR,XCR,XDR,XF0R,XF1R,     & ! Raindrop       charact.
                             XCCR,XCXR,     & !For diagnostics
             XAC,XBC,XCC,XDC,XF0C,XF2C,XC1C   ! Cloud droplet  charact.
!
!
!
!-------------------------------------------------------------------------------
!
!*       2.   MICROPHYSICAL FACTORS
!             ---------------------
!
REAL      :: XFSEDRR,XFSEDCR,                  & ! Constants for sedimentation
             XFSEDRC,XFSEDCC                     ! fluxes of R, C
!
!
REAL      :: XDIVA,                            & ! Diffusivity of water vapor
             XTHCO                               ! Thermal conductivity
REAL      :: XWMIN                               ! Min value of updraft velocity
                                                 ! to enable nucleation process
REAL      :: XTMIN                               ! Min value of
                                                 ! temperature evolution
                                                 ! to enable nucleation process
REAL      :: XCSTHEN,XCSTDCRIT                   ! Cst for HEN precalculations
INTEGER   :: NHYP                                ! Number of value of the HYP
                                                 !    functions
REAL      :: XHYPINTP1, XHYPINTP2                ! Factors defining the
                                                 ! supersaturation log scale
REAL, DIMENSION(:,:), ALLOCATABLE              & ! Tabulated HYPgeometric
          :: XHYPF12, XHYPF32                    !   functions used in HEN
INTEGER   :: NAHEN                               ! Number of value of the AHEN
                                                 !    functions
REAL      :: XAHENINTP1, XAHENINTP2              ! Factors defining the
                                                 ! temperatures in lin scale
REAL, DIMENSION(:), ALLOCATABLE                & !
          :: XAHENG,XAHENG2,XAHENG3,XPSI1, XPSI3,      & ! Twomey-CPB98 and
             XAHENF,XAHENY                       ! Feingold-Heymsfield
                                                 ! parameterization to compute Smax
REAL      :: XWCOEF_F1, XWCOEF_F2, XWCOEF_F3,  & ! COEF_F of the polynomial temp.
             XWCOEF_Y1, XWCOEF_Y2, XWCOEF_Y3     ! COEF_Y of the polynomial temp.
                                                 ! function powering W
!
!
REAL      :: XKERA1, XKERA2                      ! Constants to define the lin
                                                 ! and parabolic kernel param.
REAL      :: XSELFC                              ! Constants for cloud droplet
                                                 ! selfcollection : SELF
!
REAL      :: XAUTO1, XAUTO2, XCAUTR,           & ! Constants for cloud droplet
                 XLAUTR, XLAUTR_THRESHOLD,         & ! autoconversion : AUT
                 XITAUTR, XITAUTR_THRESHOLD, XR0     ! XR0 for KHKO autoconversion
!
REAL      :: XACCR1, XACCR2, XACCR3,           & ! Constants for the accretion
             XACCR4, XACCR5, XACCR6,           & ! process
             XACCR_CLARGE1, XACCR_CLARGE2, XACCR_RLARGE1, XACCR_RLARGE2, &
             XACCR_CSMALL1, XACCR_CSMALL2, XACCR_RSMALL1, XACCR_RSMALL2, &
             XFCACCR, XEXCACCR
!
REAL      :: XSCBU2, XSCBU3,                   & ! Constants for the raindrop
             XSCBU_EFF1, XSCBU_EFF2, XSCBUEXP1   ! breakup-selfcollection: SCBU
!
REAL      :: XSPONBUD1,XSPONBUD2,XSPONBUD3,    & ! Spontaneous Break-up
             XSPONCOEF2                          ! (drop size limiter)
!
REAL      :: X0EVAR, X1EVAR,                   & ! Constants for raindrop
             XEX0EVAR, XEX1EVAR, XEX2EVAR,     & ! evaporation: EVA
             XCEVAP                              ! for KHKO
!
REAL,DIMENSION(:,:,:,:), ALLOCATABLE :: XCONCC_INI
REAL                                 :: XCONCR_PARAM_INI
                                      ! Used to initialize the
                                      ! concentrations from mixing ratios
                                      ! (init and grid-nesting from Kessler)
!
REAL      :: X0CNDC, X2CNDC                   ! Constants for cloud droplet
                                              ! condensation/evaporation
REAL      :: XFREFFC  ! Factor to compute the cloud droplet effective radius
REAL      :: XFREFFR  ! Factor to compute the rain drop     effective radius
REAL      :: XCREC, XCRER
                      ! Factors to compute reff when cloud and rain are present
END TYPE PARAM_LIMA_WARM_t
!
TYPE(PARAM_LIMA_WARM_t), TARGET, SAVE :: PARAM_LIMA_WARM
!
REAL, POINTER :: XLBC => NULL(), &
                 XLBEXC => NULL(), &
                 XLBR => NULL(), &
                 XLBEXR => NULL(), &
                 XNR => NULL(), &
                 XAR => NULL(), &
                 XBR => NULL(), &
                 XCR => NULL(), &
                 XDR => NULL(), &
                 XF0R => NULL(), &
                 XF1R => NULL(), &
                 XCCR => NULL(), &
                 XCXR => NULL(), &
                 XAC => NULL(), &
                 XBC => NULL(), &
                 XCC => NULL(), &
                 XDC => NULL(), &
                 XF0C => NULL(), &
                 XF2C => NULL(), &
                 XC1C => NULL(), &
                 XFSEDRR => NULL(), &
                 XFSEDCR => NULL(), &
                 XFSEDRC => NULL(), &
                 XFSEDCC => NULL(), &
                 XDIVA => NULL(), &
                 XTHCO => NULL(), &
                 XWMIN => NULL(), &
                 XTMIN => NULL(), &
                 XCSTHEN => NULL(), &
                 XCSTDCRIT => NULL(), &
                 XHYPINTP1 => NULL(), &
                 XHYPINTP2 => NULL(), &
                 XAHENINTP1 => NULL(), &
                 XAHENINTP2 => NULL(), &
                 XWCOEF_F1 => NULL(), &
                 XWCOEF_F2 => NULL(), &
                 XWCOEF_F3 => NULL(), &
                 XWCOEF_Y1 => NULL(), &
                 XWCOEF_Y2 => NULL(), &
                 XWCOEF_Y3 => NULL(), &
                 XKERA1 => NULL(), &
                 XKERA2 => NULL(), &
                 XSELFC => NULL(), &
                 XAUTO1 => NULL(), &
                 XAUTO2 => NULL(), &
                 XCAUTR => NULL(), &
                 XLAUTR => NULL(), &
                 XLAUTR_THRESHOLD => NULL(), &
                 XITAUTR => NULL(), &
                 XITAUTR_THRESHOLD => NULL(), &
                 XR0 => NULL(), &
                 XACCR1 => NULL(), &
                 XACCR2 => NULL(), &
                 XACCR3 => NULL(), &
                 XACCR4 => NULL(), &
                 XACCR5 => NULL(), &
                 XACCR6 => NULL(), &
                 XACCR_CLARGE1 => NULL(), &
                 XACCR_CLARGE2 => NULL(), &
                 XACCR_RLARGE1 => NULL(), &
                 XACCR_RLARGE2 => NULL(), &
                 XACCR_CSMALL1 => NULL(), &
                 XACCR_CSMALL2 => NULL(), &
                 XACCR_RSMALL1 => NULL(), &
                 XACCR_RSMALL2 => NULL(), &
                 XFCACCR => NULL(), &
                 XEXCACCR => NULL(), &
                 XSCBU2 => NULL(), &
                 XSCBU3 => NULL(), &
                 XSCBU_EFF1 => NULL(), &
                 XSCBU_EFF2 => NULL(), &
                 XSCBUEXP1 => NULL(), &
                 XSPONBUD1 => NULL(), &
                 XSPONBUD2 => NULL(), &
                 XSPONBUD3 => NULL(), &
                 XSPONCOEF2 => NULL(), &
                 X0EVAR => NULL(), &
                 X1EVAR => NULL(), &
                 XEX0EVAR => NULL(), &
                 XEX1EVAR => NULL(), &
                 XEX2EVAR => NULL(), &
                 XCEVAP => NULL(), &
                 XCONCR_PARAM_INI => NULL(), &
                 X0CNDC => NULL(), &
                 X2CNDC => NULL(), &
                 XFREFFC => NULL(), &
                 XFREFFR => NULL(), &
                 XCREC => NULL(), &
                 XCRER => NULL()

INTEGER, POINTER :: NHYP => NULL(), &
                    NAHEN => NULL()

REAL, DIMENSION(:,:), POINTER :: XHYPF12 => NULL(), &
                                 XHYPF32 => NULL()
REAL, DIMENSION(:), POINTER :: XAHENG => NULL(), &
                               XAHENG2 => NULL(), &
                               XAHENG3 => NULL(), &
                               XPSI1 => NULL(), &
                               XPSI3 => NULL(), &
                               XAHENF => NULL(), &
                               XAHENY => NULL()
REAL,DIMENSION(:,:,:,:), POINTER :: XCONCC_INI => NULL()
!
CHARACTER(LEN=JPSVNAMELGTMAX),DIMENSION(5),PARAMETER &
                     :: CLIMA_WARM_NAMES=(/'CCLOUD  ','CRAIN   ','CCCNFREE','CCCNACTI','SPRO    '/)
                                       ! basenames of the SV articles stored
                                       ! in the binary files
CHARACTER(LEN=JPSVNAMELGTMAX),DIMENSION(5),PARAMETER &
                     :: CLIMA_WARM_CONC=(/'NC   ','NR   ','NFREE','NCCN ','SS   '/)
!                                       ! basenames of the SV articles stored
!                                       ! in the binary files for DIAG
!
!* Special issue for Below-Cloud SCAVenging of Aerosol particles
CHARACTER(LEN=JPSVNAMELGTMAX),DIMENSION(2),PARAMETER :: CAERO_MASS =(/'MASSAP', 'MAP   '/)
!
CONTAINS
SUBROUTINE PARAM_LIMA_WARM_ASSOCIATE()
IMPLICIT NONE
IF(.NOT. ASSOCIATED(XLBC)) THEN
  XLBC              => PARAM_LIMA_WARM%XLBC
  XLBEXC            => PARAM_LIMA_WARM%XLBEXC
  XLBR              => PARAM_LIMA_WARM%XLBR
  XLBEXR            => PARAM_LIMA_WARM%XLBEXR
  XNR               => PARAM_LIMA_WARM%XNR
  XAR               => PARAM_LIMA_WARM%XAR
  XBR               => PARAM_LIMA_WARM%XBR
  XCR               => PARAM_LIMA_WARM%XCR
  XDR               => PARAM_LIMA_WARM%XDR
  XF0R              => PARAM_LIMA_WARM%XF0R
  XF1R              => PARAM_LIMA_WARM%XF1R
  XCCR              => PARAM_LIMA_WARM%XCCR
  XCXR              => PARAM_LIMA_WARM%XCXR
  XAC               => PARAM_LIMA_WARM%XAC
  XBC               => PARAM_LIMA_WARM%XBC
  XCC               => PARAM_LIMA_WARM%XCC
  XDC               => PARAM_LIMA_WARM%XDC
  XF0C              => PARAM_LIMA_WARM%XF0C
  XF2C              => PARAM_LIMA_WARM%XF2C
  XC1C              => PARAM_LIMA_WARM%XC1C
  XFSEDRR           => PARAM_LIMA_WARM%XFSEDRR
  XFSEDCR           => PARAM_LIMA_WARM%XFSEDCR
  XFSEDRC           => PARAM_LIMA_WARM%XFSEDRC
  XFSEDCC           => PARAM_LIMA_WARM%XFSEDCC
  XDIVA             => PARAM_LIMA_WARM%XDIVA
  XTHCO             => PARAM_LIMA_WARM%XTHCO
  XWMIN             => PARAM_LIMA_WARM%XWMIN
  XTMIN             => PARAM_LIMA_WARM%XTMIN
  XCSTHEN           => PARAM_LIMA_WARM%XCSTHEN
  XCSTDCRIT         => PARAM_LIMA_WARM%XCSTDCRIT
  XHYPINTP1         => PARAM_LIMA_WARM%XHYPINTP1
  XHYPINTP2         => PARAM_LIMA_WARM%XHYPINTP2
  XAHENINTP1        => PARAM_LIMA_WARM%XAHENINTP1
  XAHENINTP2        => PARAM_LIMA_WARM%XAHENINTP2
  XWCOEF_F1         => PARAM_LIMA_WARM%XWCOEF_F1
  XWCOEF_F2         => PARAM_LIMA_WARM%XWCOEF_F2
  XWCOEF_F3         => PARAM_LIMA_WARM%XWCOEF_F3
  XWCOEF_Y1         => PARAM_LIMA_WARM%XWCOEF_Y1
  XWCOEF_Y2         => PARAM_LIMA_WARM%XWCOEF_Y2
  XWCOEF_Y3         => PARAM_LIMA_WARM%XWCOEF_Y3
  XKERA1            => PARAM_LIMA_WARM%XKERA1
  XKERA2            => PARAM_LIMA_WARM%XKERA2
  XSELFC            => PARAM_LIMA_WARM%XSELFC
  XAUTO1            => PARAM_LIMA_WARM%XAUTO1
  XAUTO2            => PARAM_LIMA_WARM%XAUTO2
  XCAUTR            => PARAM_LIMA_WARM%XCAUTR
  XLAUTR            => PARAM_LIMA_WARM%XLAUTR
  XLAUTR_THRESHOLD  => PARAM_LIMA_WARM%XLAUTR_THRESHOLD
  XITAUTR           => PARAM_LIMA_WARM%XITAUTR
  XITAUTR_THRESHOLD => PARAM_LIMA_WARM%XITAUTR_THRESHOLD
  XR0               => PARAM_LIMA_WARM%XR0
  XACCR1            => PARAM_LIMA_WARM%XACCR1
  XACCR2            => PARAM_LIMA_WARM%XACCR2
  XACCR3            => PARAM_LIMA_WARM%XACCR3
  XACCR4            => PARAM_LIMA_WARM%XACCR4
  XACCR5            => PARAM_LIMA_WARM%XACCR5
  XACCR6            => PARAM_LIMA_WARM%XACCR6
  XACCR_CLARGE1     => PARAM_LIMA_WARM%XACCR_CLARGE1
  XACCR_CLARGE2     => PARAM_LIMA_WARM%XACCR_CLARGE2
  XACCR_RLARGE1     => PARAM_LIMA_WARM%XACCR_RLARGE1
  XACCR_RLARGE2     => PARAM_LIMA_WARM%XACCR_RLARGE2
  XACCR_CSMALL1     => PARAM_LIMA_WARM%XACCR_CSMALL1
  XACCR_CSMALL2     => PARAM_LIMA_WARM%XACCR_CSMALL2
  XACCR_RSMALL1     => PARAM_LIMA_WARM%XACCR_RSMALL1
  XACCR_RSMALL2     => PARAM_LIMA_WARM%XACCR_RSMALL2
  XFCACCR           => PARAM_LIMA_WARM%XFCACCR
  XEXCACCR          => PARAM_LIMA_WARM%XEXCACCR
  XSCBU2            => PARAM_LIMA_WARM%XSCBU2
  XSCBU3            => PARAM_LIMA_WARM%XSCBU3
  XSCBU_EFF1        => PARAM_LIMA_WARM%XSCBU_EFF1
  XSCBU_EFF2        => PARAM_LIMA_WARM%XSCBU_EFF2
  XSCBUEXP1         => PARAM_LIMA_WARM%XSCBUEXP1
  XSPONBUD1         => PARAM_LIMA_WARM%XSPONBUD1
  XSPONBUD2         => PARAM_LIMA_WARM%XSPONBUD2
  XSPONBUD3         => PARAM_LIMA_WARM%XSPONBUD3
  XSPONCOEF2        => PARAM_LIMA_WARM%XSPONCOEF2
  X0EVAR            => PARAM_LIMA_WARM%X0EVAR
  X1EVAR            => PARAM_LIMA_WARM%X1EVAR
  XEX0EVAR          => PARAM_LIMA_WARM%XEX0EVAR
  XEX1EVAR          => PARAM_LIMA_WARM%XEX1EVAR
  XEX2EVAR          => PARAM_LIMA_WARM%XEX2EVAR
  XCEVAP            => PARAM_LIMA_WARM%XCEVAP
  XCONCR_PARAM_INI  => PARAM_LIMA_WARM%XCONCR_PARAM_INI
  X0CNDC            => PARAM_LIMA_WARM%X0CNDC
  X2CNDC            => PARAM_LIMA_WARM%X2CNDC
  XFREFFC           => PARAM_LIMA_WARM%XFREFFC
  XFREFFR           => PARAM_LIMA_WARM%XFREFFR
  XCREC             => PARAM_LIMA_WARM%XCREC
  XCRER             => PARAM_LIMA_WARM%XCRER

  NHYP              => PARAM_LIMA_WARM%NHYP
  NAHEN             => PARAM_LIMA_WARM%NAHEN
ENDIF
END SUBROUTINE PARAM_LIMA_WARM_ASSOCIATE
!
SUBROUTINE PARAM_LIMA_WARM_ALLOCATE(HNAME, KDIM1, KDIM2, KDIM3, KDIM4)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: HNAME
  INTEGER, INTENT(IN)          :: KDIM1
  INTEGER, OPTIONAL, INTENT(IN):: KDIM2
  INTEGER, OPTIONAL, INTENT(IN):: KDIM3
  INTEGER, OPTIONAL, INTENT(IN):: KDIM4

  SELECT CASE(TRIM(HNAME))
    CASE('XHYPF12')
      ALLOCATE(PARAM_LIMA_WARM%XHYPF12(KDIM1, KDIM2))
      XHYPF12 => PARAM_LIMA_WARM%XHYPF12
    CASE('XHYPF32')
      ALLOCATE(PARAM_LIMA_WARM%XHYPF32(KDIM1, KDIM2))
      XHYPF32 => PARAM_LIMA_WARM%XHYPF32
    CASE('XAHENG')
      ALLOCATE(PARAM_LIMA_WARM%XAHENG(KDIM1))
      XAHENG => PARAM_LIMA_WARM%XAHENG
    CASE('XAHENG2')
      ALLOCATE(PARAM_LIMA_WARM%XAHENG2(KDIM1))
      XAHENG2 => PARAM_LIMA_WARM%XAHENG2
    CASE('XAHENG3')
      ALLOCATE(PARAM_LIMA_WARM%XAHENG3(KDIM1))
      XAHENG3 => PARAM_LIMA_WARM%XAHENG3
    CASE('XPSI1')
      ALLOCATE(PARAM_LIMA_WARM%XPSI1(KDIM1))
      XPSI1 => PARAM_LIMA_WARM%XPSI1
    CASE('XPSI3')
      ALLOCATE(PARAM_LIMA_WARM%XPSI3(KDIM1))
      XPSI3 => PARAM_LIMA_WARM%XPSI3
    CASE('XAHENF')
      ALLOCATE(PARAM_LIMA_WARM%XAHENF(KDIM1))
      XAHENF => PARAM_LIMA_WARM%XAHENF
    CASE('XAHENY')
      ALLOCATE(PARAM_LIMA_WARM%XAHENY(KDIM1))
      XAHENY => PARAM_LIMA_WARM%XAHENY
    CASE('XCONCC_INI')
      ALLOCATE(PARAM_LIMA_WARM%XCONCC_INI(KDIM1,KDIM2,KDIM3,KDIM4))
      XCONCC_INI => PARAM_LIMA_WARM%XCONCC_INI
  END SELECT
END SUBROUTINE PARAM_LIMA_WARM_ALLOCATE
!
!-------------------------------------------------------------------------------
!
END MODULE MODD_PARAM_LIMA_WARM
