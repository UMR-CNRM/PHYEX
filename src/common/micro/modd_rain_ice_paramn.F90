!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!     ######spl
      MODULE MODD_RAIN_ICE_PARAM_n
!     ##########################
!> @file
!!****  *MODD_RAIN_ICE_PARAM_n* - declaration of some microphysical factors
!!                                extensively used in the warm and cold schemes.
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to declare some precomputed
!     microphysical paramters directly used in routine RAIN_ICE.
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (MODD_RAIN_ICE_PARAM)
!!
!!    AUTHOR
!!    ------
!!      J.-P. Pinty  *Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    04/12/95
!!       J.-P. Pinty   29/11/02 add ICE4
!!       S. Riette 11/2016: new ICE3/ICE4 processes
!!
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_PARAMETERS, ONLY: JPMODELMAX
IMPLICIT NONE
!
TYPE RAIN_ICE_PARAM_t
REAL,DIMENSION(2)      :: XFSEDC                 ! Constants for sedimentation fluxes of C
REAL      :: XFSEDR,XEXSEDR,                   & ! Constants for sedimentation
             XFSEDI,XEXCSEDI,XEXRSEDI,         & ! fluxes of R, I, S and G
             XFSEDS,XEXSEDS,                   &
             XFSEDG,XEXSEDG
!
REAL      :: XNU10,XALPHA1,XBETA1,             & ! Constants for heterogeneous
             XNU20,XALPHA2,XBETA2,             & ! ice nucleation : HEN
             XMNU0                               ! mass of nucleated ice crystal
!
REAL      :: XALPHA3,XBETA3,                   & ! Constants for homogeneous
             XHON                                ! ice nucleation : HON
!
REAL      :: XSCFAC,                           & ! Constants for raindrop
             X0EVAR,X1EVAR,XEX0EVAR,XEX1EVAR,  & ! evaporation: EVA and for
             X0DEPI,X2DEPI,                    & ! deposition : DEP on I,
             X0DEPS,X1DEPS,XEX0DEPS,XEX1DEPS,  & !                  on S and
             XRDEPSRED,&
             X0DEPG,X1DEPG,XEX0DEPG,XEX1DEPG,  & !                  on G
             XRDEPGRED
!
REAL      :: XTIMAUTI,XTEXAUTI,XCRIAUTI,       & ! Constants for pristine ice
             XT0CRIAUTI,XACRIAUTI,XBCRIAUTI      ! autoconversion : AUT
!
REAL      :: XCOLIS,XCOLEXIS,                  & ! Constants for snow
             XFIAGGS,                          & ! aggregation : AGG
             XEXIAGGS
!
REAL      :: XTIMAUTC,                         & ! Constants for cloud droplet
             XCRIAUTC                            ! autoconversion : AUT
!
REAL      :: XFCACCR,                          & ! Constants for cloud droplet
             XEXCACCR                            ! accretion on raindrops : ACC
!
REAL      :: XDCSLIM,XCOLCS,                   & ! Constants for the riming of
             XEXCRIMSS,XCRIMSS,                & ! the aggregates : RIM
             XEXCRIMSG,XCRIMSG,                & !
             XEXSRIMCG,XSRIMCG,                & !
             XEXSRIMCG2,XSRIMCG2,              & !
             XSRIMCG3,                         & !
             XGAMINC_BOUND_MIN,                & ! Min val. of Lbda_s for RIM
             XGAMINC_BOUND_MAX,                & ! Max val. of Lbda_s for RIM
             XRIMINTP1,XRIMINTP2                 ! Csts for lin. interpol. of
                                                 ! the tab. incomplete Gamma law
INTEGER      :: NGAMINC                          ! Number of tab. Lbda_s
REAL, DIMENSION(:), ALLOCATABLE          &
                       :: XGAMINC_RIM1,        & ! Tab. incomplete Gamma funct.
                          XGAMINC_RIM2,        & ! for XDS+2 and for XBS
                          XGAMINC_RIM4           ! and for 2+XDS+XBS-XBG
!
REAL      :: XFRACCSS,                         & ! Constants for the accretion
             XLBRACCS1,XLBRACCS2,XLBRACCS3,    & ! raindrops onto the aggregates
             XFSACCRG,                         & ! ACC (processes RACCSS and
             XLBSACCR1,XLBSACCR2,XLBSACCR3,    & !                SACCRG)
             XACCLBDAS_MIN,                    & ! Min val. of Lbda_s for ACC
             XACCLBDAS_MAX,                    & ! Max val. of Lbda_s for ACC
             XACCLBDAR_MIN,                    & ! Min val. of Lbda_r for ACC
             XACCLBDAR_MAX,                    & ! Max val. of Lbda_r for ACC
             XACCINTP1S,XACCINTP2S,            & ! Csts for bilin. interpol. of
             XACCINTP1R,XACCINTP2R               !   Lbda_s and Lbda_r in the
                                                 ! XKER_RACCSS and XKER_SACCRG
                                                 !            tables
INTEGER      :: NACCLBDAS,                     & ! Number of Lbda_s values and
                NACCLBDAR                        !   of Lbda_r values in the
                                                 ! XKER_RACCSS and XKER_SACCRG
                                                 !            tables
REAL,DIMENSION(:,:), ALLOCATABLE         &
                         :: XKER_RACCSS,       & ! Normalized kernel for RACCSS
                            XKER_RACCS,        & ! Normalized kernel for RACCS
                            XKER_SACCRG          ! Normalized kernel for SACCRG
REAL      :: XFSCVMG                             ! Melting-conversion factor of
                                                 ! the aggregates
!
REAL      :: XCOLIR,                           & ! Constants for rain contact
             XEXRCFRI,XRCFRI,                  & ! freezing : CFR
             XEXICFRR,XICFRR                     !
!
REAL      :: XFCDRYG,                          & ! Constants for the dry growth
             XCOLIG,XCOLEXIG,XFIDRYG,          & ! of the graupeln : DRY
             XFIDRYG2, XEXFIDRYG,              &
             XCOLSG,XCOLEXSG,XFSDRYG,          & !   processes RCDRYG
             XLBSDRYG1,XLBSDRYG2,XLBSDRYG3,    & !             RIDRYG
             XFRDRYG,                          & !             RSDRYG
             XLBRDRYG1,XLBRDRYG2,XLBRDRYG3,    & !             RRDRYG
             XDRYLBDAR_MIN,                    & ! Min val. of Lbda_r for DRY
             XDRYLBDAR_MAX,                    & ! Max val. of Lbda_r for DRY
             XDRYLBDAS_MIN,                    & ! Min val. of Lbda_s for DRY
             XDRYLBDAS_MAX,                    & ! Max val. of Lbda_s for DRY
             XDRYLBDAG_MIN,                    & ! Min val. of Lbda_g for DRY
             XDRYLBDAG_MAX,                    & ! Max val. of Lbda_g for DRY
             XDRYINTP1R,XDRYINTP2R,            & ! Csts for bilin. interpol. of
             XDRYINTP1S,XDRYINTP2S,            & ! Lbda_r, Lbda_s and Lbda_g in
             XDRYINTP1G,XDRYINTP2G               ! the XKER_SDRYG and XKER_RDRYG
                                                 !            tables
INTEGER      :: NDRYLBDAR,                     & ! Number of Lbda_r,
                NDRYLBDAS,                     & !        of Lbda_s and
                NDRYLBDAG                        !        of Lbda_g values in
                                                 ! the XKER_SDRYG and XKER_RDRYG
                                                 !            tables
REAL,DIMENSION(:,:), ALLOCATABLE         &
                         :: XKER_SDRYG,        & ! Normalized kernel for SDRYG
                            XKER_RDRYG           ! Normalized kernel for RDRYG
!
! addition of Hail category
!
REAL      :: XFSEDH,XEXSEDH                      ! Constants for sedimentation
!
!
REAL      :: X0DEPH,X1DEPH,XEX0DEPH,XEX1DEPH     ! Constants for deposition
!
REAL      :: XCOLIH, XCOLEXIH,                 & ! Constants for the dry growth
           & XCOLSH, XCOLEXSH,                 & ! of the hail
           & XCOLGH, XCOLEXGH                    !
!
REAL      :: XFWETH,XFSWETH,                   & ! Constants for the wet growth
             XLBSWETH1,XLBSWETH2,XLBSWETH3,    & ! of the hailstones : WET
             XFGWETH,                          & !   processes RSWETH
             XLBGWETH1,XLBGWETH2,XLBGWETH3,    & !             RGWETH
             XFRWETH,                          & !             RRWETH
             XLBRWETH1,XLBRWETH2,XLBRWETH3,    & !
             XWETLBDAS_MIN,                    & ! Min val. of Lbda_s for WET
             XWETLBDAS_MAX,                    & ! Max val. of Lbda_s for WET
             XWETLBDAG_MIN,                    & ! Min val. of Lbda_g for WET
             XWETLBDAG_MAX,                    & ! Max val. of Lbda_g for WET
             XWETLBDAR_MIN,                    & ! Min val. of Lbda_r for WET
             XWETLBDAR_MAX,                    & ! Max val. of Lbda_r for WET
             XWETLBDAH_MIN,                    & ! Min val. of Lbda_h for WET
             XWETLBDAH_MAX,                    & ! Max val. of Lbda_h for WET
             XWETINTP1S,XWETINTP2S,            & ! Csts for bilin. interpol. of
             XWETINTP1G,XWETINTP2G,            & ! Lbda_r, Lbda_s, Lbda_g
             XWETINTP1R,XWETINTP2R,            & ! and Lbda_h in
             XWETINTP1H,XWETINTP2H               ! the XKER_SWETH, XKER_GWETH
                                                 ! and XKER_RWETH tables
INTEGER      :: NWETLBDAS,                     & ! Number of Lbda_s,
                NWETLBDAG,                     & !        of Lbda_g,
                NWETLBDAR,                     & !        of Lbda_r and
                NWETLBDAH                        !        of Lbda_h values in
                                                 ! the XKER_SWETH, XKER_GWETH
                                                 ! and XKER_RWETH tables
REAL,DIMENSION(:,:), ALLOCATABLE         &
                         :: XKER_SWETH,        & ! Normalized kernel for SWETH
                            XKER_GWETH,        & ! Normalized kernel for GWETH
                            XKER_RWETH           ! Normalized kernel for RWETH
REAL, DIMENSION(40) :: XFRMIN                    ! Parmeters to modify melt and growth of graupels etc.
END TYPE RAIN_ICE_PARAM_t
!
TYPE(RAIN_ICE_PARAM_t), DIMENSION(JPMODELMAX), TARGET, SAVE :: RAIN_ICE_PARAM_MODEL
TYPE(RAIN_ICE_PARAM_t), POINTER, SAVE :: RAIN_ICE_PARAMN => NULL()
!
REAL,DIMENSION(:),POINTER :: XFSEDC => NULL()
REAL,DIMENSION(:),POINTER :: XFRMIN => NULL()

REAL,POINTER :: XFSEDR => NULL(), &
                XEXSEDR => NULL(), &
                XFSEDI => NULL(), &
                XEXCSEDI => NULL(), &
                XEXRSEDI => NULL(), &
                XFSEDS => NULL(), &
                XEXSEDS => NULL(), &
                XFSEDG => NULL(), &
                XEXSEDG => NULL(), &
                XNU10 => NULL(), &
                XALPHA1 => NULL(), &
                XBETA1 => NULL(), &
                XNU20 => NULL(), &
                XALPHA2 => NULL(), &
                XBETA2 => NULL(), &
                XMNU0 => NULL(), &
                XALPHA3 => NULL(), &
                XBETA3 => NULL(), &
                XHON => NULL(), &
                XSCFAC => NULL(), &
                X0EVAR => NULL(), &
                X1EVAR => NULL(), &
                XEX0EVAR => NULL(), &
                XEX1EVAR => NULL(), &
                X0DEPI => NULL(), &
                X2DEPI => NULL(), &
                X0DEPS => NULL(), &
                X1DEPS => NULL(), &
                XEX0DEPS => NULL(), &
                XEX1DEPS => NULL(), &
                XRDEPSRED => NULL(), &
                X0DEPG => NULL(), &
                X1DEPG => NULL(), &
                XEX0DEPG => NULL(), &
                XEX1DEPG => NULL(), &
                XRDEPGRED => NULL(), &
                XTIMAUTI => NULL(), &
                XTEXAUTI => NULL(), &
                XCRIAUTI => NULL(), &
                XT0CRIAUTI => NULL(), &
                XACRIAUTI => NULL(), &
                XBCRIAUTI => NULL(), &
                XCOLIS => NULL(), &
                XCOLEXIS => NULL(), &
                XFIAGGS => NULL(), &
                XEXIAGGS => NULL(), &
                XTIMAUTC => NULL(), &
                XCRIAUTC => NULL(), &
                XFCACCR => NULL(), &
                XEXCACCR => NULL(), &
                XDCSLIM => NULL(), &
                XCOLCS => NULL(), &
                XEXCRIMSS => NULL(), &
                XCRIMSS => NULL(), &
                XEXCRIMSG => NULL(), &
                XCRIMSG => NULL(), &
                XEXSRIMCG => NULL(), &
                XSRIMCG => NULL(), &
                XEXSRIMCG2 => NULL(), &
                XSRIMCG2 => NULL(), &
                XSRIMCG3 => NULL(), &
                XGAMINC_BOUND_MIN => NULL(), &
                XGAMINC_BOUND_MAX => NULL(), &
                XRIMINTP1 => NULL(), &
                XRIMINTP2 => NULL(), &
                XFRACCSS => NULL(), &
                XLBRACCS1 => NULL(), &
                XLBRACCS2 => NULL(), &
                XLBRACCS3 => NULL(), &
                XFSACCRG => NULL(), &
                XLBSACCR1 => NULL(), &
                XLBSACCR2 => NULL(), &
                XLBSACCR3 => NULL(), &
                XACCLBDAS_MIN => NULL(), &
                XACCLBDAS_MAX => NULL(), &
                XACCLBDAR_MIN => NULL(), &
                XACCLBDAR_MAX => NULL(), &
                XACCINTP1S => NULL(), &
                XACCINTP2S => NULL(), &
                XACCINTP1R => NULL(), &
                XACCINTP2R => NULL(), &
                XFSCVMG => NULL(), &
                XCOLIR => NULL(), &
                XEXRCFRI => NULL(), &
                XRCFRI => NULL(), &
                XEXICFRR => NULL(), &
                XICFRR => NULL(), &
                XFCDRYG => NULL(), &
                XCOLIG => NULL(), &
                XCOLEXIG => NULL(), &
                XFIDRYG => NULL(), &
                XFIDRYG2 => NULL(), &
                XEXFIDRYG => NULL(), &
                XCOLSG => NULL(), &
                XCOLEXSG => NULL(), &
                XFSDRYG => NULL(), &
                XLBSDRYG1 => NULL(), &
                XLBSDRYG2 => NULL(), &
                XLBSDRYG3 => NULL(), &
                XFRDRYG => NULL(), &
                XLBRDRYG1 => NULL(), &
                XLBRDRYG2 => NULL(), &
                XLBRDRYG3 => NULL(), &
                XDRYLBDAR_MIN => NULL(), &
                XDRYLBDAR_MAX => NULL(), &
                XDRYLBDAS_MIN => NULL(), &
                XDRYLBDAS_MAX => NULL(), &
                XDRYLBDAG_MIN => NULL(), &
                XDRYLBDAG_MAX => NULL(), &
                XDRYINTP1R => NULL(), &
                XDRYINTP2R => NULL(), &
                XDRYINTP1S => NULL(), &
                XDRYINTP2S => NULL(), &
                XDRYINTP1G => NULL(), &
                XDRYINTP2G => NULL(), &
                XFSEDH => NULL(), &
                XEXSEDH => NULL(), &
                X0DEPH => NULL(), &
                X1DEPH => NULL(), &
                XEX0DEPH => NULL(), &
                XEX1DEPH => NULL(), &
                XCOLIH => NULL(), &
                XCOLEXIH => NULL(), &
                XCOLSH => NULL(), &
                XCOLEXSH => NULL(), &
                XCOLGH => NULL(), &
                XCOLEXGH => NULL(), &
                XFWETH => NULL(), &
                XFSWETH => NULL(), &
                XLBSWETH1 => NULL(), &
                XLBSWETH2 => NULL(), &
                XLBSWETH3 => NULL(), &
                XFGWETH => NULL(), &
                XLBGWETH1 => NULL(), &
                XLBGWETH2 => NULL(), &
                XLBGWETH3 => NULL(), &
                XFRWETH => NULL(), &
                XLBRWETH1 => NULL(), &
                XLBRWETH2 => NULL(), &
                XLBRWETH3 => NULL(), &
                XWETLBDAS_MIN => NULL(), &
                XWETLBDAS_MAX => NULL(), &
                XWETLBDAG_MIN => NULL(), &
                XWETLBDAG_MAX => NULL(), &
                XWETLBDAR_MIN => NULL(), &
                XWETLBDAR_MAX => NULL(), &
                XWETLBDAH_MIN => NULL(), &
                XWETLBDAH_MAX => NULL(), &
                XWETINTP1S => NULL(), &
                XWETINTP2S => NULL(), &
                XWETINTP1G => NULL(), &
                XWETINTP2G => NULL(), &
                XWETINTP1R => NULL(), &
                XWETINTP2R => NULL(), &
                XWETINTP1H => NULL(), &
                XWETINTP2H => NULL()

INTEGER, POINTER :: NGAMINC => NULL(), &
                    NACCLBDAS => NULL(), &
                    NACCLBDAR => NULL(), &
                    NDRYLBDAR => NULL(), &
                    NDRYLBDAS => NULL(), &
                    NDRYLBDAG => NULL(), &
                    NWETLBDAS => NULL(), &
                    NWETLBDAG => NULL(), &
                    NWETLBDAR => NULL(), &
                    NWETLBDAH => NULL()

REAL, DIMENSION(:), POINTER :: XGAMINC_RIM1 => NULL(), &
                               XGAMINC_RIM2 => NULL(), &
                               XGAMINC_RIM4 => NULL()

REAL,DIMENSION(:,:), POINTER :: XKER_RACCSS => NULL(), &
                                XKER_RACCS => NULL(), &
                                XKER_SACCRG => NULL(), &
                                XKER_SDRYG => NULL(), &
                                XKER_RDRYG => NULL(), &
                                XKER_SWETH => NULL(), &
                                XKER_GWETH => NULL(), &
                                XKER_RWETH => NULL()
CONTAINS
SUBROUTINE RAIN_ICE_PARAM_GOTO_MODEL(KFROM, KTO)
!! This subroutine associate all the pointers to the right component of
!! the right strucuture. A value can be accessed through the structure RAIN_ICE_PARAMN
!! or through the strucuture RAIN_ICE_PARAM_MODEL(KTO) or directly through these pointers.
IMPLICIT NONE
INTEGER, INTENT(IN) :: KFROM, KTO
!
IF(.NOT. ASSOCIATED(RAIN_ICE_PARAMN, RAIN_ICE_PARAM_MODEL(KTO))) THEN
  !
  RAIN_ICE_PARAMN => RAIN_ICE_PARAM_MODEL(KTO)
  !
  XFSEDC => RAIN_ICE_PARAMN%XFSEDC
  XFRMIN => RAIN_ICE_PARAMN%XFRMIN
  !
  XFSEDR => RAIN_ICE_PARAMN%XFSEDR
  XEXSEDR => RAIN_ICE_PARAMN%XEXSEDR
  XFSEDI => RAIN_ICE_PARAMN%XFSEDI
  XEXCSEDI => RAIN_ICE_PARAMN%XEXCSEDI
  XEXRSEDI => RAIN_ICE_PARAMN%XEXRSEDI
  XFSEDS => RAIN_ICE_PARAMN%XFSEDS
  XEXSEDS => RAIN_ICE_PARAMN%XEXSEDS
  XFSEDG => RAIN_ICE_PARAMN%XFSEDG
  XEXSEDG => RAIN_ICE_PARAMN%XEXSEDG
  XNU10 => RAIN_ICE_PARAMN%XNU10
  XALPHA1 => RAIN_ICE_PARAMN%XALPHA1
  XBETA1 => RAIN_ICE_PARAMN%XBETA1
  XNU20 => RAIN_ICE_PARAMN%XNU20
  XALPHA2 => RAIN_ICE_PARAMN%XALPHA2
  XBETA2 => RAIN_ICE_PARAMN%XBETA2
  XMNU0 => RAIN_ICE_PARAMN%XMNU0
  XALPHA3 => RAIN_ICE_PARAMN%XALPHA3
  XBETA3 => RAIN_ICE_PARAMN%XBETA3
  XHON => RAIN_ICE_PARAMN%XHON
  XSCFAC => RAIN_ICE_PARAMN%XSCFAC
  X0EVAR => RAIN_ICE_PARAMN%X0EVAR
  X1EVAR => RAIN_ICE_PARAMN%X1EVAR
  XEX0EVAR => RAIN_ICE_PARAMN%XEX0EVAR
  XEX1EVAR => RAIN_ICE_PARAMN%XEX1EVAR
  X0DEPI => RAIN_ICE_PARAMN%X0DEPI
  X2DEPI => RAIN_ICE_PARAMN%X2DEPI
  X0DEPS => RAIN_ICE_PARAMN%X0DEPS
  X1DEPS => RAIN_ICE_PARAMN%X1DEPS
  XEX0DEPS => RAIN_ICE_PARAMN%XEX0DEPS
  XEX1DEPS => RAIN_ICE_PARAMN%XEX1DEPS
  XRDEPSRED => RAIN_ICE_PARAMN%XRDEPSRED
  X0DEPG => RAIN_ICE_PARAMN%X0DEPG
  X1DEPG => RAIN_ICE_PARAMN%X1DEPG
  XEX0DEPG => RAIN_ICE_PARAMN%XEX0DEPG
  XEX1DEPG => RAIN_ICE_PARAMN%XEX1DEPG
  XRDEPGRED => RAIN_ICE_PARAMN%XRDEPGRED
  XTIMAUTI => RAIN_ICE_PARAMN%XTIMAUTI
  XTEXAUTI => RAIN_ICE_PARAMN%XTEXAUTI
  XCRIAUTI => RAIN_ICE_PARAMN%XCRIAUTI
  XT0CRIAUTI => RAIN_ICE_PARAMN%XT0CRIAUTI
  XACRIAUTI => RAIN_ICE_PARAMN%XACRIAUTI
  XBCRIAUTI => RAIN_ICE_PARAMN%XBCRIAUTI
  XCOLIS => RAIN_ICE_PARAMN%XCOLIS
  XCOLEXIS => RAIN_ICE_PARAMN%XCOLEXIS
  XFIAGGS => RAIN_ICE_PARAMN%XFIAGGS
  XEXIAGGS => RAIN_ICE_PARAMN%XEXIAGGS
  XTIMAUTC => RAIN_ICE_PARAMN%XTIMAUTC
  XCRIAUTC => RAIN_ICE_PARAMN%XCRIAUTC
  XFCACCR => RAIN_ICE_PARAMN%XFCACCR
  XEXCACCR => RAIN_ICE_PARAMN%XEXCACCR
  XDCSLIM => RAIN_ICE_PARAMN%XDCSLIM
  XCOLCS => RAIN_ICE_PARAMN%XCOLCS
  XEXCRIMSS => RAIN_ICE_PARAMN%XEXCRIMSS
  XCRIMSS => RAIN_ICE_PARAMN%XCRIMSS
  XEXCRIMSG => RAIN_ICE_PARAMN%XEXCRIMSG
  XCRIMSG => RAIN_ICE_PARAMN%XCRIMSG
  XEXSRIMCG => RAIN_ICE_PARAMN%XEXSRIMCG
  XSRIMCG => RAIN_ICE_PARAMN%XSRIMCG
  XEXSRIMCG2 => RAIN_ICE_PARAMN%XEXSRIMCG2
  XSRIMCG2 => RAIN_ICE_PARAMN%XSRIMCG2
  XSRIMCG3 => RAIN_ICE_PARAMN%XSRIMCG3
  XGAMINC_BOUND_MIN => RAIN_ICE_PARAMN%XGAMINC_BOUND_MIN
  XGAMINC_BOUND_MAX => RAIN_ICE_PARAMN%XGAMINC_BOUND_MAX
  XRIMINTP1 => RAIN_ICE_PARAMN%XRIMINTP1
  XRIMINTP2 => RAIN_ICE_PARAMN%XRIMINTP2
  XFRACCSS => RAIN_ICE_PARAMN%XFRACCSS
  XLBRACCS1 => RAIN_ICE_PARAMN%XLBRACCS1
  XLBRACCS2 => RAIN_ICE_PARAMN%XLBRACCS2
  XLBRACCS3 => RAIN_ICE_PARAMN%XLBRACCS3
  XFSACCRG => RAIN_ICE_PARAMN%XFSACCRG
  XLBSACCR1 => RAIN_ICE_PARAMN%XLBSACCR1
  XLBSACCR2 => RAIN_ICE_PARAMN%XLBSACCR2
  XLBSACCR3 => RAIN_ICE_PARAMN%XLBSACCR3
  XACCLBDAS_MIN => RAIN_ICE_PARAMN%XACCLBDAS_MIN
  XACCLBDAS_MAX => RAIN_ICE_PARAMN%XACCLBDAS_MAX
  XACCLBDAR_MIN => RAIN_ICE_PARAMN%XACCLBDAR_MIN
  XACCLBDAR_MAX => RAIN_ICE_PARAMN%XACCLBDAR_MAX
  XACCINTP1S => RAIN_ICE_PARAMN%XACCINTP1S
  XACCINTP2S => RAIN_ICE_PARAMN%XACCINTP2S
  XACCINTP1R => RAIN_ICE_PARAMN%XACCINTP1R
  XACCINTP2R => RAIN_ICE_PARAMN%XACCINTP2R
  XFSCVMG => RAIN_ICE_PARAMN%XFSCVMG
  XCOLIR => RAIN_ICE_PARAMN%XCOLIR
  XEXRCFRI => RAIN_ICE_PARAMN%XEXRCFRI
  XRCFRI => RAIN_ICE_PARAMN%XRCFRI
  XEXICFRR => RAIN_ICE_PARAMN%XEXICFRR
  XICFRR => RAIN_ICE_PARAMN%XICFRR
  XFCDRYG => RAIN_ICE_PARAMN%XFCDRYG
  XCOLIG => RAIN_ICE_PARAMN%XCOLIG
  XCOLEXIG => RAIN_ICE_PARAMN%XCOLEXIG
  XFIDRYG => RAIN_ICE_PARAMN%XFIDRYG
  XFIDRYG2 => RAIN_ICE_PARAMN%XFIDRYG2
  XEXFIDRYG => RAIN_ICE_PARAMN%XEXFIDRYG
  XCOLSG => RAIN_ICE_PARAMN%XCOLSG
  XCOLEXSG => RAIN_ICE_PARAMN%XCOLEXSG
  XFSDRYG => RAIN_ICE_PARAMN%XFSDRYG
  XLBSDRYG1 => RAIN_ICE_PARAMN%XLBSDRYG1
  XLBSDRYG2 => RAIN_ICE_PARAMN%XLBSDRYG2
  XLBSDRYG3 => RAIN_ICE_PARAMN%XLBSDRYG3
  XFRDRYG => RAIN_ICE_PARAMN%XFRDRYG
  XLBRDRYG1 => RAIN_ICE_PARAMN%XLBRDRYG1
  XLBRDRYG2 => RAIN_ICE_PARAMN%XLBRDRYG2
  XLBRDRYG3 => RAIN_ICE_PARAMN%XLBRDRYG3
  XDRYLBDAR_MIN => RAIN_ICE_PARAMN%XDRYLBDAR_MIN
  XDRYLBDAR_MAX => RAIN_ICE_PARAMN%XDRYLBDAR_MAX
  XDRYLBDAS_MIN => RAIN_ICE_PARAMN%XDRYLBDAS_MIN
  XDRYLBDAS_MAX => RAIN_ICE_PARAMN%XDRYLBDAS_MAX
  XDRYLBDAG_MIN => RAIN_ICE_PARAMN%XDRYLBDAG_MIN
  XDRYLBDAG_MAX => RAIN_ICE_PARAMN%XDRYLBDAG_MAX
  XDRYINTP1R => RAIN_ICE_PARAMN%XDRYINTP1R
  XDRYINTP2R => RAIN_ICE_PARAMN%XDRYINTP2R
  XDRYINTP1S => RAIN_ICE_PARAMN%XDRYINTP1S
  XDRYINTP2S => RAIN_ICE_PARAMN%XDRYINTP2S
  XDRYINTP1G => RAIN_ICE_PARAMN%XDRYINTP1G
  XDRYINTP2G => RAIN_ICE_PARAMN%XDRYINTP2G
  XFSEDH => RAIN_ICE_PARAMN%XFSEDH
  XEXSEDH => RAIN_ICE_PARAMN%XEXSEDH
  X0DEPH => RAIN_ICE_PARAMN%X0DEPH
  X1DEPH => RAIN_ICE_PARAMN%X1DEPH
  XEX0DEPH => RAIN_ICE_PARAMN%XEX0DEPH
  XEX1DEPH => RAIN_ICE_PARAMN%XEX1DEPH
  XCOLIH => RAIN_ICE_PARAMN%XCOLIH
  XCOLEXIH => RAIN_ICE_PARAMN%XCOLEXIH
  XCOLSH => RAIN_ICE_PARAMN%XCOLSH
  XCOLEXSH => RAIN_ICE_PARAMN%XCOLEXSH
  XCOLGH => RAIN_ICE_PARAMN%XCOLGH
  XCOLEXGH => RAIN_ICE_PARAMN%XCOLEXGH
  XFWETH => RAIN_ICE_PARAMN%XFWETH
  XFSWETH => RAIN_ICE_PARAMN%XFSWETH
  XLBSWETH1 => RAIN_ICE_PARAMN%XLBSWETH1
  XLBSWETH2 => RAIN_ICE_PARAMN%XLBSWETH2
  XLBSWETH3 => RAIN_ICE_PARAMN%XLBSWETH3
  XFGWETH => RAIN_ICE_PARAMN%XFGWETH
  XLBGWETH1 => RAIN_ICE_PARAMN%XLBGWETH1
  XLBGWETH2 => RAIN_ICE_PARAMN%XLBGWETH2
  XLBGWETH3 => RAIN_ICE_PARAMN%XLBGWETH3
  XFRWETH => RAIN_ICE_PARAMN%XFRWETH
  XLBRWETH1 => RAIN_ICE_PARAMN%XLBRWETH1
  XLBRWETH2 => RAIN_ICE_PARAMN%XLBRWETH2
  XLBRWETH3 => RAIN_ICE_PARAMN%XLBRWETH3
  XWETLBDAS_MIN => RAIN_ICE_PARAMN%XWETLBDAS_MIN
  XWETLBDAS_MAX => RAIN_ICE_PARAMN%XWETLBDAS_MAX
  XWETLBDAG_MIN => RAIN_ICE_PARAMN%XWETLBDAG_MIN
  XWETLBDAG_MAX => RAIN_ICE_PARAMN%XWETLBDAG_MAX
  XWETLBDAR_MIN => RAIN_ICE_PARAMN%XWETLBDAR_MIN
  XWETLBDAR_MAX => RAIN_ICE_PARAMN%XWETLBDAR_MAX
  XWETLBDAH_MIN => RAIN_ICE_PARAMN%XWETLBDAH_MIN
  XWETLBDAH_MAX => RAIN_ICE_PARAMN%XWETLBDAH_MAX
  XWETINTP1S => RAIN_ICE_PARAMN%XWETINTP1S
  XWETINTP2S => RAIN_ICE_PARAMN%XWETINTP2S
  XWETINTP1G => RAIN_ICE_PARAMN%XWETINTP1G
  XWETINTP2G => RAIN_ICE_PARAMN%XWETINTP2G
  XWETINTP1R => RAIN_ICE_PARAMN%XWETINTP1R
  XWETINTP2R => RAIN_ICE_PARAMN%XWETINTP2R
  XWETINTP1H => RAIN_ICE_PARAMN%XWETINTP1H
  XWETINTP2H => RAIN_ICE_PARAMN%XWETINTP2H
  !
  NGAMINC => RAIN_ICE_PARAMN%NGAMINC
  NACCLBDAS => RAIN_ICE_PARAMN%NACCLBDAS
  NACCLBDAR => RAIN_ICE_PARAMN%NACCLBDAR
  NDRYLBDAR => RAIN_ICE_PARAMN%NDRYLBDAR
  NDRYLBDAS => RAIN_ICE_PARAMN%NDRYLBDAS
  NDRYLBDAG => RAIN_ICE_PARAMN%NDRYLBDAG
  NWETLBDAS => RAIN_ICE_PARAMN%NWETLBDAS
  NWETLBDAG => RAIN_ICE_PARAMN%NWETLBDAG
  NWETLBDAR => RAIN_ICE_PARAMN%NWETLBDAR
  NWETLBDAH => RAIN_ICE_PARAMN%NWETLBDAH
  !
  CALL HELPER1(XGAMINC_RIM1, RAIN_ICE_PARAMN%XGAMINC_RIM1)
  CALL HELPER1(XGAMINC_RIM2, RAIN_ICE_PARAMN%XGAMINC_RIM2)
  CALL HELPER1(XGAMINC_RIM4, RAIN_ICE_PARAMN%XGAMINC_RIM4)
  CALL HELPER2(XKER_RACCSS, RAIN_ICE_PARAMN%XKER_RACCSS)
  CALL HELPER2(XKER_RACCS, RAIN_ICE_PARAMN%XKER_RACCS)
  CALL HELPER2(XKER_SACCRG, RAIN_ICE_PARAMN%XKER_SACCRG)
  CALL HELPER2(XKER_SDRYG, RAIN_ICE_PARAMN%XKER_SDRYG)
  CALL HELPER2(XKER_RDRYG, RAIN_ICE_PARAMN%XKER_RDRYG)
  CALL HELPER2(XKER_SWETH, RAIN_ICE_PARAMN%XKER_SWETH)
  CALL HELPER2(XKER_GWETH, RAIN_ICE_PARAMN%XKER_GWETH)
  CALL HELPER2(XKER_RWETH, RAIN_ICE_PARAMN%XKER_RWETH)

ENDIF
END SUBROUTINE RAIN_ICE_PARAM_GOTO_MODEL
!
SUBROUTINE HELPER1(XPT, XARR)
  IMPLICIT NONE
  REAL, POINTER, DIMENSION(:), INTENT(INOUT) :: XPT
  REAL, DIMENSION(:), ALLOCATABLE, TARGET, INTENT(INOUT) :: XARR
  IF(ALLOCATED(XARR)) THEN
    XPT => XARR
  ELSE
    XPT => NULL()
  ENDIF
END SUBROUTINE HELPER1
!
SUBROUTINE HELPER2(XPT, XARR)
  IMPLICIT NONE
  REAL, POINTER, DIMENSION(:,:), INTENT(INOUT) :: XPT
  REAL, DIMENSION(:,:), ALLOCATABLE, TARGET, INTENT(INOUT) :: XARR
  IF(ALLOCATED(XARR)) THEN
    XPT => XARR
  ELSE
    XPT => NULL()
  ENDIF
END SUBROUTINE HELPER2
!
SUBROUTINE RAIN_ICE_PARAM_ALLOCATE(HNAME, KDIM1, KDIM2)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: HNAME
  INTEGER, INTENT(IN)           :: KDIM1
  INTEGER, OPTIONAL, INTENT(IN) :: KDIM2

  SELECT CASE(TRIM(HNAME))
    !1D arrays
    CASE('XGAMINC_RIM1')
      ALLOCATE(RAIN_ICE_PARAMN%XGAMINC_RIM1(KDIM1))
      XGAMINC_RIM1 => RAIN_ICE_PARAMN%XGAMINC_RIM1
    CASE('XGAMINC_RIM2')
      ALLOCATE(RAIN_ICE_PARAMN%XGAMINC_RIM2(KDIM1))
      XGAMINC_RIM2 => RAIN_ICE_PARAMN%XGAMINC_RIM2
    CASE('XGAMINC_RIM4')
      ALLOCATE(RAIN_ICE_PARAMN%XGAMINC_RIM4(KDIM1))
      XGAMINC_RIM4 => RAIN_ICE_PARAMN%XGAMINC_RIM4
    !
    !2D arrays
    CASE('XKER_RACCSS')
      ALLOCATE(RAIN_ICE_PARAMN%XKER_RACCSS(KDIM1, KDIM2))
      XKER_RACCSS=> RAIN_ICE_PARAMN%XKER_RACCSS
    CASE('XKER_RACCS')
      ALLOCATE(RAIN_ICE_PARAMN%XKER_RACCS(KDIM1, KDIM2))
      XKER_RACCS=> RAIN_ICE_PARAMN%XKER_RACCS
    CASE('XKER_SACCRG')
      ALLOCATE(RAIN_ICE_PARAMN%XKER_SACCRG(KDIM1, KDIM2))
      XKER_SACCRG=> RAIN_ICE_PARAMN%XKER_SACCRG
    CASE('XKER_SDRYG')
      ALLOCATE(RAIN_ICE_PARAMN%XKER_SDRYG(KDIM1, KDIM2))
      XKER_SDRYG=> RAIN_ICE_PARAMN%XKER_SDRYG
    CASE('XKER_RDRYG')
      ALLOCATE(RAIN_ICE_PARAMN%XKER_RDRYG(KDIM1, KDIM2))
      XKER_RDRYG=> RAIN_ICE_PARAMN%XKER_RDRYG
    CASE('XKER_SWETH')
      ALLOCATE(RAIN_ICE_PARAMN%XKER_SWETH(KDIM1, KDIM2))
      XKER_SWETH=> RAIN_ICE_PARAMN%XKER_SWETH
    CASE('XKER_GWETH')
      ALLOCATE(RAIN_ICE_PARAMN%XKER_GWETH(KDIM1, KDIM2))
      XKER_GWETH=> RAIN_ICE_PARAMN%XKER_GWETH
    CASE('XKER_RWETH')
      ALLOCATE(RAIN_ICE_PARAMN%XKER_RWETH(KDIM1, KDIM2))
      XKER_RWETH=> RAIN_ICE_PARAMN%XKER_RWETH
  END SELECT
END SUBROUTINE RAIN_ICE_PARAM_ALLOCATE
!
SUBROUTINE RAIN_ICE_PARAM_DEALLOCATE()
  IMPLICIT NONE
  XGAMINC_RIM1=>NULL()
  DEALLOCATE(RAIN_ICE_PARAMN%XGAMINC_RIM1)
  XGAMINC_RIM2=>NULL()
  DEALLOCATE(RAIN_ICE_PARAMN%XGAMINC_RIM2)
  XKER_RACCSS=>NULL()
  DEALLOCATE(RAIN_ICE_PARAMN%XKER_RACCSS)
  XKER_RACCS=>NULL()
  DEALLOCATE(RAIN_ICE_PARAMN%XKER_RACCS)
  XKER_SACCRG=>NULL()
  DEALLOCATE(RAIN_ICE_PARAMN%XKER_SACCRG)
  XKER_SDRYG=>NULL()
  DEALLOCATE(RAIN_ICE_PARAMN%XKER_SDRYG)
  XKER_RDRYG=>NULL()
  DEALLOCATE(RAIN_ICE_PARAMN%XKER_RDRYG)
END SUBROUTINE RAIN_ICE_PARAM_DEALLOCATE
!
END MODULE MODD_RAIN_ICE_PARAM_n
