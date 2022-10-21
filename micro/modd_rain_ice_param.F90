!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!     ######spl
      MODULE MODD_RAIN_ICE_PARAM
!     ##########################
!
!!****  *MODD_RAIN_ICE_PARAM* - declaration of some microphysical factors
!!                              extensively used in the warm and cold schemes.
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
TYPE(RAIN_ICE_PARAM_t), SAVE, TARGET :: RAIN_ICE_PARAM
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
SUBROUTINE RAIN_ICE_PARAM_ASSOCIATE()
  IMPLICIT NONE
  XFSEDC => RAIN_ICE_PARAM%XFSEDC
  XFRMIN => RAIN_ICE_PARAM%XFRMIN
  !
  XFSEDR => RAIN_ICE_PARAM%XFSEDR
  XEXSEDR => RAIN_ICE_PARAM%XEXSEDR
  XFSEDI => RAIN_ICE_PARAM%XFSEDI
  XEXCSEDI => RAIN_ICE_PARAM%XEXCSEDI
  XEXRSEDI => RAIN_ICE_PARAM%XEXRSEDI
  XFSEDS => RAIN_ICE_PARAM%XFSEDS
  XEXSEDS => RAIN_ICE_PARAM%XEXSEDS
  XFSEDG => RAIN_ICE_PARAM%XFSEDG
  XEXSEDG => RAIN_ICE_PARAM%XEXSEDG
  XNU10 => RAIN_ICE_PARAM%XNU10
  XALPHA1 => RAIN_ICE_PARAM%XALPHA1
  XBETA1 => RAIN_ICE_PARAM%XBETA1
  XNU20 => RAIN_ICE_PARAM%XNU20
  XALPHA2 => RAIN_ICE_PARAM%XALPHA2
  XBETA2 => RAIN_ICE_PARAM%XBETA2
  XMNU0 => RAIN_ICE_PARAM%XMNU0
  XALPHA3 => RAIN_ICE_PARAM%XALPHA3
  XBETA3 => RAIN_ICE_PARAM%XBETA3
  XHON => RAIN_ICE_PARAM%XHON
  XSCFAC => RAIN_ICE_PARAM%XSCFAC
  X0EVAR => RAIN_ICE_PARAM%X0EVAR
  X1EVAR => RAIN_ICE_PARAM%X1EVAR
  XEX0EVAR => RAIN_ICE_PARAM%XEX0EVAR
  XEX1EVAR => RAIN_ICE_PARAM%XEX1EVAR
  X0DEPI => RAIN_ICE_PARAM%X0DEPI
  X2DEPI => RAIN_ICE_PARAM%X2DEPI
  X0DEPS => RAIN_ICE_PARAM%X0DEPS
  X1DEPS => RAIN_ICE_PARAM%X1DEPS
  XEX0DEPS => RAIN_ICE_PARAM%XEX0DEPS
  XEX1DEPS => RAIN_ICE_PARAM%XEX1DEPS
  XRDEPSRED => RAIN_ICE_PARAM%XRDEPSRED
  X0DEPG => RAIN_ICE_PARAM%X0DEPG
  X1DEPG => RAIN_ICE_PARAM%X1DEPG
  XEX0DEPG => RAIN_ICE_PARAM%XEX0DEPG
  XEX1DEPG => RAIN_ICE_PARAM%XEX1DEPG
  XRDEPGRED => RAIN_ICE_PARAM%XRDEPGRED
  XTIMAUTI => RAIN_ICE_PARAM%XTIMAUTI
  XTEXAUTI => RAIN_ICE_PARAM%XTEXAUTI
  XCRIAUTI => RAIN_ICE_PARAM%XCRIAUTI
  XT0CRIAUTI => RAIN_ICE_PARAM%XT0CRIAUTI
  XACRIAUTI => RAIN_ICE_PARAM%XACRIAUTI
  XBCRIAUTI => RAIN_ICE_PARAM%XBCRIAUTI
  XCOLIS => RAIN_ICE_PARAM%XCOLIS
  XCOLEXIS => RAIN_ICE_PARAM%XCOLEXIS
  XFIAGGS => RAIN_ICE_PARAM%XFIAGGS
  XEXIAGGS => RAIN_ICE_PARAM%XEXIAGGS
  XTIMAUTC => RAIN_ICE_PARAM%XTIMAUTC
  XCRIAUTC => RAIN_ICE_PARAM%XCRIAUTC
  XFCACCR => RAIN_ICE_PARAM%XFCACCR
  XEXCACCR => RAIN_ICE_PARAM%XEXCACCR
  XDCSLIM => RAIN_ICE_PARAM%XDCSLIM
  XCOLCS => RAIN_ICE_PARAM%XCOLCS
  XEXCRIMSS => RAIN_ICE_PARAM%XEXCRIMSS
  XCRIMSS => RAIN_ICE_PARAM%XCRIMSS
  XEXCRIMSG => RAIN_ICE_PARAM%XEXCRIMSG
  XCRIMSG => RAIN_ICE_PARAM%XCRIMSG
  XEXSRIMCG => RAIN_ICE_PARAM%XEXSRIMCG
  XSRIMCG => RAIN_ICE_PARAM%XSRIMCG
  XEXSRIMCG2 => RAIN_ICE_PARAM%XEXSRIMCG2
  XSRIMCG2 => RAIN_ICE_PARAM%XSRIMCG2
  XSRIMCG3 => RAIN_ICE_PARAM%XSRIMCG3
  XGAMINC_BOUND_MIN => RAIN_ICE_PARAM%XGAMINC_BOUND_MIN
  XGAMINC_BOUND_MAX => RAIN_ICE_PARAM%XGAMINC_BOUND_MAX
  XRIMINTP1 => RAIN_ICE_PARAM%XRIMINTP1
  XRIMINTP2 => RAIN_ICE_PARAM%XRIMINTP2
  XFRACCSS => RAIN_ICE_PARAM%XFRACCSS
  XLBRACCS1 => RAIN_ICE_PARAM%XLBRACCS1
  XLBRACCS2 => RAIN_ICE_PARAM%XLBRACCS2
  XLBRACCS3 => RAIN_ICE_PARAM%XLBRACCS3
  XFSACCRG => RAIN_ICE_PARAM%XFSACCRG
  XLBSACCR1 => RAIN_ICE_PARAM%XLBSACCR1
  XLBSACCR2 => RAIN_ICE_PARAM%XLBSACCR2
  XLBSACCR3 => RAIN_ICE_PARAM%XLBSACCR3
  XACCLBDAS_MIN => RAIN_ICE_PARAM%XACCLBDAS_MIN
  XACCLBDAS_MAX => RAIN_ICE_PARAM%XACCLBDAS_MAX
  XACCLBDAR_MIN => RAIN_ICE_PARAM%XACCLBDAR_MIN
  XACCLBDAR_MAX => RAIN_ICE_PARAM%XACCLBDAR_MAX
  XACCINTP1S => RAIN_ICE_PARAM%XACCINTP1S
  XACCINTP2S => RAIN_ICE_PARAM%XACCINTP2S
  XACCINTP1R => RAIN_ICE_PARAM%XACCINTP1R
  XACCINTP2R => RAIN_ICE_PARAM%XACCINTP2R
  XFSCVMG => RAIN_ICE_PARAM%XFSCVMG
  XCOLIR => RAIN_ICE_PARAM%XCOLIR
  XEXRCFRI => RAIN_ICE_PARAM%XEXRCFRI
  XRCFRI => RAIN_ICE_PARAM%XRCFRI
  XEXICFRR => RAIN_ICE_PARAM%XEXICFRR
  XICFRR => RAIN_ICE_PARAM%XICFRR
  XFCDRYG => RAIN_ICE_PARAM%XFCDRYG
  XCOLIG => RAIN_ICE_PARAM%XCOLIG
  XCOLEXIG => RAIN_ICE_PARAM%XCOLEXIG
  XFIDRYG => RAIN_ICE_PARAM%XFIDRYG
  XFIDRYG2 => RAIN_ICE_PARAM%XFIDRYG2
  XEXFIDRYG => RAIN_ICE_PARAM%XEXFIDRYG
  XCOLSG => RAIN_ICE_PARAM%XCOLSG
  XCOLEXSG => RAIN_ICE_PARAM%XCOLEXSG
  XFSDRYG => RAIN_ICE_PARAM%XFSDRYG
  XLBSDRYG1 => RAIN_ICE_PARAM%XLBSDRYG1
  XLBSDRYG2 => RAIN_ICE_PARAM%XLBSDRYG2
  XLBSDRYG3 => RAIN_ICE_PARAM%XLBSDRYG3
  XFRDRYG => RAIN_ICE_PARAM%XFRDRYG
  XLBRDRYG1 => RAIN_ICE_PARAM%XLBRDRYG1
  XLBRDRYG2 => RAIN_ICE_PARAM%XLBRDRYG2
  XLBRDRYG3 => RAIN_ICE_PARAM%XLBRDRYG3
  XDRYLBDAR_MIN => RAIN_ICE_PARAM%XDRYLBDAR_MIN
  XDRYLBDAR_MAX => RAIN_ICE_PARAM%XDRYLBDAR_MAX
  XDRYLBDAS_MIN => RAIN_ICE_PARAM%XDRYLBDAS_MIN
  XDRYLBDAS_MAX => RAIN_ICE_PARAM%XDRYLBDAS_MAX
  XDRYLBDAG_MIN => RAIN_ICE_PARAM%XDRYLBDAG_MIN
  XDRYLBDAG_MAX => RAIN_ICE_PARAM%XDRYLBDAG_MAX
  XDRYINTP1R => RAIN_ICE_PARAM%XDRYINTP1R
  XDRYINTP2R => RAIN_ICE_PARAM%XDRYINTP2R
  XDRYINTP1S => RAIN_ICE_PARAM%XDRYINTP1S
  XDRYINTP2S => RAIN_ICE_PARAM%XDRYINTP2S
  XDRYINTP1G => RAIN_ICE_PARAM%XDRYINTP1G
  XDRYINTP2G => RAIN_ICE_PARAM%XDRYINTP2G
  XFSEDH => RAIN_ICE_PARAM%XFSEDH
  XEXSEDH => RAIN_ICE_PARAM%XEXSEDH
  X0DEPH => RAIN_ICE_PARAM%X0DEPH
  X1DEPH => RAIN_ICE_PARAM%X1DEPH
  XEX0DEPH => RAIN_ICE_PARAM%XEX0DEPH
  XEX1DEPH => RAIN_ICE_PARAM%XEX1DEPH
  XCOLIH => RAIN_ICE_PARAM%XCOLIH
  XCOLEXIH => RAIN_ICE_PARAM%XCOLEXIH
  XCOLSH => RAIN_ICE_PARAM%XCOLSH
  XCOLEXSH => RAIN_ICE_PARAM%XCOLEXSH
  XCOLGH => RAIN_ICE_PARAM%XCOLGH
  XCOLEXGH => RAIN_ICE_PARAM%XCOLEXGH
  XFWETH => RAIN_ICE_PARAM%XFWETH
  XFSWETH => RAIN_ICE_PARAM%XFSWETH
  XLBSWETH1 => RAIN_ICE_PARAM%XLBSWETH1
  XLBSWETH2 => RAIN_ICE_PARAM%XLBSWETH2
  XLBSWETH3 => RAIN_ICE_PARAM%XLBSWETH3
  XFGWETH => RAIN_ICE_PARAM%XFGWETH
  XLBGWETH1 => RAIN_ICE_PARAM%XLBGWETH1
  XLBGWETH2 => RAIN_ICE_PARAM%XLBGWETH2
  XLBGWETH3 => RAIN_ICE_PARAM%XLBGWETH3
  XFRWETH => RAIN_ICE_PARAM%XFRWETH
  XLBRWETH1 => RAIN_ICE_PARAM%XLBRWETH1
  XLBRWETH2 => RAIN_ICE_PARAM%XLBRWETH2
  XLBRWETH3 => RAIN_ICE_PARAM%XLBRWETH3
  XWETLBDAS_MIN => RAIN_ICE_PARAM%XWETLBDAS_MIN
  XWETLBDAS_MAX => RAIN_ICE_PARAM%XWETLBDAS_MAX
  XWETLBDAG_MIN => RAIN_ICE_PARAM%XWETLBDAG_MIN
  XWETLBDAG_MAX => RAIN_ICE_PARAM%XWETLBDAG_MAX
  XWETLBDAR_MIN => RAIN_ICE_PARAM%XWETLBDAR_MIN
  XWETLBDAR_MAX => RAIN_ICE_PARAM%XWETLBDAR_MAX
  XWETLBDAH_MIN => RAIN_ICE_PARAM%XWETLBDAH_MIN
  XWETLBDAH_MAX => RAIN_ICE_PARAM%XWETLBDAH_MAX
  XWETINTP1S => RAIN_ICE_PARAM%XWETINTP1S
  XWETINTP2S => RAIN_ICE_PARAM%XWETINTP2S
  XWETINTP1G => RAIN_ICE_PARAM%XWETINTP1G
  XWETINTP2G => RAIN_ICE_PARAM%XWETINTP2G
  XWETINTP1R => RAIN_ICE_PARAM%XWETINTP1R
  XWETINTP2R => RAIN_ICE_PARAM%XWETINTP2R
  XWETINTP1H => RAIN_ICE_PARAM%XWETINTP1H
  XWETINTP2H => RAIN_ICE_PARAM%XWETINTP2H
  !
  NGAMINC => RAIN_ICE_PARAM%NGAMINC
  NACCLBDAS => RAIN_ICE_PARAM%NACCLBDAS
  NACCLBDAR => RAIN_ICE_PARAM%NACCLBDAR
  NDRYLBDAR => RAIN_ICE_PARAM%NDRYLBDAR
  NDRYLBDAS => RAIN_ICE_PARAM%NDRYLBDAS
  NDRYLBDAG => RAIN_ICE_PARAM%NDRYLBDAG
  NWETLBDAS => RAIN_ICE_PARAM%NWETLBDAS
  NWETLBDAG => RAIN_ICE_PARAM%NWETLBDAG
  NWETLBDAR => RAIN_ICE_PARAM%NWETLBDAR
  NWETLBDAH => RAIN_ICE_PARAM%NWETLBDAH
END SUBROUTINE RAIN_ICE_PARAM_ASSOCIATE
!
SUBROUTINE RAIN_ICE_PARAM_ALLOCATE(HNAME, KDIM1, KDIM2)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: HNAME
  INTEGER, INTENT(IN)           :: KDIM1
  INTEGER, OPTIONAL, INTENT(IN) :: KDIM2

  SELECT CASE(TRIM(HNAME))
    !1D arrays
    CASE('XGAMINC_RIM1')
      ALLOCATE(RAIN_ICE_PARAM%XGAMINC_RIM1(KDIM1))
      XGAMINC_RIM1 => RAIN_ICE_PARAM%XGAMINC_RIM1
    CASE('XGAMINC_RIM2')
      ALLOCATE(RAIN_ICE_PARAM%XGAMINC_RIM2(KDIM1))
      XGAMINC_RIM2 => RAIN_ICE_PARAM%XGAMINC_RIM2
    CASE('XGAMINC_RIM4')
      ALLOCATE(RAIN_ICE_PARAM%XGAMINC_RIM4(KDIM1))
      XGAMINC_RIM4 => RAIN_ICE_PARAM%XGAMINC_RIM4
    !
    !2D arrays
    CASE('XKER_RACCSS')
      ALLOCATE(RAIN_ICE_PARAM%XKER_RACCSS(KDIM1, KDIM2))
      XKER_RACCSS=> RAIN_ICE_PARAM%XKER_RACCSS
    CASE('XKER_RACCS')
      ALLOCATE(RAIN_ICE_PARAM%XKER_RACCS(KDIM1, KDIM2))
      XKER_RACCS=> RAIN_ICE_PARAM%XKER_RACCS
    CASE('XKER_SACCRG')
      ALLOCATE(RAIN_ICE_PARAM%XKER_SACCRG(KDIM1, KDIM2))
      XKER_SACCRG=> RAIN_ICE_PARAM%XKER_SACCRG
    CASE('XKER_SDRYG')
      ALLOCATE(RAIN_ICE_PARAM%XKER_SDRYG(KDIM1, KDIM2))
      XKER_SDRYG=> RAIN_ICE_PARAM%XKER_SDRYG
    CASE('XKER_RDRYG')
      ALLOCATE(RAIN_ICE_PARAM%XKER_RDRYG(KDIM1, KDIM2))
      XKER_RDRYG=> RAIN_ICE_PARAM%XKER_RDRYG
    CASE('XKER_SWETH')
      ALLOCATE(RAIN_ICE_PARAM%XKER_SWETH(KDIM1, KDIM2))
      XKER_SWETH=> RAIN_ICE_PARAM%XKER_SWETH
    CASE('XKER_GWETH')
      ALLOCATE(RAIN_ICE_PARAM%XKER_GWETH(KDIM1, KDIM2))
      XKER_GWETH=> RAIN_ICE_PARAM%XKER_GWETH
    CASE('XKER_RWETH')
      ALLOCATE(RAIN_ICE_PARAM%XKER_RWETH(KDIM1, KDIM2))
      XKER_RWETH=> RAIN_ICE_PARAM%XKER_RWETH
  END SELECT
END SUBROUTINE RAIN_ICE_PARAM_ALLOCATE
SUBROUTINE RAIN_ICE_PARAM_DEALLOCATE()
  IMPLICIT NONE
  XGAMINC_RIM1=>NULL()
  DEALLOCATE(RAIN_ICE_PARAM%XGAMINC_RIM1)
  XGAMINC_RIM2=>NULL()
  DEALLOCATE(RAIN_ICE_PARAM%XGAMINC_RIM2)
  XKER_RACCSS=>NULL()
  DEALLOCATE(RAIN_ICE_PARAM%XKER_RACCSS)
  XKER_RACCS=>NULL()
  DEALLOCATE(RAIN_ICE_PARAM%XKER_RACCS)
  XKER_SACCRG=>NULL()
  DEALLOCATE(RAIN_ICE_PARAM%XKER_SACCRG)
  XKER_SDRYG=>NULL()
  DEALLOCATE(RAIN_ICE_PARAM%XKER_SDRYG)
  XKER_RDRYG=>NULL()
  DEALLOCATE(RAIN_ICE_PARAM%XKER_RDRYG)
END SUBROUTINE RAIN_ICE_PARAM_DEALLOCATE
END MODULE MODD_RAIN_ICE_PARAM
