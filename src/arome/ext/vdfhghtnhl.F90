!OPTIONS XOPT(HSFUN)
SUBROUTINE VDFHGHTNHL (YDVDF,YDEPHLI,YDECUMF,YDEPHY,YDPARAR,KIDIA   , KFDIA   , KLON    , KLEV   , KDRAFT, PTMST, KSTEP, &
                   & PUM1    , PVM1    , PTM1    , PQM1    , PLM1    , PIM1   , PAM1,&
                   & PAPHM1  , PAPM1   , PGEOM1  , PGEOH  , &
                   & PKMFL   , PKHFL   , PKQFL   , PMFLX  , PEXNF , PEXNH, &
! DIAGNOSTIC OUTPUT

                   & PUUH    , PVUH    , PSLGUH  , PQTUH  , PTHTVUH, PFRACB, &
                   & PZPTOP  , KPTOP   , PZPLCL  , KPLCL  , KPLZB   , &
                   &  PRICUI  , &
                   & PFPLVL  , PFPLVN  , PCLFR, &
                   & PBIR    , LDNODECP, LDRUNDRY, KPBLTYPE, &
                   & YSPP_CLDDPTH,YSPP_CLDDPTHDP, &
                   & YSPP_RFAC_TWOC,YSPP_RZC_H,YSPP_RZL_INF, &
                   & ZLENGTH_M, ZLENGTH_H,  PTKE)

!     ------------------------------------------------------------------

!**   * VDFHGHTNHL* - DETERMINES THE PBL-HEIGHT AND STRONG UPDRAFT FIELDS
!                  USING A ENTRAINING PARCEL ASCENT METHOD.

!     A.P. SIEBESMA                30/06/1999  Original (dry)
!     M. Ko"hler                    3/12/2004  Moist Version
!     Roel Neggers                 12/04/2005  Multiple updraft extension
!     Wim de Rooy/Geert Lenderink  13/06/2008 and 21/09/2010 Updates to combine TKE turbulence
!                                              with dual updraft EDMF. Lateral mixing according
!                                              to de Rooy & Siebesma MWR 2008 and QJRMS 2010
!     Wim de Rooy                  July /2015  Implementation LHARATU in Harmonie
!     Lisa Bengtsson               Feb /2017   Introduce LTOTPREC option
!     Karl-Ivar Ivarsson           Feb /2018   Code optimation
!     R. El Khatib 30-Apr-2019                 bugfix
!     Wim de Rooy                  June / 2019 Modifications among which energy
!                                              energy cascade term
!     R. El Khatib 27-Aug-2019                 Cleaning
!     Karl-Ivar Ivarsson           April /2020 Introduce LTOTPRECL option
!     U. Andrae Dec 2020                       Introduce SPP for HARMONIE-AROME
!      R. El Khatib 08-Jul-2022 Contribution to the encapsulation of YOMCST and YOETHF


!     PURPOSE
!     -------

!     DETERMINE PBL HEIGHT AND UPDRAFT FIELDS

!     INTERFACE
!     ---------

!     * VDFHGHTNHL* IS CALLED BY *VDFHGHTHL*

!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):

!     *KIDIA*        START POINT
!     *KFDIA*        END POINT
!     *KLEV*         NUMBER OF LEVELS
!     *KLON*         NUMBER OF GRID POINTS PER PACKET
!     *KDRAFT*       NUMBER OF EXPLICITLY MODELED DRAFTS - CURRENTLY 3:
!                    1: test parcel
!                    2: rising dry thermals which stop at cloud base or inversion
!                    3: rising dry thermals which become cloudy
!                    (4: downdrafts .. to be done?)

!     INPUT PARAMETERS (REAL):

!     *PTMST*        DOUBLE TIME STEP (SINGLE AT 1TH STEP)        S
!     *PUM1*         X-VELOCITY COMPONENT AT T-1                  M/S
!     *PVM1*         Y-VELOCITY COMPONENT AT T-1                  M/S
!     *PTM1*         TEMPERATURE AT T-1                           K
!     *PQM1*         SPECIFIC HUMUDITY AT T-1                     KG/KG
!     *PLM1*         SPECIFIC CLOUD LIQUID WATER AT T-1           KG/KG
!     *PIM1*         SPECIFIC CLOUD ICE AT T-1                    KG/KG
!     *PAM1*         CLOUD FRACTION AT T-1                        KG/KG
!     *PAPHM1*       PRESSURE AT HALF LEVEL AT T-1                PA
!     *PAPM1*        PRESSURE AT FULL LEVEL AT T-1                PA
!     *PGEOM1*       GEOPOTENTIAL AT T-1                          M2/S2
!     *PGEOH*        GEOPOTENTIAL AT HALF LEVEL                   M2/S2
!     *PKMFL*        SURFACE KINEMATIC MOMENTUM FLUX              M2/S2
!     *PKHFL*        SURFACE KINEMATIC HEAT FLUX                  K*M/S
!     *PKQFL*        SURFACE KINEMATIC MOISTURE FLUX              M/S
!     *PBIR*         BUOYANCY-FLUX INTEGRAL RATIO (-N/P)
!                    USED FOR DECOUPLING CRITERIA
!     *PEXNF*        ENXNER FUNCTION FOR FULL PRESSURE LEVELS (FOR OPTIMATION OF CODE)
!     *PEXNH*        ENXNER FUNCTION FOR HALF PRESSURE LEVELS (FOR OPTIMATION OF CODE),
!                    =  (P / PREF) ** R/CP

!     INPUT PARAMETERS (LOGICAL):

!     *LDNODECP*     TRUE:  NEVER DECOUPLE
!                    FALSE: MAYBE DECOUPLE
!     *LDRUNDRY*     TRUE:  RUN PARCEL WITHOUT CONDENSATION
!                    FALSE: RUN PARCEL WITH CONDENSATION

!     OUTPUT PARAMETERS (REAL):

!     *PFPLVL*       PBL PRECIPITATION FLUX AS RAIN                KG/(M**2*S)
!     *PFPLVN*       PBL PRECIPITATION FLUX AS SNOW                KG/(M**2*S)

!     *PUUH*         UPDRAFT X-MOMENTUM
!     *PVUH*         UPDRAFT Y-MOMENTUM
!     *PSLGUH*       UPDRAFT GENERALIZED LIQUID STATIC ENERGY (SLG)
!                    AT HALF LEVEL                                M2/S2
!     *PQTUH*        UPDRAFT SPECIFIC TOTAL WATER AT HALF LEVEL   KG/KG
!     *PTHTVUH*      UPDRAFT virt potential temp at HALF LEVEL (for TKE)  K
!     *PMFLX*        PBL MASS FLUX                                M/S
!     *PZPLCL*       HEIGHT OF LIFTING CONDENSATION LEVEL OF UPDRAFT          M
!     *PZPTOP*       HEIGHT OF LEVEL OF ZERO KINETIC ENERGY (W=0) OF UPDRAFT  M
!cstep/GL
!     *PBUOY_COR*    STABILITY CORRECTION PARAMETER TO BE USED FOR TKE SCHEME
!
!     *PWU*          VERTICAL VELOCITY OF SECOND UPDRAFT
!cstep/GL
!
!     OUTPUT PARAMETERS (INTEGER):

!     *KPLCL*         FIRST HALF LEVEL ABOVE REAL HEIGHT OF UPRAFT LCL
!     *KPTOP*         HIGHEST HALF LEVEL BELOW PZTOP, AND
!                       UPDRAFT TOP FULL LEVEL (PZTOP IS WITHIN THAT LAYER)
!     *KPLZB*         LEVEL OF UPRAFT ZERO BUOYANCY (LAST FULL LEVEL THAT IS POS. BUOYANT)
!     *KPBLTYPE*    -1: not defined yet
!                    0: stable PBL
!                    1: dry convective PBL (no cloud below parcel top)
!                    2: stratocumulus
!                    3: shallow cumulus
!                    4: deep cumulus

!     METHOD
!     ------

!     SEE DOCUMENTATION

!     ------------------------------------------------------------------

USE YOEPHLI  , ONLY : TEPHLI
USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE YOMCST   , ONLY : YDCST=>YRCST ! allows use of included functions. REK.
USE YOETHF   , ONLY : YDTHF=>YRTHF ! allows use of included functions. REK.
USE PARPHY   , ONLY : RKAP
USE YOECUMF  , ONLY : TECUMF
USE YOMPARAR , ONLY : TPARAR
USE YOEPHY   , ONLY : TEPHY
USE YOEVDF   , ONLY : TVDF

!for optimation
USE MODD_CST
USE MODD_RAIN_ICE_DESCR_n
USE MODD_RAIN_ICE_PARAM_n
USE MODE_TIWMX_TAB
USE MODE_TIWMX

USE SPP_MOD_TYPE, ONLY : TSPP_CONFIG_TYPE, APPLY_SPP 

IMPLICIT NONE


!*         0.1    GLOBAL VARIABLES

TYPE(TVDF)        ,INTENT(IN)    :: YDVDF
TYPE(TECUMF)      ,INTENT(IN)    :: YDECUMF
TYPE(TEPHLI)      ,INTENT(IN)    :: YDEPHLI
TYPE(TEPHY)       ,INTENT(IN)    :: YDEPHY
TYPE(TPARAR)      ,INTENT(IN)    :: YDPARAR
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV
INTEGER(KIND=JPIM),INTENT(IN)    :: KDRAFT
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTEP
INTEGER(KIND=JPIM),INTENT(INOUT) :: KPLCL(KLON,KDRAFT)
INTEGER(KIND=JPIM),INTENT(INOUT) :: KPTOP(KLON,KDRAFT)
INTEGER(KIND=JPIM),INTENT(INOUT) :: KPLZB(KLON,KDRAFT)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTMST
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUM1(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVM1(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTM1(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQM1(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLM1(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PIM1(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAM1(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHM1(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPM1(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOM1(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOH(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKMFL(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKHFL(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKQFL(KLON)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PMFLX(KLON,0:KLEV,KDRAFT)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PUUH(KLON,0:KLEV,KDRAFT)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVUH(KLON,0:KLEV,KDRAFT)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSLGUH(KLON,0:KLEV,KDRAFT)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PQTUH(KLON,0:KLEV,KDRAFT)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTHTVUH(KLON,0:KLEV,KDRAFT)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFRACB(KLON,KDRAFT)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PZPLCL(KLON,KDRAFT)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PZPTOP(KLON,KDRAFT)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFPLVL(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFPLVN(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCLFR(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PBIR(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRICUI(KLON)
TYPE(TSPP_CONFIG_TYPE),INTENT(INOUT) :: YSPP_CLDDPTH,  YSPP_CLDDPTHDP, &
                                      & YSPP_RFAC_TWOC,YSPP_RZC_H,YSPP_RZL_INF
! variables RACMO turbulence scheme
REAL(KIND=JPRB)   ,INTENT(OUT) :: ZLENGTH_M(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT) :: ZLENGTH_H(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTKE(KLON,KLEV)

! variables for optimation of code
REAL(KIND=JPRB)   ,INTENT(IN) :: PEXNF(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN) :: PEXNH(KLON,0:KLEV)

LOGICAL           ,INTENT(IN)    :: LDNODECP(KLON)
!ldrundry not used now
LOGICAL           ,INTENT(IN)    :: LDRUNDRY(KLON)
INTEGER(KIND=JPIM),INTENT(INOUT) :: KPBLTYPE(KLON)
REAL(KIND=JPRB)    :: ZENCASC(KLON,0:KLEV)

! --- variables associated with Lgeert
REAL(KIND=JPRB)    ::  PBUOY_COR (KLON,0:KLEV)
REAL(KIND=JPRB)    :: ZQCUH(KLON,0:KLEV,KDRAFT)
REAL(KIND=JPRB)    :: ZWU2H(KLON,0:KLEV,KDRAFT)

REAL(KIND=JPRB)    ::  PWU       (KLON,0:KLEV)
REAL(KIND=JPRB)    :: ZQSVAR(KLON,KLEV)
REAL(KIND=JPRB)    :: ZDQSDTEMP(KLON,KLEV)
REAL(KIND=JPRB)    :: ZFACW,ZFACI,ZESDP,ZCOR,ZLAT2CP,&
                      & ZESW,ZESI,ZES,ZQSAT,ZEL2R,ZEI2R,ZE2R
!cstep/GL ---------------------------------------------



!*         0.2    LOCAL VARIABLES

!--- mean & environmental properties ---
REAL(KIND=JPRB) ::    ZUSTAR (KLON)      , ZWSTAR(KLON)       , ZKHVFL(KLON)       , &
                    & ZWSIGMA(KLON)      , &
                    & ZSLGENH(KLON,0:KLEV),ZQLENH(KLON,0:KLEV), ZQIENH(KLON,0:KLEV), &
                    & ZQTENH(KLON,0:KLEV), ZUENH(KLON,0:KLEV) , ZVENH(KLON,0:KLEV) , &
                    & ZTVEN(KLON,KLEV),    ZQTM1 (KLON,KLEV)  , &
                    & ZSLGM1(KLON,KLEV)  ,&
                    & ZTENH(KLON,0:KLEV) , ZRHOH (KLON,0:KLEV), ZTHVEN(KLON,KLEV)

REAL(KIND=JPRB) ::    ZSTAR(KLON), ZCHICRIT(KLON,0:KLEV),ZMEANCHICRIT(KLON)   , &
                    & ZDETRSHALLOW(KLON), ZVARQ(KLON,0:KLEV), &
                    & ZFRACMB(KLON) , ZMINIMUM(KLON) ,  &
                    & ZQVENH(KLON,0:KLEV), &
                    & ZTHTVENH(KLON,0:KLEV), ZTHTLUH(KLON,0:KLEV), &
                    & ZTHTENH(KLON,0:KLEV), &
                    & ZTHTVUH(KLON,0:KLEV), ZTHTLEH(KLON,0:KLEV),&
                    & ZDQSDTU(KLON,0:KLEV),ZGAMMA(KLON,0:KLEV),ZDUMFUNC(KLON,0:KLEV), &
                    & ZPRES_0, ZKAPPA,ZQSATU(KLON,0:KLEV)


!--- updraft parameters ---
REAL(KIND=JPRB) ::    ZWUH,DZH, &
                    & ZQUH  (KLON,0:KLEV,KDRAFT), &
    & ZTUH  (KLON,0:KLEV,KDRAFT), ZEPS  (KLON,0:KLEV,KDRAFT), &
    & ZDETR (KLON,0:KLEV,KDRAFT), ZFRAC (KLON,0:KLEV,KDRAFT), &
    & ZBUOF (KLON,KLEV,KDRAFT)  , &
                    & ZDELTAMINEPS(KLON,0:KLEV), ZCAPE1(KLON)

REAL(KIND=JPRB) ::    ZQSATM, ZSATDEF, &
                    & ZUPFLXL(KLON,0:KLEV,KDRAFT), ZUPFLXN(KLON,0:KLEV,KDRAFT), &
                    & ZUPGENL(KLON,KLEV,KDRAFT), ZUPGENN(KLON,KLEV,KDRAFT), &
                    & ZDZRHO, ZPFLXTOT, ZPEVAPUP, ZFAC, ZUPMELT, ZUPMELTTEND

REAL(KIND=JPRB) ::    ZFRACB(KLON,KDRAFT), ZMFLXB(KLON,KDRAFT), ZTVEXCSURF(KLON,KDRAFT)


REAL(KIND=JPRB) ::    ZFRACMAX , ZFACMAXEXC , ZFRACTEST , ZFACTESTEXC , &
                    & ZFACEXC(KLON,KDRAFT), ZDUMFRAC, ZDUMR, TVEXCSURF, ZWT, ZWL, ZEL, ZET, &
                    & ZFACCASC(KLON,0:KLEV)

LOGICAL ::            LLDONE(KLON,KDRAFT)


INTEGER(KIND=JPIM) :: IS, JK, JL, JD, JKM

INTEGER(KIND=JPIM) :: IKSTAR(KLON)

REAL(KIND=JPRB) ::    ZQEXC   , ZTEXC   , ZDZ     , &
                    & ZCONS10 , ZTVMEAN     , &
                    & ZRG     , ZMFMAX  , ZMFS(KLON,KDRAFT)

!          REMAINING MODEL PARAMETERS

REAL(KIND=JPRB) ::    ZTAUEPS(KLON) , ZCLDDEPTH  , &
                    & ZW2THRESH         , ZSTABTHRESH       , ZBIRTHRESH , &
                    & ZCLDDEPTHDP       , ZDZCLOUD(KLON), &
                    & ZREPUST, ZGHM1
    
REAL(KIND=JPRB) ::    ZZI(KLON),ZB1(KLON),ZCOUNT(KLON), &
                    & ZFRMIN(KLON,2),ZMU,ZVAL
INTEGER(KIND=JPIM) :: ITOP, IBASE, JKO,JKE

!cstep 30082007: introduce two parcel time scales, one for the test parcel (ZTAUEPS_TEST) and one
!              : for the actual parcels (ZTAUEPS)
REAL(KIND=JPRB) :: ZTAUEPS_TEST

INTEGER(KIND=JPIM) :: IZI(KLON,KDRAFT)


REAL(KIND=JPRB) ::    ZHOOK_HANDLE


#include "surf_inq.h"

#include "vdfparcelhl.intfb.h"
#include "vdfpdftablehl.intfb.h"
#include "vdfexcuhl.intfb.h"
#include "fcttre.func.h"



!     ------------------------------------------------------------------

!*         1.     INITIALIZATION
!                 --------------

IF (LHOOK) CALL DR_HOOK('VDFHGHTNHL',0,ZHOOK_HANDLE)
ASSOCIATE(RTAUMEL=>YDECUMF%RTAUMEL, &
 & RG=>YDCST%RG, RCPD=>YDCST%RCPD, RETV=>YDCST%RETV, RLVTT=>YDCST%RLVTT, RLSTT=>YDCST%RLSTT, &
 & RATM=>YDCST%RATM, RTT=>YDCST%RTT, RLMLT=>YDCST%RLMLT, RD=>YDCST%RD, &
 & R2ES=>YDTHF%R2ES, R3LES=>YDTHF%R3LES, R3IES=>YDTHF%R3IES, R4LES=>YDTHF%R4LES, &
 & R4IES=>YDTHF%R4IES, R5LES=>YDTHF%R5LES, R5IES=>YDTHF%R5IES, R5ALVCP=>YDTHF%R5ALVCP, &
 & R5ALSCP=>YDTHF%R5ALSCP, RALVDCP=>YDTHF%RALVDCP, RALSDCP=>YDTHF%RALSDCP, &
 & RTWAT=>YDTHF%RTWAT, RTICE=>YDTHF%RTICE, RTICECU=>YDTHF%RTICECU, RTWAT_RTICE_R=>YDTHF%RTWAT_RTICE_R, &
 & RTWAT_RTICECU_R=>YDTHF%RTWAT_RTICECU_R, &
 & YSURF=>YDEPHY%YSURF,LHARATU=>YDPARAR%PHYEX%TURBN%LHARAT, &
 & LTOTPREC=>YDPARAR%LTOTPREC,LTOTPRECL=>YDPARAR%LTOTPRECL)

ZWL= 200._JPRB
ZWT= 400._JPRB
! typical entrainment values at LCL (L) and TOP (T)
ZEL= 0.002_JPRB
ZET= 0.002_JPRB
ZPRES_0 = 100000._JPRB     ! standard pressure consistent with vdfexcu
ZKAPPA = RD / RCPD         ! consistent with vdfexcu

ZFRACTEST   = 0.002_JPRB   ! top % of the PDF associated with the test parcel
CALL VDFPDFTABLEHL (ZFRACTEST, ZFACTESTEXC, ZDUMR, ZDUMR, 0) ! associated PDF scaling factor


ZFRACMAX    = 0.1_JPRB     ! total convective area fraction that is done with mass flux

CALL VDFPDFTABLEHL (ZFRACMAX, ZFACMAXEXC, ZDUMR, ZDUMR, 0) ! associated PDF scaling factor

! eddy turnover time scale used in parcel entrainment [s]  (Neggers, Siebesma & Jonker, JAS 2002)
ZTAUEPS_TEST = 400._JPRB   ! Roel's original ZTAUEPS value


!ZW2THRESH  = -1._JPRB     ! threshold parcel vertical velocity squared [m2/s2]
!CGL
ZW2THRESH   = 0.0_JPRB

ZCLDDEPTH   = 2000._JPRB   ! threshold cloud thickness for stcu/cu transition [m]

ZCLDDEPTHDP = 4000._JPRB ! threshold cloud thickness used in shallow/deep decision [m]

TVEXCSURF = 0.0_JPRB ! initialisation

IF(XFRMIN(19)>0.)ZCLDDEPTH = XFRMIN(19)
IF(XFRMIN(20)>0.)ZCLDDEPTHDP = XFRMIN(20)


ZSTABTHRESH = 20._JPRB     ! threshold stability (Klein & Hartmann criteria) [K]
ZBIRTHRESH  = 0.1_JPRB     ! threshold BIR (TKE decoupling criteria) [1]

CALL SURF_INQ(YSURF,PREPUST=ZREPUST)

! optimization
ZRG         = 1.0_JPRB/RG
ZLAT2CP     = RLVTT/RCPD
ZEL2R       = 0.62198_JPRB*RLVTT/RD
ZEI2R       = 0.62198_JPRB*RLSTT/RD

DO JL=KIDIA,KFDIA
   KPBLTYPE(JL)     = -1           ! -1 means: yet unknown

   ZZI(JL)          = 0._JPRB      ! mixed layer scalings
   ZWSTAR(JL)       = 0._JPRB

   PRICUI(JL)       = 1._JPRB      ! 1 / cumulus inversion Richardson number
   ZCAPE1(JL)       = 0._JPRB

ENDDO

DO JD=1,KDRAFT
   DO JL=KIDIA,KFDIA
      PZPLCL(JL,JD)  = -100._JPRB  ! default value: -100 (no LCL)
      PZPTOP(JL,JD)  = 0._JPRB
      KPLCL(JL,JD)   = 0           ! default value: 0 (no PBL cloud)
      KPTOP(JL,JD)   = 0
      KPLZB(JL,JD)   = 0
      LLDONE(JL,JD)  = .TRUE.      ! default: TRUE (don't launch the parcel)
      ZFRACB(JL,JD)  = 0._JPRB
      PFRACB(JL,JD)  = 0._JPRB
      ZFACEXC(JL,JD) = 0._JPRB
      ZMFLXB(JL,JD)  = 0._JPRB
      ZTVEXCSURF(JL,JD) = 0._JPRB
   ENDDO
ENDDO

DO JK=0,KLEV
   DO JL=KIDIA,KFDIA
      ZCHICRIT(JL,JK) = 0._JPRB
   ENDDO
ENDDO

!--- parcel half level parameters ---
DO JD=1,KDRAFT
   DO JK=0,KLEV
      DO JL=KIDIA,KFDIA
         PUUH(JL,JK,JD)    = 0.0_JPRB
         PVUH(JL,JK,JD)    = 0.0_JPRB
         PSLGUH(JL,JK,JD)  = 0.0_JPRB
         PQTUH(JL,JK,JD)   = 0.0_JPRB
         PTHTVUH(JL,JK,JD) = 0.0_JPRB
         PMFLX(JL,JK,JD)   = 0.0_JPRB
         ZTUH(JL,JK,JD)    = 0.0_JPRB
         ZQUH(JL,JK,JD)    = 0.0_JPRB
         ZQCUH(JL,JK,JD)   = 0.0_JPRB
         ZEPS(JL,JK,JD)    = 0.0_JPRB
         ZDETR(JL,JK,JD)   = 0.0_JPRB
         ZWU2H(JL,JK,JD)   = 0.0_JPRB
         ZFRAC(JL,JK,JD)   = 0.0_JPRB
         ZUPFLXL(JL,JK,JD) = 0.0_JPRB
         ZUPFLXN(JL,JK,JD) = 0.0_JPRB
         ZVARQ(JL,JK)       = 0.0_JPRB
         ZENCASC(JL,JK) =0.0_JPRB
      ENDDO
   ENDDO
ENDDO

!--- parcel full level parameters ---
DO JD=1,KDRAFT
   DO JK=1,KLEV
      DO JL=KIDIA,KFDIA
         ZBUOF(JL,JK,JD)    = 0.0_JPRB
         ZUPGENL(JL,JK,JD)  = 0.0_JPRB
         ZUPGENN(JL,JK,JD)  = 0.0_JPRB
      ENDDO
   ENDDO
ENDDO

! Setup SPP patterns 
IF (YSPP_CLDDPTH%LPERT) THEN
  CALL APPLY_SPP(YSPP_CLDDPTH, &
               & KLON,KIDIA,KFDIA, &
               & ZCLDDEPTH,ZFRMIN(:,1))
ELSE
  DO JL=KIDIA,KFDIA
    ZFRMIN(JL,1) = ZCLDDEPTH
  ENDDO
ENDIF

IF (YSPP_CLDDPTHDP%LPERT) THEN
  CALL APPLY_SPP(YSPP_CLDDPTHDP, &
               & KLON,KIDIA,KFDIA, &
               & ZCLDDEPTHDP,ZFRMIN(:,2))
ELSE
  DO JL=KIDIA,KFDIA
    ZFRMIN(JL,2) = ZCLDDEPTHDP
  ENDDO
ENDIF



!     -----------------------------------------------------------------

!*         2.     PREPARE FIELDS ON HALF LEVELS BY LINEAR INTERPOLATION
!*                OF CONSERVED VARIABLES
!                 -----------------------------------------------------

!*         2.1  full level cpm, slg, qt and Tv
!*

DO JK=1,KLEV
   DO JL=KIDIA,KFDIA
      ZSLGM1(JL,JK) = RCPD * PTM1(JL,JK) + PGEOM1(JL,JK)&
           & - RLVTT * PLM1(JL,JK) - RLSTT * PIM1(JL,JK)
      ZQTM1 (JL,JK) = PQM1(JL,JK) + PLM1(JL,JK) + PIM1(JL,JK)

      !          parcel goes through cloud portion of environment
      !          (added ql loading; ql,cld=ql,mean/fc; qv = qsat)
      !          safety: fc>0.1; linear interpolation between overcast
      !                  and cloudy portion for 0<fc<0.1
      !                  guaranteed to be < tv from mean conditions

      !          grid box mean virtual effect
      ZTVMEAN       = PTM1(JL,JK) * ( 1.0_JPRB + RETV * PQM1(JL,JK)&
           & - PLM1(JL,JK) - PIM1(JL,JK) )       !qli loading
      ZTVEN(JL,JK) = ZTVMEAN
      ZTHVEN(JL,JK) = ZTVEN(JL,JK) / PEXNF(JL,JK) !( PAPM1(JL,JK)/RATM )**(-RD/RCPD) * ZTVEN(JL,JK)
   ENDDO
ENDDO


!*         2.2  half-level environment interpolation (qt, ql, qi, slg)
!*              attention:  not good to interpolate everything independently
!*              better:     interpolate conserved variables and derive rest!!!
!*

DO JK=1,KLEV-1
   DO JL=KIDIA,KFDIA

      IF (JK==1) THEN
         ZGHM1 = PGEOH(JL,JK) + 50000._JPRB*RG   !avoid using top half level (=inf)
      ELSE
         ZGHM1 = PGEOH(JL,JK-1)
      ENDIF

      ZQTENH(JL,JK) = ( ZQTM1(JL,JK+1) *(ZGHM1-PGEOH(JL,JK  )) &
                  & +   ZQTM1(JL,JK)   *(PGEOH(JL,JK  )-PGEOH(JL,JK+1)) &
                  &   )                /(ZGHM1-PGEOH(JL,JK+1))
      ZQLENH(JL,JK) = ( PLM1(JL,JK+1)  *(ZGHM1-PGEOH(JL,JK  )) &
                  & +   PLM1(JL,JK)    *(PGEOH(JL,JK  )-PGEOH(JL,JK+1)) &
                  &   )                /(ZGHM1-PGEOH(JL,JK+1))
      ZQIENH(JL,JK) = ( PIM1(JL,JK+1)  *(ZGHM1-PGEOH(JL,JK  )) &
                  & +   PIM1(JL,JK)    *(PGEOH(JL,JK  )-PGEOH(JL,JK+1)) &
                  &   )                /(ZGHM1-PGEOH(JL,JK+1))

      ZQVENH(JL,JK) = ZQTENH(JL,JK) - ZQLENH(JL,JK) - ZQIENH(JL,JK)

      ZSLGENH(JL,JK)= ( ZSLGM1(JL,JK+1)*(ZGHM1-PGEOH(JL,JK  )) &
                  & +   ZSLGM1(JL,JK)  *(PGEOH(JL,JK  )-PGEOH(JL,JK+1)) &
                  &   )                /(ZGHM1-PGEOH(JL,JK+1))
      ZUENH(JL,JK)  = ( PUM1(JL,JK+1)  *(ZGHM1-PGEOH(JL,JK  )) &
                  & +   PUM1(JL,JK)    *(PGEOH(JL,JK  )-PGEOH(JL,JK+1)) &
                  &   )                /(ZGHM1-PGEOH(JL,JK+1))
      ZVENH(JL,JK)  = ( PVM1(JL,JK+1)  *(ZGHM1-PGEOH(JL,JK  )) &
                  & +   PVM1(JL,JK)    *(PGEOH(JL,JK  )-PGEOH(JL,JK+1)) &
                  &   )                /(ZGHM1-PGEOH(JL,JK+1))


      !  Calculate T at half levels from sl, for later use in density calculations
      ZTENH(JL,JK)  =  ( PTM1(JL,JK+1)  *(ZGHM1-PGEOH(JL,JK  )) &
                   & +   PTM1(JL,JK)    *(PGEOH(JL,JK  )-PGEOH(JL,JK+1)) &
                   &   )                /(ZGHM1-PGEOH(JL,JK+1))

      !  Determine thetav environment
      ZTHTENH(JL,JK) = (ZTENH(JL,JK)/PEXNH(JL,JK))
      ZTHTVENH(JL,JK) = (ZTENH(JL,JK)/PEXNH(JL,JK)) &
           &     * (1._JPRB + RETV * ZQVENH(JL,JK) &
                                ! add qice correction
           &      - ZQLENH(JL,JK) )

      ! initialize updraft thetav for dry and moist updraft with environment values
      PTHTVUH(JL,JK,1)=ZTHTVENH(JL,JK)
      PTHTVUH(JL,JK,2)=ZTHTVENH(JL,JK)
      PTHTVUH(JL,JK,3)=ZTHTVENH(JL,JK)
      PQTUH(JL,JK,1)=ZQTENH(JL,JK)
      PQTUH(JL,JK,2)=ZQTENH(JL,JK)
      PQTUH(JL,JK,3)=ZQTENH(JL,JK)
      PUUH(JL,JK,1)=ZUENH(JL,JK)
      PUUH(JL,JK,2)=ZUENH(JL,JK)
      PUUH(JL,JK,3)=ZUENH(JL,JK)
      PVUH(JL,JK,1)=ZVENH(JL,JK)
      PVUH(JL,JK,2)=ZVENH(JL,JK)
      PVUH(JL,JK,3)=ZVENH(JL,JK)
      PSLGUH(JL,JK,1)=ZSLGENH(JL,JK)
      PSLGUH(JL,JK,2)=ZSLGENH(JL,JK)
      PSLGUH(JL,JK,3)=ZSLGENH(JL,JK)


      ZRHOH(JL,JK) = PAPHM1(JL,JK)/(RD*ZTENH(JL,JK))

   ENDDO
ENDDO

!
!  Initialization lowest updraft level with environment like values
!
DO JD=1,KDRAFT
   DO JL=KIDIA,KFDIA
      PUUH(JL,KLEV,JD)  = PUUH(JL,KLEV-1,JD)
      PVUH(JL,KLEV,JD)  =  PVUH(JL,KLEV-1,JD)
      PSLGUH(JL,KLEV,JD)= PSLGUH(JL,KLEV-1,JD)
      PQTUH(JL,KLEV,JD) = PQTUH(JL,KLEV-1,JD)
      PTHTVUH(JL,KLEV,JD) =  PTHTVUH(JL,KLEV-1,JD)
   ENDDO
ENDDO



!     -----------------------------------------------------------------
!*         3.     RELEASE THE FIRST (TEST) UPDRAFT TO GET PBL HEIGHTS


!* set updraft index to 1
JD = 1

DO JL=KIDIA,KFDIA

   ZFRACB(JL,JD) = ZFRACMAX  !CGL replaced ZFRACB(JL,JD) = ZFRACTEST

   !* 3.1    Determine stability of BL using the surface buoyancy flux

   ZKHVFL(JL)  = ( 1.0_JPRB + RETV *  ZQTM1(JL,KLEV) ) * PKHFL(JL) +&
        & ( RETV * ZSLGM1(JL,KLEV) / RCPD )     * PKQFL(JL)

   IF ( ZKHVFL(JL) >= 0.0_JPRB ) THEN

      ! stable BL (no updrafts expected/needed)
      KPBLTYPE(JL)  = 0

   ELSE

      LLDONE(JL,JD) = .FALSE.  !confirm launch

      !* 3.2    Sigma-w-L60 (ignore 1-z/zi term)

      ZUSTAR  (JL)  = MAX( SQRT(PKMFL(JL)), ZREPUST )      ! u* (repust=10e-4

      ZWSIGMA(JL)      = 1.2_JPRB&
           & * ( ZUSTAR(JL)**3&
           & - 1.5_JPRB * RKAP * ZKHVFL(JL) * (PGEOH(JL,KLEV-1)-PGEOH(JL,KLEV))&
           &  / PTM1(JL,KLEV-1)&
           & ) ** ( 1.0_JPRB/3._JPRB )                         ! Kolmogorov 1/3-power



      !* 3.3    Initialize updraft

      !get the constant associated with the top ZFRACTEST % of the PDF
      ZFACEXC(JL,1) = ZFACTESTEXC
      !calculate the initial excess values
      ZWU2H(JL,KLEV-1,JD) = ( ZFACEXC(JL,1) * ZWSIGMA(JL) )**2
      ZTEXC               = - ZFACEXC(JL,1) * PKHFL(JL) / ZWSIGMA(JL)
      ZQEXC               = - ZFACEXC(JL,1) * PKQFL(JL) / ZWSIGMA(JL)
      ZTEXC            = MAX(ZTEXC, 0.0_JPRB)
      ZQEXC            = MAX(ZQEXC, 0.0_JPRB)
      PQTUH(JL,KLEV-1,JD) = ZQTENH(JL,KLEV-1) + ZQEXC
      ZQCUH(JL,KLEV-1,JD) = ZQLENH(JL,KLEV-1) + ZQIENH(JL,KLEV-1)
      ZQUH (JL,KLEV-1,JD) = PQTUH(JL,KLEV-1,JD)  - ZQCUH(JL,KLEV-1,JD)
      PSLGUH(JL,KLEV-1,JD)= ZSLGENH(JL,KLEV-1) + RCPD * ZTEXC
      ZTUH (JL,KLEV-1,JD) = ( PSLGUH (JL,KLEV-1,JD) - PGEOH(JL,KLEV-1)&
           & + RLVTT*ZQLENH(JL,KLEV-1) + RLSTT*ZQIENH(JL,KLEV-1)&
           & ) / RCPD
      PUUH(JL,KLEV-1,JD)= ZUENH(JL,KLEV-1)
      PVUH(JL,KLEV-1,JD)= ZVENH(JL,KLEV-1)

   ENDIF
ENDDO !JL


!* 3.4   Release the test updraft #1
!*          - Mainly used to get a first guess of the heights of cloud base & inversion,
!*            and to determine PBL type accordingly.

ZTAUEPS(:) = 400._JPRB

!
CALL VDFPARCELHL(YDEPHLI,YDPARAR,KIDIA,KFDIA,KLON,KLEV,KDRAFT,PGEOH,PGEOM1,PAPHM1,PUM1,PVM1,ZQTM1,ZSLGM1,ZTVEN,PUUH,PVUH,&
 & PSLGUH,PQTUH,ZWU2H,ZQCUH,ZBUOF,ZQUH,ZTUH,ZEPS,PZPLCL,KPLCL,PZPTOP,KPTOP,KPLZB,JD,ZUPGENL,ZUPGENN,ZTAUEPS,ZW2THRESH,LLDONE,KPBLTYPE)



!     -----------------------------------------------------------------
!*         4.     CLASSIFICATION OF THE CONVECTIVE PBL
!                 ------------------------------------


!* 4.1    Classify the convective PBL
!*


DO JL=KIDIA,KFDIA
   IF ( KPBLTYPE(JL)/=0 ) THEN

      !CGL loose criterium by 100 m
      IF ( PZPLCL(JL,1) > (PZPTOP(JL,1)+100._JPRB) .OR. KPLCL(JL,1) == 0 ) THEN

         !dry convective PBL
         KPBLTYPE(JL)  = 1                   !dry convective PBL
         ZDZCLOUD(JL)  = 0.0_JPRB            !cloud thickness

      ELSE

         !moist convective PBL
         ZDZCLOUD(JL)  = PZPTOP(JL,1) - PZPLCL(JL,1) !cloud thickness
         IF (ZDZCLOUD(JL)>ZFRMIN(JL,2)) THEN

            !deep convection
            KPBLTYPE(JL) = 4
         ELSE

!wc no special stratocumulus regime anymore
!           KPBLTYPE(JL) = 2   !set the type to stratocumulus for the moment
            KPBLTYPE(JL) = 3   !set the type to stratocumulus for the moment

         ENDIF

      ENDIF

   ENDIF !KPBLTYPE /=0
ENDDO !JL



!wc no special stratocumulus regime anymore
!* 4.2    Check the stratocumulus/shallow cumulus criterion (trigger function)
!*        If shallow cumulus is diagnosed, KPBLTYPE will be set to 3
!*
!CALL VDFSTCUCRITHL ( KIDIA , KFDIA  , KLON  , KLEV , KDRAFT ,&
!     &    PTM1  , ZSLGM1 , ZQTM1 , PAPM1 ,&
!     &    ZSTABTHRESH, ZCLDDEPTH, ZBIRTHRESH, ZDZCLOUD,&
!     &    KPTOP , KPBLTYPE, LDNODECP)


!     -----------------------------------------------------------------

!*         5.     CLOSURE FOR ORGANIZED UPDRAFTS (JD=2,3)
!                 ---------------------------------------


!* 5.1    Determine some mixed layer scalings
!*

DO JL=KIDIA,KFDIA

   IF ( KPBLTYPE(JL)/=0 ) THEN    !don't do this for stable PBL

      SELECT CASE (KPBLTYPE(JL))

      CASE(1)
         !Dry convective PBL - Inversion height
         ZZI(JL)   = PZPTOP(JL,1)

      CASE(2)
         !Stratocumulus - Inversion height
         !CAUTION: During decoupling in the intermediate regime (e.g. ASTEX/ATEX) the
         !   relevant ML scaling height changes from PBL inversion to level of minimum
         !   buoyancy flux. In the current setup this is not modelled yet!
         ZZI(JL)   = PZPTOP(JL,1)

      CASE(3)
         !Shallow cumulus - Level of minimum buoyancy flux
         !Assume that the moist updraft LCL is very close to this level
         ZZI(JL)   = PZPLCL(JL,1)

      CASE(4)
         !Deep cumulus - Only do a dry parcel up to cloud base
         ZZI(JL)   = PZPLCL(JL,1)

      END SELECT

      !--- Mixed layer convective velocity scale ---
      ZWSTAR(JL) = ( -ZKHVFL(JL) * RG * ZZI(JL) / ZTHVEN(JL,KLEV)  ) ** (1._JPRB/3._JPRB)
      ! CGL for the moment revert back to old constant time scale
      ZTAUEPS(JL) = 400._JPRB


   ENDIF

ENDDO


!* 5.3    Closure of updraft area fractions (JD=2,3)
!*

DO JL=KIDIA,KFDIA

   IF ( KPBLTYPE(JL)/=0 ) THEN    !don't do this for stable PBL

      SELECT CASE (KPBLTYPE(JL))

        CASE(1)
          !Dry convective PBL
          ZFRACB(JL,3) = 0._JPRB
          ZFRACB(JL,2) = ZFRACMAX - ZFRACB(JL,3)
          !          ZFRACB(JL,2) = ZFRACB(JL,2)*(1.-EXP(-ZZI(JL)/400._JPRB))
        CASE(2)
          !Stratocumulus
          ZFRACB(JL,3) = 0.1_JPRB
          ZFRACB(JL,2) = ZFRACMAX - ZFRACB(JL,3)
        CASE(3)
          !Shallow cumulus
          ZFRACB(JL,3)  = 0.03_JPRB
          ZFRACB(JL,2) = ZFRACMAX - ZFRACB(JL,3)
          !          ZFRACB(JL,2) = ZFRACB(JL,2)*(1.-EXP(-ZZI(JL)/400._JPRB) )

        CASE(4)
          !Deep cumulus
          ZFRACB(JL,3) = 0._JPRB
          ZFRACB(JL,2) = ZFRACMAX - ZFRACB(JL,3)

      END SELECT !KPBLTYPE


   ENDIF !KPBLTYPE /=0

ENDDO !JL



!     -----------------------------------------------------------------

!*         6.     CALCULATE VERTICAL PROFILES OF ALL UPDRAFTS (JD=2,3)
!                 ----------------------------------------------------


!*       6.1    Calculate the scaling factors of the updraft excess with the surface joint PDFs
!*
DO JD = 2,KDRAFT
   DO JL=KIDIA,KFDIA

      IF ( KPBLTYPE(JL)/=0 .AND. ZFRACB(JL,JD)>0._JPRB ) THEN

         !-- Get the PDF scaling factor --
         SELECT CASE (JD)

         CASE(2)
            !lower part of top ZFRACMAX %
            ZDUMFRAC = ZFRACMAX - ZFRACB(JL,2)
            CALL VDFPDFTABLEHL(ZDUMFRAC , ZFACEXC(JL,2), ZDUMR, ZDUMR, 0)
            ZFACEXC(JL,2) = ( ZFRACMAX * ZFACMAXEXC - ZDUMFRAC * ZFACEXC(JL,2) ) / ZFRACB(JL,2)
         CASE(3)
            !upper part of top ZFRACMAX %
            ZDUMFRAC = ZFRACB(JL,JD)
            CALL VDFPDFTABLEHL(ZDUMFRAC , ZFACEXC(JL,3), ZDUMR, ZDUMR, 0)

         END SELECT

      ENDIF !KPBLTYPE & ZFRACB

   ENDDO !JL
ENDDO !JD


!*       6.2    Vertical integration of dry & moist updraft budgets (JD=2,3)
!*
DO JD = 2,KDRAFT

   !-- Initialize updraft --
   DO JL=KIDIA,KFDIA

      IF ( KPBLTYPE(JL)/=0 .AND. ZFRACB(JL,JD)>0._JPRB ) THEN

         LLDONE(JL,JD) = .FALSE. !confirm launch

         ZTEXC            = - ZFACEXC(JL,JD) * PKHFL(JL) / ZWSIGMA(JL)
         ZQEXC            = - ZFACEXC(JL,JD) * PKQFL(JL) / ZWSIGMA(JL)

         ZWU2H(JL,KLEV-1,JD) = (ZWSIGMA(JL))**2

         ZTEXC            = MAX(ZTEXC, 0.0_JPRB)
         ZQEXC            = MAX(ZQEXC, 0.0_JPRB)

         ! cgl thv excess surface ; used to correct buoyancy flux
         ZTVEXCSURF(JL,JD) = ZTEXC + RETV*(ZQEXC*PTM1(JL,KLEV-1) + ZQTENH(JL,KLEV-1)*ZTEXC)

         PQTUH(JL,KLEV-1,JD) = ZQTENH(JL,KLEV-1) + ZQEXC
         ZQCUH(JL,KLEV-1,JD) = ZQLENH(JL,KLEV-1) + ZQIENH(JL,KLEV-1)
         ZQUH (JL,KLEV-1,JD) = PQTUH(JL,KLEV-1,JD)  - ZQCUH(JL,KLEV-1,JD)
         PSLGUH(JL,KLEV-1,JD)= ZSLGENH(JL,KLEV-1) + RCPD * ZTEXC
         ZTUH (JL,KLEV-1,JD) = ( PSLGUH (JL,KLEV-1,JD) - PGEOH(JL,KLEV-1) &
              & + RLVTT*ZQLENH(JL,KLEV-1) + RLSTT*ZQIENH(JL,KLEV-1) &
              & ) / RCPD
         PUUH(JL,KLEV-1,JD)= ZUENH(JL,KLEV-1)
         PVUH(JL,KLEV-1,JD)= ZVENH(JL,KLEV-1)

      ENDIF !KPBLTYPE & ZFRACB

   ENDDO !JL


   !-- Release the updraft --

   CALL VDFPARCELHL(YDEPHLI,YDPARAR,KIDIA,KFDIA,KLON,KLEV,KDRAFT,PGEOH,PGEOM1,PAPHM1,PUM1,PVM1,ZQTM1,ZSLGM1,ZTVEN,PUUH, &
     & PVUH,PSLGUH,PQTUH,ZWU2H,ZQCUH,ZBUOF,ZQUH,ZTUH,ZEPS,PZPLCL,KPLCL,PZPTOP,KPTOP,KPLZB,JD,ZUPGENL,ZUPGENN,ZTAUEPS,ZW2THRESH, &
     & LLDONE,KPBLTYPE)


ENDDO !JD


!*        6.3. In case no lcl is found in final updraft calculation, do some resque
!*
!*
! CGL      made an adjustment in vdfparcel to initialize pzplcl = -100

DO JL=KIDIA,KFDIA

   IF ( KPBLTYPE(JL)==2 .OR. KPBLTYPE(JL)==3) THEN

      IF (  PZPLCL(JL,3) < 0._JPRB .OR. KPLCL(JL,3)<KPTOP(JL,3) &
      &.OR. KPLCL(JL,3)==0  .OR. KPTOP(JL,3)==0 ) THEN
         KPBLTYPE(JL) = 1
         ZFRACB(JL,2) = ZFRACMAX
         ZFRACB(JL,3) = 0._JPRB
         KPLCL(JL,3) = 0
      ENDIF

   ENDIF

ENDDO !JL

!*        6.5  Updraft precipitation fluxes (rain and snow)
!*
DO JD = 3,KDRAFT  !moist updrafts only

   DO JK=2,KLEV
      DO JL=KIDIA,KFDIA

         ZDZRHO = ZRG * ( PAPHM1(JL,JK)-PAPHM1(JL,JK-1) )

         !-- Add precip generation to flux [kg /m2 /s: tendency * layer depth * air density] --
         ZUPFLXL(JL,JK,JD) = ZUPFLXL(JL,JK-1,JD) + ZUPGENL(JL,JK,JD) * ZDZRHO
         ZUPFLXN(JL,JK,JD) = ZUPFLXN(JL,JK-1,JD) + ZUPGENN(JL,JK,JD) * ZDZRHO

         !-- Do some melting at freezing level (snow->rain) --
         IF (ZUPFLXN(JL,JK,JD)>0._JPRB .AND. PTM1(JL,JK) > RTT) THEN
!wc
!   No melting and evaporation in the convection scheme in case LTOTPREC=TRUE
!   because this will be done inside microphysics
          IF (LTOTPREC) THEN
            ZUPMELT = 0._JPRB
          ELSE
            ZUPMELT = (1.0_JPRB+0.5_JPRB*(PTM1(JL,JK)-RTT)) * &
                 & (PTM1(JL,JK)-RTT) * RCPD/(RLMLT*RTAUMEL) * ZDZRHO
            ZUPMELT = MIN(ZUPFLXN(JL,JK,JD),ZUPMELT)
          ENDIF
            ZUPFLXL(JL,JK,JD) = ZUPFLXL(JL,JK,JD) + ZUPMELT
            ZUPFLXN(JL,JK,JD) = ZUPFLXN(JL,JK,JD) - ZUPMELT
         ENDIF

         ZPFLXTOT = ZUPFLXL(JL,JK,JD) + ZUPFLXN(JL,JK,JD)

         IF (ZPFLXTOT>0._JPRB) THEN

            !-- Saturation deficit of mean state T -
            ZESW=ESATW(PTM1(JL,JK))
            ZESI=ESATI(PTM1(JL,JK))
            ZFAC = ZUPFLXL(JL,JK,JD) / ZPFLXTOT
            ZES = ZESW*ZFAC + ZESI*(1._JPRB-ZFAC) ! Weigting according to precititation type
            ZQSAT = 0.62198_JPRB*ZES/(MAX(ZES,PAPM1(JL,JK))-0.37802_JPRB*ZES)
          ! Above boiling point for PAPHM1(JL,JK) < ZES --> no condensation. May happen in stratosphere, ZQSATU becomes 1.
            ZSATDEF=ZQSAT-PQM1(JL,JK) ! Also allow deposition = negative values

            !-- Precip evaporation tendency [kg/kg /s] (Kessler 1969, Tiedtke 1993) --
!wc
!   No melting and evaporation in the convection scheme in case LTOTPREC=TRUE
!   because this will be done inside microphysics
            IF (LTOTPREC) THEN
               ZPEVAPUP = 0._JPRB
            ELSE
               ZPEVAPUP = 0.001_JPRB * ZSATDEF * ( &             !cy32r1
                    & ( ZPFLXTOT / 0.00509_JPRB ) * &
                    & ( PAPM1(JL,JK)/PAPHM1(JL,KLEV) )**0.5_JPRB &
                    & )**0.5777_JPRB

            ENDIF
            !-- Back-partition evaporation and substract from fluxes --
            ZUPFLXL(JL,JK,JD) = ZUPFLXL(JL,JK,JD) - ZPEVAPUP * ZDZRHO * ZFAC
            ZUPFLXN(JL,JK,JD) = ZUPFLXN(JL,JK,JD) - ZPEVAPUP * ZDZRHO * (1._JPRB - ZFAC)
            ZUPFLXL(JL,JK,JD) = MAX(0._JPRB,ZUPFLXL(JL,JK,JD))
            ZUPFLXN(JL,JK,JD) = MAX(0._JPRB,ZUPFLXN(JL,JK,JD))
         ENDIF

      ENDDO
   ENDDO

   !Add contribution to total flux - weight by updraft area fraction
   !(or weighted by the mean cloud fraction in the cloud layer
   ! in case of LTOTPREC = TRUE.)

   IF (LTOTPREC .AND. (.NOT.LTOTPRECL)) THEN

      !compute mean cloud-fraction in the convective cloud layer:

      DO JL=KIDIA,KFDIA
         PCLFR(JL)=0._JPRB
      ENDDO

      DO JL=KIDIA,KFDIA
         ZB1(JL)=0.0_JPRB
         ZCOUNT(JL)=0.0_JPRB
         ITOP=KPTOP(JL,3)
         IBASE=KPLCL(JL,3)
         IF (ITOP > 0 .AND. IBASE > 0) THEN
            DO JK=ITOP,IBASE
               ZB1(JL)=ZB1(JL)+PAM1(JL,JK)
               ZCOUNT(JL)=ZCOUNT(JL)+1._JPRB
            ENDDO
         ENDIF
      ENDDO

      DO JL=KIDIA,KFDIA
         IF (ZB1(JL) > 0.0_JPRB .AND. ZCOUNT(JL) > 0.0_JPRB) THEN
            PCLFR(JL)=ZB1(JL)/ZCOUNT(JL)
         ENDIF
      ENDDO

      DO JK=0,KLEV
         DO JL=KIDIA,KFDIA
            PFPLVL(JL,JK) = PFPLVL(JL,JK) + PCLFR(JL) * ZUPFLXL(JL,JK,JD)
            PFPLVN(JL,JK) = PFPLVN(JL,JK) + PCLFR(JL) * ZUPFLXN(JL,JK,JD)
         ENDDO
      ENDDO

   ELSE
      DO JK=0,KLEV
         DO JL=KIDIA,KFDIA
            PFPLVL(JL,JK) = PFPLVL(JL,JK) + ZFRACB(JL,JD) * ZUPFLXL(JL,JK,JD)
            PFPLVN(JL,JK) = PFPLVN(JL,JK) + ZFRACB(JL,JD) * ZUPFLXN(JL,JK,JD)
         ENDDO
      ENDDO
   ENDIF

ENDDO !JD



!     -----------------------------------------------------------------

!*         7.     CONSTRUCT MASS FLUX PROFILES (JD=2,3)
!                 -------------------------------------


!*         7.1  Determine the mixed layer scaling height for JD=2,3
!*


!*         7.2  Construct subcloud / mixed layer mass fluxes
!*               - use constant area fraction, and multiply by parcel w
!*
DO JD = 2,KDRAFT

   DO JK=KLEV-1,1,-1

      DO JL=KIDIA,KFDIA

         IF ( KPBLTYPE(JL)/=0 .AND. ZFRACB(JL,JD)>0._JPRB) THEN

            ZWUH = ZWU2H(JL,JK,JD)**0.5_JPRB

            PMFLX(JL,JK,JD)  = ZFRACB(JL,JD) * ZWUH * ZRHOH(JL,JK)

            IF (JD == 3 .AND. JK >= KPLCL(JL,3) ) THEN
               ! do only below LCL !!

               PMFLX(JL,JK,JD)  = 0.35_JPRB * 0.1_JPRB * ZWSTAR(JL) *  &
                    &         ZRHOH(JL,JK)*(PGEOH(JL,JK)-PGEOH(JL,KLEV))/ &
                    &        (PGEOH(JL,KPLCL(JL,3))-PGEOH(JL,KLEV))

               !CGL cloud depth correction, if cloud thin then limit
               ZDZCLOUD(JL) = MAX(PGEOH(JL,KPTOP(JL,3))*ZRG-PGEOH(JL,KPLCL(JL,3))*ZRG,0._JPRB)
               IF( ZDZCLOUD(JL)<400._JPRB ) THEN
                  PMFLX(JL,JK,JD) = PMFLX(JL,JK,JD) * ZDZCLOUD(JL) / 400._JPRB
               ENDIF

            ENDIF

            IF (ZWU2H(JL,JK,JD)>0._JPRB) THEN
               ZFRAC(JL,JK,JD) = ZFRACB(JL,JD)
            ELSE
               ZFRAC(JL,JK,JD) = 0._JPRB
            ENDIF

         ENDIF

      ENDDO !JL

   ENDDO !JK

ENDDO !JD


! Do computation for Chi critical dpendency
!First, determine chicritmean in lower half of the
!                 cloud layer. chicritmean is used to determine
!                 the fractional detrainment coefficient according
!                 to De Rooy & Siebesma MWR 2008. In the upper
!                 half of the cloud layer a linear decreasing
!                 mass flux to 0 at cloud layer top is prescribed.

DO JL=KIDIA,KFDIA

   ! In between; determine thetavup of dry updraft (2)
   ! This can be used for the thetav flux contribution to TKE

   DO JK=KLEV-1,1,-1
      IF (KPTOP(JL,2) <= JK .AND. KPBLTYPE(JL) > 0) THEN
         PTHTVUH(JL,JK,2) = (ZTUH(JL,JK,2) / PEXNH(JL,JK))&
              &  * (1._JPRB + RETV * ZQUH(JL,JK,2) )
      ENDIF
   ENDDO

   IF (KPBLTYPE(JL) == 3   .OR. KPBLTYPE(JL) == 2 ) THEN

      ZSTAR(JL) = ((PGEOH(JL,KPTOP(JL,3))-PGEOH(JL,KLEV))*ZRG+&
      &            (PGEOH(JL,KPLCL(JL,3))-PGEOH(JL,KLEV))*ZRG)/2._JPRB
      ZMINIMUM(JL) = 10000._JPRB
      ZMEANCHICRIT(JL) = 0._JPRB

      ! Find level half way the cloud layer

      !CGL added security
      IKSTAR(JL) = KPTOP(JL,3)

      DO JK=KPLCL(JL,3),KPTOP(JL,3),-1


         IF (ABS((PGEOH(JL,JK)-PGEOH(JL,KLEV))*ZRG - ZSTAR(JL)) < ZMINIMUM(JL)) THEN
            ZMINIMUM(JL) = ABS((PGEOH(JL,JK)-PGEOH(JL,KLEV))*ZRG - ZSTAR(JL))
            IKSTAR(JL) = JK
         ENDIF


         ! determine chicrit in cloud layer

         ZTHTVENH(JL,JK) = (ZTENH(JL,JK) / PEXNH(JL,JK))&
              &   * (1._JPRB + RETV * ZQVENH(JL,JK)&
              ! add qice correction
              &    - ZQLENH(JL,JK) -ZQIENH(JL,JK) )
         ZTHTLUH(JL,JK) = (ZTUH(JL,JK,3) / PEXNH(JL,JK)) -&
              & (ZLAT2CP/PEXNH(JL,JK)*ZQCUH(JL,JK,3))
         !  N.B. ZQCUH is initialized with ql and qi from the environment!?!

         ZTHTVUH(JL,JK) = (ZTUH(JL,JK,3)/ PEXNH(JL,JK))&
              &  * (1._JPRB + RETV * ZQUH(JL,JK,3)&
              &    - ZQCUH(JL,JK,3))
         PTHTVUH(JL,JK,3)=ZTHTVUH(JL,JK)
         ZTHTLEH(JL,JK) = (ZTENH(JL,JK) / PEXNH(JL,JK)) -&
              & ZLAT2CP/ PEXNH(JL,JK)*&
              & ZQLENH(JL,JK)
         ZESW=ESATW(ZTUH(JL,JK,3))
         ZESI=ESATI(ZTUH(JL,JK,3))
         ZFAC = ZQIENH(JL,JK) / (ZQLENH(JL,JK) +ZQIENH(JL,JK) + 1.0E-20_JPRB) ! Weighting according to cloud condensate
         ZES = ZESI*ZFAC + ZESW*(1._JPRB-ZFAC)
         ZQSATU(JL,JK) = 0.62198_JPRB*ZES/(MAX(ZES,PAPHM1(JL,JK))-0.37802_JPRB*ZES)
        ! Above boiling point for PAPHM1(JL,JK) < ZES --> no condensation. May happen in stratosphere, ZQSATU becomes 1.
         ZDQSDTU(JL,JK) = ZEL2R * ZQSATU(JL,JK)/(ZTUH(JL,JK,3) * ZTUH(JL,JK,3))
         ZGAMMA(JL,JK) = ZLAT2CP*ZDQSDTU(JL,JK)
         ZDUMFUNC(JL,JK) = (1._JPRB/(1._JPRB+ZGAMMA(JL,JK)))*&
              & ((ZQTENH(JL,JK)-PQTUH(JL,JK,3)) -&
              &  PEXNH(JL,JK) * ZDQSDTU(JL,JK) *&
              &  (ZTHTLEH(JL,JK) - ZTHTLUH(JL,JK)))
         ZCHICRIT(JL,JK) = (ZTHTVENH(JL,JK)-ZTHTVUH(JL,JK))/&
              & (ZTHTLUH(JL,JK)*(RETV*(ZQTENH(JL,JK) -&
              & PQTUH(JL,JK,3)) -  (1._JPRB +RETV)*&
              & ZDUMFUNC(JL,JK)) + (ZTHTLEH(JL,JK) -&
              & ZTHTLUH(JL,JK)) * (1._JPRB + RETV*&
              & PQTUH(JL,JK,3) - (1._JPRB + RETV)*&
              & ZQCUH(JL,JK,3)) + ZLAT2CP/PEXNH(JL,JK)*&
              &  ZDUMFUNC(JL,JK))

         !CGL ; do some safety

         ZCHICRIT(JL,JK) =  MIN(ZCHICRIT(JL,JK),1._JPRB)

      ENDDO !JK

      !   Calculate weighted mean chicrit over lower half cloud layer
      !   Exclude cloud base level (see de Rooy & Siebesma MWR 2008)

      IF (KPLCL(JL,3)-IKSTAR(JL) > 1) THEN
         DO JK= KPLCL(JL,3)-1,IKSTAR(JL),-1
            ZMEANCHICRIT(JL)=ZMEANCHICRIT(JL)+ZCHICRIT(JL,JK)*&
                 & (PGEOM1(JL,JK)*ZRG-PGEOM1(JL,JK+1)*ZRG)
         ENDDO
         ZMEANCHICRIT(JL)=ZMEANCHICRIT(JL)/(PGEOM1(JL,IKSTAR(JL))*ZRG-&
              &                             PGEOM1(JL,KPLCL(JL,3))*ZRG)

         !   Fraction of mass flux that is left half way the cloud layer (ZFRACMB)
         !   is determined using LES based relation between mean chicrit and ZFRACMB
         !   Small adjustment to somewhat more active Mass flux profiles
         !   Note that the parameterization counteracts changes in the relation, therefore
         !   the ultimate differences are small

         ZFRACMB(JL)=MAX(0.05_JPRB,5.24_JPRB*ZMEANCHICRIT(JL)-0.39_JPRB)
         ZFRACMB(JL)=MIN(1.0_JPRB,ZFRACMB(JL))     !safety CGL


         ! Determine (constant) detrainment coefficient
         ! with new eps, detr should be adapted
         ! However we can also use original detr formulation as long as we also
         ! use eps=1/z in the mass flux calculation (the new eps is then only used
         ! for the updraft dilution)
         ZDETRSHALLOW(JL)=LOG(((PGEOH(JL,IKSTAR(JL))-PGEOH(JL,KLEV))*ZRG)/&
              &  ((PGEOH(JL,KPLCL(JL,3))-PGEOH(JL,KLEV))*ZRG&
              &   *ZFRACMB(JL)))/((PGEOH(JL,IKSTAR(JL))-PGEOH(JL,KLEV))*ZRG-&
              &  (PGEOH(JL,KPLCL(JL,3))-PGEOH(JL,KLEV))*ZRG)
      ELSE
         !  only for very shallow (less that 2 layers) cloud layers
         ZDETRSHALLOW(JL)=0.00275_JPRB
         ! the lower assignement only has a meaning if used in combination with
         ! extremely thin clouds
         ZFRACMB(JL)=0.3_JPRB
      ENDIF

   ENDIF ! KPBLTYPE=3

ENDDO !JL


! end computation by for chicritical


!*         7.3  Construct cloudy mass flux profile (JD=3 only)
!*


DO JK=KLEV-2,1,-1

   DO JL=KIDIA,KFDIA

      IF ( KPBLTYPE(JL)/=0 .AND. KPBLTYPE(JL)/=4 .AND. ZFRACB(JL,3)> 0._JPRB ) THEN

         IF (JK>=KPTOP(JL,3) .AND. JK<KPLCL(JL,3)) THEN

            ZWUH = ZWU2H(JL,JK,3)**0.5_JPRB

            !  Cloud layer vertical structure

            ! ZDELTAMINEPS should be named ZEPSMINDELTA

            IF (JK+1 >= IKSTAR(JL)) THEN !lower half cloud layer

               ! With changed eps formulation change zdeltamineps also

               ZDELTAMINEPS(JL,JK+1) = RG/(PGEOH(JL,JK+1)-PGEOH(JL,KLEV))-ZDETRSHALLOW(JL)

            ELSE !upper half cloud layer
               ZDELTAMINEPS(JL,JK+1) = -1._JPRB/(PGEOH(JL,KPTOP(JL,3))*ZRG -&
                    &  PGEOH(JL,JK+1)*ZRG)
            ENDIF

            ZDZ          = (PGEOH(JL,JK) - PGEOH(JL,JK+1))*ZRG

            PMFLX(JL,JK,3) = PMFLX(JL,JK+1,3) * EXP( ZDZ * ( ZDELTAMINEPS(JL,JK+1) ) ) !exact

            ! make a guess at area fraction;

            ZWUH = ZWU2H(JL,JK,3)**0.5_JPRB
            ZFRAC(JL,JK,3) = PMFLX(JL,JK,3) / ( ZWUH * ZRHOH(JL,JK) )

         ENDIF

      ENDIF

   ENDDO !JL

ENDDO !JK



!*        7.4  Compute some fluxes used in the qt-variance budget in VDFMAIN.
!*

!*        7.5  Mass flux limit according to CFL criterion
!*              (reduce M profiles uniformly by maximum excess)
!*
ZCONS10 = 1.0_JPRB/(RG*PTMST)
DO JD = 2,KDRAFT
   DO JL=KIDIA,KFDIA
      ZMFS(JL,JD) = 1.0_JPRB  ! mass flux scaling value (reduction)
   ENDDO
ENDDO

DO JD = 2,KDRAFT
   DO JK=1,KLEV-1
      DO JL=KIDIA,KFDIA
         IF ( JK >= KPTOP(JL,JD) .AND. KPTOP(JL,JD)>0) THEN
            ZMFMAX = (PAPM1(JL,JK+1)-PAPM1(JL,JK)) * ZCONS10
            IF ( PMFLX(JL,JK,JD) > ZMFMAX ) THEN
               ZMFS(JL,JD) = MIN(ZMFS(JL,JD),ZMFMAX/PMFLX(JL,JK,JD))
            ENDIF
         ENDIF
      ENDDO
   ENDDO
ENDDO

DO JD = 2,KDRAFT
   DO JK=1,KLEV
      DO JL=KIDIA,KFDIA
         PMFLX(JL,JK,JD) = PMFLX(JL,JK,JD)*ZMFS(JL,JD)
      ENDDO
   ENDDO
ENDDO




!     -----------------------------------------------------------------

!*         8.     W-SCALE USED IN CLOUD VARIANCE DISSIPATION (VDFMAIN)
!                 ----------------------------------------------------


!CGL fill this array at the end
PFRACB = ZFRACB


! ----------------------------------------------------
!
!       9. OUTPUT TO BE USED IN VDFEXCU now called from this routine
!          to calculate length scale according to Geert Lenderink
!   if LHARATU=true
!   Call vdfexcuhl to calculate length scales and TKE for use in Harmonie

IF (LHARATU) THEN

   PBUOY_COR (:,:) = 0.0_JPRB
   PWU       (:,:) = 0.0_JPRB
   DO JL=KIDIA,KFDIA

      DO JK=1,KLEV-1
         PWU (JL,JK) = (MAX( ZWU2H(JL,JK,2),0._JPRB ))**0.5_JPRB
         IF (ZWSTAR(JL)*PZPTOP(JL,2).GT.0._JPRB) THEN
            ! PBUOY_COR (JL,JK) =  ZFRAC(JL,JK,2)* ZWU2H(JL,JK,2)**(1._JPRB/2._JPRB) / ZWSTAR(JL) / PZPTOP(JL,2)
            TVEXCSURF = MAX(ZTVEXCSURF(JL,2),0.0_JPRB)
            PBUOY_COR (JL,JK) = TVEXCSURF*PWU(JL,JK)&
            &/ ZWSTAR(JL) / PZPTOP(JL,2) * RG / PTM1(JL,KLEV)
         ENDIF
      ENDDO
   ENDDO

   DO JK=1,KLEV
      JKM=MIN(JK,KLEV-1)
      DO JL=KIDIA,KFDIA
         !          dqsat/dT correction factor (1+L/cp*dqsat/dT) & alfa
         ZESW = ESATW(PTM1(JL,JK))
         ZESI = ESATI(PTM1(JL,JK))
         ZFAC = ZQIENH(JL,JKM) / (ZQLENH(JL,JKM) +ZQIENH(JL,JKM) + 1.0E-20_JPRB) ! Weighting according to cloudcondensate
         ZES = ZESI*ZFAC + ZESW*(1._JPRB-ZFAC)
         ZE2R = ZEI2R*ZFAC + ZEL2R*(1._JPRB-ZFAC)
         ZQSVAR(JL,JK) = 0.62198_JPRB*ZES/(MAX(ZES,PAPM1(JL,JK))-0.37802_JPRB*ZES)
         ! Above boiling point for PAPHM1(JL,JK) < ZES --> no condensation. May happen in stratosphere, ZQSATU becomes 1.
         ZDQSDTEMP(JL,JK) = ZE2R*ZQSVAR(JL,JK)/PTM1(JL,JK)/PTM1(JL,JK) !dqsat/dT
         !      write (913,'(i5,10f14.5)') jk, ZFAC,  ZCOR, ZDQSDTEMP(JL,JK),  ZDQSDTEMP(JL,JK)/ZQSVAR(JL,JK)
      ENDDO
   ENDDO
   DO JL=KIDIA,KFDIA
      DO JK=1,KLEV-1

!wc
!   energy cascade:   Two contributions: 1) subplume turbulence: related to eps
!                                        2) environmental turbulence: related to decrease M (detr)
!

         IF (KPBLTYPE(JL) == 3  .AND. (KPLCL(JL,3)-KPTOP(JL,3)) > 2) THEN
            ZFACCASC(JL,JK) = ZEL*(((PGEOH(JL,JK)-PGEOH(JL,KLEV))*ZRG)/((PGEOH(JL,KPLCL(JL,3))-PGEOH(JL,KLEV))*ZRG))* &
             & (1._JPRB/(1._JPRB+(((PGEOH(JL,JK)-PGEOH(JL,KPLCL(JL,3)))*ZRG)/ZWL)**2._JPRB))+ &
             & ZET*(1._JPRB/(1._JPRB+(((PGEOH(JL,JK)-PGEOH(JL,KPTOP(JL,3)))*ZRG)/ZWT)**2._JPRB))
            ZFACCASC(JL,JK) = MAX(0._JPRB,(MIN(ZFACCASC(JL,JK),0.01_JPRB)))
         ELSE
            ZFACCASC(JL,JK) = 0.002_JPRB
         ENDIF

         ZENCASC(JL,JK) = 0.5_JPRB*ZEPS(JL,JK,2)*ZWU2H(JL,JK,2)*PMFLX(JL,JK,2) + &
          & ZFACCASC(JL,JK)*ZWU2H(JL,JK,3)*PMFLX(JL,JK,3)

      ENDDO
   ENDDO
   !
   !   if LHARATU=true
   !   Call vdfexcuhl to calculate length scales and TKE for use in Harmonie
   !
   !
   !*           EXCHANGE COEFFICIENTS ABOVE THE SURFACE LAYER

   CALL VDFEXCUHL(YDVDF,YDEPHY,KIDIA,KFDIA,KLON,KLEV,PTMST,PUM1,PVM1,PTM1,PQM1,PLM1,PIM1,ZSLGM1,ZQTM1,PKMFL,PKHFL,PKQFL,PAPHM1,&
   & PAPM1,PGEOM1,PGEOH,PEXNF,ZZI,KPBLTYPE,KDRAFT,ZQSVAR,ZDQSDTEMP,PBUOY_COR,ZENCASC,PWU,&
   & YSPP_RFAC_TWOC,YSPP_RZC_H,YSPP_RZL_INF, & 
   & PTKE,PMFLX,ZLENGTH_M,ZLENGTH_H)

ENDIF ! LHARATU

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('VDFHGHTNHL',1,ZHOOK_HANDLE)

END SUBROUTINE VDFHGHTNHL
