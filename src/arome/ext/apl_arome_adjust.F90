SUBROUTINE APL_AROME_ADJUST(YDCST, YDMODEL, YDCPG_BNDS, YDCPG_OPTS, YDCPG_MISC, YDGEOMETRY, YDVARS, &
                          & YDGEOMVARS, YDMF_PHYS, YDMF_PHYS_BASE_STATE, &
                          & YSPP_ALL, YDDDH, &
                          & PRHODJM, PEXNREFM, &
                          & PPABSM, PWM, PTKEM, &
                          & PRHODREFM, &
                          & PTHM, PRM, PLIMAM, PTENDRA, &
                          & PTHS, PRS, PLIMAS, &
                          & PTENDLIMA, PZZ_F, PTENDT, &
                          & PSRCS, PNEBMNH, &
                          & PICLDFR, PWCLDFR, PSSIO, PSSIU, PIFR, &
                          & PQDM, PQVM, PQCM, PQRM, PQIM, PQSM, PQGM, PQHM, &
                          & PCPM, PRHM, PTM, &
                          & PAPHIM, PAPHIFM, &
                          & PZZ, PTENDTT, PDZZ, &
                          & PDTHRAD, PICEFR, PPRCFR, &
                          & PHLC_HRC, PHLC_HCF, PHLI_HRI, PHLI_HCF)

!**** *APL_AROME_ADJUST* - Setup and call Arome microphysics adjustment

!     Author.
!     -------
!      Rolf Heilemann Myhre *Met Norway*
!      Original : 01-09-2023
!      See APL_AROME for modification history

!-----------------------------------------------------------------------

! -   INPUT ARGUMENTS.
!     -------------------
!
!   YDCST:                Model constants
!   YDMODEL:              Model configuration settings
!   YDCPG_BNDS:           Array dimensions and bounds
!   YDCPG_OPTS:           Various dimensions
!   YDCPG_MISC:
!   YDGEOMETRY:           Geomtery settings from IFS model
!   YDVARS:               Persistent fields
!   YDMF_PHYS:            Physics fields
!   YDMF_PHYS_BASE_STATE:
!   YSPP_ALL:             SPP fields and settings
!   YDDDH:                DDH fields and settings
!
!   PTM:    Temperature.
!   PCPM:   Specific heat at constant pressure for air
!   PRM:    Moist variables at time t
!   PTKEM:  Turbulent kinetic energy
!   PTENDT: Temperature tendency
!

  USE PARKIND1, ONLY: JPRB, JPRD, JPIM

  USE YOMCST,                      ONLY: TCST
  USE TYPE_MODEL,                  ONLY: MODEL
  USE CPG_OPTS_TYPE_MOD,           ONLY: CPG_BNDS_TYPE, CPG_OPTS_TYPE
  USE CPG_TYPE_MOD,                ONLY: CPG_MISC_TYPE
  USE GEOMETRY_MOD,                ONLY: GEOMETRY
  USE FIELD_VARIABLES_MOD,         ONLY: FIELD_VARIABLES
  USE YOMGEOMVARS,                 ONLY: TGEOMVARS         

  USE MF_PHYS_TYPE_MOD,            ONLY: MF_PHYS_TYPE
  USE MF_PHYS_BASE_STATE_TYPE_MOD, ONLY: MF_PHYS_BASE_STATE_TYPE

  USE DDH_MIX,                     ONLY: TYP_DDH, NEW_ADD_FIELD_3D
  USE SPP_MOD_TYPE,                ONLY: UA_PHYS_SPP_VARS

  USE YOMHOOK,                     ONLY: LHOOK, DR_HOOK, JPHOOK

  IMPLICIT NONE

  TYPE(TCST)                    ,INTENT(IN)             :: YDCST
  TYPE(MODEL)                   ,INTENT(IN)             :: YDMODEL
  TYPE(CPG_BNDS_TYPE)           ,INTENT(IN)             :: YDCPG_BNDS
  TYPE(CPG_OPTS_TYPE)           ,INTENT(IN)             :: YDCPG_OPTS
  TYPE(CPG_MISC_TYPE)           ,INTENT(IN)             :: YDCPG_MISC
  TYPE(GEOMETRY)                ,INTENT(IN)             :: YDGEOMETRY
  TYPE(FIELD_VARIABLES)         ,INTENT(INOUT)          :: YDVARS
  TYPE(TGEOMVARS)               ,INTENT(IN)             :: YDGEOMVARS

  TYPE(MF_PHYS_TYPE)            ,INTENT(IN)             :: YDMF_PHYS
  TYPE(MF_PHYS_BASE_STATE_TYPE) ,INTENT(IN)             :: YDMF_PHYS_BASE_STATE

  TYPE(UA_PHYS_SPP_VARS)        ,INTENT(INOUT)          :: YSPP_ALL
  TYPE(TYP_DDH)                 ,INTENT(INOUT)          :: YDDDH

  REAL(KIND=JPRB)               ,INTENT(IN)             :: PRHODJM(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG+1)
  REAL(KIND=JPRB)               ,INTENT(IN)             :: PEXNREFM(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)

  REAL(KIND=JPRB)               ,INTENT(IN)             :: PPABSM(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG+1)

  REAL(KIND=JPRB)               ,INTENT(IN)             :: PWM(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG+1)
  REAL(KIND=JPRB)               ,INTENT(IN)             :: PTKEM(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG+1)

  REAL(KIND=JPRB)               ,INTENT(INOUT)          :: PRHODREFM(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG+1)

  REAL(KIND=JPRB)               ,INTENT(INOUT)          :: PTHM(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG+1)
  REAL(KIND=JPRB)               ,INTENT(INOUT)          :: PRM(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG,YDMODEL%YRML_PHY_MF%YRPARAR%NRR)
  REAL(KIND=JPRB)               ,INTENT(INOUT)          :: PLIMAM(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG,YDMODEL%YRML_GCONF%YGFL%NLIMA)
  REAL(KIND=JPRB)               ,INTENT(INOUT)          :: PTENDRA(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG,YDMODEL%YRML_PHY_MF%YRPARAR%NRR)

  REAL(KIND=JPRB)               ,INTENT(INOUT)          :: PTHS(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
  REAL(KIND=JPRB)               ,INTENT(INOUT)          :: PRS(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG,YDMODEL%YRML_PHY_MF%YRPARAR%NRR)
  REAL(KIND=JPRB)               ,INTENT(INOUT)          :: PLIMAS(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG,YDMODEL%YRML_GCONF%YGFL%NLIMA)
  REAL(KIND=JPRB)               ,INTENT(INOUT)          :: PTENDLIMA(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG,YDMODEL%YRML_GCONF%YGFL%NLIMA)
  REAL(KIND=JPRB)               ,INTENT(INOUT)          :: PZZ_F(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)

  REAL(KIND=JPRB)               ,INTENT(INOUT)          :: PTENDT(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG) ! temperature tendency

  REAL(KIND=JPRB)               ,INTENT(OUT)            :: PSRCS(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG+1)
  REAL(KIND=JPRB)               ,INTENT(OUT)            :: PNEBMNH(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)

  REAL(KIND=JPRB)               ,INTENT(OUT)            :: PICLDFR(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
  REAL(KIND=JPRB)               ,INTENT(OUT)            :: PWCLDFR(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
  REAL(KIND=JPRB)               ,INTENT(OUT)            :: PSSIO(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
  REAL(KIND=JPRB)               ,INTENT(OUT)            :: PSSIU(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
  REAL(KIND=JPRB)               ,INTENT(OUT)            :: PIFR(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)

  REAL(KIND=JPRB)               ,INTENT(OUT)            :: PQDM(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
  REAL(KIND=JPRB)               ,INTENT(OUT)            :: PQVM(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
  REAL(KIND=JPRB)               ,INTENT(OUT)            :: PQCM(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
  REAL(KIND=JPRB)               ,INTENT(OUT)            :: PQRM(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
  REAL(KIND=JPRB)               ,INTENT(OUT)            :: PQIM(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
  REAL(KIND=JPRB)               ,INTENT(OUT)            :: PQSM(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
  REAL(KIND=JPRB)               ,INTENT(OUT)            :: PQGM(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
  REAL(KIND=JPRB)               ,INTENT(OUT)            :: PQHM(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)

  REAL(KIND=JPRB)               ,INTENT(OUT)            :: PCPM(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
  REAL(KIND=JPRB)               ,INTENT(OUT)            :: PRHM(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
  REAL(KIND=JPRB)               ,INTENT(OUT)            :: PTM(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)

  REAL(KIND=JPRB)               ,INTENT(OUT)            :: PAPHIM(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)
  REAL(KIND=JPRB)               ,INTENT(OUT)            :: PAPHIFM(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)

  REAL(KIND=JPRB)               ,INTENT(OUT)            :: PZZ(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
  REAL(KIND=JPRB)               ,INTENT(OUT)            :: PTENDTT(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG) ! temperature tendency
  REAL(KIND=JPRB)               ,INTENT(OUT)            :: PDZZ(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)

  REAL(KIND=JPRB)               ,INTENT(OUT)            :: PDTHRAD(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
  REAL(KIND=JPRB)               ,INTENT(OUT)            :: PICEFR(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
  REAL(KIND=JPRB)               ,INTENT(OUT)            :: PPRCFR(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)

  REAL(KIND=JPRB)               ,INTENT(OUT)            :: PHLC_HRC(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
  REAL(KIND=JPRB)               ,INTENT(OUT)            :: PHLC_HCF(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
  REAL(KIND=JPRB)               ,INTENT(OUT)            :: PHLI_HRI(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
  REAL(KIND=JPRB)               ,INTENT(OUT)            :: PHLI_HCF(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)

  REAL(KIND=JPRB) :: ZTHSIN(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG)
  REAL(KIND=JPRB) :: ZRSIN(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG, YDMODEL%YRML_PHY_MF%YRPARAR%NRR)
  REAL(KIND=JPRB) :: ZLIMASIN(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG, YDMODEL%YRML_GCONF%YGFL%NLIMA)

  REAL(KIND=JPRB) :: ZSIGM(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG)
  REAL(KIND=JPRB) :: ZMFM(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG)

  REAL(KIND=JPRB) :: ZRC_MF(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG)
  REAL(KIND=JPRB) :: ZRI_MF(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG)
  REAL(KIND=JPRB) :: ZCF_MF(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG)
  REAL(KIND=JPRB) :: ZHLC_HRC_MF(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG)
  REAL(KIND=JPRB) :: ZHLC_HCF_MF(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG)
  REAL(KIND=JPRB) :: ZHLI_HRI_MF(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG)
  REAL(KIND=JPRB) :: ZHLI_HCF_MF(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG)
  REAL(KIND=JPRB) :: ZWEIGHT_MF_CLOUD(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG)

  REAL(KIND=JPRB) :: ZTMPAF(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG)

  REAL(KIND=JPRB) :: ZCON1(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG)
  REAL(KIND=JPRB) :: ZCON2(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG)
  REAL(KIND=JPRB) :: ZCON3(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG)

  REAL(KIND=JPRB) :: ZWNU(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)

  ! ACRANEB2 local variables
  REAL(KIND=JPRB) :: ZNEB0    (YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)  ! protected cloud fractions
  REAL(KIND=JPRB) :: ZDECRD   (YDCPG_OPTS%KLON)       ! decorrelation depth
  REAL(KIND=JPRB) :: ZCLCT_RAD(YDCPG_OPTS%KLON)       ! total cloud cover for radiation

  CHARACTER(LEN=11) :: CLNAME
  CHARACTER(LEN=2), DIMENSION(7), PARAMETER :: CLVARNAME=(/"QV","QL","QR","QI","QS","QG","QH"/)

  REAL(KIND=JPRB) :: ZQHGM(YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG)

  INTEGER(KIND=JPIM) :: JLON, JLEV, JRR, JGFL

  INTEGER(KIND=JPIM), PARAMETER :: IKL = -1
  INTEGER(KIND=JPIM), PARAMETER :: IKU = 1

  REAL(KIND=JPRB) :: ZDT, ZRHO, ZINVG, ZEPSNEB

  LOGICAL :: LLIMAINIT

  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

  !------------------------------------------------------------------

#include "aro_adjust.h"
#include "aro_adjust_lima.h"
#include "gprcp_qlirsg.intfb.h"
#include "gpgeo.intfb.h"
#include "aro_suintbudget_omp.h"
#include "acnpart.intfb.h"

    !------------------------------------------------------------------

  IF (LHOOK) CALL DR_HOOK('APL_AROME_ADJUST', 0, ZHOOK_HANDLE)

  ASSOCIATE(KLON => YDCPG_OPTS%KLON, KFLEVG => YDCPG_OPTS%KFLEVG, KIDIA => YDCPG_BNDS%KIDIA, KFDIA => YDCPG_BNDS%KFDIA)

  IF (YDMODEL%YRML_PHY_MF%YRARPHY%LMICRO) THEN

    !initialise sigma for subgrid condensation coming
    !from previous time step turbulence scheme
    IF (YDMODEL%YRML_PHY_MF%YRPARAR%PHYEX%NEBN%LSIGMAS) THEN
      ZSIGM(:,:) = YDMF_PHYS_BASE_STATE%SRC%P(:,:)
    ENDIF

    !initialise convective mas flux for subgrid condensation coming
    !from previous time step convection scheme
    IF (YDMODEL%YRML_PHY_MF%YRPARAR%PHYEX%NEBN%LSUBG_COND .AND. .NOT. &
      & YDMODEL%YRML_PHY_MF%YRPARAR%PHYEX%NEBN%LSIGMAS) THEN

      IF (YDMODEL%YRML_PHY_MF%YRARPHY%LKFBCONV) THEN
        ZMFM(:,:) = YDMF_PHYS_BASE_STATE%SRC%P(:,:)
      ELSE
        ZMFM(:,:) = 0._JPRB
      ENDIF

    ENDIF

    IF (YDMODEL%YRML_DYN%YRDYNA%LTWOTL) THEN
      ZDT = YDCPG_OPTS%ZDTPHY/2._JPRB
    ELSE
      IF (YDCPG_OPTS%NSTEP /= 0) THEN
        ZDT = YDCPG_OPTS%ZDTPHY/2._JPRB
      ELSE
        ZDT = YDCPG_OPTS%ZDTPHY
      ENDIF
    ENDIF

    ZINVG=1._JPRB/YDCST%RG

    IF (YDMODEL%YRML_PHY_MF%YRARPHY%LMFSHAL .AND. &
     & (YDMODEL%YRML_PHY_MF%YRPARAR%PHYEX%PARAM_MFSHALLN%CMF_CLOUD == 'DIRE' .OR. &
     &  YDMODEL%YRML_PHY_MF%YRPARAR%PHYEX%PARAM_MFSHALLN%CMF_CLOUD == 'BIGA')) THEN

      IF (YDCPG_OPTS%NSTEP == 0) THEN
        ZRC_MF(:,:) = 0._JPRB
        ZRI_MF(:,:) = 0._JPRB
        ZCF_MF(:,:) = 0._JPRB
        ZWEIGHT_MF_CLOUD(:,:) = 0._JPRB
        ZHLC_HRC_MF(:,:) = 0._JPRB
        ZHLC_HCF_MF(:,:) = 0._JPRB
        ZHLI_HRI_MF(:,:) = 0._JPRB
        ZHLI_HCF_MF(:,:) = 0._JPRB
      ELSE
        ZRC_MF(:,:) = YDMF_PHYS_BASE_STATE%LMF%P(:,:)
        ZCF_MF(:,:) = YDMF_PHYS_BASE_STATE%AMF%P(:,:)
        ZRI_MF(:,:) = YDMF_PHYS_BASE_STATE%IMF%P(:,:)
        ZWEIGHT_MF_CLOUD(:,:) = YDMF_PHYS_BASE_STATE%WMFC%P(:,:)
        ZHLC_HRC_MF(:,:) = YDMF_PHYS_BASE_STATE%HLMF%P(:,:)
        ZHLC_HCF_MF(:,:) =YDMF_PHYS_BASE_STATE%HLCFMF%P(:,:)
        ZHLI_HRI_MF(:,:) = YDMF_PHYS_BASE_STATE%HIMF%P(:,:)
        ZHLI_HCF_MF(:,:) = YDMF_PHYS_BASE_STATE%HICFMF%P(:,:)
      ENDIF
    ELSE
      ZRC_MF(:,:) = 0._JPRB
      ZRI_MF(:,:) = 0._JPRB
      ZCF_MF(:,:) = 0._JPRB
      ZWEIGHT_MF_CLOUD(:,:) = 0._JPRB
      ZHLC_HRC_MF(:,:) = 0._JPRB
      ZHLC_HCF_MF(:,:) = 0._JPRB
      ZHLI_HRI_MF(:,:) = 0._JPRB
      ZHLI_HCF_MF(:,:) = 0._JPRB
    ENDIF

    ! for now a copy is needed (see below, inside). I don't like than :-( REK
    ZTHSIN(:,:) = PTHS(:,:)
    ZRSIN(:,:,:) = PRS(:,:,:)

    IF (YDMODEL%YRML_PHY_MF%YRPARAR%CMICRO == 'LIMA') THEN

      !initialisation de ZZI_THRAD
      IF (YDCPG_OPTS%NSTEP==0) THEN
        LLIMAINIT = .TRUE.
        PDTHRAD(:,:) = 0._JPRB
        YDMF_PHYS_BASE_STATE%THRAD%P(:,:) = 0._JPRB
      ELSE
        PDTHRAD(:,:) = YDMF_PHYS_BASE_STATE%THRAD%P(:,:)
        LLIMAINIT = .FALSE.
      ENDIF

      ZWNU(:,:) = PWM(:,1:KFLEVG)
      IF (YDMODEL%YRML_PHY_MF%YRARPHY%LTURB) THEN
        DO JLEV=1,KFLEVG
          DO JLON=KIDIA, KFDIA
            ZWNU(JLON,JLEV) = ZWNU(JLON,JLEV) + 0.66*SQRT(PTKEM(JLON,JLEV))
          ENDDO
        ENDDO
      ENDIF

      ! for now a copy is needed (see below, inside). I don't like than :-( REK
      ZLIMASIN(:,:,:) = PLIMAS(:,:,:)

      CALL ARO_ADJUST_LIMA(YDMODEL%YRML_PHY_MF%YRPARAR%PHYEX,                          &
      & KFLEVG, IKU, IKL,                                                              &
      & KLON, KFLEVG, KIDIA, KFDIA,                                                    &
      & YDMODEL%YRML_PHY_MF%YRPARAR%NRR,                                               &
      & YDMODEL%YRML_GCONF%YGFL%NLIMA, YDCPG_OPTS%NSTEP+1,                             &
      & YDMODEL%YRML_PHY_MF%YRPARAR%PHYEX%NEBN%LSUBG_COND,                             &
      & YDMODEL%YRML_PHY_MF%YRPARAR%PHYEX%NEBN%LSIGMAS,                                &
      & ZDT, YDMODEL%YRML_PHY_MF%YRPARAR%PHYEX%NEBN%VSIGQSAT, PZZ_F,                   &
      & PRHODJM(:, 1:KFLEVG), PRHODREFM(:, 1:KFLEVG), PEXNREFM,                        &
      & PPABSM(:, 1:KFLEVG), PTHM(:, 1:KFLEVG),                                        &
      & PRM, PLIMAM, ZSIGM, ZWNU, PDTHRAD, ZMFM, ZRC_MF, ZRI_MF, ZCF_MF, ZWEIGHT_MF_CLOUD, &
      & PTHS, PRS, PLIMAS, PSRCS(:, 1:KFLEVG), PNEBMNH,                                &
      & PICEFR, PPRCFR,                                                                &
      & YDDDH, YDMODEL%YRML_DIAG%YRLDDH, YDMODEL%YRML_DIAG%YRMDDH,                     &
      & LLIMAINIT,                                                                     &
      & PHLC_HRC, PHLC_HCF, PHLI_HRI, PHLI_HCF,                                        &
      & ZHLC_HRC_MF, ZHLC_HCF_MF, ZHLI_HRI_MF, ZHLI_HCF_MF                             )

    ELSE

      CALL ARO_ADJUST(YDMODEL%YRML_PHY_MF%YRPARAR%PHYEX,           &
      & KLON, KIDIA, KFDIA, KFLEVG,                                &
      & YDMODEL%YRML_PHY_MF%YRPARAR%NRR,                           &
      & YDMODEL%YRML_PHY_MF%YRPARAR%CMICRO,                        &
      & ZDT, PZZ_F,                                                &
      & PRHODJM(:, 1:KFLEVG), PEXNREFM,                            &
      & PRHODREFM(:, 1:KFLEVG),                                    &
      & PPABSM(:, 1:KFLEVG), PTHM(:, 1:KFLEVG),                    &
      & PRM, ZSIGM, ZMFM,                                          &
      & ZRC_MF, ZRI_MF, ZCF_MF, ZWEIGHT_MF_CLOUD,                  &
      & PTHS, PRS, PSRCS(:, 1:KFLEVG),                             &
      & PNEBMNH,                                                   &
      & PICLDFR, PWCLDFR, PSSIO, PSSIU, PIFR,                      &
      & PHLC_HRC, PHLC_HCF, PHLI_HRI, PHLI_HCF,                    &
      & ZHLC_HRC_MF, ZHLC_HCF_MF, ZHLI_HRI_MF, ZHLI_HCF_MF,        &
      & YDDDH, YDMODEL%YRML_DIAG%YRLDDH, YDMODEL%YRML_DIAG%YRMDDH, &
      & YSPP_ALL%YSPP_PSIGQSAT, YSPP_ALL%YSPP_ICE_CLD_WGT)

    ENDIF

    YDVARS%A%T1(:,:) = PNEBMNH(:,:)

    !adjusted zthm and zrm
    DO JLEV = 1, KFLEVG
      DO JLON = KIDIA,KFDIA
        PTHM(JLON, JLEV) = PTHS(JLON,JLEV)*YDCPG_OPTS%ZDTPHY
      ENDDO
    ENDDO

    DO JRR=1, YDMODEL%YRML_PHY_MF%YRPARAR%NRR
      DO JLEV = 1, KFLEVG
        DO JLON = KIDIA,KFDIA
          PRM(JLON, JLEV, JRR) = PRS(JLON, JLEV, JRR)*YDCPG_OPTS%ZDTPHY
        ENDDO
      ENDDO
    ENDDO

    !initialisation de qdm utile pour
    !convertir tendance de r en tendance de q

    PQDM(:,:) = 1._JPRB
    DO JRR=1, YDMODEL%YRML_PHY_MF%YRPARAR%NRR
      DO JLEV = 1, KFLEVG
        DO JLON = KIDIA,KFDIA
          PQDM(JLON,JLEV) = PQDM(JLON,JLEV) + PRM(JLON,JLEV,JRR)
        ENDDO
      ENDDO
    ENDDO
    DO JLEV = 1, KFLEVG
      DO JLON = KIDIA,KFDIA
        PQDM(JLON,JLEV) = 1._JPRB/PQDM(JLON,JLEV)
      ENDDO
    ENDDO

    !reinitialisation des qi
    DO JLEV = 1, KFLEVG
      DO JLON = KIDIA, KFDIA
        PQVM(JLON,JLEV)=PRM(JLON,JLEV,1)*PQDM(JLON,JLEV)
        PQCM(JLON,JLEV)=PRM(JLON,JLEV,2)*PQDM(JLON,JLEV)
        PQRM(JLON,JLEV)=PRM(JLON,JLEV,3)*PQDM(JLON,JLEV)
        PQIM(JLON,JLEV)=PRM(JLON,JLEV,4)*PQDM(JLON,JLEV)
        PQSM(JLON,JLEV)=PRM(JLON,JLEV,5)*PQDM(JLON,JLEV)
        PQGM(JLON,JLEV)=PRM(JLON,JLEV,6)*PQDM(JLON,JLEV)
      ENDDO
    ENDDO

    IF (YDMODEL%YRML_PHY_MF%YRPARAR%YRTENDRA%LHAIL) THEN
      DO JLEV = 1, KFLEVG
        DO JLON = KIDIA, KFDIA
          PQHM(JLON,JLEV) = PRM(JLON,JLEV,7)*PQDM(JLON,JLEV)
        ENDDO
      ENDDO
    ELSE
      PQHM(:,:) = 0._JPRB
    ENDIF

    ! Tendances des variables LIMA
    DO JGFL = 1, YDMODEL%YRML_GCONF%YGFL%NLIMA
      DO JLEV = 1, KFLEVG
        ! RÃ©initialisation des variables LIMA
        DO JLON = KIDIA, KFDIA
          PLIMAM(JLON, JLEV, JGFL) = PLIMAS(JLON, JLEV, JGFL)*YDCPG_OPTS%ZDTPHY
          PTENDLIMA(JLON, JLEV, JGFL) = PTENDLIMA(JLON, JLEV, JGFL) &
                                    & + (PLIMAS(JLON, JLEV, JGFL) - ZLIMASIN(JLON, JLEV, JGFL))
        ENDDO
      ENDDO
    ENDDO

    !modif de R et CP
    DO JLEV = 1, KFLEVG
      DO JLON = KIDIA, KFDIA
        ZQHGM(JLON, JLEV) = PQHM(JLON, JLEV) + PQGM(JLON, JLEV)
      ENDDO
    ENDDO

    CALL GPRCP_QLIRSG(YDCST, KLON, KIDIA, KFDIA, KFLEVG, &
                    & PQ=PQVM, PQI=PQIM, PQL=PQCM, PQR=PQRM, PQS=PQSM, PQG=ZQHGM, PCP=PCPM, PR=PRHM)

    DO JLEV = 1, KFLEVG
      DO JLON = KIDIA, KFDIA
        PTM(JLON,JLEV) = PTHM(JLON,JLEV)*PEXNREFM(JLON,JLEV)
        ZRHO = YDMF_PHYS_BASE_STATE%YCPG_PHY%PREF(JLON,JLEV)/(PRHM(JLON,JLEV)*PTM(JLON,JLEV))
        PRHODREFM(JLON,JLEV) = ZRHO*PQDM(JLON,JLEV)
      ENDDO
    ENDDO

    !geopotentiel calculation
    DO JLON = KIDIA, KFDIA
      PAPHIM(JLON, KFLEVG) = YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI(JLON, KFLEVG)
    ENDDO

    CALL GPGEO(KLON, KIDIA, KFDIA, KFLEVG, PAPHIM, PAPHIFM, &
    & PTM, PRHM, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%LNPR, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%ALPH, YDGEOMETRY%YRVERT_GEOM)

    !calcul de l'altitude
    DO JLEV = 1, KFLEVG
      DO JLON = KIDIA, KFDIA
        PZZ(JLON,JLEV) = PAPHIM(JLON,JLEV)*ZINVG
        !initialisation de ZZZ_F_
        PZZ_F(JLON,JLEV) = PAPHIFM(JLON,JLEV)*ZINVG
        ! tendency of T
        PTENDT(JLON,JLEV) = PTENDT(JLON,JLEV) + (PTHS(JLON,JLEV) - ZTHSIN(JLON,JLEV))*PEXNREFM(JLON,JLEV)
        PTENDTT(JLON,JLEV) = PTHS(JLON,JLEV) - ZTHSIN(JLON,JLEV)
      ENDDO
    ENDDO

    !inversion niveaux tendances des ri et conversion en qi en multipliant par qd
    DO JRR = 1, YDMODEL%YRML_PHY_MF%YRPARAR%NRR
      DO JLEV = 1, KFLEVG
        DO JLON = KIDIA, KFDIA
          PTENDRA(JLON,JLEV,JRR) = PTENDRA(JLON,JLEV,JRR) + (PRS(JLON,JLEV,JRR) - ZRSIN(JLON,JLEV,JRR))*PQDM(JLON,JLEV)
        ENDDO
      ENDDO
    ENDDO

    !initialisation de PDZZ
    DO JLON = KIDIA, KFDIA
      PDZZ(JLON,1) = PAPHIM(JLON,0)*ZINVG - PZZ(JLON,1)
    ENDDO
    DO JLEV = 2, KFLEVG
      DO JLON = KIDIA,KFDIA
        PDZZ(JLON,JLEV) = PZZ(JLON,JLEV-1) - PZZ(JLON,JLEV)
      ENDDO
    ENDDO

  ELSE

    PTM(:,:)  = YDMF_PHYS_BASE_STATE%T%P(:,:)
    PRHM(:,:) = YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R(:,:)
    PQVM(:,:) = YDMF_PHYS_BASE_STATE%Q%P(:,:)
    PQIM(:,:) = YDMF_PHYS_BASE_STATE%I%P(:,:)
    PQCM(:,:) = YDMF_PHYS_BASE_STATE%L%P(:,:)
    PQRM(:,:) = YDMF_PHYS_BASE_STATE%R%P(:,:)
    PQSM(:,:) = YDMF_PHYS_BASE_STATE%S%P(:,:)
    PQGM(:,:) = YDMF_PHYS_BASE_STATE%G%P(:,:)

    IF (YDMODEL%YRML_PHY_MF%YRPARAR%YRTENDRA%LHAIL) THEN
      PQHM(:,:) = YDMF_PHYS_BASE_STATE%H%P(:,:)
    ELSE
      PQHM(:,:) = 0._JPRB
    ENDIF

    PCPM(:,:) = YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP(:,:)
    PAPHIM(:,:) = YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI(:,:)
    PAPHIFM(:,:) = YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF(:,:)

    DO JLEV = 1, KFLEVG
      DO JLON = KIDIA, KFDIA
        PZZ(JLON, JLEV) = YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI(JLON, JLEV)*ZINVG
      ENDDO
    ENDDO

    !initialisation of PCLFS outside LMICRO to be zero in case LMICRO=F
    YDVARS%A%T1(:,:) = 0._JPRB

  ENDIF ! ADJUSTMENT LMICRO

  IF (YDMODEL%YRML_DIAG%YRLDDH%LFLEXDIA) THEN

    DO JLEV = 1, KFLEVG
      DO JLON=KIDIA, KFDIA
        ZTMPAF(JLON,JLEV) = (PTHS(JLON,JLEV) &
                        & - ZTHSIN(JLON,JLEV))*PEXNREFM(JLON,JLEV) &
                        & * YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP(JLON,JLEV)*PCPM(JLON,JLEV)/YDCST%RG
      ENDDO
    ENDDO
    CALL NEW_ADD_FIELD_3D(YDMODEL%YRML_DIAG%YRMDDH, ZTMPAF, 'TCTADJU', YDDDH)

    DO JRR = 1, YDMODEL%YRML_PHY_MF%YRPARAR%NRR

      CLNAME='T'//CLVARNAME(JRR)//'ADJU'

      DO JLEV = 1, KFLEVG
        DO JLON=KIDIA,KFDIA
          ZTMPAF(JLON,JLEV) = (PRS(JLON,JLEV,JRR) - ZRSIN(JLON,JLEV,JRR))*PQDM(JLON,JLEV) &
                          & * YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP(JLON,JLEV)/YDCST%RG
        ENDDO
      ENDDO
      CALL NEW_ADD_FIELD_3D(YDMODEL%YRML_DIAG%YRMDDH, ZTMPAF, CLNAME, YDDDH)

      DO JLEV = 1, KFLEVG
        DO JLON=KIDIA, KFDIA
          ZTMPAF(JLON,JLEV)=YDVARS%A%T1(JLON,JLEV)*YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP(JLON,JLEV)
        ENDDO
      ENDDO
      CALL NEW_ADD_FIELD_3D(YDMODEL%YRML_DIAG%YRMDDH, ZTMPAF, 'VNT', YDDDH)

    ENDDO

    ! specific to new data flow for diagnostics
    ZCON1(:,:) = 1.0_JPRB
    ZCON2(:,:) = PQDM(:,:)

    DO JLEV = 1, KFLEVG
      DO JLON=KIDIA, KFDIA
        ZCON3(JLON, JLEV) = YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP(JLON, JLEV)*PEXNREFM(JLON, JLEV)
      ENDDO
    ENDDO

    CALL ARO_SUINTBUDGET_OMP(KLON, KFLEVG, ZCON1, ZCON2, ZCON3, YDDDH)

  ENDIF

  IF (JPRD == JPRB) THEN
    ZEPSNEB=1.E-12
  ELSE
    ZEPSNEB=1.E-06
  ENDIF

  ! protect cloudiness from being 0 or 1  (needed for ACRANEB2 and ACNPART)
  DO JLEV=YDCPG_OPTS%KTDIA, KFLEVG
    DO JLON=KIDIA, KFDIA
      ZNEB0(JLON,JLEV) = MAX(ZEPSNEB, MIN(1._JPRB-ZEPSNEB, YDVARS%A%T1(JLON, JLEV)))
    ENDDO
  ENDDO

  ! decorrelation depth for cloud overlaps
  IF (YDMODEL%YRML_PHY_MF%YRPHY%LRNUEXP) THEN
    DO JLON = KIDIA, KFDIA
      ZDECRD(JLON) = YDMODEL%YRML_PHY_MF%YRPHY0%RDECRD1 &
                 & + YDMODEL%YRML_PHY_MF%YRPHY0%RDECRD2*EXP(-((ASIN(YDGEOMVARS%GEMU(JLON)) &
                 & - YDMODEL%YRML_PHY_MF%YRPHY0%RDECRD3*YDMODEL%YRML_GCONF%YRRIP%RDECLI) &
                 & / YDMODEL%YRML_PHY_MF%YRPHY0%RDECRD4)**2)
    ENDDO
  ENDIF

  ! calculate high, medium, low and total cloud cover
  CALL ACNPART(YDCST, YDMODEL%YRML_PHY_MF, &
  & KIDIA, KFDIA, KLON, &
  & YDMODEL%YRML_PHY_MF%YRTOPH%NTNEBU, KFLEVG, &
  & YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREF, &
  & ZDECRD, ZNEB0, &
  & YDMF_PHYS%CLCH, YDMF_PHYS%CLCM, YDMF_PHYS%CLCL, &
  & YDCPG_MISC%CLCT, ZCLCT_RAD)

  END ASSOCIATE

  IF (LHOOK) CALL DR_HOOK('APL_AROME_ADJUST', 1, ZHOOK_HANDLE)

END SUBROUTINE APL_AROME_ADJUST

