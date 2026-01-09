!     ######spl
      SUBROUTINE  ARO_ADJUST_LIMA(PHYEX, &
                                  KKA,KKU,KKL,KLON,KLEV, KIDIA, KFDIA,  KRR, KSV, KTCOUNT,  &
                                  OSUBG_COND, OSIGMAS, &
                                  PTSTEP, PSIGQSAT, &
                                  PZZF, PRHODJ, PRHODREF, PEXNREF,&
                                  PPABSM, PTHT, PRT, PSVT, PSIGS, &
                                  PW_NU, PDTHRAD, &
                                  PMFCONV, PRC_MF, PRI_MF, PCF_MF, PWEIGHT_MF_CLOUD, &
                                  PTHS, PRS,  PSVS, PSRCS, PCLDFR, PICEFR, PPRCFR, &
                                  YDDDH, YDLDDH, YDMDDH, LLIMAINIT, &
                                  PHLC_HRC, PHLC_HCF, PHLI_HRI, PHLI_HCF,&
                                  PHLC_HRC_MF, PHLC_HCF_MF, PHLI_HRI_MF, PHLI_HCF_MF )
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
!     ##########################################################################
!
!!****  * -  compute the  resolved clouds and precipitation
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the  microphysical sources
!!    related to the resolved clouds and precipitation
!!
!!
!!**  METHOD
!!    ------
!!      The main actions of this routine is to call the routines computing the
!!    microphysical sources. Before that:
!!        - it computes the real absolute pressure,
!!        - negative values of the current guess of all mixing ratio are removed.
!!          This is done by a global filling algorithm based on a multiplicative
!!          method (Rood, 1987), in order to conserved the total mass in the
!!          simulation domain.
!!        - Sources are transformed in physical tendencies, by removing the
!!          multiplicative term Rhod*J.
!!        - External points values are filled owing to the use of cyclic
!!          l.b.c., in order to performe computations on the full domain.
!!      After calling to microphysical routines, the physical tendencies are
!!    switched back to prognostic variables.
!!
!!
!!    EXTERNAL
!!    --------
!!      Subroutine FMLOOK: to recover the logical unit number linked to a FMfile
!!      Subroutine SLOW_TERMS: Computes the explicit microphysical sources
!!      Subroutine FAST_TERMS: Performs the saturation adjustment for l
!!      Subroutine RAIN_ICE  : Computes the explicit microphysical sources for i
!!      Subroutine ICE_ADJUST: Performs the saturation adjustment for i+l
!!      MIN_ll,SUM3D_ll : distributed functions equivalent to MIN and SUM
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_PARAMETERS : contains declarations of parameter variables
!!         JPHEXT       : Horizontal external points number
!!         JPVEXT       : Vertical external points number
!!      Module MODD_CST
!!          XP00               ! Reference pressure
!!          XRD                ! Gaz  constant for dry air
!!          XCPD               ! Cpd (dry air)
!!
!!    REFERENCE
!!    ---------
!!
!!      Documentation AROME
!!
!!    AUTHOR
!!    ------
!!    S.Malardel and Y.Seity
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    10/03/03
!!      T. Kovacic  11-05-05, Call to budgets for NEGA1_
!!      S. Riette ice for EDKF
!!      2012-02 Y. Seity,  add possibility to run with reversed vertical levels
!!      2025-02 S. Antoine : Correction of negative values dependent on the number of moments of each species
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!USE MODD_CONF
USE MODD_PHYEX, ONLY: PHYEX_t
USE MODD_PARAMETERS
USE MODD_BUDGET, ONLY: TBUDGETDATA_PTR, NBUDGET_SV1, TBUCONF
USE MODE_BUDGET_IAL, ONLY: TBUDGETDATA_IAL
!
USE MODD_PARAM_LIMA
USE MODD_NSV
USE MODD_CH_AEROSOL, ONLY: NSP, NSOA, NCARB, LORILAM
USE MODD_SALT, ONLY: LSALT
USE MODD_DUST, ONLY: LDUST
!
USE MODI_LIMA_ADJUST_SPLIT
USE MODE_SET_CONC_LIMA
USE MODE_FILL_DIMPHYEX, ONLY: FILL_DIMPHYEX
!
USE DDH_MIX, ONLY  : TYP_DDH
USE YOMLDDH, ONLY  : TLDDH
USE YOMMDDH, ONLY  : TMDDH
!
USE MODD_DIMPHYEX,   ONLY: DIMPHYEX_t
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!

!
TYPE(PHYEX_t),            INTENT(IN)   :: PHYEX
INTEGER,                  INTENT(IN)   :: KKA    !near ground array index
INTEGER,                  INTENT(IN)   :: KKU    !uppest atmosphere array index
INTEGER,                  INTENT(IN)   :: KKL    !vert. levels type 1=MNH -1=ARO
INTEGER,                  INTENT(IN)   :: KLON     !NPROMA under CPG
INTEGER,                  INTENT(IN)   :: KLEV     !Number of vertical levels
INTEGER,                  INTENT(IN)   :: KIDIA    !
INTEGER,                  INTENT(IN)   :: KFDIA    !
INTEGER,                  INTENT(IN)   :: KRR      ! Number of moist variables
INTEGER,                  INTENT(IN)   :: KSV      ! Number of moist variables
INTEGER,                  INTENT(IN)   :: KTCOUNT  ! Temporal loop counter
LOGICAL,                  INTENT(IN)   :: OSUBG_COND ! Switch for Subgrid Cond.
LOGICAL,                  INTENT(IN)   :: OSIGMAS  ! Switch for Sigma_s:
                                        ! use values computed in CONDENSATION
                                        ! or that from turbulence scheme
REAL,                     INTENT(IN)   :: PTSTEP   ! Time step
REAL,                     INTENT(IN)   :: PSIGQSAT ! coeff applied to qsat variance contribution
!
!
REAL, DIMENSION(KLON,KLEV),   INTENT(IN)   :: PZZF     ! Height (z)
REAL, DIMENSION(KLON,KLEV),   INTENT(IN)   :: PRHODJ  !Dry density * Jacobian
REAL, DIMENSION(KLON,KLEV),   INTENT(IN)   :: PRHODREF
REAL, DIMENSION(KLON,KLEV),   INTENT(IN)   :: PEXNREF ! Reference Exner function
!
!
REAL, DIMENSION(KLON,KLEV),   INTENT(IN)   :: PPABSM  ! abs. pressure at time t-dt
REAL, DIMENSION(KLON,KLEV),   INTENT(IN)   :: PTHT    ! Theta at time t
REAL, DIMENSION(KLON,KLEV,KRR), INTENT(INOUT) :: PRT     ! Moist variables at time t
REAL, DIMENSION(KLON,KLEV,KSV), INTENT(INOUT) :: PSVT     ! Moist variables at time t
REAL, DIMENSION(KLON,KLEV),   INTENT(IN)   :: PSIGS   ! Sigma_s at time t
!
REAL, DIMENSION(KLON,KLEV),   INTENT(IN)   :: PW_NU   ! w for CCN activation
REAL, DIMENSION(KLON,KLEV),   INTENT(IN)   :: PDTHRAD ! rad theta tendency for CCN activation
!
REAL, DIMENSION(KLON,KLEV),   INTENT(IN)   :: PMFCONV ! convective mass flux
REAL, DIMENSION(KLON,KLEV),   INTENT(IN)   :: PRC_MF, PRI_MF, PCF_MF, PWEIGHT_MF_CLOUD
!
!
REAL, DIMENSION(KLON,KLEV),   INTENT(INOUT) :: PTHS  ! Theta source
REAL, DIMENSION(KLON,KLEV,KRR), INTENT(INOUT) :: PRS   ! Moist  variable sources
REAL, DIMENSION(KLON,KLEV,KSV), INTENT(INOUT) :: PSVS   ! Moist  variable sources
!
!
REAL, DIMENSION(KLON,KLEV),   INTENT(OUT)   :: PSRCS ! Second-order flux
                                                 ! s'rc'/2Sigma_s2 at time t+1
                                                 ! multiplied by Lambda_3
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT)   :: PCLDFR! Cloud fraction
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT)   :: PICEFR! Cloud fraction
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT)   :: PPRCFR! Cloud fraction
!
CHARACTER(LEN=4)    :: HACTCCN  ! kind of CCN activation
REAL, DIMENSION(KLON,1,KLEV,10)    :: ZSOLORG ![%] solubility fraction of soa
REAL, DIMENSION(KLON,1,KLEV, NSP+NCARB+NSOA) :: ZMI
!
TYPE(TYP_DDH), INTENT(INOUT), TARGET :: YDDDH
TYPE(TLDDH), INTENT(IN), TARGET :: YDLDDH
TYPE(TMDDH), INTENT(IN), TARGET :: YDMDDH
!
LOGICAL,                  INTENT(IN)    :: LLIMAINIT
!
REAL, DIMENSION(KLON,KLEV), INTENT(OUT)   :: PHLC_HRC
REAL, DIMENSION(KLON,KLEV), INTENT(OUT)   :: PHLC_HCF
REAL, DIMENSION(KLON,KLEV), INTENT(OUT)   :: PHLI_HRI
REAL, DIMENSION(KLON,KLEV), INTENT(OUT)   :: PHLI_HCF
REAL, DIMENSION(KLON,KLEV), INTENT(IN)    :: PHLC_HRC_MF
REAL, DIMENSION(KLON,KLEV), INTENT(IN)    :: PHLC_HCF_MF
REAL, DIMENSION(KLON,KLEV), INTENT(IN)    :: PHLI_HRI_MF
REAL, DIMENSION(KLON,KLEV), INTENT(IN)    :: PHLI_HCF_MF
!
!*       0.2   Declarations of local variables :

!
INTEGER :: JRR, JLON, JLEV          ! Loop index for the moist and scalar variables
!
REAL, DIMENSION(KLON,KLEV) :: ZT,ZLV,ZLS,ZCPH
REAL, DIMENSION(KLON,KLEV) :: ZCOR
                                    ! for the correction of negative rv
REAL, DIMENSION(KLON,KLEV) :: ZZZ
                                    ! model layer height
REAL  :: ZMASSTOT                   ! total mass  for one water category
                                    ! including the negative values
REAL  :: ZMASSPOS                   ! total mass  for one water category
                                    ! after removing the negative values
REAL  :: ZRATIO                     ! ZMASSTOT / ZMASSCOR
!
TYPE(TBUDGETDATA_IAL), DIMENSION(NBUDGET_SV1+NSV_LIMA-1), TARGET :: YLBUDGET
TYPE(TBUDGETDATA_PTR), DIMENSION(NBUDGET_SV1+NSV_LIMA-1) :: YLBUDGET_PTR
TYPE(DIMPHYEX_t) :: YLDIMPHYEX
REAL, DIMENSION(KLON) :: ZSIGQSAT, ZICE_CLD_WGT
!
!------------------------------------------------------------------------------
!
!*       1.     PRELIMINARY COMPUTATIONS
!               ------------------------
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('ARO_ADJUST_LIMA',0,ZHOOK_HANDLE)

CALL FILL_DIMPHYEX(YLDIMPHYEX, KLON, 1, KLEV, 0, KIDIA, KFDIA)

ZSIGQSAT(:) = PHYEX%NEBN%VSIGQSAT
HACTCCN='    '
ZMI=0.
ZSOLORG=0.
!
!*       2.     TRANSFORMATION INTO PHYSICAL TENDENCIES
!               ---------------------------------------
!
!
! personal comment:  tranfering these variables to the
!                    microphysical routines would save
!                    computing time
!
DO JLEV = 1, KLEV
  DO JLON = KIDIA, KFDIA
    ZT(JLON,JLEV) = PTHT(JLON,JLEV)*PEXNREF(JLON,JLEV)
  ENDDO
ENDDO

DO JLEV = 1, KLEV
  DO JLON = KIDIA, KFDIA
    ZLV(JLON,JLEV) = PHYEX%CST%XLVTT + (PHYEX%CST%XCPV - PHYEX%CST%XCL)*(ZT(JLON,JLEV)-PHYEX%CST%XTT)
  ENDDO
ENDDO

DO JLEV = 1, KLEV
  DO JLON = KIDIA, KFDIA
    ZLS(JLON,JLEV) = PHYEX%CST%XLSTT + (PHYEX%CST%XCPV - PHYEX%CST%XCI)*(ZT(JLON,JLEV)-PHYEX%CST%XTT)
  ENDDO
ENDDO

DO JLEV = 1, KLEV
  DO JLON = KIDIA, KFDIA
    ZCPH(JLON,JLEV) = PHYEX%CST%XCPD + PHYEX%CST%XCPV*2.*PTSTEP*PRS(JLON,JLEV,1)
  ENDDO
ENDDO

!set concentration for LIMA
DO JRR = 1, KRR
  DO JLEV = 1, KLEV
    DO JLON = KIDIA, KFDIA
      PRS(JLON, JLEV, JRR) = PRS(JLON, JLEV, JRR) * 2.*PTSTEP
    ENDDO
  ENDDO
ENDDO

DO JRR = 1, KSV
  DO JLEV = 1, KLEV
    DO JLON = KIDIA, KFDIA
      PSVS(JLON, JLEV, JRR) = PSVS(JLON, JLEV, JRR) * 2.*PTSTEP
    ENDDO
  ENDDO
ENDDO

CALL SET_CONC_LIMA (PHYEX%TNSV, YLDIMPHYEX,KRR,1,'ICE3',PRHODREF,PRT,PSVT, .NOT. LLIMAINIT)
CALL SET_CONC_LIMA (PHYEX%TNSV, YLDIMPHYEX,KRR,1,'ICE3',PRHODREF,PRS,PSVS, .NOT. LLIMAINIT)

DO JRR = 1, KRR
  DO JLEV = 1, KLEV
    DO JLON = KIDIA, KFDIA
      PRS(JLON, JLEV, JRR) = PRS(JLON, JLEV, JRR) / (2.*PTSTEP)
    ENDDO
  ENDDO
ENDDO

DO JRR = 1, KSV
  DO JLEV = 1, KLEV
    DO JLON = KIDIA, KFDIA
      PSVS(JLON, JLEV, JRR) = PSVS(JLON, JLEV, JRR) / (2.*PTSTEP)
    ENDDO
  ENDDO
ENDDO

!
!*       3.     REMOVE NEGATIVE VALUES
!               ----------------------
!
!*       3.1    Non local correction for precipitating species (Rood 87)
!
DO JRR = 3,KRR
   SELECT CASE (JRR)
   CASE(5,6,7) ! snow, graupel and hail

      IF ( MINVAL( PRS(:,:,JRR)) < 0.0 ) THEN
! For AROME, we cannot use MAX_ll so that according to JPP's advises
!  we only correct negative values but not the total mass
! compute the total water mass computation
!
!          ZMASSTOT = MAX( 0. , SUM( PRS(:,:,JRR) ))
!
! remove the negative values
!
         PRS(:,:,JRR) = MAX( 0., PRS(:,:,JRR) )
!
! compute the new total mass
!
!          ZMASSPOS = MAX(1.E-60,SUM( PRS(:,:,JRR) ))
!
! correct again in such a way to conserve the total mass
!
!          ZRATIO = ZMASSTOT / ZMASSPOS
!          PRS(:,:,JRR) = PRS(:,:,JRR) * ZRATIO

      END IF
   END SELECT
END DO
!
!*       3.2    Correct negative values
!
! Correction where rc<0
DO JLEV = 1, KLEV
  DO JLON = KIDIA, KFDIA
    IF (NMOM_C .EQ. 1 .AND. PRS(JLON,JLEV,2) < 1.E-15) THEN
      PRS(JLON, JLEV, 1) = PRS(JLON, JLEV, 1) + PRS(JLON, JLEV, 2)
      PTHS(JLON,JLEV) = PTHS(JLON, JLEV) - PRS(JLON, JLEV, 2) * ZLV(JLON, JLEV) &
                    & / ZCPH(JLON,JLEV) / PEXNREF(JLON,JLEV)
      PRS(JLON, JLEV, 2)  = 0.0
    ELSEIF (NMOM_C .GE. 2 .AND. (PRS(JLON,JLEV,2) < 1.E-15 .OR. PSVS(JLON,JLEV,NSV_LIMA_NC) < 1.E-15)) THEN
      PRS(JLON, JLEV, 1) = PRS(JLON, JLEV, 1) + PRS(JLON, JLEV, 2)
      PTHS(JLON,JLEV) = PTHS(JLON, JLEV) - PRS(JLON, JLEV, 2) * ZLV(JLON, JLEV) &
                    & / ZCPH(JLON,JLEV) / PEXNREF(JLON,JLEV)
      PRS(JLON, JLEV, 2)  = 0.0
      PSVS(JLON, JLEV, NSV_LIMA_NC) = 0.0
    ENDIF
  ENDDO
ENDDO

! Correction where rr<0
DO JLEV = 1, KLEV
  DO JLON = KIDIA, KFDIA
    IF (NMOM_R .EQ. 1 .AND. PRS(JLON,JLEV,3) < 1.E-15) THEN
      PRS(JLON, JLEV, 1) = PRS(JLON, JLEV, 1) + PRS(JLON, JLEV, 3)
      PTHS(JLON, JLEV) = PTHS(JLON, JLEV) - PRS(JLON, JLEV, 3) * ZLV(JLON, JLEV) &
                     & / ZCPH(JLON, JLEV) / PEXNREF(JLON, JLEV)
      PRS(JLON, JLEV, 3)  = 0.0
    ELSEIF (NMOM_R .GE. 2 .AND. (PRS(JLON,JLEV,3) < 1.E-15 .OR. PSVS(JLON,JLEV,NSV_LIMA_NR) < 1.E-15)) THEN
      PRS(JLON, JLEV, 1) = PRS(JLON, JLEV, 1) + PRS(JLON, JLEV, 3)
      PTHS(JLON, JLEV) = PTHS(JLON, JLEV) - PRS(JLON, JLEV, 3) * ZLV(JLON, JLEV) &
                     & / ZCPH(JLON, JLEV) / PEXNREF(JLON, JLEV)
      PRS(JLON, JLEV, 3)  = 0.0
      PSVS(JLON, JLEV, NSV_LIMA_NR) = 0.0
    ENDIF
  ENDDO
ENDDO

! Correction where ri<0
DO JLEV = 1, KLEV
  DO JLON = KIDIA, KFDIA
    IF (NMOM_I .EQ. 1 .AND. PRS(JLON, JLEV, 4) < 1.E-15) THEN
      PRS(JLON, JLEV, 1) = PRS(JLON, JLEV, 1) + PRS(JLON, JLEV, 4)
      PTHS(JLON, JLEV) = PTHS(JLON, JLEV) - PRS(JLON, JLEV, 4) * ZLS(JLON, JLEV) / ZCPH(JLON, JLEV) / PEXNREF(JLON, JLEV)
      PRS(JLON, JLEV, 4)  = 0.0
    ELSEIF (NMOM_I .GE. 2 .AND. (PRS(JLON, JLEV, 4) < 1.E-15 .AND. PSVS(JLON,JLEV,NSV_LIMA_NI) < 1.E-15)) THEN
      PRS(JLON, JLEV, 1) = PRS(JLON, JLEV, 1) + PRS(JLON, JLEV, 4)
      PTHS(JLON, JLEV) = PTHS(JLON, JLEV) - PRS(JLON, JLEV, 4) * ZLS(JLON, JLEV) / ZCPH(JLON, JLEV) / PEXNREF(JLON, JLEV)
      PRS(JLON, JLEV, 4)  = 0.0
      PSVS(JLON, JLEV, NSV_LIMA_NI) = 0.0
    ENDIF
  ENDDO
ENDDO
!
DO JRR = 1, KSV
  DO JLEV = 1, KLEV
    DO JLON = KIDIA, KFDIA
      PSVS(JLON, JLEV, JRR) = MAX(0.0, PSVS(JLON, JLEV, JRR))
    ENDDO
  ENDDO
ENDDO
!
!
!*       3.3  STORE THE BUDGET TERMS
!            ----------------------
!

DO JRR = 1, NBUDGET_SV1+NSV_LIMA-1
   YLBUDGET(JRR)%NBUDGET=JRR
   YLBUDGET(JRR)%YDDDH=>YDDDH
   YLBUDGET(JRR)%YDLDDH=>YDLDDH
   YLBUDGET(JRR)%YDMDDH=>YDMDDH
   YLBUDGET_PTR(JRR)%PTR=>YLBUDGET(JRR)
ENDDO
!
!-------------------------------------------------------------------------------
!

!*       9.     MIXED-PHASE MICROPHYSICAL SCHEME (WITH 3 ICE SPECIES)
!               -----------------------------------------------------
!
!
!*       9.2    Perform the saturation adjustment over cloud ice and cloud water
!
    ZZZ =  PZZF

    CALL LIMA_ADJUST_SPLIT(LIMAP=PHYEX%PARAM_LIMA, LIMAW=PHYEX%PARAM_LIMA_WARM, TNSV=PHYEX%TNSV, &
         D=YLDIMPHYEX, CST=PHYEX%CST, NEBN=PHYEX%NEBN, TURBN=PHYEX%TURBN, &
         BUCONF=TBUCONF, TBUDGETS=YLBUDGET_PTR, KBUDGETS=SIZE(YLBUDGET_PTR), &
         KRR=KRR, HCONDENS=PHYEX%NEBN%CCONDENS, HLAMBDA3=PHYEX%NEBN%CLAMBDA3, &
         KCARB=NCARB, KSOA=NSOA, KSP=NSP, ODUST=LDUST, OSALT=LSALT, OORILAM=LORILAM, &
         OSUBG_COND=OSUBG_COND, OSIGMAS=OSIGMAS, PTSTEP=2*PTSTEP, PSIGQSAT=ZSIGQSAT, &
         PRHODREF=PRHODREF, PRHODJ=PRHODJ, PEXNREF=PEXNREF, PSIGS=PSIGS, OMFCONV=.TRUE., PMFCONV=PMFCONV, &
         PPABST=PPABSM, PZZ=ZZZ, ODTHRAD=.TRUE., PDTHRAD=PDTHRAD, PW_NU=PW_NU, &
         PRT=PRT, PRS=PRS, PSVT=PSVT, PSVS=PSVS, &
         HACTCCN=HACTCCN, PAERO=PSVT, PSOLORG=ZSOLORG, PMI=ZMI, & 
         PTHS=PTHS, OCOMPUTE_SRC=.TRUE., PSRCS=PSRCS, PCLDFR=PCLDFR, PICEFR=PICEFR, &
         PRC_MF=PRC_MF, PRI_MF=PRI_MF, PCF_MF=PCF_MF, &
         PICE_CLD_WGT=ZICE_CLD_WGT, PWEIGHT_MF_CLOUD=PWEIGHT_MF_CLOUD, &
         PHLC_HRC=PHLC_HRC, PHLC_HCF=PHLC_HCF, PHLI_HRI=PHLI_HRI, PHLI_HCF=PHLI_HCF,            &
         PHLC_HRC_MF=PHLC_HRC_MF, PHLC_HCF_MF=PHLC_HCF_MF, PHLI_HRI_MF=PHLI_HRI_MF, PHLI_HCF_MF=PHLI_HCF_MF )
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('ARO_ADJUST_LIMA',1,ZHOOK_HANDLE)
END SUBROUTINE ARO_ADJUST_LIMA
