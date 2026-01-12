SUBROUTINE ARO_TURB_MNH(PHYEX,                                         &
                      & KKA,KKU,KKL,KIDIA,KFDIA,KLON,KLEV,             &
                      & KRR, KRRL, KRRI, KSV,                          &
                      & KGRADIENTSLEO, KGRADIENTSGOG, CMICRO, PTSTEP,  &
                      & PZZ, PZZF, PZS, PZZTOP,                        &
                      & PRHODJ, PTHVREF,                               &
                      & PSFTH, PSFRV, PSFSV, PSFU, PSFV,               &
                      & PPABSM, PUM, PVM, PWM, PTKEM, PSVM, PSRCM,     &
                      & PTHM, PRM,                                     &
                      & PRUS, PRVS, PRWS, PRTHS, PRRS,                 &
                      & PRSVS, PRTKES, PRTKES_OUT,                     &
                      & PHGRADLEO, PHGRADGOG, PDELTAX, PDELTAY, PSIGS, &
                      & PFLXZTHVMF, PFLXZUMF, PFLXZVMF, PLENGTHM, PLENGTHH, MFMOIST, &
                      & PDRUS_TURB, PDRVS_TURB,                        &
                      & PDRTHLS_TURB, PDRRTS_TURB, PDRSVS_TURB,        &
                      & PDP, PTP, PDPMF, PTPMF, PTDIFF, PTDISS, PEDR,  &
                      & YDDDH,YDLDDH,YDMDDH,&
                      & YDML_PHY_FORCING)


USE PARKIND1                 , ONLY : JPIM, JPRB
USE YOMHOOK                  , ONLY : DR_HOOK, JPHOOK, LHOOK

!     ##########################################################################
!
!!****  * -  compute the turbulence sources and the TKE evolution for Arome
!!
!!
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the turbulence sources
!!      and the TKE evolution for the Arome model
!!
!!
!!**  METHOD
!!    ------
!!      This routine calls the mesoNH turbulence scheme
!!      in its 1DIM configutation.
!!
!!
!!    EXTERNAL
!!    --------
!!      Subroutine TURB (routine de MesoNH)
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_PARAMETERS : contains declarations of parameter variables
!!         JPHEXT       : Horizontal external points number
!!         JPVEXT_TURB       : Vertical external points number
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
!!      2012-02 Y. Seity,  add possibility to run with reversed vertical levels
!!      2015-07 Wim de Rooy  possibility to run with LHARATU=TRUE (Racmo turbulence scheme)
!!      A. Marcel Jan 2025: EDMF contribution to dynamic TKE production
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!

USE MODD_PHYEX               , ONLY : PHYEX_T
USE MODD_LES                 , ONLY : TLES_T
USE MODD_PARAMETERS          , ONLY : JPVEXT_TURB
USE MODD_DIMPHYEX            , ONLY : DIMPHYEX_T
USE MODD_IO                  , ONLY : TFILEDATA
USE MODD_BUDGET              , ONLY : NBUDGET_RH, TBUCONF, TBUDGETDATA_PTR
USE MODE_BUDGET_IAL          , ONLY : TBUDGETDATA_IAL
USE MODI_TURB                , ONLY : TURB
USE MODE_FILL_DIMPHYEX       , ONLY : FILL_DIMPHYEX
USE DDH_MIX                  , ONLY : TYP_DDH
USE YOMLDDH                  , ONLY : TLDDH
USE YOMMDDH                  , ONLY : TMDDH
USE MODEL_PHYSICS_FORCING_MOD, ONLY : MODEL_PHYSICS_FORCING_TYPE
!
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
!
TYPE(PHYEX_T)                     ,INTENT(IN)             :: PHYEX                                       
INTEGER(KIND=JPIM)                ,INTENT(IN)             :: KKA                               ! Index of point near ground
INTEGER(KIND=JPIM)                ,INTENT(IN)             :: KKU                               ! Index of point near top
INTEGER(KIND=JPIM)                ,INTENT(IN)             :: KKL                               ! vert. levels type 1=MNH -1=ARO
INTEGER(KIND=JPIM)                ,INTENT(IN)             :: KIDIA                             
INTEGER(KIND=JPIM)                ,INTENT(IN)             :: KFDIA                             
INTEGER(KIND=JPIM)                ,INTENT(IN)             :: KLON                              ! KFDIA under CPG
INTEGER(KIND=JPIM)                ,INTENT(IN)             :: KLEV                              ! Number of vertical levels
INTEGER(KIND=JPIM)                ,INTENT(IN)             :: KRR                               ! Number of moist variables
INTEGER(KIND=JPIM)                ,INTENT(IN)             :: KRRL                              ! Number of liquide water variables
INTEGER(KIND=JPIM)                ,INTENT(IN)             :: KRRI                              ! Number of ice variables
INTEGER(KIND=JPIM)                ,INTENT(IN)             :: KSV                               ! Number of passive scalar
INTEGER(KIND=JPIM)                ,INTENT(IN)             :: KGRADIENTSLEO                     ! Number of stored horizontal gradients
INTEGER(KIND=JPIM)                ,INTENT(IN)             :: KGRADIENTSGOG                     ! Number of stored horizontal gradients
CHARACTER(LEN=4)                  ,INTENT(IN)             :: CMICRO                            ! Microphysics scheme
REAL(KIND=JPRB)                   ,INTENT(IN)             :: PTSTEP                            ! Time step
REAL(KIND=JPRB)                   ,INTENT(IN)             :: PZZ              (KLON,KLEV)      ! Height of layer boundaries
REAL(KIND=JPRB)                   ,INTENT(IN)             :: PZZF             (KLON,KLEV)      ! Height of level
REAL(KIND=JPRB)                   ,INTENT(IN)             :: PZS              (KLON)           ! Orography
REAL(KIND=JPRB)                   ,INTENT(IN)             :: PZZTOP           (KLON)           ! Height of highest level
REAL(KIND=JPRB)                   ,INTENT(INOUT)          :: PRHODJ           (KLON,KLEV+2)    ! Dry density * Jacobian
REAL(KIND=JPRB)                   ,INTENT(INOUT)          :: MFMOIST          (KLON,KLEV+2)    ! Moist mass flux from Dual scheme ; MFMOIST used in case LHARATU=TRUE
REAL(KIND=JPRB)                   ,INTENT(INOUT)          :: PTHVREF          (KLON,KLEV+2)    ! Virtual Potential
REAL(KIND=JPRB)                   ,INTENT(INOUT)          :: PSFTH            (KLON,1)         ! normal surface fluxes of theta and Rv
REAL(KIND=JPRB)                   ,INTENT(INOUT)          :: PSFRV            (KLON,1)         ! normal surface fluxes of theta and Rv
REAL(KIND=JPRB)                   ,INTENT(INOUT)          :: PSFU             (KLON,1)         ! normal surface fluxes of (u,v) parallel to the orography
REAL(KIND=JPRB)                   ,INTENT(INOUT)          :: PSFV             (KLON,1)         ! normal surface fluxes of (u,v) parallel to the orography
REAL(KIND=JPRB)                   ,INTENT(INOUT)          :: PSFSV            (KLON,KSV)       ! normal surface fluxes of Scalar var.
REAL(KIND=JPRB)                   ,INTENT(INOUT)          :: PPABSM           (KLON,KLEV+2)    ! Pressure at time t-1
REAL(KIND=JPRB)                   ,INTENT(INOUT)          :: PUM              (KLON,KLEV+2)    ! wind components
REAL(KIND=JPRB)                   ,INTENT(INOUT)          :: PVM              (KLON,KLEV+2)    ! wind components
REAL(KIND=JPRB)                   ,INTENT(INOUT)          :: PWM              (KLON,KLEV+2)    ! wind components
REAL(KIND=JPRB)                   ,INTENT(INOUT)          :: PTKEM            (KLON,KLEV+2)    ! TKE
REAL(KIND=JPRB)                   ,INTENT(IN)             :: PSVM             (KLON,KLEV,KSV)  ! passive scal. var.
REAL(KIND=JPRB)                   ,INTENT(INOUT)          :: PSRCM            (KLON,KLEV+2)    ! Second-order flux
REAL(KIND=JPRB)                   ,INTENT(INOUT)          :: PLENGTHM         (KLON,KLEV+2)    ! length scales vdfexcu ; PLENGTHM, PLENGTH used in case LHARATU=true
REAL(KIND=JPRB)                   ,INTENT(INOUT)          :: PLENGTHH         (KLON,KLEV+2)    ! length scales vdfexcu ; PLENGTHM, PLENGTH used in case LHARATU=true
REAL(KIND=JPRB)                   ,INTENT(INOUT)          :: PTHM             (KLON,KLEV+2)    ! pot. temp.
REAL(KIND=JPRB)                   ,INTENT(INOUT)          :: PRM              (KLON,KLEV,KRR)  ! mixing ratio
REAL(KIND=JPRB)                   ,INTENT(INOUT)          :: PRUS             (KLON,KLEV+2)    ! sources of momentum, conservative potential temperature, Turb. Kin. Energy,
REAL(KIND=JPRB)                   ,INTENT(INOUT)          :: PRVS             (KLON,KLEV+2)    ! sources of momentum, conservative potential temperature, Turb. Kin. Energy,
REAL(KIND=JPRB)                   ,INTENT(INOUT)          :: PRWS             (KLON,KLEV+2)    ! sources of momentum, conservative potential temperature, Turb. Kin. Energy,
REAL(KIND=JPRB)                   ,INTENT(INOUT)          :: PRTHS            (KLON,KLEV+2)    
REAL(KIND=JPRB)                   ,INTENT(IN)             :: PRTKES           (KLON,KLEV)      
REAL(KIND=JPRB)                   ,INTENT(OUT)            :: PRTKES_OUT       (KLON,KLEV+2)    
REAL(KIND=JPRB)                   ,INTENT(INOUT)          :: PRRS             (KLON,KLEV,KRR)  ! Source terms for all water kinds, PRRS(:,:,:,1) is used for the conservative mixing ratio
REAL(KIND=JPRB)                   ,INTENT(INOUT)          :: PRSVS            (KLON,KLEV,KSV)  ! Source terms for all passive scalar variables
REAL(KIND=JPRB)                   ,INTENT(OUT)            :: PSIGS            (KLON,KLEV+2)    ! Sigma_s at time t+1 : square root of the variance of the deviation to the saturation
REAL(KIND=JPRB)                   ,INTENT(OUT)            :: PDRUS_TURB       (KLON,KLEV+2)    ! evolution of rhoJ*U   by turbulence only
REAL(KIND=JPRB)                   ,INTENT(OUT)            :: PDRVS_TURB       (KLON,KLEV+2)    ! evolution of rhoJ*V   by turbulence only
REAL(KIND=JPRB)                   ,INTENT(OUT)            :: PDRTHLS_TURB     (KLON,KLEV+2)    ! evolution of rhoJ*thl by turbulence only
REAL(KIND=JPRB)                   ,INTENT(OUT)            :: PDRRTS_TURB      (KLON,KLEV+2)    ! evolution of rhoJ*rt  by turbulence only
REAL(KIND=JPRB)                   ,INTENT(OUT)            :: PDRSVS_TURB      (KLON,KLEV,KSV)  ! evolution of rhoJ*Sv  by turbulence only
REAL(KIND=JPRB)                   ,INTENT(INOUT)          :: PFLXZTHVMF       (KLON,KLEV+2)    
REAL(KIND=JPRB)                   ,INTENT(INOUT)          :: PFLXZUMF         (KLON,KLEV+2)    
REAL(KIND=JPRB)                   ,INTENT(INOUT)          :: PFLXZVMF         (KLON,KLEV+2)    
REAL(KIND=JPRB)                   ,INTENT(OUT)            :: PEDR             (KLON,KLEV+2)    ! EDR
REAL(KIND=JPRB)                   ,INTENT(OUT)            :: PDP              (KLON,KLEV+2)    
REAL(KIND=JPRB)                   ,INTENT(OUT)            :: PTP              (KLON,KLEV+2)    
REAL(KIND=JPRB)                   ,INTENT(OUT)            :: PDPMF            (KLON,KLEV+2)    
REAL(KIND=JPRB)                   ,INTENT(OUT)            :: PTPMF            (KLON,KLEV+2)    
REAL(KIND=JPRB)                   ,INTENT(OUT)            :: PTDIFF           (KLON,KLEV+2)    
REAL(KIND=JPRB)                   ,INTENT(OUT)            :: PTDISS           (KLON,KLEV+2)    
TYPE(TYP_DDH)                     ,INTENT(INOUT) ,TARGET  :: YDDDH                                      
TYPE(TLDDH)                       ,INTENT(IN)    ,TARGET  :: YDLDDH                                     
TYPE(TMDDH)                       ,INTENT(IN)    ,TARGET  :: YDMDDH                                     
TYPE(MODEL_PHYSICS_FORCING_TYPE)  ,INTENT(IN)             :: YDML_PHY_FORCING                           
REAL(KIND=JPRB)                   ,INTENT(INOUT)          :: PHGRADLEO        (KLON,KLEV+2,KGRADIENTSLEO)  ! Horizontal Gradients
REAL(KIND=JPRB)                   ,INTENT(INOUT)          :: PHGRADGOG        (KLON,KLEV+2,KGRADIENTSGOG)  ! Horizontal Gradients
REAL(KIND=JPRB)                   ,INTENT(IN)             :: PDELTAX                             
REAL(KIND=JPRB)                   ,INTENT(IN)             :: PDELTAY                             

#include "wrarom.intfb.h"

TYPE(TBUDGETDATA_IAL), TARGET :: YLBUDGET(NBUDGET_RH) !NBUDGET_RI is the one with the highest number needed for turb
TYPE(TBUDGETDATA_PTR) :: YLBUDGET_PTR(NBUDGET_RH)
TYPE(TFILEDATA) :: ZTFILE !I/O for MesoNH

!*       0.2   Declarations of local variables :

INTEGER(KIND=JPIM) :: JRR,JSV       ! Loop index for the moist and scalar variables
INTEGER(KIND=JPIM) :: JLEV, JLON
INTEGER(KIND=JPIM) :: II
INTEGER(KIND=JPIM) :: KSV_LGBEG, KSV_LGEND ! number of scalar variables

REAL(KIND=JPRB) :: ZDXX(KLON,KLEV+2) ! metric coefficients
REAL(KIND=JPRB) :: ZDYY(KLON,KLEV+2) ! metric coefficients
REAL(KIND=JPRB) :: ZDZZ(KLON,KLEV+2) ! metric coefficients
REAL(KIND=JPRB) :: ZDZX(KLON,KLEV+2) ! metric coefficients
REAL(KIND=JPRB) :: ZDZY(KLON,KLEV+2) ! metric coefficients

INTEGER(KIND=JPIM) :: NSV_LIMA_NR ! TODO LIMA integration : to be sent from above aro_turb_mnh
INTEGER(KIND=JPIM) :: NSV_LIMA_NS ! TODO LIMA integration : to be sent from above aro_turb_mnh
INTEGER(KIND=JPIM) :: NSV_LIMA_NG ! TODO LIMA integration : to be sent from above aro_turb_mnh
INTEGER(KIND=JPIM) :: NSV_LIMA_NH ! TODO LIMA integration : to be sent from above aro_turb_mnh

REAL(KIND=JPRB), POINTER ::  ZDIRCOSXW(:,:)
REAL(KIND=JPRB), POINTER ::  ZDIRCOSYW(:,:)
REAL(KIND=JPRB), POINTER ::  ZDIRCOSZW(:,:)

! Director Cosinus along x, y and z directions at surface w-point
REAL(KIND=JPRB), POINTER    ::  ZCOSSLOPE(:,:)   ! cosinus of the anglebetween i and the slope vector
REAL(KIND=JPRB), POINTER    ::  ZSINSLOPE(:,:)   ! sinus of the angle between i and the slope vector

REAL(KIND=JPRB)    ::  ZSEA_UCU(KLON)   ! u-speed component of sea current 
REAL(KIND=JPRB)    ::  ZSEA_VCU(KLON)    ! v-speed component of sea current

REAL(KIND=JPRB)         :: ZCEI       (KLON,KLEV+2)     
REAL(KIND=JPRB)         :: ZBL_DEPTH  (KLON,1)         
REAL(KIND=JPRB)         :: ZSBL_DEPTH (KLON,1)         
REAL(KIND=JPRB)         :: ZWTH       (KLON,KLEV+2)    ! heat flux
REAL(KIND=JPRB)         :: ZWRC       (KLON,KLEV+2)    ! cloud water flux
REAL(KIND=JPRB)         :: ZWSV       (KLON,KLEV+2,KSV)! scalar flux
REAL(KIND=JPRB)         :: ZSVM       (KLON,KLEV+2,KSV)! scalar flux
REAL(KIND=JPRB)         :: ZRSVS      (KLON,KLEV+2,KSV)! scalar flux
REAL(KIND=JPRB)         :: ZDRSVS_TURB(KLON,KLEV+2,KSV)! scalar flux
REAL(KIND=JPRB)         :: ZZZ        (KLON,KLEV+2)    ! Local value of PZZ
REAL(KIND=JPRB)         :: ZRM        (KLON,KLEV+2,KRR)
REAL(KIND=JPRB)         :: ZRRS       (KLON,KLEV+2,KRR)
REAL(KIND=JPRB), TARGET :: ZERO       (KLON,1)         
REAL(KIND=JPRB), TARGET :: ZONE       (KLON,1)         

REAL(KIND=JPRB) :: ZTWOTSTEP

TYPE(DIMPHYEX_T) :: D

TYPE(TLES_T) :: YLTLES

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
!------------------------------------------------------------------------------
!
!*       1.     PRELIMINARY COMPUTATIONS
!               ------------------------
!
IF (LHOOK) CALL DR_HOOK('ARO_TURB_MNH',0,ZHOOK_HANDLE)

CALL FILL_DIMPHYEX(D, KLON, 1, KLEV+2, JPVEXT_TURB, KIDIA, KFDIA)

YLTLES%LLES=.FALSE.
YLTLES%LLES_CALL=.FALSE.
ZTWOTSTEP=2*PTSTEP
!
!
!
!------------------------------------------------------------------------------
!
!*       2.   INITIALISATION (CAS DU MODELE 1D)
!             ---------------------------------
!
! Fichier I/O pour MesoNH (non-utilise dans AROME)
ZTFILE%LOPENED=.FALSE.

! 2D version of turbulence
KSV_LGBEG=0
KSV_LGEND=0

! tableau a recalculer a chaque pas de temps
! attention, ZDZZ est l'altitude entre deux niveaux (et pas l'Ã¯Â¿Â½paisseur de la couche)

!WRITE(20,*)'sous aro_turb_mnh PZZF', PZZF(1,58:60)
!WRITE(20,*)'sous aro_turb_mnh PZZ', PZZ(1,58:60)



DO JLEV = 1, KLEV
  DO JLON = KIDIA,KFDIA
    ZZZ(JLON,JLEV+1) = PZZ(JLON,JLEV)
  ENDDO
ENDDO

DO JLON = KIDIA,KFDIA
  ZZZ(JLON,1) = PZZTOP(JLON)
ENDDO

DO JLON = KIDIA,KFDIA
  ZDZZ(JLON,KLEV+2) = -999.
  ZDXX(JLON,KLEV+2) = -999.
  ZDYY(JLON,KLEV+2) = -999.
ENDDO

DO JLEV = 2 , KLEV
  DO JLON = KIDIA,KFDIA
    ZDZZ(JLON,JLEV)=PZZF(JLON,JLEV-1)-PZZF(JLON,JLEV)
    ZDXX(JLON,JLEV)=PDELTAX ! For LLEONARD option
    ZDYY(JLON,JLEV)=PDELTAY ! For LLEONARD option
  ENDDO
ENDDO

DO JLON = KIDIA,KFDIA
  ZZZ(JLON,KLEV+2)  = 2*PZZ(JLON,KLEV)-PZZ(JLON,KLEV-1)
  ZDZZ(JLON,1)      = ZZZ(JLON,KKU)-ZZZ(JLON,D%NKE)
  ZDZZ(JLON,KLEV+1) = PZZF(JLON,KLEV)-(1.5*ZZZ(JLON,KLEV+1)-0.5*ZZZ(JLON,KLEV))
ENDDO

! tableaux qui devront etre initialisÃ¯Â¿Â½s plus en amont dans Aladin s'il
! n'existent pas dÃ¯Â¿Â½ja. Dans le cas du 1D, il n'y a pas de relief,
!  ils ont donc des valeurs triviales.

ZERO(:,:) = 0.
ZONE(:,:) = 1.

ZDIRCOSXW=>ZONE(:,:)
ZDIRCOSYW=>ZONE(:,:)
ZDIRCOSZW=>ZONE(:,:)
ZCOSSLOPE=>ZONE(:,:)
ZSINSLOPE=>ZERO(:,:)
ZSEA_UCU=0.
ZSEA_VCU=0.

!------------------------------------------------------------------------------
!
!
!*       4.   MULTIPLICATION PAR RHODJ
!             POUR OBTENIR LES TERMES SOURCES DE MESONH
!
!             -----------------------------------------------

!        WRITE (15,*)'PRUS debut AC_TURB_MNH=',PRUS
!        WRITE (15,*)'PRVS debut AC_TURB_MNH=',PRVS
!        WRITE (15,*)'PRWS debut AC_TURB_MNH=',PRWS
!        WRITE (15,*)'PRTHS debut AC_TURB_MNH=',PRTHS
!        WRITE (15,*)'PRRS debut AC_TURB_MNH=',PRRS

DO JLEV=2,KLEV+1
  DO JLON = KIDIA,KFDIA
    PRUS(JLON,JLEV)   = PRUS(JLON,JLEV)  *PRHODJ(JLON,JLEV)
    PRVS(JLON,JLEV)   = PRVS(JLON,JLEV)  *PRHODJ(JLON,JLEV)
    PRWS(JLON,JLEV)   = PRWS(JLON,JLEV)  *PRHODJ(JLON,JLEV)
    PRTHS(JLON,JLEV)  = PRTHS(JLON,JLEV) *PRHODJ(JLON,JLEV)
    PRTKES_OUT(JLON,JLEV) = PRTKES(JLON,JLEV-1)*PRHODJ(JLON,JLEV)
  ENDDO
ENDDO
DO JRR=1,KRR
  DO JLEV=2,KLEV+1
    DO JLON = KIDIA,KFDIA
      ZRRS(JLON,JLEV,JRR) = PRRS(JLON,JLEV-1,JRR)*PRHODJ(JLON,JLEV)
    ENDDO
    DO JLON = KIDIA,KFDIA
      ZRM(JLON,JLEV,JRR) = PRM(JLON,JLEV-1,JRR)
    ENDDO
  ENDDO
  DO JLON = KIDIA,KFDIA
    ZRRS(JLON,1,JRR     )= ZRRS(JLON,2,JRR)
    ZRRS(JLON,KLEV+2,JRR)= ZRRS(JLON,KLEV+1,JRR)
    ZRM(JLON,1,JRR     )= ZRM(JLON,2,JRR)
    ZRM(JLON,KLEV+2,JRR)= ZRM(JLON,KLEV+1,JRR)
  ENDDO
ENDDO

DO JSV=1,KSV

  DO JLEV=2,KLEV+1
    DO JLON = KIDIA,KFDIA
      ZRSVS(JLON,JLEV,JSV) = PRSVS(JLON,JLEV-1,JSV)*PRHODJ(JLON,JLEV)
    ENDDO
    DO JLON = KIDIA,KFDIA
      ZSVM(JLON,JLEV,JSV) = PSVM(JLON,JLEV-1,JSV)
    ENDDO
  ENDDO

  DO JLON = KIDIA,KFDIA
    ZRSVS(JLON,1,JSV     )= ZRSVS(JLON,2,JSV)
    ZRSVS(JLON,KLEV+2,JSV)= ZRSVS(JLON,KLEV+1,JSV)
    ZSVM(JLON,1,JSV     )= ZSVM(JLON,2,JSV)
    ZSVM(JLON,KLEV+2,JSV)= ZSVM(JLON,KLEV+1,JSV)
  ENDDO
ENDDO

!------------------------------------------------------------------------------
!
!*       3.   Add 2*JPVEXT_TURB values on the vertical
!
!
CALL VERTICAL_EXTEND(PRHODJ)
CALL VERTICAL_EXTEND(PTHVREF)
CALL VERTICAL_EXTEND(PPABSM)
CALL VERTICAL_EXTEND(PUM)
CALL VERTICAL_EXTEND(PVM)
CALL VERTICAL_EXTEND(PWM)
CALL VERTICAL_EXTEND(PTKEM)
PSRCM(:,1)=0.
PSRCM(:,KLEV+2)=0.
CALL VERTICAL_EXTEND(PTHM)
CALL VERTICAL_EXTEND(PFLXZTHVMF)
CALL VERTICAL_EXTEND(PFLXZUMF)
CALL VERTICAL_EXTEND(PFLXZVMF)
IF (PHYEX%TURBN%LHARAT) THEN
  CALL VERTICAL_EXTEND(PLENGTHM)
  CALL VERTICAL_EXTEND(PLENGTHH)
ENDIF
CALL VERTICAL_EXTEND(MFMOIST)
CALL VERTICAL_EXTEND(PRUS)
CALL VERTICAL_EXTEND(PRVS)
CALL VERTICAL_EXTEND(PRWS)
CALL VERTICAL_EXTEND(PRTHS)
CALL VERTICAL_EXTEND(PRTKES_OUT)

!------------------------------------------------------------------------------
!
!
!*       5.   APPEL DE LA TURBULENCE MESONH
!
!         ---------------------------------
!pour AROME, on n'utilise pas les modifs de Mireille pour la turb au bord des nuages
ZCEI=0.0

DO JRR=1, NBUDGET_RH
  YLBUDGET(JRR)%NBUDGET=JRR
  YLBUDGET(JRR)%YDDDH=>YDDDH
  YLBUDGET(JRR)%YDLDDH=>YDLDDH
  YLBUDGET(JRR)%YDMDDH=>YDMDDH
  YLBUDGET_PTR(JRR)%PTR=>YLBUDGET(JRR)
ENDDO
CALL TURB (PHYEX%CST,PHYEX%CSTURB,TBUCONF,PHYEX%TURBN, PHYEX%NEBN, D,YLTLES,&
   & KRR, KRRL, KRRI, PHYEX%MISC%HLBCX, PHYEX%MISC%HLBCY,KGRADIENTSLEO,&
   & KGRADIENTSGOG, PHYEX%MISC%KHALO, &
   & PHYEX%TURBN%NTURBSPLIT, PHYEX%TURBN%LCLOUDMODIFLM, KSV, KSV_LGBEG, KSV_LGEND, &
   & NSV_LIMA_NR, NSV_LIMA_NS, NSV_LIMA_NG, NSV_LIMA_NH,   &
   & PHYEX%MISC%O2D, PHYEX%MISC%ONOMIXLG, PHYEX%MISC%OFLAT, PHYEX%MISC%OCOUPLES, PHYEX%MISC%OBLOWSNOW,& 
   & PHYEX%MISC%OIBM, PHYEX%MISC%OFLYER, PHYEX%MISC%OCOMPUTE_SRC, PHYEX%MISC%XRSNOW, &
   & PHYEX%MISC%OOCEAN,PHYEX%MISC%ODEEPOC, PHYEX%MISC%ODIAG_IN_RUN,   &
   & PHYEX%TURBN%CTURBLEN_CLOUD, CMICRO, PHYEX%MISC%CELEC, &
   & ZTWOTSTEP,ZTFILE,                                      &
   & ZDXX,ZDYY,ZDZZ,ZDZX,ZDZY,ZZZ,          &
   & ZDIRCOSXW,ZDIRCOSYW,ZDIRCOSZW,ZCOSSLOPE,ZSINSLOPE,    &
   & PRHODJ,PTHVREF,PHGRADLEO,PHGRADGOG,PZS,&
   & PSFTH,PSFRV,PSFSV,PSFU,PSFV,                          &
   & ZSEA_UCU,ZSEA_VCU,                                    &
   & PPABSM,PUM,PVM,PWM,PTKEM,ZSVM,PSRCM,                  &
   & PLENGTHM,PLENGTHH,MFMOIST,                            &
   & ZBL_DEPTH,ZSBL_DEPTH,                                 &
   & ZCEI, PHYEX%TURBN%XCEI_MIN, PHYEX%TURBN%XCEI_MAX, PHYEX%TURBN%XCOEF_AMPL_SAT, &
   & PTHM,ZRM, &
   & PRUS,PRVS,PRWS,PRTHS,ZRRS,ZRSVS,PRTKES_OUT,         &
   & PSIGS,                                         &
   & PFLXZTHVMF,PFLXZUMF,PFLXZVMF,ZWTH,ZWRC,ZWSV,PDP,PTP,PTDIFF,PTDISS,&
   & YLBUDGET_PTR, KBUDGETS=SIZE(YLBUDGET_PTR),PEDR=PEDR,PDPMF=PDPMF,PTPMF=PTPMF,&
   & PDRUS_TURB=PDRUS_TURB,PDRVS_TURB=PDRVS_TURB,          &
   & PDRTHLS_TURB=PDRTHLS_TURB,PDRRTS_TURB=PDRRTS_TURB,PDRSVS_TURB=ZDRSVS_TURB)
!
!
!------------------------------------------------------------------------------
!
!
!*       5.   DIVISION PAR RHODJ DES TERMES SOURCES DE MESONH
!             (on obtient des termes homogÃ¯Â¿Â½nes Ã¯Â¿Â½ des tendances)
!
!         -----------------------------------------------

DO JLEV=2,KLEV+1
  DO JLON = KIDIA,KFDIA
    PRUS(JLON,JLEV)   = PRUS(JLON,JLEV)  /PRHODJ(JLON,JLEV)
    PRVS(JLON,JLEV)   = PRVS(JLON,JLEV)  /PRHODJ(JLON,JLEV)
    PRTHS(JLON,JLEV)  = PRTHS(JLON,JLEV) /PRHODJ(JLON,JLEV)
    PRTKES_OUT(JLON,JLEV)  = PRTKES_OUT(JLON,JLEV) /PRHODJ(JLON,JLEV)
    PDRUS_TURB(JLON,JLEV)  = PDRUS_TURB(JLON,JLEV) /PRHODJ(JLON,JLEV)
    PDRVS_TURB(JLON,JLEV)  = PDRVS_TURB(JLON,JLEV) /PRHODJ(JLON,JLEV)
    PDRTHLS_TURB(JLON,JLEV)  = PDRTHLS_TURB(JLON,JLEV) /PRHODJ(JLON,JLEV)
    PDRRTS_TURB(JLON,JLEV)  = PDRRTS_TURB(JLON,JLEV) /PRHODJ(JLON,JLEV)
  ENDDO
ENDDO

DO JRR=1,KRR
  DO JLEV=2,KLEV+1
    DO JLON = KIDIA,KFDIA
      PRRS(JLON,JLEV-1,JRR) = ZRRS(JLON,JLEV,JRR)/PRHODJ(JLON,JLEV)
      PRM(JLON,JLEV-1,JRR) = ZRM(JLON,JLEV,JRR)
    ENDDO
  ENDDO
ENDDO

DO JSV=1,KSV
  DO JLEV=2,KLEV+1
    DO JLON = KIDIA,KFDIA
      PRSVS(JLON,JLEV-1,JSV) = ZRSVS(JLON,JLEV,JSV)/PRHODJ(JLON,JLEV)
      PDRSVS_TURB(JLON,JLEV-1,JSV) = ZDRSVS_TURB(JLON,JLEV,JSV)/PRHODJ(JLON,JLEV)
    ENDDO
  ENDDO
ENDDO

IF (YDML_PHY_FORCING%LMUSCLFA) THEN
   CALL WRAROM(YDML_PHY_FORCING%NMUSCLFA, 'PDP',       PDP,      KLON, D%NKT, D%NKB, D%NKE, D%NKL, .FALSE.)
   CALL WRAROM(YDML_PHY_FORCING%NMUSCLFA, 'PTP',       PTP,      KLON, D%NKT, D%NKB, D%NKE, D%NKL, .FALSE.)
   CALL WRAROM(YDML_PHY_FORCING%NMUSCLFA, 'PTDISS',    PTDISS,   KLON, D%NKT, D%NKB, D%NKE, D%NKL, .FALSE.)
   CALL WRAROM(YDML_PHY_FORCING%NMUSCLFA, 'PTDIFF',    PTDIFF,   KLON, D%NKT, D%NKB, D%NKE, D%NKL, .FALSE.)
   CALL WRAROM(YDML_PHY_FORCING%NMUSCLFA, 'WTHL_TURB', ZWTH,     KLON, D%NKT, D%NKB, D%NKE, D%NKL, .FALSE.) !Flux interpolated on mass point
   CALL WRAROM(YDML_PHY_FORCING%NMUSCLFA, 'WRT_TURB',  ZWRC,     KLON, D%NKT, D%NKB, D%NKE, D%NKL, .FALSE.) !Flux interpolated on mass point
ENDIF

IF (LHOOK) CALL DR_HOOK('ARO_TURB_MNH',1,ZHOOK_HANDLE)

CONTAINS

SUBROUTINE VERTICAL_EXTEND(PX)

! fill extra vetical levels to fit MNH interface

REAL(KIND=JPRB), INTENT(INOUT) :: PX(KLON,KLEV+2)

! NO DR_HOOK, PLEASE ! Rek

PX(:,1     )= PX(:,2)
PX(:,KLEV+2)= PX(:,KLEV+1)

END SUBROUTINE VERTICAL_EXTEND

END SUBROUTINE ARO_TURB_MNH

