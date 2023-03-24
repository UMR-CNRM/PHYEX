!     ######spl
      SUBROUTINE  ARO_TURB_MNH(PHYEX, &
                KKA,KKU,KKL,KLON,KLEV,KRR,KRRL,KRRI,KSV,              &
                KGRADIENTS, CMICRO, PTSTEP,                           &
                PZZ, PZZF, PZZTOP,                                    &
                PRHODJ, PTHVREF,                                      &
                PSFTH,PSFRV,PSFSV,PSFU,PSFV,                          &
                PPABSM,PUM,PVM,PWM,PTKEM,PEPSM,PSVM,PSRCM,            &
                PTHM,PRM,                                &
                PRUS,PRVS,PRWS,PRTHS,PRRS,PRSVSIN,PRSVS,PRTKES,PRTKES_OUT,PREPSS, &
                PHGRAD,PSIGS,                                         &
                PFLXZTHVMF,PLENGTHM,PLENGTHH,MFMOIST,                 &
                PDRUS_TURB,PDRVS_TURB,                                &
                PDRTHLS_TURB,PDRRTS_TURB,PDRSVS_TURB,                 &
                PDP,PTP,PTPMF,PTDIFF,PTDISS,PEDR,YDDDH,YDLDDH,YDMDDH)


      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
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
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PHYEX, ONLY: PHYEX_t
USE MODD_LES, ONLY: TLES_t
USE MODD_PARAMETERS,     ONLY: JPVEXT_TURB
USE MODD_DIMPHYEX,       ONLY: DIMPHYEX_t
USE MODD_IO, ONLY: TFILEDATA
USE MODD_BUDGET, ONLY: NBUDGET_RI, TBUDGETDATA, TBUCONF
!
USE MODI_TURB
!
USE MODE_FILL_DIMPHYEX, ONLY: FILL_DIMPHYEX
!
USE DDH_MIX, ONLY  : TYP_DDH
USE YOMLDDH, ONLY  : TLDDH
USE YOMMDDH, ONLY  : TMDDH

!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
!
TYPE(PHYEX_t),            INTENT(IN)   :: PHYEX
INTEGER,                  INTENT(IN)   :: KLON  !KFDIA under CPG
INTEGER,                  INTENT(IN)   :: KLEV  !Number of vertical levels
INTEGER,                  INTENT(IN)   :: KKA   !Index of point near ground
INTEGER,                  INTENT(IN)   :: KKU   !Index of point near top
INTEGER,                  INTENT(IN)   :: KKL   !vert. levels type 1=MNH -1=ARO
INTEGER,                  INTENT(IN)   :: KRR      ! Number of moist variables
INTEGER,                  INTENT(IN)   :: KRRL      ! Number of liquide water variables
INTEGER,                  INTENT(IN)   :: KRRI      ! Number of ice variables
INTEGER,                  INTENT(IN)   :: KSV     ! Number of passive scalar
INTEGER,                  INTENT(IN)   :: KGRADIENTS  ! Number of stored horizontal gradients
CHARACTER (LEN=4),        INTENT(IN)   :: CMICRO  ! Microphysics scheme
REAL,                     INTENT(IN)   :: PTSTEP   ! Time step
!
!
REAL, DIMENSION(KLON,1,KLEV),   INTENT(IN)   :: PZZ     ! Height of layer boundaries
REAL, DIMENSION(KLON,1,KLEV),   INTENT(IN)   :: PZZF    ! Height of level
REAL, DIMENSION(KLON),          INTENT(IN)   :: PZZTOP  ! Height of highest level

REAL, DIMENSION(KLON,1,KLEV+2),   INTENT(INOUT)   :: PRHODJ  !Dry density * Jacobian
! MFMOIST used in case LHARATU=TRUE
REAL, DIMENSION(KLON,1,KLEV+2),   INTENT(INOUT)   :: MFMOIST  !Moist mass flux from Dual scheme
REAL, DIMENSION(KLON,1,KLEV+2), INTENT(INOUT)     ::  PTHVREF   ! Virtual Potential
                                        ! Temperature of the reference state
!
REAL, DIMENSION(KLON,1),   INTENT(INOUT)      ::  PSFTH,PSFRV
! normal surface fluxes of theta and Rv
REAL, DIMENSION(KLON,1),   INTENT(INOUT)      ::   PSFU,PSFV
! normal surface fluxes of (u,v) parallel to the orography
REAL, DIMENSION(KLON,1,KSV), INTENT(INOUT)      ::  PSFSV
! normal surface fluxes of Scalar var.
!
!    prognostic variables at t- deltat
REAL, DIMENSION(KLON,1,KLEV+2),   INTENT(INOUT) ::  PPABSM      ! Pressure at time t-1
REAL, DIMENSION(KLON,1,KLEV+2),   INTENT(INOUT) ::  PUM,PVM,PWM ! wind components
REAL, DIMENSION(KLON,1,KLEV+2),   INTENT(INOUT) ::  PTKEM       ! TKE
REAL, DIMENSION(0,0,0),   INTENT(INOUT) ::  PEPSM       ! dissipation of TKE
REAL, DIMENSION(KLON,1,KLEV,KSV), INTENT(INOUT) ::  PSVM        ! passive scal. var.
REAL, DIMENSION(KLON,1,KLEV+2),   INTENT(INOUT) ::  PSRCM       ! Second-order flux
                      ! s'rc'/2Sigma_s2 at time t-1 multiplied by Lambda_3
!
! PLENGTHM, PLENGTH used in case LHARATU=true
REAL, DIMENSION(KLON,1,KLEV+2), INTENT(INOUT)   :: PLENGTHM, PLENGTHH  ! length scales vdfexcu

!
!   thermodynamical variables which are transformed in conservative var.
REAL, DIMENSION(KLON,1,KLEV+2),   INTENT(INOUT) ::  PTHM       ! pot. temp.
REAL, DIMENSION(KLON,1,KLEV,KRR), INTENT(INOUT) ::  PRM         ! mixing ratio
!
! sources of momentum, conservative potential temperature, Turb. Kin. Energy,
! TKE dissipation
REAL, DIMENSION(KLON,1,KLEV+2),   INTENT(INOUT) ::  PRUS,PRVS,PRWS
REAL, DIMENSION(KLON,1,KLEV+2),   INTENT(INOUT) :: PRTHS
REAL, DIMENSION(KLON,1,KLEV),     INTENT(IN)  ::  PRTKES
REAL, DIMENSION(KLON,1,KLEV+2),   INTENT(OUT) ::  PRTKES_OUT
REAL, DIMENSION(0,0,0)        , INTENT(INOUT) ::PREPSS
! Source terms for all water kinds, PRRS(:,:,:,1) is used for the conservative
! mixing ratio
REAL, DIMENSION(KLON,1,KLEV,KRR), INTENT(INOUT) ::  PRRS
! Source terms for all passive scalar variables
REAL, DIMENSION(KLON,1,KLEV,KSV), INTENT(IN)  ::  PRSVSIN
REAL, DIMENSION(KLON,1,KLEV,KSV), INTENT(OUT) ::  PRSVS
! Sigma_s at time t+1 : square root of the variance of the deviation to the
! saturation
REAL, DIMENSION(KLON,1,KLEV+2), INTENT(OUT)     ::  PSIGS
REAL, DIMENSION(KLON,1,KLEV+2), INTENT(OUT)     ::  PDRUS_TURB   ! evolution of rhoJ*U   by turbulence only
REAL, DIMENSION(KLON,1,KLEV+2), INTENT(OUT)     ::  PDRVS_TURB   ! evolution of rhoJ*V   by turbulence only
REAL, DIMENSION(KLON,1,KLEV+2), INTENT(OUT)     ::  PDRTHLS_TURB ! evolution of rhoJ*thl by turbulence only
REAL, DIMENSION(KLON,1,KLEV+2), INTENT(OUT)     ::  PDRRTS_TURB  ! evolution of rhoJ*rt  by turbulence only
REAL, DIMENSION(KLON,1,KLEV,KSV), INTENT(OUT)   ::  PDRSVS_TURB  ! evolution of rhoJ*Sv  by turbulence only
REAL, DIMENSION(KLON,1,KLEV+2), INTENT(INOUT)      ::  PFLXZTHVMF
REAL, DIMENSION(KLON,1,KLEV+2),   INTENT(OUT) ::  PEDR       ! EDR
!
REAL, DIMENSION(KLON,1,KLEV+2),  INTENT(OUT)   :: PDP, PTP, PTPMF, PTDIFF, PTDISS
!                                                !for TKE DDH budgets
!
TYPE(TYP_DDH), INTENT(INOUT), TARGET   :: YDDDH
TYPE(TLDDH),   INTENT(IN), TARGET      :: YDLDDH
TYPE(TMDDH),   INTENT(IN), TARGET      :: YDMDDH
!
!
TYPE(TBUDGETDATA), DIMENSION(NBUDGET_RI) :: YLBUDGET !NBUDGET_RI is the one with the highest number needed for turb
TYPE(TFILEDATA) :: ZTFILE !I/O for MesoNH
!*       0.2   Declarations of local variables :
!
INTEGER :: JRR,JSV       ! Loop index for the moist and scalar variables
INTEGER :: JL, JK, JLON
!
INTEGER ::II
!
!
INTEGER   :: KSV_LGBEG, KSV_LGEND ! number of scalar variables
!
REAL, DIMENSION(KLON,1,KLEV+2)   :: ZDXX,ZDYY,ZDZZ,ZDZX,ZDZY
                                        ! metric coefficients
INTEGER :: NSV_LIMA_NR, NSV_LIMA_NS, NSV_LIMA_NG, NSV_LIMA_NH ! TODO LIMA integration : to be sent from above aro_turb_mnh
REAL, DIMENSION(KLON,1) :: PZS ! TODO: to be sent from above aro_turb_mnh
REAL, POINTER ::  ZDIRCOSXW(:,:), ZDIRCOSYW(:,:), ZDIRCOSZW(:,:)
! Director Cosinus along x, y and z directions at surface w-point
REAL, POINTER    ::  ZCOSSLOPE(:,:)   ! cosinus of the anglebetween i and the slope vector
REAL, POINTER    ::  ZSINSLOPE(:,:)   ! sinus of the angle between i and the slope vector

REAL,DIMENSION(KLON,1,KLEV+2)         :: ZCEI
REAL, DIMENSION(KLON,1)               :: ZBL_DEPTH, ZSBL_DEPTH
REAL,DIMENSION(KLON,1,KLEV+2)         :: ZWTH       ! heat flux
REAL,DIMENSION(KLON,1,KLEV+2)         :: ZWRC       ! cloud water flux
REAL,DIMENSION(KLON,1,KLEV+2,KSV)     :: ZWSV,ZSVM,ZRSVS,ZDRSVS_TURB       ! scalar flux
REAL,DIMENSION(KLON,1,KLEV+2)         :: ZZZ        ! Local value of PZZ
REAL,DIMENSION(KLON,1,KLEV+2,KRR)     :: ZRM,ZRRS
REAL,DIMENSION(KLON,1,KLEV+2,KGRADIENTS) :: PHGRAD    ! Horizontal Gradients
!
REAL, DIMENSION(KLON,1), TARGET     ::  ZERO, ZONE
REAL :: ZTWOTSTEP
!
TYPE(DIMPHYEX_t) :: YLDIMPHYEX
TYPE(TLES_t) :: YLTLES
!
!------------------------------------------------------------------------------
!
!*       1.     PRELIMINARY COMPUTATIONS
!               ------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('ARO_TURB_MNH',0,ZHOOK_HANDLE)
CALL FILL_DIMPHYEX(YLDIMPHYEX, KLON, 1, KLEV+2, JPVEXT_TURB, KLON)
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
! attention, ZDZZ est l'altitude entre deux niveaux (et pas l'�paisseur de la couche)

!WRITE(20,*)'sous aro_turb_mnh PZZF', PZZF(1,1,58:60)
!WRITE(20,*)'sous aro_turb_mnh PZZ', PZZ(1,1,58:60)



ZZZ(YLDIMPHYEX%NIB:YLDIMPHYEX%NIE,1,2:KLEV+1)=PZZ(YLDIMPHYEX%NIB:YLDIMPHYEX%NIE,1,1:KLEV)
ZZZ(YLDIMPHYEX%NIB:YLDIMPHYEX%NIE,1,1) = PZZTOP(YLDIMPHYEX%NIB:YLDIMPHYEX%NIE)
ZDZZ(YLDIMPHYEX%NIB:YLDIMPHYEX%NIE,1,KLEV+2)=-999.

DO JK = 2 , KLEV
  DO JL = YLDIMPHYEX%NIB,YLDIMPHYEX%NIE
    ZDZZ(JL,1,JK)=PZZF(JL,1,JK-1)-PZZF(JL,1,JK)
  ENDDO
ENDDO

DO JL = YLDIMPHYEX%NIB,YLDIMPHYEX%NIE
  ZZZ(JL,1,KLEV+2) = 2*PZZ(JL,1,KLEV)-PZZ(JL,1,KLEV-1)
  ZDZZ(JL,1,1)=ZZZ(JL,1,KKU)-ZZZ(JL,1,YLDIMPHYEX%NKE)
  ZDZZ(JL,1,KLEV+1)=PZZF(JL,1,KLEV)-(1.5*ZZZ(JL,1,KLEV+1)-0.5*ZZZ(JL,1,KLEV))
ENDDO

! tableaux qui devront etre initialis�s plus en amont dans Aladin s'il
! n'existent pas d�ja. Dans le cas du 1D, il n'y a pas de relief,
!  ils ont donc des valeurs triviales.

ZERO(:,:) = 0.
ZONE(:,:) = 1.

ZDIRCOSXW=>ZONE(:,:)
ZDIRCOSYW=>ZONE(:,:)
ZDIRCOSZW=>ZONE(:,:)
ZCOSSLOPE=>ZONE(:,:)
ZSINSLOPE=>ZERO(:,:)

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

DO JK=2,KLEV+1
  DO JL = 1,KLON
    PRUS(JL,1,JK)   = PRUS(JL,1,JK)  *PRHODJ(JL,1,JK)
    PRVS(JL,1,JK)   = PRVS(JL,1,JK)  *PRHODJ(JL,1,JK)
    PRWS(JL,1,JK)   = PRWS(JL,1,JK)  *PRHODJ(JL,1,JK)
    PRTHS(JL,1,JK)  = PRTHS(JL,1,JK) *PRHODJ(JL,1,JK)
    PRTKES_OUT(JL,1,JK) = PRTKES(JL,1,JK-1)*PRHODJ(JL,1,JK)
  ENDDO
ENDDO
DO JRR=1,KRR
  DO JK=2,KLEV+1
    DO JL = 1,KLON
      ZRRS(JL,1,JK,JRR) = PRRS(JL,1,JK-1,JRR)*PRHODJ(JL,1,JK)
    ENDDO
    ZRM(:,1,JK,JRR) = PRM(:,1,JK-1,JRR)
  ENDDO
  ZRRS(:,1,1,JRR     )= ZRRS(:,1,2,JRR)
  ZRRS(:,1,KLEV+2,JRR)= ZRRS(:,1,KLEV+1,JRR)
  ZRM(:,1,1,JRR     )= ZRM(:,1,2,JRR)
  ZRM(:,1,KLEV+2,JRR)= ZRM(:,1,KLEV+1,JRR)
ENDDO
DO JSV=1,KSV
  DO JK=2,KLEV+1
    DO JL = 1,KLON
      ZRSVS(JL,1,JK,JSV) = PRSVSIN(JL,1,JK-1,JSV)*PRHODJ(JL,1,JK)
    ENDDO
    ZSVM(:,1,JK,JSV) = PSVM(:,1,JK-1,JSV)
  ENDDO
  ZRSVS(:,1,1,JSV     )= ZRSVS(:,1,2,JSV)
  ZRSVS(:,1,KLEV+2,JSV)= ZRSVS(:,1,KLEV+1,JSV)
  ZSVM(:,1,1,JSV     )= ZSVM(:,1,2,JSV)
  ZSVM(:,1,KLEV+2,JSV)= ZSVM(:,1,KLEV+1,JSV)
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
PSRCM(:,:,1)=0.
PSRCM(:,:,KLEV+2)=0.
CALL VERTICAL_EXTEND(PTHM)
CALL VERTICAL_EXTEND(PFLXZTHVMF)
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

DO JRR=1, NBUDGET_RI
  YLBUDGET(JRR)%NBUDGET=JRR
  YLBUDGET(JRR)%YDDDH=>YDDDH
  YLBUDGET(JRR)%YDLDDH=>YDLDDH
  YLBUDGET(JRR)%YDMDDH=>YDMDDH
ENDDO
CALL TURB (PHYEX%CST,PHYEX%CSTURB,TBUCONF,PHYEX%TURBN, PHYEX%NEBN, YLDIMPHYEX,YLTLES,&
   & PHYEX%MISC%KMI, KRR, KRRL, KRRI, PHYEX%MISC%HLBCX, PHYEX%MISC%HLBCY, KGRADIENTS, PHYEX%MISC%KHALO, &
   & PHYEX%MISC%KSPLIT,PHYEX%MISC%KMI, KSV, KSV_LGBEG, KSV_LGEND, PHYEX%MISC%CPROGRAM,&
   & NSV_LIMA_NR, NSV_LIMA_NS, NSV_LIMA_NG, NSV_LIMA_NH,   &
   & PHYEX%MISC%O2D, PHYEX%MISC%ONOMIXLG, PHYEX%MISC%OFLAT, PHYEX%MISC%OCOUPLES, PHYEX%MISC%OBLOWSNOW,& 
   & PHYEX%MISC%OIBM, PHYEX%MISC%OFLYER, PHYEX%MISC%OCOMPUTE_SRC, PHYEX%MISC%XRSNOW, &
   & PHYEX%MISC%OOCEAN,PHYEX%MISC%ODEEPOC, PHYEX%MISC%ODIAG_IN_RUN,   &
   & PHYEX%MISC%HTURBLEN_CL, CMICRO,           &
   & ZTWOTSTEP,ZTFILE,                                      &
   & ZDXX,ZDYY,ZDZZ,ZDZX,ZDZY,ZZZ,          &
   & ZDIRCOSXW,ZDIRCOSYW,ZDIRCOSZW,ZCOSSLOPE,ZSINSLOPE,    &
   & PRHODJ,PTHVREF,PHGRAD,PZS,                            &
   & PSFTH,PSFRV,PSFSV,PSFU,PSFV,                          &
   & PPABSM,PUM,PVM,PWM,PTKEM,ZSVM,PSRCM,                  &
   & PLENGTHM,PLENGTHH,MFMOIST,                            &
   & ZBL_DEPTH,ZSBL_DEPTH,                                 &
   & ZCEI, PHYEX%MISC%XCEI_MIN, PHYEX%MISC%XCEI_MAX, PHYEX%MISC%XCOEF_AMPL_SAT,    &
   & PTHM,ZRM, &
   & PRUS,PRVS,PRWS,PRTHS,ZRRS,ZRSVS,PRTKES_OUT,         &
   & PSIGS,                                         &
   & PFLXZTHVMF,ZWTH,ZWRC,ZWSV,PDP,PTP,PTDIFF,PTDISS,&
   & YLBUDGET, KBUDGETS=SIZE(YLBUDGET),PEDR=PEDR,PTPMF=PTPMF,&
   & PDRUS_TURB=PDRUS_TURB,PDRVS_TURB=PDRVS_TURB,          &
   & PDRTHLS_TURB=PDRTHLS_TURB,PDRRTS_TURB=PDRRTS_TURB,PDRSVS_TURB=ZDRSVS_TURB)
!
!
!------------------------------------------------------------------------------
!
!
!*       5.   DIVISION PAR RHODJ DES TERMES SOURCES DE MESONH
!             (on obtient des termes homog�nes � des tendances)
!
!         -----------------------------------------------

DO JK=2,KLEV+1
  DO JL = 1,KLON
    PRUS(JL,1,JK)   = PRUS(JL,1,JK)  /PRHODJ(JL,1,JK)
    PRVS(JL,1,JK)   = PRVS(JL,1,JK)  /PRHODJ(JL,1,JK)
    PRTHS(JL,1,JK)  = PRTHS(JL,1,JK) /PRHODJ(JL,1,JK)
    PRTKES_OUT(JL,1,JK)  = PRTKES_OUT(JL,1,JK) /PRHODJ(JL,1,JK)
    PDRUS_TURB(JL,1,JK)  = PDRUS_TURB(JL,1,JK) /PRHODJ(JL,1,JK)
    PDRVS_TURB(JL,1,JK)  = PDRVS_TURB(JL,1,JK) /PRHODJ(JL,1,JK)
    PDRTHLS_TURB(JL,1,JK)  = PDRTHLS_TURB(JL,1,JK) /PRHODJ(JL,1,JK)
    PDRRTS_TURB(JL,1,JK)  = PDRRTS_TURB(JL,1,JK) /PRHODJ(JL,1,JK)
  ENDDO
ENDDO

DO JRR=1,KRR
  DO JK=2,KLEV+1
    DO JL = 1,KLON
      PRRS(JL,1,JK-1,JRR) = ZRRS(JL,1,JK,JRR)/PRHODJ(JL,1,JK)
    ENDDO
    PRM(:,1,JK-1,JRR) = ZRM(:,1,JK,JRR)
  ENDDO
ENDDO

DO JSV=1,KSV
  DO JK=2,KLEV+1
    DO JL = 1,KLON
      PRSVS(JL,1,JK-1,JSV) = ZRSVS(JL,1,JK,JSV)/PRHODJ(JL,1,JK)
      PDRSVS_TURB(JL,1,JK-1,JSV) = ZDRSVS_TURB(JL,1,JK,JSV)/PRHODJ(JL,1,JK)
    ENDDO
  ENDDO
ENDDO

IF (LHOOK) CALL DR_HOOK('ARO_TURB_MNH',1,ZHOOK_HANDLE)

CONTAINS

SUBROUTINE VERTICAL_EXTEND(PX)

 ! fill extra vetical levels to fit MNH interface

REAL, DIMENSION(KLON,1,KLEV+2),   INTENT(INOUT)   :: PX
! NO DR_HOOK, PLEASE ! Rek
PX(:,1,1     )= PX(:,1,2)
PX(:,1,KLEV+2)= PX(:,1,KLEV+1)
END SUBROUTINE VERTICAL_EXTEND

END SUBROUTINE ARO_TURB_MNH
