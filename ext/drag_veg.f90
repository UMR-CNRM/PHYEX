!MNH_LIC Copyright 2009-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     #######################
       MODULE MODI_DRAG_VEG
!     #######################
!
INTERFACE

SUBROUTINE DRAG_VEG(PTSTEP,PUT,PVT,PTKET,ODEPOTREE, PVDEPOTREE, &
                    HCLOUD,PPABST,PTHT,PRT,PSVT,         &
                    PRHODJ,PZZ,PRUS, PRVS, PRTKES,       &
                    PRRS,PSVS)
!
REAL,                     INTENT(IN)    :: PTSTEP ! Time step
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PUT, PVT   ! variables
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTKET           !   at t
LOGICAL,                  INTENT(IN)    :: ODEPOTREE ! Droplet deposition on tree
REAL,                     INTENT(IN)    :: PVDEPOTREE! Velocity deposition on tree
CHARACTER (LEN=4),        INTENT(IN)    :: HCLOUD       ! Kind of microphysical scheme
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABST          !   at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT          !   at t
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRT             !   at t
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PSVT            !   at t
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ    ! dry Density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PZZ       ! Height (z)
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRUS, PRVS       ! Sources of Momentum
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRTKES           ! Sources of Tke
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRRS         
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PSVS       
!
END SUBROUTINE DRAG_VEG

END INTERFACE

END MODULE MODI_DRAG_VEG
!
!     ###################################################################
SUBROUTINE DRAG_VEG(PTSTEP,PUT,PVT,PTKET,ODEPOTREE, PVDEPOTREE, &
                    HCLOUD,PPABST,PTHT,PRT,PSVT,         &
                    PRHODJ,PZZ,PRUS, PRVS, PRTKES,       &
                    PRRS,PSVS)
!     ###################################################################
!
!!****  *DRAG_VEG_n * -
!!
!!    PURPOSE
!!    -------
!
!!**  METHOD
!!    ------
!!
!!    REFERENCE
!!    ---------
!!      
!!
!!    AUTHOR
!!    ------
!!     P. Aumond 
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    07/2009
!!       C.Lac      07/2011 : Add budgets
!!       S. Donier  06/2015 : bug surface aerosols
!!       C.Lac      07/2016 : Add droplet deposition
!!       C.Lac      10/2017 : Correction on deposition
!  C. Lac         11/2019: correction in the drag formula and application to building in addition to tree
!  P. Wautelet 28/01/2020: use the new data structures and subroutines for budgets for U
!  C. Lac         02/2020: correction missing condition for budget on RC and SV
!  P. Wautelet 04/02/2021: budgets: bugfixes for LDRAGTREE if LIMA + small optimisations and verifications
!  R. Schoetter 04/2022: bug add update halo for vegetation drag variables
!!---------------------------------------------------------------
!
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_ARGSLIST_ll, ONLY: LIST_ll
use modd_budget,      only: lbudget_u, lbudget_v, lbudget_rc, lbudget_sv,  lbudget_tke, &
                            NBUDGET_U, NBUDGET_V, NBUDGET_RC, NBUDGET_SV1, NBUDGET_TKE, &
                            tbudgets
USE MODD_CONF
USE MODD_CST
USE MODD_DYN
USE MODD_DYN_n
USE MODD_GROUND_PAR
USE MODD_NSV
USE MODD_PARAM_C2R2
USE MODD_PARAM_LIMA,  ONLY: NMOM_C
USE MODD_PARAM_n,     only: CSURF, CTURB
USE MODD_PGDFIELDS
USE MODD_VEG_n

use mode_budget,     only: Budget_store_init, Budget_store_end
use mode_msg
USE MODE_ll

USE MODI_MNHGET_SURF_PARAM_n
USE MODI_SHUMAN

IMPLICIT NONE
!  
!*       0.1   Declarations of dummy arguments :
!
REAL,                     INTENT(IN)    :: PTSTEP ! Time step
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PUT, PVT   ! variables
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTKET           !   at t
LOGICAL,                  INTENT(IN)    :: ODEPOTREE ! Droplet deposition on tree
REAL,                     INTENT(IN)    :: PVDEPOTREE! Velocity deposition on tree
CHARACTER (LEN=4),        INTENT(IN)    :: HCLOUD       ! Kind of microphysical scheme
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABST          !   at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT          !   at t
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRT             !   at t
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PSVT            !   at t
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ    ! dry Density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PZZ       ! Height (z)
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRUS, PRVS       ! Sources of Momentum
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRTKES           ! Sources of Tke
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRRS         
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PSVS       
!
!*       0.2   Declarations of local variables :
!
INTEGER ::  IIU,IJU,IKU,IKV         ! array size along the k direction
INTEGER :: JI, JJ, JK             ! loop index
INTEGER :: IINFO_ll
TYPE(LIST_ll), POINTER :: TZFIELDS_ll   ! list of fields to exchange
!
!
REAL, DIMENSION(SIZE(PUT,1),SIZE(PUT,2),SIZE(PUT,3)) ::           &
                              ZWORK1, ZWORK2, ZWORK3, ZUT_SCAL, ZVT_SCAL,   &
                              ZUS, ZVS, ZTKES, ZTKET
REAL, DIMENSION(SIZE(PUT,1),SIZE(PUT,2),SIZE(PUT,3)) ::           &
                              ZCDRAG, ZDENSITY
REAL, DIMENSION(SIZE(PUT,1),SIZE(PUT,2)) ::           &
                              ZH,ZLAI           !  LAI, Vegetation height
REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)):: ZT,ZEXN,ZLV,ZCPH                              
LOGICAL, DIMENSION(SIZE(PUT,1),SIZE(PUT,2),SIZE(PUT,3)) &
            :: GDEP
REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)):: ZWDEPR,ZWDEPS

IF ( CSURF /= 'EXTE' ) CALL PRINT_MSG( NVERB_FATAL, 'GEN', 'DRAG_VEG', 'CSURF/=EXTE not allowed' )

!Condition necessary because PTKET is used (and must be allocated)
IF ( CTURB /= 'TKEL' ) CALL PRINT_MSG( NVERB_FATAL, 'GEN', 'DRAG_VEG', 'CTURB/=TKEL not allowed' )
!
if ( lbudget_u   ) call Budget_store_init( tbudgets(NBUDGET_U  ), 'DRAG', prus  (:, :, :) )
if ( lbudget_v   ) call Budget_store_init( tbudgets(NBUDGET_V  ), 'DRAG', prvs  (:, :, :) )
if ( lbudget_tke ) call Budget_store_init( tbudgets(NBUDGET_TKE), 'DRAG', prtkes(:, :, :) )

if ( odepotree ) then
  if ( lbudget_rc ) call Budget_store_init( tbudgets(NBUDGET_RC), 'DEPOTR', prrs(:, :, :, 2) )
  if ( lbudget_sv .and. ( hcloud=='C2R2' .or. hcloud=='KHKO' ) ) &
                    call Budget_store_init( tbudgets(NBUDGET_SV1-1+(NSV_C2R2BEG+1)), 'DEPOTR', psvs(:, :, :, NSV_C2R2BEG+1) )
  if ( lbudget_sv .and.   hcloud=='LIMA' ) &
                    call Budget_store_init( tbudgets(NBUDGET_SV1-1+NSV_LIMA_NC),     'DEPOTR', psvs(:, :, :, NSV_LIMA_NC) )
end if

IIU = SIZE(PUT,1)
IJU = SIZE(PUT,2)
IKU = SIZE(PUT,3)
!
ZUS   (:,:,:) = 0.0
ZVS   (:,:,:) = 0.0
ZTKES (:,:,:) = 0.0
!
ZH  (:,:) = XUNDEF
ZLAI(:,:) = XUNDEF
!
ZCDRAG   (:,:,:) = 0.
ZDENSITY (:,:,:) = 0.
!
CALL MNHGET_SURF_PARAM_n( PH_TREE = ZH, PLAI_TREE = ZLAI )
!
WHERE ( ZH   (:,:) > (XUNDEF-1.) ) ZH   (:,:) = 0.0
WHERE ( ZLAI (:,:) > (XUNDEF-1.) ) ZLAI (:,:) = 0.0
!
!-------------------------------------------------------------------------------
!
!
!*       1.     COMPUTES THE TRUE VELOCITY COMPONENTS
!	        -------------------------------------
!
ZUT_SCAL(:,:,:) = MXF(PUT(:,:,:))
ZVT_SCAL(:,:,:) = MYF(PVT(:,:,:))
ZTKET(:,:,:)    = PTKET(:,:,:)
!
! Update halo
!
NULLIFY(TZFIELDS_ll)
CALL ADD3DFIELD_ll( TZFIELDS_ll, ZUT_SCAL, 'DRAG_VEG::ZUT_SCAL')
CALL ADD3DFIELD_ll( TZFIELDS_ll, ZVT_SCAL, 'DRAG_VEG::ZVT_SCAL')
CALL ADD3DFIELD_ll( TZFIELDS_ll, ZTKET   , 'DRAG_VEG::ZTKET'   )
CALL UPDATE_HALO_ll(TZFIELDS_ll,IINFO_ll)
CALL CLEANLIST_ll(TZFIELDS_ll)
!
!-------------------------------------------------------------------------------
!
!*      1.     Computations of wind tendency due to canopy drag
!              ------------------------------------------------
!
!
!
! Ext = - Cdrag  * u- * u- * Sv       tree canopy drag
!       - u'w'(ground)     * Sh       horizontal surfaces (ground)
!
!*      1.1    Drag coefficient by vegetation (Patton et al 2001)
!              ------------------------------
!
GDEP(:,:,:) = .FALSE.
!
DO JJ=2,(IJU-1)
   DO JI=2,(IIU-1)
      !
      ! Set density and drag coefficient for vegetation
      !
      IF (ZH(JI,JJ) /= 0) THEN
         !
         DO JK=2,(IKU-1)
            !
            IF ( (PZZ(JI,JJ,JK)-PZZ(JI,JJ,2)) .LT. ZH(JI,JJ) ) THEN
               !
               IF ((HCLOUD=='C2R2') .OR.  (HCLOUD=='KHKO')) THEN
                  IF ((PRRS(JI,JJ,JK,2) >0.) .AND. (PSVS(JI,JJ,JK,NSV_C2R2BEG+1) >0.)) &
                       GDEP(JI,JJ,JK) = .TRUE.
               ELSE IF (HCLOUD /= 'NONE' .AND. HCLOUD /= 'REVE') THEN
                  IF (PRRS(JI,JJ,JK,2) >0.) GDEP(JI,JJ,JK) = .TRUE.
               ENDIF
               !
               ZCDRAG(JI,JJ,JK)  = 0.2 !0.075
               ZDENSITY(JI,JJ,JK) = MAX((4 * (ZLAI(JI,JJ) *&
                    (PZZ(JI,JJ,JK)-PZZ(JI,JJ,2)) *&
                    (PZZ(JI,JJ,JK)-PZZ(JI,JJ,2)) *&
                    (ZH(JI,JJ)-(PZZ(JI,JJ,JK)-PZZ(JI,JJ,2)))/&
                    ZH(JI,JJ)**3)-&
                    (0.30*((ZLAI(JI,JJ) *&
                    (PZZ(JI,JJ,JK)-PZZ(JI,JJ,2)) *&
                    (PZZ(JI,JJ,JK)-PZZ(JI,JJ,2)) *&
                    (PZZ(JI,JJ,JK)-PZZ(JI,JJ,2)) /&
                    (ZH(JI,JJ)**3))-ZLAI(JI,JJ))))/&
                    ZH(JI,JJ), 0.)
               !
            ENDIF
            !
         ENDDO
      ENDIF
      !
   ENDDO
ENDDO
!
! To exclude the first vertical level already dealt in rain_ice or rain_c2r2_khko
GDEP(:,:,2) = .FALSE.
!
! Update halo
!
NULLIFY(TZFIELDS_ll)
CALL ADD3DFIELD_ll( TZFIELDS_ll, ZCDRAG  , 'DRAG_VEG::ZCDRAG')
CALL ADD3DFIELD_ll( TZFIELDS_ll, ZDENSITY, 'DRAG_VEG::ZDENSITY')
CALL UPDATE_HALO_ll(TZFIELDS_ll,IINFO_ll)
CALL CLEANLIST_ll(TZFIELDS_ll)
!
!
!*      1.2    Drag force by wall surfaces
!              ---------------------------
!
!* drag force by vertical surfaces
!
ZUS(:,:,:) = PUT(:,:,:)/( 1.0 + MXM ( ZCDRAG(:,:,:) * ZDENSITY(:,:,:) &
     * PTSTEP * SQRT(ZUT_SCAL(:,:,:)**2+ZVT_SCAL(:,:,:)**2) ) )
!
ZVS(:,:,:) = PVT(:,:,:)/( 1.0 + MYM ( ZCDRAG(:,:,:) * ZDENSITY(:,:,:) &
     * PTSTEP * SQRT(ZUT_SCAL(:,:,:)**2+ZVT_SCAL(:,:,:)**2) ) )
!
PRUS(:,:,:) = PRUS(:,:,:) + (ZUS(:,:,:)-PUT(:,:,:)) * MXM(PRHODJ(:,:,:)) / PTSTEP
!
PRVS(:,:,:) = PRVS(:,:,:) + (ZVS(:,:,:)-PVT(:,:,:)) * MYM(PRHODJ(:,:,:)) / PTSTEP
!
IF (ODEPOTREE) THEN
  IF ( HCLOUD == 'NONE' ) CALL PRINT_MSG( NVERB_FATAL, 'GEN', 'DRAG_VEG', 'LDEPOTREE=T not allowed if CCLOUD=NONE' )
  IF ( HCLOUD == 'LIMA' .AND. NMOM_C.EQ.0 ) &
    CALL PRINT_MSG( NVERB_FATAL, 'GEN', 'DRAG_VEG', 'LDEPOTREE=T not allowed if CCLOUD=LIMA and NMOM_C=0' )

  ZEXN(:,:,:)= (PPABST(:,:,:)/XP00)**(XRD/XCPD)
  ZT(:,:,:)= PTHT(:,:,:)*ZEXN(:,:,:)
  ZLV(:,:,:)=XLVTT +(XCPV-XCL) *(ZT(:,:,:)-XTT)
  ZCPH(:,:,:)=XCPD +XCPV*PRT(:,:,:,1)
  ZWDEPR(:,:,:)= 0.
  ZWDEPS(:,:,:)= 0.
  WHERE (GDEP)
   ZWDEPR(:,:,:)= PVDEPOTREE * PRT(:,:,:,2) * PRHODJ(:,:,:)
  END WHERE
  IF ( HCLOUD == 'C2R2' .OR. HCLOUD == 'KHKO' ) THEN
    WHERE (GDEP)
      ZWDEPS(:,:,:)= PVDEPOTREE * PSVT(:,:,:,NSV_C2R2BEG+1) * PRHODJ(:,:,:)
    END WHERE
  ELSE IF ( HCLOUD == 'LIMA' ) THEN
    WHERE (GDEP)
      ZWDEPS(:,:,:)= PVDEPOTREE * PSVT(:,:,:,NSV_LIMA_NC) * PRHODJ(:,:,:)
    END WHERE
  END IF
  DO JJ=2,(IJU-1)
   DO JI=2,(IIU-1)
     DO JK=2,(IKU-2) 
       IF (GDEP(JI,JJ,JK)) THEN
         PRRS(JI,JJ,JK,2) = PRRS(JI,JJ,JK,2) + (ZWDEPR(JI,JJ,JK+1)-ZWDEPR(JI,JJ,JK))/ &
                            (PZZ(JI,JJ,JK+1)-PZZ(JI,JJ,JK))
         IF ( HCLOUD == 'C2R2' .OR. HCLOUD == 'KHKO' ) THEN
           PSVS(JI,JJ,JK,NSV_C2R2BEG+1) = PSVS(JI,JJ,JK,NSV_C2R2BEG+1) + &
                                          (ZWDEPS(JI,JJ,JK+1)-ZWDEPS(JI,JJ,JK))/(PZZ(JI,JJ,JK+1)-PZZ(JI,JJ,JK))
         ELSE IF ( HCLOUD == 'LIMA' ) THEN
           PSVS(JI,JJ,JK,NSV_LIMA_NC) = PSVS(JI,JJ,JK,NSV_LIMA_NC) + &
                                        (ZWDEPS(JI,JJ,JK+1)-ZWDEPS(JI,JJ,JK))/(PZZ(JI,JJ,JK+1)-PZZ(JI,JJ,JK))
         END IF
       END IF
     END DO
    END DO
   END DO
!
!
END IF
!
!*      3.     Computations of TKE  tendency due to canopy drag
!              ------------------------------------------------

!*      3.1    Creation of TKE by wake
!              -----------------------
!
! from Kanda and Hino (1994)
!
! Ext = + Cd * u+^3  * Sv/Vair        vertical surfaces or trees             
! Ext = - Cd * e * u  * Sv        trees Destruction of TKE due to 
!   small-scale motions forced by leaves from Kanda and Hino (1994)
!
! with Vair = Vair/Vtot * Vtot = (Vair/Vtot) * Stot * Dz
! and  Sv/Vair = (Sv/Stot) * Stot/Vair = (Sv/Stot) / (Vair/Vtot) / Dz
!
ZTKES(:,:,:)=  ( ZTKET(:,:,:) + PTSTEP * ZCDRAG(:,:,:) * ZDENSITY(:,:,:) &
         * (SQRT( ZUT_SCAL(:,:,:)**2 + ZVT_SCAL(:,:,:)**2 ))**3 ) /      &
     ( 1. + PTSTEP * ZCDRAG(:,:,:) * ZDENSITY(:,:,:) * SQRT(ZUT_SCAL(:,:,:)**2+ZVT_SCAL(:,:,:)**2))
!
PRTKES(:,:,:) = PRTKES(:,:,:) + (ZTKES(:,:,:)-ZTKET(:,:,:))*PRHODJ(:,:,:)/PTSTEP

if ( lbudget_u   ) call Budget_store_end( tbudgets(NBUDGET_U  ), 'DRAG', prus  (:, :, :) )
if ( lbudget_v   ) call Budget_store_end( tbudgets(NBUDGET_V  ), 'DRAG', prvs  (:, :, :) )
if ( lbudget_tke ) call Budget_store_end( tbudgets(NBUDGET_TKE), 'DRAG', prtkes(:, :, :) )

if ( odepotree ) then
  if ( lbudget_rc ) call Budget_store_end( tbudgets(NBUDGET_RC), 'DEPOTR', prrs(:, :, :, 2) )
  if ( lbudget_sv .and. ( hcloud=='C2R2' .or. hcloud=='KHKO' ) ) &
                    call Budget_store_end( tbudgets(NBUDGET_SV1-1+(NSV_C2R2BEG+1)), 'DEPOTR', psvs(:, :, :, NSV_C2R2BEG+1) )
  if ( lbudget_sv .and.   hcloud=='LIMA' ) &
                    call Budget_store_end( tbudgets(NBUDGET_SV1-1+NSV_LIMA_NC),     'DEPOTR', psvs(:, :, :, NSV_LIMA_NC) )
end if

END SUBROUTINE DRAG_VEG
