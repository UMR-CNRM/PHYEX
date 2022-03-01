!MNH_LIC Copyright 1994-2022 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
MODULE MODE_THL_RT_FROM_TH_R_MF
IMPLICIT NONE
CONTAINS
      SUBROUTINE THL_RT_FROM_TH_R_MF( KRR,KRRL,KRRI,                  &
                                      PTH, PR, PEXN, &
                                      PTHL, PRT                      )
!     #################################################################
!
!!
!!****  *THL_RT_FROM_TH_R* - computes the conservative variables THL and RT
!!                           from TH and the non precipitating water species
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!    
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      V. Masson               * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original         20/09/02
!!      Externalisation of computations done in TURB and MF_TURB (Malardel and Pergaud, fev. 2007)
!!      V.Masson : Optimization
!!      S. Riette 2011 suppression of PLVOCPEXN and PLSOCPEXN
!!
!! --------------------------------------------------------------------------
!       
!*      0. DECLARATIONS
!          ------------
!
USE MODD_CST
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
INTEGER,                INTENT(IN)   :: KRR           ! number of moist var.
INTEGER,                INTENT(IN)   :: KRRL          ! number of liquid water var.
INTEGER,                INTENT(IN)   :: KRRI          ! number of ice water var.

REAL, DIMENSION(:,:), INTENT(IN)   :: PTH      ! theta
REAL, DIMENSION(:,:,:), INTENT(IN) :: PR       ! water species
REAL, DIMENSION(:,:), INTENT(IN)   :: PEXN    ! exner function

REAL, DIMENSION(:,:), INTENT(OUT)  :: PTHL     ! th_l
REAL, DIMENSION(:,:), INTENT(OUT)  :: PRT      ! total non precip. water
!
!-------------------------------------------------------------------------------
!
!       0.2  declaration of local variables
!

!----------------------------------------------------------------------------
REAL, DIMENSION(SIZE(PTH,1),SIZE(PTH,2)) :: ZCP, ZT
REAL, DIMENSION(SIZE(PTH,1),SIZE(PTH,2)) :: ZLVOCPEXN, ZLSOCPEXN
INTEGER :: JRR
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!----------------------------------------------------------------------------
!
!
IF (LHOOK) CALL DR_HOOK('THL_RT_FRM_TH_R_MF',0,ZHOOK_HANDLE)
!temperature
ZT(:,:) = PTH(:,:) * PEXN(:,:)

!Cp
ZCP=XCPD
IF (KRR > 0) ZCP(:,:) = ZCP(:,:) + XCPV * PR(:,:,1)
DO JRR = 2,1+KRRL  ! loop on the liquid components
  ZCP(:,:)  = ZCP(:,:) + XCL * PR(:,:,JRR)
END DO
DO JRR = 2+KRRL,1+KRRL+KRRI ! loop on the solid components
  ZCP(:,:)  = ZCP(:,:)  + XCI * PR(:,:,JRR)
END DO

IF ( KRRL >= 1 ) THEN
  IF ( KRRI >= 1 ) THEN
    !ZLVOCPEXN and ZLSOCPEXN
    ZLVOCPEXN(:,:)=(XLVTT + (XCPV-XCL) *  (ZT(:,:)-XTT) ) / ZCP(:,:) / PEXN(:,:)
    ZLSOCPEXN(:,:)=(XLSTT + (XCPV-XCI) *  (ZT(:,:)-XTT) ) / ZCP(:,:) / PEXN(:,:)
    ! Rnp 
    PRT(:,:)  = PR(:,:,1)  + PR(:,:,2)  + PR(:,:,4)
    ! Theta_l 
    PTHL(:,:)  = PTH(:,:)  - ZLVOCPEXN(:,:) * PR(:,:,2) &
                           - ZLSOCPEXN(:,:) * PR(:,:,4)
  ELSE
    !ZLVOCPEXN
    ZLVOCPEXN(:,:)=(XLVTT + (XCPV-XCL) *  (ZT(:,:)-XTT) ) / ZCP(:,:) / PEXN(:,:)
    ! Rnp
    PRT(:,:)  = PR(:,:,1)  + PR(:,:,2) 
    ! Theta_l
    PTHL(:,:) = PTH(:,:)  - ZLVOCPEXN(:,:) * PR(:,:,2)
  END IF
ELSE
  ! Rnp = rv
  PRT(:,:)  = PR(:,:,1)
  ! Theta_l = Theta
  PTHL(:,:) = PTH(:,:)
END IF
IF (LHOOK) CALL DR_HOOK('THL_RT_FRM_TH_R_MF',1,ZHOOK_HANDLE)
END SUBROUTINE THL_RT_FROM_TH_R_MF
END MODULE MODE_THL_RT_FROM_TH_R_MF
