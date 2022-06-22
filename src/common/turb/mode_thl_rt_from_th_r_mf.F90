!MNH_LIC Copyright 1994-2022 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
MODULE MODE_THL_RT_FROM_TH_R_MF
IMPLICIT NONE
CONTAINS
      SUBROUTINE THL_RT_FROM_TH_R_MF( D, CST, KRR, KRRL, KRRI,       &
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
USE MODD_DIMPHYEX,        ONLY: DIMPHYEX_t
USE MODD_CST, ONLY : CST_t
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
TYPE(DIMPHYEX_t),       INTENT(IN)   :: D
TYPE(CST_t),            INTENT(IN)   :: CST
INTEGER,                INTENT(IN)   :: KRR           ! number of moist var.
INTEGER,                INTENT(IN)   :: KRRL          ! number of liquid water var.
INTEGER,                INTENT(IN)   :: KRRI          ! number of ice water var.

REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)   :: PTH      ! theta
REAL, DIMENSION(D%NIT,D%NKT,KRR), INTENT(IN) :: PR       ! water species
REAL, DIMENSION(D%NIT,D%NKT), INTENT(IN)   :: PEXN    ! exner function

REAL, DIMENSION(D%NIT,D%NKT), INTENT(OUT)  :: PTHL     ! th_l
REAL, DIMENSION(D%NIT,D%NKT), INTENT(OUT)  :: PRT      ! total non precip. water
!
!-------------------------------------------------------------------------------
!
!       0.2  declaration of local variables
!

!----------------------------------------------------------------------------
REAL, DIMENSION(D%NIT,D%NKT) :: ZCP, ZT
REAL, DIMENSION(D%NIT,D%NKT) :: ZLVOCPEXN, ZLSOCPEXN
INTEGER :: JRR, JI, JK
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!----------------------------------------------------------------------------
!
!
IF (LHOOK) CALL DR_HOOK('THL_RT_FRM_TH_R_MF',0,ZHOOK_HANDLE)
!$mnh_expand_array(JI=D%NIB:D%NIE,JK=1:D%NKT)
!temperature
ZT(D%NIB:D%NIE,:) = PTH(D%NIB:D%NIE,:) * PEXN(D%NIB:D%NIE,:)

!Cp
ZCP(D%NIB:D%NIE,:)=CST%XCPD
IF (KRR > 0) ZCP(D%NIB:D%NIE,:) = ZCP(D%NIB:D%NIE,:) + CST%XCPV * PR(D%NIB:D%NIE,:,1)
!$mnh_end_expand_array(JI=D%NIB:D%NIE,JK=1:D%NKT)
DO JRR = 2,1+KRRL  ! loop on the liquid components
  !$mnh_expand_array(JI=D%NIB:D%NIE,JK=1:D%NKT)
  ZCP(D%NIB:D%NIE,:)  = ZCP(D%NIB:D%NIE,:) + CST%XCL * PR(D%NIB:D%NIE,:,JRR)
  !$mnh_end_expand_array(JI=D%NIB:D%NIE,JK=1:D%NKT)
END DO
DO JRR = 2+KRRL,1+KRRL+KRRI ! loop on the solid components
  !$mnh_expand_array(JI=D%NIB:D%NIE,JK=1:D%NKT)
  ZCP(D%NIB:D%NIE,:)  = ZCP(D%NIB:D%NIE,:)  + CST%XCI * PR(D%NIB:D%NIE,:,JRR)
  !$mnh_end_expand_array(JI=D%NIB:D%NIE,JK=1:D%NKT)
END DO

IF ( KRRL >= 1 ) THEN
  IF ( KRRI >= 1 ) THEN
    !$mnh_expand_array(JI=D%NIB:D%NIE,JK=1:D%NKT)
    !ZLVOCPEXN and ZLSOCPEXN
    ZLVOCPEXN(D%NIB:D%NIE,:)=(CST%XLVTT + (CST%XCPV-CST%XCL) *  (ZT(D%NIB:D%NIE,:)-CST%XTT) ) & 
                            &/ ZCP(D%NIB:D%NIE,:) / PEXN(D%NIB:D%NIE,:)
    ZLSOCPEXN(D%NIB:D%NIE,:)=(CST%XLSTT + (CST%XCPV-CST%XCI) *  (ZT(D%NIB:D%NIE,:)-CST%XTT) ) &
                            &/ ZCP(D%NIB:D%NIE,:) / PEXN(D%NIB:D%NIE,:)
    ! Rnp 
    PRT(D%NIB:D%NIE,:)  = PR(D%NIB:D%NIE,:,1)  + PR(D%NIB:D%NIE,:,2)  + PR(D%NIB:D%NIE,:,4)
    ! Theta_l 
    PTHL(D%NIB:D%NIE,:)  = PTH(D%NIB:D%NIE,:)  - ZLVOCPEXN(D%NIB:D%NIE,:) * PR(D%NIB:D%NIE,:,2) &
                           - ZLSOCPEXN(D%NIB:D%NIE,:) * PR(D%NIB:D%NIE,:,4)
    !$mnh_end_expand_array(JI=D%NIB:D%NIE,JK=1:D%NKT)
  ELSE
    !$mnh_expand_array(JI=D%NIB:D%NIE,JK=1:D%NKT)
    !ZLVOCPEXN
    ZLVOCPEXN(D%NIB:D%NIE,:)=(CST%XLVTT + (CST%XCPV-CST%XCL) *  (ZT(D%NIB:D%NIE,:)-CST%XTT) ) &
                            &/ ZCP(D%NIB:D%NIE,:) / PEXN(D%NIB:D%NIE,:)
    ! Rnp
    PRT(D%NIB:D%NIE,:)  = PR(D%NIB:D%NIE,:,1)  + PR(D%NIB:D%NIE,:,2) 
    ! Theta_l
    PTHL(D%NIB:D%NIE,:) = PTH(D%NIB:D%NIE,:)  - ZLVOCPEXN(D%NIB:D%NIE,:) * PR(D%NIB:D%NIE,:,2)
    !$mnh_end_expand_array(JI=D%NIB:D%NIE,JK=1:D%NKT)
  END IF
ELSE
  !$mnh_expand_array(JI=D%NIB:D%NIE,JK=1:D%NKT)
  ! Rnp = rv
  PRT(D%NIB:D%NIE,:)  = PR(D%NIB:D%NIE,:,1)
  ! Theta_l = Theta
  PTHL(D%NIB:D%NIE,:) = PTH(D%NIB:D%NIE,:)
  !$mnh_end_expand_array(JI=D%NIB:D%NIE,JK=1:D%NKT)
END IF
IF (LHOOK) CALL DR_HOOK('THL_RT_FRM_TH_R_MF',1,ZHOOK_HANDLE)
END SUBROUTINE THL_RT_FROM_TH_R_MF
END MODULE MODE_THL_RT_FROM_TH_R_MF
