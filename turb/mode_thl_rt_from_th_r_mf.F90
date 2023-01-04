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

REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   :: PTH      ! theta
REAL, DIMENSION(D%NIJT,D%NKT,KRR), INTENT(IN) :: PR       ! water species
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   :: PEXN    ! exner function

REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)  :: PTHL     ! th_l
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)  :: PRT      ! total non precip. water
!
!-------------------------------------------------------------------------------
!
!       0.2  declaration of local variables
!

!----------------------------------------------------------------------------
REAL, DIMENSION(D%NIJT,D%NKT) :: ZCP, ZT
REAL, DIMENSION(D%NIJT,D%NKT) :: ZLVOCPEXN, ZLSOCPEXN
INTEGER :: JRR, JIJ, JK
INTEGER :: IIJB,IIJE ! physical horizontal domain indices
INTEGER :: IKT
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!----------------------------------------------------------------------------
!
!
IF (LHOOK) CALL DR_HOOK('THL_RT_FRM_TH_R_MF',0,ZHOOK_HANDLE)
!
IIJE=D%NIJE
IIJB=D%NIJB
IKT=D%NKT
!
DO JK=1,IKT 
  DO JIJ=IIJB,IIJE 
!temperature
    ZT(JIJ,JK) = PTH(JIJ,JK) * PEXN(JIJ,JK)
    
!Cp
    ZCP(JIJ,JK)=CST%XCPD
    IF (KRR > 0) ZCP(JIJ,JK) = ZCP(JIJ,JK) + CST%XCPV * PR(JIJ,JK,1)
  ENDDO
ENDDO
DO JRR = 2,1+KRRL  ! loop on the liquid components
  DO JK=1,IKT 
    DO JIJ=IIJB,IIJE 
      ZCP(JIJ,JK)  = ZCP(JIJ,JK) + CST%XCL * PR(JIJ,JK,JRR)
    ENDDO
  ENDDO
END DO
DO JRR = 2+KRRL,1+KRRL+KRRI ! loop on the solid components
  DO JK=1,IKT 
    DO JIJ=IIJB,IIJE 
      ZCP(JIJ,JK)  = ZCP(JIJ,JK)  + CST%XCI * PR(JIJ,JK,JRR)
    ENDDO
  ENDDO
END DO

IF ( KRRL >= 1 ) THEN
  IF ( KRRI >= 1 ) THEN
    DO JK=1,IKT 
      DO JIJ=IIJB,IIJE 
    !ZLVOCPEXN and ZLSOCPEXN
        ZLVOCPEXN(JIJ,JK)=(CST%XLVTT + (CST%XCPV-CST%XCL) *  (ZT(JIJ,JK)-CST%XTT) ) & 
        &/ ZCP(JIJ,JK) / PEXN(JIJ,JK)
        ZLSOCPEXN(JIJ,JK)=(CST%XLSTT + (CST%XCPV-CST%XCI) *  (ZT(JIJ,JK)-CST%XTT) ) &
        &/ ZCP(JIJ,JK) / PEXN(JIJ,JK)
    ! Rnp 
        PRT(JIJ,JK)  = PR(JIJ,JK,1)  + PR(JIJ,JK,2)  + PR(JIJ,JK,4)
    ! Theta_l 
        PTHL(JIJ,JK)  = PTH(JIJ,JK)  - ZLVOCPEXN(JIJ,JK) * PR(JIJ,JK,2) &
        - ZLSOCPEXN(JIJ,JK) * PR(JIJ,JK,4)
      ENDDO
    ENDDO
  ELSE
    DO JK=1,IKT 
      DO JIJ=IIJB,IIJE 
    !ZLVOCPEXN
        ZLVOCPEXN(JIJ,JK)=(CST%XLVTT + (CST%XCPV-CST%XCL) *  (ZT(JIJ,JK)-CST%XTT) ) &
        &/ ZCP(JIJ,JK) / PEXN(JIJ,JK)
    ! Rnp
        PRT(JIJ,JK)  = PR(JIJ,JK,1)  + PR(JIJ,JK,2) 
    ! Theta_l
        PTHL(JIJ,JK) = PTH(JIJ,JK)  - ZLVOCPEXN(JIJ,JK) * PR(JIJ,JK,2)
      ENDDO
    ENDDO
  END IF
ELSE
  DO JK=1,IKT 
    DO JIJ=IIJB,IIJE 
  ! Rnp = rv
      PRT(JIJ,JK)  = PR(JIJ,JK,1)
  ! Theta_l = Theta
      PTHL(JIJ,JK) = PTH(JIJ,JK)
    ENDDO
  ENDDO
END IF
IF (LHOOK) CALL DR_HOOK('THL_RT_FRM_TH_R_MF',1,ZHOOK_HANDLE)
END SUBROUTINE THL_RT_FROM_TH_R_MF
END MODULE MODE_THL_RT_FROM_TH_R_MF
