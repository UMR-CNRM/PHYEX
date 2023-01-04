!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!     ######spl
     MODULE MODE_COMPUTE_FUNCTION_THERMO_MF
!    ######################################
!
IMPLICIT NONE
CONTAINS
      SUBROUTINE COMPUTE_FUNCTION_THERMO_MF(D, CST, KRR,KRRL,KRRI,OSTATNW,   &
                                       PTH, PR, PEXN, PFRAC_ICE, PPABS,      &
                                       PT,PAMOIST,PATHETA                    )
!     #################################################################
!
!!
!!****  *COMPUTE_FUNCTION_THERMO_MF* -
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
!!
!!     JP Pinty      *LA*
!!
!!    MODIFICATIONS
!!    -------------
!!     Original   24/02/03
!!     Externalisation of computations done in TURB and MF_TURB (Malardel and Pergaud, fev. 2007)
!!     Optimization : V.Masson, 09/2010
!!     S. Riette Sept 2011 : remove of unused PL?OCPEXN, use of received ice fraction
!!     Wim de Rooy June 2019: update statistical cloud scheme
!!
!! --------------------------------------------------------------------------
!
!*      0. DECLARATIONS
!          ------------
!
USE MODD_DIMPHYEX,        ONLY: DIMPHYEX_t
USE MODD_CST,             ONLY: CST_t
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
LOGICAL,                INTENT(IN)   :: OSTATNW      ! cloud scheme inclues convect. covar. contrib
INTEGER,                INTENT(IN)   :: KRR           ! number of moist var.
INTEGER,                INTENT(IN)   :: KRRL          ! number of liquid water var.
INTEGER,                INTENT(IN)   :: KRRI          ! number of ice water var.

REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   :: PTH      ! theta
REAL, DIMENSION(D%NIJT,D%NKT,KRR), INTENT(IN) :: PR       ! water species
REAL, DIMENSION(D%NIJT,D%NKT)  , INTENT(IN) :: PPABS,PEXN    ! pressure, Exner funct.
REAL, DIMENSION(D%NIJT,D%NKT)  , INTENT(IN) :: PFRAC_ICE     ! ice fraction

REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)   :: PT      ! temperature

REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)  ::  PAMOIST,PATHETA
!
!-------------------------------------------------------------------------------
!
!*       0.2   Declarations of local variables
!
REAL                :: ZEPS         ! XMV / XMD
REAL, DIMENSION(D%NIJT,D%NKT) ::      &
          ZCP,                        &  ! Cp
          ZE,                         &  ! Saturation mixing ratio
          ZDEDT,                      &  ! Saturation mixing ratio derivative
          ZAMOIST_W,                  &  ! Coefficients for s = f (Thetal,Rnp)
          ZATHETA_W,                  &  !
          ZAMOIST_I,                  &  !
          ZATHETA_I,                  &  !
          ZLVOCP,ZLSOCP

INTEGER             :: JRR, JIJ, JK
INTEGER :: IIJB,IIJE ! physical horizontal domain indices
INTEGER :: IKTB,IKTE
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('COMPUTE_FUNCTION_THERMO_MF',0,ZHOOK_HANDLE)
!
IIJE=D%NIJE
IIJB=D%NIJB
IKTB=D%NKTB
IKTE=D%NKTE
!
  ZEPS = CST%XMV / CST%XMD

!
!*       Cph
!
ZCP=CST%XCPD

IF (KRR > 0) THEN
  DO JK=IKTB,IKTE 
    DO JIJ=IIJB,IIJE 
      ZCP(JIJ,JK) = ZCP(JIJ,JK) + CST%XCPV * PR(JIJ,JK,1)
    ENDDO
  ENDDO
ENDIF

DO JRR = 2,1+KRRL  ! loop on the liquid components
  DO JK=IKTB,IKTE 
    DO JIJ=IIJB,IIJE 
      ZCP(JIJ,JK)  = ZCP(JIJ,JK) + CST%XCL * PR(JIJ,JK,JRR)
    ENDDO
  ENDDO
END DO

DO JRR = 2+KRRL,1+KRRL+KRRI ! loop on the solid components
  DO JK=IKTB,IKTE 
    DO JIJ=IIJB,IIJE 
      ZCP(JIJ,JK)  = ZCP(JIJ,JK)  + CST%XCI * PR(JIJ,JK,JRR)
    ENDDO
  ENDDO

END DO

!*      Temperature
!
DO JK=IKTB,IKTE 
  DO JIJ=IIJB,IIJE 
    PT(JIJ,JK) =  PTH(JIJ,JK) * PEXN(JIJ,JK)
  ENDDO
ENDDO
!
!
!! Liquid water
!
IF ( KRRL >= 1 ) THEN
  DO JK=IKTB,IKTE 
    DO JIJ=IIJB,IIJE 
  !
  !*       Lv/Cph
  !
      ZLVOCP(JIJ,JK) = (CST%XLVTT + (CST%XCPV-CST%XCL) *  (PT(JIJ,JK)-CST%XTT) ) / &
      & ZCP(JIJ,JK)
  !
  !*      Saturation vapor pressure with respect to water
  !
      ZE(JIJ,JK) =  EXP(CST%XALPW - CST%XBETAW/PT(JIJ,JK) - &
      &CST%XGAMW*ALOG( PT(JIJ,JK) ) )
  !
  !*      Saturation  mixing ratio with respect to water
  !
      ZE(JIJ,JK) =  ZE(JIJ,JK) * ZEPS / &
      & ( PPABS(JIJ,JK) - ZE(JIJ,JK) )
  !
  !*      Compute the saturation mixing ratio derivative (rvs')
  !
      ZDEDT(JIJ,JK) = (CST%XBETAW/PT(JIJ,JK)  - CST%XGAMW) / PT(JIJ,JK)&
      * ZE(JIJ,JK) * ( 1. + ZE(JIJ,JK) / ZEPS )
  !
  !*      Compute Amoist and Atheta
  !
      IF (OSTATNW) THEN
        ZAMOIST_W(JIJ,JK)=  1.0/( 1.0 + ZDEDT(JIJ,JK) * ZLVOCP(JIJ,JK))
        ZATHETA_W(JIJ,JK)= ZAMOIST_W(JIJ,JK) * PEXN(JIJ,JK) &
        * ZDEDT(JIJ,JK)
      ELSE
        ZAMOIST_W(JIJ,JK)= 0.5/( 1.0 + ZDEDT(JIJ,JK) * ZLVOCP(JIJ,JK) )
        ZATHETA_W(JIJ,JK)= ZAMOIST_W(JIJ,JK) * PEXN(JIJ,JK) *         &
        ( ( ZE(JIJ,JK) - PR(JIJ,JK,1) ) * ZLVOCP(JIJ,JK) /      &
        ( 1. + ZDEDT(JIJ,JK) * ZLVOCP(JIJ,JK) )           *              &
        (                                                             &
        ZE(JIJ,JK) * (1. + ZE(JIJ,JK)/ZEPS)                                &
        * ( -2.*CST%XBETAW/PT(JIJ,JK) + CST%XGAMW ) / PT(JIJ,JK)**2&
        +ZDEDT(JIJ,JK) * (1. + 2. * ZE(JIJ,JK)/ZEPS)                        &
        * ( CST%XBETAW/PT(JIJ,JK) - CST%XGAMW ) / PT(JIJ,JK)         &
        )                                                             &
        - ZDEDT(JIJ,JK)                                                   &
        )
      END IF
    ENDDO
  ENDDO
  !
  !! Solid water
  !
  IF ( KRRI >= 1 ) THEN
    DO JK=IKTB,IKTE 
      DO JIJ=IIJB,IIJE 
    !
    !*       Ls/Cph
    !
        ZLSOCP(JIJ,JK) = (CST%XLSTT + (CST%XCPV-CST%XCI) *  (PT(JIJ,JK)-CST%XTT) ) / &
        & ZCP(JIJ,JK)
    !
    !*      Saturation vapor pressure with respect to ice
    !
        ZE(JIJ,JK) =  EXP(CST%XALPI - CST%XBETAI/PT(JIJ,JK) - &
        &CST%XGAMI*ALOG( PT(JIJ,JK) ) )
    !
    !*      Saturation  mixing ratio with respect to ice
    !
        ZE(JIJ,JK) =  ZE(JIJ,JK) * ZEPS / &
        & ( PPABS(JIJ,JK) - ZE(JIJ,JK) )
    !
    !*      Compute the saturation mixing ratio derivative (rvs')
    !
        ZDEDT(JIJ,JK) = (CST%XBETAI/PT(JIJ,JK)-CST%XGAMI) /PT(JIJ,JK)&
        * ZE(JIJ,JK) * ( 1. + ZE(JIJ,JK) / ZEPS )
    !
    !*      Compute Amoist and Atheta
    !
        IF (OSTATNW) THEN
          ZAMOIST_I(JIJ,JK)= 1.0/( 1.0 + ZDEDT(JIJ,JK) *ZLVOCP(JIJ,JK))
          ZATHETA_I(JIJ,JK)= ZAMOIST_I(JIJ,JK) * PEXN(JIJ,JK) &
          * ZDEDT(JIJ,JK)
        ELSE
          ZAMOIST_I(JIJ,JK)= 0.5/(1.0 + ZDEDT(JIJ,JK) * ZLSOCP(JIJ,JK))
          ZATHETA_I(JIJ,JK)= ZAMOIST_I(JIJ,JK) * PEXN(JIJ,JK) *      &
          ( ( ZE(JIJ,JK) - PR(JIJ,JK,1) ) * ZLSOCP(JIJ,JK) /       &
          ( 1. + ZDEDT(JIJ,JK) * ZLSOCP(JIJ,JK) )           *              &
          (                                                             &
          ZE(JIJ,JK) * (1. + ZE(JIJ,JK)/ZEPS)                                &
          * ( -2.*CST%XBETAI/PT(JIJ,JK) + CST%XGAMI ) / PT(JIJ,JK)**2   &
          +ZDEDT(JIJ,JK) * (1. + 2. * ZE(JIJ,JK)/ZEPS)                        &
          * ( CST%XBETAI/PT(JIJ,JK) - CST%XGAMI ) / PT(JIJ,JK)          &
          )                                                             &
          - ZDEDT(JIJ,JK)                                                   &
          )
        END IF
      ENDDO
    ENDDO

  ELSE
    ZAMOIST_I(IIJB:IIJE,IKTB:IKTE)=0.
    ZATHETA_I(IIJB:IIJE,IKTB:IKTE)=0.
  ENDIF

  DO JK=IKTB,IKTE 
    DO JIJ=IIJB,IIJE 
      PAMOIST(JIJ,JK) = (1.0-PFRAC_ICE(JIJ,JK))*ZAMOIST_W(JIJ,JK) &
      +PFRAC_ICE(JIJ,JK) *ZAMOIST_I(JIJ,JK)
      PATHETA(JIJ,JK) = (1.0-PFRAC_ICE(JIJ,JK))*ZATHETA_W(JIJ,JK) &
      +PFRAC_ICE(JIJ,JK) *ZATHETA_I(JIJ,JK)
    ENDDO
  ENDDO
  !
ELSE
  PAMOIST(IIJB:IIJE,IKTB:IKTE) = 0.
  PATHETA(IIJB:IIJE,IKTB:IKTE) = 0.
ENDIF
IF (LHOOK) CALL DR_HOOK('COMPUTE_FUNCTION_THERMO_MF',1,ZHOOK_HANDLE)
END SUBROUTINE COMPUTE_FUNCTION_THERMO_MF
!
END MODULE MODE_COMPUTE_FUNCTION_THERMO_MF
