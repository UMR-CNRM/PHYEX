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

INTEGER             :: JRR, JI, JK
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('COMPUTE_FUNCTION_THERMO_MF',0,ZHOOK_HANDLE)
!
  ZEPS = CST%XMV / CST%XMD

!
!*       Cph
!
ZCP=CST%XCPD

IF (KRR > 0) THEN
  DO JK=D%NKTB,D%NKTE 
    DO JI=D%NIJB,D%NIJE 
      ZCP(JI,JK) = ZCP(JI,JK) + CST%XCPV * PR(JI,JK,1)
    ENDDO
  ENDDO
ENDIF

DO JRR = 2,1+KRRL  ! loop on the liquid components
  DO JK=D%NKTB,D%NKTE 
    DO JI=D%NIJB,D%NIJE 
      ZCP(JI,JK)  = ZCP(JI,JK) + CST%XCL * PR(JI,JK,JRR)
    ENDDO
  ENDDO
END DO

DO JRR = 2+KRRL,1+KRRL+KRRI ! loop on the solid components
  DO JK=D%NKTB,D%NKTE 
    DO JI=D%NIJB,D%NIJE 
      ZCP(JI,JK)  = ZCP(JI,JK)  + CST%XCI * PR(JI,JK,JRR)
    ENDDO
  ENDDO

END DO

!*      Temperature
!
DO JK=D%NKTB,D%NKTE 
  DO JI=D%NIJB,D%NIJE 
    PT(JI,JK) =  PTH(JI,JK) * PEXN(JI,JK)
  ENDDO
ENDDO
!
!
!! Liquid water
!
IF ( KRRL >= 1 ) THEN
  DO JK=D%NKTB,D%NKTE 
    DO JI=D%NIJB,D%NIJE 
  !
  !*       Lv/Cph
  !
      ZLVOCP(JI,JK) = (CST%XLVTT + (CST%XCPV-CST%XCL) *  (PT(JI,JK)-CST%XTT) ) / &
      & ZCP(JI,JK)
  !
  !*      Saturation vapor pressure with respect to water
  !
      ZE(JI,JK) =  EXP(CST%XALPW - CST%XBETAW/PT(JI,JK) - &
      &CST%XGAMW*ALOG( PT(JI,JK) ) )
  !
  !*      Saturation  mixing ratio with respect to water
  !
      ZE(JI,JK) =  ZE(JI,JK) * ZEPS / &
      & ( PPABS(JI,JK) - ZE(JI,JK) )
  !
  !*      Compute the saturation mixing ratio derivative (rvs')
  !
      ZDEDT(JI,JK) = (CST%XBETAW/PT(JI,JK)  - CST%XGAMW) / PT(JI,JK)&
      * ZE(JI,JK) * ( 1. + ZE(JI,JK) / ZEPS )
  !
  !*      Compute Amoist and Atheta
  !
      IF (OSTATNW) THEN
        ZAMOIST_W(JI,JK)=  1.0/( 1.0 + ZDEDT(JI,JK) * ZLVOCP(JI,JK))
        ZATHETA_W(JI,JK)= ZAMOIST_W(JI,JK) * PEXN(JI,JK) &
        * ZDEDT(JI,JK)
      ELSE
        ZAMOIST_W(JI,JK)= 0.5/( 1.0 + ZDEDT(JI,JK) * ZLVOCP(JI,JK) )
        ZATHETA_W(JI,JK)= ZAMOIST_W(JI,JK) * PEXN(JI,JK) *         &
        ( ( ZE(JI,JK) - PR(JI,JK,1) ) * ZLVOCP(JI,JK) /      &
        ( 1. + ZDEDT(JI,JK) * ZLVOCP(JI,JK) )           *              &
        (                                                             &
        ZE(JI,JK) * (1. + ZE(JI,JK)/ZEPS)                                &
        * ( -2.*CST%XBETAW/PT(JI,JK) + CST%XGAMW ) / PT(JI,JK)**2&
        +ZDEDT(JI,JK) * (1. + 2. * ZE(JI,JK)/ZEPS)                        &
        * ( CST%XBETAW/PT(JI,JK) - CST%XGAMW ) / PT(JI,JK)         &
        )                                                             &
        - ZDEDT(JI,JK)                                                   &
        )
      END IF
    ENDDO
  ENDDO
  !
  !! Solid water
  !
  IF ( KRRI >= 1 ) THEN
    DO JK=D%NKTB,D%NKTE 
      DO JI=D%NIJB,D%NIJE 
    !
    !*       Ls/Cph
    !
        ZLSOCP(JI,JK) = (CST%XLSTT + (CST%XCPV-CST%XCI) *  (PT(JI,JK)-CST%XTT) ) / &
        & ZCP(JI,JK)
    !
    !*      Saturation vapor pressure with respect to ice
    !
        ZE(JI,JK) =  EXP(CST%XALPI - CST%XBETAI/PT(JI,JK) - &
        &CST%XGAMI*ALOG( PT(JI,JK) ) )
    !
    !*      Saturation  mixing ratio with respect to ice
    !
        ZE(JI,JK) =  ZE(JI,JK) * ZEPS / &
        & ( PPABS(JI,JK) - ZE(JI,JK) )
    !
    !*      Compute the saturation mixing ratio derivative (rvs')
    !
        ZDEDT(JI,JK) = (CST%XBETAI/PT(JI,JK)-CST%XGAMI) /PT(JI,JK)&
        * ZE(JI,JK) * ( 1. + ZE(JI,JK) / ZEPS )
    !
    !*      Compute Amoist and Atheta
    !
        IF (OSTATNW) THEN
          ZAMOIST_I(JI,JK)= 1.0/( 1.0 + ZDEDT(JI,JK) *ZLVOCP(JI,JK))
          ZATHETA_I(JI,JK)= ZAMOIST_I(JI,JK) * PEXN(JI,JK) &
          * ZDEDT(JI,JK)
        ELSE
          ZAMOIST_I(JI,JK)= 0.5/(1.0 + ZDEDT(JI,JK) * ZLSOCP(JI,JK))
          ZATHETA_I(JI,JK)= ZAMOIST_I(JI,JK) * PEXN(JI,JK) *      &
          ( ( ZE(JI,JK) - PR(JI,JK,1) ) * ZLSOCP(JI,JK) /       &
          ( 1. + ZDEDT(JI,JK) * ZLSOCP(JI,JK) )           *              &
          (                                                             &
          ZE(JI,JK) * (1. + ZE(JI,JK)/ZEPS)                                &
          * ( -2.*CST%XBETAI/PT(JI,JK) + CST%XGAMI ) / PT(JI,JK)**2   &
          +ZDEDT(JI,JK) * (1. + 2. * ZE(JI,JK)/ZEPS)                        &
          * ( CST%XBETAI/PT(JI,JK) - CST%XGAMI ) / PT(JI,JK)          &
          )                                                             &
          - ZDEDT(JI,JK)                                                   &
          )
        END IF
      ENDDO
    ENDDO

  ELSE
    ZAMOIST_I(D%NIJB:D%NIJE,D%NKTB:D%NKTE)=0.
    ZATHETA_I(D%NIJB:D%NIJE,D%NKTB:D%NKTE)=0.
  ENDIF

  DO JK=D%NKTB,D%NKTE 
    DO JI=D%NIJB,D%NIJE 
      PAMOIST(JI,JK) = (1.0-PFRAC_ICE(JI,JK))*ZAMOIST_W(JI,JK) &
      +PFRAC_ICE(JI,JK) *ZAMOIST_I(JI,JK)
      PATHETA(JI,JK) = (1.0-PFRAC_ICE(JI,JK))*ZATHETA_W(JI,JK) &
      +PFRAC_ICE(JI,JK) *ZATHETA_I(JI,JK)
    ENDDO
  ENDDO
  !
ELSE
  PAMOIST(D%NIJB:D%NIJE,D%NKTB:D%NKTE) = 0.
  PATHETA(D%NIJB:D%NIJE,D%NKTB:D%NKTE) = 0.
ENDIF
IF (LHOOK) CALL DR_HOOK('COMPUTE_FUNCTION_THERMO_MF',1,ZHOOK_HANDLE)
END SUBROUTINE COMPUTE_FUNCTION_THERMO_MF
!
END MODULE MODE_COMPUTE_FUNCTION_THERMO_MF
