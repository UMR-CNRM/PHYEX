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
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=IKTB:IKTE)
  ZCP(IIJB:IIJE,IKTB:IKTE) = ZCP(IIJB:IIJE,IKTB:IKTE) + CST%XCPV * PR(IIJB:IIJE,IKTB:IKTE,1)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=IKTB:IKTE)
ENDIF

DO JRR = 2,1+KRRL  ! loop on the liquid components
   !$mnh_expand_array(JIJ=IIJB:IIJE,JK=IKTB:IKTE)
   ZCP(IIJB:IIJE,IKTB:IKTE)  = ZCP(IIJB:IIJE,IKTB:IKTE) + CST%XCL * PR(IIJB:IIJE,IKTB:IKTE,JRR)
   !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=IKTB:IKTE)
END DO

DO JRR = 2+KRRL,1+KRRL+KRRI ! loop on the solid components
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=IKTB:IKTE)
  ZCP(IIJB:IIJE,IKTB:IKTE)  = ZCP(IIJB:IIJE,IKTB:IKTE)  + CST%XCI * PR(IIJB:IIJE,IKTB:IKTE,JRR)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=IKTB:IKTE)

END DO

!*      Temperature
!
!$mnh_expand_array(JIJ=IIJB:IIJE,JK=IKTB:IKTE)
PT(IIJB:IIJE,IKTB:IKTE) =  PTH(IIJB:IIJE,IKTB:IKTE) * PEXN(IIJB:IIJE,IKTB:IKTE)
!$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=IKTB:IKTE)
!
!
!! Liquid water
!
IF ( KRRL >= 1 ) THEN
  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=IKTB:IKTE)
  !
  !*       Lv/Cph
  !
  ZLVOCP(IIJB:IIJE,IKTB:IKTE) = (CST%XLVTT + (CST%XCPV-CST%XCL) *  (PT(IIJB:IIJE,IKTB:IKTE)-CST%XTT) ) / &
                                      & ZCP(IIJB:IIJE,IKTB:IKTE)
  !
  !*      Saturation vapor pressure with respect to water
  !
  ZE(IIJB:IIJE,IKTB:IKTE) =  EXP(CST%XALPW - CST%XBETAW/PT(IIJB:IIJE,IKTB:IKTE) - &
                                        &CST%XGAMW*ALOG( PT(IIJB:IIJE,IKTB:IKTE) ) )
  !
  !*      Saturation  mixing ratio with respect to water
  !
  ZE(IIJB:IIJE,IKTB:IKTE) =  ZE(IIJB:IIJE,IKTB:IKTE) * ZEPS / &
                                  & ( PPABS(IIJB:IIJE,IKTB:IKTE) - ZE(IIJB:IIJE,IKTB:IKTE) )
  !
  !*      Compute the saturation mixing ratio derivative (rvs')
  !
  ZDEDT(IIJB:IIJE,IKTB:IKTE) = (CST%XBETAW/PT(IIJB:IIJE,IKTB:IKTE)  - CST%XGAMW) / PT(IIJB:IIJE,IKTB:IKTE)&
                 * ZE(IIJB:IIJE,IKTB:IKTE) * ( 1. + ZE(IIJB:IIJE,IKTB:IKTE) / ZEPS )
  !
  !*      Compute Amoist and Atheta
  !
  IF (OSTATNW) THEN
    ZAMOIST_W(IIJB:IIJE,IKTB:IKTE)=  1.0/( 1.0 + ZDEDT(IIJB:IIJE,IKTB:IKTE) * ZLVOCP(IIJB:IIJE,IKTB:IKTE))
    ZATHETA_W(IIJB:IIJE,IKTB:IKTE)= ZAMOIST_W(IIJB:IIJE,IKTB:IKTE) * PEXN(IIJB:IIJE,IKTB:IKTE) &
                                            * ZDEDT(IIJB:IIJE,IKTB:IKTE)
  ELSE
    ZAMOIST_W(IIJB:IIJE,IKTB:IKTE)= 0.5/( 1.0 + ZDEDT(IIJB:IIJE,IKTB:IKTE) * ZLVOCP(IIJB:IIJE,IKTB:IKTE) )
    ZATHETA_W(IIJB:IIJE,IKTB:IKTE)= ZAMOIST_W(IIJB:IIJE,IKTB:IKTE) * PEXN(IIJB:IIJE,IKTB:IKTE) *         &
          ( ( ZE(IIJB:IIJE,IKTB:IKTE) - PR(IIJB:IIJE,IKTB:IKTE,1) ) * ZLVOCP(IIJB:IIJE,IKTB:IKTE) /      &
            ( 1. + ZDEDT(IIJB:IIJE,IKTB:IKTE) * ZLVOCP(IIJB:IIJE,IKTB:IKTE) )           *              &
            (                                                             &
             ZE(IIJB:IIJE,IKTB:IKTE) * (1. + ZE(IIJB:IIJE,IKTB:IKTE)/ZEPS)                                &
                            * ( -2.*CST%XBETAW/PT(IIJB:IIJE,IKTB:IKTE) + CST%XGAMW ) / PT(IIJB:IIJE,IKTB:IKTE)**2&
            +ZDEDT(IIJB:IIJE,IKTB:IKTE) * (1. + 2. * ZE(IIJB:IIJE,IKTB:IKTE)/ZEPS)                        &
                          * ( CST%XBETAW/PT(IIJB:IIJE,IKTB:IKTE) - CST%XGAMW ) / PT(IIJB:IIJE,IKTB:IKTE)         &
            )                                                             &
           - ZDEDT(IIJB:IIJE,IKTB:IKTE)                                                   &
          )
  END IF
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=IKTB:IKTE)
  !
  !! Solid water
  !
  IF ( KRRI >= 1 ) THEN
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=IKTB:IKTE)
    !
    !*       Ls/Cph
    !
    ZLSOCP(IIJB:IIJE,IKTB:IKTE) = (CST%XLSTT + (CST%XCPV-CST%XCI) *  (PT(IIJB:IIJE,IKTB:IKTE)-CST%XTT) ) / &
                                        & ZCP(IIJB:IIJE,IKTB:IKTE)
    !
    !*      Saturation vapor pressure with respect to ice
    !
    ZE(IIJB:IIJE,IKTB:IKTE) =  EXP(CST%XALPI - CST%XBETAI/PT(IIJB:IIJE,IKTB:IKTE) - &
                                          &CST%XGAMI*ALOG( PT(IIJB:IIJE,IKTB:IKTE) ) )
    !
    !*      Saturation  mixing ratio with respect to ice
    !
    ZE(IIJB:IIJE,IKTB:IKTE) =  ZE(IIJB:IIJE,IKTB:IKTE) * ZEPS / &
                                    & ( PPABS(IIJB:IIJE,IKTB:IKTE) - ZE(IIJB:IIJE,IKTB:IKTE) )
    !
    !*      Compute the saturation mixing ratio derivative (rvs')
    !
    ZDEDT(IIJB:IIJE,IKTB:IKTE) = (CST%XBETAI/PT(IIJB:IIJE,IKTB:IKTE)-CST%XGAMI) /PT(IIJB:IIJE,IKTB:IKTE)&
                   * ZE(IIJB:IIJE,IKTB:IKTE) * ( 1. + ZE(IIJB:IIJE,IKTB:IKTE) / ZEPS )
    !
    !*      Compute Amoist and Atheta
    !
    IF (OSTATNW) THEN
      ZAMOIST_I(IIJB:IIJE,IKTB:IKTE)= 1.0/( 1.0 + ZDEDT(IIJB:IIJE,IKTB:IKTE) *ZLVOCP(IIJB:IIJE,IKTB:IKTE))
      ZATHETA_I(IIJB:IIJE,IKTB:IKTE)= ZAMOIST_I(IIJB:IIJE,IKTB:IKTE) * PEXN(IIJB:IIJE,IKTB:IKTE) &
                                            * ZDEDT(IIJB:IIJE,IKTB:IKTE)
    ELSE
      ZAMOIST_I(IIJB:IIJE,IKTB:IKTE)= 0.5/(1.0 + ZDEDT(IIJB:IIJE,IKTB:IKTE) * ZLSOCP(IIJB:IIJE,IKTB:IKTE))
      ZATHETA_I(IIJB:IIJE,IKTB:IKTE)= ZAMOIST_I(IIJB:IIJE,IKTB:IKTE) * PEXN(IIJB:IIJE,IKTB:IKTE) *      &
        ( ( ZE(IIJB:IIJE,IKTB:IKTE) - PR(IIJB:IIJE,IKTB:IKTE,1) ) * ZLSOCP(IIJB:IIJE,IKTB:IKTE) /       &
          ( 1. + ZDEDT(IIJB:IIJE,IKTB:IKTE) * ZLSOCP(IIJB:IIJE,IKTB:IKTE) )           *              &
          (                                                             &
           ZE(IIJB:IIJE,IKTB:IKTE) * (1. + ZE(IIJB:IIJE,IKTB:IKTE)/ZEPS)                                &
                        * ( -2.*CST%XBETAI/PT(IIJB:IIJE,IKTB:IKTE) + CST%XGAMI ) / PT(IIJB:IIJE,IKTB:IKTE)**2   &
          +ZDEDT(IIJB:IIJE,IKTB:IKTE) * (1. + 2. * ZE(IIJB:IIJE,IKTB:IKTE)/ZEPS)                        &
                        * ( CST%XBETAI/PT(IIJB:IIJE,IKTB:IKTE) - CST%XGAMI ) / PT(IIJB:IIJE,IKTB:IKTE)          &
          )                                                             &
         - ZDEDT(IIJB:IIJE,IKTB:IKTE)                                                   &
        )
    END IF
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=IKTB:IKTE)

  ELSE
    ZAMOIST_I(IIJB:IIJE,IKTB:IKTE)=0.
    ZATHETA_I(IIJB:IIJE,IKTB:IKTE)=0.
  ENDIF

  !$mnh_expand_array(JIJ=IIJB:IIJE,JK=IKTB:IKTE)
  PAMOIST(IIJB:IIJE,IKTB:IKTE) = (1.0-PFRAC_ICE(IIJB:IIJE,IKTB:IKTE))*ZAMOIST_W(IIJB:IIJE,IKTB:IKTE) &
                         +PFRAC_ICE(IIJB:IIJE,IKTB:IKTE) *ZAMOIST_I(IIJB:IIJE,IKTB:IKTE)
  PATHETA(IIJB:IIJE,IKTB:IKTE) = (1.0-PFRAC_ICE(IIJB:IIJE,IKTB:IKTE))*ZATHETA_W(IIJB:IIJE,IKTB:IKTE) &
                         +PFRAC_ICE(IIJB:IIJE,IKTB:IKTE) *ZATHETA_I(IIJB:IIJE,IKTB:IKTE)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=IKTB:IKTE)
  !
ELSE
  PAMOIST(IIJB:IIJE,IKTB:IKTE) = 0.
  PATHETA(IIJB:IIJE,IKTB:IKTE) = 0.
ENDIF
IF (LHOOK) CALL DR_HOOK('COMPUTE_FUNCTION_THERMO_MF',1,ZHOOK_HANDLE)
END SUBROUTINE COMPUTE_FUNCTION_THERMO_MF
!
END MODULE MODE_COMPUTE_FUNCTION_THERMO_MF
