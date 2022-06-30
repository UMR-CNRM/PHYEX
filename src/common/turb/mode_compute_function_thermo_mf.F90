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
      SUBROUTINE COMPUTE_FUNCTION_THERMO_MF(D, CST, KRR,KRRL,KRRI,                  &
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
REAL, DIMENSION(D%NIJT,D%NKT) ::       &
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
  !$mnh_expand_array(JI=D%NIJB:D%NIJE,JK=D%NKTB:D%NKTE)
  ZCP(D%NIJB:D%NIJE,D%NKTB:D%NKTE) = ZCP(D%NIJB:D%NIJE,D%NKTB:D%NKTE) + CST%XCPV * PR(D%NIJB:D%NIJE,D%NKTB:D%NKTE,1)
  !$mnh_end_expand_array(JI=D%NIJB:D%NIJE,JK=D%NKTB:D%NKTE)
ENDIF

DO JRR = 2,1+KRRL  ! loop on the liquid components  
   !$mnh_expand_array(JI=D%NIJB:D%NIJE,JK=D%NKTB:D%NKTE)
   ZCP(D%NIJB:D%NIJE,D%NKTB:D%NKTE)  = ZCP(D%NIJB:D%NIJE,D%NKTB:D%NKTE) + CST%XCL * PR(D%NIJB:D%NIJE,D%NKTB:D%NKTE,JRR)
   !$mnh_end_expand_array(JI=D%NIJB:D%NIJE,JK=D%NKTB:D%NKTE)
END DO

DO JRR = 2+KRRL,1+KRRL+KRRI ! loop on the solid components
  !$mnh_expand_array(JI=D%NIJB:D%NIJE,JK=D%NKTB:D%NKTE)
  ZCP(D%NIJB:D%NIJE,D%NKTB:D%NKTE)  = ZCP(D%NIJB:D%NIJE,D%NKTB:D%NKTE)  + CST%XCI * PR(D%NIJB:D%NIJE,D%NKTB:D%NKTE,JRR)
  !$mnh_end_expand_array(JI=D%NIJB:D%NIJE,JK=D%NKTB:D%NKTE)
END DO

!*      Temperature
!
!$mnh_expand_array(JI=D%NIJB:D%NIJE,JK=D%NKTB:D%NKTE)
PT(D%NIJB:D%NIJE,D%NKTB:D%NKTE) =  PTH(D%NIJB:D%NIJE,D%NKTB:D%NKTE) * PEXN(D%NIJB:D%NIJE,D%NKTB:D%NKTE)
!$mnh_end_expand_array(JI=D%NIJB:D%NIJE,JK=D%NKTB:D%NKTE)
!
!
!! Liquid water
!
IF ( KRRL >= 1 ) THEN 
  !$mnh_expand_array(JI=D%NIJB:D%NIJE,JK=D%NKTB:D%NKTE)
  !
  !*       Lv/Cph 
  !
  ZLVOCP(D%NIJB:D%NIJE,D%NKTB:D%NKTE) = (CST%XLVTT + (CST%XCPV-CST%XCL) *  (PT(D%NIJB:D%NIJE,D%NKTB:D%NKTE)-CST%XTT) ) / &
                                      & ZCP(D%NIJB:D%NIJE,D%NKTB:D%NKTE)
  !
  !*      Saturation vapor pressure with respect to water
  !
  ZE(D%NIJB:D%NIJE,D%NKTB:D%NKTE) =  EXP(CST%XALPW - CST%XBETAW/PT(D%NIJB:D%NIJE,D%NKTB:D%NKTE) - &
                                        &CST%XGAMW*ALOG( PT(D%NIJB:D%NIJE,D%NKTB:D%NKTE) ) )
  !
  !*      Saturation  mixing ratio with respect to water
  !
  ZE(D%NIJB:D%NIJE,D%NKTB:D%NKTE) =  ZE(D%NIJB:D%NIJE,D%NKTB:D%NKTE) * ZEPS / &
                                  & ( PPABS(D%NIJB:D%NIJE,D%NKTB:D%NKTE) - ZE(D%NIJB:D%NIJE,D%NKTB:D%NKTE) )
  !
  !*      Compute the saturation mixing ratio derivative (rvs')
  !
  ZDEDT(D%NIJB:D%NIJE,D%NKTB:D%NKTE) = (CST%XBETAW/PT(D%NIJB:D%NIJE,D%NKTB:D%NKTE)  - CST%XGAMW) / PT(D%NIJB:D%NIJE,D%NKTB:D%NKTE)&
                 * ZE(D%NIJB:D%NIJE,D%NKTB:D%NKTE) * ( 1. + ZE(D%NIJB:D%NIJE,D%NKTB:D%NKTE) / ZEPS )
  !
  !*      Compute Amoist
  !
  ZAMOIST_W(D%NIJB:D%NIJE,D%NKTB:D%NKTE)=  0.5 / ( 1.0 + ZDEDT(D%NIJB:D%NIJE,D%NKTB:D%NKTE) * ZLVOCP(D%NIJB:D%NIJE,D%NKTB:D%NKTE) )
  !
  !*      Compute Atheta
  !
  ZATHETA_W(D%NIJB:D%NIJE,D%NKTB:D%NKTE)= ZAMOIST_W(D%NIJB:D%NIJE,D%NKTB:D%NKTE) * PEXN(D%NIJB:D%NIJE,D%NKTB:D%NKTE) *         &
        ( ( ZE(D%NIJB:D%NIJE,D%NKTB:D%NKTE) - PR(D%NIJB:D%NIJE,D%NKTB:D%NKTE,1) ) * ZLVOCP(D%NIJB:D%NIJE,D%NKTB:D%NKTE) /      &
          ( 1. + ZDEDT(D%NIJB:D%NIJE,D%NKTB:D%NKTE) * ZLVOCP(D%NIJB:D%NIJE,D%NKTB:D%NKTE) )           *              &
          (                                                             &
           ZE(D%NIJB:D%NIJE,D%NKTB:D%NKTE) * (1. + ZE(D%NIJB:D%NIJE,D%NKTB:D%NKTE)/ZEPS)                                &
                        * ( -2.*CST%XBETAW/PT(D%NIJB:D%NIJE,D%NKTB:D%NKTE) + CST%XGAMW ) / PT(D%NIJB:D%NIJE,D%NKTB:D%NKTE)**2   &
          +ZDEDT(D%NIJB:D%NIJE,D%NKTB:D%NKTE) * (1. + 2. * ZE(D%NIJB:D%NIJE,D%NKTB:D%NKTE)/ZEPS)                        &
                        * ( CST%XBETAW/PT(D%NIJB:D%NIJE,D%NKTB:D%NKTE) - CST%XGAMW ) / PT(D%NIJB:D%NIJE,D%NKTB:D%NKTE)          &
          )                                                             &
         - ZDEDT(D%NIJB:D%NIJE,D%NKTB:D%NKTE)                                                   &
        )
  !$mnh_end_expand_array(JI=D%NIJB:D%NIJE,JK=D%NKTB:D%NKTE)
  !
  !! Solid water
  !
  IF ( KRRI >= 1 ) THEN 
    !$mnh_expand_array(JI=D%NIJB:D%NIJE,JK=D%NKTB:D%NKTE)
    !
    !*       Ls/Cph 
    !
    ZLSOCP(D%NIJB:D%NIJE,D%NKTB:D%NKTE) = (CST%XLSTT + (CST%XCPV-CST%XCI) *  (PT(D%NIJB:D%NIJE,D%NKTB:D%NKTE)-CST%XTT) ) / &
                                        & ZCP(D%NIJB:D%NIJE,D%NKTB:D%NKTE)
    !
    !*      Saturation vapor pressure with respect to ice
    !
    ZE(D%NIJB:D%NIJE,D%NKTB:D%NKTE) =  EXP(CST%XALPI - CST%XBETAI/PT(D%NIJB:D%NIJE,D%NKTB:D%NKTE) - &
                                          &CST%XGAMI*ALOG( PT(D%NIJB:D%NIJE,D%NKTB:D%NKTE) ) )
    !
    !*      Saturation  mixing ratio with respect to ice
    !
    ZE(D%NIJB:D%NIJE,D%NKTB:D%NKTE) =  ZE(D%NIJB:D%NIJE,D%NKTB:D%NKTE) * ZEPS / &
                                    & ( PPABS(D%NIJB:D%NIJE,D%NKTB:D%NKTE) - ZE(D%NIJB:D%NIJE,D%NKTB:D%NKTE) )
    !
    !*      Compute the saturation mixing ratio derivative (rvs')
    !
    ZDEDT(D%NIJB:D%NIJE,D%NKTB:D%NKTE) = (CST%XBETAI/PT(D%NIJB:D%NIJE,D%NKTB:D%NKTE) - CST%XGAMI) /PT(D%NIJB:D%NIJE,D%NKTB:D%NKTE)&
                   * ZE(D%NIJB:D%NIJE,D%NKTB:D%NKTE) * ( 1. + ZE(D%NIJB:D%NIJE,D%NKTB:D%NKTE) / ZEPS )
    !
    !*      Compute Amoist
    !
    ZAMOIST_I(D%NIJB:D%NIJE,D%NKTB:D%NKTE)=  0.5/(1.0 + ZDEDT(D%NIJB:D%NIJE,D%NKTB:D%NKTE) * ZLSOCP(D%NIJB:D%NIJE,D%NKTB:D%NKTE))
    !
    !*      Compute Atheta
    !
    ZATHETA_I(D%NIJB:D%NIJE,D%NKTB:D%NKTE)= ZAMOIST_I(D%NIJB:D%NIJE,D%NKTB:D%NKTE) * PEXN(D%NIJB:D%NIJE,D%NKTB:D%NKTE) *        &
        ( ( ZE(D%NIJB:D%NIJE,D%NKTB:D%NKTE) - PR(D%NIJB:D%NIJE,D%NKTB:D%NKTE,1) ) * ZLSOCP(D%NIJB:D%NIJE,D%NKTB:D%NKTE) /       &
          ( 1. + ZDEDT(D%NIJB:D%NIJE,D%NKTB:D%NKTE) * ZLSOCP(D%NIJB:D%NIJE,D%NKTB:D%NKTE) )           *              &
          (                                                             &
           ZE(D%NIJB:D%NIJE,D%NKTB:D%NKTE) * (1. + ZE(D%NIJB:D%NIJE,D%NKTB:D%NKTE)/ZEPS)                                &
                        * ( -2.*CST%XBETAI/PT(D%NIJB:D%NIJE,D%NKTB:D%NKTE) + CST%XGAMI ) / PT(D%NIJB:D%NIJE,D%NKTB:D%NKTE)**2   &
          +ZDEDT(D%NIJB:D%NIJE,D%NKTB:D%NKTE) * (1. + 2. * ZE(D%NIJB:D%NIJE,D%NKTB:D%NKTE)/ZEPS)                        &
                        * ( CST%XBETAI/PT(D%NIJB:D%NIJE,D%NKTB:D%NKTE) - CST%XGAMI ) / PT(D%NIJB:D%NIJE,D%NKTB:D%NKTE)          &
          )                                                             &
         - ZDEDT(D%NIJB:D%NIJE,D%NKTB:D%NKTE)                                                   &
        )
    !$mnh_end_expand_array(JI=D%NIJB:D%NIJE,JK=D%NKTB:D%NKTE)
  ELSE
    ZAMOIST_I(D%NIJB:D%NIJE,D%NKTB:D%NKTE)=0.
    ZATHETA_I(D%NIJB:D%NIJE,D%NKTB:D%NKTE)=0.
  ENDIF

  !$mnh_expand_array(JI=D%NIJB:D%NIJE,JK=D%NKTB:D%NKTE)
  PAMOIST(D%NIJB:D%NIJE,D%NKTB:D%NKTE) = (1.0-PFRAC_ICE(D%NIJB:D%NIJE,D%NKTB:D%NKTE))*ZAMOIST_W(D%NIJB:D%NIJE,D%NKTB:D%NKTE) &
                         +PFRAC_ICE(D%NIJB:D%NIJE,D%NKTB:D%NKTE) *ZAMOIST_I(D%NIJB:D%NIJE,D%NKTB:D%NKTE)
  PATHETA(D%NIJB:D%NIJE,D%NKTB:D%NKTE) = (1.0-PFRAC_ICE(D%NIJB:D%NIJE,D%NKTB:D%NKTE))*ZATHETA_W(D%NIJB:D%NIJE,D%NKTB:D%NKTE) &
                         +PFRAC_ICE(D%NIJB:D%NIJE,D%NKTB:D%NKTE) *ZATHETA_I(D%NIJB:D%NIJE,D%NKTB:D%NKTE)
  !$mnh_end_expand_array(JI=D%NIJB:D%NIJE,JK=D%NKTB:D%NKTE)
  !
ELSE
  PAMOIST(D%NIJB:D%NIJE,D%NKTB:D%NKTE) = 0.
  PATHETA(D%NIJB:D%NIJE,D%NKTB:D%NKTE) = 0.
ENDIF
IF (LHOOK) CALL DR_HOOK('COMPUTE_FUNCTION_THERMO_MF',1,ZHOOK_HANDLE)
END SUBROUTINE COMPUTE_FUNCTION_THERMO_MF
!
END MODULE MODE_COMPUTE_FUNCTION_THERMO_MF
