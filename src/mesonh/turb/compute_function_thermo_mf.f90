!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!     ######spl
     MODULE MODI_COMPUTE_FUNCTION_THERMO_MF
!    ######################################
!
INTERFACE
      
!     #################################################################
      SUBROUTINE COMPUTE_FUNCTION_THERMO_MF( KRR,KRRL,KRRI,                  &
                                       PTH, PR, PEXN, PFRAC_ICE, PPABS,      &
                                       PT, PAMOIST,PATHETA                   )
!     #################################################################

!*               1.1  Declaration of Arguments
!

INTEGER,                INTENT(IN)   :: KRR           ! number of moist var.
INTEGER,                INTENT(IN)   :: KRRL          ! number of liquid water var.
INTEGER,                INTENT(IN)   :: KRRI          ! number of ice water var.

REAL, DIMENSION(:,:), INTENT(IN)   :: PTH      ! theta
REAL, DIMENSION(:,:,:), INTENT(IN) :: PR       ! water species
REAL, DIMENSION(:,:)  , INTENT(IN) :: PPABS,PEXN    ! pressure, Exner funct.
REAL, DIMENSION(:,:)  , INTENT(IN) :: PFRAC_ICE     ! ice fraction

REAL, DIMENSION(:,:), INTENT(OUT)   :: PT      ! temperature

REAL, DIMENSION(:,:), INTENT(OUT)  ::  PAMOIST,PATHETA
!
END SUBROUTINE COMPUTE_FUNCTION_THERMO_MF

END INTERFACE
!
END MODULE MODI_COMPUTE_FUNCTION_THERMO_MF
!     ######spl
      SUBROUTINE COMPUTE_FUNCTION_THERMO_MF( KRR,KRRL,KRRI,                  &
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
USE MODD_CST
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
REAL, DIMENSION(:,:)  , INTENT(IN) :: PPABS,PEXN    ! pressure, Exner funct.
REAL, DIMENSION(:,:)  , INTENT(IN) :: PFRAC_ICE     ! ice fraction

REAL, DIMENSION(:,:), INTENT(OUT)   :: PT      ! temperature

REAL, DIMENSION(:,:), INTENT(OUT)  ::  PAMOIST,PATHETA
!
!-------------------------------------------------------------------------------
! 
!*       0.2   Declarations of local variables
!
REAL                :: ZEPS         ! XMV / XMD
REAL, DIMENSION(SIZE(PTH,1),SIZE(PTH,2)) ::     &
          ZCP,                        &  ! Cp 
          ZE,                         &  ! Saturation mixing ratio
          ZDEDT,                      &  ! Saturation mixing ratio derivative
          ZAMOIST_W,                  &  ! Coefficients for s = f (Thetal,Rnp)
          ZATHETA_W,                  &  !
          ZAMOIST_I,                  &  !
          ZATHETA_I,                  &  !
          ZLVOCP,ZLSOCP

INTEGER             :: JRR
!
!-------------------------------------------------------------------------------
!
!
  ZEPS = XMV / XMD

!
!*       Cph
!
ZCP=XCPD

IF (KRR > 0) ZCP(:,:) = ZCP(:,:) + XCPV * PR(:,:,1)

DO JRR = 2,1+KRRL  ! loop on the liquid components  
   ZCP(:,:)  = ZCP(:,:) + XCL * PR(:,:,JRR)
END DO

DO JRR = 2+KRRL,1+KRRL+KRRI ! loop on the solid components   
  ZCP(:,:)  = ZCP(:,:)  + XCI * PR(:,:,JRR)
END DO

!*      Temperature
!
PT(:,:) =  PTH(:,:) * PEXN(:,:)
!
!
!! Liquid water
!
IF ( KRRL >= 1 ) THEN 
!
!*       Lv/Cph 
!
  ZLVOCP(:,:) = (XLVTT + (XCPV-XCL) *  (PT(:,:)-XTT) ) / ZCP(:,:)
!
!*      Saturation vapor pressure with respect to water
!
  ZE(:,:) =  EXP( XALPW - XBETAW/PT(:,:) - XGAMW*ALOG( PT(:,:) ) )
!
!*      Saturation  mixing ratio with respect to water
!
  ZE(:,:) =  ZE(:,:) * ZEPS / ( PPABS(:,:) - ZE(:,:) )
!
!*      Compute the saturation mixing ratio derivative (rvs')
!
  ZDEDT(:,:) = ( XBETAW / PT(:,:)  - XGAMW ) / PT(:,:)   &
                 * ZE(:,:) * ( 1. + ZE(:,:) / ZEPS )
!
!*      Compute Amoist
!
  ZAMOIST_W(:,:)=  0.5 / ( 1.0 + ZDEDT(:,:) * ZLVOCP(:,:) )
!
!*      Compute Atheta
!
  ZATHETA_W(:,:)= ZAMOIST_W(:,:) * PEXN(:,:) *                          &
        ( ( ZE(:,:) - PR(:,:,1) ) * ZLVOCP(:,:) /                    &
          ( 1. + ZDEDT(:,:) * ZLVOCP(:,:) )           *              &
          (                                                             &
           ZE(:,:) * (1. + ZE(:,:)/ZEPS)                                &
                        * ( -2.*XBETAW/PT(:,:) + XGAMW ) / PT(:,:)**2   &
          +ZDEDT(:,:) * (1. + 2. * ZE(:,:)/ZEPS)                        &
                        * ( XBETAW/PT(:,:) - XGAMW ) / PT(:,:)          &
          )                                                             &
         - ZDEDT(:,:)                                                   &
        )

!
!! Solid water
!
  IF ( KRRI >= 1 ) THEN 

!
!*       Ls/Cph 
!
    ZLSOCP(:,:) = (XLSTT + (XCPV-XCI) *  (PT(:,:)-XTT) ) / ZCP(:,:)
!
!*      Saturation vapor pressure with respect to ice
!
    ZE(:,:) =  EXP( XALPI - XBETAI/PT(:,:) - XGAMI*ALOG( PT(:,:) ) )
!
!*      Saturation  mixing ratio with respect to ice
!
    ZE(:,:) =  ZE(:,:) * ZEPS / ( PPABS(:,:) - ZE(:,:) )
!
!*      Compute the saturation mixing ratio derivative (rvs')
!
    ZDEDT(:,:) = ( XBETAI / PT(:,:)  - XGAMI ) / PT(:,:)   &
                   * ZE(:,:) * ( 1. + ZE(:,:) / ZEPS )
!
!*      Compute Amoist
!
    ZAMOIST_I(:,:)=  0.5 / ( 1.0 + ZDEDT(:,:) * ZLSOCP(:,:) )
!
!*      Compute Atheta
!
    ZATHETA_I(:,:)= ZAMOIST_I(:,:) * PEXN(:,:) *                        &
        ( ( ZE(:,:) - PR(:,:,1) ) * ZLSOCP(:,:) /                    &
          ( 1. + ZDEDT(:,:) * ZLSOCP(:,:) )           *              &
          (                                                             &
           ZE(:,:) * (1. + ZE(:,:)/ZEPS)                                &
                        * ( -2.*XBETAI/PT(:,:) + XGAMI ) / PT(:,:)**2   &
          +ZDEDT(:,:) * (1. + 2. * ZE(:,:)/ZEPS)                        &
                        * ( XBETAI/PT(:,:) - XGAMI ) / PT(:,:)          &
          )                                                             &
         - ZDEDT(:,:)                                                   &
        )

  ELSE
    ZAMOIST_I(:,:)=0.
    ZATHETA_I(:,:)=0.
  ENDIF

  PAMOIST(:,:) = (1.0-PFRAC_ICE(:,:))*ZAMOIST_W(:,:) &
                         +PFRAC_ICE(:,:) *ZAMOIST_I(:,:)
  PATHETA(:,:) = (1.0-PFRAC_ICE(:,:))*ZATHETA_W(:,:) &
                         +PFRAC_ICE(:,:) *ZATHETA_I(:,:)

!
ELSE
  PAMOIST(:,:) = 0.
  PATHETA(:,:) = 0.
ENDIF
END SUBROUTINE COMPUTE_FUNCTION_THERMO_MF
