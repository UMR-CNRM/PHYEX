!MNH_LIC Copyright 2018-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-------------------------------------------------------------------------------
!      #################################
       MODULE MODI_LIMA_CONVERSION_MELTING_SNOW
!      #################################
!
INTERFACE
   SUBROUTINE LIMA_CONVERSION_MELTING_SNOW (LDCOMPUTE,                          &
                                            PRHODREF, PPRES, PT, PKA, PDV, PCJ, &
                                            PRVT, PRST, PLBDS,                  &
                                            P_RS_CMEL                           )
!
LOGICAL, DIMENSION(:),INTENT(IN)    :: LDCOMPUTE
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRHODREF ! Reference Exner function
REAL, DIMENSION(:),   INTENT(IN)    :: PPRES    !
REAL, DIMENSION(:),   INTENT(IN)    :: PT       !
REAL, DIMENSION(:),   INTENT(IN)    :: PKA      !
REAL, DIMENSION(:),   INTENT(IN)    :: PDV      !
REAL, DIMENSION(:),   INTENT(IN)    :: PCJ      !
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRVT    ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PRST    ! Cloud water C. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDS   ! 
!
REAL, DIMENSION(:),   INTENT(OUT)   :: P_RS_CMEL
!
END SUBROUTINE LIMA_CONVERSION_MELTING_SNOW
END INTERFACE
END MODULE MODI_LIMA_CONVERSION_MELTING_SNOW
!
!     ##############################################################################
      SUBROUTINE LIMA_CONVERSION_MELTING_SNOW (LDCOMPUTE,                          &
                                               PRHODREF, PPRES, PT, PKA, PDV, PCJ, &
                                               PRVT, PRST, PLBDS,                  &
                                               P_RS_CMEL                           )
!     ##############################################################################
!
!!    PURPOSE
!!    -------
!!      Compute the conversion-melting of snow into graupel
!!
!!
!!    AUTHOR
!!    ------
!!      J.-M. Cohard     * Laboratoire d'Aerologie*
!!      J.-P. Pinty      * Laboratoire d'Aerologie*
!!      S.    Berthet    * Laboratoire d'Aerologie*
!!      B.    Vié        * CNRM *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original             15/03/2018 
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST,              ONLY : XTT, XMV, XMD, XLVTT, XCPV, XCL, XESTT, XRV
USE MODD_PARAM_LIMA,       ONLY : XRTMIN
USE MODD_PARAM_LIMA_MIXED, ONLY : XFSCVMG
USE MODD_PARAM_LIMA_COLD,  ONLY : X0DEPS, XEX0DEPS, X1DEPS, XEX1DEPS
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
LOGICAL, DIMENSION(:),INTENT(IN)    :: LDCOMPUTE
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRHODREF ! Reference Exner function
REAL, DIMENSION(:),   INTENT(IN)    :: PPRES    !
REAL, DIMENSION(:),   INTENT(IN)    :: PT       !
REAL, DIMENSION(:),   INTENT(IN)    :: PKA      !
REAL, DIMENSION(:),   INTENT(IN)    :: PDV      !
REAL, DIMENSION(:),   INTENT(IN)    :: PCJ      !
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRVT    ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PRST    ! Cloud water C. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDS   ! 
!
REAL, DIMENSION(:),   INTENT(OUT)   :: P_RS_CMEL
!
!*       0.2   Declarations of local variables :
!
REAL, DIMENSION(SIZE(PRST)) :: ZW ! work arrays
!
!-------------------------------------------------------------------------------
!
!
!*       1.     Conversion-melting of snow
!	        --------------------------
!
!
P_RS_CMEL(:)=0.
!
ZW(:) = 0.0
WHERE( (PRST(:)>XRTMIN(5)) .AND. (PT(:)>XTT) .AND. LDCOMPUTE(:) )
   ZW(:) = PRVT(:)*PPRES(:)/((XMV/XMD)+PRVT(:)) ! Vapor pressure
   ZW(:) = PKA(:)*(XTT-PT(:)) +                                 &
              ( PDV(:)*(XLVTT + ( XCPV - XCL ) * ( PT(:) - XTT )) &
                          *(XESTT-ZW(:))/(XRV*PT(:))             )
!
! compute RSMLT
!
   ZW(:)  = XFSCVMG*MAX( 0.0,( -ZW(:) *             &
                          ( X0DEPS*           PLBDS(:)**XEX0DEPS +     &
                            X1DEPS*PCJ(:)*PLBDS(:)**XEX1DEPS ) ))!-    &
! On ne tient pas compte de la collection de pluie et gouttelettes par la neige si T>0 !!!! 
! Note that no heat is exchanged because the graupeln produced are still icy!!!
   P_RS_CMEL(:) = - ZW(:)
!
END WHERE
!
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LIMA_CONVERSION_MELTING_SNOW
