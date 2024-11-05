!MNH_LIC Copyright 2018-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-------------------------------------------------------------------------------
MODULE MODE_LIMA_CONVERSION_MELTING_SNOW
  IMPLICIT NONE
CONTAINS
!     ##############################################################################
  SUBROUTINE LIMA_CONVERSION_MELTING_SNOW (LIMAP, LIMAC, LIMAM, KSIZE, ODCOMPUTE,                   &
                                           PRHODREF, PPRES, PT, PKA, PDV, PCJ, &
                                           PRVT, PRST, PCST, PLBDS,            &
                                           P_RS_CMEL, P_CS_CMEL                )
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
!  J. Wurtz       03/2022: new snow characteristics
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST,              ONLY : XTT, XMV, XMD, XLVTT, XCPV, XCL, XESTT, XRV
USE MODD_PARAM_LIMA_MIXED, ONLY:PARAM_LIMA_MIXED_t
USE MODD_PARAM_LIMA_COLD, ONLY:PARAM_LIMA_COLD_t
USE MODD_PARAM_LIMA, ONLY:PARAM_LIMA_t
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER, INTENT(IN) :: KSIZE
LOGICAL, DIMENSION(KSIZE),INTENT(IN)    :: ODCOMPUTE
!
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRHODREF ! Reference Exner function
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PPRES    !
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PT       !
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PKA      !
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PDV      !
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PCJ      !
!
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRVT    ! 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRST    ! Snow mr at t
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PCST    ! Snow C. at t
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PLBDS   ! 
!
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_RS_CMEL
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_CS_CMEL
!
!*       0.2   Declarations of local variables :
!
TYPE(PARAM_LIMA_MIXED_t),INTENT(IN)::LIMAM
TYPE(PARAM_LIMA_COLD_t),INTENT(IN)::LIMAC
TYPE(PARAM_LIMA_t),INTENT(IN)::LIMAP
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
P_CS_CMEL(:)=0.
!
ZW(:) = 0.0
WHERE( PRST(:)>LIMAP%XRTMIN(5) .AND. PCST(:)>LIMAP%XCTMIN(5) .AND. PT(:)>XTT .AND. ODCOMPUTE(:) )
   ZW(:) = PRVT(:)*PPRES(:)/((XMV/XMD)+PRVT(:)) ! Vapor pressure
   ZW(:) = PKA(:)*(XTT-PT(:)) +                                 &
              ( PDV(:)*(XLVTT + ( XCPV - XCL ) * ( PT(:) - XTT )) &
                          *(XESTT-ZW(:))/(XRV*PT(:))             )
!
! compute RSMLT
!
   ZW(:)  = LIMAM%XFSCVMG*MAX( 0.0,( -ZW(:) * PCST(:) *                        &
                               ( LIMAC%X0DEPS*PLBDS(:)**LIMAC%XEX0DEPS +             &
                                 LIMAC%X1DEPS*PCJ(:)*PLBDS(:)**LIMAC%XEX1DEPS *      &
                                   (1+0.5*(LIMAC%XFVELOS/PLBDS(:))**LIMAP%XALPHAS)**(-LIMAP%XNUS+LIMAC%XEX1DEPS/LIMAP%XALPHAS)) ))
! On ne tient pas compte de la collection de pluie et gouttelettes par la neige si T>0 !!!! 
! Note that no heat is exchanged because the graupeln produced are still icy!!!
   P_RS_CMEL(:) = - ZW(:)
   P_CS_CMEL(:) = - ZW(:) * PCST(:) / PRST(:)
!
END WHERE
!
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LIMA_CONVERSION_MELTING_SNOW
END MODULE MODE_LIMA_CONVERSION_MELTING_SNOW
